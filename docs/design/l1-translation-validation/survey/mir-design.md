# MIR 数据结构 + Pass 6-8 + 测试框架设计

> Stage 1 Week 5 主产出（Stage 1 收尾）
> 基于：`hir-design.md`（Week 4）+ `cpp-subset-semantics.md`（Week 3）+ Mathlib/Lean stdlib 调研
> 日期：2026-04-21

---

## §0 导读

### §0.1 目标

定义 MIR（Middle IR）数据结构 + Pass 6-8 精确规格 + back-to-back 测试框架，闭合 Stage 1 设计阶段。

### §0.2 范围

- **§1 MIR 数据结构**：CFG、SSA 形式、phi 节点
- **§2 Pass 6 `ssa_build`**：HIR₄ → MIR₀（Cytron 算法）
- **§3 Pass 7 `loop_lower`**：MIR₀ → MIR₁（循环提取 + break/continue/return 下降）
- **§4 Pass 8 `codegen`**：MIR₁ → Lean 源代码
- **§5 back-to-back 测试框架**
- **§6 Pass 管道编排**
- **§7 单元测试/回归测试基建**

### §0.3 架构衔接

```
HIR₄ (来自 Week 4)
    ↓ Pass 6: ssa_build  (this doc §2)
MIR₀: SSA + phi + CFG
    ↓ Pass 7: loop_lower  (this doc §3)
MIR₁: 循环已提取为 partial def 尾递归
    ↓ Pass 8: codegen     (this doc §4)
Lean 源代码 (.lean)
    ↓ back-to-back 测试    (this doc §5)
C++ 输出 ≡ Lean #eval 输出？
```

---

## §1 MIR 数据结构

### §1.1 与 HIR 的差异

| 特性 | HIR | MIR |
|---|---|---|
| 变量版本 | 所有 `Var.version = 0` | `Var.version ≥ 0`，SSA 形式 |
| 赋值 | `AssignStmt` / `CompoundAssignStmt` 允许 | **无**（只能 `LetStmt`）|
| 结构化循环 | `WhileStmt` / `ForStmt` / `RangeForStmt` / `DoWhileStmt` | **无**（已下降为 `TailCallStmt`）|
| Break/Continue/Return | 直接节点 | **无**（已下降为 flag + TailCall）|
| Phi 节点 | **无** | `PhiStmt` 显式 |
| CFG | 隐式（结构化）| **显式**（`MIRFunc.cfg`）|
| 辅助函数 | `aux_lambdas`（lambda lifting）| `aux_defs`（循环提取）|

### §1.2 CFG 数据结构

```python
@dataclass
class BasicBlock:
    """基本块：一段无跳转的直线代码 + 末尾 terminator。"""
    bb_id: int                      # 全局唯一 id
    stmts: list['MIRStmt']          # let、phi（开头）、require 等非跳转语句
    terminator: 'Terminator'        # 末尾跳转

@dataclass
class CFG:
    """控制流图。"""
    entry: int                      # 入口 bb_id
    blocks: dict[int, BasicBlock]   # bb_id → block
    # 前驱/后继由 terminator 隐式给出，可缓存：
    preds: dict[int, list[int]] = field(default_factory=dict)
    succs: dict[int, list[int]] = field(default_factory=dict)

@dataclass
class Terminator:
    """基本块末尾的控制流。"""
    pass

@dataclass
class JumpTerm(Terminator):
    """无条件跳转。"""
    target: int                     # bb_id

@dataclass
class CondJumpTerm(Terminator):
    """条件跳转：cond ? then_bb : else_bb"""
    cond: 'ExprIR'
    then_bb: int
    else_bb: int

@dataclass
class ReturnTerm(Terminator):
    """函数返回。"""
    value: Optional['ExprIR'] = None

@dataclass
class TailCallTerm(Terminator):
    """尾递归调用（循环下降后的出口）。
    不同于 Call：表示"跳转到另一个 def 的开头，参数对应"。"""
    target_func: str                # e.g. "loop_42_ir"
    args: list['ExprIR']

# Terminator 可以扩展，但 CLPoly 不需要 switch/case（已确认 0 处）
```

### §1.3 MIR 语句

```python
@dataclass
class PhiStmt:
    """phi 节点：x_n := phi(pred1 → x_i, pred2 → x_j, ...)"""
    target: 'Var'                   # 被定义的变量（新版本）
    ty: 'TypeIR'
    sources: dict[int, 'Var']       # {pred_bb_id: source_var_version}

@dataclass
class LetStmt:  # 同 HIR，但 Var.version > 0
    var: 'Var'
    ty: 'TypeIR'
    value: 'ExprIR'                  # ExprIR 结构不变（来自 HIR）

@dataclass
class RequireStmt:  # 同 HIR
    cond: 'ExprIR'
    name: str
    source: str
    uid: int = 0

# AssignStmt / CompoundAssignStmt / IfStmt / 循环 / Break / Continue / Return
# 全部不出现在 MIR 中

MIRStmt = Union[PhiStmt, LetStmt, RequireStmt]
```

**MIR 不变量**：
- `PhiStmt` 只出现在基本块开头（连续 N 个）
- 每个 `Var(name, version)` 全局唯一定义（SSA）
- `LetStmt` 定义的 `var.version > 0`（除入口块的参数）

### §1.4 MIRFunc

```python
@dataclass
class MIRFunc:
    """MIR 层函数。
    body 和 aux_lambdas 被转换为 cfg + aux_defs。"""
    base_name: str
    instance_suffix: str
    mangled_name: str
    params: list['HIRParam']         # 同 HIRParam
    ret_ty: 'TypeIR'
    requires: list['RequireStmt']    # 签名中的 require 参数

    cfg: 'CFG'                       # 主 CFG

    # Pass 7 产出的辅助 def（循环提取、lambda 已在 HIR 层提升）
    aux_defs: list['MIRFunc'] = field(default_factory=list)

    @property
    def lean_name(self) -> str:
        if self.instance_suffix:
            return f"{self.base_name}_{self.instance_suffix}_ir"
        return f"{self.base_name}_ir"
```

### §1.5 不变量阶梯

| 阶段 | 不变量 |
|---|---|
| **MIR₀**（Pass 6 后）| 每 `Var` SSA、phi 节点已放置、CFG 完整、无 AssignStmt、仍有结构化循环边（back edge）|
| **MIR₁**（Pass 7 后）| 无结构化循环边，所有循环提取为独立 `MIRFunc.aux_defs`、break/continue/return 已下降为 TailCallTerm + flag |

---

## §2 Pass 6: `ssa_build` (HIR₄ → MIR₀)

### §2.1 职责

把 HIR₄（带结构化控制流 + AssignStmt）转为 MIR₀（SSA + phi + CFG）。核心是 **Cytron 1991 算法**。

### §2.2 算法概览（4 阶段）

```
ssa_build(hir_func):
    # 阶段 A：构造 CFG（结构化 HIR → 基本块）
    cfg = cfg_from_hir(hir_func.body)
    
    # 阶段 B：计算支配树（Dominator Tree）
    dom_tree = compute_dominator_tree(cfg)
    
    # 阶段 C：计算支配边界（Dominance Frontier）
    df = compute_dominance_frontier(cfg, dom_tree)
    
    # 阶段 D：放置 phi 节点 + 变量重命名
    place_phi_nodes(cfg, df, hir_func)
    rename_variables(cfg, dom_tree, hir_func)
    
    return MIRFunc(..., cfg=cfg)
```

### §2.3 阶段 A：CFG 构造

**输入**：结构化 HIR body（含 IfStmt、WhileStmt、ForStmt、RangeForStmt、DoWhileStmt、Break、Continue、Return）

**算法**：递归遍历，每遇到分支/循环就新建基本块 + terminator。

```python
def cfg_from_hir(body):
    bb_counter = 0
    entry_bb = new_bb()  # bb_id = 0
    current = entry_bb
    
    def process(stmts, current_bb, loop_ctx=None):
        nonlocal bb_counter
        for stmt in stmts:
            if isinstance(stmt, IfStmt):
                then_bb, else_bb, merge_bb = new_bb(), new_bb(), new_bb()
                current_bb.terminator = CondJumpTerm(stmt.cond, then_bb.id, else_bb.id)
                process(stmt.then_body, then_bb, loop_ctx)
                then_bb.terminator = JumpTerm(merge_bb.id)
                process(stmt.else_body, else_bb, loop_ctx)
                else_bb.terminator = JumpTerm(merge_bb.id)
                current_bb = merge_bb
            elif isinstance(stmt, WhileStmt):
                header_bb = new_bb()
                body_bb = new_bb()
                exit_bb = new_bb()
                current_bb.terminator = JumpTerm(header_bb.id)
                header_bb.terminator = CondJumpTerm(stmt.cond, body_bb.id, exit_bb.id)
                process(stmt.body, body_bb, loop_ctx={"header": header_bb, "exit": exit_bb})
                body_bb.terminator = JumpTerm(header_bb.id)  # back edge
                current_bb = exit_bb
            elif isinstance(stmt, ForStmt):
                # init → header → (cond ? body : exit) → body → step → header → ...
                process(stmt.init, current_bb, loop_ctx)
                # 类似 WhileStmt，但 body 尾部跳 step 再回 header
                ...
            elif isinstance(stmt, (RangeForStmt, DoWhileStmt)):
                # 展开为 ForStmt/WhileStmt 等价形式
                ...
            elif isinstance(stmt, BreakStmt):
                assert loop_ctx, "break outside loop"
                current_bb.terminator = JumpTerm(loop_ctx["exit"].id)
                current_bb = new_bb()  # unreachable after break
            elif isinstance(stmt, ContinueStmt):
                assert loop_ctx, "continue outside loop"
                # 对 while：跳 header；对 for：跳 step；对 do-while：跳条件
                current_bb.terminator = JumpTerm(loop_ctx["continue_target"].id)
                current_bb = new_bb()
            elif isinstance(stmt, ReturnStmt):
                current_bb.terminator = ReturnTerm(stmt.value)
                current_bb = new_bb()  # unreachable after return
            else:
                current_bb.stmts.append(stmt)  # let / require / expr-stmt
        return current_bb
    
    process(body, entry_bb)
    return cfg
```

**预留 terminator**：尾块没有 terminator 时补 `ReturnTerm(None)`（void 函数末尾隐式 return）或报错（非 void 函数）。

### §2.4 阶段 B：支配树（Dominator Tree）

**定义**：节点 a **支配** 节点 b ⟺ 所有从 entry 到 b 的路径都经过 a。

**算法**（简化迭代版）：

```python
def compute_dominator_tree(cfg: CFG) -> dict[int, int]:
    """返回 idom: bb_id → 直接支配者的 bb_id"""
    # 初始化
    idom = {cfg.entry: cfg.entry}
    for bb in cfg.blocks:
        if bb != cfg.entry:
            idom[bb] = None
    
    # 按 postorder 逆序迭代
    post_order = reverse_postorder(cfg)
    changed = True
    while changed:
        changed = False
        for bb in post_order:
            if bb == cfg.entry: continue
            new_idom = None
            for pred in cfg.preds[bb]:
                if idom[pred] is not None:
                    if new_idom is None:
                        new_idom = pred
                    else:
                        new_idom = intersect(new_idom, pred, idom)
            if idom[bb] != new_idom:
                idom[bb] = new_idom
                changed = True
    return idom

def intersect(b1, b2, idom):
    """CFG 中 b1 和 b2 的最近共同祖先 (ancestor)。"""
    while b1 != b2:
        while b1 < b2: b1 = idom[b1]  # 按 postorder 编号比较
        while b2 < b1: b2 = idom[b2]
    return b1
```

**参考**：Cooper, Harvey, Kennedy "A Simple, Fast Dominance Algorithm" (2001)，Cytron 原算法的简化版。

### §2.5 阶段 C：支配边界（Dominance Frontier）

**定义**：节点 a 的支配边界 DF(a) = 满足下列条件的节点 b 集合：
- a 支配 b 的某个前驱 p
- a 不严格支配 b

**算法**：

```python
def compute_dominance_frontier(cfg: CFG, idom: dict[int, int]) -> dict[int, set[int]]:
    df = {bb: set() for bb in cfg.blocks}
    for bb in cfg.blocks:
        preds = cfg.preds[bb]
        if len(preds) >= 2:
            for p in preds:
                runner = p
                while runner != idom[bb]:
                    df[runner].add(bb)
                    runner = idom[runner]
    return df
```

### §2.6 阶段 D：Phi 放置 + 变量重命名

**Phi 放置**：对每个变量 v，在所有 "v 被赋值的块的 DF" 放置 phi 节点（迭代闭包）。

**变量重命名**：DFS 支配树，每个块入口 push 新版本，出口 pop；为 AssignStmt 生成新版本；为 phi 节点填充 sources。

```python
def place_phi_nodes(cfg, df, func):
    # 对每个变量 v，收集其 def block 集合 DefBlocks(v)
    def_blocks = collect_def_blocks(cfg)
    for v_name, blocks in def_blocks.items():
        work = set(blocks)
        already_placed = set()
        while work:
            b = work.pop()
            for frontier_b in df[b]:
                if frontier_b not in already_placed:
                    # 在 frontier_b 开头放 phi
                    ty = lookup_var_type(v_name, func)
                    cfg.blocks[frontier_b].stmts.insert(
                        0, PhiStmt(
                            target=Var(v_name, version=-1),  # 暂未版本化
                            ty=ty,
                            sources={pred: Var(v_name, version=-1)
                                     for pred in cfg.preds[frontier_b]},
                        )
                    )
                    already_placed.add(frontier_b)
                    work.add(frontier_b)

def rename_variables(cfg, idom, func):
    """DFS 支配树，维护每个变量的版本栈。"""
    version_counter = {}
    stack = defaultdict(list)  # var_name → 当前版本栈
    
    # 初始化：函数参数版本 0
    for p in func.params:
        stack[p.name].append(0)
        version_counter[p.name] = 0
    
    def rename_bb(bb_id):
        bb = cfg.blocks[bb_id]
        pushes = defaultdict(int)  # 本块 push 了几次，出口 pop
        
        # 1. 重命名 phi 节点的 target（source 稍后填）
        for stmt in bb.stmts:
            if isinstance(stmt, PhiStmt):
                new_ver = next_version(stmt.target.name, version_counter)
                stmt.target.version = new_ver
                stack[stmt.target.name].append(new_ver)
                pushes[stmt.target.name] += 1
        
        # 2. 重命名非 phi 语句
        for stmt in bb.stmts:
            if isinstance(stmt, PhiStmt): continue
            # 先重命名 RHS（读 = 用栈顶版本）
            if isinstance(stmt, LetStmt):
                rename_reads_in_expr(stmt.value, stack)
                # 再定义 LHS
                new_ver = next_version(stmt.var.name, version_counter)
                stmt.var.version = new_ver
                stack[stmt.var.name].append(new_ver)
                pushes[stmt.var.name] += 1
            elif isinstance(stmt, RequireStmt):
                rename_reads_in_expr(stmt.cond, stack)
        
        # 3. 重命名 terminator 中的变量
        rename_terminator(bb.terminator, stack)
        
        # 4. 填充所有**后继块**的 phi 节点中对应当前块的 source
        for succ_id in cfg.succs[bb_id]:
            succ = cfg.blocks[succ_id]
            for stmt in succ.stmts:
                if isinstance(stmt, PhiStmt):
                    if bb_id in stmt.sources:
                        # 用栈顶版本（即当前块结尾时该变量的最新版本）
                        if stack[stmt.target.name]:
                            stmt.sources[bb_id] = Var(stmt.target.name, stack[stmt.target.name][-1])
                        else:
                            # 变量在当前路径未被定义 — 错误
                            raise Error(...)
        
        # 5. DFS 支配树子节点
        for child_id in dom_tree_children(idom, bb_id):
            rename_bb(child_id)
        
        # 6. 退出时 pop
        for name, n in pushes.items():
            for _ in range(n):
                stack[name].pop()
    
    rename_bb(cfg.entry)
```

### §2.7 复杂度

- 阶段 A：O(|Stmt|)
- 阶段 B（Cooper-Harvey-Kennedy）：O(|BB| × α(|BB|))
- 阶段 C：O(|Edges|)
- 阶段 D：O(|BB| + |Vars| × |DF|)

对 CLPoly 最复杂的 `__lll_reduce`（~150 stmts / ~30 BB），整体亚秒完成。

### §2.8 关键测试用例

Pass 6 的单元测试应覆盖：

1. **Trivial**：单块、无分支、无赋值循环
2. **Linear 赋值**：`x = 1; x = 2;` → `x_0 = 1; x_1 = 2;`
3. **If 合并**：`if (c) x = 1 else x = 2; use x;` → phi 在合并块
4. **While 循环**：`while (c) x = x+1;` → phi 在 header 块
5. **嵌套 if/while**：复杂 phi 放置
6. **Break/Continue**：多 exit 的循环
7. **完整函数**：`__squarefree_Zp` 手动验证 SSA 结果

### §2.9 实现规模

~600 行 Python（Cytron 算法 + CFG 工具）。

### §2.10 语义保持

见 `cpp-subset-semantics.md` §3 引理 L3.1（Appel 1998 引用）。

---

## §3 Pass 7: `loop_lower` (MIR₀ → MIR₁)

### §3.1 职责

把 MIR₀ 中的**循环**（通过 back edge 识别）提取为独立的 `MIRFunc`（尾递归 partial def），break/continue/return 下降为 flag + TailCall。

### §3.2 循环识别

CFG 中的**自然循环**：back edge (A → B where B dominates A) + 所有能到达 A 且被 B 支配的节点。

```python
def find_natural_loops(cfg: CFG, idom: dict[int, int]) -> list[Loop]:
    loops = []
    for a in cfg.blocks:
        for b in cfg.succs[a]:
            if dominates(b, a, idom):  # b → a 是 back edge
                header = b
                nodes = {header, a}
                work = [a]
                while work:
                    n = work.pop()
                    for p in cfg.preds[n]:
                        if p not in nodes and dominates(header, p, idom):
                            nodes.add(p)
                            work.append(p)
                loops.append(Loop(header=header, nodes=nodes))
    return loops
```

### §3.3 循环提取算法

对每个循环 L（header, nodes, exits）：

1. **识别循环变量**：在 L 中被修改的变量集 `live_vars`
2. **生成独立 `MIRFunc` `loop_N_ir`**：参数 = live_vars 的 header 版本；body = L 的 CFG（去除 back edge，改为 tail call）
3. **在外层 MIRFunc 中**：把进入 L 的 edge 替换为 `TailCallTerm("loop_N_ir", live_vars)`；L 的出口链接到外层继续

**伪代码**：

```python
def lower_loop(cfg, loop, func):
    # 1. 分析循环：live vars
    live_vars = analyze_live_vars(loop.nodes, cfg)
    
    # 2. 构造 loop body 的 sub-CFG
    sub_cfg = extract_sub_cfg(cfg, loop.nodes)
    
    # 3. back edge → TailCall 到自己
    for bb in sub_cfg.blocks:
        if isinstance(bb.terminator, JumpTerm) and bb.terminator.target == loop.header:
            bb.terminator = TailCallTerm(f"loop_{loop.id}_ir", 
                                          args=[Var(v) for v in live_vars])
    
    # 4. 循环 exit → ReturnTerm（返回当前 live_vars 状态作为 tuple）
    for exit_bb in loop.exits:
        # 该 exit 块的最后语句应该是跳出循环 — 改为 return
        new_ret = ReturnTerm(TupleExpr([Var(v) for v in live_vars]))
        exit_bb.terminator = new_ret
    
    # 5. 构造独立 MIRFunc
    loop_func = MIRFunc(
        base_name=f"loop_{loop.id}",
        instance_suffix="",
        params=[HIRParam(v, lookup_type(v), False, False, False) for v in live_vars],
        ret_ty=TupleType([lookup_type(v) for v in live_vars]),
        cfg=sub_cfg,
    )
    func.aux_defs.append(loop_func)
    
    # 6. 外层 CFG 中：进入 loop.header 的 edge → TailCall 到 loop_func
    for pred in cfg.preds[loop.header]:
        if pred not in loop.nodes:
            # pred 是循环外部，它跳到 header 应改为调用 loop_func
            bb = cfg.blocks[pred]
            # 调用后的结果解构到 live_vars 后续使用
            call_stmt = LetStmt(
                Var("_loop_result"),
                ty=loop_func.ret_ty,
                value=Call(loop_func.lean_name, [Var(v) for v in live_vars]),
            )
            bb.stmts.append(call_stmt)
            # 再把 _loop_result 解构到 live_vars
            for i, v in enumerate(live_vars):
                bb.stmts.append(LetStmt(Var(v, version=...), lookup_type(v),
                                        FieldAccess(Var("_loop_result"), str(i))))
            # terminator 指向循环外的后继
            bb.terminator = JumpTerm(loop.post_exit_bb)
    
    # 7. 从外层 CFG 中移除循环节点
    for n in loop.nodes:
        if n != loop.header:  # header 保留作为 stub？或彻底删除
            del cfg.blocks[n]
```

### §3.4 break / continue / return 下降

在循环提取时：
- **Break**：跳到循环 exit → 外层 CFG 跳到 post-loop
- **Continue**：跳到 header → sub_cfg 内部 tail call
- **Return**：在循环内 return 必须传到外层。两种方案：
  1. **Flag 方案**（`cpp-subset-semantics.md` §2.8 已定义）：循环函数返回 `(result, _ret_flag, _ret_val)`；调用方检查 flag
  2. **Exception 方案**（`Except`）：Lean 支持；CLPoly 不用（§7.3 确认）

**选 Flag 方案**：
```lean
partial def loop_42_ir (live_vars..., _ret_flag, _ret_val) : State × Bool × RetType :=
  ...
  if (原 return 点) then (curr_state, true, ret_val)
  else if (原 break) then (curr_state, false, default)  -- 正常退出
  else loop_42_ir (updated_live_vars..., _ret_flag, _ret_val)  -- continue 或 body 尾部

-- 外层:
let (state', rf, rv) := loop_42_ir (initial_state, false, default)
if rf then rv else continue_using state'
```

### §3.5 CLPoly 的循环规模

Week 1 调研数据（`cpp-construct-catalog.md` §3）：
- IfStmt 275
- WhileStmt 20
- ForStmt 145
- CXXForRangeStmt 92
- DoStmt 2（`__wang_core`、`__zassenhaus_recombine`）
- BreakStmt 26
- ContinueStmt 52

约 **260 个循环**需要被 Pass 7 提取为独立 def。加上 HIR 层的 26 个 lambda，总计 Lean 文件会有约 **350-400 个 partial def**（67 主函数 + 260 循环 + 26 lambda）。

### §3.6 实现规模

~450 行 Python。

### §3.7 语义保持

见 `cpp-subset-semantics.md` §2 引理 L2.2（Winskel）、L2.3（break/continue flag 模式）。

---

## §4 Pass 8: `codegen` (MIR₁ → Lean 源代码)

### §4.1 职责

把 MIR₁ 翻译为 Lean 4 源代码字符串。保持 1:1 对应：每个 MIR 节点一条 Lean 语句。

### §4.2 Lean 代码结构

```lean
-- Auto-generated by cpp2lean v2
import CLPoly.Model  -- clpoly_model.lean

-- Aux definitions (lambda lifts + loop extractions)
partial def _lambda__factor_Zp_1_ir (params...) : ret := ...
partial def loop_42_ir (params...) : ret := ...

-- 主函数
partial def factor_Zp_ir (f : SparsePolyZp)
    (h_prime : f.prime.toNat.Prime)
    (h_...) : Zp × Array (SparsePolyZp × UInt64) :=
  -- CFG 线性化为 let 链 + partial def 调用
  let x_0 := ...
  let x_1 := ...
  ...
  (lc, result)
```

### §4.3 CFG 线性化

从 CFG 生成 Lean 代码的挑战：Lean 无 goto，所以**必须线性化**。

**策略**：用支配关系 + phi 节点将 CFG 转为嵌套 let 链：

```python
def emit_cfg(cfg, idom):
    """按支配关系 DFS，生成 Lean 代码字符串。"""
    def emit_bb(bb_id, rename_ctx):
        bb = cfg.blocks[bb_id]
        code = []
        
        # 1. Phi 节点：合并前驱的 source → let x_n := <from pred>
        #    但 phi 的语义是"进入此块时选择"，需要调用方预先准备
        #    实际在 Lean 中，phi 体现在调用 loop_XX_ir 的参数传递
        
        # 2. 非 phi 语句
        for stmt in bb.stmts:
            if isinstance(stmt, LetStmt):
                code.append(f"let {gen_var(stmt.var)} : {gen_ty(stmt.ty)} := {gen_expr(stmt.value)}")
            elif isinstance(stmt, RequireStmt):
                # 已提升到函数签名，body 内消除
                pass
        
        # 3. Terminator
        term = bb.terminator
        if isinstance(term, JumpTerm):
            # 被支配的后继 → 嵌套；否则 error
            if dominates(bb_id, term.target, idom):
                code.extend(emit_bb(term.target, rename_ctx))
            else:
                # 汇合点：通过 phi 处理 — 应已线性化
                raise InternalError(...)
        elif isinstance(term, CondJumpTerm):
            code.append(f"if {gen_expr(term.cond)} then")
            code.append(emit_bb(term.then_bb, rename_ctx))
            code.append(f"else")
            code.append(emit_bb(term.else_bb, rename_ctx))
        elif isinstance(term, ReturnTerm):
            code.append(f"  {gen_expr(term.value) if term.value else 'Unit.unit'}")
        elif isinstance(term, TailCallTerm):
            args = " ".join(gen_expr(a) for a in term.args)
            code.append(f"  {term.target_func} {args}")
        
        return "\n".join(code)
    
    return emit_bb(cfg.entry, {})
```

### §4.4 辅助 def 的输出顺序

Lean 要求被调函数先定义。Pass 8 按依赖拓扑序输出：

1. **aux_lambdas**（HIR 层已提升的 lambda）
2. **aux_defs**（MIR 层提取的 loop）
3. **主函数**（可以互递归的一组用 `mutual`，但 CLPoly 的主函数多数不互递归）

对调用图有环的情况（如 `factorize` ↔ `__factor_multivar`），用 Lean 的 `mutual` 块：

```lean
mutual
  partial def factorize_lex_ir ... := ... __factor_multivar_ir ...
  partial def __factor_multivar_ir ... := ... factorize_lex_ir ...
end
```

### §4.5 实现规模

~350 行 Python（48% 复用 v1 `lean_codegen.py` 的 `gen_expr`/`gen_type`）。

---

## §5 Back-to-back 测试框架

### §5.1 总体架构

```
测试用例 (JSON)
    │
    ├── C++ 侧：
    │     compile cpp2lean_test_harness.cc
    │       → 读 JSON 输入 → 调 factorize/__ddf_Zp 等 → 序列化结果 → 输出 JSON
    │
    ├── Lean 侧：
    │     生成 eval_<func>.lean：
    │       import 生成的 IR.lean
    │       #eval (serialize (func_ir test_input)) 
    │     运行 `lake env lean eval_<func>.lean` → stdout 捕获 JSON
    │
    └── Diff：
          Python 工具比较两 JSON（浮点近似 tolerance）
          Diff 工具输出：通过/失败 + 差异详情
```

### §5.2 C++ 侧 harness

```cpp
// proof/cpp2lean_v2/b2b/harness.cc
#include <clpoly/polynomial_factorize.hh>
#include <json.hpp>  // nlohmann/json

using json = nlohmann::json;

json run_test(const std::string& func_name, const json& input) {
    if (func_name == "__factor_Zp") {
        SparsePolyZp f = deserialize_poly_zp(input["f"]);
        auto [lc, factors] = __factor_Zp(f);
        return {
            {"lc", serialize_zp(lc)},
            {"factors", serialize_factor_list(factors)}
        };
    }
    // ... 为 TRANSLATION_SCOPE 的每个函数写 dispatcher ...
}

int main(int argc, char** argv) {
    json test_case = json::parse(std::cin);
    json result = run_test(test_case["func"], test_case["input"]);
    std::cout << result.dump() << std::endl;
    return 0;
}
```

**手动工作**：为每个 TRANSLATION_SCOPE 函数写 deserializer + serializer（~65 个 dispatcher）。

**自动化选项**：基于 `cpp-construct-catalog.md` 的类型表自动生成 dispatcher 代码（Python 脚本生成 .cc）。

### §5.3 Lean 侧 runner

对每个函数 `foo`，生成：

```lean
-- eval_foo.lean
import CLPoly.IR
import CLPoly.Model

def read_input : IO TestInput := ...
def serialize_output (x : FooResult) : IO Unit := ...

def main : IO Unit := do
  let input ← read_input
  let result := foo_ir input.f input.h1 input.h2
  serialize_output result
```

**执行**：`lake env lean eval_foo.lean`，stdout 捕获。

### §5.4 JSON 序列化约定

| CLPoly 类型 | JSON 格式 |
|---|---|
| `Zp { val, prime }` | `{"val": 123, "prime": 4294967291}` |
| `ZZ` | 十进制字符串（保留大数）|
| `UMonomial { deg }` | `{"deg": 5}` |
| `SparsePolyZp` | `[[<umon>, <zp>], ...]` |
| `MvPolyZZ` | `{"terms": [[<mon>, <zz>], ...], "order": "lex"}` |
| `Factorization` | `{"content": ..., "factors": [[<poly>, mult], ...]}` |

### §5.5 Diff 工具

```python
def diff_results(cpp_out: dict, lean_out: dict, func_name: str):
    """比较两个 JSON，支持：
    - 浮点 tolerance
    - 数组顺序（有些结果是集合而非有序）
    - 因子重组（Z[x] 因式分解的因子可能相差单位，如 x-1 vs -(1-x)）
    """
    # 1. 精确比较标量
    # 2. 对 Factorization：先比较 content，再对 factors 做"多集合比较"
    #    允许单位差异（除 2 个因子的乘积如一致则等价）
    ...
```

### §5.6 测试用例池

**最小起步**（10 用例）：
1. 空多项式、常数多项式、一次多项式
2. 简单可分解（`(x-1)(x-2)`）
3. 带重因子（`(x-1)^2`）
4. 模素数小 prime p=3
5. 模素数大 prime p=2^63-25
6. 多变量（2 变量、3 变量）
7. Wilkinson W(10) 单变量
8. `__ddf_Zp` 的典型输入
9. `__squarefree_Zp` 的典型输入
10. `factorize` 端到端

**扩展到 100+ 用例**：从 CLPoly 现有 `test/` 目录借用已知测试向量。

### §5.7 实现规模

- C++ harness：~500 行（dispatcher + ser/des）
- Lean runner generator：~200 行 Python
- Diff 工具：~150 行 Python
- 测试数据管理：~100 行 Python

---

## §6 Pass 管道编排

### §6.1 主入口脚本

```python
# proof/cpp2lean_v2/main.py
from pass_parse import parse_pass
from pass_ref_elim import ref_elim_pass
from pass_lambda_lift import lambda_lift_pass
from pass_iter_recognize import iter_recognize_pass
from pass_operator_resolve import operator_resolve_pass
from pass_ssa_build import ssa_build_pass
from pass_loop_lower import loop_lower_pass
from codegen import emit_lean

def translate(cpp_file: str, output_file: str):
    # Frontend
    ast_jsons = parse_ast(cpp_file)  # Clang dump
    
    hir_funcs = []
    for func_name, func_ast in discover_instances(ast_jsons):
        hir0 = parse_pass(func_ast)
        assert_hir0_invariant(hir0)
        
        hir1 = ref_elim_pass(hir0)
        assert_hir1_invariant(hir1)
        
        hir2 = lambda_lift_pass(hir1)
        assert_hir2_invariant(hir2)
        
        hir3 = iter_recognize_pass(hir2)
        assert_hir3_invariant(hir3)
        
        hir4 = operator_resolve_pass(hir3)
        assert_hir4_invariant(hir4)
        
        hir_funcs.append(hir4)
    
    mir_funcs = []
    for hir in hir_funcs:
        mir0 = ssa_build_pass(hir)
        mir1 = loop_lower_pass(mir0)
        mir_funcs.append(mir1)
    
    # 全局依赖排序
    ordered = topological_sort(mir_funcs)
    
    # 生成 Lean 代码
    lean_code = emit_lean(ordered)
    with open(output_file, 'w') as f:
        f.write(lean_code)
```

### §6.2 错误处理

- 每 Pass 遇到 unsupported 构造 → `raise TranslationError(pass, func, line, reason)`
- 主入口捕获 `TranslationError`，打印上下文 + 失败函数名 + Pass 名 + 具体原因
- 不做"部分翻译"：单个函数失败 → 整体失败（除非 --allow-partial flag）

### §6.3 日志与调试

```python
# Dump 中间 HIR/MIR 状态到 _debug/ 目录（--debug flag）
def dump_hir(func: HIRFunc, stage: str):
    """dump HIR JSON（便于调试）"""
    with open(f"_debug/{func.lean_name}.{stage}.json", 'w') as f:
        json.dump(serialize_hir(func), f, indent=2)
```

### §6.4 实现规模

~200 行 Python。

---

## §7 单元测试与回归测试基建

### §7.1 测试结构

```
proof/cpp2lean_v2/
├── tests/
│   ├── test_pass1_parse.py        (~100 行)
│   ├── test_pass2_ref_elim.py     (~80 行)
│   ├── test_pass3_lambda_lift.py  (~120 行)
│   ├── test_pass4_iter_recognize.py (~150 行)
│   ├── test_pass5_operator_resolve.py (~100 行)
│   ├── test_pass6_ssa_build.py    (~200 行，多场景)
│   ├── test_pass7_loop_lower.py   (~150 行)
│   ├── test_codegen.py            (~100 行)
│   ├── test_integration.py        (端到端，~150 行)
│   └── test_regression.py         (v5 PASS 的所有函数必须持续通过)
```

### §7.2 测试数据

- 手写的小 C++ 代码片段（每 Pass 的 unit 测试）
- CLPoly 真实函数的子集（5-10 个，从 trivial 到复杂）
- v1 的 PASS 历史（v5 的 60 PASS 作为 regression 基线）

### §7.3 CI 集成

`bash proof/cpp2lean_v2/run_tests.sh` 运行全部：

```bash
python -m pytest proof/cpp2lean_v2/tests/
python proof/cpp2lean_v2/b2b/run_b2b.py  # back-to-back
```

### §7.4 实现规模

~300 行 Python（总计跨多文件）。

---

## §8 Stage 1 总体完成情况

### §8.1 5 周产出统计

| Week | 主产出 | 行数 |
|---|---|---|
| Week 1 | `cpp-construct-catalog.md` + 10 surveys | 4000+ |
| Week 2 | `type-system.md` + 9 surveys | 5400+ |
| Week 3 | `cpp-subset-semantics.md` + refs | 800 |
| Week 4 | `hir-design.md` + v1-reuse | 1350 |
| **Week 5** | **`mir-design.md`** | **本文件，目标 ~1200 行** |
| **Stage 1 合计** | | **~13000 行** |

### §8.2 Stage 2 实施规模预估（本文档 §9 更新）

从 `v1-reuse-inventory.md` + 本文档各 Pass 估计：

| Pass | 预估代码 |
|---|---|
| Pass 1 parse | 400 |
| Pass 2 ref_elim | 200 |
| Pass 3 lambda_lift | 300 |
| Pass 4 iter_recognize | 350 |
| Pass 5 operator_resolve | 250 |
| Pass 6 ssa_build | 600 |
| Pass 7 loop_lower | 450 |
| Pass 8 codegen | 350 |
| `ir_types.py`（HIR+MIR）| 250 |
| `class_map.py`（v1 复用）| 500 |
| `clang_hybrid.py`（v1 复用）| 200 |
| Orchestration (main.py)| 200 |
| back-to-back 框架 | 950 |
| 单元测试 | 300 |
| **总计** | **~5300 行** |

比 v1 5956 行减少 ~11%（代码更精简、分工更清晰）。

---

## §9 Week 5 验收

- [x] §1 MIR 数据结构完整定义（CFG、BasicBlock、Terminator、PhiStmt、MIRFunc）
- [x] §2 Pass 6 `ssa_build` 规格（Cytron 算法 4 阶段）
- [x] §3 Pass 7 `loop_lower` 规格（循环提取 + break/continue/return flag 下降）
- [x] §4 Pass 8 `codegen` 规格（CFG 线性化 + Lean 代码生成）
- [x] §5 back-to-back 测试框架设计
- [x] §6 Pass 管道编排 + 错误处理
- [x] §7 单元测试/回归测试基建
- [x] §8 Stage 1 总体完成情况 + Stage 2 代码量预估

---

## §10 Stage 1 → Stage 2 移交

**Stage 1 设计完成！** 下一步 Stage 2 实施阶段（4-5 周）的入口：

1. 创建 `proof/cpp2lean_v2/` 目录
2. 从 v1 直接复用 `clang_hybrid.py` + `class_map.py`（~723 行）
3. 基于 `hir-design.md` §1 定义 `ir_types.py`（HIR + MIR dataclass）
4. **优先实施顺序**：
   - Week 1：基础设施 + Pass 1 parse
   - Week 2：Pass 2 ref_elim + Pass 3 lambda_lift
   - Week 3：Pass 4 iter_recognize + Pass 5 operator_resolve
   - Week 4：Pass 6 ssa_build（Cytron，最关键）
   - Week 5：Pass 7 loop_lower + Pass 8 codegen
   - Week 6-8：back-to-back 框架 + 集成测试 + 修 bug

**每个 Pass 完成后立即写单元测试**（每 Pass 100-200 行测试，使 bug 早发现）。

**成功标准**：65 函数（+ factorize 3 实例 = 67 Lean def）在 back-to-back 测试中与 C++ 输出一致，对 50+ 测试向量。

---

本文件是 **Stage 1 Week 5** 的硬性产出，也是 **Stage 1 整个 5 周设计阶段**的收官文档。
Stage 2 实施可以无歧义地按本文件的 Pass 规格 + 相关的 `hir-design.md` / `type-system.md` / `cpp-subset-semantics.md` 进行。
