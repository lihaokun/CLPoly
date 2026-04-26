"""
Pass 6: ssa_build — HIR₄ → MIR₀

按 mir-design.md §2 的 Cytron 1991 算法，4 阶段：
  A. CFG 构造（结构化控制流 → 基本块 + Terminator）
  B. 支配树（Cooper-Harvey-Kennedy 简化迭代算法，2001）
  C. 支配边界 DF
  D. phi 放置 + 变量重命名（Cytron rename）

输入：HIRFunc（HIR₄，Pass 5 输出，含 IfStmt / WhileStmt / ForStmt /
      RangeForStmt / DoWhileStmt / Break / Continue / Return / AssignStmt 等
      结构化形态）
输出：MIRFunc + cfg：MIR₀（SSA + phi + CFG）

MIR₀ 不变量见 `assert_mir0_invariant`。
"""

from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR, RefType,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRFunc, HIRParam, TranslationError,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG, MIRFunc, assert_mir0_invariant,
)


# ============================================================
# §A. CFG 构造
# ============================================================

class CFGBuilder:
    """从结构化 HIR body 构造 CFG。"""

    def __init__(self, ret_ty: TypeIR):
        self.bb_counter = 0
        self.blocks: dict[int, BasicBlock] = {}
        self.ret_ty = ret_ty

    def new_bb(self) -> BasicBlock:
        bb_id = self.bb_counter
        self.bb_counter += 1
        bb = BasicBlock(bb_id=bb_id, stmts=[], terminator=None)  # type: ignore
        self.blocks[bb_id] = bb
        return bb

    def build(self, body: list[StmtIR]) -> CFG:
        entry = self.new_bb()
        end_bb = self._process(body, entry, loop_ctx=None)
        # 末尾自动补 terminator（void 函数）或报错
        if end_bb.terminator is None:
            if self.ret_ty == BaseType.UNIT:
                end_bb.terminator = ReturnTerm(value=None)
            else:
                end_bb.terminator = ReturnTerm(value=None)  # 宽松：用户应显式 return
        cfg = CFG(entry=entry.bb_id, blocks=self.blocks)
        cfg.rebuild_edges()
        return cfg

    def _process(self, stmts: list[StmtIR], current: BasicBlock,
                 loop_ctx: dict | None) -> BasicBlock:
        """递归处理 stmts，返回处理完后的"当前块"。"""
        for s in stmts:
            current = self._process_one(s, current, loop_ctx)
        return current

    def _process_one(self, s: StmtIR, current: BasicBlock,
                     loop_ctx: dict | None) -> BasicBlock:
        if current.terminator is not None:
            # 已经有 terminator（如前面 break/return 后），新建一个 unreachable 块
            current = self.new_bb()

        if isinstance(s, IfStmt):
            then_bb = self.new_bb()
            else_bb = self.new_bb()
            merge_bb = self.new_bb()
            current.terminator = CondJumpTerm(
                cond=s.cond, then_bb=then_bb.bb_id, else_bb=else_bb.bb_id)
            then_end = self._process(list(s.then_body), then_bb, loop_ctx)
            if then_end.terminator is None:
                then_end.terminator = JumpTerm(target=merge_bb.bb_id)
            else_end = self._process(list(s.else_body), else_bb, loop_ctx)
            if else_end.terminator is None:
                else_end.terminator = JumpTerm(target=merge_bb.bb_id)
            return merge_bb

        if isinstance(s, WhileStmt):
            header = self.new_bb()
            body = self.new_bb()
            exit_bb = self.new_bb()
            current.terminator = JumpTerm(target=header.bb_id)
            header.terminator = CondJumpTerm(
                cond=s.cond, then_bb=body.bb_id, else_bb=exit_bb.bb_id)
            body_ctx = {"header": header, "exit": exit_bb,
                        "continue_target": header}
            body_end = self._process(list(s.body), body, body_ctx)
            if body_end.terminator is None:
                body_end.terminator = JumpTerm(target=header.bb_id)  # back edge
            return exit_bb

        if isinstance(s, DoWhileStmt):
            # do { body } while (cond) — body 至少执行一次
            body = self.new_bb()
            check = self.new_bb()
            exit_bb = self.new_bb()
            current.terminator = JumpTerm(target=body.bb_id)
            body_ctx = {"header": check, "exit": exit_bb,
                        "continue_target": check}
            body_end = self._process(list(s.body), body, body_ctx)
            if body_end.terminator is None:
                body_end.terminator = JumpTerm(target=check.bb_id)
            cond = s.cond if s.cond is not None else Lit(True, ty=BaseType.BOOL)
            check.terminator = CondJumpTerm(
                cond=cond, then_bb=body.bb_id, else_bb=exit_bb.bb_id)
            return exit_bb

        if isinstance(s, ForStmt):
            # init → header → (cond ? body : exit); body → step → header
            current = self._process(list(s.init), current, loop_ctx)
            if current.terminator is not None:
                current = self.new_bb()
            header = self.new_bb()
            body = self.new_bb()
            step_bb = self.new_bb()
            exit_bb = self.new_bb()
            current.terminator = JumpTerm(target=header.bb_id)
            cond = s.cond if s.cond is not None else Lit(True, ty=BaseType.BOOL)
            header.terminator = CondJumpTerm(
                cond=cond, then_bb=body.bb_id, else_bb=exit_bb.bb_id)
            body_ctx = {"header": header, "exit": exit_bb,
                        "continue_target": step_bb}
            body_end = self._process(list(s.body), body, body_ctx)
            if body_end.terminator is None:
                body_end.terminator = JumpTerm(target=step_bb.bb_id)
            step_end = self._process(list(s.step), step_bb, loop_ctx)
            if step_end.terminator is None:
                step_end.terminator = JumpTerm(target=header.bb_id)
            return exit_bb

        if isinstance(s, RangeForStmt):
            # for (var : container) body — Pass 4 已 desugar decomposition；
            # 这里展开为：对容器索引循环。简化：当作 ForStmt(int i=0; i<size; ++i)
            # 的等价形式。语义上让 _i 索引访问容器（Pass 7 后续可优化）。
            # 简化建模：把 RangeFor body 包成 while-style header（无显式 step）
            header = self.new_bb()
            body = self.new_bb()
            exit_bb = self.new_bb()
            current.terminator = JumpTerm(target=header.bb_id)
            # cond：iterator 是否未到末尾——这里用 Lit(True) 占位，Pass 7
            # 详细处理（loop_lower 阶段重写）
            header.terminator = CondJumpTerm(
                cond=Lit(True, ty=BaseType.BOOL),
                then_bb=body.bb_id, else_bb=exit_bb.bb_id)
            body_ctx = {"header": header, "exit": exit_bb,
                        "continue_target": header}
            body_end = self._process(list(s.body), body, body_ctx)
            if body_end.terminator is None:
                body_end.terminator = JumpTerm(target=header.bb_id)
            return exit_bb

        if isinstance(s, BreakStmt):
            if loop_ctx is None:
                raise TranslationError(
                    pass_name="ssa_build", func_name="?",
                    reason="break outside loop")
            current.terminator = JumpTerm(target=loop_ctx["exit"].bb_id)
            return self.new_bb()  # unreachable

        if isinstance(s, ContinueStmt):
            if loop_ctx is None:
                raise TranslationError(
                    pass_name="ssa_build", func_name="?",
                    reason="continue outside loop")
            tgt = loop_ctx["continue_target"]
            current.terminator = JumpTerm(target=tgt.bb_id)
            return self.new_bb()

        if isinstance(s, ReturnStmt):
            current.terminator = ReturnTerm(value=s.value)
            return self.new_bb()

        if isinstance(s, BlockStmt):
            return self._process(list(s.stmts), current, loop_ctx)

        # let / require / assign / compound-assign / expr-stmt → 直接放当前块
        # AssignStmt 在 MIR 不允许，但 CFG 阶段保留（rename 阶段把
        # `x := expr` 转成 `let x_n := expr`）
        current.stmts.append(s)
        return current


def cfg_from_hir(body: list[StmtIR], ret_ty: TypeIR) -> CFG:
    return CFGBuilder(ret_ty).build(body)


# ============================================================
# §B. 支配树（Cooper-Harvey-Kennedy 简化迭代算法）
# ============================================================

def _reverse_postorder(cfg: CFG) -> list[int]:
    """从 entry DFS，返回逆后序（reverse postorder）—— DomTree 迭代用。"""
    visited = set()
    post: list[int] = []

    def dfs(bb_id: int):
        if bb_id in visited: return
        visited.add(bb_id)
        for succ in cfg.succs.get(bb_id, []):
            dfs(succ)
        post.append(bb_id)

    dfs(cfg.entry)
    post.reverse()
    return post


def compute_dominator_tree(cfg: CFG) -> dict[int, int]:
    """返回 idom: bb_id → 直接支配者的 bb_id（entry 自支配）。"""
    rpo = _reverse_postorder(cfg)
    rpo_index = {bb: i for i, bb in enumerate(rpo)}
    idom: dict[int, int | None] = {bb: None for bb in cfg.blocks}
    idom[cfg.entry] = cfg.entry

    def intersect(b1: int, b2: int) -> int:
        finger1, finger2 = b1, b2
        while finger1 != finger2:
            while rpo_index[finger1] > rpo_index[finger2]:
                finger1 = idom[finger1]  # type: ignore
            while rpo_index[finger2] > rpo_index[finger1]:
                finger2 = idom[finger2]  # type: ignore
        return finger1

    changed = True
    while changed:
        changed = False
        for bb in rpo:
            if bb == cfg.entry: continue
            preds = [p for p in cfg.preds.get(bb, []) if idom[p] is not None]
            if not preds: continue
            new_idom: int = preds[0]
            for p in preds[1:]:
                if idom[p] is not None:
                    new_idom = intersect(new_idom, p)
            if idom[bb] != new_idom:
                idom[bb] = new_idom
                changed = True

    # 转 dict[int, int]，未到达的块（dead code）保留 None 但断言通常不该有
    return {bb: (i if i is not None else cfg.entry) for bb, i in idom.items()}


def dom_tree_children(idom: dict[int, int]) -> dict[int, list[int]]:
    """从 idom 反向构建子节点表。"""
    children = defaultdict(list)
    for bb, i in idom.items():
        if bb != i:  # 排除 entry 自支配
            children[i].append(bb)
    return dict(children)


# ============================================================
# §C. 支配边界 DF
# ============================================================

def compute_dominance_frontier(cfg: CFG, idom: dict[int, int]
                               ) -> dict[int, set[int]]:
    df: dict[int, set[int]] = {bb: set() for bb in cfg.blocks}
    for bb in cfg.blocks:
        preds = cfg.preds.get(bb, [])
        if len(preds) >= 2:
            for p in preds:
                runner = p
                while runner != idom[bb] and runner in idom:
                    df[runner].add(bb)
                    runner = idom[runner]
                    if runner == idom[runner]:  # 到达 entry 自支配
                        if runner != idom[bb]:
                            df[runner].add(bb)
                        break
    return df


# ============================================================
# §D. Phi 放置 + 变量重命名
# ============================================================

def _collect_def_blocks_and_types(cfg: CFG, params: list[HIRParam]
                                  ) -> tuple[dict[str, set[int]], dict[str, TypeIR]]:
    """收集每个变量被赋值的块集合 + 类型。覆盖 LetStmt 和 AssignStmt。"""
    defs: dict[str, set[int]] = defaultdict(set)
    types: dict[str, TypeIR] = {}
    for p in params:
        types[p.name] = p.ty
    for bb_id, bb in cfg.blocks.items():
        for s in bb.stmts:
            if isinstance(s, LetStmt):
                defs[s.var.name].add(bb_id)
                types.setdefault(s.var.name, s.ty if s.ty is not None else UnknownType(""))
            elif isinstance(s, AssignStmt):
                tgt = s.target
                if isinstance(tgt, Var):
                    defs[tgt.name].add(bb_id)
                # FieldAccess / ArrayAccess target 不进入 SSA（视为内存写入）
            elif isinstance(s, CompoundAssignStmt):
                tgt = s.target
                if isinstance(tgt, Var):
                    defs[tgt.name].add(bb_id)
    return dict(defs), types


def place_phi_nodes(cfg: CFG, df: dict[int, set[int]],
                    params: list[HIRParam]) -> None:
    """对每个变量 v，在 def_blocks(v) 的 DF 闭包处放 phi 节点。"""
    def_blocks, types = _collect_def_blocks_and_types(cfg, params)

    for v_name, blocks in def_blocks.items():
        ty = types.get(v_name, UnknownType(""))
        already_placed: set[int] = set()
        work: list[int] = list(blocks)
        while work:
            b = work.pop()
            for fb in df.get(b, set()):
                if fb not in already_placed:
                    # 在 fb 开头插入 phi
                    phi = PhiStmt(
                        target=Var(name=v_name, version=-1, ty=ty),
                        ty=ty,
                        sources={pred: Var(name=v_name, version=-1, ty=ty)
                                 for pred in cfg.preds.get(fb, [])},
                    )
                    cfg.blocks[fb].stmts.insert(0, phi)
                    already_placed.add(fb)
                    if fb not in def_blocks.get(v_name, set()):
                        work.append(fb)


def rename_variables(cfg: CFG, idom: dict[int, int],
                     params: list[HIRParam]) -> None:
    """DFS 支配树，重命名所有 Var（LHS 新版本，RHS 用栈顶版本）。"""
    version_counter: dict[str, int] = {}
    stack: dict[str, list[int]] = defaultdict(list)
    children = dom_tree_children(idom)

    # 参数初始化为 version 0
    for p in params:
        version_counter[p.name] = 0
        stack[p.name].append(0)

    def fresh(name: str) -> int:
        version_counter[name] = version_counter.get(name, 0) + 1
        v = version_counter[name]
        stack[name].append(v)
        return v

    def cur_ver(name: str) -> int:
        if stack[name]:
            return stack[name][-1]
        return 0  # 未定义—回退到 0（参数或全局；实际上应报错）

    def rename_reads(e: ExprIR) -> ExprIR:
        if isinstance(e, Var):
            return Var(name=e.name, version=cur_ver(e.name), ty=e.ty)
        if isinstance(e, Cast):
            return Cast(expr=rename_reads(e.expr),
                        source_ty=e.source_ty, target_ty=e.target_ty,
                        cast_kind=e.cast_kind)
        if isinstance(e, BinOp):
            return BinOp(op=e.op, lhs=rename_reads(e.lhs),
                         rhs=rename_reads(e.rhs), ty=e.ty)
        if isinstance(e, UnaryOp):
            return UnaryOp(op=e.op, operand=rename_reads(e.operand), ty=e.ty)
        if isinstance(e, CondExpr):
            return CondExpr(cond=rename_reads(e.cond),
                            then_e=rename_reads(e.then_e),
                            else_e=rename_reads(e.else_e), ty=e.ty)
        if isinstance(e, FieldAccess):
            return FieldAccess(obj=rename_reads(e.obj),
                               field_name=e.field_name, ty=e.ty)
        if isinstance(e, ArrayAccess):
            return ArrayAccess(arr=rename_reads(e.arr),
                               idx=rename_reads(e.idx), ty=e.ty)
        if isinstance(e, Call):
            return Call(callee=e.callee,
                        args=[rename_reads(a) for a in e.args], ty=e.ty)
        if isinstance(e, TupleExpr):
            return TupleExpr(elems=[rename_reads(x) for x in e.elems], ty=e.ty)
        if isinstance(e, ArrayLit):
            return ArrayLit(elems=[rename_reads(x) for x in e.elems],
                            elem_ty=e.elem_ty)
        if isinstance(e, LambdaExpr):
            # lambda body 不在主 CFG 范围内（已 lifted by Pass 3）；保留原样
            return e
        return e  # Lit / UnresolvedOp / UnknownExpr

    def rename_terminator(t):
        if isinstance(t, JumpTerm) or isinstance(t, ReturnTerm) and t.value is None:
            return t
        if isinstance(t, ReturnTerm):
            return ReturnTerm(value=rename_reads(t.value))
        if isinstance(t, CondJumpTerm):
            return CondJumpTerm(cond=rename_reads(t.cond),
                                then_bb=t.then_bb, else_bb=t.else_bb)
        if isinstance(t, TailCallTerm):
            return TailCallTerm(target_func=t.target_func,
                                args=[rename_reads(a) for a in t.args])
        return t

    def rename_bb(bb_id: int):
        bb = cfg.blocks[bb_id]
        pushes: dict[str, int] = defaultdict(int)

        # 1. 重命名 phi 的 target（source 稍后填）
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                v = fresh(s.target.name)
                s.target.version = v
                pushes[s.target.name] += 1

        # 2. 处理非 phi 语句，把 AssignStmt → LetStmt 同时重命名
        new_stmts: list[StmtIR] = []
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                new_stmts.append(s)
                continue
            if isinstance(s, LetStmt):
                new_val = rename_reads(s.value) if s.value else s.value
                new_ver = fresh(s.var.name)
                pushes[s.var.name] += 1
                new_stmts.append(LetStmt(
                    var=Var(name=s.var.name, version=new_ver, ty=s.ty),
                    ty=s.ty, value=new_val))
            elif isinstance(s, AssignStmt):
                tgt = s.target
                if isinstance(tgt, Var):
                    new_val = rename_reads(s.value)
                    new_ver = fresh(tgt.name)
                    pushes[tgt.name] += 1
                    # AssignStmt → LetStmt（SSA 形态）
                    ty = tgt.ty if tgt.ty is not None else UnknownType("")
                    new_stmts.append(LetStmt(
                        var=Var(name=tgt.name, version=new_ver, ty=ty),
                        ty=ty, value=new_val))
                else:
                    # FieldAccess / ArrayAccess target — 视为副作用 stmt
                    # 包成 ExprStmt-like（这里实际不允许 ExprStmt 在 MIR；
                    # 需要 Pass 7 处理。简化：保留为 LetStmt with synthetic var name）
                    new_tgt = rename_reads(tgt)
                    new_val = rename_reads(s.value)
                    # 用合成名表示"匿名副作用"
                    syn = f"__sideeff_{bb_id}_{len(new_stmts)}"
                    fresh(syn)
                    pushes[syn] += 1
                    new_stmts.append(LetStmt(
                        var=Var(name=syn, version=cur_ver(syn),
                                ty=UnknownType("")),
                        ty=UnknownType(""),
                        value=Call(callee="__write__",
                                   args=[new_tgt, new_val],
                                   ty=BaseType.UNIT)))
            elif isinstance(s, CompoundAssignStmt):
                # MIR 不允许；Pass 5 应已展开。但保险起见：x op= e → x_n := x_(n-1) op e
                tgt = s.target
                if isinstance(tgt, Var):
                    rhs = BinOp(op=s.op,
                                lhs=Var(name=tgt.name, version=cur_ver(tgt.name)),
                                rhs=rename_reads(s.value), ty=tgt.ty)
                    new_ver = fresh(tgt.name)
                    pushes[tgt.name] += 1
                    ty = tgt.ty if tgt.ty is not None else UnknownType("")
                    new_stmts.append(LetStmt(
                        var=Var(name=tgt.name, version=new_ver, ty=ty),
                        ty=ty, value=rhs))
                # FieldAccess/ArrayAccess target 不处理（罕见）
            elif isinstance(s, RequireStmt):
                new_stmts.append(RequireStmt(
                    cond=rename_reads(s.cond), name=s.name, source=s.source))
            elif isinstance(s, ExprStmt):
                # MIR 不允许 ExprStmt；视为 side-effect Let
                syn = f"__sideeff_{bb_id}_{len(new_stmts)}"
                fresh(syn)
                pushes[syn] += 1
                new_stmts.append(LetStmt(
                    var=Var(name=syn, version=cur_ver(syn), ty=UnknownType("")),
                    ty=UnknownType(""), value=rename_reads(s.expr)))
            else:
                # 未处理类型（保险跳过；assert_mir0 会捕获）
                new_stmts.append(s)
        bb.stmts = new_stmts

        # 3. 重命名 terminator
        bb.terminator = rename_terminator(bb.terminator)

        # 4. 填后继块 phi 的 sources（按当前块的栈顶版本）
        for succ_id in cfg.succs.get(bb_id, []):
            succ = cfg.blocks[succ_id]
            for s in succ.stmts:
                if isinstance(s, PhiStmt) and bb_id in s.sources:
                    s.sources[bb_id] = Var(
                        name=s.target.name, version=cur_ver(s.target.name),
                        ty=s.ty)

        # 5. DFS 子节点
        for child in children.get(bb_id, []):
            rename_bb(child)

        # 6. 退栈
        for name, n in pushes.items():
            for _ in range(n):
                if stack[name]:
                    stack[name].pop()

    rename_bb(cfg.entry)


# ============================================================
# Pass 入口
# ============================================================

def ssa_build_pass(func: HIRFunc) -> MIRFunc:
    """HIR₄ → MIR₀。"""
    cfg = cfg_from_hir(list(func.body), func.ret_ty)
    idom = compute_dominator_tree(cfg)
    df = compute_dominance_frontier(cfg, idom)
    place_phi_nodes(cfg, df, list(func.params))
    cfg.rebuild_edges()  # phi 不影响 edges，但保险
    rename_variables(cfg, idom, list(func.params))

    return MIRFunc(
        base_name=func.base_name,
        instance_suffix=func.instance_suffix,
        mangled_name=func.mangled_name,
        qual_type=func.qual_type,
        params=list(func.params),
        ret_ty=func.ret_ty,
        requires=list(func.requires),
        cfg=cfg,
        aux_defs=[],  # Pass 7 填充
    )
