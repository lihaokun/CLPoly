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
    BaseType, NamedType, UnknownType, TypeIR, RefType, StdMapType,
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
        self._range_counter = 0  # 给 __rangefor_idx_N / __rangefor_cont_N 编号

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
        # B3 修复：从 entry 做 DFS，删除 unreachable 块
        # （break/continue/return 后 new_bb 留下的孤儿）
        reachable: set[int] = set()
        stack = [cfg.entry]
        while stack:
            bb_id = stack.pop()
            if bb_id in reachable:
                continue
            reachable.add(bb_id)
            for succ in cfg.succs.get(bb_id, []):
                if succ not in reachable:
                    stack.append(succ)
        for bb_id in list(cfg.blocks.keys()):
            if bb_id not in reachable:
                del cfg.blocks[bb_id]
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
            # B8 修复：真展开为索引循环
            #   __rangefor_cont_N := <container>
            #   __rangefor_idx_N  := 0
            # header:
            #   if __rangefor_idx_N < Array.size(__rangefor_cont_N) goto body else exit
            # body:
            #   <user var> := __rangefor_cont_N[__rangefor_idx_N]
            #   <decomposition lets if any>
            #   <s.body>
            # latch:
            #   __rangefor_idx_N := __rangefor_idx_N + 1
            #   goto header
            n = self._range_counter
            self._range_counter += 1
            cont_name = f"__rangefor_cont_{n}"
            idx_name = f"__rangefor_idx_{n}"
            # 从 s.container 推 cont_var.ty（B1 修：之前一律标 UnknownType
            # 导致下游 Pass 7 cap_param + Pass 8 codegen 类型 sorry 级联）
            cont_ty: TypeIR = UnknownType("")
            if isinstance(s.container, Var) and s.container.ty is not None \
                    and not isinstance(s.container.ty, UnknownType):
                cont_ty = s.container.ty
            cont_var = Var(name=cont_name, ty=cont_ty)
            # idx 用 Nat 与 Lean Array indexing 原生类型一致（B1 续修）
            # 之前用 Int64 → Pass 8 emit `idx < Array.size cont` 导致
            # `Int64 vs Nat` 类型不匹配批量错误。
            idx_var = Var(name=idx_name, ty=BaseType.NAT)

            # pre-block：缓存容器 + 初始化 idx
            current.stmts.append(LetStmt(
                var=cont_var, ty=cont_ty, value=s.container))
            current.stmts.append(LetStmt(
                var=idx_var, ty=BaseType.NAT,
                value=Lit(value=0, ty=BaseType.NAT)))

            header = self.new_bb()
            body = self.new_bb()
            latch = self.new_bb()
            exit_bb = self.new_bb()
            current.terminator = JumpTerm(target=header.bb_id)

            # header: idx < size(cont)
            # B1 续修：StdMap 容器用 StdMap.size 而非 Array.size
            size_callee = "StdMap.size" if isinstance(cont_ty, StdMapType) else "Array.size"
            cond = BinOp(
                op="<",
                lhs=idx_var,
                rhs=Call(callee=size_callee, args=[cont_var],
                         ty=BaseType.NAT),
                ty=BaseType.BOOL,
            )
            header.terminator = CondJumpTerm(
                cond=cond, then_bb=body.bb_id, else_bb=exit_bb.bb_id)

            # body 首：用户变量绑定
            body.stmts.append(LetStmt(
                var=s.var, ty=s.var_ty,
                value=ArrayAccess(arr=cont_var, idx=idx_var, ty=s.var_ty)))
            # decomposition：[k, v] := s.var.first / s.var.second
            if s.decomposition:
                for i, dvar in enumerate(s.decomposition):
                    field = "first" if i == 0 else "second"
                    body.stmts.append(LetStmt(
                        var=dvar, ty=dvar.ty if dvar.ty is not None else UnknownType(""),
                        value=FieldAccess(obj=s.var, field_name=field,
                                          ty=dvar.ty)))

            body_ctx = {"header": header, "exit": exit_bb,
                        "continue_target": latch}
            body_end = self._process(list(s.body), body, body_ctx)
            if body_end.terminator is None:
                body_end.terminator = JumpTerm(target=latch.bb_id)

            # P0-1 修复：is_mutable_ref（auto&）时，latch 必须把 var 回写到
            # cont[idx]，否则 body 中通过 var 的写入丢失（C++ 引用语义）
            # 必须在 idx++ 之前——读 idx 仍是当前迭代值
            # B1 续修：用 Array.set! 替代 __write__（Pass 8 codegen 端正确生成）
            if s.is_mutable_ref:
                latch.stmts.append(AssignStmt(
                    target=cont_var,
                    value=Call(callee="Array.set!",
                               args=[cont_var, idx_var, s.var],
                               ty=cont_ty),
                ))
            # latch: idx ++ then jump header（idx 是 Nat，B1 续修）
            latch.stmts.append(AssignStmt(
                target=idx_var,
                value=BinOp(op="+", lhs=idx_var,
                            rhs=Lit(value=1, ty=BaseType.NAT),
                            ty=BaseType.NAT),
            ))
            latch.terminator = JumpTerm(target=header.bb_id)

            # P0-1 关键收尾：循环退出时把 cont 写回到原 container 表达式，
            # 否则 cont 内部的修改在原变量看不见（cont 仅 loop 内部局部，
            # 原 container Var 仍读旧版本）。无论 break 还是 cond 退出都要写。
            # AssignStmt 的 rename 阶段会按 target 形态自动处理：
            #   target=Var → 直接 bump 版本
            #   target=FieldAccess/ArrayAccess → 走 B7 record-update __write__ 链
            if s.is_mutable_ref and _root_var(s.container) is not None:
                exit_bb.stmts.insert(0, AssignStmt(
                    target=s.container,
                    value=cont_var,
                ))
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

def _idx_to_nat(idx: ExprIR) -> ExprIR:
    """把 idx 转为 Nat（Array.set! / Array.get! 期望 Nat）。
    若 idx 已是 Nat（Var.ty == NAT）则原样返回，否则包 Cast 到 Nat。"""
    idx_ty = getattr(idx, 'ty', None)
    if idx_ty == BaseType.NAT:
        return idx
    src_ty = idx_ty if isinstance(idx_ty, (BaseType, NamedType)) else UnknownType("")
    return Cast(expr=idx, source_ty=src_ty, target_ty=BaseType.NAT,
                cast_kind="IntegralCast")


def _build_record_update(tgt: ExprIR, new_val: ExprIR, root_name: str) -> ExprIR:
    """从 FieldAccess/ArrayAccess 链构造 record update 表达式：
      arr[i] = v       → Array.set! arr i.toNat v
      arr[i][j] = v    → Array.set! arr i.toNat (Array.set! arr[i] j.toNat v)
      obj.field = v    → _with obj field v   （Pass 8 codegen 转 `{ obj with field := v }`）
    递归到 root 时返回最外层表达式赋给 root。
    """
    if isinstance(tgt, ArrayAccess):
        arr = tgt.arr
        idx = _idx_to_nat(tgt.idx)
        ty = tgt.ty if tgt.ty is not None else UnknownType("")
        # arr 还是 ArrayAccess 嵌套（如 m[i][j]）→ 外层 set! 内嵌 inner set!
        if isinstance(arr, ArrayAccess) \
                and _root_var(arr) is not None \
                and _root_var(arr).name == root_name:
            inner_set = Call(callee="Array.set!",
                              args=[arr, idx, new_val], ty=ty)
            outer_arr = arr.arr
            outer_idx = _idx_to_nat(arr.idx)
            return Call(callee="Array.set!",
                        args=[outer_arr, outer_idx, inner_set], ty=ty)
        return Call(callee="Array.set!", args=[arr, idx, new_val], ty=ty)
    if isinstance(tgt, FieldAccess):
        # field name 用 Lit(string) 占位（之前用 Var → 被 Pass 7 误当 free var
        # → cap_param 含 'content'/'factors' 等 field 名导致 64+ Unknown 残留）
        return Call(callee="_with",
                    args=[tgt.obj,
                          Lit(value=tgt.field_name, ty=NamedType("FieldName")),
                          new_val],
                    ty=tgt.ty if tgt.ty is not None else UnknownType(""))
    # StdMap.get!(m, k) = v → StdMap.insert m k v
    if isinstance(tgt, Call) and isinstance(tgt.callee, str) \
            and tgt.callee in ("StdMap.get!", "StdMap.find!") \
            and len(tgt.args) == 2:
        m, k = tgt.args
        return Call(callee="StdMap.insert", args=[m, k, new_val],
                    ty=getattr(m, 'ty', None) or UnknownType(""))
    return Call(callee="__write__", args=[tgt, new_val], ty=UnknownType(""))


def _root_var(e: ExprIR) -> Var | None:
    """ArrayAccess/FieldAccess 链的根 Var（用于 record-update SSA）。

    P-A 修复（TD-11 阶段 A）：识别 `Call("StdMap.get!", [m, k])` 形态——
    Pass 5 把 C++ `m[k]` 解析成此 Call。语义上等价 `m[k]` 这种 lvalue，
    root = m。后续 B7 record-update 链生成 `m := __write__(StdMap.get!(m,k), v)`。
    """
    if isinstance(e, Var):
        return e
    if isinstance(e, ArrayAccess):
        return _root_var(e.arr)
    if isinstance(e, FieldAccess):
        return _root_var(e.obj)
    if isinstance(e, Call) and isinstance(e.callee, str) \
            and e.callee == "StdMap.get!" and len(e.args) >= 1:
        return _root_var(e.args[0])
    return None


def _collect_def_blocks_and_types(cfg: CFG, params: list[HIRParam]
                                  ) -> tuple[dict[str, set[int]], dict[str, TypeIR]]:
    """收集每个变量被赋值的块集合 + 类型。覆盖 LetStmt 和 AssignStmt。

    B7 关键：`vec[i] = e` 在 SSA 形态下是 `vec` 的 record-update，必须计为
    `vec` 的 def——否则 phi 放置遗漏 → 后续读拿到旧版本（silent bug）。
    """
    defs: dict[str, set[int]] = defaultdict(set)
    types: dict[str, TypeIR] = {}
    for p in params:
        types[p.name] = p.ty

    def _record_target_def(tgt, bb_id):
        if isinstance(tgt, Var):
            defs[tgt.name].add(bb_id)
        else:
            root = _root_var(tgt)
            if root is not None:
                defs[root.name].add(bb_id)

    for bb_id, bb in cfg.blocks.items():
        for s in bb.stmts:
            if isinstance(s, LetStmt):
                defs[s.var.name].add(bb_id)
                types.setdefault(s.var.name, s.ty if s.ty is not None else UnknownType(""))
            elif isinstance(s, AssignStmt):
                _record_target_def(s.target, bb_id)
            elif isinstance(s, CompoundAssignStmt):
                _record_target_def(s.target, bb_id)
    return dict(defs), types


def _collect_var_reads(e: ExprIR, out: set[str]) -> None:
    """递归收集 e 中所有 Var 名（用于 liveness use 集）。"""
    if isinstance(e, Var):
        out.add(e.name); return
    if isinstance(e, (Lit, UnresolvedOp, UnknownExpr, LambdaExpr)):
        return
    if isinstance(e, Cast):
        _collect_var_reads(e.expr, out); return
    if isinstance(e, BinOp):
        _collect_var_reads(e.lhs, out); _collect_var_reads(e.rhs, out); return
    if isinstance(e, UnaryOp):
        _collect_var_reads(e.operand, out); return
    if isinstance(e, CondExpr):
        _collect_var_reads(e.cond, out)
        _collect_var_reads(e.then_e, out)
        _collect_var_reads(e.else_e, out); return
    if isinstance(e, FieldAccess):
        _collect_var_reads(e.obj, out); return
    if isinstance(e, ArrayAccess):
        _collect_var_reads(e.arr, out); _collect_var_reads(e.idx, out); return
    if isinstance(e, Call):
        for a in e.args: _collect_var_reads(a, out)
        return
    if isinstance(e, TupleExpr):
        for x in e.elems: _collect_var_reads(x, out)
        return
    if isinstance(e, ArrayLit):
        for x in e.elems: _collect_var_reads(x, out)
        return


def _collect_block_use_def(bb: BasicBlock) -> tuple[set[str], set[str]]:
    """计算块内 use（前向：在 def 之前被读）和 def 集。"""
    use: set[str] = set()
    defined: set[str] = set()
    for s in bb.stmts:
        if isinstance(s, LetStmt):
            tmp: set[str] = set()
            if s.value is not None:
                _collect_var_reads(s.value, tmp)
            use |= (tmp - defined)
            defined.add(s.var.name)
        elif isinstance(s, AssignStmt):
            tmp = set()
            _collect_var_reads(s.value, tmp)
            if isinstance(s.target, Var):
                use |= (tmp - defined)
                defined.add(s.target.name)
            else:
                # FieldAccess/ArrayAccess target — target 本身也是 read
                _collect_var_reads(s.target, tmp)
                use |= (tmp - defined)
        elif isinstance(s, CompoundAssignStmt):
            tmp = set()
            _collect_var_reads(s.value, tmp)
            if isinstance(s.target, Var):
                tmp.add(s.target.name)  # x op= e 读 x
                use |= (tmp - defined)
                defined.add(s.target.name)
        elif isinstance(s, RequireStmt):
            tmp = set()
            _collect_var_reads(s.cond, tmp)
            use |= (tmp - defined)
        elif isinstance(s, ExprStmt):
            tmp = set()
            _collect_var_reads(s.expr, tmp)
            use |= (tmp - defined)
    # terminator
    t = bb.terminator
    if isinstance(t, ReturnTerm) and t.value is not None:
        tmp = set(); _collect_var_reads(t.value, tmp); use |= (tmp - defined)
    elif isinstance(t, CondJumpTerm):
        tmp = set(); _collect_var_reads(t.cond, tmp); use |= (tmp - defined)
    elif isinstance(t, TailCallTerm):
        for a in t.args:
            tmp = set(); _collect_var_reads(a, tmp); use |= (tmp - defined)
    return use, defined


def _compute_live_in(cfg: CFG) -> dict[int, set[str]]:
    """后向数据流：live_in = use ∪ (live_out − def)；用于 pruned SSA。

    在 phi 放置 *之前* 调用，所以不必处理 phi 的 source/target。
    """
    use_def: dict[int, tuple[set[str], set[str]]] = {
        bb_id: _collect_block_use_def(bb) for bb_id, bb in cfg.blocks.items()
    }
    live_in: dict[int, set[str]] = {bb_id: set() for bb_id in cfg.blocks}
    changed = True
    while changed:
        changed = False
        for bb_id, bb in cfg.blocks.items():
            use, defined = use_def[bb_id]
            live_out: set[str] = set()
            for succ in cfg.succs.get(bb_id, []):
                live_out |= live_in[succ]
            new_in = use | (live_out - defined)
            if new_in != live_in[bb_id]:
                live_in[bb_id] = new_in
                changed = True
    return live_in


def place_phi_nodes(cfg: CFG, df: dict[int, set[int]],
                    params: list[HIRParam]) -> None:
    """对每个变量 v，在 def_blocks(v) 的 DF 闭包处放 phi 节点。

    B2 修复：pruned SSA — 只在 v ∈ live_in(fb) 时放置 phi，避免在
    分支汇合处生成永不被读的 phi（其 sources 全为 ver=0 的悬空引用）。
    """
    def_blocks, types = _collect_def_blocks_and_types(cfg, params)
    live_in = _compute_live_in(cfg)

    for v_name, blocks in def_blocks.items():
        ty = types.get(v_name, UnknownType(""))
        already_placed: set[int] = set()
        work: list[int] = list(blocks)
        while work:
            b = work.pop()
            for fb in df.get(b, set()):
                if fb in already_placed:
                    continue
                # 仅在 v 在 fb live-in 时放 phi（pruned SSA）
                if v_name not in live_in.get(fb, set()):
                    continue
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

    # 收集所有被赋值的名字 + 参数名（用于 B7 record-update root 判定）
    def_blocks_local, _ = _collect_def_blocks_and_types(cfg, params)
    def_names: set[str] = set(def_blocks_local.keys()) | {p.name for p in params}

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
                    # B7 修复：FieldAccess/ArrayAccess target — record update：
                    # `vec[i] = e`  →  `vec_(n+1) := vec_n.set! i e`（Lean Array.set!）
                    # `obj.field = e` → `obj_(n+1) := { obj_n with field := e }`
                    # FieldAccess+ArrayAccess 嵌套（如 m[i][j]=v）→ 递归构造
                    new_tgt = rename_reads(tgt)
                    new_val = rename_reads(s.value)
                    root = _root_var(tgt)
                    if root is not None and root.name in def_names:
                        ty = root.ty if root.ty is not None else UnknownType("")
                        new_ver = fresh(root.name)
                        pushes[root.name] += 1
                        # 递归从 tgt（FieldAccess/ArrayAccess 链）构造 set/with 表达式
                        new_value_expr = _build_record_update(new_tgt, new_val, root.name)
                        new_stmts.append(LetStmt(
                            var=Var(name=root.name, version=new_ver, ty=ty),
                            ty=ty,
                            value=new_value_expr))
                    else:
                        # root 是 param 或外部引用：仍以合成 sideeff 保留
                        # （完整 fix 需要 param 视为可变 → 跨 Pass 工程，留 TODO）
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

def _collect_iter_origins(stmts: list[StmtIR]) -> dict[str, tuple[ExprIR, ExprIR]]:
    """P-B 修复（TD-11 阶段 B）：扫描 body 找 iterator 创建点，建 it→(map, key) 映射。

    识别形态：
      let it : Iterator := StdMap.find(m, k)    → it_origins[it] = (m, k)

    递归进入 IfStmt / WhileStmt / ForStmt / RangeForStmt / DoWhileStmt / BlockStmt 子 body。
    返回 dict：iter 名 → (source map expr, key expr)。

    Limitation：仅 LetStmt 直接绑定的 find 来源；alias chain（`let it2 := it1`）不跨。
    跨循环边界的 alias 不处理（保守）。
    """
    origins: dict[str, tuple[ExprIR, ExprIR]] = {}

    def walk(ss: list[StmtIR]):
        for s in ss:
            if isinstance(s, LetStmt) and isinstance(s.value, Call):
                callee = s.value.callee
                if isinstance(callee, str) and callee == "StdMap.find" \
                        and len(s.value.args) == 2:
                    origins[s.var.name] = (s.value.args[0], s.value.args[1])
            if isinstance(s, IfStmt):
                walk(s.then_body); walk(s.else_body or [])
            elif isinstance(s, (WhileStmt, RangeForStmt, DoWhileStmt)):
                walk(s.body)
            elif isinstance(s, ForStmt):
                walk(s.init); walk(s.body)
            elif isinstance(s, BlockStmt):
                walk(s.stmts)
    walk(stmts)
    return origins


def _rewrite_iter_deref_in_expr(e: ExprIR,
                                origins: dict[str, tuple[ExprIR, ExprIR]]) -> ExprIR:
    """递归把 `Call("__deref__", [Var(it)])` 改写为 `Call("StdMap.get!", [m, k])`
    （仅当 it 在 origins 中）。其余表达式 passthrough。

    `__deref__(it).second` 形态：FieldAccess(obj=Call("__deref__", [it]), field="second")
    → FieldAccess(obj=Call("StdMap.get!", [m, k]), ...)? 不，更直接：
    Pass 5 把 `it->second` 解析为 FieldAccess(__deref__(it), "second")；语义上
    等价 `m[k].second`，但 StdMap 的 .second 就是 value 本身——所以
    `__deref__(it).second = v` 应改写为 `StdMap.get!(m, k) = v`（去掉 .second）。

    我们采用：仅当上游 AssignStmt.target 是 `FieldAccess(__deref__(it), "second")`
    时把 target 整体替换为 `Call("StdMap.get!", [m, k])`。本函数只处理 read 端
    的递归改写——target 替换在 _rewrite_iter_deref_assign 直接做。
    """
    if isinstance(e, Call) and isinstance(e.callee, str) \
            and e.callee == "__deref__" and len(e.args) == 1:
        it_arg = e.args[0]
        if isinstance(it_arg, Var) and it_arg.name in origins:
            m_expr, k_expr = origins[it_arg.name]
            # `__deref__(it)` 的语义在 StdMap 视角是 (key, value) 对；外部通常是
            # `__deref__(it).first` 或 `__deref__(it).second`。如果在 expr context
            # 我们不能简单替换为 StdMap.get!（因为 StdMap.get! 只返回 value，
            # 不含 first）。所以这里仅在 .second FieldAccess 包裹时整体重写——
            # 由外层 FieldAccess 检测处理。本函数对 __deref__ 本身保持原样。
            return e
    if isinstance(e, FieldAccess):
        # __deref__(it).second 改为 StdMap.get!(m, k)
        if isinstance(e.obj, Call) and isinstance(e.obj.callee, str) \
                and e.obj.callee == "__deref__" and len(e.obj.args) == 1:
            it_arg = e.obj.args[0]
            if isinstance(it_arg, Var) and it_arg.name in origins \
                    and e.field_name == "second":
                m_expr, k_expr = origins[it_arg.name]
                return Call(callee="StdMap.get!",
                            args=[m_expr, k_expr],
                            ty=e.ty)
        return FieldAccess(
            obj=_rewrite_iter_deref_in_expr(e.obj, origins),
            field_name=e.field_name, ty=e.ty)
    if isinstance(e, BinOp):
        return BinOp(op=e.op,
                     lhs=_rewrite_iter_deref_in_expr(e.lhs, origins),
                     rhs=_rewrite_iter_deref_in_expr(e.rhs, origins),
                     ty=e.ty)
    if isinstance(e, UnaryOp):
        return UnaryOp(op=e.op,
                       operand=_rewrite_iter_deref_in_expr(e.operand, origins),
                       ty=e.ty)
    if isinstance(e, CondExpr):
        return CondExpr(cond=_rewrite_iter_deref_in_expr(e.cond, origins),
                        then_e=_rewrite_iter_deref_in_expr(e.then_e, origins),
                        else_e=_rewrite_iter_deref_in_expr(e.else_e, origins),
                        ty=e.ty)
    if isinstance(e, Cast):
        return Cast(expr=_rewrite_iter_deref_in_expr(e.expr, origins),
                    source_ty=e.source_ty, target_ty=e.target_ty,
                    cast_kind=e.cast_kind)
    if isinstance(e, ArrayAccess):
        return ArrayAccess(arr=_rewrite_iter_deref_in_expr(e.arr, origins),
                           idx=_rewrite_iter_deref_in_expr(e.idx, origins),
                           ty=e.ty)
    if isinstance(e, Call):
        return Call(callee=e.callee,
                    args=[_rewrite_iter_deref_in_expr(a, origins) for a in e.args],
                    ty=e.ty)
    if isinstance(e, TupleExpr):
        return TupleExpr(elems=[_rewrite_iter_deref_in_expr(x, origins)
                                for x in e.elems], ty=e.ty)
    if isinstance(e, ArrayLit):
        return ArrayLit(elems=[_rewrite_iter_deref_in_expr(x, origins)
                               for x in e.elems], elem_ty=e.elem_ty)
    return e


def _rewrite_iter_deref_in_stmts(stmts: list[StmtIR],
                                 origins: dict[str, tuple[ExprIR, ExprIR]]
                                 ) -> list[StmtIR]:
    """递归把 stmts 中所有 `__deref__(it).second` 改写为 `StdMap.get!(m, k)`。"""
    out: list[StmtIR] = []
    for s in stmts:
        if isinstance(s, LetStmt):
            new_v = _rewrite_iter_deref_in_expr(s.value, origins) if s.value else s.value
            out.append(LetStmt(var=s.var, ty=s.ty, value=new_v))
        elif isinstance(s, AssignStmt):
            new_target = _rewrite_iter_deref_in_expr(s.target, origins)
            new_value = _rewrite_iter_deref_in_expr(s.value, origins)
            out.append(AssignStmt(target=new_target, value=new_value))
        elif isinstance(s, CompoundAssignStmt):
            new_target = _rewrite_iter_deref_in_expr(s.target, origins)
            new_value = _rewrite_iter_deref_in_expr(s.value, origins)
            out.append(CompoundAssignStmt(target=new_target, op=s.op, value=new_value))
        elif isinstance(s, RequireStmt):
            out.append(RequireStmt(cond=_rewrite_iter_deref_in_expr(s.cond, origins),
                                   name=s.name, source=s.source))
        elif isinstance(s, ExprStmt):
            out.append(ExprStmt(expr=_rewrite_iter_deref_in_expr(s.expr, origins)))
        elif isinstance(s, ReturnStmt):
            out.append(ReturnStmt(
                value=_rewrite_iter_deref_in_expr(s.value, origins) if s.value else None))
        elif isinstance(s, IfStmt):
            out.append(IfStmt(
                cond=_rewrite_iter_deref_in_expr(s.cond, origins),
                then_body=_rewrite_iter_deref_in_stmts(s.then_body, origins),
                else_body=_rewrite_iter_deref_in_stmts(s.else_body or [], origins)))
        elif isinstance(s, WhileStmt):
            out.append(WhileStmt(
                cond=_rewrite_iter_deref_in_expr(s.cond, origins),
                body=_rewrite_iter_deref_in_stmts(s.body, origins)))
        elif isinstance(s, DoWhileStmt):
            out.append(DoWhileStmt(
                cond=_rewrite_iter_deref_in_expr(s.cond, origins),
                body=_rewrite_iter_deref_in_stmts(s.body, origins)))
        elif isinstance(s, ForStmt):
            out.append(ForStmt(
                init=_rewrite_iter_deref_in_stmts(s.init, origins),
                cond=_rewrite_iter_deref_in_expr(s.cond, origins),
                step=_rewrite_iter_deref_in_stmts(s.step, origins),
                body=_rewrite_iter_deref_in_stmts(s.body, origins)))
        elif isinstance(s, RangeForStmt):
            out.append(RangeForStmt(
                var=s.var, var_ty=s.var_ty,
                container=_rewrite_iter_deref_in_expr(s.container, origins),
                body=_rewrite_iter_deref_in_stmts(s.body, origins),
                decomposition=s.decomposition,
                is_mutable_ref=s.is_mutable_ref))
        elif isinstance(s, BlockStmt):
            out.append(BlockStmt(stmts=_rewrite_iter_deref_in_stmts(s.stmts, origins)))
        else:
            out.append(s)
    return out


def _collect_seq_iter_origins(stmts: list[StmtIR]
                              ) -> dict[str, tuple[ExprIR, str]]:
    """CF-2 阶段 2：收集 sequence iterator origin。

    识别形态：
      let it : Iterator := SparsePolyZZ.toList(c)    # begin 等价
      let it : Iterator := <method>.begin(c)
      let it : Iterator := <method>.toList(c)
    返回 dict：iter 名 → (container expr, idx_var_name)。
    idx_var_name 是为该 iter 生成的索引名（`__iter_idx_<原 it 名>`）。
    """
    origins: dict[str, tuple[ExprIR, str]] = {}

    def is_seq_iter_creation(value: ExprIR) -> ExprIR | None:
        """识别 `<X>.toList(c)` / `<X>.begin(c)` 形态（含模板参数标注），返回 c。"""
        if isinstance(value, Call) and isinstance(value.callee, str) \
                and len(value.args) == 1:
            cs = value.callee.split()[0]  # 剥 `{a0}` 等模板参数标注
            if cs.endswith(".toList") or cs.endswith(".begin"):
                return value.args[0]
        return None

    def walk(ss: list[StmtIR]):
        for s in ss:
            if isinstance(s, LetStmt) and isinstance(s.value, Call):
                # 仅识别 ty 为 Iterator 的（避免误抓普通 toList）
                ty_name = ""
                if isinstance(s.ty, NamedType):
                    ty_name = s.ty.name
                if ty_name in ("Iterator", "ConstIterator"):
                    cont = is_seq_iter_creation(s.value)
                    if cont is not None:
                        origins[s.var.name] = (cont, f"__iter_idx_{s.var.name}")
            if isinstance(s, IfStmt):
                walk(s.then_body); walk(s.else_body or [])
            elif isinstance(s, (WhileStmt, RangeForStmt, DoWhileStmt)):
                walk(s.body)
            elif isinstance(s, ForStmt):
                walk(s.init); walk(s.body)
            elif isinstance(s, BlockStmt):
                walk(s.stmts)
    walk(stmts)
    return origins


def _rewrite_seq_iter_in_expr(e: ExprIR,
                              seq_origins: dict[str, tuple[ExprIR, str]]) -> ExprIR:
    """CF-2 阶段 2：递归把 sequence iterator 操作改写为索引形态。
      `__deref__(Var(it))`            → `c[idx]`
      `__deref__(Var(it)).field`      → `c[idx].field`（FieldAccess 链）
      `Iterator.deref!(Var(it))`      → `c[idx]`
      `Iterator.advance(Var(it))`     → `idx + 1`
      `Var(it) != <toList>(c)`        → `idx < Array.size(c)`
      `Var(it) == <toList>(c)`        → `idx >= Array.size(c)`
    """
    if isinstance(e, Call):
        cs = e.callee if isinstance(e.callee, str) else ""
        # __deref__(it) / Iterator.deref!(it) → c[idx]
        if cs in ("__deref__", "Iterator.deref!") and len(e.args) == 1:
            inner = e.args[0]
            if isinstance(inner, Var) and inner.name in seq_origins:
                cont, idx_name = seq_origins[inner.name]
                return ArrayAccess(arr=cont,
                                   idx=Var(name=idx_name, ty=BaseType.INT64),
                                   ty=e.ty)
        # Iterator.advance(it) → idx + 1（注意：在 AssignStmt rhs 才有意义）
        if cs == "Iterator.advance" and len(e.args) == 1:
            inner = e.args[0]
            if isinstance(inner, Var) and inner.name in seq_origins:
                _, idx_name = seq_origins[inner.name]
                return BinOp(op="+",
                             lhs=Var(name=idx_name, ty=BaseType.INT64),
                             rhs=Lit(value=1, ty=BaseType.INT64),
                             ty=BaseType.INT64)
        # CF-2 补：`__ctor__Iterator.fromList(<X>.toList(c))` 是 begin 位置的字面值
        # → 整数 0（用于 Array.erase 单实参场景如 `r.data().erase(r.data().begin())`）
        # callee 可能带模板参数标注（如 `__ctor__Iterator.fromList {a0}`），用 startswith
        if cs.startswith("__ctor__Iterator.fromList") and len(e.args) == 1:
            inner = e.args[0]
            if isinstance(inner, Call) and isinstance(inner.callee, str):
                # 同样剥模板参数：toList 可能形态 `SparsePolyZZ.toList`
                inner_cs = inner.callee.split()[0]  # 去掉 `{a0}` 等
                if inner_cs.endswith(".toList") and len(inner.args) == 1:
                    return Lit(value=0, ty=BaseType.INT64)
        return Call(callee=e.callee,
                    args=[_rewrite_seq_iter_in_expr(a, seq_origins) for a in e.args],
                    ty=e.ty)
    if isinstance(e, FieldAccess):
        # Recurse first (handles __deref__(it).field)
        new_obj = _rewrite_seq_iter_in_expr(e.obj, seq_origins)
        return FieldAccess(obj=new_obj, field_name=e.field_name, ty=e.ty)
    if isinstance(e, BinOp):
        # 关键：`it != toList(c)` / `it == toList(c)` 比较改写
        if e.op in ("==", "!="):
            lhs = e.lhs; rhs = e.rhs
            # 让 it 在左侧（标准化）
            if isinstance(rhs, Var) and rhs.name in seq_origins:
                lhs, rhs = rhs, lhs
            if isinstance(lhs, Var) and lhs.name in seq_origins:
                cont, idx_name = seq_origins[lhs.name]
                # rhs 必须是 toList(c) / begin(c) 形态指向相同 container
                rhs_cont = None
                if isinstance(rhs, Call) and isinstance(rhs.callee, str) \
                        and len(rhs.args) == 1:
                    rcs = rhs.callee.split()[0]  # 剥模板参数
                    if rcs.endswith(".toList") or rcs.endswith(".begin") \
                            or rcs.endswith(".end"):
                        rhs_cont = rhs.args[0]
                if rhs_cont is not None and _exprs_equal(cont, rhs_cont):
                    cont_ty_local = getattr(cont, 'ty', None)
                    size_call_callee = "StdMap.size" if isinstance(cont_ty_local, StdMapType) else "Array.size"
                    size_call = Call(callee=size_call_callee, args=[cont],
                                     ty=BaseType.INT64)
                    if e.op == "!=":
                        return BinOp(op="<",
                                     lhs=Var(name=idx_name, ty=BaseType.INT64),
                                     rhs=size_call, ty=BaseType.BOOL)
                    else:
                        return BinOp(op=">=",
                                     lhs=Var(name=idx_name, ty=BaseType.INT64),
                                     rhs=size_call, ty=BaseType.BOOL)
        return BinOp(op=e.op,
                     lhs=_rewrite_seq_iter_in_expr(e.lhs, seq_origins),
                     rhs=_rewrite_seq_iter_in_expr(e.rhs, seq_origins),
                     ty=e.ty)
    if isinstance(e, UnaryOp):
        return UnaryOp(op=e.op,
                       operand=_rewrite_seq_iter_in_expr(e.operand, seq_origins),
                       ty=e.ty)
    if isinstance(e, CondExpr):
        return CondExpr(cond=_rewrite_seq_iter_in_expr(e.cond, seq_origins),
                        then_e=_rewrite_seq_iter_in_expr(e.then_e, seq_origins),
                        else_e=_rewrite_seq_iter_in_expr(e.else_e, seq_origins),
                        ty=e.ty)
    if isinstance(e, Cast):
        return Cast(expr=_rewrite_seq_iter_in_expr(e.expr, seq_origins),
                    source_ty=e.source_ty, target_ty=e.target_ty,
                    cast_kind=e.cast_kind)
    if isinstance(e, ArrayAccess):
        return ArrayAccess(arr=_rewrite_seq_iter_in_expr(e.arr, seq_origins),
                           idx=_rewrite_seq_iter_in_expr(e.idx, seq_origins),
                           ty=e.ty)
    if isinstance(e, TupleExpr):
        return TupleExpr(elems=[_rewrite_seq_iter_in_expr(x, seq_origins)
                                for x in e.elems], ty=e.ty)
    if isinstance(e, ArrayLit):
        return ArrayLit(elems=[_rewrite_seq_iter_in_expr(x, seq_origins)
                               for x in e.elems], elem_ty=e.elem_ty)
    return e


def _exprs_equal(a: ExprIR, b: ExprIR) -> bool:
    """简化表达式相等（仅 Var/FieldAccess/Cast 链）。"""
    if isinstance(a, Var) and isinstance(b, Var):
        return a.name == b.name
    if isinstance(a, FieldAccess) and isinstance(b, FieldAccess):
        return a.field_name == b.field_name and _exprs_equal(a.obj, b.obj)
    if isinstance(a, Cast) and isinstance(b, Cast):
        return _exprs_equal(a.expr, b.expr)
    return False


def _rewrite_seq_iter_in_stmts(stmts: list[StmtIR],
                               seq_origins: dict[str, tuple[ExprIR, str]]
                               ) -> list[StmtIR]:
    """CF-2 阶段 2：递归 stmts 改写。

    特别处理：
      LetStmt(var=it, ty=Iterator, value=Call("<X>.toList", [c]))
        → LetStmt(var=__iter_idx_it, ty=INT64, value=Lit(0))
      AssignStmt(target=Var(it), value=Iterator.advance(it))
        → AssignStmt(target=Var(__iter_idx_it), value=BinOp("+", idx, 1))
      AssignStmt(target=Var(it), value=Var(it2)) (alias)
        → AssignStmt(target=__iter_idx_it, value=Var(__iter_idx_it2))（若两边都是 seq iter）
    """
    out: list[StmtIR] = []
    for s in stmts:
        if isinstance(s, LetStmt) and s.var.name in seq_origins:
            # 替换为 idx := 0
            cont, idx_name = seq_origins[s.var.name]
            out.append(LetStmt(
                var=Var(name=idx_name, ty=BaseType.INT64),
                ty=BaseType.INT64,
                value=Lit(value=0, ty=BaseType.INT64),
            ))
            continue
        if isinstance(s, AssignStmt) and isinstance(s.target, Var) \
                and s.target.name in seq_origins:
            # `it = Iterator.advance(it)` 或 `it = it2` (alias)
            tgt_idx_name = seq_origins[s.target.name][1]
            new_value = _rewrite_seq_iter_in_expr(s.value, seq_origins)
            # 若 value 仍是 Var(other_it) 且 other_it 也是 seq iter（alias）
            if isinstance(new_value, Var) and new_value.name in seq_origins:
                src_idx_name = seq_origins[new_value.name][1]
                new_value = Var(name=src_idx_name, ty=BaseType.INT64)
            out.append(AssignStmt(
                target=Var(name=tgt_idx_name, ty=BaseType.INT64),
                value=new_value,
            ))
            continue
        # 其它 stmt：递归
        if isinstance(s, LetStmt):
            new_v = _rewrite_seq_iter_in_expr(s.value, seq_origins) if s.value else s.value
            out.append(LetStmt(var=s.var, ty=s.ty, value=new_v))
        elif isinstance(s, AssignStmt):
            out.append(AssignStmt(
                target=_rewrite_seq_iter_in_expr(s.target, seq_origins),
                value=_rewrite_seq_iter_in_expr(s.value, seq_origins)))
        elif isinstance(s, CompoundAssignStmt):
            out.append(CompoundAssignStmt(
                target=_rewrite_seq_iter_in_expr(s.target, seq_origins),
                op=s.op,
                value=_rewrite_seq_iter_in_expr(s.value, seq_origins)))
        elif isinstance(s, RequireStmt):
            out.append(RequireStmt(
                cond=_rewrite_seq_iter_in_expr(s.cond, seq_origins),
                name=s.name, source=s.source))
        elif isinstance(s, ExprStmt):
            out.append(ExprStmt(
                expr=_rewrite_seq_iter_in_expr(s.expr, seq_origins)))
        elif isinstance(s, ReturnStmt):
            out.append(ReturnStmt(
                value=_rewrite_seq_iter_in_expr(s.value, seq_origins)
                if s.value else None))
        elif isinstance(s, IfStmt):
            out.append(IfStmt(
                cond=_rewrite_seq_iter_in_expr(s.cond, seq_origins),
                then_body=_rewrite_seq_iter_in_stmts(s.then_body, seq_origins),
                else_body=_rewrite_seq_iter_in_stmts(s.else_body or [], seq_origins)))
        elif isinstance(s, WhileStmt):
            out.append(WhileStmt(
                cond=_rewrite_seq_iter_in_expr(s.cond, seq_origins),
                body=_rewrite_seq_iter_in_stmts(s.body, seq_origins)))
        elif isinstance(s, DoWhileStmt):
            out.append(DoWhileStmt(
                cond=_rewrite_seq_iter_in_expr(s.cond, seq_origins),
                body=_rewrite_seq_iter_in_stmts(s.body, seq_origins)))
        elif isinstance(s, ForStmt):
            out.append(ForStmt(
                init=_rewrite_seq_iter_in_stmts(s.init, seq_origins),
                cond=_rewrite_seq_iter_in_expr(s.cond, seq_origins),
                step=_rewrite_seq_iter_in_stmts(s.step, seq_origins),
                body=_rewrite_seq_iter_in_stmts(s.body, seq_origins)))
        elif isinstance(s, RangeForStmt):
            out.append(RangeForStmt(
                var=s.var, var_ty=s.var_ty,
                container=_rewrite_seq_iter_in_expr(s.container, seq_origins),
                body=_rewrite_seq_iter_in_stmts(s.body, seq_origins),
                decomposition=s.decomposition,
                is_mutable_ref=s.is_mutable_ref))
        elif isinstance(s, BlockStmt):
            out.append(BlockStmt(
                stmts=_rewrite_seq_iter_in_stmts(s.stmts, seq_origins)))
        else:
            out.append(s)
    return out


def ssa_build_pass(func: HIRFunc) -> MIRFunc:
    """HIR₄ → MIR₀。"""
    # P-B 修复（TD-11 阶段 B）：在 SSA build 之前，把 `__deref__(it).second`
    # 改写为 `StdMap.get!(m, k)`，让 _root_var 能识别 root = m。
    iter_origins = _collect_iter_origins(list(func.body))
    body = _rewrite_iter_deref_in_stmts(list(func.body), iter_origins) \
        if iter_origins else list(func.body)
    # CF-2 阶段 2：sequence iterator (begin/toList) → 索引化（idx : Int64）
    seq_origins = _collect_seq_iter_origins(body)
    if seq_origins:
        body = _rewrite_seq_iter_in_stmts(body, seq_origins)
    cfg = cfg_from_hir(body, func.ret_ty)
    idom = compute_dominator_tree(cfg)
    df = compute_dominance_frontier(cfg, idom)
    place_phi_nodes(cfg, df, list(func.params))
    cfg.rebuild_edges()  # phi 不影响 edges，但保险
    rename_variables(cfg, idom, list(func.params))

    # B6 修复：aux_lambdas 必须递归 SSA-build，否则 Pass 8 codegen 会丢这些 def
    aux_defs = [ssa_build_pass(aux) for aux in func.aux_lambdas]

    return MIRFunc(
        base_name=func.base_name,
        instance_suffix=func.instance_suffix,
        mangled_name=func.mangled_name,
        qual_type=func.qual_type,
        params=list(func.params),
        ret_ty=func.ret_ty,
        requires=list(func.requires),
        cfg=cfg,
        aux_defs=aux_defs,
    )
