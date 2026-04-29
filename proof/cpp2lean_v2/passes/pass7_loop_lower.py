"""
Pass 7: loop_lower — MIR₀ → MIR₁

把 SSA CFG 中的 natural loops 提取为 partial def aux MIRFunc，主 cfg 变 DAG。

输入：MIRFunc，cfg 含 back edges + phi at loop headers
输出：MIRFunc，cfg 是 DAG，每个 loop 提取为 aux_defs 内独立 MIRFunc，
      用 TailCallTerm 实现尾递归（Lean 端对应 partial def）。

算法（按 loop-extraction-design.md §2 + mir-design.md MIR₁ 不变量）：

1. DomTree 找 back edges：edge (e, h) 中 h dominates e
2. 按"内层先"顺序处理（reverse postorder + nesting depth）
3. 每个 natural loop（header h + back edge sources）：
   a) 收集 loop body BBs：h ∪ {p | p reaches some back-edge source via h-dominated paths}
   b) loop func params = phi targets at h（loop-carried SSA vars）
      + free vars（loop 内读但定义在外）
   c) loop func body = copy of loop BBs，with：
      - h 块的 phi 移除（phi targets 已成 params）
      - back edge `JumpTerm/CondJumpTerm → h` → TailCallTerm(loop_name, [phi sources at h from e])
      - exit edge `→ exit_bb (not in loop)` → ReturnTerm(TupleExpr(phi targets at h))
   d) 主 cfg：h 块替换为"call block"（LetStmt(__loop_ret := loop_N(args)) +
      destructure stmts），跳到 loop exit 后块；删除 loop body BBs
   e) 重建 cfg edges
4. 重复直到无 back edge

与 lifted lambda 的关系：Pass 3b 已用 aux_defs 装 lifted lambda；Pass 7 复用此机制
装 loop functions。Lean codegen 端两者都是 partial def。

简化（当前 corpus）：
- 98% loops 单 back edge — 主路径
- 多 back edge / multi-exit loops（5 个）：合并多 exit 为 single ret tuple
- 嵌套：reverse postorder 自然处理（内层先）
"""

from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, ExprIR,
    LetStmt, RequireStmt, StmtIR,
    HIRParam, MIRStmt,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG, MIRFunc, TranslationError,
    assert_mir1_invariant,
)


# ============================================================
# §A. DomTree + Back Edge 检测（复用 Pass 6 的 idom 算法）
# ============================================================

def _reverse_postorder(cfg: CFG) -> list[int]:
    """从 entry DFS，返回逆后序。"""
    visited: set[int] = set()
    post: list[int] = []

    def dfs(bb_id: int):
        if bb_id in visited:
            return
        visited.add(bb_id)
        for s in cfg.succs.get(bb_id, []):
            dfs(s)
        post.append(bb_id)

    dfs(cfg.entry)
    return list(reversed(post))


def _compute_idom(cfg: CFG) -> dict[int, int]:
    """Cooper-Harvey-Kennedy 简化迭代 DomTree。"""
    rpo = _reverse_postorder(cfg)
    rpo_index = {bb: i for i, bb in enumerate(rpo)}
    idom: dict[int, int | None] = {bb: None for bb in cfg.blocks}
    idom[cfg.entry] = cfg.entry

    def intersect(b1: int, b2: int) -> int:
        f1, f2 = b1, b2
        while f1 != f2:
            while rpo_index[f1] > rpo_index[f2]:
                f1 = idom[f1]  # type: ignore
            while rpo_index[f2] > rpo_index[f1]:
                f2 = idom[f2]  # type: ignore
        return f1

    changed = True
    while changed:
        changed = False
        for bb in rpo:
            if bb == cfg.entry:
                continue
            preds = [p for p in cfg.preds.get(bb, []) if idom[p] is not None]
            if not preds:
                continue
            new_idom = preds[0]
            for p in preds[1:]:
                new_idom = intersect(p, new_idom)
            if idom[bb] != new_idom:
                idom[bb] = new_idom
                changed = True
    # type: ignore[return-value]
    return {bb: d for bb, d in idom.items() if d is not None}


def _dominates(a: int, b: int, idom: dict[int, int]) -> bool:
    """a 是否 dominate b（含 a == b）。"""
    cur = b
    while True:
        if cur == a:
            return True
        nxt = idom.get(cur, cur)
        if nxt == cur:
            return False
        cur = nxt


def _find_back_edges(cfg: CFG, idom: dict[int, int]) -> list[tuple[int, int]]:
    """back edges：edge (e, h) 中 h dominate e。"""
    out: list[tuple[int, int]] = []
    for e in cfg.blocks:
        for s in cfg.succs.get(e, []):
            if _dominates(s, e, idom):
                out.append((e, s))
    return out


# ============================================================
# §B. Natural Loop 收集
# ============================================================

def _natural_loop_body(cfg: CFG, header: int, back_edge_sources: list[int],
                       idom: dict[int, int]) -> set[int]:
    """Natural loop = header ∪ {所有反向可达 back_edge_sources 且 dominated by header
    的块}。"""
    body: set[int] = {header}
    work = list(back_edge_sources)
    while work:
        bb = work.pop()
        if bb in body:
            continue
        if not _dominates(header, bb, idom):
            continue
        body.add(bb)
        for p in cfg.preds.get(bb, []):
            if p not in body:
                work.append(p)
    return body


def _find_innermost_loop(cfg: CFG, idom: dict[int, int],
                          back_edges: list[tuple[int, int]]
                          ) -> tuple[int, list[int], set[int]] | None:
    """选择"内层先"的 loop：
       - 按 header 分组 back edges
       - 选 loop body 不含其它 loop header 的（即最内层）
    返回 (header, back_edge_sources, body_bbs)；无 loop 返回 None。
    """
    # 按 header 分组
    by_header: dict[int, list[int]] = defaultdict(list)
    for e, h in back_edges:
        by_header[h].append(e)

    # 计算每个 loop 的 body
    loop_bodies: dict[int, set[int]] = {}
    for h, es in by_header.items():
        loop_bodies[h] = _natural_loop_body(cfg, h, es, idom)

    # 选内层：body 不含其它 loop header（除自身）
    for h, body in loop_bodies.items():
        is_innermost = True
        for h2 in by_header:
            if h2 != h and h2 in body:
                is_innermost = False
                break
        if is_innermost:
            return h, by_header[h], body
    # 不应到这（至少有一个最内层 loop）；保险返回任意
    h = next(iter(by_header))
    return h, by_header[h], loop_bodies[h]


# ============================================================
# §C. Loop 提取
# ============================================================

def _collect_phi_targets_at(bb: BasicBlock) -> list[tuple[Var, TypeIR]]:
    """收集块的所有 phi targets（loop-carried vars）。"""
    out: list[tuple[Var, TypeIR]] = []
    for s in bb.stmts:
        if isinstance(s, PhiStmt):
            out.append((s.target, s.ty))
    return out


def _collect_phi_sources_from(bb: BasicBlock, src_bb: int) -> list[Var]:
    """收集块的 phi sources，按 phi 顺序，从 src_bb 边过来的 source。"""
    out: list[Var] = []
    for s in bb.stmts:
        if isinstance(s, PhiStmt):
            src = s.sources.get(src_bb)
            if src is None:
                # 不应发生（phi.sources 已与 preds 一致，invariant 保证）
                src = Var(name=s.target.name, version=0, ty=s.ty)
            out.append(src)
    return out


def _collect_var_reads_in_expr(e: ExprIR, out: set[str]) -> None:
    """递归收集 expr 中的 Var.name。"""
    if isinstance(e, Var):
        out.add(e.name); return
    if isinstance(e, (Lit, UnresolvedOp)):
        return
    if isinstance(e, Cast):
        _collect_var_reads_in_expr(e.expr, out); return
    if isinstance(e, BinOp):
        _collect_var_reads_in_expr(e.lhs, out)
        _collect_var_reads_in_expr(e.rhs, out); return
    if isinstance(e, UnaryOp):
        _collect_var_reads_in_expr(e.operand, out); return
    if isinstance(e, CondExpr):
        _collect_var_reads_in_expr(e.cond, out)
        _collect_var_reads_in_expr(e.then_e, out)
        _collect_var_reads_in_expr(e.else_e, out); return
    if isinstance(e, FieldAccess):
        _collect_var_reads_in_expr(e.obj, out); return
    if isinstance(e, ArrayAccess):
        _collect_var_reads_in_expr(e.arr, out)
        _collect_var_reads_in_expr(e.idx, out); return
    if isinstance(e, Call):
        for a in e.args:
            _collect_var_reads_in_expr(a, out)
        return
    from ir_types import TupleExpr, ArrayLit, BlockExpr, LambdaExpr
    if isinstance(e, TupleExpr):
        for x in e.elems: _collect_var_reads_in_expr(x, out)
        return
    if isinstance(e, ArrayLit):
        for x in e.elems: _collect_var_reads_in_expr(x, out)
        return


def _collect_var_reads_in_expr_versioned(e: ExprIR, out: set[tuple[str, int]],
                                          tys: dict[tuple[str, int], TypeIR]) -> None:
    """收集 expr 中所有 Var 引用，含版本号 + 类型。"""
    if isinstance(e, Var):
        out.add((e.name, e.version))
        if (e.name, e.version) not in tys and e.ty is not None:
            tys[(e.name, e.version)] = e.ty
        return
    if isinstance(e, (Lit, UnresolvedOp)):
        return
    if isinstance(e, Cast):
        _collect_var_reads_in_expr_versioned(e.expr, out, tys); return
    if isinstance(e, BinOp):
        _collect_var_reads_in_expr_versioned(e.lhs, out, tys)
        _collect_var_reads_in_expr_versioned(e.rhs, out, tys); return
    if isinstance(e, UnaryOp):
        _collect_var_reads_in_expr_versioned(e.operand, out, tys); return
    if isinstance(e, CondExpr):
        _collect_var_reads_in_expr_versioned(e.cond, out, tys)
        _collect_var_reads_in_expr_versioned(e.then_e, out, tys)
        _collect_var_reads_in_expr_versioned(e.else_e, out, tys); return
    if isinstance(e, FieldAccess):
        _collect_var_reads_in_expr_versioned(e.obj, out, tys); return
    if isinstance(e, ArrayAccess):
        _collect_var_reads_in_expr_versioned(e.arr, out, tys)
        _collect_var_reads_in_expr_versioned(e.idx, out, tys); return
    if isinstance(e, Call):
        for a in e.args:
            _collect_var_reads_in_expr_versioned(a, out, tys)
        return
    from ir_types import TupleExpr, ArrayLit
    if isinstance(e, TupleExpr):
        for x in e.elems: _collect_var_reads_in_expr_versioned(x, out, tys)
        return
    if isinstance(e, ArrayLit):
        for x in e.elems: _collect_var_reads_in_expr_versioned(x, out, tys)
        return


def _collect_loop_free_vars(cfg: CFG, body_bbs: set[int],
                             phi_target_at_header: list[tuple[Var, TypeIR]]
                             ) -> list[Var]:
    """收集 loop body 内 free Var（带版本）—— body 内 read 但 def 在外（也不是
    header phi target）。返回带版本的 Var 列表。
    """
    # defs in loop = phi target names + body 内 LetStmt vars + body 内 phi target
    defs_in_loop: set[str] = {v.name for v, _ in phi_target_at_header}
    reads_in_loop: set[tuple[str, int]] = set()
    read_tys: dict[tuple[str, int], TypeIR] = {}
    for bb_id in body_bbs:
        bb = cfg.blocks[bb_id]
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                defs_in_loop.add(s.target.name)
                # phi sources from outside loop = init args (not free vars in body sense)
                # phi sources from inside loop = loop-defined
                # 这里跳过 phi.sources 的收集（init args 是 caller 的事）
            elif isinstance(s, LetStmt):
                if s.value is not None:
                    _collect_var_reads_in_expr_versioned(s.value, reads_in_loop, read_tys)
                defs_in_loop.add(s.var.name)
            elif isinstance(s, RequireStmt):
                _collect_var_reads_in_expr_versioned(s.cond, reads_in_loop, read_tys)
        # terminator
        t = bb.terminator
        if isinstance(t, ReturnTerm) and t.value is not None:
            _collect_var_reads_in_expr_versioned(t.value, reads_in_loop, read_tys)
        elif isinstance(t, CondJumpTerm):
            _collect_var_reads_in_expr_versioned(t.cond, reads_in_loop, read_tys)
        elif isinstance(t, TailCallTerm):
            for a in t.args:
                _collect_var_reads_in_expr_versioned(a, reads_in_loop, read_tys)

    # 过滤出 free vars: name 不在 defs_in_loop（含 phi target name）
    free: list[Var] = []
    for (name, version) in sorted(reads_in_loop):
        if name in defs_in_loop:
            continue
        ty = read_tys.get((name, version), UnknownType(""))
        free.append(Var(name=name, version=version, ty=ty))
    return free


def _var_param_name(v: Var) -> str:
    """SSA Var → HIRParam name（含版本编码）。"""
    return f"{v.name}_{v.version}" if v.version > 0 else v.name


def _is_exit_edge(succ_bb: int, body_bbs: set[int]) -> bool:
    """succ 不在 loop body → exit edge。"""
    return succ_bb not in body_bbs


def _replace_terminator_for_back_edge(t, header: int, loop_name: str,
                                       header_block: BasicBlock,
                                       free_vars: list[str]) -> object:
    """back edge → TailCallTerm（递归）。
    args 顺序：phi sources at header from src_bb + free vars。
    src_bb 由 caller 传（这里需要重新 design — 改为 caller 处理）。
    """
    raise NotImplementedError("inline below")


def _build_loop_func(cfg: CFG, header: int, body_bbs: set[int],
                     idom: dict[int, int], loop_id: int
                     ) -> tuple[MIRFunc, list[str], list[int], list[Var]]:
    """从 cfg 提取 loop 为独立 MIRFunc。

    返回 (loop_func, phi_target_names, exit_targets, free_vars_versioned)。
    free_vars_versioned 用于 caller 在 splice 时传 init args（含版本）。
    """
    loop_name = f"_loop_{loop_id}"
    header_block = cfg.blocks[header]

    # 1. params: phi targets at header + free vars（含版本号编入 name）
    phi_targets: list[tuple[Var, TypeIR]] = _collect_phi_targets_at(header_block)
    phi_target_names = [t[0].name for t in phi_targets]
    free_vars: list[Var] = _collect_loop_free_vars(cfg, body_bbs, phi_targets)
    # cap_params: 每个 free var 编码 version 进 name（如 'k_1'）
    cap_params = [HIRParam(name=_var_param_name(v), ty=v.ty or UnknownType(""),
                            is_ref=False, is_const_ref=False, is_output=False)
                  for v in free_vars]
    # phi_params: phi.target 同样编码 version（如 'i_2'）
    phi_params = [HIRParam(name=_var_param_name(t[0]), ty=t[1], is_ref=False,
                            is_const_ref=False, is_output=False)
                  for t in phi_targets]
    new_params = phi_params + cap_params

    # 2. 收集 exit targets（按 first-encountered 顺序分配 kind）
    exit_targets: list[int] = []
    for bb_id in sorted(body_bbs):  # deterministic order
        for s in cfg.succs.get(bb_id, []):
            if s not in body_bbs and s not in exit_targets:
                exit_targets.append(s)
    n_exits = len(exit_targets)
    exit_kind: dict[int, int] = {t: i for i, t in enumerate(exit_targets)}

    # 3. ret_ty: TupleType(exit_kind, phi target tys)
    ret_tys: list[TypeIR] = [BaseType.INT64]  # exit_kind first
    ret_tys.extend(t[1] for t in phi_targets)
    if len(ret_tys) == 1:
        ret_ty: TypeIR = ret_tys[0]
    elif len(ret_tys) == 2:
        ret_ty = PairType(ret_tys[0], ret_tys[1])
    else:
        ret_ty = TupleType(tuple(ret_tys))

    # 4. body cfg = copy of body BBs，改写 terminator
    new_blocks: dict[int, BasicBlock] = {}
    for bb_id in body_bbs:
        bb = cfg.blocks[bb_id]
        new_stmts: list[MIRStmt] = []
        for s in bb.stmts:
            # header 块的 phi 全部移除（target 已成 params）
            if bb_id == header and isinstance(s, PhiStmt):
                continue
            new_stmts.append(s)
        new_term = _rewrite_loop_terminator(
            bb, body_bbs, header, header_block, loop_name, free_vars,
            phi_target_names, exit_kind, n_exits)
        new_blocks[bb_id] = BasicBlock(bb_id=bb_id, stmts=new_stmts,
                                        terminator=new_term)

    # CondJumpTerm 含 back/exit 边的块拆分（生成 dummy dispatch BBs）
    fresh_bb_id = [max(cfg.blocks.keys(), default=0) + 100 + loop_id * 1000]
    _split_cond_to_exit_or_back(new_blocks, body_bbs, header, header_block,
                                  loop_name, free_vars, phi_target_names,
                                  exit_kind, n_exits, fresh_bb_id)

    new_cfg = CFG(entry=header, blocks=new_blocks)
    new_cfg.rebuild_edges()

    loop_func = MIRFunc(
        base_name=loop_name,
        instance_suffix="",
        mangled_name="",
        qual_type=f"loop in <outer> | id={loop_id}",
        params=new_params,
        ret_ty=ret_ty,
        requires=[],
        cfg=new_cfg,
        aux_defs=[],
    )
    return loop_func, phi_target_names, exit_targets, free_vars


def _rewrite_loop_terminator(bb: BasicBlock, body_bbs: set[int],
                              header: int, header_block: BasicBlock,
                              loop_name: str, free_vars: list[Var],
                              phi_target_names: list[str],
                              exit_kind: dict[int, int],
                              n_exits: int):
    """改写 loop body 内某 block 的 terminator：
    - JumpTerm(t):
        t == header → TailCallTerm
        t in body   → keep
        t in exit   → ReturnTerm((kind=exit_kind[t], phi_targets...))
    - CondJumpTerm(then, else):
        若 then/else 不全在 body：拆出 dispatch BBs（**简化：调用方负责**，本函数
        处理"两边都同类"或"一边在 body 一边非 body"的简单 case）
    """
    src_bb_id = bb.bb_id
    t = bb.terminator
    if isinstance(t, JumpTerm):
        if t.target == header:
            return _make_tail_call(header_block, src_bb_id, loop_name, free_vars)
        if t.target not in body_bbs:
            return _make_exit_return(t.target, exit_kind, header_block,
                                      src_bb_id, phi_target_names, n_exits)
        return t  # in-loop jump unchanged
    if isinstance(t, CondJumpTerm):
        # 处理：两边任意 mix（back/exit/inner）→ 用辅助 BB 拆分
        # 但当前函数无法新增 BB（_build_loop_func 持有 dict）。
        # 简化策略：CondJumpTerm 只能两边都 in-loop（保持），
        # 或两边都 exit/back（上层负责）
        # 实际 Pass 6 RangeFor 真展开：header 块 CondJumpTerm(idx<size, body, exit)
        # → 两边 inner+exit 混合（body in body_bbs, exit not in body_bbs）
        # 这里我们需要拆分。返回一个 "needs split" 标记？
        # 简化做法：直接构造 new bbs in caller; here we return original term, caller
        # post-processes.
        return t
    if isinstance(t, ReturnTerm):
        # body 内有 return - 这等价于函数级 return。loop func 内的 return 应该
        # 是 outer function 的 return；用特殊 exit_kind=-1 表示"forward return"，
        # caller 在 dispatch 时把这个 case 对应到 outer ReturnTerm。
        # 简化：暂不支持 loop body 内的 ReturnTerm（罕见，非 break/continue 标准 case）
        return t
    return t


def _make_exit_return(target: int, exit_kind: dict[int, int],
                       header_block: BasicBlock, src_bb_id: int,
                       phi_target_names: list[str], n_exits: int) -> ReturnTerm:
    """exit edge → ReturnTerm((kind=k, phi_target_0, phi_target_1, ...))"""
    from ir_types import TupleExpr
    k = exit_kind[target]
    elems: list[ExprIR] = [Lit(value=k, ty=BaseType.INT64)]
    # phi target latest values when exiting through this edge
    # 用 src_bb_id 的视角：phi.sources at header from src_bb_id（如果存在）
    phi_targets_at_h: list[tuple[Var, TypeIR]] = _collect_phi_targets_at(header_block)
    for v, ty in phi_targets_at_h:
        src_var = None
        for s in header_block.stmts:
            if isinstance(s, PhiStmt) and s.target.name == v.name:
                src_var = s.sources.get(src_bb_id)
                break
        if src_var is None:
            # 直接用 phi target var（可能 src 不是 phi pred — 这种 case 罕见，
            # 用最新 SSA 版本作 fallback）
            src_var = v
        elems.append(src_var)
    if len(elems) == 1:
        return ReturnTerm(value=elems[0])
    return ReturnTerm(value=TupleExpr(elems=elems, ty=None))


def _make_tail_call(header_block: BasicBlock, src_bb_id: int, loop_name: str,
                    free_vars: list[Var]) -> TailCallTerm:
    """tail call args = phi sources at h from src_bb (back edge) + 带版本的 free vars。"""
    phi_srcs = _collect_phi_sources_from(header_block, src_bb_id)
    args: list[ExprIR] = list(phi_srcs)
    args.extend(free_vars)  # 带版本的 Var；body 内 free var 引用版本固定
    return TailCallTerm(target_func=loop_name, args=args)


def _split_cond_to_exit_or_back(loop_blocks: dict[int, BasicBlock],
                                 body_bbs: set[int], header: int,
                                 header_block: BasicBlock,
                                 loop_name: str, free_vars: list[Var],
                                 phi_target_names: list[str],
                                 exit_kind: dict[int, int], n_exits: int,
                                 fresh_bb_id: list[int]) -> None:
    """对 loop body 内含 CondJumpTerm 且至少一边是 back/exit 的块，拆分：
    把 back/exit 边替换为指向新 dummy BB（含 TailCallTerm/ReturnTerm）。
    """
    new_bb_to_add: dict[int, BasicBlock] = {}
    for bb_id, bb in list(loop_blocks.items()):
        t = bb.terminator
        if not isinstance(t, CondJumpTerm):
            continue
        new_then = t.then_bb
        new_else = t.else_bb
        for side in ("then", "else"):
            tgt = t.then_bb if side == "then" else t.else_bb
            if tgt in body_bbs and tgt != header:
                continue  # in-body, keep
            # back edge or exit edge — 创建新 dummy BB
            new_id = fresh_bb_id[0]
            fresh_bb_id[0] += 1
            if tgt == header:
                term = _make_tail_call(header_block, bb_id, loop_name, free_vars)
            else:
                term = _make_exit_return(tgt, exit_kind, header_block, bb_id,
                                          phi_target_names, n_exits)
            new_bb_to_add[new_id] = BasicBlock(
                bb_id=new_id, stmts=[], terminator=term)
            if side == "then":
                new_then = new_id
            else:
                new_else = new_id
        if new_then != t.then_bb or new_else != t.else_bb:
            loop_blocks[bb_id] = BasicBlock(
                bb_id=bb_id, stmts=bb.stmts,
                terminator=CondJumpTerm(cond=t.cond, then_bb=new_then,
                                          else_bb=new_else),
            )
    loop_blocks.update(new_bb_to_add)


# ============================================================
# §D. Pass 入口
# ============================================================

def loop_lower_pass(func: MIRFunc) -> MIRFunc:
    """MIR₀ → MIR₁。递归处理 aux_defs；对自身 cfg 提取所有 loops。"""
    # 递归先处理已有 aux_defs（lifted lambda）
    new_aux: list[MIRFunc] = [loop_lower_pass(a) for a in func.aux_defs]

    if func.cfg is None:
        return replace(func, aux_defs=new_aux)

    cfg = func.cfg
    loop_id_counter = 0

    while True:
        idom = _compute_idom(cfg)
        back_edges = _find_back_edges(cfg, idom)
        if not back_edges:
            break

        loop_info = _find_innermost_loop(cfg, idom, back_edges)
        if loop_info is None:
            break
        header, back_srcs, body_bbs = loop_info

        # 提取 loop 函数
        loop_func, phi_target_names, exit_targets, free_vars = _build_loop_func(
            cfg, header, body_bbs, idom, loop_id_counter)
        loop_id_counter += 1
        new_aux.append(loop_func)

        # 主 cfg 中替换 loop：header 改为 call block + dispatch BBs
        cfg = _splice_loop_call(cfg, header, body_bbs, loop_func,
                                  phi_target_names, exit_targets, free_vars)

    return replace(func, cfg=cfg, aux_defs=new_aux)


def _splice_loop_call(cfg: CFG, header: int, body_bbs: set[int],
                       loop_func: MIRFunc,
                       phi_target_names: list[str],
                       exit_targets: list[int],
                       free_vars: list[Var]) -> CFG:
    """主 cfg 中：header 块替换为 call block；删 loop body 块；分发到对应 exit。

    Loop func 返回 `(exit_kind: Int64, phi_target_0, ...)`。
    Caller:
       header block:
         let __loop_ret_N := loop_N(init_args)
         let __exit_kind := __loop_ret_N.elem0  (或 fst if PairType)
         let phi_v0 := __loop_ret_N.elem1  ...
         goto dispatch_0  # 或单 exit 时直接 goto exit
       dispatch_k: if __exit_kind == k then exit_k else dispatch_{k+1}
    """
    header_block = cfg.blocks[header]
    phi_targets = _collect_phi_targets_at(header_block)

    # init args
    outer_preds = [p for p in cfg.preds.get(header, []) if p not in body_bbs]
    if not outer_preds:
        raise TranslationError("loop_lower", loop_func.base_name,
                                "loop header has no outer pred")
    init_outer_pred = outer_preds[0]
    init_phi_srcs = _collect_phi_sources_from(header_block, init_outer_pred)
    n_phi = len(phi_targets)
    init_args: list[ExprIR] = list(init_phi_srcs)
    init_args.extend(free_vars)  # 带版本的 free var Var 直接传

    n_exits = len(exit_targets)
    if n_exits == 0:
        raise TranslationError("loop_lower", loop_func.base_name,
                                "loop has no exit edges")

    # tmp var receiving loop func result
    tmp_name = f"__loop_ret_{loop_func.base_name.split('_')[-1]}"
    tmp_var = Var(name=tmp_name, version=1, ty=loop_func.ret_ty)

    # 字段命名：n_total = 1 + n_phi（kind + phi_targets）
    n_total = 1 + n_phi
    if n_total == 1:
        # 没 phi target，仅 kind，ret_ty = INT64
        kind_field = None  # 直接用 tmp_var 即 kind
        phi_fields: list[str] = []
    elif n_total == 2:
        kind_field = "fst"
        phi_fields = ["snd"]
    else:
        kind_field = "elem0"
        phi_fields = [f"elem{k+1}" for k in range(n_phi)]

    new_header_stmts: list[MIRStmt] = []
    new_header_stmts.append(LetStmt(
        var=tmp_var, ty=loop_func.ret_ty,
        value=Call(callee=loop_func.base_name, args=init_args,
                    ty=loop_func.ret_ty),
    ))
    # destructure phi targets
    for k, (v, ty) in enumerate(phi_targets):
        new_header_stmts.append(LetStmt(
            var=v, ty=ty,
            value=FieldAccess(obj=tmp_var, field_name=phi_fields[k], ty=ty)
                  if phi_fields else tmp_var,
        ))

    # 决定 exit dispatch
    new_blocks: dict[int, BasicBlock] = {}
    fresh_bb_id = max(cfg.blocks.keys(), default=0) + 1

    if n_exits == 1:
        # 单 exit：直接 jump
        new_header_term = JumpTerm(target=exit_targets[0])
    else:
        # 多 exit：构造 dispatch chain
        # let kind := tmp_var.<kind_field>
        kind_var = Var(name=f"{tmp_name}__kind", version=1, ty=BaseType.INT64)
        if kind_field is None:
            new_header_stmts.append(LetStmt(
                var=kind_var, ty=BaseType.INT64, value=tmp_var,
            ))
        else:
            new_header_stmts.append(LetStmt(
                var=kind_var, ty=BaseType.INT64,
                value=FieldAccess(obj=tmp_var, field_name=kind_field,
                                   ty=BaseType.INT64),
            ))
        # dispatch chain：dispatch_0 → dispatch_1 → ... → dispatch_{n-2}
        # 每个 dispatch_k：if kind == k then exit_k else dispatch_{k+1}
        dispatch_ids = [fresh_bb_id + i for i in range(n_exits - 1)]
        fresh_bb_id += n_exits - 1
        new_header_term = JumpTerm(target=dispatch_ids[0])
        for k, did in enumerate(dispatch_ids):
            cond = BinOp(op="==",
                          lhs=kind_var,
                          rhs=Lit(value=k, ty=BaseType.INT64),
                          ty=BaseType.BOOL)
            then_bb = exit_targets[k]
            # else: 下一个 dispatch，最后 dispatch 的 else = exit_targets[-1]
            else_bb = dispatch_ids[k + 1] if k + 1 < len(dispatch_ids) \
                else exit_targets[-1]
            new_blocks[did] = BasicBlock(
                bb_id=did, stmts=[],
                terminator=CondJumpTerm(cond=cond, then_bb=then_bb,
                                          else_bb=else_bb),
            )

    new_header = BasicBlock(
        bb_id=header, stmts=new_header_stmts, terminator=new_header_term)
    new_blocks[header] = new_header

    # exit_pred_map：进入每个 exit_target 的新 caller-side pred BB
    # （用于重写 outer phi 的 source key）
    exit_pred_map: dict[int, int] = {}
    if n_exits == 1:
        exit_pred_map[exit_targets[0]] = header
    else:
        # dispatch_ids 已在上面创建（new_blocks 中已加）
        # 重建 dispatch_ids 列表以便引用
        # （因为 dispatch_ids 在 dispatch chain 创建处是局部变量；这里复读 new_blocks）
        # 简便：找 new_blocks 中刚加入的 dispatch BBs（按 bb_id 顺序）
        # 实际我们之前生成 dispatch_ids = [fresh_bb_id + i for i in 0..n_exits-2]
        # fresh_bb_id 已 += n_exits-1 后丢失。改用 reverse-engineer：
        first_dispatch = max(cfg.blocks.keys(), default=0) + 1
        dispatch_ids_local = [first_dispatch + i for i in range(n_exits - 1)]
        for k, et in enumerate(exit_targets):
            if k < n_exits - 1:
                exit_pred_map[et] = dispatch_ids_local[k]
            else:
                exit_pred_map[et] = dispatch_ids_local[-1]

    # 保留所有非 loop body 块；删 loop body 块（除 header 已替换）
    # 同时重写其中 phi.sources：把 body_bb 来源合并为 exit_pred_map[本 bb]
    for bb_id, bb in cfg.blocks.items():
        if bb_id == header:
            continue  # 已加
        if bb_id in body_bbs:
            continue  # 删
        # 检查该 BB 是否有 phi 引用 body_bb 作 source
        new_stmts: list[MIRStmt] = []
        any_changed = False
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                new_sources: dict[int, Var] = {}
                merged_body_var: Var | None = None
                for src_bb, src_var in s.sources.items():
                    if src_bb in body_bbs:
                        merged_body_var = src_var  # 取最后的（多 path 应同 var）
                        any_changed = True
                    else:
                        new_sources[src_bb] = src_var
                if merged_body_var is not None:
                    # body 来源合并为单 source from exit_pred
                    new_pred = exit_pred_map.get(bb_id, header)
                    new_sources[new_pred] = merged_body_var
                new_stmts.append(PhiStmt(target=s.target, ty=s.ty,
                                          sources=new_sources))
            else:
                new_stmts.append(s)
        if any_changed:
            new_blocks[bb_id] = BasicBlock(
                bb_id=bb_id, stmts=new_stmts, terminator=bb.terminator)
        else:
            new_blocks[bb_id] = bb

    new_cfg = CFG(entry=cfg.entry, blocks=new_blocks)
    new_cfg.rebuild_edges()
    return new_cfg
