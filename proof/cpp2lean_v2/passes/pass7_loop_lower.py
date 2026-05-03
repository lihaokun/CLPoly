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
    ArrayType, OptionType, StdMapType, RefType,
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
    """收集块的 phi sources，按 phi 顺序，从 src_bb 边过来的 source。

    阶段 F #1 修复：phi target 在 loop body 内部新生（如 find-loop 的
    `term_1`、EDF 的 `g_8`），outer pred 没 source，旧版返回
    `Var(target.name, v0)` 导致下游 emit 时找不到名字 → Unknown identifier。
    改为返回 Lean `default` 关键字，由 Inhabited 实例提供占位初值。
    """
    out: list[Var] = []
    for s in bb.stmts:
        if isinstance(s, PhiStmt):
            src = s.sources.get(src_bb)
            if src is None:
                src = Var(name="default", version=0, ty=s.ty)
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
        # 优先非-Unknown ty（同一 SSA Var 多处引用时，先记的可能是 Unknown 占位，
        # 后遇到的可能是正确类型 — 后者覆盖前者；与 setdefault 反过来）
        if e.ty is not None and not isinstance(e.ty, UnknownType):
            tys[(e.name, e.version)] = e.ty
        elif (e.name, e.version) not in tys and e.ty is not None:
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
        # 阶段 G-A：callee 可能是 Var（local var lambda 调用），递归收集
        # 让 cap_param 自动覆盖。其他形态（str / UnresolvedOp）跳过。
        if isinstance(e.callee, Var):
            _collect_var_reads_in_expr_versioned(e.callee, out, tys)
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
                             phi_target_at_header: list[tuple[Var, TypeIR]],
                             extra_aux_for_ty: list[MIRFunc] | None = None,
                             func_param_tys: dict[str, TypeIR] | None = None
                             ) -> list[Var]:
    """收集 loop body 内 free Var（带版本）—— body 内 read 但 def 在外（也不是
    header phi target）。返回带版本的 Var 列表。
    """
    # defs in loop = phi target names + body 内 LetStmt vars + body 内 phi target
    defs_in_loop: set[str] = {v.name for v, _ in phi_target_at_header}
    reads_in_loop: set[tuple[str, int]] = set()
    read_tys: dict[tuple[str, int], TypeIR] = {}
    # 阶段 A 续修：扫全 cfg 收集 (name, version) → ty 表（caller-scope LetStmt 的
    # ty 比 read 处 Var.ty 更可靠；当 Var.ty=None / UnknownType 时 fallback）
    cfg_def_tys: dict[tuple[str, int], TypeIR] = {}
    def _scan_cfg(c: CFG):
        for bb_g in c.blocks.values():
            for s_g in bb_g.stmts:
                if isinstance(s_g, LetStmt):
                    if s_g.ty is not None and not isinstance(s_g.ty, UnknownType):
                        cfg_def_tys.setdefault((s_g.var.name, s_g.var.version), s_g.ty)
                elif isinstance(s_g, PhiStmt):
                    if s_g.ty is not None and not isinstance(s_g.ty, UnknownType):
                        cfg_def_tys.setdefault((s_g.target.name, s_g.target.version), s_g.ty)
    _scan_cfg(cfg)
    # 阶段 A 续修：扫已提取的 aux_defs（new_aux）的 cfg —— nested loop 时外层
    # cap_param 引用的 var 可能定义在内层（已提取的）loop 的 cfg 中
    if extra_aux_for_ty:
        for aux in extra_aux_for_ty:
            if aux.cfg is not None:
                _scan_cfg(aux.cfg)
    # header BB 的 phi 是 loop-carried — sources 要么在循环外（init arg），要么
    # 在循环内某 body BB（loop-defined）。body 内非-header BB 的 phi（Pass 6 在
    # 嵌套合流点也可能插 phi）—— 其 sources 可能引用 outer-scope Var。
    # 策略：对所有 phi.sources 都扫；下面的 filter（name ∈ defs_in_loop）会自然
    # 过滤掉 loop-defined 的 source，外层引用即落入 free_vars。
    for bb_id in body_bbs:
        bb = cfg.blocks[bb_id]
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                defs_in_loop.add(s.target.name)
                # P1-rescan-phi（第十二轮）：扫 phi.sources 把 outer-scope 引用纳入；
                # 若 source 是 loop-defined（含 header phi 的循环内分支），下面 filter
                # 阶段 `name in defs_in_loop` 会过滤掉，无需在此区分。
                for src in s.sources.values():
                    if isinstance(src, Var):
                        reads_in_loop.add((src.name, src.version))
                        if (src.name, src.version) not in read_tys \
                                and src.ty is not None:
                            read_tys[(src.name, src.version)] = src.ty
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
    # P1-cap-blacklist（第十二轮）：跳过合成名前缀（lambda/loop ref / 临时返回值），
    # 与 rescan 段（line ~538）保持过滤一致——这些是函数引用 / Pass 7 内部生成
    # 名，不应作为 cap_param 进 loop_func.params。
    free: list[Var] = []
    for (name, version) in sorted(reads_in_loop):
        if name in defs_in_loop:
            continue
        if name.startswith(("_lambda_", "_loop_", "__loop_ret_",
                              "__refret_", "__hoist_lam")):
            continue
        # 阶段 A/C 续修：read_tys → cfg_def_tys → func_param_tys (version=0
        # 才匹配) → UnknownType 四级 fallback
        # 阶段 G-A 修复：read_tys 若是 NamedType("Lambda")/"LambdaRef" 占位，
        # 优先 cfg_def_tys（LetStmt 实际 ty，如 PolyZp / FuncType）
        ty = read_tys.get((name, version))
        is_lambda_placeholder = (isinstance(ty, NamedType)
                                   and ty.name in ("Lambda", "LambdaRef"))
        if ty is None or isinstance(ty, UnknownType) or is_lambda_placeholder:
            cfg_ty = cfg_def_tys.get((name, version))
            if cfg_ty is not None:
                ty = cfg_ty
        if (ty is None or isinstance(ty, UnknownType)) and version == 0 \
                and func_param_tys and name in func_param_tys:
            param_ty = func_param_tys[name]
            if param_ty is not None and not isinstance(param_ty, UnknownType):
                ty = param_ty
        if ty is None:
            ty = UnknownType("")
        free.append(Var(name=name, version=version, ty=ty))
    return free


def _var_param_name(v: Var) -> str:
    """SSA Var → HIRParam name（含版本编码）。"""
    return f"{v.name}_{v.version}" if v.version > 0 else v.name


def _strip_ref_ty(ty: TypeIR | None) -> TypeIR:
    """剥 RefType 外壳；递归处理 Pair / Tuple / Array / Option / StdMap 内层。
    阶段 A 修复：Pass 7 emit 的 cap_param / live_out / ret_ty 不应含 RefType。
    """
    if ty is None:
        return UnknownType("")
    if isinstance(ty, RefType):
        return _strip_ref_ty(ty.inner)
    if isinstance(ty, ArrayType):
        inner = _strip_ref_ty(ty.elem)
        return ArrayType(elem=inner) if inner is not ty.elem else ty
    if isinstance(ty, PairType):
        f = _strip_ref_ty(ty.fst); s = _strip_ref_ty(ty.snd)
        return PairType(fst=f, snd=s) if (f is not ty.fst or s is not ty.snd) else ty
    if isinstance(ty, TupleType):
        new_elems = tuple(_strip_ref_ty(e) for e in ty.elems)
        return TupleType(elems=new_elems) if new_elems != ty.elems else ty
    if isinstance(ty, OptionType):
        inner = _strip_ref_ty(ty.inner)
        return OptionType(inner=inner) if inner is not ty.inner else ty
    if isinstance(ty, StdMapType):
        k = _strip_ref_ty(ty.key); v = _strip_ref_ty(ty.value)
        return StdMapType(key=k, value=v) if (k is not ty.key or v is not ty.value) else ty
    return ty


def _reaching_def_at(cfg: CFG, name: str, src_bb: int,
                     idom: dict[int, int]) -> Var | None:
    """Find latest def of `name` reaching src_bb（沿 dominator 链）。

    返回 Var 或 None（无 def 可达——通常是 SSA 漏 def 的征兆）。
    """
    cur = src_bb
    visited: set[int] = set()
    while True:
        if cur in visited:
            return None
        visited.add(cur)
        bb = cfg.blocks.get(cur)
        if bb is None:
            return None
        for s in reversed(bb.stmts):
            if isinstance(s, LetStmt) and s.var.name == name:
                return s.var
            if isinstance(s, PhiStmt) and s.target.name == name:
                return s.target
        nxt = idom.get(cur, cur)
        if nxt == cur:
            return None
        cur = nxt


def _collect_loop_live_outs(cfg: CFG, body_bbs: set[int]) -> list[Var]:
    """收集 loop body 内 def 且在 caller 端 (非-body BB) 被 read 的 SSA Var。

    包括：
    - body 内 LetStmt.var
    - body 内 PhiStmt.target（含 header phi）
    满足"在非-body BB 中被引用"。

    返回按 (name, version) 排序的 Var 列表。
    """
    body_def_vars: dict[tuple[str, int], Var] = {}
    for bb_id in body_bbs:
        for s in cfg.blocks[bb_id].stmts:
            if isinstance(s, LetStmt):
                body_def_vars[(s.var.name, s.var.version)] = s.var
            elif isinstance(s, PhiStmt):
                body_def_vars[(s.target.name, s.target.version)] = s.target

    non_body_reads: set[tuple[str, int]] = set()
    read_tys: dict[tuple[str, int], TypeIR] = {}
    for bb_id, bb in cfg.blocks.items():
        if bb_id in body_bbs:
            continue
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.value is not None:
                _collect_var_reads_in_expr_versioned(s.value, non_body_reads, read_tys)
            elif isinstance(s, PhiStmt):
                for src_var in s.sources.values():
                    if isinstance(src_var, Var):
                        non_body_reads.add((src_var.name, src_var.version))
            elif isinstance(s, RequireStmt):
                _collect_var_reads_in_expr_versioned(s.cond, non_body_reads, read_tys)
        t = bb.terminator
        if isinstance(t, ReturnTerm) and t.value is not None:
            _collect_var_reads_in_expr_versioned(t.value, non_body_reads, read_tys)
        elif isinstance(t, CondJumpTerm):
            _collect_var_reads_in_expr_versioned(t.cond, non_body_reads, read_tys)
        elif isinstance(t, TailCallTerm):
            for a in t.args:
                _collect_var_reads_in_expr_versioned(a, non_body_reads, read_tys)

    live_outs: list[Var] = []
    for key, var in body_def_vars.items():
        if key in non_body_reads:
            live_outs.append(var)
    return sorted(live_outs, key=lambda v: (v.name, v.version))


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
                     idom: dict[int, int], loop_id: int,
                     host_base: str,
                     extra_aux_for_ty: list[MIRFunc] | None = None,
                     func_param_tys: dict[str, TypeIR] | None = None
                     ) -> tuple[MIRFunc, list[Var], list[int], list[Var], list[Var]]:
    """从 cfg 提取 loop 为独立 MIRFunc。

    返回 (loop_func, phi_target_vars, exit_targets, free_vars, live_outs)。
    - phi_target_vars: caller 端做 init args 用
    - live_outs: caller 端做 destructure 用（含 phi targets 若 post-loop 用）

    P0-A 修复（第十一轮）：loop 名 `_loop_{host_base}_{id}` 避免全局命名冲突
    （多个 host 都从 id=0 开始 → 旧的 `_loop_0` 在 aux_defs 跨 host 重名，导致
    codegen 时同名碰撞 + arity 校验歧义）。
    """
    loop_name = f"_loop_{host_base}_{loop_id}"
    header_block = cfg.blocks[header]

    # 1. params: phi targets at header + free vars（含版本号编入 name）
    phi_targets: list[tuple[Var, TypeIR]] = _collect_phi_targets_at(header_block)
    phi_target_vars = [Var(name=t[0].name, version=t[0].version, ty=t[1])
                        for t in phi_targets]
    free_vars: list[Var] = _collect_loop_free_vars(cfg, body_bbs, phi_targets,
                                                     extra_aux_for_ty=extra_aux_for_ty,
                                                     func_param_tys=func_param_tys)
    # 阶段 A：cap_param ty 剥 RefType（Pass 7 silent regression 修复）
    cap_params = [HIRParam(name=_var_param_name(v), ty=_strip_ref_ty(v.ty),
                            is_ref=False, is_const_ref=False, is_output=False)
                  for v in free_vars]
    phi_params = [HIRParam(name=_var_param_name(t[0]), ty=t[1], is_ref=False,
                            is_const_ref=False, is_output=False)
                  for t in phi_targets]
    new_params = phi_params + cap_params

    # 2. exit targets（按 first-encountered 顺序分配 kind）
    exit_targets: list[int] = []
    for bb_id in sorted(body_bbs):
        for s in cfg.succs.get(bb_id, []):
            if s not in body_bbs and s not in exit_targets:
                exit_targets.append(s)
    n_exits = len(exit_targets)
    exit_kind: dict[int, int] = {t: i for i, t in enumerate(exit_targets)}

    # 3. P7-6 修复：live_outs = body 内 def 且 caller 端 read 的 SSA Var。
    # 这替代之前"仅 phi targets 作 ret"的简化模型——含 body 内 def 的非-phi var
    # 也作为 ret tuple 元素（如 __edf_Zp 的 g_8）。
    live_outs: list[Var] = _collect_loop_live_outs(cfg, body_bbs)
    # ret_ty: (exit_kind, ...live_outs)
    ret_tys: list[TypeIR] = [BaseType.INT64]
    # 阶段 A：live_out ret_ty 也剥 RefType
    ret_tys.extend(_strip_ref_ty(v.ty) for v in live_outs)
    if len(ret_tys) == 1:
        ret_ty: TypeIR = ret_tys[0]  # 仅 kind（无 live_out）
    elif len(ret_tys) == 2:
        ret_ty = PairType(ret_tys[0], ret_tys[1])
    else:
        ret_ty = TupleType(tuple(ret_tys))

    # 4. body cfg copy + 改写 terminator（用 reaching-def 修 P7-7）
    new_blocks: dict[int, BasicBlock] = {}
    for bb_id in body_bbs:
        bb = cfg.blocks[bb_id]
        new_stmts: list[MIRStmt] = []
        for s in bb.stmts:
            if bb_id == header and isinstance(s, PhiStmt):
                continue
            new_stmts.append(s)
        new_term = _rewrite_loop_terminator(
            bb, body_bbs, header, header_block, loop_name, free_vars,
            exit_kind, n_exits, live_outs, idom, cfg)
        new_blocks[bb_id] = BasicBlock(bb_id=bb_id, stmts=new_stmts,
                                        terminator=new_term)

    fresh_bb_id = [max(cfg.blocks.keys(), default=0) + 100 + loop_id * 1000]
    _split_cond_to_exit_or_back(new_blocks, body_bbs, header, header_block,
                                  loop_name, free_vars,
                                  exit_kind, n_exits, fresh_bb_id,
                                  live_outs, idom, cfg)

    # P7-6 续：terminator rewrite 后，ReturnTerm 通过 reaching-def 引入新 Var
    # 引用（如 outer-scope 的 vp_1）。重扫 new_blocks 找漏的 free vars，加入
    # params。
    all_def_keys: set[tuple[str, int]] = set()
    all_def_keys.update((p.name, 0) for p in cap_params)  # free vars by name
    for v, _ty in phi_targets:
        all_def_keys.add((v.name, v.version))
    for v in free_vars:
        all_def_keys.add((v.name, v.version))
    for bb in new_blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt):
                all_def_keys.add((s.var.name, s.var.version))
            elif isinstance(s, PhiStmt):
                all_def_keys.add((s.target.name, s.target.version))
    extra_reads: set[tuple[str, int]] = set()
    extra_read_tys: dict[tuple[str, int], TypeIR] = {}
    for bb in new_blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.value is not None:
                _collect_var_reads_in_expr_versioned(s.value, extra_reads, extra_read_tys)
            elif isinstance(s, RequireStmt):
                _collect_var_reads_in_expr_versioned(s.cond, extra_reads, extra_read_tys)
            elif isinstance(s, PhiStmt):
                # P1-rescan-phi（第十二轮）：body 内非-header BB 的 PhiStmt
                # （header phi 已在 482 行剔除），其 sources 可能引用 reaching-def
                # 后的 outer-scope Var。漏扫导致 loop_func 内 undef 引用。
                for src in s.sources.values():
                    if isinstance(src, Var):
                        extra_reads.add((src.name, src.version))
                        if (src.name, src.version) not in extra_read_tys \
                                and src.ty is not None:
                            extra_read_tys[(src.name, src.version)] = src.ty
        t = bb.terminator
        if isinstance(t, ReturnTerm) and t.value is not None:
            _collect_var_reads_in_expr_versioned(t.value, extra_reads, extra_read_tys)
        elif isinstance(t, CondJumpTerm):
            _collect_var_reads_in_expr_versioned(t.cond, extra_reads, extra_read_tys)
        elif isinstance(t, TailCallTerm):
            for a in t.args:
                _collect_var_reads_in_expr_versioned(a, extra_reads, extra_read_tys)
    extra_free: list[Var] = []
    for (n, v) in sorted(extra_reads):
        if (n, v) in all_def_keys:
            continue
        if v == 0 and (n, 0) in all_def_keys:
            continue
        # 跳过函数引用 / loop ret / refret
        if n.startswith(("_lambda_", "_loop_", "__loop_ret_",
                          "__refret_", "__hoist_lam")):
            continue
        ty = extra_read_tys.get((n, v), UnknownType(""))
        extra_free.append(Var(name=n, version=v, ty=ty))
    if extra_free:
        free_vars = list(free_vars) + extra_free
        for ev in extra_free:
            cap_params.append(HIRParam(name=_var_param_name(ev),
                                         ty=_strip_ref_ty(ev.ty),
                                         is_ref=False, is_const_ref=False,
                                         is_output=False))
        new_params = phi_params + cap_params
        # P0-A 修复（第十一轮审视）：terminator 改写发生在重扫之前，已构造的
        # TailCallTerm.args 用的是初始 free_vars。这里追加 extra_free 进所有
        # 自调用的 TailCallTerm.args，否则 args 数 < params 数 → silent arity 错位。
        for bb_id, bb in list(new_blocks.items()):
            t = bb.terminator
            if isinstance(t, TailCallTerm) and t.target_func == loop_name:
                new_args = list(t.args) + list(extra_free)
                new_blocks[bb_id] = BasicBlock(
                    bb_id=bb_id, stmts=bb.stmts,
                    terminator=TailCallTerm(target_func=loop_name,
                                             args=new_args),
                )

    # P0-A 修复（续）：自检——所有自调用 TailCallTerm.args 数必须 == new_params 数。
    expected_arity = len(new_params)
    for bb_id, bb in new_blocks.items():
        t = bb.terminator
        if isinstance(t, TailCallTerm) and t.target_func == loop_name:
            if len(t.args) != expected_arity:
                raise TranslationError(
                    pass_name="loop_lower",
                    func_name=loop_name,
                    reason=(f"bb[{bb_id}]: TailCallTerm arity mismatch — "
                            f"args={len(t.args)} expected={expected_arity}"),
                )

    new_cfg = CFG(entry=header, blocks=new_blocks)
    new_cfg.rebuild_edges()

    loop_func = MIRFunc(
        base_name=loop_name,
        instance_suffix="",
        mangled_name="",
        qual_type=f"loop in {host_base} | id={loop_id}",
        params=new_params,
        ret_ty=ret_ty,
        requires=[],
        cfg=new_cfg,
        aux_defs=[],
    )
    return loop_func, phi_target_vars, exit_targets, free_vars, live_outs


def _rewrite_loop_terminator(bb: BasicBlock, body_bbs: set[int],
                              header: int, header_block: BasicBlock,
                              loop_name: str, free_vars: list[Var],
                              exit_kind: dict[int, int],
                              n_exits: int,
                              live_outs: list[Var],
                              idom: dict[int, int],
                              orig_cfg: CFG):
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
            return _make_exit_return(t.target, exit_kind, src_bb_id,
                                      live_outs, idom, orig_cfg, n_exits)
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
                       src_bb_id: int, live_outs: list[Var],
                       idom: dict[int, int], orig_cfg: CFG,
                       n_exits: int) -> ReturnTerm:
    """exit edge → ReturnTerm((kind=k, ...live_out_values_at_src_bb))。

    P7-7 修复：用 reaching-def 找 src_bb 处每个 live_out.name 的最新版本，
    替代之前用 phi.target 自身的错误 fallback。
    """
    from ir_types import TupleExpr
    k = exit_kind[target]
    elems: list[ExprIR] = [Lit(value=k, ty=BaseType.INT64)]
    for lo in live_outs:
        rd = _reaching_def_at(orig_cfg, lo.name, src_bb_id, idom)
        if rd is None:
            # 阶段 F #1 修复：find-loop 模式中 live_out（如 `term_1`）只在循环
            # body 内部 def，循环 0 次进入时 reaching-def 不存在 → 用 Lean
            # `default` 关键字作占位（要求 Var.ty 有 Inhabited 实例）。
            rd = Var(name="default", version=0, ty=lo.ty)
        elems.append(rd)
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
                                 exit_kind: dict[int, int], n_exits: int,
                                 fresh_bb_id: list[int],
                                 live_outs: list[Var],
                                 idom: dict[int, int],
                                 orig_cfg: CFG) -> None:
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
                term = _make_exit_return(tgt, exit_kind, bb_id, live_outs,
                                          idom, orig_cfg, n_exits)
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
    # 阶段 C：func.params 作为 caller-scope ty 来源（free_var name == param name 时
    # 使用 param.ty）
    func_param_tys = {p.name: p.ty for p in func.params}

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
        # 跨实例命名冲突修复（与 Pass 3、Pass 4 保持一致）：含 instance_suffix
        host_base = func.base_name
        if func.instance_suffix:
            host_base = f"{host_base}_{func.instance_suffix}"
        # 阶段 A 续修：把 new_aux（已提取的 loops + lifted lambdas）的 cfg 也
        # 纳入 cfg_def_tys 来源——nested loop 提取时 outer cap_param 引用的 var
        # 可能在 inner loop 的 cfg 中定义（已被提取掉，main cfg 找不到）
        loop_func, phi_target_vars, exit_targets, free_vars, live_outs = \
            _build_loop_func(cfg, header, body_bbs, idom, loop_id_counter,
                             host_base=host_base, extra_aux_for_ty=new_aux,
                             func_param_tys=func_param_tys)
        loop_id_counter += 1
        new_aux.append(loop_func)

        # 主 cfg 中替换 loop：header 改为 call block + dispatch BBs
        cfg = _splice_loop_call(cfg, header, body_bbs, loop_func,
                                  phi_target_vars, exit_targets, free_vars,
                                  live_outs)

    return replace(func, cfg=cfg, aux_defs=new_aux)


def _splice_loop_call(cfg: CFG, header: int, body_bbs: set[int],
                       loop_func: MIRFunc,
                       phi_target_vars: list[Var],
                       exit_targets: list[int],
                       free_vars: list[Var],
                       live_outs: list[Var]) -> CFG:
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

    # init args: phi sources at h from outer pred + free vars
    outer_preds = [p for p in cfg.preds.get(header, []) if p not in body_bbs]
    if not outer_preds:
        raise TranslationError("loop_lower", loop_func.base_name,
                                "loop header has no outer pred")
    init_outer_pred = outer_preds[0]
    init_phi_srcs = _collect_phi_sources_from(header_block, init_outer_pred)
    init_args: list[ExprIR] = list(init_phi_srcs)
    init_args.extend(free_vars)

    n_exits = len(exit_targets)
    if n_exits == 0:
        raise TranslationError("loop_lower", loop_func.base_name,
                                "loop has no exit edges")

    # tmp var receiving loop func result
    # P0-A 修复（续）：含 host 的 loop_name → tmp 名也唯一
    tmp_name = f"__loop_ret_{loop_func.base_name.removeprefix('_loop_')}"
    tmp_var = Var(name=tmp_name, version=1, ty=loop_func.ret_ty)

    # P7-6 修复：destructure 字段以 live_outs 为准（替代之前的 phi_targets）
    # B1 续修（2026-05-02）：n>2 tuple 直接生成 Lean 投影路径（"1", "2.1",
    # "2.2.1", ..., "2.2...2"）—— Pass 8 emit FieldAccess 直接输出 obj.<path>
    n_lo = len(live_outs)
    n_total = 1 + n_lo  # kind + live_outs
    if n_total == 1:
        kind_field: str | None = None  # 仅 kind
        lo_fields: list[str] = []
    elif n_total == 2:
        kind_field = "fst"
        lo_fields = ["snd"]
    else:
        # n > 2: nested Prod (a, (b, (c, d)))
        #   elem0 (kind) → .1
        #   elem1 → .2.1
        #   elem2 → .2.2.1 (middle) or .2.2 (last)
        kind_field = "1"
        lo_fields = []
        for k in range(n_lo):
            # k 是 lo 内索引；对应 tuple 总位置 k+1（kind 在位置 0）
            n_dots = k + 1  # 多少个 .2
            if (k + 1) == n_total - 1:
                # 最后一个：纯 .2 链
                lo_fields.append("2" * 1 if n_dots == 1 else ".".join(["2"] * n_dots))
            else:
                # 中间：.2 链 + .1
                lo_fields.append(".".join(["2"] * n_dots) + ".1")

    new_header_stmts: list[MIRStmt] = []
    new_header_stmts.append(LetStmt(
        var=tmp_var, ty=loop_func.ret_ty,
        value=Call(callee=loop_func.base_name, args=init_args,
                    ty=loop_func.ret_ty),
    ))
    # destructure live_outs：每个 live_out 在 caller 端被恢复为其原 SSA Var
    # （body 的 def 已删，destructure LetStmt 即为新 def，不破坏 SSA 单赋值）
    for k, lo in enumerate(live_outs):
        # 阶段 A：destructure LetStmt ty 剥 RefType
        lo_ty = _strip_ref_ty(lo.ty)
        new_header_stmts.append(LetStmt(
            var=lo, ty=lo_ty,
            value=FieldAccess(obj=tmp_var, field_name=lo_fields[k], ty=lo_ty)
                  if lo_fields else tmp_var,
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
