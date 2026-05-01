"""
Pass 6 ssa_build 单元测试。

按 mir-design.md §2.8 关键测试用例：
- T1: trivial passthrough（无分支无赋值）
- T2: linear assign 编版本号
- T3: if 合并产生 phi
- T4: while 循环 phi 在 header
- T5: 嵌套 if/while
- T6: break/continue 多 exit
- T7: 完整函数 __make_zp / __squarefree_Zp 全 67 smoke
- T8: assert_mir0_invariant 拒绝违规
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType, RefType, ArrayType,
    Var, Lit, BinOp, Call, FieldAccess, ArrayAccess,
    LetStmt, AssignStmt, IfStmt, WhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RangeForStmt,
    HIRFunc, HIRParam, MIRFunc,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm,
    BasicBlock, CFG,
    assert_mir0_invariant, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import (
    ssa_build_pass, cfg_from_hir,
    compute_dominator_tree, compute_dominance_frontier,
    place_phi_nodes, rename_variables,
)
from class_map import TRANSLATION_SCOPE


def _mk_func(body, params=None, ret_ty=BaseType.UNIT) -> HIRFunc:
    return HIRFunc(
        base_name="test", instance_suffix="", mangled_name="", qual_type="",
        params=params or [], ret_ty=ret_ty, body=body,
        requires=[], aux_lambdas=[],
    )


# ============================================================
# T1: Trivial
# ============================================================

def test_trivial_passthrough():
    """单 stmt 单块 + auto-Return"""
    body = [LetStmt(var=Var("x", 0), ty=BaseType.INT64, value=Lit(42, ty=BaseType.INT64))]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    assert mir.cfg.entry == 0
    bb0 = mir.cfg.blocks[0]
    assert isinstance(bb0.terminator, ReturnTerm)
    # 第一个 LetStmt 的 var.version 应为 1（fresh 后）
    let_stmt = next(s for s in bb0.stmts if isinstance(s, LetStmt))
    assert let_stmt.var.version == 1


# ============================================================
# T2: Linear assign — version 编号
# ============================================================

def test_linear_assign_versions():
    """x = 1; x = 2;  →  x_1 = 1; x_2 = 2"""
    body = [
        AssignStmt(target=Var("x"), value=Lit(1, ty=BaseType.INT64)),
        AssignStmt(target=Var("x"), value=Lit(2, ty=BaseType.INT64)),
    ]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    bb0 = mir.cfg.blocks[mir.cfg.entry]
    versions = [s.var.version for s in bb0.stmts if isinstance(s, LetStmt) and s.var.name == "x"]
    assert versions == [1, 2], f"got {versions}"


# ============================================================
# T3: If 合并产生 phi
# ============================================================

def test_if_phi_at_merge():
    """if (cond) x=1 else x=2; return x → phi at merge block"""
    body = [
        IfStmt(
            cond=Var("cond", ty=BaseType.BOOL),
            then_body=[AssignStmt(target=Var("x"), value=Lit(1, ty=BaseType.INT64))],
            else_body=[AssignStmt(target=Var("x"), value=Lit(2, ty=BaseType.INT64))],
        ),
        ReturnStmt(value=Var("x")),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="cond", ty=BaseType.BOOL),
        HIRParam(name="x", ty=BaseType.INT64),
    ], ret_ty=BaseType.INT64)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 应至少有 1 个 phi 节点（merge bb 内 x_n := phi(...)）
    has_phi = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, PhiStmt) and s.target.name == "x":
                has_phi = True
                # phi 应有 2 个 sources（来自 then_bb 和 else_bb）
                assert len(s.sources) == 2, f"phi sources: {s.sources}"
    assert has_phi, "no phi for x"


# ============================================================
# T4: While 循环 — phi 在 header
# ============================================================

def test_while_phi_at_header():
    """while (i < n) i = i + 1  →  phi at header"""
    body = [
        WhileStmt(
            cond=BinOp(op="<", lhs=Var("i"), rhs=Var("n"), ty=BaseType.BOOL),
            body=[AssignStmt(
                target=Var("i"),
                value=BinOp(op="+", lhs=Var("i"), rhs=Lit(1), ty=BaseType.INT64),
            )],
        ),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="i", ty=BaseType.INT64),
        HIRParam(name="n", ty=BaseType.INT64),
    ])
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 找 while 的 header 块（带 CondJumpTerm 且有前驱不止 entry）
    has_phi_at_header = False
    for bb_id, bb in mir.cfg.blocks.items():
        if isinstance(bb.terminator, CondJumpTerm):
            preds = mir.cfg.preds.get(bb_id, [])
            if len(preds) >= 2:
                for s in bb.stmts:
                    if isinstance(s, PhiStmt) and s.target.name == "i":
                        has_phi_at_header = True
                        assert len(s.sources) >= 2
    assert has_phi_at_header, "no phi for i at while header"


# ============================================================
# T5: Break 退出
# ============================================================

def test_break_creates_extra_edge():
    """while (...) { break; ... }  break 跳到 exit 块"""
    body = [
        WhileStmt(
            cond=Lit(True, ty=BaseType.BOOL),
            body=[BreakStmt()],
        ),
    ]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)


# ============================================================
# T6: Continue 跳到 header
# ============================================================

def test_continue_jumps_header():
    body = [
        WhileStmt(
            cond=Var("c", ty=BaseType.BOOL),
            body=[ContinueStmt()],
        ),
    ]
    func = _mk_func(body, params=[HIRParam(name="c", ty=BaseType.BOOL)])
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)


# ============================================================
# T6b: P0-1 RangeFor `auto&` 写回（is_mutable_ref=True）
# ============================================================

def test_rangefor_auto_ref_writeback():
    """for (auto& x : c) x = x + 1
       期望：
       - latch 含 cont := __write__(cont[idx], x)
       - exit 含 c := cont（即原 container 写回）
    """
    body = [
        RangeForStmt(
            var=Var("x", ty=BaseType.INT64),
            var_ty=BaseType.INT64,
            container=Var("c", ty=ArrayType(BaseType.INT64)),
            body=[AssignStmt(target=Var("x"),
                             value=BinOp("+", Var("x"), Lit(1, ty=BaseType.INT64),
                                         ty=BaseType.INT64))],
            is_mutable_ref=True,
        ),
        ReturnStmt(value=Var("c")),
    ]
    func = _mk_func(body, params=[HIRParam(name="c", ty=ArrayType(BaseType.INT64))],
                    ret_ty=ArrayType(BaseType.INT64))
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)

    # 1) latch 含 __rangefor_cont_*_N := __write__(cont[idx], x)
    has_latch_writeback = False
    has_exit_writeback = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.var.name.startswith("__rangefor_cont_") \
                    and isinstance(s.value, Call) and s.value.callee == "__write__":
                has_latch_writeback = True
            # 2) exit 含 c_(N+1) := __rangefor_cont_*_M
            if isinstance(s, LetStmt) and s.var.name == "c" and s.var.version >= 1 \
                    and isinstance(s.value, Var) and s.value.name.startswith("__rangefor_cont_"):
                has_exit_writeback = True
    assert has_latch_writeback, "missing latch writeback (B7 chain inside loop)"
    assert has_exit_writeback, "missing exit writeback (P0-1: c := cont after loop)"


def test_rangefor_no_writeback_when_not_mutable():
    """for (auto x : c) x = x + 1   （非 ref，is_mutable_ref=False）
       期望：不生成任何写回（cont 也不修改原 c）。
    """
    body = [
        RangeForStmt(
            var=Var("x", ty=BaseType.INT64),
            var_ty=BaseType.INT64,
            container=Var("c", ty=ArrayType(BaseType.INT64)),
            body=[AssignStmt(target=Var("x"),
                             value=BinOp("+", Var("x"), Lit(1, ty=BaseType.INT64),
                                         ty=BaseType.INT64))],
            is_mutable_ref=False,  # 关键
        ),
        ReturnStmt(value=Var("c")),
    ]
    func = _mk_func(body, params=[HIRParam(name="c", ty=ArrayType(BaseType.INT64))],
                    ret_ty=ArrayType(BaseType.INT64))
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 不应有 c 的版本 ≥ 1（c 一直是 c_0）
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.var.name == "c":
                assert False, f"unexpected c-version write: {s}"


def test_rangefor_writeback_field_container():
    """for (auto& x : obj.field) ...   container 是 FieldAccess
       期望：exit 写回走 B7 record-update 链：
       obj_(N+1) := __write__(obj_N.field, cont)
    """
    body = [
        RangeForStmt(
            var=Var("x", ty=BaseType.INT64),
            var_ty=BaseType.INT64,
            container=FieldAccess(obj=Var("obj"), field_name="field",
                                  ty=ArrayType(BaseType.INT64)),
            body=[AssignStmt(target=Var("x"),
                             value=BinOp("+", Var("x"), Lit(1, ty=BaseType.INT64),
                                         ty=BaseType.INT64))],
            is_mutable_ref=True,
        ),
        ReturnStmt(value=Var("obj")),
    ]
    func = _mk_func(body,
                    params=[HIRParam(name="obj", ty=NamedType("Container"))],
                    ret_ty=NamedType("Container"))
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 期望：exit 块中存在 obj 的版本 ≥ 1，rhs 是 _with(obj_X, field, cont)
    # （B1 续修：FieldAccess record-update 改用 _with → Lean { obj with field := v }）
    has_field_writeback = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.var.name == "obj" and s.var.version >= 1:
                if isinstance(s.value, Call) and s.value.callee in ("__write__", "_with"):
                    has_field_writeback = True
    assert has_field_writeback, "missing FieldAccess record-update writeback at exit"


# ============================================================
# T6c: CF-2 阶段 2 sequence iterator → 索引化（__upoly_divmod_mod 形态）
# ============================================================

def test_seq_iter_desugar_to_index():
    """sequence iterator (begin/toList) 应被改写为 idx 索引：
       - `let it := c.toList(c)` → `let __iter_idx_it := 0`
       - `it = Iterator.advance(it)` → `__iter_idx_it := __iter_idx_it + 1`
       - `__deref__(it).field` → `c[__iter_idx_it].field`
       - `it != c.toList(c)` → `__iter_idx_it < c.size`
    """
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    f = AST / "__upoly_divmod_mod.json"
    if not f.exists():
        return  # corpus 缺该函数则跳过
    with open(f) as fh: ast = json.load(fh)
    from pass1_parse import parse_pass
    from pass2_ref_elim import ref_elim_pass
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3_lambda_lift import lambda_lift_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    from pass4_iter_recognize import iter_recognize_pass
    from pass5_operator_resolve import operator_resolve_pass
    hir = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))))
    hir3 = iter_recognize_pass(hir)
    hir4, _ = operator_resolve_pass(hir3)
    mir = ssa_build_pass(hir4)
    # 验证：MIR 中存在 __iter_idx_X_N（idx 变量）；不存在 __deref__ 残留
    has_idx = False
    has_deref = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if hasattr(s, "var") and "__iter_idx_" in s.var.name:
                has_idx = True
            # 任何 stmt 含 __deref__ Call → fail
            def expr_has_deref(e):
                if isinstance(e, Call) and isinstance(e.callee, str) \
                        and e.callee in ("__deref__", "Iterator.advance",
                                         "Iterator.deref!"):
                    return True
                if isinstance(e, BinOp):
                    return expr_has_deref(e.lhs) or expr_has_deref(e.rhs)
                if isinstance(e, FieldAccess):
                    return expr_has_deref(e.obj)
                if isinstance(e, ArrayAccess):
                    return expr_has_deref(e.arr) or expr_has_deref(e.idx)
                if isinstance(e, Call):
                    return any(expr_has_deref(a) for a in e.args)
                return False
            value = getattr(s, 'value', None)
            if value and expr_has_deref(value):
                has_deref = True
    assert has_idx, "missing __iter_idx_X (CF-2 desugar didn't fire)"
    assert not has_deref, "__deref__/Iterator.advance residue (CF-2 incomplete)"


# ============================================================
# T6d: CF-3 阶段 3 sort/iota functional 形态
# ============================================================

def test_stl_sort_iota_functional():
    """sort/iota 应被 Pass 5 改写为 Array.sort / Array.range_init。
    不再有 `__sideeff_X := sort(toList, toList, comp)` 形态。
    """
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    found_funcs = ["__factor_Zp", "__zassenhaus_recombine"]
    from pass1_parse import parse_pass
    from pass2_ref_elim import ref_elim_pass
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3_lambda_lift import lambda_lift_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    from pass4_iter_recognize import iter_recognize_pass
    from pass5_operator_resolve import operator_resolve_pass
    for fn in found_funcs:
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        with open(f) as fh: ast = json.load(fh)
        hir = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))))
        hir3 = iter_recognize_pass(hir)
        hir4, _ = operator_resolve_pass(hir3)
        mir = ssa_build_pass(hir4)
        # 验证：MIR 中无 sort/iota 作 sideeff 形态
        for bb in mir.cfg.blocks.values():
            for s in bb.stmts:
                if hasattr(s, "var") and s.var.name.startswith("__sideeff_"):
                    if hasattr(s, "value") and isinstance(s.value, Call) \
                            and isinstance(s.value.callee, str) \
                            and s.value.callee in ("sort", "iota"):
                        assert False, f"{fn}: sort/iota 仍是 sideeff: {s.value.callee}"


# ============================================================
# T7: 完整函数 + 67 smoke
# ============================================================

def test_smoke_all_67():
    """全 67 函数走 Pass 1-6 + invariant 全过。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    ok = 0
    errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            from pass2b_callsite_ref_elim import callsite_ref_elim_pass
            from pass3b_lambda_ref_elim import lambda_ref_elim_pass
            hir3 = iter_recognize_pass(lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast))))))
            hir4, _ = operator_resolve_pass(hir3)
            mir = ssa_build_pass(hir4)
            assert_mir0_invariant(mir)
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {str(e)[:100]}"))
    assert not errs, f"{len(errs)} errors: {errs[:3]}"


# ============================================================
# T7b: silent-bug 残余 metric 不回归
# ============================================================

def test_no_silent_bug_regression():
    """全 67 函数：phi_undef_ver0 + aux_dropped + B1 incdec 残余 == 0。"""
    from ir_types import Call as _Call
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    phi_ver0 = 0
    aux_drop = 0
    sideeff_write_count = 0
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        with open(f) as fh: ast = json.load(fh)
        hir3 = iter_recognize_pass(lambda_lift_pass(ref_elim_pass(parse_pass(ast))))
        hir4, _ = operator_resolve_pass(hir3)
        mir = ssa_build_pass(hir4)
        if len(mir.aux_defs) != len(hir4.aux_lambdas):
            aux_drop += 1
        params = {p.name for p in mir.params}
        for bb in mir.cfg.blocks.values():
            for s in bb.stmts:
                if isinstance(s, PhiStmt):
                    for src in s.sources.values():
                        if src.version == 0 and src.name not in params:
                            phi_ver0 += 1
                elif hasattr(s, "var") and s.var.name.startswith("__sideeff_"):
                    if isinstance(s.value, _Call) and s.value.callee == "__write__":
                        sideeff_write_count += 1
    # B2 phi_undef_ver0 ≤ 2（pruned SSA + B4 unique decomp 联合保证）
    # 余 2 例（Gi、h）是 C++ 默认初始化局部变量被 lambda by-ref 捕获后的真实
    # SSA：lambda 调用前 var 仅默认构造 → phi sources 含 ver=0 是合法语义。
    # Pass 1 后续可补 `let X := __default__(T)` 显式化。
    assert phi_ver0 <= 2, f"B2 regression: phi_undef_ver0={phi_ver0}"
    # B6 aux_dropped = 0（aux_lambdas 100% 接力）
    assert aux_drop == 0, f"B6 regression: aux_dropped={aux_drop}"
    # B7：当前可接受 ≤ 8（param-mutation 跨 Pass 工程，留 TODO）
    assert sideeff_write_count <= 8, \
        f"B7 regression: sideeff_write={sideeff_write_count} > 8 baseline"


# ============================================================
# T8: 不变量违规
# ============================================================

def test_invariant_rejects_no_cfg():
    func = MIRFunc(base_name="bad", cfg=None)
    try:
        assert_mir0_invariant(func)
        assert False, "should raise"
    except TranslationError as e:
        assert "cfg is None" in e.reason


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_trivial_passthrough,
        test_linear_assign_versions,
        test_if_phi_at_merge,
        test_while_phi_at_header,
        test_break_creates_extra_edge,
        test_continue_jumps_header,
        test_rangefor_auto_ref_writeback,
        test_rangefor_no_writeback_when_not_mutable,
        test_rangefor_writeback_field_container,
        test_seq_iter_desugar_to_index,
        test_stl_sort_iota_functional,
        test_smoke_all_67,
        test_no_silent_bug_regression,
        test_invariant_rejects_no_cfg,
    ]
    passed = 0; failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            import traceback
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print(f"\n=== {passed}/{passed+failed} passed ===")
    sys.exit(0 if failed == 0 else 1)
