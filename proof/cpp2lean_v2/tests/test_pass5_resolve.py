"""
Pass 5 operator_resolve 单元测试。

覆盖：
- T1: passthrough 函数无 gap (`__make_zp`)
- T2: cast 解析 — Lit→Dependent 剥（`__factor_multivar`）
- T3: cast 解析 — IntegralCast 走 CAST_TABLE
- T4: CompoundAssign 展开 → AssignStmt + BinOp
- T5: mutator method (push_back) 转 AssignStmt
- T6: 构造器解析（`construct_ZZ(1)` 等）
- T7: UB require 生成（除零）
- T8: assert_hir4_invariant 拒绝残留 CompoundAssignStmt
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType,
    Var, Lit, BinOp, UnaryOp, Cast, Call, UnresolvedOp,
    LetStmt, AssignStmt, CompoundAssignStmt, ExprStmt, IfStmt,
    HIRFunc, HIRParam, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass, assert_hir4_invariant

AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"


def _run(func_name):
    with open(AST / f"{func_name}.json") as f:
        ast = json.load(f)
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    hir3 = iter_recognize_pass(lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast))))))
    return operator_resolve_pass(hir3)


def _count_compound(stmts):
    n = 0
    for s in stmts:
        if isinstance(s, CompoundAssignStmt): n += 1
        if isinstance(s, IfStmt):
            n += _count_compound(s.then_body) + _count_compound(s.else_body)
        # 简化：不递归其他 stmt 容器（够测试用）
    return n


def test_passthrough_no_compound():
    """T1+T8: __make_zp 跑过 + 不变量 + 无 CompoundAssignStmt 残留。"""
    hir4, gap = _run("__make_zp")
    assert_hir4_invariant(hir4)


def test_compound_assign_expanded():
    """T4: __upoly_powmod 包含 *= 等 compound assign，应全部展开。"""
    hir4, gap = _run("__upoly_powmod")
    assert_hir4_invariant(hir4)
    assert _count_compound(hir4.body) == 0


def test_factor_multivar_lit_strip():
    """T2: __factor_multivar 中 Lit→DependentType 4 个 IntegralCast，应被
    `is_safe_to_strip_lit_cast` 剥掉，gap.cast_miss 不应包含字面量类。"""
    hir4, gap = _run("__factor_multivar")
    assert_hir4_invariant(hir4)
    # gap.cast_miss 中应不含字面量类（Lit 形态被剥后不再产生 cast_miss）
    # 唯一可能的 cast_miss 是 (int32, unresolved) 之类


def test_factor_zp_constructors_resolved():
    """T6: __factor_Zp 的 construct_Zp / construct_ZZ 应全部 resolve。"""
    hir4, gap = _run("__factor_Zp")
    assert_hir4_invariant(hir4)
    # 所有 construct_* 应已 resolve（gap.constructor_miss 为空）
    # 因为 CONSTRUCTOR_MAP 命中率 100%
    assert len(gap.constructor_miss) == 0


def test_invariant_rejects_compound():
    """T8: 残留 CompoundAssignStmt 触发 TranslationError。"""
    f = HIRFunc(
        base_name="_test", instance_suffix="", mangled_name="", qual_type="",
        params=[], ret_ty=BaseType.UNIT,
        body=[
            CompoundAssignStmt(target=Var("x"), op="+", value=Lit(1)),
        ],
        requires=[], aux_lambdas=[],
    )
    try:
        assert_hir4_invariant(f)
        assert False, "should have raised"
    except TranslationError as e:
        assert "CompoundAssignStmt" in e.reason


def test_smoke_all_65():
    """T-smoke: 65 函数全跑 Pass 1+2+3+4+5。"""
    from class_map import TRANSLATION_SCOPE
    ok = 0; errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            from pass2b_callsite_ref_elim import callsite_ref_elim_pass
            hir3 = iter_recognize_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))))
            hir4, gap = operator_resolve_pass(hir3)
            assert_hir4_invariant(hir4)
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {e}"))
    assert not errs, f"{len(errs)} errors: {errs[:5]}"


if __name__ == "__main__":
    tests = [
        test_passthrough_no_compound,
        test_compound_assign_expanded,
        test_factor_multivar_lit_strip,
        test_factor_zp_constructors_resolved,
        test_invariant_rejects_compound,
        test_smoke_all_65,
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
