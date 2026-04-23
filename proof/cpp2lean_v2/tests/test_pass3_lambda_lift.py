"""
Pass 3 lambda_lift 单元测试 + 全量烟测。

典型场景：
- T1: 无 lambda 函数（__make_zp）passthrough
- T2: 无捕获 lambda (__factor_Zp sort comparator)
- T3: [&] 捕获 (__mtshl_multi_bdp compute_error)
- T4: 多 lambda 宿主 (__lll_reduce 5 lambdas)
- T5: 65 函数全跑 → HIR₂ 不变量
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, Var, LambdaExpr,
    HIRFunc, HIRParam,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import (
    lambda_lift_pass, assert_hir2_invariant,
    _collect_free_vars, _collect_modified,
)
from class_map import TRANSLATION_SCOPE

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"


def _load_hir1(func_name: str) -> HIRFunc:
    with open(AST_CACHE_DIR / f"{func_name}.json") as f:
        ast = json.load(f)
    return ref_elim_pass(parse_pass(ast))


# ============================================================
# 单元测试
# ============================================================

def test_no_lambda_passthrough():
    """T1: __make_zp 无 lambda → aux_lambdas 空，body 不变。"""
    hir1 = _load_hir1("__make_zp")
    hir2 = lambda_lift_pass(hir1)
    assert hir2.aux_lambdas == [], f"expected empty aux_lambdas"
    assert len(hir2.body) == len(hir1.body)
    assert_hir2_invariant(hir2)


def test_sort_comparator_no_capture():
    """T2: __factor_Zp 有 sort lambda（无捕获）→ 提升 1 个 _lambda__factor_Zp_1。"""
    hir1 = _load_hir1("__factor_Zp")
    hir2 = lambda_lift_pass(hir1)
    assert len(hir2.aux_lambdas) == 1, f"expected 1 lambda, got {len(hir2.aux_lambdas)}"
    aux = hir2.aux_lambdas[0]
    assert aux.base_name == "_lambda___factor_Zp_1"
    # 无捕获 → 参数数等于 lambda 自身的 params（2 个）
    assert len(aux.params) == 2, f"expected 2 params (no captures), got {len(aux.params)}: {[p.name for p in aux.params]}"
    assert_hir2_invariant(hir2)


def test_capture_with_ref():
    """T3: __mtshl_multi_bdp 有 `compute_error` lambda 捕获 `result`、`F`、`c` 等。"""
    hir1 = _load_hir1("__mtshl_multi_bdp")
    hir2 = lambda_lift_pass(hir1)
    assert len(hir2.aux_lambdas) >= 1
    aux = hir2.aux_lambdas[0]
    # 应有捕获（来自外层 func 的 result、F、c 等参数）
    cap_names = [p.name for p in aux.params if p.name in {"result", "F", "c", "x1", "x2", "alpha2"}]
    assert len(cap_names) >= 2, f"expected some captures of outer vars, got {[p.name for p in aux.params]}"
    assert_hir2_invariant(hir2)


def test_multi_lambda_host():
    """T4: __lll_reduce 有 5 lambdas → 提升 5 个 _lambda__lll_reduce_{1..5}。"""
    hir1 = _load_hir1("__lll_reduce")
    hir2 = lambda_lift_pass(hir1)
    # __lll_reduce 实际 lambda 数 = 5（dot/round_qq/row_sub/row_swap/sort 比较器）
    assert len(hir2.aux_lambdas) == 5, \
        f"expected 5 lambdas, got {len(hir2.aux_lambdas)}: {[l.base_name for l in hir2.aux_lambdas]}"
    names = [l.base_name for l in hir2.aux_lambdas]
    for i in range(1, 6):
        assert f"_lambda___lll_reduce_{i}" in names
    assert_hir2_invariant(hir2)


# ============================================================
# 辅助函数单元测试
# ============================================================

def test_free_vars_basic():
    """自由变量分析基础测试。"""
    from ir_types import (Var, LetStmt, ExprStmt, Call, UnresolvedOp,
                          BaseType, NamedType, UnknownType)
    # body: let x := 1; use y (y free, x local)
    body = [
        LetStmt(Var("x"), BaseType.INT32, Call("f", [])),
        ExprStmt(Var("y")),
    ]
    free = _collect_free_vars(body, set())
    assert "y" in free, f"expected 'y' free, got {free}"
    assert "x" not in free, f"x should be local, not free"


def test_modified_detection():
    """修改检测：x = 1 vs x += 1 vs ++x 都算。"""
    from ir_types import (Var, AssignStmt, CompoundAssignStmt, ExprStmt,
                          Call, UnresolvedOp, Lit, BaseType)
    body = [
        AssignStmt(Var("a"), Lit(1)),
        CompoundAssignStmt(Var("b"), "+", Lit(1)),
        ExprStmt(Call(UnresolvedOp("operator++"), [Var("c")])),
    ]
    mod = _collect_modified(body)
    assert mod == {"a", "b", "c"}, f"got {mod}"


# ============================================================
# 全量烟测
# ============================================================

def test_smoke_all_65():
    """T-smoke: 65 函数全跑 Pass 1+2+3，验证 HIR₂ 不变量 + 统计 lambda 数。"""
    n_lambdas_total = 0
    n_funcs_with_lambdas = 0
    errors = []
    for func_name in sorted(TRANSLATION_SCOPE):
        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            continue
        try:
            hir1 = _load_hir1(func_name)
            hir2 = lambda_lift_pass(hir1)
            assert_hir2_invariant(hir2)
            if hir2.aux_lambdas:
                n_funcs_with_lambdas += 1
                n_lambdas_total += len(hir2.aux_lambdas)
        except Exception as e:
            errors.append((func_name, f"{type(e).__name__}: {e}"))

    print(f"\n  Functions with lambdas: {n_funcs_with_lambdas}")
    print(f"  Total lambdas lifted:   {n_lambdas_total}")
    # Week 1 scan_lambdas 记录 26 个 in-scope lambda，14 个宿主
    assert n_lambdas_total >= 20, f"expected >=20 lambdas, got {n_lambdas_total}"
    assert not errors, f"{len(errors)} errors: {errors[:5]}"


if __name__ == "__main__":
    tests = [
        test_free_vars_basic,
        test_modified_detection,
        test_no_lambda_passthrough,
        test_sort_comparator_no_capture,
        test_capture_with_ref,
        test_multi_lambda_host,
        test_smoke_all_65,
    ]
    passed = 0
    failed = 0
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
