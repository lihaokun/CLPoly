"""
Pass 2 ref_elim 单元测试 + 全量烟测。

典型场景：
- T1: void + 2 ref 参数 (__upoly_divmod) → 返回 PairType(SparsePolyZp, SparsePolyZp)
- T2: Zp + 1 ref 参数 (__upoly_make_monic) → 返回 PairType(Zp, SparsePolyZp)
- T3: 无 ref 参数 (__make_zp) → 不变
- T4: 65 函数烟测 → 全部满足 HIR₁ 不变量
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, PairType, TupleType, UnknownType,
    TupleExpr, Var, Lit, FieldAccess,
    ReturnStmt, HIRFunc, HIRParam, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass, assert_hir1_invariant
from class_map import TRANSLATION_SCOPE

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"


def _load_hir0(func_name: str) -> HIRFunc:
    with open(AST_CACHE_DIR / f"{func_name}.json") as f:
        ast = json.load(f)
    return parse_pass(ast)


# ============================================================
# 单元测试
# ============================================================

def test_no_ref_params_unchanged():
    """T3: __make_zp 无 ref 参数 → Pass 2 后签名不变。"""
    hir0 = _load_hir0("__make_zp")
    hir1 = ref_elim_pass(hir0)
    assert len(hir1.params) == 2
    assert hir1.params[0].name == "val"
    assert not hir1.params[0].is_ref
    assert not hir1.params[0].is_output
    assert hir1.ret_ty == NamedType("Zp")  # 不变
    # body 也不变
    assert len(hir1.body) == len(hir0.body)
    assert_hir1_invariant(hir1)


def test_two_ref_void_to_pair():
    """T1: __upoly_divmod (void + 2 ref) → ret_ty = PairType(SparsePolyZp, SparsePolyZp)."""
    hir0 = _load_hir0("__upoly_divmod")
    # Pass 1 应该把 q, r 标为 is_ref
    assert hir0.params[0].is_ref and hir0.params[1].is_ref
    assert hir0.ret_ty == BaseType.UNIT

    hir1 = ref_elim_pass(hir0)
    assert len(hir1.params) == 4  # 参数数不变
    # 所有 params 的 is_ref/is_output 都应为 False
    assert all(not p.is_ref for p in hir1.params)
    assert all(not p.is_output for p in hir1.params)
    # ret_ty 应是 PairType(q.ty, r.ty) 即 PairType(SparsePolyZp, SparsePolyZp)
    assert isinstance(hir1.ret_ty, PairType), f"ret_ty should be PairType, got {hir1.ret_ty}"
    assert hir1.ret_ty.fst == NamedType("SparsePolyZp")
    assert hir1.ret_ty.snd == NamedType("SparsePolyZp")
    # body 末尾应追加 ReturnStmt((q, r))
    assert len(hir1.body) > 0
    last = hir1.body[-1]
    assert isinstance(last, ReturnStmt)
    # value 应是 TupleExpr 或 PairType 的等价 — 我们用 TupleExpr 表示 2 元
    assert isinstance(last.value, TupleExpr)
    assert len(last.value.elems) == 2
    assert isinstance(last.value.elems[0], Var) and last.value.elems[0].name == "q"
    assert isinstance(last.value.elems[1], Var) and last.value.elems[1].name == "r"
    assert_hir1_invariant(hir1)


def test_nonvoid_one_ref_to_pair():
    """T2: __upoly_make_monic (Zp + 1 ref `f`) → ret_ty = PairType(Zp, SparsePolyZp)."""
    hir0 = _load_hir0("__upoly_make_monic")
    # Pass 1 应把 f 标为 is_ref
    assert hir0.params[0].is_ref
    assert hir0.ret_ty == NamedType("Zp")

    hir1 = ref_elim_pass(hir0)
    assert len(hir1.params) == 1
    assert not hir1.params[0].is_ref
    # 返回类型应是 PairType(Zp, SparsePolyZp)
    assert isinstance(hir1.ret_ty, PairType)
    assert hir1.ret_ty.fst == NamedType("Zp")
    assert hir1.ret_ty.snd == NamedType("SparsePolyZp")
    # 原函数有多个 return；所有 return 值都应包装为 tuple
    from ir_types import IfStmt
    def count_returns(body):
        n = 0
        for s in body:
            if isinstance(s, ReturnStmt):
                assert isinstance(s.value, TupleExpr), f"Return value not TupleExpr: {s.value}"
                assert len(s.value.elems) == 2
                n += 1
            elif isinstance(s, IfStmt):
                n += count_returns(s.then_body)
                n += count_returns(s.else_body)
            elif hasattr(s, "body") and isinstance(s.body, list):
                n += count_returns(s.body)
        return n
    n_ret = count_returns(hir1.body)
    assert n_ret >= 2, f"expected ≥2 Return, got {n_ret}"
    assert_hir1_invariant(hir1)


def test_four_ref_mtshl_wmds():
    """T4: __mtshl_step_j 含 F: Array<MvPolyZp> [REF] + lc_tau: Array<...> [REF]
    → ret_ty 含多 ref 的 tuple（Wang 模块内最复杂的输出参数组合）。"""
    hir0 = _load_hir0("__mtshl_step_j")
    ref_params = [p for p in hir0.params if p.is_ref]
    assert len(ref_params) >= 1, f"expected ref params, got {[p.name for p in hir0.params]}"

    hir1 = ref_elim_pass(hir0)
    assert all(not p.is_ref for p in hir1.params)
    assert_hir1_invariant(hir1)


# ============================================================
# 全量烟测
# ============================================================

def test_smoke_all_65():
    """T-smoke: 65 函数全跑 Pass 2，全部应满足 HIR₁ 不变量。"""
    n_with_refs = 0
    n_without_refs = 0
    errors = []
    for func_name in sorted(TRANSLATION_SCOPE):
        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            continue
        try:
            hir0 = _load_hir0(func_name)
            has_ref = any(p.is_ref for p in hir0.params)
            hir1 = ref_elim_pass(hir0)
            assert_hir1_invariant(hir1)
            if has_ref:
                n_with_refs += 1
            else:
                n_without_refs += 1
        except Exception as e:
            errors.append((func_name, f"{type(e).__name__}: {e}"))

    print(f"\n  {n_with_refs} functions had ref params (transformed to tuple)")
    print(f"  {n_without_refs} functions had no ref params (passthrough)")
    print(f"  {len(errors)} errors")
    for fn, err in errors:
        print(f"    FAIL {fn}: {err}")
    assert not errors, f"{len(errors)} errors: {errors}"


if __name__ == "__main__":
    tests = [
        test_no_ref_params_unchanged,
        test_two_ref_void_to_pair,
        test_nonvoid_one_ref_to_pair,
        test_four_ref_mtshl_wmds,
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
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\n=== {passed}/{passed+failed} passed ===")
    sys.exit(0 if failed == 0 else 1)
