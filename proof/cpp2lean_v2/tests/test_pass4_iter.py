"""
Pass 4 iter_recognize 单元测试。

覆盖：
- T1: structured-binding desugar（decomposition 非空 → body 头部 LetStmt）
- T2: A 形态 compact-erase（LetStmt 前缀）— __upoly_mod_coeff
- T3: A 形态 compact-erase（ExprStmt(operator=) 前缀，变量复用）— __hensel_step 第 2 处
- T4: B 形态 classic-ForStmt — __extract_monomial_content body[5]
- T5: B 形态 classic-WhileStmt — __extract_monomial_content 内嵌
- T6: passthrough（无 filter-loop 且无 decomposition） — __make_zp
- T7: assert_hir3_invariant 拒绝残留 decomposition
- T8: assert_hir3_invariant 拒绝未识别的 A 形态
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType, PairType,
    Var, Lit, Call, UnresolvedOp, FieldAccess,
    LetStmt, AssignStmt, ExprStmt, IfStmt, WhileStmt, ForStmt, DoWhileStmt,
    BlockStmt, RangeForStmt, ReturnStmt,
    HIRFunc, HIRParam, LambdaExpr, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import (
    iter_recognize_pass, assert_hir3_invariant,
    _desugar_decomposition, _match_filter_loop_A, _match_filter_loop_B,
    _stmt_rhs,
)

AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"


def _load_hir2(func_name: str) -> HIRFunc:
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    with open(AST / f"{func_name}.json") as f:
        ast = json.load(f)
    return lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))))


def _count_filter_calls(stmts):
    """递归统计 Call("Array.filter", ...)。"""
    from ir_types import BlockStmt, DoWhileStmt
    n = 0
    for s in stmts:
        if isinstance(s, AssignStmt) and isinstance(s.value, Call):
            if isinstance(s.value.callee, str) and s.value.callee == "Array.filter":
                n += 1
        if isinstance(s, RangeForStmt): n += _count_filter_calls(s.body)
        elif isinstance(s, IfStmt):
            n += _count_filter_calls(s.then_body) + _count_filter_calls(s.else_body)
        elif isinstance(s, ForStmt):
            n += _count_filter_calls(s.init) + _count_filter_calls(s.step) + _count_filter_calls(s.body)
        elif isinstance(s, (WhileStmt, DoWhileStmt)): n += _count_filter_calls(s.body)
        elif isinstance(s, BlockStmt): n += _count_filter_calls(s.stmts)
    return n


def _count_decomposition(stmts):
    n = 0
    for s in stmts:
        if isinstance(s, RangeForStmt):
            if s.decomposition: n += 1
            n += _count_decomposition(s.body)
        elif isinstance(s, IfStmt):
            n += _count_decomposition(s.then_body) + _count_decomposition(s.else_body)
        elif isinstance(s, ForStmt):
            n += _count_decomposition(s.init) + _count_decomposition(s.step) + _count_decomposition(s.body)
        elif isinstance(s, WhileStmt): n += _count_decomposition(s.body)
    return n


# ============================================================
# T1: structured-binding desugar
# ============================================================

def test_structured_binding_desugar():
    """T1: RangeForStmt.decomposition = [k, v] → body 头部有 LetStmt(k, _x.fst), LetStmt(v, _x.snd)。"""
    # 人工构造：无容器，直接 HIRFunc
    x = Var(name="__x", version=0, ty=PairType(BaseType.INT32, BaseType.INT32))
    rf = RangeForStmt(
        var=x,
        var_ty=PairType(BaseType.INT32, BaseType.INT32),
        container=Var(name="mymap", version=0, ty=UnknownType("")),
        body=[ReturnStmt(value=Var(name="k", version=0))],
        decomposition=[Var(name="k", version=0), Var(name="v", version=0)],
    )
    new_rf = _desugar_decomposition(rf)
    assert new_rf.decomposition == [], f"decomposition should be empty, got {new_rf.decomposition}"
    assert len(new_rf.body) == 3, f"expected 3 stmts (2 let + 1 return), got {len(new_rf.body)}"
    l0, l1, _ = new_rf.body
    assert isinstance(l0, LetStmt) and l0.var.name == "k"
    assert isinstance(l0.value, FieldAccess) and l0.value.field_name == "fst"
    assert isinstance(l1, LetStmt) and l1.var.name == "v"
    assert isinstance(l1.value, FieldAccess) and l1.value.field_name == "snd"


# ============================================================
# T2-T5: 识别 6 处 filter-loop
# ============================================================

def test_compact_erase_pure_body():
    """T2: __hensel_step — 2 处 pure compact-erase（1 LetStmt 前缀 + 1 ExprStmt(operator=) 前缀）。

    fdiv_r 在独立前置 for-range 里，不在 compact-erase body 内，所以是纯筛选。
    """
    hir3 = iter_recognize_pass(_load_hir2("__hensel_step"))
    assert_hir3_invariant(hir3)
    n = _count_filter_calls(hir3.body)
    assert n == 2, f"expected 2 Array.filter for __hensel_step, got {n}"


def _count_filtermap_calls(stmts):
    """统计 Call("Array.filterMap", ...) (CF-1 mutate-then-filter)。"""
    n = 0
    for s in stmts:
        if isinstance(s, AssignStmt) and isinstance(s.value, Call):
            if isinstance(s.value.callee, str) and s.value.callee == "Array.filterMap":
                n += 1
        if isinstance(s, RangeForStmt): n += _count_filtermap_calls(s.body)
        elif isinstance(s, IfStmt):
            n += _count_filtermap_calls(s.then_body) + _count_filtermap_calls(s.else_body)
        elif isinstance(s, ForStmt):
            n += _count_filtermap_calls(s.init) + _count_filtermap_calls(s.step) + _count_filtermap_calls(s.body)
        elif isinstance(s, (WhileStmt, DoWhileStmt)): n += _count_filtermap_calls(s.body)
        elif isinstance(s, BlockStmt): n += _count_filtermap_calls(s.stmts)
    return n


def test_upoly_mod_coeff_mutate_filter():
    """T3 (CF-1 阶段 1): __upoly_mod_coeff — body 头部有 fdiv_r(it->second,...) mutator
    + 后续 if-write-advance；Pass 4 应识别为 mutate-then-filter，emit Array.filterMap。"""
    hir3 = iter_recognize_pass(_load_hir2("__upoly_mod_coeff"))
    assert_hir3_invariant(hir3)
    n_filter = _count_filter_calls(hir3.body)
    n_filtermap = _count_filtermap_calls(hir3.body)
    assert n_filter == 0, f"expected 0 Array.filter, got {n_filter}"
    assert n_filtermap >= 1, f"expected >=1 Array.filterMap (mutate-then-filter), got {n_filtermap}"


def test_hensel_step_linear_mutate_filter():
    """T4 (CF-1 阶段 1): __hensel_step_linear — body 头部 fdiv_q + fdiv_r mutators
    + if-write-advance；Pass 4 应识别 mutate-then-filter。"""
    hir3 = iter_recognize_pass(_load_hir2("__hensel_step_linear"))
    assert_hir3_invariant(hir3)
    n_filtermap = _count_filtermap_calls(hir3.body)
    assert n_filtermap >= 1, f"expected >=1 Array.filterMap, got {n_filtermap}"


def test_classic_both_containers_with_pred_inversion():
    """T5: __extract_monomial_content — B-ForStmt + B-WhileStmt 两处，pred 都应反转
    （C++ `if (EraseCond) erase` 保留 `!EraseCond` 的元素）。"""
    hir3 = iter_recognize_pass(_load_hir2("__extract_monomial_content"))
    assert_hir3_invariant(hir3)
    n = _count_filter_calls(hir3.body)
    assert n == 2, f"expected 2 Array.filter (B-For + B-While), got {n}"

    # 深扫两个 filter 的 pred，验证都是 UnaryOp("!")
    preds = []
    # P7-8 修复后：Array.filter 的 lambda 已被 lifted 到 aux_lambdas，
    # Call args[1] 是 Var(lifted_name) 不是 LambdaExpr。从 aux_lambdas 找
    # 对应 lifted lambda 的 body 提取 pred。
    from ir_types import AssignStmt, Call, Var, LambdaExpr, ReturnStmt
    lifted_names: list[str] = []
    def walk(stmts):
        for s in stmts:
            if isinstance(s, AssignStmt) and isinstance(s.value, Call):
                if isinstance(s.value.callee, str) and s.value.callee == "Array.filter":
                    arg = s.value.args[1]
                    if isinstance(arg, Var):
                        lifted_names.append(arg.name)
            if isinstance(s, RangeForStmt): walk(s.body)
            elif isinstance(s, IfStmt): walk(s.then_body); walk(s.else_body)
            elif isinstance(s, ForStmt): walk(s.body)
            elif isinstance(s, WhileStmt): walk(s.body)
            elif hasattr(s, 'stmts'): walk(s.stmts)
    walk(hir3.body)
    # 从 aux_lambdas 找对应 lifted func，提取 pred
    aux_by_name = {a.base_name: a for a in hir3.aux_lambdas}
    from ir_types import UnaryOp
    for i, name in enumerate(lifted_names):
        aux = aux_by_name.get(name)
        assert aux is not None, f"lifted {name} not in aux_lambdas"
        # body 末尾的 ReturnStmt.value 即 pred
        ret = aux.body[-1]
        assert isinstance(ret, ReturnStmt) and ret.value is not None
        p = ret.value
        assert isinstance(p, UnaryOp) and p.op == "!", \
            f"pred[{i}] should be UnaryOp('!', ...) for form B inversion, got {type(p).__name__}"
    assert len(lifted_names) == 2, f"expected 2 Array.filter calls, got {len(lifted_names)}"


# ============================================================
# T7: passthrough
# ============================================================

def test_make_zp_passthrough():
    """T7: __make_zp 无 filter-loop 无 decomposition → HIR₃ 与 HIR₂ 结构一致。"""
    hir2 = _load_hir2("__make_zp")
    hir3 = iter_recognize_pass(hir2)
    assert_hir3_invariant(hir3)
    assert _count_filter_calls(hir3.body) == 0
    assert _count_decomposition(hir3.body) == 0


# ============================================================
# T8: invariant 违反检测
# ============================================================

def test_invariant_rejects_leftover_decomposition():
    """T8a: 残留 decomposition → assert_hir3_invariant 抛 TranslationError。"""
    bad_rf = RangeForStmt(
        var=Var(name="__x", version=0),
        var_ty=PairType(BaseType.INT32, BaseType.INT32),
        container=Var(name="m", version=0, ty=UnknownType("")),
        body=[],
        decomposition=[Var(name="k", version=0)],  # 故意不空
    )
    func = HIRFunc(
        base_name="_test", instance_suffix="", mangled_name="", qual_type="",
        params=[], ret_ty=BaseType.UNIT, body=[bad_rf], requires=[], aux_lambdas=[],
    )
    try:
        assert_hir3_invariant(func)
        assert False, "should have raised"
    except TranslationError as e:
        assert "decomposition" in e.reason


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_structured_binding_desugar,
        test_compact_erase_pure_body,
        test_upoly_mod_coeff_mutate_filter,
        test_hensel_step_linear_mutate_filter,
        test_classic_both_containers_with_pred_inversion,
        test_make_zp_passthrough,
        test_invariant_rejects_leftover_decomposition,
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
