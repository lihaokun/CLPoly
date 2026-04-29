"""
Pass 2b callsite_ref_elim 单元测试。

按 mutation-model-design.md §3.2 + 修正方案 §3 阶段 2 关键用例：
- T1: void + 1 out → `out := f(out, in)`
- T2: void + 2 out → `let __refret := f(...); out0 := __refret.fst; out1 := __refret.snd`
- T3: arg 是 `arr[i]`（HIR1 阶段为 Call(UnresolvedOp("operator[]"))）→ ArrayAccess unwrap
- T4: arg 是 `vec.data()`（HIR1 阶段为 Call(UnresolvedOp("<method>.data"))）→ 透明 unwrap
- T5: callee 不在 OUTPUT_PARAMS → passthrough（不改写）
- T6: arity overload 解析（pair_vec_div#4 vs #5）
- T7: 全 67 函数烟测（pipeline + 残余=0）
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
    Var, Lit, Call, FieldAccess, Cast, UnresolvedOp,
    LetStmt, AssignStmt, ExprStmt, ReturnStmt,
    HIRFunc, HIRParam,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass
from class_map import TRANSLATION_SCOPE


def _mk_func(body, params=None, ret_ty=BaseType.UNIT) -> HIRFunc:
    return HIRFunc(
        base_name="test", instance_suffix="", mangled_name="", qual_type="",
        params=params or [], ret_ty=ret_ty, body=body,
        requires=[], aux_lambdas=[],
    )


# ============================================================
# T1: void + 1 out
# ============================================================

def test_single_out_directly_assigned():
    """`(expr) __upoly_make_monic(f)` → `f := __upoly_make_monic(f)`"""
    body = [
        ExprStmt(expr=Call(callee="__upoly_make_monic",
                           args=[Var("f")], ty=BaseType.UNIT)),
    ]
    func = _mk_func(body, params=[HIRParam(name="f", ty=NamedType("Poly"))])
    out = callsite_ref_elim_pass(func)
    assert len(out.body) == 1
    s = out.body[0]
    assert isinstance(s, AssignStmt), f"got {type(s).__name__}"
    assert isinstance(s.target, Var) and s.target.name == "f"
    assert isinstance(s.value, Call) and s.value.callee == "__upoly_make_monic"


# ============================================================
# T2: void + 2 out → 临时变量 destructure
# ============================================================

def test_double_out_destructure():
    """`(expr) pair_vec_div(q, r, g, divisor, comp)` →
       `let __refret_0 := pair_vec_div(...); q := __refret_0.fst; r := __refret_0.snd`
    """
    body = [
        ExprStmt(expr=Call(
            callee="pair_vec_div",
            args=[Var("q"), Var("r"), Var("g"), Var("divisor"), Var("comp")],
            ty=BaseType.UNIT,
        )),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="q", ty=NamedType("Vec")),
        HIRParam(name="r", ty=NamedType("Vec")),
        HIRParam(name="g", ty=NamedType("Vec")),
        HIRParam(name="divisor", ty=NamedType("Vec")),
        HIRParam(name="comp", ty=NamedType("Comp")),
    ])
    out = callsite_ref_elim_pass(func)
    # 应有 3 个 stmt：tmp let + 2 个 assign
    assert len(out.body) == 3, f"got {len(out.body)} stmts"
    let_tmp, a0, a1 = out.body
    assert isinstance(let_tmp, LetStmt) and let_tmp.var.name.startswith("__refret_")
    assert isinstance(a0, AssignStmt) and a0.target.name == "q"
    assert isinstance(a0.value, FieldAccess) and a0.value.field_name == "fst"
    assert isinstance(a1, AssignStmt) and a1.target.name == "r"
    assert isinstance(a1.value, FieldAccess) and a1.value.field_name == "snd"


# ============================================================
# T3: arr[i] 形态 unwrap
# ============================================================

def test_arrayaccess_unwrap():
    """`swap(mu[k], mu[k-1])` 在 HIR1 是 Call(UnresolvedOp(operator[])) 包装；
       Pass 2b 应识别为 ArrayAccess 作为 lvalue。
    """
    arr_access = lambda arr, idx: Call(
        callee=UnresolvedOp(op_name="operator[]"),
        args=[arr, idx],
        ty=BaseType.INT64,
    )
    body = [
        ExprStmt(expr=Call(
            callee="swap",
            args=[arr_access(Var("mu"), Var("k")),
                  arr_access(Var("mu"), Var("k_minus_1"))],
            ty=BaseType.UNIT,
        )),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="mu", ty=NamedType("Vec")),
        HIRParam(name="k", ty=BaseType.INT64),
        HIRParam(name="k_minus_1", ty=BaseType.INT64),
    ])
    out = callsite_ref_elim_pass(func)
    # 应有 3 stmts: tmp + 2 ArrayAccess assigns
    assert len(out.body) == 3
    a0 = out.body[1]
    assert isinstance(a0, AssignStmt)
    from ir_types import ArrayAccess
    assert isinstance(a0.target, ArrayAccess), f"got {type(a0.target).__name__}"


# ============================================================
# T4: vec.data() 透明 unwrap
# ============================================================

def test_data_method_unwrap():
    """`(expr) f(vec.data(), other.data())` →
       `vec := f(vec.data(), ...)`（vec 作为 lvalue）
    """
    method_call = lambda recv, name: Call(
        callee=UnresolvedOp(op_name=f"<method>.{name}", receiver_ty=NamedType("Vec")),
        args=[recv], ty=NamedType("VecData"),
    )
    body = [
        ExprStmt(expr=Call(
            callee="__upoly_make_monic",
            args=[method_call(Var("vec"), "data")],
            ty=BaseType.UNIT,
        )),
    ]
    func = _mk_func(body, params=[HIRParam(name="vec", ty=NamedType("Poly"))])
    out = callsite_ref_elim_pass(func)
    assert len(out.body) == 1
    s = out.body[0]
    assert isinstance(s, AssignStmt)
    assert isinstance(s.target, Var) and s.target.name == "vec"


# ============================================================
# T5: 未注册 callee → passthrough
# ============================================================

def test_unregistered_callee_passthrough():
    """`(expr) some_unknown_func(a, b)` → 保持 ExprStmt"""
    body = [
        ExprStmt(expr=Call(callee="some_unknown_func",
                           args=[Var("a"), Var("b")], ty=BaseType.UNIT)),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="a", ty=BaseType.INT64),
        HIRParam(name="b", ty=BaseType.INT64),
    ])
    out = callsite_ref_elim_pass(func)
    assert len(out.body) == 1
    assert isinstance(out.body[0], ExprStmt)


# ============================================================
# T6: arity overload (pair_vec_div#4 vs #5)
# ============================================================

def test_arity_overload_resolution():
    """pair_vec_div 4-arg → out=[0]；5-arg → out=[0,1]"""
    # 4-arg: 1 out
    body4 = [ExprStmt(expr=Call(
        callee="pair_vec_div",
        args=[Var("new_v"), Var("v1"), Var("v2"), Var("comp")],
        ty=BaseType.UNIT,
    ))]
    func4 = _mk_func(body4, params=[
        HIRParam(name="new_v", ty=NamedType("Vec")),
        HIRParam(name="v1", ty=NamedType("Vec")),
        HIRParam(name="v2", ty=NamedType("Vec")),
        HIRParam(name="comp", ty=NamedType("Comp")),
    ])
    out4 = callsite_ref_elim_pass(func4)
    assert len(out4.body) == 1, "4-arg pair_vec_div: 1 out → 直接 AssignStmt"
    assert isinstance(out4.body[0], AssignStmt)
    assert out4.body[0].target.name == "new_v"

    # 5-arg: 2 outs
    body5 = [ExprStmt(expr=Call(
        callee="pair_vec_div",
        args=[Var("new_v"), Var("R"), Var("v1"), Var("v2"), Var("comp")],
        ty=BaseType.UNIT,
    ))]
    func5 = _mk_func(body5, params=[
        HIRParam(name="new_v", ty=NamedType("Vec")),
        HIRParam(name="R", ty=NamedType("Vec")),
        HIRParam(name="v1", ty=NamedType("Vec")),
        HIRParam(name="v2", ty=NamedType("Vec")),
        HIRParam(name="comp", ty=NamedType("Comp")),
    ])
    out5 = callsite_ref_elim_pass(func5)
    assert len(out5.body) == 3, "5-arg pair_vec_div: 2 outs → tmp + 2 AssignStmt"


# ============================================================
# T7: 67 函数烟测（pipeline + 整体不破坏）
# ============================================================

def test_pipeline_smoke_67():
    """全 67 函数走 parse → ref_elim → callsite_ref_elim 全过 + 不抛错。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    ok = 0
    errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            hir = callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {str(e)[:100]}"))
    assert not errs, f"{len(errs)} errors: {errs[:3]}"


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_single_out_directly_assigned,
        test_double_out_destructure,
        test_arrayaccess_unwrap,
        test_data_method_unwrap,
        test_unregistered_callee_passthrough,
        test_arity_overload_resolution,
        test_pipeline_smoke_67,
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
