"""
Pass 3b lambda_ref_elim 单元测试。

按 mutation-model-design.md §5 G3 + 修正方案
docs/fixes/cpp2lean-v2-lambda-by-ref-capture.md 第三部分关键用例：

- T1: 单 modified capture（`[&x]` lambda 修改 x）
- T2: 双 modified captures（`[&M, &U]` 修改两个）
- T3: mixed modified + read-only capture
- T4: 0 modified captures（read-only `[&]`）
- T5: 67 函数烟测（含 6 个 by-ref capture 函数）
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType, ArrayType, PairType,
    Var, Lit, Call, FieldAccess, UnresolvedOp,
    LetStmt, AssignStmt, ExprStmt, ReturnStmt, IfStmt,
    HIRFunc, HIRParam,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass3b_lambda_ref_elim import (
    lambda_ref_elim_pass, _build_alias_map, _extract_cap_info,
)
from class_map import TRANSLATION_SCOPE


def _mk_lifted(name: str, params: list[HIRParam], body=None,
               ret_ty=BaseType.UNIT, n_caps: int = None) -> HIRFunc:
    """构造测试用 lifted lambda。
    n_caps 默认 = 所有 leading is_ref/is_const_ref 参数的个数（旧测试假设）。
    """
    if n_caps is None:
        n_caps = 0
        for p in params:
            if p.is_ref or p.is_const_ref:
                n_caps += 1
            else:
                break
    return HIRFunc(
        base_name=name, instance_suffix="", mangled_name="",
        qual_type=f"lambda lifted (test) | n_caps={n_caps} | modified_captures=[]",
        params=params, ret_ty=ret_ty,
        body=body or [],
        requires=[], aux_lambdas=[],
    )


def _mk_outer(body, params=None, ret_ty=BaseType.UNIT,
              aux_lambdas=None) -> HIRFunc:
    return HIRFunc(
        base_name="test_outer", instance_suffix="", mangled_name="", qual_type="",
        params=params or [], ret_ty=ret_ty, body=body,
        requires=[], aux_lambdas=aux_lambdas or [],
    )


# ============================================================
# T1: 单 modified capture
# ============================================================

def test_single_modified_capture():
    """
    [&x] lambda body 改 x → lifted: params=[x[REF]], 1 modified.
    Outer body 中 `lam(arg)` (operator() 形态) → `x := lam(x, arg)`
    """
    lifted = _mk_lifted(
        "_lambda_test_1",
        params=[
            HIRParam(name="x", ty=BaseType.INT64, is_ref=True),
            HIRParam(name="i", ty=BaseType.INT64),
        ],
        body=[AssignStmt(target=Var("x"),
                         value=Lit(42, ty=BaseType.INT64))],
    )
    outer_body = [
        LetStmt(var=Var("alias"), ty=NamedType("LambdaRef"),
                value=Var("_lambda_test_1")),
        # 调用：Call(UnresolvedOp("operator()"), [Var("alias"), Var("i")])
        ExprStmt(expr=Call(
            callee=UnresolvedOp(op_name="operator()"),
            args=[Var("alias"), Var("i")],
            ty=BaseType.UNIT,
        )),
    ]
    outer = _mk_outer(outer_body,
                      params=[HIRParam(name="x", ty=BaseType.INT64),
                              HIRParam(name="i", ty=BaseType.INT64)],
                      aux_lambdas=[lifted])
    out = lambda_ref_elim_pass(outer)
    # 应有 LetStmt(alias) + AssignStmt(x := lam(x, i))
    assert len(out.body) == 2
    assign = out.body[1]
    assert isinstance(assign, AssignStmt), f"got {type(assign).__name__}"
    assert isinstance(assign.target, Var) and assign.target.name == "x"
    assert isinstance(assign.value, Call) and assign.value.callee == "_lambda_test_1"
    assert len(assign.value.args) == 2  # cap x + arg i
    assert assign.value.args[0].name == "x"
    assert assign.value.args[1].name == "i"


# ============================================================
# T2: 双 modified captures
# ============================================================

def test_double_modified_captures():
    """
    [&M, &U] lambda → lifted: params=[M[REF], U[REF], i, j], 2 modified.
    Outer 调用：tmp + 2 destructure
    """
    lifted = _mk_lifted(
        "_lambda_test_2",
        params=[
            HIRParam(name="M", ty=NamedType("Matrix"), is_ref=True),
            HIRParam(name="U", ty=NamedType("Matrix"), is_ref=True),
            HIRParam(name="i", ty=BaseType.INT64),
            HIRParam(name="j", ty=BaseType.INT64),
        ],
    )
    outer_body = [
        LetStmt(var=Var("row_swap"), ty=NamedType("LambdaRef"),
                value=Var("_lambda_test_2")),
        ExprStmt(expr=Call(
            callee=UnresolvedOp(op_name="operator()"),
            args=[Var("row_swap"), Var("i"), Var("j")],
            ty=BaseType.UNIT,
        )),
    ]
    outer = _mk_outer(outer_body,
                      params=[HIRParam(name="M", ty=NamedType("Matrix")),
                              HIRParam(name="U", ty=NamedType("Matrix")),
                              HIRParam(name="i", ty=BaseType.INT64),
                              HIRParam(name="j", ty=BaseType.INT64)],
                      aux_lambdas=[lifted])
    out = lambda_ref_elim_pass(outer)
    # 应有 LetStmt(row_swap) + LetStmt(__refret_lam_0) + 2 AssignStmts
    assert len(out.body) == 4, f"got {len(out.body)} stmts"
    let_ref = out.body[1]
    assert isinstance(let_ref, LetStmt)
    assert let_ref.var.name.startswith("__refret_lam_")
    assert isinstance(let_ref.value, Call)
    assert len(let_ref.value.args) == 4  # M, U, i, j
    assert [a.name for a in let_ref.value.args] == ["M", "U", "i", "j"]
    # M := tmp.fst, U := tmp.snd
    a0, a1 = out.body[2], out.body[3]
    assert isinstance(a0, AssignStmt) and a0.target.name == "M"
    assert isinstance(a0.value, FieldAccess) and a0.value.field_name == "fst"
    assert isinstance(a1, AssignStmt) and a1.target.name == "U"
    assert isinstance(a1.value, FieldAccess) and a1.value.field_name == "snd"


# ============================================================
# T3: mixed modified + read-only capture
# ============================================================

def test_mixed_modified_and_const_caps():
    """
    [&out, &in_] (in_ const) 仅 out 修改 → lifted: params=[out[REF], in_[CONST-REF], k]
    Outer 调用：out := lam(out, in_, k)（仅 out destructure）
    """
    lifted = _mk_lifted(
        "_lambda_test_3",
        params=[
            HIRParam(name="out", ty=BaseType.INT64, is_ref=True),
            HIRParam(name="in_", ty=BaseType.INT64, is_const_ref=True),
            HIRParam(name="k", ty=BaseType.INT64),
        ],
    )
    outer_body = [
        LetStmt(var=Var("compute"), ty=NamedType("LambdaRef"),
                value=Var("_lambda_test_3")),
        ExprStmt(expr=Call(
            callee=UnresolvedOp(op_name="operator()"),
            args=[Var("compute"), Var("k")],
            ty=BaseType.UNIT,
        )),
    ]
    outer = _mk_outer(outer_body,
                      params=[HIRParam(name="out", ty=BaseType.INT64),
                              HIRParam(name="in_", ty=BaseType.INT64),
                              HIRParam(name="k", ty=BaseType.INT64)],
                      aux_lambdas=[lifted])
    out = lambda_ref_elim_pass(outer)
    # 单 modified → 直接 AssignStmt（无临时变量）
    assert len(out.body) == 2
    assign = out.body[1]
    assert isinstance(assign, AssignStmt) and assign.target.name == "out"
    # cap 全部 prepend：[out, in_, k]
    assert [a.name for a in assign.value.args] == ["out", "in_", "k"]


# ============================================================
# T4: 0 modified captures（仅 read-only）
# ============================================================

def test_no_modified_captures_only_prepend():
    """
    [&n, &dot] 全 const-ref → lifted: params=[n[CONST-REF], dot[CONST-REF], a, b]
    Outer 调用：保留 ExprStmt（不 destructure）但 prepend caps
    """
    lifted = _mk_lifted(
        "_lambda_test_4",
        params=[
            HIRParam(name="n", ty=BaseType.INT64, is_const_ref=True),
            HIRParam(name="dot", ty=NamedType("Func"), is_const_ref=True),
            HIRParam(name="a", ty=BaseType.INT64),
            HIRParam(name="b", ty=BaseType.INT64),
        ],
    )
    outer_body = [
        LetStmt(var=Var("cmp"), ty=NamedType("LambdaRef"),
                value=Var("_lambda_test_4")),
        ExprStmt(expr=Call(
            callee=UnresolvedOp(op_name="operator()"),
            args=[Var("cmp"), Var("a"), Var("b")],
            ty=BaseType.UNIT,
        )),
    ]
    outer = _mk_outer(outer_body,
                      params=[HIRParam(name="n", ty=BaseType.INT64),
                              HIRParam(name="dot", ty=NamedType("Func")),
                              HIRParam(name="a", ty=BaseType.INT64),
                              HIRParam(name="b", ty=BaseType.INT64)],
                      aux_lambdas=[lifted])
    out = lambda_ref_elim_pass(outer)
    # ExprStmt(Call(_lambda_test_4, [n, dot, a, b])) — 4 args, 仍是 ExprStmt
    assert len(out.body) == 2
    es = out.body[1]
    assert isinstance(es, ExprStmt) and isinstance(es.expr, Call)
    assert es.expr.callee == "_lambda_test_4"
    assert [a.name for a in es.expr.args] == ["n", "dot", "a", "b"]


# ============================================================
# T5: 67 函数烟测
# ============================================================

def test_pipeline_smoke_67():
    """全 67 函数走 parse → ref_elim → 2b → lift → 3b 全过 + 不抛错。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    ok = 0
    errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            hir = lambda_ref_elim_pass(lambda_lift_pass(
                callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast)))))
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {str(e)[:100]}"))
    assert not errs, f"{len(errs)} errors: {errs[:3]}"


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_single_modified_capture,
        test_double_modified_captures,
        test_mixed_modified_and_const_caps,
        test_no_modified_captures_only_prepend,
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
