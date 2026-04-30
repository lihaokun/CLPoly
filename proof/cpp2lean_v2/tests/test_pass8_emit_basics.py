"""Pass 8 S1 单测：emit_type / emit_expr / emit_stmt 基础。"""

from __future__ import annotations
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, ArrayType, PairType, TupleType, OptionType,
    StdMapType, RefType, UnknownType,
    Var, Lit, BinOp, UnaryOp, CondExpr, Call, ArrayAccess, FieldAccess, Cast,
    TupleExpr, ArrayLit, UnknownExpr, UnresolvedOp,
    LetStmt, RequireStmt, HIRParam,
)
from pass8_codegen import (
    emit_type, emit_lit, emit_expr, emit_stmt, emit_param, emit_var_name,
    EmitCtx,
)


def _ctx() -> EmitCtx:
    return EmitCtx()


# ============================================================
# 类型 emit
# ============================================================

def test_emit_type_basetype():
    assert emit_type(BaseType.UINT64) == "UInt64"
    assert emit_type(BaseType.INT64) == "Int64"
    assert emit_type(BaseType.NAT) == "Nat"
    assert emit_type(BaseType.BOOL) == "Bool"
    assert emit_type(BaseType.UNIT) == "Unit"
    print("PASS test_emit_type_basetype")


def test_emit_type_named():
    assert emit_type(NamedType("Zp")) == "Zp"
    assert emit_type(NamedType("SparsePolyZp")) == "SparsePolyZp"
    # Lambda/LambdaRef 残留 → sorry
    assert "sorry" in emit_type(NamedType("Lambda"))
    assert "sorry" in emit_type(NamedType("LambdaRef"))
    print("PASS test_emit_type_named")


def test_emit_type_composite():
    # Array
    assert emit_type(ArrayType(BaseType.UINT64)) == "Array UInt64"
    assert emit_type(ArrayType(NamedType("Zp"))) == "Array Zp"
    # Pair
    assert emit_type(PairType(BaseType.UINT64, BaseType.INT64)) == "(UInt64 × Int64)"
    # Tuple
    assert emit_type(TupleType((BaseType.UINT64, BaseType.INT64, BaseType.BOOL))) \
        == "(UInt64 × Int64 × Bool)"
    # Option
    assert emit_type(OptionType(BaseType.NAT)) == "Option Nat"
    # StdMap
    assert emit_type(StdMapType(BaseType.UINT64, NamedType("Zp"))) == "StdMap UInt64 Zp"
    # 嵌套：Array (UInt64 × Int64)
    assert emit_type(ArrayType(PairType(BaseType.UINT64, BaseType.INT64))) \
        == "Array (UInt64 × Int64)"
    print("PASS test_emit_type_composite")


def test_emit_type_residual_sorry():
    # RefType 残留 → inner + 注释
    assert "/- ref residual -/" in emit_type(RefType(BaseType.UINT64))
    # UnknownType → sorry
    assert "sorry" in emit_type(UnknownType(""))
    assert "std::size_t" in emit_type(UnknownType("std::size_t"))
    print("PASS test_emit_type_residual_sorry")


# ============================================================
# 字面量 emit
# ============================================================

def test_emit_lit():
    # 整数字面量不带类型标注，让 Lean 推断
    assert emit_lit(Lit(42, BaseType.INT64)) == "42"
    assert emit_lit(Lit(0, BaseType.UINT64)) == "0"
    assert emit_lit(Lit(True, BaseType.BOOL)) == "true"
    assert emit_lit(Lit(False, BaseType.BOOL)) == "false"
    assert emit_lit(Lit(0, BaseType.UNIT)) == "()"
    print("PASS test_emit_lit")


# ============================================================
# 表达式 emit
# ============================================================

def test_emit_expr_var():
    # version=0 → 无后缀
    v0 = Var(name="m", version=0, ty=BaseType.UINT64)
    assert emit_expr(v0, _ctx()) == "m"
    # version=2 → m_2
    v2 = Var(name="m", version=2, ty=BaseType.UINT64)
    assert emit_expr(v2, _ctx()) == "m_2"
    # 关键字撞名 → «match»
    vk = Var(name="match", version=0)
    assert emit_expr(vk, _ctx()) == "«match»"
    print("PASS test_emit_expr_var")


def test_emit_expr_binop():
    a = Var(name="a", version=0, ty=BaseType.UINT64)
    b = Var(name="b", version=0, ty=BaseType.UINT64)
    e_add = BinOp("+", a, b, ty=BaseType.UINT64)
    assert emit_expr(e_add, _ctx()) == "(a + b)"
    e_shr = BinOp(">>", a, Lit(1, BaseType.UINT64), ty=BaseType.UINT64)
    assert emit_expr(e_shr, _ctx()) == "(a >>> 1)"
    e_eq = BinOp("==", a, Lit(0, BaseType.UINT64), ty=BaseType.BOOL)
    assert emit_expr(e_eq, _ctx()) == "(a == 0)"
    print("PASS test_emit_expr_binop")


def test_emit_expr_condexpr():
    cond = Var(name="flag", version=0, ty=BaseType.BOOL)
    a = Var(name="a", version=0)
    b = Var(name="b", version=0)
    e = CondExpr(cond=cond, then_e=a, else_e=b, ty=BaseType.UINT64)
    assert emit_expr(e, _ctx()) == "(if flag then a else b)"
    print("PASS test_emit_expr_condexpr")


def test_emit_expr_field_array_tuple():
    obj = Var(name="p", version=0, ty=NamedType("Pair"))
    fa = FieldAccess(obj=obj, field_name="fst", ty=BaseType.UINT64)
    assert emit_expr(fa, _ctx()) == "p.fst"

    arr = Var(name="xs", version=0)
    idx = Var(name="i", version=0, ty=BaseType.UINT64)
    aa = ArrayAccess(arr=arr, idx=idx, ty=BaseType.UINT64)
    # idx 强制 () 包裹（防裸字面量 0.toNat decimal 歧义）
    assert emit_expr(aa, _ctx()) == "(xs[(i).toNat]!)"

    te = TupleExpr(elems=[Var("a"), Var("b"), Var("c")])
    assert emit_expr(te, _ctx()) == "(a, b, c)"

    al = ArrayLit(elems=[Lit(1, BaseType.UINT64), Lit(2, BaseType.UINT64)],
                   elem_ty=BaseType.UINT64)
    assert emit_expr(al, _ctx()) == "#[1, 2]"
    print("PASS test_emit_expr_field_array_tuple")


def test_emit_expr_cast():
    a = Var(name="a", version=0, ty=BaseType.UINT64)
    # UINT64 → UINT32 截断
    c1 = Cast(expr=a, source_ty=BaseType.UINT64, target_ty=BaseType.UINT32,
                cast_kind="IntegralCast")
    assert emit_expr(c1, _ctx()) == "(a).toUInt32"
    # UINT64 → BOOL
    c2 = Cast(expr=a, source_ty=BaseType.UINT64, target_ty=BaseType.BOOL,
                cast_kind="IntegralToBoolean")
    assert emit_expr(c2, _ctx()) == "(a != 0)"
    # 同类型 NoOp
    c3 = Cast(expr=a, source_ty=BaseType.UINT64, target_ty=BaseType.UINT64,
                cast_kind="NoOp")
    assert emit_expr(c3, _ctx()) == "a"
    # NamedType 兜底
    c4 = Cast(expr=Lit(0, BaseType.UINT64), source_ty=BaseType.UINT64,
                target_ty=NamedType("Zp"), cast_kind="ConstructorConversion")
    assert "Zp" in emit_expr(c4, _ctx())
    print("PASS test_emit_expr_cast")


def test_emit_expr_call():
    # Call 已知 LEAN_BUILTINS（class_map 已注册）→ 直接输出名字
    c1 = Call(callee="Array.empty", args=[])
    assert emit_expr(c1, _ctx()) == "Array.empty"
    # Call _lambda_ → 加 _ir 后缀（lean_name 规则）
    c2 = Call(callee="_lambda_foo_filt1",
                args=[Var("xs"), Var("m")])
    assert emit_expr(c2, _ctx()) == "(_lambda_foo_filt1_ir xs m)"
    # Call _loop_ → 加 _ir 后缀
    c3 = Call(callee="_loop_bar_0", args=[Var("k")])
    assert emit_expr(c3, _ctx()) == "(_loop_bar_0_ir k)"
    # 已带 _ir 后缀
    c4 = Call(callee="__make_zp_ir", args=[Var("p")])
    assert emit_expr(c4, _ctx()) == "(__make_zp_ir p)"
    # 未知 callee 但是合法 ident → 直接输出
    c5 = Call(callee="some_var", args=[])
    assert emit_expr(c5, _ctx()) == "some_var"
    print("PASS test_emit_expr_call")


def test_emit_expr_unknown():
    # UnknownExpr / UnresolvedOp → sorry
    assert "sorry" in emit_expr(UnknownExpr(kind="weird"), _ctx())
    assert "sorry" in emit_expr(UnresolvedOp(op_name="operator?"), _ctx())
    print("PASS test_emit_expr_unknown")


# ============================================================
# 语句 emit
# ============================================================

def test_emit_stmt_let():
    ctx = EmitCtx(indent=1)
    s = LetStmt(var=Var("x", version=2, ty=BaseType.UINT64), ty=BaseType.UINT64,
                  value=Lit(42, BaseType.UINT64))
    out = emit_stmt(s, ctx)
    assert out == "  let x_2 : UInt64 := 42"
    # sorry type → 省略标注
    s2 = LetStmt(var=Var("y"), ty=UnknownType(""), value=Var("rhs"))
    out2 = emit_stmt(s2, EmitCtx(indent=0))
    assert out2 == "let y := rhs"
    print("PASS test_emit_stmt_let")


def test_emit_stmt_require():
    ctx = EmitCtx(indent=0)
    s = RequireStmt(cond=BinOp("!=", Var("p"), Lit(0, BaseType.UINT64),
                                 ty=BaseType.BOOL),
                     name="hp", source="div_by_zero")
    out = emit_stmt(s, ctx)
    assert out.startswith("-- require")
    assert "hp" in out
    print("PASS test_emit_stmt_require")


# ============================================================
# Param emit
# ============================================================

def test_emit_param():
    p = HIRParam(name="m", ty=BaseType.UINT64, is_ref=False,
                  is_const_ref=False, is_output=False)
    assert emit_param(p) == "(m : UInt64)"
    # 关键字撞名
    pk = HIRParam(name="match", ty=BaseType.BOOL, is_ref=False,
                    is_const_ref=False, is_output=False)
    assert emit_param(pk) == "(«match» : Bool)"
    print("PASS test_emit_param")


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_emit_type_basetype,
        test_emit_type_named,
        test_emit_type_composite,
        test_emit_type_residual_sorry,
        test_emit_lit,
        test_emit_expr_var,
        test_emit_expr_binop,
        test_emit_expr_condexpr,
        test_emit_expr_field_array_tuple,
        test_emit_expr_cast,
        test_emit_expr_call,
        test_emit_expr_unknown,
        test_emit_stmt_let,
        test_emit_stmt_require,
        test_emit_param,
    ]
    passed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n=== {passed}/{len(tests)} passed ===")
    sys.exit(0 if passed == len(tests) else 1)
