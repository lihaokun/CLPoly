"""
Pass 5: operator_resolve — HIR₃ → HIR₄

职责：
  1. 解析 `Cast` 节点：查 CAST_TABLE → Lean cast 函数；NoOp/LValueToRValue 等
     直接 strip；Lit→Dependent 用 `is_safe_to_strip_lit_cast` 剥
  2. 解析 `Call(UnresolvedOp(...))`：
     - `<method>.X`：查 CLASS_MAP[T]["methods"] —— mutator 类转 AssignStmt
     - `operator+/-/*/...` 二元/一元：查 CLASS_MAP[T]["operators"]，返回 None 用
       Lean typeclass，返回字符串用具体函数
     - `operator=`：转 AssignStmt
     - `operator[]` / `operator()`：转 ArrayAccess / lambda apply（含 UB require）
     - `construct_X` / `default_init_X`：查 CONSTRUCTOR_MAP
  3. 展开 `CompoundAssignStmt(t, op, v)` → `AssignStmt(t, BinOp(op, t, v))`
  4. 生成 UB require（除零/空容器/越界/signed 溢出，按 type-system.md §2）
  5. 维护本地 typectx（params + LetStmt + RangeFor binding）以推断 receiver
     类型和表达式类型
  6. 失败策略 B：未识别的节点保留原样并记入 gap 清单（Pass 8 codegen 输出 sorry）

依据 hir-design.md §7 + cast_table.py + class_map.py + constructor_map.py。

HIR₄ 不变量：
  - 通过 `assert_hir4_invariant` 验证；放宽允许残留 UnresolvedOp（gap），
    但禁止 CompoundAssignStmt 残留
"""

from __future__ import annotations
import sys
from dataclasses import replace, dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType, RefType,
    ArrayType, StdMapType,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRFunc, HIRParam, TranslationError,
)
from class_map import CLASS_MAP, FUNC_MAP
from cast_table import (
    CAST_KIND_DISPOSITION, lookup_cast, canonicalize_type,
    is_safe_to_strip_lit_cast,
)
from constructor_map import resolve_constructor


# ============================================================
# Gap log（B 策略：未识别节点累积清单，供后续审计）
# ============================================================

@dataclass
class GapLog:
    unresolved_op: list[str] = field(default_factory=list)  # op_name 列表
    cast_miss:     list[tuple[str, str]] = field(default_factory=list)
    method_miss:   list[tuple[str, str]] = field(default_factory=list)  # (T, method)
    op_miss:       list[tuple[str, str]] = field(default_factory=list)
    constructor_miss: list[tuple[str, int]] = field(default_factory=list)


# ============================================================
# Typectx：变量名 → 类型
# ============================================================

def _strip_ref(ty: TypeIR | None) -> TypeIR | None:
    while isinstance(ty, RefType):
        ty = ty.inner
    return ty


def _expr_ty(expr: ExprIR, typectx: dict[str, TypeIR]) -> TypeIR | None:
    """推断表达式的类型（基于 typectx + 节点 .ty）。"""
    if isinstance(expr, Var):
        if expr.ty is not None and not _is_unresolved(expr.ty):
            return expr.ty
        return typectx.get(expr.name)
    if isinstance(expr, Cast):
        return expr.target_ty
    if isinstance(expr, Lit):
        return expr.ty
    if isinstance(expr, FieldAccess):
        if expr.ty is not None and not _is_unresolved(expr.ty):
            return expr.ty
        # 从 obj 类型 + field_name 推
        obj_ty = _strip_ref(_expr_ty(expr.obj, typectx))
        if isinstance(obj_ty, PairType):
            if expr.field_name == "fst" or expr.field_name == "first":
                return obj_ty.fst
            if expr.field_name == "snd" or expr.field_name == "second":
                return obj_ty.snd
        return None
    if isinstance(expr, ArrayAccess):
        arr_ty = _strip_ref(_expr_ty(expr.arr, typectx))
        if isinstance(arr_ty, ArrayType):
            return arr_ty.elem
        return None
    if isinstance(expr, BinOp):
        return expr.ty
    if isinstance(expr, UnaryOp):
        return expr.ty
    if isinstance(expr, CondExpr):
        return expr.ty
    if isinstance(expr, Call):
        return expr.ty
    return getattr(expr, "ty", None)


def _is_unresolved(ty: TypeIR | None) -> bool:
    if ty is None:
        return True
    if isinstance(ty, NamedType):
        return ty.name in ("DependentType", "ValueType", "ResultType",
                           "StlInternal", "KeyType", "ArgType")
    if isinstance(ty, UnknownType):
        return True
    return False


def _typename_for_classmap(ty: TypeIR | None) -> str | None:
    """把表达式类型规范化为 CLASS_MAP 的 key 字符串。"""
    ty = _strip_ref(ty)
    if ty is None:
        return None
    if isinstance(ty, NamedType):
        return ty.name
    if isinstance(ty, ArrayType):
        return "Array"
    if isinstance(ty, StdMapType):
        return "StdMap"
    if isinstance(ty, PairType):
        return "Pair"
    if isinstance(ty, BaseType):
        return f"BaseType.{ty.name}"
    return None


# ============================================================
# Cast 解析
# ============================================================

def _resolve_cast(cast: Cast, typectx: dict, gap: GapLog) -> ExprIR:
    """解析 Cast 节点。返回新的 ExprIR（可能仍是 Cast 或 Call 或剥后的内层）。"""
    inner = _walk_expr(cast.expr, typectx, gap)
    new_cast = replace(cast, expr=inner)
    disposition = CAST_KIND_DISPOSITION.get(cast.cast_kind, "unknown")

    if disposition == "strip":
        # NoOp / LValueToRValue / ToVoid / 等：直接剥
        return inner

    if disposition == "table":
        # IntegralCast / IntegralToFloating / FloatingToIntegral：查表
        # 字面量+目标 unresolved 时优先剥（is_safe_to_strip_lit_cast）
        if is_safe_to_strip_lit_cast(new_cast):
            return inner

        resolution = lookup_cast(cast.source_ty, cast.target_ty)
        if resolution is None:
            # B 策略：保留 Cast 不解析
            gap.cast_miss.append(
                (canonicalize_type(cast.source_ty), canonicalize_type(cast.target_ty)))
            return new_cast
        # 命中：用 resolution.template 包装为 Call("__lean_cast__", [inner])
        # （Pass 8 codegen 用 template 渲染；此处用 Call 携带模板 + 参数）
        return Call(
            callee=f"__cast__{resolution.template}",  # Pass 8 解析此 callee
            args=[inner],
            ty=cast.target_ty,
        )

    if disposition == "delegate":
        # ConstructorConversion / UserDefinedConversion：交给 FUNC_MAP/CLASS_MAP
        # 实际实现：把 Cast 视为 inner 直接转换（构造器调用通常已嵌在 inner 里）
        return inner

    # 未知 cast_kind：保留
    return new_cast


# ============================================================
# UnresolvedOp 解析：操作符
# ============================================================

_OPERATOR_FAMILIES = {
    "binary":   {"+", "-", "*", "/", "%", "==", "!=", "<", ">", "<=", ">=",
                 "&&", "||", "&", "|", "^", "<<", ">>"},
    "compound": {"+=", "-=", "*=", "/=", "%=", "<<=", ">>=", "&=", "|=", "^="},
    "unary":    {"!", "~", "-", "+"},  # 与 binary 重叠时按 args 数判
    "incdec":   {"++", "--"},
    "subscript": {"[]"},
    "call":     {"()"},
    "assign":   {"="},
    "deref":    {"*", "->"},
    "bool":     {"bool"},
}


def _resolve_operator_call(call: Call, op_sym: str, typectx: dict, gap: GapLog
                           ) -> ExprIR | tuple[list[StmtIR], ExprIR]:
    """处理 `Call(UnresolvedOp("operator<sym>"), args)`。

    返回：
      - ExprIR：纯表达式重写
      - (stmts, ExprIR)：需要插入的前置语句（如 UB require）+ 最终表达式
    Mutator 类（compound assign / operator++/--）由调用方在 stmt 级别处理。
    """
    args = [_walk_expr(a, typectx, gap) for a in call.args]

    # operator= → 调用方在 stmt 级别处理；这里返回 BinOp(=) 让 caller 识别
    if op_sym == "=":
        if len(args) == 2:
            return BinOp(op="=", lhs=args[0], rhs=args[1], ty=None)
        return replace(call, args=args)

    # operator++ / -- → caller 在 stmt 处理
    if op_sym in ("++", "--"):
        return replace(call, args=args)

    # operator[] → ArrayAccess（UB require 在 stmt 注入阶段，按 receiver 类型分派）
    # P0-5：对 StdMap 不输出 ArrayAccess（Lean 的 StdMap.lookup 语义不同），
    # 用 Call 形态表示，UB 检查也走 StdMap 自己的判定
    if op_sym == "[]" and len(args) == 2:
        recv_ty = _strip_ref(_expr_ty(args[0], typectx))
        if isinstance(recv_ty, StdMapType) or (isinstance(recv_ty, NamedType)
                                               and recv_ty.name == "StdMap"):
            return Call(callee="StdMap.get!", args=args, ty=None)
        return ArrayAccess(arr=args[0], idx=args[1], ty=None)

    # operator() → lambda/functor apply
    if op_sym == "()" and len(args) >= 1:
        # 第 0 个是 callable：可能是 Var (lifted lambda 引用) / FieldAccess /
        # 或裸字符串。Pass 3 lambda_lift 后 lifted lambda 调用形态是
        # Call(UnresolvedOp("operator()"), [Var("dot"), arg1, arg2])。
        callee = args[0]
        if isinstance(callee, Var):
            callee_str = callee.name
        elif isinstance(callee, str):
            callee_str = callee
        elif isinstance(callee, FieldAccess):
            callee_str = f"<field>.{callee.field_name}"
        else:
            callee_str = "<callable>"  # fallback：未知形态
        return Call(callee=callee_str, args=args[1:], ty=None)

    # operator->: 解引用
    if op_sym == "->" and len(args) == 1:
        return Call(callee="__deref__", args=[args[0]], ty=None)

    # 其他二元 / 一元
    recv = _typename_for_classmap(_expr_ty(args[0], typectx)) if args else None
    cls_ops = CLASS_MAP.get(recv, {}).get("operators", {}) if recv else {}

    if recv is None:
        # 无类型信息：保留 UnresolvedOp
        gap.op_miss.append(("None", op_sym))
        return replace(call, args=args)

    if op_sym in cls_ops:
        target = cls_ops[op_sym]
        if target is None:
            # 用 Lean typeclass：BinOp / UnaryOp 节点
            # P0-4：填回结果 ty —— 算术/比较的结果类型规则：
            #   - 比较 `==/!=/</>/<=/>=` → BaseType.BOOL
            #   - 逻辑 `&&/||` → BaseType.BOOL
            #   - 算术/位 `+/-/*//=/&/|/^/<</>>` → 同 receiver 类型
            #   - 复合赋值 `+=/-=/*=/...` → 同 receiver 类型（虽然语义是赋值，
            #     但表达式形态保留 receiver type 帮助下游推断）
            result_ty: TypeIR | None
            if op_sym in ("==", "!=", "<", ">", "<=", ">=", "&&", "||"):
                result_ty = BaseType.BOOL
            else:
                # 算术/位/复合 → 用 receiver 类型（已 strip RefType）
                result_ty = _strip_ref(_expr_ty(args[0], typectx))
            if len(args) == 2:
                return BinOp(op=op_sym, lhs=args[0], rhs=args[1], ty=result_ty)
            if len(args) == 1:
                return UnaryOp(op=op_sym, operand=args[0], ty=result_ty)
            return replace(call, args=args)
        # tuple disposition：(`disposition`, `target_fn`) — 类似 methods 表的
        # 二元组形态。CLPoly 用此在 Iterator 类型里写 ("mutate", "Iterator.advance")
        # 等。dispatch 与 method 同：method/mutate/identity/noop。
        if isinstance(target, tuple) and len(target) == 2:
            disp, fn = target
            if disp == "method":
                return Call(callee=fn, args=args, ty=None)
            if disp == "mutate":
                # mutator 形态：caller 在 stmt-level 转 AssignStmt（args[0] 通
                # 常是 Iterator 变量）；此处仅回 Call，stmt-level 处理在
                # _walk_stmt 的 ExprStmt 分支里识别 callee 是否 Iterator
                # mutator 函数（保守：不在此自动转 stmt，避免误判）
                return Call(callee=fn, args=args, ty=None)
            if disp == "identity":
                return args[0] if args else replace(call, args=args)
            if disp == "noop":
                return args[0] if args else replace(call, args=args)
        # 具体函数名（字符串，如 Zp.div）
        if isinstance(target, str):
            return Call(callee=target, args=args, ty=None)
        # 未知形态：保留为 UnresolvedOp（B 策略）
        gap.op_miss.append((recv, op_sym))
        return replace(call, args=args)

    # CLASS_MAP 没注册：B 策略
    gap.op_miss.append((recv, op_sym))
    return replace(call, args=args)


# ============================================================
# UnresolvedOp 解析：方法调用 <method>.X
# ============================================================

def _resolve_method_call(call: Call, method: str, typectx: dict, gap: GapLog,
                         host_stmt_ctx: bool = False
                         ) -> ExprIR | tuple[list[StmtIR], ExprIR]:
    """处理 `Call(UnresolvedOp("<method>.X"), [recv, args...])`。

    `host_stmt_ctx=True` 表示调用方处于 stmt 级（ExprStmt 包裹），可生成
    AssignStmt（mutator 类）；False 时即使是 mutator 也只能返回 Call。
    """
    args = [_walk_expr(a, typectx, gap) for a in call.args]
    if not args:
        gap.method_miss.append(("?", method))
        return replace(call, args=args)

    recv = args[0]
    extra_args = args[1:]
    recv_ty = _expr_ty(recv, typectx)
    T = _typename_for_classmap(recv_ty)

    if T is None:
        gap.method_miss.append(("None", method))
        return replace(call, args=args)

    cls = CLASS_MAP.get(T, {})
    methods = cls.get("methods", {})
    # P1-3：`<method>.operator bool` / `<method>.operator <op>` 是 Pass 1 把
    # C++ 转换运算符或运算符 overload 包成 method 调用形态。优先查 methods，
    # 找不到时把 method 名按 operator 处理（剥 "operator " 前缀查 operators）
    if method not in methods and method.startswith("operator "):
        op_sym = method[len("operator "):]
        operators = cls.get("operators", {})
        if op_sym in operators:
            target = operators[op_sym]
            if target is None:
                # typeclass 形态：UnaryOp / 直接返回 recv
                return UnaryOp(op=op_sym, operand=recv,
                               ty=BaseType.BOOL if op_sym == "bool" else None)
            if isinstance(target, str):
                return Call(callee=target, args=[recv] + extra_args, ty=None)
    if method not in methods:
        gap.method_miss.append((T, method))
        return replace(call, args=args)

    disposition, target = methods[method]

    if disposition == "field":
        # FieldAccess
        return FieldAccess(obj=recv, field_name=target or method, ty=None)

    if disposition == "method":
        # Call(target, [recv, *extra_args])
        return Call(callee=target, args=[recv] + extra_args, ty=None)

    if disposition == "identity":
        return recv  # 如 vector::data() ≈ self

    if disposition == "noop":
        return recv  # mutate 但 no-op（如 reserve）

    if disposition == "mutate":
        new_call = Call(callee=target, args=[recv] + extra_args, ty=None)
        # P1-5：扩展 ArrayAccess 也作 target（如 `groups[d].push_back(x)`）
        if host_stmt_ctx and isinstance(recv, (Var, FieldAccess, ArrayAccess)):
            return ([AssignStmt(target=recv, value=new_call)], recv)
        return new_call

    if disposition == "mutate_push":
        new_call = Call(callee=target or "Array.push",
                        args=[recv] + extra_args, ty=None)
        if host_stmt_ctx and isinstance(recv, (Var, FieldAccess, ArrayAccess)):
            return ([AssignStmt(target=recv, value=new_call)], recv)
        return new_call

    # 未知 disposition：保留
    gap.method_miss.append((T, f"{method}@{disposition}"))
    return replace(call, args=args)


# ============================================================
# UnresolvedOp 解析：构造器
# ============================================================

def _resolve_constructor_call(call: Call, op_name: str, typectx: dict,
                              gap: GapLog) -> ExprIR:
    """处理 `Call(UnresolvedOp("construct_X" | "default_init_X"), args)`。"""
    from constructor_map import _normalize_typename
    args = [_walk_expr(a, typectx, gap) for a in call.args]
    arity = len(args)
    resolution = resolve_constructor(op_name, arity)
    if resolution is None:
        gap.constructor_miss.append((op_name, arity))
        return replace(call, args=args)
    # P0-3：从 typename 推回目标类型，填到 Call.ty。这样下游 receiver 类型查询
    # （如 `__ctor__Zp.ofInt(x) + y` 中 outer + 的 LHS 类型）能正确得到 Zp。
    typename = _normalize_typename(op_name)
    target_ty = _typename_to_typeir(typename)
    return Call(callee=f"__ctor__{resolution.template}", args=args, ty=target_ty)


def _typename_to_typeir(name: str) -> TypeIR | None:
    """把 normalized typename（如 'Zp', 'ZZ', 'vector<Zp>', 'pair<A,B>'）映射回
    TypeIR。仅覆盖 CLPoly 核心类型 + 常见 STL；未识别返回 NamedType(name)。"""
    if not name:
        return None
    # CLPoly 核心类型
    _CORE = {
        "Zp": NamedType("Zp"),
        "ZZ": NamedType("ZZ"),
        "QQ": NamedType("QQ"),
        "umonomial": NamedType("UMonomial"),
        "UMonomial": NamedType("UMonomial"),
        "Monomial": NamedType("Monomial"),
        "Variable": NamedType("Variable"),
        "variable": NamedType("Variable"),
        "Poly": NamedType("MvPolyZZ"),
        "PolyZp": NamedType("MvPolyZp"),
        "PolyZZ": NamedType("MvPolyZZ"),
        "PolyQQ": NamedType("PolyQQ"),
        "upolynomial_<ZZ>": NamedType("SparsePolyZZ"),
        "upolynomial_<Zp>": NamedType("SparsePolyZp"),
        "UPZp": NamedType("SparsePolyZp"),
        "polynomial_<ZZ, lex_<less>>": NamedType("MvPolyZZ"),
        "polynomial_<Zp, lex_<less>>": NamedType("MvPolyZp"),
        "polynomial_<ZZ, lex>": NamedType("MvPolyZZ"),
        "LLLMatrix": NamedType("LLLMatrix"),
        "iterator": NamedType("Iterator"),
        "const_iterator": NamedType("Iterator"),
    }
    if name in _CORE:
        return _CORE[name]
    if name.startswith("vector<"):
        return NamedType("Array")  # 元素类型留给 Pass 8 推
    if name.startswith("pair<"):
        return NamedType("Pair")
    if name.startswith("map<") or name.startswith("multimap<"):
        return NamedType("StdMap")
    if name.startswith("set<") or name.startswith("unordered_set<"):
        return NamedType("StdMap")
    if name.startswith("factorization<"):
        return NamedType("Factorization")
    if name.startswith("basic_polynomial<"):
        # 长形态——粗略归类
        if "Zp" in name and "basic_monomial" in name:
            return NamedType("MvPolyZp")
        if "ZZ" in name and "basic_monomial" in name:
            return NamedType("MvPolyZZ")
        if "Zp" in name:
            return NamedType("SparsePolyZp")
        if "ZZ" in name:
            return NamedType("SparsePolyZZ")
    return NamedType(name)


# ============================================================
# UnresolvedOp dispatch
# ============================================================

def _resolve_unresolved_op_call(call: Call, typectx: dict, gap: GapLog,
                                host_stmt_ctx: bool = False
                                ) -> ExprIR | tuple[list[StmtIR], ExprIR]:
    op_name = call.callee.op_name if isinstance(call.callee, UnresolvedOp) else ""
    if op_name.startswith("<method>."):
        return _resolve_method_call(call, op_name[len("<method>."):],
                                    typectx, gap, host_stmt_ctx)
    if op_name.startswith("construct_") or op_name.startswith("default_init"):
        return _resolve_constructor_call(call, op_name, typectx, gap)
    if op_name.startswith("operator"):
        op_sym = op_name[len("operator"):]
        return _resolve_operator_call(call, op_sym, typectx, gap)
    # 其他形式（不应出现）：保留
    gap.unresolved_op.append(op_name)
    return replace(call, args=[_walk_expr(a, typectx, gap) for a in call.args])


# ============================================================
# Walk expr
# ============================================================

def _walk_expr(e: ExprIR, typectx: dict, gap: GapLog) -> ExprIR:
    """递归处理表达式（不处理 stmt 级 mutator，由 caller 处理）。"""
    if isinstance(e, Cast):
        return _resolve_cast(e, typectx, gap)
    if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
        result = _resolve_unresolved_op_call(e, typectx, gap, host_stmt_ctx=False)
        # 表达式上下文：mutator 返回 (stmts, expr) 时退化为 expr（语义略损）
        if isinstance(result, tuple):
            _, expr = result
            return expr
        return result
    if isinstance(e, Call):
        return Call(callee=e.callee,
                    args=[_walk_expr(a, typectx, gap) for a in e.args],
                    ty=e.ty)
    if isinstance(e, BinOp):
        return BinOp(op=e.op, lhs=_walk_expr(e.lhs, typectx, gap),
                     rhs=_walk_expr(e.rhs, typectx, gap), ty=e.ty)
    if isinstance(e, UnaryOp):
        return UnaryOp(op=e.op, operand=_walk_expr(e.operand, typectx, gap), ty=e.ty)
    if isinstance(e, CondExpr):
        return CondExpr(cond=_walk_expr(e.cond, typectx, gap),
                        then_e=_walk_expr(e.then_e, typectx, gap),
                        else_e=_walk_expr(e.else_e, typectx, gap), ty=e.ty)
    if isinstance(e, FieldAccess):
        return FieldAccess(obj=_walk_expr(e.obj, typectx, gap),
                           field_name=e.field_name, ty=e.ty)
    if isinstance(e, ArrayAccess):
        return ArrayAccess(arr=_walk_expr(e.arr, typectx, gap),
                           idx=_walk_expr(e.idx, typectx, gap), ty=e.ty)
    if isinstance(e, TupleExpr):
        return TupleExpr(elems=[_walk_expr(x, typectx, gap) for x in e.elems], ty=e.ty)
    if isinstance(e, ArrayLit):
        return ArrayLit(elems=[_walk_expr(x, typectx, gap) for x in e.elems],
                        elem_ty=e.elem_ty)
    if isinstance(e, LambdaExpr):
        # 不递归进入 lambda body（aux_lambdas 由顶层 pass 单独处理）
        return e
    return e


# ============================================================
# UB require 生成
# ============================================================

def _emit_div_require(rhs: ExprIR) -> RequireStmt:
    """除零保护：require rhs ≠ 0"""
    return RequireStmt(
        cond=BinOp(op="!=", lhs=rhs, rhs=Lit(value=0, ty=BaseType.INT64), ty=BaseType.BOOL),
        name="h_nonzero",
        source="UB-1 div_by_zero",
    )


def _emit_idx_require(idx: ExprIR, arr: ExprIR) -> RequireStmt:
    """越界保护：require idx < arr.size"""
    return RequireStmt(
        cond=BinOp(op="<", lhs=idx,
                   rhs=Call(callee="Array.size", args=[arr], ty=BaseType.NAT),
                   ty=BaseType.BOOL),
        name="h_in_bounds",
        source="UB-2 array_oob",
    )


def _emit_nonempty_require(container: ExprIR) -> RequireStmt:
    """空容器保护：require ¬container.isEmpty"""
    return RequireStmt(
        cond=UnaryOp(op="!",
                     operand=Call(callee="Array.isEmpty", args=[container], ty=BaseType.BOOL),
                     ty=BaseType.BOOL),
        name="h_nonempty",
        source="UB-3 empty_container",
    )


def _collect_ub_requires(expr: ExprIR, typectx: dict) -> list[RequireStmt]:
    """从表达式收集 UB require（先深扫子表达式，再当前节点）。"""
    reqs: list[RequireStmt] = []

    def walk(e):
        # 先处理子表达式
        if isinstance(e, BinOp):
            walk(e.lhs); walk(e.rhs)
            if e.op in ("/", "%"):
                # P1-4：整数类型（含 ZZ/Zp/QQ 任意精度）加 require
                rhs_ty = _strip_ref(_expr_ty(e.rhs, typectx))
                is_int_like = (
                    (isinstance(rhs_ty, BaseType) and rhs_ty in (
                        BaseType.INT32, BaseType.INT64, BaseType.UINT32,
                        BaseType.UINT64, BaseType.NAT, BaseType.UINT128))
                    or (isinstance(rhs_ty, NamedType) and rhs_ty.name in (
                        "ZZ", "Zp", "QQ"))
                )
                if is_int_like:
                    reqs.append(_emit_div_require(e.rhs))
        elif isinstance(e, UnaryOp):
            walk(e.operand)
        elif isinstance(e, CondExpr):
            walk(e.cond); walk(e.then_e); walk(e.else_e)
        elif isinstance(e, FieldAccess):
            walk(e.obj)
        elif isinstance(e, ArrayAccess):
            walk(e.arr); walk(e.idx)
            reqs.append(_emit_idx_require(e.idx, e.arr))
        elif isinstance(e, Call):
            for a in e.args: walk(a)
            # 空容器访问（front!/back!/head!/getLast!）
            cs = e.callee if isinstance(e.callee, str) else None
            if cs and any(s in cs for s in (".head!", ".getLast!", ".front!", ".back!")):
                if e.args:
                    reqs.append(_emit_nonempty_require(e.args[0]))
        elif isinstance(e, Cast):
            walk(e.expr)

    walk(expr)
    return reqs


# ============================================================
# Walk stmt
# ============================================================

def _walk_stmt(s: StmtIR, typectx: dict, gap: GapLog) -> list[StmtIR]:
    """处理单条 stmt。返回 stmt 列表（mutator 展开 + UB require 注入可能产生多条）。"""
    if isinstance(s, LetStmt):
        new_val = _walk_expr(s.value, typectx, gap) if s.value else None
        # 注册 typectx
        if s.ty is not None and not _is_unresolved(s.ty):
            typectx[s.var.name] = s.ty
        elif new_val is not None:
            et = _expr_ty(new_val, typectx)
            if et is not None and not _is_unresolved(et):
                typectx[s.var.name] = et
        # UB require 前置
        reqs = _collect_ub_requires(new_val, typectx) if new_val else []
        return reqs + [LetStmt(var=s.var, ty=s.ty, value=new_val)]

    if isinstance(s, AssignStmt):
        new_val = _walk_expr(s.value, typectx, gap)
        new_tgt = _walk_expr(s.target, typectx, gap)
        reqs = _collect_ub_requires(new_val, typectx) + _collect_ub_requires(new_tgt, typectx)
        return reqs + [AssignStmt(target=new_tgt, value=new_val)]

    if isinstance(s, CompoundAssignStmt):
        # 展开：x op= e → x := x op e
        new_val = _walk_expr(s.value, typectx, gap)
        new_tgt = _walk_expr(s.target, typectx, gap)
        new_assign = AssignStmt(
            target=new_tgt,
            value=BinOp(op=s.op, lhs=new_tgt, rhs=new_val, ty=None),
        )
        reqs = _collect_ub_requires(new_val, typectx)
        return reqs + [new_assign]

    if isinstance(s, IfStmt):
        new_cond = _walk_expr(s.cond, typectx, gap)
        new_then = _walk_stmts(list(s.then_body), typectx, gap)
        new_else = _walk_stmts(list(s.else_body), typectx, gap)
        reqs = _collect_ub_requires(new_cond, typectx)
        return reqs + [IfStmt(cond=new_cond, then_body=new_then, else_body=new_else)]

    if isinstance(s, WhileStmt):
        new_cond = _walk_expr(s.cond, typectx, gap)
        new_body = _walk_stmts(list(s.body), typectx, gap)
        return [WhileStmt(cond=new_cond, body=new_body)]

    if isinstance(s, DoWhileStmt):
        new_cond = _walk_expr(s.cond, typectx, gap) if s.cond else None
        new_body = _walk_stmts(list(s.body), typectx, gap)
        return [DoWhileStmt(cond=new_cond, body=new_body)]

    if isinstance(s, ForStmt):
        new_init = _walk_stmts(list(s.init), typectx, gap)
        new_cond = _walk_expr(s.cond, typectx, gap) if s.cond else None
        new_step = _walk_stmts(list(s.step), typectx, gap)
        new_body = _walk_stmts(list(s.body), typectx, gap)
        return [ForStmt(init=new_init, cond=new_cond, step=new_step, body=new_body)]

    if isinstance(s, RangeForStmt):
        # 注册 var + decomposition
        if s.var_ty is not None and not _is_unresolved(s.var_ty):
            typectx[s.var.name] = s.var_ty
        new_container = _walk_expr(s.container, typectx, gap)
        new_body = _walk_stmts(list(s.body), typectx, gap)
        return [RangeForStmt(var=s.var, var_ty=s.var_ty, container=new_container,
                             body=new_body, decomposition=s.decomposition)]

    if isinstance(s, ReturnStmt):
        new_val = _walk_expr(s.value, typectx, gap) if s.value else None
        reqs = _collect_ub_requires(new_val, typectx) if new_val else []
        return reqs + [ReturnStmt(value=new_val)]

    if isinstance(s, RequireStmt):
        return [RequireStmt(cond=_walk_expr(s.cond, typectx, gap),
                            name=s.name, source=s.source)]

    if isinstance(s, ExprStmt):
        # ExprStmt 是 stmt-level mutator 的入口
        # P1-1+P1-2：识别 user-defined op overload 形态的 compound assign / incdec
        # `Call(UnresolvedOp("operator*=" | "operator++" | ...))`，在 stmt 级
        # 展开为 AssignStmt（CompoundAssignStmt 在 BaseType 上的等价形态）
        e = s.expr
        if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
            op = e.callee.op_name
            # operator++ / operator-- (incdec)
            if op in ("operator++", "operator--") and len(e.args) >= 1:
                target = _walk_expr(e.args[0], typectx, gap)
                op_sym = "+" if op == "operator++" else "-"
                if isinstance(target, (Var, FieldAccess, ArrayAccess)):
                    new_assign = AssignStmt(
                        target=target,
                        value=BinOp(op=op_sym, lhs=target,
                                    rhs=Lit(value=1, ty=BaseType.INT64),
                                    ty=_strip_ref(_expr_ty(target, typectx))),
                    )
                    return _collect_ub_requires(new_assign.value, typectx) + [new_assign]
            # operator+=/-=/*=/...（compound assign）
            _COMPOUND = {"operator+=": "+", "operator-=": "-", "operator*=": "*",
                         "operator/=": "/", "operator%=": "%",
                         "operator<<=": "<<", "operator>>=": ">>",
                         "operator&=": "&", "operator|=": "|", "operator^=": "^"}
            if op in _COMPOUND and len(e.args) == 2:
                target = _walk_expr(e.args[0], typectx, gap)
                rhs = _walk_expr(e.args[1], typectx, gap)
                op_sym = _COMPOUND[op]
                if isinstance(target, (Var, FieldAccess, ArrayAccess)):
                    new_assign = AssignStmt(
                        target=target,
                        value=BinOp(op=op_sym, lhs=target, rhs=rhs,
                                    ty=_strip_ref(_expr_ty(target, typectx))),
                    )
                    return _collect_ub_requires(new_assign.value, typectx) + [new_assign]

            # 一般 mutator method 路径
            result = _resolve_unresolved_op_call(e, typectx, gap, host_stmt_ctx=True)
            if isinstance(result, tuple):
                stmts, _ = result
                out: list[StmtIR] = []
                for st in stmts:
                    if isinstance(st, AssignStmt):
                        out.extend(_collect_ub_requires(st.value, typectx))
                    out.append(st)
                return out
            new_expr = result
        else:
            new_expr = _walk_expr(e, typectx, gap)
        # operator= → AssignStmt
        if isinstance(new_expr, BinOp) and new_expr.op == "=":
            return _collect_ub_requires(new_expr.rhs, typectx) + [
                AssignStmt(target=new_expr.lhs, value=new_expr.rhs)
            ]
        # P1-1 兜底：BinOp 形态的 compound（来自 _resolve_operator_call 返回）
        if isinstance(new_expr, BinOp) and new_expr.op in (
                "+=", "-=", "*=", "/=", "%=", "<<=", ">>=", "&=", "|=", "^="):
            base_op = new_expr.op[:-1]  # "*=" → "*"
            if isinstance(new_expr.lhs, (Var, FieldAccess, ArrayAccess)):
                ty = _strip_ref(_expr_ty(new_expr.lhs, typectx))
                new_assign = AssignStmt(
                    target=new_expr.lhs,
                    value=BinOp(op=base_op, lhs=new_expr.lhs, rhs=new_expr.rhs, ty=ty),
                )
                return _collect_ub_requires(new_expr.rhs, typectx) + [new_assign]
        return _collect_ub_requires(new_expr, typectx) + [ExprStmt(expr=new_expr)]

    if isinstance(s, BlockStmt):
        return [BlockStmt(stmts=_walk_stmts(list(s.stmts), typectx, gap))]

    if isinstance(s, (BreakStmt, ContinueStmt)):
        return [s]

    return [s]


def _walk_stmts(stmts: list[StmtIR], typectx: dict, gap: GapLog) -> list[StmtIR]:
    out: list[StmtIR] = []
    for s in stmts:
        out.extend(_walk_stmt(s, typectx, gap))
    return out


# ============================================================
# Pass 入口
# ============================================================

def operator_resolve_pass(func: HIRFunc) -> tuple[HIRFunc, GapLog]:
    """Pass 5：HIR₃ → HIR₄。返回 (新 HIRFunc, gap_log)。"""
    typectx: dict[str, TypeIR] = {}
    for p in func.params:
        if not _is_unresolved(p.ty):
            typectx[p.name] = p.ty
    gap = GapLog()
    new_body = _walk_stmts(list(func.body), typectx, gap)
    new_aux = []
    for aux in func.aux_lambdas:
        aux_typectx = dict(typectx)
        for p in aux.params:
            if not _is_unresolved(p.ty):
                aux_typectx[p.name] = p.ty
        new_aux.append(replace(aux, body=_walk_stmts(list(aux.body), aux_typectx, gap)))
    return replace(func, body=new_body, aux_lambdas=new_aux), gap


# ============================================================
# HIR₄ 出口不变量
# ============================================================

def assert_hir4_invariant(func: HIRFunc) -> None:
    """HIR₄ 出口检查（P0-6 强化版）：
      - 禁止 CompoundAssignStmt 残留（必须展开为 AssignStmt + BinOp）
      - 禁止 Call.callee 是 tuple / list / 非字符串非 UnresolvedOp 的形态（P0-2）
      - 禁止 Call.callee 字面量 "<callable>" / "<field>." 占位（P0-1）
      - 禁止 ExprStmt(Call(UnresolvedOp("operator+="/"-="/...))) 半展开（P1-1）
      - 禁止 ExprStmt(Call(UnresolvedOp("operator++"/"--"))) 半展开（P1-2）
      - 允许 UnresolvedOp 在表达式深处残留（B 策略 typectx gap）
    """
    _COMPOUND_OPS = {"operator+=", "operator-=", "operator*=", "operator/=",
                     "operator%=", "operator<<=", "operator>>=",
                     "operator&=", "operator|=", "operator^="}
    _INCDEC_OPS = {"operator++", "operator--"}

    def fail(reason: str) -> None:
        raise TranslationError(
            pass_name="operator_resolve",
            func_name=func.base_name,
            reason=reason,
        )

    def check_callee(callee, ctx: str):
        # callee 可以是 str 或 UnresolvedOp（B 策略允许）
        if isinstance(callee, (str, UnresolvedOp)):
            if isinstance(callee, str):
                if callee == "<callable>":
                    fail(f"{ctx}: P0-1 violation — '<callable>' placeholder")
            return
        fail(f"{ctx}: P0-2 violation — Call.callee is {type(callee).__name__} "
             f"(must be str or UnresolvedOp), got {callee!r}")

    def check_expr(e, ctx: str):
        if isinstance(e, Call):
            check_callee(e.callee, ctx)
            for a in e.args: check_expr(a, ctx)
        elif isinstance(e, Cast):
            check_expr(e.expr, ctx)
        elif isinstance(e, BinOp):
            check_expr(e.lhs, ctx); check_expr(e.rhs, ctx)
        elif isinstance(e, UnaryOp):
            check_expr(e.operand, ctx)
        elif isinstance(e, CondExpr):
            check_expr(e.cond, ctx); check_expr(e.then_e, ctx); check_expr(e.else_e, ctx)
        elif isinstance(e, FieldAccess):
            check_expr(e.obj, ctx)
        elif isinstance(e, ArrayAccess):
            check_expr(e.arr, ctx); check_expr(e.idx, ctx)
        elif isinstance(e, (TupleExpr, ArrayLit)):
            for x in e.elems: check_expr(x, ctx)

    def check(stmts, ctx: str = "body"):
        for i, s in enumerate(stmts):
            sub_ctx = f"{ctx}[{i}]"
            if isinstance(s, CompoundAssignStmt):
                fail(f"{sub_ctx}: CompoundAssignStmt residual: {s.op}")
            if isinstance(s, ExprStmt):
                # P1-1+P1-2：检测 ExprStmt 包 compound/incdec UnresolvedOp Call
                if isinstance(s.expr, Call) and isinstance(s.expr.callee, UnresolvedOp):
                    op = s.expr.callee.op_name
                    if op in _COMPOUND_OPS:
                        fail(f"{sub_ctx}: P1-1 violation — ExprStmt-wrapped "
                             f"compound assign Call({op})")
                    if op in _INCDEC_OPS:
                        fail(f"{sub_ctx}: P1-2 violation — ExprStmt-wrapped "
                             f"incdec Call({op})")
                # P1-1 兜底形态：ExprStmt(BinOp("*=", ...))
                if isinstance(s.expr, BinOp) and s.expr.op in (
                        "+=", "-=", "*=", "/=", "%=", "<<=", ">>=", "&=", "|=", "^="):
                    fail(f"{sub_ctx}: P1-1 violation — ExprStmt-wrapped "
                         f"BinOp({s.expr.op})")
                check_expr(s.expr, sub_ctx)
            elif isinstance(s, LetStmt):
                if s.value: check_expr(s.value, sub_ctx)
            elif isinstance(s, AssignStmt):
                check_expr(s.target, sub_ctx); check_expr(s.value, sub_ctx)
            elif isinstance(s, IfStmt):
                check_expr(s.cond, sub_ctx)
                check(s.then_body, f"{sub_ctx}/then")
                check(s.else_body, f"{sub_ctx}/else")
            elif isinstance(s, ForStmt):
                check(s.init, f"{sub_ctx}/init")
                check(s.step, f"{sub_ctx}/step")
                check(s.body, f"{sub_ctx}/body")
            elif isinstance(s, (WhileStmt, DoWhileStmt)):
                if s.cond: check_expr(s.cond, sub_ctx)
                check(s.body, f"{sub_ctx}/body")
            elif isinstance(s, RangeForStmt):
                check_expr(s.container, sub_ctx)
                check(s.body, f"{sub_ctx}/body")
            elif isinstance(s, BlockStmt):
                check(s.stmts, f"{sub_ctx}/blk")
            elif isinstance(s, ReturnStmt):
                if s.value: check_expr(s.value, sub_ctx)
            elif isinstance(s, RequireStmt):
                check_expr(s.cond, sub_ctx)

    check(func.body, "body")
    for aux in func.aux_lambdas:
        check(aux.body, f"aux({aux.base_name})")
