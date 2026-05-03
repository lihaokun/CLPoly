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
# Iterator-offset pattern 识别（A 类设计 gap 修复）
# ============================================================

def _try_extract_iter_offset(e: ExprIR
                             ) -> tuple[ExprIR, ExprIR | None, str] | None:
    """识别 iterator + offset 形态，返回 (容器, offset_expr_or_None, op)：

    - `Call("Array.toList", [v])`                          → (v, None, "begin")
    - `Call(UnresolvedOp("operator+"), [Call(toList,[v]), o])`→ (v, o, "+")
    - `Call(UnresolvedOp("operator-"), [Call(toList,[v]), o])`→ (v, o, "-")
    - `Call("__ctor__Iterator.fromList {a0}", [inner])`     → 透传 inner

    其他 → None
    """
    # 透传 Iterator wrapper（已 resolve 后形态）
    if isinstance(e, Call):
        cs = e.callee if isinstance(e.callee, str) else None
        if cs and cs.startswith("__ctor__Iterator.fromList"):
            if e.args:
                return _try_extract_iter_offset(e.args[0])
    # 透传 raw Iterator ctor: Call(UnresolvedOp("construct_iterator"/"construct_const_iterator"))
    if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
        if e.callee.op_name in ("construct_iterator", "construct_const_iterator"):
            if e.args:
                return _try_extract_iter_offset(e.args[0])
    # 剥 Cast
    if isinstance(e, Cast):
        return _try_extract_iter_offset(e.expr)
    # toList(v) 单独形式 = begin（已 resolve 后形态）
    if isinstance(e, Call):
        cs = e.callee if isinstance(e.callee, str) else None
        if cs and (cs.endswith(".toList") or cs == "Array.toList"):
            if e.args:
                return (e.args[0], None, "begin")
    # raw 形态：Call(UnresolvedOp("<method>.begin/.end/.cbegin/.cend"))
    if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
        op_name = e.callee.op_name
        if op_name in ("<method>.begin", "<method>.end",
                       "<method>.cbegin", "<method>.cend"):
            if e.args:
                return (e.args[0], None, "begin")
    # operator+ / operator- 嵌套
    if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
        op = e.callee.op_name
        if op in ("operator+", "operator-") and len(e.args) == 2:
            inner_lhs = _try_extract_iter_offset(e.args[0])
            if inner_lhs is not None and inner_lhs[1] is None:
                # lhs 是 begin，offset 是 args[1]
                return (inner_lhs[0], e.args[1], "+" if op == "operator+" else "-")
    # BinOp 形态（万一已被 Lean-推出）
    if isinstance(e, BinOp) and e.op in ("+", "-"):
        inner = _try_extract_iter_offset(e.lhs)
        if inner is not None and inner[1] is None:
            return (inner[0], e.rhs, e.op)
    return None


def _is_lvalue_target(e: ExprIR) -> bool:
    """检查 e 是否合法 AssignStmt target。
    Var / FieldAccess / ArrayAccess 是显式 lvalue；
    Call("StdMap.get!", [m, k]) 是 map subscript 等价 lvalue（Pass 6 _root_var 阶段 A 处理）。
    """
    if isinstance(e, (Var, FieldAccess, ArrayAccess)):
        return True
    if isinstance(e, Call) and isinstance(e.callee, str) \
            and e.callee == "StdMap.get!" and len(e.args) == 2:
        return True
    return False


def _exprs_eq(a: ExprIR, b: ExprIR) -> bool:
    """简化的表达式相等判断（用于 Pattern B 检验 begin 和 end 的容器是同一个）。"""
    if isinstance(a, Var) and isinstance(b, Var):
        return a.name == b.name
    if isinstance(a, FieldAccess) and isinstance(b, FieldAccess):
        return a.field_name == b.field_name and _exprs_eq(a.obj, b.obj)
    if isinstance(a, Cast) and isinstance(b, Cast):
        return _exprs_eq(a.expr, b.expr)
    return False


def _stl_extract_toList_container(arg: ExprIR) -> ExprIR | None:
    """CF-3 阶段 3: 从 `Array.toList(c)` / `<method>.begin(c)` / `<method>.end(c)` /
    Cast(NoOp, ...) 等包装中抽出底层容器 expr。返回 None 表示不匹配。

    HIR3 阶段（Pass 4 输入未到 Pass 5 时）begin/end 仍可能是 UnresolvedOp；
    HIR4（Pass 5 输出）已被 CLASS_MAP 解析成 `Array.toList(c)`。本函数支持
    两种形态：直接 toList Call、或仍未解析的 UnresolvedOp(<method>.begin/end)。
    """
    cur = arg
    if isinstance(cur, Cast):
        cur = cur.expr
    if isinstance(cur, Call):
        cs = cur.callee
        # Pass 5 已解析后的 Array.toList(c)
        if isinstance(cs, str) and cs == "Array.toList" and len(cur.args) == 1:
            return cur.args[0]
        # 未解析 UnresolvedOp("<method>.begin"|"<method>.end") — 极少见
        if isinstance(cs, UnresolvedOp) \
                and cs.op_name in ("<method>.begin", "<method>.end") \
                and len(cur.args) == 1:
            return cur.args[0]
    return None


# ============================================================
# Cast 解析
# ============================================================

_CAST_PENDING_REQUIRES: list = []  # P0-5: 临时收集 Cast 触发的 UB require


def _resolve_cast(cast: Cast, typectx: dict, gap: GapLog) -> ExprIR:
    """解析 Cast 节点。返回新的 ExprIR（可能仍是 Cast 或 Call 或剥后的内层）。

    P0-5（轮 2 修复）：CAST_TABLE 命中后若有 ub_kind，生成对应 RequireStmt
    挂到 `_CAST_PENDING_REQUIRES` 全局列表，由 `_collect_ub_requires` 统一发射。
    """
    inner = _walk_expr(cast.expr, typectx, gap)
    new_cast = replace(cast, expr=inner)
    disposition = CAST_KIND_DISPOSITION.get(cast.cast_kind, "unknown")

    if disposition == "strip":
        return inner

    if disposition == "table":
        if is_safe_to_strip_lit_cast(new_cast):
            return inner

        resolution = lookup_cast(cast.source_ty, cast.target_ty)
        if resolution is None:
            gap.cast_miss.append(
                (canonicalize_type(cast.source_ty), canonicalize_type(cast.target_ty)))
            return new_cast
        # P0-5：发射 ub_kind 对应的 RequireStmt
        if resolution.ub_kind:
            _emit_cast_ub_require(inner, resolution.ub_kind)
        return Call(
            callee=f"__cast__{resolution.template}",
            args=[inner],
            ty=cast.target_ty,
        )

    if disposition == "delegate":
        # ConstructorConversion / UserDefinedConversion：剥外层（构造器调用已在 inner）
        # P0-1 副作用：若 source==target 实质是 noop，inner 已是目标类型
        return inner

    return new_cast


def _emit_cast_ub_require(inner: ExprIR, ub_kind: str) -> None:
    """根据 cast 的 ub_kind 生成 RequireStmt 挂到 pending list。"""
    cond: ExprIR | None = None
    name = "h_cast"
    if ub_kind == "nonneg":
        cond = BinOp(op=">=", lhs=inner,
                     rhs=Lit(value=0, ty=BaseType.INT64), ty=BaseType.BOOL)
        name = "h_nonneg"
    elif ub_kind == "fits_int32":
        cond = BinOp(
            op="&&",
            lhs=BinOp(op=">=", lhs=inner,
                      rhs=Lit(value=-(2**31), ty=BaseType.INT64), ty=BaseType.BOOL),
            rhs=BinOp(op="<=", lhs=inner,
                      rhs=Lit(value=2**31 - 1, ty=BaseType.INT64), ty=BaseType.BOOL),
            ty=BaseType.BOOL,
        )
        name = "h_fits_int32"
    elif ub_kind == "fits_int64":
        cond = BinOp(
            op="&&",
            lhs=BinOp(op=">=", lhs=inner,
                      rhs=Lit(value=-(2**63), ty=BaseType.INT64), ty=BaseType.BOOL),
            rhs=BinOp(op="<=", lhs=inner,
                      rhs=Lit(value=2**63 - 1, ty=BaseType.INT64), ty=BaseType.BOOL),
            ty=BaseType.BOOL,
        )
        name = "h_fits_int64"
    if cond is not None:
        _CAST_PENDING_REQUIRES.append(
            RequireStmt(cond=cond, name=name, source=f"UB cast: {ub_kind}")
        )


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

    # P2-A（轮 2 修复）：识别 STL set::find == set::end 惯用法 →
    # `(c.find? k).isNone` / `.isSome`。形态：
    #   args[0] = Call("Array.find?" / "StdMap.find", [c, k])
    #   args[1] = Call(toList-ish, [c])  （end iterator 被映射到 toList 形态）
    if op_sym in ("==", "!=") and len(args) == 2:
        def _is_find_call(e):
            return (isinstance(e, Call) and isinstance(e.callee, str)
                    and e.callee in ("Array.find?", "Array.findVal",
                                       "StdMap.find"))
        def _is_end_marker(e):
            return (isinstance(e, Call) and isinstance(e.callee, str)
                    and (e.callee == "Array.toList" or e.callee.endswith(".toList")
                         or e.callee in ("StdMap.end", "Array.end")))
        find_call = end_call = None
        if _is_find_call(args[0]) and _is_end_marker(args[1]):
            find_call, end_call = args[0], args[1]
        elif _is_find_call(args[1]) and _is_end_marker(args[0]):
            find_call, end_call = args[1], args[0]
        if find_call is not None and end_call is not None:
            # 容器一致性检查：find?(c, k) 的 c 应等于 toList(c) 的 c
            if (len(find_call.args) >= 1 and len(end_call.args) >= 1
                    and _exprs_eq(find_call.args[0], end_call.args[0])):
                # `find == end` ⇔ `.isNone`；`find != end` ⇔ `.isSome`
                method = "Option.isNone" if op_sym == "==" else "Option.isSome"
                # 重新构造 find?：从 find_call 直接复用
                return Call(callee=method, args=[find_call], ty=BaseType.BOOL)

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
    # P0-10：StdMap.get! 返回 value type，填到 Call.ty 让下游 method 调用能查
    # CLASS_MAP（如 `groups[d].push_back(idx)` 中 push_back 的 receiver type）
    if op_sym == "[]" and len(args) == 2:
        recv_ty = _strip_ref(_expr_ty(args[0], typectx))
        if isinstance(recv_ty, StdMapType):
            return Call(callee="StdMap.get!", args=args, ty=recv_ty.value)
        if isinstance(recv_ty, NamedType) and recv_ty.name == "StdMap":
            return Call(callee="StdMap.get!", args=args, ty=None)
        # ArrayAccess for vector/Array
        elem_ty = None
        if isinstance(recv_ty, ArrayType):
            elem_ty = recv_ty.elem
        return ArrayAccess(arr=args[0], idx=args[1], ty=elem_ty)

    # operator() → lambda/functor apply
    if op_sym == "()" and len(args) >= 1:
        # 第 0 个是 callable：可能是 Var (lifted lambda 引用) / FieldAccess /
        # 或裸字符串。Pass 3 lambda_lift 后 lifted lambda 调用形态是
        # Call(UnresolvedOp("operator()"), [Var("dot"), arg1, arg2])。
        callee = args[0]
        if isinstance(callee, Var):
            # P0-3（轮 2 修复）：检查是否是 lifted lambda 别名，重映射到 mangled name
            alias_ty = typectx.get(callee.name)
            if isinstance(alias_ty, NamedType) and alias_ty.name.startswith("__lambda_alias__:"):
                callee_str = alias_ty.name[len("__lambda_alias__:"):]
            else:
                callee_str = callee.name
        elif isinstance(callee, str):
            callee_str = callee
        elif isinstance(callee, FieldAccess):
            callee_str = f"<field>.{callee.field_name}"
        else:
            callee_str = "<callable>"  # fallback：未知形态
        # P0-4（轮 2 修复）：CALL_OPERATOR_MAP fallback——若 callable 是 RNG/distribution
        # 等 functor，重映射到对应 functor.next/Rng.next 调用形态
        recv_ty = _strip_ref(_expr_ty(callee, typectx))
        if isinstance(recv_ty, NamedType):
            if recv_ty.name in ("UniformIntDist", "UniformRealDist", "Distribution"):
                # dist(rng) → Rng.next rng dist
                if len(args) == 2:
                    return Call(callee="Rng.next", args=[args[1], callee], ty=None)
        return Call(callee=callee_str, args=args[1:], ty=None)

    # operator->: 解引用
    if op_sym == "->" and len(args) == 1:
        return Call(callee="__deref__", args=[args[0]], ty=None)

    # 其他二元 / 一元
    recv = _typename_for_classmap(_expr_ty(args[0], typectx)) if args else None
    cls_ops = CLASS_MAP.get(recv, {}).get("operators", {}) if recv else {}

    if recv is None:
        # P0-9：args[0] 类型未知。对**数值/比较二元 op**，输出 BinOp 让 Lean
        # typeclass (HAdd/HMul/HSub/Eq/LT) 自动推断；对 args[0] 是 List/Array
        # 等容器形态的，保留 UnresolvedOp（Lean 没 List+Int typeclass，会真编译失败，
        # 标 gap 让人知道是 iterator-offset 模式 —— A 类设计 issue）
        _NUMERIC_OPS = {"+", "-", "*", "/", "%", "==", "!=", "<", ">", "<=", ">=",
                        "&&", "||", "&", "|", "^", "<<", ">>"}
        if op_sym in _NUMERIC_OPS and len(args) == 2:
            # 检查 args[0] 是不是 List 形态（来自 Array.toList）
            def _is_list_form(e: ExprIR) -> bool:
                if isinstance(e, Call):
                    cs = e.callee if isinstance(e.callee, str) else None
                    return bool(cs and ("toList" in cs or cs == "Array.toList"))
                return False
            if not _is_list_form(args[0]) and not _is_list_form(args[1]):
                # 真正的数值算术，让 Lean typeclass 推
                result_ty = (BaseType.BOOL
                             if op_sym in ("==", "!=", "<", ">", "<=", ">=", "&&", "||")
                             else None)
                return BinOp(op=op_sym, lhs=args[0], rhs=args[1], ty=result_ty)
        # iterator-offset / list-offset 等：保留 gap
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
    # Pattern A：`vec.erase(begin + j)` 早识别，避免 walk operator+ 误记 gap
    if method == "erase" and len(call.args) == 2:
        iter_info = _try_extract_iter_offset(call.args[1])
        if iter_info is not None:
            container, offset, op = iter_info
            if op == "+" and offset is not None:
                recv_raw = call.args[0]
                if _exprs_eq(container, recv_raw):
                    recv_w = _walk_expr(recv_raw, typectx, gap)
                    off_w = _walk_expr(offset, typectx, gap)
                    new_call = Call(callee="Array.eraseIdx",
                                    args=[recv_w, off_w], ty=None)
                    if (host_stmt_ctx
                            and isinstance(recv_w, (Var, FieldAccess, ArrayAccess))):
                        return ([AssignStmt(target=recv_w, value=new_call)], recv_w)
                    return new_call

    args = [_walk_expr(a, typectx, gap) for a in call.args]
    if not args:
        gap.method_miss.append(("?", method))
        return replace(call, args=args)

    recv = args[0]
    extra_args = args[1:]
    recv_ty = _expr_ty(recv, typectx)
    T = _typename_for_classmap(recv_ty)

    # B1 stage 4: 剥 "<method>." 前缀（Pass 1 把 C++ 转换运算符标记为
    # "<method>.operator bool" 形态，干扰 Pass 5 的 startswith 检测）
    if method.startswith("<method>."):
        method = method[len("<method>."):]
    if T is None:
        # T=None 且 method 是 unary 转换 / binary operator → fallback
        if method.startswith("operator "):
            op_sym = method[len("operator "):]
            _BINARY_OPS = {"==", "!=", "<", ">", "<=", ">=", "+", "-", "*",
                           "/", "%", "&&", "||", "&", "|", "^", "<<", ">>"}
            if op_sym in _BINARY_OPS and len(args) >= 2:
                return BinOp(op=op_sym, lhs=args[0], rhs=args[1],
                             ty=BaseType.BOOL if op_sym in ("==", "!=", "<",
                                 ">", "<=", ">=", "&&", "||") else None)
            # operator bool / ! 等 unary：单 arg 时 fallback UnaryOp
            if len(args) >= 1 and op_sym in ("bool", "!"):
                return UnaryOp(op=op_sym, operand=args[0],
                               ty=BaseType.BOOL if op_sym == "bool" else None)
        gap.method_miss.append(("None", method))
        return replace(call, args=args)

    cls = CLASS_MAP.get(T, {})
    methods = cls.get("methods", {})
    # P0-7：对 BaseType.BOOL 上的 `operator bool`，identity 剥（bool→bool 是 noop）
    if T == "BaseType.BOOL" and method == "operator bool":
        return recv
    # P1-3：`<method>.operator bool` / `<method>.operator <op>` 是 Pass 1 把
    # C++ 转换运算符或运算符 overload 包成 method 调用形态。优先查 methods，
    # 找不到时把 method 名按 operator 处理（剥 "operator " 前缀查 operators）
    if method not in methods and method.startswith("operator "):
        op_sym = method[len("operator "):]
        operators = cls.get("operators", {})
        if op_sym in operators:
            target = operators[op_sym]
            if target is None:
                # typeclass 形态：unary (bool/!) → UnaryOp / 直接返回 recv；
                # binary (==/!=/<...) → BinOp(op, recv, extra_args[0])
                _BINARY_OPS = {"==", "!=", "<", ">", "<=", ">=", "+", "-",
                               "*", "/", "%", "&&", "||", "&", "|", "^",
                               "<<", ">>"}
                if op_sym in _BINARY_OPS and len(extra_args) >= 1:
                    return BinOp(op=op_sym, lhs=recv, rhs=extra_args[0],
                                 ty=BaseType.BOOL if op_sym in ("==", "!=", "<",
                                     ">", "<=", ">=", "&&", "||") else None)
                return UnaryOp(op=op_sym, operand=recv,
                               ty=BaseType.BOOL if op_sym == "bool" else None)
            if isinstance(target, str):
                return Call(callee=target, args=[recv] + extra_args, ty=None)
    if method not in methods:
        # B1 stage 4：methods 找不到 + method 是 operator → fallback Bin/UnaryOp
        if method.startswith("operator "):
            op_sym = method[len("operator "):]
            _BINARY_OPS = {"==", "!=", "<", ">", "<=", ">=", "+", "-", "*",
                           "/", "%", "&&", "||", "&", "|", "^", "<<", ">>"}
            if op_sym in _BINARY_OPS and len(extra_args) >= 1:
                return BinOp(op=op_sym, lhs=recv, rhs=extra_args[0],
                             ty=BaseType.BOOL if op_sym in ("==", "!=", "<",
                                 ">", "<=", ">=", "&&", "||") else None)
            if op_sym in ("bool", "!"):
                return UnaryOp(op=op_sym, operand=recv,
                               ty=BaseType.BOOL if op_sym == "bool" else None)
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
        # Pattern A：识别 `vec.erase(begin + j)` → `Array.eraseIdx(vec, j)`
        # 形态：method == "erase" / target == "Array.erase" / 1 个 extra arg
        # 是 Iterator.fromList(begin+offset) 形态
        if (method == "erase" and len(extra_args) == 1):
            iter_info = _try_extract_iter_offset(extra_args[0])
            if iter_info is not None:
                container, offset, op = iter_info
                if op == "+" and offset is not None and _exprs_eq(container, recv):
                    new_call = Call(callee="Array.eraseIdx",
                                    args=[recv, offset], ty=None)
                    if host_stmt_ctx and isinstance(recv, (Var, FieldAccess, ArrayAccess)):
                        return ([AssignStmt(target=recv, value=new_call)], recv)
                    return new_call

        new_call = Call(callee=target, args=[recv] + extra_args, ty=None)
        # P1-5：扩展 ArrayAccess 也作 target（如 `groups[d].push_back(x)`）
        if host_stmt_ctx and _is_lvalue_target(recv):
            return ([AssignStmt(target=recv, value=new_call)], recv)
        return new_call

    if disposition == "mutate_push":
        new_call = Call(callee=target or "Array.push",
                        args=[recv] + extra_args, ty=None)
        if host_stmt_ctx and _is_lvalue_target(recv):
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
    typename = _normalize_typename(op_name)
    arity = len(call.args)

    # Pattern B：识别 `vector<T>(begin±i, end±j)` 形态。**raw args 检测**（在 walk
    # 之前），匹配时只 walk container + offset 子表达式，避免 walk outer
    # operator+/- 时误记 None-receiver gap。
    if ("vector<" in typename and arity >= 2
            and op_name.startswith("construct_")):
        a0_info = _try_extract_iter_offset(call.args[0])
        a1_info = _try_extract_iter_offset(call.args[1])
        if a0_info is not None and a1_info is not None:
            c0, off0, op0 = a0_info
            c1, off1, op1 = a1_info
            if _exprs_eq(c0, c1):
                cw = _walk_expr(c0, typectx, gap)
                off0_w = _walk_expr(off0, typectx, gap) if off0 is not None else None
                off1_w = _walk_expr(off1, typectx, gap) if off1 is not None else None
                if op0 == "+" and op1 == "begin" and off0_w is not None:
                    # P0-10（轮 2 修复）：Lean Array.drop 第二参 Nat。包 .toNat
                    nat_off = Call(callee="__cast__({x}).toNat",
                                   args=[off0_w], ty=BaseType.NAT)
                    return Call(callee="Array.drop", args=[cw, nat_off],
                                ty=NamedType("Array"))
                if op0 == "begin" and op1 == "-" and off1_w is not None:
                    # P0-9（轮 2 修复）：Lean Array.dropLast 是 0-arg "去掉最后 1 个"，
                    # 没有 2-arg 形态。`vector(begin, end-k)` = `arr.take (size-k)`
                    size_call = Call(callee="Array.size", args=[cw], ty=BaseType.NAT)
                    take_count = BinOp(op="-", lhs=size_call, rhs=off1_w, ty=BaseType.NAT)
                    return Call(callee="Array.take", args=[cw, take_count],
                                ty=NamedType("Array"))
                if op0 == "begin" and op1 == "begin":
                    return Call(callee="id", args=[cw], ty=NamedType("Array"))
                if op0 == "begin" and op1 == "+" and off1_w is not None:
                    nat_off = Call(callee="__cast__({x}).toNat",
                                   args=[off1_w], ty=BaseType.NAT)
                    return Call(callee="Array.take", args=[cw, nat_off],
                                ty=NamedType("Array"))
                if op0 == "+" and op1 == "+" and off0_w and off1_w:
                    inner = Call(callee="Array.take", args=[cw, off1_w],
                                 ty=NamedType("Array"))
                    return Call(callee="Array.drop", args=[inner, off0_w],
                                ty=NamedType("Array"))

    args = [_walk_expr(a, typectx, gap) for a in call.args]

    # P0-1（轮 2 修复）：copy ctor 检测 — 1-arg 构造且 arg 类型 == 目标类型时
    # 是 copy/identity，直接返回 arg。避免 `<construct_Zp>(Zp_value)` 被错误
    # 包成 `Zp.ofInt(Zp_value)` 等数据表 1-arg fallback。
    target_ty = _typename_to_typeir(typename)
    if arity == 1 and target_ty is not None:
        arg_ty = _strip_ref(_expr_ty(args[0], typectx))
        if (arg_ty is not None
                and canonicalize_type(target_ty) == canonicalize_type(arg_ty)
                and not _is_unresolved(arg_ty)):
            return args[0]

    resolution = resolve_constructor(op_name, arity)
    if resolution is None:
        gap.constructor_miss.append((op_name, arity))
        return replace(call, args=args)
    # P0-3：从 typename 推回目标类型，填到 Call.ty。这样下游 receiver 类型查询
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
        # A 方案：参数化 Factorization。inner 解析参考 Pass 1。
        inner = name[len("factorization<"):].rstrip(">").strip()
        if inner.startswith("clpoly::"):
            inner = inner[len("clpoly::"):]
        if "upolynomial_<ZZ>" in inner:
            return NamedType("Factorization SparsePolyZZ")
        if "upolynomial_<Zp>" in inner:
            return NamedType("Factorization SparsePolyZp")
        if "Zp" in inner and ("polynomial_" in inner or inner == "PolyZp"):
            return NamedType("Factorization MvPolyZp")
        return NamedType("Factorization MvPolyZZ")
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
        # P0-7（轮 2 修复）：内联 LambdaExpr 来自 Pass 4 filter-loop pred（Pass 3
        # lambda_lift 已经把所有 outer-level lambda 提升到 aux_lambdas，但 Pass 4
        # 在 _build_pred_lambda 里**新生成**内联 LambdaExpr 用作 Array.filter 的
        # pred，body 内的 UnresolvedOp 需要 Pass 5 处理）
        inner_typectx = dict(typectx)
        for p in e.params:
            if p.ty is not None and not _is_unresolved(p.ty):
                inner_typectx[p.name] = p.ty
        new_body = _walk_stmts(list(e.body), inner_typectx, gap)
        return replace(e, body=new_body)
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
    """从表达式收集 UB require（先深扫子表达式，再当前节点）。
    P0-5：吸收 _CAST_PENDING_REQUIRES（_resolve_cast 留下的 cast UB）"""
    reqs: list[RequireStmt] = []
    # 先吸收 cast pending requires
    if _CAST_PENDING_REQUIRES:
        reqs.extend(_CAST_PENDING_REQUIRES)
        _CAST_PENDING_REQUIRES.clear()

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
        # P0-3（轮 2 修复）：检测 `let X := Var("_lambda_<host>_<N>")` 形态——
        # X 是 lifted lambda 的本地别名（如 `auto dot = ...; dot(a, b)`）。
        # 标记 X 类型为 `NamedType("__lambda_alias__:<mangled>")`，op_sym=="()"
        # 解析时用 mangled name 作 callee 而非 X 自身
        if (new_val is not None and isinstance(new_val, Var)
                and new_val.name.startswith("_lambda_")):
            typectx[s.var.name] = NamedType(f"__lambda_alias__:{new_val.name}")
        elif s.ty is not None and not _is_unresolved(s.ty):
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
        # P0-8（轮 2 修复）：`StdMap.get!(m, k) = v` 形态——Lean 不允许 Call 作
        # 左值。识别 LHS 是 `StdMap.get!` 的 Call → 转 `m = StdMap.insert m k v`
        if (isinstance(new_tgt, Call) and isinstance(new_tgt.callee, str)
                and new_tgt.callee == "StdMap.get!" and len(new_tgt.args) == 2):
            m, k = new_tgt.args
            insert_call = Call(callee="StdMap.insert", args=[m, k, new_val], ty=None)
            reqs = (_collect_ub_requires(new_val, typectx)
                    + _collect_ub_requires(insert_call, typectx))
            if isinstance(m, (Var, FieldAccess, ArrayAccess)):
                return reqs + [AssignStmt(target=m, value=insert_call)]
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
                             body=new_body, decomposition=s.decomposition,
                             is_mutable_ref=s.is_mutable_ref)]

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
        # CF-3 阶段 3：STL 算法 sort/iota 取 begin/end iterator 范围 →
        # 转为 functional Array.sort / Array.range_init 形态
        # 形态：sort(Array.toList(c), Array.toList(c), comp) → c := Array.sort(c, comp)
        #       iota(Array.toList(c), Array.toList(c), v)   → c := Array.range_init(c.size, v)
        if isinstance(e, Call) and isinstance(e.callee, str) \
                and e.callee in ("sort", "iota") and len(e.args) == 3:
            container_arg = _stl_extract_toList_container(e.args[0])
            container_arg2 = _stl_extract_toList_container(e.args[1])
            if container_arg is not None and container_arg2 is not None \
                    and _exprs_eq(container_arg, container_arg2):
                container = _walk_expr(container_arg, typectx, gap)
                if isinstance(container, (Var, FieldAccess, ArrayAccess)):
                    third_arg = _walk_expr(e.args[2], typectx, gap)
                    if e.callee == "sort":
                        new_assign = AssignStmt(
                            target=container,
                            value=Call(callee="Array.sort",
                                       args=[container, third_arg],
                                       ty=_expr_ty(container, typectx)),
                        )
                    else:  # iota
                        new_assign = AssignStmt(
                            target=container,
                            value=Call(callee="Array.range_init",
                                       args=[container, third_arg],
                                       ty=_expr_ty(container, typectx)),
                        )
                    return [new_assign]
        # B1 修复：内建数值类型的 ++/--，Clang AST 直接吐出 UnaryOp，
        # 不会包装成 Call(UnresolvedOp("operator++"))。在 stmt 级也展开为 AssignStmt
        if isinstance(e, UnaryOp) and e.op in ("++", "--"):
            target = _walk_expr(e.operand, typectx, gap)
            if isinstance(target, (Var, FieldAccess, ArrayAccess)):
                target_ty = _strip_ref(_expr_ty(target, typectx))
                op_sym = "+" if e.op == "++" else "-"
                new_assign = AssignStmt(
                    target=target,
                    value=BinOp(op=op_sym, lhs=target,
                                rhs=Lit(value=1, ty=BaseType.INT64),
                                ty=target_ty),
                )
                return _collect_ub_requires(new_assign.value, typectx) + [new_assign]
        if isinstance(e, Call) and isinstance(e.callee, UnresolvedOp):
            op = e.callee.op_name
            # operator++ / operator-- (incdec)
            if op in ("operator++", "operator--") and len(e.args) >= 1:
                target = _walk_expr(e.args[0], typectx, gap)
                if isinstance(target, (Var, FieldAccess, ArrayAccess)):
                    target_ty = _strip_ref(_expr_ty(target, typectx))
                    # P0-2（轮 2 修复）：先查 CLASS_MAP[T]["operators"]——
                    # Iterator 等 user-defined 类型的 ++/-- 应走 mutate 路径
                    # （如 `it = Iterator.advance(it)`）而非硬编码 `it = it + 1`
                    T = _typename_for_classmap(target_ty)
                    cls_ops = CLASS_MAP.get(T, {}).get("operators", {}) if T else {}
                    op_short = "++" if op == "operator++" else "--"
                    if op_short in cls_ops:
                        target_def = cls_ops[op_short]
                        if isinstance(target_def, tuple) and len(target_def) == 2:
                            disp, fn = target_def
                            if disp in ("mutate", "method"):
                                new_call = Call(callee=fn, args=[target], ty=target_ty)
                                return _collect_ub_requires(new_call, typectx) + [
                                    AssignStmt(target=target, value=new_call)
                                ]
                        elif isinstance(target_def, str):
                            new_call = Call(callee=target_def, args=[target], ty=target_ty)
                            return _collect_ub_requires(new_call, typectx) + [
                                AssignStmt(target=target, value=new_call)
                            ]
                    # 默认（数值类型）：x = x ± 1
                    op_sym = "+" if op == "operator++" else "-"
                    new_assign = AssignStmt(
                        target=target,
                        value=BinOp(op=op_sym, lhs=target,
                                    rhs=Lit(value=1, ty=BaseType.INT64),
                                    ty=target_ty),
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
        # 形式清理：ExprStmt 主体是裸 Var/Lit（C++ `(void)x;` 或 reserve()
        # 等 noop disposition 退化）— 真合法 dead code，丢弃避免 sideeff 噪声
        if isinstance(new_expr, (Var, Lit)):
            return []
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
