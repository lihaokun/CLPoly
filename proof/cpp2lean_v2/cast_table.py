"""
CAST_TABLE：C++ Cast 到 Lean 表达式的映射（Pass 5 operator_resolve 用）。

依据 `docs/design/l1-translation-validation/survey/type-system.md §8`。

覆盖范围（5670 cast 全部）：
  - NoOp / LValueToRValue / ToVoid / BuiltinFnToFnPtr / ArrayToPointerDecay /
    FunctionToPointerDecay — 恒等剥（无 require）
  - IntegralCast（656 处实测，19 种 src/dst 对）
  - IntegralToFloating（4 处，全在 __heuristic_starting_precision）
  - FloatingToIntegral（1 处）
  - ConstructorConversion / UserDefinedConversion → 不在此查表，由
    Pass 5 查 FUNC_MAP / CLASS_MAP 构造器

查询约定：`_resolve_cast(cast) -> CastResolution | None`
  - 命中 → `CastResolution(template="({x}).toInt64", ub_kind="fits_int64" | None)`
  - 未命中 → None（Pass 5 走 B 策略：保留原 Cast，记 gap）

`template` 约定：`{x}` 占位源表达式；Pass 5 替换为实际 Lean 源码
`ub_kind`：若非 None，Pass 5 生成对应 `RequireStmt`（见 UB_REQUIRE_BUILDERS）
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from ir_types import BaseType, NamedType, UnknownType, TypeIR


# ============================================================
# CastResolution
# ============================================================

@dataclass(frozen=True)
class CastResolution:
    """Cast 解析结果。"""
    template: str         # Lean 表达式模板，`{x}` 是源值占位符
    ub_kind: Optional[str]  # None 表示无 UB；否则为 UB_REQUIRE_BUILDERS 的键


# ============================================================
# 类型规范化：TypeIR → 规范化字符串键
# ============================================================

def canonicalize_type(ty: TypeIR | None) -> str:
    """规范化 TypeIR 为 CAST_TABLE 键字符串。

    - BaseType.X → 'int32' / 'int64' / 'uint32' / 'uint64' / 'nat' / 'bool' / 'float64' / ...
    - NamedType("ValueType"/"DependentType"/"ResultType"/"StlInternal") → 'unresolved'
      （Pass 1 未穿透的 C++ 模板依赖名；Pass 5 标 gap）
    - NamedType(其他) → 'named:<name>'
    - UnknownType → 'unknown'
    - None → 'none'
    """
    if ty is None:
        return "none"
    if isinstance(ty, BaseType):
        mapping = {
            BaseType.INT32: "int32",
            BaseType.INT64: "int64",
            BaseType.UINT32: "uint32",
            BaseType.UINT64: "uint64",
            BaseType.UINT128: "uint128",
            BaseType.NAT: "nat",
            BaseType.BOOL: "bool",
            BaseType.FLOAT: "float64",
            BaseType.UNIT: "unit",
        }
        return mapping.get(ty, f"base_other:{ty.name}")
    if isinstance(ty, NamedType):
        if ty.name in ("ValueType", "DependentType", "ResultType", "StlInternal"):
            return "unresolved"
        return f"named:{ty.name}"
    if isinstance(ty, UnknownType):
        # Pass 1 对 `std::size_t` / `std::size_type` 有时保留为 UnknownType（
        # typedef 未穿透）。CLPoly 里它们一律是 unsigned 容器索引，规范化为 'nat'。
        raw = getattr(ty, "raw", "") or ""
        if raw in ("std::size_t", "size_t", "std::size_type", "size_type"):
            return "nat"
        return "unknown"
    return f"other:{type(ty).__name__}"


# ============================================================
# CAST_TABLE：(src_canon, dst_canon) → CastResolution
# ============================================================

# 术语：
#   - 'nat'   = Lean 的 `Nat`（Pass 1 把 C++ `size_t`/`size_type` 已映射到
#               BaseType.NAT，简化了 Pass 5 的转换；L1 严格语义是 `USize`
#               但 Pass 1 做了化简，视为 L1→L2 过渡）
#   - 'int32' = Lean `Int32`（L1 wrapping）
#   - 'int64' = Lean `Int64`
#   - 'uint64' = Lean `UInt64`
CAST_TABLE: dict[tuple[str, str], CastResolution] = {
    # ====================================================================
    # IntegralCast — 19 条实测 + 设计补全
    # ====================================================================

    # 1. 同型 noop（4 处 uint64→uint64，1 处 ResultType→ResultType 等）
    ("int32",  "int32"):   CastResolution("{x}", None),
    ("int64",  "int64"):   CastResolution("{x}", None),
    ("uint32", "uint32"):  CastResolution("{x}", None),
    ("uint64", "uint64"):  CastResolution("{x}", None),
    ("nat",    "nat"):     CastResolution("{x}", None),
    ("bool",   "bool"):    CastResolution("{x}", None),

    # 2. Widening（无 UB）
    ("int32",  "int64"):   CastResolution("({x}).toInt64", None),            # 70 次
    ("uint32", "uint64"):  CastResolution("({x}).toUInt64", None),           # 2 次
    ("int32",  "nat"):     CastResolution("({x}).toNat", "nonneg"),          # 452 次（require x ≥ 0）
    ("int32",  "uint64"):  CastResolution("({x}).toNat.toUInt64", "nonneg"), # 41 次
    # int → size_t 走和 int → nat 相同路径（CLPoly Pass 1 size_t 常是 NAT；
    # UnknownType("std::size_t") 出现 4 次 — fallback 到 int32→nat 含义）

    # 3. Narrowing（带 UB-7 / fits require）
    ("nat",    "int32"):   CastResolution("({x}).toInt32", "fits_int32"),    # 46 次 — require x ≤ INT32_MAX
    ("int64",  "int32"):   CastResolution("({x}).toInt32", "fits_int32"),    # 6 次  — require -2^31 ≤ x < 2^31
    ("int64",  "uint64"):  CastResolution("({x}).toNat.toUInt64", "nonneg"), # 5 次  — require x ≥ 0
    ("uint64", "int64"):   CastResolution("({x}).toInt64", "fits_int64"),    # 2 次  — UB-7
    ("uint64", "int32"):   CastResolution("({x}).toInt32", "fits_int32"),    # 1 次  — UB-7

    # 4. bool 参与（设计需要，实测 0 次但 §8.2 列出）
    ("bool",   "int32"):   CastResolution("(if {x} then 1 else 0)", None),
    ("int32",  "bool"):    CastResolution("({x} != 0)", None),

    # ====================================================================
    # IntegralToFloating（4 处实测，__heuristic_starting_precision）
    # ====================================================================
    ("int32",  "float64"): CastResolution("Int.toFloat ({x})", None),
    ("uint64", "float64"): CastResolution("Nat.toFloat (UInt64.toNat {x})", None),

    # ====================================================================
    # FloatingToIntegral（1 处实测）
    # ====================================================================
    ("float64", "int32"):  CastResolution("Float.toInt32 ({x})", None),  # C++ truncate 语义
}


# ============================================================
# UB require 生成模板（Pass 5 用）
# ============================================================

# 从 ub_kind 到 "Lean 前置条件表达式模板"；Pass 5 用实际值填 `{x}` 或 `{dst}`
UB_REQUIRE_TEMPLATES: dict[str, str] = {
    "nonneg":      "{x} ≥ 0",                    # 非负（int → unsigned 类）
    "fits_int32":  "-2147483648 ≤ {x} ∧ {x} ≤ 2147483647",  # Int32 范围
    "fits_int64":  "-9223372036854775808 ≤ {x} ∧ {x} ≤ 9223372036854775807",  # Int64 范围
    "fits_uint64": "0 ≤ {x} ∧ {x} ≤ 18446744073709551615",  # UInt64 范围
}


# ============================================================
# Cast kind → disposition
# ============================================================

# Pass 5 按 cast_kind 分派：
#   "strip":    直接剥 Cast，返回内部 expr（不查 CAST_TABLE）
#   "table":    查 CAST_TABLE[(src, dst)]
#   "delegate": 委托 FUNC_MAP/CLASS_MAP 构造器（ConstructorConversion、
#               UserDefinedConversion）
#   "unknown":  其他 kind —— Pass 5 按 B 策略保留原 Cast
CAST_KIND_DISPOSITION: dict[str, str] = {
    "NoOp":                   "strip",
    "LValueToRValue":         "strip",
    "ToVoid":                 "strip",
    "BuiltinFnToFnPtr":       "strip",
    "ArrayToPointerDecay":    "strip",
    "FunctionToPointerDecay": "strip",

    "IntegralCast":           "table",
    "IntegralToFloating":     "table",
    "FloatingToIntegral":     "table",

    "ConstructorConversion":  "delegate",
    "UserDefinedConversion":  "delegate",
}


def lookup_cast(src: TypeIR | None, dst: TypeIR | None) -> CastResolution | None:
    """在 CAST_TABLE 中查找，返回 None 表示未命中（含类型规范化失败）。"""
    s = canonicalize_type(src)
    d = canonicalize_type(dst)
    return CAST_TABLE.get((s, d))


# ============================================================
# 字面量整数类型范围（用于 Lit cast 简化的安全约束）
# ============================================================

_INT_TYPE_RANGES: dict[str, tuple[int, int]] = {
    "int8":   (-(2**7),  2**7 - 1),
    "uint8":  (0,        2**8 - 1),
    "int16":  (-(2**15), 2**15 - 1),
    "uint16": (0,        2**16 - 1),
    "int32":  (-(2**31), 2**31 - 1),
    "uint32": (0,        2**32 - 1),
    "int64":  (-(2**63), 2**63 - 1),
    "uint64": (0,        2**64 - 1),
}

# int16 范围：保守"通用安全"区间。任何字面量在此范围内，必然能容纳到
# 主流整数类型（uint8 略小但 CLPoly 不用；int16/+ 全部容纳；负数被 uint8/16 拒
# 但本约束不允许负数走 uint8/16）。
_INT16_RANGE = _INT_TYPE_RANGES["int16"]


def is_safe_to_strip_lit_cast(cast) -> bool:
    """判定 `Cast(Lit(int_val), source, target, IntegralCast)` 是否能安全剥掉。

    语义保证：剥 Cast 后的 Lit 在 Lean 端通过 OfNat typeclass 自动 elaborate
    到目标类型，得到与 C++ implicit conversion 相同的数值——前提条件：
      1. 内部表达式是 `Lit`，且 `value` 是 Python `int`（不是 float/bool/string）
      2. **`value ∈ int16 范围 [-32768, 32767]`**（"足够小"的通用安全区间）
      3. 若 `target_ty` 已知（不是 unresolved），`value` 必须在 target 类型范围内；
         若 target 未知（DependentType），仅靠条件 2 保证

    覆盖场景：CLPoly 的 4 个 `Lit(0|1) → DependentType` miss——0 和 1 远在 int16
    内，且 target 通常是 UInt64/Int64/Int32（任何主流整数都容纳 0/1）。

    超出条件的 cast 保留为 Pass 5 的 fallback（生成 sorry / 保留 UnresolvedOp）。
    """
    # 延迟导入避免循环
    from ir_types import Cast, Lit
    if not (isinstance(cast, Cast) and isinstance(cast.expr, Lit)):
        return False
    if cast.cast_kind != "IntegralCast":
        return False
    val = cast.expr.value
    if not isinstance(val, int) or isinstance(val, bool):
        return False

    # 条件 1+2：int16 范围保险阀
    if not (_INT16_RANGE[0] <= val <= _INT16_RANGE[1]):
        return False

    # 条件 3：若 target 已知，必须在 target 范围内
    target_canon = canonicalize_type(cast.target_ty)
    if target_canon in _INT_TYPE_RANGES:
        lo, hi = _INT_TYPE_RANGES[target_canon]
        if not (lo <= val <= hi):
            return False
    elif target_canon == "nat":
        if val < 0:
            return False
    # target 是 unresolved/unknown/named/...：靠条件 2 兜底

    return True
