"""
CONSTRUCTOR_MAP：C++ 构造器调用 (`Call(UnresolvedOp("construct_X"), args)`) 到
Lean 表达式的映射（Pass 5 operator_resolve 用）。

依据 `docs/design/l1-translation-validation/survey/type-system.md §3 / §8.5`
+ HIR₃ 实证扫描（631 calls / 108 (name, arity) 对）。

**两层结构**：

1. **`CLPOLY_CONSTRUCTORS`** — 按 typename 显式枚举 CLPoly 核心类型 +
   aux struct + Iterator。每条 entry 按 arity 索引到 Lean 模板。
2. **`STL_CONSTRUCTOR_PATTERNS`** — 规则匹配 STL 容器 `std::vector<T>` /
   `std::pair<A,B>` / `std::map<K,V>` / `std::set<T>` / RNG 三件套。规则
   覆盖所有泛型实例，避免 60+ 条冗余枚举。

**resolve API**：`resolve_constructor(op_name, arity) -> ConstructorResolution | None`

未命中 → 返回 None，Pass 5 走 B 策略保留 UnresolvedOp。
"""

from __future__ import annotations
import re
from dataclasses import dataclass
from typing import Optional


# ============================================================
# ConstructorResolution
# ============================================================

@dataclass(frozen=True)
class ConstructorResolution:
    """构造器解析结果。"""
    template: str          # Lean 表达式模板，{a0}/{a1}/... 是位置参数占位
    is_default: bool = False  # arity=0 的"空构造"


# ============================================================
# 类型名规范化（剥 const、clpoly::、std::、模板内嵌套前缀）
# ============================================================

def _normalize_typename(name: str) -> str:
    """剥 'construct_' 前缀 + 'const ' / 'const_' / 'clpoly::' / 'std::' 等装饰。

    返回 normalized typename string（保留模板括号 `<...>`）。
    """
    s = name
    if s.startswith("construct_"):
        s = s[len("construct_"):]
    # 剥 'const '
    if s.startswith("const "):
        s = s[len("const "):]
    if s.startswith("const_"):
        s = s[len("const_"):]
    # 剥 'clpoly::'
    s = s.replace("clpoly::", "").replace("std::", "")
    # 剥 'typename '（如 `typename Poly::monomial_type`）
    if s.startswith("typename "):
        s = s[len("typename "):]
    return s.strip()


# ============================================================
# CLPoly 核心类型构造器（按 arity 索引）
# ============================================================

CLPOLY_CONSTRUCTORS: dict[str, dict[int, ConstructorResolution]] = {
    # ZZ：任意精度整数
    "ZZ": {
        0: ConstructorResolution("(0 : Int)", is_default=True),
        1: ConstructorResolution("(({a0}) : Int)"),  # ZZ(int_lit) / ZZ(uint64) etc.
    },
    # Zp：模 p 整数
    "Zp": {
        0: ConstructorResolution("Zp.default", is_default=True),
        1: ConstructorResolution("Zp.ofInt {a0}"),
        # C++ Zp(int64_t v, uint64_t p) 做模归约 → 用 Zp.ofInt（Lean Model 提供）
        # ofInt 签名 (Int, UInt64)，Int64 → Int 用 .toInt
        2: ConstructorResolution("Zp.ofInt ({a0}).toInt {a1}"),
    },
    # QQ：有理数
    "QQ": {
        0: ConstructorResolution("(0 : Rat)", is_default=True),
        1: ConstructorResolution("QQ.ofInt {a0}"),
        2: ConstructorResolution("QQ.mk {a0} {a1}"),  # (num, den)
    },

    # umonomial：一元单项式（degree 类型）
    "umonomial": {
        0: ConstructorResolution("UMonomial.default", is_default=True),
        1: ConstructorResolution("UMonomial.mk {a0}"),
    },
    "UMonomial": {
        0: ConstructorResolution("UMonomial.default", is_default=True),
        1: ConstructorResolution("UMonomial.mk {a0}"),
    },

    # basic_monomial<lex_<...>>：多元单项式
    "basic_monomial<lex_<less>>": {
        0: ConstructorResolution("MvMonomial.empty", is_default=True),
        1: ConstructorResolution("MvMonomial.mk {a0}"),
    },
    # 别名（typename Poly::monomial_type）
    "Monomial": {
        0: ConstructorResolution("MvMonomial.empty", is_default=True),
        1: ConstructorResolution("MvMonomial.mk {a0}"),
    },

    # variable / Variable
    "variable": {
        0: ConstructorResolution("Variable.default", is_default=True),
        1: ConstructorResolution("Variable.mk {a0}"),
    },
    "Variable": {
        0: ConstructorResolution("Variable.default", is_default=True),
        1: ConstructorResolution("Variable.mk {a0}"),
    },

    # upolynomial_<T>（一元多项式 typedef）
    # upolynomial_<ZZ> 在 HIR 里通常是 SparsePolyZZ；upolynomial_<Zp> 是 SparsePolyZp
    "upolynomial_<ZZ>": {
        0: ConstructorResolution("SparsePolyZZ.empty", is_default=True),
        1: ConstructorResolution("({a0} : SparsePolyZZ)"),  # copy / convert
    },
    "upolynomial_<Zp>": {
        0: ConstructorResolution("SparsePolyZp.empty", is_default=True),
        1: ConstructorResolution("({a0} : SparsePolyZp)"),
    },
    "UPZp": {
        0: ConstructorResolution("SparsePolyZp.empty", is_default=True),
        1: ConstructorResolution("({a0} : SparsePolyZp)"),
    },

    # polynomial_<T, lex_<...>>（多元多项式）
    "polynomial_<ZZ, lex_<less>>": {
        0: ConstructorResolution("MvPolyZZ.empty", is_default=True),
        1: ConstructorResolution("MvPolyZZ.mk {a0}"),  # (comp_ptr) 这里假设 1-arg 是 comp_ptr
    },
    "polynomial_<Zp, lex_<less>>": {
        0: ConstructorResolution("MvPolyZp.empty", is_default=True),
        1: ConstructorResolution("MvPolyZp.mk {a0}"),
    },
    "polynomial_<ZZ, lex>": {  # 不带 less 后缀的实例化
        0: ConstructorResolution("MvPolyZZ.empty", is_default=True),
        1: ConstructorResolution("MvPolyZZ.mk {a0}"),
    },

    # CLPoly 别名
    "Poly": {
        0: ConstructorResolution("MvPolyZZ.empty", is_default=True),
        1: ConstructorResolution("MvPolyZZ.mk {a0}"),
    },
    "PolyZp": {
        0: ConstructorResolution("MvPolyZp.empty", is_default=True),
        1: ConstructorResolution("MvPolyZp.mk {a0}"),
    },
    "PolyZZ": {
        0: ConstructorResolution("MvPolyZZ.empty", is_default=True),
        1: ConstructorResolution("MvPolyZZ.mk {a0}"),
    },

    # value_type 占位（`construct_value_type(a, b)` —— 通常是 pair-like 初始化）
    "value_type": {
        2: ConstructorResolution("({a0}, {a1})"),  # pair 形式默认
    },

    # factorization<T>
    "factorization<Poly>": {
        0: ConstructorResolution("Factorization.empty", is_default=True),
        1: ConstructorResolution("({a0} : Factorization)"),
    },
    "factorization<polynomial_<ZZ, lex_<less>>>": {
        0: ConstructorResolution("Factorization.empty", is_default=True),
        1: ConstructorResolution("({a0} : Factorization)"),
    },
    "factorization<polynomial_<ZZ, grlex_<less>>>": {
        0: ConstructorResolution("Factorization.empty", is_default=True),
        1: ConstructorResolution("({a0} : Factorization)"),
    },
    "factorization<upolynomial_<ZZ>>": {
        0: ConstructorResolution("Factorization.empty", is_default=True),
        1: ConstructorResolution("({a0} : Factorization)"),
    },

    # Iterator（compact-erase 残留 + Pass 4 漏识别兜底）
    "iterator": {
        1: ConstructorResolution("Iterator.fromList {a0}"),
    },
    "const_iterator": {
        1: ConstructorResolution("Iterator.fromList {a0}"),
    },

    # LLLMatrix
    "LLLMatrix": {
        0: ConstructorResolution("LLLMatrix.empty", is_default=True),
        1: ConstructorResolution("LLLMatrix.mk {a0}"),
        3: ConstructorResolution("LLLMatrix.replicate {a0} {a1} {a2}"),  # (n, m, default)
    },

    # Aux structs（CLPoly 内部）
    "__wang_lc_result<less>": {
        0: ConstructorResolution("WangLcResult.default", is_default=True),
        1: ConstructorResolution("WangLcResult.mk {a0}"),
    },
    "__prime_selection_result": {
        0: ConstructorResolution("PrimeSelectionResult.default", is_default=True),
        1: ConstructorResolution("PrimeSelectionResult.mk {a0}"),
    },
    "__hensel_node": {
        0: ConstructorResolution("HenselNode.default", is_default=True),
        1: ConstructorResolution("HenselNode.mk {a0}"),
    },

    # 复合 const/uless/grlex 长形态（HIR 里偶现）
    "basic_polynomial<umonomial, ZZ, uless>": {
        1: ConstructorResolution("({a0} : SparsePolyZZ)"),
    },
    "basic_polynomial<umonomial, Zp, uless>": {
        1: ConstructorResolution("({a0} : SparsePolyZp)"),
    },

    # Dependent nested typedefs: `typename T::monomial_type` 等。Pass 1 剥
    # `typename ` 前缀后剩 `T::monomial_type`，typedef 语义已知映射到具体
    # CLPoly 类型。
    "Poly::monomial_type": {
        0: ConstructorResolution("MvMonomial.empty", is_default=True),
    },
    "PolyZp::monomial_type": {
        0: ConstructorResolution("MvMonomial.empty", is_default=True),
    },
    "PolyZZ::monomial_type": {
        0: ConstructorResolution("MvMonomial.empty", is_default=True),
    },

    # std::string —— 仅用于 assert 消息字符串（Variable name 字面量等）
    # 2 args: (char_array, allocator) → String.mk
    "basic_string<char>": {
        2: ConstructorResolution("String.mk {a0}"),
    },
    "string": {
        0: ConstructorResolution("\"\""),
        1: ConstructorResolution("({a0} : String)"),
        2: ConstructorResolution("String.mk {a0}"),
    },
}


# ============================================================
# STL 容器构造器规则（pattern matching）
# ============================================================

# `std::vector<T>(args...)`：
# - arity 0: 空 array
# - arity 1: copy 或 convert（视 arg 类型）
# - arity 2: (n, default_val) — n 个 default_val 的填充
# - arity 3: (n, default_val, allocator) — alloc 忽略
_VECTOR_PATTERNS = {
    0: ConstructorResolution("#[]", is_default=True),
    1: ConstructorResolution("({a0} : Array _)"),       # copy
    2: ConstructorResolution("Array.replicate (({a0}).toNat) {a1}"),  # size → Nat
    3: ConstructorResolution("Array.replicate (({a0}).toNat) {a1}"),
}

# `std::pair<A, B>(a, b)`
_PAIR_PATTERNS = {
    2: ConstructorResolution("({a0}, {a1})"),
}

# `std::map<K, V>` / `std::set<T>`
_MAP_PATTERNS = {
    0: ConstructorResolution("StdMap.empty", is_default=True),
    1: ConstructorResolution("({a0} : StdMap _ _)"),
}
_SET_PATTERNS = {
    0: ConstructorResolution("StdMap.empty", is_default=True),  # set 视为 map<T, Unit>
    1: ConstructorResolution("({a0} : StdMap _ Unit)"),
}

# RNG 三件套
_RNG_PATTERNS = {
    "mt19937": {
        0: ConstructorResolution("Rng.default", is_default=True),
        1: ConstructorResolution("Rng.mk {a0}"),
    },
    "random_device": {
        0: ConstructorResolution("Rng.default", is_default=True),
    },
    "uniform_int_distribution": {
        2: ConstructorResolution("UniformIntDist.mk {a0} {a1}"),
    },
}


def _try_stl_pattern(typename: str, arity: int) -> ConstructorResolution | None:
    """规则匹配 STL 容器构造器。typename 已经 normalize（剥 std::/clpoly:: 等）。"""
    # vector<T>
    if typename.startswith("vector<"):
        return _VECTOR_PATTERNS.get(arity)
    # pair<A, B>
    if typename.startswith("pair<"):
        return _PAIR_PATTERNS.get(arity)
    # map<K, V> / multimap<K, V>（CLPoly 不用 multimap，但保险）
    if typename.startswith("map<") or typename.startswith("multimap<"):
        return _MAP_PATTERNS.get(arity)
    # set<T>
    if typename.startswith("set<") or typename.startswith("unordered_set<"):
        return _SET_PATTERNS.get(arity)
    # mt19937 / mersenne_twister_engine
    if typename.startswith("mt19937") or "mersenne_twister" in typename:
        return _RNG_PATTERNS["mt19937"].get(arity)
    if typename == "random_device":
        return _RNG_PATTERNS["random_device"].get(arity)
    if typename.startswith("uniform_int_distribution"):
        return _RNG_PATTERNS["uniform_int_distribution"].get(arity)
    if typename.startswith("uniform_real_distribution"):
        return _RNG_PATTERNS["uniform_int_distribution"].get(arity)  # alias 复用
    return None


# ============================================================
# default_init
# ============================================================

DEFAULT_INIT_MAP: dict[str, ConstructorResolution] = {
    "default_init_const allocator_type": ConstructorResolution("()", is_default=True),
    "default_init_const std::allocator<char>": ConstructorResolution("()", is_default=True),
    "default_init_int":  ConstructorResolution("0", is_default=True),
    "default_init_bool": ConstructorResolution("false", is_default=True),
    # 其他 default_init_* 走规则
}


# ============================================================
# resolve API
# ============================================================

def resolve_constructor(op_name: str, arity: int) -> ConstructorResolution | None:
    """对 `Call(UnresolvedOp("construct_X" | "default_init_X"), args)` 解析。

    返回 ConstructorResolution（template + is_default 标志），或 None
    （未命中 → Pass 5 走 B 策略保留为 UnresolvedOp）。
    """
    # default_init_*
    if op_name.startswith("default_init"):
        # 精确匹配
        if op_name in DEFAULT_INIT_MAP:
            return DEFAULT_INIT_MAP[op_name]
        # default_init_<TypeName> → 用 TypeName 的 0-arity 构造器
        rest = op_name[len("default_init_"):].strip() if op_name.startswith("default_init_") else ""
        if rest:
            tn = _normalize_typename(rest)
            if tn in CLPOLY_CONSTRUCTORS and 0 in CLPOLY_CONSTRUCTORS[tn]:
                return CLPOLY_CONSTRUCTORS[tn][0]
            stl = _try_stl_pattern(tn, 0)
            if stl is not None:
                return stl
        return None

    # construct_*
    if not op_name.startswith("construct_"):
        return None

    typename = _normalize_typename(op_name)

    # 1. CLPoly 核心 / aux struct / Iterator
    if typename in CLPOLY_CONSTRUCTORS:
        entry = CLPOLY_CONSTRUCTORS[typename].get(arity)
        if entry is not None:
            return entry

    # 2. STL 规则
    stl = _try_stl_pattern(typename, arity)
    if stl is not None:
        return stl

    # 3. 未命中
    return None
