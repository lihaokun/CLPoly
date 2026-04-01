"""
C++ → Lean 映射表

翻译器查此表。不在表中 → sorry。
"""

from ir_types import BaseType, StructType

# ============================================================
# CLASS_MAP: 类方法映射
# ============================================================

CLASS_MAP = {
    "Zp": {
        "lean_type": StructType("Zp", []),
        "constructors": {
            (BaseType.INT64, BaseType.UINT64): "Zp.ofInt",
            (BaseType.UINT64, BaseType.UINT64): "Zp.ofUInt64",
            (BaseType.INT64,): "Zp.ofInt",  # 单参数时用默认 prime
            (): "default",
        },
        "methods": {
            "number": ("field", "val"),
            "prime": ("field", "prime"),
            "val": ("field", "val"),
            "inv": ("method", "Zp.inv"),
        },
        "operators": {
            "+": None, "-": None, "*": None,
            "/": "Zp.div",
            "==": None, "!=": None,
        },
    },

    "ZZ": {
        "lean_type": BaseType.INT64,
        "constructors": {
            (BaseType.INT64,): "id",
            (BaseType.UINT64,): "Int.ofNat",
        },
        "methods": {},
        "operators": {
            "+": None, "-": None, "*": None,
            "/": None, "%": None,
            "==": None, "!=": None, "<": None, ">": None,
        },
    },

    "UMonomial": {
        "lean_type": StructType("UMonomial", []),
        "constructors": {
            (BaseType.UINT64,): "UMonomial.mk",
            (): "default",
        },
        "methods": {
            "deg": ("field", "deg"),
        },
    },

    "SparsePolyZp": {
        "lean_type": StructType("SparsePolyZp", []),
        "constructors": {
            (): "SparsePolyZp.empty",
        },
        "methods": {
            "empty": ("method", "Array.isEmpty"),
            "front": ("method", "SparsePolyZp.front!"),
            "back": ("method", "SparsePolyZp.back!"),
            "size": ("method", "SparsePolyZp.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "SparsePolyZp.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
            "comp": ("method", "SparsePolyZp.comp"),
        },
    },
}

# ============================================================
# FUNC_MAP: 独立函数映射
# ============================================================

FUNC_MAP = {
    # 函数名 → (Lean 函数名, 参数规则)
    # 参数规则:
    #   "direct" = 直接传参
    #   "identity" = 返回第一个实参（如 std::move）
    #   "make_pair" = Prod.mk
    #   "out2_drop1(a,b,c,d,e)" = 输出前2个，输入中间2个，丢弃最后1个

    "derivative": ("SparsePolyZp.derivative", "direct"),
    "polynomial_GCD": ("SparsePolyZp.gcd", "direct"),
    "pair_vec_div": ("SparsePolyZp.divmod", "out2_drop1"),
    "get_deg": ("SparsePolyZp.getDeg", "direct"),
    "move": ("id", "identity"),
    "make_pair": ("Prod.mk", "make_pair"),
}

# ============================================================
# FIELD_MAP: 成员字段映射（C++ → Lean）
# ============================================================

FIELD_MAP = {
    "first": "fst",
    "second": "snd",
    "val": "val",
    "deg": "deg",
    "prime": "prime",
}

# ============================================================
# LEAN_BUILTINS: Lean 标准库函数（不加 _ir 后缀）
# ============================================================

LEAN_BUILTINS = {
    "Prod.mk", "Array.empty", "Array.mk", "Array.set!",
    "Array.push", "Array.size_u64", "Array.isEmpty",
    "Zp.ofInt", "Zp.ofUInt64", "Zp.div", "Zp.inv",
    "UMonomial.mk",
    "SparsePolyZp.empty", "SparsePolyZp.front!",
    "SparsePolyZp.back!", "SparsePolyZp.getDeg",
    "SparsePolyZp.size_u64", "SparsePolyZp.normalization",
    "SparsePolyZp.derivative", "SparsePolyZp.gcd",
    "SparsePolyZp.divmod", "SparsePolyZp.comp",
    "Int.ofNat", "id",
}

# ============================================================
# NOOP_METHODS: 调用后返回 self（无操作）
# ============================================================

NOOP_METHODS = {"reserve", "data"}
