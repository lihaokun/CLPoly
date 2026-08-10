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
            "+=": None, "-=": None, "*=": None, "/=": "Zp.div",
        },
    },

    "ZZ": {
        "lean_type": StructType("ZZ", []),
        "constructors": {
            (BaseType.INT64,): "id",
            (BaseType.UINT64,): "Int.ofNat",
            (): "((0 : Int))",
        },
        "methods": {
            "sizeinbase": ("method", "ZZ.sizeinbase"),
            "fdiv_ui": ("method", "ZZ.fdiv_ui"),
        },
        "operators": {
            # ZZ → bool: 非零为 true
            "bool": "ZZ.toBool",
            # ZZ comparison
            "==": "BEq.beq",
            "!=": "(fun a b => a != b)",
        },
        "operators": {
            "+": None, "-": None, "*": None,
            "/": None, "%": None,
            "==": None, "!=": None, "<": None, ">": None,
            "<=": None, ">=": None,
            "<<": None, ">>": None,  # P1-6: ZZ << k / ZZ >> k 用 typeclass
            "&": None, "|": None, "^": None,  # P1-6: 位运算
            "bool": "ZZ.toBool",  # operator bool: ZZ → Bool (nonzero check)
            # compound assignment：Pass 5 将 `x += y` 展开为 `x := x + y`；
            # 这里注册占位让 CLASS_MAP 查表命中，Pass 5 走 compound 展开分支
            "+=": None, "-=": None, "*=": None, "/=": None, "%=": None,
            "<<=": None, ">>=": None, "|=": None, "&=": None, "^=": None,
        },
        "static_methods": {
            "fdiv_q": "ZZ.fdiv_q",
            "fdiv_r": "ZZ.fdiv_r",
            "invert": "ZZ.invert",
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
            "empty": ("method", "UMonomial.isEmpty"),
            "normalization": ("noop", None),
        },
    },

    "DenseUPolyZp": {
        "lean_type": StructType("DenseUPolyZp", []),
        # Constructor bodies are part of the strict generated GCD closure;
        # these names must never be mapped to SparsePolyZp.gcd or another L2
        # implementation.
        "constructors": {},
        "methods": {
            "lead": ("method", "dense_upoly_zp_lead_ir"),
            "empty": ("method", "dense_upoly_zp_empty_ir"),
            "__strip": ("mutate", "dense_upoly_zp___strip_ir"),
            "nmod_inv": ("method", "dense_upoly_zp_nmod_inv_ir"),
            "nmod_mul": ("method", "dense_upoly_zp_nmod_mul_ir"),
            "nmod_add": ("method", "dense_upoly_zp_nmod_add_ir"),
            "nmod_sub": ("method", "dense_upoly_zp_nmod_sub_ir"),
            "__precompute": ("mutate", "dense_upoly_zp___precompute_ir"),
            "scalar_mul": ("mutate", "dense_upoly_zp_scalar_mul_ir"),
            "deg": ("method", "dense_upoly_zp_deg_ir"),
            "to_upoly": ("method", "dense_upoly_zp_to_upoly_ir"),
            "_gcd_hgcd": ("method", "dense_upoly_zp__gcd_hgcd_ir"),
            "_poly_divrem": ("method", "dense_upoly_zp__poly_divrem_ir"),
            "_hgcd_recursive": ("method", "dense_upoly_zp__hgcd_recursive_ir"),
            "_gcd_euclid": ("method", "dense_upoly_zp__gcd_euclid_ir"),
            "_hgcd_iter": ("method", "dense_upoly_zp__hgcd_iter_ir"),
            "_mul": ("method", "dense_upoly_zp__mul_ir"),
            "_poly_sub": ("method", "dense_upoly_zp__poly_sub_ir"),
            "_poly_add": ("method", "dense_upoly_zp__poly_add_ir"),
            "_mat_mul": ("method", "dense_upoly_zp__mat_mul_ir"),
            "_mat_row_update": ("method", "dense_upoly_zp__mat_row_update_ir"),
            "_classical_mul": ("method", "dense_upoly_zp__classical_mul_ir"),
            "_kar_mul": ("method", "dense_upoly_zp__kar_mul_ir"),
            "_mat_mul_entry": ("method", "dense_upoly_zp__mat_mul_entry_ir"),
        },
    },

    "SparsePolyZp": {
        "lean_type": StructType("SparsePolyZp", []),
        "constructors": {
            (): "SparsePolyZp.empty",
        },
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None,
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
            "begin": ("method", "SparsePolyZp.toList"),
            "end": ("method", "SparsePolyZp.toList"),
            "clear": ("mutate", "Array.clearVec"),
            "assign": ("mutate", "Array.replicateMut"),  # vec.assign(n, val) → Array.replicateMut recv n val（receiver 忽略，复用 mutate 模板）
            "resize": ("noop", None),
            "erase": ("mutate", "Array.pop"),
            "at": ("method", "Array.get!"),
        },
    },
    "SparsePolyZZ": {
        "lean_type": StructType("SparsePolyZZ", []),
        "constructors": {
            (): "SparsePolyZZ.empty",
        },
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None,
        },
        "methods": {
            "empty": ("method", "Array.isEmpty"),
            "front": ("method", "SparsePolyZZ.front!"),
            "back": ("method", "SparsePolyZZ.back!"),
            "size": ("method", "SparsePolyZZ.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "SparsePolyZZ.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
            "comp": ("method", "SparsePolyZZ.comp"),
            "begin": ("method", "SparsePolyZZ.toList"),
            "end": ("method", "SparsePolyZZ.toList"),
            "clear": ("mutate", "Array.clearVec"),
            "assign": ("mutate", "Array.replicateMut"),  # vec.assign(n, val) → Array.replicateMut recv n val（receiver 忽略，复用 mutate 模板）
            "erase": ("mutate", "Array.eraseAny"),  # P1-7: range/iter erase（兜底，接受任意 idx 类型）
            "at":    ("method", "Array.get!"),
        },
    },

    "QQ": {
        "lean_type": StructType("QQ", []),
        "constructors": {
            (BaseType.INT64,): "QQ.ofInt",
            (BaseType.INT64, BaseType.INT64): "QQ.mk",
            (): "((0 : Rat))",
        },
        "methods": {
            "get_num": ("method", "QQ.num"),
            "get_den": ("method", "QQ.den"),
        },
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None, "<": None, ">": None,
            "<=": None, ">=": None,
            "+=": None, "-=": None, "*=": None, "/=": None,
        },
    },

    "MvPolyZZ": {
        "lean_type": StructType("MvPolyZZ", []),
        "constructors": {
            (): "MvPolyZZ.empty",
        },
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None,
        },
        "methods": {
            "empty": ("method", "MvPolyZZ.isEmpty"),
            "front": ("method", "MvPolyZZ.front!"),
            "back": ("method", "MvPolyZZ.back!"),
            "size": ("method", "MvPolyZZ.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "MvPolyZZ.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
            "comp": ("method", "MvPolyZZ.comp"),
            "comp_ptr": ("method", "MvPolyZZ.comp"),
            "begin": ("method", "MvPolyZZ.toList"),
            "end": ("method", "MvPolyZZ.toList"),
            "clear": ("mutate", "Array.clearVec"),
        },
    },

    "MvPolyZp": {
        "lean_type": StructType("MvPolyZp", []),
        "constructors": {
            (): "MvPolyZp.empty",
        },
        "operators": {
            "+": None, "-": None, "*": None,
            "==": None, "!=": None,
            "+=": None, "-=": None, "*=": None,
        },
        "methods": {
            "empty": ("method", "MvPolyZp.isEmpty"),
            "front": ("method", "MvPolyZp.front!"),
            "back": ("method", "MvPolyZp.back!"),
            "size": ("method", "MvPolyZp.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "MvPolyZp.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
            "comp": ("method", "MvPolyZp.comp"),
            "comp_ptr": ("method", "MvPolyZp.comp"),
        },
    },

    "MvMonomial": {
        "lean_type": StructType("MvMonomial", []),
        "constructors": {
            (): "MvMonomial.empty",
        },
        "methods": {
            "deg": ("field", "deg"),
            "empty": ("method", "MvMonomial.isEmpty"),
            "push_back": ("mutate_push", "Array.push"),
            "begin": ("method", "MvMonomial.toList"),
            "end": ("method", "MvMonomial.toList"),
        },
    },

    "Variable": {
        "lean_type": StructType("Variable", []),
        "constructors": {},
        "operators": {
            "==": None, "!=": None, "<": None, ">": None,
            "<=": None, ">=": None,
        },
        "methods": {
            "insert": ("mutate", "Variable.insert"),
            "find": ("method", "Variable.find"),
            "end": ("method", "Variable.end"),
            "begin": ("method", "Variable.begin"),
            "size": ("method", "Variable.size"),
            "empty": ("method", "Variable.isEmpty"),
            "at": ("method", "Variable.get!"),
            "erase": ("mutate", "Variable.erase"),
        },
    },

    "LLLMatrix": {
        "lean_type": StructType("LLLMatrix", []),
        "constructors": {
            (): "LLLMatrix.empty",
        },
        "methods": {
            "size": ("method", "LLLMatrix.size"),
            "assign": ("mutate", "Array.replicateMut"),  # vec.assign(n, val) → Array.replicateMut recv n val（receiver 忽略，复用 mutate 模板）
            "push_back": ("mutate_push", "Array.push"),
            "reserve": ("noop", None),
            "swap": ("mutate", "LLLMatrix.swap"),
        },
    },

    "HenselNode": {
        "lean_type": StructType("HenselNode", []),
        "constructors": {
            (): "HenselNode.default",
        },
        "methods": {},
    },

    "PrimeSelectionResult": {
        "lean_type": StructType("PrimeSelectionResult", []),
        "constructors": {},
        "methods": {},
    },

    "WangLcResult": {
        "lean_type": StructType("WangLcResult", []),
        "constructors": {},
        "methods": {},
    },

    # Iterator — Pass 4 已识别绝大多数 iter 循环；本表用于 impure compact-erase
    # 等 Pass 4 拒识场景的兜底（保留为通用 Lean 操作）
    "Iterator": {
        "lean_type": StructType("Iterator", []),
        "constructors": {},
        "operators": {
            "++": ("mutate", "Iterator.advance"),
            "--": ("mutate", "Iterator.retreat"),
            "*":  ("method", "Iterator.deref!"),
            "->": ("method", "Iterator.deref!"),
            "==": None, "!=": None,
        },
        "methods": {},
    },

    "StdMap": {
        "lean_type": StructType("StdMap", []),
        "constructors": {
            (): "StdMap.empty",
        },
        "methods": {
            "find": ("method", "StdMap.find"),
            "end": ("method", "StdMap.end"),
            "begin": ("method", "StdMap.begin"),
            "at": ("method", "StdMap.get!"),
            "erase": ("mutate", "StdMap.erase"),
            "insert": ("mutate", "StdMap.insert"),
            "size": ("method", "StdMap.size"),
            "empty": ("method", "StdMap.isEmpty"),
        },
    },

    "Factorization": {
        "lean_type": StructType("Factorization", []),
        "constructors": {
            (): "Factorization.empty",
        },
        "methods": {
            "content": ("field", "content"),
            "factors": ("field", "factors"),
        },
    },

    # ====================================================================
    # Array — 通用 std::vector<T> 容器（Pass 1 把 std::vector 解析为
    # ArrayType；Pass 5 用 "Array" 作为 receiver 键查方法）
    # ====================================================================
    "Array": {
        "lean_type": StructType("Array", []),
        "constructors": {
            (): "#[]",
        },
        "methods": {
            "size":         ("method",       "Array.size"),
            "empty":        ("method",       "Array.isEmpty"),
            "front":        ("method",       "Array.head!"),       # UB-3
            "back":         ("method",       "Array.getLast!"),    # UB-3
            "data":         ("identity",     None),                # vector::data() ≈ self
            "begin":        ("method",       "Array.toList"),      # 配合 Pass 4 已识别 iter
            "end":          ("method",       "Array.toList"),
            "cbegin":       ("method",       "Array.toList"),
            "cend":         ("method",       "Array.toList"),
            "push_back":    ("mutate_push",  "Array.push"),
            "pop_back":     ("mutate",       "Array.pop"),
            "emplace_back": ("mutate_push",  "Array.push"),
            "resize":       ("mutate",       "Array.resize"),
            "reserve":      ("noop",         None),
            "clear":        ("mutate",       "Array.clearVec"),
            "assign":       ("mutate",       "Array.replicateMut"),
            "erase":        ("mutate",       "Array.eraseAny"),       # Pass 4 漏识别的兜底
            "at":           ("method",       "Array.get!"),        # UB-2
            "insert":       ("mutate",       "Array.insert"),
            "swap":         ("mutate",       "Array.swap"),
            "find":         ("method",       "Array.findVal"),  # 按值找（不是 predicate）
            "contains":     ("method",       "Array.contains"),
        },
        "operators": {
            "[]":  "Array.get!",   # UB-2
            "()":  None,           # 不应出现
            "==":  None, "!=": None,
            "=":   None,           # operator= → AssignStmt（特殊处理）
        },
    },

    # ====================================================================
    # Monomial — basic_monomial<...> typedef 别名；alias 到 MvMonomial
    # ====================================================================
    "Monomial": {
        "lean_type": StructType("MvMonomial", []),
        "constructors": {
            (): "MvMonomial.empty",
        },
        "methods": {
            "deg":          ("field",        "deg"),
            "empty":        ("method",       "MvMonomial.isEmpty"),
            "size":         ("method",       "MvMonomial.size_u64"),
            "front":        ("method",       "MvMonomial.front!"),
            "back":         ("method",       "MvMonomial.back!"),
            "data":         ("identity",     None),
            "push_back":    ("mutate_push",  "Array.push"),
            "pop_back":     ("mutate",       "Array.pop"),
            "begin":        ("method",       "MvMonomial.toList"),
            "end":          ("method",       "MvMonomial.toList"),
            "normalization": ("mutate",      "MvMonomial.normalization"),
            "clear":        ("mutate",       "Array.clearVec"),
            "reserve":      ("noop",         None),
        },
    },

    # ====================================================================
    # Poly / PolyZp / PolyQQ / PolyZZ — CLPoly 别名 typedef 到 MvPoly*
    # （Pass 1 line 217 _TYPEDEF_ALIASES 里映射到独立 NamedType；这里
    # 注册等同方法集，避免 Pass 5 漏命中）
    # ====================================================================
    "Poly": {
        "lean_type": StructType("MvPolyZZ", []),
        "constructors": {(): "MvPolyZZ.empty"},
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None,
        },
        "methods": {
            # Poly == polynomial_<ZZ, lex_<var_order>> == MvPolyZZ alias
            "size":          ("method",       "MvPolyZZ.size_u64"),
            "empty":         ("method",       "MvPolyZZ.isEmpty"),
            "front":         ("method",       "MvPolyZZ.front!"),
            "back":          ("method",       "MvPolyZZ.back!"),
            "data":          ("identity",     None),
            "comp":          ("method",       "MvPolyZZ.comp"),
            "comp_ptr":      ("method",       "MvPolyZZ.comp"),
            "begin":         ("method",       "MvPolyZZ.toList"),
            "end":           ("method",       "MvPolyZZ.toList"),
            "push_back":     ("mutate_push",  "Array.push"),
            "pop_back":      ("mutate",       "Array.pop"),
            "normalization": ("mutate",       "MvPolyZZ.normalization"),
            "reserve":       ("noop",         None),
            "clear":         ("mutate",       "Array.clearVec"),
            "assign":        ("mutate",       "Array.replicateMut"),
        },
    },

    "PolyZp": {
        "lean_type": StructType("MvPolyZp", []),
        "constructors": {(): "MvPolyZp.empty"},
        "operators": {
            "+": None, "-": None, "*": None, "/": None,
            "==": None, "!=": None,
        },
        "methods": {
            "size":          ("method",       "MvPolyZp.size_u64"),
            "empty":         ("method",       "MvPolyZp.isEmpty"),
            "front":         ("method",       "MvPolyZp.front!"),
            "back":          ("method",       "MvPolyZp.back!"),
            "data":          ("identity",     None),
            "comp":          ("method",       "MvPolyZp.comp"),
            "comp_ptr":      ("method",       "MvPolyZp.comp"),
            "begin":         ("method",       "MvPolyZp.toList"),
            "end":           ("method",       "MvPolyZp.toList"),
            "push_back":     ("mutate_push",  "Array.push"),
            "pop_back":      ("mutate",       "Array.pop"),
            "normalization": ("mutate",       "MvPolyZp.normalization"),
            "reserve":       ("noop",         None),
            "clear":         ("mutate",       "Array.clearVec"),
            "assign":        ("mutate",       "Array.replicateMut"),
        },
    },

    "PolyZZ": {
        "lean_type": StructType("MvPolyZZ", []),
        "constructors": {(): "MvPolyZZ.empty"},
        "methods": {
            # 与 Poly 等同
            "size":          ("method",       "MvPolyZZ.size_u64"),
            "empty":         ("method",       "MvPolyZZ.isEmpty"),
            "front":         ("method",       "MvPolyZZ.front!"),
            "back":          ("method",       "MvPolyZZ.back!"),
            "data":          ("identity",     None),
            "comp":          ("method",       "MvPolyZZ.comp"),
            "comp_ptr":      ("method",       "MvPolyZZ.comp"),
            "begin":         ("method",       "MvPolyZZ.toList"),
            "end":           ("method",       "MvPolyZZ.toList"),
            "push_back":     ("mutate_push",  "Array.push"),
            "pop_back":      ("mutate",       "Array.pop"),
            "normalization": ("mutate",       "MvPolyZZ.normalization"),
            "reserve":       ("noop",         None),
            "clear":         ("mutate",       "Array.clearVec"),
        },
    },

    "PolyQQ": {
        "lean_type": StructType("PolyQQ", []),
        "constructors": {(): "PolyQQ.empty"},
        "methods": {
            "size":          ("method",       "PolyQQ.size_u64"),
            "empty":         ("method",       "PolyQQ.isEmpty"),
            "front":         ("method",       "PolyQQ.front!"),
            "back":          ("method",       "PolyQQ.back!"),
            "data":          ("identity",     None),
            "comp":          ("method",       "PolyQQ.comp"),
            "begin":         ("method",       "PolyQQ.toList"),
            "end":           ("method",       "PolyQQ.toList"),
            "push_back":     ("mutate_push",  "Array.push"),
            "normalization": ("mutate",       "PolyQQ.normalization"),
            "reserve":       ("noop",         None),
            "clear":         ("mutate",       "Array.clearVec"),
        },
    },
}

# ============================================================
# FUNC_MAP: 独立函数映射
# ============================================================

FUNC_MAP = {
    # 函数名 → (Lean 函数名, 参数规则, 输出参数索引)
    # 参数规则:
    #   "direct" = 直接传参
    #   "identity" = 返回第一个实参（如 std::move）
    #   "make_pair" = Prod.mk
    #   "output" = 输出参数模式：output_indices 指定的参数从 args 移到返回值
    # 输出参数索引：从 0 开始的参数位置列表（C++ 调用中的位置）

    "derivative": ("derivative", "direct"),
    "inv_prime": ("inv_prime_ir", "direct"),
    "__preinvert_limb": ("dense_upoly_zp___preinvert_limb_ir", "direct"),
    "_umul128": ("dense_upoly_zp__umul128_ir", "direct"),
    "_add_carry3": ("dense_upoly_zp__add_carry3_ir", "direct"),
    "_lll_mod_preinv": ("dense_upoly_zp__lll_mod_preinv_ir", "direct"),
    "divrem": ("dense_upoly_zp_divrem_nonalias_ir", "direct"),
    "__builtin_clzll": ("uint64_clz", "direct"),
    "polynomial_GCD": ("polynomial_GCD", "direct"),
    "polynomial_GCD_eea": ("polynomial_GCD_eea", "direct"),
    "pair_vec_div5": ("pair_vec_div5", "direct"),
    "pair_vec_div": ("pair_vec_div", "direct"),
    "pair_vec_multiplies": ("pair_vec_multiplies", "direct"),
    "get_deg": ("get_deg", "direct"),
    "move": ("id", "identity"),
    "make_pair": ("Prod.mk", "make_pair"),
    "pow": ("HPow.hPow", "direct"),
    # 素数
    "next_prime_64": ("next_prime_64", "direct"),
    "prev_prime_64": ("prev_prime_64", "direct"),
    # 多项式转换
    "poly_convert": ("poly_convert", "direct"),
    # squarefreefactorize 移到 TRANSLATION_SCOPE（参 L820 附近），由 Pass 8
    # 模板实例后缀派发 emit `squarefreefactorize_lex_ir`；不再在 FUNC_MAP 里。
    "degree": ("degree", "direct"),
    "is_number": ("is_number", "direct"),
    # 排序（翻译为 identity，排序不影响正确性）
    "sort": ("id", "identity"),
    # swap
    "swap": ("swap", "direct"),
    # abs
    "abs": ("Int.natAbs", "direct"),
    # ZZ 静态方法
    "fdiv_q": ("ZZ.fdiv_q", "direct"),
    "fdiv_r": ("ZZ.fdiv_r", "direct"),
    "invert": ("ZZ.invert", "direct"),
    # std 函数
    "max": ("max", "direct"),
    "min": ("min", "direct"),
    "iota": ("List.range", "direct"),
    # 多项式操作
    "get_variables": ("get_variables", "direct"),
    "assign": ("assign", "direct"),
    # 多项式查询
    "leadcoeff": ("leadcoeff", "direct"),
    "get_variables": ("get_variables", "direct"),
    "content": ("content", "direct"),
    "upoly_prem": ("upoly_prem", "direct"),
    "polynomial_mod": ("polynomial_mod", "direct"),
    "cont": ("cont", "direct"),
    # C++ std::log(double) → Lean Float.log（自然对数；Float 内置）
    "log": ("Float.log", "direct"),
    # C++ ceil(double) / floor(double) → 返回 Float（Lean 同名）
    "ceil": ("Float.ceil", "direct"),
    "floor": ("Float.floor", "direct"),
    "get_first_deg": ("get_first_deg", "direct"),
    "gcd": ("gcd", "direct"),
    "pp": ("pp", "direct"),
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
    "Prod.mk", "Array.clearVec", "Array.mk", "Array.set!",
    "Array.push", "Array.size_u64", "Array.isEmpty",
    "Zp.ofInt", "Zp.ofUInt64", "Zp.div", "Zp.inv",
    "UMonomial.mk",
    "SparsePolyZp.empty", "SparsePolyZp.front!",
    "SparsePolyZp.back!", "SparsePolyZp.getDeg",
    "SparsePolyZp.size_u64", "SparsePolyZp.normalization",
    "SparsePolyZp.derivative", "SparsePolyZp.gcd",
    "SparsePolyZp.divmod", "SparsePolyZp.comp",
    "Int.ofNat", "id",
    "Rng.next", "Rng.step",
    "#[]", "default",
    "HPow.hPow",
}

# ============================================================
# NOOP_METHODS: 调用后返回 self（无操作）
# ============================================================

# ============================================================
# NOOP_METHODS / MUTATING_METHODS：从 CLASS_MAP 自动派生
# ============================================================

def _derive_method_sets():
    """从 CLASS_MAP 派生 noop/identity/mutate 方法集合。"""
    noop = set()
    mutating = set()
    for cls_info in CLASS_MAP.values():
        for method_name, (category, _lean_name) in cls_info.get("methods", {}).items():
            if category in ("noop", "identity"):
                noop.add(method_name)
            elif category in ("mutate", "mutate_push"):
                mutating.add(method_name)
    return noop, mutating

NOOP_METHODS, MUTATING_METHODS = _derive_method_sets()

# ============================================================
# OPERATOR_MAP: C++ operator overload → Lean 运算符
# key = C++ operator 名（Clang referencedDecl.name 中的子串）
# value = (最少参数数, Lean 运算符字符串)
# ============================================================

OPERATOR_MAP = {
    "operator/":  (3, "/"),
    "operator*":  (3, "*"),
    "operator+":  (3, "+"),
    "operator-":  (3, "-"),
    "operator==": (3, "=="),
    "operator!=": (3, "!="),
    "operator<":  (3, "<"),
    "operator>":  (3, ">"),
    "operator<=": (3, "<="),
    "operator>=": (3, ">="),
    "operator%":  (3, "%"),
    "operator<<": (3, "<<"),
    "operator>>": (3, ">>"),
}

# ============================================================
# STRUCT_COERCE_MAP: StructType → BaseType 隐式转换
# key = (struct_name, BaseType)
# value = 字段访问表达式模板，{e} 是内部表达式
# ============================================================

STRUCT_COERCE_MAP = {
    ("Zp", BaseType.UINT64): "{e}.val",
    ("Zp", BaseType.INT64):  "({e}.val.toNat : Int)",
    ("UMonomial", BaseType.UINT64): "{e}.deg",
    ("ZZ", BaseType.INT64): "{e}",           # ZZ = Int，同类型
    ("ZZ", BaseType.UINT64): "{e}.toNat.toUInt64",
    ("MvMonomial", BaseType.UINT64): "{e}.deg",
}

# ============================================================
# UNSAFE_CAST_PAIRS: unsigned→signed 转换需要证明值在范围内
# (source_type, target_type) → lean_prop 模板，{e} 是内部表达式
# ============================================================

UNSAFE_CAST_PAIRS = {
    (BaseType.UINT64, BaseType.INT64):  "{e}.val ≤ Int64.max",
    (BaseType.UINT128, BaseType.INT64): "{e} ≤ Int64.max",
}

# ============================================================
# EMPTY_CONTAINER_METHODS: 空容器上调用是 UB 的方法
# 从 CLASS_MAP 中 method 名含 "!" 的自动派生
# ============================================================

def _derive_empty_container_methods():
    """从 CLASS_MAP 派生需要非空检查的方法。"""
    result = set()
    for cls_info in CLASS_MAP.values():
        for _method_name, (category, lean_name) in cls_info.get("methods", {}).items():
            if lean_name and "!" in lean_name:
                result.add(lean_name)
    return result

EMPTY_CONTAINER_METHODS = _derive_empty_container_methods()

# ============================================================
# ASSERT_FAIL_NAMES: assert 宏展开后的函数名
# ============================================================

ASSERT_FAIL_NAMES = {"__assert_fail", "__assert_rtn", "__assert"}

# ============================================================
# CALL_OPERATOR_MAP: C++ 函数调用运算符 operator() → Lean 函数
# 用于 std::uniform_int_distribution(rng) 等模式
# value = (lean_func, arg_order)
#   arg_order: "reverse" = operator()(this, arg) → lean_func arg this
#              "direct"  = operator()(this, arg) → lean_func this arg
# ============================================================

CALL_OPERATOR_MAP = {
    "operator()": ("Rng.next", "reverse"),  # dist(rng) → Rng.next rng dist
}

# ============================================================
# 翻译范围：本次翻译的 C++ 函数
# 范围内函数互相调用加 _ir 后缀
# ============================================================

# M2：TRANSLATION_SCOPE 中函数的输出参数（从 C++ 源码扫描）
# func_name → [输出参数索引]（非 const 引用参数位置）
TRANSLATION_SCOPE_OUTPUT_PARAMS = {
    "_umul128": [0, 1],                 # hi, lo
    "_add_carry3": [0],                 # word3 accumulator
    "divrem": [0, 1],                   # dense quotient and remainder
    "__upoly_make_monic": [0],           # f
    "__upoly_mod_coeff": [0],            # f
    "__upoly_divmod_mod": [0, 1],        # q, r
    "__edf_Zp": [0, 3],                  # result, rng
    "__upoly_random": [2],               # (max_deg, p, rng&) — rng output param
    "__hensel_tree_build_recursive": [0], # nodes
    "__hensel_step": [0],                # node
    "__hensel_extract_factors": [2],      # factors
    "__hensel_lift_recursive": [0],       # nodes
    "__hensel_step_linear": [0],         # node
    "__hensel_lift_linear_recursive": [0], # nodes
    "__build_cld_matrix": [0],           # M
    "__lll_reduce": [0, 1],              # M, U
    "__si_vandermonde_solve": [2],       # coeffs
    "__si_theta_array_eval": [5],        # images（C++ 第 6 参数 = 索引 5；之前误填 6 越界）
    # MTSHL 函数原注册全部错索引（误填"第 N 个/第 M 个"形如 [3, 9]）。
    # 重新对照 polynomial_factorize_wang.hh 真实参数列表（agent 第七轮发现）：
    "__mtshl_zp_univar_mdp": [2],        # (F, c, sigma&) — 3 args，sigma at idx 2
    "__mtshl_multi_bdp": [5],            # (F, c, x1, x2, alpha2, result&) — 6 args
    "__mtshl_sparse_int": [6],           # (F, c, forms, x1, aux_vars, p, result&) — 7 args
    "__mtshl_wmds": [5],                 # (F, c, x1, aux_vars, ideal_alphas_zp, result&) — 6 args
    "__mtshl_step_j": [1],               # (aj, F&, lc_tau, xj, alpha_j, x1, aux_vars, ideal_alphas_zp, p) — 9 args，F at idx 1
    # P0-2 同源补全（修正方案 §4 同类问题梳理 表）：
    # 这些 callee 不在 TRANSLATION_SCOPE 内（外部 / lib 函数），但作为
    # mutating call 调用点也需要 SSA bump。Pass 2b 用此表改写调用点。
    "poly_convert": [1],                 # p_out（upolynomial.hh）
    # pair_vec_div 有 2 个 overload — 按 arity 区分：
    "pair_vec_div#4": [0],               # (new_v&, v1, v2, comp)（basic.hh:568）
    "pair_vec_div#5": [0, 1],            # (new_v&, R&, v1, v2, comp)（basic.hh:698）
    "fdiv_r": [0],                       # r（ZZ.hh:810）
    "fdiv_q": [0],                       # q（ZZ.hh:731）
    # ZZ.invert: 模逆元，bool 返回成功与否 + 写到 out 的 inv 值
    # Pass 2b 把 out 当 ref-out → 返回 (Bool, ZZ) tuple，Pass 2b destructure
    "invert": [0],                       # out (inv 值)
    "ZZ.invert": [0],                    # 同上（Pass 5 已重命名后再次匹配）
    "swap#2": [0, 1],                    # std::swap 双向（mu/M 行交换）。
                                         # `#2` 区分 interval.hh 的 1-arg 成员 swap
    # polynomial_GCD 4-arg Bezout 形式：return GCD + 修改 s, t (Bezout coefficients)
    "polynomial_GCD#4": [2, 3],          # F, G, s&, t&（s,t 为 out）
    # __extract_monomial_content(f, var_factors&) — var_factors 输出
    # Pass 8 emit lean_name = __extract_monomial_content_lex_ir，但 Pass 2b 看 callee
    # 是 Pass 1 的 base_name `__extract_monomial_content`
    "__extract_monomial_content": [1],   # var_factors&
    # 阶段 G+：dist(rng) ref-out 修
    # Pass 5 把 `dist(rng)` translate 为 Rng.next_advance(rng, dist)
    # Pass 2b（在 Pass 5 之后重跑）按 idx=0 把 rng 当 ref-out destructure
    "Rng.next_advance": [0],             # rng&
}


def get_output_params(callee: str, num_args: int) -> list[int]:
    """查询函数的输出参数索引列表。
    优先查 `name#arity` 形式（多 overload 区分），fall back 到 `name`。
    返回空列表表示该函数不是 ref-out（无需调用点改写）。
    """
    key_arity = f"{callee}#{num_args}"
    if key_arity in TRANSLATION_SCOPE_OUTPUT_PARAMS:
        return TRANSLATION_SCOPE_OUTPUT_PARAMS[key_arity]
    return TRANSLATION_SCOPE_OUTPUT_PARAMS.get(callee, [])


# 第八轮 P0 修复：哪些 ref-out 函数是非 void（含 orig_ret 作为返回 tuple 第一个元素）。
# Pass 2b 据此决定 destructure 字段命名约定：
#   void + 1 out  → 直接 AssignStmt(out := f(...))
#   void + 2 out  → fst, snd (refs)
#   void + N≥3 out → elem0..N-1 (refs)
#   non-void + 1 out → fst (orig_ret), snd (ref) ← 必须解构
#   non-void + N≥2 out → elem0 (orig_ret), elem1..N (refs)
TRANSLATION_SCOPE_NONVOID_REFOUT: set[str] = {
    "__upoly_make_monic",       # Zp
    "__build_cld_matrix",       # int
    "__lll_reduce",             # std::vector<int>
    "__si_vandermonde_solve",   # bool
    "__mtshl_zp_univar_mdp",    # bool
    "__mtshl_multi_bdp",        # bool
    "__mtshl_sparse_int",       # bool
    "__mtshl_wmds",             # bool
    "__mtshl_step_j",           # bool
    "polynomial_GCD#4",         # polynomial_<ZZ,...>
    "invert",                   # bool（ZZ::invert 返回 success bit）
    "ZZ.invert",                # 同上（Pass 5 重命名后形式）
}


def is_callee_nonvoid(callee: str, num_args: int) -> bool:
    """callee 是否非 void（影响 destructure 字段命名）。"""
    return (f"{callee}#{num_args}" in TRANSLATION_SCOPE_NONVOID_REFOUT
            or callee in TRANSLATION_SCOPE_NONVOID_REFOUT)

# ============================================================
# REFINEMENT_MAP: L1 translated function → L2 algorithm model
#
# Each entry maps a TRANSLATION_SCOPE function (by base_name) to its
# L2 mathematical model. The generator (Pass 9) uses this to emit
# refinement theorem skeletons.
#
# Fields:
#   l2_name        — L2 function name (in CLPoly.Algorithm.* or Model)
#   l2_import      — Lean import path for the L2 module
#   refinement_file — Output file name (without .lean)
#   bridge         — Type bridge: "toPoly" | "toZMod" | "identity" | "toMathlib"
#   result_kind    — Theorem shape (see pass9 docstring)
#   cpp_source     — Original C++ source file (for docstring provenance)
#   doc            — Human-readable description
# ============================================================

REFINEMENT_MAP = {
    "__squarefree_Zp": {
        "l1_name": "__squarefree_Zp_ir",
        "l2_name": "sqfZp",
        "l2_call": "sqfZp (SparsePolyZp.toPoly p f)",
        "l2_import": "CLPoly.Algorithm.SquarefreeZp",
        "refinement_file": "SquarefreeZp",
        "result_kind": "map_eq",
        "cpp_source": "clpoly/polynomial_factorize_zp.hh",
        "doc": "无平方因子分解（Yun 算法）",
        # Pass 9 emits this public, proved contract after the strict semantic
        # implementation has been discharged in Refinement/SquarefreeZp.lean.
        # Keep the original C++ spelling (including its leading `__`) in the
        # theorem name so generated contracts remain mechanically traceable.
        "verified_contract": {
            "theorem_name": "__squarefree_Zp_raw_ir_refines_sqfZp",
            "proof_import": "CLPoly.Refinement.StrictSquarefreeGenerated",
            "proof_theorem": (
                "Refinement.StrictSquarefreeGenerated."
                "__squarefree_Zp_raw_ir_refines_sqfZp"
            ),
            "kind": "strict_squarefree_zp",
        },
    },
    "__ddf_Zp": {
        "l1_name": "__ddf_Zp_ir",
        "l2_name": "ddf",
        "l2_call": "ddf (SparsePolyZp.toPoly p f)",
        "l2_import": "CLPoly.Algorithm.DDF",
        "refinement_file": "DDF",
        "result_kind": "map_eq",
        "cpp_source": "clpoly/polynomial_factorize_zp.hh",
        "doc": "不同度数因子分解",
        "verified_contract": {
            "theorem_name": "__ddf_Zp_raw_ir_refines_ddf",
            "proof_import": "CLPoly.Refinement.DDF",
            "proof_theorem": "Refinement.StrictDDF.strictDDFEntryIR_refines_ddf",
            "kind": "strict_ddf_zp",
        },
    },
    "__make_zp": {
        "l1_name": "__make_zp_ir",
        "l2_name": "Zp.ofInt",
        "l2_call": "Zp.ofInt val.toInt modulus",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZpArith",
        "result_kind": "direct_eq",
        "cpp_source": "clpoly/number/ZZ.hh",
        "doc": "Zp 构造函数",
    },
    "__upoly_make_monic": {
        "l1_name": "__upoly_make_monic_ir",
        "l2_name": "SparsePolyZp.makeMonic",
        "l2_call": "SparsePolyZp.makeMonic f",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZpArith",
        "result_kind": "map_eq_pair",
        "cpp_source": "clpoly/upolynomial.hh",
        "doc": "首一化：multiply by inv(lc)",
    },
    "__upoly_divmod": {
        "l1_name": "__upoly_divmod_ir",
        "l2_name": "SparsePolyZp.divmod",
        "l2_call": "SparsePolyZp.divmod f g",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZpArith",
        "result_kind": "pair_eq",
        "cpp_source": "clpoly/upolynomial.hh",
        "doc": "多项式长除法：q = f / g, r = f mod g",
    },
    "__symmetric_mod": {
        "l1_name": "__symmetric_mod_ir",
        "l2_name": "ZZ.symmetricMod",
        "l2_call": "ZZ.symmetricMod a m",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZZArith",
        "result_kind": "direct_eq",
        "cpp_source": "clpoly/number/ZZ.hh",
        "doc": "对称模运算：a mod m → [-m/2, m/2]",
    },
    "__binomial": {
        "l1_name": "__binomial_ir",
        "l2_name": "ZZ.binomial",
        "l2_call": "ZZ.binomial n k",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZZArith",
        "result_kind": "direct_eq",
        "cpp_source": "clpoly/polynomial_factorize_univar.hh",
        "doc": "二项式系数 n choose k",
    },
    "__isqrt_ceil": {
        "l1_name": "__isqrt_ceil_ir",
        "l2_name": "ZZ.isqrtCeil",
        "l2_call": "ZZ.isqrtCeil n",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZZArith",
        "result_kind": "direct_eq",
        "cpp_source": "clpoly/number/ZZ.hh",
        "doc": "向上取整平方根 ceil(sqrt(n))",
    },
    "__upoly_mod": {
        "l1_name": "__upoly_mod_ir",
        "l2_name": "SparsePolyZp.divmod",
        "l2_call": "(SparsePolyZp.divmod f g).snd",
        "l2_import": "CLPoly.Model",
        "refinement_file": "ZpArith",
        "result_kind": "direct_eq",
        "cpp_source": "clpoly/upolynomial.hh",
        "doc": "多项式取模 remainder = f mod g",
    },
}

TRANSLATION_SCOPE = {
    # Zp 模块 (13)
    "__make_zp", "__upoly_make_monic", "__upoly_mod", "__upoly_divmod",
    "__upoly_powmod", "__upoly_random", "__extract_pth_root",
    "__squarefree_Zp", "__upoly_subtract_x", "__upoly_subtract_one",
    "__ddf_Zp", "__edf_Zp", "__factor_Zp",
    # Univar 模块 (33)
    "__symmetric_mod", "__upoly_symmetric_mod", "__upoly_norm_l2_sq",
    "__binomial", "__isqrt_ceil", "__mignotte_bound",
    "__upoly_mod_coeff", "__upoly_divmod_mod", "__upoly_mul_mod",
    "__hensel_tree_build_recursive", "__hensel_tree_build",
    "__hensel_step", "__hensel_extract_factors",
    "__hensel_lift_recursive", "__hensel_lift",
    "__heuristic_starting_precision", "__hensel_step_linear",
    "__hensel_lift_linear_recursive", "__upoly_norm_l1",
    "__upoly_primitive", "__subset_product_mod", "__upoly_const_term",
    "__zassenhaus_recombine", "__cld_polys", "__build_cld_matrix",
    "__lll_reduce", "__extract_candidates", "__vanhoeij_recombine",
    "__factor_recombine", "__select_prime", "__lll_factorize",
    "__factor_squarefree_primitive_ZZ", "__upoly_to_poly",
    "factorize",
    # Wang 模块 (18) — 不含经典 Wang Hensel 死代码
    # 已移除：__multivar_hensel_lift, __hensel_lift_one_var,
    #         __hensel_lc_correct, __multivar_diophantine, __pseudo_remainder_x1
    #   （经典 Wang Hensel 路线，已被 MTSHL 替代，无调用者）
    # 已移除：__taylor_coeff（ZZ 版，仅被已排除的死代码
    #         __multivar_diophantine / __hensel_lift_one_var 调用，无实例化）
    "__si_vandermonde_solve", "__mtshl_zp_univar_mdp",
    "__assign_partial_zp", "__extract_monomial_content",
    "__factor_multivar",
    "__mtshl_coeff_bound", "__mtshl_lift", "__mtshl_multi_bdp",
    "__mtshl_sparse_int", "__mtshl_step_j", "__mtshl_wmds",
    "__polynomial_to_zp",
    "__select_eval_point", "__si_theta_array_eval",
    "__symmetric_mod_poly", "__taylor_coeff_zp",
    "__wang_core", "__wang_leading_coeff",
    # GCD 模块：稀疏 SQF 入口仍单独提取；polynomial_GCD 必须经
    # dense_upoly_zp 的显式 raw heap 闭包（StrictGCD/StrictHGCD）接回，禁止
    # 映射到数学规格或 L1 占位实现。详 strict GCD/HGCD source gates。
    "squarefreefactorize",
}
