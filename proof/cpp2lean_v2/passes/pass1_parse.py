"""
Pass 1: parse — Clang AST JSON → HIR₀

依据 `hir-design.md` §3。
输入：单个 FunctionDecl 节点（已经是实例化 mangledName≠None）
输出：HIRFunc（可能含 UnknownStmt/UnknownExpr）

设计要点：
- P1 原则：类型严格从 AST qualType 读取，不推断
- 不认识的 AST kind → UnknownStmt/UnknownExpr（不报错）
- `__assert_fail` 模式识别为 RequireStmt
- 运算符 → UnresolvedOp（后续 Pass 5 resolve）
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, ArrayType, PairType, TupleType, OptionType,
    StdMapType, RefType, UnknownType, TypeIR,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr, IteratorExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRParam, HIRFunc, TranslationError,
)


# ============================================================
# 类型解析（type-system.md §1-3）
# ============================================================

def parse_type(qt: str) -> TypeIR:
    """C++ qualType 字符串 → TypeIR。

    严格按 type-system.md 的映射表；未识别的返回 UnknownType（不报错）。
    """
    if not qt:
        return UnknownType("")
    qt = qt.strip()

    # 剥离 const / volatile 修饰（前导或后置）
    is_const = False
    if qt.startswith("const "):
        is_const = True
        qt = qt[len("const "):].strip()
    if qt.endswith(" const"):
        is_const = True
        qt = qt[:-len(" const")].strip()
    if qt.startswith("volatile "):
        qt = qt[len("volatile "):].strip()

    # 引用 / 指针
    if qt.endswith("&&"):
        inner = parse_type(qt[:-2].strip())
        return RefType(inner, is_const, is_rvalue=True)
    if qt.endswith("&"):
        inner = parse_type(qt[:-1].strip())
        return RefType(inner, is_const, is_rvalue=False)
    if qt.endswith("*"):
        # 指针（如 `const lex_<var_order>* comp_ptr`，CLPoly 仅 1 处）
        inner = parse_type(qt[:-1].strip())
        return RefType(inner, is_const, is_rvalue=False, is_pointer=True)

    # C 数组 char[N]（字符串字面量）
    import re as _re
    m = _re.match(r"^(.+)\[(\d+)\]$", qt)
    if m:
        elem_qt = m.group(1).strip()
        if elem_qt in ("char", "const char"):
            return NamedType("String")
        return ArrayType(parse_type(elem_qt))

    # 基础数值
    _BASE_MAP = {
        "uint64_t": BaseType.UINT64,
        "unsigned long": BaseType.UINT64,
        "unsigned long long": BaseType.UINT64,
        "int64_t": BaseType.INT64,
        "long": BaseType.INT64,
        "long long": BaseType.INT64,
        "uint32_t": BaseType.UINT32,
        "unsigned int": BaseType.UINT32,
        "unsigned": BaseType.UINT32,
        "int32_t": BaseType.INT32,
        "int": BaseType.INT32,
        "short": BaseType.INT32,
        "size_t": BaseType.NAT,
        "ptrdiff_t": BaseType.INT64,
        "bool": BaseType.BOOL,
        "_Bool": BaseType.BOOL,
        "char": BaseType.UINT32,
        "unsigned char": BaseType.UINT32,
        "double": BaseType.FLOAT,
        "float": BaseType.FLOAT,
        "void": BaseType.UNIT,
    }
    if qt in _BASE_MAP:
        return _BASE_MAP[qt]

    # STL 容器
    if qt.startswith("std::vector<") and qt.endswith(">"):
        elem = parse_type(_extract_template_arg(qt, "std::vector<"))
        return ArrayType(elem)
    if qt.startswith("std::pair<") and qt.endswith(">"):
        args = _extract_template_args(qt, "std::pair<")
        if len(args) == 2:
            return PairType(parse_type(args[0]), parse_type(args[1]))
    if qt.startswith("std::tuple<") and qt.endswith(">"):
        args = _extract_template_args(qt, "std::tuple<")
        return TupleType(tuple(parse_type(a) for a in args))
    if qt.startswith("std::map<") and qt.endswith(">"):
        args = _extract_template_args(qt, "std::map<")
        if len(args) >= 2:
            return StdMapType(parse_type(args[0]), parse_type(args[1]))
    if qt.startswith("std::list<") and qt.endswith(">"):
        # list → Array（CLPoly 不需要 linked list）
        elem = parse_type(_extract_template_arg(qt, "std::list<"))
        return ArrayType(elem)
    if qt.startswith("std::set<") and qt.endswith(">"):
        # set → Array（CLPoly 仅 __extract_monomial_content 用）
        elem = parse_type(_extract_template_arg(qt, "std::set<"))
        return ArrayType(elem)
    if qt.startswith("std::unordered_map<") and qt.endswith(">"):
        args = _extract_template_args(qt, "std::unordered_map<")
        if len(args) >= 2:
            return StdMapType(parse_type(args[0]), parse_type(args[1]))
    if qt.startswith("std::unordered_set<") and qt.endswith(">"):
        elem = parse_type(_extract_template_arg(qt, "std::unordered_set<"))
        return ArrayType(elem)

    # vector<...> / pair<...> / map<...> 无 std:: 前缀的情况（已去除 clpoly:: 后）
    if qt.startswith("vector<") and qt.endswith(">"):
        elem = parse_type(_extract_template_arg(qt, "vector<"))
        return ArrayType(elem)
    if qt.startswith("map<") and qt.endswith(">"):
        args = _extract_template_args(qt, "map<")
        if len(args) >= 2:
            return StdMapType(parse_type(args[0]), parse_type(args[1]))
    if qt.startswith("set<") and qt.endswith(">"):
        elem = parse_type(_extract_template_arg(qt, "set<"))
        return ArrayType(elem)

    # CLPoly 类型
    # 带完整 clpoly:: 前缀
    if qt.startswith("clpoly::"):
        qt = qt[len("clpoly::"):]

    _CLPOLY_MAP = {
        "ZZ": NamedType("ZZ"),
        "QQ": NamedType("QQ"),
        "Zp": NamedType("Zp"),
        "variable": NamedType("Variable"),
        "umonomial": NamedType("UMonomial"),
        "less": NamedType("Less"),  # MonomialOrder tag
        "uless": NamedType("ULess"),
    }
    if qt in _CLPOLY_MAP:
        return _CLPOLY_MAP[qt]

    # 模板化 CLPoly 类型
    if qt.startswith("upolynomial_<"):
        inner = _extract_template_arg(qt, "upolynomial_<")
        # upolynomial_<Zp> → SparsePolyZp
        if inner == "Zp":
            return NamedType("SparsePolyZp")
        if inner == "ZZ":
            return NamedType("SparsePolyZZ")
        return NamedType(f"SparsePoly_{inner}")

    if qt.startswith("polynomial_<"):
        # polynomial_<ZZ, lex_<less>> → MvPolyZZ
        args = _extract_template_args(qt, "polynomial_<")
        if args and args[0] == "ZZ":
            return NamedType("MvPolyZZ")
        if args and args[0] == "Zp":
            return NamedType("MvPolyZp")
        return NamedType(qt)

    if qt.startswith("basic_monomial<"):
        return NamedType("Monomial")

    if qt.startswith("basic_polynomial<"):
        # 底层表示；根据内部类型映射
        args = _extract_template_args(qt, "basic_polynomial<")
        if len(args) >= 2 and args[1] == "Zp":
            return NamedType("MvPolyZp") if "basic_monomial" in args[0] else NamedType("SparsePolyZp")
        if len(args) >= 2 and args[1] == "ZZ":
            return NamedType("MvPolyZZ") if "basic_monomial" in args[0] else NamedType("SparsePolyZZ")
        return NamedType(qt)

    if qt.startswith("factorization<"):
        return NamedType("Factorization")

    if qt.startswith("lex_<"):
        return NamedType("Lex")
    if qt.startswith("grlex_<"):
        return NamedType("Grlex")

    # CLPoly 辅助 struct（type-system.md §3.3.9）—— 含模板版本
    _CLPOLY_STRUCTS_PREFIX = {
        "__prime_selection_result": NamedType("PrimeSelectionResult"),
        "__wang_lc_result": NamedType("WangLcResult"),
        "__hensel_node": NamedType("HenselNode"),
        "__profile_stats": NamedType("ProfileStats"),
    }
    for prefix, ty in _CLPOLY_STRUCTS_PREFIX.items():
        if qt == prefix or qt.startswith(prefix + "<"):
            return ty

    # CLPoly 常用 typedef（出现在函数体内的 `using Poly = ...`）
    _TYPEDEF_ALIASES = {
        "Poly": NamedType("Poly"),          # 上下文依赖：通常是 polynomial_<ZZ, lex_<...>>
        "PolyZp": NamedType("PolyZp"),      # polynomial_<Zp, ...>
        "PolyQQ": NamedType("PolyQQ"),
        "PolyZZ": NamedType("PolyZZ"),
        "UPZp": NamedType("SparsePolyZp"),  # upolynomial_<Zp>
        "UPZZ": NamedType("SparsePolyZZ"),
        "LLLMatrix": NamedType("LLLMatrix"),
        "Matrix": NamedType("Matrix"),
    }
    if qt in _TYPEDEF_ALIASES:
        return _TYPEDEF_ALIASES[qt]

    # STL 迭代器 typedef（嵌套在容器内：value_type / size_type / iterator / const_iterator）
    _STL_NESTED_TYPES = {
        "size_type": BaseType.NAT,
        "difference_type": BaseType.INT64,
        "value_type": NamedType("ValueType"),       # 依容器具体化；HIR 阶段用占位
        "reference": NamedType("Reference"),
        "const_reference": NamedType("Reference"),
        "pointer": NamedType("Pointer"),
        "const_pointer": NamedType("Pointer"),
        "iterator": NamedType("Iterator"),
        "const_iterator": NamedType("Iterator"),
        "reverse_iterator": NamedType("ReverseIterator"),
        "allocator_type": NamedType("Allocator"),
        "allocator": NamedType("Allocator"),
        "mapped_type": NamedType("ValueType"),
        "key_type": NamedType("KeyType"),
        "result_type": NamedType("ResultType"),
        "argument_type": NamedType("ArgType"),
    }
    if qt in _STL_NESTED_TYPES:
        return _STL_NESTED_TYPES[qt]

    # STL 随机数三件套
    if qt in ("std::mt19937", "std::mersenne_twister_engine", "std::minstd_rand",
              "std::mt19937_64", "std::random_device"):
        return NamedType("Rng")
    if qt.startswith("std::uniform_int_distribution<"):
        return NamedType("UniformIntDist")
    if qt.startswith("std::uniform_real_distribution<"):
        return NamedType("UniformRealDist")
    if qt.startswith("std::bernoulli_distribution") or \
       qt.startswith("std::normal_distribution") or \
       qt.startswith("std::discrete_distribution"):
        return NamedType("Distribution")

    # std::pair<...> / std::tuple_element<...> 等 SFINAE 元编程的剩余物
    # 这些 qualType 带 typename / ::__type 后缀，统一归为 DependentType
    if (qt.startswith("std::pair<")) and ("typename" in qt or "::" in qt[-20:]):
        return NamedType("DependentPair")
    if qt.startswith("std::tuple_element<"):
        return NamedType("DependentType")
    if qt.startswith("pair<") and not qt.startswith("pair<std::"):  # 去前缀后的 pair<umonomial, ZZ>
        # 尝试解析两个参数
        args = _extract_template_args(qt, "pair<")
        if len(args) == 2:
            return PairType(parse_type(args[0]), parse_type(args[1]))

    # basic_monomial / basic_polynomial 的裸名（templates 已剥离的情况）
    if qt.startswith("basic_polynomial<basic_monomial<"):
        args = _extract_template_args(qt, "basic_polynomial<")
        if len(args) >= 2 and args[1] == "Zp":
            return NamedType("MvPolyZp")
        if len(args) >= 2 and args[1] == "ZZ":
            return NamedType("MvPolyZZ")
        return NamedType("MvPoly")

    # umonomial 裸名
    if qt == "umonomial":
        return NamedType("UMonomial")

    # 空 qualType — Clang JSON 偶尔产出，作为 unknown 默认
    if qt == "":
        return UnknownType("")

    # Lambda closure 类型：`(lambda at .../file.hh:NNN:CC)` — 匿名
    if qt.startswith("(lambda at ") and qt.endswith(")"):
        return NamedType("Lambda")

    # __int128
    if qt in ("__int128", "unsigned __int128"):
        return BaseType.UINT128

    # basic_string 等（仅 assert 消息用）
    if qt.startswith("std::basic_string") or qt.startswith("const char"):
        return NamedType("String")

    # typename T::... 形式（SFINAE 元编程剩余物）
    if qt.startswith("typename "):
        return NamedType("DependentType")

    # __gnu_cxx::__normal_iterator / std::__normal_iterator — C++ STL 内部实现
    if "__normal_iterator" in qt or "_Rb_tree" in qt or "_List_iterator" in qt:
        return NamedType("Iterator")

    # STL 内部 typedef（_Self / _Bit_reference / _Node_ptr 等实现细节）
    if qt.startswith("_") and len(qt) > 1:
        return NamedType("StlInternal")

    # STL 具体类型已在上面处理；其他未识别
    return UnknownType(qt)


def _extract_template_arg(qt: str, prefix: str) -> str:
    """提取单参数模板的唯一参数，如 `std::vector<T>` → `T`."""
    body = qt[len(prefix):-1].strip()
    # 取逗号前（若无逗号，整个）
    return _split_template_args(body)[0]


def _extract_template_args(qt: str, prefix: str) -> list[str]:
    """提取多参数模板的所有参数。"""
    body = qt[len(prefix):-1].strip()
    return _split_template_args(body)


def _split_template_args(s: str) -> list[str]:
    """按逗号分割，尊重嵌套 < > 对。"""
    result = []
    depth = 0
    start = 0
    for i, c in enumerate(s):
        if c == "<":
            depth += 1
        elif c == ">":
            depth -= 1
        elif c == "," and depth == 0:
            result.append(s[start:i].strip())
            start = i + 1
    if start < len(s):
        result.append(s[start:].strip())
    return result


def parse_param(parm_json: dict) -> HIRParam:
    """解析 ParmVarDecl → HIRParam。"""
    name = parm_json.get("name", "")
    qt = parm_json.get("type", {}).get("qualType", "") or \
         parm_json.get("type", {}).get("desugaredQualType", "")
    ty = parse_type(qt)

    # is_ref / is_const_ref 判断（指针单独处理）
    is_const_ref = False
    is_ref = False
    is_pointer = False
    if isinstance(ty, RefType):
        if ty.is_pointer:
            is_pointer = True
        elif ty.is_const:
            is_const_ref = True
        else:
            is_ref = True
        ty = ty.inner

    return HIRParam(
        name=name,
        ty=ty,
        is_ref=is_ref,
        is_const_ref=is_const_ref,
        is_output=is_ref,  # 指针不视作 output（1 处是 const 指针）
    )


def parse_return_type(qt: str) -> TypeIR:
    """从函数 qualType 提取返回类型。
    qualType 形如 `RetType (Param1, Param2)`。"""
    if "(" not in qt:
        return parse_type(qt)
    ret = qt[:qt.index("(")].strip()
    return parse_type(ret)


# ============================================================
# 表达式解析
# ============================================================

def parse_expr(node: Any) -> ExprIR:
    """AST 节点 → ExprIR。不认识的返回 UnknownExpr。"""
    if not isinstance(node, dict):
        return UnknownExpr("not_a_dict", raw=str(node)[:100])

    kind = node.get("kind", "")
    qt = node.get("type", {}).get("qualType", "") or \
         node.get("type", {}).get("desugaredQualType", "")
    ty = parse_type(qt) if qt else None

    if kind == "DeclRefExpr":
        rd = node.get("referencedDecl", {}) or {}
        name = rd.get("name", "<unknown>")
        return Var(name=name, version=0, ty=ty)

    if kind == "IntegerLiteral":
        raw_value = node.get("value", "0")
        try:
            value = int(raw_value)
        except (ValueError, TypeError):
            value = 0
        return Lit(value=value, ty=ty if ty else BaseType.INT32)

    if kind == "CXXBoolLiteralExpr":
        value = node.get("value", False)
        return Lit(value=bool(value), ty=BaseType.BOOL)

    if kind == "FloatingLiteral":
        raw_value = node.get("value", "0.0")
        try:
            value = float(raw_value)
        except (ValueError, TypeError):
            value = 0.0
        return Lit(value=value, ty=BaseType.FLOAT)

    if kind == "StringLiteral":
        return Lit(value=node.get("value", ""), ty=NamedType("String"))

    if kind == "BinaryOperator":
        op = node.get("opcode", "?")
        inner = node.get("inner", [])
        if len(inner) >= 2:
            lhs = parse_expr(inner[0])
            rhs = parse_expr(inner[1])
            return BinOp(op=op, lhs=lhs, rhs=rhs, ty=ty)
        return UnknownExpr(kind, children=[parse_expr(c) for c in inner])

    if kind == "UnaryOperator":
        op = node.get("opcode", "?")
        inner = node.get("inner", [])
        operand = parse_expr(inner[0]) if inner else UnknownExpr("empty_unary")
        return UnaryOp(op=op, operand=operand, ty=ty)

    if kind == "ConditionalOperator":
        inner = node.get("inner", [])
        if len(inner) >= 3:
            return CondExpr(
                cond=parse_expr(inner[0]),
                then_e=parse_expr(inner[1]),
                else_e=parse_expr(inner[2]),
                ty=ty,
            )
        return UnknownExpr(kind, children=[parse_expr(c) for c in inner])

    if kind == "CallExpr":
        inner = node.get("inner", [])
        if not inner:
            return UnknownExpr(kind)
        callee_node = inner[0]
        callee_name = _extract_callee_name(callee_node)
        args = [parse_expr(a) for a in inner[1:]]
        return Call(callee=callee_name, args=args, ty=ty)

    if kind == "CXXOperatorCallExpr":
        inner = node.get("inner", [])
        if not inner:
            return UnknownExpr(kind)
        # 第一个 inner 是 operator function ref，提取 op 名字
        op_name = _extract_callee_name(inner[0])
        if op_name is None or not op_name.startswith("operator"):
            op_name = op_name or "<unknown-op>"
        args = [parse_expr(a) for a in inner[1:]]
        receiver_ty = args[0].ty if args and hasattr(args[0], "ty") else None
        return Call(
            callee=UnresolvedOp(op_name=op_name, receiver_ty=receiver_ty),
            args=args,
            ty=ty,
        )

    if kind == "CXXMemberCallExpr":
        inner = node.get("inner", [])
        if not inner:
            return UnknownExpr(kind)
        # 第一个 inner 是 MemberExpr(obj, method_name)
        mem = inner[0]
        method_name = _extract_member_name(mem)
        obj_expr = _extract_member_object(mem)
        obj = parse_expr(obj_expr) if obj_expr else UnknownExpr("empty-obj")
        args = [parse_expr(a) for a in inner[1:]]
        # 构造 receiver type
        receiver_ty = obj.ty if hasattr(obj, "ty") else None
        # 编码为 Call(UnresolvedOp("TypeName.method"), [obj, *args])
        op_name = f"<method>.{method_name}"
        return Call(
            callee=UnresolvedOp(op_name=op_name, receiver_ty=receiver_ty),
            args=[obj] + args,
            ty=ty,
        )

    if kind == "ArraySubscriptExpr":
        inner = node.get("inner", [])
        if len(inner) >= 2:
            return ArrayAccess(
                arr=parse_expr(inner[0]),
                idx=parse_expr(inner[1]),
                ty=ty,
            )
        return UnknownExpr(kind)

    if kind == "MemberExpr":
        inner = node.get("inner", [])
        obj = parse_expr(inner[0]) if inner else UnknownExpr("empty-obj")
        field_name = node.get("name", "")
        return FieldAccess(obj=obj, field_name=field_name, ty=ty)

    if kind in ("ImplicitCastExpr", "CStyleCastExpr",
                "CXXStaticCastExpr", "CXXFunctionalCastExpr"):
        inner = node.get("inner", [])
        if not inner:
            return UnknownExpr(kind)
        sub = parse_expr(inner[0])
        source_qt = inner[0].get("type", {}).get("qualType", "") if isinstance(inner[0], dict) else ""
        source_ty = parse_type(source_qt) if source_qt else UnknownType("")
        cast_kind = node.get("castKind", "?")
        target_ty = ty if ty else UnknownType("")
        return Cast(expr=sub, source_ty=source_ty, target_ty=target_ty, cast_kind=cast_kind)

    if kind in ("CXXConstructExpr", "CXXTemporaryObjectExpr"):
        # T(args...) 形式（临时对象用 CXXTemporaryObjectExpr）
        args = [parse_expr(c) for c in node.get("inner", [])]
        return Call(
            callee=UnresolvedOp(op_name=f"construct_{qt}", receiver_ty=ty),
            args=args,
            ty=ty,
        )

    if kind == "ImplicitValueInitExpr":
        # 默认初始化：T() 隐式
        # 返回 Lit(0) 或 Call(construct)，依上下文；简化为 construct 调用
        return Call(
            callee=UnresolvedOp(op_name=f"default_init_{qt}", receiver_ty=ty),
            args=[],
            ty=ty,
        )

    if kind == "CXXBindTemporaryExpr" or kind == "MaterializeTemporaryExpr" or kind == "ExprWithCleanups":
        # 透明包装：取 inner[0]
        inner = node.get("inner", [])
        if inner:
            return parse_expr(inner[0])
        return UnknownExpr(kind)

    if kind == "ParenExpr":
        inner = node.get("inner", [])
        if inner:
            return parse_expr(inner[0])
        return UnknownExpr(kind)

    if kind == "CXXThisExpr":
        return Var(name="this", version=0, ty=ty)

    if kind == "CXXDefaultArgExpr":
        # 默认参数：inner[0] 是原表达式，若无则用类型默认值
        inner = node.get("inner", [])
        if inner:
            return parse_expr(inner[0])
        # 无 inner：回退为 default_init
        return Call(
            callee=UnresolvedOp(op_name=f"default_init_{qt}", receiver_ty=ty),
            args=[],
            ty=ty,
        )

    if kind == "LambdaExpr":
        # Clang JSON LambdaExpr 的 inner：
        #   [0] CXXRecordDecl (closure class) — 包含 operator() CXXMethodDecl
        #   [1] CompoundStmt (lambda body)
        # captures 字段 Clang JSON 略过，由 Pass 3 lambda_lift 从源码/AST 上下文推导
        lam_params: list[HIRParam] = []
        lam_body: list[StmtIR] = []
        for c in node.get("inner", []):
            if not isinstance(c, dict):
                continue
            if c.get("kind") == "CXXRecordDecl":
                # 找 operator() 方法
                for cc in c.get("inner", []):
                    if not isinstance(cc, dict):
                        continue
                    cc_kind = cc.get("kind")
                    # 具体 operator()
                    if cc_kind == "CXXMethodDecl" and cc.get("name") == "operator()":
                        for ccc in cc.get("inner", []):
                            if isinstance(ccc, dict) and ccc.get("kind") == "ParmVarDecl":
                                lam_params.append(parse_param(ccc))
                    # 模板 operator()（generic lambda，CLPoly 已消除；兜底）
                    elif cc_kind == "FunctionTemplateDecl":
                        for ccc in cc.get("inner", []):
                            if isinstance(ccc, dict) and ccc.get("kind") == "CXXMethodDecl" \
                                    and ccc.get("name") == "operator()":
                                for cccc in ccc.get("inner", []):
                                    if isinstance(cccc, dict) and cccc.get("kind") == "ParmVarDecl":
                                        lam_params.append(parse_param(cccc))
            elif c.get("kind") == "CompoundStmt":
                lam_body = _parse_stmts(c)
        return LambdaExpr(
            captures=[],  # 由 Pass 3 lambda_lift 从源码正则 + 自由变量分析推导
            params=lam_params,
            body=lam_body,
            ty=ty,
        )

    if kind == "InitListExpr":
        elems = [parse_expr(c) for c in node.get("inner", [])]
        elem_ty = elems[0].ty if elems and hasattr(elems[0], "ty") else UnknownType("")
        return ArrayLit(elems=elems, elem_ty=elem_ty)

    if kind == "CXXStdInitializerListExpr":
        inner = node.get("inner", [])
        if inner:
            return parse_expr(inner[0])
        return UnknownExpr(kind)

    if kind == "UnresolvedLookupExpr":
        name = node.get("name", "<unresolved>")
        return Var(name=name, version=0, ty=ty)

    # 忽略调试/元信息节点
    if kind in ("SourceLocExpr", "PredefinedExpr"):
        return Lit(value="", ty=NamedType("String"))

    # 未识别
    return UnknownExpr(
        kind=kind,
        children=[parse_expr(c) for c in node.get("inner", [])],
        raw=str(node.get("name", ""))[:50],
    )


def _extract_callee_name(node: Any) -> str | None:
    """从 CallExpr 的 callee 节点提取函数名。"""
    if not isinstance(node, dict):
        return None
    kind = node.get("kind", "")
    if kind == "DeclRefExpr":
        rd = node.get("referencedDecl", {}) or {}
        return rd.get("name")
    if kind == "ImplicitCastExpr":
        for c in node.get("inner", []):
            name = _extract_callee_name(c)
            if name:
                return name
    if kind == "MemberExpr":
        return node.get("name")
    if kind == "UnresolvedLookupExpr":
        return node.get("name")
    return None


def _extract_member_name(node: Any) -> str:
    """从 MemberExpr 节点提取方法/字段名。"""
    if isinstance(node, dict):
        return node.get("name", "")
    return ""


def _extract_member_object(node: Any) -> Any:
    """从 MemberExpr 提取 object 部分（inner[0]）。"""
    if isinstance(node, dict):
        inner = node.get("inner", [])
        if inner:
            return inner[0]
    return None


# ============================================================
# 语句解析
# ============================================================

def parse_stmt(node: Any) -> StmtIR | list[StmtIR]:
    """AST 节点 → StmtIR 或 list[StmtIR]（DeclStmt 可能产生多个 LetStmt）。

    不认识的 → UnknownStmt（不报错）。
    """
    if not isinstance(node, dict):
        return UnknownStmt("not_a_dict")

    kind = node.get("kind", "")

    if kind == "CompoundStmt":
        stmts = []
        for c in node.get("inner", []):
            s = parse_stmt(c)
            if isinstance(s, list):
                stmts.extend(s)
            else:
                stmts.append(s)
        return BlockStmt(stmts=stmts)

    if kind == "DeclStmt":
        result = []
        for c in node.get("inner", []):
            if not isinstance(c, dict):
                continue
            child_kind = c.get("kind")
            if child_kind == "VarDecl":
                name = c.get("name", "")
                qt = c.get("type", {}).get("qualType", "") or \
                     c.get("type", {}).get("desugaredQualType", "")
                ty = parse_type(qt)
                # 初始值：inner[0] if present（跳过 TypeLoc 等元信息）
                init = None
                for cc in c.get("inner", []):
                    if isinstance(cc, dict) and cc.get("kind") not in ("TypeLoc",):
                        init = parse_expr(cc)
                        break
                if init is None:
                    # 无初始化：用 default_init
                    init = Call(
                        callee=UnresolvedOp(op_name=f"default_init_{qt}", receiver_ty=ty),
                        args=[],
                        ty=ty,
                    )
                result.append(LetStmt(var=Var(name), ty=ty, value=init))
            elif child_kind == "DecompositionDecl":
                # 结构化绑定：auto [a, b] = pair_expr;
                # 拆为：let _dec := expr; let a := _dec.fst; let b := _dec.snd
                qt = c.get("type", {}).get("qualType", "")
                ty = parse_type(qt)

                # 推导 binding 类型（从 pair/tuple 的 template 参数拆出）
                binding_types: list[TypeIR] = []
                if isinstance(ty, PairType):
                    binding_types = [ty.fst, ty.snd]
                elif isinstance(ty, TupleType):
                    binding_types = list(ty.elems)
                # RefType(PairType) 也要剥一层
                elif isinstance(ty, RefType) and isinstance(ty.inner, PairType):
                    binding_types = [ty.inner.fst, ty.inner.snd]
                elif isinstance(ty, RefType) and isinstance(ty.inner, TupleType):
                    binding_types = list(ty.inner.elems)

                init = None
                bindings = []  # list[(name, TypeIR)]
                idx = 0
                for cc in c.get("inner", []):
                    if isinstance(cc, dict):
                        if cc.get("kind") == "BindingDecl":
                            bname = cc.get("name", "")
                            btype_qt = ""
                            t_node = cc.get("type")
                            if isinstance(t_node, dict):
                                btype_qt = t_node.get("qualType", "") or t_node.get("desugaredQualType", "")
                            if btype_qt:
                                bty = parse_type(btype_qt)
                            elif idx < len(binding_types):
                                bty = binding_types[idx]
                            else:
                                bty = UnknownType("")
                            bindings.append((bname, bty))
                            idx += 1
                        elif init is None and cc.get("kind") not in ("TypeLoc",):
                            init = parse_expr(cc)
                if init is None:
                    init = UnknownExpr("decomp_uninit")
                # 先 let 整个解构目标
                decomp_var = f"__decomp_{id(c) & 0xfffff}"
                result.append(LetStmt(
                    var=Var(decomp_var),
                    ty=ty,
                    value=init,
                ))
                # 再按 binding 序分别 let
                for i, (bname, bty) in enumerate(bindings):
                    result.append(LetStmt(
                        var=Var(bname),
                        ty=bty,
                        value=FieldAccess(
                            obj=Var(decomp_var),
                            field_name=f"fst" if i == 0 else f"snd" if i == 1 else f"_{i}",
                        ),
                    ))
            elif child_kind in ("UsingDecl", "TypeAliasDecl"):
                # 忽略：这些在翻译阶段不产生运行时代码
                pass
            else:
                result.append(UnknownStmt(f"decl_{child_kind}"))
        return result if len(result) != 1 else result[0]

    if kind == "IfStmt":
        # 特殊：识别 `if (!cond) __assert_fail(...)` 模式 → RequireStmt
        inner = node.get("inner", [])
        if len(inner) >= 2:
            cond_node = inner[0]
            then_node = inner[1]
            require_stmt = _try_parse_assert_pattern(cond_node, then_node)
            if require_stmt is not None:
                return require_stmt
            cond = parse_expr(cond_node)
            then_body = _parse_stmts(then_node)
            else_body = _parse_stmts(inner[2]) if len(inner) >= 3 else []
            return IfStmt(cond=cond, then_body=then_body, else_body=else_body)
        return UnknownStmt(kind)

    if kind == "WhileStmt":
        inner = node.get("inner", [])
        if len(inner) >= 2:
            cond = parse_expr(inner[0])
            body = _parse_stmts(inner[1])
            return WhileStmt(cond=cond, body=body)
        return UnknownStmt(kind)

    if kind == "ForStmt":
        # Clang JSON: ForStmt inner = [init, null_placeholder?, cond, step, body]
        inner = node.get("inner", [])
        # 标准位置：inner[0]=init, [1]=null/placeholder, [2]=cond, [3]=step, [4]=body
        init_node = inner[0] if len(inner) > 0 else None
        cond_node = inner[2] if len(inner) > 2 else None
        step_node = inner[3] if len(inner) > 3 else None
        body_node = inner[4] if len(inner) > 4 else None
        init = _parse_stmts(init_node) if init_node else []
        cond = parse_expr(cond_node) if cond_node else Lit(True, BaseType.BOOL)
        step = [ExprStmt(parse_expr(step_node))] if step_node else []
        body = _parse_stmts(body_node) if body_node else []
        return ForStmt(init=init, cond=cond, step=step, body=body)

    if kind == "CXXForRangeStmt":
        # Clang AST inner 约定：
        #   [0] init stmt（通常是 null placeholder）
        #   [1] DeclStmt __range1（container init）
        #   [2] DeclStmt __begin1
        #   [3] DeclStmt __end1
        #   [4] 循环条件（__begin1 != __end1）
        #   [5] 循环步进（++__begin1）
        #   [6] DeclStmt 循环变量 var（或 DecompositionDecl 结构化绑定）
        #   [7] body（可能是 CompoundStmt 或 单个表达式/语句）
        inner = node.get("inner", [])
        var = Var("__range_var")
        var_ty: TypeIR = UnknownType("")
        container: ExprIR = UnknownExpr("container")
        body: list[StmtIR] = []
        decomposition: list[Var] | None = None

        def extract_varDecl_or_decomp(declstmt: dict):
            """从 DeclStmt 里找 VarDecl 或 DecompositionDecl。"""
            for cc in declstmt.get("inner", []):
                if isinstance(cc, dict) and cc.get("kind") in ("VarDecl", "DecompositionDecl"):
                    return cc
            return None

        def find_init_expr(vardecl: dict):
            """从 VarDecl 的 inner 里找初始化表达式（第一个非 child Decl）。"""
            for ccc in vardecl.get("inner", []):
                if isinstance(ccc, dict):
                    # 跳过 TypeLoc 等元信息
                    if ccc.get("kind") not in ("TypeLoc", "NamespaceSpec"):
                        return ccc
            return None

        # inner[1]: __range1 DeclStmt → container
        if len(inner) > 1 and isinstance(inner[1], dict) and inner[1].get("kind") == "DeclStmt":
            vd = extract_varDecl_or_decomp(inner[1])
            if vd and vd.get("kind") == "VarDecl":
                init = find_init_expr(vd)
                if init:
                    container = parse_expr(init)

        # inner[6]: 循环变量或结构化绑定
        if len(inner) > 6 and isinstance(inner[6], dict) and inner[6].get("kind") == "DeclStmt":
            vd = extract_varDecl_or_decomp(inner[6])
            if vd:
                if vd.get("kind") == "DecompositionDecl":
                    decomposition = [
                        Var(b.get("name", ""))
                        for b in vd.get("inner", [])
                        if isinstance(b, dict) and b.get("kind") == "BindingDecl"
                    ]
                    # DecompositionDecl 自身的 qualType 是 pair 类型
                    qt = vd.get("type", {}).get("qualType", "")
                    var_ty = parse_type(qt)
                    var = Var("__decomp")  # 占位，解构由 Pass 4 完成
                elif vd.get("kind") == "VarDecl":
                    nm = vd.get("name", "")
                    qt = vd.get("type", {}).get("qualType", "")
                    var = Var(nm)
                    var_ty = parse_type(qt)

        # inner[7]: body（可能是 CompoundStmt 或单语句）
        if len(inner) > 7 and isinstance(inner[7], dict):
            body = _parse_stmts(inner[7])

        return RangeForStmt(
            var=var, var_ty=var_ty,
            container=container, body=body,
            decomposition=decomposition,
        )

    if kind == "DoStmt":
        inner = node.get("inner", [])
        if len(inner) >= 2:
            body = _parse_stmts(inner[0])
            cond = parse_expr(inner[1])
            return DoWhileStmt(body=body, cond=cond)
        return UnknownStmt(kind)

    if kind == "BreakStmt":
        return BreakStmt()

    if kind == "ContinueStmt":
        return ContinueStmt()

    if kind == "ReturnStmt":
        inner = node.get("inner", [])
        value = parse_expr(inner[0]) if inner else None
        return ReturnStmt(value=value)

    if kind == "NullStmt":
        return BlockStmt(stmts=[])

    # assert 的三元表达式形式：
    #   (cond ? (void)0 : __assert_fail(...))
    # Clang 预处理把 `assert(cond)` 展开成这种形式，作为 ExprStmt 出现。
    # 可能被 CStyleCastExpr(void) / ExprWithCleanups / CXXBindTemporaryExpr / ParenExpr 等包裹。
    unwrapped = _unwrap_expr(node)
    if isinstance(unwrapped, dict) and unwrapped.get("kind") == "ConditionalOperator":
        assert_req = _try_parse_assert_ternary(unwrapped)
        if assert_req is not None:
            return assert_req

    # 表达式语句
    if kind in ("BinaryOperator", "UnaryOperator", "CallExpr",
                "CXXOperatorCallExpr", "CXXMemberCallExpr", "CompoundAssignOperator"):
        # 作为顶层语句 → 识别为 AssignStmt/CompoundAssignStmt 或 ExprStmt
        if kind == "BinaryOperator" and node.get("opcode") == "=":
            inner = node.get("inner", [])
            if len(inner) >= 2:
                lhs = parse_expr(inner[0])
                rhs = parse_expr(inner[1])
                return AssignStmt(target=lhs, value=rhs)
        if kind == "CompoundAssignOperator":
            op = node.get("opcode", "?")  # "+=", "-=", ...
            base_op = op[:-1] if len(op) > 1 else op  # "+=" → "+"
            inner = node.get("inner", [])
            if len(inner) >= 2:
                return CompoundAssignStmt(
                    target=parse_expr(inner[0]),
                    op=base_op,
                    value=parse_expr(inner[1]),
                )
        # 其他表达式 → ExprStmt
        return ExprStmt(expr=parse_expr(node))

    # 其他表达式（作为语句）
    if kind in ("ImplicitCastExpr", "CStyleCastExpr",
                "CXXStaticCastExpr", "CXXBindTemporaryExpr",
                "MaterializeTemporaryExpr", "ExprWithCleanups",
                "ParenExpr", "CXXConstructExpr", "ConditionalOperator",
                "DeclRefExpr", "IntegerLiteral", "StringLiteral",
                "MemberExpr", "ArraySubscriptExpr"):
        return ExprStmt(expr=parse_expr(node))

    # 未识别
    children = []
    for c in node.get("inner", []):
        s = parse_stmt(c)
        if isinstance(s, list):
            children.extend(s)
        else:
            children.append(s)
    return UnknownStmt(kind=kind, children=children)


def _parse_stmts(node: Any) -> list[StmtIR]:
    """解析一个 node（通常是 CompoundStmt）为 stmt list。
    单语句 body（不在 CompoundStmt 内）也处理。"""
    if not isinstance(node, dict):
        return []
    if node.get("kind") == "CompoundStmt":
        result = []
        for c in node.get("inner", []):
            s = parse_stmt(c)
            if isinstance(s, list):
                result.extend(s)
            else:
                result.append(s)
        return result
    # 单语句
    s = parse_stmt(node)
    return s if isinstance(s, list) else [s]


# ============================================================
# Assert 模式识别
# ============================================================

def _try_parse_assert_pattern(cond_node: Any, then_node: Any) -> RequireStmt | None:
    """识别 `if (!cond) __assert_fail(...)` → RequireStmt(cond, ...).

    C++ assert(cond) 在 debug 编译时展开为 if (!cond) __assert_fail(msg,...)。
    """
    if not isinstance(then_node, dict):
        return None
    # then_node 可能是 CompoundStmt 包装单个 __assert_fail 调用，或直接 CallExpr
    assert_call = _find_assert_fail_call(then_node)
    if assert_call is None:
        return None

    # cond 应该是 UnaryOperator(!, cond_inner)
    if not isinstance(cond_node, dict):
        return None
    if cond_node.get("kind") != "UnaryOperator" or cond_node.get("opcode") != "!":
        return None
    inner = cond_node.get("inner", [])
    if not inner:
        return None
    true_cond = parse_expr(inner[0])
    return RequireStmt(
        cond=true_cond,
        name="h_assert",
        source="assert",
    )


def _find_assert_fail_call(node: Any) -> dict | None:
    """递归查找 `__assert_fail(...)` CallExpr。"""
    if not isinstance(node, dict):
        return None
    kind = node.get("kind", "")
    if kind == "CallExpr":
        inner = node.get("inner", [])
        if inner:
            callee_name = _extract_callee_name(inner[0])
            if callee_name == "__assert_fail":
                return node
    if kind in ("CompoundStmt", "BlockStmt",
                # 经常包装 __assert_fail 的节点
                "ImplicitCastExpr", "CStyleCastExpr",
                "ExprWithCleanups", "CXXBindTemporaryExpr",
                "MaterializeTemporaryExpr", "ParenExpr"):
        for c in node.get("inner", []):
            found = _find_assert_fail_call(c)
            if found:
                return found
    return None


_UNWRAP_KINDS = {
    "CStyleCastExpr", "CXXBindTemporaryExpr", "MaterializeTemporaryExpr",
    "ExprWithCleanups", "ParenExpr", "ImplicitCastExpr",
}


def _unwrap_expr(node: Any) -> Any:
    """剥离对语义无贡献的 wrapper（cast 到 void、临时对象绑定等），
    返回最内层的有意义节点。用于识别嵌套 wrapper 下的 assert 三元。"""
    current = node
    while isinstance(current, dict) and current.get("kind") in _UNWRAP_KINDS:
        inner = current.get("inner", [])
        if not inner:
            break
        current = inner[0]
    return current


def _try_parse_assert_ternary(cond_op_node: dict) -> RequireStmt | None:
    """识别 `cond ? (void)0 : __assert_fail(...)` ConditionalOperator
    作为 assert 展开的等价形式 → RequireStmt(cond)。

    ConditionalOperator 的 inner 有 3 个子节点：[cond, then_e, else_e]。
    模式：else_e 含 `__assert_fail(...)` 调用。
    """
    if not isinstance(cond_op_node, dict):
        return None
    if cond_op_node.get("kind") != "ConditionalOperator":
        return None
    inner = cond_op_node.get("inner", [])
    if len(inner) < 3:
        return None
    # else 分支递归查 __assert_fail
    else_branch = inner[2]
    assert_call = _find_assert_fail_call(else_branch)
    if assert_call is None:
        return None
    # then 分支应该是 (void)0 或类似 void 操作；不强制检查
    cond_expr = parse_expr(inner[0])
    return RequireStmt(
        cond=cond_expr,
        name="h_assert",
        source="assert",
    )


# ============================================================
# Pass 1 入口
# ============================================================

def parse_pass(ast_json: dict) -> HIRFunc:
    """AST FunctionDecl 节点 → HIRFunc（HIR₀ 阶段）。

    前提：ast_json 是一个实例化的 FunctionDecl（mangledName 非空）。
    """
    if not isinstance(ast_json, dict):
        raise TranslationError("parse", "<unknown>", "AST not a dict")
    if ast_json.get("kind") != "FunctionDecl":
        raise TranslationError(
            "parse", ast_json.get("name", "<unknown>"),
            f"expected FunctionDecl, got {ast_json.get('kind')}",
        )

    base_name = ast_json.get("name", "")
    mangled_name = ast_json.get("mangledName", "")
    qual_type = ast_json.get("type", {}).get("qualType", "")

    # 解析参数
    params = []
    body_node: Any = None
    for c in ast_json.get("inner", []):
        if not isinstance(c, dict):
            continue
        k = c.get("kind")
        if k == "ParmVarDecl":
            params.append(parse_param(c))
        elif k == "CompoundStmt":
            body_node = c

    # 返回类型
    ret_ty = parse_return_type(qual_type)

    # body
    if body_node:
        body_stmt = parse_stmt(body_node)
        if isinstance(body_stmt, BlockStmt):
            body = body_stmt.stmts
        elif isinstance(body_stmt, list):
            body = body_stmt
        else:
            body = [body_stmt]
    else:
        body = []

    # 实例化 suffix（基于 qual_type）
    instance_suffix = _infer_instance_suffix(qual_type)

    return HIRFunc(
        base_name=base_name,
        instance_suffix=instance_suffix,
        mangled_name=mangled_name,
        qual_type=qual_type,
        params=params,
        ret_ty=ret_ty,
        body=body,
    )


def _infer_instance_suffix(qual_type: str) -> str:
    """根据 qualType 推断 factorize 3 实例的后缀。其他函数无后缀。"""
    if not qual_type:
        return ""
    if "upolynomial_<ZZ>" in qual_type:
        return "upoly"
    # polynomial_<ZZ, grlex_<less>> 先检查（比 lex_<less> 更具体）
    if "grlex_<less>" in qual_type:
        return "grlex"
    if "lex_<less>" in qual_type:
        return "lex"
    return ""


# ============================================================
# 不变量检查（HIR₀ — 宽松）
# ============================================================

# IR 节点的联合类型（仅 HIR Stmt/Expr，不递归到 TypeIR 防止循环）
_IR_NODE_TYPES = (
    # Stmt 家族
    LetStmt, AssignStmt, CompoundAssignStmt,
    IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
    BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt,
    # Expr 家族
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr, IteratorExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr,
)


def _walk_ir(node, visitor):
    """只递归 IR Stmt/Expr 节点（不追 TypeIR 防循环）。"""
    if not isinstance(node, _IR_NODE_TYPES):
        return
    visitor(node)
    # 通过反射收集子 IR 节点
    for fname in node.__dataclass_fields__:
        v = getattr(node, fname)
        if isinstance(v, _IR_NODE_TYPES):
            _walk_ir(v, visitor)
        elif isinstance(v, list):
            for x in v:
                _walk_ir(x, visitor)
        elif isinstance(v, dict):
            for x in v.values():
                _walk_ir(x, visitor)


def assert_hir0_invariant(func: HIRFunc):
    """HIR₀ 出口 assert。
    允许 UnknownStmt/UnknownExpr，但打印数量警告。"""
    unknown_count = 0

    def visit(node):
        nonlocal unknown_count
        if isinstance(node, (UnknownStmt, UnknownExpr)):
            unknown_count += 1

    for s in func.body:
        _walk_ir(s, visit)

    if unknown_count > 0:
        import sys as _sys
        print(
            f"[parse] WARNING: {func.base_name} has {unknown_count} Unknown nodes",
            file=_sys.stderr,
        )
