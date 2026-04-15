"""
CLPoly C++ → Lean IR 翻译器：Clang AST JSON 解析器

输入：clang++ -Xclang -ast-dump=json -fsyntax-only 产出的 JSON
输出：list[FuncIR]
"""

from __future__ import annotations
import json
import re
from ir_types import *


def _get_qual_type(node: dict) -> str:
    """从 AST 节点获取类型字符串。优先 desugaredQualType（展开 auto/typedef/value_type）。"""
    type_info = node.get("type", {})
    desugar = type_info.get("desugaredQualType", "")
    qual = type_info.get("qualType", "")
    # 优先用 desugaredQualType：展开 auto、typedef、value_type 等别名
    if desugar and desugar != qual:
        return desugar
    return qual or desugar


def _get_class_name(node: dict) -> str | None:
    """从 Clang AST 节点提取 CLPoly 类名。parse_type 统一处理所有类型。"""
    qual = _get_qual_type(node)
    typ = parse_type(qual)
    if isinstance(typ, StructType):
        return typ.name
    if isinstance(typ, ArrayType):
        # ArrayType 的 elem 决定了对应的 CLASS_MAP 类名
        from class_map import CLASS_MAP
        for cls_name, cls_info in CLASS_MAP.items():
            if isinstance(cls_info.get("lean_type"), StructType):
                # SparsePolyZp = Array (UMonomial × Zp)
                if cls_name == "SparsePolyZp":
                    return cls_name
        return None
    return None


def _types_equal(a: TypeIR, b: TypeIR) -> bool:
    """判断两个 TypeIR 是否等价（用于 copy constructor 检测）。"""
    if type(a) != type(b):
        return False
    if isinstance(a, BaseType):
        return a == b
    if isinstance(a, StructType):
        return a.name == b.name
    if isinstance(a, ArrayType):
        return _types_equal(a.elem, b.elem)
    if isinstance(a, PairType):
        return _types_equal(a.fst, b.fst) and _types_equal(a.snd, b.snd)
    return a == b


def _match_ctor(ctors: dict, arg_types: tuple) -> str | None:
    """匹配构造函数：先按参数类型精确匹配，再按参数数量匹配。"""
    # 精确类型匹配
    for key, val in ctors.items():
        if isinstance(key, tuple) and len(key) == len(arg_types):
            if all(_types_equal(k, a) if not isinstance(k, str) else False
                   for k, a in zip(key, arg_types)):
                return val
    # 按数量匹配（fallback）
    for key, val in ctors.items():
        if isinstance(key, tuple) and len(key) == len(arg_types):
            return val
    return None


def _extract_callee_name(node: dict) -> str:
    """从 AST 节点（可能是 ImplicitCastExpr 链）中提取被引用的函数/运算符名。"""
    # 直接有 referencedDecl
    ref = node.get("referencedDecl", {})
    if ref.get("name"):
        return ref["name"]
    # 展开 ImplicitCastExpr / CXXBindTemporaryExpr 等包裹层
    for child in node.get("inner", []):
        result = _extract_callee_name(child)
        if result:
            return result
    return ""


def _parse_output_params_from_sig(sig: str) -> list[int]:
    """从函数签名字符串解析输出参数位置（非 const 引用 = 输出）。

    例：'void (ZZ &, const ZZ &, const ZZ &)' → [0]
    """
    m = re.search(r'\((.+)\)\s*$', sig)
    if not m:
        return []
    params_str = m.group(1)
    # 按顶层逗号分割（不分割 <> 内的逗号）
    params = []
    depth = 0
    current = ''
    for ch in params_str:
        if ch in '<(':
            depth += 1
        elif ch in '>)':
            depth -= 1
        if ch == ',' and depth == 0:
            params.append(current.strip())
            current = ''
        else:
            current += ch
    if current.strip():
        params.append(current.strip())
    output = []
    for i, p in enumerate(params):
        if '&' in p:
            before_amp = p[:p.rindex('&')]
            if 'const' not in before_amp:
                output.append(i)
    return output


def _get_callee_sig(callee_node: dict) -> str:
    """从 CallExpr 的 callee 节点获取被调函数签名。"""
    ref = callee_node
    while ref.get("kind") in ("ImplicitCastExpr",):
        inner = ref.get("inner", [])
        ref = inner[0] if inner else ref
        break
    rd = ref.get("referencedDecl", {})
    return rd.get("type", {}).get("qualType", "")


# Lambda body 全局注册表：AST 层填充，SSA 层消费
_LAMBDA_REGISTRY: dict[int, tuple[list[str], list, list]] = {}
_LAMBDA_COUNTER = 0


def _extract_root_var(expr) -> Var | None:
    """从 Cast/FieldAccess 链中提取根变量。"""
    if isinstance(expr, Var):
        return expr
    if isinstance(expr, Cast):
        return _extract_root_var(expr.expr)
    if isinstance(expr, FieldAccess):
        return _extract_root_var(expr.obj)
    return None


def _extract_func_name(expr) -> str:
    """从表达式中提取函数名。处理 Var、Cast(Var)、等。"""
    if isinstance(expr, Var):
        return expr.name
    if isinstance(expr, Cast):
        return _extract_func_name(expr.expr)
    if isinstance(expr, UnknownExpr):
        return expr.kind
    return "unknown_func"


def _extract_template_arg(s: str, prefix: str) -> str | None:
    """从 'prefix<arg>' 中提取 arg，正确处理嵌套 <>。"""
    start = s.find(prefix + "<")
    if start < 0:
        # 试 std::prefix<
        start = s.find("std::" + prefix + "<")
        if start < 0:
            return None
        start += 5  # skip "std::"
    start += len(prefix) + 1  # skip "prefix<"
    depth = 1
    i = start
    while i < len(s) and depth > 0:
        if s[i] == '<':
            depth += 1
        elif s[i] == '>':
            depth -= 1
        i += 1
    if depth == 0:
        return s[start:i - 1].strip()
    return None


def _extract_pair_args(s: str) -> tuple[str, str] | None:
    """从 'pair<A, B>' 中提取 A 和 B，正确处理嵌套 <>。"""
    inner = _extract_template_arg(s, "pair")
    if inner is None:
        return None
    # 在最外层找逗号分隔点
    depth = 0
    for i, ch in enumerate(inner):
        if ch == '<':
            depth += 1
        elif ch == '>':
            depth -= 1
        elif ch == ',' and depth == 0:
            return (inner[:i].strip(), inner[i + 1:].strip())
    return None


def parse_translation_unit(json_str: str) -> list[FuncIR]:
    ast = json.loads(json_str)
    funcs = []
    for node in ast.get("inner", []):
        if node.get("kind") != "FunctionDecl":
            continue
        # 跳过无函数体的声明
        if not any(c.get("kind") == "CompoundStmt" for c in node.get("inner", [])):
            continue
        # 跳过隐式/系统声明
        if node.get("isImplicit"):
            continue
        loc = node.get("loc", {})
        if loc.get("includedFrom"):
            continue
        funcs.append(parse_func_decl(node))
    return funcs


def parse_func_decl(node: dict) -> FuncIR:
    name = node.get("name", "?")
    # 解析返回类型：从 type.qualType 中提取 "ret_type(params)" 的 ret_type 部分
    qual_type = node.get("type", {}).get("qualType", "")
    ret_str = qual_type.split("(")[0].strip()
    ret_type = parse_type(ret_str)

    params = []
    body_stmts = []
    source_line = node.get("loc", {}).get("line", 0)
    source_file = node.get("loc", {}).get("file", "")

    for child in node.get("inner", []):
        kind = child.get("kind")
        if kind == "ParmVarDecl":
            params.append(parse_param(child))
        elif kind == "CompoundStmt":
            body_stmts = parse_compound(child)

    return FuncIR(
        name=name, params=params, ret_type=ret_type,
        body=body_stmts, source_file=source_file, source_line=source_line
    )


def parse_param(node: dict) -> ParamIR:
    name = node.get("name", "_")
    qual = _get_qual_type(node)
    is_output = "&" in qual and "const" not in qual
    is_const_ref = "const" in qual and "&" in qual
    # 去掉引用和 const 标记后解析基础类型
    base = qual.replace("const", "").replace("&", "").strip()
    return ParamIR(
        name=name, typ=parse_type(base),
        is_output=is_output, is_const_ref=is_const_ref
    )


TYPE_MAP = {
    "uint64_t": BaseType.UINT64,
    "unsigned long": BaseType.UINT64,
    "unsigned long long": BaseType.UINT64,
    "int64_t": BaseType.INT64,
    "long": BaseType.INT64,
    "long long": BaseType.INT64,
    "int": BaseType.INT64,
    "uint32_t": BaseType.UINT32,
    "unsigned int": BaseType.UINT32,
    "unsigned __int128": BaseType.UINT128,
    "__int128 unsigned": BaseType.UINT128,
    "bool": BaseType.BOOL,
    "_Bool": BaseType.BOOL,
    "void": BaseType.VOID,
    "size_t": BaseType.UINT64,
    "size_type": BaseType.UINT64,
    "std::size_t": BaseType.UINT64,
    "auto": "auto",
    "const auto &": "auto",
    "auto &": "auto",
    "double": BaseType.INT64,  # 近似：float → Int（用于 log/sqrt 估计）
}

# CLPoly 特有类型映射
# 注意：子串匹配按顺序执行，长的放前面（"upolynomial_<ZZ>" 在 "ZZ" 之前）
CLPOLY_TYPE_MAP = {
    # 多项式类型（长匹配优先）
    "upolynomial_<Zp>": StructType("SparsePolyZp", []),
    "upolynomial_<ZZ>": StructType("SparsePolyZZ", []),
    "polynomial_<ZZ": StructType("MvPolyZZ", []),     # polynomial_<ZZ, lex_<...>>
    "polynomial_<Zp": StructType("MvPolyZp", []),     # polynomial_<Zp, lex_<...>>
    # Clang 展开后的 basic_polynomial 形式
    "basic_polynomial": StructType("MvPolyZZ", []),   # fallback（下面会精确匹配系数）
    "SparsePolyZp": StructType("SparsePolyZp", []),
    "SparsePolyZZ": StructType("SparsePolyZZ", []),
    # 模板别名（Clang 可能输出 PolyZp/PolyZZ 而非完整类型）
    "PolyZp": StructType("MvPolyZp", []),
    "PolyZZ": StructType("MvPolyZZ", []),
    # 单项式
    "umonomial": StructType("UMonomial", []),
    "basic_monomial": StructType("MvMonomial", []),    # 多变量单项式
    # 数值类型
    "Zp": StructType("Zp", []),
    "QQ": StructType("QQ", []),
    "ZZ": StructType("ZZ", []),                        # CLPoly 大整数 → Lean Int（通过 abbrev）
    "variable": StructType("Variable", []),
    # factorization 结构
    "factorization": StructType("Factorization", []),
    # 内部结构
    "__hensel_node": StructType("HenselNode", []),
    "__prime_selection_result": StructType("PrimeSelectionResult", []),
    "__wang_lc_result": StructType("WangLcResult", []),
    "LLLMatrix": StructType("LLLMatrix", []),
    # RNG
    "std::mt19937": BaseType.UINT64,
    "std::uniform_int_distribution": BaseType.UINT64,
    "result_type": BaseType.UINT64,
    # 迭代器 → 忽略（翻译为 sorry 或 identity）
    "iterator": BaseType.UINT64,
    "const_iterator": BaseType.UINT64,
}

# 注：CLASS_MAP (class_map.py) 统一管理所有类方法映射，此处不再重复


def parse_type(qual: str) -> TypeIR:
    qual = qual.strip()
    clean = qual.replace("const", "").replace("&", "").replace("*", "").strip()

    # 1. 基本类型（精确匹配）
    if clean in TYPE_MAP:
        return TYPE_MAP[clean]

    # 2. CLPoly 叶子类型（不含 < 的，如 "Zp", "ZZ"）——按 key 长度降序匹配（长的优先）
    if "<" not in clean:
        for key, val in sorted(CLPOLY_TYPE_MAP.items(), key=lambda x: -len(x[0])):
            if key in clean:
                return val

    # 3. CLPoly 模板类型（含 < 的，如 "upolynomial_<ZZ>"）——子串匹配
    #    在 vector/pair 拆解之前，避免被拆成 vector<pair<UMonomial, ZZ>>
    for key, val in CLPOLY_TYPE_MAP.items():
        if "<" in key and key in clean:
            return val

    # 4. 复合类型：vector 在 pair 之前（vector<pair<...>> 外层是 vector）
    inner_str = _extract_template_arg(clean, "vector")
    if inner_str:
        return ArrayType(parse_type(inner_str))
    pair_args = _extract_pair_args(clean)
    if pair_args:
        return PairType(parse_type(pair_args[0]), parse_type(pair_args[1]))

    # 5. Lambda 类型 → auto（Lean 会推断）
    if "lambda" in clean:
        return "auto"

    # 6. std 容器和工具类型
    if "mersenne_twister" in clean:
        return BaseType.UINT64
    if "std::map" in clean or "map<" in clean:
        return StructType("StdMap", [])
    if "_Self" in clean:
        return BaseType.UINT64  # STL 内部 typedef，占位

    # 6. basic_polynomial 精细匹配（Clang 展开后的 polynomial_ 别名）
    if "basic_polynomial" in clean:
        if "Zp" in clean:
            return StructType("MvPolyZp", [])
        return StructType("MvPolyZZ", [])

    # 6. CLPoly 叶子类型 fallback（含 < 但不在模板映射中的）
    for key, val in sorted(CLPOLY_TYPE_MAP.items(), key=lambda x: -len(x[0])):
        if key in clean:
            return val

    # reference (*)(size_type) → 未知但非致命类型
    if "reference" in clean or "pointer" in clean:
        return "auto"
    # lex_<...> * 指针 → 忽略（单项式序，翻译中不用）
    if "lex_<" in clean or "grlex" in clean:
        return BaseType.UINT64  # 占位
    # Poly / PolyZp 等模板别名
    if clean in ("Poly", "PolyZp", "PolyZZ"):
        return StructType("MvPolyZZ", [])
    # BMono 等模板别名
    if "BMono" in clean or "monomial_type" in clean:
        return StructType("MvMonomial", [])

    return qual  # 未知类型保留原始字符串


def parse_compound(node: dict) -> list[StmtIR]:
    return [parse_stmt(c) for c in node.get("inner", [])]


def parse_stmt(node: dict) -> StmtIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    # 1. 透传
    if kind in TRANSPARENT_NODES:
        if inner:
            return parse_stmt(inner[0])
        return ExprStmt(Lit(0))

    # 2. 元信息
    if kind in META_NODES:
        return ExprStmt(Lit(0))

    # 3. 查语句注册表
    handler = STMT_HANDLERS.get(kind)
    if handler:
        return handler(node)

    # 4. DeclStmt：内部分发 VarDecl / DecompositionDecl
    if kind == "DeclStmt":
        stmts = []
        for decl in inner:
            if decl.get("kind") == "VarDecl":
                stmts.append(parse_var_decl(decl))
            elif decl.get("kind") == "DecompositionDecl":
                stmts.extend(_handle_decomposition_decl(decl))
        return stmts[0] if len(stmts) == 1 else UnknownStmt("MultiDecl", stmts) if stmts else ExprStmt(Lit(0))

    if kind == "CompoundStmt":
        # 扁平化：不嵌套 compound
        return UnknownStmt("CompoundStmt", [parse_stmt(c) for c in inner])

    if kind == "IfStmt":
        return parse_if(node)

    if kind == "ForStmt":
        return parse_for(node)

    if kind == "WhileStmt":
        return parse_while(node)

    if kind == "DoStmt":
        # do { body } while (cond) → 翻译为 while（首次无条件执行由 SSA 处理）
        return parse_while(node)

    if kind == "ReturnStmt":
        if inner:
            return ReturnStmt(parse_expr(inner[0]))
        return ReturnStmt()

    if kind == "BreakStmt":
        return UnknownStmt("BreakStmt")

    if kind == "ContinueStmt":
        return UnknownStmt("ContinueStmt")

    if kind == "CXXForRangeStmt":
        # range-for: for (auto& x : collection) { body }
        # inner 通常有 [range_decl, begin_expr, end_expr, cond, step, loop_var, body]
        # 简化：提取 collection 和 body，翻译为 Array.foldl
        return parse_range_for(node)

    # assert(cond) 在 Clang AST 中通常是 CallExpr 到 __assert_fail 或类似
    # 或者是一个包含 __builtin_expect 的表达式
    # 简化处理：检测 assert 模式
    if kind == "BinaryOperator" and node.get("opcode") == "=":
        lhs_expr = parse_expr(inner[0]) if inner else UnknownExpr("empty")
        rhs_expr = parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("empty")
        if isinstance(lhs_expr, Var):
            return AssignStmt(lhs_expr, rhs_expr)
        return ExprStmt(BinOp("=", lhs_expr, rhs_expr))

    if kind == "CompoundAssignOperator":
        opcode = node.get("opcode", "?=")
        base_op = opcode.replace("=", "")
        if len(inner) >= 2:
            lhs = parse_expr(inner[0])
            rhs = parse_expr(inner[1])
            if isinstance(lhs, Var):
                return AssignStmt(lhs, BinOp(base_op, lhs, rhs))
        return ExprStmt(UnknownExpr(kind))

    if kind == "CXXThrowExpr":
        return Throw("runtime_error")

    if kind == "ExprWithCleanups":
        # Clang 包裹 throw 等需要清理的表达式
        if inner:
            return parse_stmt(inner[0])

    if kind == "CXXOperatorCallExpr":
        inner = node.get("inner", [])
        if len(inner) >= 3:
            callee_name = _extract_callee_name(inner[0])
            # 区分 operator= vs operator*=/+=/etc（复合赋值）
            is_assign = callee_name == "operator="
            compound_op = None
            if not is_assign and "operator" in callee_name and "=" in callee_name:
                # operator*= → *, operator+= → +, operator-= → -, etc
                op_part = callee_name.replace("operator", "").replace("=", "")
                if op_part in ("+", "-", "*", "/", "%", "<<", ">>", "&", "|", "^"):
                    compound_op = op_part
                    is_assign = True
            if is_assign:
                # inner[1] = LHS, inner[2] = RHS
                rhs_node = inner[2]
                while rhs_node.get("kind") in ("ImplicitCastExpr", "CStyleCastExpr", "MaterializeTemporaryExpr", "CXXBindTemporaryExpr"):
                    sub = rhs_node.get("inner", [])
                    if sub:
                        rhs_node = sub[0]
                    else:
                        break
                lhs = parse_expr(inner[1])
                rhs = parse_expr(rhs_node)
                # 复合赋值 a *= b → a = a * b
                if compound_op is not None:
                    rhs = BinOp(compound_op, lhs, rhs)
                if isinstance(lhs, Var):
                    return AssignStmt(lhs, rhs)
                # LHS 是 field access → 保留字段信息给 SSA transform_member_assign
                if isinstance(lhs, FieldAccess):
                    return ExprStmt(BinOp("=", lhs, rhs))
                if isinstance(lhs, Cast):
                    lhs_var = _extract_root_var(lhs)
                    if lhs_var:
                        return AssignStmt(lhs_var, rhs)
        return ExprStmt(parse_expr(node))

    # 表达式语句（含语句位置的 ParenExpr、CStyleCastExpr 等）
    if kind in ("CallExpr", "CXXMemberCallExpr", "UnaryOperator",
                "BinaryOperator", "ParenExpr", "CStyleCastExpr",
                "CXXFunctionalCastExpr", "CXXBindTemporaryExpr"):
        return ExprStmt(parse_expr(node))

    # assert 宏展开检测：(cond ? void() : __assert_fail(...))
    if kind == "ParenExpr" and inner:
        child = inner[0]
        if child.get("kind") == "ConditionalOperator":
            cond_inner = child.get("inner", [])
            # __assert_fail 在 false 分支（第三个子节点）
            if len(cond_inner) >= 3:
                false_branch = cond_inner[2]
                if _is_assert_fail(false_branch):
                    cond_expr = parse_expr(cond_inner[0])
                    return Require(cond_expr, _gen_require_name(cond_expr), "assert")

    # 5. 查表达式注册表（表达式出现在语句位置 → ExprStmt）
    expr_handler = EXPR_HANDLERS.get(kind)
    if expr_handler:
        return ExprStmt(expr_handler(node))

    return UnknownStmt(kind, [parse_stmt(c) for c in inner] if inner else [])


def parse_var_decl(node: dict) -> LetStmt:
    name = node.get("name", "_")
    qual = _get_qual_type(node)
    typ = parse_type(qual)
    init_inner = node.get("inner", [])
    value = parse_expr(init_inner[0]) if init_inner else Lit(0)
    return LetStmt(Var(name), typ, value)


def parse_if(node: dict) -> IfStmt:
    inner = node.get("inner", [])
    # Clang IfStmt inner: [cond, then, else?]
    # 可能有额外的 null 节点
    parts = [c for c in inner if c is not None]
    cond = parse_expr(parts[0]) if parts else UnknownExpr("empty_cond")
    then_body = _wrap_body(parse_stmt(parts[1])) if len(parts) > 1 else []
    else_body = _wrap_body(parse_stmt(parts[2])) if len(parts) > 2 else []
    return IfStmt(cond, then_body, else_body)


def parse_for(node: dict) -> StmtIR:
    """ForStmt → UnknownStmt("ForLoop", [init, cond, step, body])。

    Clang ForStmt inner 固定 5 个子节点：[init, condVar, cond, step, body]
    condVar 通常是 NULL（kind=""）。
    任何子节点可能是 NULL。
    """
    inner = node.get("inner", [])

    def get_child(idx: int) -> dict | None:
        if idx < len(inner):
            child = inner[idx]
            if child.get("kind", ""):
                return child
        return None

    # 固定位置：[0]=init, [1]=condVar(null), [2]=cond, [3]=step, [4]=body
    init_node = get_child(0)
    # [1] = condVar, 通常 null，跳过
    cond_node = get_child(2)
    step_node = get_child(3)
    body_node = get_child(4)

    return UnknownStmt("ForLoop", [
        parse_stmt(init_node) if init_node else UnknownStmt("null"),
        ExprStmt(parse_expr(cond_node)) if cond_node else ExprStmt(Lit(1, BaseType.BOOL)),
        ExprStmt(parse_expr(step_node)) if step_node else UnknownStmt("null"),
        parse_stmt(body_node) if body_node else UnknownStmt("null"),
    ])


def parse_range_for(node: dict) -> StmtIR:
    """CXXForRangeStmt → ForLoop over array indices。

    简化翻译：for (auto& x : arr) { body }
    → for (int __i = 0; __i < arr.size; __i++) { auto& x = arr[__i]; body }
    """
    inner = node.get("inner", [])
    # 找循环变量和 body
    loop_var_name = "elem"
    loop_var_type = "auto"
    collection = UnknownExpr("collection")
    body_stmt = UnknownStmt("empty")

    for child in inner:
        k = child.get("kind", "")
        if k == "DeclStmt":
            for decl in child.get("inner", []):
                if decl.get("kind") == "VarDecl":
                    vname = decl.get("name", "")
                    if vname.startswith("__range"):
                        # __range1 = collection
                        range_inner = decl.get("inner", [])
                        if range_inner:
                            collection = parse_expr(range_inner[0])
                    elif not vname.startswith("__"):
                        # 用户的循环变量（如 term）
                        loop_var_name = vname
                        loop_var_type = parse_type(decl.get("type", {}).get("qualType", "auto"))
        elif k == "CompoundStmt":
            body_stmt = parse_stmt(child)
        elif k in ("ExprWithCleanups", "CXXOperatorCallExpr", "BinaryOperator",
                   "CXXMemberCallExpr", "CallExpr", "UnaryOperator"):
            # 单语句 range-for body（无 CompoundStmt 包裹）
            body_stmt = parse_stmt(child)

    # 翻译为 ForLoop
    idx_var = Var("__idx")
    return UnknownStmt("RangeForLoop", [
        ExprStmt(collection),                    # 被遍历的集合
        LetStmt(Var(loop_var_name), loop_var_type, Lit(0)),  # 循环变量名+类型
        body_stmt,                               # 循环体
    ])


def parse_while(node: dict) -> StmtIR:
    inner = node.get("inner", [])
    cond = parse_expr(inner[0]) if inner else UnknownExpr("empty")
    body = parse_stmt(inner[1]) if len(inner) > 1 else UnknownStmt("empty")
    return UnknownStmt("WhileLoop", [ExprStmt(cond), body])


# ============================================================
# 注册表
# ============================================================

# 透传节点：不改变语义的包装层，递归第一个子节点
TRANSPARENT_NODES = {
    "ExprWithCleanups", "MaterializeTemporaryExpr", "CXXBindTemporaryExpr",
    "SubstNonTypeTemplateParmExpr",
    # 注：ParenExpr 不在此列——它在 parse_stmt 中可能包裹 assert 宏，需要特殊处理
}

# 元信息节点：调试信息，对算法翻译无语义价值
META_NODES = {
    "SourceLocExpr", "PredefinedExpr", "TypeAliasDecl",
    "FullComment", "GCCAsmStmt", "NullStmt",
}

# 表达式处理器（新增/覆盖用，缺的走 _parse_expr_inner fallback）
EXPR_HANDLERS: dict[str, callable] = {}

# 语句处理器（新增/覆盖用，缺的走 legacy parse_stmt fallback）
STMT_HANDLERS: dict[str, callable] = {}


def _handle_operator_subscript(node: dict) -> ExprIR:
    """operator[] → ArrayAccess"""
    inner = node.get("inner", [])
    if len(inner) >= 3:
        return ArrayAccess(parse_expr(inner[1]), parse_expr(inner[2]))
    return UnknownExpr("operator[]_bad_args")


def _handle_decomposition_decl(node: dict):
    """auto [a, b] = expr → let pair := expr; let a := pair.1; let b := pair.2"""
    inner = node.get("inner", [])
    init_expr = None
    binding_names = []
    for child in inner:
        if child.get("kind") == "BindingDecl":
            binding_names.append(child.get("name", "_"))
        elif init_expr is None and child.get("kind") not in ("BindingDecl",):
            init_expr = parse_expr(child)
    if init_expr is None:
        return [UnknownStmt("DecompositionDecl_no_init")]
    result = []
    pair_var = Var("_decomp")
    result.append(LetStmt(pair_var, "auto", init_expr))
    for i, name in enumerate(binding_names):
        if len(binding_names) == 1:
            accessor = pair_var
        elif i == 0:
            accessor = FieldAccess(pair_var, "1")
        elif i == len(binding_names) - 1:
            accessor = pair_var
            for _ in range(i):
                accessor = FieldAccess(accessor, "2")
        else:
            accessor = pair_var
            for _ in range(i):
                accessor = FieldAccess(accessor, "2")
            accessor = FieldAccess(accessor, "1")
        result.append(LetStmt(Var(name), "auto", accessor))
    return result


def parse_expr(node: dict) -> ExprIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    # 1. 透传
    if kind in TRANSPARENT_NODES:
        if inner:
            return parse_expr(inner[0])
        return Lit(0)

    # 2. 元信息
    if kind in META_NODES:
        return Lit(0)

    # 3. 查注册表（新增/覆盖的 handler）
    handler = EXPR_HANDLERS.get(kind)
    if handler:
        result = handler(node)
        ast_qual = _get_qual_type(node)
        if ast_qual:
            result._ast_type = parse_type(ast_qual)
        return result

    # 4. legacy fallback（现有 _parse_expr_inner）
    result = _parse_expr_inner(node)
    ast_qual = _get_qual_type(node)
    if ast_qual:
        result._ast_type = parse_type(ast_qual)
    return result


def _parse_expr_inner(node: dict) -> ExprIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    if kind == "IntegerLiteral":
        return Lit(int(node.get("value", "0")))

    if kind == "CXXBoolLiteralExpr":
        return Lit(1 if node.get("value", False) else 0, BaseType.BOOL)

    if kind == "StringLiteral":
        val = node.get("value", "")
        return Lit(0)  # 字符串字面量在数学算法中不影响正确性

    if kind == "ImplicitValueInitExpr":
        # 零值初始化（如 ZZ() → 0, vector() → #[]）
        return Lit(0)

    if kind == "FloatingLiteral":
        # 浮点数 → 取整（用于 log/ceil 估计，精确值不重要）
        val = node.get("value", "0")
        try:
            return Lit(int(float(val)))
        except (ValueError, OverflowError):
            return Lit(0)

    if kind == "DeclRefExpr":
        ref = node.get("referencedDecl", {})
        return Var(ref.get("name", "?"))

    if kind == "BinaryOperator":
        op = node.get("opcode", "?")
        if len(inner) >= 2:
            return BinOp(op, parse_expr(inner[0]), parse_expr(inner[1]))
        return UnknownExpr(f"BinOp_{op}")

    if kind == "UnaryOperator":
        op = node.get("opcode", "?")
        if inner:
            operand = parse_expr(inner[0])
            # 后缀 ++/-- 也用 UnaryOp 表示
            return UnaryOp(op, operand)
        return UnknownExpr(f"UnaryOp_{op}")

    if kind == "ConditionalOperator":
        if len(inner) >= 3:
            return CondExpr(parse_expr(inner[0]), parse_expr(inner[1]), parse_expr(inner[2]))

    # CallExpr 在后面统一处理（含 std::move, make_pair 检测）

    if kind == "ArraySubscriptExpr":
        if len(inner) >= 2:
            return ArrayAccess(parse_expr(inner[0]), parse_expr(inner[1]))

    if kind == "MemberExpr":
        fname = node.get("name", "?")
        from class_map import FIELD_MAP
        lean_fname = FIELD_MAP.get(fname, fname)
        if inner:
            return FieldAccess(parse_expr(inner[0]), lean_fname)

    if kind in ("CStyleCastExpr", "ImplicitCastExpr"):
        if inner:
            target_qual = node.get("type", {}).get("qualType", "")
            target_type = parse_type(target_qual)
            # 直接取子节点（不追踪链），保留每一层 Cast
            source_qual = inner[0].get("type", {}).get("qualType", "")
            source_type = parse_type(source_qual)
            child_expr = parse_expr(inner[0])
            if target_type != source_type:
                # lambda → const lambda 等无害转换：直接透传
                if isinstance(source_type, str) and isinstance(target_type, str):
                    if "lambda" in source_type or "lambda" in target_type:
                        return child_expr
                return Cast(child_expr, target_type, source_type)
            return child_expr

    if kind == "ParenExpr":
        if inner:
            return parse_expr(inner[0])

    if kind == "CXXMemberCallExpr":
        if inner:
            member = inner[0]
            method_name = member.get("name", "?")
            # 提取 this 对象：跳过 ImplicitCastExpr 链
            obj_node = member.get("inner", [{}])[0] if member.get("inner") else {}
            while obj_node.get("kind") in ("ImplicitCastExpr", "CStyleCastExpr"):
                sub = obj_node.get("inner", [])
                if sub:
                    obj_node = sub[0]
                else:
                    break
            obj = parse_expr(obj_node) if obj_node else UnknownExpr("no_obj")
            # 查 CLASS_MAP
            from class_map import CLASS_MAP, FIELD_MAP
            obj_type = _get_class_name(obj_node)
            class_info = CLASS_MAP.get(obj_type)
            if class_info and method_name in class_info.get("methods", {}):
                category, lean_name = class_info["methods"][method_name]
                if category == "field":
                    return FieldAccess(obj, lean_name)
                elif category == "method":
                    return Call(lean_name, [obj])
                elif category == "mutate":
                    return Call("_mutate_" + lean_name, [obj])
                elif category == "mutate_push":
                    elem = parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("no_elem")
                    return ArrayPush(obj, elem)
                elif category == "noop":
                    return obj
                elif category == "identity":
                    return obj
            # operator bool → obj != 0
            if method_name == "operator bool":
                return BinOp("!=", obj, Lit(0))
            # 常见容器方法 fallback（类型未知但方法名明确）
            _COMMON_METHODS = {
                "size": ("method", "Array.size_u64"),
                "empty": ("method", "Array.isEmpty"),
                "front": ("method", "Array.front!"),
                "back": ("method", "Array.back!"),
                "push_back": ("mutate_push", "Array.push"),
                "normalization": ("noop", None),
                "reserve": ("noop", None),
                "data": ("identity", None),
                "prime": ("field", "prime"),
                "comp_ptr": ("noop", None),
            }
            if method_name in _COMMON_METHODS:
                cat, ln = _COMMON_METHODS[method_name]
                if cat == "field": return FieldAccess(obj, ln)
                if cat == "method": return Call(ln, [obj])
                if cat == "mutate": return Call("_mutate_" + ln, [obj])
                if cat == "mutate_push":
                    elem = parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("no_elem")
                    return ArrayPush(obj, elem)
                if cat in ("noop", "identity"): return obj
            # 不在 CLASS_MAP → sorry
            args = [parse_expr(a) for a in inner[1:]]
            return UnknownExpr(f"unmapped:{obj_type}.{method_name}")

    if kind in ("CXXConstructExpr", "CXXTemporaryObjectExpr"):
        qual_type = node.get("type", {}).get("qualType", "")
        typ = parse_type(qual_type)
        from class_map import CLASS_MAP

        # 检测 copy/move 构造：1 个参数且参数类型 == 构造类型 → identity
        if len(inner) == 1:
            arg_qual = inner[0].get("type", {}).get("qualType", "")
            arg_typ = parse_type(arg_qual)
            if _types_equal(typ, arg_typ):
                return parse_expr(inner[0])

        # 查 CLASS_MAP 构造函数
        class_name = typ.name if isinstance(typ, StructType) else None
        class_info = CLASS_MAP.get(class_name) if class_name else None
        if class_info:
            arg_types = tuple(parse_type(a.get("type", {}).get("qualType", ""))
                              for a in inner if a.get("type"))
            ctors = class_info.get("constructors", {})
            # 按参数类型匹配（优先精确匹配，再按数量匹配）
            ctor = _match_ctor(ctors, arg_types)
            if ctor is None and () in ctors and len(inner) == 0:
                ctor = ctors[()]
            if ctor == "default":
                return Call("default", [])  # Lean 的 default 实例
            if ctor:
                # 字面量构造函数（如 "((0 : Int))"）→ 直接作为表达式
                if ctor.startswith("(") and len(inner) == 0:
                    return Var(ctor)  # codegen 会原样输出
                args = [parse_expr(a) for a in inner]
                return Call(ctor, args)
        # pair 构造
        if isinstance(typ, PairType) and len(inner) == 2:
            return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
        # ArrayType 空构造 → #[]（通用空数组）
        if isinstance(typ, ArrayType) and len(inner) == 0:
            return Call("#[]", [])
        # BaseType 构造：不可能是 Prod.mk（非 PairType）
        if isinstance(typ, BaseType):
            if inner:
                return parse_expr(inner[-1])  # 取最后一个参数
            return Lit(0)
        # fallback：PairType 已在上方处理，其他未知类型尝试 Prod.mk
        if inner and len(inner) == 2:
            return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
        if inner and len(inner) == 1:
            return parse_expr(inner[0])
        return Lit(0)

    if kind == "CXXStdInitializerListExpr":
        if inner:
            return Call("Array.mk", [parse_expr(c) for c in inner])

    if kind == "InitListExpr":
        # {a, b} → Prod.mk a b (pair) 或 Array.mk [a, b] (array)
        if len(inner) == 2:
            return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
        return Call("Array.mk", [parse_expr(c) for c in inner])

    if kind == "CallExpr":
        if inner:
            # 检测模板中的成员调用：CallExpr → CXXDependentScopeMemberExpr
            if inner[0].get("kind") == "CXXDependentScopeMemberExpr":
                member_name = inner[0].get("member", "?")
                obj_inner = inner[0].get("inner", [])
                obj = parse_expr(obj_inner[0]) if obj_inner else UnknownExpr("no_obj")
                obj_type = _get_class_name(obj_inner[0]) if obj_inner else None
                extra_args = [parse_expr(a) for a in inner[1:]]
                # 查 CLASS_MAP
                from class_map import CLASS_MAP
                class_info = CLASS_MAP.get(obj_type) if obj_type else None
                if class_info and member_name in class_info.get("methods", {}):
                    category, lean_name = class_info["methods"][member_name]
                    if category == "field":
                        return FieldAccess(obj, lean_name)
                    elif category == "method":
                        return Call(lean_name, [obj] + extra_args)
                    elif category == "mutate":
                        return Call("_mutate_" + lean_name, [obj] + extra_args)
                    elif category == "mutate_push":
                        elem = extra_args[0] if extra_args else UnknownExpr("no_elem")
                        return ArrayPush(obj, elem)
                    elif category in ("noop", "identity"):
                        return obj
                # fallback
                return Call(f"{obj_type or 'Unknown'}.{member_name}", [obj] + extra_args)

            func_expr = parse_expr(inner[0])
            func_name = _extract_func_name(func_expr)
            all_args = [parse_expr(a) for a in inner[1:]]

            # 1. FUNC_MAP 特殊规则（identity, make_pair）
            from class_map import FUNC_MAP
            if func_name in FUNC_MAP:
                entry = FUNC_MAP[func_name]
                lean_name = entry[0]
                rule = entry[1]
                if rule == "identity":
                    return parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("identity_empty")
                elif rule == "make_pair":
                    return Call("Prod.mk", all_args)
                elif rule == "direct":
                    func_name = lean_name  # 替换为 Lean 名，后续走统一路径

            # 2. 运算符解析（G2）：func_name 以 operator 开头 → UnaryOp / BinOp
            if func_name.startswith("operator"):
                op = func_name.replace("operator", "").strip()
                if len(all_args) == 1:
                    return UnaryOp(op, all_args[0])
                elif len(all_args) == 2:
                    return BinOp(op, all_args[0], all_args[1])

            # 3. 自动检测输出参数（从被调函数签名，对所有 CallExpr 通用）
            #    设计 §3：所有参数都传入，返回值 = Prod(输出参数最终值)
            callee_sig = _get_callee_sig(inner[0])
            output_indices = _parse_output_params_from_sig(callee_sig) if callee_sig else []
            if output_indices:
                call = Call(func_name, all_args)  # 全部参数传入
                out_vars = [all_args[i] for i in output_indices if i < len(all_args)]
                if len(out_vars) == 1:
                    # 单输出 → BinOp("=", out_var, call)
                    return BinOp("=", out_vars[0], call)
                else:
                    # 多输出 → BinOp("=", Prod(out_vars), call)
                    return BinOp("=", Call("Prod.mk", out_vars), call)

            # 3. 无输出参数 → 普通调用
            return Call(func_name, all_args)

    if kind == "CXXOperatorCallExpr":
        # 运算符调用（赋值、算术、比较、函数调用运算符）
        if inner:
            ref_name = _extract_callee_name(inner[0])
            from class_map import OPERATOR_MAP, CALL_OPERATOR_MAP
            # 0a. 下标运算符 operator[] → ArrayAccess
            if "operator[]" in ref_name and len(inner) >= 3:
                return ArrayAccess(parse_expr(inner[1]), parse_expr(inner[2]))
            # 0b. 赋值运算符 operator= → BinOp("=", lhs, rhs)（SSA 变换会处理）
            if ref_name == "operator=" and len(inner) >= 3:
                return BinOp("=", parse_expr(inner[1]), parse_expr(inner[2]))
            # 1a. 一元运算符（callee + 1 operand = 2 inner）
            if len(inner) == 2 and "operator" in ref_name:
                op = ref_name.replace("operator", "").strip()
                if op in ("->", "*"):
                    # 解引用/成员访问 → identity（纯函数式无指针）
                    return parse_expr(inner[1])
                if op in ("++", "--"):
                    # 自增/自减 → identity（SSA 中已处理）
                    return parse_expr(inner[1])
                return UnaryOp(op, parse_expr(inner[1]))
            # 1b. 二元运算符：按长度降序匹配（"operator==" 优先于 "operator="）
            for op_key, (min_args, lean_op) in sorted(
                    OPERATOR_MAP.items(), key=lambda x: -len(x[0])):
                if op_key in ref_name and len(inner) >= min_args:
                    return BinOp(lean_op, parse_expr(inner[1]), parse_expr(inner[2]))
            # 2. 函数调用运算符 operator()
            if "operator()" in ref_name and len(inner) >= 2:
                # M3：区分 RNG distribution vs lambda
                callee_type = _get_qual_type(inner[1]) if len(inner) > 1 else ""
                if "uniform_int_distribution" in callee_type or "mersenne_twister" in callee_type:
                    # RNG 路径：dist(rng) → Rng.next rng upper
                    args = [parse_expr(a) for a in inner[1:]]
                    if len(args) == 2:
                        args = [args[1], args[0]]  # reverse: Rng.next rng upper
                    return Call("Rng.next", args)
                else:
                    # Lambda / 一般函数对象：直接调用 callee
                    callee = parse_expr(inner[1])
                    args = [parse_expr(a) for a in inner[2:]]
                    return Call(callee if isinstance(callee, str) else
                               (callee.name if isinstance(callee, Var) else "unknown_callee"),
                               args)
            # 3. Fallback：ref_name 提取失败时，扫描子节点找运算符名
            if len(inner) >= 3:
                for c in inner:
                    name = c.get("name", "")
                    if name.startswith("operator"):
                        op = name.replace("operator", "").strip()
                        return BinOp(op, parse_expr(inner[1]), parse_expr(inner[2]))
                # 完全无法识别的二元运算 → 用 inner[0] 的类型推断
                return BinOp("?", parse_expr(inner[1]), parse_expr(inner[2]))
            if len(inner) == 2:
                for c in inner:
                    name = c.get("name", "")
                    if name.startswith("operator"):
                        op = name.replace("operator", "").strip()
                        if op not in ("->", "++", "--"):
                            return UnaryOp(op, parse_expr(inner[1]))

    if kind == "CompoundAssignOperator":
        op = node.get("opcode", "?=").replace("=", "")
        if len(inner) >= 2:
            return BinOp(op, parse_expr(inner[0]), parse_expr(inner[1]))

    # CXXDependentScopeMemberExpr：模板中的成员调用（Gi.empty()、f.comp_ptr() 等）
    # 有 .member 字段（方法名）和 inner[0]（this 对象）
    if kind == "CXXDependentScopeMemberExpr":
        member_name = node.get("member", "?")
        obj = parse_expr(inner[0]) if inner else UnknownExpr("no_obj")
        # 先尝试作为字段访问
        from class_map import FIELD_MAP
        if member_name in FIELD_MAP:
            return FieldAccess(obj, FIELD_MAP[member_name])
        # 尝试查 CLASS_MAP（用 obj 的类型）
        obj_type = _get_class_name(inner[0]) if inner else None
        if obj_type:
            from class_map import CLASS_MAP
            class_info = CLASS_MAP.get(obj_type)
            if class_info and member_name in class_info.get("methods", {}):
                category, lean_name = class_info["methods"][member_name]
                if category == "field":
                    return FieldAccess(obj, lean_name)
                elif category == "method":
                    return Call(lean_name, [obj])
                elif category == "mutate":
                    return Call("_mutate_" + lean_name, [obj])
                elif category == "mutate_push":
                    elem = parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("no_elem")
                    return ArrayPush(obj, elem)
                elif category in ("noop", "identity"):
                    return obj
        # fallback：生成方法调用（codegen 的 . 规则会处理）
        return Call(f"{obj_type or 'Unknown'}.{member_name}", [obj])

    # UnresolvedLookupExpr：模板中的 ADL 函数调用（leadcoeff(f)、degree(f) 等）
    if kind == "UnresolvedLookupExpr":
        func_name = node.get("name", "?")
        # 查 FUNC_MAP
        from class_map import FUNC_MAP, TRANSLATION_SCOPE
        if func_name in FUNC_MAP:
            lean_name = FUNC_MAP[func_name][0]
            return Var(lean_name)
        # 翻译范围内的函数
        if func_name in TRANSLATION_SCOPE:
            return Var(func_name)
        # assert
        from class_map import ASSERT_FAIL_NAMES
        if func_name in ASSERT_FAIL_NAMES:
            return Lit(0)  # assert 由 UB collector 处理
        # 运算符：返回运算符名
        if func_name.startswith("operator"):
            return Var(func_name)
        return UnknownExpr(f"UnresolvedLookupExpr:{func_name}")

    # CXXUnresolvedConstructExpr：模板中的构造函数调用
    if kind == "CXXUnresolvedConstructExpr":
        qual_type = node.get("type", {}).get("qualType", "")
        typ = parse_type(qual_type)
        if inner:
            return parse_expr(inner[0])
        return Lit(0)

    # CXXDefaultArgExpr：默认参数 — 直接用默认值
    # CXXBindTemporaryExpr：临时对象绑定 — 直接透传子节点
    if kind == "CXXBindTemporaryExpr":
        if inner:
            return parse_expr(inner[0])
        return Lit(0)

    if kind == "CXXDefaultArgExpr":
        if inner:
            return parse_expr(inner[0])
        # 无子节点时，Clang 可能内联了值
        return Lit(0)  # fallback: 大多数默认参数是 0

    # ParmVarDecl 出现在表达式位置（lambda capture 等）
    if kind == "ParmVarDecl":
        return Var(node.get("name", "_param"))

    # CXXStaticCastExpr: static_cast<T>(expr) — 同 CStyleCastExpr
    if kind == "CXXStaticCastExpr":
        if inner:
            target_qual = node.get("type", {}).get("qualType", "")
            target_type = parse_type(target_qual)
            source_qual = inner[0].get("type", {}).get("qualType", "")
            source_type = parse_type(source_qual)
            child_expr = parse_expr(inner[0])
            if target_type != source_type:
                return Cast(child_expr, target_type, source_type)
            return child_expr

    # LambdaExpr: 简单 lambda → 内联，复杂 → 提取为辅助函数（设计 §4.4）
    if kind == "LambdaExpr":
        # 尝试简单 lambda（单 return）→ 内联
        for child in inner:
            if child.get("kind") == "CompoundStmt":
                body_inner = child.get("inner", [])
                if len(body_inner) == 1 and body_inner[0].get("kind") == "ReturnStmt":
                    ret_inner = body_inner[0].get("inner", [])
                    if ret_inner:
                        return parse_expr(ret_inner[0])
            # 找 CXXMethodDecl 中的简单 body
            if child.get("kind") == "CXXMethodDecl":
                for sub in child.get("inner", []):
                    if sub.get("kind") == "CompoundStmt":
                        method_body = sub.get("inner", [])
                        if len(method_body) == 1 and method_body[0].get("kind") == "ReturnStmt":
                            ret_inner = method_body[0].get("inner", [])
                            if ret_inner:
                                return parse_expr(ret_inner[0])

        # 复杂 lambda → 提取 capture 信息 + body
        captures = []
        lambda_params = []
        lambda_body = []
        for child in inner:
            if child.get("kind") == "CXXRecordDecl":
                for field in child.get("inner", []):
                    if field.get("kind") == "FieldDecl":
                        cap_name = field.get("name", "_cap")
                        captures.append(cap_name)
                    # 参数在 CXXMethodDecl(operator()) 内
                    if field.get("kind") == "CXXMethodDecl" and field.get("name") == "operator()":
                        for sub in field.get("inner", []):
                            if sub.get("kind") == "ParmVarDecl":
                                lambda_params.append(parse_param(sub))
            # body 是 LambdaExpr 的直接子节点 CompoundStmt
            if child.get("kind") == "CompoundStmt":
                lambda_body = parse_compound(child)

        # 复杂 lambda → 注册 body，SSA 层创建辅助函数
        global _LAMBDA_COUNTER
        _LAMBDA_COUNTER += 1
        lambda_id = _LAMBDA_COUNTER
        _LAMBDA_REGISTRY[lambda_id] = (captures, lambda_params, lambda_body)
        # 返回占位调用（SSA 层会替换为带 capture 的辅助函数调用）
        return Call(f"_lambda_{lambda_id}", [Var(c) for c in captures] if captures else [])

    # ParenListExpr: 括号列表 → 取最后一个（C++ 逗号表达式语义）
    if kind == "ParenListExpr":
        if inner:
            return parse_expr(inner[-1])
        return Lit(0)

    # 兜底：递归第一个子节点
    if inner:
        return parse_expr(inner[0])

    return UnknownExpr(kind)


def _is_assert_fail(node: dict) -> bool:
    """检测 assert 宏展开后的失败函数调用。查 ASSERT_FAIL_NAMES。"""
    from class_map import ASSERT_FAIL_NAMES
    if node.get("kind") == "CallExpr":
        inner = node.get("inner", [])
        if inner:
            callee = inner[0]
            ref = callee.get("referencedDecl", {})
            callee_name = ref.get("name", "")
            if any(name in callee_name for name in ASSERT_FAIL_NAMES):
                return True
            # 递归检查（可能有隐式转换包裹）
            if callee.get("inner"):
                return _is_assert_fail({"kind": "CallExpr", "inner": callee.get("inner")})
    # 递归检查子节点
    for child in node.get("inner", []):
        if _is_assert_fail(child):
            return True
    return False


def _gen_require_name(cond: ExprIR) -> str:
    """从条件表达式生成 require 参数名。"""
    if isinstance(cond, BinOp) and cond.op == "!=":
        if isinstance(cond.lhs, Var):
            return f"h_{cond.lhs.name}_ne"
        if isinstance(cond.rhs, Var):
            return f"h_{cond.rhs.name}_ne"
    if isinstance(cond, Var):
        return f"h_{cond.name}"
    return "h_assert"


def _wrap_body(stmt: StmtIR) -> list[StmtIR]:
    """将单条语句包装为列表；CompoundStmt 展开。"""
    if isinstance(stmt, UnknownStmt) and stmt.kind == "CompoundStmt":
        return stmt.children
    return [stmt]
