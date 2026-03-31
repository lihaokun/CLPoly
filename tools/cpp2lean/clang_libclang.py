"""
CLPoly C++ → Lean IR 翻译器：libclang 版 AST 解析器

替换 clang_ast.py（JSON 版）。使用 libclang Python 绑定直接遍历 AST，
无需 773MB JSON dump。

用法：
    from clang_libclang import parse_file
    funcs = parse_file("clpoly/polynomial_factorize_zp.hh",
                       args=["-std=c++17", "-I./"],
                       target_file="factorize_zp")
"""

from __future__ import annotations
from clang.cindex import Index, Cursor, CursorKind, TypeKind
from ir_types import *


# ============================================================
# 系统头文件路径（自动检测）
# ============================================================

def _get_system_includes() -> list[str]:
    import subprocess
    try:
        result = subprocess.run(
            ["g++", "-E", "-Wp,-v", "-x", "c++", "-"],
            input="", capture_output=True, text=True
        )
        paths = []
        capture = False
        for line in result.stderr.splitlines():
            if "#include <...>" in line:
                capture = True
                continue
            if line.startswith("End of search"):
                break
            if capture and line.startswith(" "):
                paths.append(f"-I{line.strip()}")
        return paths
    except Exception:
        return ["-I/usr/include", "-I/usr/lib/gcc/x86_64-linux-gnu/13/include"]


SYSTEM_INCLUDES = _get_system_includes()


# ============================================================
# 类型映射
# ============================================================

def translate_type(clang_type) -> TypeIR:
    """Clang Type → TypeIR。"""
    spelling = clang_type.spelling.replace("const", "").replace("&", "").replace("*", "").strip()

    # 基本类型
    type_map = {
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
        "bool": BaseType.BOOL,
        "void": BaseType.VOID,
    }
    if spelling in type_map:
        return type_map[spelling]

    # CLPoly 特有类型
    clpoly_map = {
        "upolynomial_<Zp>": StructType("SparsePolyZp", []),
        "umonomial": StructType("UMonomial", []),
        "Zp": StructType("Zp", []),
    }
    for key, val in clpoly_map.items():
        if key in spelling:
            return val

    # vector<T>
    if "vector" in spelling:
        # 尝试从模板参数获取元素类型
        if clang_type.get_num_template_arguments() > 0:
            elem = translate_type(clang_type.get_template_argument_type(0))
            return ArrayType(elem)
        return ArrayType(BaseType.UINT64)  # fallback

    # pair<A,B>
    if "pair" in spelling:
        if clang_type.get_num_template_arguments() >= 2:
            fst = translate_type(clang_type.get_template_argument_type(0))
            snd = translate_type(clang_type.get_template_argument_type(1))
            return PairType(fst, snd)

    return spelling  # 未知类型保留字符串


# ============================================================
# 表达式翻译
# ============================================================

def translate_expr(cursor: Cursor) -> ExprIR:
    """Clang Cursor → ExprIR。"""
    kind = cursor.kind
    children = list(cursor.get_children())

    if kind == CursorKind.INTEGER_LITERAL:
        tokens = list(cursor.get_tokens())
        val = int(tokens[0].spelling) if tokens else 0
        return Lit(val)

    if kind == CursorKind.DECL_REF_EXPR:
        return Var(cursor.spelling)

    if kind == CursorKind.BINARY_OPERATOR:
        tokens = list(cursor.get_tokens())
        # 从 token 中提取运算符
        op = _extract_binop(tokens, children)
        if len(children) >= 2:
            return BinOp(op, translate_expr(children[0]), translate_expr(children[1]))

    if kind == CursorKind.UNARY_OPERATOR:
        tokens = list(cursor.get_tokens())
        op = tokens[0].spelling if tokens else "?"
        # 后缀 ++ 的 token 在最后
        if tokens and tokens[-1].spelling in ("++", "--"):
            op = tokens[-1].spelling
        if children:
            return UnaryOp(op, translate_expr(children[0]))

    if kind == CursorKind.CONDITIONAL_OPERATOR:
        if len(children) >= 3:
            return CondExpr(translate_expr(children[0]),
                          translate_expr(children[1]),
                          translate_expr(children[2]))

    if kind == CursorKind.CALL_EXPR:
        if children:
            func_name = children[0].spelling or _get_func_name(children[0])
            # std::move → 忽略
            if "move" in func_name and len(children) >= 2:
                return translate_expr(children[1])
            # std::make_pair
            if "make_pair" in func_name and len(children) >= 3:
                return Call("Prod.mk", [translate_expr(children[1]),
                                        translate_expr(children[2])])
            args = [translate_expr(c) for c in children[1:]]
            return Call(func_name, args)

    if kind == CursorKind.CXX_MEMBER_CALL_EXPR:
        if children:
            member = children[0]
            method_name = member.spelling
            # 获取 this 对象
            member_children = list(member.get_children())
            obj = translate_expr(member_children[0]) if member_children else UnknownExpr("no_obj")
            # 无参 getter → 字段访问
            getter_map = {"empty": "isEmpty", "size": "size", "front": "front!",
                         "back": "back!", "deg": "deg", "prime": "prime", "val": "val"}
            if method_name in getter_map:
                return FieldAccess(obj, getter_map[method_name])
            # push_back → ArrayPush
            if method_name == "push_back" and len(children) >= 2:
                return ArrayPush(obj, translate_expr(children[1]))
            # reserve → no-op
            if method_name in ("reserve", "normalization"):
                return Call(method_name, [obj])
            args = [translate_expr(c) for c in children[1:]]
            return Call(method_name, [obj] + args)

    if kind == CursorKind.ARRAY_SUBSCRIPT_EXPR:
        if len(children) >= 2:
            return ArrayAccess(translate_expr(children[0]), translate_expr(children[1]))

    if kind == CursorKind.MEMBER_REF_EXPR:
        field = cursor.spelling
        # C++ pair .first/.second → Lean .fst/.snd
        field_map = {"first": "fst", "second": "snd"}
        lean_field = field_map.get(field, field)
        if children:
            return FieldAccess(translate_expr(children[0]), lean_field)

    if kind in (CursorKind.CSTYLE_CAST_EXPR, CursorKind.CXX_STATIC_CAST_EXPR,
                CursorKind.CXX_FUNCTIONAL_CAST_EXPR):
        target = translate_type(cursor.type)
        if children:
            source = translate_type(children[0].type)
            return Cast(translate_expr(children[0]), target, source)

    if kind == CursorKind.IMPLICIT_CAST_EXPR:
        target = translate_type(cursor.type)
        if children:
            source = translate_type(children[0].type)
            if target != source:
                return Cast(translate_expr(children[0]), target, source)
            return translate_expr(children[0])

    if kind == CursorKind.PAREN_EXPR:
        if children:
            return translate_expr(children[0])

    if kind == CursorKind.CXX_CONSTRUCT_EXPR:
        typ = translate_type(cursor.type)
        if isinstance(typ, StructType) and typ.name == "SparsePolyZp":
            return Call("Array.empty", [])
        if isinstance(typ, ArrayType):
            return Call("Array.empty", [])
        if children and len(children) == 1:
            return translate_expr(children[0])
        return Lit(0)

    if kind == CursorKind.MATERIALIZE_TEMPORARY_EXPR:
        if children:
            return translate_expr(children[0])

    if kind == CursorKind.CXX_BIND_TEMPORARY_EXPR:
        if children:
            return translate_expr(children[0])

    if kind == CursorKind.EXPR_WITH_CLEANUPS:
        if children:
            return translate_expr(children[0])

    # 兜底
    if children:
        return translate_expr(children[0])
    return UnknownExpr(kind.name)


# ============================================================
# 语句翻译
# ============================================================

def translate_stmt(cursor: Cursor) -> StmtIR:
    """Clang Cursor → StmtIR。"""
    kind = cursor.kind
    children = list(cursor.get_children())

    if kind == CursorKind.DECL_STMT:
        stmts = []
        for child in children:
            if child.kind == CursorKind.VAR_DECL:
                stmts.append(translate_var_decl(child))
        return stmts[0] if len(stmts) == 1 else UnknownStmt("MultiDecl", stmts)

    if kind == CursorKind.COMPOUND_STMT:
        return UnknownStmt("CompoundStmt", [translate_stmt(c) for c in children])

    if kind == CursorKind.IF_STMT:
        cond = translate_expr(children[0]) if children else UnknownExpr("empty")
        then_body = _wrap_body(translate_stmt(children[1])) if len(children) > 1 else []
        else_body = _wrap_body(translate_stmt(children[2])) if len(children) > 2 else []
        return IfStmt(cond, then_body, else_body)

    if kind == CursorKind.FOR_STMT:
        return translate_for(cursor, children)

    if kind == CursorKind.WHILE_STMT:
        cond = translate_expr(children[0]) if children else Lit(1, BaseType.BOOL)
        body = translate_stmt(children[1]) if len(children) > 1 else UnknownStmt("empty")
        return UnknownStmt("WhileLoop", [ExprStmt(cond), body])

    if kind == CursorKind.CXX_FOR_RANGE_STMT:
        return translate_range_for(cursor, children)

    if kind == CursorKind.RETURN_STMT:
        if children:
            return ReturnStmt(translate_expr(children[0]))
        return ReturnStmt()

    if kind == CursorKind.BREAK_STMT:
        return UnknownStmt("BreakStmt")

    if kind == CursorKind.CONTINUE_STMT:
        return UnknownStmt("ContinueStmt")

    if kind == CursorKind.CXX_THROW_EXPR:
        return Throw("runtime_error")

    # assert 检测（通过 token 扫描）
    tokens = list(cursor.get_tokens())
    if tokens and tokens[0].spelling == "assert":
        # assert(cond) — 提取条件
        if children:
            cond = translate_expr(children[0])
            return Require(cond, _gen_require_name(cond), "assert")

    # 赋值
    if kind == CursorKind.BINARY_OPERATOR:
        op_tokens = list(cursor.get_tokens())
        op = _extract_binop(op_tokens, children)
        if op == "=" and len(children) >= 2:
            lhs = translate_expr(children[0])
            rhs = translate_expr(children[1])
            if isinstance(lhs, Var):
                return AssignStmt(lhs, rhs)
            return ExprStmt(BinOp("=", lhs, rhs))

    if kind == CursorKind.COMPOUND_ASSIGNMENT_OPERATOR:
        tokens = list(cursor.get_tokens())
        op = _extract_compound_op(tokens)
        if len(children) >= 2:
            lhs = translate_expr(children[0])
            rhs = translate_expr(children[1])
            if isinstance(lhs, Var):
                return AssignStmt(lhs, BinOp(op, lhs, rhs))
        return ExprStmt(UnknownExpr(kind.name))

    if kind == CursorKind.CXX_OPERATOR_CALL_EXPR:
        # operator= (move/copy)
        if len(children) >= 3:
            lhs = translate_expr(children[1])
            rhs = translate_expr(children[2])
            if isinstance(lhs, Var):
                return AssignStmt(lhs, rhs)

    # 表达式语句
    if kind in (CursorKind.CALL_EXPR, CursorKind.CXX_MEMBER_CALL_EXPR,
                CursorKind.UNARY_OPERATOR, CursorKind.EXPR_WITH_CLEANUPS):
        return ExprStmt(translate_expr(cursor))

    return UnknownStmt(kind.name, [translate_stmt(c) for c in children] if children else [])


def translate_var_decl(cursor: Cursor) -> LetStmt:
    name = cursor.spelling
    typ = translate_type(cursor.type)
    children = list(cursor.get_children())
    value = translate_expr(children[0]) if children else Lit(0)
    return LetStmt(Var(name), typ, value)


def translate_for(cursor: Cursor, children: list) -> StmtIR:
    """ForStmt → UnknownStmt("ForLoop", [init, cond, step, body])。"""
    # Clang ForStmt children: [init, cond, step, body] (some may be None/missing)
    init = translate_stmt(children[0]) if len(children) > 0 and children[0].kind != CursorKind.NULL_STMT else UnknownStmt("null")
    cond = ExprStmt(translate_expr(children[1])) if len(children) > 1 and children[1].kind != CursorKind.NULL_STMT else ExprStmt(Lit(1, BaseType.BOOL))
    step = ExprStmt(translate_expr(children[2])) if len(children) > 2 and children[2].kind != CursorKind.NULL_STMT else UnknownStmt("null")
    body = translate_stmt(children[3]) if len(children) > 3 else UnknownStmt("empty")
    return UnknownStmt("ForLoop", [init, cond, step, body])


def translate_range_for(cursor: Cursor, children: list) -> StmtIR:
    """CXXForRangeStmt → UnknownStmt("RangeForLoop", ...)。"""
    loop_var_name = "elem"
    loop_var_type = "auto"
    collection = UnknownExpr("collection")
    body_stmt = UnknownStmt("empty")

    for child in children:
        if child.kind == CursorKind.VAR_DECL:
            name = child.spelling
            if name.startswith("__range"):
                range_children = list(child.get_children())
                if range_children:
                    collection = translate_expr(range_children[0])
            elif not name.startswith("__"):
                loop_var_name = name
                loop_var_type = translate_type(child.type)
        elif child.kind == CursorKind.COMPOUND_STMT:
            body_stmt = translate_stmt(child)
        elif child.kind == CursorKind.DECL_STMT:
            for sub in child.get_children():
                if sub.kind == CursorKind.VAR_DECL:
                    name = sub.spelling
                    if name.startswith("__range"):
                        sub_children = list(sub.get_children())
                        if sub_children:
                            collection = translate_expr(sub_children[0])
                    elif not name.startswith("__"):
                        loop_var_name = name
                        loop_var_type = translate_type(sub.type)

    return UnknownStmt("RangeForLoop", [
        ExprStmt(collection),
        LetStmt(Var(loop_var_name), loop_var_type, Lit(0)),
        body_stmt,
    ])


# ============================================================
# 函数翻译
# ============================================================

def translate_function(cursor: Cursor) -> FuncIR:
    """FunctionDecl → FuncIR。"""
    name = cursor.spelling
    ret_type = translate_type(cursor.result_type)

    params = []
    body_stmts = []

    for child in cursor.get_children():
        if child.kind == CursorKind.PARM_DECL:
            pname = child.spelling or f"_arg{len(params)}"
            ptyp = translate_type(child.type)
            is_ref = "&" in child.type.spelling
            is_const = "const" in child.type.spelling
            params.append(ParamIR(
                name=pname, typ=ptyp,
                is_output=is_ref and not is_const,
                is_const_ref=is_ref and is_const
            ))
        elif child.kind == CursorKind.COMPOUND_STMT:
            for stmt_cursor in child.get_children():
                body_stmts.append(translate_stmt(stmt_cursor))

    loc = cursor.location
    return FuncIR(
        name=name, params=params, ret_type=ret_type,
        body=body_stmts,
        source_file=loc.file.name if loc.file else "",
        source_line=loc.line,
    )


# ============================================================
# 入口
# ============================================================

def parse_file(filename: str, args: list[str] = None,
               target_file: str = None) -> list[FuncIR]:
    """解析 C++ 文件，返回函数列表。

    target_file: 只提取文件名包含此子串的函数（过滤系统头文件）。
    """
    if args is None:
        args = ["-std=c++17", "-I./"]

    index = Index.create()
    tu = index.parse(filename, args=args + SYSTEM_INCLUDES)

    errors = [d for d in tu.diagnostics if d.severity >= 3]
    if errors:
        print(f"Warning: {len(errors)} compilation errors")
        for e in errors[:5]:
            print(f"  {e}")

    funcs = []
    for cursor in tu.cursor.walk_preorder():
        if cursor.kind == CursorKind.FUNCTION_DECL and cursor.is_definition():
            loc = cursor.location
            if loc.file:
                if target_file and target_file not in loc.file.name:
                    continue
                funcs.append(translate_function(cursor))

    return funcs


# ============================================================
# 辅助
# ============================================================

def _extract_binop(tokens, children) -> str:
    """从 token 列表中提取二元运算符。"""
    if len(children) >= 2:
        # 运算符在两个子表达式的 token 之间
        ops = {"+", "-", "*", "/", "%", "<<", ">>", "<", ">", "<=", ">=",
               "==", "!=", "&&", "||", "=", "+=", "-=", "*=", "/=", "%="}
        for t in tokens:
            if t.spelling in ops:
                return t.spelling
    return "?"


def _extract_compound_op(tokens) -> str:
    """从 compound assignment token 提取基础运算符。"""
    for t in tokens:
        if t.spelling in ("+=", "-=", "*=", "/=", "%=", "<<=", ">>="):
            return t.spelling[:-1]  # 去掉 =
    return "?"


def _get_func_name(cursor: Cursor) -> str:
    """从函数引用 cursor 提取函数名。"""
    if cursor.kind == CursorKind.DECL_REF_EXPR:
        return cursor.spelling
    for child in cursor.get_children():
        name = _get_func_name(child)
        if name != "?":
            return name
    return cursor.spelling or "?"


def _gen_require_name(cond: ExprIR) -> str:
    if isinstance(cond, BinOp) and cond.op == "!=":
        if isinstance(cond.lhs, Var):
            return f"h_{cond.lhs.name}_ne"
    return "h_assert"


def _wrap_body(stmt: StmtIR) -> list[StmtIR]:
    if isinstance(stmt, UnknownStmt) and stmt.kind == "CompoundStmt":
        return stmt.children
    return [stmt]
