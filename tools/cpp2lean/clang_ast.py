"""
CLPoly C++ → Lean IR 翻译器：Clang AST JSON 解析器

输入：clang++ -Xclang -ast-dump=json -fsyntax-only 产出的 JSON
输出：list[FuncIR]
"""

from __future__ import annotations
import json
import re
from ir_types import *


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
    qual = node.get("type", {}).get("qualType", "")
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
}

# CLPoly 特有类型映射
CLPOLY_TYPE_MAP = {
    "SparsePolyZp": StructType("SparsePolyZp", []),
    "upolynomial_<Zp>": StructType("SparsePolyZp", []),
    "umonomial": StructType("UMonomial", []),
    "basic_monomial": StructType("UMonomial", []),
    "Zp": StructType("Zp", []),
    "ZZ": BaseType.INT64,              # CLPoly 大整数 → Int（Lean 任意精度）
    "std::mt19937": BaseType.UINT64,   # 随机数生成器 → 参数化种子
    "std::uniform_int_distribution": BaseType.UINT64,
}

# CLPoly 方法调用映射
METHOD_MAP = {
    "empty": "isEmpty",
    "front": "front!",
    "back": "back!",
    "size": "size",
    "push_back": "push",
    "reserve": "reserve",     # 在 Lean 中是 no-op
    "normalization": "normalize",
    "deg": "deg",
    "prime": "prime",
    "comp": "comp",
}


def parse_type(qual: str) -> TypeIR:
    qual = qual.strip()
    # 去掉 const 和引用
    clean = qual.replace("const", "").replace("&", "").replace("*", "").strip()
    if clean in TYPE_MAP:
        return TYPE_MAP[clean]
    # CLPoly 类型
    for key, val in CLPOLY_TYPE_MAP.items():
        if key in clean:
            return val
    # vector<T> — 需要匹配嵌套 <>
    m = re.match(r"(?:std::)?vector<(.+)>$", clean)
    if m:
        inner_str = m.group(1)
        # 检查嵌套 <>：只在最外层的 > 处截断
        inner_str = _extract_template_arg(clean, "vector")
        if inner_str:
            return ArrayType(parse_type(inner_str))
    # pair<A,B>
    pair_args = _extract_pair_args(clean)
    if pair_args:
        return PairType(parse_type(pair_args[0]), parse_type(pair_args[1]))
    return qual  # 未知类型保留原始字符串


def parse_compound(node: dict) -> list[StmtIR]:
    return [parse_stmt(c) for c in node.get("inner", [])]


def parse_stmt(node: dict) -> StmtIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    if kind == "DeclStmt":
        stmts = []
        for decl in inner:
            if decl.get("kind") == "VarDecl":
                stmts.append(parse_var_decl(decl))
        return stmts[0] if len(stmts) == 1 else UnknownStmt("MultiDecl", stmts)

    if kind == "CompoundStmt":
        # 扁平化：不嵌套 compound
        return UnknownStmt("CompoundStmt", [parse_stmt(c) for c in inner])

    if kind == "IfStmt":
        return parse_if(node)

    if kind == "ForStmt":
        return parse_for(node)

    if kind == "WhileStmt":
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
        # 检测 operator=（赋值/move 赋值）
        if len(inner) >= 3:
            callee = inner[0]
            callee_name = callee.get("referencedDecl", {}).get("name", "")
            if "operator=" in callee_name or callee.get("kind") == "ImplicitCastExpr":
                # inner[1] = LHS, inner[2] = RHS
                # 跳过 RHS 的 ImplicitCastExpr 链（move 语义无关紧要）
                rhs_node = inner[2]
                while rhs_node.get("kind") in ("ImplicitCastExpr", "CStyleCastExpr", "MaterializeTemporaryExpr", "CXXBindTemporaryExpr"):
                    sub = rhs_node.get("inner", [])
                    if sub:
                        rhs_node = sub[0]
                    else:
                        break
                lhs = parse_expr(inner[1])
                rhs = parse_expr(rhs_node)
                if isinstance(lhs, Var):
                    return AssignStmt(lhs, rhs)
                # LHS 是 field access（如 c = expr）→ 也当赋值
                if isinstance(lhs, (FieldAccess, Cast)):
                    lhs_var = _extract_root_var(lhs)
                    if lhs_var:
                        return AssignStmt(lhs_var, rhs)
        return ExprStmt(parse_expr(node))

    # 表达式语句
    if kind in ("CallExpr", "CXXMemberCallExpr", "UnaryOperator",
                "BinaryOperator"):
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

    return UnknownStmt(kind, [parse_stmt(c) for c in inner] if inner else [])


def parse_var_decl(node: dict) -> LetStmt:
    name = node.get("name", "_")
    qual = node.get("type", {}).get("qualType", "")
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


def parse_expr(node: dict) -> ExprIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    if kind == "IntegerLiteral":
        return Lit(int(node.get("value", "0")))

    if kind == "CXXBoolLiteralExpr":
        return Lit(1 if node.get("value", False) else 0, BaseType.BOOL)

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
        if inner:
            return FieldAccess(parse_expr(inner[0]), fname)

    if kind in ("CStyleCastExpr", "ImplicitCastExpr"):
        if inner:
            # 提取目标类型
            target_qual = node.get("type", {}).get("qualType", "")
            target_type = parse_type(target_qual)
            # 追踪到最内层找真正的源类型（跳过中间的 ImplicitCastExpr 链）
            source_node = inner[0]
            while source_node.get("kind") in ("ImplicitCastExpr", "CStyleCastExpr"):
                sub = source_node.get("inner", [])
                if sub:
                    source_node = sub[0]
                else:
                    break
            source_qual = source_node.get("type", {}).get("qualType", "")
            source_type = parse_type(source_qual)
            child_expr = parse_expr(source_node)
            # 如果源和目标不同 → 保留 Cast
            if target_type != source_type:
                return Cast(child_expr, target_type, source_type)
            return child_expr

    if kind == "ParenExpr":
        if inner:
            return parse_expr(inner[0])

    if kind == "CXXMemberCallExpr":
        if inner:
            member = inner[0]
            method_name = member.get("name", "?")
            # 提取 this 对象：跳过 ImplicitCastExpr 链到最内层
            obj_node = member.get("inner", [{}])[0] if member.get("inner") else {}
            while obj_node.get("kind") in ("ImplicitCastExpr", "CStyleCastExpr"):
                sub = obj_node.get("inner", [])
                if sub:
                    obj_node = sub[0]
                else:
                    break
            obj = parse_expr(obj_node) if obj_node else UnknownExpr("no_obj")
            # 方法名映射
            lean_method = METHOD_MAP.get(method_name, method_name)
            if method_name == "push_back" and len(inner) > 1:
                elem = parse_expr(inner[1])
                return ArrayPush(obj, elem)
            # 无参 getter → 字段访问
            if method_name in ("empty", "size", "front", "back", "deg", "prime", "val"):
                return FieldAccess(obj, lean_method)
            if method_name in ("reserve",):
                return Call("reserve", [obj])  # no-op
            if method_name in ("normalization", "normalize"):
                return Call("_mutate_normalize", [obj])  # mutation method
            if method_name == "data":
                return obj  # .data() → 返回 self（Array 本身）
            args = [parse_expr(a) for a in inner[1:]]
            return Call(lean_method, [obj] + args)

    if kind == "CXXConstructExpr":
        qual_type = node.get("type", {}).get("qualType", "")
        typ = parse_type(qual_type)
        if isinstance(typ, StructType) and typ.name == "SparsePolyZp":
            return Call("Array.empty", [])
        if isinstance(typ, ArrayType):
            return Call("Array.empty", [])
        # 双参数构造（pair(a, b)）
        if inner and len(inner) == 2:
            return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
        # 单参数构造（如 umonomial(deg)）
        if inner and len(inner) == 1:
            return parse_expr(inner[0])
        # 零参数（默认构造）
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
            func_expr = parse_expr(inner[0])
            func_name = _extract_func_name(func_expr)
            # std::move → 忽略（纯函数式不需要）
            if "move" in func_name:
                return parse_expr(inner[1]) if len(inner) > 1 else UnknownExpr("move_empty")
            # std::make_pair → 直接构造 pair
            if "make_pair" in func_name and len(inner) >= 3:
                return Call("Prod.mk", [parse_expr(inner[1]), parse_expr(inner[2])])
            args = [parse_expr(a) for a in inner[1:]]
            return Call(func_name, args)

    if kind == "CXXOperatorCallExpr":
        # 非赋值的运算符调用（如 ZZ::operator/, operator==）
        if inner:
            callee = inner[0]
            ref_name = callee.get("referencedDecl", {}).get("name", "")
            if "operator/" in ref_name and len(inner) >= 3:
                return BinOp("/", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator*" in ref_name and len(inner) >= 3:
                return BinOp("*", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator+" in ref_name and len(inner) >= 3:
                return BinOp("+", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator-" in ref_name and len(inner) >= 3:
                return BinOp("-", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator==" in ref_name and len(inner) >= 3:
                return BinOp("==", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator!=" in ref_name and len(inner) >= 3:
                return BinOp("!=", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator<" in ref_name and len(inner) >= 3:
                return BinOp("<", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator>" in ref_name and len(inner) >= 3:
                return BinOp(">", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator<<" in ref_name and len(inner) >= 3:
                return BinOp("<<", parse_expr(inner[1]), parse_expr(inner[2]))
            if "operator>>" in ref_name and len(inner) >= 3:
                return BinOp(">>", parse_expr(inner[1]), parse_expr(inner[2]))

    if kind == "CompoundAssignOperator":
        op = node.get("opcode", "?=").replace("=", "")
        if len(inner) >= 2:
            return BinOp(op, parse_expr(inner[0]), parse_expr(inner[1]))

    # 兜底：递归第一个子节点
    if inner:
        return parse_expr(inner[0])

    return UnknownExpr(kind)


def _is_assert_fail(node: dict) -> bool:
    """检测 __assert_fail 调用。"""
    if node.get("kind") == "CallExpr":
        inner = node.get("inner", [])
        if inner:
            callee = inner[0]
            ref = callee.get("referencedDecl", {})
            if "__assert_fail" in ref.get("name", ""):
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
