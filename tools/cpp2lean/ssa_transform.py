"""
CLPoly C++ → Lean IR 翻译器：SSA 变换引擎

输入：FuncIR（clang_ast 产出）
输出：SSAFunc（变量单赋值，循环→TailRec，assert→Require）
"""

from __future__ import annotations
from copy import deepcopy
from ir_types import *


class VarEnv:
    """变量版本环境 + 类型环境。"""

    def __init__(self):
        self.versions: dict[str, int] = {}
        self.types: dict[str, TypeIR] = {}    # 变量名 → 类型

    def current(self, name: str) -> Var:
        return Var(name, self.versions.get(name, 0))

    def bump(self, name: str) -> Var:
        v = self.versions.get(name, 0) + 1
        self.versions[name] = v
        return Var(name, v)

    def set_type(self, name: str, typ: TypeIR):
        if typ != "auto":
            self.types[name] = typ

    def get_type(self, name: str) -> TypeIR:
        return self.types.get(name, BaseType.UINT64)  # 默认 UInt64

    def fork(self) -> VarEnv:
        new = VarEnv()
        new.versions = dict(self.versions)
        new.types = dict(self.types)
        return new

    def merge(self, other: VarEnv) -> dict[str, tuple[Var, Var]]:
        diff = {}
        all_names = set(self.versions.keys()) | set(other.versions.keys())
        for name in all_names:
            v1 = self.versions.get(name, 0)
            v2 = other.versions.get(name, 0)
            if v1 != v2:
                diff[name] = (Var(name, v1), Var(name, v2))
        return diff


# ============================================================
# 主入口
# ============================================================

def transform_func(func: FuncIR) -> SSAFunc:
    env = VarEnv()

    # 初始化参数版本 + 类型
    for p in func.params:
        env.versions[p.name] = 0
        env.set_type(p.name, p.typ)

    # 分离输出参数
    new_params, out_names = transform_params(func.params)

    # 变换函数体
    new_body = transform_body(func.body, env)

    # 收集 require + 检测 throw
    requires = collect_requires(new_body)
    has_throw = detect_throws(new_body)

    # 移除已提升到签名的 Require
    clean_body = [s for s in new_body if not isinstance(s, Require)]

    # 如果有输出参数，在末尾添加返回
    ret_type = func.ret_type
    if out_names:
        if func.ret_type == BaseType.VOID:
            # void + 输出参数 → 返回输出参数的最终版本
            if len(out_names) == 1:
                out_var = env.current(out_names[0])
                clean_body.append(ReturnStmt(out_var))
            else:
                # 多输出 → tuple（简化：暂不处理）
                pass

    if has_throw:
        ret_type = ExceptType(ret_type)

    return SSAFunc(
        name=func.name,
        params=new_params,
        ret_type=ret_type,
        requires=requires,
        body=clean_body,
        has_throw=has_throw,
    )


# ============================================================
# 参数变换
# ============================================================

def transform_params(params: list[ParamIR]) -> tuple[list[ParamIR], list[str]]:
    """分离输出参数 → 返回值。返回 (新参数列表, 输出参数名列表)。"""
    new_params = []
    out_names = []
    for p in params:
        if p.is_output:
            out_names.append(p.name)
            # 输出参数变为普通输入参数（初始值）
            new_params.append(ParamIR(p.name, p.typ, is_output=False))
        else:
            new_params.append(p)
    return new_params, out_names


# ============================================================
# 函数体变换
# ============================================================

def transform_body(stmts: list[StmtIR], env: VarEnv) -> list[StmtIR]:
    result = []
    i = 0
    while i < len(stmts):
        stmt = stmts[i]
        # 检测 early return 模式：if (cond) return; → if cond then return else {rest}
        if isinstance(stmt, IfStmt) and not stmt.else_body:
            then_has_return = any(isinstance(s, ReturnStmt) for s in stmt.then_body)
            if then_has_return and i + 1 < len(stmts):
                # 剩余语句变成 else 分支
                cond = rename_expr(stmt.cond, env)
                env_then = env.fork()
                then_body = transform_body(stmt.then_body, env_then)
                env_else = env.fork()
                else_body = transform_body(stmts[i + 1:], env_else)
                result.append(IfStmt(cond, then_body, else_body))
                env.versions = env_else.versions
                env.types.update(env_else.types)
                return result  # 剩余语句已在 else 中处理

        transformed = transform_stmt(stmt, env)
        if isinstance(transformed, list):
            result.extend(transformed)
            # 如果最后生成的是 TailRec，后面的 return 被循环替代
            if any(isinstance(t, TailRec) for t in transformed):
                break  # 循环返回值就是函数返回值
        else:
            result.append(transformed)
            if isinstance(transformed, TailRec):
                break
        i += 1
    return result


def transform_stmt(stmt: StmtIR, env: VarEnv) -> StmtIR | list[StmtIR]:
    if isinstance(stmt, LetStmt):
        # let x : T := e → let x_1 : T := rename(e)
        new_value = rename_expr(stmt.value, env)
        var = env.bump(stmt.var.name) if stmt.var.name != "_" else stmt.var
        env.set_type(stmt.var.name, stmt.typ)
        return LetStmt(var, stmt.typ, new_value)

    if isinstance(stmt, AssignStmt):
        # x = e → let x_{n+1} : T := rename(e)（T 从之前的声明获取）
        new_value = rename_expr(stmt.value, env)
        typ = env.get_type(stmt.target.name)
        new_var = env.bump(stmt.target.name)
        return LetStmt(new_var, typ, new_value)

    if isinstance(stmt, IfStmt):
        return transform_if(stmt, env)

    if isinstance(stmt, ReturnStmt):
        if stmt.value is not None:
            return ReturnStmt(rename_expr(stmt.value, env))
        return stmt

    if isinstance(stmt, Require):
        return Require(rename_expr(stmt.cond, env), stmt.name, stmt.source)

    if isinstance(stmt, Throw):
        return stmt

    if isinstance(stmt, ExprStmt):
        expr = stmt.expr
        # i++ / ++i → let i_{n+1} := i_n + 1
        if isinstance(expr, UnaryOp) and expr.op == "++" and isinstance(expr.operand, Var):
            old = env.current(expr.operand.name)
            new_var = env.bump(expr.operand.name)
            typ = env.get_type(expr.operand.name)
            return LetStmt(new_var, typ, BinOp("+", old, Lit(1)))
        if isinstance(expr, UnaryOp) and expr.op == "--" and isinstance(expr.operand, Var):
            old = env.current(expr.operand.name)
            new_var = env.bump(expr.operand.name)
            typ = env.get_type(expr.operand.name)
            return LetStmt(new_var, typ, BinOp("-", old, Lit(1)))
        # _mutate_normalize(obj) → let obj_{n+1} := normalize obj_n
        if isinstance(expr, Call) and expr.func == "_mutate_normalize" and expr.args:
            if isinstance(expr.args[0], Var):
                old = env.current(expr.args[0].name)
                new_var = env.bump(expr.args[0].name)
                typ = env.get_type(expr.args[0].name)
                return LetStmt(new_var, typ, Call("normalize", [old]))
        # __upoly_make_monic(obj) → let obj_{n+1} := upoly_make_monic obj_n
        if isinstance(expr, Call) and expr.func == "__upoly_make_monic" and expr.args:
            if isinstance(expr.args[0], Var):
                old = env.current(expr.args[0].name)
                new_var = env.bump(expr.args[0].name)
                typ = env.get_type(expr.args[0].name)
                return LetStmt(new_var, typ, Call("__upoly_make_monic", [old]))
        # arr.push(x) → let arr_{n+1} := arr_n.push(x)
        if isinstance(expr, ArrayPush) and isinstance(expr.arr, Var):
            old = env.current(expr.arr.name)
            new_var = env.bump(expr.arr.name)
            typ = env.get_type(expr.arr.name)
            return LetStmt(new_var, typ, ArrayPush(old, rename_expr(expr.elem, env)))
        # field = expr（成员赋值）→ functional update + Array.set
        if isinstance(expr, BinOp) and expr.op == "=":
            if isinstance(expr.lhs, FieldAccess):
                return transform_member_assign(expr.lhs, expr.rhs, env)
        return ExprStmt(rename_expr(expr, env))

    if isinstance(stmt, UnknownStmt):
        if stmt.kind == "ForLoop":
            return transform_for_loop(stmt, env)
        if stmt.kind == "WhileLoop":
            return transform_while_loop(stmt, env)
        if stmt.kind == "RangeForLoop":
            return transform_range_for(stmt, env)
        if stmt.kind == "CompoundStmt":
            return transform_body(stmt.children, env)
        if stmt.kind == "MultiDecl":
            return transform_body(stmt.children, env)
        if stmt.kind == "BreakStmt":
            return stmt  # 由 transform_for_loop 处理
        if stmt.kind == "ContinueStmt":
            return stmt  # 由 transform_for_loop 处理
        return stmt

    return stmt


# ============================================================
# if/else 变换（phi 节点）
# ============================================================

def transform_if(stmt: IfStmt, env: VarEnv) -> list[StmtIR]:
    cond = rename_expr(stmt.cond, env)

    # fork 两个分支
    env_then = env.fork()
    env_else = env.fork()

    then_body = transform_body(stmt.then_body, env_then)
    else_body = transform_body(stmt.else_body, env_else) if stmt.else_body else []

    # 检测两分支中是否有 return/break（提前退出）
    then_returns = any(isinstance(s, ReturnStmt) for s in then_body)
    else_returns = any(isinstance(s, ReturnStmt) for s in else_body)

    # 如果 then 分支 return，else 分支继续 → 不需要 phi
    if then_returns and not else_returns:
        result = [IfStmt(cond, then_body, else_body)]
        # else 分支的 env 成为后续的 env
        env.versions = env_else.versions
        return result

    if else_returns and not then_returns:
        result = [IfStmt(cond, then_body, else_body)]
        env.versions = env_then.versions
        return result

    if then_returns and else_returns:
        return [IfStmt(cond, then_body, else_body)]

    # 两分支都不 return → 需要 phi
    diff = env_then.merge(env_else)

    if not diff:
        return [IfStmt(cond, then_body, else_body)]

    # 有变量分歧 → 提取各分支对变量的赋值表达式，生成 phi
    then_exprs = _extract_last_assigns(then_body)
    else_exprs = _extract_last_assigns(else_body)

    # 按依赖顺序排列 phi 节点：如果 A 的表达式引用了 B，B 先输出
    result: list[StmtIR] = []
    remaining = dict(diff)
    emitted = set()
    max_iter = len(remaining) + 1
    while remaining and max_iter > 0:
        max_iter -= 1
        for name in list(remaining.keys()):
            then_var, else_var = remaining[name]
            then_expr = then_exprs.get(name, then_var)
            else_expr = else_exprs.get(name, else_var)
            # 检查依赖：表达式中引用的 diff 变量是否已 emit
            deps = _extract_vars(then_expr) | _extract_vars(else_expr)
            unresolved = deps & set(remaining.keys()) - {name} - emitted
            if not unresolved:
                new_var = env.bump(name)
                typ = env.get_type(name)
                phi_expr = CondExpr(cond, then_expr, else_expr)
                result.append(LetStmt(new_var, typ, phi_expr))
                emitted.add(name)
                del remaining[name]

    # 剩余的（循环依赖）直接输出
    for name, (then_var, else_var) in remaining.items():
        new_var = env.bump(name)
        typ = env.get_type(name)
        then_expr = then_exprs.get(name, then_var)
        else_expr = else_exprs.get(name, else_var)
        result.append(LetStmt(new_var, typ, CondExpr(cond, then_expr, else_expr)))

    return result


def _extract_last_assigns(stmts: list[StmtIR]) -> dict[str, ExprIR]:
    """从语句列表中提取每个变量最后一次赋值的 RHS 表达式。"""
    result = {}
    for s in stmts:
        if isinstance(s, LetStmt):
            result[s.var.name] = s.value
        elif isinstance(s, AssignStmt):
            result[s.target.name] = s.value
    return result


# ============================================================
# 循环变换
# ============================================================

def transform_for_loop(stmt: UnknownStmt, env: VarEnv) -> StmtIR | list[StmtIR]:
    """ForLoop → init_stmts + TailRec。

    stmt.children = [init, cond_expr, step, body]
    """
    children = stmt.children
    if len(children) < 4:
        return stmt

    init_stmt = children[0]
    cond_stmt = children[1]
    step_stmt = children[2]
    body_stmt = children[3]

    # 提取条件表达式
    cond_expr = Lit(1, BaseType.BOOL)
    if isinstance(cond_stmt, ExprStmt):
        cond_expr = cond_stmt.expr

    # 识别循环变量
    loop_vars = identify_loop_vars(body_stmt)

    # 提取 init（在 env 中注册）
    init_stmts = []
    if isinstance(init_stmt, LetStmt):
        transformed_init = transform_stmt(init_stmt, env)
        init_stmts.append(transformed_init)
        loop_vars.add(init_stmt.var.name)
    elif isinstance(init_stmt, AssignStmt):
        transformed_init = transform_stmt(init_stmt, env)
        init_stmts.append(transformed_init)
        loop_vars.add(init_stmt.target.name)

    # 检测 step 中修改的变量（如 i++）
    step_vars = identify_loop_vars_from_step(step_stmt)
    loop_vars |= step_vars

    # 构建 TailRec 参数：循环变量的当前版本 + 类型
    params = []
    for name in sorted(loop_vars):
        var = env.current(name)
        typ = env.get_type(name)
        params.append((var, typ))

    # rename 条件表达式中的变量
    cond_renamed = rename_expr(cond_expr, env)
    exit_cond = negate_expr(cond_renamed)

    # 检测 break 条件（在 rename 前）
    break_cond_raw = extract_break_cond(body_stmt)
    break_cond = rename_expr(break_cond_raw, env) if break_cond_raw else None

    # 循环体 SSA 变换
    loop_env = env.fork()
    body_stmts_raw = flatten_stmts(body_stmt)
    body_stmts = transform_body(body_stmts_raw, loop_env)

    # 步进 SSA 变换
    step_stmts = transform_step(step_stmt, loop_env)

    func_name = f"loop_{abs(hash(str(stmt))) % 100000}"

    tailrec = TailRec(
        func_name=func_name,
        params=params,
        exit_cond=exit_cond,
        break_cond=break_cond,
        body=body_stmts,
        step=step_stmts,
    )

    result = init_stmts + [tailrec]
    return result if len(result) > 1 else result[0]


def identify_loop_vars_from_step(step_stmt: StmtIR) -> set[str]:
    """从步进语句（如 i++）中提取变量。"""
    result = set()
    if isinstance(step_stmt, ExprStmt):
        expr = step_stmt.expr
        # i++ / ++i
        if isinstance(expr, UnaryOp) and expr.op in ("++", "--"):
            if isinstance(expr.operand, Var):
                result.add(expr.operand.name)
        # i += 1
        if isinstance(expr, BinOp) and expr.op in ("+", "-"):
            if isinstance(expr.lhs, Var):
                result.add(expr.lhs.name)
    if isinstance(step_stmt, AssignStmt):
        result.add(step_stmt.target.name)
    return result


def transform_step(step_stmt: StmtIR, env: VarEnv) -> list[StmtIR]:
    """翻译步进语句（i++ → let i_{n+1} := i_n + 1）。"""
    if isinstance(step_stmt, UnknownStmt) and step_stmt.kind == "null":
        return []

    if isinstance(step_stmt, ExprStmt):
        expr = step_stmt.expr
        # i++ / ++i → let i_{n+1} := i_n + 1
        if isinstance(expr, UnaryOp) and expr.op in ("++",):
            if isinstance(expr.operand, Var):
                old_var = env.current(expr.operand.name)
                new_var = env.bump(expr.operand.name)
                typ = env.get_type(expr.operand.name)
                return [LetStmt(new_var, typ, BinOp("+", old_var, Lit(1)))]
        if isinstance(expr, UnaryOp) and expr.op in ("--",):
            if isinstance(expr.operand, Var):
                old_var = env.current(expr.operand.name)
                new_var = env.bump(expr.operand.name)
                typ = env.get_type(expr.operand.name)
                return [LetStmt(new_var, typ, BinOp("-", old_var, Lit(1)))]

    # 通用：当作普通语句处理
    result = transform_stmt(step_stmt, env)
    if isinstance(result, list):
        return result
    return [result]


def transform_while_loop(stmt: UnknownStmt, env: VarEnv) -> StmtIR:
    """WhileLoop → TailRec。"""
    children = stmt.children
    if len(children) < 2:
        return stmt

    cond_expr = Lit(1, BaseType.BOOL)
    if isinstance(children[0], ExprStmt):
        cond_expr = children[0].expr

    body_stmt = children[1]

    # 识别循环变量：只包含被赋值的变量（不包含只被读取的）
    loop_vars = identify_loop_vars(body_stmt)

    # rename 条件
    cond_renamed = rename_expr(cond_expr, env)
    break_cond = extract_break_cond(body_stmt)
    break_cond_renamed = rename_expr(break_cond, env) if break_cond else None

    # 从循环体中推断未知变量的类型
    _infer_types_from_body(body_stmt, env)
    params = [(env.current(name), env.get_type(name)) for name in sorted(loop_vars)]

    # SSA 变换循环体
    loop_env = env.fork()
    body_stmts = transform_body(flatten_stmts(body_stmt), loop_env)

    return TailRec(
        func_name=f"while_{abs(hash(str(stmt))) % 100000}",
        params=params,
        exit_cond=negate_expr(cond_renamed),
        break_cond=break_cond_renamed,
        body=body_stmts,
        step=[],
    )


def transform_range_for(stmt: UnknownStmt, env: VarEnv) -> list[StmtIR]:
    """RangeForLoop → partial def 遍历 Array。

    for (auto& elem : arr) { body }
    → let rec range_loop (idx : Nat) (acc : ...) :=
        if idx >= arr.size then acc
        else let elem := arr[idx]!; body; range_loop (idx+1) acc'
    """
    children = stmt.children
    if len(children) < 3:
        return [stmt]

    collection = children[0].expr if isinstance(children[0], ExprStmt) else UnknownExpr("?")
    # children[1] 是 LetStmt(var, type, _) — 循环变量名+类型
    if isinstance(children[1], LetStmt):
        loop_var_name = children[1].var.name
        loop_var_type = children[1].typ
    else:
        loop_var_name = "elem"
        loop_var_type = "auto"
    body_stmt = children[2]

    collection_renamed = rename_expr(collection, env)

    # 识别 body 中被修改的变量
    body_vars = identify_loop_vars(body_stmt)

    idx_name = "__idx"
    env.versions[idx_name] = 0
    env.set_type(idx_name, BaseType.UINT64)
    idx_var = env.bump(idx_name)

    # 给集合一个 SSA 名称（避免 UNKNOWN）
    coll_var = Var("__coll")
    env.versions["__coll"] = 0
    coll_var = env.bump("__coll")

    # 循环体前添加 let elem := __coll[idx]
    elem_let = LetStmt(Var(loop_var_name), loop_var_type,
                        ArrayAccess(coll_var, idx_var))
    env.set_type(loop_var_name, loop_var_type)

    # __coll 类型
    coll_type = env.get_type(collection.name) if isinstance(collection, Var) else "auto"
    env.set_type("__coll", coll_type if coll_type != "auto" else BaseType.UINT64)

    # 循环体 SSA 变换（可能修改 __coll）
    loop_env = env.fork()
    loop_env.versions[loop_var_name] = 1
    body_stmts = [elem_let] + transform_body(flatten_stmts(body_stmt), loop_env)

    # 循环参数 = idx + body 中修改的变量（含 __coll 如果被 Array.set 修改）
    all_loop_vars = set(body_vars)
    # 扫描 body_stmts 中实际产生的 LetStmt 变量
    for s in body_stmts:
        if isinstance(s, LetStmt) and s.var.name in ("__coll",):
            all_loop_vars.add(s.var.name)

    params = [(idx_var, BaseType.UINT64)]
    for name in sorted(all_loop_vars):
        params.append((env.current(name), env.get_type(name)))

    func_name = f"range_{abs(hash(str(stmt))) % 100000}"

    tailrec = TailRec(
        func_name=func_name,
        params=params,
        exit_cond=BinOp(">=", idx_var, FieldAccess(coll_var, "size")),
        break_cond=None,
        body=body_stmts,
        step=[LetStmt(Var(idx_name, idx_var.version + 1), BaseType.UINT64,
                       BinOp("+", idx_var, Lit(1)))],
    )

    coll_init = LetStmt(coll_var, coll_type, collection_renamed)
    idx_init = LetStmt(idx_var, BaseType.UINT64, Lit(0))
    return [coll_init, idx_init, tailrec]


def transform_member_assign(lhs: FieldAccess, rhs: ExprIR, env: VarEnv) -> StmtIR | list[StmtIR]:
    """成员赋值 elem.field1.field2 = expr → functional update。

    检测模式：
    - term.second.val = expr → { term with snd := { term.snd with val := expr } }
    - term.field = expr → { term with field := expr }

    如果 term 来自 __coll[__idx]（range-for），追加 Array.set。
    """
    renamed_rhs = rename_expr(rhs, env)

    # 收集字段路径：term.second.val → [("second","snd"), ("val","val")]
    path = []
    node = lhs
    while isinstance(node, FieldAccess):
        # C++ field → Lean field 映射
        lean_field = {
            "first": "fst", "second": "snd",
            "val": "val", "deg": "deg", "p": "p",
        }.get(node.field_name, node.field_name)
        path.append((node.field_name, lean_field))
        node = node.obj
    path.reverse()

    # node 现在是根变量（如 term）
    root_var = rename_expr(node, env) if isinstance(node, Var) else node
    root_name = node.name if isinstance(node, Var) else None

    if not root_name or not path:
        return ExprStmt(BinOp(":=", rename_expr(lhs, env), renamed_rhs))

    # 构建 functional update 表达式
    update_expr = _build_functional_update(root_var, path, renamed_rhs)

    # 创建新版本的 root 变量
    new_root = env.bump(root_name)
    typ = env.get_type(root_name)
    result = [LetStmt(new_root, typ, update_expr)]

    # 如果 root 来自 __coll[__idx]（range-for 模式），更新数组
    if "__coll" in env.versions and "__idx" in env.versions:
        coll_var = env.current("__coll")
        idx_var = env.current("__idx")
        new_coll = env.bump("__coll")
        coll_typ = env.get_type("__coll")
        set_expr = Call("Array.set!", [coll_var, idx_var, new_root])
        result.append(LetStmt(new_coll, coll_typ, set_expr))

    return result


def _build_functional_update(root: ExprIR, path: list[tuple[str, str]], value: ExprIR) -> ExprIR:
    """构建嵌套的 functional update。

    path = [("second", "snd"), ("val", "val")]
    → Call("struct_update", [root, "snd", Call("struct_update", [root.snd, "val", value])])

    在 Lean 中表示为 { root with snd := { root.snd with val := value } }
    用 Call("_with", ...) 表示，gen_expr 特殊处理。
    """
    if len(path) == 1:
        _, lean_field = path[0]
        return Call("_with", [root, Var(lean_field), value])
    else:
        _, lean_field = path[0]
        inner_obj = FieldAccess(root, lean_field)
        inner_update = _build_functional_update(inner_obj, path[1:], value)
        return Call("_with", [root, Var(lean_field), inner_update])


def _infer_types_from_body(stmt: StmtIR, env: VarEnv):
    """从循环体中的 LetStmt 提取变量类型到 env。"""
    stmts = flatten_stmts(stmt) if not isinstance(stmt, list) else stmt
    for s in stmts:
        if isinstance(s, LetStmt) and s.typ != "auto":
            env.set_type(s.var.name, s.typ)
        if isinstance(s, IfStmt):
            _infer_types_from_body_list(s.then_body, env)
            _infer_types_from_body_list(s.else_body, env)


def _infer_types_from_body_list(stmts: list[StmtIR], env: VarEnv):
    for s in stmts:
        if isinstance(s, LetStmt) and s.typ != "auto":
            env.set_type(s.var.name, s.typ)


def _extract_vars(expr: ExprIR) -> set[str]:
    """提取表达式中的变量名。"""
    result = set()
    if isinstance(expr, Var):
        result.add(expr.name)
    elif isinstance(expr, BinOp):
        result |= _extract_vars(expr.lhs) | _extract_vars(expr.rhs)
    elif isinstance(expr, UnaryOp):
        result |= _extract_vars(expr.operand)
    elif isinstance(expr, CondExpr):
        result |= _extract_vars(expr.cond) | _extract_vars(expr.then_e) | _extract_vars(expr.else_e)
    elif isinstance(expr, Call):
        for a in expr.args:
            result |= _extract_vars(a)
    elif isinstance(expr, ArrayAccess):
        result |= _extract_vars(expr.arr) | _extract_vars(expr.idx)
    return result


def identify_loop_vars(stmt: StmtIR) -> set[str]:
    """识别语句中被**修改**（而非新声明）的变量名。

    AssignStmt = 修改已有变量 → 加入
    LetStmt = 新声明 → 不加入（循环体内的局部变量不是循环参数）
    """
    result = set()
    if isinstance(stmt, AssignStmt):
        result.add(stmt.target.name)
    elif isinstance(stmt, (IfStmt,)):
        for s in stmt.then_body + stmt.else_body:
            result |= identify_loop_vars(s)
    elif isinstance(stmt, UnknownStmt):
        for c in stmt.children:
            result |= identify_loop_vars(c)
    elif isinstance(stmt, ExprStmt):
        expr = stmt.expr
        if isinstance(expr, BinOp) and expr.op == "=":
            if isinstance(expr.lhs, Var):
                result.add(expr.lhs.name)
        # i++ / i--
        if isinstance(expr, UnaryOp) and expr.op in ("++", "--"):
            if isinstance(expr.operand, Var):
                result.add(expr.operand.name)
        # arr.push(x) → arr 被修改
        if isinstance(expr, ArrayPush) and isinstance(expr.arr, Var):
            result.add(expr.arr.name)
        # method call 可能修改 this 对象
        if isinstance(expr, Call) and expr.args and isinstance(expr.args[0], Var):
            # normalize, reserve 等方法修改第一个参数
            if expr.func in ("normalize", "reserve", "normalization"):
                result.add(expr.args[0].name)
    return result


def extract_break_cond(stmt: StmtIR) -> ExprIR | None:
    """提取循环体中第一个 if(cond) break 的条件。"""
    stmts = flatten_stmts(stmt)
    for s in stmts:
        if isinstance(s, IfStmt) and s.then_body:
            if any(isinstance(t, UnknownStmt) and t.kind == "BreakStmt"
                   for t in s.then_body):
                return s.cond
    return None


# ============================================================
# 辅助
# ============================================================

def rename_expr(expr: ExprIR, env: VarEnv) -> ExprIR:
    """将表达式中的变量引用替换为当前版本。"""
    if isinstance(expr, Var):
        return env.current(expr.name)
    if isinstance(expr, BinOp):
        return BinOp(expr.op, rename_expr(expr.lhs, env), rename_expr(expr.rhs, env))
    if isinstance(expr, UnaryOp):
        return UnaryOp(expr.op, rename_expr(expr.operand, env))
    if isinstance(expr, CondExpr):
        return CondExpr(rename_expr(expr.cond, env),
                        rename_expr(expr.then_e, env),
                        rename_expr(expr.else_e, env))
    if isinstance(expr, Call):
        return Call(expr.func, [rename_expr(a, env) for a in expr.args])
    if isinstance(expr, ArrayAccess):
        return ArrayAccess(rename_expr(expr.arr, env), rename_expr(expr.idx, env))
    if isinstance(expr, FieldAccess):
        return FieldAccess(rename_expr(expr.obj, env), expr.field_name)
    if isinstance(expr, Cast):
        return Cast(rename_expr(expr.expr, env), expr.target_type, expr.source_type)
    if isinstance(expr, ArrayPush):
        return ArrayPush(rename_expr(expr.arr, env), rename_expr(expr.elem, env))
    return expr


def negate_expr(expr: ExprIR) -> ExprIR:
    """取反表达式。"""
    if isinstance(expr, BinOp):
        neg_map = {"<": ">=", ">": "<=", "<=": ">", ">=": "<", "==": "!=", "!=": "=="}
        if expr.op in neg_map:
            return BinOp(neg_map[expr.op], expr.lhs, expr.rhs)
    return UnaryOp("!", expr)


def flatten_stmts(stmt: StmtIR) -> list[StmtIR]:
    """将 CompoundStmt 展开为扁平列表。"""
    if isinstance(stmt, UnknownStmt) and stmt.kind == "CompoundStmt":
        result = []
        for c in stmt.children:
            result.extend(flatten_stmts(c))
        return result
    return [stmt]


def collect_requires(stmts: list[StmtIR]) -> list[Require]:
    """提取所有 Require 节点。"""
    result = []
    for s in stmts:
        if isinstance(s, Require):
            result.append(s)
    return result


def detect_throws(stmts: list[StmtIR]) -> bool:
    """检测是否包含 Throw。"""
    for s in stmts:
        if isinstance(s, Throw):
            return True
        if isinstance(s, IfStmt):
            if detect_throws(s.then_body) or detect_throws(s.else_body):
                return True
        if isinstance(s, UnknownStmt):
            if detect_throws(s.children):
                return True
    return False
