"""
CLPoly C++ → Lean IR 翻译器：SSA 变换引擎

输入：FuncIR（clang_ast 产出）
输出：SSAFunc（变量单赋值，循环→TailRec，assert→Require）
"""

from __future__ import annotations
from copy import deepcopy
from ir_types import *


class FuncCtx:
    """函数级上下文：收集提取的循环函数。"""
    def __init__(self):
        self.aux_defs: list[SSAFunc] = []


# === ExprStmt 黑名单：已知无副作用的表达式 ===
_KNOWN_NOOP_FUNCS = {
    # CLASS_MAP noop 方法
    "normalization", "comp_ptr", "reserve", "comp",
    "MvPolyZp.normalization", "MvPolyZZ.normalization",
    "MvMonomial.normalization", "SparsePolyZp.normalization",
    "SparsePolyZZ.normalization",
    # sort 占位
    "SparsePolyZp.toList", "SparsePolyZZ.toList",
    # 副作用已由其他机制处理
    "poly_convert", "SparsePolyZZ.compactNonzero",
    # iota init → result 不需要（初始化数组已由构造函数处理）
    "List.range",
}

def _is_known_noop(expr) -> bool:
    """判断 ExprStmt 中的表达式是否已知无副作用。"""
    if isinstance(expr, (Var, Lit)):
        return True
    if isinstance(expr, CondExpr):
        return True  # assert 检查
    if isinstance(expr, (FieldAccess, ArrayAccess, Cast)):
        return True  # 纯读取
    if isinstance(expr, Call):
        func = expr.func if isinstance(expr.func, str) else ""
        if func in _KNOWN_NOOP_FUNCS:
            return True
        # identity 函数（如 id, std::move 的翻译结果）
        if func == "id":
            return True
    if isinstance(expr, UnaryOp) and expr.op in ("*", "->"):
        return True  # 指针解引用作为语句 = noop
    return False


def _has_iterator_ops(stmts: list) -> bool:
    """检测语句列表中是否含迭代器模式标志。

    标志：变量名含 'it' 或 'out' + 成员赋值（BinOp(":=")）
    或 ExprStmt 含 Call("Array.pop") 等 erase 操作。
    """
    def check_expr(expr):
        if isinstance(expr, BinOp) and expr.op == ":=":
            return True  # 成员赋值 = 迭代器写操作
        if isinstance(expr, UnaryOp) and expr.op in ("->", "++"):
            return True
        return False

    for s in stmts:
        if isinstance(s, LetStmt):
            if check_expr(s.value):
                return True
        if isinstance(s, ExprStmt):
            if isinstance(s.expr, BinOp) and s.expr.op == ":=":
                return True
            if check_expr(s.expr):
                return True
        if isinstance(s, IfStmt):
            if _has_iterator_ops(s.then_body) or _has_iterator_ops(s.else_body):
                return True
        if isinstance(s, ReturnStmt) and s.value:
            if check_expr(s.value):
                return True
    return False


def _replace_iterator_body(loop_func: 'SSAFunc', outer_env: 'VarEnv') -> None:
    """将含迭代器操作的循环函数体替换为 filter/filterMap。

    迭代器循环的语义：遍历数组，对每个元素做变换+过滤。
    替换为：arr.filter (fun term => term.snd != 0) 或 arr.filterMap (fun term => ...)
    """
    # 找到数组参数（通常名为 e、ep、f 等，类型含 Array 或 SparsePolyZZ）
    arr_param = None
    for p in loop_func.params:
        if "Array" in str(p.typ) or "Sparse" in str(p.typ) or "MvPoly" in str(p.typ):
            arr_param = p
            break
    if arr_param is None and loop_func.params:
        arr_param = loop_func.params[0]  # fallback

    if arr_param is None:
        return

    # 替换 body 为 compactNonzero（过滤零系数）
    arr_var = Var(arr_param.name.split("_")[0],
                  int(arr_param.name.split("_")[-1]) if "_" in arr_param.name else 0)
    # 尝试解析参数名
    try:
        parts = arr_param.name.rsplit("_", 1)
        arr_var = Var(parts[0], int(parts[1]))
    except (ValueError, IndexError):
        arr_var = Var(arr_param.name, 0)

    loop_func.body = [ReturnStmt(Call("SparsePolyZZ.compactNonzero", [arr_var]))]
    loop_func.ret_type = outer_env.get_type(arr_var.name) if arr_var.name in outer_env.versions else "auto"


def finalize_loop_func(loop_func: 'SSAFunc', outer_env: 'VarEnv',
                       call_stmts: list) -> None:
    """后处理 pass：确保循环函数无自由变量（设计 §5.1）。

    扫描 loop_func.body 的全部自由变量，将遗漏的追加为参数。
    同时更新 call_stmts 中的调用（添加对应参数）。
    检测并替换含迭代器操作的循环体。
    """
    # 迭代器检测：仅当参数名含 it/out 双指针模式时替换为 compactNonzero
    # 更保守的检测：避免误替换非 compact 循环（如 __hensel_step Part 2）
    param_name_set = {p.name.split("_")[0] for p in loop_func.params}
    is_dual_pointer = ("it" in param_name_set and "out" in param_name_set)
    if is_dual_pointer and _has_iterator_ops(loop_func.body):
        _replace_iterator_body(loop_func, outer_env)
        return  # 替换后不需要闭包检测

    param_names = {p.name for p in loop_func.params}
    free = collect_free_vars(loop_func.body)

    missing = []
    for vname, vver in sorted(free):
        lean_name = Var(vname, vver).lean_name()
        if lean_name not in param_names:
            # 只追加 outer_env 中已知的变量（排除迭代器裸标记等）
            if vname in outer_env.versions:
                missing.append((vname, vver, lean_name))

    if not missing:
        return

    for vname, vver, lean_name in missing:
        typ = outer_env.get_type(vname)
        loop_func.params.append(ParamIR(lean_name, typ))

    # 更新调用点：在 call_stmts 中找到 Call(loop_func.name + "_ir", ...) 并追加参数
    func_call_name = loop_func.name + "_ir"
    for stmt in call_stmts:
        _add_missing_args_to_calls(stmt, func_call_name,
                                    [Var(vn, vv) for vn, vv, _ in missing])
    # 更新递归调用（body 内）
    for stmt in loop_func.body:
        _add_missing_args_to_calls(stmt, func_call_name,
                                    [Var(vn, vv) for vn, vv, _ in missing])


def _add_missing_args_to_calls(stmt, func_name: str, extra_args: list) -> None:
    """递归查找 Call(func_name, ...) 并追加 extra_args。"""
    if isinstance(stmt, LetStmt):
        _add_missing_args_to_call_expr(stmt.value, func_name, extra_args)
    elif isinstance(stmt, ReturnStmt) and stmt.value:
        _add_missing_args_to_call_expr(stmt.value, func_name, extra_args)
    elif isinstance(stmt, IfStmt):
        for s in stmt.then_body + stmt.else_body:
            _add_missing_args_to_calls(s, func_name, extra_args)
    elif isinstance(stmt, ExprStmt):
        _add_missing_args_to_call_expr(stmt.expr, func_name, extra_args)
    elif isinstance(stmt, list):
        for s in stmt:
            _add_missing_args_to_calls(s, func_name, extra_args)


def _add_missing_args_to_call_expr(expr, func_name: str, extra_args: list) -> None:
    """在表达式树中查找 Call(func_name) 并追加参数。"""
    if isinstance(expr, Call) and expr.func == func_name:
        expr.args.extend(extra_args)
    elif isinstance(expr, Call):
        for a in expr.args:
            _add_missing_args_to_call_expr(a, func_name, extra_args)
    elif isinstance(expr, BinOp):
        _add_missing_args_to_call_expr(expr.lhs, func_name, extra_args)
        _add_missing_args_to_call_expr(expr.rhs, func_name, extra_args)
    elif isinstance(expr, CondExpr):
        _add_missing_args_to_call_expr(expr.cond, func_name, extra_args)
        _add_missing_args_to_call_expr(expr.then_e, func_name, extra_args)
        _add_missing_args_to_call_expr(expr.else_e, func_name, extra_args)


class LoopCtx:
    """循环函数内部上下文：break/continue 语义替换用。"""
    def __init__(self, func_name: str, modified_vars: list[str],
                 step_stmts: list = None, env: 'VarEnv' = None,
                 all_param_names: list[str] = None,
                 has_parent_return: bool = False):
        self.func_name = func_name
        self.modified_vars = modified_vars
        self.step_stmts = step_stmts or []
        self.env = env
        self.all_param_names = all_param_names or []  # 全部参数名（含 idx、coll、闭包）
        self.has_parent_return = has_parent_return  # 循环体含父函数 return


class VarEnv:
    """变量版本环境 + 类型环境。"""

    def __init__(self):
        self.versions: dict[str, int] = {}
        self.types: dict[str, TypeIR] = {}    # 变量名 → 类型
        self.lambda_modified: dict[str, list[str]] = {}  # S4: lambda 变量名 → 修改的 capture 列表

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
        new.lambda_modified = dict(self.lambda_modified)
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

# 函数体整体替换表（设计 §6.5）
# 对整个函数体是不可翻译模式的函数，用预定义的语义等价 Lean body 替换
def _override_upoly_mod_coeff(params, env):
    f_cur = env.current("f")
    m_cur = env.current("m")
    f_new = env.bump("f")
    return [LetStmt(f_new, env.get_type("f"),
                    Call("SparsePolyZZ.modCoeff", [f_cur, m_cur]))]

def _override_upoly_divmod_mod(params, env):
    q_cur = env.current("q")
    r_cur = env.current("r")
    f_cur = env.current("f")
    g_cur = env.current("g")
    m_cur = env.current("m")
    stmts = [ExprStmt(Call("pair_vec_div", [q_cur, r_cur, f_cur, g_cur]))]
    q_new = env.bump("q")
    stmts.append(LetStmt(q_new, env.get_type("q"),
                          Call("SparsePolyZZ.modCoeff", [q_cur, m_cur])))
    r_new = env.bump("r")
    stmts.append(LetStmt(r_new, env.get_type("r"),
                          Call("SparsePolyZZ.modCoeff", [r_cur, m_cur])))
    return stmts

FUNC_BODY_OVERRIDE = {
    "__upoly_mod_coeff": _override_upoly_mod_coeff,
    "__upoly_divmod_mod": _override_upoly_divmod_mod,
}


def transform_func(func: FuncIR) -> SSAFunc:
    env = VarEnv()

    # 初始化参数版本 + 类型
    for p in func.params:
        env.versions[p.name] = 0
        env.set_type(p.name, p.typ)

    # 分离输出参数
    new_params, out_names = transform_params(func.params)

    # 检查函数体整体替换（设计 §6.5）
    if func.name in FUNC_BODY_OVERRIDE:
        override_fn = FUNC_BODY_OVERRIDE[func.name]
        new_body = override_fn(func.params, env)
        ctx = FuncCtx()
    else:
        # 正常变换函数体
        ctx = FuncCtx()
        new_body = transform_body(func.body, env, ctx)

    # 收集 require + 检测 throw
    requires = collect_requires(new_body)
    has_throw = detect_throws(new_body)

    # 移除已提升到签名的 Require
    clean_body = [s for s in new_body if not isinstance(s, Require)]

    # 输出参数签名改造（设计 §3.4）
    ret_type = func.ret_type
    if out_names:
        out_vars = [env.current(n) for n in out_names]
        if func.ret_type == BaseType.VOID:
            # void + 输出参数 → 返回输出参数
            if len(out_names) == 1:
                ret_type = env.get_type(out_names[0])
                clean_body.append(ReturnStmt(out_vars[0]))
            else:
                # 多输出 → Prod
                ret_type = env.get_type(out_names[0])  # codegen 用 auto
                clean_body.append(ReturnStmt(Call("Prod.mk", out_vars)))
        else:
            # 非 void + 输出参数 → 包装所有 return 路径为 Prod(原返回值, 输出参数...)
            def _wrap_all_returns(stmts, ovs):
                result = []
                for s in stmts:
                    if isinstance(s, ReturnStmt) and s.value:
                        result.append(ReturnStmt(Call("Prod.mk", [s.value] + ovs)))
                    elif isinstance(s, IfStmt):
                        result.append(IfStmt(s.cond,
                                             _wrap_all_returns(s.then_body, ovs),
                                             _wrap_all_returns(s.else_body, ovs)))
                    else:
                        result.append(s)
                return result
            clean_body = _wrap_all_returns(clean_body, out_vars)
            # 如果没有任何 return 被包装，追加 fallback
            has_any_return = any(isinstance(s, ReturnStmt) for s in clean_body)
            if not has_any_return:
                clean_body.append(ReturnStmt(Call("Prod.mk", [Lit(0)] + out_vars)))
            ret_type = env.get_type(out_names[0])  # codegen 用 auto

    if has_throw:
        ret_type = ExceptType(ret_type)

    return SSAFunc(
        name=func.name,
        params=new_params,
        ret_type=ret_type,
        requires=requires,
        body=clean_body,
        has_throw=has_throw,
        aux_defs=ctx.aux_defs,
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

def _detect_iterator_compact(stmts: list[StmtIR], i: int) -> tuple[str, str, int] | None:
    """检测迭代器压缩模式（设计 §6.4）。

    模式：it = begin; out = begin; for(...) { transform; if(pred) compact }; erase
    返回 (arr_var_name, mod_var_name, 消耗的语句数) 或 None。
    """
    if i + 3 > len(stmts):
        return None

    # stmt[i]: LetStmt(it, _, Call("*.toList", [f]))
    s0 = stmts[i]
    if not (isinstance(s0, LetStmt) and isinstance(s0.value, Call)):
        return None
    if "toList" not in s0.value.func and "begin" not in s0.value.func:
        return None

    # stmt[i+1]: LetStmt(out, _, ...) — out 初始化（可以是 it 的拷贝或另一个 begin）
    s1 = stmts[i + 1] if i + 1 < len(stmts) else None
    if not isinstance(s1, LetStmt):
        return None

    # stmt[i+2]: ForLoop 或 WhileLoop（循环体含 operator->）
    s2 = stmts[i + 2] if i + 2 < len(stmts) else None
    if not isinstance(s2, UnknownStmt):
        return None
    if s2.kind not in ("ForLoop",):
        return None

    # 检测循环体是否含 operator->（迭代器解引用）
    def has_operator_arrow(stmt):
        if isinstance(stmt, LetStmt) and isinstance(stmt.value, FieldAccess):
            return True
        if isinstance(stmt, ExprStmt) and isinstance(stmt.expr, Call):
            if "operator" in stmt.expr.func:
                return True
        if isinstance(stmt, UnknownStmt):
            return any(has_operator_arrow(c) for c in stmt.children)
        if isinstance(stmt, IfStmt):
            return any(has_operator_arrow(c) for c in stmt.then_body + stmt.else_body)
        return False

    # 找数组变量名
    arr_name = None
    if s0.value.args:
        arg = s0.value.args[0]
        if isinstance(arg, Var):
            arr_name = arg.name

    if arr_name is None:
        return None

    # stmt[i+3]: erase 调用（ExprStmt 或 LetStmt 含 pop/erase）
    consumed = 3
    if i + 3 < len(stmts):
        s3 = stmts[i + 3]
        if isinstance(s3, ExprStmt) and isinstance(s3.expr, Call):
            if "pop" in s3.expr.func or "erase" in s3.expr.func:
                consumed = 4
        elif isinstance(s3, LetStmt) and isinstance(s3.value, Call):
            if "pop" in s3.value.func or "erase" in s3.value.func:
                consumed = 4

    return (arr_name, s0.var.name, consumed)


def transform_body(stmts: list[StmtIR], env: VarEnv, ctx: FuncCtx = None,
                    loop_ctx: LoopCtx = None) -> list[StmtIR]:
    if ctx is None:
        ctx = FuncCtx()
    result = []
    i = 0
    while i < len(stmts):
        stmt = stmts[i]

        # --- G4: 迭代器压缩模式 → filterMap / filter ---
        compact = _detect_iterator_compact(stmts, i)
        if compact:
            arr_name, _, consumed = compact
            arr_var = env.current(arr_name) if arr_name in env.versions else Var(arr_name)
            new_var = env.bump(arr_name)
            typ = env.get_type(arr_name)
            # 检测是 modCoeff（有 fdiv_r）还是纯 compact（无 transform）
            # 简化：如果 stmts[i+2] 循环体含 fdiv_r 调用 → modCoeff，否则 → compactNonzero
            loop_body_str = str(stmts[i + 2]) if i + 2 < len(stmts) else ""
            if "fdiv_r" in loop_body_str or "fdiv_q" in loop_body_str:
                # modCoeff 模式：需要 m 参数
                # m 通常是循环体中 fdiv_r 的最后一个参数——从函数参数中获取
                m_var = Var("m")  # 占位：精化证明时确定
                for p in env.versions:
                    if p == "m":
                        m_var = env.current("m")
                        break
                result.append(LetStmt(new_var, typ,
                              Call("SparsePolyZZ.modCoeff", [arr_var, m_var])))
            else:
                result.append(LetStmt(new_var, typ,
                              Call("SparsePolyZZ.compactNonzero", [arr_var])))
            i += consumed
            continue

        # --- break/continue 处理（仅在循环函数体内） ---
        if loop_ctx and isinstance(stmt, UnknownStmt):
            if stmt.kind == "BreakStmt":
                # break → return 当前状态
                exit_expr = _make_loop_return(loop_ctx.modified_vars, env)
                result.append(ReturnStmt(exit_expr))
                return result  # break 后的代码不可达
            if stmt.kind == "ContinueStmt":
                # continue → 先执行 step（for-loop），再递归
                if loop_ctx.step_stmts:
                    result.extend(loop_ctx.step_stmts)
                recurse = _make_loop_recurse(loop_ctx, env)
                result.append(ReturnStmt(recurse))
                return result  # continue 后的代码不可达

        # --- 检测 early return/break/continue 模式 ---
        # if (cond) return/break/continue; → if cond then <action> else {rest}
        if isinstance(stmt, IfStmt) and not stmt.else_body:
            then_has_exit = any(
                isinstance(s, ReturnStmt) or
                (isinstance(s, UnknownStmt) and s.kind in ("BreakStmt", "ContinueStmt"))
                for s in stmt.then_body
            )
            if then_has_exit and i + 1 < len(stmts):
                cond = rename_expr(stmt.cond, env)
                env_then = env.fork()
                then_body = transform_body(stmt.then_body, env_then, ctx, loop_ctx)
                env_else = env.fork()
                else_body = transform_body(stmts[i + 1:], env_else, ctx, loop_ctx)
                result.append(IfStmt(cond, then_body, else_body))
                env.versions = env_else.versions
                env.types.update(env_else.types)
                return result

        transformed = transform_stmt(stmt, env, ctx, loop_ctx)
        if isinstance(transformed, list):
            result.extend(transformed)
        else:
            result.append(transformed)
        i += 1
    return result


def _make_loop_return(modified_vars: list[str], env: VarEnv) -> ExprIR:
    """构造循环退出返回值：当前修改变量的值。"""
    if len(modified_vars) == 0:
        return Lit(0)
    elif len(modified_vars) == 1:
        return env.current(modified_vars[0])
    else:
        return Call("Prod.mk", [env.current(n) for n in modified_vars])


def _make_loop_recurse(loop_ctx: LoopCtx, env: VarEnv) -> ExprIR:
    """构造 continue 的递归调用。传全部循环参数。"""
    args = []
    for name in loop_ctx.all_param_names:
        if name in env.versions:
            args.append(env.current(name))
        else:
            args.append(Var(name))
    # finalize_loop_func 会补充遗漏的闭包变量
    return Call(loop_ctx.func_name + "_ir", args)


def transform_stmt(stmt: StmtIR, env: VarEnv, ctx: FuncCtx = None,
                    loop_ctx: LoopCtx = None) -> StmtIR | list[StmtIR]:
    if ctx is None:
        ctx = FuncCtx()
    if isinstance(stmt, LetStmt):
        # Lambda 辅助函数创建（设计 §4.4）
        if isinstance(stmt.value, Call) and stmt.value.func.startswith("_lambda_"):
            from clang_ast import _LAMBDA_REGISTRY
            lambda_id_str = stmt.value.func.replace("_lambda_", "")
            try:
                lambda_id = int(lambda_id_str)
            except ValueError:
                lambda_id = None
            if lambda_id is not None and lambda_id in _LAMBDA_REGISTRY:
                captures_from_ast, lam_params, lam_body = _LAMBDA_REGISTRY[lambda_id]
                func_name = f"_lambda_{lambda_id}"
                # 用 clang_ast 提供的 capture 列表作为 primary source
                # 再用 collect_free_vars 补充遗漏
                param_names = {lp.name for lp in lam_params}
                capture_names = []
                seen = set()
                for cname in captures_from_ast:
                    if cname in param_names or cname in seen:
                        continue
                    if cname in env.versions:
                        seen.add(cname)
                        capture_names.append(cname)
                # fallback: collect_free_vars 补充 AST 未捕获的变量
                raw_free = collect_free_vars(lam_body)
                for vname, vver in sorted(raw_free):
                    if vname in param_names or vname in seen:
                        continue
                    if vname in env.versions:
                        seen.add(vname)
                        capture_names.append(vname)
                # SSA 变换 lambda body
                lam_env = env.fork()
                for lp in lam_params:
                    lam_env.versions[lp.name] = 0
                    lam_env.set_type(lp.name, lp.typ)
                lam_body_ssa = transform_body(lam_body, lam_env, ctx)
                # 构造参数列表：capture + formal params
                capture_params = []
                for cn in capture_names:
                    capture_params.append(ParamIR(env.current(cn).lean_name(),
                                                   env.get_type(cn)))
                all_params = capture_params + list(lam_params)
                # S4: 如果 lambda 修改了 capture，追加返回修改后的值
                # M4: 保守策略 — 所有 capture 都假设可能被修改
                # 未修改的 capture 写回 = noop，语义安全
                modified_caps = list(capture_names)
                if modified_caps and not any(isinstance(s, ReturnStmt) for s in lam_body_ssa):
                    if len(modified_caps) == 1:
                        lam_body_ssa.append(ReturnStmt(lam_env.current(modified_caps[0])))
                    else:
                        lam_body_ssa.append(ReturnStmt(
                            Call("Prod.mk", [lam_env.current(c) for c in modified_caps])))
                # 记录 lambda 变量名 → 修改的 capture 列表
                if modified_caps:
                    env.lambda_modified[stmt.var.name] = modified_caps
                lam_func = SSAFunc(
                    name=func_name,
                    params=all_params,
                    ret_type="auto",
                    requires=[],
                    body=lam_body_ssa,
                )
                finalize_loop_func(lam_func, env, [])
                ctx.aux_defs.append(lam_func)

        # N1: LetStmt 中的 output-param 函数调用 → 解构返回值
        # 解包 Cast 层（type annotation）和 BinOp("=") 层（引用参数赋值）
        call_value = stmt.value
        if isinstance(call_value, Cast):
            call_value = call_value.expr
        # BinOp("=", Var("f"), Call("make_monic", [f])) → 解包到 Call
        if isinstance(call_value, BinOp) and call_value.op == "=" and isinstance(call_value.rhs, (Call, Cast)):
            call_value = call_value.rhs
            if isinstance(call_value, Cast):
                call_value = call_value.expr
        if isinstance(call_value, Call):
            func_base = call_value.func.replace("_ir", "") if isinstance(call_value.func, str) and call_value.func.endswith("_ir") else (call_value.func if isinstance(call_value.func, str) else "")
            from class_map import TRANSLATION_SCOPE_OUTPUT_PARAMS
            if func_base in TRANSLATION_SCOPE_OUTPUT_PARAMS:
                out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[func_base]
                out_vars_orig = [call_value.args[i] for i in out_indices if i < len(call_value.args)]
                new_value = rename_expr(call_value, env)
                stmts = []
                # _out := call(...)
                pair_var = env.bump("_out")
                env.set_type("_out", "auto")
                stmts.append(LetStmt(pair_var, "auto", new_value))
                # ret_val := _out.1 (原返回值，绑定到 stmt.var)
                ret_var = env.bump(stmt.var.name) if stmt.var.name != "_" else Var("_")
                env.set_type(stmt.var.name, stmt.typ)
                stmts.append(LetStmt(ret_var, stmt.typ, FieldAccess(pair_var, "1")))
                # 输出参数 := _out.2, _out.2.2, ...
                for i, ov in enumerate(out_vars_orig):
                    if isinstance(ov, Var):
                        new_var = env.bump(ov.name)
                        if i == 0 and len(out_vars_orig) == 1:
                            accessor = FieldAccess(pair_var, "2")
                        elif i == 0:
                            accessor = FieldAccess(pair_var, "2")
                            accessor = FieldAccess(accessor, "1")
                        elif i == len(out_vars_orig) - 1:
                            accessor = FieldAccess(pair_var, "2")
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                        else:
                            accessor = FieldAccess(pair_var, "2")
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                            accessor = FieldAccess(accessor, "1")
                        stmts.append(LetStmt(new_var, env.get_type(ov.name), accessor))
                return stmts

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
        return transform_if(stmt, env, ctx, loop_ctx)

    if isinstance(stmt, ReturnStmt):
        # 循环体内的 ReturnStmt：转换为 _ret_flag=true + _ret_val=val + break
        if loop_ctx and loop_ctx.has_parent_return and stmt.value is not None:
            stmts = []
            ret_val = rename_expr(stmt.value, env)
            flag_var = env.bump("_ret_flag")
            stmts.append(LetStmt(flag_var, BaseType.BOOL, Lit(True)))
            val_var = env.bump("_ret_val")
            stmts.append(LetStmt(val_var, "auto", ret_val))
            # break → return modified_vars tuple
            exit_expr = _make_loop_return(loop_ctx.modified_vars, env)
            stmts.append(ReturnStmt(exit_expr))
            return stmts
        if stmt.value is not None:
            return ReturnStmt(rename_expr(stmt.value, env))
        return stmt

    if isinstance(stmt, Require):
        return Require(rename_expr(stmt.cond, env), stmt.name, stmt.source)

    if isinstance(stmt, Throw):
        return stmt

    if isinstance(stmt, ExprStmt):
        expr = stmt.expr
        # 解包 Cast 层（type annotation 不影响语义）
        while isinstance(expr, Cast):
            expr = expr.expr
        # M2：TRANSLATION_SCOPE 函数调用的输出参数捕获
        if isinstance(expr, Call):
            func_base = expr.func.replace("_ir", "") if expr.func.endswith("_ir") else expr.func
            from class_map import TRANSLATION_SCOPE_OUTPUT_PARAMS
            if func_base in TRANSLATION_SCOPE_OUTPUT_PARAMS:
                out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[func_base]
                all_args = [rename_expr(a, env) for a in expr.args]
                call = Call(expr.func, all_args)
                out_vars = [expr.args[i] for i in out_indices if i < len(expr.args)]
                # 单输出
                if len(out_vars) == 1 and isinstance(out_vars[0], Var):
                    new_var = env.bump(out_vars[0].name)
                    return LetStmt(new_var, env.get_type(out_vars[0].name), call)
                # 多输出
                if len(out_vars) > 1:
                    stmts = []
                    pair_var = env.bump("_out")
                    env.set_type("_out", "auto")
                    stmts.append(LetStmt(pair_var, "auto", call))
                    for i, ov in enumerate(out_vars):
                        if isinstance(ov, Var):
                            new_var = env.bump(ov.name)
                            if i == 0:
                                accessor = FieldAccess(pair_var, "1")
                            elif i == len(out_vars) - 1:
                                accessor = pair_var
                                for _ in range(i):
                                    accessor = FieldAccess(accessor, "2")
                            else:
                                accessor = pair_var
                                for _ in range(i):
                                    accessor = FieldAccess(accessor, "2")
                                accessor = FieldAccess(accessor, "1")
                            stmts.append(LetStmt(new_var, env.get_type(ov.name), accessor))
                    return stmts
        # i++ / ++i → let i_{n+1} := i_n + 1
        if isinstance(expr, UnaryOp) and expr.op in ("++", "--"):
            delta = Lit(1) if expr.op == "++" else Lit(-1)
            op = "+" if expr.op == "++" else "-"
            if isinstance(expr.operand, Var):
                old = env.current(expr.operand.name)
                new_var = env.bump(expr.operand.name)
                typ = env.get_type(expr.operand.name)
                return LetStmt(new_var, typ, BinOp(op, old, Lit(1)))
            # arr[i]++ → Array.set! arr i (arr[i] + 1)
            if isinstance(expr.operand, ArrayAccess) and isinstance(expr.operand.arr, Var):
                arr_var = env.current(expr.operand.arr.name)
                idx = rename_expr(expr.operand.idx, env)
                old_val = ArrayAccess(arr_var, idx)
                new_arr = env.bump(expr.operand.arr.name)
                return LetStmt(new_arr, env.get_type(expr.operand.arr.name),
                               Call("Array.set!", [arr_var, idx, BinOp(op, old_val, Lit(1))]))
        # _mutate_* 方法（来自 CLASS_MAP 的 mutate 类别）
        # _mutate_X(obj) → let obj_{n+1} := X obj_n
        if isinstance(expr, Call) and isinstance(expr.func, str) and expr.func.startswith("_mutate_") and expr.args:
            lean_name = expr.func[len("_mutate_"):]
            if isinstance(expr.args[0], Var):
                old = env.current(expr.args[0].name)
                new_var = env.bump(expr.args[0].name)
                typ = env.get_type(expr.args[0].name)
                return LetStmt(new_var, typ, Call(lean_name, [old]))
            # _mutate_X(arr[i]) → noop（normalization 等对数组元素原地操作，
            # 在 Lean 模型中 normalization 是 noop）
            if lean_name in ("MvPolyZp.normalization", "MvPolyZZ.normalization",
                             "MvMonomial.normalization", "SparsePolyZp.normalization",
                             "SparsePolyZZ.normalization", "Array.empty"):
                return ExprStmt(rename_expr(expr.args[0], env))
        # compound assignment: BinOp("+"/"-"/"*", ArrayAccess/Var, rhs) as ExprStmt
        # = operator-=, +=, *= 翻译结果（AST 只保留运算，丢失赋值）
        # 还原为：lhs = lhs op rhs
        if isinstance(expr, BinOp) and expr.op in ("+", "-", "*", "/", "%"):
            lhs = expr.lhs
            # 1D: arr[i] -= val → arr[i] = arr[i] - val → Array.set!
            if isinstance(lhs, ArrayAccess) and isinstance(lhs.arr, Var):
                arr_var = env.current(lhs.arr.name)
                idx = rename_expr(lhs.idx, env)
                rhs = rename_expr(expr.rhs, env)
                old_val = ArrayAccess(arr_var, idx)
                new_val = BinOp(expr.op, old_val, rhs)
                new_arr = env.bump(lhs.arr.name)
                return LetStmt(new_arr, env.get_type(lhs.arr.name),
                               Call("Array.set!", [arr_var, idx, new_val]))
            # 2D: arr[i][k] -= val → 取行、改元素、放回
            if isinstance(lhs, ArrayAccess) and isinstance(lhs.arr, ArrayAccess):
                outer = lhs.arr
                if isinstance(outer.arr, Var):
                    arr_var = env.current(outer.arr.name)
                    i_idx = rename_expr(outer.idx, env)
                    k_idx = rename_expr(lhs.idx, env)
                    rhs = rename_expr(expr.rhs, env)
                    row = ArrayAccess(arr_var, i_idx)
                    old_val = ArrayAccess(row, k_idx)
                    new_val = BinOp(expr.op, old_val, rhs)
                    row_var = Var("_row")
                    row_upd = Var("_row_upd")
                    stmts = [
                        LetStmt(row_var, "auto", row),
                        LetStmt(row_upd, "auto", Call("Array.set!", [row_var, k_idx, new_val])),
                    ]
                    new_arr = env.bump(outer.arr.name)
                    stmts.append(LetStmt(new_arr, env.get_type(outer.arr.name),
                                         Call("Array.set!", [arr_var, i_idx, row_upd])))
                    return stmts
            # x -= val → x = x - val
            if isinstance(lhs, Var):
                old = env.current(lhs.name)
                rhs = rename_expr(expr.rhs, env)
                new_var = env.bump(lhs.name)
                return LetStmt(new_var, env.get_type(lhs.name), BinOp(expr.op, old, rhs))
        # arr.push(x) → let arr_{n+1} := arr_n.push(x)
        if isinstance(expr, ArrayPush) and isinstance(expr.arr, Var):
            old = env.current(expr.arr.name)
            new_var = env.bump(expr.arr.name)
            typ = env.get_type(expr.arr.name)
            return LetStmt(new_var, typ, ArrayPush(old, rename_expr(expr.elem, env)))
        # obj.field.push(x) → let obj_{n+1} := obj_n.push x（Factorization 等容器字段）
        if isinstance(expr, ArrayPush) and isinstance(expr.arr, FieldAccess):
            fa = expr.arr
            if isinstance(fa.obj, Var):
                old = env.current(fa.obj.name)
                new_var = env.bump(fa.obj.name)
                typ = env.get_type(fa.obj.name)
                elem = rename_expr(expr.elem, env)
                return LetStmt(new_var, typ, ArrayPush(old, elem))
        # arr[i].push(x) → let arr_{n+1} := Array.set! arr_n i (arr_n[i]!.push x)
        if isinstance(expr, ArrayPush) and isinstance(expr.arr, ArrayAccess):
            aa = expr.arr
            if isinstance(aa.arr, Var):
                arr_old = env.current(aa.arr.name)
                idx = rename_expr(aa.idx, env)
                elem = rename_expr(expr.elem, env)
                row = ArrayAccess(arr_old, idx)
                pushed = ArrayPush(row, elem)
                new_arr = env.bump(aa.arr.name)
                return LetStmt(new_arr, env.get_type(aa.arr.name),
                               Call("Array.set!", [arr_old, idx, pushed]))
        # x = expr（变量赋值）→ let x_{n+1} := expr
        # ":=" 是迭代器赋值（it->second = val），语义与 "=" 相同
        if isinstance(expr, BinOp) and expr.op in ("=", ":="):
            if isinstance(expr.lhs, Var):
                # N1: x = f(args) 且 f 有 output params → 解构
                rhs = expr.rhs
                if isinstance(rhs, Cast):
                    rhs = rhs.expr
                if isinstance(rhs, Call) and isinstance(rhs.func, str):
                    func_base = rhs.func.replace("_ir", "") if rhs.func.endswith("_ir") else rhs.func
                    from class_map import TRANSLATION_SCOPE_OUTPUT_PARAMS
                    if func_base in TRANSLATION_SCOPE_OUTPUT_PARAMS:
                        out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[func_base]
                        out_vars_orig = [rhs.args[i] for i in out_indices if i < len(rhs.args)]
                        call_renamed = rename_expr(rhs, env)
                        stmts = []
                        pair_var = env.bump("_out")
                        env.set_type("_out", "auto")
                        stmts.append(LetStmt(pair_var, "auto", call_renamed))
                        # x := _out.1 (原返回值)
                        new_var = env.bump(expr.lhs.name)
                        stmts.append(LetStmt(new_var, env.get_type(expr.lhs.name),
                                             FieldAccess(pair_var, "1")))
                        # 输出参数 := _out.2
                        for i, ov in enumerate(out_vars_orig):
                            if isinstance(ov, Var):
                                new_ov = env.bump(ov.name)
                                accessor = FieldAccess(pair_var, "2")
                                if len(out_vars_orig) > 1 and i < len(out_vars_orig) - 1:
                                    for _ in range(i):
                                        accessor = FieldAccess(accessor, "2")
                                    accessor = FieldAccess(accessor, "1")
                                elif len(out_vars_orig) > 1:
                                    for _ in range(i):
                                        accessor = FieldAccess(accessor, "2")
                                stmts.append(LetStmt(new_ov, env.get_type(ov.name), accessor))
                        return stmts
                new_value = rename_expr(expr.rhs, env)
                typ = env.get_type(expr.lhs.name)
                new_var = env.bump(expr.lhs.name)
                return LetStmt(new_var, typ, new_value)
            # 多输出赋值：Prod.mk(out₁, out₂) = call → 解构
            if isinstance(expr.lhs, Call) and expr.lhs.func == "Prod.mk":
                out_vars = expr.lhs.args  # [Var("q"), Var("r")]
                new_value = rename_expr(expr.rhs, env)
                stmts = []
                if len(out_vars) == 1 and isinstance(out_vars[0], Var):
                    new_var = env.bump(out_vars[0].name)
                    return LetStmt(new_var, env.get_type(out_vars[0].name), new_value)
                # 多输出：let _result := call; let q' := _result.1; let r' := _result.2
                pair_var = Var("_out")
                pair_ver = env.bump("_out")
                env.set_type("_out", "auto")
                stmts.append(LetStmt(pair_ver, "auto", new_value))
                for i, ov in enumerate(out_vars):
                    if isinstance(ov, Var):
                        new_var = env.bump(ov.name)
                        # Prod 嵌套：A × (B × C) → .1, .2.1, .2.2, ...
                        if i == 0:
                            accessor = FieldAccess(pair_ver, "1")
                        elif i == len(out_vars) - 1:
                            accessor = pair_ver
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                        else:
                            accessor = pair_ver
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                            accessor = FieldAccess(accessor, "1")
                        stmts.append(LetStmt(new_var, env.get_type(ov.name), accessor))
                return stmts
            # M4：1D 数组赋值 arr[i] = val → Array.set!
            if isinstance(expr.lhs, ArrayAccess) and isinstance(expr.lhs.arr, Var):
                arr_var = env.current(expr.lhs.arr.name)
                idx = rename_expr(expr.lhs.idx, env)
                val = rename_expr(expr.rhs, env)
                new_arr = env.bump(expr.lhs.arr.name)
                return LetStmt(new_arr, env.get_type(expr.lhs.arr.name),
                               Call("Array.set!", [arr_var, idx, val]))
            # M4：2D 数组赋值 arr[i][j] = val → 取行、改行、放回
            if isinstance(expr.lhs, ArrayAccess) and isinstance(expr.lhs.arr, ArrayAccess):
                outer = expr.lhs.arr  # arr[i]
                if isinstance(outer.arr, Var):
                    arr_var = env.current(outer.arr.name)
                    i_idx = rename_expr(outer.idx, env)
                    j_idx = rename_expr(expr.lhs.idx, env)
                    val = rename_expr(expr.rhs, env)
                    row_var = Var("_row")
                    row_updated = Var("_row_upd")
                    stmts = [
                        LetStmt(row_var, "auto", ArrayAccess(arr_var, i_idx)),
                        LetStmt(row_updated, "auto",
                                Call("Array.set!", [row_var, j_idx, val])),
                    ]
                    new_arr = env.bump(outer.arr.name)
                    stmts.append(LetStmt(new_arr, env.get_type(outer.arr.name),
                                 Call("Array.set!", [arr_var, i_idx, row_updated])))
                    return stmts
            if isinstance(expr.lhs, FieldAccess):
                return transform_member_assign(expr.lhs, expr.rhs, env)
            # 解包 Cast 层再检查 FieldAccess（迭代器 it->second = val）
            lhs_unwrap = expr.lhs
            while isinstance(lhs_unwrap, Cast):
                lhs_unwrap = lhs_unwrap.expr
            if isinstance(lhs_unwrap, FieldAccess):
                return transform_member_assign(lhs_unwrap, expr.rhs, env)
            # catch-all: 未识别的 BinOp(=) lhs → 尝试提取根变量
            root_vars = _extract_all_root_var_names(expr.lhs)
            if root_vars:
                # 有可识别的根变量 → 当作变量赋值
                var_name = sorted(root_vars)[0]
                new_value = rename_expr(expr.rhs, env)
                new_var = env.bump(var_name)
                return LetStmt(new_var, env.get_type(var_name), new_value)
        # S4: Lambda 调用结果写回 — lambda 修改了 capture 变量
        # 检测：Call 的函数名是一个 Var，且该 Var 在 env.lambda_modified 中
        func_name_raw = ""
        if isinstance(expr, Call):
            func_name_raw = expr.func if isinstance(expr.func, str) else \
                            (expr.func.name if isinstance(expr.func, Var) else "")
        if func_name_raw in env.lambda_modified:
            modified_captures = env.lambda_modified[func_name_raw]
            call_renamed = rename_expr(expr, env)
            stmts = []
            if len(modified_captures) == 1:
                new_var = env.bump(modified_captures[0])
                stmts.append(LetStmt(new_var, env.get_type(modified_captures[0]),
                                     call_renamed))
            else:
                pair_var = env.bump("_lam_out")
                env.set_type("_lam_out", "auto")
                stmts.append(LetStmt(pair_var, "auto", call_renamed))
                for i, cname in enumerate(modified_captures):
                    new_var = env.bump(cname)
                    if i == 0:
                        accessor = FieldAccess(pair_var, "1")
                    elif i == len(modified_captures) - 1:
                        accessor = pair_var
                        for _ in range(i):
                            accessor = FieldAccess(accessor, "2")
                    else:
                        accessor = pair_var
                        for _ in range(i):
                            accessor = FieldAccess(accessor, "2")
                        accessor = FieldAccess(accessor, "1")
                    stmts.append(LetStmt(new_var, env.get_type(cname), accessor))
            return stmts
        # 通用 Call fallback：未被上面 handler 匹配的函数调用
        # 如果在 lambda_modified 中 → M4 已处理（上面）
        # 否则：丢弃返回值但保留调用（可能有副作用通过其他机制传播）
        if isinstance(expr, Call) and isinstance(expr.func, str):
            func = expr.func
            # pair_vec_div 在 ExprStmt 位置 = 输出参数已由 M2 handler 处理
            # next_combination / mv_next_combination / normalize_factor = 副作用通过循环机制传播
            if func in ("pair_vec_div", "next_combination", "mv_next_combination",
                        "normalize_factor", "all_div", "lift_vars"):
                return ExprStmt(rename_expr(expr, env))  # 保留调用，副作用由其他机制处理
        # BinOp(=) lhs=Call → operator[] 赋值（如 map[key] = val）
        # 排除 Prod.mk（已由 line 876 处理，如果到这里说明是嵌套场景）
        if isinstance(expr, BinOp) and expr.op == "=" and isinstance(expr.lhs, Call):
            if isinstance(expr.lhs.func, str) and expr.lhs.func != "Prod.mk":
                return ExprStmt(rename_expr(expr, env))
            # Prod.mk 赋值 fallback：解构
            if isinstance(expr.lhs.func, str) and expr.lhs.func == "Prod.mk":
                out_vars = expr.lhs.args
                new_value = rename_expr(expr.rhs, env)
                stmts = []
                pair_var = env.bump("_out")
                env.set_type("_out", "auto")
                stmts.append(LetStmt(pair_var, "auto", new_value))
                for i, ov in enumerate(out_vars):
                    if isinstance(ov, Var):
                        new_var = env.bump(ov.name)
                        if i == 0:
                            accessor = FieldAccess(pair_var, "1")
                        elif i == len(out_vars) - 1:
                            accessor = pair_var
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                        else:
                            accessor = pair_var
                            for _ in range(i):
                                accessor = FieldAccess(accessor, "2")
                            accessor = FieldAccess(accessor, "1")
                        stmts.append(LetStmt(new_var, env.get_type(ov.name), accessor))
                return stmts if stmts else ExprStmt(rename_expr(expr, env))
        # === 黑名单策略：分类 ExprStmt，拒绝静默丢弃未知变异 ===
        renamed = rename_expr(expr, env)
        if _is_known_noop(expr):
            return ExprStmt(renamed)
        # 未识别的 ExprStmt → 报警（不再静默 let _ :=）
        import sys
        func_desc = ""
        if isinstance(expr, Call):
            func_desc = expr.func if isinstance(expr.func, str) else str(expr.func)
        elif isinstance(expr, BinOp):
            func_desc = f"BinOp({expr.op})"
        else:
            func_desc = type(expr).__name__
        print(f"WARNING: unhandled ExprStmt (possible mutation lost): {func_desc}",
              file=sys.stderr)
        return ExprStmt(renamed)

    if isinstance(stmt, UnknownStmt):
        if stmt.kind == "ForLoop":
            return transform_for_loop(stmt, env, ctx)
        if stmt.kind == "WhileLoop":
            return transform_while_loop(stmt, env, ctx)
        if stmt.kind == "RangeForLoop":
            return transform_range_for(stmt, env, ctx)
        if stmt.kind == "CompoundStmt":
            return transform_body(stmt.children, env, ctx, loop_ctx)
        if stmt.kind == "MultiDecl":
            return transform_body(stmt.children, env, ctx, loop_ctx)
        if stmt.kind == "BreakStmt":
            return stmt  # 由 transform_for_loop 处理
        if stmt.kind == "ContinueStmt":
            return stmt  # 由 transform_for_loop 处理
        return stmt

    return stmt


# ============================================================
# if/else 变换（phi 节点）
# ============================================================

def transform_if(stmt: IfStmt, env: VarEnv, ctx: FuncCtx = None,
                  loop_ctx: LoopCtx = None) -> list[StmtIR]:
    if ctx is None:
        ctx = FuncCtx()
    cond = rename_expr(stmt.cond, env)

    # fork 两个分支
    env_then = env.fork()
    env_else = env.fork()

    then_body = transform_body(stmt.then_body, env_then, ctx, loop_ctx)
    else_body = transform_body(stmt.else_body, env_else, ctx, loop_ctx) if stmt.else_body else []

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

    # 两分支都不 return → 需要 phi（或同步版本）
    diff = env_then.merge(env_else)

    # 只处理在 if 之前已存在的变量（排除 if 内首次声明的局部变量）
    pre_existing = set(env.versions.keys())
    diff = {k: v for k, v in diff.items() if k in pre_existing}

    # 收集两分支一致推进的变量（版本相同但高于外层）
    agreed = {}
    for name in pre_existing:
        v_then = env_then.versions.get(name, 0)
        v_else = env_else.versions.get(name, 0)
        v_outer = env.versions.get(name, 0)
        if v_then == v_else and v_then > v_outer:
            agreed[name] = (Var(name, v_then), Var(name, v_then))

    if not diff and not agreed:
        return [IfStmt(cond, then_body, else_body)]

    # 有变量分歧 → 用 BlockExpr 封装分支内的 let 链 + 最终值
    # 这确保分支内的中间变量在 Lean 的 if-then 作用域内可见
    result: list[StmtIR] = []

    # 收集分支中与特定变量相关的所有 let 语句链 + 最终值
    def _extract_stmts_for_var(stmts: list[StmtIR], var_name: str) -> tuple[list[StmtIR], ExprIR]:
        """提取分支中与 var_name 相关的所有 let 语句链（含跨变量依赖）+ 最终值。"""
        # 建索引：lean_name → LetStmt
        stmt_by_name = {}
        for s in stmts:
            if isinstance(s, LetStmt):
                stmt_by_name[s.var.lean_name()] = s

        # 找该变量的最后赋值
        final_stmt = None
        for s in reversed(stmts):
            if isinstance(s, LetStmt) and s.var.name == var_name:
                final_stmt = s
                break
        if not final_stmt:
            return [], None

        # 从最终赋值反向追踪所有依赖（用 lean_name 精确匹配）
        needed = set()
        queue = [final_stmt.var.lean_name()]
        while queue:
            vn = queue.pop()
            if vn in needed:
                continue
            needed.add(vn)
            s = stmt_by_name.get(vn)
            if s:
                for dep_lean_name in _extract_lean_names(s.value):
                    if dep_lean_name in stmt_by_name and dep_lean_name not in needed:
                        queue.append(dep_lean_name)

        # 按原始顺序输出
        relevant = []
        for s in stmts:
            if isinstance(s, LetStmt) and s.var.lean_name() in needed:
                relevant.append(s)

        return relevant, final_stmt.var

    # 合并 diff（版本不同）和 agreed（版本一致但需 phi）
    all_phi = {}
    all_phi.update(diff)
    all_phi.update(agreed)
    for name in sorted(all_phi.keys()):
        then_var, else_var = all_phi[name]
        # then 分支：提取相关 let 链
        then_stmts, then_final = _extract_stmts_for_var(then_body, name)
        # else 分支：提取相关 let 链
        else_stmts, else_final = _extract_stmts_for_var(else_body, name)

        then_expr = BlockExpr(then_stmts, then_final) if then_stmts and then_final else then_var
        else_expr = BlockExpr(else_stmts, else_final) if else_stmts and else_final else else_var

        new_var = env.bump(name)
        typ = env.get_type(name)
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

def transform_for_loop(stmt: UnknownStmt, env: VarEnv, ctx: FuncCtx = None) -> list[StmtIR]:
    """ForLoop → init_stmts + 提取为独立 partial def + 调用。

    stmt.children = [init, cond_expr, step, body]
    """
    if ctx is None:
        ctx = FuncCtx()
    children = stmt.children
    if len(children) < 4:
        return [stmt]

    init_stmt = children[0]
    cond_stmt = children[1]
    step_stmt = children[2]
    body_stmt = children[3]

    # 条件表达式
    cond_expr = Lit(1, BaseType.BOOL)
    if isinstance(cond_stmt, ExprStmt):
        cond_expr = cond_stmt.expr

    # 识别循环变量
    loop_vars = identify_loop_vars(body_stmt) - _collect_local_decls(body_stmt)

    # init（在 env 中注册）
    init_stmts = []
    if isinstance(init_stmt, LetStmt):
        transformed_init = transform_stmt(init_stmt, env, ctx)
        init_stmts.append(transformed_init) if not isinstance(transformed_init, list) else init_stmts.extend(transformed_init)
        loop_vars.add(init_stmt.var.name)
    elif isinstance(init_stmt, AssignStmt):
        transformed_init = transform_stmt(init_stmt, env, ctx)
        init_stmts.append(transformed_init) if not isinstance(transformed_init, list) else init_stmts.extend(transformed_init)
        loop_vars.add(init_stmt.target.name)

    # step 变量
    step_vars = identify_loop_vars_from_step(step_stmt)
    loop_vars |= step_vars
    modified_names = sorted(loop_vars)

    loop_id = abs(hash(str(stmt))) % 100000
    func_name = f"loop_{loop_id}"

    # S5: 检测循环体内是否有父函数 return
    has_parent_ret = _body_has_return(body_stmt)
    if has_parent_ret:
        # 注入 _ret_flag, _ret_val 作为循环状态
        if "_ret_flag" not in modified_names:
            modified_names.append("_ret_flag")
        if "_ret_val" not in modified_names:
            modified_names.append("_ret_val")
        modified_names = sorted(modified_names)
        # 在 env 中初始化
        env.versions["_ret_flag"] = 0
        env.set_type("_ret_flag", BaseType.BOOL)
        env.versions["_ret_val"] = 0
        env.set_type("_ret_val", "auto")

    # N2: 先用 fork env 计算 step（给 continue 用），再变换 body
    loop_env = env.fork()
    step_env = loop_env.fork()
    pre_step_stmts = transform_step(step_stmt, step_env)

    lctx = LoopCtx(func_name, modified_names, env=loop_env,
                    all_param_names=list(modified_names),
                    has_parent_return=has_parent_ret,
                    step_stmts=pre_step_stmts)  # continue 时可用
    body_stmts = transform_body(flatten_stmts(body_stmt), loop_env, ctx, loop_ctx=lctx)
    step_stmts = transform_step(step_stmt, loop_env)  # 正常路径的 step
    lctx.step_stmts = step_stmts  # 更新为 body 后版本

    # 闭包变量
    all_loop_stmts = body_stmts + step_stmts
    free = collect_free_vars(all_loop_stmts)
    loop_declared = set(modified_names)
    closure_vars = []
    seen_closure = set()
    for vname, vver in sorted(free):
        if vname in loop_declared or vname in seen_closure:
            continue
        if vname in env.versions:
            seen_closure.add(vname)
            closure_vars.append(env.current(vname))

    # 更新 LoopCtx 参数列表（闭包变量在 body 变换之后才可知）
    lctx.all_param_names = list(modified_names) + [cv.name for cv in closure_vars]

    # 循环函数参数
    loop_params = []
    for name in modified_names:
        loop_params.append(ParamIR(env.current(name).lean_name(), env.get_type(name)))
    for cv in closure_vars:
        loop_params.append(ParamIR(cv.lean_name(), env.get_type(cv.name)))

    # 退出条件
    cond_renamed = rename_expr(cond_expr, env)
    exit_cond = negate_expr(cond_renamed)

    # 退出返回值
    param_versions = {name: env.current(name) for name in modified_names}
    if len(modified_names) == 0:
        ret_type = BaseType.VOID
        exit_expr = Lit(0)
    elif len(modified_names) == 1:
        ret_type = env.get_type(modified_names[0])
        exit_expr = param_versions[modified_names[0]]
    else:
        ret_type = env.get_type(modified_names[0])
        exit_expr = Call("Prod.mk", [param_versions[n] for n in modified_names])

    # 递归调用（body + step 后的最新版本）
    recurse_args = []
    for name in modified_names:
        recurse_args.append(loop_env.current(name))
    for cv in closure_vars:
        recurse_args.append(cv)
    recurse_call = Call(func_name + "_ir", recurse_args)

    # 组装循环函数体
    loop_body_full = [
        IfStmt(exit_cond,
               [ReturnStmt(exit_expr)],
               body_stmts + step_stmts + [ReturnStmt(recurse_call)])
    ]

    loop_func = SSAFunc(
        name=func_name,
        params=loop_params,
        ret_type=ret_type,
        requires=[],
        body=loop_body_full,
    )
    ctx.aux_defs.append(loop_func)

    # 主函数中的调用
    call_args = []
    for name in modified_names:
        call_args.append(env.current(name))
    for cv in closure_vars:
        call_args.append(cv)
    call_expr = Call(func_name + "_ir", call_args)

    result_stmts = list(init_stmts)
    # S5: 初始化 _ret_flag/_ret_val
    if has_parent_ret:
        flag_init = env.bump("_ret_flag") if env.versions.get("_ret_flag", 0) == 0 else env.current("_ret_flag")
        if env.versions.get("_ret_flag", 0) <= 1:
            result_stmts.append(LetStmt(env.current("_ret_flag"), BaseType.BOOL, Lit(False)))
            result_stmts.append(LetStmt(env.current("_ret_val"), "auto", Lit(0)))
    if len(modified_names) == 0:
        result_stmts.append(ExprStmt(call_expr))
    elif len(modified_names) == 1:
        new_var = env.bump(modified_names[0])
        result_stmts.append(LetStmt(new_var, env.get_type(modified_names[0]), call_expr))
    else:
        pair_var = env.bump(f"_loop_{loop_id}")
        env.set_type(f"_loop_{loop_id}", "auto")
        result_stmts.append(LetStmt(pair_var, "auto", call_expr))
        for i, name in enumerate(modified_names):
            new_var = env.bump(name)
            if i == 0:
                accessor = FieldAccess(pair_var, "1")
            elif i == len(modified_names) - 1:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
            else:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
                accessor = FieldAccess(accessor, "1")
            result_stmts.append(LetStmt(new_var, env.get_type(name), accessor))

    finalize_loop_func(loop_func, env, result_stmts)

    # S5: 循环后检查 _ret_flag → 提前 return _ret_val
    if has_parent_ret:
        ret_flag = env.current("_ret_flag")
        ret_val = env.current("_ret_val")
        result_stmts.append(IfStmt(ret_flag, [ReturnStmt(ret_val)], []))

    return result_stmts


def identify_loop_vars_from_step(step_stmt: StmtIR) -> set[str]:
    """从步进语句（如 i++, p = next_p(p)）中提取变量。"""
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
        # p = next_p(p) 赋值形式
        if isinstance(expr, BinOp) and expr.op == "=":
            result |= _extract_all_root_var_names(expr.lhs)
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


def transform_while_loop(stmt: UnknownStmt, env: VarEnv, ctx: FuncCtx = None) -> list[StmtIR]:
    """WhileLoop → 提取为独立 partial def + 调用。"""
    if ctx is None:
        ctx = FuncCtx()
    children = stmt.children
    if len(children) < 2:
        return [stmt]

    cond_expr = Lit(1, BaseType.BOOL)
    if isinstance(children[0], ExprStmt):
        cond_expr = children[0].expr

    body_stmt = children[1]

    # 识别循环变量
    loop_vars = identify_loop_vars(body_stmt) - _collect_local_decls(body_stmt)
    modified_names = sorted(loop_vars)

    loop_id = abs(hash(str(stmt))) % 100000
    func_name = f"while_{loop_id}"

    # S5: 检测循环体内是否有父函数 return
    has_parent_ret = _body_has_return(body_stmt)
    if has_parent_ret:
        if "_ret_flag" not in modified_names:
            modified_names.append("_ret_flag")
        if "_ret_val" not in modified_names:
            modified_names.append("_ret_val")
        modified_names = sorted(modified_names)
        env.versions["_ret_flag"] = 0
        env.set_type("_ret_flag", BaseType.BOOL)
        env.versions["_ret_val"] = 0
        env.set_type("_ret_val", "auto")

    # 从循环体推断类型
    _infer_types_from_body(body_stmt, env)

    # SSA 变换循环体（传入 loop_ctx 支持 break/continue）
    loop_env = env.fork()
    lctx = LoopCtx(func_name, modified_names, env=loop_env,
                    all_param_names=list(modified_names),
                    has_parent_return=has_parent_ret)
    body_stmts = transform_body(flatten_stmts(body_stmt), loop_env, ctx, loop_ctx=lctx)

    # 闭包变量
    free = collect_free_vars(body_stmts)
    loop_declared = set(modified_names)
    closure_vars = []
    seen_closure = set()
    for vname, vver in sorted(free):
        if vname in loop_declared or vname in seen_closure:
            continue
        if vname in env.versions:
            seen_closure.add(vname)
            closure_vars.append(env.current(vname))

    # 更新 LoopCtx 参数列表（闭包变量在 body 变换之后才可知）
    lctx.all_param_names = list(modified_names) + [cv.name for cv in closure_vars]

    # 循环函数参数
    loop_params = []
    for name in modified_names:
        loop_params.append(ParamIR(env.current(name).lean_name(), env.get_type(name)))
    for cv in closure_vars:
        loop_params.append(ParamIR(cv.lean_name(), env.get_type(cv.name)))

    # 退出条件（while 条件取反）
    cond_renamed = rename_expr(cond_expr, env)
    exit_cond = negate_expr(cond_renamed)

    # 退出返回值
    param_versions = {name: env.current(name) for name in modified_names}
    if len(modified_names) == 0:
        ret_type = BaseType.VOID
        exit_expr = Lit(0)
    elif len(modified_names) == 1:
        ret_type = env.get_type(modified_names[0])
        exit_expr = param_versions[modified_names[0]]
    else:
        ret_type = env.get_type(modified_names[0])
        exit_expr = Call("Prod.mk", [param_versions[n] for n in modified_names])

    # 递归调用
    recurse_args = []
    for name in modified_names:
        recurse_args.append(loop_env.current(name))
    for cv in closure_vars:
        recurse_args.append(cv)
    recurse_call = Call(func_name + "_ir", recurse_args)

    # 组装循环函数体
    loop_body_full = [
        IfStmt(exit_cond,
               [ReturnStmt(exit_expr)],
               body_stmts + [ReturnStmt(recurse_call)])
    ]

    loop_func = SSAFunc(
        name=func_name,
        params=loop_params,
        ret_type=ret_type,
        requires=[],
        body=loop_body_full,
    )
    ctx.aux_defs.append(loop_func)

    # 主函数中的调用
    call_args = []
    for name in modified_names:
        call_args.append(env.current(name))
    for cv in closure_vars:
        call_args.append(cv)
    call_expr = Call(func_name + "_ir", call_args)

    result_stmts = []
    # M3: while-loop 的 _ret_flag/_ret_val 初始化
    if has_parent_ret:
        result_stmts.append(LetStmt(env.current("_ret_flag"), BaseType.BOOL, Lit(False)))
        result_stmts.append(LetStmt(env.current("_ret_val"), "auto", Lit(0)))
    if len(modified_names) == 0:
        result_stmts.append(ExprStmt(call_expr))
    elif len(modified_names) == 1:
        new_var = env.bump(modified_names[0])
        result_stmts.append(LetStmt(new_var, env.get_type(modified_names[0]), call_expr))
    else:
        pair_var = env.bump(f"_loop_{loop_id}")
        env.set_type(f"_loop_{loop_id}", "auto")
        result_stmts.append(LetStmt(pair_var, "auto", call_expr))
        for i, name in enumerate(modified_names):
            new_var = env.bump(name)
            if i == 0:
                accessor = FieldAccess(pair_var, "1")
            elif i == len(modified_names) - 1:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
            else:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
                accessor = FieldAccess(accessor, "1")
            result_stmts.append(LetStmt(new_var, env.get_type(name), accessor))

    finalize_loop_func(loop_func, env, result_stmts)

    # S5: 循环后检查 _ret_flag → 提前 return _ret_val
    if has_parent_ret:
        ret_flag = env.current("_ret_flag")
        ret_val = env.current("_ret_val")
        result_stmts.append(IfStmt(ret_flag, [ReturnStmt(ret_val)], []))

    return result_stmts


def collect_free_vars(stmts: list[StmtIR]) -> set[tuple[str, int]]:
    """收集语句列表中的自由变量（被引用但未在本作用域声明）。"""
    referenced: set[tuple[str, int]] = set()
    declared: set[tuple[str, int]] = set()

    def scan_expr(expr):
        if isinstance(expr, Var):
            referenced.add((expr.name, expr.version))
        elif isinstance(expr, BinOp):
            scan_expr(expr.lhs); scan_expr(expr.rhs)
        elif isinstance(expr, UnaryOp):
            scan_expr(expr.operand)
        elif isinstance(expr, Call):
            for a in expr.args:
                scan_expr(a)
        elif isinstance(expr, ArrayAccess):
            scan_expr(expr.arr); scan_expr(expr.idx)
        elif isinstance(expr, FieldAccess):
            scan_expr(expr.obj)
        elif isinstance(expr, ArrayPush):
            scan_expr(expr.arr); scan_expr(expr.elem)
        elif isinstance(expr, Cast):
            scan_expr(expr.expr)
        elif isinstance(expr, CondExpr):
            scan_expr(expr.cond); scan_expr(expr.then_e); scan_expr(expr.else_e)
        elif isinstance(expr, BlockExpr):
            for s in expr.stmts:
                scan_stmt(s)
            scan_expr(expr.value)

    def scan_stmt(s):
        if isinstance(s, LetStmt):
            declared.add((s.var.name, s.var.version))
            scan_expr(s.value)
        elif isinstance(s, AssignStmt):
            scan_expr(s.value)
        elif isinstance(s, IfStmt):
            scan_expr(s.cond)
            for sub in s.then_body + s.else_body:
                scan_stmt(sub)
        elif isinstance(s, ReturnStmt) and s.value:
            scan_expr(s.value)
        elif isinstance(s, ExprStmt):
            scan_expr(s.expr)
        elif isinstance(s, Require):
            scan_expr(s.cond)
        elif isinstance(s, UnknownStmt):
            for c in s.children:
                scan_stmt(c)

    for s in stmts:
        scan_stmt(s)
    return referenced - declared


def transform_range_for(stmt: UnknownStmt, env: VarEnv, ctx: FuncCtx = None) -> list[StmtIR]:
    """RangeForLoop → 提取为独立 partial def + 调用。

    for (auto& elem : arr) { body }
    →
    partial def range_N (idx : UInt64) (coll : Array T) (vars...) (closure...) : RetType :=
      if idx >= coll.size then <return modified vars>
      else let elem := coll[idx]!; <body>; range_N (idx+1) coll vars' closure

    -- 在主函数中：
    let coll := arr
    let (var1, var2, ...) := range_N 0 coll var1 var2 ... closure1 closure2
    """
    if ctx is None:
        ctx = FuncCtx()

    children = stmt.children
    if len(children) < 3:
        return [stmt]

    collection = children[0].expr if isinstance(children[0], ExprStmt) else UnknownExpr("?")
    if isinstance(children[1], LetStmt):
        loop_var_name = children[1].var.name
        loop_var_type = children[1].typ
    else:
        loop_var_name = "elem"
        loop_var_type = "auto"
    body_stmt = children[2]

    # 结构化绑定（children[3] 存绑定名列表）
    decomp_bindings = []
    if len(children) > 3 and isinstance(children[3], ExprStmt):
        lit = children[3].expr
        if isinstance(lit, Lit) and isinstance(lit.value, list):
            decomp_bindings = lit.value
            loop_var_name = "elem"  # 统一用 elem，然后解构

    collection_renamed = rename_expr(collection, env)

    # auto 类型推导（从集合类型）
    if loop_var_type == "auto" or (isinstance(loop_var_type, str) and "auto" in str(loop_var_type)):
        coll_type_raw = env.get_type(collection.name) if isinstance(collection, Var) else None
        if isinstance(coll_type_raw, ArrayType):
            loop_var_type = coll_type_raw.elem
        elif isinstance(coll_type_raw, StructType):
            from ir_types import PairType
            elem_map = {
                "SparsePolyZp": PairType(StructType("UMonomial", []), StructType("Zp", [])),
                "SparsePolyZZ": PairType(StructType("UMonomial", []), StructType("ZZ", [])),
                "MvPolyZZ": PairType(StructType("MvMonomial", []), StructType("ZZ", [])),
                "MvPolyZp": PairType(StructType("MvMonomial", []), StructType("Zp", [])),
            }
            loop_var_type = elem_map.get(coll_type_raw.name, loop_var_type)

    # 识别被修改的变量
    body_vars = identify_loop_vars(body_stmt) - _collect_local_decls(body_stmt)

    loop_id = abs(hash(str(stmt))) % 100000
    func_name = f"range_{loop_id}"

    # idx 和 coll 变量
    idx_name = f"__idx"
    coll_name = f"__coll"

    # 在循环函数的 env 中初始化变量
    loop_env = env.fork()
    loop_env.versions[idx_name] = 0
    loop_env.set_type(idx_name, BaseType.UINT64)
    idx_var = loop_env.bump(idx_name)  # __idx_1

    loop_env.versions[coll_name] = 0
    coll_type = env.get_type(collection.name) if isinstance(collection, Var) else "auto"
    loop_env.set_type(coll_name, coll_type if coll_type != "auto" else BaseType.UINT64)
    coll_var = loop_env.bump(coll_name)  # __coll_1

    # 循环体前：let elem := coll[idx]
    elem_let = LetStmt(Var(loop_var_name), loop_var_type,
                        ArrayAccess(coll_var, idx_var))
    loop_env.set_type(loop_var_name, loop_var_type)
    loop_env.versions[loop_var_name] = 0

    # 结构化绑定解构：在 elem_let 之后插入 let var := elem.fst, let deg := elem.snd
    decomp_stmts = []
    if decomp_bindings:
        elem_var = Var(loop_var_name)
        if len(decomp_bindings) == 2:
            decomp_stmts.append(LetStmt(Var(decomp_bindings[0]), "auto",
                                         FieldAccess(elem_var, "fst")))
            decomp_stmts.append(LetStmt(Var(decomp_bindings[1]), "auto",
                                         FieldAccess(elem_var, "snd")))
        else:
            # 3+ 绑定：嵌套 Prod 投影
            for i, bname in enumerate(decomp_bindings):
                if i < len(decomp_bindings) - 1:
                    accessor = FieldAccess(elem_var, "fst") if i == 0 else \
                               FieldAccess(Var(f"_decomp_{i}"), "fst")
                    decomp_stmts.append(LetStmt(Var(bname), "auto", accessor))
                    if i + 1 < len(decomp_bindings) - 1:
                        rest = FieldAccess(elem_var, "snd") if i == 0 else \
                               FieldAccess(Var(f"_decomp_{i}"), "snd")
                        decomp_stmts.append(LetStmt(Var(f"_decomp_{i+1}"), "auto", rest))
                else:
                    accessor = FieldAccess(elem_var, "snd") if len(decomp_bindings) == 2 else \
                               FieldAccess(Var(f"_decomp_{i-1}"), "snd")
                    decomp_stmts.append(LetStmt(Var(bname), "auto", accessor))
        # 注册解构变量到环境
        for bname in decomp_bindings:
            loop_env.versions[bname] = 0
            loop_env.set_type(bname, "auto")

    # 修改变量列表（循环携带的状态）
    modified_names = sorted(body_vars)

    # SSA 变换循环体（传入 loop_ctx 支持 break/continue）
    # all_param_names 包含全部循环参数（idx、coll、modified、闭包——闭包后面补）
    all_loop_params = [idx_name, coll_name] + modified_names
    lctx = LoopCtx(func_name, modified_names, step_stmts=[], env=loop_env,
                    all_param_names=all_loop_params)
    all_body_raw = [elem_let] + decomp_stmts + flatten_stmts(body_stmt)
    body_stmts = transform_body(all_body_raw, loop_env, ctx, loop_ctx=lctx)

    # 引用写回：如果循环变量（elem/term）被修改了，追加 Array.set! 写回集合
    # C++ range-for 给的是引用：修改 term = 修改 f[i]
    # Lean 给的是拷贝：需要显式写回到 __coll（遍历集合），
    # 使下次迭代从修改后的数组取元素
    loop_var_ver = loop_env.versions.get(loop_var_name, 0)
    coll_needs_writeback = loop_var_ver > 0
    if coll_needs_writeback:
        # 追加：let __coll_2 := Array.set! __coll_1 __idx_1 term_latest
        # 注意：直接 bump coll_var（遍历集合参数），不创建额外参数
        new_coll = loop_env.bump(coll_name)
        body_stmts.append(LetStmt(new_coll, coll_type if coll_type != "auto" else BaseType.UINT64,
                                   Call("Array.set!",
                                        [coll_var,  # __coll_1（当前遍历集合）
                                         idx_var,
                                         loop_env.current(loop_var_name)])))

    # M1 修复（遗留）：检测 body 中是否有 Array.set! 修改了 __coll → 加入 modified_names
    def _body_modifies_coll(stmts, coll_n):
        for s in stmts:
            if isinstance(s, LetStmt) and isinstance(s.value, Call):
                if "Array.set!" in s.value.func and s.var.name == coll_n:
                    return True
            if isinstance(s, IfStmt):
                if _body_modifies_coll(s.then_body, coll_n) or _body_modifies_coll(s.else_body, coll_n):
                    return True
        return False

    if _body_modifies_coll(body_stmts, coll_name):
        if coll_name not in modified_names:
            modified_names = sorted(set(modified_names) | {coll_name})
            # 重新 SSA 变换（modified_names 变了影响 loop_ctx）
            loop_env2 = env.fork()
            loop_env2.versions[idx_name] = 0
            loop_env2.set_type(idx_name, BaseType.UINT64)
            loop_env2.bump(idx_name)
            loop_env2.versions[coll_name] = 0
            loop_env2.set_type(coll_name, coll_type if coll_type != "auto" else BaseType.UINT64)
            loop_env2.bump(coll_name)
            loop_env2.set_type(loop_var_name, loop_var_type)
            loop_env2.versions[loop_var_name] = 0
            lctx2 = LoopCtx(func_name, modified_names, step_stmts=[], env=loop_env2)
            body_stmts = transform_body(all_body_raw, loop_env2, ctx, loop_ctx=lctx2)
            loop_env = loop_env2

    # 闭包变量：body 中引用但不是 idx/coll/modified/loop_var 的变量
    free = collect_free_vars(body_stmts)
    loop_declared = {idx_name, coll_name, loop_var_name} | set(modified_names)
    # 过滤：只保留外部 env 中存在的变量
    closure_vars = []
    seen_closure = set()
    for vname, vver in sorted(free):
        if vname in loop_declared:
            continue
        if vname in seen_closure:
            continue
        if vname in env.versions:
            seen_closure.add(vname)
            closure_vars.append(env.current(vname))

    # --- 构造循环函数 ---

    # 循环函数参数：idx, coll, modified_vars, closure_vars
    loop_params = [ParamIR(idx_var.lean_name(), BaseType.UINT64)]
    loop_params.append(ParamIR(coll_var.lean_name(), coll_type if coll_type != "auto" else BaseType.UINT64))
    for name in modified_names:
        # 参数名 = 进入循环前的版本（和退出表达式、调用者一致）
        loop_params.append(ParamIR(env.current(name).lean_name(), env.get_type(name)))
    for cv in closure_vars:
        loop_params.append(ParamIR(cv.lean_name(), env.get_type(cv.name)))

    # 循环函数体：if exit then <return> else <body; step; recurse>
    # step: idx + 1
    idx_next = Var(idx_name, idx_var.version + 1)
    step_stmt = LetStmt(idx_next, BaseType.UINT64, BinOp("+", idx_var, Lit(1)))

    # 递归调用参数（如果有写回，用修改后的 __coll）
    coll_for_recurse = loop_env.current(coll_name) if coll_needs_writeback else coll_var
    recurse_args = [idx_next, coll_for_recurse]
    for name in modified_names:
        recurse_args.append(loop_env.current(name))
    for cv in closure_vars:
        recurse_args.append(cv)
    recurse_call = Call(func_name + "_ir", recurse_args)

    # 退出返回值：循环函数参数中修改变量的版本
    # 这些是循环开始时的版本——如果循环条件一开始就不满足，直接返回初始值
    param_versions = {}  # name → 参数中的 Var
    for name in modified_names:
        # 循环函数参数中的版本 = 外部 env 的当前版本（进入循环前）
        param_versions[name] = env.current(name)

    if len(modified_names) == 0:
        ret_type = BaseType.VOID
        exit_expr = Lit(0)
    elif len(modified_names) == 1:
        name = modified_names[0]
        ret_type = env.get_type(name)
        exit_expr = param_versions[name]
    else:
        types = [env.get_type(n) for n in modified_names]
        ret_type = types[0]  # 简化：codegen 用 auto
        exit_expr = Call("Prod.mk", [param_versions[n] for n in modified_names])

    # 组装循环函数体
    exit_cond = BinOp(">=", idx_var, Call("Array.size_u64", [coll_var]))
    loop_body_full = [
        IfStmt(exit_cond,
               [ReturnStmt(exit_expr)],
               body_stmts + [step_stmt, ReturnStmt(recurse_call)])
    ]

    # 注册循环函数
    loop_func = SSAFunc(
        name=func_name,
        params=loop_params,
        ret_type=ret_type,
        requires=[],
        body=loop_body_full,
    )
    ctx.aux_defs.append(loop_func)

    # --- 主函数中的调用 ---
    result_stmts = []

    # let coll := collection
    coll_init_var = env.bump(f"__coll_{loop_id}")
    env.set_type(f"__coll_{loop_id}", coll_type if coll_type != "auto" else BaseType.UINT64)
    result_stmts.append(LetStmt(coll_init_var, coll_type, collection_renamed))

    # 调用参数
    call_args = [Lit(0), coll_init_var]
    for name in modified_names:
        call_args.append(env.current(name))
    for cv in closure_vars:
        call_args.append(cv)
    call_expr = Call(func_name + "_ir", call_args)

    # 绑定返回值
    if len(modified_names) == 0:
        result_stmts.append(ExprStmt(call_expr))
    elif len(modified_names) == 1:
        name = modified_names[0]
        new_var = env.bump(name)
        result_stmts.append(LetStmt(new_var, env.get_type(name), call_expr))
    else:
        # 多变量：解构 Prod
        pair_var = env.bump(f"_loop_{loop_id}")
        env.set_type(f"_loop_{loop_id}", "auto")
        result_stmts.append(LetStmt(pair_var, "auto", call_expr))
        for i, name in enumerate(modified_names):
            new_var = env.bump(name)
            if i == 0:
                accessor = FieldAccess(pair_var, "1")
            elif i == len(modified_names) - 1:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
            else:
                accessor = pair_var
                for _ in range(i):
                    accessor = FieldAccess(accessor, "2")
                accessor = FieldAccess(accessor, "1")
            result_stmts.append(LetStmt(new_var, env.get_type(name), accessor))

    # 写回映射：__coll → 原集合变量（如 f）
    if coll_needs_writeback and isinstance(collection, Var):
        orig_name = collection.name
        if orig_name in env.versions and orig_name != coll_name:
            coll_latest = env.current(coll_name) if coll_name in env.versions else None
            if coll_latest:
                new_orig = env.bump(orig_name)
                result_stmts.append(LetStmt(new_orig, env.get_type(orig_name), coll_latest))

    # G3 后处理：确保循环函数无自由变量
    finalize_loop_func(loop_func, env, result_stmts)

    return result_stmts


def transform_member_assign(lhs: FieldAccess, rhs: ExprIR, env: VarEnv) -> StmtIR | list[StmtIR]:
    """成员赋值 elem.field1.field2 = expr → functional update。

    检测模式：
    - term.second.val = expr → { term with snd := { term.snd with val := expr } }
    - term.field = expr → { term with field := expr }

    如果 term 来自 __coll[__idx]（range-for），追加 Array.set。
    """
    renamed_rhs = rename_expr(rhs, env)

    # 收集字段路径：term.second.val → [("second","snd"), ("val","val")]
    from class_map import FIELD_MAP
    path = []
    node = lhs
    while isinstance(node, FieldAccess):
        lean_field = FIELD_MAP.get(node.field_name, node.field_name)
        path.append((node.field_name, lean_field))
        node = node.obj
    path.reverse()

    # node 现在是根变量（如 term）
    root_var = rename_expr(node, env) if isinstance(node, Var) else node
    root_name = node.name if isinstance(node, Var) else None

    # 根是 ArrayAccess（nodes[idx].field = val）→ Array.set! 模式
    if isinstance(node, ArrayAccess) and path:
        arr = rename_expr(node.arr, env)
        idx = rename_expr(node.idx, env)
        elem_var = Var("_elem")
        update_expr = _build_functional_update(elem_var, path, renamed_rhs)
        # let _elem := arr[idx]!; let _elem' := { _elem with field := val }; let arr' := arr.set! idx _elem'
        arr_name = node.arr.name if isinstance(node.arr, Var) else None
        if arr_name:
            stmts = []
            stmts.append(LetStmt(elem_var, "auto", ArrayAccess(arr, idx)))
            updated_elem = Var("_elem_updated")
            stmts.append(LetStmt(updated_elem, "auto", update_expr))
            new_arr = env.bump(arr_name)
            stmts.append(LetStmt(new_arr, env.get_type(arr_name),
                                 Call("Array.set!", [arr, idx, updated_elem])))
            return stmts
        return ExprStmt(BinOp(":=", rename_expr(lhs, env), renamed_rhs))

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


def _extract_lean_names(expr: ExprIR) -> set[str]:
    """提取表达式中的变量 lean_name（含版本号）。"""
    result = set()
    if isinstance(expr, Var):
        result.add(expr.lean_name())
    elif isinstance(expr, BinOp):
        result |= _extract_lean_names(expr.lhs) | _extract_lean_names(expr.rhs)
    elif isinstance(expr, UnaryOp):
        result |= _extract_lean_names(expr.operand)
    elif isinstance(expr, CondExpr):
        result |= _extract_lean_names(expr.cond) | _extract_lean_names(expr.then_e) | _extract_lean_names(expr.else_e)
    elif isinstance(expr, Call):
        for a in expr.args:
            result |= _extract_lean_names(a)
    elif isinstance(expr, ArrayAccess):
        result |= _extract_lean_names(expr.arr) | _extract_lean_names(expr.idx)
    elif isinstance(expr, FieldAccess):
        result |= _extract_lean_names(expr.obj)
    elif isinstance(expr, Cast):
        result |= _extract_lean_names(expr.expr)
    elif isinstance(expr, ArrayPush):
        result |= _extract_lean_names(expr.arr) | _extract_lean_names(expr.elem)
    return result


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


def _collect_local_decls(stmt: StmtIR) -> set[str]:
    """收集语句中所有 LetStmt 声明的变量名（循环体内的局部变量）。"""
    result = set()
    if isinstance(stmt, LetStmt):
        result.add(stmt.var.name)
    elif isinstance(stmt, IfStmt):
        for s in stmt.then_body + stmt.else_body:
            result |= _collect_local_decls(s)
    elif isinstance(stmt, UnknownStmt):
        for c in stmt.children:
            result |= _collect_local_decls(c)
    return result


def _extract_root_var_name(expr) -> str | None:
    """从赋值 LHS 提取根变量名。统一处理所有 LHS 模式。

    Var("x") → "x"
    ArrayAccess(Var("arr"), idx) → "arr"
    ArrayAccess(ArrayAccess(Var("M"), i), j) → "M"
    FieldAccess(Var("term"), "snd") → "term"
    FieldAccess(ArrayAccess(Var("arr"), i), "snd") → "arr"
    Call("Prod.mk", [Var("q"), Var("r")]) → "q"（+ "r" 由调用者处理）
    """
    if isinstance(expr, Var):
        return expr.name
    if isinstance(expr, ArrayAccess):
        return _extract_root_var_name(expr.arr)
    if isinstance(expr, FieldAccess):
        return _extract_root_var_name(expr.obj)
    if isinstance(expr, Call) and expr.func == "Prod.mk":
        for a in expr.args:
            name = _extract_root_var_name(a)
            if name:
                return name
    return None


def _extract_all_root_var_names(expr) -> set[str]:
    """从赋值 LHS 提取所有根变量名。处理 Prod.mk 多输出。"""
    if isinstance(expr, Call) and expr.func == "Prod.mk":
        names = set()
        for a in expr.args:
            name = _extract_root_var_name(a)
            if name:
                names.add(name)
        return names
    name = _extract_root_var_name(expr)
    return {name} if name else set()


def identify_loop_vars(stmt: StmtIR) -> set[str]:
    """识别语句中被**修改**（而非新声明）的变量名。

    统一检测所有修改模式：
    - x = e（简单赋值）
    - arr[i] = e, arr[i][j] = e（数组赋值）
    - term.field = e（成员赋值）
    - (q, r) = f(...)（多输出赋值）
    - arr.push(x)（追加）
    - f.normalization()（原地修改方法）
    - i++, i--（自增自减）
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
        # 统一赋值检测：从 LHS 提取所有根变量
        if isinstance(expr, BinOp) and expr.op == "=":
            result |= _extract_all_root_var_names(expr.lhs)
        # i++ / i--
        if isinstance(expr, UnaryOp) and expr.op in ("++", "--"):
            if isinstance(expr.operand, Var):
                result.add(expr.operand.name)
        # arr.push(x) → arr 被修改（含嵌套 arr[i].push / arr.field.push）
        if isinstance(expr, ArrayPush):
            result |= _extract_all_root_var_names(expr.arr)
        # N3: compound assignment operator (operator-=, +=, *=)
        if isinstance(expr, Call) and isinstance(expr.func, str):
            if any(op in expr.func for op in ["-=", "+=", "*=", "/=", "%="]):
                if expr.args:
                    result |= _extract_all_root_var_names(expr.args[0])
        # method call 可能修改 this 对象
        if isinstance(expr, Call) and expr.args and isinstance(expr.args[0], Var):
            from class_map import MUTATING_METHODS, NOOP_METHODS
            func = expr.func if isinstance(expr.func, str) else ""
            if func.startswith("_mutate_"):
                result.add(expr.args[0].name)
            elif func in MUTATING_METHODS or func in NOOP_METHODS:
                result.add(expr.args[0].name)
    return result


def _body_has_return(stmt: StmtIR) -> bool:
    """检测循环体中是否含 ReturnStmt（含嵌套循环内的 return）。

    嵌套循环内的 return 也需要外层循环传播 _ret_flag，
    因为内层循环设置 flag 后外层需要检查并 break。
    """
    if isinstance(stmt, ReturnStmt):
        return True
    if isinstance(stmt, IfStmt):
        return any(_body_has_return(s) for s in stmt.then_body + stmt.else_body)
    if isinstance(stmt, UnknownStmt):
        return any(_body_has_return(c) for c in stmt.children)
    return False


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
    """将表达式中的变量引用替换为当前版本。传播 _ast_type。"""
    t = getattr(expr, '_ast_type', None)

    if isinstance(expr, Var):
        r = env.current(expr.name)
        r._ast_type = t or expr._ast_type
        return r
    if isinstance(expr, BinOp):
        r = BinOp(expr.op, rename_expr(expr.lhs, env), rename_expr(expr.rhs, env))
        r._ast_type = t
        return r
    if isinstance(expr, UnaryOp):
        r = UnaryOp(expr.op, rename_expr(expr.operand, env))
        r._ast_type = t
        return r
    if isinstance(expr, CondExpr):
        r = CondExpr(rename_expr(expr.cond, env),
                     rename_expr(expr.then_e, env),
                     rename_expr(expr.else_e, env))
        r._ast_type = t
        return r
    if isinstance(expr, Call):
        r = Call(expr.func, [rename_expr(a, env) for a in expr.args])
        r._ast_type = t
        return r
    if isinstance(expr, ArrayAccess):
        r = ArrayAccess(rename_expr(expr.arr, env), rename_expr(expr.idx, env))
        r._ast_type = t
        return r
    if isinstance(expr, FieldAccess):
        r = FieldAccess(rename_expr(expr.obj, env), expr.field_name)
        r._ast_type = t
        return r
    if isinstance(expr, Cast):
        r = Cast(rename_expr(expr.expr, env), expr.target_type, expr.source_type)
        r._ast_type = t
        return r
    if isinstance(expr, ArrayPush):
        r = ArrayPush(rename_expr(expr.arr, env), rename_expr(expr.elem, env))
        r._ast_type = t
        return r
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
