"""
CLPoly C++ → Lean IR 翻译器：UB 证明目标收集器

遍历 SSAFunc，识别所有可能的未定义行为（UB）点，
生成对应的 Lean 证明目标。

原则：不做推断。只标注 C++ 标准定义的 UB。
- unsigned +/-/* → 无 UB（mod 2^64）
- unsigned / 和 % → 除以零 UB
- signed +/-/* → 溢出 UB
- shift → 移位量 ≥ 位宽 UB
- 数组访问 → 越界 UB
- assert → 断言失败（已在 SSA 阶段提升为 require）
"""

from __future__ import annotations
from ir_types import *
from lean_codegen import gen_expr, gen_type


def collect_all(func: SSAFunc) -> list[UBObligation]:
    """收集函数中所有 UB 证明目标。"""
    obligations = []
    _scan_stmts(func.body, func.name, obligations, [])
    return obligations


def _scan_stmts(stmts: list[StmtIR], func_name: str,
                obligations: list[UBObligation], context: list[str]):
    for stmt in stmts:
        _scan_stmt(stmt, func_name, obligations, context)


def _scan_stmt(stmt: StmtIR, func_name: str,
               obligations: list[UBObligation], context: list[str]):
    if isinstance(stmt, LetStmt):
        ctx = context + [f"let {stmt.var.lean_name()} := {gen_expr(stmt.value)}"]
        _scan_expr(stmt.value, func_name, obligations, ctx)

    elif isinstance(stmt, IfStmt):
        _scan_expr(stmt.cond, func_name, obligations, context)
        _scan_stmts(stmt.then_body, func_name, obligations, context)
        _scan_stmts(stmt.else_body, func_name, obligations, context)

    elif isinstance(stmt, ReturnStmt):
        if stmt.value:
            _scan_expr(stmt.value, func_name, obligations, context)

    elif isinstance(stmt, ExprStmt):
        _scan_expr(stmt.expr, func_name, obligations, context)

    elif isinstance(stmt, TailRec):
        _scan_expr(stmt.exit_cond, func_name, obligations, context)
        if stmt.break_cond:
            _scan_expr(stmt.break_cond, func_name, obligations, context)
        _scan_stmts(stmt.body, func_name, obligations, context)
        _scan_stmts(stmt.step, func_name, obligations, context)

    elif isinstance(stmt, UnknownStmt):
        _scan_stmts(stmt.children, func_name, obligations, context)


def _scan_expr(expr: ExprIR, func_name: str,
               obligations: list[UBObligation], context: list[str]):
    """扫描表达式中的 UB 点。"""

    if isinstance(expr, BinOp):
        _scan_expr(expr.lhs, func_name, obligations, context)
        _scan_expr(expr.rhs, func_name, obligations, context)

        # 除以零：a / b, a % b
        if expr.op in ("/", "%"):
            rhs_str = gen_expr(expr.rhs)
            obligations.append(UBObligation(
                func_name=func_name,
                source_line=0,
                ub_type=UBType.DIV_BY_ZERO,
                lean_prop=f"{rhs_str} ≠ 0",
                context=list(context),
            ))

        # 移位越界：a << n, a >> n
        if expr.op in ("<<", ">>"):
            rhs_str = gen_expr(expr.rhs)
            obligations.append(UBObligation(
                func_name=func_name,
                source_line=0,
                ub_type=UBType.SHIFT_OOB,
                lean_prop=f"{rhs_str} < 64",
                context=list(context),
            ))

        # 有符号溢出：signed int 的 +, -, *
        # unsigned 无 UB（mod 2^64），不检查。
        # signed 溢出是 UB（C++ [expr]/4），必须检查。
        # 检测方法：从 Clang AST 的类型信息判断操作数是否 signed。
        # 当前简化：检查表达式中是否有 Cast 到 INT64 的操作数。
        if expr.op in ("+", "-", "*"):
            if _is_signed_context(expr):
                lhs_str = gen_expr(expr.lhs)
                rhs_str = gen_expr(expr.rhs)
                obligations.append(UBObligation(
                    func_name=func_name,
                    source_line=0,
                    ub_type=UBType.SIGNED_OVERFLOW,
                    lean_prop=f"Int.noOverflow ({lhs_str} {expr.op} {rhs_str})",
                    context=list(context),
                ))

    elif isinstance(expr, ArrayAccess):
        _scan_expr(expr.arr, func_name, obligations, context)
        _scan_expr(expr.idx, func_name, obligations, context)
        idx_str = gen_expr(expr.idx)
        arr_str = gen_expr(expr.arr)
        obligations.append(UBObligation(
            func_name=func_name,
            source_line=0,
            ub_type=UBType.ARRAY_OOB,
            lean_prop=f"{idx_str} < {arr_str}.size",
            context=list(context),
        ))

    elif isinstance(expr, FieldAccess):
        _scan_expr(expr.obj, func_name, obligations, context)
        # .front! / .back! 在空容器上是 UB
        if expr.field_name in ("front!", "back!"):
            obj_str = gen_expr(expr.obj)
            obligations.append(UBObligation(
                func_name=func_name,
                source_line=0,
                ub_type=UBType.ARRAY_OOB,
                lean_prop=f"¬ {obj_str}.isEmpty",
                context=list(context),
            ))

    elif isinstance(expr, Call):
        for arg in expr.args:
            _scan_expr(arg, func_name, obligations, context)

    elif isinstance(expr, Cast):
        _scan_expr(expr.expr, func_name, obligations, context)
        # unsigned → signed 转换：值可能超出 signed 范围
        # C++17: 实现定义（非 UB），但 C++20 前可能导致错误行为
        # 保守处理：生成证明目标
        if (isinstance(expr.source_type, BaseType) and
            isinstance(expr.target_type, BaseType)):
            # UINT64 → INT64：值 > INT64_MAX 时行为实现定义
            if (expr.source_type == BaseType.UINT64 and
                expr.target_type == BaseType.INT64):
                inner_str = gen_expr(expr.expr)
                obligations.append(UBObligation(
                    func_name=func_name,
                    source_line=0,
                    ub_type=UBType.SIGNED_OVERFLOW,
                    lean_prop=f"{inner_str}.val ≤ Int64.max",
                    context=list(context),
                ))
            # UINT128 → INT64 同理
            if (expr.source_type == BaseType.UINT128 and
                expr.target_type == BaseType.INT64):
                inner_str = gen_expr(expr.expr)
                obligations.append(UBObligation(
                    func_name=func_name,
                    source_line=0,
                    ub_type=UBType.SIGNED_OVERFLOW,
                    lean_prop=f"{inner_str} ≤ Int64.max",
                    context=list(context),
                ))

    elif isinstance(expr, CondExpr):
        _scan_expr(expr.cond, func_name, obligations, context)
        _scan_expr(expr.then_e, func_name, obligations, context)
        _scan_expr(expr.else_e, func_name, obligations, context)

    elif isinstance(expr, UnaryOp):
        _scan_expr(expr.operand, func_name, obligations, context)

    elif isinstance(expr, FieldAccess):
        _scan_expr(expr.obj, func_name, obligations, context)

    elif isinstance(expr, ArrayPush):
        _scan_expr(expr.arr, func_name, obligations, context)
        _scan_expr(expr.elem, func_name, obligations, context)

    elif isinstance(expr, BinOp) and expr.op == ":=":
        _scan_expr(expr.lhs, func_name, obligations, context)
        _scan_expr(expr.rhs, func_name, obligations, context)


# ============================================================
# 输出格式化
# ============================================================

def _is_trivial(ob: UBObligation) -> bool:
    """过滤明显无意义的 UB 目标。"""
    # 常数移位量：如果 RHS 是字面量且 < 64，不需要证明
    if ob.ub_type == UBType.SHIFT_OOB:
        # "(64 : UInt64) < 64" 这种是 UInt128 >>> 64 产生的，不是真 UB
        if "64) < 64" in ob.lean_prop or "(0 : " in ob.lean_prop:
            return True
    # signed overflow 对 unsigned 运算的误报
    if ob.ub_type == UBType.SIGNED_OVERFLOW:
        # 如果 prop 中只有 UInt64 字面量，说明是 unsigned 运算
        if ".toNat.toUInt64" in ob.lean_prop:
            return True
    return False


def _is_signed_context(expr: BinOp) -> bool:
    """判断二元运算是否在 signed 上下文中。

    只有当**两个操作数的最终类型都是 signed**时才算 signed 运算。
    如果任一操作数是 UInt64（或 Cast 到 UInt64），则是 unsigned 运算。
    """
    def _final_type(e: ExprIR) -> BaseType | None:
        if isinstance(e, Cast):
            return e.target_type if isinstance(e.target_type, BaseType) else None
        if isinstance(e, Lit):
            return e.typ
        return None

    lhs_type = _final_type(expr.lhs)
    rhs_type = _final_type(expr.rhs)

    # 只有两边都明确是 INT64 才算 signed
    if lhs_type == BaseType.INT64 and rhs_type == BaseType.INT64:
        return True
    # 一边是 INT64 且另一边未知（Var 等），保守认为 signed
    if (lhs_type == BaseType.INT64 and rhs_type is None) or \
       (rhs_type == BaseType.INT64 and lhs_type is None):
        return True
    return False


def inject_obligations(func: SSAFunc, obligations: list[UBObligation]):
    """将 UB 证明目标注入到函数的 requires 列表中。"""
    seen = set()
    for ob in obligations:
        # 过滤无意义/不可用的目标
        if _is_trivial(ob):
            continue
        # 引用循环内部变量（__idx, __coll）的目标不能提升到函数签名
        if "__idx" in ob.lean_prop or "__coll" in ob.lean_prop:
            continue
        if ob.lean_prop in seen:
            continue
        seen.add(ob.lean_prop)
        # 生成 require 名称
        tag = ob.ub_type.name.lower()
        idx = len(func.requires) + 1
        name = f"h_{tag}_{idx}"
        # 构造 Require 节点（cond 用 Var 包装 lean_prop 字符串）
        req = Require(
            cond=Var(ob.lean_prop),  # lean_prop 作为 Prop 表达式
            name=name,
            source=tag,
        )
        func.requires.append(req)


def format_obligations(obligations: list[UBObligation]) -> str:
    """格式化为 Lean 注释 + theorem 骨架。"""
    if not obligations:
        return "-- No UB proof obligations.\n"

    lines = [f"-- UB proof obligations: {len(obligations)} total", ""]

    by_func = {}
    for ob in obligations:
        by_func.setdefault(ob.func_name, []).append(ob)

    for func_name, obs in by_func.items():
        lines.append(f"-- Function: {func_name}")
        for i, ob in enumerate(obs):
            tag = ob.ub_type.name.lower()
            lines.append(f"--   [{i+1}] {ob.ub_type.name}: {ob.lean_prop}")
            if ob.context:
                lines.append(f"--       context: {ob.context[-1]}")
        lines.append("")

    # 生成 theorem 骨架
    lines.append("-- Proof obligation theorems (fill sorry to complete):")
    lines.append("")
    for func_name, obs in by_func.items():
        for i, ob in enumerate(obs):
            tag = ob.ub_type.name.lower()
            name = f"ub_{func_name}_{tag}_{i+1}"
            lines.append(f"theorem {name} : {ob.lean_prop} := by")
            lines.append(f"  sorry")
            lines.append("")

    return "\n".join(lines)


def format_summary(obligations: list[UBObligation]) -> str:
    """简要统计。"""
    by_type = {}
    for ob in obligations:
        by_type[ob.ub_type] = by_type.get(ob.ub_type, 0) + 1

    lines = [f"UB obligations: {len(obligations)} total"]
    for typ, count in sorted(by_type.items(), key=lambda x: x[0].name):
        lines.append(f"  {typ.name}: {count}")
    return "\n".join(lines)
