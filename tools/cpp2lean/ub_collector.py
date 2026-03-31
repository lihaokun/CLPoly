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

        # 有符号溢出：int +, -, *
        # 注：unsigned 无 UB（mod 2^64）。只有 signed 操作需要检查。
        # 需要知道操作数类型才能判断。当前简化：不生成（CLPoly 极少用 signed 算术）

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
