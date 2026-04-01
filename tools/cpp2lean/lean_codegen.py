"""
CLPoly C++ → Lean IR 翻译器：Lean 4 代码生成器

输入：list[SSAFunc]（或 list[FuncIR] for Phase 0 简化模式）
输出：Lean 4 源代码字符串
"""

from __future__ import annotations
from ir_types import *


# ============================================================
# 类型生成
# ============================================================

def gen_type(t: TypeIR) -> str:
    if isinstance(t, BaseType):
        return t.value
    if isinstance(t, ArrayType):
        return f"Array ({gen_type(t.elem)})"
    if isinstance(t, PairType):
        return f"({gen_type(t.fst)} × {gen_type(t.snd)})"
    if isinstance(t, StructType):
        return t.name
    if isinstance(t, ExceptType):
        return f"Except Error ({gen_type(t.inner)})"
    if isinstance(t, str):
        return f"/- unknown: {t} -/ sorry"
    return "sorry"


# ============================================================
# 表达式生成
# ============================================================

def gen_coercion(expr_str: str, source, target) -> str:
    """统一类型转换：source 类型表达式 → target 类型。"""
    if source == target:
        return expr_str
    # BaseType → BaseType
    if isinstance(source, BaseType) and isinstance(target, BaseType):
        key = (source, target)
        if key in CAST_TABLE:
            return CAST_TABLE[key].format(e=expr_str)
    # StructType → BaseType
    if isinstance(source, StructType) and isinstance(target, BaseType):
        if source.name == "Zp" and target == BaseType.UINT64:
            return f"{expr_str}.val"
        if source.name == "Zp" and target == BaseType.INT64:
            return f"({expr_str}.val.toNat : Int)"
        if source.name == "UMonomial" and target == BaseType.UINT64:
            return f"{expr_str}.deg"
    # BaseType → StructType
    if isinstance(source, BaseType) and isinstance(target, StructType):
        return f"({expr_str} : {gen_type(target)})"
    # str types (unrecognized)
    if isinstance(source, str) or isinstance(target, str):
        return f"({expr_str} : {gen_type(target)})"
    # fallback
    return f"({expr_str} : {gen_type(target)})"


def gen_require_prop(req) -> str:
    """Require → Lean Prop。"""
    cond = req.cond
    # UB obligations 注入的 Var：lean_prop 字符串直接作为 Prop
    if isinstance(cond, Var) and any(c in cond.name for c in ("≠", "<", "≤", "¬", "noOverflow")):
        return cond.name
    if isinstance(cond, BinOp) and cond.op == "!=":
        return f"{gen_expr(cond.lhs)} ≠ {gen_expr(cond.rhs)}"
    if isinstance(cond, BinOp) and cond.op == "==":
        return f"{gen_expr(cond.lhs)} = {gen_expr(cond.rhs)}"
    if isinstance(cond, BinOp) and cond.op == "<":
        return f"{gen_expr(cond.lhs)} < {gen_expr(cond.rhs)}"
    return gen_expr(cond)


# ============================================================
# 统一类型转换表：(source, target) → Lean 表达式模板
# {e} 是内部表达式。翻译器不推断类型，只查此表。
# ============================================================

CAST_TABLE = {
    # 截断
    (BaseType.UINT128, BaseType.UINT64): "(uint128_lo {e})",
    (BaseType.UINT64, BaseType.UINT32):  "({e}).toUInt32",
    (BaseType.INT64, BaseType.UINT32):   "({e}).toNat.toUInt32",
    # 扩展
    (BaseType.UINT64, BaseType.UINT128): "({e} : UInt128)",
    (BaseType.UINT32, BaseType.UINT64):  "({e}).toUInt64",
    (BaseType.UINT32, BaseType.UINT128): "(({e}).toUInt64 : UInt128)",
    # 有符号 ↔ 无符号
    (BaseType.INT64, BaseType.UINT64):   "({e}).toNat.toUInt64",
    (BaseType.UINT64, BaseType.INT64):   "(({e}).toNat : Int)",
    (BaseType.INT64, BaseType.UINT128):  "(({e}).toNat.toUInt64 : UInt128)",
    (BaseType.UINT32, BaseType.INT64):   "(({e}).toNat : Int)",
    (BaseType.UINT128, BaseType.INT64):  "((uint128_lo {e}).toNat : Int)",
    # Bool
    (BaseType.BOOL, BaseType.UINT64):    "(if {e} then (1 : UInt64) else 0)",
    (BaseType.BOOL, BaseType.INT64):     "(if {e} then (1 : Int) else 0)",
    (BaseType.UINT64, BaseType.BOOL):    "({e} != 0)",
    (BaseType.INT64, BaseType.BOOL):     "({e} != 0)",
}

# Lean 标准库函数（不加 _ir 后缀）
LEAN_STDLIB = {
    "Prod.mk": "Prod.mk",
    "Array.empty": "#[]",
    "Array.mk": "Array.mk",
    "Array.set!": "Array.set!",
    "Array.size_u64": "Array.size_u64",
    "Zp.mk": "Zp.mk",
    "UMonomial.mk": "UMonomial.mk",
}

# no-op 函数（C++ 内存管理等，翻译时丢弃）
NOOP_FUNCS = {"reserve"}

BINOP_MAP = {
    "+": "+", "-": "-", "*": "*", "/": "/", "%": "%",
    "<<": "<<<", ">>": ">>>",
    "<": "<", ">": ">", "<=": "<=", ">=": ">=",
    "==": "==", "!=": "!=",
    "&&": "&&", "||": "||",
}


def gen_expr(expr: ExprIR) -> str:
    if isinstance(expr, Var):
        return expr.lean_name()

    if isinstance(expr, Lit):
        if expr.typ == BaseType.BOOL:
            return "true" if expr.value else "false"
        return f"({expr.value} : {expr.typ.value})"

    if isinstance(expr, BinOp):
        # 成员赋值 := → Lean 注释（纯函数式中需要 Array.set，暂标注）
        if expr.op == ":=":
            return f"/- {gen_expr(expr.lhs)} := {gen_expr(expr.rhs)} -/ sorry"
        lean_op = BINOP_MAP.get(expr.op, f"/- {expr.op} -/ sorry")
        return f"({gen_expr(expr.lhs)} {lean_op} {gen_expr(expr.rhs)})"

    if isinstance(expr, UnaryOp):
        if expr.op == "!":
            return f"(! {gen_expr(expr.operand)})"
        if expr.op == "-":
            return f"(- {gen_expr(expr.operand)})"
        return f"/- unary {expr.op} -/ {gen_expr(expr.operand)}"

    if isinstance(expr, CondExpr):
        return f"(if {gen_expr(expr.cond)} then {gen_expr(expr.then_e)} else {gen_expr(expr.else_e)})"

    if isinstance(expr, Call):
        func_name = expr.func
        args = " ".join(gen_expr(a) for a in expr.args)
        # _mutate_ 前缀 → 去掉前缀（SSA 处理在 ssa_transform）
        if func_name.startswith("_mutate_"):
            func_name = func_name[len("_mutate_"):]
        # functional update: { obj with field := value }
        if func_name == "_with" and len(expr.args) == 3:
            obj = gen_expr(expr.args[0])
            field = expr.args[1].name if isinstance(expr.args[1], Var) else gen_expr(expr.args[1])
            val = gen_expr(expr.args[2])
            return f"{{ {obj} with {field} := {val} }}"
        # 查找顺序：LEAN_BUILTINS → LEAN_STDLIB → NOOP → TRANSLATION_SCOPE → sorry
        from class_map import LEAN_BUILTINS, TRANSLATION_SCOPE
        if func_name in LEAN_BUILTINS:
            return f"({func_name} {args})" if args else func_name
        if func_name in LEAN_STDLIB:
            lean_name = LEAN_STDLIB[func_name]
            return f"({lean_name} {args})" if args else lean_name
        if func_name in NOOP_FUNCS:
            return args.split()[0] if args else "()"
        # 翻译范围内的函数 → 加 _ir
        if func_name in TRANSLATION_SCOPE:
            return f"({func_name}_ir {args})" if args else f"{func_name}_ir"
        # 都不在 → sorry
        return f"/- unknown func: {func_name} -/ sorry"

    if isinstance(expr, ArrayAccess):
        idx = gen_expr(expr.idx)
        return f"({gen_expr(expr.arr)}[{idx}.toNat]!)"

    if isinstance(expr, FieldAccess):
        return f"{gen_expr(expr.obj)}.{expr.field_name}"

    if isinstance(expr, Cast):
        source = expr.source_type
        target = expr.target_type
        # 同类型 → 跳过
        if source == target:
            return gen_expr(expr.expr)
        # 字面量 Cast → 直接生成目标类型的字面量
        if isinstance(expr.expr, Lit) and isinstance(target, BaseType):
            return f"({expr.expr.value} : {target.value})"
        inner_str = gen_expr(expr.expr)
        # BaseType 查表
        if isinstance(source, BaseType) and isinstance(target, BaseType):
            key = (source, target)
            if key in CAST_TABLE:
                return CAST_TABLE[key].format(e=inner_str)
        # StructType → BaseType（如 Zp → UInt64：取 .val）
        if isinstance(source, StructType) and isinstance(target, BaseType):
            if source.name == "Zp" and target == BaseType.UINT64:
                return f"{inner_str}.val"
            return f"({inner_str} : {gen_type(target)})"
        # BaseType → StructType（如 Int → Zp：构造）
        if isinstance(source, BaseType) and isinstance(target, StructType):
            return f"({inner_str} : {gen_type(target)})"
        # Nat → UInt64（Array.size 返回 Nat）
        if source == "unsigned long" or (isinstance(source, str) and "size" in str(source)):
            return f"({inner_str} : {gen_type(target)})"
        # 同类型或兼容 → 类型标注
        if isinstance(target, (BaseType, StructType)):
            return f"({inner_str} : {gen_type(target)})"
        return f"/- CAST: {source} → {target} -/ {inner_str}"

    if isinstance(expr, ArrayPush):
        return f"({gen_expr(expr.arr)}.push {gen_expr(expr.elem)})"

    if isinstance(expr, UnknownExpr):
        return f"/- UNKNOWN: {expr.kind} -/ sorry"

    return "sorry"


# ============================================================
# 语句生成
# ============================================================

def gen_stmt(stmt: StmtIR, indent: int = 1) -> list[str]:
    pad = "  " * indent
    lines = []

    if isinstance(stmt, LetStmt):
        typ = gen_type(stmt.typ)
        val = gen_expr(stmt.value)
        # 自动 coercion：如果表达式类型和声明类型不同，插入转换
        actual = getattr(stmt.value, '_ast_type', None)
        if actual and actual != stmt.typ:
            val = gen_coercion(val, actual, stmt.typ)
        lines.append(f"{pad}let {stmt.var.lean_name()} : {typ} := {val}")

    elif isinstance(stmt, AssignStmt):
        val = gen_expr(stmt.value)
        lines.append(f"{pad}let {stmt.target.lean_name()} := {val}")

    elif isinstance(stmt, IfStmt):
        lines.append(f"{pad}if {gen_expr(stmt.cond)} then")
        for s in stmt.then_body:
            lines.extend(gen_stmt(s, indent + 1))
        if stmt.else_body:
            lines.append(f"{pad}else")
            for s in stmt.else_body:
                lines.extend(gen_stmt(s, indent + 1))

    elif isinstance(stmt, ReturnStmt):
        if stmt.value is not None:
            val = gen_expr(stmt.value)
            # 返回值如有 _ast_type 则标注类型
            actual = getattr(stmt.value, '_ast_type', None)
            if actual:
                val = f"({val} : {gen_type(actual)})"
            lines.append(f"{pad}{val}")
        else:
            lines.append(f"{pad}()")

    elif isinstance(stmt, Require):
        lines.append(f"{pad}-- require {stmt.name} : {gen_expr(stmt.cond)} [{stmt.source}]")

    elif isinstance(stmt, Throw):
        lines.append(f"{pad}.error .{stmt.error_tag}")

    elif isinstance(stmt, ExprStmt):
        lines.append(f"{pad}let _ := {gen_expr(stmt.expr)}")

    elif isinstance(stmt, TailRec):
        lines.extend(gen_tailrec(stmt, indent))

    elif isinstance(stmt, UnknownStmt):
        if stmt.kind == "CompoundStmt":
            for child in stmt.children:
                lines.extend(gen_stmt(child, indent))
        elif stmt.kind in ("BreakStmt", "ContinueStmt"):
            lines.append(f"{pad}/- {stmt.kind}: handled by TailRec -/")
        else:
            lines.append(f"{pad}/- UNKNOWN_STMT: {stmt.kind} -/")
            lines.append(f"{pad}sorry")

    return lines


def gen_tailrec(tr: TailRec, indent: int) -> list[str]:
    """TailRec → partial def 尾递归。"""
    pad = "  " * indent
    lines = []

    # 参数列表
    params = " ".join(f"({v.lean_name()} : {gen_type(t)})"
                      for v, t in tr.params)

    lines.append(f"{pad}-- loop: {tr.func_name}")
    lines.append(f"{pad}let rec {tr.func_name} {params} :=")

    inner_pad = "  " * (indent + 1)

    # 退出条件
    lines.append(f"{inner_pad}if {gen_expr(tr.exit_cond)} then")
    # 返回累积状态：优先选名为 result/acc/sum 的参数，否则最后一个
    if tr.params:
        acc_name = _pick_accumulator(tr.params)
        lines.append(f"{inner_pad}  {acc_name}")
    else:
        lines.append(f"{inner_pad}  ()")

    # break 条件
    if tr.break_cond is not None:
        lines.append(f"{inner_pad}else if {gen_expr(tr.break_cond)} then")
        if tr.params:
            lines.append(f"{inner_pad}  {tr.params[-1][0].lean_name()}")
        else:
            lines.append(f"{inner_pad}  ()")

    lines.append(f"{inner_pad}else")

    # 预计算 latest（用于 continue 递归调用）
    latest: dict[str, str] = {}
    for v, _ in tr.params:
        latest[v.name] = v.lean_name()
    for s in tr.body + tr.step:
        if isinstance(s, LetStmt) and s.var.name in latest:
            latest[s.var.name] = s.var.lean_name()
        if isinstance(s, IfStmt):
            for sub in s.then_body + s.else_body:
                if isinstance(sub, LetStmt) and sub.var.name in latest:
                    latest[sub.var.name] = sub.var.lean_name()

    # 循环体
    body_pad = "  " * (indent + 2)
    for s in tr.body:
        if isinstance(s, UnknownStmt) and s.kind in ("BreakStmt", "ContinueStmt"):
            continue
        if isinstance(s, IfStmt) and _is_break_if(s):
            continue
        if isinstance(s, IfStmt) and _is_continue_if(s):
            # continue = 跳过剩余 body，执行 step，用当前（非更新）的累积变量递归
            lines.append(f"{body_pad}if {gen_expr(s.cond)} then")
            cont_pad = "  " * (indent + 3)
            # step（如 i++）
            for step_s in tr.step:
                lines.extend(gen_stmt(step_s, indent + 3))
            # 递归：step 变量更新，其他变量保持当前版本（不是 latest）
            cont_args = {}
            for v, _ in tr.params:
                cont_args[v.name] = v.lean_name()  # 保持进入循环体时的版本
            for step_s in tr.step:
                if isinstance(step_s, LetStmt) and step_s.var.name in cont_args:
                    cont_args[step_s.var.name] = step_s.var.lean_name()
            if tr.params:
                c_args = " ".join(cont_args[v.name] for v, _ in tr.params)
                lines.append(f"{cont_pad}{tr.func_name} {c_args}")
            lines.append(f"{body_pad}else")
            continue
        body_lines = gen_stmt(s, indent + 2)
        lines.extend(body_lines)

    # 步进 + 递归调用
    for s in tr.step:
        lines.extend(gen_stmt(s, indent + 2))

    # 递归调用使用 latest 中的最新版本
    if tr.params:
        args = " ".join(latest[v.name] for v, _ in tr.params)
        lines.append(f"{body_pad}{tr.func_name} {args}")
    else:
        lines.append(f"{body_pad}{tr.func_name}")

    # 初始调用
    if tr.params:
        init_args = " ".join(v.lean_name() for v, _ in tr.params)
        lines.append(f"{pad}{tr.func_name} {init_args}")

    return lines


def _pick_accumulator(params: list) -> str:
    """选择循环的累积变量作为退出返回值。"""
    # 优先选名字含 result/acc/sum/out 的
    for v, _ in params:
        name = v.lean_name()
        if any(kw in v.name for kw in ("result", "acc", "sum", "out", "coll")):
            return name
    # 否则返回最后一个参数
    return params[-1][0].lean_name()


def _is_break_if(stmt: IfStmt) -> bool:
    return any(isinstance(s, UnknownStmt) and s.kind == "BreakStmt"
               for s in stmt.then_body)

def _is_continue_if(stmt: IfStmt) -> bool:
    return any(isinstance(s, UnknownStmt) and s.kind == "ContinueStmt"
               for s in stmt.then_body)


# ============================================================
# 函数生成
# ============================================================

def gen_func_ir(func: FuncIR) -> str:
    """Phase 0 简化模式：直接从 FuncIR 生成（无 SSA 变换）。"""
    # 参数
    params = []
    for p in func.params:
        typ = gen_type(p.typ)
        params.append(f"({p.name} : {typ})")
    param_str = " ".join(params)

    # 返回类型
    ret = gen_type(func.ret_type)

    # 函数体
    body_lines = []
    for stmt in func.body:
        body_lines.extend(gen_stmt(stmt))

    body_str = "\n".join(body_lines)

    return f"""partial def {func.name}_ir {param_str} : {ret} :=
{body_str}"""


def gen_ssa_func(func: SSAFunc) -> str:
    """完整模式：从 SSAFunc 生成（Phase 1+）。"""
    # 值参数（先于 require，因为 require 可能引用参数名）
    val_params = []
    for p in func.params:
        typ = gen_type(p.typ)
        val_params.append(f"({p.name} : {typ})")

    # require 参数（在值参数之后）
    req_params = []
    for r in func.requires:
        req_params.append(f"({r.name} : {gen_require_prop(r)})")

    all_params = " ".join(val_params + req_params)
    ret = gen_type(func.ret_type)

    body_lines = []
    for stmt in func.body:
        body_lines.extend(gen_stmt(stmt))

    body_str = "\n".join(body_lines)

    return f"""partial def {func.name}_ir {all_params} : {ret} :=
{body_str}"""


# ============================================================
# 文件生成
# ============================================================

PRELUDE = """-- UInt128 模型
structure UInt128 where
  hi : UInt64
  lo : UInt64
deriving Repr

-- 128 位算术（简化实现，精化证明时替换为正确版本）
instance : Mul UInt128 where
  mul a b := { hi := 0, lo := a.lo * b.lo }  -- placeholder

instance : HShiftRight UInt128 UInt64 UInt64 where
  hShiftRight a n := if n == 64 then a.hi else a.lo  -- placeholder

instance : Coe UInt64 UInt128 where
  coe x := { hi := 0, lo := x }

def uint128_lo (x : UInt128) : UInt64 := x.lo

-- C++ 移位允许 RHS 为不同整数类型；Lean 需要显式实例
instance : HShiftLeft UInt64 UInt32 UInt64 where
  hShiftLeft a n := a <<< n.toUInt64

instance : HShiftRight UInt64 UInt32 UInt64 where
  hShiftRight a n := a >>> n.toUInt64
"""


def generate_file(funcs: list[FuncIR], header: str = "") -> str:
    parts = [
        "-- Auto-generated Lean IR from C++ via cpp2lean",
        "-- DO NOT EDIT: regenerate from C++ source",
        "",
    ]
    if header:
        parts.append(header)
        parts.append("")

    # 检测是否需要 UInt128
    needs_uint128 = any(_uses_uint128(f) for f in funcs)
    if needs_uint128:
        parts.append(PRELUDE)
        parts.append("")

    for func in funcs:
        parts.append(gen_func_ir(func))
        parts.append("")

    return "\n".join(parts)


def _uses_uint128(func) -> bool:
    """检查函数是否使用 UInt128。"""
    from ir_types import BaseType
    for p in func.params:
        if p.typ == BaseType.UINT128:
            return True
    if func.ret_type == BaseType.UINT128:
        return True
    # 检查函数体中的类型
    def check_stmts(stmts):
        for s in stmts:
            if hasattr(s, 'typ') and s.typ == BaseType.UINT128:
                return True
            for attr in ('then_body', 'else_body', 'body', 'children'):
                sub = getattr(s, attr, None)
                if sub and isinstance(sub, list) and check_stmts(sub):
                    return True
        return False
    return check_stmts(func.body)
