"""
Pass 8: codegen — MIR₁ → Lean 4 源码字符串

输入：MIRFunc（含 cfg + aux_defs，无 back edge，所有循环已 lower 为 partial def）
输出：完整 .lean 源代码（含 import + namespace）

阶段（按工作流细化方案）：
  S1（本提交）：emit_type / emit_expr / emit_stmt + 单测
  S2: emit_cfg + 终止子 + merge BB lambda
  S3: emit_mirfunc + codegen_pass + 全 67 函数 smoke
  S4: lake build 端到端

设计决策（已确认）：
  (a) 命名      — MIRFunc.lean_name（既存）
  (b) 类型      — emit_type 扩展 v1 gen_type 4 个新分支
  (c) loop      — aux_defs 顶层 partial def，TailCallTerm → 普通调用
  (d-A2) merge  — 局部 `let bb_N := fun phi_targets => ...`
  (e-B2) phi    — merge BB lambda params
  (f) tailcall  — `target_func_ir args...`
  (g) 关键字    — Lean 4 «$name» 包裹保险
  (h) 解构      — `let (a, b, c) := tmp`

Sorry 降级（残留容错，不阻塞 Pass 8）：
  RefType            → sorry /- ref residual: <inner> -/
  UnknownType("")    → sorry /- unknown -/
  UnknownType(raw)   → sorry /- unknown: raw -/
  Lambda/LambdaRef   → sorry /- lambda residual -/
  UnknownExpr        → sorry /- unknown expr: kind -/
"""

from __future__ import annotations
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType,
    OptionType, StdMapType, ArrayType, RefType,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, RequireStmt, PhiStmt, MIRStmt,
    JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG, HIRParam, MIRFunc,
)


# ============================================================
# 常量表（与 v1 lean_codegen 对齐 + v2 新增）
# ============================================================

# Lean 4 关键字（不全，只列 SSA 名常见冲突；如撞需 «»）
LEAN_KEYWORDS = frozenset({
    "let", "fun", "match", "with", "do", "if", "then", "else",
    "where", "by", "and", "or", "not", "in", "of", "from",
    "structure", "inductive", "class", "instance", "def", "theorem",
    "axiom", "variable", "namespace", "section", "end", "open",
    "import", "private", "protected", "partial", "noncomputable",
    "mutual", "deriving", "syntax", "macro", "elab", "set_option",
    "abbrev", "example", "Type", "Prop", "Sort",
})

BINOP_MAP = {
    "+": "+", "-": "-", "*": "*", "/": "/", "%": "%",
    "<<": "<<<", ">>": ">>>",
    "<": "<", ">": ">", "<=": "<=", ">=": ">=",
    "==": "==", "!=": "!=",
    "&&": "&&", "||": "||",
    "&": "&&&", "|": "|||", "^": "^^^",  # bitwise（Lean 4 BitVec/UInt 用三字符）
}

CAST_TABLE = {
    # 截断
    (BaseType.UINT128, BaseType.UINT64): "(uint128_lo {e})",
    (BaseType.UINT64, BaseType.UINT32):  "({e}).toUInt32",
    (BaseType.INT64, BaseType.UINT32):   "({e}).toNat.toUInt32",
    # 扩展
    (BaseType.UINT64, BaseType.UINT128): "({e} : UInt128)",
    (BaseType.UINT32, BaseType.UINT64):  "({e}).toUInt64",
    (BaseType.UINT32, BaseType.UINT128): "(({e}).toUInt64 : UInt128)",
    (BaseType.NAT,    BaseType.UINT64):  "({e}).toUInt64",
    (BaseType.NAT,    BaseType.INT64):   "(({e} : Int))",
    (BaseType.UINT64, BaseType.NAT):     "({e}).toNat",
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


# ============================================================
# EmitCtx — emit 阶段共享状态
# ============================================================

@dataclass
class EmitCtx:
    """共享 emit 状态。

    indent — 当前 emit 缩进层级（每级 2 空格）
    cfg / merge_bbs — S2 emit_cfg 注入；emit_jump_to 在 ctx 内查询
    """
    indent: int = 0
    cfg: Optional[CFG] = None
    merge_bbs: set[int] = field(default_factory=set)


# ============================================================
# 工具：Sorry 降级 + 名字保护
# ============================================================

def _sorry(reason: str) -> str:
    return f"sorry /- {reason} -/"


def _safe_ident(name: str) -> str:
    """Lean 关键字撞名 → 用 «» 包裹。"""
    if name in LEAN_KEYWORDS:
        return f"«{name}»"
    return name


def emit_var_name(v: Var) -> str:
    """SSA Var → Lean 标识符（含版本 + «» 保护）。"""
    return _safe_ident(v.lean_name())


# ============================================================
# 类型 emit
# ============================================================

def emit_type(ty: Optional[TypeIR]) -> str:
    """TypeIR → Lean 类型表达式。"""
    if ty is None:
        return _sorry("type is None")
    if isinstance(ty, BaseType):
        return ty.value
    if isinstance(ty, NamedType):
        # Lambda/LambdaRef 残留：Pass 3 lambda_lift 漏过的 case
        if ty.name in ("Lambda", "LambdaRef"):
            return _sorry(f"lambda residual type: {ty.name}")
        # Option 单名（无 type arg）—— Pass 4 上游漏指定 inner type；
        # 降级为 Option Unit 占位（语义不重要，B2B 测试时再细化）
        if ty.name == "Option":
            return "Option Unit"
        return ty.name
    if isinstance(ty, ArrayType):
        return f"Array {_paren(emit_type(ty.elem))}"
    if isinstance(ty, PairType):
        return f"({emit_type(ty.fst)} × {emit_type(ty.snd)})"
    if isinstance(ty, TupleType):
        if len(ty.elems) == 0:
            return "Unit"
        if len(ty.elems) == 1:
            return emit_type(ty.elems[0])
        # Lean 4 元组 = 嵌套 Prod，糖语法 (T1 × T2 × T3) 右结合
        return "(" + " × ".join(emit_type(e) for e in ty.elems) + ")"
    if isinstance(ty, OptionType):
        return f"Option {_paren(emit_type(ty.inner))}"
    if isinstance(ty, StdMapType):
        return f"StdMap {_paren(emit_type(ty.key))} {_paren(emit_type(ty.value))}"
    if isinstance(ty, RefType):
        # MIR 阶段不应有 RefType 残留——降级为 inner 类型 + 注释
        return f"{emit_type(ty.inner)} /- ref residual -/"
    if isinstance(ty, UnknownType):
        if ty.raw:
            return _sorry(f"unknown type: {ty.raw}")
        return _sorry("unknown")
    return _sorry(f"unhandled type: {type(ty).__name__}")


def _paren(s: str) -> str:
    """若 s 包含空格且未被包裹，加圆括号。用于嵌套 Lean 类型表达式。"""
    s = s.strip()
    if not s:
        return s
    if " " in s and not (s.startswith("(") and s.endswith(")")):
        return f"({s})"
    return s


# ============================================================
# 字面量 emit
# ============================================================

def emit_lit(l: Lit) -> str:
    """Lit → Lean 字面量。

    整数字面量不带类型标注（让 Lean 上下文推断）—— Pass 1 的 Lit.ty 常常过于
    具体（C++ `0` 默认 Int32，但 Lean 上下文可能期望 Int / Nat / Int64）。
    Pass 5 的 Cast 链处理显式转换；本函数仅给值 + 必要时显式标注。

    Bool / Unit / Float / String 保留类型语义（符号不同）。
    """
    if l.ty == BaseType.BOOL:
        return "true" if l.value else "false"
    if l.ty == BaseType.UNIT:
        return "()"
    if l.ty == BaseType.FLOAT:
        return f"({l.value} : Float)"
    if isinstance(l.value, str):
        # 字符串字面量
        return f"\"{l.value}\""
    # 整数字面量：不标注类型，让 Lean 推断
    return f"{l.value}"


# ============================================================
# 表达式 emit
# ============================================================

def emit_expr(e: ExprIR, ctx: EmitCtx) -> str:
    """ExprIR → Lean 表达式字符串。"""
    if isinstance(e, Var):
        return emit_var_name(e)

    if isinstance(e, Lit):
        return emit_lit(e)

    if isinstance(e, BinOp):
        lean_op = BINOP_MAP.get(e.op, None)
        if lean_op is None:
            return _sorry(f"unmapped binop: {e.op}")
        return f"({emit_expr(e.lhs, ctx)} {lean_op} {emit_expr(e.rhs, ctx)})"

    if isinstance(e, UnaryOp):
        if e.op == "!":
            return f"(! {emit_expr(e.operand, ctx)})"
        if e.op == "-":
            return f"(- {emit_expr(e.operand, ctx)})"
        if e.op == "~":
            return f"(~~~ {emit_expr(e.operand, ctx)})"
        # ++ / -- / * / -> 在 MIR 中不应残留——Pass 5 已展开
        return _sorry(f"unmapped unary: {e.op}")

    if isinstance(e, CondExpr):
        return (f"(if {emit_expr(e.cond, ctx)} "
                f"then {emit_expr(e.then_e, ctx)} "
                f"else {emit_expr(e.else_e, ctx)})")

    if isinstance(e, Call):
        return emit_call(e, ctx)

    if isinstance(e, ArrayAccess):
        idx = emit_expr(e.idx, ctx)
        arr = emit_expr(e.arr, ctx)
        # Lean Array.get! 需 Nat。按 idx 类型分派：
        #   Nat                → 直接用
        #   UInt32/64/128 / Int → .toNat
        #   Int32/Int64        → .toNatClampNeg（无 .toNat 方法；负数夹零）
        # idx 强制 (...)  包裹（即使是裸字面量），避免 `0.toNat` 被解析 decimal。
        idx_p = idx if (idx.startswith("(") and idx.endswith(")")) else f"({idx})"
        idx_ty = getattr(e.idx, 'ty', None)
        if idx_ty == BaseType.NAT:
            return f"({arr}[{idx_p}]!)"
        if idx_ty in (BaseType.INT32, BaseType.INT64):
            return f"({arr}[{idx_p}.toNatClampNeg]!)"
        return f"({arr}[{idx_p}.toNat]!)"

    if isinstance(e, FieldAccess):
        obj = emit_expr(e.obj, ctx)
        # C++ std::pair 字段名 first/second → Lean Prod fst/snd
        field = e.field_name
        if field == "first": field = "fst"
        elif field == "second": field = "snd"
        return f"{obj}.{field}"

    if isinstance(e, Cast):
        return emit_cast(e, ctx)

    if isinstance(e, TupleExpr):
        if len(e.elems) == 0:
            return "()"
        if len(e.elems) == 1:
            return emit_expr(e.elems[0], ctx)
        return "(" + ", ".join(emit_expr(x, ctx) for x in e.elems) + ")"

    if isinstance(e, ArrayLit):
        elems = ", ".join(emit_expr(x, ctx) for x in e.elems)
        return f"#[{elems}]"

    if isinstance(e, BlockExpr):
        # Pass 6+ 不应有 BlockExpr 残留，但容错
        lines = []
        for s in e.stmts:
            lines.append(emit_stmt(s, ctx))
        lines.append(emit_expr(e.value, ctx))
        return "\n".join(lines)

    if isinstance(e, LambdaExpr):
        # Pass 3 lambda_lift 应已消除——降级
        return _sorry("lambda residual in expr")

    if isinstance(e, UnknownExpr):
        return _sorry(f"unknown expr: {e.kind}")

    if isinstance(e, UnresolvedOp):
        # Pass 5 operator_resolve 应已替换——降级
        return _sorry(f"unresolved op: {e.op_name}")

    return _sorry(f"unhandled expr: {type(e).__name__}")


def emit_call(e: Call, ctx: EmitCtx) -> str:
    """Call → Lean 调用表达式。callee 解析顺序：FUNC_MAP / LEAN_BUILTINS /
    LEAN_STDLIB / NOOP_FUNCS / TRANSLATION_SCOPE（加 _ir）/ 显式带 . 方法 /
    _lambda_/_loop_/已 _ir 后缀 / 兜底 sorry。"""
    callee = e.callee
    if isinstance(callee, UnresolvedOp):
        return _sorry(f"unresolved call: {callee.op_name}")
    if not isinstance(callee, str):
        return _sorry(f"non-str callee: {type(callee).__name__}")

    args_str = " ".join(_paren(emit_expr(a, ctx)) for a in e.args)
    no_args = not args_str.strip()

    # 查表（顺序与 v1 一致）
    try:
        from class_map import (
            LEAN_BUILTINS, TRANSLATION_SCOPE, FUNC_MAP,
            ASSERT_FAIL_NAMES,
        )
    except ImportError:
        LEAN_BUILTINS = set(); TRANSLATION_SCOPE = set()
        FUNC_MAP = {}; ASSERT_FAIL_NAMES = set()

    LEAN_STDLIB = {
        "Prod.mk": "Prod.mk",
        "Array.empty": "#[]",
        "Array.mk": "Array.mk",
        "Array.set!": "Array.set!",
        "Array.size_u64": "Array.size_u64",
        "Array.filter": "Array.filter",
        "Array.filterMap": "Array.filterMap",
        "Zp.mk": "Zp.mk",
        "UMonomial.mk": "UMonomial.mk",
        "Option.some": "Option.some",
        "Option.none": "Option.none",
    }
    NOOP_FUNCS = {"reserve"}

    # assert / unknown_func → ()
    if callee in ASSERT_FAIL_NAMES or callee == "unknown_func":
        return "()"

    # __ctor__<template> — Pass 5 constructor 解析后形式
    # 模板含 {a0}/{a1}/... 占位符，按位置 args 替换为实际表达式
    if callee.startswith("__ctor__"):
        template = callee[len("__ctor__"):]
        arg_strs = [emit_expr(a, ctx) for a in e.args]
        formatted = template
        for i, s in enumerate(arg_strs):
            formatted = formatted.replace(f"{{a{i}}}", s)
        # 模板可能是 `Zp.mk {a0} {a1}` 这种空格分隔形式，包圆括号
        return f"({formatted})"

    # __cast__<template> — Pass 5 cast 解析后形式
    # 模板含 {x} 占位符（仅 1 个 arg = inner expr）
    if callee.startswith("__cast__"):
        template = callee[len("__cast__"):]
        if len(e.args) == 1:
            inner = emit_expr(e.args[0], ctx)
            formatted = template.replace("{x}", inner)
            return f"({formatted})"
        # 多 arg 退化：保留 callee 名（不应发生）
        return _sorry(f"__cast__ with !=1 arg: {len(e.args)}")

    # _mutate_ 前缀 → 去掉（v1 兼容）
    raw = callee
    if callee.startswith("_mutate_"):
        callee = callee[len("_mutate_"):]

    # functional update {obj with field := value}
    if callee == "_with" and len(e.args) == 3:
        obj = emit_expr(e.args[0], ctx)
        field_arg = e.args[1]
        field = field_arg.name if isinstance(field_arg, Var) else emit_expr(field_arg, ctx)
        val = emit_expr(e.args[2], ctx)
        return f"{{ {obj} with {field} := {val} }}"

    if callee in LEAN_BUILTINS:
        return callee if no_args else f"({callee} {args_str})"
    if callee in FUNC_MAP:
        lean_name = FUNC_MAP[callee][0]
        return lean_name if no_args else f"({lean_name} {args_str})"
    _func_map_targets = {v[0] for v in FUNC_MAP.values()}
    if callee in _func_map_targets:
        return callee if no_args else f"({callee} {args_str})"
    if callee in LEAN_STDLIB:
        lean_name = LEAN_STDLIB[callee]
        return lean_name if no_args else f"({lean_name} {args_str})"
    if callee in NOOP_FUNCS:
        return args_str.split()[0] if args_str else "()"
    if callee in TRANSLATION_SCOPE:
        return f"{callee}_ir" if no_args else f"({callee}_ir {args_str})"
    # 已带 _ir 后缀（loop / lifted lambda）
    if callee.endswith("_ir"):
        return callee if no_args else f"({callee} {args_str})"
    # _lambda_<host>_N / _loop_<host>_<id> — Call.callee 来自 base_name（无 _ir），
    # MIRFunc.lean_name 给它加了 _ir，所以这里也要加。
    if callee.startswith(("_lambda_", "_loop_")):
        ir_name = f"{callee}_ir"
        return ir_name if no_args else f"({ir_name} {args_str})"
    # 含 "." 的方法名（CLASS_MAP 已映射）
    if "." in callee:
        return callee if no_args else f"({callee} {args_str})"
    # 普通 ident
    if callee.isidentifier():
        return callee if no_args else f"({callee} {args_str})"
    return _sorry(f"unknown callee: {raw}")


def emit_cast(e: Cast, ctx: EmitCtx) -> str:
    """Cast → Lean 类型转换表达式（查 CAST_TABLE，否则兜底 (e : T) 标注）。"""
    src, tgt = e.source_ty, e.target_ty
    inner = emit_expr(e.expr, ctx)
    # 同类型 / NoOp
    if src == tgt:
        return inner
    # 字面量 cast → 直接生成目标类型的字面量
    if isinstance(e.expr, Lit) and isinstance(tgt, BaseType):
        return f"({e.expr.value} : {tgt.value})"
    # BaseType → BaseType 查表
    if isinstance(src, BaseType) and isinstance(tgt, BaseType):
        tmpl = CAST_TABLE.get((src, tgt))
        if tmpl is not None:
            return tmpl.format(e=inner)
    # NamedType ↔ BaseType: 兜底类型标注（Lean Coe 实例自动处理）
    if isinstance(tgt, (BaseType, NamedType)):
        return f"({inner} : {emit_type(tgt)})"
    # 其他兜底
    return f"({inner} : {emit_type(tgt)})"


# ============================================================
# 语句 emit（仅 LetStmt / RequireStmt / PhiStmt——MIR 唯一允许）
# ============================================================

def emit_stmt(s: MIRStmt, ctx: EmitCtx) -> str:
    """MIRStmt → Lean let-binding 行（含 indent）。"""
    pad = "  " * ctx.indent

    if isinstance(s, LetStmt):
        ty_str = emit_type(s.ty)
        val_str = emit_expr(s.value, ctx)
        var_str = emit_var_name(s.var)
        # 若类型是 sorry 或 unknown 残留，省略类型标注让 Lean 推断
        if "sorry" in ty_str:
            return f"{pad}let {var_str} := {val_str}"
        return f"{pad}let {var_str} : {ty_str} := {val_str}"

    if isinstance(s, RequireStmt):
        # MIR 中 require 通常已被 hoist 到 func sig；body 内残留作注释
        cond = emit_expr(s.cond, ctx)
        return f"{pad}-- require ({s.name}): {cond}"

    if isinstance(s, PhiStmt):
        # PhiStmt 在 emit_cfg 中通过 merge BB lambda params 处理；
        # 不应单独 emit（S2 实现）。容错：注释占位
        var_str = emit_var_name(s.target)
        return f"{pad}-- phi {var_str} (handled at merge BB lambda)"

    return f"{pad}-- unhandled MIRStmt: {type(s).__name__}"


# ============================================================
# Param emit（用于 S3 emit_mirfunc）
# ============================================================

_anon_param_counter = [0]


def emit_param(p: HIRParam) -> str:
    """HIRParam → Lean 函数参数 `(name : T)`。空名给匿名占位（防 parse 错误）。"""
    name = p.name if p.name else f"_anon_{_anon_param_counter[0]}"
    if not p.name:
        _anon_param_counter[0] += 1
    return f"({_safe_ident(name)} : {emit_type(p.ty)})"


def emit_params(params: list[HIRParam]) -> str:
    """params → Lean 参数列表（空格分隔）。"""
    return " ".join(emit_param(p) for p in params) if params else ""


# ============================================================
# S2: CFG / 终止子 emit
# ============================================================

def emit_cfg(cfg: CFG, base_indent: int = 1) -> str:
    """MIRFunc.cfg → Lean 函数体表达式（不含 def 头/尾）。

    算法（见 pass8_codegen 文档 S2）：
      1. 找 merge BBs（preds > 1）
      2. 后序排序 merge BBs（被调用者 → 调用者）
      3. emit `let bb_<id> := fun phi_targets => ...` 链
      4. emit entry BB inline（CondJumpTerm 拆 if-then-else，跳到 merge → bb_<id> args）
    """
    # 1. find merge BBs
    merge_bbs: set[int] = {
        bb_id for bb_id in cfg.blocks
        if bb_id != cfg.entry and len(cfg.preds.get(bb_id, [])) > 1
    }

    # 2. postorder DFS（successor first，自身后；这样 callee 先于 caller emit）
    visited: set[int] = set()
    postorder: list[int] = []

    def dfs(bb_id: int):
        if bb_id in visited:
            return
        visited.add(bb_id)
        for s in cfg.succs.get(bb_id, []):
            dfs(s)
        postorder.append(bb_id)

    dfs(cfg.entry)

    # 3 + 4. 构造 ctx 给所有 emit 函数共享
    ctx = EmitCtx(indent=base_indent, cfg=cfg, merge_bbs=merge_bbs)

    # 3. emit merge BB lambdas（后序：被调用者先 emit）
    merge_lines: list[str] = []
    for bb_id in postorder:
        if bb_id in merge_bbs:
            merge_lines.append(emit_merge_lambda(bb_id, ctx))

    # 4. emit entry inline（最后一个表达式）
    entry_lines = emit_bb_inline(cfg.entry, ctx)

    # 拼接
    parts: list[str] = []
    for ml in merge_lines:
        parts.append(ml)
    parts.append(entry_lines)
    return "\n".join(parts)


def emit_merge_lambda(bb_id: int, ctx: EmitCtx) -> str:
    """合并块 → `let bb_<id> := fun phi_targets => <body>`。

    phi_targets 取 BB 开头连续 PhiStmt 的 target Var；body 跳过 phi 后 inline。
    """
    assert ctx.cfg is not None
    bb = ctx.cfg.blocks[bb_id]
    phi_targets: list[Var] = []
    for s in bb.stmts:
        if isinstance(s, PhiStmt):
            phi_targets.append(s.target)
        else:
            break
    pad = "  " * ctx.indent
    inner_pad = "  " * (ctx.indent + 1)
    if phi_targets:
        params_str = " ".join(emit_var_name(v) for v in phi_targets)
        header = f"{pad}let bb_{bb_id} := fun {params_str} =>"
    else:
        # 0 phi（理论可能：merge BB 无 SSA 冲突）→ 0 元 lambda
        header = f"{pad}let bb_{bb_id} := (fun _ : Unit =>"
    inner_ctx = EmitCtx(indent=ctx.indent + 1, cfg=ctx.cfg,
                          merge_bbs=ctx.merge_bbs)
    body = emit_bb_inline(bb_id, inner_ctx)
    if not phi_targets:
        return f"{header}\n{body}\n{pad}) ()"
    return f"{header}\n{body}"


def emit_bb_inline(bb_id: int, ctx: EmitCtx) -> str:
    """BB stmts（跳过开头 PhiStmt）+ terminator → Lean 表达式（多行字符串）。"""
    assert ctx.cfg is not None
    bb = ctx.cfg.blocks[bb_id]
    pad = "  " * ctx.indent
    lines: list[str] = []
    seen_non_phi = False
    for s in bb.stmts:
        if isinstance(s, PhiStmt) and not seen_non_phi:
            continue  # phi 在 merge lambda params 中处理
        seen_non_phi = True
        if isinstance(s, RequireStmt):
            # require 当前作注释（语义由调用方前置 prop 处理；S3 时再决定）
            lines.append(emit_stmt(s, ctx))
            continue
        if isinstance(s, LetStmt):
            lines.append(emit_stmt(s, ctx))
            continue
        # 不应到这（PhiStmt 在 phi 段后已 break；其它 MIRStmt 不允许）
        lines.append(f"{pad}-- unhandled stmt in MIR body: {type(s).__name__}")

    # terminator
    term_str = emit_terminator(bb.terminator, bb_id, ctx)
    lines.append(term_str)
    return "\n".join(lines)


def emit_terminator(t, src_bb_id: int, ctx: EmitCtx) -> str:
    """Terminator → Lean 表达式片段（多行；首行已含 ctx.indent 缩进）。"""
    pad = "  " * ctx.indent
    if isinstance(t, JumpTerm):
        return emit_jump_to(t.target, src_bb_id, ctx)
    if isinstance(t, CondJumpTerm):
        cond_str = emit_expr(t.cond, ctx)
        # 子分支缩进 +1
        inner_ctx = EmitCtx(indent=ctx.indent + 1, cfg=ctx.cfg,
                              merge_bbs=ctx.merge_bbs)
        then_str = emit_jump_to(t.then_bb, src_bb_id, inner_ctx)
        else_str = emit_jump_to(t.else_bb, src_bb_id, inner_ctx)
        return (f"{pad}if {cond_str} then\n"
                f"{then_str}\n"
                f"{pad}else\n"
                f"{else_str}")
    if isinstance(t, ReturnTerm):
        if t.value is None:
            return f"{pad}()"
        return f"{pad}{emit_expr(t.value, ctx)}"
    if isinstance(t, TailCallTerm):
        # target_func 是 base_name（如 `_loop_<host>_<id>`），需加 `_ir` 后缀
        callee = f"{t.target_func}_ir"
        if not t.args:
            return f"{pad}{callee}"
        args_str = " ".join(_paren(emit_expr(a, ctx)) for a in t.args)
        return f"{pad}{callee} {args_str}"
    return f"{pad}{_sorry(f'unhandled terminator: {type(t).__name__}')}"


# ============================================================
# S3: MIRFunc / 文件级 emit
# ============================================================

def emit_mirfunc(f: MIRFunc) -> str:
    """单个 MIRFunc → 顶层 `partial def {lean_name} (params) : ret := <body>`。"""
    sig_params = emit_params(f.params)
    sig_params_str = f" {sig_params}" if sig_params else ""
    ret_ty = emit_type(f.ret_ty)
    sig = f"partial def {f.lean_name}{sig_params_str} : {ret_ty} :="
    if f.cfg is None:
        # 罕见：无 cfg 的 MIRFunc（如纯 sorry stub）
        return f"{sig}\n  sorry /- no cfg -/"
    body = emit_cfg(f.cfg, base_indent=1)
    return f"{sig}\n{body}"


def codegen_pass(top: MIRFunc) -> str:
    """完整 mir 树 → Lean 4 .lean 源码（单 top）。

    布局：
      import CLPoly.Model
      namespace Generated
      <aux_defs[*]> 前置（按拓扑序：lifted lambda + loops 先于宿主）
      <top> 主 def
      end Generated
    """
    out: list[str] = []
    out.append("-- Auto-generated by cpp2lean v2 Pass 8")
    out.append("import CLPoly.Model")
    out.append("")
    out.append("namespace Generated")
    out.append("")
    for f in _topo_collect_funcs([top]):
        out.append(emit_mirfunc(f))
        out.append("")
    out.append("end Generated")
    return "\n".join(out)


def codegen_corpus(top_funcs: list[MIRFunc]) -> str:
    """全 corpus → 单一 .lean 源码（aggregate，含一个全局 `mutual ... end`）。

    所有 top 函数 + 它们的 aux_defs 全部在同一 mutual 块中——允许任意调用顺序，
    简化文件分块策略（B1 阶段不做精细 SCC 拆分）。

    输出：单个 .lean 文件字符串。
    """
    out: list[str] = []
    out.append("-- Auto-generated by cpp2lean v2 Pass 8 (corpus aggregate)")
    out.append("import CLPoly.Model")
    out.append("")
    out.append("namespace Generated")
    out.append("")
    out.append("mutual")
    funcs = _topo_collect_funcs(top_funcs)
    for f in funcs:
        out.append(emit_mirfunc(f))
        out.append("")
    out.append("end")  # mutual end
    out.append("")
    out.append("end Generated")
    return "\n".join(out)


def _topo_collect_funcs(roots: list[MIRFunc]) -> list[MIRFunc]:
    """从 roots 列表收集所有 MIRFunc，按"被依赖者先"顺序返回。

    对每个 root：DFS 走 aux_defs（先递归到叶子），后序追加。跨 root 共享
    visited 集，按 id() 去重。
    """
    emit_order: list[MIRFunc] = []
    visited: set[int] = set()

    def _walk(f: MIRFunc):
        if id(f) in visited:
            return
        visited.add(id(f))
        for a in f.aux_defs:
            _walk(a)
        emit_order.append(f)

    for r in roots:
        _walk(r)
    return emit_order


def emit_jump_to(target_bb: int, src_bb_id: int, ctx: EmitCtx) -> str:
    """跳转到 target_bb：
      - target ∈ merge_bbs → `bb_<target> <phi_args from src>`
      - target 单前驱 → inline emit_bb_inline(target)
    """
    assert ctx.cfg is not None
    pad = "  " * ctx.indent
    if target_bb in ctx.merge_bbs:
        target_block = ctx.cfg.blocks[target_bb]
        phi_args: list[str] = []
        for s in target_block.stmts:
            if isinstance(s, PhiStmt):
                src_var = s.sources.get(src_bb_id)
                if src_var is None:
                    phi_args.append(_sorry(f"missing phi src from bb[{src_bb_id}]"))
                else:
                    phi_args.append(_paren(emit_expr(src_var, ctx)))
            else:
                break
        if phi_args:
            return f"{pad}bb_{target_bb} " + " ".join(phi_args)
        # 0-phi merge：调零元 lambda（与 emit_merge_lambda 的 `(fun _ => ...) ()` 对应）
        return f"{pad}bb_{target_bb}"
    # 单前驱 → inline
    return emit_bb_inline(target_bb, ctx)

