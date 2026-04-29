"""
Pass 2b: callsite_ref_elim — HIR₁ → HIR₁'

调用点改写：把 `f(out, in)` 形式的 ExprStmt / LetStmt / AssignStmt 中的
mutating call 改写为 destructure，使输出参数的修改通过 SSA 链传播。

`mutation-model-design.md` §3.2 已设计变换规则；本 Pass 落实之。
查询表是 `class_map.TRANSLATION_SCOPE_OUTPUT_PARAMS` + `get_output_params`。

变换规则：
  void + 1 out:
    `(expr) f(out, in)`
      ↓
    `out := f(out, in)`
  void + 2 out:
    `(expr) f(out1, in, out2)`
      ↓
    `let __refret_N := f(out1, in, out2)
     out1 := __refret_N.fst
     out2 := __refret_N.snd`
  void + 3+ out:
    `let __refret_N := f(args)
     out_k := __refret_N.elem<k>` (k = 0,1,...)
  非 void + N out:
    *暂不支持*（ExprStmt 形式 = 调用结果丢弃；LetStmt/AssignStmt 捕获返回
    较罕见，留 follow-up）

工程：仅处理 ExprStmt(Call(callee=str, args=...)) 形式。callee 必须是
str（自由函数）；method/operator/UnresolvedOp 不在此 Pass 范围。
"""

from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR,
    PairType, TupleType,
    Var, Lit, Call, FieldAccess, ArrayAccess, Cast, UnresolvedOp,
    LetStmt, AssignStmt, IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
    BreakStmt, ContinueStmt, ReturnStmt, RequireStmt, ExprStmt, BlockStmt,
    HIRFunc, StmtIR, ExprIR, TranslationError,
)
from class_map import get_output_params


_FIELD_NAMES = ["fst", "snd"]
_ELEM_NAMES = lambda i: f"elem{i}"  # 3+ 元素时改用 elem0/elem1/...

# 透明 lvalue 包装（C++ 中返回内部引用的 method）：
#   `vec.data()` 返回内部数组的引用，写到该引用 = 写 vec 自己
#   `obj.first` / `obj.second` 已经是 FieldAccess，不在此列
# 这些 method 在 HIR1 阶段未被 Pass 5 解析，仍是 Call(UnresolvedOp("<method>.X"))
_TRANSPARENT_METHODS = {"data", "ref"}


def _unwrap_to_lvalue(e: ExprIR) -> ExprIR | None:
    """剥离 Cast(NoOp) + 透明 method call + UnresolvedOp("operator[]")，
    返回底层 lvalue（Var/FieldAccess/ArrayAccess）。

    HIR1 阶段（Pass 5 之前）大部分 lvalue 表达式仍是 Call/UnresolvedOp 形式：
      `vec.data()`     → Call(UnresolvedOp("<method>.data"), [vec])
      `arr[i]`         → Call(UnresolvedOp("operator[]"), [arr, i])
      `(T&)x`          → Cast(NoOp, x)
    需要识别这些形态并转换为目标 IR 中的对应 lvalue。
    """
    cur = e
    while True:
        if isinstance(cur, Cast) and cur.cast_kind in ("NoOp", "LValueToRValue", None):
            cur = cur.expr
            continue
        if isinstance(cur, Call) and isinstance(cur.callee, UnresolvedOp):
            op = cur.callee.op_name
            # `obj.data()` 等透明 method
            if op.startswith("<method>.") and len(cur.args) == 1:
                method_name = op[len("<method>."):]
                if method_name in _TRANSPARENT_METHODS:
                    cur = cur.args[0]
                    continue
            # `arr[i]` —— 重塑为 ArrayAccess
            if op == "operator[]" and len(cur.args) == 2:
                arr_inner = _unwrap_to_lvalue(cur.args[0])
                if arr_inner is None:
                    return None
                return ArrayAccess(arr=arr_inner, idx=cur.args[1], ty=cur.ty)
        break
    if isinstance(cur, (Var, FieldAccess, ArrayAccess)):
        return cur
    return None


def _is_lvalue(e: ExprIR) -> bool:
    """LHS 友好检查：能 unwrap 到 Var/FieldAccess/ArrayAccess 即合法。"""
    return _unwrap_to_lvalue(e) is not None


def _ref_elim_call_at_stmt(call: Call, counter: list[int]) -> list[StmtIR] | None:
    """对 Call 节点尝试做 ref-elim 改写。返回 None 表示不需改写。

    counter: 单元素 list 用作 mutable counter（生成 __refret_N 唯一名）。
    """
    callee = call.callee
    if not isinstance(callee, str):
        return None
    out_indices = get_output_params(callee, len(call.args))
    if not out_indices:
        return None

    # 提取 out 参数列表（必须是 lvalue 才能赋值）
    # 透明包装如 `f_new.data()` 应剥离到底层 Var
    out_args: list[ExprIR] = []
    for idx in out_indices:
        if idx >= len(call.args):
            return None  # 索引越界——保险跳过
        unwrapped = _unwrap_to_lvalue(call.args[idx])
        if unwrapped is None:
            # rvalue out（如字面量）— 应该是 bug，跳过此 callsite
            return None
        out_args.append(unwrapped)

    n_out = len(out_indices)
    if n_out == 1:
        # 单输出：直接 `out := f(args)`，不需要临时变量
        return [AssignStmt(target=out_args[0], value=call)]

    # 多输出：临时变量 + 字段提取
    tmp_id = counter[0]
    counter[0] += 1
    tmp_name = f"__refret_{tmp_id}"
    tmp_var = Var(name=tmp_name, ty=UnknownType(""))
    out: list[StmtIR] = [LetStmt(var=tmp_var, ty=UnknownType(""), value=call)]
    field_names = _FIELD_NAMES if n_out == 2 else [_ELEM_NAMES(k) for k in range(n_out)]
    for k, target in enumerate(out_args):
        fld = field_names[k] if k < len(field_names) else _ELEM_NAMES(k)
        out.append(AssignStmt(
            target=target,
            value=FieldAccess(obj=tmp_var, field_name=fld, ty=UnknownType("")),
        ))
    return out


def _rewrite_stmt(s: StmtIR, counter: list[int]) -> list[StmtIR]:
    """处理一条 stmt，返回替换后的 stmt 列表。递归 if/while/for/block。"""
    if isinstance(s, ExprStmt):
        if isinstance(s.expr, Call):
            replaced = _ref_elim_call_at_stmt(s.expr, counter)
            if replaced is not None:
                return replaced
        return [s]

    if isinstance(s, IfStmt):
        return [replace(s,
                        then_body=_rewrite_stmts(list(s.then_body), counter),
                        else_body=_rewrite_stmts(list(s.else_body or []), counter))]
    if isinstance(s, WhileStmt):
        return [replace(s, body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, ForStmt):
        return [replace(s,
                        init=_rewrite_stmts(list(s.init), counter),
                        step=_rewrite_stmts(list(s.step), counter),
                        body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, RangeForStmt):
        return [replace(s, body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, DoWhileStmt):
        return [replace(s, body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, BlockStmt):
        return [replace(s, stmts=_rewrite_stmts(list(s.stmts), counter))]

    # LetStmt / AssignStmt / RequireStmt / ReturnStmt / Break / Continue：
    # 当前不改写其内嵌 Call（参见 docstring 设计：留 follow-up）
    return [s]


def _rewrite_stmts(stmts: list[StmtIR], counter: list[int]) -> list[StmtIR]:
    out: list[StmtIR] = []
    for s in stmts:
        out.extend(_rewrite_stmt(s, counter))
    return out


def callsite_ref_elim_pass(func: HIRFunc) -> HIRFunc:
    """HIR₁ → HIR₁'：调用点 ref-elim 改写。"""
    counter = [0]
    new_body = _rewrite_stmts(list(func.body), counter)
    new_aux = []
    for aux in func.aux_lambdas:
        new_aux.append(callsite_ref_elim_pass(aux))
    return replace(func, body=new_body, aux_lambdas=new_aux)
