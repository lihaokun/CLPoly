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
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType, RefType,
    PairType, TupleType,
    Var, Lit, Call, FieldAccess, ArrayAccess, Cast, UnresolvedOp, UnaryOp,
    LetStmt, AssignStmt, IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
    BreakStmt, ContinueStmt, ReturnStmt, RequireStmt, ExprStmt, BlockStmt,
    HIRFunc, StmtIR, ExprIR, TranslationError,
)
from class_map import get_output_params, is_callee_nonvoid


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
    nonvoid = is_callee_nonvoid(callee, len(call.args))

    # callee 字段命名约定（与 Pass 2 ref_elim 配套）：
    #   void + 1 out → ret_ty = ref_type (无 wrap), 单变量直接接收
    #   void + 2 out → ret_ty = PairType(refs)  → fst, snd
    #   void + 3+ out → ret_ty = TupleType(refs) → elem0..N-1
    #   non-void + 1 out → ret_ty = PairType(orig_ret, ref) → fst (ret), snd (ref)
    #   non-void + 2+ out → ret_ty = TupleType(orig_ret, refs) → elem0 (ret), elem1..N
    if not nonvoid and n_out == 1:
        # void + 1 out：直接 `out := f(args)`
        return [AssignStmt(target=out_args[0], value=call)]

    # 否则用 tmp + destructure
    tmp_id = counter[0]
    counter[0] += 1
    tmp_name = f"__refret_{tmp_id}"
    # 阶段 B：tmp_var ty 从 call.ty + out_args[i].ty 推导
    out_arg_tys = [_strip_ref_ty_b2(getattr(a, 'ty', None)) for a in out_args]
    tmp_ty = _build_refret_tuple_ty(call.ty if nonvoid else None, out_arg_tys)
    tmp_var = Var(name=tmp_name, ty=tmp_ty)
    # 阶段 F #8：polynomial_GCD#4 (Bezout 形式) 与 2-arg 形式同名，Lean 端
    # 区分名为 polynomial_GCD_eea。在 callsite refret 转换时同步 rename。
    refret_call = call
    if callee == "polynomial_GCD" and len(call.args) == 4:
        refret_call = Call(callee="polynomial_GCD_eea", args=call.args, ty=call.ty)
    elif callee == "pair_vec_div" and len(call.args) == 5:
        # 5-arg overload: (new_v, R&, v1, v2, comp) — Lean 端 pair_vec_div5 返回 (q, R) tuple
        refret_call = Call(callee="pair_vec_div5", args=call.args, ty=call.ty)
    elif callee == "poly_convert" and len(call.args) == 3:
        refret_call = Call(callee="poly_convert3", args=call.args, ty=call.ty)
    out: list[StmtIR] = [LetStmt(var=tmp_var, ty=tmp_ty, value=refret_call)]

    if nonvoid:
        n_total = 1 + n_out
        if n_total == 2:
            # fst = orig_ret (discarded), snd = ref
            out_fields = ["snd"]
        else:
            # 阶段 F #8 修复：n>2 tuple 用 Lean 嵌套 Prod 投影路径
            # tuple = (orig_ret, out0, out1, ..., outN-1) = (a, (b, (c, ..., d)))
            #   orig_ret → .1
            #   out_k (k=0..n_out-1) → tuple position k+1
            out_fields = []
            for k in range(n_out):
                pos = k + 1  # 在 tuple 中的位置（0=orig_ret）
                n_dots = pos
                if pos == n_total - 1:
                    # 末尾元素：纯 .2 链
                    out_fields.append(".".join(["2"] * n_dots))
                else:
                    # 中间元素：.2 链 + .1
                    out_fields.append(".".join(["2"] * n_dots) + ".1")
    else:
        if n_out == 2:
            out_fields = list(_FIELD_NAMES)
        elif n_out >= 3:
            # 阶段 F #8 修复：void + n_out>=3 tuple 同样用嵌套 Prod 投影
            out_fields = []
            for k in range(n_out):
                n_dots = k
                if k == n_out - 1:
                    out_fields.append(".".join(["2"] * n_dots) if n_dots > 0 else "1")
                elif k == 0:
                    out_fields.append("1")
                else:
                    out_fields.append(".".join(["2"] * n_dots) + ".1")
        else:
            out_fields = [_ELEM_NAMES(k) for k in range(n_out)]

    for k, target in enumerate(out_args):
        out.append(AssignStmt(
            target=target,
            value=FieldAccess(obj=tmp_var, field_name=out_fields[k],
                              ty=out_arg_tys[k] if out_arg_tys[k] else UnknownType("")),
        ))
    return out


def _strip_ref_ty_b2(ty):
    """剥 RefType 外壳（Pass 2b 内部使用，避免 forward dep on Pass 7）。"""
    if ty is None:
        return None
    if isinstance(ty, RefType):
        return _strip_ref_ty_b2(ty.inner)
    return ty


def _build_refret_tuple_ty(orig_ret_ty, out_arg_tys):
    """构造 __refret tmp 的 tuple 类型。
    nonvoid + 1 out: PairType(orig_ret, out0)
    nonvoid + 2+ out: TupleType(orig_ret, out0, out1, ...)
    void + 1 out: out0
    void + 2 out: PairType(out0, out1)
    void + 3+ out: TupleType(out0, out1, ...)
    None 元素 → UnknownType
    """
    def _safe(ty):
        return ty if ty is not None else UnknownType("")
    n_out = len(out_arg_tys)
    if orig_ret_ty is not None:
        # nonvoid
        elems = [_safe(orig_ret_ty)] + [_safe(t) for t in out_arg_tys]
    else:
        # void
        if n_out == 1:
            return _safe(out_arg_tys[0])
        elems = [_safe(t) for t in out_arg_tys]
    if len(elems) == 2:
        return PairType(fst=elems[0], snd=elems[1])
    if len(elems) >= 3:
        return TupleType(elems=tuple(elems))
    return UnknownType("")


def _hoist_ref_call(expr: ExprIR, counter: list[int]
                    ) -> tuple[ExprIR, list[StmtIR]] | None:
    """若 expr 是 ref-out Call（含 UnaryOp("!", Call) 包装），生成 hoist：
       前置 stmts: `let __refret_N := f(args); <out_i> := __refret_N.<field>`
       替换 expr：`__refret_N.fst`（或对 UnaryOp("!", ...) 包装：`!__refret_N.fst`）
    返回 None 表示 expr 不是需 hoist 的形态。
    """
    # 剥 UnaryOp("!") 包装
    not_wrap = False
    inner = expr
    if isinstance(expr, UnaryOp) and expr.op == "!":
        inner = expr.operand
        not_wrap = True
    if not isinstance(inner, Call):
        return None
    callee = inner.callee
    if not isinstance(callee, str):
        return None
    out_indices = get_output_params(callee, len(inner.args))
    if not out_indices:
        return None

    out_args: list[ExprIR] = []
    for idx in out_indices:
        if idx >= len(inner.args):
            return None
        unwrapped = _unwrap_to_lvalue(inner.args[idx])
        if unwrapped is None:
            return None
        out_args.append(unwrapped)

    tmp_id = counter[0]
    counter[0] += 1
    tmp_name = f"__refret_{tmp_id}"
    # 阶段 B：tmp_var ty 推导（_hoist_ref_call 总是 nonvoid，因 expr 即 inner ret）
    out_arg_tys = [_strip_ref_ty_b2(getattr(a, 'ty', None)) for a in out_args]
    tmp_ty = _build_refret_tuple_ty(inner.ty, out_arg_tys)
    tmp_var = Var(name=tmp_name, ty=tmp_ty)

    n_out = len(out_indices)
    # 字段命名：返回 (orig_ret, out_0, out_1, ...)
    n_total = 1 + n_out  # orig_ret + outs
    if n_total == 2:
        ret_field = "fst"
        out_fields = ["snd"]
    else:
        # 阶段 F #8 修复：n>2 tuple 用 Lean 嵌套 Prod 投影路径
        ret_field = "1"
        out_fields = []
        for k in range(n_out):
            pos = k + 1
            n_dots = pos
            if pos == n_total - 1:
                out_fields.append(".".join(["2"] * n_dots))
            else:
                out_fields.append(".".join(["2"] * n_dots) + ".1")

    # 阶段 F #8：polynomial_GCD#4 → polynomial_GCD_eea（同 _ref_elim_call_at_stmt）
    refret_inner = inner
    if isinstance(inner.callee, str) and inner.callee == "polynomial_GCD" \
            and len(inner.args) == 4:
        refret_inner = Call(callee="polynomial_GCD_eea", args=inner.args, ty=inner.ty)
    elif isinstance(inner.callee, str) and inner.callee == "pair_vec_div" \
            and len(inner.args) == 5:
        refret_inner = Call(callee="pair_vec_div5", args=inner.args, ty=inner.ty)
    elif isinstance(inner.callee, str) and inner.callee == "poly_convert" \
            and len(inner.args) == 3:
        refret_inner = Call(callee="poly_convert3", args=inner.args, ty=inner.ty)
    pre_stmts: list[StmtIR] = [
        LetStmt(var=tmp_var, ty=tmp_ty, value=refret_inner),
    ]
    for k, target in enumerate(out_args):
        pre_stmts.append(AssignStmt(
            target=target,
            value=FieldAccess(obj=tmp_var, field_name=out_fields[k],
                              ty=out_arg_tys[k] if out_arg_tys[k] else UnknownType("")),
        ))
    new_inner = FieldAccess(obj=tmp_var, field_name=ret_field, ty=inner.ty)
    new_expr = UnaryOp(op="!", operand=new_inner, ty=expr.ty) if not_wrap else new_inner
    return new_expr, pre_stmts


def _rewrite_stmt(s: StmtIR, counter: list[int]) -> list[StmtIR]:
    """处理一条 stmt，返回替换后的 stmt 列表。递归 if/while/for/block。"""
    if isinstance(s, ExprStmt):
        if isinstance(s.expr, Call):
            replaced = _ref_elim_call_at_stmt(s.expr, counter)
            if replaced is not None:
                return replaced
        return [s]

    # P0 修复（agent 第七轮发现）：嵌入在 LetStmt/AssignStmt RHS 与 IfStmt/While
    # cond 内的 ref-out Call 也要 hoist + destructure，否则 out 参数 SSA 不 bump
    # + Pair/Tuple 返回值被当 Bool/Int 处理。
    if isinstance(s, LetStmt) and isinstance(s.value, Call):
        hr = _hoist_ref_call(s.value, counter)
        if hr is not None:
            new_value, pre = hr
            return pre + [LetStmt(var=s.var, ty=s.ty, value=new_value)]

    if isinstance(s, AssignStmt) and isinstance(s.value, Call):
        hr = _hoist_ref_call(s.value, counter)
        if hr is not None:
            new_value, pre = hr
            return pre + [AssignStmt(target=s.target, value=new_value)]

    if isinstance(s, IfStmt):
        new_then = _rewrite_stmts(list(s.then_body), counter)
        new_else = _rewrite_stmts(list(s.else_body or []), counter)
        hr = _hoist_ref_call(s.cond, counter)
        if hr is not None:
            new_cond, pre = hr
            return pre + [replace(s, cond=new_cond, then_body=new_then,
                                  else_body=new_else)]
        return [replace(s, then_body=new_then, else_body=new_else)]
    if isinstance(s, WhileStmt):
        new_body = _rewrite_stmts(list(s.body), counter)
        hr = _hoist_ref_call(s.cond, counter)
        if hr is not None:
            new_cond, pre = hr
            return pre + [replace(s, cond=new_cond, body=new_body)]
        return [replace(s, body=new_body)]
    if isinstance(s, DoWhileStmt):
        new_body = _rewrite_stmts(list(s.body), counter)
        hr = _hoist_ref_call(s.cond, counter)
        if hr is not None:
            new_cond, pre = hr
            # DoWhile cond 在 body 之后求值；hoist 放 body 末尾 + continue 路径
            # （此 Pass 不识别 continue —— 简化：放 body 末尾，假设无 continue）
            return [replace(s, body=new_body + pre, cond=new_cond)]
        return [replace(s, body=new_body)]
    if isinstance(s, ForStmt):
        return [replace(s,
                        init=_rewrite_stmts(list(s.init), counter),
                        step=_rewrite_stmts(list(s.step), counter),
                        body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, RangeForStmt):
        return [replace(s, body=_rewrite_stmts(list(s.body), counter))]
    if isinstance(s, BlockStmt):
        return [replace(s, stmts=_rewrite_stmts(list(s.stmts), counter))]

    # 其它 stmt（RequireStmt/ReturnStmt/Break/Continue/CompoundAssignStmt）：
    # 嵌套 Call hoist 暂不覆盖（corpus 未观测；P3-dormant follow-up）
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
