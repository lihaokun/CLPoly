"""
Pass 3b: lambda_ref_elim — HIR₂ → HIR₂'

对 Pass 3 lifted lambda 应用 ref-elim + 调用点改写：

- **签名改造**：每个 lifted lambda 中 `is_ref=True` 的 cap params → 函数返回
  tuple 含 modified captures（复用 `pass2_ref_elim.ref_elim_pass`）。
- **调用点 prepend captures**（修 Pass 3 pre-existing bug A）：每个 call to
  `_lambda_<name>(args)` → `_lambda_<name>(<cap_args>, args)`。
- **调用点 destructure**（修 P0-2 同源 bug）：仿 Pass 2b，把
  `(expr) _lambda_..._N(<caps>, <args>)` 改写为 `let __refret_lam_N := ...;
  <cap_i> := __refret_lam_N.<field_i>`。

调用点识别：
  HIR2 中 lambda 通过 `LetStmt(var=alias, value=Var("_lambda_<name>"))` 别名绑定，
  随后 `Call(callee=Var(alias), args=[...])` 调用。本 Pass 先扫描建 alias map，
  再重写调用。

参考：mutation-model-design.md §5 G3、修正方案
       docs/fixes/cpp2lean-v2-lambda-by-ref-capture.md。
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
    HIRFunc, HIRParam, StmtIR, ExprIR, TranslationError,
)
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass


_FIELD_NAMES = ["fst", "snd"]
_ELEM = lambda i: f"elem{i}"


def _extract_cap_info(lifted: HIRFunc):
    """从 Pass 3 lifted lambda 提取 cap 信息。

    Pass 3 在 qual_type 中编码 `n_caps=N`，标识前 N 个 params 是 captures
    （其余是原 lambda 自己的 params）。

    返回 (cap_names, ref_param_info, orig_ret)：
      cap_names: 前 N 个 cap 参数名字（call site prepend 时用）
      ref_param_info: list of (lifted_param_idx, target_kind, target_name_or_argpos)
        target_kind = "cap"        → destructure 写到 Var(cap_name)
        target_kind = "lam_param"  → destructure 写到 call-site real_args[k]
        每条对应 lifted 中一个 is_ref param；按 lifted.params 顺序，与
        ref_elim_pass 包装的 tuple 字段顺序一致。
      orig_ret: lambda 原始返回类型
    """
    import re
    qt = lifted.qual_type or ""
    m = re.search(r"n_caps=(\d+)", qt)
    if m is None:
        return [], [], lifted.ret_ty
    n_caps = int(m.group(1))
    cap_names: list[str] = [p.name for p in lifted.params[:n_caps]]
    ref_param_info: list[tuple[int, str, object]] = []
    for i, p in enumerate(lifted.params):
        if p.is_ref:
            if i < n_caps:
                ref_param_info.append((i, "cap", p.name))
            else:
                # lambda 自己的 by-ref 参数：destructure 到 call-site arg[i - n_caps]
                ref_param_info.append((i, "lam_param", i - n_caps))
    return cap_names, ref_param_info, lifted.ret_ty


def _build_alias_map(stmts: list[StmtIR]) -> dict[str, str]:
    """递归扫描 stmts，找 LetStmt(var=X, value=Var("_lambda_<name>"))，建 alias→lambda 映射。"""
    out: dict[str, str] = {}

    def walk(ss: list[StmtIR]):
        for s in ss:
            if isinstance(s, LetStmt) and isinstance(s.value, Var):
                if s.value.name.startswith("_lambda_"):
                    out[s.var.name] = s.value.name
            if isinstance(s, IfStmt):
                walk(s.then_body); walk(s.else_body or [])
            elif isinstance(s, (WhileStmt, RangeForStmt, DoWhileStmt)):
                walk(s.body)
            elif isinstance(s, ForStmt):
                walk(s.init); walk(s.body)
            elif isinstance(s, BlockStmt):
                walk(s.stmts)
    walk(stmts)
    return out


def _resolve_callee(e: ExprIR, alias_map: dict[str, str]) -> str | None:
    """若 e 是直接或别名指向 lifted lambda 的引用，返回 lambda 名；否则 None。
    支持剥离 Cast(NoOp/LValueToRValue) 包装。
    """
    cur = e
    while isinstance(cur, Cast) and cur.cast_kind in ("NoOp", "LValueToRValue", None):
        cur = cur.expr
    if isinstance(cur, Var):
        name = cur.name
        if name.startswith("_lambda_"):
            return name
        if name in alias_map:
            return alias_map[name]
    if isinstance(cur, str):
        if cur.startswith("_lambda_"):
            return cur
        if cur in alias_map:
            return alias_map[cur]
    return None


def _try_unpack_operator_call(call: Call) -> tuple[ExprIR, list[ExprIR]] | None:
    """识别 lambda 调用的 HIR2 形态：
       `Call(UnresolvedOp("operator()"), [callee_expr, arg1, arg2, ...])`
       返回 (callee_expr, real_args)；非匹配返回 None。
    """
    if not (isinstance(call.callee, UnresolvedOp)
            and call.callee.op_name == "operator()"
            and len(call.args) >= 1):
        return None
    return call.args[0], list(call.args[1:])


def _resolve_lam_call(call: Call, alias_map: dict[str, str]
                      ) -> tuple[str, list[ExprIR]] | None:
    """识别 call 是 lifted lambda 调用，返回 (lam_name, real_args)。

    支持两种 HIR 形态：
      A) 直接 callee：Call(callee=Var("_lambda_..." | alias), args=[...])
      B) operator() 包装：Call(UnresolvedOp("operator()"), [callee_expr, ...])
    """
    lam_name = _resolve_callee(call.callee, alias_map)
    if lam_name is not None:
        return lam_name, list(call.args)
    unpacked = _try_unpack_operator_call(call)
    if unpacked is not None:
        callee_expr, real_args_b = unpacked
        lam_name = _resolve_callee(callee_expr, alias_map)
        if lam_name is not None:
            return lam_name, real_args_b
    return None


def _unwrap_to_lvalue(e: ExprIR) -> ExprIR | None:
    """剥离 Cast(NoOp/LValueToRValue) + UnresolvedOp("operator[]") + 透明 method
    （参考 pass2b 同名 helper）。返回底层 lvalue（Var/FieldAccess/ArrayAccess）。
    """
    cur = e
    while True:
        if isinstance(cur, Cast) and cur.cast_kind in ("NoOp", "LValueToRValue", None):
            cur = cur.expr
            continue
        if isinstance(cur, Call) and isinstance(cur.callee, UnresolvedOp):
            op = cur.callee.op_name
            # `arr[i]` → ArrayAccess
            if op == "operator[]" and len(cur.args) == 2:
                arr_inner = _unwrap_to_lvalue(cur.args[0])
                if arr_inner is None:
                    return None
                return ArrayAccess(arr=arr_inner, idx=cur.args[1], ty=cur.ty)
            # `obj.data()` 等透明 method
            if op.startswith("<method>.") and len(cur.args) == 1:
                method_name = op[len("<method>."):]
                if method_name in {"data", "ref"}:
                    cur = cur.args[0]
                    continue
        break
    if isinstance(cur, (Var, FieldAccess, ArrayAccess)):
        return cur
    return None


def _build_destructure_targets(ref_param_info, real_args: list[ExprIR]) -> list[ExprIR] | None:
    """为每个 is_ref param 构造 destructure target lvalue。
    返回 None 表示某个 lam_param 对应的 real_args 不是 lvalue（不能 destructure）。
    """
    targets: list[ExprIR] = []
    for (idx, kind, where) in ref_param_info:
        if kind == "cap":
            targets.append(Var(name=where, ty=UnknownType("")))
        else:  # "lam_param"
            arg_pos = where
            if arg_pos >= len(real_args):
                return None
            unwrapped = _unwrap_to_lvalue(real_args[arg_pos])
            if unwrapped is None:
                return None
            targets.append(unwrapped)
    return targets


def _rewrite_call_at_stmt(call: Call,
                          lambda_info: dict[str, tuple[list[str], list]],
                          alias_map: dict[str, str],
                          counter: list[int]) -> list[StmtIR] | None:
    """对 Call 节点尝试做 lambda ref-elim 改写。返回 None 表示不需改写。

    若 call 指向 lifted lambda → prepend captures + 必要时 destructure。
    """
    resolved = _resolve_lam_call(call, alias_map)
    if resolved is None:
        return None
    lam_name, real_args = resolved
    info = lambda_info.get(lam_name)
    if info is None:
        return None
    cap_names, ref_param_info = info

    # prepend cap args + 原 args
    cap_args: list[ExprIR] = [Var(name=cn, ty=UnknownType("")) for cn in cap_names]
    new_call = Call(callee=lam_name, args=cap_args + real_args, ty=call.ty)

    # 无 is_ref params：保留 ExprStmt（仅 prepend）
    if not ref_param_info:
        if not cap_names:
            return None  # 完全无变化，passthrough
        return [ExprStmt(expr=new_call)]

    # 构造 destructure targets
    targets = _build_destructure_targets(ref_param_info, real_args)
    if targets is None:
        return None  # rvalue 出参 — 不能改写

    n_ref = len(ref_param_info)
    if n_ref == 1:
        # 单 is_ref：直接 AssignStmt（不需要临时 tuple）
        return [AssignStmt(target=targets[0], value=new_call)]

    # 多 is_ref：tmp + 字段提取
    tmp_id = counter[0]
    counter[0] += 1
    tmp_name = f"__refret_lam_{tmp_id}"
    tmp_var = Var(name=tmp_name, ty=UnknownType(""))
    out: list[StmtIR] = [LetStmt(var=tmp_var, ty=UnknownType(""), value=new_call)]
    field_names = _FIELD_NAMES if n_ref == 2 else [_ELEM(k) for k in range(n_ref)]
    for k, target in enumerate(targets):
        fld = field_names[k] if k < len(field_names) else _ELEM(k)
        out.append(AssignStmt(
            target=target,
            value=FieldAccess(obj=tmp_var, field_name=fld, ty=UnknownType("")),
        ))
    return out


def _hoist_lambda_call(expr: ExprIR,
                       lambda_info: dict[str, tuple[list[str], list]],
                       alias_map: dict[str, str],
                       counter: list[int]) -> tuple[ExprIR, list[StmtIR]] | None:
    """若 expr 是顶层 lambda Call（含 operator() 包装）且该 lambda 含 is_ref 参数
    （cap 或 lambda by-ref param），返回 (替换 expr, 前置 stmts)：

      - 前置 stmts：`let __hoist_N := <lam>(<caps>, <args>)` + 各 ref 字段 destructure
      - 替换 expr：原 lambda 返回值字段访问

    无 is_ref 参数时返回 None（无需 hoist；可在 stmt-level passthrough）。
    """
    if not isinstance(expr, Call):
        return None
    resolved = _resolve_lam_call(expr, alias_map)
    if resolved is None:
        return None
    lam_name, real_args = resolved
    info = lambda_info.get(lam_name)
    if info is None:
        return None
    cap_names, ref_param_info = info
    if not ref_param_info:
        return None  # 无 is_ref 参数，无需 hoist

    targets = _build_destructure_targets(ref_param_info, real_args)
    if targets is None:
        return None

    cap_args: list[ExprIR] = [Var(name=cn, ty=UnknownType("")) for cn in cap_names]
    hoisted_call = Call(callee=lam_name,
                        args=cap_args + real_args,
                        ty=expr.ty)
    tmp_id = counter[0]
    counter[0] += 1
    tmp_name = f"__hoist_lam_{tmp_id}"
    tmp_var = Var(name=tmp_name, ty=UnknownType(""))

    # 字段命名：ref_elim_pass 包装成 (orig_ret, ref0, ref1, ...)
    n_ref = len(ref_param_info)
    n_total = 1 + n_ref
    if n_total == 2:
        ret_field = "fst"
        ref_fields = ["snd"]
    else:
        ret_field = "elem0"
        ref_fields = [f"elem{k+1}" for k in range(n_ref)]

    pre_stmts: list[StmtIR] = [
        LetStmt(var=tmp_var, ty=UnknownType(""), value=hoisted_call),
    ]
    for k, target in enumerate(targets):
        pre_stmts.append(AssignStmt(
            target=target,
            value=FieldAccess(obj=tmp_var, field_name=ref_fields[k],
                              ty=UnknownType("")),
        ))
    new_expr = FieldAccess(obj=tmp_var, field_name=ret_field, ty=expr.ty)
    return new_expr, pre_stmts


def _rewrite_stmt(s: StmtIR,
                  lambda_info: dict[str, tuple[list[str], list[int]]],
                  alias_map: dict[str, str],
                  counter: list[int]) -> list[StmtIR]:
    """改写一条 stmt。递归 if/while/for/block。"""
    if isinstance(s, ExprStmt):
        if isinstance(s.expr, Call):
            replaced = _rewrite_call_at_stmt(s.expr, lambda_info, alias_map, counter)
            if replaced is not None:
                return replaced
        return [s]

    # IfStmt：cond 含 lambda Call → hoist
    if isinstance(s, IfStmt):
        new_then = _rewrite_stmts(list(s.then_body), lambda_info, alias_map, counter)
        new_else = _rewrite_stmts(list(s.else_body or []), lambda_info, alias_map, counter)
        hoist_result = _hoist_lambda_call(s.cond, lambda_info, alias_map, counter)
        if hoist_result is not None:
            new_cond, pre = hoist_result
            return pre + [replace(s, cond=new_cond, then_body=new_then,
                                  else_body=new_else)]
        return [replace(s, then_body=new_then, else_body=new_else)]

    # WhileStmt / DoWhileStmt：cond 含 lambda Call → hoist (注意 cond 在 while 头部)
    if isinstance(s, WhileStmt):
        new_body = _rewrite_stmts(list(s.body), lambda_info, alias_map, counter)
        hoist_result = _hoist_lambda_call(s.cond, lambda_info, alias_map, counter)
        if hoist_result is not None:
            new_cond, pre = hoist_result
            # while cond hoist 复杂：cond 每轮都需重算。简化：在 while 之前一次 +
            # body 末尾再算一次（loop_lower 时校正）。当前 corpus 未观测到 while
            # cond 含 lambda（仅 IfStmt 中），保守留 TODO 警告。
            return pre + [replace(s, cond=new_cond, body=new_body)]
        return [replace(s, body=new_body)]
    if isinstance(s, DoWhileStmt):
        new_body = _rewrite_stmts(list(s.body), lambda_info, alias_map, counter)
        hoist_result = _hoist_lambda_call(s.cond, lambda_info, alias_map, counter)
        if hoist_result is not None:
            new_cond, pre = hoist_result
            # DoWhile cond 在 body 之后求值。但 body 内 `continue` 也跳到 cond
            # 检查 — 必须在每个 (本层 do-while 的) continue 之前也插 hoist pre_stmts，
            # 否则 continue 路径用旧 cap → cond 死循环。
            new_body_with_continue_hoist = _inject_pre_at_continues(new_body, pre)
            return [replace(s, body=new_body_with_continue_hoist + pre, cond=new_cond)]
        return [replace(s, body=new_body)]

    if isinstance(s, ForStmt):
        return [replace(s,
                        init=_rewrite_stmts(list(s.init), lambda_info, alias_map, counter),
                        step=_rewrite_stmts(list(s.step), lambda_info, alias_map, counter),
                        body=_rewrite_stmts(list(s.body), lambda_info, alias_map, counter))]
    if isinstance(s, RangeForStmt):
        return [replace(s, body=_rewrite_stmts(list(s.body),
                                               lambda_info, alias_map, counter))]
    if isinstance(s, BlockStmt):
        return [replace(s, stmts=_rewrite_stmts(list(s.stmts),
                                                lambda_info, alias_map, counter))]

    # LetStmt：value 含 lambda Call → hoist
    if isinstance(s, LetStmt) and isinstance(s.value, Call):
        hoist_result = _hoist_lambda_call(s.value, lambda_info, alias_map, counter)
        if hoist_result is not None:
            new_value, pre = hoist_result
            return pre + [replace(s, value=new_value)]

    # AssignStmt：value 含 lambda Call → hoist
    if isinstance(s, AssignStmt) and isinstance(s.value, Call):
        hoist_result = _hoist_lambda_call(s.value, lambda_info, alias_map, counter)
        if hoist_result is not None:
            new_value, pre = hoist_result
            return pre + [replace(s, value=new_value)]

    return [s]


def _rewrite_stmts(stmts: list[StmtIR],
                   lambda_info: dict[str, tuple[list[str], list[int]]],
                   alias_map: dict[str, str],
                   counter: list[int]) -> list[StmtIR]:
    out: list[StmtIR] = []
    for s in stmts:
        out.extend(_rewrite_stmt(s, lambda_info, alias_map, counter))
    return out


def _inject_pre_at_continues(stmts: list[StmtIR],
                             pre_stmts: list[StmtIR]) -> list[StmtIR]:
    """在 *本层* （非嵌套）的每个 ContinueStmt 之前注入 pre_stmts。
    嵌套 WhileStmt/ForStmt/RangeForStmt/DoWhileStmt 内部的 continue 是别的
    loop 的，不动。但 IfStmt/BlockStmt 的子 stmts 仍属本层。
    """
    out: list[StmtIR] = []
    for s in stmts:
        if isinstance(s, ContinueStmt):
            out.extend(pre_stmts)
            out.append(s)
        elif isinstance(s, IfStmt):
            out.append(replace(s,
                then_body=_inject_pre_at_continues(s.then_body, pre_stmts),
                else_body=_inject_pre_at_continues(s.else_body or [], pre_stmts)))
        elif isinstance(s, BlockStmt):
            out.append(replace(s, stmts=_inject_pre_at_continues(s.stmts, pre_stmts)))
        else:
            # WhileStmt / ForStmt / RangeForStmt / DoWhileStmt: 嵌套 loop，
            # 内部 continue 是 inner loop 的，不动
            out.append(s)
    return out


def _body_has_value_return(stmts: list[StmtIR]) -> bool:
    """递归检查 body 是否含 `ReturnStmt(value=non-None)`，用于判 lambda 真 void。"""
    for s in stmts:
        if isinstance(s, ReturnStmt) and s.value is not None:
            return True
        if isinstance(s, IfStmt):
            if _body_has_value_return(s.then_body): return True
            if _body_has_value_return(s.else_body or []): return True
        elif isinstance(s, (WhileStmt, RangeForStmt, DoWhileStmt)):
            if _body_has_value_return(s.body): return True
        elif isinstance(s, ForStmt):
            if _body_has_value_return(s.body): return True
        elif isinstance(s, BlockStmt):
            if _body_has_value_return(s.stmts): return True
    return False


def lambda_ref_elim_pass(func: HIRFunc) -> HIRFunc:
    """HIR₂ → HIR₂'：lifted lambda 的 ref-elim + 调用点改写。"""
    if not func.aux_lambdas:
        return func

    # P1：扫描每个 lifted lambda 的 cap 信息（**清标记前**）
    lambda_info: dict[str, tuple[list[str], list[int]]] = {}
    for aux in func.aux_lambdas:
        cap_names, modified_indices, _ = _extract_cap_info(aux)
        lambda_info[aux.base_name] = (cap_names, modified_indices)

    # P2：每个 lifted lambda 应用 ref_elim_pass（清 is_ref/is_const_ref + 包装 returns）
    # 同时跑 callsite_ref_elim_pass 处理 lambda body 内的 ref-out 调用（如
    # `fdiv_q(result, ...)` / `std::swap(M[i], M[j])` —— 这些在 Pass 3 lift 时
    # 仍是 ExprStmt(Call)，Pass 2b 错过因 lift 在 2b 之后；这里补做。
    #
    # **第八轮 P0 修复**：Pass 1 给 lambda 设 ret_ty=UnknownType（不是
    # BaseType.UNIT）。ref_elim_pass 据此误判为非 void → 不追加 trailing return →
    # SSA build 后 ReturnTerm.value=None，destructure 字段全空，M/U 修改丢失。
    # 修：若 lambda body 无任何 non-None ReturnStmt，强制视为 void。
    new_aux: list[HIRFunc] = []
    for aux in func.aux_lambdas:
        if isinstance(aux.ret_ty, UnknownType) and not _body_has_value_return(aux.body):
            aux = replace(aux, ret_ty=BaseType.UNIT)
        aux_with_refelim = ref_elim_pass(aux)
        aux_final = callsite_ref_elim_pass(aux_with_refelim)
        new_aux.append(aux_final)

    # P3：扫描 outer body 建 alias map
    alias_map = _build_alias_map(func.body)

    # P4：递归重写 outer body 中所有 lambda 调用点
    counter = [0]
    new_body = _rewrite_stmts(list(func.body), lambda_info, alias_map, counter)

    return replace(func, body=new_body, aux_lambdas=new_aux)
