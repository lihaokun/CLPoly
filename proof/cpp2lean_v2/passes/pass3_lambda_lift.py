"""
Pass 3: lambda_lift — HIR₁ → HIR₂

职责：把所有内联 LambdaExpr 提升为独立 HIRFunc（存入 aux_lambdas）；
调用点的 LambdaExpr 节点替换为对新函数的 Var 引用。

依据 hir-design.md §5 + cpp-subset-semantics.md §5（引理 L5.1）。

HIR₂ 不变量：func body + aux_lambdas body 内无 LambdaExpr 节点。

## 算法

1. 遍历 func.body 找 LambdaExpr
2. 对每个 LambdaExpr：
   - 推导 captures：自由变量分析（body 中的 Var.name ∉ params ∪ local）
   - 检测 modified_captures（body 中 AssignStmt/CompoundAssignStmt 对 capture 的修改）
   - 构造新 HIRFunc `_lambda_<host>_<N>_ir`
   - 参数 = captures + 原 lambda params
   - body 原样复制
   - 替换原 LambdaExpr 为 Var("_lambda_<host>_<N>")（Pass 5 operator_resolve 处理调用点）
3. 子 lambda 递归处理（但本项目无 nested lambda）

## 简化说明

- Pass 3 不处理 modified capture 的 call-site 写回（留给 Pass 5 或专用 sub-pass）
- 但在 aux HIRFunc 上标记 `modified_captures` 供后续使用
- 对 CLPoly 现状：26 个 lambda，其中约 10 个 `[&]` 默认捕获（保守全视作可能修改）
"""

from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRFunc, HIRParam, TranslationError,
)


# ============================================================
# 自由变量与 local 变量分析
# ============================================================

def _collect_vars_expr(e: ExprIR, result: set[str]):
    """收集表达式中引用的变量名（Var.name），加入 result。"""
    if isinstance(e, Var):
        result.add(e.name)
    elif isinstance(e, BinOp):
        _collect_vars_expr(e.lhs, result); _collect_vars_expr(e.rhs, result)
    elif isinstance(e, UnaryOp):
        _collect_vars_expr(e.operand, result)
    elif isinstance(e, CondExpr):
        _collect_vars_expr(e.cond, result)
        _collect_vars_expr(e.then_e, result)
        _collect_vars_expr(e.else_e, result)
    elif isinstance(e, Call):
        if isinstance(e.callee, UnresolvedOp):
            pass  # UnresolvedOp 不是 Var
        # args
        for a in e.args:
            _collect_vars_expr(a, result)
    elif isinstance(e, ArrayAccess):
        _collect_vars_expr(e.arr, result); _collect_vars_expr(e.idx, result)
    elif isinstance(e, FieldAccess):
        _collect_vars_expr(e.obj, result)
    elif isinstance(e, Cast):
        _collect_vars_expr(e.expr, result)
    elif isinstance(e, BlockExpr):
        _collect_vars_stmts(e.stmts, result, locals_set=None)
        _collect_vars_expr(e.value, result)
    elif isinstance(e, TupleExpr):
        for x in e.elems: _collect_vars_expr(x, result)
    elif isinstance(e, ArrayLit):
        for x in e.elems: _collect_vars_expr(x, result)
    elif isinstance(e, LambdaExpr):
        # 嵌套 lambda：它的自由变量也算父 lambda 的使用
        # （CLPoly 目前没有嵌套 lambda，保守收集）
        for p in e.params:
            pass  # lambda 的 param 是内层局部
        # 收集其 body 的自由变量（不包括其自身 params / locals）
        nested_vars = _collect_free_vars(e.body, set(p.name for p in e.params))
        result |= nested_vars
    # Lit, UnresolvedOp, UnknownExpr: 无变量引用


def _collect_vars_stmt(s: StmtIR, used: set[str], locals_set: set[str] | None):
    """遍历语句，收集 used vars；同时扩充 locals_set（若提供）。"""
    if isinstance(s, LetStmt):
        _collect_vars_expr(s.value, used)
        if locals_set is not None:
            locals_set.add(s.var.name)
    elif isinstance(s, AssignStmt):
        _collect_vars_expr(s.target, used)
        _collect_vars_expr(s.value, used)
    elif isinstance(s, CompoundAssignStmt):
        _collect_vars_expr(s.target, used)
        _collect_vars_expr(s.value, used)
    elif isinstance(s, IfStmt):
        _collect_vars_expr(s.cond, used)
        _collect_vars_stmts(s.then_body, used, locals_set)
        _collect_vars_stmts(s.else_body, used, locals_set)
    elif isinstance(s, WhileStmt):
        _collect_vars_expr(s.cond, used)
        _collect_vars_stmts(s.body, used, locals_set)
    elif isinstance(s, ForStmt):
        _collect_vars_stmts(s.init, used, locals_set)
        _collect_vars_expr(s.cond, used)
        _collect_vars_stmts(s.step, used, locals_set)
        _collect_vars_stmts(s.body, used, locals_set)
    elif isinstance(s, RangeForStmt):
        _collect_vars_expr(s.container, used)
        if locals_set is not None:
            locals_set.add(s.var.name)
            if s.decomposition:
                for v in s.decomposition:
                    locals_set.add(v.name)
        _collect_vars_stmts(s.body, used, locals_set)
    elif isinstance(s, DoWhileStmt):
        _collect_vars_stmts(s.body, used, locals_set)
        _collect_vars_expr(s.cond, used)
    elif isinstance(s, BlockStmt):
        _collect_vars_stmts(s.stmts, used, locals_set)
    elif isinstance(s, ReturnStmt):
        if s.value is not None:
            _collect_vars_expr(s.value, used)
    elif isinstance(s, RequireStmt):
        _collect_vars_expr(s.cond, used)
    elif isinstance(s, ExprStmt):
        _collect_vars_expr(s.expr, used)
    # BreakStmt / ContinueStmt / UnknownStmt: 无


def _collect_vars_stmts(stmts: list[StmtIR], used: set[str],
                        locals_set: set[str] | None):
    for s in stmts:
        _collect_vars_stmt(s, used, locals_set)


def _collect_free_vars(body: list[StmtIR], params: set[str]) -> set[str]:
    """收集 body 里的自由变量 = used - params - locals。"""
    used: set[str] = set()
    locals_set: set[str] = set()
    _collect_vars_stmts(body, used, locals_set)
    return used - params - locals_set


# ============================================================
# 修改 capture 检测（AssignStmt / CompoundAssignStmt on capture var）
# ============================================================

def _collect_modified(stmts: list[StmtIR]) -> set[str]:
    """返回 body 中通过赋值被修改的顶层变量名集合（Var 目标）。
    对 ArrayAccess/FieldAccess 递归找根 Var。"""
    modified: set[str] = set()

    def target_root(expr: ExprIR) -> str | None:
        if isinstance(expr, Var):
            return expr.name
        if isinstance(expr, ArrayAccess):
            return target_root(expr.arr)
        if isinstance(expr, FieldAccess):
            return target_root(expr.obj)
        if isinstance(expr, Cast):
            return target_root(expr.expr)
        # Pass 1 对 a[i] 可能生成 Call(UnresolvedOp("operator[]"), [arr, idx])
        # 而非 ArrayAccess，需要递归到 arr 上。
        if (isinstance(expr, Call)
                and isinstance(expr.callee, UnresolvedOp)
                and expr.callee.op_name == "operator[]"
                and expr.args):
            return target_root(expr.args[0])
        return None

    def walk(s: StmtIR):
        if isinstance(s, AssignStmt):
            n = target_root(s.target)
            if n: modified.add(n)
        elif isinstance(s, CompoundAssignStmt):
            n = target_root(s.target)
            if n: modified.add(n)
        elif isinstance(s, IfStmt):
            for t in s.then_body: walk(t)
            for t in s.else_body: walk(t)
        elif isinstance(s, WhileStmt):
            for t in s.body: walk(t)
        elif isinstance(s, ForStmt):
            for t in s.init: walk(t)
            for t in s.step: walk(t)
            for t in s.body: walk(t)
        elif isinstance(s, RangeForStmt):
            for t in s.body: walk(t)
        elif isinstance(s, DoWhileStmt):
            for t in s.body: walk(t)
        elif isinstance(s, BlockStmt):
            for t in s.stmts: walk(t)
        # UnresolvedOp("operator++/--/+=/-=/*=/ 等") 作为 ExprStmt 也是修改；
        # `swap(a, b)` / `std::swap(a, b)` 修改两个参数。
        elif isinstance(s, ExprStmt):
            e = s.expr
            if isinstance(e, Call):
                callee = e.callee
                if isinstance(callee, UnresolvedOp):
                    op = callee.op_name
                    if any(m in op for m in ("operator++", "operator--", "operator+=",
                                              "operator-=", "operator*=", "operator/=",
                                              "operator%=", "operator<<=", "operator>>=",
                                              "operator|=", "operator&=", "operator^=")):
                        if e.args:
                            n = target_root(e.args[0])
                            if n: modified.add(n)
                elif isinstance(callee, str) and callee in ("swap", "std::swap"):
                    for a in e.args[:2]:
                        n = target_root(a)
                        if n: modified.add(n)

    for s in stmts:
        walk(s)
    return modified


# ============================================================
# Lambda 提升
# ============================================================

def _lift_lambda(lam: LambdaExpr, host_name: str, counter: list[int],
                 outer_params: list[HIRParam],
                 outer_typectx: dict[str, TypeIR]) -> tuple[Var, HIRFunc]:
    """把 LambdaExpr 提升为独立 HIRFunc；返回 (Var 替换节点, 新 HIRFunc)。

    outer_params: 宿主 HIRFunc 的参数（决定 capture 的类型）
    outer_typectx: {var_name → type}，从宿主 body 的 LetStmt 推导
    """
    counter[0] += 1
    lam_name = f"_lambda_{host_name}_{counter[0]}"

    # 1. 分析自由变量 = captures
    param_names = {p.name for p in lam.params}
    free_vars = _collect_free_vars(lam.body, param_names)

    # 2. 过滤掉 free_vars 中不是外层可见的（可能是未在 outer_typectx 的东西，
    #    如已知 class/function 名 —— 保守只取外层已定义的）
    capture_names = sorted(free_vars & set(outer_typectx.keys()))

    # 3. 检测 modified captures
    modified = _collect_modified(lam.body)
    modified_captures = sorted(set(capture_names) & modified)

    # 4. 构造新 HIRFunc 的 params = captures + 原 lambda params
    #    被修改的 capture 用可变 [REF]；只读的用 [CONST-REF]。
    cap_params = [
        HIRParam(
            name=cn,
            ty=outer_typectx.get(cn, UnknownType("")),
            is_ref=(cn in modified),
            is_const_ref=(cn not in modified),
            is_output=False,
        )
        for cn in capture_names
    ]
    new_params = cap_params + list(lam.params)

    # 5. 返回类型：若有 modified capture → (原返回, modified_cap...) tuple
    #    简化：Pass 3 不改返回类型（Pass 5 或后续 pass 可加）。
    # Lambda 原始返回类型从 ty 取；若无则用 UnknownType
    orig_ret = lam.ty if lam.ty else UnknownType("")
    # 若 Lambda 的 ty 是 NamedType("Lambda") 或其他，取 UnknownType("")
    if isinstance(orig_ret, NamedType) and orig_ret.name == "Lambda":
        orig_ret = UnknownType("")

    # 6. 生成新 HIRFunc
    # qual_type 编码 cap 信息（Pass 3b 用于区分 cap_params vs lambda by-ref params，
    # 因为 cap_params 与 by-ref lambda 参数都用 is_ref=True 标记，无法仅靠
    # HIRParam 字段区分）：
    #   - n_caps：cap_params 的个数（即前 N 个 params 是 captures）
    #   - modified_captures：被修改的 capture 名字列表
    qual_type = (f"lambda in {host_name} | n_caps={len(cap_params)} "
                 f"| modified_captures={modified_captures}")
    new_func = HIRFunc(
        base_name=lam_name,
        instance_suffix="",
        mangled_name="",
        qual_type=qual_type,
        params=new_params,
        ret_ty=orig_ret,
        body=lam.body,
        requires=[],
        aux_lambdas=[],
    )

    # 7. 调用点替换：Var 指向新函数
    replacement = Var(name=lam_name, version=0, ty=NamedType("LambdaRef"))
    return replacement, new_func


# ============================================================
# 遍历 + 替换
# ============================================================

def _walk_and_lift(stmts: list[StmtIR], host_name: str, counter: list[int],
                   outer_params: list[HIRParam],
                   typectx: dict[str, TypeIR],
                   aux_collect: list[HIRFunc]) -> list[StmtIR]:
    """递归遍历 stmts，把所有 LambdaExpr 提升，替换为 Var。"""
    result: list[StmtIR] = []
    for s in stmts:
        result.append(_rewrite_stmt(s, host_name, counter, outer_params,
                                    typectx, aux_collect))
    return result


def _rewrite_expr(e: ExprIR, host_name: str, counter: list[int],
                  outer_params: list[HIRParam],
                  typectx: dict[str, TypeIR],
                  aux_collect: list[HIRFunc]) -> ExprIR:
    if isinstance(e, LambdaExpr):
        replacement, new_func = _lift_lambda(e, host_name, counter,
                                              outer_params, typectx)
        aux_collect.append(new_func)
        return replacement
    # 递归子表达式
    if isinstance(e, BinOp):
        return replace(e,
            lhs=_rewrite_expr(e.lhs, host_name, counter, outer_params, typectx, aux_collect),
            rhs=_rewrite_expr(e.rhs, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, UnaryOp):
        return replace(e,
            operand=_rewrite_expr(e.operand, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, CondExpr):
        return replace(e,
            cond=_rewrite_expr(e.cond, host_name, counter, outer_params, typectx, aux_collect),
            then_e=_rewrite_expr(e.then_e, host_name, counter, outer_params, typectx, aux_collect),
            else_e=_rewrite_expr(e.else_e, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, Call):
        new_args = [_rewrite_expr(a, host_name, counter, outer_params, typectx, aux_collect)
                    for a in e.args]
        return replace(e, args=new_args)
    if isinstance(e, ArrayAccess):
        return replace(e,
            arr=_rewrite_expr(e.arr, host_name, counter, outer_params, typectx, aux_collect),
            idx=_rewrite_expr(e.idx, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, FieldAccess):
        return replace(e,
            obj=_rewrite_expr(e.obj, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, Cast):
        return replace(e,
            expr=_rewrite_expr(e.expr, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, BlockExpr):
        return replace(e,
            stmts=_walk_and_lift(e.stmts, host_name, counter, outer_params, typectx, aux_collect),
            value=_rewrite_expr(e.value, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(e, TupleExpr):
        return replace(e, elems=[_rewrite_expr(x, host_name, counter, outer_params, typectx, aux_collect)
                                 for x in e.elems])
    if isinstance(e, ArrayLit):
        return replace(e, elems=[_rewrite_expr(x, host_name, counter, outer_params, typectx, aux_collect)
                                 for x in e.elems])
    # Var / Lit / UnresolvedOp / UnknownExpr 原样
    return e


def _rewrite_stmt(s: StmtIR, host_name: str, counter: list[int],
                  outer_params: list[HIRParam],
                  typectx: dict[str, TypeIR],
                  aux_collect: list[HIRFunc]) -> StmtIR:
    if isinstance(s, LetStmt):
        new_val = _rewrite_expr(s.value, host_name, counter, outer_params, typectx, aux_collect)
        # 把这个 let 的 var 加入 typectx（后续语句里是 local）
        typectx[s.var.name] = s.ty
        return replace(s, value=new_val)
    if isinstance(s, AssignStmt):
        return replace(s,
            target=_rewrite_expr(s.target, host_name, counter, outer_params, typectx, aux_collect),
            value=_rewrite_expr(s.value, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, CompoundAssignStmt):
        return replace(s,
            target=_rewrite_expr(s.target, host_name, counter, outer_params, typectx, aux_collect),
            value=_rewrite_expr(s.value, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, IfStmt):
        return replace(s,
            cond=_rewrite_expr(s.cond, host_name, counter, outer_params, typectx, aux_collect),
            then_body=_walk_and_lift(s.then_body, host_name, counter, outer_params, typectx, aux_collect),
            else_body=_walk_and_lift(s.else_body, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, WhileStmt):
        return replace(s,
            cond=_rewrite_expr(s.cond, host_name, counter, outer_params, typectx, aux_collect),
            body=_walk_and_lift(s.body, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, ForStmt):
        return replace(s,
            init=_walk_and_lift(s.init, host_name, counter, outer_params, typectx, aux_collect),
            cond=_rewrite_expr(s.cond, host_name, counter, outer_params, typectx, aux_collect),
            step=_walk_and_lift(s.step, host_name, counter, outer_params, typectx, aux_collect),
            body=_walk_and_lift(s.body, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, RangeForStmt):
        return replace(s,
            container=_rewrite_expr(s.container, host_name, counter, outer_params, typectx, aux_collect),
            body=_walk_and_lift(s.body, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, DoWhileStmt):
        return replace(s,
            body=_walk_and_lift(s.body, host_name, counter, outer_params, typectx, aux_collect),
            cond=_rewrite_expr(s.cond, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, BlockStmt):
        return replace(s,
            stmts=_walk_and_lift(s.stmts, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, ReturnStmt):
        if s.value is not None:
            return replace(s,
                value=_rewrite_expr(s.value, host_name, counter, outer_params, typectx, aux_collect),
            )
        return s
    if isinstance(s, RequireStmt):
        return replace(s,
            cond=_rewrite_expr(s.cond, host_name, counter, outer_params, typectx, aux_collect),
        )
    if isinstance(s, ExprStmt):
        return replace(s,
            expr=_rewrite_expr(s.expr, host_name, counter, outer_params, typectx, aux_collect),
        )
    # BreakStmt / ContinueStmt / UnknownStmt 原样
    return s


# ============================================================
# Pass 3 入口
# ============================================================

def lambda_lift_pass(func: HIRFunc) -> HIRFunc:
    """HIR₁ → HIR₂。提升所有 LambdaExpr 为 aux_lambdas。"""
    counter = [0]  # mutable box，方便嵌套计数
    aux_collect: list[HIRFunc] = []

    # 初始 typectx = 参数
    typectx: dict[str, TypeIR] = {p.name: p.ty for p in func.params}

    new_body = _walk_and_lift(
        func.body, func.base_name, counter,
        outer_params=func.params,
        typectx=typectx,
        aux_collect=aux_collect,
    )

    return HIRFunc(
        base_name=func.base_name,
        instance_suffix=func.instance_suffix,
        mangled_name=func.mangled_name,
        qual_type=func.qual_type,
        params=func.params,
        ret_ty=func.ret_ty,
        body=new_body,
        requires=func.requires,
        aux_lambdas=list(func.aux_lambdas) + aux_collect,
    )


# ============================================================
# 不变量检查
# ============================================================

def assert_hir2_invariant(func: HIRFunc):
    """HIR₂ 出口 assert：func.body 及 aux_lambdas[*].body 内均无 LambdaExpr。"""
    def _check_no_lambda_expr(e: ExprIR):
        if isinstance(e, LambdaExpr):
            raise TranslationError(
                pass_name="lambda_lift",
                func_name=func.base_name,
                reason=f"inline LambdaExpr still present in body",
            )
        # 递归
        for attr in dir(e):
            if attr.startswith("_"): continue
            try:
                v = getattr(e, attr)
            except Exception:
                continue
            if isinstance(v, list):
                for x in v:
                    if isinstance(x, tuple(_EXPR_TYPES)):
                        _check_no_lambda_expr(x)
                    elif isinstance(x, tuple(_STMT_TYPES)):
                        _check_no_lambda_stmt(x)

    def _check_no_lambda_stmt(s):
        for attr in dir(s):
            if attr.startswith("_"): continue
            try:
                v = getattr(s, attr)
            except Exception:
                continue
            if isinstance(v, list):
                for x in v:
                    if isinstance(x, tuple(_EXPR_TYPES)):
                        _check_no_lambda_expr(x)
                    elif isinstance(x, tuple(_STMT_TYPES)):
                        _check_no_lambda_stmt(x)
            elif isinstance(v, tuple(_EXPR_TYPES)):
                _check_no_lambda_expr(v)

    for s in func.body:
        _check_no_lambda_stmt(s)
    for aux in func.aux_lambdas:
        for s in aux.body:
            _check_no_lambda_stmt(s)


_EXPR_TYPES = (Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
               ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
               BlockExpr, TupleExpr, ArrayLit, UnknownExpr)

_STMT_TYPES = (LetStmt, AssignStmt, CompoundAssignStmt,
               IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
               BreakStmt, ContinueStmt, ReturnStmt,
               RequireStmt, ExprStmt, BlockStmt, UnknownStmt)
