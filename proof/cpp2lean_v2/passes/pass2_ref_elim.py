"""
Pass 2: ref_elim — HIR₀ → HIR₁

职责：消除 `T&` / `T&&` 参数，转换为 tuple 返回值。
依据 hir-design.md §4 + cpp-subset-semantics.md §4.3-4.4（引理 L4.1）。

HIR₁ 不变量：func.params 所有 is_ref == False（const ref 保留）。

## 算法

1. 扫参数表，收集 ref 参数列表。若无 ref，仅清理标志返回。
2. 否则：
   - 每个 ref 参数在参数列表中保留（作为"初值输入"），但标记 is_ref=False
   - 返回类型改为 tuple：
     * void 原返回 + N 个 ref → 单 ref 时直接用 ref 类型；多 ref 时 PairType/TupleType
     * 非 void 原返回 + N 个 ref → PairType/TupleType[原返回, ref1, ..., refN]
   - body 内每个 ReturnStmt 的值包装为 tuple
   - void 函数若 body 末尾非 ReturnStmt，追加 ReturnStmt(tuple)
"""

from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, ArrayType, PairType, TupleType, TypeIR,
    Var, TupleExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt,
    IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
    BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRFunc, HIRParam, TranslationError,
)


def ref_elim_pass(func: HIRFunc) -> HIRFunc:
    """HIR₀ → HIR₁。消除 non-const ref 参数，转为 tuple 返回值。"""

    ref_params = [p for p in func.params if p.is_ref]

    if not ref_params:
        # 无 ref 参数：仅清理 is_output 标志（保持 is_ref=False 已是）
        new_params = [
            HIRParam(
                name=p.name,
                ty=p.ty,
                is_ref=False,
                is_const_ref=p.is_const_ref,
                is_output=False,
            )
            for p in func.params
        ]
        return HIRFunc(
            base_name=func.base_name,
            instance_suffix=func.instance_suffix,
            mangled_name=func.mangled_name,
            qual_type=func.qual_type,
            params=new_params,
            ret_ty=func.ret_ty,
            body=func.body,
            requires=func.requires,
            aux_lambdas=func.aux_lambdas,
        )

    # 1. 改写参数列表：ref 参数保留为"初值参数"，但清除 is_ref 标志
    new_params: list[HIRParam] = []
    for p in func.params:
        if p.is_ref:
            new_params.append(HIRParam(
                name=p.name,
                ty=p.ty,
                is_ref=False,
                is_const_ref=False,
                is_output=False,  # Pass 1 的 is_output 在 HIR₁ 不再使用
            ))
        else:
            new_params.append(HIRParam(
                name=p.name,
                ty=p.ty,
                is_ref=False,
                is_const_ref=p.is_const_ref,
                is_output=False,
            ))

    # 2. 计算新返回类型：
    #    void + N ref → N==1: ref 类型；N>1: TupleType
    #    非 void + N ref → TupleType / PairType(orig, ref1, ..., refN)
    ref_types = [p.ty for p in ref_params]
    original_is_void = func.ret_ty == BaseType.UNIT

    if original_is_void:
        new_ret_ty = _tuple_of(ref_types)
    else:
        new_ret_ty = _tuple_of([func.ret_ty] + ref_types)

    # 3. 改写 body：每个 ReturnStmt 的值包装为 tuple
    ref_names = [p.name for p in ref_params]
    new_body = _rewrite_returns(func.body, ref_names, original_is_void)

    # 4. void 函数若末尾无 ReturnStmt，追加一条
    if original_is_void and not _body_always_returns(new_body):
        trailing_return = ReturnStmt(
            value=_make_tuple_expr([Var(n) for n in ref_names])
        )
        new_body = new_body + [trailing_return]

    return HIRFunc(
        base_name=func.base_name,
        instance_suffix=func.instance_suffix,
        mangled_name=func.mangled_name,
        qual_type=func.qual_type,
        params=new_params,
        ret_ty=new_ret_ty,
        body=new_body,
        requires=func.requires,
        aux_lambdas=func.aux_lambdas,
    )


def _tuple_of(tys: list[TypeIR]) -> TypeIR:
    """把类型列表聚合为 tuple（1 个直接返回；2 个用 PairType；3+ 用 TupleType）。"""
    if len(tys) == 1:
        return tys[0]
    if len(tys) == 2:
        return PairType(tys[0], tys[1])
    return TupleType(tuple(tys))


def _make_tuple_expr(elems: list[ExprIR]) -> ExprIR:
    """把 ExprIR 列表聚合为 tuple 表达式。"""
    if len(elems) == 1:
        return elems[0]
    return TupleExpr(elems=elems)


def _rewrite_returns(stmts: list[StmtIR], ref_names: list[str],
                     original_is_void: bool) -> list[StmtIR]:
    """递归改写所有 ReturnStmt。
    - void 函数：return; → return (ref_tuple)
    - 非 void 函数：return v; → return (v, ref_tuple...)
    """
    result = []
    for s in stmts:
        result.append(_rewrite_stmt_returns(s, ref_names, original_is_void))
    return result


def _rewrite_stmt_returns(s: StmtIR, ref_names: list[str],
                          original_is_void: bool) -> StmtIR:
    """递归改写单条语句内的 ReturnStmt。"""
    if isinstance(s, ReturnStmt):
        ref_vars = [Var(n) for n in ref_names]
        if original_is_void:
            # void return; → return (ref_tuple)
            new_val = _make_tuple_expr(ref_vars)
        else:
            # return v; → return (v, ref1, ..., refN)
            orig = s.value
            if orig is None:
                # 理论上非 void 不应裸 return，容错处理
                new_val = _make_tuple_expr(ref_vars)
            else:
                new_val = _make_tuple_expr([orig] + ref_vars)
        return ReturnStmt(value=new_val)

    if isinstance(s, IfStmt):
        return IfStmt(
            cond=s.cond,
            then_body=_rewrite_returns(s.then_body, ref_names, original_is_void),
            else_body=_rewrite_returns(s.else_body, ref_names, original_is_void),
        )

    if isinstance(s, WhileStmt):
        return WhileStmt(
            cond=s.cond,
            body=_rewrite_returns(s.body, ref_names, original_is_void),
        )

    if isinstance(s, ForStmt):
        return ForStmt(
            init=_rewrite_returns(s.init, ref_names, original_is_void),
            cond=s.cond,
            step=_rewrite_returns(s.step, ref_names, original_is_void),
            body=_rewrite_returns(s.body, ref_names, original_is_void),
        )

    if isinstance(s, RangeForStmt):
        return RangeForStmt(
            var=s.var,
            var_ty=s.var_ty,
            container=s.container,
            body=_rewrite_returns(s.body, ref_names, original_is_void),
            decomposition=s.decomposition,
        )

    if isinstance(s, DoWhileStmt):
        return DoWhileStmt(
            body=_rewrite_returns(s.body, ref_names, original_is_void),
            cond=s.cond,
        )

    if isinstance(s, BlockStmt):
        return BlockStmt(
            stmts=_rewrite_returns(s.stmts, ref_names, original_is_void),
        )

    # 其他语句不含 ReturnStmt，原样返回
    return s


def _body_always_returns(body: list[StmtIR]) -> bool:
    """判断 body 的末尾路径是否总会 return（用于决定是否追加 trailing return）。

    保守策略：只在**显式末尾 ReturnStmt** 时认为 "always returns"。
    if/while/for 内部的 return 不算（因为可能走 else 或 0 次循环）。
    """
    if not body:
        return False
    last = body[-1]
    if isinstance(last, ReturnStmt):
        return True
    if isinstance(last, IfStmt):
        # 仅 then+else 都以 return 结尾时才保证 return
        return (_body_always_returns(last.then_body)
                and _body_always_returns(last.else_body))
    if isinstance(last, BlockStmt):
        return _body_always_returns(last.stmts)
    return False


# ============================================================
# 不变量检查
# ============================================================

def assert_hir1_invariant(func: HIRFunc):
    """HIR₁ 出口 assert：
    - 所有 params.is_ref == False
    - 所有 params.is_output == False
    """
    for i, p in enumerate(func.params):
        if p.is_ref:
            raise TranslationError(
                pass_name="ref_elim",
                func_name=func.base_name,
                reason=f"param[{i}] {p.name} still has is_ref=True after ref_elim",
            )
        if p.is_output:
            raise TranslationError(
                pass_name="ref_elim",
                func_name=func.base_name,
                reason=f"param[{i}] {p.name} still has is_output=True after ref_elim",
            )
