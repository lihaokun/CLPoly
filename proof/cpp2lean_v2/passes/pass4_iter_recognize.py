"""
Pass 4: iter_recognize — HIR₂ → HIR₃

职责：
  1. Structured-binding desugar（23 处）：`RangeForStmt.decomposition` 非空时，
     把 `[k, v]` 绑定展开为 body 头部的 LetStmt(k, _x.fst), LetStmt(v, _x.snd)；
     decomposition 字段清空。
  2. Filter-loop 识别（6 处）：两种源形态都转为 `AssignStmt(v, Array.filter(v, pred))`：
       - 形态 A compact-erase（4 处）：4 件套
         [s0 赋 begin; s1 赋 alias of it; s2 For/WhileStmt; s3 erase ExprStmt]
       - 形态 B classic iter-loop（2 处）：For/WhileStmt 单体，body 内深嵌 erase
  3. 并行双迭代器（1 处 zip-walk，`__upoly_divmod_mod`）不处理，原样透传。

依据 hir-design.md §6（已更新 2026-04-23 实证）。

HIR₃ 不变量：
  - `RangeForStmt.decomposition == None` 对所有 RangeForStmt
  - 不存在未识别的 filter-loop 4 件套（A 形态）序列
  - 允许残留 `Call(UnresolvedOp("<method>.xxx"))`；Pass 5 作为通用 Call 处理
"""

from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ir_types import (
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType, RefType, ArrayType,
    StdMapType,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr, ExprIR,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt, StmtIR,
    HIRFunc, HIRParam, TranslationError,
)


# ============================================================
# Helper：callee / expr 判别
# ============================================================

def _callee_str(e: ExprIR) -> str | None:
    """Call 的 callee 归一为字符串（str 或 UnresolvedOp.op_name）。"""
    if isinstance(e, Call):
        c = e.callee
        if isinstance(c, str): return c
        if isinstance(c, UnresolvedOp): return c.op_name
    return None


def _strip_cast(e: ExprIR) -> ExprIR:
    """剥 Cast 外壳（NoOp / LValueToRValue / ConstructorConversion 等都剥）。"""
    while isinstance(e, Cast):
        e = e.expr
    return e


def _expr_has_call_matching(e: ExprIR, substr: str) -> bool:
    """深度搜索：e 内是否含 callee 含 substr 的 Call。"""
    if isinstance(e, Call):
        cs = _callee_str(e)
        if cs and substr in cs: return True
        for a in e.args:
            if _expr_has_call_matching(a, substr): return True
    if isinstance(e, Cast):
        return _expr_has_call_matching(e.expr, substr)
    if isinstance(e, BinOp):
        return _expr_has_call_matching(e.lhs, substr) or _expr_has_call_matching(e.rhs, substr)
    if isinstance(e, UnaryOp):
        return _expr_has_call_matching(e.operand, substr)
    if isinstance(e, CondExpr):
        return (_expr_has_call_matching(e.cond, substr)
                or _expr_has_call_matching(e.then_e, substr)
                or _expr_has_call_matching(e.else_e, substr))
    if isinstance(e, (ArrayAccess,)):
        return _expr_has_call_matching(e.arr, substr) or _expr_has_call_matching(e.idx, substr)
    if isinstance(e, FieldAccess):
        return _expr_has_call_matching(e.obj, substr)
    if isinstance(e, (TupleExpr, ArrayLit)):
        return any(_expr_has_call_matching(x, substr) for x in e.elems)
    return False


_ITER_METHODS = (".begin", ".end", ".rbegin", ".rend", ".cbegin", ".cend")


def _expr_has_iter_method(e: ExprIR) -> bool:
    return any(_expr_has_call_matching(e, m) for m in _ITER_METHODS)


def _expr_has_erase(e: ExprIR) -> bool:
    return _expr_has_call_matching(e, ".erase")


def _stmt_rhs(s: StmtIR) -> ExprIR | None:
    """抽 stmt 的 'value' 端（赋值 stmt 的右侧），否则 None。
    覆盖 LetStmt / AssignStmt / ExprStmt(Call("operator=", [lhs, rhs]))。"""
    if isinstance(s, LetStmt): return s.value
    if isinstance(s, AssignStmt): return s.value
    if isinstance(s, ExprStmt) and isinstance(s.expr, Call):
        if _callee_str(s.expr) == "operator=" and len(s.expr.args) == 2:
            return s.expr.args[1]
    return None


def _is_erase_stmt(s: StmtIR) -> bool:
    """top-level erase 调用的 ExprStmt（A 形态的尾部）。"""
    if isinstance(s, ExprStmt) and isinstance(s.expr, Call):
        cs = _callee_str(s.expr)
        return bool(cs and ".erase" in cs)
    return False


def _deep_has_erase(stmts: list[StmtIR]) -> bool:
    """body 内（递归）是否含 erase 调用——可嵌在 operator= / 普通 Call.args 里。"""
    for s in stmts:
        if isinstance(s, ExprStmt) and _expr_has_erase(s.expr): return True
        if isinstance(s, AssignStmt) and _expr_has_erase(s.value): return True
        if isinstance(s, IfStmt):
            if _deep_has_erase(s.then_body) or _deep_has_erase(s.else_body): return True
        elif isinstance(s, ForStmt):
            if _deep_has_erase(s.init) or _deep_has_erase(s.step) or _deep_has_erase(s.body): return True
        elif isinstance(s, (WhileStmt, DoWhileStmt, RangeForStmt)):
            if _deep_has_erase(s.body): return True
        elif isinstance(s, BlockStmt):
            if _deep_has_erase(s.stmts): return True
    return False


# ============================================================
# Filter-loop 识别
# ============================================================

class FilterLoopMatch:
    """识别到的 filter-loop 描述。

    kind:
      "A" / "B-For" / "B-While": pure filter → Array.filter
      "A-mut": mutate-then-filter → Array.filterMap（CF-1 阶段 1）
    """
    __slots__ = ("start_idx", "end_idx", "container", "pred_lambda", "kind")

    def __init__(self, start_idx, end_idx, container, pred_lambda, kind):
        self.start_idx = start_idx
        self.end_idx = end_idx  # exclusive
        self.container = container
        self.pred_lambda = pred_lambda
        self.kind = kind


def _extract_container_from_call(e: ExprIR) -> ExprIR | None:
    """从 `<method>.begin(X)` / `<method>.end(X)` 抽出容器 X。"""
    stripped = _strip_cast(e)
    if isinstance(stripped, Call):
        cs = _callee_str(stripped)
        if cs and any(m in cs for m in _ITER_METHODS) and stripped.args:
            return stripped.args[0]
    return None


def _build_pred_lambda(pred_cond: ExprIR, loop_body_it_name: str,
                       element_ty: TypeIR | None) -> LambdaExpr:
    """从 IfStmt.cond（在 it 上 deref）构造 lambda(x: elem_ty) -> Bool。

    把 pred_cond 里所有对 `*it` / `it->field` / `Call("operator*", [Var(it)])`
    的引用替换为 `Var("x")` / `FieldAccess(Var("x"), ...)`。
    """
    x_name = "__x"
    elem_ty = element_ty if element_ty else UnknownType("")

    def rewrite(e: ExprIR) -> ExprIR:
        if isinstance(e, Var):
            return e
        if isinstance(e, Call):
            cs = _callee_str(e)
            # `*it` → x
            if cs == "operator*" and len(e.args) == 1:
                inner = _strip_cast(e.args[0])
                if isinstance(inner, Var) and inner.name == loop_body_it_name:
                    return Var(name=x_name, version=0, ty=elem_ty)
            # `it->field` 常表示为 FieldAccess(Call("operator->", [Var(it)]), field) 或
            # 直接 `operator->(it).field` — 见 HIR₂ dump `<operator->>(it).second`
            # 这种 FieldAccess(Call("operator->", ...), ...) 会走到 FieldAccess 分支
            return replace(e, args=[rewrite(a) for a in e.args])
        if isinstance(e, FieldAccess):
            obj = _strip_cast(e.obj)
            # it->field：obj = Call(UnresolvedOp("operator->"), [Var(it)])
            if isinstance(obj, Call) and _callee_str(obj) == "operator->":
                inner = _strip_cast(obj.args[0]) if obj.args else None
                if isinstance(inner, Var) and inner.name == loop_body_it_name:
                    return FieldAccess(obj=Var(name=x_name, version=0, ty=elem_ty),
                                       field_name=e.field_name, ty=e.ty)
            return replace(e, obj=rewrite(e.obj))
        if isinstance(e, BinOp):
            return replace(e, lhs=rewrite(e.lhs), rhs=rewrite(e.rhs))
        if isinstance(e, UnaryOp):
            return replace(e, operand=rewrite(e.operand))
        if isinstance(e, CondExpr):
            return replace(e, cond=rewrite(e.cond),
                           then_e=rewrite(e.then_e), else_e=rewrite(e.else_e))
        if isinstance(e, Cast):
            return replace(e, expr=rewrite(e.expr))
        if isinstance(e, ArrayAccess):
            return replace(e, arr=rewrite(e.arr), idx=rewrite(e.idx))
        if isinstance(e, (TupleExpr, ArrayLit)):
            return replace(e, elems=[rewrite(x) for x in e.elems])
        return e

    new_body_expr = rewrite(pred_cond)
    pred_param = HIRParam(name=x_name, ty=elem_ty, is_ref=False,
                          is_const_ref=True, is_output=False)
    return LambdaExpr(
        captures=[],
        params=[pred_param],
        body=[ReturnStmt(value=new_body_expr)],
        ty=NamedType("Bool"),
    )


def _get_it_name_from_stmt(s: StmtIR) -> str | None:
    """从赋值-like stmt 抽取左手变量名。"""
    if isinstance(s, LetStmt): return s.var.name
    if isinstance(s, AssignStmt) and isinstance(s.target, Var): return s.target.name
    if isinstance(s, ExprStmt) and isinstance(s.expr, Call):
        if _callee_str(s.expr) == "operator=" and len(s.expr.args) >= 1:
            lhs = _strip_cast(s.expr.args[0])
            if isinstance(lhs, Var): return lhs.name
    return None


def _find_pred_ifstmt_in_body(body: list[StmtIR], it_name: str) -> IfStmt | None:
    """在 body 中找第一个 IfStmt，其 cond 涉及 `*it` / `it->field`。
    返回 IfStmt 对象（caller 依 form A/B 决定 pred 是否需反转），找不到返回 None。"""
    for s in body:
        if isinstance(s, IfStmt):
            if _refers_to_deref(s.cond, it_name):
                return s
        elif isinstance(s, BlockStmt):
            r = _find_pred_ifstmt_in_body(s.stmts, it_name)
            if r is not None: return r
    return None


def _is_field_assign_to_it(s: StmtIR, it_name: str
                           ) -> tuple[str, ExprIR] | None:
    """检查 stmt 是否为 `__deref__(it).field = expr` 形态的 mutator。
    返回 (field_name, value_expr)；非该形态返回 None。

    AssignStmt(target=FieldAccess(obj=Call("__deref__", [Var(it)]), field=...))
    或 AssignStmt(target=Call("StdMap.get!", [m, k]), ...) 不在此处理（独立模式）。
    """
    if not isinstance(s, AssignStmt):
        return None
    tgt = _strip_cast(s.target)
    if not isinstance(tgt, FieldAccess):
        return None
    obj = _strip_cast(tgt.obj)
    # 形态 A：FieldAccess(obj=Call("__deref__", [Var(it)]), field)
    if isinstance(obj, Call) and _callee_str(obj) == "__deref__" \
            and len(obj.args) == 1:
        inner = _strip_cast(obj.args[0])
        if isinstance(inner, Var) and inner.name == it_name:
            return tgt.field_name, s.value
    # 形态 B：FieldAccess(obj=Call("operator->", [Var(it)]), field)
    if isinstance(obj, Call) and _callee_str(obj) == "operator->" \
            and len(obj.args) == 1:
        inner = _strip_cast(obj.args[0])
        if isinstance(inner, Var) and inner.name == it_name:
            return tgt.field_name, s.value
    return None


def _extract_head_mutators(body: list[StmtIR], it_name: str
                           ) -> tuple[list[tuple[str, ExprIR]], int]:
    """从 body 头部连续提取 `it->field = expr` mutators。
    返回 (mutators, consumed_count)：mutators 列表 + 消耗的 stmt 数（即 IfStmt 之前）。
    """
    mutators: list[tuple[str, ExprIR]] = []
    i = 0
    while i < len(body):
        m = _is_field_assign_to_it(body[i], it_name)
        if m is None:
            break
        mutators.append(m)
        i += 1
    return mutators, i


def _is_pure_filter_body(body: list[StmtIR], it_name: str, out_name: str | None) -> bool:
    """检查 compact-erase 循环 body 是否"纯筛选"——除 IfStmt / 迭代器推进外无其他 side-effect。

    纯筛选 body 允许的 stmt：
    - IfStmt（包含 pred 判断 + copy-and-advance）
    - ExprStmt(Call("operator++/--", [it or out]))
    - AssignStmt 给 it / out 的
    - ExprStmt(Call("operator=", [it or out, ...]))
    禁止的 stmt（非纯化信号）：
    - 任何目标是 `*it` / `it->field` 的 ExprStmt（如 fdiv_r, fdiv_q 对 it->second 的原地修改）
    """
    def is_allowed_stmt(s: StmtIR) -> bool:
        if isinstance(s, IfStmt): return True
        # operator++/-- / operator= 对 it/out 自身（AssignStmt 或 ExprStmt(Call)）
        target: ExprIR | None = None
        call: Call | None = None
        if isinstance(s, ExprStmt) and isinstance(s.expr, Call):
            call = s.expr
            cs = _callee_str(call)
            if cs and any(cs.startswith(op) for op in ("operator++", "operator--", "operator=")):
                if call.args:
                    target = _strip_cast(call.args[0])
        elif isinstance(s, AssignStmt):
            target = _strip_cast(s.target)
        if target is None: return False
        # target 必须是 it 或 out 变量自身（不是 deref/field）
        if isinstance(target, Var):
            if target.name == it_name: return True
            if out_name is not None and target.name == out_name: return True
        return False

    return all(is_allowed_stmt(s) for s in body)


def _refers_to_deref(e: ExprIR, it_name: str) -> bool:
    """e 内是否引用 `*it_name` 或 `it_name->field`。"""
    if isinstance(e, Call):
        cs = _callee_str(e)
        if cs in ("operator*", "operator->"):
            for a in e.args:
                inner = _strip_cast(a)
                if isinstance(inner, Var) and inner.name == it_name: return True
        for a in e.args:
            if _refers_to_deref(a, it_name): return True
    if isinstance(e, Cast): return _refers_to_deref(e.expr, it_name)
    if isinstance(e, BinOp):
        return _refers_to_deref(e.lhs, it_name) or _refers_to_deref(e.rhs, it_name)
    if isinstance(e, UnaryOp): return _refers_to_deref(e.operand, it_name)
    if isinstance(e, CondExpr):
        return (_refers_to_deref(e.cond, it_name)
                or _refers_to_deref(e.then_e, it_name)
                or _refers_to_deref(e.else_e, it_name))
    if isinstance(e, FieldAccess): return _refers_to_deref(e.obj, it_name)
    if isinstance(e, ArrayAccess):
        return _refers_to_deref(e.arr, it_name) or _refers_to_deref(e.idx, it_name)
    if isinstance(e, (TupleExpr, ArrayLit)):
        return any(_refers_to_deref(x, it_name) for x in e.elems)
    return False


def _match_filter_loop_A(stmts: list[StmtIR], idx: int) -> FilterLoopMatch | None:
    """匹配形态 A（compact-erase 双指针 4 件套）。

    特征：
      s0: rhs 为 Call("<method>.begin", [container])
      s1: rhs 为任意（Var / Cast / Call 皆可，常见 Var(it) 或 Call("construct_iterator", [it])）
      s2: ForStmt or WhileStmt（loop 主体）
      s3: ExprStmt(Call("<method>.erase", ...))  — top-level
    """
    if idx + 3 >= len(stmts): return None
    s0, s1, s2, s3 = stmts[idx:idx+4]

    rhs0 = _stmt_rhs(s0)
    rhs1 = _stmt_rhs(s1)
    if rhs0 is None or rhs1 is None: return None

    # s0.rhs 必须是 Call("<method>.begin")
    rhs0_stripped = _strip_cast(rhs0)
    if not isinstance(rhs0_stripped, Call): return None
    cs0 = _callee_str(rhs0_stripped)
    if not (cs0 and ".begin" in cs0): return None
    container = rhs0_stripped.args[0] if rhs0_stripped.args else None
    if container is None: return None

    # s2 是 loop
    if not isinstance(s2, (ForStmt, WhileStmt)): return None

    # s3 是 top-level erase
    if not _is_erase_stmt(s3): return None

    it_name = _get_it_name_from_stmt(s0)
    out_name = _get_it_name_from_stmt(s1)
    if it_name is None: return None

    # CF-1 阶段 1（mutate-then-filter）：先尝试识别 body 头部的 mutator stmts，
    # 余下应是 pure filter pattern。若 mutator 非空则走 Array.filterMap 路径。
    body_list = list(s2.body)
    mutators, head_consumed = _extract_head_mutators(body_list, it_name)
    rest_body = body_list[head_consumed:]

    if mutators:
        # mutate-then-filter：rest_body 应是纯 if-write-advance 形态
        if not _is_pure_filter_body(rest_body, it_name, out_name):
            return None
        if_stmt = _find_pred_ifstmt_in_body(rest_body, it_name)
        if if_stmt is None: return None
        # 推 elem_ty：从 container 类型剥（Array → elem_ty）
        elem_ty = None
        if isinstance(container, Var):
            ct = container.ty
            if isinstance(ct, ArrayType):
                elem_ty = ct.elem
            elif isinstance(ct, StdMapType):
                # std::map<K, V> 的 elem 是 (K, V) pair
                elem_ty = PairType(ct.key, ct.value)
            elif isinstance(ct, NamedType):
                # 已知 CLPoly 容器 NamedType → elem 类型查表
                _CONTAINER_ELEM = {
                    "SparsePolyZp": PairType(NamedType("UMonomial"), NamedType("Zp")),
                    "SparsePolyZZ": PairType(NamedType("UMonomial"), NamedType("ZZ")),
                    "MvPolyZp": PairType(NamedType("MvMonomial"), NamedType("Zp")),
                    "MvPolyZZ": PairType(NamedType("MvMonomial"), NamedType("ZZ")),
                }
                elem_ty = _CONTAINER_ELEM.get(ct.name)
        pred_lambda = _build_filter_map_lambda(
            mutators, if_stmt.cond, it_name, elem_ty)
        return FilterLoopMatch(
            start_idx=idx, end_idx=idx + 4,
            container=container, pred_lambda=pred_lambda, kind="A-mut",
        )

    # 纯筛选约束（无 mutator）：C++ compact-erase body 内若混有
    # side-effect（如 fdiv_r(it->second, ...)），丢弃识别以避免语义错失。
    if not _is_pure_filter_body(body_list, it_name, out_name):
        return None

    # 找 body 内 IfStmt 作 pred（form A 的 then 分支是 "copy-and-advance"，
    # pred 不反转——保留 cond=true 的元素）
    if_stmt = _find_pred_ifstmt_in_body(body_list, it_name)
    if if_stmt is None: return None

    # B1 续修：从 container 类型推 elem_ty（与 mut 路径一致）
    pure_elem_ty = None
    if isinstance(container, Var):
        ct = container.ty
        if isinstance(ct, ArrayType):
            pure_elem_ty = ct.elem
        elif isinstance(ct, StdMapType):
            pure_elem_ty = PairType(ct.key, ct.value)
        elif isinstance(ct, NamedType):
            _CONTAINER_ELEM_PURE = {
                "SparsePolyZp": PairType(NamedType("UMonomial"), NamedType("Zp")),
                "SparsePolyZZ": PairType(NamedType("UMonomial"), NamedType("ZZ")),
                "MvPolyZp": PairType(NamedType("MvMonomial"), NamedType("Zp")),
                "MvPolyZZ": PairType(NamedType("MvMonomial"), NamedType("ZZ")),
            }
            pure_elem_ty = _CONTAINER_ELEM_PURE.get(ct.name)
    pred_lambda = _build_pred_lambda(if_stmt.cond, it_name, element_ty=pure_elem_ty)
    return FilterLoopMatch(
        start_idx=idx, end_idx=idx + 4,
        container=container, pred_lambda=pred_lambda, kind="A",
    )


def _match_filter_loop_B(stmts: list[StmtIR], idx: int) -> FilterLoopMatch | None:
    """匹配形态 B（classic iter-loop + 条件 erase，单 stmt）。

    特征：
      s = ForStmt with init 含 LetStmt(it, Call("<method>.begin")) 且 body 含 erase，或
      s = WhileStmt with cond 含 Call("<method>.end") 且 body 含 erase。
    """
    if idx >= len(stmts): return None
    s = stmts[idx]
    if not isinstance(s, (ForStmt, WhileStmt)): return None

    container = None
    it_name = None

    if isinstance(s, ForStmt):
        # init 找 LetStmt(it, Call("<method>.begin", [container]))
        for init_s in s.init:
            if isinstance(init_s, LetStmt):
                v = _strip_cast(init_s.value) if init_s.value else None
                if isinstance(v, Call):
                    cs = _callee_str(v)
                    if cs and ".begin" in cs and v.args:
                        container = v.args[0]
                        it_name = init_s.var.name
                        break
    if container is None:
        # 尝试从 cond 找 `.end(container)`
        if s.cond:
            def find_end(e):
                nonlocal container, it_name
                stripped = _strip_cast(e)
                if isinstance(stripped, Call):
                    cs = _callee_str(stripped)
                    if cs and ".end" in cs and stripped.args:
                        container = stripped.args[0]
                        return
                    for a in stripped.args:
                        if container is None: find_end(a)
                if isinstance(stripped, Cast):
                    find_end(stripped.expr)
                if isinstance(stripped, BinOp):
                    find_end(stripped.lhs); find_end(stripped.rhs)
            find_end(s.cond)
        # it_name 从 cond 里第一个 Var（operator!= 的左参）找
        if container is not None and it_name is None:
            def find_it_var(e):
                nonlocal it_name
                stripped = _strip_cast(e)
                if isinstance(stripped, Var):
                    it_name = stripped.name
                    return
                if isinstance(stripped, Call):
                    cs = _callee_str(stripped)
                    if cs == "operator!=" and stripped.args:
                        find_it_var(stripped.args[0])
                        return
                    for a in stripped.args:
                        if it_name is None: find_it_var(a)
                if isinstance(stripped, Cast): find_it_var(stripped.expr)
            if s.cond: find_it_var(s.cond)

    if container is None or it_name is None: return None
    if not _deep_has_erase(s.body): return None

    # 找 body 内 IfStmt
    if_stmt = _find_pred_ifstmt_in_body(list(s.body), it_name)
    if if_stmt is None: return None

    # P0-a：form B 的 IfStmt 结构为 `if (cond) {erase} else {++it}` 或
    # `if (cond) {++it} else {erase}`。filter 保留的是 non-erase 分支的元素：
    #   - erase 在 then → pred = NOT cond
    #   - erase 在 else → pred = cond
    erase_in_then = _deep_has_erase(list(if_stmt.then_body))
    erase_in_else = _deep_has_erase(list(if_stmt.else_body))
    if erase_in_then and not erase_in_else:
        pred_expr = UnaryOp(op="!", operand=if_stmt.cond, ty=NamedType("Bool"))
    elif erase_in_else and not erase_in_then:
        pred_expr = if_stmt.cond
    else:
        # 模糊情况（两侧都或都不含 erase）— 不识别
        return None

    # B1 续修：从 container 类型推 elem_ty
    b_elem_ty = None
    if isinstance(container, Var):
        ct = container.ty
        if isinstance(ct, ArrayType):
            b_elem_ty = ct.elem
        elif isinstance(ct, StdMapType):
            b_elem_ty = PairType(ct.key, ct.value)
        elif isinstance(ct, NamedType):
            _CONTAINER_ELEM_B = {
                "SparsePolyZp": PairType(NamedType("UMonomial"), NamedType("Zp")),
                "SparsePolyZZ": PairType(NamedType("UMonomial"), NamedType("ZZ")),
                "MvPolyZp": PairType(NamedType("MvMonomial"), NamedType("Zp")),
                "MvPolyZZ": PairType(NamedType("MvMonomial"), NamedType("ZZ")),
            }
            b_elem_ty = _CONTAINER_ELEM_B.get(ct.name)
    pred_lambda = _build_pred_lambda(pred_expr, it_name, element_ty=b_elem_ty)
    kind = f"B-{type(s).__name__}"
    return FilterLoopMatch(
        start_idx=idx, end_idx=idx + 1,
        container=container, pred_lambda=pred_lambda, kind=kind,
    )


def _match_filter_loop(stmts: list[StmtIR], idx: int) -> FilterLoopMatch | None:
    """先尝试形态 A（消费 4 条），再尝试形态 B（消费 1 条）。"""
    m = _match_filter_loop_A(stmts, idx)
    if m is not None: return m
    return _match_filter_loop_B(stmts, idx)


def _filter_loop_to_assign(m: FilterLoopMatch) -> AssignStmt:
    """把 FilterLoopMatch 转为 AssignStmt(v, Call("Array.filter|filterMap", [v, pred]))。
    StdMap 容器用 StdMap.filter（B1 续修）。"""
    is_mut = m.kind == "A-mut"
    cont_ty = getattr(m.container, 'ty', None) if isinstance(m.container, Var) else None
    if isinstance(cont_ty, StdMapType):
        callee = "StdMap.filterMap" if is_mut else "StdMap.filter"
    else:
        callee = "Array.filterMap" if is_mut else "Array.filter"
    return AssignStmt(
        target=m.container,
        value=Call(
            callee=callee,
            args=[m.container, m.pred_lambda],
            ty=None,
        ),
    )


def _build_filter_map_lambda(
    mutators: list[tuple[str, ExprIR]],
    filter_cond: ExprIR,
    it_name: str,
    elem_ty: TypeIR | None,
) -> LambdaExpr:
    """构造 mutate-then-filter 的 lambda（用于 Array.filterMap）。

    生成形态：
      fun (__x: ElemTy) =>
        let __m_<f0> := <expr0_subst>     -- 每个 mutator 一个 let
        let __m_<f1> := <expr1_subst>
        let __x_mut := { __x with field0 := __m_f0, field1 := __m_f1, ... }
        if <filter_cond_subst> then Some __x_mut else None

    其中 expr_subst 是把 `__deref__(it).field` 替换为 `__x.field` 或 `__m_<field>`
    （后者用于"mutator 后续读自身已修改的 field"）。

    elem_ty 用于决定 record-update 形态：PairType → TupleExpr(fst, snd)；
    其它 → Call("__ctor__pair", ...) 退化（Pass 8 codegen 兜底）。
    """
    x_name = "__x"
    elem_ty = elem_ty if elem_ty else UnknownType("")

    # 把 `__deref__(it).field` / `it->field` / `*it` → `__x.field` 或 `__x` 的 rewrite
    # mutator_field_to_var：mutated 字段名 → 该字段当前最新 Var（让序列内引用拿到 mutated 后的值）
    mutator_field_to_var: dict[str, Var] = {}

    def rewrite_iter_ref(e: ExprIR) -> ExprIR:
        if isinstance(e, Call):
            cs = _callee_str(e)
            # `*it` / `Iterator.deref!(it)` → __x（整体）
            if cs in ("operator*", "Iterator.deref!", "__deref__") and len(e.args) == 1:
                inner = _strip_cast(e.args[0])
                if isinstance(inner, Var) and inner.name == it_name:
                    return Var(name=x_name, version=0, ty=elem_ty)
            return Call(callee=e.callee,
                        args=[rewrite_iter_ref(a) for a in e.args], ty=e.ty)
        if isinstance(e, FieldAccess):
            obj = _strip_cast(e.obj)
            # `it->field` 或 `__deref__(it).field` 或 `*it.field`
            if isinstance(obj, Call):
                cs = _callee_str(obj)
                if cs in ("operator->", "__deref__", "Iterator.deref!") \
                        and len(obj.args) == 1:
                    inner = _strip_cast(obj.args[0])
                    if isinstance(inner, Var) and inner.name == it_name:
                        # 若该字段已被 mutator 修改过，用最新版本
                        if e.field_name in mutator_field_to_var:
                            return mutator_field_to_var[e.field_name]
                        return FieldAccess(
                            obj=Var(name=x_name, version=0, ty=elem_ty),
                            field_name=e.field_name, ty=e.ty)
            return FieldAccess(obj=rewrite_iter_ref(e.obj),
                               field_name=e.field_name, ty=e.ty)
        if isinstance(e, BinOp):
            return BinOp(op=e.op, lhs=rewrite_iter_ref(e.lhs),
                         rhs=rewrite_iter_ref(e.rhs), ty=e.ty)
        if isinstance(e, UnaryOp):
            return UnaryOp(op=e.op, operand=rewrite_iter_ref(e.operand), ty=e.ty)
        if isinstance(e, CondExpr):
            return CondExpr(cond=rewrite_iter_ref(e.cond),
                            then_e=rewrite_iter_ref(e.then_e),
                            else_e=rewrite_iter_ref(e.else_e), ty=e.ty)
        if isinstance(e, Cast):
            return Cast(expr=rewrite_iter_ref(e.expr),
                        source_ty=e.source_ty, target_ty=e.target_ty,
                        cast_kind=e.cast_kind)
        if isinstance(e, ArrayAccess):
            return ArrayAccess(arr=rewrite_iter_ref(e.arr),
                               idx=rewrite_iter_ref(e.idx), ty=e.ty)
        if isinstance(e, (TupleExpr, ArrayLit)):
            return type(e)(elems=[rewrite_iter_ref(x) for x in e.elems],
                           **({"ty": e.ty} if hasattr(e, "ty") else {})
                                if not isinstance(e, ArrayLit) else
                          {"elem_ty": e.elem_ty})
        return e

    # 构造 lambda body
    body: list[StmtIR] = []
    for field_name, value_expr in mutators:
        # `__m_<field>_<idx>` 用全局 idx 防同 field 多次 mutate
        m_var_name = f"__m_{field_name}_{len(mutator_field_to_var)}"
        m_var = Var(name=m_var_name, version=0, ty=UnknownType(""))
        body.append(LetStmt(
            var=m_var, ty=UnknownType(""),
            value=rewrite_iter_ref(value_expr),
        ))
        mutator_field_to_var[field_name] = m_var

    # 构造 mutated x（PairType → TupleExpr(fst, snd)）
    if mutators and isinstance(elem_ty, PairType):
        # PairType: 字段 fst / snd
        new_fst = mutator_field_to_var.get("fst") or mutator_field_to_var.get("first") or \
            FieldAccess(obj=Var(name=x_name, version=0, ty=elem_ty),
                        field_name="fst", ty=elem_ty.fst)
        new_snd = mutator_field_to_var.get("snd") or mutator_field_to_var.get("second") or \
            FieldAccess(obj=Var(name=x_name, version=0, ty=elem_ty),
                        field_name="snd", ty=elem_ty.snd)
        x_mut_expr: ExprIR = TupleExpr(elems=[new_fst, new_snd], ty=elem_ty)
    else:
        # 退化：直接用原 __x（mutator 修改通过 mutator_field_to_var 反映在 cond，
        # 但 x_mut 本身保留原结构 — codegen 阶段处理）
        x_mut_expr = Var(name=x_name, version=0, ty=elem_ty)

    x_mut_var = Var(name="__x_mut", version=0, ty=elem_ty)
    body.append(LetStmt(var=x_mut_var, ty=elem_ty, value=x_mut_expr))

    # filter cond 改写后用于决定 some / none
    cond_subst = rewrite_iter_ref(filter_cond)
    body.append(ReturnStmt(
        value=CondExpr(
            cond=cond_subst,
            then_e=Call(callee="Option.some",
                        args=[x_mut_var], ty=None),
            else_e=Call(callee="Option.none", args=[], ty=None),
            ty=None,
        )
    ))

    pred_param = HIRParam(name=x_name, ty=elem_ty, is_ref=False,
                          is_const_ref=True, is_output=False)
    return LambdaExpr(
        captures=[],
        params=[pred_param],
        body=body,
        ty=NamedType("Option"),
    )


# ============================================================
# Structured-binding desugar
# ============================================================

def _desugar_decomposition(rf: RangeForStmt) -> RangeForStmt:
    """若 decomposition 非空，在 body 头部注入 LetStmt 展开。"""
    if not rf.decomposition: return rf

    # 元素类型（_x 的类型）来自 rf.var_ty；剥 RefType 外壳
    # （range-for 的 C++ 元素类型常为 Pair<A,B>& / const Pair<A,B>&，
    # 原始 TypeIR 是 RefType(inner=PairType(...))；穿透至 PairType 用 .fst/.snd）
    elem_ty = rf.var_ty
    while isinstance(elem_ty, RefType):
        elem_ty = elem_ty.inner
    x_var = rf.var

    prelude: list[StmtIR] = []

    if isinstance(elem_ty, PairType):
        # 2 元素：.fst / .snd
        field_names = ["fst", "snd"]
        field_tys = [elem_ty.fst, elem_ty.snd]
    elif isinstance(elem_ty, TupleType):
        field_names = [f"elem{i}" for i in range(len(elem_ty.elems))]
        field_tys = elem_ty.elems
    else:
        # fallback：直接用 decomposition 里的名字与未知类型
        field_names = [f"elem{i}" for i in range(len(rf.decomposition))]
        field_tys = [UnknownType("")] * len(rf.decomposition)

    for i, dv in enumerate(rf.decomposition):
        fn = field_names[i] if i < len(field_names) else f"elem{i}"
        fty = field_tys[i] if i < len(field_tys) else UnknownType("")
        prelude.append(LetStmt(
            var=Var(name=dv.name, version=0, ty=fty),
            ty=fty,
            value=FieldAccess(obj=x_var, field_name=fn, ty=fty),
        ))

    # P0a 修复（嵌套 destructure-bound ranged-for + inner mutation）：
    # `for (auto& [fac, mult] : c)` 在 C++ 中 fac/mult 是引用绑定——body
    # 中 fac=... 应等价于 c[idx].fst=...。Pass 4 把 destructure 展开为
    # `let fac := __decomp.fst` 后，body 内的 fac mutation 不会反向写回
    # __decomp，导致 P0-1 的 cont latch 用旧 __decomp 写回 → 修改丢失。
    # 解决：is_mutable_ref + decomposition 时，body 尾自动 repack：
    #   `__decomp := __ctor__pair(fac_latest, mult_latest)`（PairType 用 TupleExpr）
    postlude: list[StmtIR] = []
    if rf.is_mutable_ref:
        elems = [Var(name=dv.name, ty=field_tys[i] if i < len(field_tys)
                                       else UnknownType(""))
                 for i, dv in enumerate(rf.decomposition)]
        postlude.append(AssignStmt(
            target=x_var,
            value=TupleExpr(elems=elems, ty=elem_ty),
        ))

    return replace(rf, body=prelude + list(rf.body) + postlude,
                   decomposition=[])


# ============================================================
# 主 walk 递归
# ============================================================

def _walk_and_rewrite(stmts: list[StmtIR]) -> list[StmtIR]:
    """穿透所有容器 stmt，识别 filter-loop + desugar decomposition。"""
    out: list[StmtIR] = []
    i = 0
    while i < len(stmts):
        s = stmts[i]

        # 1) filter-loop 优先匹配（只在当前 stmts 层；不进入 s 本身）
        m = _match_filter_loop(stmts, i)
        if m is not None:
            out.append(_filter_loop_to_assign(m))
            i = m.end_idx
            continue

        # 2) 递归处理子 stmt
        if isinstance(s, RangeForStmt):
            # 先 desugar decomposition，body 再递归
            s2 = _desugar_decomposition(s)
            s2 = replace(s2, body=_walk_and_rewrite(list(s2.body)))
            out.append(s2)
        elif isinstance(s, IfStmt):
            out.append(replace(s,
                then_body=_walk_and_rewrite(list(s.then_body)),
                else_body=_walk_and_rewrite(list(s.else_body)),
            ))
        elif isinstance(s, ForStmt):
            out.append(replace(s,
                init=_walk_and_rewrite(list(s.init)),
                step=_walk_and_rewrite(list(s.step)),
                body=_walk_and_rewrite(list(s.body)),
            ))
        elif isinstance(s, WhileStmt):
            out.append(replace(s, body=_walk_and_rewrite(list(s.body))))
        elif isinstance(s, DoWhileStmt):
            out.append(replace(s, body=_walk_and_rewrite(list(s.body))))
        elif isinstance(s, BlockStmt):
            out.append(replace(s, stmts=_walk_and_rewrite(list(s.stmts))))
        else:
            out.append(s)
        i += 1
    return out


# ============================================================
# Pass 入口
# ============================================================

def iter_recognize_pass(func: HIRFunc) -> HIRFunc:
    """Pass 4：HIR₂ → HIR₃。处理主 body + aux_lambdas[*].body。

    P7-8 修复 + P0-B（第十一轮）：Pass 4 CF-1 生成的 LambdaExpr (Array.filter/
    filterMap pred) 在 Pass 3 lambda_lift 之后才生成 → 没被 lift。这里在 Pass 4
    末尾自己 lift：把 Call 中的 LambdaExpr 提到 aux_lambdas，并**正确处理 capture**：
    复用 Pass 3 的 `_collect_free_vars` + outer typectx 计算捕获，把 captures
    作为 lifted func 前置 params；Call args 仍替换为单个 Var(lifted_name)（与
    Pass 3 行为一致，由 Pass 8 codegen 端基于 qual_type.n_caps 做 partial-app）。
    """
    new_body = _walk_and_rewrite(list(func.body))
    new_aux = [
        replace(aux, body=_walk_and_rewrite(list(aux.body)))
        for aux in func.aux_lambdas
    ]
    # 提取 CF-1 生成的 filter/filterMap lambdas
    extra_aux: list[HIRFunc] = []
    counter = [0]
    main_typectx: dict[str, TypeIR] = {p.name: p.ty for p in func.params}
    # 跨实例命名冲突修复（与 Pass 3、Pass 7 保持一致）：含 instance_suffix
    main_host = func.base_name
    if func.instance_suffix:
        main_host = f"{main_host}_{func.instance_suffix}"
    new_body = _lift_filtermap_lambdas(new_body, main_host,
                                          main_typectx, extra_aux, counter)
    new_aux_2: list[HIRFunc] = []
    for aux in new_aux:
        aux_extra: list[HIRFunc] = []
        aux_typectx: dict[str, TypeIR] = {p.name: p.ty for p in aux.params}
        new_aux_body = _lift_filtermap_lambdas(list(aux.body), aux.base_name,
                                                aux_typectx, aux_extra, counter)
        new_aux_2.append(replace(aux, body=new_aux_body))
        new_aux_2.extend(aux_extra)
    return replace(func, body=new_body,
                    aux_lambdas=new_aux_2 + extra_aux)


# 引入 Pass 3 的 free-var 工具（外部依赖；Pass 3 已 1:1 验证过）
from pass3_lambda_lift import _collect_free_vars as _p3_collect_free_vars
from pass3_lambda_lift import _collect_modified as _p3_collect_modified


def _lift_filtermap_lambdas(stmts: list[StmtIR], host_name: str,
                              typectx: dict[str, TypeIR],
                              extra_aux: list[HIRFunc],
                              counter: list[int]) -> list[StmtIR]:
    """递归扫描 stmts，找 Array.filter/filterMap Call 中嵌入的 LambdaExpr，提取
    为 HIRFunc，加入 extra_aux，Call args 中替换为 Var(lifted_name)。

    捕获处理（P0-B 第十一轮）：lifted func 的 params = capture_params + lambda
    原 params；qual_type 编码 n_caps + modified_captures，与 Pass 3 一致。

    typectx 在递归中累积 LetStmt 引入的 local var 类型，用于决定哪些 free var
    属于 host 可见 capture（vs 全局函数引用）。
    """
    def _try_lift_call(c: Call) -> Call | None:
        """若 c 是 Array.filter/filterMap 且 args[1] 是 LambdaExpr，lift 之；
        返回替换后的 Call；否则返回 None。"""
        if not (isinstance(c.callee, str)
                and c.callee in ("Array.filter", "Array.filterMap",
                                   "StdMap.filter", "StdMap.filterMap")
                and len(c.args) == 2 and isinstance(c.args[1], LambdaExpr)):
            return None
        lam: LambdaExpr = c.args[1]
        counter[0] += 1
        lam_name = f"_lambda_{host_name}_filt{counter[0]}"
        param_names = {p.name for p in lam.params}
        free_vars = _p3_collect_free_vars(lam.body, param_names)
        capture_names = sorted(free_vars & set(typectx.keys()))
        modified = _p3_collect_modified(lam.body)
        modified_captures = sorted(set(capture_names) & modified)
        cap_params = [
            HIRParam(
                name=cn,
                ty=typectx.get(cn, UnknownType("")),
                is_ref=(cn in modified),
                is_const_ref=(cn not in modified),
                is_output=False,
            )
            for cn in capture_names
        ]
        ret_ty = lam.ty if lam.ty else UnknownType("")
        lifted = HIRFunc(
            base_name=lam_name,
            instance_suffix="", mangled_name="",
            qual_type=(f"lambda from filterMap in {host_name} "
                        f"| n_caps={len(cap_params)} "
                        f"| modified_captures={modified_captures}"),
            params=cap_params + list(lam.params),
            ret_ty=ret_ty,
            body=list(lam.body),
            requires=[], aux_lambdas=[],
        )
        extra_aux.append(lifted)
        return Call(callee=c.callee,
                      args=[c.args[0],
                             Var(name=lam_name, version=0,
                                 ty=NamedType("LambdaRef"))],
                      ty=c.ty)

    out: list[StmtIR] = []
    for s in stmts:
        if isinstance(s, AssignStmt) and isinstance(s.value, Call):
            replaced = _try_lift_call(s.value)
            if replaced is not None:
                out.append(AssignStmt(target=s.target, value=replaced))
                continue
        if isinstance(s, LetStmt):
            # P1-filter-letstmt（第十二轮）：LetStmt.value 也可能含 Array.filter
            # （CF-1 当前只生成 AssignStmt 形态，但手写 HIR / 未来 CF-2 可能生成
            # `let r := Array.filter(...)`，防御性 lift）。
            new_value = s.value
            if isinstance(s.value, Call):
                replaced = _try_lift_call(s.value)
                if replaced is not None:
                    new_value = replaced
            # 累积 typectx 后继续（注意：先 lift RHS，后入 typectx——RHS 内 lambda
            # 不应视 s.var 为 capture，符合 C++ 作用域语义）
            typectx[s.var.name] = s.ty
            out.append(replace(s, value=new_value) if new_value is not s.value else s)
            continue
        if isinstance(s, IfStmt):
            out.append(replace(s,
                then_body=_lift_filtermap_lambdas(list(s.then_body),
                                                    host_name, typectx, extra_aux, counter),
                else_body=_lift_filtermap_lambdas(list(s.else_body or []),
                                                    host_name, typectx, extra_aux, counter)))
            continue
        if isinstance(s, WhileStmt):
            out.append(replace(s, body=_lift_filtermap_lambdas(list(s.body),
                host_name, typectx, extra_aux, counter)))
            continue
        if isinstance(s, DoWhileStmt):
            out.append(replace(s, body=_lift_filtermap_lambdas(list(s.body),
                host_name, typectx, extra_aux, counter)))
            continue
        if isinstance(s, ForStmt):
            out.append(replace(s,
                init=_lift_filtermap_lambdas(list(s.init), host_name, typectx,
                                                extra_aux, counter),
                step=_lift_filtermap_lambdas(list(s.step), host_name, typectx,
                                                extra_aux, counter),
                body=_lift_filtermap_lambdas(list(s.body), host_name, typectx,
                                                extra_aux, counter)))
            continue
        if isinstance(s, RangeForStmt):
            # range-for 引入循环变量 + decomposition
            typectx[s.var.name] = s.var_ty
            if s.decomposition:
                for dv in s.decomposition:
                    typectx[dv.name] = UnknownType("")
            out.append(replace(s, body=_lift_filtermap_lambdas(list(s.body),
                host_name, typectx, extra_aux, counter)))
            continue
        if isinstance(s, BlockStmt):
            out.append(replace(s, stmts=_lift_filtermap_lambdas(list(s.stmts),
                host_name, typectx, extra_aux, counter)))
            continue
        out.append(s)
    return out


# ============================================================
# HIR₃ 出口不变量
# ============================================================

def assert_hir3_invariant(func: HIRFunc) -> None:
    """- 所有 RangeForStmt.decomposition 为空（或 None）
    - 不存在未识别的 filter-loop A 形态序列
    """
    def check_stmts(stmts: list[StmtIR], ctx: str):
        i = 0
        while i < len(stmts):
            s = stmts[i]
            if isinstance(s, RangeForStmt):
                if s.decomposition:
                    raise TranslationError(
                        pass_name="iter_recognize",
                        func_name=func.base_name,
                        reason=f"{ctx}: RangeForStmt with non-empty decomposition at stmt[{i}]",
                    )
                check_stmts(list(s.body), f"{ctx}/rf[{i}]")
            elif isinstance(s, IfStmt):
                check_stmts(list(s.then_body), f"{ctx}/if.t[{i}]")
                check_stmts(list(s.else_body), f"{ctx}/if.e[{i}]")
            elif isinstance(s, ForStmt):
                check_stmts(list(s.init), f"{ctx}/for.init[{i}]")
                check_stmts(list(s.step), f"{ctx}/for.step[{i}]")
                check_stmts(list(s.body), f"{ctx}/for.body[{i}]")
            elif isinstance(s, (WhileStmt, DoWhileStmt)):
                check_stmts(list(s.body), f"{ctx}/loop[{i}]")
            elif isinstance(s, BlockStmt):
                check_stmts(list(s.stmts), f"{ctx}/blk[{i}]")
            # 未识别的 A 形态 4 件套：只有 Pass 4 在顶层 stmts 消费；若出现 →
            # 说明识别器漏了。检查：i..i+3 四件套 + _match_filter_loop_A 成功
            m = _match_filter_loop_A(stmts, i)
            if m is not None:
                raise TranslationError(
                    pass_name="iter_recognize",
                    func_name=func.base_name,
                    reason=f"{ctx}: unrecognized filter-loop A at stmts[{i}..{i+3}]",
                )
            i += 1

    check_stmts(list(func.body), "body")
    for aux in func.aux_lambdas:
        check_stmts(list(aux.body), f"aux({aux.base_name})")
