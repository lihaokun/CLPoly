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
    BaseType, NamedType, UnknownType, TypeIR, PairType, TupleType, RefType,
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
    """识别到的 filter-loop 描述。"""
    __slots__ = ("start_idx", "end_idx", "container", "pred_lambda", "kind")

    def __init__(self, start_idx, end_idx, container, pred_lambda, kind):
        self.start_idx = start_idx
        self.end_idx = end_idx  # exclusive
        self.container = container
        self.pred_lambda = pred_lambda
        self.kind = kind  # "A" or "B-For" or "B-While"


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

    # P0-b：纯筛选约束。C++ compact-erase body 内若混有
    # side-effect（如 fdiv_r(it->second, ...)），丢弃识别以避免语义错失。
    if not _is_pure_filter_body(list(s2.body), it_name, out_name):
        return None

    # 找 body 内 IfStmt 作 pred（form A 的 then 分支是 "copy-and-advance"，
    # pred 不反转——保留 cond=true 的元素）
    if_stmt = _find_pred_ifstmt_in_body(list(s2.body), it_name)
    if if_stmt is None: return None

    pred_lambda = _build_pred_lambda(if_stmt.cond, it_name, element_ty=None)
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

    pred_lambda = _build_pred_lambda(pred_expr, it_name, element_ty=None)
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
    """把 FilterLoopMatch 转为 AssignStmt(v, Call("Array.filter", [v, pred]))。"""
    return AssignStmt(
        target=m.container,
        value=Call(
            callee="Array.filter",
            args=[m.container, m.pred_lambda],
            ty=None,
        ),
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

    return replace(rf, body=prelude + list(rf.body), decomposition=[])


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
    """Pass 4：HIR₂ → HIR₃。处理主 body + aux_lambdas[*].body。"""
    new_body = _walk_and_rewrite(list(func.body))
    new_aux = [
        replace(aux, body=_walk_and_rewrite(list(aux.body)))
        for aux in func.aux_lambdas
    ]
    return replace(func, body=new_body, aux_lambdas=new_aux)


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
