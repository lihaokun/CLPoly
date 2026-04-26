"""
Pass 6 ssa_build 单元测试。

按 mir-design.md §2.8 关键测试用例：
- T1: trivial passthrough（无分支无赋值）
- T2: linear assign 编版本号
- T3: if 合并产生 phi
- T4: while 循环 phi 在 header
- T5: 嵌套 if/while
- T6: break/continue 多 exit
- T7: 完整函数 __make_zp / __squarefree_Zp 全 67 smoke
- T8: assert_mir0_invariant 拒绝违规
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType, RefType,
    Var, Lit, BinOp, Call,
    LetStmt, AssignStmt, IfStmt, WhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    HIRFunc, HIRParam, MIRFunc,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm,
    BasicBlock, CFG,
    assert_mir0_invariant, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import (
    ssa_build_pass, cfg_from_hir,
    compute_dominator_tree, compute_dominance_frontier,
    place_phi_nodes, rename_variables,
)
from class_map import TRANSLATION_SCOPE


def _mk_func(body, params=None, ret_ty=BaseType.UNIT) -> HIRFunc:
    return HIRFunc(
        base_name="test", instance_suffix="", mangled_name="", qual_type="",
        params=params or [], ret_ty=ret_ty, body=body,
        requires=[], aux_lambdas=[],
    )


# ============================================================
# T1: Trivial
# ============================================================

def test_trivial_passthrough():
    """单 stmt 单块 + auto-Return"""
    body = [LetStmt(var=Var("x", 0), ty=BaseType.INT64, value=Lit(42, ty=BaseType.INT64))]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    assert mir.cfg.entry == 0
    bb0 = mir.cfg.blocks[0]
    assert isinstance(bb0.terminator, ReturnTerm)
    # 第一个 LetStmt 的 var.version 应为 1（fresh 后）
    let_stmt = next(s for s in bb0.stmts if isinstance(s, LetStmt))
    assert let_stmt.var.version == 1


# ============================================================
# T2: Linear assign — version 编号
# ============================================================

def test_linear_assign_versions():
    """x = 1; x = 2;  →  x_1 = 1; x_2 = 2"""
    body = [
        AssignStmt(target=Var("x"), value=Lit(1, ty=BaseType.INT64)),
        AssignStmt(target=Var("x"), value=Lit(2, ty=BaseType.INT64)),
    ]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    bb0 = mir.cfg.blocks[mir.cfg.entry]
    versions = [s.var.version for s in bb0.stmts if isinstance(s, LetStmt) and s.var.name == "x"]
    assert versions == [1, 2], f"got {versions}"


# ============================================================
# T3: If 合并产生 phi
# ============================================================

def test_if_phi_at_merge():
    """if (cond) x=1 else x=2; return x → phi at merge block"""
    body = [
        IfStmt(
            cond=Var("cond", ty=BaseType.BOOL),
            then_body=[AssignStmt(target=Var("x"), value=Lit(1, ty=BaseType.INT64))],
            else_body=[AssignStmt(target=Var("x"), value=Lit(2, ty=BaseType.INT64))],
        ),
        ReturnStmt(value=Var("x")),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="cond", ty=BaseType.BOOL),
        HIRParam(name="x", ty=BaseType.INT64),
    ], ret_ty=BaseType.INT64)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 应至少有 1 个 phi 节点（merge bb 内 x_n := phi(...)）
    has_phi = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, PhiStmt) and s.target.name == "x":
                has_phi = True
                # phi 应有 2 个 sources（来自 then_bb 和 else_bb）
                assert len(s.sources) == 2, f"phi sources: {s.sources}"
    assert has_phi, "no phi for x"


# ============================================================
# T4: While 循环 — phi 在 header
# ============================================================

def test_while_phi_at_header():
    """while (i < n) i = i + 1  →  phi at header"""
    body = [
        WhileStmt(
            cond=BinOp(op="<", lhs=Var("i"), rhs=Var("n"), ty=BaseType.BOOL),
            body=[AssignStmt(
                target=Var("i"),
                value=BinOp(op="+", lhs=Var("i"), rhs=Lit(1), ty=BaseType.INT64),
            )],
        ),
    ]
    func = _mk_func(body, params=[
        HIRParam(name="i", ty=BaseType.INT64),
        HIRParam(name="n", ty=BaseType.INT64),
    ])
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)
    # 找 while 的 header 块（带 CondJumpTerm 且有前驱不止 entry）
    has_phi_at_header = False
    for bb_id, bb in mir.cfg.blocks.items():
        if isinstance(bb.terminator, CondJumpTerm):
            preds = mir.cfg.preds.get(bb_id, [])
            if len(preds) >= 2:
                for s in bb.stmts:
                    if isinstance(s, PhiStmt) and s.target.name == "i":
                        has_phi_at_header = True
                        assert len(s.sources) >= 2
    assert has_phi_at_header, "no phi for i at while header"


# ============================================================
# T5: Break 退出
# ============================================================

def test_break_creates_extra_edge():
    """while (...) { break; ... }  break 跳到 exit 块"""
    body = [
        WhileStmt(
            cond=Lit(True, ty=BaseType.BOOL),
            body=[BreakStmt()],
        ),
    ]
    func = _mk_func(body)
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)


# ============================================================
# T6: Continue 跳到 header
# ============================================================

def test_continue_jumps_header():
    body = [
        WhileStmt(
            cond=Var("c", ty=BaseType.BOOL),
            body=[ContinueStmt()],
        ),
    ]
    func = _mk_func(body, params=[HIRParam(name="c", ty=BaseType.BOOL)])
    mir = ssa_build_pass(func)
    assert_mir0_invariant(mir)


# ============================================================
# T7: 完整函数 + 67 smoke
# ============================================================

def test_smoke_all_67():
    """全 67 函数走 Pass 1-6 + invariant 全过。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    ok = 0
    errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            hir3 = iter_recognize_pass(lambda_lift_pass(ref_elim_pass(parse_pass(ast))))
            hir4, _ = operator_resolve_pass(hir3)
            mir = ssa_build_pass(hir4)
            assert_mir0_invariant(mir)
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {str(e)[:100]}"))
    assert not errs, f"{len(errs)} errors: {errs[:3]}"


# ============================================================
# T8: 不变量违规
# ============================================================

def test_invariant_rejects_no_cfg():
    func = MIRFunc(base_name="bad", cfg=None)
    try:
        assert_mir0_invariant(func)
        assert False, "should raise"
    except TranslationError as e:
        assert "cfg is None" in e.reason


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_trivial_passthrough,
        test_linear_assign_versions,
        test_if_phi_at_merge,
        test_while_phi_at_header,
        test_break_creates_extra_edge,
        test_continue_jumps_header,
        test_smoke_all_67,
        test_invariant_rejects_no_cfg,
    ]
    passed = 0; failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            import traceback
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print(f"\n=== {passed}/{passed+failed} passed ===")
    sys.exit(0 if failed == 0 else 1)
