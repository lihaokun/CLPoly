"""
MIR ir_types 数据结构单元测试。

覆盖：
- T1: PhiStmt / LetStmt / RequireStmt 构造
- T2: Terminator 4 种（Jump / CondJump / Return / TailCall）
- T3: BasicBlock 组装 + stmts + terminator
- T4: CFG entry + blocks + rebuild_edges 计算 preds/succs
- T5: MIRFunc 完整字段 + lean_name 属性
- T6: assert_mir0_invariant 通过最小合法 CFG（1 block + ReturnTerm）
- T7: assert_mir0_invariant 拒绝 CFG=None
- T8: assert_mir0_invariant 拒绝 phi 在 non-phi 之后
- T9: assert_mir0_invariant 拒绝 SSA 重复定义
- T10: assert_mir0_invariant 拒绝残留 AssignStmt
"""

from __future__ import annotations
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))

from ir_types import (
    BaseType, NamedType,
    Var, Lit,
    LetStmt, AssignStmt, RequireStmt,
    HIRParam,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG, MIRFunc,
    assert_mir0_invariant,
    TranslationError,
)


def _v(name: str, version: int = 0, ty=BaseType.INT64) -> Var:
    return Var(name=name, version=version, ty=ty)


# ============================================================
# T1: 基础节点构造
# ============================================================

def test_phi_letreq_construct():
    phi = PhiStmt(
        target=_v("x", 1),
        ty=BaseType.INT64,
        sources={0: _v("x", 0), 1: _v("x", 0)},
    )
    assert phi.target.version == 1
    assert phi.sources[0].name == "x"

    lt = LetStmt(var=_v("y", 1), ty=BaseType.INT64, value=Lit(0, ty=BaseType.INT64))
    assert isinstance(lt.value, Lit)

    rq = RequireStmt(cond=Lit(True, ty=BaseType.BOOL), name="h", source="UB-1")
    assert rq.name == "h"


# ============================================================
# T2: Terminator 4 种
# ============================================================

def test_terminator_kinds():
    j = JumpTerm(target=2)
    assert j.target == 2
    cj = CondJumpTerm(cond=Lit(True, ty=BaseType.BOOL), then_bb=1, else_bb=2)
    assert cj.then_bb == 1 and cj.else_bb == 2
    rt = ReturnTerm(value=Lit(0, ty=BaseType.INT64))
    assert rt.value is not None
    rt_void = ReturnTerm()
    assert rt_void.value is None
    tc = TailCallTerm(target_func="loop_42_ir", args=[Lit(1), Lit(2)])
    assert tc.target_func == "loop_42_ir" and len(tc.args) == 2


# ============================================================
# T3-T4: BasicBlock + CFG + rebuild_edges
# ============================================================

def test_cfg_rebuild_edges():
    # 三块 CFG: entry=0; 0 → 1/2 (CondJump); 1 → return; 2 → return
    bb0 = BasicBlock(bb_id=0, stmts=[],
                     terminator=CondJumpTerm(cond=Lit(True), then_bb=1, else_bb=2))
    bb1 = BasicBlock(bb_id=1, stmts=[], terminator=ReturnTerm(value=Lit(1)))
    bb2 = BasicBlock(bb_id=2, stmts=[], terminator=ReturnTerm(value=Lit(2)))
    cfg = CFG(entry=0, blocks={0: bb0, 1: bb1, 2: bb2})
    cfg.rebuild_edges()
    assert cfg.succs[0] == [1, 2]
    assert cfg.preds[1] == [0]
    assert cfg.preds[2] == [0]
    assert cfg.succs[1] == [] and cfg.succs[2] == []


# ============================================================
# T5: MIRFunc + lean_name
# ============================================================

def test_mirfunc_lean_name():
    bb = BasicBlock(bb_id=0, stmts=[], terminator=ReturnTerm())
    cfg = CFG(entry=0, blocks={0: bb})

    f1 = MIRFunc(base_name="__make_zp", cfg=cfg)
    assert f1.lean_name == "__make_zp_ir"

    f2 = MIRFunc(base_name="factorize", instance_suffix="lex", cfg=cfg)
    assert f2.lean_name == "factorize_lex_ir"


# ============================================================
# T6: 不变量通过
# ============================================================

def test_mir0_invariant_minimal_pass():
    bb = BasicBlock(bb_id=0,
                    stmts=[LetStmt(var=_v("x", 1), ty=BaseType.INT64, value=Lit(42))],
                    terminator=ReturnTerm(value=_v("x", 1)))
    cfg = CFG(entry=0, blocks={0: bb})
    func = MIRFunc(base_name="test_fn", cfg=cfg, ret_ty=BaseType.INT64)
    assert_mir0_invariant(func)  # 不抛错


# ============================================================
# T7-T10: 不变量拒绝违规
# ============================================================

def test_mir0_invariant_no_cfg():
    func = MIRFunc(base_name="bad", cfg=None)
    try:
        assert_mir0_invariant(func)
        assert False, "should have raised"
    except TranslationError as e:
        assert "cfg is None" in e.reason


def test_mir0_invariant_phi_after_nonphi():
    bb = BasicBlock(
        bb_id=0,
        stmts=[
            LetStmt(var=_v("x", 1), ty=BaseType.INT64, value=Lit(0)),
            PhiStmt(target=_v("y", 1), ty=BaseType.INT64,
                    sources={0: _v("y", 0)}),  # 错位！
        ],
        terminator=ReturnTerm(),
    )
    cfg = CFG(entry=0, blocks={0: bb})
    func = MIRFunc(base_name="bad", cfg=cfg)
    try:
        assert_mir0_invariant(func)
        assert False, "should have raised"
    except TranslationError as e:
        assert "PhiStmt after non-phi" in e.reason


def test_mir0_invariant_ssa_redef():
    bb = BasicBlock(
        bb_id=0,
        stmts=[
            LetStmt(var=_v("x", 1), ty=BaseType.INT64, value=Lit(1)),
            LetStmt(var=_v("x", 1), ty=BaseType.INT64, value=Lit(2)),  # 重复定义
        ],
        terminator=ReturnTerm(),
    )
    cfg = CFG(entry=0, blocks={0: bb})
    func = MIRFunc(base_name="bad", cfg=cfg)
    try:
        assert_mir0_invariant(func)
        assert False, "should have raised"
    except TranslationError as e:
        assert "SSA violation" in e.reason


def test_mir0_invariant_disallowed_stmt():
    bb = BasicBlock(
        bb_id=0,
        stmts=[
            AssignStmt(target=_v("x", 0), value=Lit(0)),  # MIR 不允许
        ],
        terminator=ReturnTerm(),
    )
    cfg = CFG(entry=0, blocks={0: bb})
    func = MIRFunc(base_name="bad", cfg=cfg)
    try:
        assert_mir0_invariant(func)
        assert False, "should have raised"
    except TranslationError as e:
        assert "disallowed MIR stmt" in e.reason


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_phi_letreq_construct,
        test_terminator_kinds,
        test_cfg_rebuild_edges,
        test_mirfunc_lean_name,
        test_mir0_invariant_minimal_pass,
        test_mir0_invariant_no_cfg,
        test_mir0_invariant_phi_after_nonphi,
        test_mir0_invariant_ssa_redef,
        test_mir0_invariant_disallowed_stmt,
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
