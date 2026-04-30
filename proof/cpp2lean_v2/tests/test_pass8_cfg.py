"""Pass 8 S2 单测：emit_cfg + 终止子 + merge BB lambda。"""

from __future__ import annotations
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, UnknownType,
    Var, Lit, BinOp, Call,
    LetStmt, RequireStmt, PhiStmt,
    JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG,
)
from pass8_codegen import emit_cfg, EmitCtx


def _make_cfg(entry: int, blocks_dict: dict[int, BasicBlock]) -> CFG:
    cfg = CFG(entry=entry, blocks=blocks_dict)
    cfg.rebuild_edges()
    return cfg


# ============================================================
# Case 1: 单 BB 直接 return
# ============================================================

def test_single_bb_return():
    """entry → return 42"""
    bb0 = BasicBlock(
        bb_id=0, stmts=[
            LetStmt(var=Var("x", 0, BaseType.INT64), ty=BaseType.INT64,
                      value=Lit(42, BaseType.INT64)),
        ],
        terminator=ReturnTerm(value=Var("x", 0, BaseType.INT64)),
    )
    cfg = _make_cfg(0, {0: bb0})
    out = emit_cfg(cfg)
    print("--- single_bb ---")
    print(out)
    assert "let x : Int64 := (42 : Int64)" in out
    assert out.rstrip().endswith("x")  # 最后一行是 return value
    print("PASS test_single_bb_return")


# ============================================================
# Case 2: 菱形 if-else（merge BB lambda）
# ============================================================

def test_diamond_merge():
    """
    BB0: cond_jump (a > 0) → BB1, BB2
    BB1: jump BB3 (phi src x_1 = 1)
    BB2: jump BB3 (phi src x_2 = 2)
    BB3: phi x_3 := {1: x_1, 2: x_2}; return x_3
    """
    a = Var("a", 0, BaseType.INT64)
    x1 = Var("x", 1, BaseType.INT64)
    x2 = Var("x", 2, BaseType.INT64)
    x3 = Var("x", 3, BaseType.INT64)
    bb0 = BasicBlock(
        bb_id=0, stmts=[],
        terminator=CondJumpTerm(
            cond=BinOp(">", a, Lit(0, BaseType.INT64), ty=BaseType.BOOL),
            then_bb=1, else_bb=2),
    )
    bb1 = BasicBlock(
        bb_id=1, stmts=[
            LetStmt(var=x1, ty=BaseType.INT64, value=Lit(1, BaseType.INT64)),
        ],
        terminator=JumpTerm(target=3),
    )
    bb2 = BasicBlock(
        bb_id=2, stmts=[
            LetStmt(var=x2, ty=BaseType.INT64, value=Lit(2, BaseType.INT64)),
        ],
        terminator=JumpTerm(target=3),
    )
    bb3 = BasicBlock(
        bb_id=3, stmts=[
            PhiStmt(target=x3, ty=BaseType.INT64,
                      sources={1: x1, 2: x2}),
        ],
        terminator=ReturnTerm(value=x3),
    )
    cfg = _make_cfg(0, {0: bb0, 1: bb1, 2: bb2, 3: bb3})
    out = emit_cfg(cfg)
    print("--- diamond_merge ---")
    print(out)
    assert "let bb_3 := fun x_3 =>" in out
    # then 分支调用 bb_3 with x_1
    assert "bb_3 x_1" in out
    # else 分支调用 bb_3 with x_2
    assert "bb_3 x_2" in out
    # if cond
    assert "if (a > (0 : Int64)) then" in out
    # merge body: x_3 直接 return
    assert "x_3" in out.split("=>")[1] if "=>" in out else False
    print("PASS test_diamond_merge")


# ============================================================
# Case 3: BB 含 TailCallTerm（自调用）
# ============================================================

def test_tailcall():
    """
    BB0: cond → BB1 (back, tail call) | BB2 (return)
    """
    k = Var("k", 0, BaseType.INT64)
    bb0 = BasicBlock(
        bb_id=0, stmts=[],
        terminator=CondJumpTerm(
            cond=BinOp(">", k, Lit(0, BaseType.INT64), ty=BaseType.BOOL),
            then_bb=1, else_bb=2),
    )
    bb1 = BasicBlock(
        bb_id=1, stmts=[],
        terminator=TailCallTerm(
            target_func="_loop_foo_0",
            args=[BinOp("-", k, Lit(1, BaseType.INT64), ty=BaseType.INT64)]),
    )
    bb2 = BasicBlock(
        bb_id=2, stmts=[],
        terminator=ReturnTerm(value=Lit(0, BaseType.INT64)),
    )
    cfg = _make_cfg(0, {0: bb0, 1: bb1, 2: bb2})
    out = emit_cfg(cfg)
    print("--- tailcall ---")
    print(out)
    # tail call 添加 _ir 后缀
    assert "_loop_foo_0_ir" in out
    # 返回 0
    assert "(0 : Int64)" in out
    print("PASS test_tailcall")


# ============================================================
# Case 4: 嵌套 merge（菱形里再嵌菱形）
# ============================================================

def test_nested_merges():
    """
    BB0 → BB1 / BB2
    BB1 → BB3 / BB4
    BB2 → BB5 (merge of 0→2 single + 5 may also merge from 4)
    BB3 → BB5
    BB4 → BB5
    BB5: merge → return
    """
    cond1 = Var("c1", 0, BaseType.BOOL)
    cond2 = Var("c2", 0, BaseType.BOOL)
    x_a = Var("x", 1, BaseType.INT64)  # from BB3
    x_b = Var("x", 2, BaseType.INT64)  # from BB4
    x_c = Var("x", 3, BaseType.INT64)  # from BB2
    x_m = Var("x", 4, BaseType.INT64)  # merge in BB5

    bb0 = BasicBlock(0, [], CondJumpTerm(cond=cond1, then_bb=1, else_bb=2))
    bb1 = BasicBlock(1, [], CondJumpTerm(cond=cond2, then_bb=3, else_bb=4))
    bb2 = BasicBlock(2, [
        LetStmt(var=x_c, ty=BaseType.INT64, value=Lit(30, BaseType.INT64)),
    ], JumpTerm(target=5))
    bb3 = BasicBlock(3, [
        LetStmt(var=x_a, ty=BaseType.INT64, value=Lit(10, BaseType.INT64)),
    ], JumpTerm(target=5))
    bb4 = BasicBlock(4, [
        LetStmt(var=x_b, ty=BaseType.INT64, value=Lit(20, BaseType.INT64)),
    ], JumpTerm(target=5))
    bb5 = BasicBlock(5, [
        PhiStmt(target=x_m, ty=BaseType.INT64,
                  sources={2: x_c, 3: x_a, 4: x_b}),
    ], ReturnTerm(value=x_m))
    cfg = _make_cfg(0, {0: bb0, 1: bb1, 2: bb2, 3: bb3, 4: bb4, 5: bb5})
    out = emit_cfg(cfg)
    print("--- nested_merges ---")
    print(out)
    # 仅 BB5 是 merge（preds = {2,3,4}），其余皆 single-pred → inline
    assert "let bb_5 := fun x_4 =>" in out
    # 三个调用点：
    assert "bb_5 x_3" in out  # from BB2
    assert "bb_5 x_1" in out  # from BB3
    assert "bb_5 x_2" in out  # from BB4
    # BB1 的 if-else 嵌入在 BB0 的 then 分支
    assert out.count("if ") == 2  # outer + inner
    print("PASS test_nested_merges")


# ============================================================
# Runner
# ============================================================

if __name__ == "__main__":
    tests = [
        test_single_bb_return,
        test_diamond_merge,
        test_tailcall,
        test_nested_merges,
    ]
    passed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except AssertionError as e:
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:
            import traceback
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
    print(f"\n=== {passed}/{len(tests)} passed ===")
    sys.exit(0 if passed == len(tests) else 1)
