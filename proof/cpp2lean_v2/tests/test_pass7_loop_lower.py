"""Pass 7 loop_lower 单测。

按 mir-design.md MIR₁ 不变量 + loop-extraction-design.md 关键测试用例：
- T1: 无循环 passthrough
- T2: 单 while 循环（单 phi、单 exit）
- T3: 嵌套 while 循环
- T4: while + break（多 exit dispatch）
- T5: 67 函数烟测 + invariant
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, Var, Lit, BinOp, Call,
    LetStmt, AssignStmt, WhileStmt, IfStmt, BreakStmt, ReturnStmt,
    HIRFunc, HIRParam, MIRFunc, TailCallTerm, ReturnTerm, JumpTerm,
    CondJumpTerm, assert_mir0_invariant, assert_mir1_invariant,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass3b_lambda_ref_elim import lambda_ref_elim_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import ssa_build_pass
from pass7_loop_lower import loop_lower_pass, _find_back_edges, _compute_idom
from class_map import TRANSLATION_SCOPE


def _full_pipeline(ast: dict) -> MIRFunc:
    h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
        ref_elim_pass(parse_pass(ast)))))
    h3 = iter_recognize_pass(h)
    h4, _ = operator_resolve_pass(h3)
    return ssa_build_pass(h4)


# T1: 无循环 passthrough
def test_no_loop_passthrough():
    """无循环函数：MIR₁ 与 MIR₀ 结构相同（aux_defs 仅含 lifted lambda）。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    f = AST / "__make_zp.json"
    with open(f) as fh: ast = json.load(fh)
    mir0 = _full_pipeline(ast)
    mir1 = loop_lower_pass(mir0)
    # __make_zp 无循环
    n_loops = sum(1 for aux in mir1.aux_defs
                   if aux.base_name.startswith("_loop_"))
    assert n_loops == 0, f"expected 0 loops, got {n_loops}"
    # MIR₁ invariant 通过
    assert_mir1_invariant(mir1)


# T2: 含简单 while 循环的函数
def test_simple_while_loop():
    """`__binomial` 含简单 while；提取后 invariant 通过 + 0 back edge。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    f = AST / "__binomial.json"
    if not f.exists():
        return
    with open(f) as fh: ast = json.load(fh)
    mir0 = _full_pipeline(ast)
    mir1 = loop_lower_pass(mir0)
    # 至少 1 loop 提取
    loops = [a for a in mir1.aux_defs if a.base_name.startswith("_loop_")]
    assert len(loops) >= 1, f"expected ≥1 loop, got {len(loops)}"
    # 主 cfg 0 back edge
    if mir1.cfg is not None:
        idom = _compute_idom(mir1.cfg)
        be = _find_back_edges(mir1.cfg, idom)
        assert len(be) == 0, f"main cfg has {len(be)} back edges"
    # 每个 loop func 内部也无 back edge（loop body 末尾是 TailCallTerm 不是 back edge）
    for lp in loops:
        if lp.cfg:
            idom = _compute_idom(lp.cfg)
            be = _find_back_edges(lp.cfg, idom)
            assert len(be) == 0, f"{lp.base_name}: {len(be)} back edges"
    assert_mir1_invariant(mir1)


# T3: 嵌套 while 循环（多 loops in single func）
def test_nested_loops():
    """`__lll_reduce` 含多重嵌套循环；测试嵌套提取。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    f = AST / "__lll_reduce.json"
    if not f.exists():
        return
    with open(f) as fh: ast = json.load(fh)
    mir0 = _full_pipeline(ast)
    mir1 = loop_lower_pass(mir0)
    loops = [a for a in mir1.aux_defs if a.base_name.startswith("_loop_")]
    assert len(loops) >= 2, f"expected ≥2 loops in __lll_reduce, got {len(loops)}"
    assert_mir1_invariant(mir1)


# T4: while + break（多 exit）
def test_multi_exit_loop():
    """`__edf_Zp` 含 while-true + break，多 exit dispatch。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    f = AST / "__edf_Zp.json"
    if not f.exists():
        return
    with open(f) as fh: ast = json.load(fh)
    mir0 = _full_pipeline(ast)
    mir1 = loop_lower_pass(mir0)
    # 在 loop 提取后，主 cfg 仍要满足 invariant
    assert_mir1_invariant(mir1)
    # tail call 节点应至少 1 个（loop 体内）
    has_tail_call = False
    for aux in mir1.aux_defs:
        if not aux.base_name.startswith("_loop_"): continue
        for bb in aux.cfg.blocks.values():
            if isinstance(bb.terminator, TailCallTerm):
                has_tail_call = True; break
    assert has_tail_call, "no TailCallTerm found in any loop func"


# T5: 67 函数烟测
def test_smoke_all_67():
    """67 函数全过 Pass 7 + invariant。"""
    AST = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
    ok = 0
    errs = []
    for fn in sorted(TRANSLATION_SCOPE):
        f = AST / f"{fn}.json"
        if not f.exists(): continue
        try:
            with open(f) as fh: ast = json.load(fh)
            mir0 = _full_pipeline(ast)
            mir1 = loop_lower_pass(mir0)
            assert_mir1_invariant(mir1)
            ok += 1
        except Exception as e:
            errs.append((fn, f"{type(e).__name__}: {str(e)[:100]}"))
    assert not errs, f"{len(errs)} errors: {errs[:3]}"


if __name__ == "__main__":
    tests = [
        test_no_loop_passthrough,
        test_simple_while_loop,
        test_nested_loops,
        test_multi_exit_loop,
        test_smoke_all_67,
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
