"""Pass 7 loop_lower 全量烟测（65 函数 + factorize x3）。"""

from __future__ import annotations
import sys
import json
import traceback
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from class_map import TRANSLATION_SCOPE
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass3b_lambda_ref_elim import lambda_ref_elim_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import ssa_build_pass
from pass7_loop_lower import loop_lower_pass, _compute_idom, _find_back_edges
from ir_types import TailCallTerm, MIRFunc, assert_mir1_invariant
from smoke_pass2_full import _get_factorize_instances

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass7-smoke.md"


def main():
    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1-7 smoke on {len(targets)} functions...", file=sys.stderr)

    ok = 0; fail = 0; fail_details = []
    total_loops = 0
    total_tailcalls = 0
    total_blocks = 0

    def count_recursive(mir: MIRFunc) -> tuple[int, int, int]:
        n_loops = 0; n_tc = 0; n_bb = 0
        if mir.cfg is not None:
            n_bb += len(mir.cfg.blocks)
            for bb in mir.cfg.blocks.values():
                if isinstance(bb.terminator, TailCallTerm):
                    n_tc += 1
        for aux in mir.aux_defs:
            if aux.base_name.startswith("_loop_"):
                n_loops += 1
            sub_loops, sub_tc, sub_bb = count_recursive(aux)
            n_loops += sub_loops; n_tc += sub_tc; n_bb += sub_bb
        return n_loops, n_tc, n_bb

    def process(tag, ast):
        nonlocal ok, fail, total_loops, total_tailcalls, total_blocks
        try:
            h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                ref_elim_pass(parse_pass(ast)))))
            h3 = iter_recognize_pass(h)
            h4, _ = operator_resolve_pass(h3)
            mir0 = ssa_build_pass(h4)
            mir1 = loop_lower_pass(mir0)
            assert_mir1_invariant(mir1)
            # 验证：主 cfg 0 back edge
            if mir1.cfg is not None:
                idom = _compute_idom(mir1.cfg)
                be = _find_back_edges(mir1.cfg, idom)
                if be:
                    raise RuntimeError(f"residual back edges in main cfg: {be}")
            n_loops, n_tc, n_bb = count_recursive(mir1)
            total_loops += n_loops
            total_tailcalls += n_tc
            total_blocks += n_bb
            ok += 1
        except Exception as e:
            fail += 1
            fail_details.append((tag, f"{type(e).__name__}: {str(e)[:150]}"))

    for fn in targets:
        if fn == "factorize":
            for suffix, ast in _get_factorize_instances():
                process(f"factorize_{suffix}", ast)
            continue
        cache = AST_CACHE_DIR / f"{fn}.json"
        if not cache.exists():
            fail += 1; fail_details.append((fn, "no AST cache")); continue
        with open(cache) as f: ast = json.load(f)
        process(fn, ast)

    n_total = ok + fail
    avg_loops = total_loops / max(ok, 1)

    lines = []
    lines.append("# Pass 1-7 全量烟测")
    lines.append("")
    lines.append(f"- 目标 HIRs：**{n_total}**（factorize 展开 3 实例）")
    lines.append(f"- OK: **{ok}** / FAIL: **{fail}**")
    lines.append(f"- 平均 loop 数: **{avg_loops:.1f}** / 函数")
    lines.append(f"- 总 loops 提取: {total_loops}")
    lines.append(f"- 总 TailCallTerm: {total_tailcalls}")
    lines.append(f"- 总 BasicBlock 数（含 aux_defs）: {total_blocks}")
    lines.append("")
    if fail_details:
        lines.append("## FAIL")
        for fn, r in fail_details[:20]:
            lines.append(f"- `{fn}`: {r}")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    print(f"\n=== {ok} OK / {fail} FAIL ===", file=sys.stderr)
    print(f"  total loops: {total_loops} (avg {avg_loops:.1f}/func)", file=sys.stderr)
    print(f"  total tail calls: {total_tailcalls}", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
