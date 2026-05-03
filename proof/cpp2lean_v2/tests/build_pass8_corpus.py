"""Pass 8 S4 续：跑全 corpus → 聚合到 Generated/Corpus.lean → 准备 lake build。

仅做生成，不跑 lake（避免 Python 命令行 PATH 问题）。
"""

from __future__ import annotations
import sys
import json
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
from pass7_loop_lower import loop_lower_pass
from pass8_codegen import codegen_corpus
from smoke_pass2_full import _get_factorize_instances

AST_CACHE = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "Corpus.lean"


def main():
    targets = sorted(TRANSLATION_SCOPE)
    mirs = []
    skipped = []
    for fn in targets:
        if fn == "factorize":
            for suf, ast in _get_factorize_instances():
                tag = f"factorize_{suf}"
                try:
                    h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                        ref_elim_pass(parse_pass(ast)))))
                    h3 = iter_recognize_pass(h)
                    h4, _ = operator_resolve_pass(h3)
                    # 阶段 G+：Pass 5 emit 的 Rng.next_advance / 类似 ref-out
                    # 调用，需 Pass 2b 再处理一遍（output_params 对它们生效）
                    h4 = callsite_ref_elim_pass(h4, callee_filter={"Rng.next_advance"})
                    mir0 = ssa_build_pass(h4)
                    mir1 = loop_lower_pass(mir0)
                    mirs.append(mir1)
                except Exception as e:
                    skipped.append((tag, type(e).__name__))
            continue
        cache = AST_CACHE / f"{fn}.json"
        if not cache.exists():
            skipped.append((fn, "no AST cache"))
            continue
        try:
            with open(cache) as f: ast = json.load(f)
            h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                ref_elim_pass(parse_pass(ast)))))
            h3 = iter_recognize_pass(h)
            h4, _ = operator_resolve_pass(h3)
            # 阶段 G+：Pass 5 emit 的 Rng.next_advance 等 ref-out 调用
            # 需 Pass 2b 再处理一遍（output_params 对它们生效）
            h4 = callsite_ref_elim_pass(h4, callee_filter={"Rng.next_advance"})
            mir0 = ssa_build_pass(h4)
            mir1 = loop_lower_pass(mir0)
            mirs.append(mir1)
        except Exception as e:
            skipped.append((fn, type(e).__name__))

    src = codegen_corpus(mirs)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(src)

    n_lines = src.count("\n") + 1
    n_sorry = src.count("sorry")
    print(f"OK: {len(mirs)} top funcs aggregated", file=sys.stderr)
    print(f"  output: {OUT}", file=sys.stderr)
    print(f"  lines:  {n_lines}", file=sys.stderr)
    print(f"  sorry:  {n_sorry}", file=sys.stderr)
    if skipped:
        print(f"  skipped: {len(skipped)}", file=sys.stderr)
        for s in skipped: print(f"    {s}", file=sys.stderr)


if __name__ == "__main__":
    main()
