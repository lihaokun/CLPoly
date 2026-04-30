"""Pass 8 codegen 全 corpus smoke：65 + factorize x3 = 67 函数。

不验证 lake build（S4 阶段做），仅确保 codegen 不抛异常 + 输出可写。
按 sorry 数量与代码大小做 quick stats。
"""

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
from pass7_loop_lower import loop_lower_pass
from pass8_codegen import codegen_pass
from smoke_pass2_full import _get_factorize_instances

AST_CACHE = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_DIR = Path("/tmp/v2_lean_dump")
OUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT = (V2_ROOT.parent.parent / "docs" / "design" /
            "l1-translation-validation" / "survey" / "pass8-smoke.md")


def main():
    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1-8 codegen smoke on {len(targets)} functions", file=sys.stderr)

    ok = 0; fail = 0
    fail_details: list[tuple[str, str]] = []
    sorry_total = 0
    line_total = 0
    sorry_per_func: list[tuple[str, int, int]] = []  # (name, sorry, lines)

    def process(tag: str, ast):
        nonlocal ok, fail, sorry_total, line_total
        try:
            h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                ref_elim_pass(parse_pass(ast)))))
            h3 = iter_recognize_pass(h)
            h4, _ = operator_resolve_pass(h3)
            mir0 = ssa_build_pass(h4)
            mir1 = loop_lower_pass(mir0)
            lean = codegen_pass(mir1)
            (OUT_DIR / f"{tag}.lean").write_text(lean)
            n_sorry = lean.count("sorry")
            n_lines = lean.count("\n") + 1
            sorry_total += n_sorry
            line_total += n_lines
            sorry_per_func.append((tag, n_sorry, n_lines))
            ok += 1
        except Exception as e:
            fail += 1
            fail_details.append((tag, f"{type(e).__name__}: {str(e)[:200]}"))

    for fn in targets:
        if fn == "factorize":
            for suf, ast in _get_factorize_instances():
                process(f"factorize_{suf}", ast)
            continue
        cache = AST_CACHE / f"{fn}.json"
        if not cache.exists():
            fail += 1; fail_details.append((fn, "no AST cache")); continue
        with open(cache) as f: ast = json.load(f)
        process(fn, ast)

    n_total = ok + fail
    avg_sorry = sorry_total / max(ok, 1)
    avg_lines = line_total / max(ok, 1)
    sorry_per_func.sort(key=lambda x: -x[1])

    lines: list[str] = []
    lines.append("# Pass 1-8 codegen 全量烟测")
    lines.append("")
    lines.append(f"- 目标 mir1：**{n_total}**（factorize 展开 3 实例）")
    lines.append(f"- OK: **{ok}** / FAIL: **{fail}**")
    lines.append(f"- 总 sorry：**{sorry_total}**（avg {avg_sorry:.1f}/func）")
    lines.append(f"- 总输出行数：{line_total}（avg {avg_lines:.1f}/func）")
    lines.append(f"- 输出目录：`/tmp/v2_lean_dump/`")
    lines.append("")
    if fail_details:
        lines.append("## FAIL")
        for tag, r in fail_details[:20]:
            lines.append(f"- `{tag}`: {r}")
        lines.append("")
    lines.append("## sorry 数量 Top 15（残留 hot-spot）")
    for tag, n_s, n_l in sorry_per_func[:15]:
        lines.append(f"- `{tag}`: sorry={n_s}, lines={n_l}")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines))

    print(f"\n=== {ok} OK / {fail} FAIL ===", file=sys.stderr)
    print(f"  total sorry: {sorry_total} (avg {avg_sorry:.1f}/func)", file=sys.stderr)
    print(f"  total lines: {line_total} (avg {avg_lines:.1f}/func)", file=sys.stderr)
    print(f"\nReport: {REPORT}", file=sys.stderr)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
