"""
Pass 5 operator_resolve 全量烟测 — 65 函数（factorize 展开 3 实例）。

流程：Pass 1+2+3+4+5；验证 HIR₄ 不变量 + 统计 gap。

输出：
- stdout: OK / FAIL + gap 统计
- `docs/.../survey/pass5-smoke.md`: 详细报告 + Top gap 清单（B 策略残余）
- `/tmp/hir4_dump/`: 每函数 HIR₄ dump
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path
from collections import Counter

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import HIRFunc, TranslationError
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import (
    operator_resolve_pass, assert_hir4_invariant, GapLog,
)
from class_map import TRANSLATION_SCOPE

sys.path.insert(0, str(V2_ROOT / "tests"))
from dump_all_hir import fmt_func
from smoke_pass2_full import _get_factorize_instances

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass5-smoke.md"
DUMP_DIR = Path("/tmp/hir4_dump")


def fmt_hir4(func: HIRFunc) -> str:
    parts = [fmt_func(func)]
    if func.aux_lambdas:
        parts.append("")
        parts.append(f"## Lifted lambdas ({len(func.aux_lambdas)})")
        for aux in func.aux_lambdas:
            parts.append("")
            parts.append(fmt_func(aux))
    return "\n".join(parts)


def main():
    DUMP_DIR.mkdir(parents=True, exist_ok=True)
    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1+2+3+4+5 smoke on {len(targets)} functions...", file=sys.stderr)

    ok = 0; fail = 0; fail_details = []
    gap_global = Counter()
    gap_by_func = {}

    def process(tag, ast, dump_name):
        nonlocal ok, fail
        try:
            from pass2b_callsite_ref_elim import callsite_ref_elim_pass
            from pass3b_lambda_ref_elim import lambda_ref_elim_pass
            hir3 = iter_recognize_pass(lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast))))))
            hir4, gap = operator_resolve_pass(hir3)
            assert_hir4_invariant(hir4)
            (DUMP_DIR / f"{dump_name}.txt").write_text(fmt_hir4(hir4))
            n_gap = (len(gap.unresolved_op) + len(gap.cast_miss) +
                     len(gap.method_miss) + len(gap.op_miss) +
                     len(gap.constructor_miss))
            gap_by_func[tag] = n_gap
            for x in gap.unresolved_op:    gap_global[("unresolved_op", x)] += 1
            for x in gap.cast_miss:        gap_global[("cast_miss", x)] += 1
            for x in gap.method_miss:      gap_global[("method_miss", x)] += 1
            for x in gap.op_miss:          gap_global[("op_miss", x)] += 1
            for x in gap.constructor_miss: gap_global[("ctor_miss", x)] += 1
            ok += 1
        except TranslationError as e:
            fail += 1; fail_details.append((tag, f"TranslationError: {e.reason}"))
        except Exception as e:
            fail += 1; fail_details.append((tag, f"{type(e).__name__}: {e}"))

    for fn in targets:
        if fn == "factorize":
            instances = _get_factorize_instances()
            if not instances:
                fail += 1
                fail_details.append((fn, "no factorize instances"))
                continue
            for suffix, ast in instances:
                tag = f"factorize_{suffix}"
                process(tag, ast, tag)
            continue
        cache = AST_CACHE_DIR / f"{fn}.json"
        if not cache.exists():
            fail += 1; fail_details.append((fn, "no AST cache")); continue
        with open(cache) as f: ast = json.load(f)
        process(fn, ast, fn)

    # 类别聚合
    gap_by_kind = Counter()
    for (k, _), n in gap_global.items():
        gap_by_kind[k] += n
    total_gap = sum(gap_by_kind.values())

    # 报告
    lines = []
    lines.append("# Pass 1 + 2 + 3 + 4 + 5 全量烟测")
    lines.append("")
    lines.append(f"- 目标函数数：**{len(targets)}**（factorize 3 实例 → {len(gap_by_func)} HIRs）")
    lines.append(f"- OK: **{ok}** / FAIL: **{fail}**")
    lines.append(f"- Gap 总数: **{total_gap}**（B 策略残余，Pass 8 codegen 时输出 sorry）")
    lines.append("")

    if fail_details:
        lines.append("## FAIL 列表")
        lines.append("")
        for fn, reason in fail_details:
            lines.append(f"- `{fn}`: {reason}")
        lines.append("")

    lines.append("## Gap 分类")
    lines.append("")
    lines.append("| 类别 | 数量 |")
    lines.append("|---|---|")
    for k, n in gap_by_kind.most_common():
        lines.append(f"| `{k}` | {n} |")
    lines.append("")

    lines.append("## Top 20 Gap 详情")
    lines.append("")
    lines.append("| 次数 | 类别 | 详情 |")
    lines.append("|---|---|---|")
    for (k, x), n in sorted(gap_global.items(), key=lambda v: -v[1])[:20]:
        lines.append(f"| {n} | {k} | `{x}` |")
    lines.append("")

    lines.append("## Per-function Gap 数")
    lines.append("")
    funcs_with_gap = [(t, n) for t, n in gap_by_func.items() if n > 0]
    funcs_with_gap.sort(key=lambda v: -v[1])
    if funcs_with_gap:
        lines.append("| 函数 | gap |")
        lines.append("|---|---|")
        for t, n in funcs_with_gap[:30]:
            lines.append(f"| `{t}` | {n} |")
        if len(funcs_with_gap) > 30:
            lines.append(f"| ... +{len(funcs_with_gap)-30} more | |")
    else:
        lines.append("✅ 0 gap，所有节点都已 resolve")
    lines.append("")

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    print(f"\n=== {ok} OK / {fail} FAIL ===", file=sys.stderr)
    print(f"  total gap (B 策略残余): {total_gap}", file=sys.stderr)
    print(f"  funcs with gap: {len(funcs_with_gap)}/{len(gap_by_func)}", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    print(f"HIR₄ dumps: {DUMP_DIR}", file=sys.stderr)
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
