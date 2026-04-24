"""
Pass 4 iter_recognize 全量烟测 — 65 函数（factorize 展开 3 实例）。

流程：Pass 1+2+3+4；验证 HIR₃ 不变量 + 统计 filter-loop 识别数 + 统计
structured-binding desugar 数。

输出：
- stdout: OK / FAIL
- `docs/.../survey/pass4-smoke.md`: 详细报告
- `/tmp/hir3_dump/`: 每函数 HIR₃ dump（主体 + aux_lambdas）
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    AssignStmt, Call, IfStmt, ForStmt, WhileStmt, DoWhileStmt,
    RangeForStmt, BlockStmt,
    HIRFunc, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass, assert_hir3_invariant
from class_map import TRANSLATION_SCOPE

sys.path.insert(0, str(V2_ROOT / "tests"))
from dump_all_hir import fmt_func
from smoke_pass2_full import _get_factorize_instances

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass4-smoke.md"
DUMP_DIR = Path("/tmp/hir3_dump")


def _count_filter_and_decomp(stmts):
    n_filter = 0
    n_decomp = 0
    for s in stmts:
        if isinstance(s, AssignStmt) and isinstance(s.value, Call):
            if isinstance(s.value.callee, str) and s.value.callee == "Array.filter":
                n_filter += 1
        if isinstance(s, RangeForStmt):
            if s.decomposition: n_decomp += 1
            f2, d2 = _count_filter_and_decomp(s.body); n_filter += f2; n_decomp += d2
        elif isinstance(s, IfStmt):
            f2, d2 = _count_filter_and_decomp(s.then_body); n_filter += f2; n_decomp += d2
            f2, d2 = _count_filter_and_decomp(s.else_body); n_filter += f2; n_decomp += d2
        elif isinstance(s, ForStmt):
            for xs in (s.init, s.step, s.body):
                f2, d2 = _count_filter_and_decomp(xs); n_filter += f2; n_decomp += d2
        elif isinstance(s, (WhileStmt, DoWhileStmt)):
            f2, d2 = _count_filter_and_decomp(s.body); n_filter += f2; n_decomp += d2
        elif isinstance(s, BlockStmt):
            f2, d2 = _count_filter_and_decomp(s.stmts); n_filter += f2; n_decomp += d2
    return n_filter, n_decomp


def fmt_hir3(func: HIRFunc) -> str:
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
    print(f"Pass 1+2+3+4 smoke on {len(targets)} functions...", file=sys.stderr)

    ok = 0; fail = 0; fail_details = []
    per_func = {}
    total_filter = 0
    total_decomp_remaining = 0

    def process(tag, ast, dump_name):
        nonlocal ok, fail, total_filter, total_decomp_remaining
        try:
            hir2 = lambda_lift_pass(ref_elim_pass(parse_pass(ast)))
            hir3 = iter_recognize_pass(hir2)
            assert_hir3_invariant(hir3)

            nf_main, nd_main = _count_filter_and_decomp(hir3.body)
            nf_aux = sum(_count_filter_and_decomp(a.body)[0] for a in hir3.aux_lambdas)
            nd_aux = sum(_count_filter_and_decomp(a.body)[1] for a in hir3.aux_lambdas)

            total_filter += nf_main + nf_aux
            total_decomp_remaining += nd_main + nd_aux

            per_func[tag] = {"status": "OK", "n_filter": nf_main + nf_aux}
            (DUMP_DIR / f"{dump_name}.txt").write_text(fmt_hir3(hir3))
            ok += 1
        except TranslationError as e:
            fail += 1
            fail_details.append((tag, f"TranslationError: {e.reason}"))
            per_func[tag] = {"status": "FAIL", "reason": e.reason}
        except Exception as e:
            fail += 1
            tb = traceback.format_exc()[:500]
            fail_details.append((tag, f"{type(e).__name__}: {e}"))
            per_func[tag] = {"status": "FAIL", "reason": f"{type(e).__name__}: {e}"}

    for fn in targets:
        if fn == "factorize":
            instances = _get_factorize_instances()
            if not instances:
                fail += 1
                fail_details.append((fn, "clang returned 0 factorize instances"))
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

    # ======== 报告 ========
    EXPECTED_FILTER = 4  # P0-b 修复后：__upoly_mod_coeff + __hensel_step_linear 的 compact-erase 因 body impure（内含 fdiv_r/fdiv_q side-effect）被拒识，原 4 件套保留于 HIR₃（技术债，留 Pass 5+ 处理）。纯 filter-loop：__hensel_step ×2 + __extract_monomial_content ×2 = 4
    hosts_with_filter = [(t, d) for t, d in per_func.items()
                         if d.get("status") == "OK" and d["n_filter"] > 0]

    lines = []
    lines.append("# Pass 1 + 2 + 3 + 4 全量烟测")
    lines.append("")
    lines.append(f"- 目标函数数（TRANSLATION_SCOPE）：**{len(targets)}**（factorize 展开 3 实例 → {len(per_func)} HIRs）")
    lines.append(f"- OK: **{ok}** / FAIL: **{fail}**")
    lines.append("")
    lines.append(f"- 识别 filter-loop 总数: **{total_filter}**（预期 {EXPECTED_FILTER}）")
    lines.append(f"- 残留 `RangeForStmt.decomposition`: **{total_decomp_remaining}**（预期 0）")
    lines.append(f"- 有 filter-loop 的宿主: **{len(hosts_with_filter)}**")
    lines.append("")

    if fail_details:
        lines.append("## FAIL 列表")
        lines.append("")
        for fn, reason in fail_details:
            lines.append(f"- `{fn}`: {reason}")
        lines.append("")

    if hosts_with_filter:
        lines.append("## Filter-loop 识别详情")
        lines.append("")
        lines.append("| 宿主函数 | 识别数 |")
        lines.append("|---|---|")
        for t, d in sorted(hosts_with_filter):
            lines.append(f"| `{t}` | {d['n_filter']} |")
        lines.append("")

    lines.append("## 核对")
    lines.append("")
    if total_filter == EXPECTED_FILTER and total_decomp_remaining == 0:
        lines.append("✅ filter-loop = 预期 & decomposition 全清空")
    else:
        lines.append(f"⚠️ 差异：filter {total_filter} vs 预期 {EXPECTED_FILTER}；残留 decomp {total_decomp_remaining}")

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    print(f"\n=== {ok} OK / {fail} FAIL ===", file=sys.stderr)
    print(f"  filter-loops recognized: {total_filter} (expected {EXPECTED_FILTER})", file=sys.stderr)
    print(f"  remaining decomposition: {total_decomp_remaining} (expected 0)", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    print(f"HIR₃ dumps: {DUMP_DIR}", file=sys.stderr)

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
