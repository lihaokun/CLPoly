"""
Pass 3 lambda_lift 全量烟测 — 65 函数（factorize 展开为 3 实例）。

流程：Pass 1 parse → Pass 2 ref_elim → Pass 3 lambda_lift；
验证 HIR₂ 不变量（body/aux_lambdas 中无 LambdaExpr）+ 统计 lambda 提升效果。

输出：
- stdout: OK / FAIL 统计
- `docs/.../survey/pass3-smoke.md`: 详细报告
- `/tmp/hir2_dump/`: 每函数的 HIR₂ dump（主体 + aux_lambdas）

可以随时重跑作为回归测试。
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path
from collections import defaultdict

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import HIRFunc, TranslationError
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass, assert_hir2_invariant
from class_map import TRANSLATION_SCOPE

sys.path.insert(0, str(V2_ROOT / "tests"))
from dump_all_hir import fmt_func
from smoke_pass2_full import _get_factorize_instances

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass3-smoke.md"
DUMP_DIR = Path("/tmp/hir2_dump")


def fmt_hir2(func: HIRFunc) -> str:
    """格式化 HIR₂（主函数 + aux_lambdas）。"""
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
    print(f"Pass 1+2+3 smoke on {len(targets)} functions...", file=sys.stderr)

    ok_count = 0
    fail_count = 0
    fail_details = []
    per_func = {}

    def process_ast(tag: str, ast: dict, dump_name: str):
        nonlocal ok_count, fail_count
        try:
            hir0 = parse_pass(ast)
            hir1 = ref_elim_pass(hir0)
            hir2 = lambda_lift_pass(hir1)
            assert_hir2_invariant(hir2)

            ok_count += 1
            per_func[tag] = {
                "status": "OK",
                "n_lambdas": len(hir2.aux_lambdas),
                "lambda_names": [l.base_name for l in hir2.aux_lambdas],
                "lambda_param_counts": [len(l.params) for l in hir2.aux_lambdas],
            }
            (DUMP_DIR / f"{dump_name}.txt").write_text(fmt_hir2(hir2))

        except TranslationError as e:
            fail_count += 1
            fail_details.append((tag, f"TranslationError: {e.reason}"))
            per_func[tag] = {"status": "FAIL", "reason": e.reason}
        except AssertionError as e:
            fail_count += 1
            fail_details.append((tag, f"AssertionError: {e}"))
            per_func[tag] = {"status": "FAIL", "reason": str(e)}
        except Exception as e:
            fail_count += 1
            tb = traceback.format_exc()[:500]
            fail_details.append((tag, f"{type(e).__name__}: {e}"))
            per_func[tag] = {"status": "FAIL", "reason": f"{type(e).__name__}: {e}", "traceback": tb}

    for func_name in targets:
        if func_name == "factorize":
            instances = _get_factorize_instances()
            if not instances:
                fail_count += 1
                fail_details.append((func_name, "clang returned 0 factorize instances"))
                per_func[func_name] = {"status": "FAIL", "reason": "0 instances"}
                continue
            for suffix, ast in instances:
                tag = f"factorize_{suffix}"
                process_ast(tag, ast, tag)
            continue

        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            fail_count += 1
            fail_details.append((func_name, "no AST cache"))
            continue
        with open(cache_file) as f:
            ast = json.load(f)
        process_ast(func_name, ast, func_name)

    # ======== 汇总 ========
    hosts_with_lambdas = [(tag, d) for tag, d in per_func.items()
                          if d.get("status") == "OK" and d["n_lambdas"] > 0]
    total_lambdas = sum(d["n_lambdas"] for _, d in hosts_with_lambdas)

    # ======== 输出报告 ========
    lines = []
    lines.append("# Pass 1 + Pass 2 + Pass 3 全量烟测")
    lines.append("")
    lines.append(f"- 目标函数数（TRANSLATION_SCOPE）：**{len(targets)}**（factorize 展开 3 实例 → {len(per_func)} HIRs）")
    lines.append(f"- OK: **{ok_count}**")
    lines.append(f"- FAIL: **{fail_count}**")
    lines.append("")
    lines.append(f"- 有 lambda 的宿主 HIR 数: **{len(hosts_with_lambdas)}**")
    lines.append(f"- 提升的 lambda 总数: **{total_lambdas}**")
    lines.append("")

    if fail_details:
        lines.append("## FAIL 列表")
        lines.append("")
        for fn, reason in fail_details:
            lines.append(f"- `{fn}`: {reason}")
        lines.append("")

    # 详情：有 lambda 的宿主
    if hosts_with_lambdas:
        lines.append("## 宿主函数 lambda 提升详情")
        lines.append("")
        lines.append("| 宿主函数 | 提升数 | 生成的 aux_lambdas |")
        lines.append("|---|---|---|")
        for tag, d in sorted(hosts_with_lambdas):
            names = ", ".join(f"`{n}`({pc})" for n, pc in
                              zip(d["lambda_names"], d["lambda_param_counts"]))
            lines.append(f"| `{tag}` | {d['n_lambdas']} | {names} |")
        lines.append("")
        lines.append("_(params count 包含 captures + lambda 自身参数)_")
        lines.append("")

    # 核验：与 scan_lambdas 预期值对照
    # lambdas.md 报告 26 in-scope lambdas (14 hosts)；其中 factorize.hh:108 的 lambda
    # 被错误归属 __factor_multivar（实际在 factorize 多变量 primary 中），
    # 因此 TRANSLATION_SCOPE sweep（factorize 展开 3 实例）预期：
    #   - __factor_multivar: 1 (wang.hh:2537)
    #   - factorize_lex:    1 (factorize.hh:108)
    #   - factorize_upoly:  1 (upoly.hh:1628 sort 比较器)
    #   - factorize_grlex:  0 (wrapper, 无 lambda)
    # 其他 12 个宿主按 scan_lambdas 计 23
    # 总计: 1 + 23 + 1 + 1 + 0 = 26
    EXPECTED_LAMBDAS = 26
    lines.append("## 与 `lambdas.md` 预期值核对")
    lines.append("")
    lines.append(f"- 本轮测得 lambda 总数: **{total_lambdas}**")
    lines.append(f"- `lambdas.md` in-scope 预期: **{EXPECTED_LAMBDAS}**")
    lines.append("")
    if total_lambdas == EXPECTED_LAMBDAS:
        lines.append("✅ 与预期值完全一致")
    else:
        lines.append(f"⚠️ 不一致，差异 {abs(total_lambdas - EXPECTED_LAMBDAS)}")
    lines.append("")

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    # stdout 总结
    print(f"\n=== {ok_count} OK / {fail_count} FAIL ===", file=sys.stderr)
    print(f"  hosts with lambdas: {len(hosts_with_lambdas)}", file=sys.stderr)
    print(f"  total lambdas:      {total_lambdas} (expected: {EXPECTED_LAMBDAS})", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    print(f"HIR₂ dumps: {DUMP_DIR}", file=sys.stderr)

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
