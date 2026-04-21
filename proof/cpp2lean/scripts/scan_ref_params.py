#!/usr/bin/env python3
"""scan_ref_params.py — 分析每个函数的参数模式。

对每个 TRANSLATION_SCOPE 函数：
- 列出每个参数的类型 + 传递方式（const ref / non-const ref / 指针 / 值）
- 标识输出参数（non-const ref 且被函数修改）
- 输出汇总：有多少函数用值传递、ref、pointer

为 HIR ref_elim Pass 设计提供输入。

输出：`docs/design/l1-translation-validation/survey/ref_params.md`
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
from collections import defaultdict, Counter

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "ref_params.md"

sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE, TRANSLATION_SCOPE_OUTPUT_PARAMS


def classify_param_type(qt: str) -> str:
    """根据 qualType 判断传递方式。"""
    if not qt:
        return "unknown"
    # 注意顺序：const T& 包含 const，必须先检查
    if "const " in qt and qt.endswith("&"):
        return "const-ref"
    elif qt.endswith("&&"):
        return "rvalue-ref"
    elif qt.endswith("&"):
        return "ref"
    elif qt.endswith("*"):
        return "pointer"
    else:
        return "value"


def extract_params(ast_func: dict) -> list[dict]:
    params = []
    for c in ast_func.get("inner", []):
        if isinstance(c, dict) and c.get("kind") == "ParmVarDecl":
            name = c.get("name")
            qt = c.get("type", {}).get("qualType", "")
            kind = classify_param_type(qt)
            params.append({"name": name, "type": qt, "kind": kind})
    return params


def main():
    per_func = {}
    for f in sorted(AST_CACHE_DIR.glob("*.json")):
        fn = f.stem
        with open(f) as fp:
            ast = json.load(fp)
        # ast 顶层应是 FunctionDecl
        if ast.get("kind") != "FunctionDecl":
            continue
        params = extract_params(ast)
        per_func[fn] = params

    # 汇总统计
    kind_counts = Counter()
    for params in per_func.values():
        for p in params:
            kind_counts[p["kind"]] += 1

    # 每函数至少 1 个 non-const ref 的函数数
    has_out_param = [fn for fn, ps in per_func.items() if any(p["kind"] == "ref" for p in ps)]
    has_ptr = [fn for fn, ps in per_func.items() if any(p["kind"] == "pointer" for p in ps)]
    has_rvalue = [fn for fn, ps in per_func.items() if any(p["kind"] == "rvalue-ref" for p in ps)]

    lines = []
    lines.append("# 参数传递模式扫描")
    lines.append("")
    lines.append("为 HIR `ref_elim` Pass 设计提供输入。")
    lines.append("")
    lines.append("## 全局统计")
    lines.append("")
    total_params = sum(kind_counts.values())
    lines.append(f"- 总参数数: **{total_params}**")
    for k in ["value", "const-ref", "ref", "pointer", "rvalue-ref", "unknown"]:
        c = kind_counts.get(k, 0)
        pct = 100 * c / total_params if total_params else 0
        lines.append(f"- `{k}`: {c} ({pct:.0f}%)")
    lines.append("")
    lines.append(f"- 含 non-const ref 参数的函数（疑似输出参数）: **{len(has_out_param)}**")
    lines.append(f"- 含指针参数的函数: **{len(has_ptr)}**")
    lines.append(f"- 含 rvalue-ref 参数的函数: **{len(has_rvalue)}**")
    lines.append("")

    # 对比 TRANSLATION_SCOPE_OUTPUT_PARAMS（已配置的 output param 清单）
    lines.append("## 输出参数对比（AST 扫描 vs TRANSLATION_SCOPE_OUTPUT_PARAMS）")
    lines.append("")
    lines.append("| 函数 | AST 扫出的 non-const ref 位置 | 已配置的 OUTPUT_PARAMS | 是否一致 |")
    lines.append("|---|---|---|---|")
    for fn in sorted(has_out_param):
        params = per_func[fn]
        ast_ref_indices = [i for i, p in enumerate(params) if p["kind"] == "ref"]
        configured = TRANSLATION_SCOPE_OUTPUT_PARAMS.get(fn, [])
        match = "✅" if sorted(ast_ref_indices) == sorted(configured) else "❌"
        lines.append(f"| `{fn}` | {ast_ref_indices} | {configured} | {match} |")
    lines.append("")

    # 每函数的参数详情
    lines.append("## 按函数参数详情")
    lines.append("")
    for fn in sorted(per_func.keys()):
        params = per_func[fn]
        if not params:
            lines.append(f"### `{fn}` (0 params)")
            lines.append("")
            continue
        lines.append(f"### `{fn}` ({len(params)} params)")
        for i, p in enumerate(params):
            tag = {
                "value": "VAL",
                "const-ref": "CREF",
                "ref": "**REF (output?)**",
                "pointer": "PTR",
                "rvalue-ref": "RVAL",
            }.get(p["kind"], p["kind"])
            lines.append(f"- [{i}] `{p.get('name') or '<unnamed>'}` : `{(p['type'] or '')[:80]}` — {tag}")
        lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Total params: {total_params}")
    for k in ["value", "const-ref", "ref", "pointer", "rvalue-ref", "unknown"]:
        c = kind_counts.get(k, 0)
        print(f"  {k}: {c}")
    print(f"Functions with non-const ref: {len(has_out_param)}")
    print(f"  configured OUTPUT_PARAMS match: {sum(1 for fn in has_out_param if sorted([i for i,p in enumerate(per_func[fn]) if p['kind']=='ref']) == sorted(TRANSLATION_SCOPE_OUTPUT_PARAMS.get(fn, [])))}/{len(has_out_param)}")


if __name__ == "__main__":
    main()
