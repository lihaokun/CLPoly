#!/usr/bin/env python3
"""
verify_lambda_hypothesis.py — 验证所有 dependent AST 节点是否都在 generic lambda 内部。

对每个 CXXDependentScopeMemberExpr / UnresolvedLookupExpr / TemplateTypeParmDecl
/ FunctionTemplateDecl 节点：
- 记录它的最近 LambdaExpr 祖先（有/无）
- 若无，打印其完整祖先链（说明它来自非 lambda 的模板上下文）
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"

TARGETS = {
    "CXXDependentScopeMemberExpr",
    "UnresolvedLookupExpr",
    "TemplateTypeParmDecl",
    "FunctionTemplateDecl",
}


def walk(node, target, ancestors, out):
    if isinstance(node, dict):
        kind = node.get("kind")
        if kind in target:
            # 检查祖先链中是否有 LambdaExpr
            lambda_ancestor_idx = -1
            for i, a in enumerate(ancestors):
                if a == "LambdaExpr":
                    lambda_ancestor_idx = i
            out.append({
                "kind": kind,
                "ancestors": ancestors[:],
                "has_lambda_ancestor": lambda_ancestor_idx >= 0,
                "lambda_depth": len(ancestors) - lambda_ancestor_idx if lambda_ancestor_idx >= 0 else None,
                "node": {k: v for k, v in node.items() if k not in ("inner", "range", "loc")}
            })
        new_anc = ancestors + [kind] if kind else ancestors
        for c in node.get("inner", []):
            walk(c, target, new_anc, out)
    elif isinstance(node, list):
        for x in node:
            walk(x, target, ancestors, out)


def main():
    per_func_summary = {}  # func → {"in_lambda": count, "outside_lambda": count, "samples_outside": [...]}

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        findings = []
        walk(ast, TARGETS, [], findings)
        if not findings:
            continue

        in_lambda = sum(1 for f in findings if f["has_lambda_ancestor"])
        outside = [f for f in findings if not f["has_lambda_ancestor"]]
        per_func_summary[func_name] = {
            "total": len(findings),
            "in_lambda": in_lambda,
            "outside_lambda": len(outside),
            "samples_outside": outside[:3],
        }

    # 报告
    print("# Lambda Hypothesis Verification")
    print()
    print("For each function with dependent AST nodes, report how many are inside")
    print("a LambdaExpr vs outside. If any are outside, Option 1 (replace `auto` with")
    print("explicit types in lambdas) will NOT fully resolve the issue.")
    print()

    grand_in = grand_out = 0
    for func, d in sorted(per_func_summary.items()):
        print(f"## `{func}`")
        print(f"- total: {d['total']}")
        print(f"- in lambda: **{d['in_lambda']}**")
        print(f"- outside lambda: **{d['outside_lambda']}**")
        grand_in += d["in_lambda"]
        grand_out += d["outside_lambda"]
        if d["samples_outside"]:
            print(f"- samples outside lambda:")
            for s in d["samples_outside"]:
                print(f"  - kind: `{s['kind']}`, ancestors: `{' > '.join(s['ancestors'][-6:])}`")
                print(f"    node: `{json.dumps(s['node'], ensure_ascii=False)[:200]}`")
        print()

    print(f"## Grand Total")
    print(f"- in lambda: **{grand_in}**")
    print(f"- outside lambda: **{grand_out}**")
    if grand_out == 0:
        print()
        print("✅ **Hypothesis verified**: ALL dependent AST nodes are inside generic lambdas.")
        print("Option 1 (use explicit types in lambda parameters) will fully resolve the issue.")
    else:
        print()
        print(f"❌ **Hypothesis rejected**: {grand_out} nodes are outside lambda context.")
        print("Option 1 alone is insufficient.")


if __name__ == "__main__":
    main()
