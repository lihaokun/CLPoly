#!/usr/bin/env python3
"""
drill_dependent.py — 深挖 UnresolvedLookupExpr / CXXDependentScopeMemberExpr /
TemplateTypeParmDecl / FunctionTemplateDecl 出现在哪些函数、什么上下文。

用于验证 AST 是否完全单态化的假设。
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"

TARGETS = {
    "UnresolvedLookupExpr",
    "CXXDependentScopeMemberExpr",
    "TemplateTypeParmDecl",
    "FunctionTemplateDecl",
    "DependentScopeDeclRefExpr",
}


def find_nodes(node, target_kind, path=None, out=None):
    if path is None:
        path = []
    if out is None:
        out = []
    if isinstance(node, dict):
        k = node.get("kind")
        if k == target_kind:
            out.append((path[:], node))
        new_path = path + [k] if k else path
        for c in node.get("inner", []):
            find_nodes(c, target_kind, new_path, out)
    elif isinstance(node, list):
        for x in node:
            find_nodes(x, target_kind, path, out)
    return out


def main():
    per_kind_per_func = defaultdict(lambda: defaultdict(int))
    per_kind_samples = defaultdict(list)  # kind → list of (func_name, path, node)

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        for target in TARGETS:
            matches = find_nodes(ast, target)
            if matches:
                per_kind_per_func[target][func_name] = len(matches)
                # 存前 2 个样例
                for path, node in matches[:2]:
                    per_kind_samples[target].append((func_name, path, node))

    print("# Dependent/unresolved AST nodes drill-down")
    print()
    for target in sorted(TARGETS):
        funcs = per_kind_per_func[target]
        total = sum(funcs.values())
        print(f"## `{target}` (total {total} across {len(funcs)} functions)")
        print()
        if not funcs:
            print("— none —")
            print()
            continue
        # 按出现次数排序
        for fn, cnt in sorted(funcs.items(), key=lambda x: -x[1]):
            print(f"- `{fn}`: {cnt}")
        print()
        # 展示 2 个样例
        samples = per_kind_samples[target][:3]
        for func_name, path, node in samples:
            print(f"### sample from `{func_name}`")
            print(f"- parent path: `{ ' > '.join(path[-5:]) }`")
            slim = {k: v for k, v in node.items() if k != "inner"}
            print("```json")
            print(json.dumps(slim, ensure_ascii=False, indent=2)[:600])
            print("```")
            # 打印 inner 的 kinds（第一层）
            inners = node.get("inner", [])
            if inners:
                kinds = [i.get("kind", "?") for i in inners if isinstance(i, dict)]
                print(f"- inner kinds: {kinds}")
            print()


if __name__ == "__main__":
    main()
