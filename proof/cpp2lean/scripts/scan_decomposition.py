#!/usr/bin/env python3
"""scan_decomposition.py — 扫描结构化绑定（DecompositionDecl）详情。

从 AST 找每个 DecompositionDecl：
- 所在函数、行号
- 绑定名列表（如 `[k, v]`、`[gk, mk]`）
- 被解构的容器类型（pair / tuple / map entry / custom struct）
- 是否出现在 range-for 中
- 是 `auto&` / `auto` / `const auto&`

输出：`docs/design/l1-translation-validation/survey/decompositions.md`
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
from collections import defaultdict, Counter

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "decompositions.md"


def walk(node, func_name, parent_stack, out):
    if isinstance(node, dict):
        if node.get("kind") == "DecompositionDecl":
            # 分析此节点
            # 绑定名：inner 里的 BindingDecl
            bindings = []
            for c in node.get("inner", []):
                if isinstance(c, dict) and c.get("kind") == "BindingDecl":
                    bindings.append({
                        "name": c.get("name"),
                        "type": c.get("type", {}).get("qualType"),
                    })

            # 自身的 type（绑定左边的 `auto&` 类型）
            self_type = node.get("type", {}).get("qualType")

            # 被解构的值：inner 中的 init 表达式（通常第一个非 BindingDecl）
            init_kind = None
            init_type = None
            for c in node.get("inner", []):
                if isinstance(c, dict) and c.get("kind") != "BindingDecl":
                    init_kind = c.get("kind")
                    init_type = c.get("type", {}).get("qualType")
                    break

            # 是否在 range-for 中（parent_stack 往上有 CXXForRangeStmt）
            in_range_for = any(p.get("kind") == "CXXForRangeStmt"
                               for p in parent_stack)

            out.append({
                "func": func_name,
                "bindings": bindings,
                "self_type": self_type,
                "init_kind": init_kind,
                "init_type": init_type,
                "in_range_for": in_range_for,
            })
        for c in node.get("inner", []):
            walk(c, func_name, parent_stack + [node], out)
    elif isinstance(node, list):
        for x in node:
            walk(x, func_name, parent_stack, out)


def main():
    all_decomps = []
    for f in sorted(AST_CACHE_DIR.glob("*.json")):
        fn = f.stem
        with open(f) as fp:
            ast = json.load(fp)
        walk(ast, fn, [], all_decomps)

    lines = []
    lines.append("# 结构化绑定（DecompositionDecl）扫描")
    lines.append("")
    lines.append(f"总计 **{len(all_decomps)}** 处，分布于 **{len(set(d['func'] for d in all_decomps))}** 个函数。")
    lines.append("")

    # 分类统计
    in_for = sum(1 for d in all_decomps if d["in_range_for"])
    outside = len(all_decomps) - in_for
    lines.append(f"- 在 range-for 中: **{in_for}**")
    lines.append(f"- 独立 DeclStmt: **{outside}**")
    lines.append("")

    # 按 container 类型分组
    type_groups = Counter()
    for d in all_decomps:
        # 提取 init_type 的顶层模板名
        it = d["init_type"] or ""
        if "pair" in it:
            type_groups["pair"] += 1
        elif "tuple" in it:
            type_groups["tuple"] += 1
        elif "map" in it:
            type_groups["map entry"] += 1
        elif "basic_monomial" in it:
            type_groups["monomial (var, deg)"] += 1
        else:
            type_groups[f"其他: {it[:50]}"] += 1
    lines.append("## 绑定 container 类型")
    lines.append("")
    for t, c in type_groups.most_common():
        lines.append(f"- {t}: {c}")
    lines.append("")

    # 绑定名模式
    name_patterns = Counter()
    for d in all_decomps:
        names = tuple(b.get("name") or "?" for b in d["bindings"])
        name_patterns[names] += 1
    lines.append("## 常见绑定名模式（前 20）")
    lines.append("")
    for names, c in name_patterns.most_common(20):
        lines.append(f"- `[{', '.join(names)}]`: {c}")
    lines.append("")

    # 按函数详情
    by_func = defaultdict(list)
    for d in all_decomps:
        by_func[d["func"]].append(d)
    lines.append("## 按函数详情")
    lines.append("")
    for fn in sorted(by_func.keys()):
        decomps = by_func[fn]
        lines.append(f"### `{fn}` ({len(decomps)})")
        for i, d in enumerate(decomps, 1):
            names = ", ".join(b.get("name") or "?" for b in d["bindings"])
            in_for = " [range-for]" if d["in_range_for"] else ""
            lines.append(f"- `[{names}]`{in_for} — container: `{(d['init_type'] or '')[:80]}`")
        lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Total decompositions: {len(all_decomps)} ({in_for} in range-for, {outside} standalone)")


if __name__ == "__main__":
    main()
