#!/usr/bin/env python3
"""scan_iterators.py — 扫描 CLPoly 源码中的迭代器模式。

为 HIR iter_recognize Pass 设计输入。六种典型模式：

1. Range-for：`for (auto& x : v)` — 从 AST 的 CXXForRangeStmt 统计
2. Classic iterator loop：`for (auto it = v.begin(); it != v.end(); ++it)`
3. 双指针 compact：`auto it, out = v.begin(); ...; v.erase(out, v.end());`
4. std::find* / std::remove_if / std::count_if 等 STL 算法调用
5. StdMap 遍历：`for (const auto& [k, v] : m)` — 结构化绑定在 range-for
6. begin/end/++/* 的其他零散使用

输出 `docs/design/l1-translation-validation/survey/iterators.md`。
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
PROJECT_ROOT = HERE.parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = PROJECT_ROOT / "docs" / "design" / "l1-translation-validation" / "survey" / "iterators.md"
CLPOLY_DIR = PROJECT_ROOT / "clpoly"

sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE


# === 1. AST 侧：Range-for 统计 ===

def find_range_for_hosts(cache_dir: Path) -> dict[str, int]:
    """扫每个 AST，统计 CXXForRangeStmt + 是否有 DecompositionDecl。"""
    hosts = defaultdict(lambda: {"range_for": 0, "range_for_with_binding": 0})
    for f in sorted(cache_dir.glob("*.json")):
        func_name = f.stem
        with open(f) as fp:
            ast = json.load(fp)
        def walk(n, inside_range_for=False):
            if isinstance(n, dict):
                k = n.get("kind")
                new_inside = inside_range_for
                if k == "CXXForRangeStmt":
                    hosts[func_name]["range_for"] += 1
                    new_inside = True
                if k in ("DecompositionDecl", "BindingDecl") and inside_range_for:
                    # 仅计 range-for 的第一个 DecompositionDecl
                    pass  # 交给下面的检测
                for c in n.get("inner", []):
                    walk(c, new_inside)
            elif isinstance(n, list):
                for x in n:
                    walk(x, inside_range_for)
        # 另一轮专门检测 "range-for + DecompositionDecl" 组合
        def walk2(n):
            if isinstance(n, dict):
                if n.get("kind") == "CXXForRangeStmt":
                    # 检查 inner 是否含 DecompositionDecl
                    for c in n.get("inner", []):
                        if isinstance(c, dict):
                            # 递归找第一层
                            def has_decomp(node):
                                if not isinstance(node, dict):
                                    return False
                                if node.get("kind") == "DecompositionDecl":
                                    return True
                                for cc in node.get("inner", []):
                                    if has_decomp(cc):
                                        return True
                                return False
                            if has_decomp(c):
                                hosts[func_name]["range_for_with_binding"] += 1
                                break
                for c in n.get("inner", []):
                    walk2(c)
            elif isinstance(n, list):
                for x in n:
                    walk2(x)
        walk(ast)
        walk2(ast)
    return hosts


# === 2. 源码侧：Classic iterator / double-pointer / STL algorithms ===

# 匹配 `for (auto it = X.begin(); ... ++it)` 或类似
CLASSIC_ITER_RE = re.compile(
    r"for\s*\(\s*auto\s+\w+\s*=\s*[^;]+\.begin\s*\(\s*\)",
    re.MULTILINE
)

# 迭代器初始化：`auto <name> = <expr>.begin()` 或紧随 `auto <name2> = <name1>`
# 用来定位潜在的双指针 compact 开始点
ITER_INIT_RE = re.compile(
    r"auto\s+(\w+)\s*=\s*([^;]+?)\.begin\s*\(\s*\)",
    re.MULTILINE
)

# 双指针在一行：`auto it = X.begin(), out = it;` / `auto it = X.begin(), out = X.begin();`
COMPACT_IT_OUT_ONELINE_RE = re.compile(
    r"auto\s+\w+\s*=\s*[^;,]+?\.begin\s*\(\s*\)\s*,\s*\w+\s*=",
    re.MULTILINE
)

# Compact erase（松散）：任何 `.erase(X, Y.end())` 形式
COMPACT_ERASE_RE = re.compile(
    r"(\w+(?:\.\w+\s*\(\s*\))*?)\.erase\s*\(\s*\w+\s*,\s*\1\.end\s*\(\s*\)\s*\)",
    re.MULTILINE
)

# 任一 .erase(...) 调用（含位置删除、范围删除）
ANY_ERASE_RE = re.compile(
    r"\w+(?:\.\w+\s*\(\s*\))*?\.erase\s*\(",
    re.MULTILINE
)

# STL 算法
STL_ALG_RE = re.compile(
    r"std::(sort|stable_sort|partial_sort|find|find_if|find_if_not|remove_if|"
    r"unique|reverse|rotate|count|count_if|any_of|all_of|none_of|for_each|"
    r"transform|accumulate|copy|copy_if|move|swap|lower_bound|upper_bound|"
    r"binary_search|iota|fill|fill_n|generate|partition|nth_element|"
    r"min_element|max_element|max|min|minmax|adjacent_find|mismatch|equal|"
    r"includes|set_union|set_intersection|set_difference|merge|distance|advance|"
    r"next|prev)\b",
    re.MULTILINE
)

# 结构化绑定 range-for
STRUCT_BIND_RANGE_RE = re.compile(
    r"for\s*\(\s*(?:const\s+)?auto\s*&?\s*\[\s*[\w,\s]+\s*\]\s*:",
    re.MULTILINE
)


def _find_host_function(file_path: Path, line: int, known_hosts: set[str]) -> str | None:
    """从 class_map.TRANSLATION_SCOPE + factorize 中找 line 对应的宿主函数。"""
    lines = file_path.read_text().splitlines()
    best = None
    best_line = -1
    for i, src_line in enumerate(lines):
        if i + 1 > line:
            break
        for fn in known_hosts:
            idx = src_line.find(fn + "(")
            if idx < 0:
                continue
            if idx > 0 and (src_line[idx - 1].isalnum() or src_line[idx - 1] == "_"):
                continue
            # 粗略判定是否是定义：行首全空白或已有返回类型
            pre = src_line[:idx].strip()
            prev = lines[i - 1].strip() if i > 0 else ""
            def likely_def(s):
                if not s: return False
                if s.startswith("inline") or s.startswith("template"): return True
                if s.endswith(">") or s.endswith("&"): return True
                return s[-1].isalpha() or s[-1] == "_"
            is_def = (not pre and likely_def(prev)) or likely_def(pre)
            if is_def and (i + 1) > best_line:
                best = fn
                best_line = i + 1
    return best


def scan_source_patterns():
    files = sorted(CLPOLY_DIR.glob("polynomial_factorize*.hh"))
    known_hosts = set(TRANSLATION_SCOPE) | {"factorize", "squarefreefactorize"}

    classic_iter = defaultdict(list)  # func → [(file, line, text)]
    compact_it_out = defaultdict(list)
    compact_erase = defaultdict(list)
    any_erase = defaultdict(list)
    iter_init = defaultdict(list)
    stl_algs = Counter()  # alg_name → total count
    stl_alg_hosts = defaultdict(lambda: Counter())  # alg_name → {host: count}
    struct_bind_range = defaultdict(list)

    for f in files:
        content = f.read_text()
        # 1) Classic iterator loops
        for m in CLASSIC_ITER_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                classic_iter[host].append((f.name, line, m.group(0)[:80]))
        # 2a) Iterator 初始化（auto X = Y.begin()）
        for m in ITER_INIT_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                iter_init[host].append((f.name, line, m.group(0)[:80]))
        # 2b) 双指针 one-liner
        for m in COMPACT_IT_OUT_ONELINE_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                compact_it_out[host].append((f.name, line, m.group(0)[:80]))
        # 3a) Compact erase (严格：X.erase(Y, X.end()))
        for m in COMPACT_ERASE_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                compact_erase[host].append((f.name, line, m.group(0)[:80]))
        # 3b) 所有 .erase(...) 调用
        for m in ANY_ERASE_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                any_erase[host].append((f.name, line, m.group(0)[:80]))
        # 4) STL algorithms
        for m in STL_ALG_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            alg = m.group(1)
            stl_algs[alg] += 1
            if host:
                stl_alg_hosts[alg][host] += 1
        # 5) Structured-binding range-for
        for m in STRUCT_BIND_RANGE_RE.finditer(content):
            line = content[: m.start()].count("\n") + 1
            host = _find_host_function(f, line, known_hosts)
            if host:
                struct_bind_range[host].append((f.name, line, m.group(0)[:80]))

    return {
        "classic_iter": dict(classic_iter),
        "iter_init": dict(iter_init),
        "compact_it_out": dict(compact_it_out),
        "compact_erase": dict(compact_erase),
        "any_erase": dict(any_erase),
        "stl_algs": stl_algs,
        "stl_alg_hosts": dict(stl_alg_hosts),
        "struct_bind_range": dict(struct_bind_range),
    }


def main():
    range_for_hosts = find_range_for_hosts(AST_CACHE_DIR)
    src_patterns = scan_source_patterns()

    total_range_for = sum(h["range_for"] for h in range_for_hosts.values())
    total_range_for_binding = sum(h["range_for_with_binding"] for h in range_for_hosts.values())
    total_classic = sum(len(v) for v in src_patterns["classic_iter"].values())
    total_iter_init = sum(len(v) for v in src_patterns["iter_init"].values())
    total_compact_it_out = sum(len(v) for v in src_patterns["compact_it_out"].values())
    total_compact_erase = sum(len(v) for v in src_patterns["compact_erase"].values())
    total_any_erase = sum(len(v) for v in src_patterns["any_erase"].values())
    total_struct_bind = sum(len(v) for v in src_patterns["struct_bind_range"].values())

    lines = []
    lines.append("# 迭代器模式扫描")
    lines.append("")
    lines.append("为 HIR `iter_recognize` Pass 设计输入。")
    lines.append("")
    lines.append("## 全局统计")
    lines.append("")
    lines.append(f"| 模式 | 数量 | 备注 |")
    lines.append(f"|---|---|---|")
    lines.append(f"| Range-for `for (auto& x : v)` | **{total_range_for}** | AST `CXXForRangeStmt` 统计 |")
    lines.append(f"| Range-for 含结构化绑定 `for (auto& [k,v] : m)` | **{total_range_for_binding}** | 需要 HIR 特殊处理 |")
    lines.append(f"| Classic iterator loop `for (auto it = X.begin();...;++it)` | **{total_classic}** | 源码正则 |")
    lines.append(f"| Iterator 初始化 `auto X = Y.begin()` | **{total_iter_init}** | 源码正则（潜在双指针 compact 信号）|")
    lines.append(f"| One-liner 双指针 (`auto it = X.begin(), out = ...`) | **{total_compact_it_out}** | 源码正则 |")
    lines.append(f"| Compact erase (`X.erase(out, X.end())`) | **{total_compact_erase}** | 源码正则（严格匹配）|")
    lines.append(f"| 任意 `.erase(...)` 调用 | **{total_any_erase}** | 源码正则（含位置删除、range 删除）|")
    lines.append(f"| Structured-binding range-for | **{total_struct_bind}** | 源码正则 |")
    lines.append("")

    # STL 算法调用
    lines.append("## STL 算法调用直方图")
    lines.append("")
    lines.append("| 算法 | 总次数 | 首次出现函数 |")
    lines.append("|---|---|---|")
    for alg, cnt in src_patterns["stl_algs"].most_common():
        hosts = src_patterns["stl_alg_hosts"].get(alg, Counter())
        host_list = ", ".join(f"`{h}`({c})" for h, c in hosts.most_common(3))
        if len(hosts) > 3:
            host_list += f", ... (共 {len(hosts)})"
        lines.append(f"| `std::{alg}` | {cnt} | {host_list or '(无 in-scope 宿主)'} |")
    lines.append("")

    # Range-for 分布
    lines.append("## Range-for 分布（按宿主函数）")
    lines.append("")
    sorted_hosts = sorted(range_for_hosts.items(), key=lambda kv: -kv[1]["range_for"])
    lines.append("| 宿主 | range-for | 其中含结构化绑定 |")
    lines.append("|---|---|---|")
    for fn, h in sorted_hosts:
        if h["range_for"] == 0:
            continue
        lines.append(f"| `{fn}` | {h['range_for']} | {h['range_for_with_binding']} |")
    lines.append("")

    # Classic iterator 详情
    if src_patterns["classic_iter"]:
        lines.append("## Classic iterator loop 详情")
        lines.append("")
        for fn, items in sorted(src_patterns["classic_iter"].items()):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text}`")
            lines.append("")

    # Double-pointer compact 详情
    if src_patterns["compact_it_out"]:
        lines.append("## 双指针 compact 详情")
        lines.append("")
        for fn, items in sorted(src_patterns["compact_it_out"].items()):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text}`")
            lines.append("")

    # Compact erase 详情
    if src_patterns["compact_erase"]:
        lines.append("## Compact erase 详情（X.erase(out, X.end()) 严格匹配）")
        lines.append("")
        for fn, items in sorted(src_patterns["compact_erase"].items()):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text}`")
            lines.append("")

    # Iterator 初始化详情（所有 auto X = Y.begin()）
    if src_patterns["iter_init"]:
        lines.append("## Iterator 初始化详情（auto X = Y.begin()，潜在 compact 信号）")
        lines.append("")
        for fn, items in sorted(src_patterns["iter_init"].items(), key=lambda kv: -len(kv[1])):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text}`")
            lines.append("")

    # 任意 erase 详情
    if src_patterns["any_erase"]:
        lines.append("## `.erase(...)` 全部调用")
        lines.append("")
        for fn, items in sorted(src_patterns["any_erase"].items(), key=lambda kv: -len(kv[1])):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text}`")
            lines.append("")

    # 结构化绑定 range-for 详情
    if src_patterns["struct_bind_range"]:
        lines.append("## Structured-binding range-for 详情")
        lines.append("")
        for fn, items in sorted(src_patterns["struct_bind_range"].items()):
            lines.append(f"### `{fn}` ({len(items)})")
            for file_name, line, text in items:
                lines.append(f"- {file_name}:{line} — `{text[:100]}`")
            lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Range-for: {total_range_for} ({total_range_for_binding} with binding)")
    print(f"Classic iter: {total_classic}, iter-init: {total_iter_init}, compact-it/out one-liner: {total_compact_it_out}")
    print(f"Compact-erase strict: {total_compact_erase}, any-erase: {total_any_erase}")
    print(f"Struct-bind range-for: {total_struct_bind}")
    print(f"STL algs: {sum(src_patterns['stl_algs'].values())} across {len(src_patterns['stl_algs'])} unique")


if __name__ == "__main__":
    main()
