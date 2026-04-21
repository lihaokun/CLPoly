#!/usr/bin/env python3
"""scan_stl.py — 全量 STL 依赖扫描。

从缓存的 AST 中找所有 std::* 符号（referencedDecl / qualType），分类输出：

1. **容器类型**: std::vector / std::pair / std::map / std::set / std::tuple 等
2. **算法**: std::sort / std::find / std::move / ... (扩展 scan_iterators 结果)
3. **随机数**: std::mt19937 / std::uniform_int_distribution / std::random_device
4. **工具**: std::string / std::numeric_limits / std::initializer_list 等
5. **智能指针/引用包装**: std::unique_ptr / std::shared_ptr / std::ref（如有）
6. **异常/错误**: std::runtime_error / std::exception / std::logic_error（如有）

目的：为 clpoly_model.lean 补全 Lean 端 shim 提供完整清单。

输出：`docs/design/l1-translation-validation/survey/stl.md`
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "stl.md"


# === 分类字典 ===
CATEGORIES = {
    "容器": [
        "std::vector", "std::array", "std::pair", "std::tuple",
        "std::map", "std::unordered_map", "std::set", "std::unordered_set",
        "std::deque", "std::list", "std::stack", "std::queue",
        "std::__detail::_Node_iterator", "std::__cxx11::basic_string",
    ],
    "迭代器": [
        "std::back_insert_iterator", "std::move_iterator",
        "std::reverse_iterator", "std::__normal_iterator",
    ],
    "算法": [
        "std::sort", "std::stable_sort", "std::partial_sort",
        "std::find", "std::find_if", "std::find_if_not",
        "std::remove", "std::remove_if", "std::unique", "std::reverse",
        "std::copy", "std::copy_if", "std::move",
        "std::fill", "std::fill_n", "std::generate",
        "std::transform", "std::accumulate",
        "std::any_of", "std::all_of", "std::none_of", "std::for_each",
        "std::count", "std::count_if",
        "std::min", "std::max", "std::minmax",
        "std::min_element", "std::max_element",
        "std::lower_bound", "std::upper_bound", "std::binary_search",
        "std::iota", "std::partition", "std::nth_element",
        "std::adjacent_find", "std::mismatch", "std::equal", "std::includes",
        "std::set_union", "std::set_intersection", "std::set_difference",
        "std::merge", "std::rotate",
        "std::distance", "std::advance", "std::next", "std::prev",
        "std::swap",
    ],
    "数学": [
        "std::log", "std::log2", "std::log10", "std::exp",
        "std::sqrt", "std::cbrt", "std::pow",
        "std::ceil", "std::floor", "std::round", "std::trunc",
        "std::abs", "std::fabs", "std::sin", "std::cos", "std::tan",
        "std::atan2", "std::atan", "std::asin", "std::acos",
        "std::hypot", "std::fmod", "std::fma",
        "std::numeric_limits",
    ],
    "随机": [
        "std::mt19937", "std::mt19937_64", "std::minstd_rand",
        "std::random_device",
        "std::uniform_int_distribution", "std::uniform_real_distribution",
        "std::bernoulli_distribution", "std::normal_distribution",
        "std::discrete_distribution",
    ],
    "工具": [
        "std::string", "std::basic_string", "std::string_view",
        "std::to_string", "std::stoi", "std::stoll", "std::stoul", "std::stoull",
        "std::forward", "std::move", "std::make_pair", "std::make_tuple",
        "std::get", "std::initializer_list",
        "std::size_t", "std::ptrdiff_t", "std::nullptr_t",
        "std::begin", "std::end", "std::size", "std::empty",
    ],
    "智能指针": [
        "std::unique_ptr", "std::shared_ptr", "std::weak_ptr",
        "std::make_unique", "std::make_shared",
        "std::ref", "std::cref", "std::reference_wrapper",
    ],
    "异常": [
        "std::exception", "std::runtime_error", "std::logic_error",
        "std::out_of_range", "std::invalid_argument", "std::range_error",
        "std::overflow_error", "std::underflow_error",
        "std::bad_alloc", "std::bad_cast",
    ],
    "IO/流": [
        "std::cout", "std::cerr", "std::cin", "std::endl",
        "std::ostream", "std::istream", "std::stringstream",
        "std::ofstream", "std::ifstream", "std::fstream",
    ],
    "时间/并发": [
        "std::chrono", "std::thread", "std::mutex", "std::atomic",
        "std::lock_guard", "std::unique_lock", "std::condition_variable",
    ],
}


def classify(symbol: str) -> str | None:
    """根据符号前缀判断类别，返回类别名或 None（未分类）。"""
    for cat, keywords in CATEGORIES.items():
        for kw in keywords:
            if symbol.startswith(kw):
                return cat
    return None


def walk_ast(node, func_name: str, type_counts: Counter, type_hosts: defaultdict,
             ref_counts: Counter, ref_hosts: defaultdict):
    """遍历 AST，收集 std::* 的 qualType 和 referencedDecl.name。"""
    if isinstance(node, dict):
        # qualType
        t = node.get("type")
        if isinstance(t, dict):
            qt = t.get("qualType", "")
            # 提取所有 std::... 前缀
            for m in re.finditer(r"std::[\w:<>,\s]*\w", qt):
                sym = m.group(0)
                # 截短到第一个 `<` 之前
                clean = re.sub(r"<.*", "", sym).strip()
                type_counts[clean] += 1
                type_hosts[clean].add(func_name)

        # referencedDecl（函数调用等）
        rd = node.get("referencedDecl")
        if isinstance(rd, dict):
            name = rd.get("name", "")
            if name.startswith("std::") or (
                name and node.get("referencedDecl", {}).get("qualType", "").startswith("std::")
            ):
                ref_counts[name] += 1
                ref_hosts[name].add(func_name)
            else:
                # 看 mangledName 是否以 _ZSt 开头（std 符号）
                mangled = rd.get("mangledName", "")
                if mangled.startswith("_ZSt") and name:
                    # std::name
                    key = "std::" + name
                    ref_counts[key] += 1
                    ref_hosts[key].add(func_name)

        for c in node.get("inner", []):
            walk_ast(c, func_name, type_counts, type_hosts, ref_counts, ref_hosts)
    elif isinstance(node, list):
        for x in node:
            walk_ast(x, func_name, type_counts, type_hosts, ref_counts, ref_hosts)


def main():
    type_counts: Counter[str] = Counter()
    type_hosts: dict[str, set[str]] = defaultdict(set)
    ref_counts: Counter[str] = Counter()
    ref_hosts: dict[str, set[str]] = defaultdict(set)

    for f in sorted(AST_CACHE_DIR.glob("*.json")):
        fn = f.stem
        with open(f) as fp:
            ast = json.load(fp)
        walk_ast(ast, fn, type_counts, type_hosts, ref_counts, ref_hosts)

    # 合并 type + ref 到一个 unified view
    all_symbols = set(type_counts.keys()) | set(ref_counts.keys())

    # 按类别分组
    by_category: dict[str, list] = defaultdict(list)
    unclassified = []
    for sym in all_symbols:
        cat = classify(sym)
        total = type_counts[sym] + ref_counts[sym]
        hosts = type_hosts[sym] | ref_hosts[sym]
        entry = {
            "symbol": sym,
            "type_count": type_counts[sym],
            "ref_count": ref_counts[sym],
            "total": total,
            "hosts": hosts,
        }
        if cat:
            by_category[cat].append(entry)
        else:
            unclassified.append(entry)

    # 输出
    lines = []
    lines.append("# STL 依赖调研")
    lines.append("")
    lines.append("为 `clpoly_model.lean` 端 Lean shim 设计提供完整清单。")
    lines.append("")
    lines.append(f"总计 **{len(all_symbols)}** 个 std::* 符号。")
    lines.append("")

    for cat in ["容器", "迭代器", "算法", "数学", "随机", "工具",
                "智能指针", "异常", "IO/流", "时间/并发"]:
        entries = sorted(by_category.get(cat, []), key=lambda e: -e["total"])
        if not entries:
            continue
        lines.append(f"## {cat}")
        lines.append("")
        lines.append("| 符号 | 类型出现 | 调用出现 | 合计 | 宿主函数 |")
        lines.append("|---|---|---|---|---|")
        for e in entries:
            hosts = sorted(e["hosts"])
            if len(hosts) <= 3:
                host_str = ", ".join(f"`{h}`" for h in hosts)
            else:
                host_str = ", ".join(f"`{h}`" for h in hosts[:3]) + f" … (共 {len(hosts)})"
            lines.append(f"| `{e['symbol']}` | {e['type_count']} | {e['ref_count']} | {e['total']} | {host_str} |")
        lines.append("")

    if unclassified:
        lines.append("## 未分类")
        lines.append("")
        lines.append("| 符号 | 类型出现 | 调用出现 | 合计 | 宿主 |")
        lines.append("|---|---|---|---|---|")
        for e in sorted(unclassified, key=lambda x: -x["total"])[:50]:
            hosts = sorted(e["hosts"])
            host_str = ", ".join(f"`{h}`" for h in hosts[:3])
            if len(hosts) > 3:
                host_str += f" … ({len(hosts)})"
            lines.append(f"| `{e['symbol']}` | {e['type_count']} | {e['ref_count']} | {e['total']} | {host_str} |")
        lines.append("")

    # Lean shim 建议
    lines.append("## Lean shim 设计建议")
    lines.append("")
    lines.append("对每类 STL 符号，推荐的 Lean 端实现方式：")
    lines.append("")
    lines.append("| 类别 | Lean 对应 | 备注 |")
    lines.append("|---|---|---|")
    lines.append("| `std::vector<T>` | `Array T` | 直接等价；`.push_back` → `.push`；`.erase(it, end)` → `.take n` |")
    lines.append("| `std::pair<A,B>` | `A × B` | 原生 product type |")
    lines.append("| `std::tuple<A,B,C>` | `A × B × C` / 结构体 | 三元 product 或自定义 |")
    lines.append("| `std::map<K,V>` | `StdMap K V`（自定义）| 按 key 排序的 Array；CLPoly 用它做线性探针 |")
    lines.append("| `std::sort` + comparator | `Array.qsortWith` | lean 4.x 有 `Array.qsort` |")
    lines.append("| `std::iota(b, e, v)` | `Array.range` | 生成 [v, v+1, ..., ) |")
    lines.append("| `std::move` | `id` | Lean 值语义，move 是 no-op |")
    lines.append("| `std::swap(a, b)` | let-rebind | Lean 无 in-place swap，用 `(b, a) := (a, b)` |")
    lines.append("| `std::max`/`std::min` | `max` / `min` | Mathlib/stdlib |")
    lines.append("| `std::mt19937` + `uniform_int_distribution` | `Rng` 结构体 + `Rng.next` | 公理化（EDF 已用）|")
    lines.append("| `std::random_device` | `axiom Rng.seed : Nat` | 仅作 seed 源 |")
    lines.append("| `std::log` | `Float.log` / `Nat.log` | 依上下文选具体数值类型 |")
    lines.append("| `std::ceil` / `std::floor` | `Nat.ceilDiv` / 截断 | 注意 double 语义 |")
    lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    for cat, entries in by_category.items():
        print(f"  {cat}: {len(entries)} symbols, {sum(e['total'] for e in entries)} total occurrences")
    print(f"  未分类: {len(unclassified)}")


if __name__ == "__main__":
    main()
