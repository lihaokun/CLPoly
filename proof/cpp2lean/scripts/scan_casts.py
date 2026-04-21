#!/usr/bin/env python3
"""scan_casts.py — 枚举所有 (source_type, target_type) 类型转换对。

从缓存的 AST 中遍历 ImplicitCastExpr / CStyleCastExpr / CXXStaticCastExpr /
CXXFunctionalCastExpr / CXXConstructExpr，提取:
- castKind（如 IntegralCast / NoOp / LValueToRValue / ...）
- source type（从 inner[0] 的 type.qualType 推断）
- target type（从自身 type.qualType）
- 宿主函数 + 出现次数

为 type-system.md §8 (CAST_TABLE) 提供完整枚举。

输出：`docs/design/l1-translation-validation/survey/casts.md`
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "casts.md"

CAST_KINDS = {
    "ImplicitCastExpr",
    "CStyleCastExpr",
    "CXXStaticCastExpr",
    "CXXFunctionalCastExpr",
    # CXXConstructExpr 在某些场景也是"转换"，但多数是对象构造；暂不归入 cast
}


def extract_source_type(node: dict) -> str | None:
    """从 CastExpr 的 inner[0] 里找 source type。"""
    for c in node.get("inner", []):
        if isinstance(c, dict):
            t = c.get("type")
            if isinstance(t, dict):
                return t.get("qualType") or t.get("desugaredQualType")
    return None


def simplify_type(qt: str) -> str:
    """压缩类型名：去掉命名空间 `clpoly::`、去掉模板参数里的具体类型（简化为类名）。
    用于直方图分组。"""
    if not qt:
        return "?"
    s = qt.strip()
    # 去掉 const/volatile 修饰
    s = s.replace("const ", "").replace("volatile ", "")
    # 去掉引用/指针
    s = s.rstrip("&*").strip()
    # 压缩 namespace
    s = s.replace("clpoly::", "")
    s = s.replace("std::__cxx11::", "std::")
    s = s.replace("std::_", "std::")
    # 模板参数压缩：保留类名 + "<...>"
    if "<" in s:
        base = s[: s.index("<")]
        # 保留简短的模板名
        # polynomial_<ZZ, lex_<less>> -> polynomial_<...>
        s = base + "<...>"
    return s


def walk(node, func_name, cast_pair_counts, cast_kind_counts, pair_hosts):
    if isinstance(node, dict):
        k = node.get("kind")
        if k in CAST_KINDS:
            target_qt = node.get("type", {}).get("qualType", "") if isinstance(node.get("type"), dict) else ""
            source_qt = extract_source_type(node) or ""
            target_s = simplify_type(target_qt)
            source_s = simplify_type(source_qt)
            cast_kind = node.get("castKind", "?")
            key = (source_s, target_s, cast_kind)
            cast_pair_counts[key] += 1
            pair_hosts[key].add(func_name)
            cast_kind_counts[cast_kind] += 1
        for c in node.get("inner", []):
            walk(c, func_name, cast_pair_counts, cast_kind_counts, pair_hosts)
    elif isinstance(node, list):
        for x in node:
            walk(x, func_name, cast_pair_counts, cast_kind_counts, pair_hosts)


def main():
    cast_pair_counts: Counter[tuple[str, str, str]] = Counter()
    cast_kind_counts: Counter[str] = Counter()
    pair_hosts: dict[tuple, set[str]] = defaultdict(set)

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        walk(ast, func_name, cast_pair_counts, cast_kind_counts, pair_hosts)

    total = sum(cast_pair_counts.values())

    lines = []
    lines.append("# 类型转换（Cast）枚举")
    lines.append("")
    lines.append(f"总计 **{total}** 次 cast（4 种 CastExpr kind）。按 (source, target, castKind) 三元组去重统计。")
    lines.append("")

    # CastKind 直方图
    lines.append("## castKind 直方图")
    lines.append("")
    lines.append("| castKind | 次数 | 含义 |")
    lines.append("|---|---|---|")
    castkind_hints = {
        "NoOp": "无操作（类型同质，可消除）",
        "LValueToRValue": "左值转右值（读取）",
        "IntegralCast": "整数类型转换（需检查位宽/符号）",
        "IntegralToBoolean": "整数转 bool",
        "BooleanToSignedIntegral": "bool 转整数",
        "FloatingToIntegral": "浮点转整数（截断）",
        "IntegralToFloating": "整数转浮点",
        "ConstructorConversion": "调构造函数",
        "UserDefinedConversion": "用户自定义转换",
        "DerivedToBase": "派生类到基类",
        "FunctionToPointerDecay": "函数名到函数指针",
        "ArrayToPointerDecay": "数组到指针",
        "NullToPointer": "nullptr 转指针",
    }
    for kind, cnt in cast_kind_counts.most_common():
        hint = castkind_hints.get(kind, "")
        lines.append(f"| `{kind}` | {cnt} | {hint} |")
    lines.append("")

    # (source, target) 对直方图
    lines.append("## (source, target, castKind) 三元组 — 前 80")
    lines.append("")
    lines.append("| Source | Target | CastKind | 次数 | 宿主数 |")
    lines.append("|---|---|---|---|---|")
    for (s, t, k), cnt in cast_pair_counts.most_common(80):
        hosts = pair_hosts[(s, t, k)]
        lines.append(f"| `{s}` | `{t}` | `{k}` | {cnt} | {len(hosts)} |")
    lines.append("")

    # 基础数值间的转换（重点分析）
    lines.append("## 基础数值类型间的转换（UB 分析重点）")
    lines.append("")
    NUMERIC = {"uint64_t", "int64_t", "int", "unsigned int", "unsigned",
               "long", "unsigned long", "long long", "unsigned long long",
               "bool", "char", "unsigned char", "size_t",
               "double", "float", "short", "unsigned short",
               "int32_t", "int16_t", "int8_t",
               "uint32_t", "uint16_t", "uint8_t"}
    lines.append("| Source | Target | Kind | 次数 |")
    lines.append("|---|---|---|---|")
    numeric_pairs = [(k, v) for k, v in cast_pair_counts.items()
                     if any(n in k[0] for n in NUMERIC) and any(n in k[1] for n in NUMERIC)]
    numeric_pairs.sort(key=lambda x: -x[1])
    for (s, t, k), cnt in numeric_pairs[:50]:
        lines.append(f"| `{s}` | `{t}` | `{k}` | {cnt} |")
    lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Total casts: {total}")
    print(f"Unique (source, target, kind) triples: {len(cast_pair_counts)}")
    print(f"Top 5 castKinds:")
    for k, c in cast_kind_counts.most_common(5):
        print(f"  {k}: {c}")


if __name__ == "__main__":
    main()
