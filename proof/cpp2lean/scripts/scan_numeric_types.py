#!/usr/bin/env python3
"""scan_numeric_types.py — 基础数值类型使用模式扫描。

从缓存的 AST 中统计每个基础数值类型（uint64_t / int64_t / int / bool / double /
size_t 等）的：
- 作为变量类型出现次数
- 作为函数参数出现次数
- 作为函数返回类型出现次数
- 参与的运算符（从同一表达式的 BinaryOperator / UnaryOperator 附近推断）

为 type-system.md §1 (基础数值) 和 §2 (UB 分析) 提供输入。

输出：`docs/design/l1-translation-validation/survey/numeric-types.md`
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "numeric-types.md"


# 已知的基础数值类型（qualType 可能的形式）
NUMERIC_TYPES = {
    "uint64_t": ["uint64_t", "unsigned long", "unsigned long long"],
    "int64_t": ["int64_t", "long", "long long"],
    "uint32_t": ["uint32_t", "unsigned int", "unsigned"],
    "int32_t": ["int32_t", "int"],
    "uint16_t": ["uint16_t", "unsigned short"],
    "int16_t": ["int16_t", "short"],
    "uint8_t": ["uint8_t", "unsigned char"],
    "int8_t": ["int8_t", "signed char"],
    "size_t": ["size_t", "unsigned long"],  # 与 uint64 重叠
    "ptrdiff_t": ["ptrdiff_t"],
    "bool": ["bool", "_Bool"],
    "char": ["char"],
    "double": ["double"],
    "float": ["float"],
}


def classify(qt: str) -> str | None:
    """根据 qualType 的前导部分判断基础数值类型。"""
    if not qt:
        return None
    qt = qt.strip()
    # 去掉 const/volatile/& 修饰
    qt_clean = qt.replace("const ", "").replace("volatile ", "").rstrip("&*")
    qt_clean = qt_clean.strip()
    # 精确匹配
    for canonical, aliases in NUMERIC_TYPES.items():
        if qt_clean in aliases:
            return canonical
    return None


def walk(node, func_name, counts, var_types, ret_types, param_types):
    if isinstance(node, dict):
        kind = node.get("kind")
        t = node.get("type", {}) if isinstance(node.get("type"), dict) else {}
        qt = t.get("qualType") or t.get("desugaredQualType", "")
        canonical = classify(qt)
        if canonical:
            counts[canonical] += 1
            if kind == "VarDecl":
                var_types[canonical] += 1
            elif kind == "ParmVarDecl":
                param_types[canonical] += 1
        # FunctionDecl 返回类型：从 qualType 的开头部分
        if kind == "FunctionDecl":
            # qualType: "RetType (Param1, Param2)"
            if "(" in qt:
                ret = qt[: qt.index("(")].strip()
                # 去掉前导修饰
                ret_canonical = classify(ret)
                if ret_canonical:
                    ret_types[ret_canonical] += 1

        for c in node.get("inner", []):
            walk(c, func_name, counts, var_types, ret_types, param_types)
    elif isinstance(node, list):
        for x in node:
            walk(x, func_name, counts, var_types, ret_types, param_types)


def walk_operator_types(node, op_types):
    """统计每个 BinaryOperator/UnaryOperator/CompoundAssignOperator 的结果类型。"""
    if isinstance(node, dict):
        kind = node.get("kind")
        if kind in ("BinaryOperator", "UnaryOperator", "CompoundAssignOperator"):
            t = node.get("type", {})
            if isinstance(t, dict):
                qt = t.get("qualType") or t.get("desugaredQualType", "")
                canonical = classify(qt)
                if canonical:
                    op = node.get("opcode", "?")
                    op_types[(kind, op, canonical)] += 1
        for c in node.get("inner", []):
            walk_operator_types(c, op_types)
    elif isinstance(node, list):
        for x in node:
            walk_operator_types(x, op_types)


def main():
    total_counts = Counter()
    var_types = Counter()
    ret_types = Counter()
    param_types = Counter()
    op_types = Counter()
    per_func_types = defaultdict(Counter)

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        walk(ast, func_name, total_counts, var_types, ret_types, param_types)
        walk_operator_types(ast, op_types)

    lines = []
    lines.append("# 基础数值类型使用扫描")
    lines.append("")
    lines.append("## 全局频次")
    lines.append("")
    lines.append("| 类型 | 总出现 | 作为变量 | 作为参数 | 作为返回类型 |")
    lines.append("|---|---|---|---|---|")
    for t, _ in total_counts.most_common():
        lines.append(f"| `{t}` | {total_counts[t]} | {var_types[t]} | {param_types[t]} | {ret_types[t]} |")
    lines.append("")

    # 运算符结果类型分布
    lines.append("## 运算符结果类型（基础数值类型）")
    lines.append("")
    lines.append("| Operator | 结果类型 | 次数 |")
    lines.append("|---|---|---|")
    for (kind, op, t), cnt in sorted(op_types.items(), key=lambda x: -x[1]):
        lines.append(f"| `{kind}::{op}` | `{t}` | {cnt} |")
    lines.append("")

    # 每类型的运算符组合
    lines.append("## 按类型聚合的运算符")
    lines.append("")
    by_type = defaultdict(Counter)
    for (kind, op, t), cnt in op_types.items():
        by_type[t][(kind, op)] += cnt
    for t in sorted(by_type.keys(), key=lambda x: -sum(by_type[x].values())):
        lines.append(f"### `{t}` (合计 {sum(by_type[t].values())})")
        lines.append("")
        for (kind, op), cnt in by_type[t].most_common():
            lines.append(f"- `{kind}::{op}`: {cnt}")
        lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Types seen:")
    for t, c in total_counts.most_common():
        print(f"  {t}: total={c} (var={var_types[t]} param={param_types[t]} ret={ret_types[t]})")


if __name__ == "__main__":
    main()
