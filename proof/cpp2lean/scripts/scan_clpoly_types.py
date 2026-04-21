#!/usr/bin/env python3
"""scan_clpoly_types.py — 扫描 CLPoly 核心类型在 AST 中的使用模式。

为每个 CLPoly 类型（ZZ、QQ、Zp、Variable、polynomial_、upolynomial_、
factorization、basic_monomial、umonomial 等）统计：
- 作为变量/参数/返回类型的频次
- 被调用的方法列表（及频次）
- 被访问的字段列表
- 构造方式（默认构造、带参构造、拷贝、移动、poly_convert）
- 在哪些宿主函数中使用

为 type-system.md §3 CLPoly 核心类型设计提供输入。

输出：`docs/design/l1-translation-validation/survey/clpoly-types-usage.md`
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "clpoly-types-usage.md"

# 关心的 CLPoly 类型（前缀匹配）
CLPOLY_TYPES = [
    "clpoly::ZZ",
    "clpoly::QQ",
    "clpoly::Zp",
    "clpoly::variable",
    "clpoly::umonomial",
    "clpoly::basic_monomial",
    "clpoly::polynomial_",
    "clpoly::upolynomial_",
    "clpoly::factorization",
    "clpoly::lex_",
    "clpoly::grlex_",
    "clpoly::__hensel_node",
    "clpoly::__prime_selection_result",
    "clpoly::__wang_lc_result",
]


def classify(qt: str) -> str | None:
    """根据 qualType 返回归属的 CLPoly 类型（返回基础类型名，去掉模板参数）。"""
    if not qt:
        return None
    # 去修饰
    s = qt.replace("const ", "").rstrip("&*").strip()
    for t in CLPOLY_TYPES:
        if s.startswith(t):
            return t
    return None


def walk(node, func_name, type_stats):
    """遍历收集。type_stats 是 defaultdict，key 为类型名。"""
    if isinstance(node, dict):
        kind = node.get("kind")
        t = node.get("type", {}) if isinstance(node.get("type"), dict) else {}
        qt = t.get("qualType") or t.get("desugaredQualType", "")
        canonical = classify(qt)

        if canonical:
            stats = type_stats[canonical]
            stats["total"] += 1
            stats["hosts"].add(func_name)
            if kind == "VarDecl":
                stats["var"] += 1
            elif kind == "ParmVarDecl":
                stats["param"] += 1

        # 方法调用
        if kind == "CXXMemberCallExpr":
            # 找调用对象类型 + 方法名
            method_name = None
            obj_canonical = None
            for c in node.get("inner", []):
                if isinstance(c, dict) and c.get("kind") == "MemberExpr":
                    method_name = c.get("name", "")
                    # 对象类型
                    obj_t = c.get("type", {})
                    if isinstance(obj_t, dict):
                        obj_qt = obj_t.get("qualType", "")
                        obj_canonical = classify(obj_qt)
                    # 再深入 inner 找 obj
                    for cc in c.get("inner", []):
                        if isinstance(cc, dict):
                            ot = cc.get("type", {})
                            if isinstance(ot, dict):
                                ot_qt = ot.get("qualType", "")
                                oc = classify(ot_qt)
                                if oc:
                                    obj_canonical = oc
                                    break
                    break
            if obj_canonical and method_name:
                type_stats[obj_canonical]["methods"][method_name] += 1

        # 字段访问
        if kind == "MemberExpr":
            field_name = node.get("name", "")
            # 找 container type
            for c in node.get("inner", []):
                if isinstance(c, dict):
                    ct = c.get("type", {})
                    if isinstance(ct, dict):
                        cqt = ct.get("qualType", "")
                        oc = classify(cqt)
                        if oc and field_name:
                            type_stats[oc]["fields"][field_name] += 1
                    break

        # 构造函数
        if kind == "CXXConstructExpr":
            tqt = node.get("type", {}).get("qualType", "") if isinstance(node.get("type"), dict) else ""
            canonical_ctor = classify(tqt)
            if canonical_ctor:
                n_args = len([c for c in node.get("inner", [])
                             if isinstance(c, dict) and c.get("kind") not in ("CXXTemporaryExpr",)])
                ctor_kind = (
                    "default" if n_args == 0 else
                    "copy/move" if n_args == 1 else
                    f"{n_args}-arg"
                )
                type_stats[canonical_ctor]["ctors"][ctor_kind] += 1

        for c in node.get("inner", []):
            walk(c, func_name, type_stats)
    elif isinstance(node, list):
        for x in node:
            walk(x, func_name, type_stats)


def main():
    # 每类型的统计
    type_stats = defaultdict(lambda: {
        "total": 0,
        "var": 0,
        "param": 0,
        "hosts": set(),
        "methods": Counter(),
        "fields": Counter(),
        "ctors": Counter(),
    })

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        walk(ast, func_name, type_stats)

    lines = []
    lines.append("# CLPoly 核心类型使用模式扫描")
    lines.append("")
    lines.append("为 `type-system.md §3 CLPoly 核心类型` 设计提供输入。")
    lines.append("")

    # 类型汇总表
    lines.append("## 使用频次（总览）")
    lines.append("")
    lines.append("| 类型 | 总出现 | 变量 | 参数 | 方法数 | 字段数 | 宿主数 |")
    lines.append("|---|---|---|---|---|---|---|")
    sorted_types = sorted(type_stats.items(), key=lambda kv: -kv[1]["total"])
    for t, s in sorted_types:
        lines.append(f"| `{t}` | {s['total']} | {s['var']} | {s['param']} | {len(s['methods'])} | {len(s['fields'])} | {len(s['hosts'])} |")
    lines.append("")

    # 每类型详情
    for t, s in sorted_types:
        if s["total"] < 5:
            continue  # 跳过几乎没用到的
        lines.append(f"## `{t}`")
        lines.append("")
        lines.append(f"- 总出现: **{s['total']}**（var {s['var']} + param {s['param']}）")
        lines.append(f"- 宿主函数: **{len(s['hosts'])}** 个")
        lines.append("")
        if s["methods"]:
            lines.append("### 方法调用")
            lines.append("")
            lines.append("| 方法 | 调用次数 |")
            lines.append("|---|---|")
            for m, c in s["methods"].most_common(20):
                lines.append(f"| `.{m}()` | {c} |")
            lines.append("")
        if s["fields"]:
            lines.append("### 字段访问")
            lines.append("")
            lines.append("| 字段 | 访问次数 |")
            lines.append("|---|---|")
            for f, c in s["fields"].most_common(20):
                lines.append(f"| `.{f}` | {c} |")
            lines.append("")
        if s["ctors"]:
            lines.append("### 构造方式")
            lines.append("")
            lines.append("| 构造种类 | 次数 |")
            lines.append("|---|---|")
            for k, c in s["ctors"].most_common():
                lines.append(f"| `{k}` | {c} |")
            lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Types seen:")
    for t, s in sorted_types[:15]:
        top_methods = [m for m, _ in s["methods"].most_common(3)]
        print(f"  {t}: total={s['total']} hosts={len(s['hosts'])} methods={top_methods}")


if __name__ == "__main__":
    main()
