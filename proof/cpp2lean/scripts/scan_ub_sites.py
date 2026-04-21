#!/usr/bin/env python3
"""scan_ub_sites.py — 全 65 函数的 UB 站点枚举。

对每个缓存的 AST，扫描所有可能触发 UB 的操作，生成完整目录：

- UB-1 除以零：`BinaryOperator::/` 或 `%`（uint/int）
- UB-2 数组越界：`CXXOperatorCallExpr::operator[]` + `.at`（如有）
- UB-3 空容器 front/back：`CXXMemberCallExpr` 方法名 == `front`/`back`/`front!`/`back!`
- UB-4 移位越界：`BinaryOperator::<<` 或 `>>`（非流输出，即整数类型）
- UB-6 有符号溢出：`BinaryOperator::+`/`-`/`*`（signed int 类型）
- UB-7 unsigned→signed 截断：Cast 中 uint → int（限比特数）
- UB-8 assert：`CallExpr` callee 名 == `__assert_fail` 或 `assert`

为每站点记录：
- 宿主函数
- 源码行号（从 AST range.begin.line）
- UB 类型
- 可能还需要的 require 条件摘要

输出：`docs/design/l1-translation-validation/survey/ub-sites.md`
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent.parent
AST_CACHE_DIR = HERE / "_ast_cache"
OUT_MD = HERE.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "ub-sites.md"


def is_integer_type(qt: str) -> bool:
    """判断 qualType 是否整数类型。"""
    if not qt:
        return False
    qt = qt.strip().replace("const ", "").rstrip("&*").strip()
    return any(qt == n for n in [
        "uint64_t", "int64_t", "uint32_t", "int32_t",
        "uint16_t", "int16_t", "uint8_t", "int8_t",
        "int", "unsigned int", "unsigned", "long", "unsigned long",
        "long long", "unsigned long long", "short", "unsigned short",
        "size_t", "ptrdiff_t", "char", "unsigned char", "signed char",
    ])


def is_signed_int_type(qt: str) -> bool:
    qt = (qt or "").strip().replace("const ", "").rstrip("&*").strip()
    return qt in {"int64_t", "int32_t", "int16_t", "int8_t",
                  "int", "long", "long long", "short", "signed char",
                  "ptrdiff_t"}


def get_line(node) -> int | None:
    r = node.get("range") if isinstance(node, dict) else None
    if isinstance(r, dict):
        b = r.get("begin") or r.get("line") or {}
        if isinstance(b, dict):
            return b.get("line")
        elif isinstance(r, dict):
            return r.get("line")
    return None


def walk(node, func_name, sites):
    """遍历 AST，把 UB 站点收集到 sites（按 UB 类型分组）。"""
    if isinstance(node, dict):
        kind = node.get("kind")
        opcode = node.get("opcode")
        result_qt = node.get("type", {}).get("qualType", "") if isinstance(node.get("type"), dict) else ""
        line = get_line(node)

        # UB-1: 除以零
        if kind == "BinaryOperator" and opcode in ("/", "%"):
            if is_integer_type(result_qt):
                sites["UB-1 Div0"].append({
                    "func": func_name,
                    "line": line,
                    "op": opcode,
                    "type": result_qt,
                })
        # UB-1 (operator overload): CXXOperatorCallExpr::operator/ or operator%
        if kind == "CXXOperatorCallExpr":
            # 取 referenced operator 名
            for c in node.get("inner", []):
                if isinstance(c, dict):
                    rd = None
                    if c.get("kind") == "ImplicitCastExpr":
                        for cc in c.get("inner", []):
                            if isinstance(cc, dict) and cc.get("kind") == "DeclRefExpr":
                                rd = cc.get("referencedDecl")
                    elif c.get("kind") == "DeclRefExpr":
                        rd = c.get("referencedDecl")
                    if isinstance(rd, dict):
                        name = rd.get("name", "")
                        if name in ("operator/", "operator%"):
                            sites["UB-1 Div0"].append({
                                "func": func_name,
                                "line": line,
                                "op": name,
                                "type": result_qt,
                            })
                        elif name in ("operator[]",):
                            sites["UB-2 OOB"].append({
                                "func": func_name,
                                "line": line,
                                "op": name,
                                "type": result_qt,
                            })
                        elif name in ("operator<<", "operator>>"):
                            # 仅整数类型（排除流操作）
                            if is_integer_type(result_qt):
                                sites["UB-4 Shift"].append({
                                    "func": func_name,
                                    "line": line,
                                    "op": name,
                                    "type": result_qt,
                                })
                        break

        # UB-2: 数组越界（BinaryOperator 里没有，在 CXXOperatorCallExpr[] 或 ArraySubscriptExpr）
        if kind == "ArraySubscriptExpr":
            sites["UB-2 OOB"].append({
                "func": func_name,
                "line": line,
                "op": "[]",
                "type": result_qt,
            })

        # UB-3: 空容器 front/back / front!/back!
        if kind == "CXXMemberCallExpr":
            # 找方法名
            for c in node.get("inner", []):
                if isinstance(c, dict) and c.get("kind") == "MemberExpr":
                    name = c.get("name", "")
                    if name in ("front", "back", "front!", "back!"):
                        sites["UB-3 EmptyContainer"].append({
                            "func": func_name,
                            "line": line,
                            "method": name,
                            "type": result_qt,
                        })
                    break

        # UB-4: 移位越界（基础类型）
        if kind == "BinaryOperator" and opcode in ("<<", ">>"):
            if is_integer_type(result_qt):
                sites["UB-4 Shift"].append({
                    "func": func_name,
                    "line": line,
                    "op": opcode,
                    "type": result_qt,
                })

        # UB-6: 有符号溢出
        if kind == "BinaryOperator" and opcode in ("+", "-", "*"):
            if is_signed_int_type(result_qt):
                sites["UB-6 SignedOverflow"].append({
                    "func": func_name,
                    "line": line,
                    "op": opcode,
                    "type": result_qt,
                })
        if kind == "CompoundAssignOperator" and opcode in ("+=", "-=", "*="):
            if is_signed_int_type(result_qt):
                sites["UB-6 SignedOverflow"].append({
                    "func": func_name,
                    "line": line,
                    "op": opcode,
                    "type": result_qt,
                })
        if kind == "UnaryOperator" and opcode == "-":
            if is_signed_int_type(result_qt):
                sites["UB-6 SignedOverflow"].append({
                    "func": func_name,
                    "line": line,
                    "op": "unary-",
                    "type": result_qt,
                })

        # UB-7: unsigned→signed cast（IntegralCast）
        if kind in ("ImplicitCastExpr", "CStyleCastExpr",
                    "CXXStaticCastExpr", "CXXFunctionalCastExpr"):
            if node.get("castKind") == "IntegralCast":
                target_qt = result_qt
                source_qt = None
                for c in node.get("inner", []):
                    if isinstance(c, dict):
                        t = c.get("type", {})
                        if isinstance(t, dict):
                            source_qt = t.get("qualType") or t.get("desugaredQualType")
                        break
                # unsigned → signed?
                UNSIGNED = {"uint64_t", "uint32_t", "uint16_t", "uint8_t",
                            "unsigned int", "unsigned", "unsigned long",
                            "unsigned long long", "unsigned short",
                            "size_t", "unsigned char"}
                src_clean = (source_qt or "").replace("const ", "").rstrip("&*").strip()
                tgt_clean = target_qt.replace("const ", "").rstrip("&*").strip()
                if src_clean in UNSIGNED and is_signed_int_type(tgt_clean):
                    sites["UB-7 UnsignedToSigned"].append({
                        "func": func_name,
                        "line": line,
                        "source": src_clean,
                        "target": tgt_clean,
                    })

        # UB-8: assert
        if kind == "CallExpr":
            # 找 callee
            for c in node.get("inner", []):
                if isinstance(c, dict):
                    rd = None
                    if c.get("kind") == "ImplicitCastExpr":
                        for cc in c.get("inner", []):
                            if isinstance(cc, dict) and cc.get("kind") == "DeclRefExpr":
                                rd = cc.get("referencedDecl")
                    elif c.get("kind") == "DeclRefExpr":
                        rd = c.get("referencedDecl")
                    if isinstance(rd, dict):
                        name = rd.get("name", "")
                        if name == "__assert_fail":
                            sites["UB-8 Assert"].append({
                                "func": func_name,
                                "line": line,
                            })
                        break

        for c in node.get("inner", []):
            walk(c, func_name, sites)
    elif isinstance(node, list):
        for x in node:
            walk(x, func_name, sites)


def main():
    all_sites = defaultdict(list)  # UB 类型 → [sites]

    for cache_file in sorted(AST_CACHE_DIR.glob("*.json")):
        func_name = cache_file.stem
        with open(cache_file) as f:
            ast = json.load(f)
        walk(ast, func_name, all_sites)

    total = sum(len(v) for v in all_sites.values())

    # 分组统计
    by_func = defaultdict(lambda: Counter())
    for ub_type, sites in all_sites.items():
        for s in sites:
            by_func[s["func"]][ub_type] += 1

    lines = []
    lines.append("# UB 站点全量枚举")
    lines.append("")
    lines.append("> 参考 `ub-catalog.md` 的 UB-1 到 UB-8 分类定义。")
    lines.append("> 本文件是 65 函数的**实测**UB 站点清单。")
    lines.append("")
    lines.append("## 全局统计")
    lines.append("")
    lines.append("| UB 类型 | 站点数 |")
    lines.append("|---|---|")
    for ub_type in sorted(all_sites.keys()):
        lines.append(f"| {ub_type} | {len(all_sites[ub_type])} |")
    lines.append(f"| **合计** | **{total}** |")
    lines.append("")

    # 按宿主函数汇总
    lines.append("## 按函数汇总（UB 站点数）")
    lines.append("")
    lines.append("| 函数 | UB-1 Div | UB-2 OOB | UB-3 Empty | UB-4 Shift | UB-6 Signed | UB-7 U→S | UB-8 Assert | 合计 |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    sorted_funcs = sorted(by_func.items(), key=lambda kv: -sum(kv[1].values()))
    for fn, cnts in sorted_funcs:
        row_total = sum(cnts.values())
        lines.append(
            f"| `{fn}` | {cnts.get('UB-1 Div0', 0)} "
            f"| {cnts.get('UB-2 OOB', 0)} "
            f"| {cnts.get('UB-3 EmptyContainer', 0)} "
            f"| {cnts.get('UB-4 Shift', 0)} "
            f"| {cnts.get('UB-6 SignedOverflow', 0)} "
            f"| {cnts.get('UB-7 UnsignedToSigned', 0)} "
            f"| {cnts.get('UB-8 Assert', 0)} "
            f"| **{row_total}** |"
        )
    lines.append("")

    # 每 UB 类型的典型站点（前 10）
    for ub_type in sorted(all_sites.keys()):
        sites = all_sites[ub_type]
        lines.append(f"## {ub_type} 站点样本（前 10）")
        lines.append("")
        by_host = defaultdict(int)
        for s in sites:
            by_host[s["func"]] += 1
        lines.append(f"总计 {len(sites)}，分布于 {len(by_host)} 个函数。")
        lines.append("")
        lines.append("| 宿主 | 行 | 细节 |")
        lines.append("|---|---|---|")
        for s in sites[:10]:
            detail_keys = [k for k in s if k not in ("func", "line")]
            detail = " / ".join(f"{k}=`{s[k]}`" for k in detail_keys)
            lines.append(f"| `{s['func']}` | {s.get('line')} | {detail} |")
        lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Total UB sites: {total}")
    for ub_type in sorted(all_sites.keys()):
        print(f"  {ub_type}: {len(all_sites[ub_type])}")


if __name__ == "__main__":
    main()
