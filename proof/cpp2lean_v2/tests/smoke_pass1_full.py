"""
Pass 1 parse 全量烟测 — 65 函数。

使用 Week 1 Day 1 survey 阶段产生的 AST cache（proof/cpp2lean/_ast_cache/）。
每 AST 对应一个已实例化的 FunctionDecl（mangledName ≠ None）。

输出：
  - 每函数的 Pass 1 结果（OK / FAIL: 原因 / WARN: Unknown 数）
  - UnknownExpr/UnknownStmt 的 kind 直方图
  - UnknownType 的 qualType 直方图
  - 引用参数自动判定 vs scan_ref_params.md 的对比
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path
from collections import Counter, defaultdict

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, RefType, UnknownType, TypeIR,
    UnknownStmt, UnknownExpr,
    HIRFunc, HIRParam,
    LambdaExpr,
    AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, LetStmt,
    Var, Lit, BinOp, UnaryOp, Call, UnresolvedOp,
    ArrayAccess, FieldAccess, Cast, BlockExpr, TupleExpr, ArrayLit, Capture,
    CondExpr, TranslationError,
)
from pass1_parse import parse_pass, assert_hir0_invariant, _walk_ir
from class_map import TRANSLATION_SCOPE, TRANSLATION_SCOPE_OUTPUT_PARAMS

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass1-smoke.md"


# ============================================================
# 统计辅助
# ============================================================

def collect_unknowns(func: HIRFunc) -> dict:
    """遍历 HIRFunc，收集 Unknown kind 直方图 + UnknownType 统计。"""
    stmt_unknowns = Counter()  # kind → count
    expr_unknowns = Counter()
    types_with_unknown = []  # (ctx_description, raw_qt)
    all_type_refs = 0
    unknown_type_refs = 0

    # 扫类型（在 params + body 里出现的 TypeIR）
    def check_type(ty: TypeIR, ctx: str):
        nonlocal all_type_refs, unknown_type_refs
        all_type_refs += 1
        if isinstance(ty, UnknownType):
            unknown_type_refs += 1
            types_with_unknown.append((ctx, ty.raw))
        elif isinstance(ty, RefType):
            check_type(ty.inner, ctx + ".ref")
        # 不深入 ArrayType/PairType 等（简化）

    # 参数类型
    for p in func.params:
        check_type(p.ty, f"param:{p.name}")
    # 返回类型
    check_type(func.ret_ty, "ret")

    # 扫 IR 节点找 Unknown
    def visit(node):
        nonlocal unknown_type_refs, all_type_refs
        if isinstance(node, UnknownStmt):
            stmt_unknowns[node.kind] += 1
        elif isinstance(node, UnknownExpr):
            expr_unknowns[node.kind] += 1
        # 检查节点的 ty 字段（若有）
        ty = getattr(node, "ty", None)
        if isinstance(ty, TypeIR.__args__):
            check_type(ty, f"expr:{type(node).__name__}")

    for s in func.body:
        _walk_ir(s, visit)

    return {
        "stmt_unknowns": dict(stmt_unknowns),
        "expr_unknowns": dict(expr_unknowns),
        "all_type_refs": all_type_refs,
        "unknown_type_refs": unknown_type_refs,
        "unknown_type_samples": types_with_unknown[:5],
    }


def detect_ref_params(func: HIRFunc) -> list[int]:
    """返回被标为 is_ref 的参数索引。"""
    return [i for i, p in enumerate(func.params) if p.is_ref]


# ============================================================
# 主入口
# ============================================================

def main():
    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1 smoke test on {len(targets)} functions...", file=sys.stderr)

    ok_count = 0
    warn_count = 0
    fail_count = 0

    per_func_stats = {}  # func → {stmts, exprs, unknown_types, ...}
    global_stmt_unknowns = Counter()
    global_expr_unknowns = Counter()
    global_unknown_types = Counter()
    ref_mismatches = []  # (func, auto_detected, configured)

    for func_name in targets:
        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            fail_count += 1
            per_func_stats[func_name] = {
                "status": "FAIL",
                "reason": "no AST cache",
            }
            continue

        try:
            with open(cache_file) as f:
                ast = json.load(f)
            hir = parse_pass(ast)
        except TranslationError as e:
            fail_count += 1
            per_func_stats[func_name] = {
                "status": "FAIL",
                "reason": f"TranslationError: {e.reason}",
            }
            continue
        except Exception as e:
            fail_count += 1
            tb = traceback.format_exc()
            per_func_stats[func_name] = {
                "status": "FAIL",
                "reason": f"{type(e).__name__}: {e}",
                "traceback": tb[:500],
            }
            continue

        # 成功 parse → 收集 Unknown 统计
        stats = collect_unknowns(hir)
        # 更新全局
        for k, c in stats["stmt_unknowns"].items():
            global_stmt_unknowns[k] += c
        for k, c in stats["expr_unknowns"].items():
            global_expr_unknowns[k] += c
        for _, raw in stats["unknown_type_samples"]:
            global_unknown_types[raw[:60]] += 1

        total_unknowns = (sum(stats["stmt_unknowns"].values())
                          + sum(stats["expr_unknowns"].values())
                          + stats["unknown_type_refs"])

        # ref 参数对比
        auto_refs = detect_ref_params(hir)
        configured = TRANSLATION_SCOPE_OUTPUT_PARAMS.get(func_name, [])
        if sorted(auto_refs) != sorted(configured):
            ref_mismatches.append((func_name, auto_refs, configured))

        status = "OK" if total_unknowns == 0 else "WARN"
        if status == "OK":
            ok_count += 1
        else:
            warn_count += 1

        per_func_stats[func_name] = {
            "status": status,
            "stmt_unknowns": stats["stmt_unknowns"],
            "expr_unknowns": stats["expr_unknowns"],
            "unknown_type_refs": stats["unknown_type_refs"],
            "total_unknowns": total_unknowns,
            "n_params": len(hir.params),
            "n_body_stmts": len(hir.body),
            "ret_ty": str(hir.ret_ty),
            "auto_refs": auto_refs,
            "configured_refs": configured,
            "ref_match": sorted(auto_refs) == sorted(configured),
        }

    # ======== 输出报告 ========
    lines = []
    lines.append("# Pass 1 parse 全量烟测 (65 函数)")
    lines.append("")
    lines.append(f"- OK (0 Unknown): **{ok_count}**")
    lines.append(f"- WARN (有 Unknown 节点): **{warn_count}**")
    lines.append(f"- FAIL (parse 崩溃): **{fail_count}**")
    lines.append("")

    # 全局 Unknown 直方图
    lines.append("## 全局 Unknown 节点直方图")
    lines.append("")
    lines.append("### Unknown Stmt kinds")
    lines.append("| kind | count |")
    lines.append("|---|---|")
    if global_stmt_unknowns:
        for k, c in global_stmt_unknowns.most_common():
            lines.append(f"| `{k}` | {c} |")
    else:
        lines.append("| (none) | 0 |")
    lines.append("")
    lines.append("### Unknown Expr kinds")
    lines.append("| kind | count |")
    lines.append("|---|---|")
    if global_expr_unknowns:
        for k, c in global_expr_unknowns.most_common():
            lines.append(f"| `{k}` | {c} |")
    else:
        lines.append("| (none) | 0 |")
    lines.append("")
    lines.append("### Unknown Type qualTypes (top 20)")
    lines.append("| qualType | count |")
    lines.append("|---|---|")
    if global_unknown_types:
        for qt, c in global_unknown_types.most_common(20):
            lines.append(f"| `{qt}` | {c} |")
    else:
        lines.append("| (none) | 0 |")
    lines.append("")

    # ref 参数一致性
    lines.append("## Ref 参数自动检测 vs TRANSLATION_SCOPE_OUTPUT_PARAMS")
    lines.append("")
    if not ref_mismatches:
        lines.append("✅ 全部一致（auto_detect 与 configured 完全匹配）")
    else:
        lines.append(f"⚠️  有 **{len(ref_mismatches)}** 个函数不匹配：")
        lines.append("")
        lines.append("| 函数 | auto | configured |")
        lines.append("|---|---|---|")
        for fn, auto, cfg in ref_mismatches:
            lines.append(f"| `{fn}` | {auto} | {cfg} |")
    lines.append("")

    # 按状态分组列出函数
    lines.append("## 函数明细")
    lines.append("")
    lines.append("### FAIL 列表")
    lines.append("")
    fail_list = [(fn, s) for fn, s in per_func_stats.items() if s["status"] == "FAIL"]
    if fail_list:
        for fn, s in fail_list:
            lines.append(f"- `{fn}`: {s['reason']}")
    else:
        lines.append("(无)")
    lines.append("")

    lines.append("### WARN 列表（按 total_unknowns 降序）")
    lines.append("")
    warn_list = [(fn, s) for fn, s in per_func_stats.items() if s["status"] == "WARN"]
    warn_list.sort(key=lambda x: -x[1]["total_unknowns"])
    if warn_list:
        lines.append("| 函数 | total Unknown | stmt kinds | expr kinds | unknown type refs |")
        lines.append("|---|---|---|---|---|")
        for fn, s in warn_list[:30]:  # 前 30
            stmt_kinds = ",".join(sorted(s["stmt_unknowns"].keys()))[:40]
            expr_kinds = ",".join(sorted(s["expr_unknowns"].keys()))[:40]
            lines.append(f"| `{fn}` | {s['total_unknowns']} | {stmt_kinds or '-'} | {expr_kinds or '-'} | {s['unknown_type_refs']} |")
    else:
        lines.append("(无)")
    lines.append("")

    lines.append("### OK 列表")
    lines.append("")
    ok_list = [fn for fn, s in per_func_stats.items() if s["status"] == "OK"]
    if ok_list:
        for fn in ok_list:
            lines.append(f"- `{fn}`")
    else:
        lines.append("(无)")

    OUT_MD.write_text("\n".join(lines))
    print(f"\n=== {ok_count} OK / {warn_count} WARN / {fail_count} FAIL ===", file=sys.stderr)
    print(f"\nGlobal unknowns:", file=sys.stderr)
    print(f"  Stmt: {dict(global_stmt_unknowns.most_common(10))}", file=sys.stderr)
    print(f"  Expr: {dict(global_expr_unknowns.most_common(10))}", file=sys.stderr)
    print(f"  Types: {dict(global_unknown_types.most_common(10))}", file=sys.stderr)
    print(f"  Ref mismatches: {len(ref_mismatches)}", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
