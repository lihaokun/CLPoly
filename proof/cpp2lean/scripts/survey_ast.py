#!/usr/bin/env python3
"""
survey_ast.py — 调研 67 个目标函数的 Clang AST，输出构造统计。

用途：为 translator v2 重构的 Stage 1 Week 1 提供全量 AST 构造清单。

输出：
  docs/design/l1-translation-validation/survey/
    ├── ast-kinds.md       — 所有 AST node kind 的直方图 + 首次出现
    ├── operators.md        — 所有运算符的直方图（CXXOperatorCallExpr / UnaryOperator
    │                         / BinaryOperator / CompoundAssignOperator）
    ├── types.md            — 所有类型（type.qualType）的直方图
    └── summary.md          — 总览 + 潜在难点提示

缓存：
  proof/cpp2lean/_ast_cache/<func_name>.json — 每函数的 Clang AST（去除
                                               base64 binary 字段后的精简版）

用法：
  cd proof/cpp2lean
  python3 scripts/survey_ast.py              # 增量（跳过已缓存函数）
  python3 scripts/survey_ast.py --refresh    # 强制重 dump
"""

from __future__ import annotations

import os
import sys
import json
import time
import argparse
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


# ----------------------------------------------------------------------
# 路径配置
# ----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent.parent
PROJECT_ROOT = HERE.parent.parent

INSTANTIATE_CC = HERE / "instantiate.cc"
AST_CACHE_DIR = HERE / "_ast_cache"
SURVEY_DIR = PROJECT_ROOT / "docs" / "design" / "l1-translation-validation" / "survey"


# ----------------------------------------------------------------------
# 目标函数清单（从 class_map.py 复用）
# ----------------------------------------------------------------------
sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE  # noqa: E402


# ----------------------------------------------------------------------
# Clang 调用
# ----------------------------------------------------------------------
def _system_includes() -> list[str]:
    """从 g++ 拿系统 include 路径 + libclang 内建头文件。"""
    try:
        result = subprocess.run(
            ["g++", "-E", "-Wp,-v", "-x", "c++", "-"],
            input="", capture_output=True, text=True
        )
        paths = []
        capture = False
        for line in result.stderr.splitlines():
            if "#include <...>" in line:
                capture = True
                continue
            if line.startswith("End of search"):
                break
            if capture and line.startswith(" "):
                paths.append(f"-I{line.strip()}")
        import glob
        for d in glob.glob("/usr/lib/llvm-*/lib/clang/*/include"):
            paths.append(f"-I{d}")
        return paths
    except Exception:
        return []


SYS_INCLUDES = _system_includes()


def dump_ast_json(func_name: str, timeout: int = 30) -> dict | None:
    """通过 clang++ -ast-dump-filter 拿单个函数的 AST JSON。

    filter 是子串匹配，可能命中多个；我们选 name 精确匹配的。
    对模板函数，优先选有 CompoundStmt body 的实例化版本。
    """
    cmd = [
        "clang++",
        "-std=c++17",
        f"-I{PROJECT_ROOT}",
    ] + SYS_INCLUDES + [
        "-Xclang", "-ast-dump=json",
        "-Xclang", f"-ast-dump-filter={func_name}",
        "-fsyntax-only",
        str(INSTANTIATE_CC),
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            cwd=HERE
        )
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT", file=sys.stderr)
        return None

    if not result.stdout.strip():
        print(f"  empty stdout", file=sys.stderr)
        return None

    # filter 子串匹配可能输出多个顶层 JSON 对象，逐个解析，选精确匹配
    decoder = json.JSONDecoder()
    raw = result.stdout.strip()
    pos = 0
    candidates = []
    while pos < len(raw):
        try:
            obj, end = decoder.raw_decode(raw[pos:])
            candidates.append(obj)
            pos += end
            while pos < len(raw) and raw[pos] in " \t\n\r":
                pos += 1
        except json.JSONDecodeError:
            break

    # has_body: 检查 FunctionDecl 是否有 CompoundStmt body（非仅声明）
    def has_body(node: dict) -> bool:
        if node.get("kind") == "FunctionDecl":
            for c in node.get("inner", []):
                if isinstance(c, dict) and c.get("kind") == "CompoundStmt":
                    return True
            return False
        if node.get("kind") == "FunctionTemplateDecl":
            for child in node.get("inner", []):
                if isinstance(child, dict) and has_body(child):
                    return True
        return False

    # 展开 FunctionTemplateDecl 的所有内嵌 FunctionDecl（含模板定义与实例化）
    def expand(node: dict) -> list[dict]:
        if node.get("kind") == "FunctionTemplateDecl":
            out = []
            for child in node.get("inner", []):
                if isinstance(child, dict) and child.get("kind") == "FunctionDecl":
                    out.append(child)
            return out
        return [node]

    # 先挑精确名的候选
    exact_named = [c for c in candidates if c.get("name") == func_name]

    # 展开模板实例化
    expanded = []
    for c in exact_named:
        expanded.extend(expand(c))

    # 优先级：
    # 1. 有 mangledName + 有 body（= 已实例化且带实现）
    # 2. 有 mangledName（即使无 body，意味着它至少是具体声明，不是模板定义）
    # 3. 有 body（模板定义 fallback）
    # 4. 任一
    mangled_with_body = [n for n in expanded
                         if n.get("mangledName") and has_body(n)]
    if mangled_with_body:
        return mangled_with_body[-1]  # 若有多个实例化，取最后一个

    mangled_only = [n for n in expanded if n.get("mangledName")]
    if mangled_only:
        return mangled_only[-1]

    with_body = [n for n in expanded if has_body(n)]
    if with_body:
        return with_body[-1]

    if expanded:
        return expanded[0]
    return None


# ----------------------------------------------------------------------
# 缓存
# ----------------------------------------------------------------------
def _strip_large(obj: object, max_string: int = 200) -> object:
    """缩减 AST：去掉 range/loc 细节和超长字符串，方便 git 存储。"""
    if isinstance(obj, dict):
        # 保留关键字段，去掉位置细节里的 offset/tokLen 等噪声
        out = {}
        for k, v in obj.items():
            if k in ("range", "loc"):
                # 只保留 begin/end 的 line 和 col 基本信息
                if isinstance(v, dict):
                    if "begin" in v:
                        b = v["begin"]
                        if isinstance(b, dict):
                            out[k] = {
                                "line": b.get("line"),
                                "col": b.get("col"),
                            }
                continue
            out[k] = _strip_large(v, max_string)
        return out
    if isinstance(obj, list):
        return [_strip_large(x, max_string) for x in obj]
    if isinstance(obj, str) and len(obj) > max_string:
        return obj[:max_string] + "...<truncated>"
    return obj


def cache_path(func_name: str) -> Path:
    # 清理文件名
    safe = func_name.replace("/", "_")
    return AST_CACHE_DIR / f"{safe}.json"


def get_ast(func_name: str, refresh: bool = False) -> dict | None:
    p = cache_path(func_name)
    if not refresh and p.exists():
        with open(p) as f:
            return json.load(f)

    print(f"  dumping {func_name}...", end="", flush=True, file=sys.stderr)
    t0 = time.time()
    ast = dump_ast_json(func_name)
    dt = time.time() - t0
    if ast is None:
        print(f" FAILED ({dt:.1f}s)", file=sys.stderr)
        return None
    slim = _strip_large(ast)
    with open(p, "w") as f:
        json.dump(slim, f, indent=2, ensure_ascii=False)
    print(f" OK ({dt:.1f}s, {len(json.dumps(slim))//1024}KB)", file=sys.stderr)
    return slim


# ----------------------------------------------------------------------
# 统计
# ----------------------------------------------------------------------
class Stats:
    def __init__(self):
        self.kind_counts: Counter[str] = Counter()
        self.kind_first: dict[str, str] = {}  # kind → 首次出现的函数名
        self.kind_example: dict[str, str] = {}  # kind → 样例 JSON 片段
        self.operator_counts: Counter[str] = Counter()  # "CXXOperatorCallExpr::+"
        self.operator_first: dict[str, str] = {}
        self.type_counts: Counter[str] = Counter()
        self.type_first: dict[str, str] = {}
        self.func_kinds: dict[str, set[str]] = defaultdict(set)  # 每函数用到的 kinds

    def walk(self, node: object, func_name: str, depth: int = 0):
        if isinstance(node, dict):
            kind = node.get("kind")
            if kind:
                self.kind_counts[kind] += 1
                self.func_kinds[func_name].add(kind)
                if kind not in self.kind_first:
                    self.kind_first[kind] = func_name
                    # 存一个样例（去掉 inner 避免太大）
                    sample = {k: v for k, v in node.items() if k != "inner"}
                    self.kind_example[kind] = json.dumps(sample, ensure_ascii=False)[:400]

                # 运算符专项统计
                if kind == "CXXOperatorCallExpr":
                    op = self._extract_operator_name(node)
                    key = f"CXXOperatorCallExpr::{op}"
                    self.operator_counts[key] += 1
                    if key not in self.operator_first:
                        self.operator_first[key] = func_name
                elif kind in ("UnaryOperator", "BinaryOperator",
                              "CompoundAssignOperator"):
                    op = node.get("opcode", "?")
                    key = f"{kind}::{op}"
                    self.operator_counts[key] += 1
                    if key not in self.operator_first:
                        self.operator_first[key] = func_name

            # 类型统计
            t = node.get("type")
            if isinstance(t, dict):
                qt = t.get("qualType")
                if qt:
                    self.type_counts[qt] += 1
                    if qt not in self.type_first:
                        self.type_first[qt] = func_name

            # 递归 inner
            for c in node.get("inner", []):
                self.walk(c, func_name, depth + 1)

        elif isinstance(node, list):
            for x in node:
                self.walk(x, func_name, depth + 1)

    def _extract_operator_name(self, node: dict) -> str:
        """从 CXXOperatorCallExpr 提取运算符名（operator+, operator[], operator() 等）。"""
        # 第一个 inner 通常是 DeclRefExpr 到 operator 函数
        for c in node.get("inner", []):
            if isinstance(c, dict):
                rd = c.get("referencedDecl")
                if isinstance(rd, dict):
                    name = rd.get("name", "")
                    if name.startswith("operator"):
                        return name
                # ImplicitCastExpr 包 DeclRefExpr
                for cc in c.get("inner", []):
                    if isinstance(cc, dict):
                        rd2 = cc.get("referencedDecl")
                        if isinstance(rd2, dict):
                            name = rd2.get("name", "")
                            if name.startswith("operator"):
                                return name
        return "?"


# ----------------------------------------------------------------------
# 输出
# ----------------------------------------------------------------------
def _md_table(header: list[str], rows: list[list[str]]) -> str:
    out = []
    out.append("| " + " | ".join(header) + " |")
    out.append("|" + "|".join(["---"] * len(header)) + "|")
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def write_ast_kinds(stats: Stats, out_path: Path):
    lines = []
    lines.append("# AST Kind 直方图")
    lines.append("")
    lines.append(f"> 生成时间：{time.strftime('%Y-%m-%d %H:%M')}")
    lines.append(f"> 覆盖 {len(stats.func_kinds)} 个函数")
    lines.append("")
    lines.append("按出现次数降序。每 kind 列首次出现的函数名（样例）。")
    lines.append("")

    rows = []
    for kind, count in stats.kind_counts.most_common():
        rows.append([kind, count, stats.kind_first.get(kind, "?")])
    lines.append(_md_table(["Kind", "Count", "First Seen In"], rows))
    lines.append("")

    # 罕见 kind 的样例（出现 <= 3 次的）
    rare = [k for k, c in stats.kind_counts.items() if c <= 3]
    if rare:
        lines.append("## 罕见 Kind 样例（≤3 次出现）")
        lines.append("")
        for k in sorted(rare):
            lines.append(f"### `{k}` （{stats.kind_counts[k]} 次，首现 `{stats.kind_first.get(k)}`）")
            lines.append("")
            lines.append("```json")
            lines.append(stats.kind_example.get(k, "<no sample>"))
            lines.append("```")
            lines.append("")

    out_path.write_text("\n".join(lines))


def write_operators(stats: Stats, out_path: Path):
    lines = []
    lines.append("# 运算符直方图")
    lines.append("")
    lines.append(f"> 生成时间：{time.strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("包括 `CXXOperatorCallExpr`（用户类型重载）、`UnaryOperator`（基本类型 `!` `-` `~` `++` 等）、`BinaryOperator`（`+ - * / == ` 等）、`CompoundAssignOperator`（`+= -= *=` 等）。")
    lines.append("")

    rows = []
    for key, count in stats.operator_counts.most_common():
        rows.append([key, count, stats.operator_first.get(key, "?")])
    lines.append(_md_table(["Operator", "Count", "First Seen In"], rows))

    out_path.write_text("\n".join(lines))


def write_types(stats: Stats, out_path: Path):
    lines = []
    lines.append("# 类型（qualType）直方图")
    lines.append("")
    lines.append(f"> 生成时间：{time.strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("按出现次数降序。截取前 200 条。")
    lines.append("")

    rows = []
    for qt, count in stats.type_counts.most_common(200):
        rows.append([f"`{qt[:80]}`" + ("..." if len(qt) > 80 else ""),
                     count, stats.type_first.get(qt, "?")])
    lines.append(_md_table(["Type", "Count", "First Seen In"], rows))

    # 基础数值/容器类型分组
    lines.append("")
    lines.append("## 类型分组")
    lines.append("")
    groups = {
        "基础数值": ["int", "long", "unsigned", "short", "char", "bool",
                 "float", "double", "size_t", "ptrdiff_t",
                 "int64_t", "uint64_t", "int32_t", "uint32_t"],
        "CLPoly 数域": ["clpoly::ZZ", "clpoly::QQ", "clpoly::Zp", "clpoly::Zm"],
        "CLPoly 多项式": ["clpoly::polynomial_", "clpoly::upolynomial_",
                      "clpoly::basic_monomial", "clpoly::variable"],
        "STL 容器": ["std::vector", "std::pair", "std::map", "std::set",
                   "std::tuple", "std::array", "std::list"],
        "STL 迭代器/引用": ["std::__", "iterator"],
        "STL 算法/随机": ["std::sort", "std::swap", "std::max", "std::min",
                      "std::random", "std::uniform"],
        "GMP": ["mpz_", "__mpz_"],
        "__int128": ["__int128"],
    }
    for group_name, keywords in groups.items():
        matched = [(qt, c) for qt, c in stats.type_counts.items()
                   if any(k in qt for k in keywords)]
        matched.sort(key=lambda x: -x[1])
        if not matched:
            continue
        lines.append(f"### {group_name}")
        lines.append("")
        for qt, c in matched[:30]:
            lines.append(f"- `{qt}` — {c}")
        lines.append("")

    out_path.write_text("\n".join(lines))


def write_summary(stats: Stats, missing: list[str], out_path: Path):
    lines = []
    lines.append("# AST 调研总览")
    lines.append("")
    lines.append(f"> 生成时间：{time.strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("## 覆盖情况")
    lines.append("")
    lines.append(f"- 目标函数总数（TRANSLATION_SCOPE）: **{len(TRANSLATION_SCOPE)}**")
    lines.append(f"- 成功 dump AST 的函数数: **{len(stats.func_kinds)}**")
    lines.append(f"- 失败/缺失的函数数: **{len(missing)}**")
    if missing:
        lines.append("")
        lines.append("### 失败/缺失函数")
        for m in sorted(missing):
            lines.append(f"- `{m}`")
    lines.append("")
    lines.append("## 全局统计")
    lines.append("")
    lines.append(f"- 不同的 AST kind 种数: **{len(stats.kind_counts)}**")
    lines.append(f"- 不同的运算符种数: **{len(stats.operator_counts)}**")
    lines.append(f"- 不同的类型（qualType）种数: **{len(stats.type_counts)}**")
    lines.append(f"- AST 节点总数（累计）: **{sum(stats.kind_counts.values())}**")
    lines.append("")
    lines.append("## 潜在难点提示")
    lines.append("")

    lines.append("### 控制流构造")
    cf_kinds = ["IfStmt", "WhileStmt", "ForStmt", "CXXForRangeStmt",
                "DoStmt", "BreakStmt", "ContinueStmt", "ReturnStmt",
                "SwitchStmt", "CaseStmt", "DefaultStmt",
                "CXXTryStmt", "CXXThrowExpr", "CXXCatchStmt", "GotoStmt"]
    for k in cf_kinds:
        c = stats.kind_counts.get(k, 0)
        mark = "✅" if c > 0 else "—"
        lines.append(f"- {mark} `{k}`: {c}")
    lines.append("")

    lines.append("### Lambda / 闭包")
    lam_kinds = ["LambdaExpr", "CXXRecordDecl", "CallExpr", "CXXOperatorCallExpr"]
    for k in lam_kinds:
        c = stats.kind_counts.get(k, 0)
        mark = "✅" if c > 0 else "—"
        lines.append(f"- {mark} `{k}`: {c}")
    lines.append("")

    lines.append("### 结构化绑定 / 迭代器")
    sb_kinds = ["DecompositionDecl", "BindingDecl", "CXXBindTemporaryExpr",
                "MaterializeTemporaryExpr", "CXXConstructExpr",
                "CXXDependentScopeMemberExpr"]
    for k in sb_kinds:
        c = stats.kind_counts.get(k, 0)
        mark = "✅" if c > 0 else "—"
        lines.append(f"- {mark} `{k}`: {c}")
    lines.append("")

    lines.append("### 模板相关")
    tpl_kinds = ["FunctionTemplateDecl", "ClassTemplateDecl",
                 "ClassTemplateSpecializationDecl",
                 "TemplateTypeParmDecl", "UnresolvedLookupExpr",
                 "CXXDependentScopeMemberExpr", "DependentScopeDeclRefExpr"]
    for k in tpl_kinds:
        c = stats.kind_counts.get(k, 0)
        mark = "✅" if c > 0 else "—"
        lines.append(f"- {mark} `{k}`: {c}")
    lines.append("")

    lines.append("### 内存/生命周期")
    mem_kinds = ["CXXNewExpr", "CXXDeleteExpr", "ExprWithCleanups",
                 "CXXThisExpr", "UnresolvedMemberExpr"]
    for k in mem_kinds:
        c = stats.kind_counts.get(k, 0)
        mark = "✅" if c > 0 else "—"
        lines.append(f"- {mark} `{k}`: {c}")
    lines.append("")

    lines.append("## 每函数使用的 kind 数量")
    lines.append("")
    rows = [(fn, len(ks)) for fn, ks in stats.func_kinds.items()]
    rows.sort(key=lambda x: -x[1])
    lines.append(_md_table(["Function", "Distinct Kinds"], rows))

    out_path.write_text("\n".join(lines))


# ----------------------------------------------------------------------
# 主入口
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true",
                        help="强制重 dump 所有函数（忽略缓存）")
    parser.add_argument("--only", nargs="*",
                        help="只处理指定函数（调试用）")
    args = parser.parse_args()

    AST_CACHE_DIR.mkdir(exist_ok=True)
    SURVEY_DIR.mkdir(parents=True, exist_ok=True)

    targets = sorted(TRANSLATION_SCOPE)
    if args.only:
        targets = [t for t in targets if t in args.only]

    print(f"Target: {len(targets)} functions", file=sys.stderr)

    stats = Stats()
    missing = []

    for i, func_name in enumerate(targets, 1):
        print(f"[{i}/{len(targets)}] {func_name}", file=sys.stderr)
        ast = get_ast(func_name, refresh=args.refresh)
        if ast is None:
            missing.append(func_name)
            continue
        stats.walk(ast, func_name)

    # 输出
    write_ast_kinds(stats, SURVEY_DIR / "ast-kinds.md")
    write_operators(stats, SURVEY_DIR / "operators.md")
    write_types(stats, SURVEY_DIR / "types.md")
    write_summary(stats, missing, SURVEY_DIR / "summary.md")

    print(f"\nOutput in: {SURVEY_DIR}", file=sys.stderr)
    print(f"Cache in: {AST_CACHE_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
