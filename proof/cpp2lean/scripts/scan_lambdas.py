#!/usr/bin/env python3
"""scan_lambdas.py — 全量统计 TRANSLATION_SCOPE 里每个 Lambda 的特征。

Clang JSON 对 LambdaExpr 的 captures 字段做了省略，不直接可用。
改用 **源码正则扫描 + AST 交叉验证** 的混合策略。

输出 `docs/design/l1-translation-validation/survey/lambdas.md`：
- 每 Lambda 的宿主函数 + 源位置
- 捕获模式（[&] / [=] / 具名 / 混合 / 空）+ 具体捕获清单
- 是否 generic（auto 参数）
- body 统计（顶层语句数）+ 行数
- 调用方式分类（std::sort 比较器 / stored / 其他）
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
OUT_MD = PROJECT_ROOT / "docs" / "design" / "l1-translation-validation" / "survey" / "lambdas.md"

CLPOLY_DIR = PROJECT_ROOT / "clpoly"

# --- 步骤 1: 从 C++ 源码提取所有 lambda ---

# Lambda 语法：`[CAPTURE](PARAMS) { BODY }` 或 `[CAPTURE](PARAMS) -> RET { BODY }`
# CAPTURE:
#   []   — 无捕获
#   [&]  — 默认引用捕获
#   [=]  — 默认值捕获
#   [&, x, y]  — 默认引用，y 值
#   [=, &x]    — 默认值，x 引用
#   [x, y, &z] — 具名混合
#   [this] — 捕获 this
# 我们用松散正则：找 `[...]` 后紧跟 `(` 的模式

LAMBDA_HEAD_RE = re.compile(
    r"""
    \[                       # 开始方括号
    ([^\]\[]*)               # 捕获内容（不含嵌套方括号）
    \]
    \s*                      # 可选空白
    \(                       # 参数列表开始
    """, re.VERBOSE
)


def parse_capture(cap_str: str) -> dict:
    """解析捕获字符串，如 '&', '=', '&x, y', '', '&, x'."""
    cap = cap_str.strip()
    if cap == "":
        return {"mode": "none", "default": None, "items": []}
    parts = [p.strip() for p in cap.split(",")]
    default = None
    items = []
    if parts[0] in ("&", "="):
        default = parts[0]
        parts = parts[1:]
    for p in parts:
        p = p.strip()
        if not p:
            continue
        if p.startswith("&"):
            items.append({"name": p[1:].strip(), "by_ref": True})
        elif p == "this":
            items.append({"name": "this", "by_ref": False, "is_this": True})
        elif p.startswith("="):
            items.append({"name": p[1:].strip(), "by_ref": False})
        else:
            # 默认按 default 捕获类型
            items.append({"name": p, "by_ref": default == "&" if default else False})
    mode = (
        "implicit-ref" if default == "&" and not items else
        "implicit-val" if default == "=" and not items else
        "explicit" if not default else
        f"{default}+explicit"
    )
    return {"mode": mode, "default": default, "items": items}


sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE as _TRANSLATION_SCOPE

# 扩展要识别的"宿主函数"集合：包含 TRANSLATION_SCOPE + factorize + 已知顶层入口
_KNOWN_HOSTS = set(_TRANSLATION_SCOPE) | {"factorize", "squarefreefactorize"}

# 每文件缓存"函数名 → 起始行号"
_FUNC_START_CACHE: dict[Path, dict[str, int]] = {}


def _index_function_starts(file_path: Path) -> dict[str, int]:
    """扫一遍文件，找每个已知函数名的首个定义行。

    匹配模式：函数名出现在 `fn_name\s*\(` 形式，且前一行（若是返回类型单独一行）或同一行
    以合理的返回类型开头（`inline` / `template` / `auto` / `void` / `bool` / `int` /
    `std::` / `ZZ` / `Zp` / `QQ` / `polynomial_` / `upolynomial_` / `factorization` / `__` / `std::pair`）。
    """
    if file_path in _FUNC_START_CACHE:
        return _FUNC_START_CACHE[file_path]

    lines = file_path.read_text().splitlines()
    starts: dict[str, int] = {}

    for i, line in enumerate(lines):
        # 检查所有已知函数名是否以 `fn_name(` 的形式出现
        for fn in _KNOWN_HOSTS:
            idx = line.find(fn + "(")
            if idx < 0:
                continue
            # 确认 fn 前是非 word char（避免匹配 `__fn_name_suffix`）
            if idx > 0 and (line[idx - 1].isalnum() or line[idx - 1] == "_"):
                continue
            # 确认这是"定义"而不是"调用"——简单判据：所在行或上一行含
            # 返回类型关键词（inline/auto/template 或大写开头的标识）
            # 或 line 开头就是 fn 名（前面全是空白，表示多行返回类型声明）
            line_stripped = line[:idx].strip()
            prev_line = lines[i - 1].strip() if i > 0 else ""

            def looks_like_return_type(s: str) -> bool:
                if not s:
                    return False
                # 以 `>` 结尾（模板结束）
                if s.endswith(">") or s.endswith("&"):
                    return True
                # 以 `inline` 开头
                if s.startswith("inline") or s.startswith("template"):
                    return True
                # 其他可能的返回类型片段
                return s[-1].isalpha() or s[-1] == "_"

            is_def = False
            if not line_stripped:
                # fn 名在行首（排除空白），上一行是返回类型
                if looks_like_return_type(prev_line):
                    is_def = True
            elif looks_like_return_type(line_stripped):
                # 返回类型和 fn 名在同一行
                is_def = True

            if is_def and fn not in starts:
                starts[fn] = i + 1  # 1-indexed

    _FUNC_START_CACHE[file_path] = starts
    return starts


def find_function_at(file_path: Path, target_line: int) -> str | None:
    """从已知函数名表中找包含 target_line 的函数名（最大起始行 ≤ target_line）。"""
    starts = _index_function_starts(file_path)
    best = None
    best_line = -1
    for fn, sl in starts.items():
        if sl <= target_line and sl > best_line:
            best = fn
            best_line = sl
    return best


def scan_source_file(file_path: Path) -> list[dict]:
    """扫源码找所有 lambda，返回列表。"""
    content = file_path.read_text()
    results = []
    offset = 0
    while True:
        m = LAMBDA_HEAD_RE.search(content, offset)
        if not m:
            break
        # 确认这是一个 lambda（不是数组访问等）：前面应该是表达式上下文
        # 简单启发式：前面不紧邻 `)` / `]` / identifier（那些意味着 operator[] 或索引）
        start_pos = m.start()
        # 找前一个非空白字符
        pre_idx = start_pos - 1
        while pre_idx >= 0 and content[pre_idx] in " \t\n":
            pre_idx -= 1
        if pre_idx >= 0:
            c = content[pre_idx]
            # 如果紧接 identifier char 或 )] 或 .>  则是下标/调用，不是 lambda
            if c.isalnum() or c == "_":
                offset = m.end()
                continue
        # 从 pre_idx 回推：如果是 "return ", "(", ",", "=", "&&", "||" 等则是 lambda 上下文
        # 但简单起见，不再细化，接受这里为 lambda
        cap_str = m.group(1)
        # 计算行号
        line_no = content[: start_pos].count("\n") + 1
        # 找对应的 `)` 然后 `{`（body 开始）
        body_start = content.find("{", m.end())
        if body_start < 0:
            offset = m.end()
            continue
        # body_end: 匹配大括号
        depth = 1
        j = body_start + 1
        while j < len(content) and depth > 0:
            if content[j] == "{": depth += 1
            elif content[j] == "}": depth -= 1
            j += 1
        body = content[body_start:j]
        body_line_count = body.count("\n") + 1

        # 参数字符串
        params_end = body_start - 1  # 略过空白
        while params_end > m.end() and content[params_end] in " \t\n>-":
            params_end -= 1
        # 参数 ( ... )
        # 从 m.end()-1 开始
        paren_start = m.end() - 1
        depth = 1
        k = paren_start + 1
        while k < len(content) and depth > 0:
            if content[k] == "(": depth += 1
            elif content[k] == ")": depth -= 1
            k += 1
        params_str = content[paren_start + 1: k - 1]

        # 是否 generic（含 auto）
        is_generic = bool(re.search(r"\bauto\b", params_str))

        results.append({
            "file": file_path.name,
            "line": line_no,
            "host_function": find_function_at(file_path, line_no),
            "capture_raw": cap_str,
            "capture_parsed": parse_capture(cap_str),
            "params": params_str.strip(),
            "is_generic": is_generic,
            "body_lines": body_line_count,
            "body_preview": body[:200].replace("\n", " | "),
        })

        offset = j
    return results


# --- 步骤 2: 按 TRANSLATION_SCOPE 过滤 ---
sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE


def main():
    all_lambdas = []
    for f in sorted(CLPOLY_DIR.glob("polynomial_factorize*.hh")):
        all_lambdas.extend(scan_source_file(f))

    # 过滤：host_function 在 TRANSLATION_SCOPE 内（或是 `factorize` 的 dispatcher 之类）
    in_scope = [l for l in all_lambdas if l["host_function"] in TRANSLATION_SCOPE]
    out_of_scope = [l for l in all_lambdas if l["host_function"] not in TRANSLATION_SCOPE]

    # 输出报告
    lines = []
    lines.append("# Lambda 扫描")
    lines.append("")
    lines.append(f"总计 **{len(all_lambdas)}** 个 lambda，其中 **{len(in_scope)}** 个在 TRANSLATION_SCOPE 内。")
    lines.append("")
    lines.append("> 源码正则扫描（Clang JSON 对 LambdaExpr 的 captures 字段省略，无法用 AST）。")
    lines.append("")

    # 全局统计
    lines.append("## 全局统计（仅 in-scope）")
    lines.append("")
    n_generic = sum(1 for l in in_scope if l["is_generic"])
    n_no_cap = sum(1 for l in in_scope if l["capture_parsed"]["mode"] == "none")
    n_implicit_ref = sum(1 for l in in_scope if l["capture_parsed"]["mode"] == "implicit-ref")
    n_implicit_val = sum(1 for l in in_scope if l["capture_parsed"]["mode"] == "implicit-val")
    n_explicit = sum(1 for l in in_scope if l["capture_parsed"]["mode"] == "explicit")
    n_mixed = sum(1 for l in in_scope if l["capture_parsed"]["mode"].endswith("+explicit"))

    lines.append(f"- Generic（`auto` 参数）: **{n_generic}** — 若仍有请再修 Clang 单态化问题")
    lines.append(f"- 捕获模式 `[]`（无捕获）: **{n_no_cap}**")
    lines.append(f"- 捕获模式 `[&]`（默认引用）: **{n_implicit_ref}**")
    lines.append(f"- 捕获模式 `[=]`（默认值）: **{n_implicit_val}**")
    lines.append(f"- 捕获模式 `[x, y, ...]`（具名）: **{n_explicit}**")
    lines.append(f"- 捕获模式 `[&, x, ...]` / `[=, &x, ...]`（混合）: **{n_mixed}**")
    lines.append("")

    # body 行数分布
    lines.append("## Body 行数直方图（in-scope）")
    lines.append("")
    lines_buckets = Counter()
    for l in in_scope:
        bl = l["body_lines"]
        if bl <= 1: lines_buckets["1 行"] += 1
        elif bl <= 3: lines_buckets["2-3 行"] += 1
        elif bl <= 10: lines_buckets["4-10 行"] += 1
        elif bl <= 30: lines_buckets["11-30 行"] += 1
        else: lines_buckets["30+ 行"] += 1
    for k in ["1 行", "2-3 行", "4-10 行", "11-30 行", "30+ 行"]:
        if k in lines_buckets:
            lines.append(f"- {k}: {lines_buckets[k]}")
    lines.append("")

    # 按宿主函数详情
    lines.append("## 按宿主函数详情")
    lines.append("")
    by_func = defaultdict(list)
    for l in in_scope:
        by_func[l["host_function"]].append(l)
    for fn in sorted(by_func.keys()):
        lambdas = by_func[fn]
        lines.append(f"### `{fn}` ({len(lambdas)} lambda{'s' if len(lambdas)>1 else ''})")
        lines.append("")
        for i, l in enumerate(lambdas, 1):
            cap = l["capture_parsed"]
            cap_desc = f"`[{l['capture_raw']}]`" if l['capture_raw'] else "`[]`"
            tags = [cap["mode"]]
            if l["is_generic"]:
                tags.append("GENERIC")
            lines.append(f"**Lambda {i}** — {l['file']}:{l['line']} — {cap_desc} — {' / '.join(tags)}")
            lines.append(f"- params: `{l['params'][:100]}`")
            lines.append(f"- body: {l['body_lines']} 行")
            if cap["items"]:
                caps_str = ", ".join(
                    f"{'&' if c['by_ref'] else '='}{c['name']}"
                    for c in cap["items"]
                )
                lines.append(f"- 具名捕获: {caps_str}")
            lines.append(f"- body preview: `{l['body_preview'][:150]}`")
            lines.append("")

    # out-of-scope 列表
    if out_of_scope:
        lines.append("## Out-of-scope lambda（宿主函数不在 TRANSLATION_SCOPE）")
        lines.append("")
        for l in out_of_scope:
            lines.append(f"- `{l['host_function']}` @ {l['file']}:{l['line']} — `[{l['capture_raw']}]`")
        lines.append("")

    OUT_MD.write_text("\n".join(lines))
    print(f"Written: {OUT_MD}")
    print(f"Total lambdas: {len(all_lambdas)}, in-scope: {len(in_scope)}")
    print(f"  Generic: {n_generic}")
    print(f"  Capture modes: no-cap={n_no_cap}, [&]={n_implicit_ref}, [=]={n_implicit_val}, explicit={n_explicit}, mixed={n_mixed}")


if __name__ == "__main__":
    main()
