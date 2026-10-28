#!/usr/bin/env python3
"""Per-function elaboration check for Generated/Corpus.lean.

把 lake build 的所有 errors 按函数归类：
1. 解析 Corpus.lean 抽取每个 `partial def NAME` 的行号范围（直到下一个 def 或文件末）
2. 解析 lake build error log（stdin or 给定文件）提取每个 error 行号
3. 对每个函数计算 error count
4. 输出 per-function 状态：通过 / 错数

用法：
  python3 per_function_check.py <corpus.lean> [<lake_errors.txt>]

如果不给 lake_errors.txt 就跑 lake build 自己抓输出。
"""

from __future__ import annotations

import re
import sys
import subprocess
from pathlib import Path
from collections import defaultdict


def parse_corpus_funcs(corpus_path: Path) -> list[tuple[str, int, int]]:
    """返回 [(name, start_line, end_line), ...]，行号 1-based 包含。"""
    text = corpus_path.read_text()
    lines = text.splitlines()
    funcs: list[tuple[str, int]] = []  # (name, line_no)
    pat = re.compile(r"^partial def (\S+)")
    for i, line in enumerate(lines, 1):
        m = pat.match(line)
        if m:
            funcs.append((m.group(1), i))
    # end line = next func's start - 1（最后一个 = 文件末）
    out: list[tuple[str, int, int]] = []
    for k, (name, start) in enumerate(funcs):
        end = funcs[k + 1][1] - 1 if k + 1 < len(funcs) else len(lines)
        out.append((name, start, end))
    return out


def parse_errors(error_log: str) -> list[tuple[int, str]]:
    """从 lake build 输出提取 (line_no, first_50_chars_of_msg)。"""
    out: list[tuple[int, str]] = []
    pat = re.compile(r"^error: CLPoly/Generated/Corpus\.lean:(\d+):\d+: (.+)$")
    for line in error_log.splitlines():
        m = pat.match(line)
        if m:
            line_no = int(m.group(1))
            msg = m.group(2)[:80]
            out.append((line_no, msg))
    return out


def assign_errors_to_funcs(
    funcs: list[tuple[str, int, int]],
    errors: list[tuple[int, str]],
) -> dict[str, list[tuple[int, str]]]:
    """每个 error 归到包含其行号的函数。返回 {func_name: [(line, msg), ...]}."""
    by_func: dict[str, list[tuple[int, str]]] = defaultdict(list)
    # funcs 已按 start 顺序；用二分或线性都行
    for line_no, msg in errors:
        for name, s, e in funcs:
            if s <= line_no <= e:
                by_func[name].append((line_no, msg))
                break
    return by_func


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    corpus_path = Path(sys.argv[1])
    if not corpus_path.exists():
        print(f"corpus not found: {corpus_path}", file=sys.stderr)
        sys.exit(1)

    if len(sys.argv) >= 3:
        error_log = Path(sys.argv[2]).read_text()
    else:
        # 跑 lake build 自己抓
        lean_dir = corpus_path.parent.parent.parent  # proof/lean
        proc = subprocess.run(
            ["lake", "build", "CLPoly.Generated.Corpus"],
            cwd=lean_dir, capture_output=True, text=True,
        )
        error_log = proc.stdout + "\n" + proc.stderr

    funcs = parse_corpus_funcs(corpus_path)
    errors = parse_errors(error_log)
    by_func = assign_errors_to_funcs(funcs, errors)

    # 报告
    n_funcs = len(funcs)
    n_passed = sum(1 for f in funcs if f[0] not in by_func)
    n_failed = n_funcs - n_passed

    print(f"=== Per-function 状态 ({corpus_path.name}) ===")
    print(f"总函数数：{n_funcs}")
    print(f"通过：{n_passed} ({100*n_passed/n_funcs:.1f}%)")
    print(f"失败：{n_failed}")
    print(f"总错误：{len(errors)}")
    print()

    # 失败函数按错数降序
    failed = [(name, s, e) for (name, s, e) in funcs if name in by_func]
    failed.sort(key=lambda f: -len(by_func[f[0]]))

    print("=== 失败函数（按错数降序）===")
    for name, s, e in failed:
        errs = by_func[name]
        print(f"  {name} (lines {s}-{e}): {len(errs)} errors")
        for ln, msg in errs[:3]:  # show top 3
            print(f"    L{ln}: {msg[:70]}")
        if len(errs) > 3:
            print(f"    ... +{len(errs)-3} more")
    print()

    # 按错数分布
    bucket_counts: dict[str, int] = defaultdict(int)
    for name, s, e in failed:
        n = len(by_func[name])
        if n == 1: bucket_counts["1"] += 1
        elif n <= 3: bucket_counts["2-3"] += 1
        elif n <= 5: bucket_counts["4-5"] += 1
        elif n <= 10: bucket_counts["6-10"] += 1
        else: bucket_counts[">10"] += 1
    print("=== 失败函数错数分布 ===")
    for bucket in ["1", "2-3", "4-5", "6-10", ">10"]:
        if bucket in bucket_counts:
            print(f"  {bucket} 错: {bucket_counts[bucket]} 函数")


if __name__ == "__main__":
    main()
