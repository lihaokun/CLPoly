#!/usr/bin/env python3
"""检查 `#print axioms` 审计输出只依赖 Lean/Mathlib 的标准公理。

用法：check_axiom_audit.py AUDIT.lean OUTPUT.txt

AUDIT.lean 中每条 `#print axioms NAME` 都必须在 OUTPUT.txt（`lake env lean AUDIT.lean`
的输出）里有对应报告，且依赖的公理只能是 propext / Classical.choice / Quot.sound。
出现 sorryAx、Lean.ofReduceBool（native_decide）等其他公理、缺少报告或有 error 时失败。
"""
import re
import sys

ALLOWED = {"propext", "Classical.choice", "Quot.sound"}

PRINT_RE = re.compile(r"^\s*#print\s+axioms\s+(\S+)", re.MULTILINE)
DEPENDS_RE = re.compile(r"'([^']+)' depends on axioms: \[([^\]]*)\]")
NONE_RE = re.compile(r"'([^']+)' does not depend on any axioms")
# Lean 的错误消息格式："<file>:<line>:<col>: error: ..." 或行首 "error: ..."
ERROR_RE = re.compile(r"(^|:\s?)error:", re.MULTILINE)


def normalize(name: str) -> str:
    return name.replace("«", "").replace("»", "")


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    audit_path, output_path = sys.argv[1], sys.argv[2]
    with open(audit_path, encoding="utf-8") as f:
        audited = [normalize(n) for n in PRINT_RE.findall(f.read())]
    with open(output_path, encoding="utf-8") as f:
        output = f.read()

    problems = []
    error_lines = [line for line in output.splitlines() if ERROR_RE.search(line)]
    if error_lines:
        problems.append("输出中有 error：\n" + "\n".join(error_lines)[:4000])

    reports = {}
    for name, axioms in DEPENDS_RE.findall(output):
        reports[normalize(name)] = {a.strip() for a in axioms.split(",") if a.strip()}
    for name in NONE_RE.findall(output):
        reports[normalize(name)] = set()

    missing = [n for n in audited if n not in reports]
    if missing:
        problems.append(f"{len(missing)} 条 #print axioms 没有对应报告，例如：" + ", ".join(missing[:5]))
    for name in audited:
        extra = reports.get(name, set()) - ALLOWED
        if extra:
            problems.append(f"{name} 依赖非标准公理：{', '.join(sorted(extra))}")

    if not audited:
        problems.append(f"{audit_path} 里没有 #print axioms")

    if problems:
        print(f"FAIL {audit_path}")
        for p in problems:
            print("  - " + p)
        return 1
    print(f"OK   {audit_path}: {len(audited)} 个定理，只依赖 {', '.join(sorted(ALLOWED))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
