#!/usr/bin/env python3
"""enumerate_instances.py — 枚举 TRANSLATION_SCOPE 中每个函数的全部实例化。

通过一次性 dump 完整 TU AST，找所有 FunctionDecl，按名字分组，输出实例化表。
"""
from __future__ import annotations
import json
import sys
import subprocess
import glob
from pathlib import Path
from collections import defaultdict

HERE = Path(__file__).resolve().parent.parent
PROJECT_ROOT = HERE.parent.parent

sys.path.insert(0, str(HERE))
from class_map import TRANSLATION_SCOPE


def system_includes():
    try:
        r = subprocess.run(["g++", "-E", "-Wp,-v", "-x", "c++", "-"],
                           input="", capture_output=True, text=True)
        paths = []
        cap = False
        for line in r.stderr.splitlines():
            if "#include <...>" in line:
                cap = True; continue
            if line.startswith("End of search"):
                break
            if cap and line.startswith(" "):
                paths.append(f"-I{line.strip()}")
        for d in glob.glob("/usr/lib/llvm-*/lib/clang/*/include"):
            paths.append(f"-I{d}")
        return paths
    except Exception:
        return []


def dump_full_ast():
    cmd = ["clang++", "-std=c++17", f"-I{PROJECT_ROOT}"] + system_includes() + [
        "-Xclang", "-ast-dump=json", "-fsyntax-only",
        str(HERE / "instantiate.cc"),
    ]
    print(f"Dumping full AST (this takes a while, no filter)...", file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE, timeout=120)
    print(f"stdout size: {len(r.stdout)//1024//1024}MB", file=sys.stderr)
    return json.loads(r.stdout)


def walk(node, collect):
    if isinstance(node, dict):
        kind = node.get("kind")
        name = node.get("name", "")
        if kind == "FunctionDecl" and name in TRANSLATION_SCOPE:
            has_body = any(c.get("kind") == "CompoundStmt"
                           for c in node.get("inner", []) if isinstance(c, dict))
            collect[name].append({
                "mangled": node.get("mangledName"),
                "qualType": node.get("type", {}).get("qualType"),
                "has_body": has_body,
            })
        for c in node.get("inner", []):
            walk(c, collect)
    elif isinstance(node, list):
        for x in node:
            walk(x, collect)


def main():
    ast = dump_full_ast()
    collect = defaultdict(list)
    walk(ast, collect)

    print("# Instantiation enumeration")
    print()
    print(f"TRANSLATION_SCOPE size: {len(TRANSLATION_SCOPE)}")
    print(f"Functions found: {len(collect)}")
    print()
    print("## Functions with multiple instantiations")
    print()
    multi_inst = 0
    total_with_body = 0
    for name, instances in sorted(collect.items()):
        with_body = [i for i in instances if i["has_body"]]
        mangled_wb = [i for i in with_body if i["mangled"]]
        if len(mangled_wb) > 1:
            multi_inst += 1
            print(f"### `{name}` — {len(mangled_wb)} instantiations")
            for inst in mangled_wb:
                qt = (inst["qualType"] or "")[:120]
                print(f"- `{qt}`")
            print()
        total_with_body += 1 if with_body else 0

    print(f"\n## Summary")
    print(f"- Functions with exactly 1 mangled+body instance: {total_with_body - multi_inst}")
    print(f"- Functions with 2+ mangled+body instances: {multi_inst}")

    # 检查 TRANSLATION_SCOPE 中未找到任何实例化的
    missing = set(TRANSLATION_SCOPE) - set(collect.keys())
    no_body = [n for n, insts in collect.items()
               if not any(i["has_body"] for i in insts)]
    no_mangled = [n for n, insts in collect.items()
                  if any(i["has_body"] for i in insts) and
                  not any(i["has_body"] and i["mangled"] for i in insts)]
    if missing:
        print(f"\n## Missing from AST ({len(missing)})")
        for m in sorted(missing):
            print(f"- `{m}`")
    if no_body:
        print(f"\n## Found but no body ({len(no_body)})")
        for n in sorted(no_body):
            print(f"- `{n}`")
    if no_mangled:
        print(f"\n## Found with body but no mangled (= template def only, no instantiation)")
        for n in sorted(no_mangled):
            insts = collect[n]
            for i in insts:
                if i["has_body"]:
                    qt = (i["qualType"] or "")[:100]
                    print(f"- `{n}`: qt=`{qt}`")


if __name__ == "__main__":
    main()
