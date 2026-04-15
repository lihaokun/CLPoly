#!/usr/bin/env python3
"""
CLPoly C++ → Lean IR 翻译器

用法:
  clang++ -Xclang -ast-dump=json -fsyntax-only file.cpp 2>/dev/null | python3 cpp2lean.py
  # 或
  python3 cpp2lean.py < ast.json
"""

import sys
from clang_ast import parse_translation_unit
from ssa_transform import transform_func
from lean_codegen import generate_file, gen_ssa_func


def main():
    use_ssa = "--ssa" in sys.argv

    json_str = sys.stdin.read()
    funcs = parse_translation_unit(json_str)

    # 统计
    unknown_count = 0
    for f in funcs:
        unknown_count += count_unknowns(f.body)

    if use_ssa:
        from lean_codegen import PRELUDE, _uses_uint128
        from ub_collector import collect_all, inject_obligations
        ssa_funcs = [transform_func(f) for f in funcs]
        # UB 证明目标注入到函数签名
        for sf in ssa_funcs:
            obs = collect_all(sf)
            inject_obligations(sf, obs)
        parts = ["-- Auto-generated Lean IR from C++ via cpp2lean (SSA mode)", ""]
        if any(_uses_uint128(f) for f in funcs):
            parts.append(PRELUDE)
        for sf in ssa_funcs:
            parts.append(gen_ssa_func(sf))
            parts.append("")
        lean_code = "\n".join(parts)
    else:
        lean_code = generate_file(funcs)
    print(lean_code)

    # 报告
    print(f"\n-- Translation summary:", file=sys.stderr)
    print(f"--   Functions: {len(funcs)}", file=sys.stderr)
    for f in funcs:
        print(f"--     {f.name} ({len(f.params)} params, {len(f.body)} stmts)", file=sys.stderr)
    if unknown_count > 0:
        print(f"--   Unknown nodes: {unknown_count} (marked with sorry)", file=sys.stderr)


def count_unknowns(stmts) -> int:
    count = 0
    for s in stmts:
        if hasattr(s, 'kind') and hasattr(s, 'children'):
            if s.kind not in ("CompoundStmt",):
                count += 1
            count += count_unknowns(getattr(s, 'children', []))
        for attr in ('then_body', 'else_body', 'body'):
            if hasattr(s, attr):
                count += count_unknowns(getattr(s, attr))
        if hasattr(s, 'value') and hasattr(s.value, 'kind') and isinstance(s.value, type) and hasattr(s.value, 'children'):
            count += 1
    return count


if __name__ == "__main__":
    main()
