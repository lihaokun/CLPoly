"""Generate the strict C++ member-function core used by polynomial_GCD.

Every definition in this file starts from a concrete Clang AST body in
``dense_upoly_zp.hh``.  The outer constructors, divrem, Euclid/HGCD branch and
polynomial_GCD wrappers are added to ``METHOD_ROOTS`` only after their own
dependency closures are placeholder-free.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from scripts.survey_ast import dump_ast_json
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass2b_callsite_ref_elim import callsite_ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass3b_lambda_ref_elim import lambda_ref_elim_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import ssa_build_pass
from pass7_loop_lower import loop_lower_pass
from pass8_codegen import codegen_corpus


OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictGCD.lean"
METHOD_ROOTS = (
    "inv_prime",
    "deg",
    "lead",
    "nmod_mul",
    "nmod_inv",
    "scalar_mul",
    "to_upoly",
)


def lower_method(name: str):
    ast = dump_ast_json(name, timeout=120)
    expected_kind = "FunctionDecl" if name == "inv_prime" else "CXXMethodDecl"
    if ast is None or ast.get("kind") != expected_kind:
        raise RuntimeError(f"missing concrete C++ AST for {name}: {expected_kind}")
    hir = parse_pass(ast)
    hir = ref_elim_pass(hir)
    hir = callsite_ref_elim_pass(hir)
    hir = lambda_lift_pass(hir)
    hir = lambda_ref_elim_pass(hir)
    hir = iter_recognize_pass(hir)
    hir, gaps = operator_resolve_pass(hir)
    if any((gaps.unresolved_op, gaps.cast_miss, gaps.method_miss,
            gaps.op_miss, gaps.constructor_miss)):
        raise RuntimeError(f"unresolved {name} translation gaps: {gaps}")
    mir = loop_lower_pass(ssa_build_pass(hir))
    mir.base_name = name if name == "inv_prime" else f"dense_upoly_zp_{name}"
    return mir


def generate_strict_gcd() -> str:
    source = codegen_corpus(
        [lower_method(name) for name in METHOD_ROOTS],
        namespace="Generated.StrictGCD")
    forbidden = ("sorry", "partial def", "unresolved call")
    found = [token for token in forbidden if token in source]
    if found:
        raise RuntimeError(f"strict GCD core contains placeholders: {found}")
    return source


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_gcd()
    if args.check:
        if not args.output.exists() or args.output.read_text() != source:
            print(f"FAIL: {args.output} is not reproducible", file=sys.stderr)
            return 1
        print(f"PASS: {args.output} is reproducible and placeholder-free")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source)
    print(f"generated {args.output} ({source.count(chr(10)) + 1} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
