"""Generate the admitted strict C++ foundations used by polynomial_GCD.

Every definition in this file starts from a concrete Clang AST body in
``dense_upoly_zp.hh``.  Buffer-manipulating routines, divrem, Euclid/HGCD and
polynomial_GCD are added to ``METHOD_ROOTS`` only after they use explicit raw
heap semantics and their dependency closures are placeholder-free.  In
particular, an Array ``get!``/``set!`` model is not admissible: its default
out-of-bounds behaviour is not C++ semantics.
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
from ir_types import BinOp, Cast, IfStmt, UnaryOp, Var


OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictGCD.lean"
METHOD_ROOTS = (
    "inv_prime",
    "__preinvert_limb",
    "_umul128",
    "_add_carry3",
    "_lll_mod_preinv",
    "nmod_mul",
    "nmod_inv",
)

CONSTRUCTOR_ARITIES = {
    "dense_upoly_zp_default": 0,
    "dense_upoly_zp_of_prime": 1,
    "dense_upoly_zp_of_sparse": 2,
}


def dense_constructor_ast(name: str):
    record = dump_ast_json("dense_upoly_zp", timeout=120)
    if record is None or record.get("kind") != "CXXRecordDecl":
        raise RuntimeError("missing dense_upoly_zp record AST")
    arity = CONSTRUCTOR_ARITIES[name]
    matches = []
    for node in record.get("inner", []):
        if node.get("kind") != "CXXConstructorDecl":
            continue
        params = [x for x in node.get("inner", [])
                  if x.get("kind") == "ParmVarDecl"]
        param_types = [x.get("type", {}).get("qualType", "") for x in params]
        type_match = (
            name == "dense_upoly_zp_default" or
            name == "dense_upoly_zp_of_prime" and param_types == ["uint64_t"] or
            name == "dense_upoly_zp_of_sparse" and
              len(param_types) == 2 and "upolynomial_<Zp>" in param_types[0]
        )
        if len(params) == arity and type_match:
            matches.append(node)
    if len(matches) != 1:
        raise RuntimeError(f"ambiguous dense constructor {name}: {len(matches)}")
    return matches[0]


def dense_method_ast(name: str):
    """Select an exact method from the `dense_upoly_zp` record.

    This is required for names such as `gcd` that collide with free/template
    functions elsewhere in the translation unit.
    """
    record = dump_ast_json("dense_upoly_zp", timeout=120)
    if record is None or record.get("kind") != "CXXRecordDecl":
        raise RuntimeError("missing dense_upoly_zp record AST")
    matches = [node for node in record.get("inner", [])
               if node.get("kind") == "CXXMethodDecl" and
               node.get("name") == name and
               any(child.get("kind") == "CompoundStmt"
                   for child in node.get("inner", [])
                   if isinstance(child, dict))]
    if len(matches) != 1:
        raise RuntimeError(f"ambiguous dense method {name}: {len(matches)}")
    return matches[0]


def lower_method(name: str):
    ast = (dense_constructor_ast(name) if name in CONSTRUCTOR_ARITIES else
           dense_method_ast(name) if name == "gcd" else
           dump_ast_json(name, timeout=120))
    expected_kind = ("CXXConstructorDecl" if name in CONSTRUCTOR_ARITIES else
                     "FunctionDecl" if name == "inv_prime" else "CXXMethodDecl")
    if ast is None or ast.get("kind") != expected_kind:
        raise RuntimeError(f"missing concrete C++ AST for {name}: {expected_kind}")
    hir = parse_pass(ast)
    if name == "divrem":
        hir.body = _specialize_divrem_nonalias(hir.body)
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
    mir.base_name = (
        "dense_upoly_zp_divrem_nonalias" if name == "divrem" else
        name if name == "inv_prime" or name in CONSTRUCTOR_ARITIES else
        f"dense_upoly_zp_{name}")
    return mir


def _strip_casts(expr):
    while isinstance(expr, Cast):
        expr = expr.expr
    return expr


def _address_pair(expr):
    expr = _strip_casts(expr)
    if not isinstance(expr, BinOp) or expr.op not in ("==", "="):
        return None
    names = []
    for side in (expr.lhs, expr.rhs):
        side = _strip_casts(side)
        if not isinstance(side, UnaryOp) or side.op != "&":
            return None
        operand = _strip_casts(side.operand)
        if not isinstance(operand, Var):
            return None
        names.append(operand.name)
    return tuple(names)


def _flatten_alias_disjunction(expr):
    expr = _strip_casts(expr)
    if isinstance(expr, BinOp) and expr.op == "||":
        return (_flatten_alias_disjunction(expr.lhs) +
                _flatten_alias_disjunction(expr.rhs))
    pair = _address_pair(expr)
    return [pair] if pair is not None else []


def _specialize_divrem_nonalias(body):
    """Remove only the verified alias-protection branch.

    The sole production caller is `gcd` with four distinct local objects
    `q`, `r`, `a`, and `b`.  The later GCD closure gate checks those declaration
    identities before this specialization is admitted.  Here we additionally
    require the source guard to be exactly the four documented comparisons;
    any source drift fails generation instead of silently pruning a branch.
    """
    expected = {("Q", "A"), ("Q", "B"), ("R", "A"), ("R", "B")}
    matches = []
    for index, stmt in enumerate(body):
        if isinstance(stmt, IfStmt):
            pairs = set(_flatten_alias_disjunction(stmt.cond))
            if pairs == expected:
                matches.append(index)
    if len(matches) != 1:
        raise RuntimeError(
            f"divrem nonalias specialization guard mismatch: {matches}")
    return [stmt for index, stmt in enumerate(body) if index != matches[0]]


def generate_strict_gcd() -> str:
    source = codegen_corpus(
        [lower_method(name) for name in METHOD_ROOTS],
        namespace="Generated.StrictGCD")
    forbidden = (
        "sorry", "partial def", "unresolved call", "Array.get!", "Array.set!",
        "get!", "set!", " default", "fuel", "Safe", "oracle", "fallback",
        "axiom",
    )
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
