"""Admission gate for the raw-heap `_poly_divrem` specialization.

The generic backend cannot yet thread `RawHeap` through arbitrary pointer
expressions.  Until that pass replaces this specialization, this gate binds
the Lean control-flow lowering to the exact Clang AST and rejects source or
artifact drift.  It is deliberately not presented as a completed refinement
proof: raw-to-safe invariants remain a separate obligation.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from build_strict_gcd import dense_method_ast
from ir_types import ForStmt, PairType, PtrType
from pass1_parse import parse_pass


OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictDivrem.lean"
REFINEMENT_ROOT = V2_ROOT.parent / "lean" / "CLPoly" / "Refinement"
AST_SHA256 = "6ed59f9a349887be5e6c3bc49f103a64d7df8d34ad0730f7a727ec9cb54f5c15"
LEAN_SHA256 = "2a555b0c01958ebb35d2b4518ecde4ea1e7bed847daf91d63a6ee67bbc154cf0"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_ast(node):
    """Remove Clang process-local declaration identities before hashing."""
    if isinstance(node, dict):
        return {key: stable_ast(value) for key, value in node.items()
                if key != "id" and not key.lower().endswith("id") and
                key != "referencedMemberDecl"}
    if isinstance(node, list):
        return [stable_ast(value) for value in node]
    return node


def check_ast() -> None:
    ast = dense_method_ast("_poly_divrem")
    actual = digest(json.dumps(
        stable_ast(ast), sort_keys=True, separators=(",", ":")).encode())
    if actual != AST_SHA256:
        raise RuntimeError(
            f"_poly_divrem Clang AST drift: expected {AST_SHA256}, got {actual}")

    hir = parse_pass(ast)
    if len(hir.params) != 8 or not isinstance(hir.ret_ty, PairType):
        raise RuntimeError("_poly_divrem signature no longer matches specialization")
    pointer_positions = (1, 2, 3, 5, 7)
    if any(not isinstance(hir.params[i].ty, PtrType) for i in pointer_positions):
        raise RuntimeError("_poly_divrem raw-pointer parameters were erased")

    def count_loops(stmts) -> int:
        total = 0
        for stmt in stmts:
            if isinstance(stmt, ForStmt):
                total += 1 + count_loops(stmt.body)
            for attr in ("then_body", "else_body", "body"):
                child = getattr(stmt, attr, None)
                if isinstance(child, list) and not isinstance(stmt, ForStmt):
                    total += count_loops(child)
        return total

    if count_loops(hir.body) != 4:
        raise RuntimeError("_poly_divrem no longer has the four admitted source loops")


def main() -> int:
    check_ast()
    source = OUT.read_bytes()
    actual = digest(source)
    if actual != LEAN_SHA256:
        raise RuntimeError(
            f"StrictDivrem lowering drift: expected {LEAN_SHA256}, got {actual}")
    forbidden = (
        b"sorry", b"partial def", b"fuel", b"Safe", b"oracle", b"fallback",
        b"axiom", b"get!", b"set!", b" default",
    )
    found = [token.decode() for token in forbidden if token in source]
    if found:
        raise RuntimeError(f"StrictDivrem contains forbidden constructs: {found}")
    withdrawn = {
        "ZpArith.lean": (b"__upoly_divmod_ir_refines", b"__upoly_mod_ir_refines"),
        "DDF.lean": (b"strict_upoly_powmod_refines", b"strict_upoly_subtract_x_refines"),
        "SquarefreeZp.lean": (b"squarefree_Zp_ir_refines", b"polynomial_GCD_refines"),
    }
    for name, claims in withdrawn.items():
        legacy_source = (REFINEMENT_ROOT / name).read_bytes()
        present = [claim.decode() for claim in claims if claim in legacy_source]
        if present:
            raise RuntimeError(
                f"{name} re-exports dispatch-based refinement claims: {present}")
    print("PASS: StrictDivrem is bound to the exact four-loop C++ AST and has no fallback")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
