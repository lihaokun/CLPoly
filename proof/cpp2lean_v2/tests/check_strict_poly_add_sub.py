#!/usr/bin/env python3
"""Bind strict raw add/sub lowering to the exact C++ method bodies."""

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
from check_strict_divrem import stable_ast

LEAN = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictPolyAddSub.lean"
REFINEMENT = (V2_ROOT.parent / "lean" / "CLPoly" / "Impl" /
              "StrictPolyAddSubRefinement.lean")

EXPECTED = {
    "nmod_add": "fafffbacebad32ddb1ce8b2a6b02e3d1b054bde99618f7cfdddfa8e3f5360b61",
    "nmod_sub": "f41724ba4b20bc8e3ee33221e2c3b9681aac98ceea2cf9de9b08afc94e67fb50",
    "_poly_add": "c530288b9a28f2bcfd1d2c50eb911eb139bae29067df7c152739fa1216c0d3ae",
    "_poly_sub": "850b0f46fb5b6368ae9c8f7c9f35e356c5cd0bdf8b8351179a16f97be174863e",
}


def digest(name: str) -> str:
    ast = dense_method_ast(name)
    return hashlib.sha256(json.dumps(
        stable_ast(ast), sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def main() -> None:
    actual = {name: digest(name) for name in EXPECTED}
    missing = [name for name, value in EXPECTED.items() if value == "TO_FILL"]
    if missing:
        print(json.dumps(actual, indent=2))
        raise SystemExit("record the reviewed AST hashes")
    for name, expected in EXPECTED.items():
        if actual[name] != expected:
            raise SystemExit(f"dense_upoly_zp::{name} AST drift")
    source = LEAN.read_text()
    forbidden = ("sorry", "partial def", "fuel", "oracle", "fallback",
                 "Generated.Corpus", "Array.get!", "Array.set!")
    found = [token for token in forbidden if token in source]
    if found:
        raise SystemExit(f"strict add/sub contains forbidden constructs: {found}")
    for fragment in ("addCommonLoop", "subCommonLoop", "subNegTailLoop",
                     "RawExec", "copyU64", "normaliseU64"):
        if fragment not in source:
            raise SystemExit(f"strict add/sub drift: missing {fragment}")
    refinement = REFINEMENT.read_text()
    found = [token for token in forbidden if token in refinement]
    if found:
        raise SystemExit(f"strict add/sub refinement contains forbidden constructs: {found}")
    for fragment in ("addCommonLoop_ok", "subCommonLoop_ok",
                     "subNegTailLoop_ok", "polyAdd_ok", "polySub_ok",
                     "nmodAdd_toNat", "nmodSub_toNat", "nmodAdd_lt",
                     "nmodSub_lt", "nmodAdd_cast", "nmodSub_cast",
                     "ExactOrDisjoint", "addCommonLoop_preserves_outside",
                     "subCommonLoop_preserves_outside",
                     "subNegTailLoop_preserves_outside",
                     "addCommonLoop_value", "subCommonLoop_value",
                     "nmodNeg", "nmodNeg_lt", "nmodNeg_cast",
                     "subNegTailLoop_value",
                     "sameAddress_eq_true_iff",
                     "region_ne_of_exactOrDisjoint_not_sameAddress",
                     "copyTail_preserves_prefix",
                     "addCommonLoop_preserves_input_tail",
                     "subCommonLoop_preserves_input_tail",
                     "addLeftLongTail", "addRightLongTail",
                     "polyAdd_equalLength_refines",
                     "polySub_equalLength_refines",
                     "polyAdd_leftLong_refines",
                     "polyAdd_rightLong_refines", "polyAdd_refines",
                     "RawDensePolyRep",
                     "RawHeap.SameLayout"):
        if fragment not in refinement:
            raise SystemExit(f"strict add/sub refinement drift: missing {fragment}")
    print("PASS: raw polynomial add/sub are pinned to exact C++ ASTs")


if __name__ == "__main__":
    main()
