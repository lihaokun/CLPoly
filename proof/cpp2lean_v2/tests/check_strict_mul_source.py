#!/usr/bin/env python3
"""Bind strict raw multiplication work to the exact C++ source family."""

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

LEAN = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictMul.lean"
REFINEMENT = (V2_ROOT.parent / "lean" / "CLPoly" / "Impl" /
              "StrictMulRefinement.lean")

EXPECTED = {
    "_classical_mul": "0d701d8406109ef552a49a4538ed63c1776f63ae03bf3156e3f0083b76a47f06",
    "_kar_mul": "9d7b5c319affbb66f8b53edd73ea0d4c3ee73f8d68c7248e978cead304183d51",
    "_mul": "01ffd3f114721110d16e48864c9ae390bd17fbe1c1984184ddc0825e0bb7b24d",
}


def digest(name: str) -> str:
    return hashlib.sha256(json.dumps(
        stable_ast(dense_method_ast(name)), sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def main() -> None:
    actual = {name: digest(name) for name in EXPECTED}
    if any(value == "TO_FILL" for value in EXPECTED.values()):
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
        raise SystemExit(f"strict multiplication contains forbidden constructs: {found}")
    for fragment in ("classicalDotLoop", "classicalOuterLoop",
                     "dense_upoly_zp__classical_mul_ir", "RawExec",
                     "karAddHalvesLoop", "karSubLoop", "karAssembleLoop",
                     "karOddTail",
                     "karPrepareHalves",
                     "dense_upoly_zp__kar_mul_ir", "copyU64",
                     "readU64", "writeU64", "termination_by"):
        if fragment not in source:
            raise SystemExit(f"strict multiplication drift: missing {fragment}")
    refinement = REFINEMENT.read_text()
    found = [token for token in forbidden if token in refinement]
    if found:
        raise SystemExit(f"strict multiplication refinement contains forbidden constructs: {found}")
    for fragment in ("classicalDotLoop_ok", "classical_index_bounds",
                     "slicePolyRep_prefix_exists",
                     "karAddHalvesLoop_ok", "karSubLoop_ok",
                     "karAssembleLoop_ok",
                     "karOddTail_ok", "karOddTail_preserves_region_ne",
                     "karPrepareHalves_ok",
                     "karOddTail_values", "karOddTail_coeffs",
                     "karOddTail_preserves_own_prefixes",
                     "karOddTail_refines_slices_odd",
                     "karAddHalvesLoop_preserves_outside",
                     "karAddHalvesLoop_current_values",
                     "karAddHalvesLoop_current_coeffs",
                     "karAddHalvesLoop_coeffs",
                     "karHalfSumPoly", "coeff_karHalfSumPoly",
                     "karAddHalvesLoop_refines_slices",
                     "karSubLoop_preserves_outside",
                     "karAssembleLoop_preserves_outside",
                     "classicalDotNat", "classicalDotNat_ok",
                     "classicalDotNat_bound", "classicalDotLoop_modEq",
                     "classicalDotLoop_raw_sum",
                     "classicalDotLoop_exact_zero",
                     "classicalDotReduced_toNat",
                     "classicalDotReduced_cast",
                     "classicalDotPoly",
                     "classicalDotPoly_eq_sum_Icc",
                     "classicalDotPoly_source_eq_coeff",
                     "classicalDotNat_cast_eq_poly",
                     "classicalReduced_source_eq_coeff",
                     "classicalOuterLoop_preserves_outside",
                     "classicalOuterLoop_same_prefix_region_ne",
                     "canonicalU64Prefix_of_same_prefix",
                     "ClassicalCoeffPrefix",
                     "slicePolyRep_of_classicalCoeffPrefix",
                     "canonicalU64Prefix_of_classicalCoeffPrefix",
                     "normaliseU64_eq_length_of_classicalCoeffPrefix",
                     "mul_coeff_zero_of_slice_lengths",
                     "mul_last_coeff_ne_zero_of_rawDense",
                     "classicalCoeffPrefix_succ_of_write",
                     "classicalOuterLoop_preserves_coeff_prefix",
                     "classicalOuterLoop_refines_coeff_prefix",
                     "classicalOuterLoop_ok", "classicalMul_ok",
                     "classicalMul_refines_slice",
                     "classicalMul_refines",
                     "RawHeap.SameLayout"):
        if fragment not in refinement:
            raise SystemExit(f"strict multiplication refinement drift: missing {fragment}")
    print("PASS: multiplication source family is pinned and schoolbook raw lowering is strict")


if __name__ == "__main__":
    main()
