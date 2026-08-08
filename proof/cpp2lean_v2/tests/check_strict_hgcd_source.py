#!/usr/bin/env python3
"""Bind the strict HGCD invariant to the exact C++ source family.

This gate does not claim raw HGCD refinement is complete.  It freezes the
four source bodies whose pointer/matrix operations must subsequently discharge
the determinant and transform obligations in ``StrictHGCDRefinement``.
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
from check_strict_divrem import stable_ast
from pass1_parse import parse_pass


EXPECTED = {
    "_mat_one": "82467677117ac0f62183eb2d4ff8879c7297787a0ba8a5523150027b6aff03dc",
    "_mat_row_update": "193d713664e9c1c8c08989e0c113a1c876a77be8ca8892ebe11cd972d03da8bf",
    "_hgcd_iter": "b4bdca04357c87b7f56b7ace1824b6bcd1c244b66f696638359134cce29a9a4d",
    "_hgcd_recursive": "a39dc9dc390042e7873087b671ff4721b1dd0ba1d9c7b25dbd4a12b6acf29191",
    "_gcd_hgcd": "7eabada4b3f368249f26de6945b34b8bfa37c0f5066a1aec5161feba072bdb9f",
    "_gcd_euclid": "c7bfe30f2d2da9ba4eb217092546166f10b830afdbaba44e3640424f71e810f2",
}

REQUIRED_CALLS = {
    "_mat_one": (),
    "_mat_row_update": ("_mul", "_poly_add"),
    "_hgcd_iter": ("_poly_divrem", "_mat_row_update"),
    "_hgcd_recursive": ("_hgcd_iter", "_hgcd_recursive", "_poly_divrem"),
    "_gcd_hgcd": ("_poly_divrem", "_hgcd_recursive", "_gcd_euclid"),
    "_gcd_euclid": ("_poly_divrem",),
}

LEAN_FILES = (
    V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictHGCD.lean",
    V2_ROOT.parent / "lean" / "CLPoly" / "Impl" / "StrictHGCDRawRefinement.lean",
    V2_ROOT.parent / "lean" / "CLPoly" / "Impl" / "StrictHGCDRefinement.lean",
)


def main() -> None:
    for name, expected in EXPECTED.items():
        ast = dense_method_ast(name)
        digest = hashlib.sha256(json.dumps(
            stable_ast(ast), sort_keys=True, separators=(",", ":")
        ).encode()).hexdigest()
        if digest != expected:
            raise SystemExit(
                f"dense_upoly_zp::{name} AST drift: expected {expected}, got {digest}")
        rendered = repr(parse_pass(ast).body)
        missing = [callee for callee in REQUIRED_CALLS[name]
                   if callee not in rendered]
        if missing:
            raise SystemExit(f"{name} call structure drift: missing {missing}")

    source = "\n".join(path.read_text() for path in LEAN_FILES)
    forbidden = ("sorry", "partial def", "fuel", "oracle", "fallback",
                 "Generated.Corpus", "Array.get!", "Array.set!")
    found = [token for token in forbidden if token in source]
    if found:
        raise SystemExit(f"strict HGCD invariant contains forbidden constructs: {found}")
    required = (
        "dense_upoly_zp__mat_one_ir",
        "matOne_refines",
        "matOne_preserves_prefix",
        "copyU64_refines_rawDense",
        "hgcdMatPolyRep_of_same_prefixes",
        "copyU64_preserves_hgcdMatPolyRep",
        "HgcdMatRawDenseRep",
        "hgcdMatRawDenseRep_of_same_prefixes",
        "copyU64_preserves_hgcdMatRawDenseRep",
        "rawDensePolyRep_zero_length",
        "rawDensePolyRep_one_of_read_one",
        "dense_upoly_zp__mat_row_update_ir",
        "HgcdIterState",
        "hgcdIterLoop",
        "hgcdIterInit",
        "hgcdIterInit_refines",
        "hgcdIterLoop_stop",
        "hgcdIterLoop_step_shape",
        "dense_upoly_zp__hgcd_iter_ir",
        "polyDivrem_remainder_lt",
        "matRowUpdate_zero_exec",
        "matRowUpdate_zero_refines",
        "matRowUpdate_nonzero_success_shape",
        "matRowUpdate_mul_result",
        "matRowUpdate_mul_preserves_entry1",
        "matRowUpdate_mul_preserves_guard",
        "matRowUpdate_add_preserves_entry0",
        "matRowUpdate_nonzero_sum_rep",
        "matRowUpdate_nonzero_refines",
        "matRowUpdate_refines",
        "matRowUpdate_preserves_guard",
        "normalize_gcd_eq_of_hgcd_transform",
        "HgcdTransform",
        "hgcdTransform_euclid_step",
        "hgcdRowUpdate_determinant",
        "normalize_gcd_eq_of_det_one_transform",
        "normalize_gcd_eq_of_det_neg_one_transform",
    )
    missing = [fragment for fragment in required if fragment not in source]
    if missing:
        raise SystemExit(f"strict HGCD invariant drift: missing {missing}")
    print("PASS: HGCD source family is pinned to determinant/GCD invariants")


if __name__ == "__main__":
    main()
