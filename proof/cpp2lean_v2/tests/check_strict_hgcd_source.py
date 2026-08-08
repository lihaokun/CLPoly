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
    "_hgcd_iter": "b4bdca04357c87b7f56b7ace1824b6bcd1c244b66f696638359134cce29a9a4d",
    "_hgcd_recursive": "a39dc9dc390042e7873087b671ff4721b1dd0ba1d9c7b25dbd4a12b6acf29191",
    "_gcd_hgcd": "7eabada4b3f368249f26de6945b34b8bfa37c0f5066a1aec5161feba072bdb9f",
    "_gcd_euclid": "c7bfe30f2d2da9ba4eb217092546166f10b830afdbaba44e3640424f71e810f2",
}

REQUIRED_CALLS = {
    "_hgcd_iter": ("_poly_divrem", "_mat_row_update"),
    "_hgcd_recursive": ("_hgcd_iter", "_hgcd_recursive", "_poly_divrem"),
    "_gcd_hgcd": ("_poly_divrem", "_hgcd_recursive", "_gcd_euclid"),
    "_gcd_euclid": ("_poly_divrem",),
}

LEAN = (V2_ROOT.parent / "lean" / "CLPoly" / "Impl" /
        "StrictHGCDRefinement.lean")


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

    source = LEAN.read_text()
    forbidden = ("sorry", "partial def", "fuel", "oracle", "fallback",
                 "Generated.Corpus", "Array.get!", "Array.set!")
    found = [token for token in forbidden if token in source]
    if found:
        raise SystemExit(f"strict HGCD invariant contains forbidden constructs: {found}")
    required = (
        "normalize_gcd_eq_of_hgcd_transform",
        "normalize_gcd_eq_of_det_one_transform",
        "normalize_gcd_eq_of_det_neg_one_transform",
    )
    missing = [fragment for fragment in required if fragment not in source]
    if missing:
        raise SystemExit(f"strict HGCD invariant drift: missing {missing}")
    print("PASS: HGCD source family is pinned to determinant/GCD invariants")


if __name__ == "__main__":
    main()
