"""Known-factor outer entry cases, including irreducibles without inert primes."""
import argparse
import json
from pathlib import Path

from generate_vanhoeij_entry_vectors import encode, multiply

OUT = Path(__file__).resolve().parent / "vectors/__factor_squarefree_primitive_ZZ_vanhoeij.json"


def generate():
    families = [
        ("two_quadratics", [[1, 0, 1], [2, 0, 1]]),
        ("three_quadratics", [[1, 0, 1], [2, 0, 1], [3, 0, 1]]),
        ("nonmonic_quadratics", [[1, 0, 2], [1, 0, 3]]),
        ("mixed_degrees", [[-1, 1], [2, 1], [1, 0, 1]]),
        ("phi8", [[1, 0, 0, 0, 1]]),
        ("phi12", [[1, 0, -1, 0, 1]]),
        ("phi24", [[1, 0, 0, 0, -1, 0, 0, 0, 1]]),
        # The second quartic is 16*Phi_8(x/2), so is also irreducible over Q.
        ("two_quartics", [[1, 0, 0, 0, 1], [16, 0, 0, 0, 1]]),
        ("large_height", [[1 + 17 * 2**48, 0, 1], [2 + 17 * 2**48, 0, 1]]),
    ]
    cases = []
    for name, factors in families:
        f = [1]
        for factor in factors:
            f = multiply(f, factor)
        for large_prime in [False, True]:
            # Both wrappers check the unordered known factors, then return the
            # unchanged physical order for the exact C++/Lean comparison.
            cases.append({"name": f"{name}_{'large' if large_prime else 'small'}_prime",
                          "args": [encode(f), {"type": "Bool", "val": large_prime},
                                   [encode(factor) for factor in factors]]})
    return {"fn": "__factor_squarefree_primitive_ZZ", "cases": cases}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    vectors = generate()
    content = json.dumps(vectors, indent=2) + "\n"
    if args.check:
        if not OUT.exists() or OUT.read_text() != content:
            raise SystemExit(f"stale vectors: {OUT}")
    else:
        OUT.write_text(content)
    print(f"{len(vectors['cases'])} outer entry cases")
