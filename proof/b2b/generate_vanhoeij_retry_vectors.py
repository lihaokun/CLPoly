"""Full retry inputs with independently known irreducible factors."""
import argparse
import json
import random
from pathlib import Path

from generate_vanhoeij_entry_vectors import encode, evaluate, multiply

OUT = Path(__file__).resolve().parent / "vectors/__lll_factorize.json"


def encode_mod(poly, prime):
    return {"type": "SparsePolyZp", "val": [
        [i, [str(poly[i] % prime), str(prime)]]
        for i in reversed(range(len(poly))) if poly[i] % prime]}


def generate():
    families = [
        ("two_quadratics", [[1, 0, 1], [2, 0, 1]], 17),
        ("irreducible_quartic", [[1, 0, 0, 0, 1]], 17),
        ("irreducible_quadratic", [[1, 0, 1]], 17),
        ("three_quadratics", [[1, 0, 1], [2, 0, 1], [3, 0, 1]], 97),
        ("nonmonic_quadratics", [[1, 0, 2], [1, 0, 3]], 97),
        ("linear_fast_path", [[-1, 1], [2, 1], [3, 2]], 17),
    ]
    # Positive quadratics are irreducible over Q; perturbing their constants
    # by multiples of p preserves the distinct simple modular roots.
    for bits in [16, 48, 96]:
        families.append((f"large_height_{bits}",
                         [[1 + 17 * 2**bits, 0, 1], [2 + 17 * 2**bits, 0, 1]], 17))
    rng = random.Random(20261016)
    cases = []
    for name, factors, prime in families:
        f = [1]
        for factor in factors:
            f = multiply(f, factor)
        derivative = [i * f[i] for i in range(1, len(f))]
        roots = [r for r in range(prime) if evaluate(f, r) % prime == 0]
        assert len(roots) == len(f) - 1
        assert all(evaluate(derivative, root) % prime for root in roots)
        modular = [[-root, 1] for root in roots]
        product = [f[-1]]
        for factor in modular:
            product = multiply(product, factor)
        assert all((a - b) % prime == 0 for a, b in zip(product, f))
        expected = sorted([encode(factor) for factor in factors],
                          key=lambda value: json.dumps(value, separators=(",", ":")))
        for ordering in range(3):
            if ordering == 1:
                modular.reverse()
            elif ordering == 2:
                rng.shuffle(modular)
            cases.append({"name": f"{name}_p{prime}_order{ordering}",
                          "args": [encode(f), [encode_mod(poly, prime) for poly in modular], prime],
                          "expected": expected})
    return {"fn": "__lll_factorize", "cases": cases}


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
    print(f"{len(vectors['cases'])} full retry cases")
