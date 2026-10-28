"""Valid split Hensel inputs, generated independently by lifting simple roots."""
import argparse
import json
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent / "vectors/__vanhoeij_try_recombine.json"


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def evaluate(poly, root):
    result = 0
    for coefficient in reversed(poly):
        result = result * root + coefficient
    return result


def encode(poly):
    return {"type": "SparsePolyZZ", "val": [[i, str(poly[i])] for i in reversed(range(len(poly))) if poly[i]]}


def generate():
    # Positive quadratics have no rational root. The quartic is Phi_8.
    families = [
        ("two_quadratics", [[1, 0, 1], [2, 0, 1]], 17, False),
        ("irreducible_quartic", [[1, 0, 0, 0, 1]], 17, False),
        ("irreducible_quadratic", [[1, 0, 1]], 17, False),
        ("three_quadratics", [[1, 0, 1], [2, 0, 1], [3, 0, 1]], 97, False),
        ("nonmonic_quadratics", [[1, 0, 2], [1, 0, 3]], 97, False),
        ("linear_fast_path", [[-1, 1], [2, 1], [3, 2]], 17, True),
    ]
    cases = []
    rng = random.Random(20261016)
    for name, factors, prime, singleton in families:
        f = [1]
        for factor in factors:
            f = multiply(f, factor)
        derivative = [i * f[i] for i in range(1, len(f))]
        roots = [r for r in range(prime) if evaluate(f, r) % prime == 0]
        assert len(roots) == len(f) - 1
        assert all(evaluate(derivative, r) % prime for r in roots)
        for exponent in ([8, 16] if singleton else [1, 8, 16]):
            modulus = prime ** exponent
            lifted_roots = []
            for root in roots:
                precision = prime
                for _ in range(1, exponent):
                    assert evaluate(f, root) % precision == 0
                    digit = (-evaluate(f, root) // precision *
                             pow(evaluate(derivative, root), -1, prime)) % prime
                    root += precision * digit
                    precision *= prime
                assert evaluate(f, root) % modulus == 0
                lifted_roots.append(root)
            product = [f[-1]]
            for root in lifted_roots:
                product = multiply(product, [-root, 1])
            assert all((a - b) % modulus == 0 for a, b in zip(product, f))
            complete = singleton or exponent > 1
            expected_factors = [encode(factor) for factor in factors] if complete else []
            for ordering in range(3):
                if ordering == 1:
                    lifted_roots.reverse()
                elif ordering == 2:
                    rng.shuffle(lifted_roots)
                lifted = [encode([(-root) % modulus, 1]) for root in lifted_roots]
                case = {"name": f"{name}_p{prime}_a{exponent}_order{ordering}",
                        "args": [encode(f), lifted, {"type": "ZZ", "val": str(modulus)},
                                 complete, expected_factors]}
                # The wrapper checks known factors/completion independently;
                # B2B additionally compares exact lattice statistics.
                if singleton or not complete:
                    ordered = sorted(expected_factors, key=lambda x: json.dumps(x, separators=(",", ":")))
                    case["expected"] = [complete, ordered, 0, 0]
                cases.append(case)
    return {"fn": "__vanhoeij_try_recombine", "cases": cases}


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
    print(f"{len(vectors['cases'])} complete-entry cases")
