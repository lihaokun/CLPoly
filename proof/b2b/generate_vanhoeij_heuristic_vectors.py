"""Independent float-bit and integer-power oracles for the checked prefix."""
import argparse
import json
import math
import random
import struct
from pathlib import Path

from generate_vanhoeij_entry_vectors import encode

ROOT = Path(__file__).resolve().parent / "vectors"
LIMIT = 2**31 - 1


def bits(value):
    return struct.unpack("<Q", struct.pack("<d", value))[0]


def clamp(value):
    if math.isnan(value) or value >= LIMIT:
        return LIMIT
    return 1 if value <= 1 else max(1, int(value))


def exponent(prime, target):
    power, count = 1, 0
    while power <= target:
        power *= prime
        count += 1
    return count


def generate():
    values = [float("nan"), float("inf"), -float("inf"), -0.0, 0.0,
              -1.0, 0.5, 1.0, 1.5, 2.0, 17.9, float(LIMIT),
              math.nextafter(float(LIMIT), 0.0), float(LIMIT - 1),
              float(LIMIT + 1), float(2**32), -float(2**32)]
    cases = [{"args": [bits(v)], "expected": clamp(v)} for v in values]
    rng = random.Random(20261016)
    for _ in range(200):
        raw = rng.getrandbits(64)
        value = struct.unpack("<d", struct.pack("<Q", raw))[0]
        cases.append({"args": [raw], "expected": clamp(value)})
    out = {"__vanhoeij_clamp_precision": cases}
    cases = []
    for prime in [2, 3, 17, 97, 65537, 18446744073709551557]:
        targets = [-3, -1, 0, 1, 2, prime - 1, prime, prime + 1]
        for a in [2, 3, 8, 16, 32]:
            targets += [prime**a - 1, prime**a, prime**a + 1]
        for target in targets:
            cases.append({"args": [prime, {"type": "ZZ", "val": str(target)}],
                          "expected": exponent(prime, target)})
    out["__vanhoeij_full_precision"] = cases
    cases = []
    for poly in [[1, 1], [1, 0, 1], [3, -4, 5], [2**96 + 1, 0, -7, 0, 3]]:
        degree = len(poly) - 1
        square = sum(c*c for c in poly)
        norm = math.isqrt(square)
        norm += norm*norm < square
        bound = math.comb(degree, degree // 2) * norm
        target = 2 * abs(poly[-1]) * bound
        terms = sum(c != 0 for c in poly)
        for prime in [2, 17, 97, 18446744073709551557]:
            for count in [2, 24, 32, LIMIT]:
                logp = math.log(prime)
                value = math.ceil((2.5 * count + prime.bit_length()) * math.log(2.0) / logp
                                  + math.log(terms) / (2.0 * logp))
                full = exponent(prime, target)
                cases.append({"args": [encode(poly), count, prime],
                              "expected": [min(full, clamp(value)), full]})
    out["__heuristic_starting_precision"] = cases
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for name, cases in generate().items():
        path = ROOT / (name + ".json")
        content = json.dumps({"fn": name, "cases": cases}, indent=2) + "\n"
        if args.check:
            if not path.exists() or path.read_text() != content:
                raise SystemExit(f"stale vectors: {path}")
        else:
            path.write_text(content)
        print(f"{name}: {len(cases)} cases")
