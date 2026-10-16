"""Boundary vectors with an independent integer/bit-length parameter oracle."""
import argparse
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "vectors"


def zz(value):
    return {"type": "ZZ", "val": str(value)}


def expected_parameters(count, height, modulus):
    # Closed-form least power of two, independent of the source doubling loop.
    required = max(1, (height + 31) // 32)
    base_exponent = (required - 1).bit_length()
    divisor = 1 << min(2 * base_exponent, max(base_exponent, modulus.bit_length() - count - 16))
    reduced = modulus // divisor
    residual = (height + divisor - 1) // divisor + 2 * count
    if modulus <= height or reduced <= 1 or reduced <= 8 * count * residual:
        # A non-null expected value makes the runner check rejection against
        # the independent oracle too (JSON null means "no oracle" there).
        return []
    return [zz(divisor), zz(reduced), zz(residual)]


def generate():
    cases = []
    seen = set()

    def add(count, height, modulus):
        if modulus <= 1 or (count, height, modulus) in seen:
            return
        seen.add((count, height, modulus))
        cases.append({"name": f"parameters_{len(cases)}", "args": [count, zz(height), zz(modulus)],
                      "expected": expected_parameters(count, height, modulus)})

    for count in [2, 3, 8, 24, 32, 127, 2147483647]:
        for height in [0, 1, 31, 32, 33, 63, 64, 65, (1 << 80) - 1, 1 << 80,
                       (1 << 80) + 1, (1 << 260) + 7]:
            divisor = 1 << (max(1, (height + 31) // 32) - 1).bit_length()
            threshold = 8 * count * ((height + divisor - 1) // divisor + 2 * count) * divisor
            for boundary in [height, divisor, 2 * divisor, threshold]:
                for offset in [-1, 0, 1]:
                    add(count, height, boundary + offset)
        for bits in [17, 18, 19, 39, 40, 41, 47, 48, 49, 64, 256, 512]:
            for offset in [-1, 0, 1]:
                add(count, 33, (1 << bits) + offset)
        # Cross the input-fixed divisor ceiling, including its transition.
        for base_exponent in [0, 1, 5, 20, 80]:
            height = 32 * (1 << base_exponent)
            for bits in [count + 16 + 2 * base_exponent - 1,
                         count + 16 + 2 * base_exponent,
                         count + 16 + 2 * base_exponent + 1]:
                if bits > 1024:
                    continue
                for offset in [-1, 0, 1]:
                    add(count, height, (1 << bits) + offset)
    rng = random.Random(20261015)
    for _ in range(200):
        add(rng.choice([2, 3, 24, 32, 127]), rng.getrandbits(rng.randrange(1, 300)),
            max(2, rng.getrandbits(rng.randrange(2, 600))))

    def poly(degrees):
        return {"type": "SparsePolyZZ", "val": [[d, "1"] for d in sorted(set(degrees), reverse=True)]}

    widths = []
    degree_lists = [[], [[]], [[], []], [[0]], [[], [0], []], [[1, 0], [4]],
                    [[2147483646]], [[], [2147483645, 0], [3], []]]
    for _ in range(100):
        degree_lists.append([rng.sample(range(500), rng.randrange(6)) for _ in range(rng.randrange(10))])
    for i, polynomials in enumerate(degree_lists):
        widths.append({"name": f"width_{i}", "args": [[poly(ds) for ds in polynomials]],
                       "expected": max([d + 1 for ds in polynomials for d in ds] + [0])})
    return {"__vanhoeij_choose_parameters": cases, "__vanhoeij_cld_width": widths}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for name, cases in generate().items():
        path = ROOT / f"{name}.json"
        content = json.dumps({"fn": name, "cases": cases}, ensure_ascii=False, indent=2) + "\n"
        if args.check:
            if not path.exists() or path.read_text() != content:
                raise SystemExit(f"stale vectors: {path}")
        else:
            path.write_text(content)
        print(f"{name}: {len(cases)} cases")
