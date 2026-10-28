"""Independent integer oracle for checked precision growth, including overflow."""
import argparse
import json
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent / "vectors/__vanhoeij_next_precision.json"
LIMIT = 2**31 - 1


def generate():
    points = [0, 1, 2, 3, 17, LIMIT // 2 - 1, LIMIT // 2,
              LIMIT // 2 + 1, LIMIT - 1, LIMIT]
    pairs = [(a, b) for a in points for b in points]
    rng = random.Random(20261016)
    pairs += [(rng.randrange(LIMIT + 1), rng.randrange(LIMIT + 1)) for _ in range(100)]
    cases = []
    for precision, full in pairs:
        if precision == LIMIT:
            expected = {"overflow": True}
        else:
            expected = {"value": min(LIMIT, max(full, 1, 2 * precision))}
        cases.append({"args": [precision, full], "expected": expected})
    return {"fn": "__vanhoeij_next_precision", "cases": cases}


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
    print(f"{len(vectors['cases'])} checked precision cases")
