"""Generate only the currently honest strict C++ L1 dependencies for DDF.

The DDF entry itself must not be emitted until its C++ ``polynomial_GCD``
dependency is translated as L1.  Mapping that call to ``CLPoly.Model`` would
silently turn the purported L1 proof into an L2 fallback.

The same restriction applies to ``__upoly_mod`` and ``__upoly_powmod``: their
current translation reaches ``pair_vec_div5`` and hence the hand-written L2
``SparsePolyZp.divmod`` instance.  They are deliberately excluded until the
raw four-loop division implementation is connected to their source call.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from build_pass8_corpus import generate_corpus
from pass8_codegen import codegen_corpus


OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictDDF.lean"
STRICT_DDF_ROOTS = {
    "__make_zp",
    "__upoly_subtract_x",
}


def generate_strict_ddf() -> str:
    _, skipped, roots = generate_corpus()
    if skipped:
        details = ", ".join(f"{name}:{reason}" for name, reason in skipped)
        raise RuntimeError(f"full MIR generation skipped functions: {details}")
    selected = [f for f in roots if f.base_name in STRICT_DDF_ROOTS]
    found = {f.base_name for f in selected}
    missing = STRICT_DDF_ROOTS - found
    if missing:
        raise RuntimeError(f"missing strict DDF roots: {sorted(missing)}")
    source = codegen_corpus(selected, namespace="Generated.StrictDDF")
    if "sorry" in source or "partial def" in source or "      default" in source:
        raise RuntimeError("strict DDF output contains an opaque placeholder")
    if "polynomial_GCD" in source:
        raise RuntimeError("strict DDF output contains the untranslated C++ GCD boundary")
    forbidden_dispatch = (
        "pair_vec_div", "HasPolyDivmod", "SparsePolyZp.divmod",
        "Array.get!", "Array.set!", "front!", "back!",
    )
    if any(token in source for token in forbidden_dispatch):
        raise RuntimeError("strict DDF output contains an L2 polynomial-division dispatch")
    return source


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_ddf()
    if args.check:
        if not args.output.exists() or args.output.read_text() != source:
            print(f"FAIL: {args.output} is not reproducible", file=sys.stderr)
            return 1
        print(f"PASS: {args.output} is reproducible and placeholder-free")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source)
    print(f"generated {args.output} ({source.count(chr(10)) + 1} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
