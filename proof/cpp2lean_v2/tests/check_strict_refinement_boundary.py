#!/usr/bin/env python3
"""Reject legacy L2 dispatch from the strict L1→L2 refinement boundary."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FILES = [
    ROOT / "proof/lean/CLPoly/Pipeline/L1.lean",
    ROOT / "proof/lean/CLPoly/Refinement/Basic.lean",
    ROOT / "proof/lean/CLPoly/Refinement/ZpArith.lean",
    ROOT / "proof/lean/CLPoly/Refinement/SquarefreeZp.lean",
    ROOT / "proof/lean/CLPoly/Refinement/DDF.lean",
    ROOT / "proof/lean/CLPoly/Refinement/EDF.lean",
    ROOT / "proof/lean/CLPoly/Refinement/Hensel.lean",
    ROOT / "proof/lean/CLPoly/Refinement/Recombine.lean",
]

FORBIDDEN = (
    "import CLPoly.Generated.Corpus",
    "HasPolyDivmod.polyDivmod",
    "def henselGeneratedCandidate",
    "def henselCandidateToPk",
)


def main() -> None:
    failures: list[str] = []
    for path in FILES:
        source = path.read_text()
        for token in FORBIDDEN:
            if token in source:
                failures.append(f"{path.relative_to(ROOT)}: forbidden {token!r}")
    if failures:
        raise SystemExit("strict refinement boundary violation:\n" + "\n".join(failures))
    print("PASS: strict refinement boundary cannot import legacy L2 dispatch")


if __name__ == "__main__":
    main()
