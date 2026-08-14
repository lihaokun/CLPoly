#!/usr/bin/env python3
"""Reject legacy L2 dispatch from the strict L1→L2 refinement boundary."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FILES = [ROOT / "proof/lean/CLPoly/Pipeline/L1.lean"]
FILES.extend(sorted((ROOT / "proof/lean/CLPoly/Refinement").glob("*.lean")))
FILES.extend(sorted((ROOT / "proof/lean/CLPoly/Generated").glob("Strict*.lean")))
FILES.extend(sorted((ROOT / "proof/lean/CLPoly/Impl").glob("Strict*.lean")))

FORBIDDEN = (
    "import CLPoly.Generated.Corpus",
    "HasPolyDivmod.polyDivmod",
    "partial def",
    "sorry",
    "fuel",
    "def henselGeneratedCandidate",
    "def henselCandidateToPk",
    "multiplyNormalizeMod :",
    "ops.multiplyNormalizeMod",
    "  zassenhaus :",
    "ops.zassenhaus ",
    "  prepareCandidates :",
)


def strip_lean_comments(source: str) -> str:
    """Remove nested Lean block comments and end-of-line comments."""
    output: list[str] = []
    index = 0
    depth = 0
    while index < len(source):
        pair = source[index:index + 2]
        if pair == "/-":
            depth += 1
            index += 2
        elif pair == "-/" and depth:
            depth -= 1
            index += 2
        elif pair == "--" and depth == 0:
            newline = source.find("\n", index)
            index = len(source) if newline < 0 else newline
        else:
            if depth == 0:
                output.append(source[index])
            index += 1
    return "".join(output)


def main() -> None:
    failures: list[str] = []
    for path in FILES:
        source = path.read_text()
        # Policy tokens in explanatory comments are harmless; audit the Lean
        # declarations themselves so phrases such as "without a counter" do
        # not make the check fail merely for documenting the guarantee.
        code = strip_lean_comments(source)
        for token in FORBIDDEN:
            if token in code:
                failures.append(f"{path.relative_to(ROOT)}: forbidden {token!r}")
    if failures:
        raise SystemExit("strict refinement boundary violation:\n" + "\n".join(failures))
    print("PASS: strict refinement boundary cannot import legacy L2 dispatch")


if __name__ == "__main__":
    main()
