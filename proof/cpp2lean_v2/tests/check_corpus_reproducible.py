"""Check that the checked-in Lean corpus is reproducible from C++ AST/MIR.

This deliberately writes only to a temporary directory.  A non-zero diff is a
strict L1 provenance failure: content living in Generated/Corpus.lean is not by
itself evidence that cpp2lean emitted it.
"""

from __future__ import annotations

import difflib
import sys
import tempfile
from pathlib import Path

from build_pass8_corpus import OUT, generate_corpus


def main() -> int:
    generated, skipped, _ = generate_corpus()
    if skipped:
        for name, reason in skipped:
            print(f"SKIPPED {name}: {reason}", file=sys.stderr)
        return 2

    checked_in = OUT.read_text()
    if checked_in == generated:
        print("PASS: Generated/Corpus.lean is reproducible (zero diff)")
        return 0

    with tempfile.TemporaryDirectory(prefix="clpoly-corpus-repro-") as tmp:
        candidate = Path(tmp) / "Corpus.lean"
        candidate.write_text(generated)
        diff = difflib.unified_diff(
            checked_in.splitlines(), generated.splitlines(),
            fromfile=str(OUT), tofile=str(candidate), lineterm="")
        print("FAIL: Generated/Corpus.lean is not reproducible", file=sys.stderr)
        for line in list(diff)[:200]:
            print(line, file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
