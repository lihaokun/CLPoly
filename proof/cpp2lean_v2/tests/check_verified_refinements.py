"""Check that the centralized proved refinement API matches Pass 9 output."""

from __future__ import annotations

import sys
from pathlib import Path


V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from class_map import REFINEMENT_MAP
from pass9_refine_gen import REFINEMENT_DIR, _emit_verified_contracts_file


def main() -> None:
    output = REFINEMENT_DIR / "Generated.lean"
    expected = _emit_verified_contracts_file(list(REFINEMENT_MAP.values()))
    if not output.exists() or output.read_text() != expected:
        raise SystemExit(f"generated verified refinement output is stale: {output}")


if __name__ == "__main__":
    main()
