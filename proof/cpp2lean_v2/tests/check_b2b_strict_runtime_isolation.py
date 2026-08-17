#!/usr/bin/env python3
"""Keep unsafe physical B2B witnesses out of the formal proof boundary."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BOUNDARY_ROOTS = (
    ROOT / "proof/lean/CLPoly/Generated",
    ROOT / "proof/lean/CLPoly/Refinement",
    ROOT / "proof/lean/CLPoly/Pipeline",
)
FORBIDDEN = (
    "B2B.StrictRuntime",
    "StrictRuntime.erasedValue",
    "import B2B",
)


def main() -> None:
    failures: list[str] = []
    for boundary in BOUNDARY_ROOTS:
        for path in sorted(boundary.rglob("*.lean")):
            source = path.read_text()
            for token in FORBIDDEN:
                if token in source:
                    failures.append(
                        f"{path.relative_to(ROOT)}: forbidden B2B dependency {token!r}"
                    )
    if failures:
        raise SystemExit("B2B strict-runtime isolation violation:\n" + "\n".join(failures))
    print("PASS: B2B physical runtime cannot flow into Generated/Refinement/Pipeline")


if __name__ == "__main__":
    main()
