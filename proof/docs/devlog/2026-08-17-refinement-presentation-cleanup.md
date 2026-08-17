# Refinement presentation cleanup

## Goal

Make the completed strict C++ L1 to Lean L2 result easy to find and remove
obsolete presentation paths without weakening the proof boundary.

## Changes

- Added `CLPoly/Refinement/README.md` as the human-facing index of the two
  end-to-end results and their component contracts.
- Ordered `CLPoly/Refinement/Generated.lean` so the public `__factor_Zp` and
  `__factor_squarefree_primitive_ZZ` theorems appear first.
- Renamed `StrictSquarefreeGenerated.lean` to `SquarefreeZpEntry.lean`: the
  module is a checked entry-point refinement, not generated proof text.
- Renamed `AxiomAudit.lean` to `FinalAxiomAudit.lean` and documented its role.
- Removed the stale `Generated/Corpus.lean.bak` snapshot.
- Removed Pass 9's unreachable legacy skeleton emitter, including every path
  that could write `sorry`, placeholder types, or incomplete theorem files.
  Pass 9 now emits only centralized wrappers around completed strict proofs.
- Simplified the public `CLPoly.lean` umbrella to expose the stable
  specification, L2 algorithms, strict refinements, and L1 pipeline layers.
  Experimental `E1`/`E2`/`E3` modules remain in the repository but are no
  longer part of the public import surface.

No production C++ source was changed.

## Verification

- `python3 passes/pass9_refine_gen.py`: regenerated the public contracts.
- all strict component generators with `--check`: reproducible and
  placeholder-free.
- `python3 tests/check_verified_refinements.py`: passed.
- `python3 tests/check_strict_refinement_boundary.py`: passed.
- `python3 tests/check_b2b_strict_runtime_isolation.py`: passed.
- `bash tests/run_all_smoke.sh`: all cpp2lean v2 unit and Pass 1--7 smoke
  tests passed (68 translated targets, zero failures in each full pass).
- `lake build`: passed, 3556 jobs.
- full C++/Lean B2B suite: 158/158 passed, including all 6 FactorZp and all 5
  FactorZZ final-entry cases.
- `lake env lean CLPoly/Refinement/FinalAxiomAudit.lean`: passed. Both public
  final theorems depend only on `propext`, `Classical.choice`, and
  `Quot.sound`.

The cleanup changes discoverability only: the public results still execute
the generated strict, well-founded L1 entries and prove `FactorZpCorrect` and
`FactorZZCorrect`; no fuel model or semantic oracle was introduced.
