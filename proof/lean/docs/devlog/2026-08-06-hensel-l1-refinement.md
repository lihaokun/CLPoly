# Hensel L1 refinement boundary

Date: 2026-08-06

## Corrections to the old skeleton

- Removed the placeholder `hensel_l2 := []` and the false equality against a
  list of integer polynomials.
- The generated routine returns integer representatives and its actual lifting
  modulus. The L2 observable candidate now maps those representatives into
  `ZMod (p^k)`, the codomain required by `HenselCorrect`.
- The explicit exponent is guarded by the `Int32` narrowing bound.
- The generated tree entry is guarded by its documented requirement that at
  least two sparse factors are present; singleton inputs use the certified L2
  branch.

## Implementation

- Added total helpers for the explicit target `p^k - 1` and both coefficient
  scaling loops in `Generated.HenselSafe`.
- Proved `henselTargetPowSafe_eq` and `henselTargetSafe_eq`.
- Added `henselGeneratedCandidate`, `henselZpSafe`, and
  `henselZpSafe_correct`.
- Added `henselZpSafe_of_generated_correct`: a generated candidate that passes
  the full contract is adopted verbatim.
- Added `toSparsePolyZZ` and switched `Pipeline.hensel_l1` to the refinement
  entry point.

## Proof boundary

The generated implementation lifts by repeated modulus squaring, so its final
integer modulus need not equal `p^k`. Correctness is observed after mapping the
integer representatives to `ZMod (p^k)`. A generated candidate is used only
when the complete `HenselCorrect` proposition holds; all unsupported or rejected
cases use the already-proved multifactor lift.

## Audit

- Hensel-scoped scan found no `sorry`, `admit`, or `native_decide` in
  `Generated/HenselSafe.lean`, `Refinement/Hensel.lean`, or the Pipeline
  integration.
- `henselTargetSafe_eq` depends only on `propext`.
- `henselZpSafe_correct`, generated-candidate adoption, `hensel_l1_correct`, and
  `factor_ZZ_cpp_correct` depend only on
  `[propext, Classical.choice, Quot.sound]`.
- Targeted Pipeline build succeeded with 3308 jobs.
- Full `lake build` succeeded with 3466 jobs.
