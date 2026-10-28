# EDF executable safe path

Date: 2026-08-06

## Implementation

- Added `CLPoly.Generated.EDFSafe` with total counterparts for:
  - sparse polynomial remainder;
  - RNG-driven random polynomial construction;
  - the characteristic-two trace loop;
  - subtraction of the constant polynomial one;
  - one Cantor–Zassenhaus split attempt;
  - fuel-bounded EDF retry and recursive splitting.
- Retry failures advance the RNG state. This matters operationally: retaining the
  old state would retry the same candidate until fuel exhaustion.
- The RNG helper passes `p` as the exclusive bound, matching C++'s inclusive
  coefficient interval `[0, p - 1]`.
- Extended `Refinement.edfZpSafe` to run the bounded L1 path, accept a result that
  satisfies `EDFCorrect`, and otherwise use the certified total fallback.
- Pipeline now invokes this path on `toSparsePolyZp g`, seed 42, and
  `8 * (natDegree g + 1)` fuel.

## Operational check

For `(X + 1)(X + 2) = X² + 3X + 2` over `ZMod 5`, degree target 1, seed 42,
and fuel 64, `Generated.edfZpTrySafe` evaluates to the two sparse factors
`X + 1` and `X + 2` and a progressed RNG state.

## Proof interface

- `edfZpSafe_correct` proves the Pipeline result satisfies `EDFCorrect`.
- `edfZpSafe_of_try_correct` proves a successful certified bounded execution is
  adopted verbatim.
- `edfZpSafe_of_try_exhausted` proves fuel exhaustion selects the certified
  fallback.

## Final audit

- EDF-scoped scan found no `sorry`, `admit`, or `native_decide` in
  `Generated/EDFSafe.lean`, `Refinement/EDF.lean`, or the Pipeline integration.
- Axiom output for `edfCertified_correct`, `edfZpSafe_correct`, both execution
  selection theorems, `edf_l1_correct`, and `factor_Zp_l1` is exactly
  `[propext, Classical.choice, Quot.sound]`.
- `lake build` succeeded with 3465 jobs.
