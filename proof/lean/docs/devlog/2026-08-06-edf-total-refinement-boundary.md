# EDF total refinement boundary

Date: 2026-08-06

## Outcome

- Rejected the old `__edf_Zp_ir_refines ... (edf ... [])` skeleton because the
  proposition is false for successful random splits.
- Documented the concrete termination obstruction: the translated xorshift RNG
  has `0` as a fixed point, so the generated `while (true)` EDF can keep producing
  the zero polynomial forever.
- Added `Refinement.edfZpSafe`, a total specification boundary for EDF.
- Proved `Refinement.edfZpSafe_correct` under the monic, squarefree, equal-degree
  hypotheses supplied by DDF.
- Switched `Pipeline.edf_l1` to the verified refinement-layer entry point and
  reduced `edf_l1_correct` to the new theorem.

## Audit

- EDF scoped placeholder scan: no `sorry`, `admit`, or `native_decide` in
  `CLPoly/Refinement/EDF.lean` or its Pipeline implementation.
- Axiom audit for `edfZpSafe_correct`, `edf_l1_correct`, and `factor_Zp_l1`:
  `[propext, Classical.choice, Quot.sound]`.
- Targeted build: `lake build CLPoly.Refinement.EDF CLPoly.Pipeline.L1`
  succeeded (3306 jobs).
- Full build: `lake build` succeeded (3464 jobs).

## Proof boundary

This closes the total, correctness-preserving EDF interface used by Pipeline. It
does not claim that the generated partial RNG loop terminates for every seed.
Such a claim is refuted by seed zero in the current RNG model. A future executable
refinement can add bounded retries and an explicit certified fallback, or state a
conditional trace theorem for terminating generated executions.
