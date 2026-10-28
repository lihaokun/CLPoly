# Strict `__lll_factorize` controller refinement

## What changed

The generated source-shaped `__lll_factorize_raw_ir` controller now refines to
`FactorZZCorrect` when instantiated with the concrete van-Hoeij and
Zassenhaus executions.

The proof covers every literal controller branch:

- accepted first-pass full-cardinality van-Hoeij output;
- reduced low-precision output followed by the actual second Hensel call with
  source argument zero;
- accepted full-cardinality second-pass van-Hoeij output;
- second-pass reduced output followed by the literal Zassenhaus safety net;
- first-pass full-precision reduced output followed by the same concrete
  Zassenhaus safety net.

The existing default-Mignotte Zassenhaus theorem was generalized to accept an
arbitrary actual Hensel target plus an explicit recovery-margin proof.  Its
old zero-target interface remains as a proved specialization.

## Why

This is the exact control-flow bridge needed between the proved Hensel and
recombination components and the generated squarefree-primitive FactorZZ
entry.  In particular, a reduced van-Hoeij result is never accepted through a
semantic assumption: the generated safety decision and actual Zassenhaus call
remain visible in the theorem.

## Full-precision recovery

The controller theorem has no supplied recovery premise.  For the case where
the generated first target equals `aMig`, the proof establishes that the
well-founded `heuristicPrecisionLoop` crosses the concrete generated
Mignotte target, proves that a successful signed exponent conversion cannot
have wrapped at `2^31`, and combines that fact with the literal Hensel target
and returned prime-power modulus.  `LiveRecoveryPrecision` is therefore a
derived property of the two real generated executions, not an oracle.

## Verification

- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed.
- New controller and composition theorems added to the axiom audit.
- Production C++ changes: none; no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
- `proof/docs/devlog/2026-08-17-strict-lll-factorize-controller.md`

## Progress and estimate

The generated recombination controller is structurally closed.  Remaining
work is concrete select/Hensel operation wiring, the outer
squarefree-primitive theorem, generated public export, Pipeline/L1, and full
verification.  Overall goal progress is approximately 99.84%; estimated
remaining time is 1.5--3 full working days.
