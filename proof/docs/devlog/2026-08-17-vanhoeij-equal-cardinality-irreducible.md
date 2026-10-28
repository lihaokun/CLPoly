# Full-cardinality van-Hoeij irreducibility

## What changed

The FactorZZ composition layer now proves that a full-cardinality result from
the literal concrete `__vanhoeij_recombine_raw_ir` execution consists entirely
of irreducible integer polynomials.

The theorem combines:

- the associated integer product returned by that execution;
- the physical nonunit/modular-nonunit/leading-coefficient facts returned by
  that execution;
- the actual `LiveHenselProduct` certificate;
- modular irreducibility of the physically lifted factors; and
- equality between output and lifted cardinalities.

## Why

`__lll_factorize` accepts a van-Hoeij result without the exhaustive safety net
only when its physical cardinality has not decreased.  This theorem closes
the mathematical obligation for precisely that controller branch.  A
reduced-cardinality result remains assigned to the literal full-precision
Zassenhaus branch.

## Key decisions

- The recombination product helper was made public within the refinement
  namespace so the cross-module FactorZZ theorem can state the exact generated
  product without duplicating it.
- Output primitivity is derived from associatedness with the primitive source,
  not postulated separately.
- No result array or irreducibility witness is supplied by an operation field.

## Verification

- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed.
- The composition theorem was added to the axiom audit.
- Production C++ changes: none; no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
- `proof/docs/devlog/2026-08-17-vanhoeij-equal-cardinality-irreducible.md`

## Progress and estimate

The accepted van-Hoeij controller branch now has its missing irreducibility
argument.  Next is the source-shaped `__lll_factorize_raw_ir` case split and
the final squarefree-primitive entry composition.  Overall goal progress is
approximately 99.77%; estimated remaining time is 2.5--4 full working days.
