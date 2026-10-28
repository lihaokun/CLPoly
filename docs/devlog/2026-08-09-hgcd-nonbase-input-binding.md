# HGCD non-base input binding

## Date

2026-08-09

## What changed

- Extended `HgcdRecursiveNonBasePackage` with the complete input polynomials.
- Bound each complete input to the same low/high polynomial witnesses used by
  the physical first-child workspace through the exact source split at
  `lenA / 2`.

## Why

The upcoming unified body invariant must not combine a child theorem and raw
workspace from one input with reconstruction equations for another input.
Keeping both split equalities in the same dependent package makes that
cross-call substitution impossible and lets the early/non-early branch proofs
consume one source-bound object.

## Key decisions

- The fields are proof/mathematical representation data and do not affect
  `RawExec`.
- The split uses the exact half-length and `X` shift already used by generated
  reconstruction.
- No recursive semantic theorem was reintroduced into the package.

## Verification

- `lake env lean CLPoly/Impl/StrictHGCDRawRefinement.lean` passed.
- No generated computation or admitted axiom changed.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-nonbase-input-binding.md`
