# HGCD non-base induction boundary

## Date

2026-08-09

## What changed

- Removed the recursive semantic theorem from
  `HgcdRecursiveNonBasePackage`; it now contains only the mathematical
  low/high split and physical first-call workspace.
- Added `HgcdRecursiveNonBaseCallbackRefines` as the separate proposition that
  the top-level strong induction must establish for the actual smaller call.
- Updated the branch-admissible body wrapper and its full-execution erasure
  theorem to accept that induction fact separately.
- Added a bridge that combines the physical package with the supplied child
  theorem into first-call admissibility.

## Why

The well-founded recursive theorem must derive child semantics from strict
length decrease.  Storing that semantic result inside caller-provided physical
data would leave a circular assumption at the public boundary.  The new API
makes the induction injection point explicit and prevents the non-base package
from carrying the theorem being proved.

## Key decisions

- The body still executes the exact generated `RawExec`; proof arguments do
  not select a result or add an executable branch.
- Polynomial split witnesses remain in the package because they describe raw
  input representation, while callback refinement is kept outside.
- The base branch continues to request neither package nor callback theorem.

## Verification

- Direct Lean type-check passed.
- `lake build CLPoly.Impl.StrictHGCDRawRefinement`: 3302 jobs succeeded.
- `git diff --check` passed.
- Axiom audit reports only `propext`, `Classical.choice`, and `Quot.sound`.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-nonbase-ih-boundary.md`
