# HGCD middle and second-child induction boundary

## Date

2026-08-09

## What changed

- Added `HgcdRecursiveMiddleWorkspace`, bundling the concrete allocation,
  capacity, and non-aliasing facts used by the generated middle divrem.
- Added `HgcdRecursiveSecondCallbackRefines`, the semantic proposition supplied
  by induction for the actual second recursive input.
- Changed the complete non-early body theorem to accept the physical second
  workspace and second-child induction theorem separately, assembling
  admissibility internally.

## Why

The top-level strong recursion must distinguish invariant memory-layout inputs
from semantic results derived for smaller calls.  The previous non-early API
still accepted a preassembled second-call admissibility record, obscuring that
boundary.  The revised theorem now exposes the exact induction injection point
at `middle.lenC0 < bound`.

## Key decisions

- All twelve middle divrem premises remain unchanged in meaning; only their
  packaging changed.
- The second callback proposition is indexed by the actual `middle` record and
  both source-derived ordering/decrease proofs.
- No result, branch, or raw heap operation was replaced.

## Verification

- `lake env lean CLPoly/Impl/StrictHGCDRawRefinement.lean` passed.
- Existing complete non-early proof type-checks through the new boundary.
- No proof placeholder or bounded recursion was introduced.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-middle-second-ih-boundary.md`
