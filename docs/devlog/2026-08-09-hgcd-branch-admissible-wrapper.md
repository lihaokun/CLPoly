# HGCD branch-admissible body wrapper

## Date

2026-08-09

## What changed

- Added `HgcdRecursiveNonBasePackage`, which packages the polynomial split and
  actual first-child admissibility for a non-base recursive invocation.
- Added `hgcdRecursiveBodyBranchAdmissible`.  It requests that package only
  after the generated half-length base guard fails.
- Proved `hgcdRecursiveBodyBranchAdmissible_eq_body`, showing that both guard
  branches erase exactly to the source-shaped `hgcdRecursiveBody`.

## Why

The older admissible wrapper requested first-child data unconditionally even
though the generated base branch performs no recursive call.  The lower-level
body had already been corrected, so this change carries the same exact control
flow boundary through the public body wrapper needed by the forthcoming
well-founded fixpoint.

## Key decisions

- The non-base package is proof/input representation data; it does not select
  or replace the `RawExec` result.
- The base path constructs only contradiction eliminators for unreachable
  proof arguments and executes the same generated base body.
- The erasure theorem compares full `RawExec` values, retaining both faults and
  every successful computational field.

## Verification

- `lake env lean CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `lake build CLPoly.Impl.StrictHGCDRawRefinement`: 3302 jobs succeeded.
- `git diff --check`
- Relevant implementation search found no proof placeholders, fuel recursion,
  fallback implementation, or `partial def`.
- The new equality theorem depends only on `propext`, `Classical.choice`, and
  `Quot.sound`.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-branch-admissible-wrapper.md`
