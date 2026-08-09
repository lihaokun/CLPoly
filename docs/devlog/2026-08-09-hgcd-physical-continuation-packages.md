# HGCD physical continuation packages

## Date

2026-08-09

## What changed

- Added separate physical packages for the successful early and non-early
  continuations of the generated HGCD recursive body.
- Recorded the actual first dispatch, reconstruction, middle division,
  second dispatch, and finish executions in the appropriate package.
- Added wrapper theorems that attach the first and second strictly-smaller
  recursive semantic hypotheses to those physical packages and derive the
  common raw HGCD invariant.

## Why

The top-level well-founded proof must case-split on the generated execution
and recover branch-local memory evidence without accepting the desired
polynomial semantics as an oracle.  These packages provide that physical
boundary while keeping recursive semantics exclusively in the induction
hypotheses.

## Key decisions

- The packages contain no `HgcdRecursiveCallbackRefinesAt` field.
- The early package has only the first-child semantic hypothesis at its
  wrapper theorem.
- The non-early package requires exactly two semantic hypotheses, each
  indexed by the actual smaller input length used by its generated call.
- No fuel parameter or L2 fallback was introduced.

## Problems and resolution

The first-child decrease proof also needs the non-base guard.  The guard is
therefore an index of both continuation packages rather than an untracked
assumption reconstructed later.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-physical-continuation-packages.md`
