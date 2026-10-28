# HGCD recursive providers made branch-local

## Date

2026-08-09

## What changed

- Changed `hgcdRecursiveBodyBelow` so the first-recursion length provider and
  first-reconstruction bound provider are demanded only after the generated
  non-base guard has been taken.
- Migrated the strict HGCD raw-refinement proofs and execution-erasure theorem
  to the branch-local provider interface.
- Removed both recursive-provider parameters from the base-case raw-invariant
  theorem.  Its execution now supplies the impossible non-base providers by
  eliminating the contradiction with the actual base guard.

## Why

The generated C++ base path does not execute a recursive call or reconstruct a
recursive result.  Requiring those workspaces at the top-level interface made
the proof boundary stronger than the real execution and risked forcing callers
to manufacture evidence for code that is not run.  Making the providers
conditional mirrors the source control flow and is necessary for a genuine
well-founded raw-to-safe wrapper.

## Key decisions

- The computational `RawExec` result was not changed.  Only proof-valued
  arguments were moved under the actual non-base guard.
- The early and non-early proofs still derive their providers jointly from the
  same admissible first dispatch; no caller-selected semantic result was added.
- The base theorem has no fallback provider or specification oracle.

## Verification

- `lake env lean CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `lake build`: 3470 jobs completed successfully.
- `git diff --check`
- Relevant-file search found no `sorry`, `admit`, `sorryAx`, oracle, fallback,
  fuel, or `partial def` token.
- Axiom audit for the base, early, non-early, body-erasure, and admissible-body
  theorems reports only `propext`, `Classical.choice`, and `Quot.sound`.

## Files

- `proof/lean/CLPoly/Generated/StrictHGCD.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-branch-local-providers.md`
