# HGCD child workspace and induction separation

## Date

2026-08-09

## What changed

- Added `HgcdRecursiveSecondCallWorkspace`, containing only the matrix-prefix,
  iterator workspace, frame, and finish workspace for the actual second child.
- Refactored `HgcdRecursiveSecondCallAdmissible` into a physical workspace plus
  the smaller-call semantic theorem.
- Added constructors that attach a well-founded child theorem to an existing
  first- or second-call physical workspace.
- Migrated the complete non-early body proof to the separated projections.

## Why

The top-level HGCD fixpoint must obtain recursive polynomial semantics from its
strictly smaller induction hypotheses.  Keeping those semantics mixed with
memory-layout premises would make it possible for a caller to prepackage the
very theorem the fixpoint must prove.  Both recursive call sites now have the
same boundary: physical facts are inputs, semantic facts are attached from the
well-founded induction result.

## Key decisions

- No execution definition or computational result changed.
- The second-child package remains tied to the actual `middle` and `second`
  records and to their concrete strict-decrease proofs.
- The constructor theorems do not synthesize results; they only combine an
  existing physical workspace with the exact callback refinement theorem.

## Verification

- `lake env lean CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `lake build CLPoly.Impl.StrictHGCDRawRefinement`: 3302 jobs succeeded.
- Relevant implementation search found no placeholders, bounded recursion,
  alternate implementation, or `partial def`.
- Axiom audit reports only `propext`, `Classical.choice`, and `Quot.sound`.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `docs/devlog/2026-08-09-hgcd-child-workspace-separation.md`
