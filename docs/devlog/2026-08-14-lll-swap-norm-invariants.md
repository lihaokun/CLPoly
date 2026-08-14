# LLL swap norm invariants

## Date

2026-08-14

## What changed

Proved two exact properties of the generated `normsAfterLovaszSwap` operation:

- every updated Gram–Schmidt norm remains strictly positive when the incoming
  norms are positive;
- the product of the two norms affected by an adjacent Lovász swap is exactly
  preserved.

## Why

The strict well-founded implementation of the generated C++ LLL loop needs a
discrete determinant rank.  Establishing validity of the state after the real
swap requires the generated rational norm update to stay positive and to agree
with the unchanged higher Gram determinants.  These lemmas discharge the
norm-only part of that argument without assuming the recursive result.

## Key decisions

The proofs unfold the actual generated array update and use its exact division
formula.  They do not introduce fuel, a semantic oracle, or an abstract result
provider.  The remaining boundary-prefix argument will be connected to the
matrix through an explicit Gram–Schmidt execution invariant.

## Problems and resolution

Lean distinguishes bounded array reads from `getElem!` reads.  The positivity
proof therefore first converts the generated total-read theorem to the bounded
read at the proven index.  For the product identity, cancellation is justified
by a separately proved strictly-positive denominator.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-14-lll-swap-norm-invariants.md`
