# Strict Hensel tree-builder execution bridge

## What changed

- Connected the generated tree builder's multiplication callback to the
  verified strict sparse-to-dense generated multiplication path.
- Connected its EEA callback to the generated, degree-measured polynomial EEA.
- Proved the two concrete half-open product loops terminate and denote the
  corresponding products of input factors.
- Proved every checked node write succeeds without changing the array size.
- Proved the generated recursive builder succeeds by well-founded induction on
  `stop - start`; recursive calls use the exact source midpoint and append
  children in left-before-right order.
- Proved the public raw entry allocates its real root node and successfully
  executes the recursive builder for every canonical, nonempty factor array of
  size at least two.

## Refinement boundary

This milestone is an execution and termination bridge, not yet the final tree
semantic theorem.  Every polynomial result is obtained from the generated L1
operations.  There is no expected tree argument, specification oracle, L2
fallback, or fuel counter.  The next theorem will strengthen the recursive
postcondition with the stored interval-product and Bezout invariants needed by
the full Hensel pipeline.

The follow-up node-local bridge now also proves that `poly_convert` preserves
the represented polynomial, that the six concrete stores record the exact two
interval products and EEA certificate, and that coprime interval products turn
the monic EEA gcd into the unit Bezout certificate consumed by lifting.  The
remaining strengthening is recursive: preserve these algebraic fields while
child links are installed and establish the invariant at every appended node.

The recursive frame strengthening is now proved as well.  A child build may
change its own root and append descendants, but it leaves every other
pre-existing array entry identical.  Consequently the root's concrete
interval-product and monic-gcd certificate survives both recursive calls.  A
semantic entry theorem exposes that certificate, and a coprime-root corollary
turns it into the unit Bezout invariant required by the first lift.  The next
remaining step is to collect the recursively established child certificates
into the complete generated tree topology invariant, rather than exposing only
the root certificate.

The coprimality premise is also no longer required separately for each split.
Two additional `stop - index` well-founded proofs show first that one factor
is coprime to a later interval product, and then that products over any two
adjacent intervals are coprime.  Therefore an index-level pairwise-coprime
property of the concrete factor array supplies the unit-gcd fact at the root
and, in the forthcoming recursive topology theorem, at every child split.

## Verification

```text
lake env lean CLPoly/Refinement/Hensel.lean
git diff --check
placeholder scan: no sorry, admit, axiom, implemented_by, or opaque declaration
```

## Files

- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `docs/devlog/2026-08-11-strict-hensel-tree-builder-execution.md`
