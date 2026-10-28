# Derive Hensel adjusted coprimality and nonemptiness

## What changed

The FactorZZ refinement now derives two Hensel tree-entry properties from the
actual prime-selection and first-factor-adjustment certificates:

- every adjusted modular factor is physically nonempty;
- every two distinct adjusted modular factors are `IsCoprime`.

The pairwise proof first establishes a generic PID lemma: a squarefree list
product has pairwise coprime entries.  The concrete adjusted-factor product is
then rewritten, using the already proved literal adjustment theorem, to the
source polynomial mapped into the selected finite field.  Its squarefreeness
is exactly `SelectionCorrect.goodPrime.sqfree` from the actual select-prime
execution.

Nonemptiness follows from the existing proof that adjustment preserves the
irreducibility of every selected factor; an irreducible decoded polynomial is
nonzero, and the sparse representation size lemma turns that into a positive
physical array size.

## Why

These facts were fields of the public Hensel readiness argument even though
they are consequences of the concrete upstream execution.  Proving them here
removes any need for a caller to assert pairwise coprimality or nonemptiness
when the generated FactorZZ wrapper is narrowed.

No abstract factor result, tree, modulus, or output witness is supplied by
these theorems.

## Verification

- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed.
- Production C++ changes: none; there is no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/docs/devlog/2026-08-17-derived-hensel-adjusted-pairwise.md`

## Next step

Derive physical canonicality of the adjusted array, assemble the complete
algebraic invariant from the actual selection/input/successful control flow,
and expose only runtime readiness in the generated final FactorZZ theorem.
