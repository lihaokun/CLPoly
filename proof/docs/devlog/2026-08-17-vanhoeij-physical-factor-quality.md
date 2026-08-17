# Van-Hoeij physical factor quality

## What changed

The physical factor-quality invariant proved for concrete Zassenhaus execution
is now propagated through the complete generated van-Hoeij loop and its public
`__vanhoeij_recombine_raw_ir` entry.

The proof follows all generated branches:

- candidate preparation and literal `validateCandidates` extraction;
- active-index removal and lattice reset;
- bounded precision retry with the same lexicographic well-founded measure;
- actual concrete Zassenhaus fallback and physical array append; and
- final positive-degree append followed by degree sorting.

For every returned factor it proves integer nonunit, selected-prime modular
nonunit, and nonzero leading coefficient modulo the selected prime.

## Why

The full-cardinality van-Hoeij acceptance branch in `__lll_factorize` must use
the same modular-cardinality irreducibility argument as the concrete
Zassenhaus/component proofs.  That argument was previously missing facts for
the actual van-Hoeij output members.

## Key decisions

- The theorem is specialized to `concreteVanHoeijRawOps` and
  `concreteVanHoeijTermination`.  This is intentional: the fallback proof must
  know that the generated operation record executes the concrete legal
  combination recursion, not an arbitrary abstract termination record.
- Gathered active factors inherit irreducibility by physical membership in the
  original lifted array.
- Successful validation obtains the next `fStar`, accumulated-result quality,
  leading coefficient, canonicality, and primitivity from the same raw run.

## Verification

- `lake env lean CLPoly/Refinement/Recombine.lean`: passed.
- New entry and loop theorems added to the axiom audit.
- Production C++ changes: none; no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
- `proof/docs/devlog/2026-08-17-vanhoeij-physical-factor-quality.md`

## Progress and estimate

Both concrete recombination engines now expose the physical member facts
needed by final composition.  The next proof combines van-Hoeij size/product/
quality with the selected Hensel certificate in the full-cardinality branch,
and uses the existing terminal Zassenhaus theorem in the safety branch.
Overall goal progress is approximately 99.73%; estimated remaining time is
3--4 full working days for the final FactorZZ controller theorem, concrete
operation wiring, generated public theorem, Pipeline/L1 integration, and all
full verification gates.
