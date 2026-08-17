# Derive the Hensel entry invariant from concrete FactorZZ execution

## What changed

The generated public refinement theorem for
`__factor_squarefree_primitive_ZZ` now asks only for
`HenselLiftRuntimeReadiness`.  This residual contract contains the two
low-level execution obligations for the recursive lift and normalization
stages.  It does not contain algebraic facts about the input, selected
factors, target exponent, or eventual result.

`concreteHenselLift_success_of_selection` reconstructs the complete
`HenselLiftEntryInvariant` before applying the existing concrete Hensel
refinement.  In particular, it derives:

- the source leading term from the nonempty physical input array;
- exponent admissibility from successful execution of the generated raw
  Hensel controller (the negative-exponent branch faults);
- nonnegativity of the generated Mignotte or explicit target;
- canonical and nonempty adjusted factors from the literal first-factor
  adjustment and the physical prime-selection result;
- monicity from `SelectionCorrect.quality`;
- pairwise coprimality from the selected squarefree modular product.

The generic composition theorem also passes the actual
`SelectionCorrect` and `SelectionPhysical` evidence to its Hensel callback,
so these properties cannot be supplied for an unrelated or chosen factor
array.  Pass 9, `Refinement/Generated.lean`, and `Pipeline/L1.lean` expose the
same reduced public boundary.

## Why

The old `HenselLiftEntryReadiness` public premise bundled facts that are
already consequences of the real C++ selection, adjustment, target
calculation, and successful Hensel execution.  Leaving them at the final
boundary weakened the end-to-end refinement statement by asking callers to
certify source-derived algebra instead of proving it along the executed path.

The remaining runtime readiness is universally quantified over the actual
intermediate results of the generated low-level calls.  It cannot choose the
Hensel output, modulus, factorization, or recombination result and therefore
is not an existence or semantic oracle.

## Verification

- Pass 9 centralized refinement output check: passed.
- Pass 9 Python syntax check: passed.
- `lake build CLPoly.Refinement.FactorZZ CLPoly.Refinement.Generated
  CLPoly.Pipeline.L1`: passed (3551 jobs).
- strict refinement-boundary audit: passed.
- permanent axiom audit extended with the three new derivation theorems;
  only `propext`, `Classical.choice`, and `Quot.sound` occur.
- `git diff --check`: passed.

No production C++ file changed in this proof-only step, so it creates no new
native or C++/Lean B2B change surface.  Full entry-level FactorZp and FactorZZ
B2B remains the final acceptance gate.
