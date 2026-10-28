# Remove artificial word bounds from strict Hensel refinement

## What changed

The strict Hensel and FactorZZ refinement chain no longer assumes either
`2 * p <= UInt64.size` or `p * p <= UInt64.size` for the selected prime.
Those assumptions came from older polynomial bridge lemmas, even though the
current `Zp` implementation forms additions and products in `Nat`, reduces
them modulo `p`, and only then converts the reduced value back to `UInt64`.

`CLPoly.Math.Univariate` now proves the `Zp -> ZMod` addition and
multiplication homomorphisms directly from that actual implementation.  The
corresponding sparse-polynomial addition, subtraction, scaling,
multiplication, well-formedness, canonicality, and ring bridge lemmas require
only a positive modulus.  The strict Hensel EEA, recursive tree builder, and
full generated Hensel entry use primality to obtain that fact.

The obsolete bounds were then removed from:

- `StrictHenselEEAAlgebraicInvariant.step` and the generated EEA refinement;
- the recursive and root-level generated Hensel-tree refinements;
- `__hensel_lift_upoly_raw_ir_refines`;
- `concreteHenselLift_success`;
- the FactorZZ Hensel-product, van-Hoeij, Zassenhaus, LLL-controller, and
  outer-controller theorem chain;
- `CandidatePhysical.twicePrimeFits` and
  `SelectedPrimePhysical.primeSquareFits`.

## Why

The large-prime C++ path can select primes for which `2 * p` or `p * p` does
not fit in one machine word.  Retaining either premise at the public Hensel or
FactorZZ boundary would therefore turn the theorem into a small-prime-only
statement and would not refine the real C++ large-prime execution.  Proving
the modular semantics after reduction removes the false restriction at its
source rather than hiding it in a higher-level certificate.

## Verification

- `lake env lean CLPoly/Math/Univariate.lean`: passed;
- `lake build CLPoly.Math.Univariate`: passed;
- `lake env lean CLPoly/Refinement/Hensel.lean`: passed;
- `lake build CLPoly.Refinement.Hensel`: passed;
- `lake build CLPoly.Refinement.SelectPrime`: passed;
- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed;
- residual search over the Hensel/FactorZZ/select-prime chain found no
  `hp2`, `h2p`, `primeSquareFits`, or `twicePrimeFits` occurrence;
- `git diff --check`: passed.

No production C++ source changed in this proof step, so it introduces no new
native or C++/Lean B2B change surface.  Full FactorZZ B2B remains a required
final pipeline gate.

## Files

- `proof/lean/CLPoly/Math/Univariate.lean`
- `proof/lean/CLPoly/Refinement/EDF.lean`
- `proof/lean/CLPoly/Refinement/SelectPrime.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
