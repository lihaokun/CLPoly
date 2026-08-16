# Preserve physical nonempty van-Hoeij validation results

## What changed

`exactDivmodRaw_divisor_nonempty_of_success` unfolds the generated sparse
long-division entry and proves that any successful execution on a nonempty
dividend must have taken its positive-size divisor branch.  Thus a candidate
factor accepted by literal exact division cannot be an empty sparse array.

`validateCandidatesLoop_result_nonempty` follows all branches of the generated
candidate-validation loop.  Rejected, unavailable, empty, and trivial
candidates preserve the existing result invariant.  The successful branch
uses the actual trial product, symmetric recovery, primitive call, exact
division, quotient primitive call, consumed marking, and recursive push; the
new exact-division lemma proves that the physically pushed factor is nonempty.
`validateCandidates_result_nonempty` packages the invariant for the complete
generated validation entry.

## Why

The van-Hoeij cardinality closure requires every returned polynomial to be
nonzero and nonunit before equality with the modular irreducible-atom count can
force one normalized factor per output.  This step establishes the first of
those facts at the physical sparse-array level without replacing execution by
an abstract factorization witness.

## Verification and remaining work

- `CLPoly/Refinement/Recombine.lean` compiles directly.
- The theorem is included in `CLPoly.Refinement.AxiomAudit`.
- C++ changes: none, so there is no new C++ regression or B2B change surface.
- Next: preserve canonicality for every pushed factor, derive polynomial
  nonzeroness, and use the nonempty candidate's positive modular degree to
  rule out units.

## Recover canonical and nonzero L2 factors

`validationRecoveredFactor_canonical` now packages the exact successful
validation prefix.  Starting from the current canonical nonempty `fStar`, it
proves the physical leading-coefficient singleton canonical, transports this
through the actual trial-product loop, obtains positivity from the successful
`symmetricModRaw` branch, and then follows symmetric recovery and the literal
primitive call to the returned factor.

`sparsePolyZZ_toPoly_ne_zero_of_canonical_nonempty` proves the general sparse
representation bridge by comparing the physical head with the L2 leading
coefficient.  `validationRecoveredFactor_ne_zero` combines it with the actual
successful exact-division trace, whose divisor-size check supplies physical
nonemptiness.  Thus the precise candidate factor that the generated loop can
push is nonzero as an L2 polynomial.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: thread this certificate through the complete result array and prove
  positive modular degree/nonunit status.

## Thread L2 nonzeroness through complete validation

`validateCandidatesLoop_result_ne_zero` strengthens the live validation
invariant to carry both the canonical current quotient and pointwise L2
nonzeroness of the physical result array.  On the successful branch it applies
`validationRecoveredFactor_ne_zero` to exactly the factor being pushed,
constructs the pushed-array invariant occurrence-sensitively, and derives the
next quotient's canonical representation from the same exact-division and
primitive-normalization traces.  Every non-success branch reuses the unchanged
physical state.

`validateCandidates_result_ne_zero` specializes the invariant to the complete
generated entry.  The final result array needed by the normalized-factor count
therefore has a direct execution-derived nonzero certificate for every member.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: establish nonunit status using the positive degree of nonempty products
  of the selected modular irreducible factors.

## Prove an accepted physical validation factor is nonunit

`candidateAvailableLoop_true_valid_suffix` reverses the literal checked-index
scan: a `true` result proves that every remaining `Int32` index passed both
nonnegativity and physical-array bounds.  Its entry specialization produces
the exact `CandidateIndicesValid` certificate consumed by the trial-product
refinement.

For a nonempty candidate, `selectedProductMod_not_isUnit_of_nonempty_irreducible`
uses the actual first checked array occurrence.  That occurrence denotes an
irreducible modular factor and divides the exact selected list product, so the
product cannot be a unit.  `validationRecoveredFactor_mod_eq_selected` then
follows the same trial-product, symmetric-recovery, and primitive-content
trace at the prime divisor of the Hensel modulus.

Finally, `validationRecoveredFactor_not_isUnit` proves the physical recovered
integer factor is nonunit.  Nonzero source leading coefficient and the
nonzero selected irreducible product force the actual primitive content to be
nonzero modulo the prime.  Removing the two resulting constant units gives an
association between the mapped recovered factor and exact selected product;
an integer unit would map to a modular unit, contradicting that association.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: thread this nonunit invariant through all successful validation pushes
  and combine it with the result-count equality theorem.

## Preserve the selected-prime leading coefficient after validation

`validationQuotient_leading_mod_ne_zero` proves that the exact primitive
quotient installed after a successful validation extraction retains a nonzero
leading coefficient modulo the selected prime.  It uses
`successfulTrialExtraction_toPoly`, hence the same returned exact quotient and
the same literal quotient-primitive call.  Taking leading coefficients and
mapping the resulting integer equality modulo the prime shows that a zero next
leading coefficient would force the prior live leading coefficient to zero.

This closes the recursive state premise needed to reapply
`validationRecoveredFactor_not_isUnit` at every later successful candidate,
rather than only at the first validation push.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: package canonicality, nonempty live quotient, nonzero modular leading
  coefficient, consumed-array size, and pointwise nonunit results into the
  complete generated validation-loop invariant.

## Preserve nonunit results through the generated validation loop

`validateCandidatesLoop_result_nonunit` follows the literal generated
candidate-index recursion and carries the complete physical state needed by
the final van-Hoeij cardinality argument.  At each accepted exact division it
uses the actual trial product, symmetric recovery, primitive normalization,
exact quotient, quotient normalization, and consumed-marker update.  It proves
simultaneously that the live quotient remains canonical and physically
nonempty, its leading coefficient remains nonzero modulo the selected prime,
the consumed bitmap retains the active lifted-factor length, and every factor
actually pushed into the result array is nonunit over `ZZ`.

`validateCandidates_result_nonunit` instantiates that invariant at the public
generated `validateCandidates` entry, including its literal false-filled
consumed bitmap and active-factor count.  Neither theorem replaces execution
with an abstract factorization or existence oracle; both invert and recurse
over the generated `RawExec` control flow.

- Termination measure: `candidates.size - candidateIndex`, inherited directly
  from the generated well-founded validation loop.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: apply the invariant to the validation call inside the generated
  van-Hoeij iteration, combine nonzero/nonunit output membership with exact
  output cardinality and associated products, and close irreducibility through
  the literal retry loop.
