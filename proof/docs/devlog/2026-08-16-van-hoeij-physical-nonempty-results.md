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
