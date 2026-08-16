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
