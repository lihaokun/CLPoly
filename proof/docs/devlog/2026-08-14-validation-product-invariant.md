# Candidate-validation product invariant

The local successful trial-extraction identity is now lifted through the full
generated candidate-validation loop.

For every successful execution, including all skipped-candidate, failed-trial,
and repeated-extraction branches, the proof produces the concrete accumulated
integer scalar `s` such that

`toPoly fStar * product(existing result)`

equals

`C(s) * (toPoly fStar' * product(updated result))`.

The scalar is the product of quotient contents returned by the actual
generated primitive-normalization calls.  It is kept explicit rather than
silently assumed to be a unit.  This isolates the next mathematical obligation:
derive unit content from the primitive top-level input and the concrete
primitive outputs.

The proof is by the same well-founded candidate-index measure as the generated
execution and handles every raw callee result explicitly.  It uses no semantic
validation callback or factorization witness.

Estimated FactorZZ progress: about 94.5%.  Remaining work is the content-unit
bridge, irreducibility of emitted recombination groups, Zassenhaus fallback and
main-loop closure, then generated FactorZZ/Pipeline integration and full audit.
Estimated remaining focused time: 2--4 workdays.
