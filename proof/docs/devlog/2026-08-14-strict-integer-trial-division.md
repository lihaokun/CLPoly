# Strict integer trial division for recombination

Audit found that the legacy aggregate corpus had no specialized
`HasPolyDivmod SparsePolyZZ` instance, so its generic `pair_vec_div5` dispatch
could evaluate to `(default, default)`.  The strict FactorZZ dependency closure
now bypasses that path entirely.

The recombination generator emits a concrete sparse integer long division:

- append the shifted, negatively scaled divisor terms to the current
  remainder;
- normalize the concrete term array;
- append the computed leading quotient term;
- recurse only after checking that the actual normalized remainder degree
  strictly decreased.

The recursion measure is the current remainder degree (`divisionRank`), not a
fuel parameter.  Nonexact leading-coefficient division returns the current
quotient and nonempty remainder, matching the trial-division use site;
malformed zero divisors or nondecreasing intermediates are raw faults.

The refinement proof establishes the scaled-subtraction polynomial semantics
and then the complete successful-execution identity

`toPoly dividend = toPoly divisor * toPoly quotient + toPoly remainder`.

Candidate validation now calls this generated function directly and has no
remaining divmod callback.  The strict-boundary audit continues to reject
`HasPolyDivmod.polyDivmod`, the legacy corpus, partial definitions, fuel, and
admitted proofs.

Estimated FactorZZ progress: about 94%.  Remaining work is the full successful
candidate extraction/product invariant, Zassenhaus fallback semantics,
van-Hoeij result closure, and generated FactorZZ/Pipeline final theorem and
full audit.  Estimated focused time: 1--2 workdays.
