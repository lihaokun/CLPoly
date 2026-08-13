# Strict `__select_prime` entry refinement

The generated well-founded outer loop is now instantiated with an actual
runtime callback.  For every machine word proved prime by the generated GMP
enumerator, the callback constructs only the corresponding dense arithmetic
state, raw workspaces, and the EDF well-founded termination certificate.  It
then executes the strict coefficient reduction, sparse polynomial reduction,
derivative, heap GCD, make-monic, DDF, and EDF body already proved in
`tryCandidateRaw_factored_refines`.

The entry theorem unfolds the generated original `__select_prime` guard,
checks both concrete initial primes (`2` and `2^64 - 59`), and invokes the
outer loop invariant with the exact C++ initial `UINT64_MAX` best count.  Its
result satisfies `SelectionCorrect`: the stored prime is good, the stored
modular factors reconstruct the reduced integer input up to association, and
every stored factor is irreducible and monic.

No factor witness or `GoodPrime` proposition is present in the runtime
provider.  No fuel, partial recursion, semantic oracle, or unproved axiom is
used.  A subsequent representation lemma will derive the currently explicit
`front!`/L2-leading-coefficient premise from sparse ZZ canonicality.

Pass 9 now emits the public theorem
`__select_prime_raw_ir_refines_SelectionCorrect` in the centralized
`Refinement/Generated.lean` module.  The checked output is derived from the
`REFINEMENT_MAP` entry and is reproducible with
`check_verified_refinements.py`; it is not a manually maintained wrapper.
The strict-boundary audit was also corrected to strip nested Lean comments
before looking for banned constructs, so documentation saying that an
implementation has no execution counter no longer creates a false failure.
