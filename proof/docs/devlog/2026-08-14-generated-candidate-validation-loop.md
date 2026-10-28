# Generated van-Hoeij candidate validation loop

This checkpoint removes the former whole-loop `validateCandidates` callback
from the strict recombination boundary.

The generator now emits the actual C++ candidate-validation control flow:

- skip empty and trivial candidates;
- scan candidate indices for already-consumed factors;
- build the selected lifted-factor product;
- apply symmetric modular recovery and primitive-part computation;
- execute trial division and inspect the concrete remainder;
- append a successful factor, primitive-normalize the quotient, mark every
  consumed index, and continue with the updated state.

The loop is total and well-founded on `candidates.size - candidateIndex`.  Its
remaining raw-operation fields correspond only to concrete C++ callees and
return computed data, never propositions or factorization witnesses.

The refinement layer proves directly from successful executions that:

- `markConsumedLoop` preserves the consumed-array size;
- the complete validation loop preserves that size;
- successful `gatherActive` produces exactly one lifted entry per active index;
- therefore a validation result with any marked bit makes reverse deletion
  strictly shrink the active array.

Consequently the van-Hoeij main-loop termination certificate no longer needs
an external consumed-size contract.

Estimated FactorZZ progress after this checkpoint: about 92%.  Remaining work
is to instantiate and prove the concrete symmetric-mod, primitive-part, and
integer trial-division callees, lift their product invariant through the full
validation and van-Hoeij/Zassenhaus loops, then generate the FactorZZ entry,
wire Pipeline/L1, and run the complete build/test/axiom audit.  Estimated
remaining focused time: 1.5--3 workdays.
