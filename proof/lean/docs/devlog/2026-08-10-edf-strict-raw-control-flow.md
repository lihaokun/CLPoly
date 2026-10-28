# EDF strict raw control flow

The generated EDF model now replaces cpp2lean's `partial def` loops with total,
well-founded recursion while preserving the C++ `__edf_Zp` control flow.

- The characteristic-two trace loop decreases on `d - i`.
- The randomized retry loop decreases on an explicit rank of the real RNG state;
  it is not a fuel parameter.
- The recursive factor split decreases on polynomial degree.
- Successful splits retain equations for both the concrete random call and the
  complete candidate-producing `mod`/trace-or-`pow`/subtract/`gcd` execution.
- Failed retries may decrease only after proving that same concrete candidate
  execution failed the proper-factor guard.
- The generator rejects placeholders, partial definitions, fuel, and an L2 EDF
  execution fallback.

This commit establishes the strict executable L1 side. It does not yet claim the
final L1-to-L2 EDF refinement theorem; the next step is to instantiate the raw
operations and discharge their semantic invariants in `CLPoly/Refinement/EDF.lean`,
then publish the generated contract under `CLPoly/Refinement/Generated/EDF.lean`.

## Semantic base branch

The first direct generated-execution bridge is now in
`CLPoly/Refinement/EDF.lean`.  It proves that the exact generated degree-equals-
`d` branch executes its concrete `makeMonic` call and pushes that exact sparse
factor.  The accompanying L2 lemma proves the pushed monic equal-degree factor
is irreducible and satisfies `EDFCorrect`.  The public theorem remains absent:
retry and recursive split simulation still have to be proved before exporting
the generated contract.

## Concrete subtract-one primitive

The odd-characteristic candidate path no longer obtains `subtractOne` from an
abstract operation record.  The generator now emits the exact
`__upoly_subtract_one` range-for control flow as a total function, decreasing on
the unprocessed input-array suffix.  `candidateRun` calls this generated entry
directly, so a future operation-record instance cannot substitute an oracle or
an L2 subtraction for this C++ step.

The first structural loop invariant is proved in
`CLPoly/Refinement/EDFSubtractOne.lean`: after the unique constant term has
been processed, the exact generated loop preserves `found = true` and copies
the remaining concrete array suffix verbatim.  The proof unfolds the generated
well-founded recursion and tracks `Array.toList`; it assumes no polynomial
semantic result.

The complementary `found = false` invariant is also closed for inputs whose
suffix contains no constant term.  From the exact loop trace and exact
normalization epilogue, the proof now derives an actual successful raw-entry
equation and the L2 identity `toPoly output = toPoly input - 1`.  This covers
the source branch that appends the concrete `p - 1` coefficient; the existing
finite-field coefficient lemma was made reusable rather than reproved or
replaced by an oracle.
