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
