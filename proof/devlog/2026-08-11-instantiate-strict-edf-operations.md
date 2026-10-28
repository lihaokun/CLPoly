# Instantiate the strict EDF operation record

`traceLoop_terminates_canonical` proves the generated characteristic-two loop
terminates and preserves canonical representation by induction on the exact
source distance `d-i`, provided its concrete square/add/mod call has the
already-proved execution contract.

`strictEDFRawOps` is now the unique computational record supplied to the
generated EDF shell. It connects random generation, raw modular reduction,
square/add/mod, powmod, GCD, exact division, and make-monic to their strict L1
implementations. The recursive invariant is `EDFEntryInvariant`; no semantic
result or split proof is stored in this computational record.

Verification:

```text
lake build CLPoly.Refinement.EDF
```
