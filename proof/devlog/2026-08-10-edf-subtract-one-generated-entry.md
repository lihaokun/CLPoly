# Close the generated EDF subtract-one semantic entry

The theorem
`Refinement.StrictEDF.__upoly_subtract_one_raw_ir_refines` now proves that the
exact cpp2lean-generated C++ `__upoly_subtract_one` entry terminates with an
`ok` result and that decoding that result gives the input polynomial minus
one.

The proof covers both concrete epilogue branches. It uses the generated
well-founded loop and its structural execution trace; it does not use an L2
implementation as an oracle or fallback.

Verification:

```text
lake build CLPoly.Refinement.EDFSubtractOne
#print axioms Refinement.StrictEDF.__upoly_subtract_one_raw_ir_refines
```

The axiom report contains only Lean/Mathlib foundations: `propext`,
`Classical.choice`, and `Quot.sound`.

Files changed:

- `proof/lean/CLPoly/Refinement/EDFSubtractOne.lean`
