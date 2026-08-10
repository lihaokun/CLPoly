# Prove canonical output for generated EDF subtract-one

The exact cpp2lean-generated `__upoly_subtract_one_raw_ir` now has a combined
execution contract proving:

- successful termination;
- decoded polynomial semantics equal to input minus one;
- canonical sparse output representation.

Canonical preservation is proved structurally from the generated execution
trace: strict degree ordering, reduced coefficients, and nonzero stored
coefficients are each established. It is not inferred from polynomial
extensional equality.

The public `zp_sub_reduced` representation lemma is shared from the strict DDF
subtract-X proof rather than duplicating the same word-level arithmetic proof.

Files changed:

- `proof/lean/CLPoly/Refinement/DDFSubtractX.lean`
- `proof/lean/CLPoly/Refinement/EDFSubtractOne.lean`
