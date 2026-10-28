# Prove EDF subtract-one polynomial semantics

The exact cpp2lean-generated `__upoly_subtract_one_raw_ir` loop is now related
to a structural list trace. The trace preserves the generated branch behavior:
it decrements an existing constant coefficient, removes it when the result is
zero, or appends `-1` when no constant term exists.

The trace's decoded polynomial is proved equal to the input polynomial minus
one. This extends the earlier no-constant result to all reduced term lists and
does not replace the generated execution with an L2 oracle.

Files changed:

- `proof/lean/CLPoly/Refinement/DDFSubtractX.lean`
- `proof/lean/CLPoly/Refinement/EDFSubtractOne.lean`
