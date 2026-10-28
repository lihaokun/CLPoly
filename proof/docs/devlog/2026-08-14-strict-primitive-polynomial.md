# Strict integer primitive-part refinement

The strict van-Hoeij candidate validator no longer accepts primitive-part
normalization as a raw callback.

The generator now emits C++ `__upoly_primitive` as three concrete pieces:

- a well-founded fold accumulating the nonnegative gcd of coefficient
  absolute values;
- the source leading-coefficient sign adjustment;
- a well-founded coefficient loop which checks nonzero/exact divisibility and
  performs integer division before appending each output term.

From a successful execution, the refinement proof establishes the exact
integer-polynomial identity

`toPoly input = C(content) * toPoly primitive`.

The empty-input branch and negative-leading-coefficient branch are included.
No primitive-factor proposition, existence witness, or semantic operation is
returned by a callback; exactness follows from the raw loop's divisibility
checks and `Int` division laws.

Estimated FactorZZ progress: about 94%.  The main remaining recombination
callee is concrete sparse integer trial division, after which the extraction
product invariant can be carried through van Hoeij/Zassenhaus and into the
generated FactorZZ/Pipeline theorem.  Estimated remaining focused time:
1--2 workdays.
