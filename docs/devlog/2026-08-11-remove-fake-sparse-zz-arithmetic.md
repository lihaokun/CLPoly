# Remove fake SparsePolyZZ arithmetic

## Problem

The global `SparsePolyZZ` `HMul`, `HAdd`, and `HSub` instances all returned
array concatenation.  In particular, multiplication did not multiply
monomials or coefficients.  Although strict Hensel currently takes its raw
operations through an explicit `HenselStepRawOps` interface, these instances
were a dangerous residual path back to non-semantic arithmetic.

## Change

- Removed the three concatenation instances.
- Reused the existing normalized sparse addition and subtraction functions.
- Added convolution multiplication followed by sparse normalization.
- Installed the three global instances only after those real operations are
  defined.

This is a cleanup prerequisite, not the final C++ refinement bridge.  The
strict Hensel raw operations still need to be instantiated from the generated
`pair_vec_add`, `pair_vec_sub`, and `pair_vec_multiplies` semantics and related
to these mathematical results before the public Hensel contract is complete.

## Verification

```text
lake build CLPoly.Model CLPoly.Refinement.Hensel
Build completed successfully (1923 jobs).
```
