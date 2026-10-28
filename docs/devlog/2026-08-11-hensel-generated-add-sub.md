# Generate and refine Hensel sparse addition and subtraction

## Generated L1 semantics

The strict Hensel generator now emits direct lowerings of the C++
`pair_vec_add` and `pair_vec_sub` iterator merges used by
`basic_polynomial::operator+` and `operator-`.

Each loop carries both concrete source indices and the concrete output array.
Its six control-flow cases match iterator exhaustion, degree comparison,
equal-degree coefficient combination, zero elimination, and termination.
Recursion is well founded on the sum of the two remaining suffix lengths; no
fuel parameter is present.

## Refinement

`pairVecAddLoop_toPolyMod` and `pairVecSubLoop_toPolyMod` prove invariants over
the exact unconsumed array suffixes.  The public internal bridges
`__upoly_add_raw_ir_refines` and `__upoly_sub_raw_ir_refines` establish that
the successful raw outputs decode to polynomial addition and subtraction in
`ZMod m`.

The proof does not execute `SparsePolyZZ.addReal` or `subReal`, and it does not
replace a failing raw path with an L2 result.

## Verification

```text
python3 proof/cpp2lean_v2/tests/build_strict_hensel.py --check
lake build CLPoly.Refinement.Hensel
Build completed successfully (1923 jobs).
```
