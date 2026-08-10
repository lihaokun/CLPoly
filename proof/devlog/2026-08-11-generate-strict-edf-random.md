# Generate the strict EDF random-polynomial entry

The strict EDF generator now emits a total, well-founded lowering of the C++
`__upoly_random` descending-degree loop.

Each iteration calls `Rng.next_advance rng p`, giving the same `[0,p)` range
as C++ `uniform_int_distribution(0,p-1)`, advances the RNG state, omits zero
coefficients, and appends nonzero coefficients in descending degree order.

This deliberately does not reuse the old corpus-level
`UniformIntDist.mk = 0` placeholder, which would make every generated
polynomial empty and make retry progress impossible. A direct theorem proves
the new raw entry always returns its computed polynomial and advanced RNG
state.

The follow-up invariant proof inspects those concrete draws. It proves that
every retained coefficient is reduced modulo `p`, every retained coefficient
is nonzero, and the descending `Array.push` order preserves
`SparsePolyZp.Canonical`. Consequently
`__upoly_random_raw_ir_canonical` certifies the exact generated output without
post-normalization or an L2 replacement.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Generated.StrictEDF CLPoly.Refinement.EDFRandom
lake build CLPoly.Refinement.EDF
```

Files changed:

- `proof/cpp2lean_v2/tests/build_strict_edf.py`
- `proof/lean/CLPoly/Generated/StrictEDF.lean`
- `proof/lean/CLPoly/Refinement/EDFRandom.lean`
- `proof/lean/CLPoly/Refinement/EDF.lean`
