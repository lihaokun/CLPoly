# Parameterize the strict EDF RNG boundary

An audit against `clpoly/polynomial_factorize_zp.hh` found that the shared
Lean `Rng` model is a `UInt64` xorshift state, whereas the C++ EDF entry takes
`std::mt19937&`. Binding the strict EDF proof to that test model would prove
the right coefficient range but the wrong concrete RNG trajectory.

The generated strict EDF state machine is now polymorphic in its RNG state.
`RandomEngine State` exposes the source draw-and-advance boundary together
with the required `[0, upper)` range proof. The generated random loop,
`EDFRawOps`, retry state, split state, and recursive EDF entry all preserve
that exact state type. Consequently a future public theorem must supply a
real C++ RNG adapter; the xorshift compatibility model can no longer be used
implicitly.

The existing random-output proof was generalized over this engine and still
proves termination, canonical output, and the strict degree bound without
depending on a particular pseudo-random algorithm.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Generated.StrictEDF CLPoly.Refinement.EDFRandom
lake build CLPoly.Refinement.EDF
```
