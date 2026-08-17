# Execute the strict FactorZZ final entry back to back

## What changed

The B2B registry now exposes the original C++ entry
`__factor_squarefree_primitive_ZZ`.  Its Lean side executes the generated
strict outer controller and supplies concrete physical implementations for
every major callee:

- prime selection, including actual candidate iteration and modular checks;
- polynomial reduction, derivative, raw GCD, make-monic, DDF, and EDF;
- the exact seeded `std::mt19937`-compatible random stream and concrete EDF
  retry traces;
- multifactor Hensel lifting with concrete raw division and multiplication;
- concrete van-Hoeij and exhaustive Zassenhaus recombination.

No L2 factorization routine, factor witness, candidate-prime table, fuel
parameter, or partial definition supplies the observed result.  The runtime
is test-only: the reverse-dependency audit prevents its erased `Prop` fields
from entering `Generated`, `Refinement`, or `Pipeline`.  The formal final
entry theorem remains the safe, kernel-checked theorem in the refinement
boundary.

## Large-prime execution

Direct evaluation of Lean's mathematical `Nat.Prime` decision procedure made
the descending UInt64 prime iterator impractically slow.  The B2B runtime now
computes the same candidate search with deterministic seven-base UInt64
Miller--Rabin and logarithmic modular exponentiation.  Every candidate is
tested; there is no hard-coded selected prime or factorization result.  All
loops use their actual source variants (`exponent`, square count, or machine
candidate rank), not an artificial fuel limit.

The runtime tests cover ordinary primes, the top-of-UInt64 initial prime, and
two strong pseudoprimes.  The complete large-prime FactorZZ request also
matches the GMP-backed C++ entry.

## Observable result ordering

C++ and Lean may enumerate the same irreducible factors in different orders.
Both B2B serializers now sort the complete encoded `ArraySPZZ`, as the existing
FactorZp serializers already do.  This changes only factor ordering: it does
not remove, merge, normalize, or replace any returned factor.  Explicit
expected results independently pin every factor in every final-entry vector.

## Verification

- strict final FactorZZ C++/Lean B2B: 5/5 passed;
- covered an irreducible result, a two-factor split, three factors, coefficients
  requiring nontrivial Hensel precision, and `useLargePrime = true`;
- exact expected factor arrays are checked in addition to cross-driver equality;
- deterministic UInt64 primality/runtime tests: 8/8 passed;
- existing MT19937/uniform tests: 8/8 passed;
- C++ B2B codec round trips: all passed;
- `B2B.StrictRuntime`, `B2B.Registry`, and `B2B.Driver`: built successfully;
- B2B reverse-dependency isolation audit: passed;
- `git diff --check`: passed.

No production C++ source changed in this step, so it creates no new native
regression surface.  The B2B driver executes the current production C++
FactorZZ entry directly.

## Files

- `proof/lean/B2B/StrictRuntime.lean`
- `proof/lean/B2B/Registry.lean`
- `proof/lean/B2B/Types.lean`
- `proof/lean/B2B/TestStrictRuntime.lean`
- `proof/b2b/types/b2b_types.cc`
- `proof/b2b/vectors/__factor_squarefree_primitive_ZZ.json`
- `proof/docs/devlog/2026-08-17-factor-zz-final-entry-b2b.md`
