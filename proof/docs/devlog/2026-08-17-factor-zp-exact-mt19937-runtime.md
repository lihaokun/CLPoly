# Match the concrete C++ FactorZp random engine

## What changed

The B2B runtime now implements the source `std::mt19937(42)` state: the
624-word seed recurrence, twist, tempering and state index.  It also mirrors
the libc++ `uniform_int_distribution<uint64_t>(0, p - 1)` behavior used by
the current C++ driver, including low-bit rejection for small ranges and the
two-word independent-bits construction for ranges wider than 32 bits.

This engine feeds the generated `__upoly_random_raw_ir`; its output and
advanced state are stored in the concrete EDF `RetryTrace`.  Consequently the
strict FactorZp runtime now follows the same seeded random choices as the C++
entry, rather than merely choosing another semantically valid split.

## Coverage

The final-entry corpus now covers:

- an irreducible quadratic and the EDF early-return branch;
- a split linear DDF component and recursive EDF;
- a non-monic polynomial with multiplicity two;
- the derivative-zero p-th-root branch with multiplicity five;
- a degree-two EDF component split into two irreducible quadratics;
- `x^3 - x` over the prime `2^64 - 59`, exercising the 64-bit distribution
  and large-prime arithmetic path.

## Verification

- C++-probed uniform sequences and Lean sequences match for upper bounds
  `2`, `3`, `5`, `7`, `11`, `4294967311`, `9223372036854775783`, and
  `18446744073709551557`.
- Strict FactorZp final-entry C++/Lean B2B: 6/6 passed.
- The small and large ranges include rejection cases, so the comparison
  checks state advancement as well as individual mapped values.
- No production C++ source changed in this step.

With this runtime evidence the final FactorZp entry is executable across its
SQF, DDF and EDF branches.  The remaining end-to-end runtime task is the
FactorZZ outer entry, which reuses this exact FactorZp engine inside prime
selection before Hensel lifting and recombination.
