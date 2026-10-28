# Add the observable boundary for full factorization B2B

## What changed

The C++ B2B registry now exposes the two original source entries that form the
final refinement targets:

- `__factor_Zp`, including its leading coefficient, complete factor array,
  and multiplicities;
- `__factor_squarefree_primitive_ZZ`, including the source large-prime mode
  and complete integer-factor array.

Two shared observable encodings were added on both sides:
`FactorZpResult` and `ArraySPZZ`.  They preserve the existing typed nested
encoding for every polynomial, coefficient, prime, and multiplicity instead
of comparing a product checksum or a Boolean success flag.

The C++ FactorZZ registry branch temporarily sets the same global
`__g_use_large_prime` flag read by the real source entry and restores its
previous value on both normal and exceptional exits.  It calls the production
entry directly; no B2B-only factorization wrapper is involved.

The B2B Makefile now obtains GMP flags from `pkg-config` and, when available,
the Homebrew include prefix needed by the repository's Boost dependency.

## Verification

- forced clean-equivalent rebuild with `make -B test_types cpp_driver`:
  passed;
- C++ codec tests: all passed, including both complete result encodings;
- Lean codec tests via `lake env lean --run B2B/TestTypes.lean`: all passed;
- direct C++ driver execution of both final entries on `x^2 - 1`: passed and
  returned both physical factors in the complete encodings;
- `lake build B2B.Types B2B.TestTypes`: passed;
- `git diff --check`: passed.

## Remaining gate

This change establishes only the necessary observable C++ entry and shared
protocol.  It deliberately does not register a Lean-side implementation yet:
the latter must execute `Generated.StrictFactorZp.__factor_Zp_raw_ir` and
`Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir` with real
physical workspace and RNG implementations.  Calling the legacy generated
corpus, a high-level L2 factorizer, or replaying C++ outputs would not test the
strict L1 and is therefore not accepted as full pipeline B2B.
