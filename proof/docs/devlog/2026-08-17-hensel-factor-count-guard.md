# Hensel factor-count narrowing guard

## What changed

The strict Hensel refinement exposed an unchecked conversion from the
`size_t` modular-factor count to the `int` indices used by the C++ Hensel
tree.  Both production entry points now reject counts above `INT_MAX` before
the conversion.  The generated strict Lean L1 executes the same guard and
returns `RawFault.arithmeticOverflow` on the rejected branch.

The Hensel refinement proves that any successful strict execution passed the
guard.  FactorZZ no longer receives `factorCountFits` as a fact attached to
prime selection; it derives the fact from the actual successful Hensel run
before using machine-index comparison lemmas.

## Why

Without the check, a factor count above `INT_MAX` could wrap or truncate when
passed to `__hensel_tree_build_recursive`, invalidating its half-open ranges
and node indices.  More importantly for source refinement, Lean could not
justify the source narrowing from successful execution because the C++
program did not enforce the required range.

## Verification

- `test_recombine`: 121/121 native tests pass, including direct checks at
  `INT_MAX` and `INT_MAX + 1`.
- C++/Lean B2B `__hensel_factor_count_fits`: 5/5 cases pass at `0`, `2`,
  `INT_MAX`, `INT_MAX + 1`, and `UINT64_MAX`.
- `build_strict_hensel.py --check` passes and now validates the new C++
  source anchors.
- The Hensel and FactorZZ Lean modules are compiled as part of this change's
  proof gate.

## Key decision

The regression tests call the range predicate directly, so the exact
overflow boundary is tested without attempting to allocate an impossible
array.  End-to-end Hensel behavior remains covered by the existing native
recombination suite.

## Files

- `clpoly/polynomial_factorize_univar.hh`
- `test/test_recombine.cc`
- `proof/b2b/registry/cpp_registry.cc`
- `proof/b2b/vectors/__hensel_factor_count_fits.json`
- `proof/lean/B2B/Registry.lean`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
