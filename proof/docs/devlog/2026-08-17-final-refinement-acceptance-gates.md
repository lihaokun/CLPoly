# Close the final strict refinement acceptance gates

## What changed

The final audit found and corrected three test-infrastructure defects that
could otherwise hide regressions:

1. The DDF strict generator still expected the retired per-component
   `Refinement/Generated/DDF.lean` file.  DDF now follows every other strict
   component: it generates only its strict L1 artifact, while Pass 9 emits the
   public theorem into the single `Refinement/Generated.lean` module.
2. The B2B pipeline smoke request intentionally calls an unknown function,
   but the runner treated matching C++ and Lean errors as a failure.  Expected
   errors are now asserted explicitly, so the full suite can be genuinely
   green rather than carrying an explained red case.
3. The cpp2lean smoke fixtures did not add the active Homebrew include prefix
   on macOS.  Clang therefore recovered from missing Boost headers with
   incorrect `int` types.  The fixture include surface was narrowed and the
   package include path was added.  The Pass 6 report formatter was also made
   aware that `Call.callee` may be an SSA `Var`, as required by the translator
   IR.  Finally, `pipefail` prevents a failing smoke stage from being hidden by
   its output-filtering `tail` process.

These changes affect generators, tests, and generated survey reports only.
No production C++ algorithm changed.

## Final proof-boundary evidence

- Pass 9 reproducibility for `CLPoly/Refinement/Generated.lean`: passed.
- Strict SQF, DDF, EDF, GCD, FactorZp, select-prime, Hensel, recombination,
  and FactorZZ generator checks: passed and placeholder-free.
- Strict-boundary audit: passed; the formal boundary contains no legacy
  Corpus dispatch, `partial def`, `sorry`, fuel interface, or semantic
  recombination fallback.
- B2B runtime reverse-dependency audit: passed; unsafe test-only witnesses
  cannot enter `Generated`, `Refinement`, or `Pipeline`.
- Public final theorem axiom audit: both FactorZp and FactorZZ depend only on
  `propext`, `Classical.choice`, and `Quot.sound`.
- Full `lake build`: 3559 jobs, zero errors.

The centralized public theorem names are:

- `Refinement.__factor_Zp_raw_ir_refines_FactorZpCorrect`;
- `Refinement.__factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect`.

`CLPoly/Pipeline/L1.lean` re-exports both strict generated entries directly.

## Executable regression evidence

- all C++/Lean B2B vectors: 158/158 passed, including 6 final FactorZp and 5
  final FactorZZ cases;
- cpp2lean v2 unit and Pass 1--7 smoke suite: passed with failure propagation
  enabled;
- B2B codec tests: passed;
- deterministic MT19937 and UInt64 primality runtime tests: passed;
- freshly compiled native `test_recombine`: 121/121 passed.

The repository's aggregate native Make target remains affected by its
pre-existing macOS `realroot.hh` full-library build issue.  The relevant
factorization regression was therefore freshly compiled from
`test/test_recombine.cc`, `clpoly/upolynomial.cc`, and `clpoly/variable.cc`
with C++17 and the active Homebrew Boost/GMP include and link paths; it did not
reuse a stale binary.

## Files

- `proof/cpp2lean_v2/tests/build_strict_ddf.py`
- `proof/cpp2lean_v2/tests/test_pass1_parse.py`
- `proof/cpp2lean_v2/tests/smoke_pass2_full.py`
- `proof/cpp2lean_v2/tests/smoke_pass6_full.py`
- `proof/cpp2lean_v2/tests/run_all_smoke.sh`
- `proof/cpp2lean_v2/tests/fixtures/upoly_mod.cc`
- `proof/cpp2lean_v2/tests/fixtures/upoly_divmod.cc`
- `proof/cpp2lean_v2/tests/fixtures/upoly_const_term.cc`
- `proof/b2b/runner/diff.py`
- `proof/b2b/vectors/_pipeline_smoke.json`
- `docs/design/l1-translation-validation/survey/pass1-smoke.md`
- `docs/design/l1-translation-validation/survey/pass2-smoke.md`
- `docs/design/l1-translation-validation/survey/pass3-smoke.md`
- `docs/design/l1-translation-validation/survey/pass4-smoke.md`
- `docs/design/l1-translation-validation/survey/pass5-smoke.md`
- `docs/design/l1-translation-validation/survey/pass6-smoke.md`
- `docs/design/l1-translation-validation/survey/pass7-smoke.md`
- `proof/docs/devlog/2026-08-17-final-refinement-acceptance-gates.md`
