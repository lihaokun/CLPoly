# Execute the strict FactorZp pipeline at the B2B boundary

## What changed

The Lean B2B driver can now execute the generated
`Generated.StrictFactorZp.__factor_Zp_raw_ir` entry for inputs whose EDF
components take the source degree-equals-`d` branch.  Its operation record is
assembled directly from the checked raw make-monic, SQF, DDF, EDF, arithmetic,
GCD, division and modular-reduction implementations.  It does not call an L2
factorizer, replay a C++ result, or use the legacy generated corpus.

The initial cross-language vector is the irreducible polynomial
`x^2 + 2` over `F_5`.  It executes the complete outer sequence
make-monic → SQF → DDF → EDF → multiplicity attachment → degree sort.  Both
drivers return the same complete `FactorZpResult` encoding.

## Physical runtime boundary

The formal refinement continues to quantify over checked physical workspace
contracts.  The test runtime supplies concrete raw heaps, pointer regions,
capacities and operation records.  It uses an unsafe erased inhabitant only
for compiler-erased proof fields.  A new reverse-dependency audit rejects any
import or reference from `Generated`, `Refinement`, or `Pipeline` back into
this B2B-only module.

The executor currently rejects an EDF component that requires randomized
splitting.  This is intentional and observable as `assertionFailure`; it is
not a fallback.  General FactorZp B2B still requires the exact C++ MT19937
draw stream and a concrete finite `RetryTrace`.  Therefore this commit is an
execution checkpoint, not completion of full FactorZp B2B.

## Verification

- `lake build B2B.Registry B2B.Driver`: passed.
- C++/Lean B2B for `__factor_Zp` irreducible quadratic: 1/1 passed.
- All effective B2B vectors: 147/147 passed.  The only runner failure is the
  deliberate `_pipeline_smoke.json` unknown-function case, rejected by both
  drivers.
- strict refinement boundary audit: passed.
- B2B strict-runtime reverse-dependency audit: passed.
- `git diff --check`: passed.
- Production C++ source changes: none, so this checkpoint creates no new
  native regression surface.

## Remaining work

Implement the exact C++ `std::mt19937(42)`/uniform-distribution execution and
construct the corresponding finite retry traces from actual candidate runs;
then add split/retry FactorZp vectors before using the same physical runtime
for the full FactorZZ entry.
