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

At this checkpoint the executor rejected an EDF component that required
randomized splitting.  This was intentional and observable as
`assertionFailure`; it was not a fallback.  The following runtime-trace commit
removes that limitation.

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

Implement runtime retry traces and then the exact C++
`std::mt19937(42)`/uniform-distribution stream before using the same physical
runtime for the full FactorZZ entry.
