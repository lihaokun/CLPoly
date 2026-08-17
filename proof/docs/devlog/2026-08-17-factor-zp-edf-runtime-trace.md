# Execute randomized EDF through concrete finite traces

## What changed

The FactorZp B2B runtime now executes EDF components that require splitting.
For each source retry it runs the generated random-polynomial operation and
the generated candidate operation, inspects the actual GCD result, and builds
the matching `RetryTrace.empty`, `RetryTrace.failed`, or
`RetryTrace.success` constructor.  The generated EDF then consumes that
finite trace structurally.  There is no fuel, retry bound, L2 factorizer, or
factor candidate supplied independently of the raw execution.

The new vector factors `x^2 - 1` over `F_5`.  This forces the DDF degree-one
component through the EDF retry and recursive split branches.  Both linear
factors and their multiplicities are compared with the C++ entry.

Because C++ performs only a degree comparison in its final `std::sort`, the
relative order of equal-degree factors is not part of the mathematical
result and can differ with a valid random trace.  Both B2B codecs now
canonicalize the complete encoded factor records before comparison.  No
factor or multiplicity is discarded.

## Verification

- forced C++ codec/driver rebuild: passed;
- C++ codec tests: passed;
- Lean B2B types, strict runtime, registry and driver builds: passed;
- final `__factor_Zp` B2B: 2/2 passed, including the EDF split case;
- the trace is assembled from executed raw random/candidate calls and is
  structurally consumed by `Generated.StrictEDF.__edf_Zp_raw_ir`.

## Follow-up

The following MT19937 commit replaces the temporary range-correct engine with
the literal source engine and platform distribution behavior.
