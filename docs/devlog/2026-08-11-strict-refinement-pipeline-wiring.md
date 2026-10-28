# Wire strict generated refinement contracts into Pipeline

## Outcome

`CLPoly/Pipeline/L1.lean` now imports the centralized generated refinement
module and exposes checked Pipeline stage boundaries for the exact generated
C++ entries:

- `__squarefree_Zp_raw_ir` to `sqfZp`;
- `__ddf_Zp_raw_ir` to `ddf`;
- `__edf_Zp_raw_ir` to `EDFCorrect`, including the RNG transition;
- `__hensel_step_raw_ir` to `HenselStepCorrect`.

These Pipeline theorems delegate only to the public contracts generated in
`CLPoly.Refinement.Generated`; they do not call the L2 algorithms as an
execution fallback or manufacture an existence witness.

The strict SQF proof no longer imports the legacy `Generated.Corpus` module.
Two unused bridge lemmas tied to old Corpus entries were removed, and the
existing invariant library now aliases its measure directly to the strict
generated SQF control-flow measure.

## Verification

- `lake build` completed successfully with 3498 jobs.
- A zero-placeholder scan found no executable `sorry`, `admit`, or declared
  `axiom` anywhere under `proof/lean/CLPoly`.
- No CLPoly module imports `CLPoly.Generated.Corpus`.
- The strict generated/refinement/Pipeline closure contains no definition or
  parameter named as a fuel, oracle, or fallback implementation.
- `#print axioms` for all four generated contracts and all four Pipeline
  boundaries reports only `propext`, `Classical.choice`, and `Quot.sound`.

## Remaining work

The strict stages are now individually connected to Pipeline, but an
end-to-end generated execution composition is not yet exported.  That wrapper
must transport SQF output representation invariants into DDF, DDF components
into EDF, and the relevant integer-factor state into repeated Hensel steps;
it must not reuse the L2-only existence wrappers in the instantiate modules.
