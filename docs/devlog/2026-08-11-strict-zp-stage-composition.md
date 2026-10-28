# Compose strict SQF, DDF, and EDF stages

## Outcome

The strict Zp stages now transport concrete representation invariants across
their real generated outputs:

- the SQF refinement returns canonicality for every sparse factor stored in
  its actual output array;
- `strictSQFStage_preparesDDF` proves the prime tag, signed degree bound,
  monicity, and squarefreeness required by DDF for every such array member;
- `strictSQFStage_runsDDF` executes the generated DDF entry on each unchanged
  SQF output factor;
- the DDF loop invariant now records canonicality, positive polynomial
  degree, and positive degree-class labels for all accumulated raw outputs;
- `strictDDFStage_preparesEDF` constructs the full `EDFEntryInvariant` for
  every concrete DDF output component;
- `strictDDFStage_runsEDF` executes the generated EDF entry on each unchanged
  DDF component and preserves the concrete RNG transition.

All output witnesses come from successful generated raw executions.  The
proofs use L2 correctness only to establish mathematical properties of those
same concrete arrays; they do not replace an array, select a specification
witness, or call an L2 algorithm as executable fallback.

## Verification

The modified strict SQF, DDF, generated public-contract, and Pipeline modules
all pass direct Lean checking.  The centralized generator reproduces the
strengthened SQF and DDF public contracts.

## Remaining work

The Zp path is now composed through EDF at the per-component execution level.
The remaining Pipeline work is the repeated Hensel-step composition under a
well-founded precision measure, followed by the final repository-wide build,
zero-placeholder, degeneration, and axiom audits.
