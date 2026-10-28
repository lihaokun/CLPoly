# Generate the strict Hensel public refinement contract

The quadratic `__hensel_step` proof now closes the exact generated raw C++
entry against `HenselStepCorrect`.  The proof composes the generated factor
and Bezout phases and requires an explicit refinement invariant containing
the source safety, degree-bound, divisibility, and input Hensel facts.  It
does not use a fuel interpreter, a specification oracle, or an L2 fallback.

Pass 9 now treats this checked theorem like the completed SQF, DDF, and EDF
contracts.  It emits the public theorem
`Refinement.__hensel_step_raw_ir_refines` into
`CLPoly/Refinement/Generated.lean`, while the proof implementation remains in
`CLPoly/Refinement/Hensel.lean`.  This makes the generated module the single
index of completed public L1-to-L2 refinement contracts.
