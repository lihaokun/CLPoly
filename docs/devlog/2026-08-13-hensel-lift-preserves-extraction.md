# Hensel lift preserves the builder extraction certificate

The strict generated Hensel step updates only polynomial payloads.  A theorem
derived from its successful raw execution equation now proves preservation of
`left`, `right`, `leaf_start`, and `leaf_end`.  Array-level topology equality
then composes the root `set!`, left subtree run, and right subtree run.

Induction over the actual semantic traces proves topology preservation for:

- every `__hensel_lift_recursive_raw_ir` traversal; and
- every iteration of the well-founded `__hensel_lift_loop_raw_ir` loop.

The constructor-shaped builder/extraction certificate is transported across
that topology equality.  `henselLiftLoopThenExtractRawIR_refines` consequently
executes the real generated extraction traversal on the real lifted array and
returns its `HenselExtractCorrect` trace.  No replacement tree, factor array,
fuel counter, specification oracle, or L2 fallback is involved.

Verification:

- Direct Lean checking reports no errors or `sorry`.
- `lake build CLPoly.Refinement.Hensel` completed successfully (3326 jobs).
- `#print axioms` for the raw-step preservation, recursive and loop
  preservation, certificate transport, and lift-to-extract composition reports
  only `propext`, `Classical.choice`, and `Quot.sound` as applicable; no
  `sorryAx` or custom axiom appears.
