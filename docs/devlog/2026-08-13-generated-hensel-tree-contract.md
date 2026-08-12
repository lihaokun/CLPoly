# Generated Hensel tree-builder contract

Pass 9 now emits the public theorem
`Refinement.__hensel_tree_build_raw_ir_refines`, preserving the original C++
double-underscore spelling.  The wrapper exposes the proved strict raw builder
contract: exact preorder size/root, concrete child pointers, the whole-tree
`HenselExtractInvariant`, and the initial algebraic root invariant.

`CLPoly.Pipeline.L1.strictHenselTreeBuildStage` consumes that generated public
contract.  It deliberately does not claim completion of the full
`__hensel_lift_upoly_ir` entry: the strict top-level lowering and composition
of build, lift loop, extraction, and normalization remain a separate proof
obligation.

Verification:

- Pass 9 regenerated `CLPoly/Refinement/Generated.lean`.
- The Python generator sources pass `py_compile` with their cache redirected
  to `/tmp`.
- Direct Lean checks of `Generated.lean` and `Pipeline/L1.lean` report no
  errors or `sorry`.
- `lake build CLPoly.Refinement.Generated CLPoly.Pipeline.L1` completed all
  3490 jobs successfully.
