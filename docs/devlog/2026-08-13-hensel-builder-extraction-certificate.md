# Hensel builder-to-extraction certificate

The generated C++-derived `__hensel_tree_build_recursive_raw_ir` proof now
returns a constructor-shaped certificate for every node allocated by the real
preorder builder.  The certificate is preserved across later subtree growth
and disjoint parent updates, so it covers the final concrete array rather than
an independently supplied mathematical tree.

`HenselTreeRawBuildCertificate.toExtractCertificate` interprets the stored raw
`Int32` child indices under the source `h_fits_int32` bound and derives the
exact recursive certificate consumed by the generated extraction traversal.
Consequently `strictHenselTreeBuildRawIR_refines_topology_root` now exports
`HenselExtractInvariant (henselTreeBuildTopology 0 factors.size 0) output` for
the actual result of `__hensel_tree_build_raw_ir`.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: no errors and no `sorry`.
- `lake build CLPoly.Refinement.Hensel`: completed successfully (3326 jobs).
- `git diff --check`: clean.
- `#print axioms` for the conversion, recursive builder theorem, and root
  topology theorem reports only `propext`, `Classical.choice`, and
  `Quot.sound`; there is no `sorryAx` or custom axiom.

This bridge uses neither fuel nor a specification oracle or L2 fallback.
