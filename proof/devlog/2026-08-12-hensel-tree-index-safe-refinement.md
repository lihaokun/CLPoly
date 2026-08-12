# Hensel tree: index-safe generated refinement

The generated Hensel tree builder now has a build-checked refinement proof for
its concrete recursive construction path. The proof follows the generated L1
program: it computes interval products, executes the generated polynomial EEA,
stores the six node fields, appends child nodes, writes their raw `Int32`
pointers, and recurses left before right.

The recursive proof uses the source interval length `stop - start` as its
well-founded measure. The midpoint lemmas establish strict decrease for both
recursive calls; no fuel parameter or L2 fallback is used.

This checkpoint strengthens the array part of the proof in particular:

- every `getElem` result carries an explicit bound for the array version it
  indexes;
- every checked node setter preserves array size and frames all other nodes;
- the six stores for `g`, `h`, `s`, `t`, `leaf_start`, and `leaf_end` are
  composed in their actual generated execution order;
- recursive construction preserves prior nodes, establishes the exact node
  count, and identifies the root range and child pointers;
- pairwise-coprime factor hypotheses are connected to the polynomial GCD and
  Bezout invariant stored at every constructed internal node.

Verification:

- `lake build CLPoly.Refinement.Hensel` completed successfully (3326 jobs);
- the modified source contains no `sorry`, `admit`, or `axiom` declarations;
- the frame/store theorems use only `propext` and `Quot.sound`; the recursive
  existence theorem additionally uses `Classical.choice` and no project axiom;
- `git diff --check` passes.

The remaining Hensel work is the semantic topology bridge from the constructed
array to the lift/extraction traversal, followed by composition of the complete
generated Hensel entry theorem and its Pipeline contract.

## Follow-up: constructor-shaped extraction certificates

The topology bridge now has an explicit inductive certificate whose node
constructor records the actual `nodes[index]?` result, both raw child-pointer
matches, and certificates for every present child.  It converts directly to
the existing `HenselExtractInvariant` and is preserved only across genuine
array-prefix growth that leaves every old successful lookup unchanged.  This
is the transport mechanism needed when the generated builder finishes a left
subtree and then appends the disjoint right subtree.

The module again passes `lake build CLPoly.Refinement.Hensel`.  The certificate
conversion is axiom-free; prefix transport uses only `propext` and
`Quot.sound` and introduces no project axiom or placeholder.
