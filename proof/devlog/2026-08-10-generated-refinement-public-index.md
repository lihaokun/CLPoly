# Add a public index for completed generated refinement contracts

The completed cpp2lean L1→L2 refinement contracts now have one import entry:
`CLPoly.Refinement.Generated`.

The index currently exposes SQF and DDF. A stage is added only after its strict
top-level refinement theorem has been proved; in particular, unfinished EDF
intermediate lemmas are not presented as a final contract.

The directory roles remain deliberately separate:

- `CLPoly/Generated`: generated executable L1 semantics;
- `CLPoly/Refinement`: detailed proof implementation;
- `CLPoly/Refinement/Generated`: generated, public final theorem contracts.

Files changed:

- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/lean/CLPoly.lean`
- `proof/devlog/2026-08-10-generated-refinement-public-index.md`
