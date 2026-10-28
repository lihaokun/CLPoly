# Centralize generated refinement contracts

Generated executable models remain under `CLPoly/Generated`, while generated
public L1→L2 theorem contracts now live under `CLPoly/Refinement/Generated`.
This separates executable artifacts from proof-facing API files.

The DDF contract moved from
`CLPoly/Generated/StrictDDFRefinement.lean` to
`CLPoly/Refinement/Generated/DDF.lean`.  Its public name is now the uniform
root-namespace declaration
`Refinement.__ddf_Zp_raw_ir_refines_ddf`, matching the existing generated SQF
contract layout.  The strict implementation lemmas remain organized under
`Refinement.StrictDDF`.

The DDF generator writes and reproducibly checks the new location.  Both
centralized public contracts build together:

```text
python3 ../cpp2lean_v2/tests/build_strict_ddf.py --check
lake build CLPoly.Refinement.Generated.DDF \
  CLPoly.Refinement.Generated.SquarefreeZp
```

Future EDF and Hensel public contracts will use the same directory and root
namespace convention.
