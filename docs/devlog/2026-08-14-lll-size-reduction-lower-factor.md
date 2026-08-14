# LLL size reduction as a concrete lower-factor row operation

## Date

2026-08-14

## What changed

Defined the exact mu result of one generated `sizeReduceAt` coefficient and
proved that every successful generated execution returns that concrete array.
Its target row and all non-target rows were characterized elementwise.

Lifted those array results to the finite-matrix identity
`L' = E(k,j,-q) * L` for every prefix containing row `k`, where `E` is the
same transvection used by the actual integer basis-row subtraction.  Prefixes
before row `k` are proved unchanged.

## Why

The generated LLL loop may be used in a well-founded recursion only after the
relationship between its actual basis, mu coefficients, and squared norms is
preserved by every concrete step.  This identity supplies the mu side of the
`G = L D Lᵀ` preservation proof for size reduction.

## Key decisions

The coefficient is obtained from the successful generated `roundQQ` call, and
the result is unfolded through the generated `reduceMuPrefixLoop`.  The proof
does not assume an arbitrary updated lower factor or use a semantic oracle.

## Problems and resolution

The diagonal case `l=j` is distinct from the strict-prefix update because the
unit lower factor has `L[j,j]=1`.  The proof handles the five concrete column
positions separately and shows that the explicit final subtraction by `q`
exactly implements that diagonal contribution.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-14-lll-size-reduction-lower-factor.md`
