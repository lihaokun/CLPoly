# Concrete Gram–Schmidt execution invariant

## Date

2026-08-14

## What changed

Added an explicit rational `G = L D Lᵀ` invariant for the generated LLL
state.  `G` is the Gram matrix of the actual integer basis rows, `L` is the
unit lower-triangular matrix read from the generated `mu` array, and `D` is the
diagonal matrix read from the generated squared norms.

Proved that this invariant implies the exact prefix Gram determinant identity
used by the natural-number LLL termination rank.

Added `ConcreteLLLExecutionValid` as the loop-facing invariant.  It combines
the concrete storage/positivity conditions with the factorization and derives
the older determinant-validity view through
`toConcreteLLLValid`.  Changing only the loop index is proved to preserve this
full invariant.

The execution invariant now also records that the generated `mu` storage is
exactly square with the basis dimension.  This reflects `zeroQQMatrix size
size` in `initializeLLL` and supplies every source bound needed by the concrete
size-reduction loop; no bounds are inferred from an abstract successful result.

## Why

Prefix determinants alone do not determine the boundary prefix created by an
adjacent row swap.  The stronger factorization records the missing relationship
between the real basis matrix and the generated Gram–Schmidt arrays, allowing
the swap-validity proof to be derived instead of assumed.

## Key decisions

The invariant is an equality over concrete finite matrices and generated array
reads.  It contains neither an output witness nor an execution oracle.  The
unit lower-triangular determinant is proved equal to one, so taking determinants
of the factorization yields exactly the existing prefix norm product.

## Problems and resolution

The determinant bridge required mapping the integer Gram matrix into the
rationals and using `Int.cast_det`.  The lower-triangular determinant was
discharged with Mathlib's block-triangular determinant theorem over `Fin`.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-14-concrete-gram-schmidt-invariant.md`
