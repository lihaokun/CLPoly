# Exact execution of the LLL mu-prefix update

## Date

2026-08-14

## What changed

Defined the pure array value produced by the generated
`reduceMuPrefixLoop` and proved that the generated `RawExec` loop returns
exactly that value whenever its concrete source bounds hold.

Also proved its closed elementwise behavior: outer and row sizes are
preserved, all non-target rows are unchanged, and the target row is updated
exactly on `[index, limit)` by `mu[k,l] -= q * mu[source,l]`.  A direct
corollary exposes this formula for every successful generated execution.

## Why

`sizeReduceAt` changes an integer basis row and the corresponding row of the
Gram–Schmidt coefficient matrix.  Proving that `G = L D Lᵀ` is preserved
requires the exact generated mu update, not merely preservation of array sizes
or norms.

## Key decisions

The proof-side array function uses the same `limit - index` well-founded
measure as the generated loop.  It is only an exact value characterization;
the generated loop remains the executable L1 path.  No fuel, partial
definition, result witness, or semantic oracle is introduced.

## Problems and resolution

Each recursive update changes the target row while the source row must remain
available.  The induction therefore explicitly proves outer-array size, target
row size, and source row size preservation after every bounded `Array.set`.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-14-reduce-mu-prefix-execution.md`
