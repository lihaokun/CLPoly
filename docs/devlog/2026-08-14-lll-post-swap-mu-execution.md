# Exact execution of the LLL post-swap mu correction

## Date

2026-08-14

## What changed

Defined and proved the exact array semantics of the generated
`updateMuAfterSwapLoop`.  The proof covers the per-row update, preservation of
outer and row sizes, unchanged rows before the starting index, and the exact
updated row returned for every row at or after that index.

A direct theorem now equates every successful generated `RawExec` result with
the strictly well-founded pure array computation and exposes the corrected row
formula.

## Why

The failed-Lovász branch swaps two adjacent basis rows and then repairs the two
corresponding mu columns below the pair.  Proving the concrete
`G = L D Lᵀ` invariant after that branch requires these actual generated array
assignments, not an assumed post-swap lower factor.

## Key decisions

The proof helper uses the same `mu.size - row` well-founded measure as the
generated loop.  It characterizes the generated result only; execution remains
on the generated L1 function.  No fuel, partial definition, or semantic oracle
is introduced.

## Problems and resolution

The two assignments overlap when expressed through total array reads, and the
second assignment must read the old value at `k` while using the newly computed
value for `k-1`.  A dedicated exact row function records this evaluation order,
after which the outer-loop induction updates each row once.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-14-lll-post-swap-mu-execution.md`
