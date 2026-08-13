# Strict van-Hoeij recombination bootstrap

`Refinement/Recombine.lean` is no longer an empty placeholder.  The first two
source range loops from `__vanhoeij_recombine` now have generated,
well-founded raw semantics and direct refinement theorems:

- `gatherActiveLoop` maps the current active-index array back into the
  original Hensel-lifted factor array.  Negative and out-of-range indices are
  explicit raw faults rather than `get!` defaults.
- `appendFallbackLoop` appends every factor returned by the actual
  Zassenhaus fallback to the accumulated result.

Both loops terminate by `array.size - index`.  Their proofs identify the
successful concrete arrays with the corresponding mapped/appended L2 lists.
They do not assume a recombination result, irreducibility, or an abstract UFD
factor witness.  The remaining work is the candidate-validation loop,
reverse active deletion, precision escalation, and the surrounding main
loop.

The reverse active-deletion loop is now generated as
`removeConsumedLoop`.  It follows the source order from `active.size - 1`
down to zero and performs the actual `eraseIdxIfInBounds` mutation whenever
the matching consumed bit is true.  Its recursion measure is the remaining
prefix length.  The refinement proof shows that equal active/consumed sizes
make every access safe, successful execution exists, and the resulting active
array cannot grow.  The next strengthening will use `found_any` to prove the
strict size decrease needed by the outer loop measure.

That strengthening is now proved: if at least one in-range consumed bit is
true, every successful execution of the actual reverse-erasure loop returns
an array with strictly smaller size.  The proof follows the same reverse
recursion, separating deletion of the current last index from preservation of
a marked bit in the remaining prefix.  This is the successful-extraction arm
of the outer van-Hoeij well-founded measure.
