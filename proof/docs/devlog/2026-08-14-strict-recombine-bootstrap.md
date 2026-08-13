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
