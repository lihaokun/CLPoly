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

The no-factor arm now has an exact `nextPrecision` transition and a derived
`precisionRank`.  A retry (`0 → J0`, then `J → 2J`) strictly decreases this
rank while the source bound `J ≤ J_max` holds; exceeding the bound returns the
fallback action.  The rank is erased termination evidence computed from the
real precision state, not an execution parameter or counter supplied to the
algorithm.

The two decreasing arms are now combined in generated `vanHoeijLoop`, whose
measure is `(active.size, precisionRank target initial maximum)`.  The loop
executes the real active gather, raw candidate preparation and validation,
reverse removal, precision transition, Zassenhaus fallback, and fallback
append operations.  The successful-extraction decrease is carried as erased
proof data around the very same `removeConsumed` result; no extra runtime
assertion was added to the C++ semantics.  `removalTermination` derives this
certificate from the concrete validator's one-consumed-bit-per-active-entry
size invariant and the already proved strict removal theorem.

Candidate validation now contains generated strict loops for its two index
scans.  `candidateAvailableLoop` rejects a candidate as soon as one of its
active-relative indices was consumed by an earlier successful candidate;
negative or out-of-range indices are raw faults.  `markConsumedLoop` performs
the actual indexed Boolean updates after successful trial division.  Their
refinement proofs establish successful execution for valid candidate indices
and preservation of the consumed-array size, which is the representation
fact needed to instantiate the main-loop removal certificate.

The source trial-product range loop is also generated.  For each candidate
index it calls a raw multiply/normalize/mod operation and threads the concrete
sparse product into the next iteration.  Its refinement theorem proves, by
the generated loop's `candidate.size - index` measure, that every successful
run decodes modulo `m` to the starting product times the ordered product of
the selected active lifted factors.  The remaining single-step contract must
now be discharged by strict integer polynomial multiplication together with
the existing strict coefficient-reduction theorem.
