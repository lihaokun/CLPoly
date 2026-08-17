# Add the strict EDF odd-candidate raw execution bridge

The odd-characteristic branch of generated `candidateRun` now has an exact
raw-to-safe execution theorem. The proof follows the generated match chain in
source order:

1. `ops.powmod` succeeds with `hpow`;
2. generated `__upoly_subtract_one_raw_ir` succeeds with `hminus`;
3. `ops.gcd` succeeds with `factor`;
4. generated `candidateRun` returns that same `factor`.

Every intermediate value is fixed by a concrete `RawExec = .ok ...` equation;
there is no L2 witness, specification oracle, or fallback result.

Files changed:

- `proof/lean/CLPoly/Refinement/EDF.lean`
