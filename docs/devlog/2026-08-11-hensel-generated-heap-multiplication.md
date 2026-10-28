# Generate and refine Hensel heap multiplication

## Source boundary

`basic_polynomial::operator*` delegates to `pair_vec_multiplies`.  Its
observable heap loop repeatedly:

1. extracts the maximum monomial key;
2. executes every `addmul` contribution chained at that key;
3. appends the accumulated coefficient only when nonzero;
4. reinserts the advanced frontier nodes.

The strict generator now validates these anchors directly against
`clpoly/basic.hh`, as well as the two modular-division calls and multiplication
order in the actual `__hensel_step` source.  Generation fails if those source
anchors change.

The Lean L1 state erases allocation addresses and raw pointer storage, but it
retains every pending monomial/coefficient contribution and the observable
maximum-key extraction, `addmul`, zero-elimination, and frontier-removal
steps.  It is well founded on the number of pending contributions and has no
fuel argument.

## Refinement

The proof establishes:

- each heap-key coefficient is the sum of exactly its pending contributions;
- extracting a key partitions and preserves the complete frontier polynomial;
- the heap loop preserves the accumulator plus all remaining products;
- the initial Cartesian product frontier decodes to polynomial
  multiplication.

Consequently `__upoly_mul_raw_ir_refines` proves the generated raw output
decodes to `toPolyMod m a * toPolyMod m b`.

`strictHenselRawOps` now fixes the Hensel entry to the generated multiplication,
addition, and subtraction entries.  It does not use the `SparsePolyZZ` model
instances.

## Verification

```text
python3 proof/cpp2lean_v2/tests/build_strict_hensel.py --check
lake build CLPoly.Refinement.Hensel
Build completed successfully (1923 jobs).
```
