# Executable well-founded LLL rank

## What changed

The concrete determinant-based rank used by the generated strict LLL loop is
now computational.  The `noncomputable` markers were removed from the integer
Gram matrix, determinant potential, concrete LLL termination/execution bundle,
and the concrete van-Hoeij operation bundle.

The affected Gram--Schmidt norm proofs were also rebuilt from source.  Their
quadratic-form statement now ranges over every physical lattice column
(`matrix.size`), rather than incorrectly truncating the represented vector to
the number of selected rows.

## Why

`FactorZZRawOps` must execute the actual generated well-founded L1 during B2B
testing.  A noncomputable termination rank would make that execution
impossible, even though its proof fields are erased.  Lean's finite integer
matrix determinant is executable, so no fuel counter, partial definition, or
semantic replacement is necessary.

## Key decisions

- Retain the proved determinant/index lexicographic termination measure.
- Use Mathlib's executable finite determinant directly.
- Repair the physical-column dimension in the norm theorem instead of adding
  an assumption that would conceal the mismatch.
- Validate with a real `#eval` probe of both `concreteLLLRank` and the rank
  reached through `concreteVanHoeijRawOps`; the probe returned normally and was
  not retained as production source.

## Problems encountered

A direct source rebuild exposed errors in the recently added Gram--Schmidt
proof block that cached `.olean` files had hidden.  The block was repaired and
the target was then rebuilt through Lake successfully.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `docs/devlog/2026-08-17-executable-well-founded-lll-rank.md`
