# Strict `__factor_squarefree_primitive_ZZ` generated entry

The generated module `CLPoly.Generated.StrictFactorZZ` replaces the legacy
`partial` aggregate-corpus definitions for the two top-level C++ functions
used by squarefree primitive integer factorization.

The generated `__lll_factorize_raw_ir` preserves both C++ phases:

1. compute heuristic and full Mignotte precisions;
2. Hensel-lift to the heuristic precision and run van Hoeij recombination;
3. when the result is incomplete and the heuristic precision was lower, run
   full-precision Hensel lifting followed by a second recombination.

The generated `__factor_squarefree_primitive_ZZ_raw_ir` preserves the source
precondition guard, prime-selection call, irreducible/single-factor branch,
and the call to the strict LLL factorization entry.  Its operation interface
contains only effectful raw calls.  In particular, it contains no L2
factorization function, existence witness, fuel counter, or semantic oracle.

The next proof layer must instantiate the four raw boundaries with concrete
strict refinements for prime selection, heuristic precision, Hensel lifting,
and van Hoeij recombination before proving `FactorZZCorrect`.
