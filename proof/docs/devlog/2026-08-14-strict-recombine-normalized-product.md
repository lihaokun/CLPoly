# Strict recombination: normalized integer products

Commit scope: the concrete candidate-product multiplication used by strict
van Hoeij recombination.

- Replaced the model's opaque `Array.qsort` normalization tail with executable
  `List.mergeSort` followed by `toArray`.  This retains descending-degree
  sorting while exposing the standard library theorem `List.mergeSort_perm`.
- Proved that updating the coefficient of an equal-degree term preserves the
  represented polynomial after adding the new monomial.
- Proved the fold that groups equal degrees preserves the full polynomial sum.
- Proved filtering zero coefficients and sorting by degree preserve that sum.
- Combined those results with the already proved nested multiplication loops
  to show `multiplyNormalizeRaw` returns exactly the product of its two input
  integer polynomials.

The theorem is computational: it refers to the generated nested term loops and
the concrete normalization implementation.  It introduces no fuel, partial
definition, semantic callback, existence witness, axiom, or admitted proof.

Estimated progress after this commit: about 91% of the current FactorZZ goal.
Remaining work is candidate divisibility/validation and quotient semantics,
then the full van Hoeij/Zassenhaus result invariant, public generated FactorZZ
entry theorem, pipeline wiring, and final no-axiom/no-degeneration audit.
Estimated remaining focused time: 1.5--3 workdays; the divisibility/quotient
bridge is the main uncertainty.
