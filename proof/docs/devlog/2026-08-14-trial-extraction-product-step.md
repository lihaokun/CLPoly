# Trial extraction product step

The successful integer trial-division branch is now summarized by a proved
algebraic step derived only from concrete raw executions.

If generated exact division returns an empty remainder and generated
primitive normalization of the quotient succeeds, then the old remaining
polynomial is exactly the product of the extracted factor, the new primitive
remainder, and the concrete quotient content.  This is the local product
identity required for the candidate-validation loop invariant.

The lemma deliberately makes no irreducibility claim.  Irreducibility must be
derived later from the Hensel-factor grouping and lattice/recombination
conditions rather than being hidden in a validation callback.

Estimated FactorZZ progress remains about 94%.  Product preservation through
the validation/main loops and the nontrivial irreducibility/fallback closure
remain before the generated FactorZZ theorem.  Current remaining estimate is
2--4 focused workdays; this supersedes the earlier shorter estimate because
the irreducibility bridge has now been separated from mere successful trial
division instead of being treated as implicit.
