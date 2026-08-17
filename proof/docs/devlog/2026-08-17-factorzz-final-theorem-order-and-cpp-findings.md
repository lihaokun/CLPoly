# FactorZZ final theorem ordering and C++ findings

## Goal

Make the strongest concrete FactorZZ refinement the conclusion of its
implementation module and expose the production C++ defects found during the
formal refinement effort from the Refinement landing page.

## Changes

- Moved
  `concreteSelectHensel___factor_squarefree_primitive_ZZ_raw_ir_refines` after
  every supporting lemma in `Refinement/FactorZZ.lean`.  It is now the final
  theorem before the namespace closes.
- Expanded its documentation to state that Pass 9 uses it as the proof behind
  the generated public FactorZZ contract.
- Added a documented C++ findings section to `Refinement/README.md`, covering
  the large-prime DDF/EDF narrowing error, Karatsuba workspace underallocation,
  van-Hoeij coordinate/safety defects, and unchecked Hensel factor-count
  narrowing, together with links to their detailed records.

This change only reorders an existing theorem and adds documentation.  It does
not alter the theorem statement, proof term, generated public wrapper, strict
L1 semantics, or production C++ source.
