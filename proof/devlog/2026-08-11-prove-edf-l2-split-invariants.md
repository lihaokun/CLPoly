# Prove the L2 EDF split invariants

`edfPolynomialSplit_of_properDivisor` proves the mathematical core of the
recursive C++ split step. Given the actual proper monic divisor selected by
the candidate GCD, it identifies the normalized exact quotient and proves:

- their product is the original polynomial;
- both recursive inputs are monic, squarefree, and positive-degree;
- both degrees strictly decrease;
- every irreducible divisor still has the target equal degree.

The proof derives the quotient from `divByMonic`, its zero remainder from the
actual divisibility premise, and its monicity from the leading-coefficient
equation. It does not select factors through the L2 existence theorem.

Verification:

```text
lake build CLPoly.Refinement.EDF
```
