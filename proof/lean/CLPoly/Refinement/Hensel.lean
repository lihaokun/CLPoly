/-
  CLPoly/Refinement/Hensel.lean — generated Hensel candidate and certified L2 boundary.

  The generated routine returns integer representatives together with an actual
  lifting modulus.  The observable L2 result is obtained by mapping those
  representatives to `ZMod (p^k)`, not by equating two integer-polynomial lists.
-/
import CLPoly.Generated.HenselSafe
import CLPoly.Pipeline.FactorZZInstantiate
import CLPoly.Refinement.Basic
import CLPoly.Math.Univariate

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- Map the integer representatives returned by generated Hensel lifting into
the ring required by `HenselCorrect`. -/
noncomputable def henselCandidateToPk (p k : ℕ)
    (candidate : Array SparsePolyZZ × Int) :
    List (Polynomial (ZMod (p ^ k))) :=
  (candidate.1.map fun factor =>
    Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (SparsePolyZZ.toPoly factor)).toList

/-- The actual generated Hensel candidate for explicit exponent `k`.  The caller
guards the narrowing conversion to `Int32`. -/
noncomputable def henselGeneratedCandidate (p k : ℕ) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) : List (Polynomial (ZMod (p ^ k))) :=
  henselCandidateToPk p k
    (Generated.__hensel_lift_upoly_ir f factors (UInt64.ofNat p) (k : Int32))

/-- Total certified Hensel boundary.  A generated candidate is adopted exactly
when it satisfies the complete L2 contract; narrowing overflow or a rejected
candidate uses the already-proved multifactor lift. -/
noncomputable def henselZpSafe (p : ℕ) [Fact (Nat.Prime p)]
    (k : ℕ) (hk : 0 < k) (f : Polynomial ℤ) (fSparse : SparsePolyZZ)
    (factors : List (Polynomial (ZMod p))) (factorsSparse : Array SparsePolyZp) :
    List (Polynomial (ZMod (p ^ k))) := by
  classical
  if hfit : k ≤ 2147483647 then
    if hcount : 2 ≤ factorsSparse.size then
      let candidate := henselGeneratedCandidate p k fSparse factorsSparse
      exact if HenselCorrect f k factors candidate then candidate
        else hensel_lift p k hk f factors
    else
      exact hensel_lift p k hk f factors
  else
    exact hensel_lift p k hk f factors

/-- Correctness of the generated-candidate/verified-fallback boundary. -/
theorem henselZpSafe_correct (p : ℕ) [Fact (Nat.Prime p)]
    (k : ℕ) (hk : 0 < k) (f : Polynomial ℤ) (fSparse : SparsePolyZZ)
    (factors : List (Polynomial (ZMod p))) (factorsSparse : Array SparsePolyZp)
    (hne : factors ≠ [])
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = factors.prod)
    (hcop : factors.Pairwise fun a b => IsCoprime a b) :
    HenselCorrect f k factors
      (henselZpSafe p k hk f fSparse factors factorsSparse) := by
  classical
  unfold henselZpSafe
  split
  · split
    · dsimp only
      split
      · assumption
      · exact hensel_lift_correct p k hk f factors hne hprod hcop
    · exact hensel_lift_correct p k hk f factors hne hprod hcop
  · exact hensel_lift_correct p k hk f factors hne hprod hcop

/-- A generated candidate that passes the contract is adopted verbatim. -/
theorem henselZpSafe_of_generated_correct (p : ℕ) [Fact (Nat.Prime p)]
    (k : ℕ) (hk : 0 < k) (f : Polynomial ℤ) (fSparse : SparsePolyZZ)
    (factors : List (Polynomial (ZMod p))) (factorsSparse : Array SparsePolyZp)
    (hfit : k ≤ 2147483647) (hcount : 2 ≤ factorsSparse.size)
    (hcandidate : HenselCorrect f k factors
      (henselGeneratedCandidate p k fSparse factorsSparse)) :
    henselZpSafe p k hk f fSparse factors factorsSparse =
      henselGeneratedCandidate p k fSparse factorsSparse := by
  classical
  unfold henselZpSafe
  rw [dif_pos hfit, dif_pos hcount]
  dsimp only
  exact if_pos hcandidate

end Refinement
