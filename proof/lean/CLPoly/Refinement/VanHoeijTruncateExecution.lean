import CLPoly.Refinement.VanHoeijTruncation

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeij

theorem truncateTermsLoop_succeeds (poly : SparsePolyZZ) (divisor : ZZ)
    (index : Nat) (hd : divisor ≠ 0) :
    ∃ output, truncateTermsLoop poly divisor index = .ok output := by
  induction hmeasure : poly.size - index using Nat.strong_induction_on
      generalizing poly index with
  | h measure ih =>
    rw [truncateTermsLoop]
    split
    next hi =>
      exact ih _ (by simp only [Array.size_set]; omega) _ (index + 1) rfl
    next hi => exact ⟨poly, rfl⟩

theorem truncateCldLoop_succeeds (cld : Array SparsePolyZZ) (divisor : ZZ)
    (index : Nat) (hd : divisor ≠ 0) :
    ∃ output, truncateCldLoop cld divisor index = .ok output := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on
      generalizing cld index with
  | h measure ih =>
    rw [truncateCldLoop]
    split
    next hi =>
      obtain ⟨poly, hp⟩ := truncateTermsLoop_succeeds cld[index] divisor 0 hd
      rw [hp]
      exact ih _ (by simp only [Array.size_set]; omega) _ (index + 1) rfl
    next hi => exact ⟨cld, rfl⟩

theorem truncateCld_succeeds (cld : Array SparsePolyZZ) (divisor : ZZ)
    (hd : 0 < divisor) : ∃ output, truncateCld cld divisor = .ok output := by
  unfold truncateCld
  rw [if_pos hd]
  exact truncateCldLoop_succeeds cld divisor 0 (ne_of_gt hd)

end Refinement.StrictRecombine
