import CLPoly.Refinement.VanHoeijExecution

set_option autoImplicit false

namespace Refinement.StrictHensel
open _root_.Generated.StrictHensel

theorem pushNonzero_fdiv_coefficientsReduced (m : Nat) (hm : 0 < m)
    (result : SparsePolyZZ) (degree : Nat) (value : Int)
    (hr : CoefficientsReduced m result) :
    CoefficientsReduced m (pushNonzero result degree (ZZ.fdiv_r value value (m : Int))) := by
  intro term ht
  unfold pushNonzero at ht
  split at ht
  · simp only [Array.toList_push, List.mem_append, List.mem_singleton] at ht
    rcases ht with ht | rfl
    · exact hr term ht
    · exact ⟨ZZ.fdiv_r_nonneg value value (m : Int) (by exact_mod_cast hm),
        ZZ.fdiv_r_lt value value (m : Int) (by exact_mod_cast hm)⟩
  · exact hr term ht

theorem divmodMergeLoop_coefficientsReduced (m : Nat) (hm : 0 < m)
    (r g : SparsePolyZZ) (coefficient : Int) (degreeShift rIndex gIndex : Nat)
    (result : SparsePolyZZ) (hr : CoefficientsReduced m result) :
    CoefficientsReduced m (divmodMergeLoop r g coefficient (m : Int) degreeShift rIndex gIndex result) := by
  induction hmeasure : (r.size - rIndex) + (g.size - gIndex) using Nat.strong_induction_on
      generalizing rIndex gIndex result with
  | h measure ih =>
    rw [divmodMergeLoop]
    split
    next hmore =>
      split
      · exact ih _ (by omega) _ _ _ (pushNonzero_fdiv_coefficientsReduced m hm _ _ _ hr) rfl
      · split
        · exact ih _ (by omega) _ _ _ (pushNonzero_fdiv_coefficientsReduced m hm _ _ _ hr) rfl
        · dsimp only
          split
          · exact ih _ (by omega) _ _ _ (pushNonzero_fdiv_coefficientsReduced m hm _ _ _ hr) rfl
          · split <;>
              exact ih _ (by omega) _ _ _ (pushNonzero_fdiv_coefficientsReduced m hm _ _ _ hr) rfl
    next hdone => exact hr

theorem divmodLoop_remainder_coefficientsReduced (m : Nat) (hm : 0 < m)
    (g : SparsePolyZZ) (inverse : Int) {r q : SparsePolyZZ}
    (trace : DivmodTrace g inverse (m : Int) r q) (hr : CoefficientsReduced m r) :
    CoefficientsReduced m (divmodLoop g inverse (m : Int) trace).2 := by
  induction trace with
  | done r q _ => exact hr
  | vanished r q _ _ next ih =>
    apply ih
    intro term ht
    rw [Array.toList_eraseIdxIfInBounds] at ht
    exact hr term (List.mem_of_mem_eraseIdx ht)
  | subtract r q _ _ next ih =>
    apply ih
    exact divmodMergeLoop_coefficientsReduced m hm r g _ _ 1 1 #[]
      (by intro term ht; simp at ht)

theorem divmod_raw_remainder_coefficientsReduced
    (termination : DivmodTermination) (f g q r : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (hrun : __upoly_divmod_mod_raw_ir termination f g (m : Int) = .ok (q, r)) :
    CoefficientsReduced m r := by
  simp only [__upoly_divmod_mod_raw_ir, __upoly_mod_coeff_raw_ir] at hrun
  split at hrun <;> try contradiction
  next hnotEmpty =>
    split at hrun <;> try contradiction
    next hvalid =>
      split at hrun <;> try contradiction
      next hinvertible =>
        have hr := divmodLoop_remainder_coefficientsReduced m hm g
          (ZZ.invert 0 g[0]!.2 (m : Int)).2
          (termination.trace f g (modCoeffOutput f (m : Int)) (m : Int) rfl hvalid)
          (modCoeffOutput_coefficientsReduced f m hm)
        have he := Except.ok.inj hrun
        rw [he] at hr
        exact hr

end Refinement.StrictHensel
