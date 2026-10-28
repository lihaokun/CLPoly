import CLPoly.Generated.StrictVanHoeij

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem nextCldTarget_progress (target width : Int32)
    (ht : 0 < target) (hw : target < width) :
    target.toInt < (Generated.StrictVanHoeij.nextCldTarget target width).toInt ∧
      (Generated.StrictVanHoeij.nextCldTarget target width).toInt ≤ width.toInt := by
  have ht' := Int32.lt_iff_toInt_lt.mp ht
  change (0 : Int) < target.toInt at ht'
  have hw' := Int32.lt_iff_toInt_lt.mp hw
  have hwpos : 0 ≤ width.toInt := by omega
  have hhalf : (width / 2).toInt = width.toInt / 2 := by
    rw [Int32.toInt_div_of_ne_right width 2 (by decide)]
    change width.toInt.tdiv 2 = width.toInt / 2
    exact Int.tdiv_eq_ediv_of_nonneg hwpos
  unfold Generated.StrictVanHoeij.nextCldTarget
  split
  · exact ⟨hw', by omega⟩
  · rename_i hnot
    have hhalfBound : target.toInt ≤ width.toInt / 2 := by
      change ¬ width / 2 < target at hnot
      rw [Int32.lt_iff_toInt_lt, hhalf] at hnot
      omega
    have hlimit := Int32.toInt_lt width
    have hdouble : (2 * target).toInt = 2 * target.toInt := by
      rw [Int32.toInt_mul]
      change (2 * target.toInt).bmod (2 ^ 32) = 2 * target.toInt
      apply Int.bmod_eq_of_le <;> omega
    rw [hdouble]
    constructor <;> omega

theorem nextCldTarget_decreases_remaining (target width : Int32)
    (ht : 0 < target) (hw : target < width) :
    (width.toInt - (Generated.StrictVanHoeij.nextCldTarget target width).toInt).toNat <
      (width.toInt - target.toInt).toNat := by
  obtain ⟨hprogress, hbound⟩ := nextCldTarget_progress target width ht hw
  omega

end Refinement.StrictRecombine
