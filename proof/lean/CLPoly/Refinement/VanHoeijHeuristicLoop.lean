import CLPoly.Generated.StrictVanHoeijHeuristic
import CLPoly.Refinement.VanHoeijWidth

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeijHeuristic

theorem clampPrecision_positive (value : Float) : 0 < clampPrecision value := by
  unfold clampPrecision
  split
  · decide
  · split
    · decide
    · rw [Int32.lt_iff_toInt_lt, Int32.toInt_zero, int32_max_toInt]
      change (0 : Int) < max 1 (Float.toInt32 value).toInt
      omega

theorem heuristic_increment_exact (exponent : Int32) (hnotMax : exponent ≠ 2147483647) :
    (exponent + 1).toInt = exponent.toInt + 1 := by
  apply Generated.StrictVanHoeij.int32_increment_toInt_of_lt exponent 2147483647
  apply Int32.lt_iff_toInt_lt.mpr
  have hu := Int32.toInt_lt exponent
  have hn : exponent.toInt ≠ 2147483647 := by
    intro he
    exact hnotMax (Int32.toInt_inj.mp he)
  change exponent.toInt < (2147483647 : Int)
  omega

theorem fullPrecisionLoop_rank_decreases (p : UInt64) (hp : 2 ≤ p.toNat)
    (target pa : ZZ) (hpa : 0 < pa) (hcontinue : pa ≤ target) :
    (target + 1 - pa * (p.toNat : Int)).toNat < (target + 1 - pa).toNat := by
  have hpInt : (2 : Int) ≤ p.toNat := by exact_mod_cast hp
  have hd := Int.mul_le_mul_of_nonneg_left (c := pa) hpInt (Int.le_of_lt hpa)
  change pa * 2 ≤ pa * (p.toNat : Int) at hd
  change (0 : Int) < pa at hpa
  change (pa : Int) ≤ target at hcontinue
  simp only [ZZ] at *
  rw [Int.toNat_lt_toNat (by omega : (0 : Int) < target + 1 - pa)]
  omega

theorem fullPrecisionLoop_output_ge (p : UInt64) (hp : 2 ≤ p.toNat)
    (target pa : ZZ) (hpa : 0 < pa) (exponent output : Int32)
    (hrun : fullPrecisionLoop p hp target pa hpa exponent = .ok output) :
    exponent.toInt ≤ output.toInt := by
  rw [fullPrecisionLoop] at hrun
  split at hrun
  next hcontinue =>
    split at hrun
    next hmax => contradiction
    next hnotMax =>
      have hn := fullPrecisionLoop_output_ge p hp target (pa * (p.toNat : Int))
        (Int.mul_pos hpa (by omega)) (exponent + 1) output hrun
      rw [heuristic_increment_exact exponent hnotMax] at hn
      omega
  next hstop =>
    have he := Except.ok.inj hrun
    rw [he]
termination_by (target + 1 - pa).toNat
decreasing_by exact fullPrecisionLoop_rank_decreases p hp target pa hpa (by assumption)

theorem fullPrecisionLoop_succeeds (p : UInt64) (hp : 2 ≤ p.toNat)
    (target pa : ZZ) (hpa : 0 < pa) (exponent : Int32)
    (hnonnegative : 0 ≤ exponent.toInt)
    (hpower : pa = (p.toNat : Int) ^ exponent.toNatClampNeg)
    (hfeasible : target < (p.toNat : Int) ^ 2147483647) :
    ∃ output, fullPrecisionLoop p hp target pa hpa exponent = .ok output := by
  rw [fullPrecisionLoop]
  split
  next hcontinue =>
    have hnotMax : exponent ≠ 2147483647 := by
      intro he
      rw [he] at hpower
      change pa = (p.toNat : Int) ^ 2147483647 at hpower
      simp only [ZZ] at hcontinue hfeasible hpower
      omega
    rw [if_neg hnotMax]
    have hi := heuristic_increment_exact exponent hnotMax
    have hn : (exponent + 1).toNatClampNeg = exponent.toNatClampNeg + 1 := by
      change (exponent + 1).toInt.toNat = exponent.toInt.toNat + 1
      rw [hi]
      omega
    apply fullPrecisionLoop_succeeds p hp target (pa * (p.toNat : Int)) (Int.mul_pos hpa (by omega))
      (exponent + 1) (by omega) ?_ hfeasible
    rw [hpower, hn, pow_succ]
  next hstop => exact ⟨exponent, rfl⟩
termination_by (target + 1 - pa).toNat
decreasing_by exact fullPrecisionLoop_rank_decreases p hp target pa hpa (by assumption)

theorem fullPrecision_succeeds (p : UInt64) (hp : 2 ≤ p.toNat) (target : ZZ)
    (hfeasible : target < (p.toNat : Int) ^ 2147483647) :
    ∃ output, fullPrecision p target = .ok output := by
  unfold fullPrecision
  rw [dif_pos hp]
  exact fullPrecisionLoop_succeeds p hp target 1 (by decide) 0 (by decide) (by simp) hfeasible

theorem fullPrecision_positive (p : UInt64) (hp : 2 ≤ p.toNat) (target : ZZ)
    (htarget : 1 ≤ target) (output : Int32)
    (hrun : fullPrecision p target = .ok output) : 0 < output := by
  unfold fullPrecision at hrun
  rw [dif_pos hp, fullPrecisionLoop, dif_pos htarget] at hrun
  have hzero : (0 : Int32) ≠ 2147483647 := by decide
  rw [if_neg hzero] at hrun
  have hg := fullPrecisionLoop_output_ge p hp target (1 * (p.toNat : Int))
    (Int.mul_pos (by decide) (by omega)) (0 + 1) output hrun
  apply Int32.lt_iff_toInt_lt.mpr
  change (0 : Int) < output.toInt
  change (1 : Int) ≤ output.toInt at hg
  omega

end Refinement.StrictRecombine
