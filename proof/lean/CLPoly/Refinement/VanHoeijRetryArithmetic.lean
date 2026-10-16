import CLPoly.Generated.StrictVanHoeijRetry
import CLPoly.Refinement.VanHoeijWidth

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Normal growth has the exact source value, including its saturation branch. -/
theorem nextPrecision_success_exact (precision full next : Int32)
    (hpositive : 0 ≤ precision)
    (hrun : Generated.StrictVanHoeijRetry.nextPrecision precision full = .ok next) :
    next.toInt = (if 1073741823 < precision.toInt then 2147483647
      else max full.toInt (max 1 (2 * precision.toInt))) ∧
      precision.toInt < next.toInt ∧ full.toInt ≤ next.toInt := by
  have hg := Generated.StrictVanHoeijRetry.nextPrecision_increases precision full next hpositive hrun
  have hp := Int32.le_iff_toInt_le.mp hpositive
  change (0 : Int) ≤ precision.toInt at hp
  unfold Generated.StrictVanHoeijRetry.nextPrecision at hrun
  split at hrun
  next hmax => contradiction
  next hnotMax =>
    split at hrun
    next hlarge =>
      have hl : 1073741823 < precision.toInt := by
        change (2147483647 : Int32) / 2 < precision at hlarge
        exact Int32.lt_iff_toInt_lt.mp hlarge
      have he := Except.ok.inj hrun
      have hu := Int32.toInt_lt full
      rw [← he] at hg ⊢
      change (2147483647 : Int) = _ ∧ precision.toInt < 2147483647 ∧ full.toInt ≤ 2147483647
      rw [if_pos hl]
      exact ⟨rfl, hg, by omega⟩
    next hsmall =>
      have hs : ¬1073741823 < precision.toInt := by
        simpa only [Int32.lt_iff_toInt_lt] using hsmall
      have hdouble : (2 * precision).toInt = 2 * precision.toInt := by
        rw [Int32.toInt_mul]
        change (2 * precision.toInt).bmod (2 ^ 32) = 2 * precision.toInt
        apply Int.bmod_eq_of_le <;> omega
      have he := Except.ok.inj hrun
      rw [← he, int32_max_toInt, int32_max_toInt, hdouble, if_neg hs]
      exact ⟨rfl, by simpa only [← he, int32_max_toInt, hdouble] using hg, le_max_left _ _⟩

theorem nextPrecision_overflow (precision full : Int32)
    (hmax : precision = 2147483647) :
    Generated.StrictVanHoeijRetry.nextPrecision precision full = .error .arithmeticOverflow := by
  simp only [Generated.StrictVanHoeijRetry.nextPrecision, hmax, ↓reduceIte]

theorem nextPrecision_succeeds (precision full : Int32)
    (hmax : precision ≠ 2147483647) :
    ∃ next, Generated.StrictVanHoeijRetry.nextPrecision precision full = .ok next := by
  unfold Generated.StrictVanHoeijRetry.nextPrecision
  rw [if_neg hmax]
  split <;> exact ⟨_, rfl⟩

end Refinement.StrictRecombine
