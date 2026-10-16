import CLPoly.Refinement.VanHoeijPartitionProductExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem primitiveRaw_content_of_positive_leading
    (input factor : SparsePolyZZ) (content : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hne : 0 < input.size) (hpositive : (0 : ZZ) < input[0].2)
    (hrun : Generated.StrictRecombine.primitiveRaw input = .ok (content, factor)) :
    content = (SparsePolyZZ.toPoly input).content := by
  have hempty : input.isEmpty = false := by simp [Array.isEmpty, Nat.ne_of_gt hne]
  unfold Generated.StrictRecombine.primitiveRaw at hrun
  simp only [hempty, Bool.false_eq_true, ↓reduceDIte] at hrun
  rw [if_neg (not_lt_of_ge hpositive.le)] at hrun
  cases hd : Generated.StrictRecombine.primitiveDivideLoop input
      (Generated.StrictRecombine.contentLoop input 0 0 : ZZ) 0 #[] with
  | error fault => simp [hd] at hrun
  | ok output =>
    rw [hd] at hrun
    have he := congrArg Prod.fst (Except.ok.inj hrun)
    exact he.symm.trans (contentLoop_zero_eq_content input hcanonical)

/-- Positive scaling of a positive-leading primitive target is removed by
the actual gcd/sign/primitive-division implementation. -/
theorem primitiveRaw_recovers_positive_scaled_target
    (input recovered : SparsePolyZZ) (content scalar : ZZ) (target : Polynomial Int)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hscalar : 0 < scalar) (hprimitive : target.IsPrimitive) (hleading : 0 < target.leadingCoeff)
    (hinput : SparsePolyZZ.toPoly input = Polynomial.C scalar * target)
    (hrun : Generated.StrictRecombine.primitiveRaw input = .ok (content, recovered)) :
    content = scalar ∧ SparsePolyZZ.toPoly recovered = target := by
  have ht : Polynomial.C scalar * target ≠ 0 :=
    mul_ne_zero (Polynomial.C_ne_zero.mpr (ne_of_gt hscalar)) hprimitive.ne_zero
  have hne : 0 < input.size := by
    by_contra hn
    have he : input = #[] := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hn)
    simp only [he, SparsePolyZZ.toPoly, Array.toList_empty, List.map_nil, List.sum_nil] at hinput
    exact ht hinput.symm
  have hpositive : (0 : ZZ) < input[0].2 := by
    rw [← sparsePolyZZ_leadingCoeff_eq_head input hcanonical hne, hinput,
      Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C]
    exact mul_pos hscalar hleading
  have hc := primitiveRaw_content_of_positive_leading input recovered content hcanonical hne hpositive hrun
  rw [hinput, Polynomial.content_C_mul, hprimitive.content_eq_one, mul_one,
    Int.normalize_of_nonneg hscalar.le] at hc
  refine ⟨hc, ?_⟩
  have he := primitiveRaw_toPoly input recovered content hrun
  rw [hinput, hc] at he
  exact (mul_left_cancel₀ (Polynomial.C_ne_zero.mpr (ne_of_gt hscalar)) he).symm

end Refinement.StrictRecombine
