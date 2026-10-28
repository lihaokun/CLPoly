import CLPoly.Refinement.VanHoeijRowSeparation

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem truncateCld_zero_width_coeff_bound
    (original truncated : Array SparsePolyZZ) (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : Generated.StrictRecombine.cldSpiralWidth truncated = 0)
    (hcanonical : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (i : Nat) (hi : i < original.size) (degree : Nat) :
    |(SparsePolyZZ.toPoly original[i]!).coeff degree| < divisor := by
  have ht := truncateCld_spec original truncated divisor htruncate
  have htc : StrictPolynomialMod.SparsePolyZZCanonical truncated[i]! := by
    rw [ht.2.2 i hi]
    exact normalization_canonical _
  have hz := cldSpiralWidth_sparseCoeff_zero truncated i (by rw [ht.2.1]; exact hi)
    htc degree (by rw [hwidth]; omega)
  rw [truncateCld_sparseCoeff original truncated divisor htruncate i hi (hcanonical i hi),
    sparseCoeff_eq_coeff _ (hcanonical i hi)] at hz
  exact floor_zero_abs_lt _ divisor ht.1 hz

theorem rowResidualPolynomial_zero_width_unit (original : Array SparsePolyZZ)
    (i : Fin original.size) (modulus : Int) :
    rowResidualPolynomial original original.size 0
      (fun j => if j = i.val then 1 else 0) modulus = SparsePolyZZ.toPoly original[i.val]! := by
  simp only [rowResidualPolynomial, rowCarryPolynomial, Finset.range_zero,
    Finset.sum_empty, mul_zero, add_zero, rowCldCombination]
  classical
  simp only [Fin.is_lt, ↓reduceIte]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hji
    have hv : j.val ≠ i.val := fun h => hji (Fin.ext h)
    simp [hv]
  · simp

/-- Zero truncated width forces the full modular support of each genuine
irreducible divisor to contain at most one index at sufficient precision.
This uses no hypothetical LLL execution on the zero-width source branch. -/
theorem CldSeparationData.zero_width_local_indices_equal
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : Generated.StrictRecombine.cldSpiralWidth truncated = 0)
    (factor : Polynomial Int) (hfactor : factor ∣ SparsePolyZZ.toPoly f)
    (hprimitive : factor.IsPrimitive)
    (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : factor.natDegree ≠ 0)
    (left right : Fin original.size)
    (hleft : StrictHensel.toPolyMod modulus active[left.val]! ∣
      factor.map (Int.castRingHom (ZMod modulus)))
    (hright : StrictHensel.toPolyMod modulus active[right.val]! ∣
      factor.map (Int.castRingHom (ZMod modulus)))
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) active.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    left = right := by
  have hd := (truncateCld_spec original truncated divisor htruncate).1
  have hK : 1 ≤ separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree := by
    unfold separationCoordinateBound
    have := separationCutoff_nonneg active.size (SparsePolyZZ.toPoly f).natDegree
    have : 0 ≤ ((active.size + (SparsePolyZZ.toPoly f).natDegree : Nat) : Int) *
        2 ^ (active.size + (SparsePolyZZ.toPoly f).natDegree) *
        separationCutoff active.size (SparsePolyZZ.toPoly f).natDegree := by positivity
    omega
  have hD : divisor ≤ separationResidualHeight active.size (SparsePolyZZ.toPoly f).natDegree divisor := by
    unfold separationResidualHeight
    have hdim : (1 : Int) ≤ 2 + 3 * ((active.size + (SparsePolyZZ.toPoly f).natDegree : Nat) : Int) := by
      omega
    have hm := mul_le_mul hdim hK (by norm_num) (by positivity)
    nlinarith
  have he := data.scalars_equal_of_sufficient_precision truncated divisor htruncate
    original.size (fun j => if j = left.val then 1 else 0)
    (by intro i _; dsimp only; split <;> simp only [abs_one, abs_zero] <;> omega)
    (by
      intro degree
      rw [hwidth, rowResidualPolynomial_zero_width_unit original left]
      exact (truncateCld_zero_width_coeff_bound original truncated divisor htruncate hwidth
        (fun i hi => (data.values i hi).canonical) left.val left.isLt degree).le.trans hD)
    factor (Polynomial.ne_zero_of_natDegree_gt data.source_degree_pos) hfactor hprimitive
    hirreducible hfactorDegree left right left.isLt right.isLt hleft hright hresultant hcoordinates
  by_contra hne
  have hv : right.val ≠ left.val := fun h => hne (Fin.ext h.symm)
  simp [hv] at he

end Refinement.StrictRecombine
