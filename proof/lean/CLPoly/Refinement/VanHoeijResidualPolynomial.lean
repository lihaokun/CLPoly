import CLPoly.Refinement.VanHoeijColumnCoverage
import CLPoly.Refinement.VanHoeijCld
import Mathlib.Data.ZMod.Basic

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Original (untruncated) CLD combination, using exactly the physical
transform's factor coordinates; modulus-row coordinates contribute zero. -/
noncomputable def rowCldCombination (original : Array SparsePolyZZ) (size : Nat)
    (coefficients : Nat → Int) : Polynomial Int :=
  ∑ i : Fin size, if i.val < original.size then
    Polynomial.C (coefficients i.val) * SparsePolyZZ.toPoly original[i.val]! else 0

/-- Polynomial carrying the integral modulus-row coefficients in the
actual alternating degree order. -/
noncomputable def rowCarryPolynomial (count width : Nat)
    (coefficients : Nat → Int) : Polynomial Int :=
  ∑ degree ∈ Finset.range width,
    Polynomial.monomial degree (coefficients (count + cldSpiralIndex width degree))

noncomputable def rowResidualPolynomial (original : Array SparsePolyZZ)
    (size width : Nat) (coefficients : Nat → Int) (modulus : Int) : Polynomial Int :=
  rowCldCombination original size coefficients +
    Polynomial.C modulus * rowCarryPolynomial original.size width coefficients

theorem rowCldCombination_coeff (original : Array SparsePolyZZ) (size : Nat)
    (coefficients : Nat → Int)
    (hcanonical : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (degree : Nat) :
    (rowCldCombination original size coefficients).coeff degree =
      ∑ i : Fin size, if i.val < original.size then coefficients i.val *
        Generated.StrictRecombine.sparseCoeff original[i.val]! degree else 0 := by
  unfold rowCldCombination
  rw [Polynomial.finset_sum_coeff]
  apply Finset.sum_congr rfl
  intro i _
  by_cases hi : i.val < original.size
  · simp only [hi, ↓reduceIte, Polynomial.coeff_C_mul]
    rw [sparseCoeff_eq_coeff _ (hcanonical i.val hi)]
  · simp [hi]

theorem rowCarryPolynomial_coeff (count width : Nat) (coefficients : Nat → Int)
    (degree : Nat) :
    (rowCarryPolynomial count width coefficients).coeff degree =
      if degree < width then coefficients (count + cldSpiralIndex width degree) else 0 := by
  simp [rowCarryPolynomial, Polynomial.finset_sum_coeff, Polynomial.coeff_monomial]

theorem rowResidualPolynomial_coeff (original : Array SparsePolyZZ) (size width : Nat)
    (coefficients : Nat → Int) (modulus : Int)
    (hcanonical : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (degree : Nat) :
    (rowResidualPolynomial original size width coefficients modulus).coeff degree =
      (∑ i : Fin size, if i.val < original.size then coefficients i.val *
        Generated.StrictRecombine.sparseCoeff original[i.val]! degree else 0) +
      modulus * (if degree < width then coefficients (original.size + cldSpiralIndex width degree) else 0) := by
  rw [rowResidualPolynomial, Polynomial.coeff_add, Polynomial.coeff_C_mul,
    rowCldCombination_coeff original size coefficients hcanonical, rowCarryPolynomial_coeff]

/-- One polynomial, with a physical carry at every covered degree, has
exactly the original CLD combination modulo the original modulus. -/
theorem rowResidualPolynomial_modular (original : Array SparsePolyZZ) (size width : Nat)
    (coefficients : Nat → Int) (modulus : Nat) :
    Polynomial.map (Int.castRingHom (ZMod modulus))
      (rowResidualPolynomial original size width coefficients (modulus : Int)) =
    Polynomial.map (Int.castRingHom (ZMod modulus))
      (rowCldCombination original size coefficients) := by
  simp [rowResidualPolynomial, Polynomial.map_add, Polynomial.map_mul]

/-- Full coefficient height, including degrees whose floor coefficients
vanished before the width scan. No bound on the original degree is needed
for this height statement. -/
theorem CanonicalCldBasis.rowResidualPolynomial_bound
    (original truncated : Array SparsePolyZZ) (width : Nat)
    (scale modulus divisor K : Int) (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : width = Generated.StrictRecombine.cldSpiralWidth truncated)
    (hcanonical : CanonicalCldBasis truncated width scale (modulus / divisor) width matrix)
    (horiginal : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (hscale : 1 ≤ scale) (hm : 1 ≤ modulus / divisor) (hK : 0 ≤ K)
    (coefficients : Nat → Int)
    (hphysical : ∀ column, |Matrix.vecMul (fun i : Fin matrix.size => coefficients i.val)
      (basisPrefixMatrix matrix matrix.size matrix.size) column| ≤ K)
    (hvalues : ∀ i, i < original.size → ∀ degree,
      |Generated.StrictRecombine.sparseCoeff original[i]! degree| ≤ modulus)
    (degree : Nat) :
    |(rowResidualPolynomial original matrix.size width coefficients modulus).coeff degree| ≤
      divisor * (2 + 3 * (matrix.size : Int)) * K := by
  have ht := truncateCld_spec original truncated divisor htruncate
  have hsize := hcanonical.size_eq
  rw [ht.2.1] at hsize
  rw [rowResidualPolynomial_coeff original matrix.size width coefficients modulus horiginal]
  by_cases hd : degree < width
  · let column : Fin matrix.size :=
      ⟨original.size + cldSpiralIndex width degree, by
        have hi := cldSpiralIndex_lt width degree hd
        omega⟩
    have hb := CanonicalCldBasis.truncated_physical_row_residual_bound
      original truncated width width scale modulus divisor K matrix htruncate hcanonical
      horiginal hscale hm hK (fun i => coefficients i.val) hphysical column
      (by dsimp [column]; omega) (fun i hi => hvalues i hi _)
    dsimp only [column] at hb
    simp only [Nat.add_sub_cancel_left, cldSpiralDegree_index width degree hd] at hb
    simp only [hd, ↓reduceIte]
    simpa only [mul_comm] using hb.1
  · have hzero : ∀ i, i < original.size →
        Generated.StrictRecombine.sparseCoeff original[i]! degree / divisor = 0 := by
      intro i hi
      have htc : StrictPolynomialMod.SparsePolyZZCanonical truncated[i]! := by
        rw [ht.2.2 i hi]
        exact normalization_canonical _
      have hz := cldSpiralWidth_sparseCoeff_zero truncated i (by rw [ht.2.1]; exact hi)
        htc degree (by rw [← hwidth]; omega)
      rw [truncateCld_sparseCoeff original truncated divisor htruncate i hi (horiginal i hi)] at hz
      exact hz
    let weights : Fin matrix.size → Int := fun i =>
      if i.val < original.size then coefficients i.val else 0
    let values : Fin matrix.size → Int := fun i =>
      if i.val < original.size then Generated.StrictRecombine.sparseCoeff original[i.val]! degree else 0
    have he : ∀ i ∈ (Finset.univ : Finset (Fin matrix.size)), |weights i| ≤ K := by
      intro i _
      dsimp only [weights]
      split
      next hi =>
        have hc := hcanonical.factor_coordinate (fun j => coefficients j.val) i
          (by rw [ht.2.1]; exact hi)
        have hp := hphysical i
        rw [hc, abs_mul, abs_of_nonneg (by omega : 0 ≤ scale)] at hp
        nlinarith only [hp, hscale, abs_nonneg (coefficients i.val)]
      next hi => simpa using hK
    have hv : ∀ i ∈ (Finset.univ : Finset (Fin matrix.size)), |values i| ≤ divisor := by
      intro i _
      dsimp only [values]
      split
      next hi => exact le_of_lt (floor_zero_abs_lt _ divisor ht.1 (hzero i.val hi))
      next hi => simpa using le_of_lt ht.1
    have hs := weighted_sum_abs_le Finset.univ weights values matrix.size K divisor
      (by simp) hK (le_of_lt ht.1) he hv
    have hsum : (∑ i : Fin matrix.size, weights i * values i) =
        ∑ i : Fin matrix.size, if i.val < original.size then coefficients i.val *
          Generated.StrictRecombine.sparseCoeff original[i.val]! degree else 0 := by
      apply Finset.sum_congr rfl
      intro i _
      by_cases hi : i.val < original.size <;> simp [weights, values, hi]
    rw [hsum] at hs
    simp only [hd, ↓reduceIte, mul_zero, add_zero]
    calc
      _ ≤ (matrix.size : Int) * K * divisor := hs
      _ ≤ divisor * (2 + 3 * (matrix.size : Int)) * K := by
        have hn : (0 : Int) ≤ matrix.size := by positivity
        nlinarith only [mul_nonneg (le_of_lt ht.1) hK,
          mul_nonneg hn (mul_nonneg (le_of_lt ht.1) hK)]

/-- An input bound on the original CLD degrees also bounds the complete
residual, including its carry polynomial. The zero polynomial is covered. -/
theorem rowResidualPolynomial_natDegree_lt
    (original : Array SparsePolyZZ) (size width limit : Nat)
    (coefficients : Nat → Int) (modulus : Int) (hlimit : 0 < limit)
    (hwidth : width ≤ limit)
    (hcanonical : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (hdegree : ∀ i, i < original.size →
      (SparsePolyZZ.toPoly original[i]!).natDegree < limit) :
    (rowResidualPolynomial original size width coefficients modulus).natDegree < limit := by
  have hb : (rowResidualPolynomial original size width coefficients modulus).natDegree ≤ limit - 1 := by
    apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
    intro d hd
    rw [rowResidualPolynomial_coeff original size width coefficients modulus hcanonical]
    have hn : ¬d < width := by omega
    simp only [hn, ↓reduceIte, mul_zero, add_zero]
    apply Finset.sum_eq_zero
    intro i _
    split
    next hi =>
      rw [sparseCoeff_eq_coeff _ (hcanonical i.val hi),
        Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_lt_of_le (hdegree i.val hi) (by omega))]
      simp
    next hi => rfl
  omega

/-- A complete retained row from the actual LLL run supplies a single
bounded residual polynomial. Original CLD canonicality and coefficient
bounds come from the existing modular CLD certificate. -/
theorem lllReduceWithNorms_residualPolynomial
    (f : SparsePolyZZ) (active original truncated : Array SparsePolyZZ) (width : Nat)
    (scale divisor bound : ZZ) (modulus : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hcld : ∀ i, i < original.size → ModularCldValue f active[i]! original[i]! modulus)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : width = Generated.StrictRecombine.cldSpiralWidth truncated)
    (hcanonical : CanonicalCldBasis truncated width scale ((modulus : Int) / divisor) width matrix)
    (hscale : 1 ≤ scale) (hm : 1 ≤ (modulus : Int) / divisor) (hB : 0 ≤ bound)
    (hinput : ConcreteLLLInputValid matrix)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (row : Fin matrix.size)
    (hkeep : row.val < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    let coefficients := fun i => (output.2.1[row.val]!)[i]!
    let K := (output.1.size : Int) * 2 ^ output.1.size * bound + 1
    let H := rowResidualPolynomial original matrix.size width coefficients (modulus : Int)
    Polynomial.map (Int.castRingHom (ZMod modulus)) H =
      Polynomial.map (Int.castRingHom (ZMod modulus)) (rowCldCombination original matrix.size coefficients) ∧
      ∀ degree, |H.coeff degree| ≤ divisor * (2 + 3 * (matrix.size : Int)) * K := by
  dsimp only
  refine ⟨rowResidualPolynomial_modular original matrix.size width _ modulus, ?_⟩
  let K : Int := (output.1.size : Int) * 2 ^ output.1.size * bound + 1
  have hK : 0 ≤ K := by dsimp [K]; positivity
  obtain ⟨hs, hrepr⟩ := lllReduceWithNorms_row_representation matrix hinput bound output hrun
  have hphysical : ∀ c : Fin matrix.size,
      |Matrix.vecMul (fun i : Fin matrix.size => (output.2.1[row.val]!)[i.val]!)
        (basisPrefixMatrix matrix matrix.size matrix.size) c| ≤ K := by
    intro c
    rw [hrepr row c]
    exact lllReduceWithNorms_retained_coordinate_abs_le matrix hinput bound hB output hrun
      row.val c.val (by rw [hs]; exact row.isLt) (by rw [hs]; exact c.isLt) hkeep
  intro degree
  exact CanonicalCldBasis.rowResidualPolynomial_bound original truncated width
    scale (modulus : Int) divisor K matrix htruncate hwidth hcanonical
    (fun i hi => (hcld i hi).canonical) hscale hm hK
    (fun i => (output.2.1[row.val]!)[i]!) hphysical
    (fun i hi d => (hcld i hi).sparseCoeff_abs_le d) degree

end Refinement.StrictRecombine
