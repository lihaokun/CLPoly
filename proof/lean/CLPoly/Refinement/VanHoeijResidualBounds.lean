import CLPoly.Refinement.VanHoeijTruncation
import CLPoly.Refinement.VanHoeijRowBounds

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Flooring a negative centered coefficient costs at most one additional
unit. This uses Euclidean division, as does the physical truncation loop. -/
theorem floor_abs_le_modulus_quotient_add_one (value modulus divisor : Int)
    (hd : 0 < divisor) (hv : |value| ≤ modulus) :
    |value / divisor| ≤ modulus / divisor + 1 := by
  have hvdiv := Int.mul_ediv_add_emod value divisor
  have hmdiv := Int.mul_ediv_add_emod modulus divisor
  have hvr0 := Int.emod_nonneg value (ne_of_gt hd)
  have hvrD := Int.emod_lt_of_pos value hd
  have hmr0 := Int.emod_nonneg modulus (ne_of_gt hd)
  have hmrD := Int.emod_lt_of_pos modulus hd
  obtain ⟨hvl, hvu⟩ := abs_le.mp hv
  apply abs_le.mpr
  constructor
  · by_contra h
    have hq : value / divisor ≤ -(modulus / divisor + 1) - 1 := by omega
    nlinarith
  · by_contra h
    have hq : modulus / divisor + 2 ≤ value / divisor := by omega
    nlinarith

theorem weighted_sum_abs_le {ι : Type*} (support : Finset ι)
    (weights values : ι → Int) (radius : Nat) (E V : Int)
    (hcard : support.card ≤ radius) (hE : 0 ≤ E) (hV : 0 ≤ V)
    (he : ∀ i ∈ support, |weights i| ≤ E)
    (hv : ∀ i ∈ support, |values i| ≤ V) :
    |∑ i ∈ support, weights i * values i| ≤ (radius : Int) * E * V := by
  have hc : (support.card : Int) ≤ radius := by exact_mod_cast hcard
  calc
    _ ≤ ∑ i ∈ support, |weights i * values i| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i ∈ support, E * V := by
      apply Finset.sum_le_sum
      intro i hi
      rw [abs_mul]
      exact mul_le_mul (he i hi) (hv i hi) (abs_nonneg _) hE
    _ = (support.card : Int) * E * V := by simp [mul_assoc]
    _ ≤ _ := by
      simpa only [mul_assoc] using mul_le_mul_of_nonneg_right hc (mul_nonneg hE hV)

/-- The integral coefficient of the modulus row is controlled by the
physical data coordinate and factor coefficients, independently of m'. -/
theorem truncated_row_carry_abs_le {ι : Type*} (support : Finset ι)
    (weights values : ι → Int) (radius : Nat)
    (modulus divisor E Y carry data : Int)
    (hcard : support.card ≤ radius) (hd : 0 < divisor)
    (hm : 1 ≤ modulus / divisor) (hE : 0 ≤ E) (hY : 0 ≤ Y)
    (he : ∀ i ∈ support, |weights i| ≤ E)
    (hv : ∀ i ∈ support, |values i| ≤ modulus)
    (hy : |data| ≤ Y)
    (hdata : data = (∑ i ∈ support, weights i * (values i / divisor)) +
      carry * (modulus / divisor)) :
    |carry| ≤ Y + 2 * (radius : Int) * E := by
  have hs := weighted_sum_abs_le support weights (fun i => values i / divisor)
    radius E (modulus / divisor + 1) hcard hE (by omega) he
    (fun i hi => floor_abs_le_modulus_quotient_add_one _ _ _ hd (hv i hi))
  have heq : carry * (modulus / divisor) =
      data - ∑ i ∈ support, weights i * (values i / divisor) := by omega
  have hc : |carry| * (modulus / divisor) ≤
      Y + (radius : Int) * E * (modulus / divisor + 1) := by
    calc
      _ = |carry * (modulus / divisor)| := by rw [abs_mul, abs_of_nonneg (by omega : 0 ≤ modulus / divisor)]
      _ = _ := congrArg abs heq
      _ ≤ |data| + |∑ i ∈ support, weights i * (values i / divisor)| := abs_sub _ _
      _ ≤ _ := add_le_add hy hs
  have hrE : 0 ≤ (radius : Int) * E := mul_nonneg (by positivity) hE
  have hc' : |carry| * (modulus / divisor) ≤
      (Y + 2 * (radius : Int) * E) * (modulus / divisor) := by
    nlinarith [mul_nonneg hY (show 0 ≤ modulus / divisor - 1 by omega),
      mul_nonneg hrE (show 0 ≤ modulus / divisor - 1 by omega)]
  exact le_of_mul_le_mul_right hc' (show 0 < modulus / divisor by omega)

/-- Undoing the common floor division produces this exact residual, with
the same integral modulus-row coefficient as in the truncated lattice. -/
theorem truncated_row_residual_identity {ι : Type*} (support : Finset ι)
    (weights values : ι → Int) (modulus divisor carry : Int) :
    (∑ i ∈ support, weights i * values i) + carry * modulus =
      divisor * ((∑ i ∈ support, weights i * (values i / divisor)) +
        carry * (modulus / divisor)) +
      (∑ i ∈ support, weights i * (values i % divisor)) + carry * (modulus % divisor) := by
  have hs : (∑ i ∈ support, weights i * values i) =
      divisor * (∑ i ∈ support, weights i * (values i / divisor)) +
        ∑ i ∈ support, weights i * (values i % divisor) := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _
    have hi := congrArg (weights i * ·) (Int.mul_ediv_add_emod (values i) divisor)
    nlinarith only [hi]
  rw [hs]
  have hm := congrArg (carry * ·) (Int.mul_ediv_add_emod modulus divisor)
  nlinarith only [hm]

/-- A single explicitly constructed residual has both the required
congruence and a height bound for arbitrary signed row coefficients. -/
theorem truncated_row_residual_bound {ι : Type*} (support : Finset ι)
    (weights values : ι → Int) (radius : Nat)
    (modulus divisor E Y carry data : Int)
    (hcard : support.card ≤ radius) (hd : 0 < divisor)
    (hm : 1 ≤ modulus / divisor) (hE : 0 ≤ E) (hY : 0 ≤ Y)
    (he : ∀ i ∈ support, |weights i| ≤ E)
    (hv : ∀ i ∈ support, |values i| ≤ modulus)
    (hy : |data| ≤ Y)
    (hdata : data = (∑ i ∈ support, weights i * (values i / divisor)) +
      carry * (modulus / divisor)) :
    |(∑ i ∈ support, weights i * values i) + carry * modulus| ≤
        divisor * (2 * Y + 3 * (radius : Int) * E) ∧
      modulus ∣ ((∑ i ∈ support, weights i * values i) + carry * modulus) -
        (∑ i ∈ support, weights i * values i) := by
  have hc := truncated_row_carry_abs_le support weights values radius modulus divisor
    E Y carry data hcard hd hm hE hY he hv hy hdata
  have hr (v : Int) : |v % divisor| ≤ divisor := by
    rw [abs_of_nonneg (Int.emod_nonneg v (ne_of_gt hd))]
    exact le_of_lt (Int.emod_lt_of_pos v hd)
  have hs := weighted_sum_abs_le support weights (fun i => values i % divisor)
    radius E divisor hcard hE (le_of_lt hd) he (fun i _ => hr (values i))
  constructor
  · rw [truncated_row_residual_identity, ← hdata]
    calc
      |divisor * data + (∑ i ∈ support, weights i * (values i % divisor)) +
          carry * (modulus % divisor)| ≤
          |divisor * data + ∑ i ∈ support, weights i * (values i % divisor)| +
          |carry * (modulus % divisor)| := abs_add_le _ _
      _ ≤ (|divisor * data| + |∑ i ∈ support, weights i * (values i % divisor)|) +
          |carry * (modulus % divisor)| := add_le_add
            (abs_add_le (divisor * data) (∑ i ∈ support, weights i * (values i % divisor)))
            le_rfl
      _ ≤ (divisor * Y + (radius : Int) * E * divisor) +
          (Y + 2 * (radius : Int) * E) * divisor := by
        apply add_le_add
        · apply add_le_add _ hs
          rw [abs_mul, abs_of_nonneg (le_of_lt hd)]
          exact mul_le_mul_of_nonneg_left hy (le_of_lt hd)
        · rw [abs_mul]
          exact mul_le_mul hc (hr modulus) (abs_nonneg _) (by positivity)
      _ = _ := by ring
  · refine ⟨carry, ?_⟩
    ring

/-- The same physical transform returned by LLL expresses each reduced
row as an integral combination of the original canonical rows. -/
theorem lllReduceWithNorms_row_representation
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    output.1.size = matrix.size ∧
      ∀ row column : Fin matrix.size,
        Matrix.vecMul (fun i : Fin matrix.size => (output.2.1[row.val]!)[i.val]!)
          (basisPrefixMatrix matrix matrix.size matrix.size) column =
          (output.1[row.val]!)[column.val]! := by
  obtain ⟨state, _, hm, ht, _, hrel, _⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  refine ⟨by rw [← hm]; exact hrel.1, ?_⟩
  intro row column
  have he := congrFun (congrFun hrel.2.2.2 row) column
  rw [hm, ht] at he
  simpa only [basisPrefixMatrix, Matrix.mul_apply, Matrix.vecMul, dotProduct] using he.symm

/-- Residual height for one actually added CLD column. The coefficient
bounds are derived from the physical vector, not postulated for the transform. -/
theorem CanonicalCldBasis.truncated_physical_row_residual_bound
    (original truncated : Array SparsePolyZZ) (width columns : Nat)
    (scale modulus divisor K : ZZ) (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : CanonicalCldBasis truncated width scale (modulus / divisor) columns matrix)
    (horiginal : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (hscale : 1 ≤ scale) (hm : 1 ≤ modulus / divisor) (hK : 0 ≤ K)
    (coefficients : Fin matrix.size → ZZ)
    (hphysical : ∀ column, |Matrix.vecMul coefficients
      (basisPrefixMatrix matrix matrix.size matrix.size) column| ≤ K)
    (column : Fin matrix.size) (hcolumn : original.size ≤ column.val)
    (hvalues : ∀ i, i < original.size →
      |Generated.StrictRecombine.sparseCoeff original[i]!
        (cldSpiralDegree width (column.val - original.size))| ≤ modulus) :
    let total := ∑ i : Fin matrix.size, if i.val < original.size then
      coefficients i * Generated.StrictRecombine.sparseCoeff original[i.val]!
        (cldSpiralDegree width (column.val - original.size)) else 0
    |total + coefficients column * modulus| ≤ divisor * (2 + 3 * (matrix.size : ZZ)) * K ∧
      modulus ∣ (total + coefficients column * modulus) - total := by
  have ht := truncateCld_spec original truncated divisor htruncate
  let weights : Fin matrix.size → Int := fun i => if i.val < original.size then coefficients i else 0
  let values : Fin matrix.size → Int := fun i => if i.val < original.size then
    Generated.StrictRecombine.sparseCoeff original[i.val]!
      (cldSpiralDegree width (column.val - original.size)) else 0
  have he : ∀ i ∈ (Finset.univ : Finset (Fin matrix.size)), |weights i| ≤ K := by
    intro i _
    dsimp only [weights]
    split
    next hi =>
      have hc := hcanonical.factor_coordinate coefficients i (by rw [ht.2.1]; exact hi)
      have hp := hphysical i
      rw [hc, abs_mul, abs_of_nonneg (le_trans (by decide : (0 : ZZ) ≤ 1) hscale)] at hp
      dsimp only [ZZ] at hscale hp ⊢
      nlinarith only [hp, hscale, abs_nonneg (coefficients i)]
    next hi => simpa using hK
  have hmod : 0 ≤ modulus := by
    have heq := Int.mul_ediv_add_emod modulus divisor
    have hr := Int.emod_nonneg modulus (ne_of_gt ht.1)
    nlinarith only [heq, hr, mul_nonneg (le_of_lt ht.1) (le_trans (by decide : (0 : Int) ≤ 1) hm)]
  have hv : ∀ i ∈ (Finset.univ : Finset (Fin matrix.size)), |values i| ≤ modulus := by
    intro i _
    dsimp only [values]
    split
    next hi => exact hvalues i.val hi
    next hi => simpa using hmod
  have hdata := hcanonical.data_coordinate coefficients column (by rw [ht.2.1]; exact hcolumn)
  rw [ht.2.1] at hdata
  have hsum : (∑ i : Fin matrix.size, if i.val < original.size then
      coefficients i * Generated.StrictRecombine.sparseCoeff truncated[i.val]!
        (cldSpiralDegree width (column.val - original.size)) else 0) =
      ∑ i : Fin matrix.size, weights i * (values i / divisor) := by
    apply Finset.sum_congr rfl
    intro i _
    by_cases hi : i.val < original.size
    · simp only [hi, ↓reduceIte, weights, values]
      rw [truncateCld_sparseCoeff original truncated divisor htruncate i.val hi (horiginal i.val hi)]
    · simp [hi, weights, values]
  rw [hsum] at hdata
  have hb := truncated_row_residual_bound Finset.univ weights values matrix.size
    modulus divisor K K (coefficients column)
    (Matrix.vecMul coefficients (basisPrefixMatrix matrix matrix.size matrix.size) column)
    (by simp) ht.1 hm hK hK he hv (hphysical column) hdata
  have htotal : (∑ i : Fin matrix.size, weights i * values i) =
      ∑ i : Fin matrix.size, if i.val < original.size then
        coefficients i * Generated.StrictRecombine.sparseCoeff original[i.val]!
          (cldSpiralDegree width (column.val - original.size)) else 0 := by
    apply Finset.sum_congr rfl
    intro i _
    by_cases hi : i.val < original.size <;> simp [weights, values, hi]
  rw [htotal] at hb
  constructor
  · calc
      _ ≤ divisor * (2 * K + 3 * (matrix.size : ZZ) * K) := hb.1
      _ = _ := by ring
  · exact hb.2

/-- The residual bound instantiated with the actual LLL output and actual
retention test. K depends only on dimension and the short-vector cutoff;
the remaining dependence on the truncation divisor is kept explicit. -/
theorem lllReduceWithNorms_truncated_residual_bound
    (original truncated : Array SparsePolyZZ) (width columns : Nat)
    (scale modulus divisor bound : ZZ) (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : CanonicalCldBasis truncated width scale (modulus / divisor) columns matrix)
    (horiginal : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (hscale : 1 ≤ scale) (hm : 1 ≤ modulus / divisor) (hB : 0 ≤ bound)
    (hinput : ConcreteLLLInputValid matrix)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (row column : Fin matrix.size)
    (hkeep : row.val < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size)
    (hcolumn : original.size ≤ column.val)
    (hvalues : ∀ i, i < original.size →
      |Generated.StrictRecombine.sparseCoeff original[i]!
        (cldSpiralDegree width (column.val - original.size))| ≤ modulus) :
    let K := (output.1.size : ZZ) * 2 ^ output.1.size * bound + 1
    let total := ∑ i : Fin matrix.size, if i.val < original.size then
      (output.2.1[row.val]!)[i.val]! * Generated.StrictRecombine.sparseCoeff original[i.val]!
        (cldSpiralDegree width (column.val - original.size)) else 0
    |total + (output.2.1[row.val]!)[column.val]! * modulus| ≤
        divisor * (2 + 3 * (matrix.size : ZZ)) * K ∧
      modulus ∣ (total + (output.2.1[row.val]!)[column.val]! * modulus) - total := by
  let K : ZZ := (output.1.size : ZZ) * 2 ^ output.1.size * bound + 1
  have hK : 0 ≤ K := by dsimp [K]; positivity
  obtain ⟨hs, hrepr⟩ := lllReduceWithNorms_row_representation matrix hinput bound output hrun
  have hphysical : ∀ c : Fin matrix.size,
      |Matrix.vecMul (fun i : Fin matrix.size => (output.2.1[row.val]!)[i.val]!)
        (basisPrefixMatrix matrix matrix.size matrix.size) c| ≤ K := by
    intro c
    rw [hrepr row c]
    exact lllReduceWithNorms_retained_coordinate_abs_le matrix hinput bound hB output hrun
      row.val c.val (by rw [hs]; exact row.isLt) (by rw [hs]; exact c.isLt) hkeep
  exact CanonicalCldBasis.truncated_physical_row_residual_bound original truncated
    width columns scale modulus divisor K matrix htruncate hcanonical horiginal hscale hm hK
    (fun i => (output.2.1[row.val]!)[i.val]!) hphysical column hcolumn hvalues

end Refinement.StrictRecombine
