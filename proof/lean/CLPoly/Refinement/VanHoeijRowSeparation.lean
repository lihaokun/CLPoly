import CLPoly.Refinement.VanHoeijSeparationHeight

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem residual_difference_natDegree_le (f H : Polynomial Int) (scalar : Int)
    (hH : H.natDegree < f.natDegree) :
    (H - Polynomial.C scalar * Polynomial.derivative f).natDegree ≤ f.natDegree := by
  apply (Polynomial.natDegree_sub_le H (Polynomial.C scalar * Polynomial.derivative f)).trans
  apply max_le (le_of_lt hH)
  exact (Polynomial.natDegree_C_mul_le scalar (Polynomial.derivative f)).trans
    ((Polynomial.natDegree_derivative_le f).trans (Nat.sub_le _ _))

/-- The resultant threshold forces equality of bounded row coefficients
at any two local factors belonging to the same genuine irreducible divisor.
All local CLD and coprimality conditions come from the execution data. -/
theorem CldSeparationData.scalars_equal_of_sufficient_precision
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (size : Nat) (coefficients : Nat → Int)
    (hcoefficients : ∀ i : Fin size, i.val < original.size →
      |coefficients i.val| ≤ separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree)
    (hresidual : ∀ k, |(rowResidualPolynomial original size
      (Generated.StrictRecombine.cldSpiralWidth truncated) coefficients (modulus : Int)).coeff k| ≤
        separationResidualHeight active.size (SparsePolyZZ.toPoly f).natDegree divisor)
    (factor : Polynomial Int) (hf : SparsePolyZZ.toPoly f ≠ 0)
    (hfactor : factor ∣ SparsePolyZZ.toPoly f) (hprimitive : factor.IsPrimitive)
    (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : factor.natDegree ≠ 0)
    (left right : Fin size) (hl : left.val < original.size) (hr : right.val < original.size)
    (hleft : StrictHensel.toPolyMod modulus active[left.val]! ∣
      factor.map (Int.castRingHom (ZMod modulus)))
    (hright : StrictHensel.toPolyMod modulus active[right.val]! ∣
      factor.map (Int.castRingHom (ZMod modulus)))
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) active.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    coefficients left.val = coefficients right.val := by
  haveI : NeZero modulus := ⟨by have := data.modulus_gt_one; omega⟩
  let H := rowResidualPolynomial original size (Generated.StrictRecombine.cldSpiralWidth truncated)
    coefficients (modulus : Int)
  let G := H - Polynomial.C (coefficients left.val) * Polynomial.derivative (SparsePolyZZ.toPoly f)
  let B := separationPolynomialHeight (SparsePolyZZ.toPoly f) active.size divisor
  have hd := (truncateCld_spec original truncated divisor htruncate).1
  have hB : 1 ≤ B := separationPolynomialHeight_ge_one _ _ _ (le_of_lt hd)
  have hp := data.residual_properties truncated divisor htruncate size coefficients
  have hGdegree : G.natDegree ≤ (SparsePolyZZ.toPoly f).natDegree :=
    residual_difference_natDegree_le _ _ _ hp.2.1
  have hFdegree : factor.natDegree ≤ (SparsePolyZZ.toPoly f).natDegree :=
    Polynomial.natDegree_le_of_dvd hfactor hf
  have hdimension : factor.natDegree + G.natDegree ≤ 2 * (SparsePolyZZ.toPoly f).natDegree := by omega
  have hGB : ∀ k, |G.coeff k| ≤ B :=
    residual_difference_coeff_bound _ H active.size divisor (coefficients left.val)
      (le_of_lt hd) (hcoefficients left hl) hresidual
  have hFB : ∀ k, |factor.coeff k| ≤ B :=
    source_factor_coeff_le_separationPolynomialHeight _ factor active.size divisor (le_of_lt hd) hf hfactor
  have hlocal : StrictHensel.toPolyMod modulus active[left.val]! ∣
      G.map (Int.castRingHom (ZMod modulus)) := by
    simpa only [G, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_C,
      ← Polynomial.derivative_map, StrictHensel.toPolyMod] using (hp.2.2 left hl).1
  have hlarge : ((factor.natDegree + G.natDegree).factorial : Int) *
      B ^ (factor.natDegree + G.natDegree) < (modulus : Int) :=
    (resultant_coefficient_bound_mono _ _ B hdimension hB).trans_lt hresultant
  have hdiv : factor ∣ G := CLPoly.VanHoeij.primitive_irreducible_dvd_of_sufficient_modulus
    modulus data.modulus_gt_one factor G (StrictHensel.toPolyMod modulus active[left.val]!)
    hprimitive hirreducible hfactorDegree (data.monic left.val hl)
    (data.local_degree_pos left.val hl) hleft hlocal B (by omega) hFB hGB hlarge
  have hrightG := hright.trans (Polynomial.map_dvd (Int.castRingHom (ZMod modulus)) hdiv)
  apply CLPoly.VanHoeij.bounded_local_cld_scalars_equal modulus data.modulus_gt_one
    (H.map (Int.castRingHom (ZMod modulus))) (Polynomial.derivative (StrictHensel.toPolyMod modulus f))
    (StrictHensel.toPolyMod modulus active[right.val]!) (coefficients left.val) (coefficients right.val)
    (separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree)
    (hcoefficients left hl) (hcoefficients right hr) hcoordinates
    (data.monic right.val hr) (data.local_degree_pos right.val hr)
    (data.derivative_coprime right.val hr)
  · simpa only [G, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_C,
      ← Polynomial.derivative_map, StrictHensel.toPolyMod] using hrightG
  · exact (hp.2.2 right hr).1

theorem local_factor_dvd_of_support (active : Array SparsePolyZZ) (modulus : Nat)
    (factor : Polynomial Int) (indices : List Nat) (scalar : ZMod modulus)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) = Polynomial.C scalar *
      (indices.map fun i => StrictHensel.toPolyMod modulus active[i]!).prod)
    (index : Nat) (hi : index ∈ indices) :
    StrictHensel.toPolyMod modulus active[index]! ∣ factor.map (Int.castRingHom (ZMod modulus)) := by
  rw [hsupport]
  exact dvd_mul_of_dvd_right (List.dvd_prod (List.mem_map.mpr ⟨index, hi, rfl⟩)) _

/-- Actual retained transform rows are constant on any genuine irreducible
factor's Hensel support once the explicit input thresholds are exceeded. -/
theorem CldSeparationData.retained_row_equal_on_support
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (scale divisor bound : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : CanonicalCldBasis truncated (Generated.StrictRecombine.cldSpiralWidth truncated)
      scale ((modulus : Int) / divisor) (Generated.StrictRecombine.cldSpiralWidth truncated) matrix)
    (hscale : 1 ≤ scale) (hm : 1 ≤ (modulus : Int) / divisor)
    (hinput : ConcreteLLLInputValid matrix) (hB : 0 ≤ bound)
    (hsize : matrix.size ≤ active.size + (SparsePolyZZ.toPoly f).natDegree)
    (hbound : bound ≤ separationCutoff active.size (SparsePolyZZ.toPoly f).natDegree)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output)
    (row : Fin matrix.size)
    (hkeep : row.val < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size)
    (factor : Polynomial Int) (hfactor : factor ∣ SparsePolyZZ.toPoly f)
    (hprimitive : factor.IsPrimitive) (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : factor.natDegree ≠ 0)
    (indices : List Nat) (scalar : ZMod modulus)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) = Polynomial.C scalar *
      (indices.map fun i => StrictHensel.toPolyMod modulus active[i]!).prod)
    (left right : Fin matrix.size) (hl : left.val < original.size) (hr : right.val < original.size)
    (hlSupport : left.val ∈ indices) (hrSupport : right.val ∈ indices)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) active.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    (output.2.1[row.val]!)[left.val]! = (output.2.1[row.val]!)[right.val]! := by
  have hb := data.retained_uniform_bounds truncated scale divisor bound matrix htruncate hcanonical
    hscale hm hinput hB hsize hbound output hrun row hkeep
  exact data.scalars_equal_of_sufficient_precision truncated divisor htruncate matrix.size
    (fun i => (output.2.1[row.val]!)[i]!) hb.1 hb.2 factor
    (Polynomial.ne_zero_of_natDegree_gt data.source_degree_pos) hfactor hprimitive hirreducible
    hfactorDegree left right hl hr
    (local_factor_dvd_of_support active modulus factor indices scalar hsupport left.val hlSupport)
    (local_factor_dvd_of_support active modulus factor indices scalar hsupport right.val hrSupport)
    hresultant hcoordinates

end Refinement.StrictRecombine
