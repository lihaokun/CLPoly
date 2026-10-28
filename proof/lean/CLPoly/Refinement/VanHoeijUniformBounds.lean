import CLPoly.Refinement.VanHoeijSeparationInput
import CLPoly.Refinement.VanHoeijParameters

set_option autoImplicit false

namespace Refinement.StrictRecombine

def separationCutoff (count degree : Nat) : Int :=
  (count : Int) * 32 * 32 + (degree : Int) * (32 + 2 * (count : Int)) ^ 2

def separationCoordinateBound (count degree : Nat) : Int :=
  ((count + degree : Nat) : Int) * 2 ^ (count + degree) * separationCutoff count degree + 1

def separationResidualHeight (count degree : Nat) (divisor : Int) : Int :=
  divisor * (2 + 3 * ((count + degree : Nat) : Int)) * separationCoordinateBound count degree

theorem separationCutoff_nonneg (count degree : Nat) : 0 ≤ separationCutoff count degree := by
  unfold separationCutoff
  positivity

theorem separationCoordinateBound_nonneg (count degree : Nat) :
    0 ≤ separationCoordinateBound count degree := by
  unfold separationCoordinateBound
  have := separationCutoff_nonneg count degree
  positivity

theorem retained_coordinate_bound_mono (size limit : Nat) (bound upper : Int)
    (hsize : size ≤ limit) (hb : 0 ≤ bound) (hu : bound ≤ upper) :
    (size : Int) * 2 ^ size * bound + 1 ≤ (limit : Int) * 2 ^ limit * upper + 1 := by
  have hs : (size : Int) ≤ limit := by exact_mod_cast hsize
  have hp : (2 : Int) ^ size ≤ 2 ^ limit := pow_le_pow_right₀ (by norm_num) hsize
  have hm := mul_le_mul hs hp (by positivity) (by positivity : (0 : Int) ≤ limit)
  have ht := mul_le_mul hm hu hb (by positivity : (0 : Int) ≤ (limit : Int) * 2 ^ limit)
  omega

/-- The actual rebuilt dimension and number of added columns are bounded
by the factor count and source degree from the execution certificate. -/
theorem CldSeparationData.rebuild_size_bounds
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (divisor scale : ZZ) (target : Nat)
    (initial matrix : Generated.StrictRecombine.LLLMatrix) (columns : Nat)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hscale : scale ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix truncated.size scale = .ok initial)
    (htarget : target ≤ Generated.StrictRecombine.cldSpiralWidth truncated)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix initial truncated 0 target
      ((modulus : Int) / divisor) = .ok (matrix, columns)) :
    columns ≤ (SparsePolyZZ.toPoly f).natDegree ∧
      matrix.size ≤ active.size + (SparsePolyZZ.toPoly f).natDegree := by
  have ht := truncateCld_spec original truncated divisor htruncate
  have hw := truncateCld_width_le original truncated divisor htruncate
    (fun i hi => (data.values i hi).canonical) _ data.degrees
  have hc := buildCldMatrix_added_le truncated 0 target ((modulus : Int) / divisor) initial
    (matrix, columns) hbuild
  have hb := buildCldMatrix_canonical truncated 0 target scale ((modulus : Int) / divisor)
    initial (matrix, columns)
    (makeInitialMatrix_canonical truncated _ scale ((modulus : Int) / divisor) initial hscale hinit) hbuild
  have hs := hb.size_eq
  dsimp only at hc hs
  rw [ht.2.1, data.size_eq] at hs
  constructor <;> omega

/-- Instantiation of the actual parameter-selection cutoff with the
physical factor count. Its bound contains no modulus. -/
theorem selected_cutoff_le_separationCutoff
    (count : Int32) (radius degree columns : Nat) (height modulus : ZZ)
    (parameters : Generated.StrictVanHoeijParameters.Parameters)
    (hcount : count.toInt = (radius : Int)) (hheight : 0 ≤ height)
    (hrun : Generated.StrictVanHoeijParameters.selectParameters count height modulus = some parameters)
    (hcolumns : columns ≤ degree) :
    0 ≤ (radius : Int) * 32 * 32 + (columns : Int) * parameters.residual * parameters.residual ∧
      (radius : Int) * 32 * 32 + (columns : Int) * parameters.residual * parameters.residual ≤
        separationCutoff radius degree := by
  have hn : 0 ≤ count.toInt := by rw [hcount]; positivity
  have hb := selectParameters_squared_cutoff_le count height modulus parameters hn hheight
    hrun columns degree hcolumns
  rw [hcount] at hb
  refine ⟨?_, hb⟩
  have hs := sq_nonneg parameters.residual
  have hmul := mul_nonneg (show (0 : Int) ≤ columns by positivity) hs
  nlinarith only [hmul, show (0 : Int) ≤ (radius : Int) * 32 * 32 by positivity]

theorem lllReduceWithNorms_uniform_coordinate_bound
    (count degree : Nat) (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (hB : 0 ≤ bound)
    (hsize : matrix.size ≤ count + degree) (hbound : bound ≤ separationCutoff count degree)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output)
    (row column : Nat) (hrow : row < output.1.size) (hcolumn : column < output.1.size)
    (hkeep : row < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    |(output.1[row]!)[column]!| ≤ separationCoordinateBound count degree := by
  have hs := (lllReduceWithNorms_row_representation matrix hinput bound output hrun).1
  have hb := lllReduceWithNorms_retained_coordinate_abs_le matrix hinput bound hB output hrun
    row column hrow hcolumn hkeep
  exact hb.trans (retained_coordinate_bound_mono output.1.size (count + degree)
    bound (separationCutoff count degree) (by rwa [hs]) hB hbound)

/-- Uniform factor-coordinate and full residual height bounds for the same
actual retained transform row. D remains explicit in the residual height. -/
theorem CldSeparationData.retained_uniform_bounds
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
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    let coefficients := fun i => (output.2.1[row.val]!)[i]!
    let K := separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree
    (∀ i : Fin matrix.size, i.val < original.size → |coefficients i.val| ≤ K) ∧
      ∀ d, |(rowResidualPolynomial original matrix.size
        (Generated.StrictRecombine.cldSpiralWidth truncated) coefficients (modulus : Int)).coeff d| ≤
        separationResidualHeight active.size (SparsePolyZZ.toPoly f).natDegree divisor := by
  let K := separationCoordinateBound active.size (SparsePolyZZ.toPoly f).natDegree
  have hK : 0 ≤ K := separationCoordinateBound_nonneg _ _
  have ht := truncateCld_spec original truncated divisor htruncate
  obtain ⟨hs, hrepr⟩ := lllReduceWithNorms_row_representation matrix hinput bound output hrun
  have hphysical : ∀ c : Fin matrix.size,
      |Matrix.vecMul (fun i : Fin matrix.size => (output.2.1[row.val]!)[i.val]!)
        (basisPrefixMatrix matrix matrix.size matrix.size) c| ≤ K := by
    intro c
    rw [hrepr row c]
    exact lllReduceWithNorms_uniform_coordinate_bound active.size (SparsePolyZZ.toPoly f).natDegree
      matrix hinput bound hB hsize hbound output hrun row.val c.val
      (by rw [hs]; exact row.isLt) (by rw [hs]; exact c.isLt) hkeep
  constructor
  · intro i hi
    have hc := hcanonical.factor_coordinate (fun j => (output.2.1[row.val]!)[j.val]!) i
      (by rw [ht.2.1]; exact hi)
    have hp := hphysical i
    rw [hc, abs_mul, abs_of_nonneg (le_trans (by decide : (0 : ZZ) ≤ 1) hscale)] at hp
    dsimp only [ZZ] at hscale hp ⊢
    nlinarith only [hp, hscale, abs_nonneg ((output.2.1[row.val]!)[i.val]!)]
  · intro d
    have hb := CanonicalCldBasis.rowResidualPolynomial_bound original truncated _ scale
      (modulus : Int) divisor K matrix htruncate rfl hcanonical
      (fun i hi => (data.values i hi).canonical) hscale hm hK
      (fun i => (output.2.1[row.val]!)[i]!) hphysical
      (fun i hi degree => (data.values i hi).sparseCoeff_abs_le degree) d
    apply hb.trans
    unfold separationResidualHeight
    have hn : (matrix.size : Int) ≤ ((active.size + (SparsePolyZZ.toPoly f).natDegree : Nat) : Int) :=
      by exact_mod_cast hsize
    apply mul_le_mul_of_nonneg_right _ hK
    apply mul_le_mul_of_nonneg_left _ (le_of_lt ht.1)
    exact add_le_add le_rfl (mul_le_mul_of_nonneg_left hn (by norm_num))

end Refinement.StrictRecombine
