import CLPoly.Refinement.VanHoeijCertificate
import Mathlib.LinearAlgebra.Matrix.Block

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem arraySetBang_getElemBang_ne {α : Type*} [Inhabited α]
    (input : Array α) (index position : Nat) (value : α)
    (hposition : position < input.size) (hne : index ≠ position) :
    (input.set! index value)[position]! = input[position]! := by
  simp only [Array.set!_eq_setIfInBounds]
  rw [getElem!_pos _ position (by simpa using hposition),
    Array.getElem_setIfInBounds_ne hposition hne,
    ← getElem!_pos input position hposition]

/-- Entry formula for the repaired canonical row basis [sI | C; 0 | mI]. -/
def canonicalCldEntry (cld : Array SparsePolyZZ) (width : Nat)
    (scale modulus : ZZ) (row column : Nat) : ZZ :=
  if column < cld.size then (if row = column then scale else 0)
  else if row < cld.size then
    Generated.StrictRecombine.sparseCoeff cld[row]!
      (cldSpiralDegree width (column - cld.size))
  else if row = column then modulus else 0

structure CanonicalCldBasis (cld : Array SparsePolyZZ) (width : Nat)
    (scale modulus : ZZ) (columns : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) : Prop where
  size_eq : matrix.size = cld.size + columns
  row_size : ∀ row, row < matrix.size → matrix[row]!.size = matrix.size
  entry : ∀ row column, row < matrix.size → column < matrix.size →
    (matrix[row]!)[column]! = canonicalCldEntry cld width scale modulus row column

theorem makeInitialMatrix_canonical
    (cld : Array SparsePolyZZ) (width : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hscale : scale ≠ 0)
    (hrun : Generated.StrictRecombine.makeInitialMatrix cld.size scale = .ok matrix) :
    CanonicalCldBasis cld width scale modulus 0 matrix := by
  obtain ⟨hvalid, hsize⟩ := makeInitialMatrix_input_valid cld.size scale matrix hscale hrun
  have hprefix := (makeInitialMatrix_basisPrefix cld.size scale matrix hrun).2
  refine ⟨by simpa using hsize, ?_, ?_⟩
  · intro row hrow
    simpa only [getElem!_pos matrix row hrow] using hvalid.rows_square row hrow
  · intro row column hrow hcolumn
    have hr : row < cld.size := by omega
    have hc : column < cld.size := by omega
    have he := congrFun (congrFun hprefix ⟨row, hr⟩) ⟨column, hc⟩
    simpa [basisPrefixMatrix, Matrix.diagonal_apply, canonicalCldEntry, hc,
      Fin.ext_iff] using he

theorem appendCldColumn_canonical
    (cld : Array SparsePolyZZ) (width columns : Nat) (scale modulus : ZZ)
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (hcanonical : CanonicalCldBasis cld width scale modulus columns matrix)
    (hrun : Generated.StrictVanHoeij.appendCldColumn matrix cld columns
      (cldSpiralDegree width columns) modulus = .ok output) :
    CanonicalCldBasis cld width scale modulus (columns + 1) output := by
  obtain ⟨hsize, hrows, hlast⟩ := appendCldColumn_rows matrix output cld columns
    (cldSpiralDegree width columns) modulus hrun
  have hn := hcanonical.size_eq
  refine ⟨by omega, ?_, ?_⟩
  · intro row hr
    by_cases hold : row < matrix.size
    · rw [hrows row hold]
      split <;> simp [hcanonical.row_size row hold, hsize]
    · have heq : row = matrix.size := by omega
      subst row
      rw [hlast]
      simp [zeroMatrixRow_size, hsize, hn]
  · intro row column hr hc
    by_cases hrold : row < matrix.size
    · rw [hrows row hrold]
      have hrs := hcanonical.row_size row hrold
      by_cases hcold : column < matrix.size
      · have hne : cld.size + columns ≠ column := by omega
        have hentry := hcanonical.entry row column hrold hcold
        split
        · rw [arraySetBang_getElemBang_ne _ _ column _ (by simp [hrs]; omega) hne,
            arrayPush_getElem!_lt _ _ column (by omega)]
          exact hentry
        · rw [arrayPush_getElem!_lt _ _ column (by omega)]
          exact hentry
      · have heq : column = matrix.size := by omega
        subst column
        by_cases hfactor : row < cld.size
        · simp [hfactor, getElem!_pos, hrs, hn, canonicalCldEntry]
        · have hne : row ≠ cld.size + columns := by omega
          have hpush : (matrix[row]!.push 0)[matrix.size]! = 0 := by
            simpa only [hrs] using arrayPush_getElem!_last matrix[row]! (0 : ZZ)
          rw [if_neg hfactor, hpush]
          simp [canonicalCldEntry, hn, hfactor, hne]
    · have heq : row = matrix.size := by omega
      subst row
      rw [hlast]
      have hz := zeroMatrixRow_get (cld.size + columns + 1) column (by omega)
      by_cases heq : column = cld.size + columns
      · subst column
        simp [getElem!_pos, zeroMatrixRow_size, canonicalCldEntry, hn]
      · have hne : cld.size + columns ≠ column := Ne.symm heq
        have hrnot : ¬cld.size + columns < cld.size := by omega
        rw [arraySetBang_getElemBang_ne _ _ column _
          (by rw [zeroMatrixRow_size]; omega) hne, hz]
        simp [canonicalCldEntry, hn, hrnot, hne]

theorem buildCldMatrixLoop_canonical
    (cld : Array SparsePolyZZ) (width current target added : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (output : Generated.StrictRecombine.LLLMatrix × Nat)
    (hcanonical : CanonicalCldBasis cld width scale modulus (current + added) matrix)
    (hrun : Generated.StrictVanHoeij.buildCldMatrixLoop cld current target width modulus
      added matrix = .ok output) :
    CanonicalCldBasis cld width scale modulus (current + output.2) output.1 := by
  induction hmeasure : target - added using Nat.strong_induction_on
      generalizing added matrix with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.buildCldMatrixLoop] at hrun
    split at hrun
    next htarget =>
      dsimp only at hrun
      split at hrun
      next hwidth =>
        split at hrun
        next fault hstep => contradiction
        next nextMatrix hstep =>
          have hnext := appendCldColumn_canonical cld width (current + added)
            scale modulus matrix nextMatrix hcanonical (by
              simpa only [cldSpiralDegree] using hstep)
          apply ih (target - (added + 1)) (by omega) (added + 1) nextMatrix
            (by simpa only [Nat.add_assoc] using hnext) hrun rfl
      next hwidth =>
        cases Except.ok.inj hrun
        exact hcanonical
    next htarget =>
      cases Except.ok.inj hrun
      exact hcanonical

theorem buildCldMatrix_canonical
    (cld : Array SparsePolyZZ) (current target : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (output : Generated.StrictRecombine.LLLMatrix × Nat)
    (hcanonical : CanonicalCldBasis cld (Generated.StrictRecombine.cldSpiralWidth cld)
      scale modulus current matrix)
    (hrun : Generated.StrictVanHoeij.buildCldMatrix matrix cld current target modulus =
      .ok output) :
    CanonicalCldBasis cld (Generated.StrictRecombine.cldSpiralWidth cld)
      scale modulus (current + output.2) output.1 := by
  exact buildCldMatrixLoop_canonical cld _ current target 0 scale modulus matrix output
    (by simpa using hcanonical) hrun

/-- The concrete canonical matrix is a valid full-rank integer LLL input.
Neither a nonsingularity oracle nor the obsolete identity-modulus matrix is used. -/
theorem CanonicalCldBasis.input_valid
    {cld : Array SparsePolyZZ} {width columns : Nat} {scale modulus : ZZ}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hcanonical : CanonicalCldBasis cld width scale modulus columns matrix)
    (hscale : scale ≠ 0) (hmodulus : modulus ≠ 0) :
    ConcreteLLLInputValid matrix := by
  refine ⟨?_, ?_⟩
  · intro row hrow
    simpa only [getElem!_pos matrix row hrow] using hcanonical.row_size row hrow
  · have htri : (basisPrefixMatrix matrix matrix.size matrix.size).BlockTriangular id := by
      intro row column hlt
      change column < row at hlt
      change (matrix[row.val]!)[column.val]! = 0
      rw [hcanonical.entry row.val column.val row.isLt column.isLt]
      have hne : row.val ≠ column.val := by omega
      by_cases hc : column.val < cld.size
      · simp [canonicalCldEntry, hc, hne]
      · have hr : ¬row.val < cld.size := by omega
        simp [canonicalCldEntry, hc, hr, hne]
    rw [Matrix.det_of_upperTriangular htri]
    apply Finset.prod_ne_zero_iff.mpr
    intro index _
    change (matrix[index.val]!)[index.val]! ≠ 0
    rw [hcanonical.entry index.val index.val index.isLt index.isLt]
    by_cases hi : index.val < cld.size <;>
      simp [canonicalCldEntry, hi, hscale, hmodulus]

/-- Composition of the actual initialization and column-building executions. -/
theorem initialized_buildCldMatrix_input_valid
    (cld : Array SparsePolyZZ) (target : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (output : Generated.StrictRecombine.LLLMatrix × Nat)
    (hscale : scale ≠ 0) (hmodulus : modulus ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix cld.size scale = .ok matrix)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix matrix cld 0 target modulus =
      .ok output) :
    ConcreteLLLInputValid output.1 := by
  exact (buildCldMatrix_canonical cld 0 target scale modulus matrix output
    (makeInitialMatrix_canonical cld _ scale modulus matrix hscale hinit) hbuild).input_valid
      hscale hmodulus

theorem CanonicalCldBasis.factor_coordinate
    {cld : Array SparsePolyZZ} {width columns : Nat} {scale modulus : ZZ}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hcanonical : CanonicalCldBasis cld width scale modulus columns matrix)
    (coefficients : Fin matrix.size → ZZ) (column : Fin matrix.size)
    (hcolumn : column.val < cld.size) :
    Matrix.vecMul coefficients (basisPrefixMatrix matrix matrix.size matrix.size) column =
      coefficients column * scale := by
  simp only [Matrix.vecMul, dotProduct]
  calc
    (∑ row, coefficients row * basisPrefixMatrix matrix matrix.size matrix.size row column) =
        ∑ row, if row = column then coefficients row * scale else 0 := by
      apply Finset.sum_congr rfl
      intro row _
      change coefficients row * (matrix[row.val]!)[column.val]! = _
      rw [hcanonical.entry row.val column.val row.isLt column.isLt]
      simp [canonicalCldEntry, hcolumn, Fin.ext_iff]
    _ = coefficients column * scale := by simp

theorem split_dotProduct_bound
    (size firstCount lastCount : Nat) (hsize : size = firstCount + lastCount)
    (value : Fin size → QQ) (firstBound lastBound : QQ)
    (hfirst : ∀ index, index.val < firstCount → value index ^ 2 ≤ firstBound ^ 2)
    (hlast : ∀ index, firstCount ≤ index.val → value index ^ 2 ≤ lastBound ^ 2) :
    dotProduct value value ≤
      (firstCount : QQ) * firstBound ^ 2 + (lastCount : QQ) * lastBound ^ 2 := by
  subst size
  simp only [dotProduct, ← pow_two]
  rw [Fin.sum_univ_add]
  calc
    _ ≤ (∑ _ : Fin firstCount, firstBound ^ 2) +
        ∑ _ : Fin lastCount, lastBound ^ 2 := by
      apply add_le_add
      · apply Finset.sum_le_sum
        intro index _
        exact hfirst (Fin.castAdd lastCount index) index.isLt
      · apply Finset.sum_le_sum
        intro index _
        exact hlast (Fin.natAdd firstCount index) (by simp)
    _ = _ := by simp

/-- Data coordinates are the selected CLD sum plus the integral modulus
carry from the corresponding added basis row. -/
theorem CanonicalCldBasis.data_coordinate
    {cld : Array SparsePolyZZ} {width columns : Nat} {scale modulus : ZZ}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hcanonical : CanonicalCldBasis cld width scale modulus columns matrix)
    (coefficients : Fin matrix.size → ZZ) (column : Fin matrix.size)
    (hcolumn : cld.size ≤ column.val) :
    Matrix.vecMul coefficients (basisPrefixMatrix matrix matrix.size matrix.size) column =
      (∑ row : Fin matrix.size, if row.val < cld.size then
        coefficients row * Generated.StrictRecombine.sparseCoeff cld[row.val]!
          (cldSpiralDegree width (column.val - cld.size)) else 0) +
        coefficients column * modulus := by
  simp only [Matrix.vecMul, dotProduct]
  calc
    (∑ row, coefficients row * basisPrefixMatrix matrix matrix.size matrix.size row column) =
        ∑ row, ((if row.val < cld.size then coefficients row *
          Generated.StrictRecombine.sparseCoeff cld[row.val]!
            (cldSpiralDegree width (column.val - cld.size)) else 0) +
          if row = column then coefficients row * modulus else 0) := by
      apply Finset.sum_congr rfl
      intro row _
      change coefficients row * (matrix[row.val]!)[column.val]! = _
      rw [hcanonical.entry row.val column.val row.isLt column.isLt]
      have hc : ¬column.val < cld.size := by omega
      by_cases hr : row.val < cld.size
      · have hne : row ≠ column := by intro heq; subst row; omega
        simp [canonicalCldEntry, hc, hr, hne]
      · simp [canonicalCldEntry, hc, hr, Fin.ext_iff]
    _ = _ := by rw [Finset.sum_add_distrib]; simp

/-- The bound passed to LLL follows from actual factor-coordinate entries
and the bounded residuals of the actual data coordinates. -/
theorem CanonicalCldBasis.vector_norm_bound
    {cld : Array SparsePolyZZ} {width columns : Nat} {scale modulus : ZZ}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hcanonical : CanonicalCldBasis cld width scale modulus columns matrix)
    (coefficients : Fin matrix.size → ZZ) (residual : ZZ)
    (hindicator : ∀ index, index.val < cld.size →
      coefficients index = 0 ∨ coefficients index = 1)
    (hresidual : ∀ index : Fin matrix.size, cld.size ≤ index.val →
      |Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index| ≤ residual) :
    dotProduct
      (fun index => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
      (fun index => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ)) ≤
      ((cld.size : ZZ) * scale * scale + (columns : ZZ) * residual * residual : ZZ) := by
  have hbound := split_dotProduct_bound matrix.size cld.size columns hcanonical.size_eq
    (fun index => ((Matrix.vecMul coefficients
      (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
    (scale : QQ) (residual : QQ) ?_ ?_
  · simpa [pow_two, mul_assoc] using hbound
  · intro index hi
    dsimp only
    rw [hcanonical.factor_coordinate coefficients index hi]
    rcases hindicator index hi with hz | ho
    · simp [hz, sq_nonneg]
    · simp [ho]
  · intro index hi
    dsimp only
    have hr := abs_le.mp (hresidual index hi)
    have hsquare : (Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index) ^ 2 ≤ residual ^ 2 := by
      nlinarith [sq_nonneg (Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index + residual)]
    exact_mod_cast hsquare

/-- The two concrete preconditions needed by the downstream LLL/grouping
proof, derived from the actual matrix-building execution. -/
theorem initialized_buildCldMatrix_short_vector
    (cld : Array SparsePolyZZ) (target : Nat) (scale modulus residual : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (output : Generated.StrictRecombine.LLLMatrix × Nat)
    (hscale : scale ≠ 0) (hmodulus : modulus ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix cld.size scale = .ok matrix)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix matrix cld 0 target modulus =
      .ok output)
    (coefficients : Fin output.1.size → ZZ)
    (hindicator : ∀ index, index.val < cld.size →
      coefficients index = 0 ∨ coefficients index = 1)
    (hresidual : ∀ index : Fin output.1.size, cld.size ≤ index.val →
      |Matrix.vecMul coefficients
        (basisPrefixMatrix output.1 output.1.size output.1.size) index| ≤ residual) :
    ConcreteLLLInputValid output.1 ∧
      dotProduct
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix output.1 output.1.size output.1.size) index : ZZ) : QQ))
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix output.1 output.1.size output.1.size) index : ZZ) : QQ)) ≤
        ((cld.size : ZZ) * scale * scale + (output.2 : ZZ) * residual * residual : ZZ) := by
  have hcanonical := buildCldMatrix_canonical cld 0 target scale modulus matrix output
    (makeInitialMatrix_canonical cld _ scale modulus matrix hscale hinit) hbuild
  refine ⟨hcanonical.input_valid hscale hmodulus, ?_⟩
  simpa only [Nat.zero_add] using
    hcanonical.vector_norm_bound coefficients residual hindicator hresidual

/-- The modular carry is bounded using only |C_i| ≤ m and H < m.
This also covers either convention at a symmetric-remainder endpoint. -/
theorem modular_carry_bound (radius : Nat) (modulus height actual total carry : Int)
    (hm : 0 < modulus) (hh : height < modulus)
    (ha : |actual| ≤ height) (ht : |total| ≤ (radius : Int) * modulus)
    (heq : total = actual + carry * modulus) :
    -(radius : Int) ≤ carry ∧ carry ≤ (radius : Int) := by
  obtain ⟨haLow, haHigh⟩ := abs_le.mp ha
  obtain ⟨htLow, htHigh⟩ := abs_le.mp ht
  constructor
  · by_contra h
    have hc : carry ≤ -(radius : Int) - 1 := by omega
    nlinarith
  · by_contra h
    have hc : (radius : Int) + 1 ≤ carry := by omega
    nlinarith

theorem ceil_div_mul_ge (height divisor : Int) (hd : 0 < divisor) :
    height ≤ divisor * ((height + divisor - 1) / divisor) := by
  have heq := Int.mul_ediv_add_emod (height + divisor - 1) divisor
  have hrem := Int.emod_lt_of_pos (height + divisor - 1) hd
  omega

/-- The error bound for the actual floor-divided coefficients and modulus,
derived from a modular CLD identity rather than an assumed residual. -/
theorem truncated_subset_residual_bound {ι : Type*}
    (support : Finset ι) (values : ι → Int) (radius : Nat)
    (modulus divisor height actual carry : Int)
    (hcard : support.card ≤ radius) (hm : 0 < modulus) (hd : 0 < divisor)
    (hh : height < modulus) (ha : |actual| ≤ height)
    (hvalues : ∀ index ∈ support, |values index| ≤ modulus)
    (heq : ∑ index ∈ support, values index = actual + carry * modulus) :
    |(∑ index ∈ support, values index / divisor) - carry * (modulus / divisor)| ≤
      (height + divisor - 1) / divisor + 2 * (radius : Int) := by
  have hcardInt : (support.card : Int) ≤ radius := by exact_mod_cast hcard
  have htotal : |∑ index ∈ support, values index| ≤ (radius : Int) * modulus := by
    calc
      _ ≤ ∑ index ∈ support, |values index| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ _index ∈ support, modulus := Finset.sum_le_sum hvalues
      _ = (support.card : Int) * modulus := by simp
      _ ≤ (radius : Int) * modulus := mul_le_mul_of_nonneg_right hcardInt (le_of_lt hm)
  have hc := modular_carry_bound radius modulus height actual _ carry hm hh ha htotal heq
  have hremNonneg : 0 ≤ ∑ index ∈ support, values index % divisor := by
    apply Finset.sum_nonneg
    intro index _
    exact Int.emod_nonneg _ (ne_of_gt hd)
  have hremBound : (∑ index ∈ support, values index % divisor) ≤ (radius : Int) * divisor := by
    calc
      _ ≤ ∑ _index ∈ support, divisor := by
        apply Finset.sum_le_sum
        intro index _
        exact le_of_lt (Int.emod_lt_of_pos _ hd)
      _ = (support.card : Int) * divisor := by simp
      _ ≤ (radius : Int) * divisor := mul_le_mul_of_nonneg_right hcardInt (le_of_lt hd)
  have hsum : divisor * (∑ index ∈ support, values index / divisor) +
      (∑ index ∈ support, values index % divisor) = ∑ index ∈ support, values index := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro index _
    exact Int.mul_ediv_add_emod (values index) divisor
  have hmod := Int.mul_ediv_add_emod modulus divisor
  have hidentity : divisor * ((∑ index ∈ support, values index / divisor) -
      carry * (modulus / divisor)) = actual - (∑ index ∈ support, values index % divisor) +
        carry * (modulus % divisor) := by
    rw [heq] at hsum
    nlinarith
  apply abs_le.mpr
  exact truncated_cld_coordinate_bound divisor radius height
    ((height + divisor - 1) / divisor) actual carry
    (∑ index ∈ support, values index % divisor) (modulus % divisor)
    ((∑ index ∈ support, values index / divisor) - carry * (modulus / divisor))
    hd (by exact_mod_cast Nat.zero_le radius) (ceil_div_mul_ge height divisor hd)
    (abs_le.mp ha) hc ⟨hremNonneg, hremBound⟩
    ⟨Int.emod_nonneg _ (ne_of_gt hd), Int.emod_lt_of_pos _ hd⟩ hidentity

/-- Explicit short lattice vector from a modular support and its CLD carries.
Its first coordinates are indicators; the added coordinates are negative
carries. The norm bound is derived, not a caller premise. -/
theorem CanonicalCldBasis.truncated_indicator_short
    {cld : Array SparsePolyZZ} {width columns : Nat} {scale : ZZ}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (modulus divisor height : Int)
    (hcanonical : CanonicalCldBasis cld width scale (modulus / divisor) columns matrix)
    (support : Finset (Fin matrix.size))
    (raw : Fin matrix.size → Fin matrix.size → Int)
    (actual carry : Fin matrix.size → Int)
    (hcard : support.card ≤ cld.size)
    (hsupport : ∀ index ∈ support, index.val < cld.size)
    (hm : 0 < modulus) (hd : 0 < divisor) (hh : height < modulus)
    (hactual : ∀ column, cld.size ≤ column.val → |actual column| ≤ height)
    (hraw : ∀ column, cld.size ≤ column.val →
      ∀ row ∈ support, |raw row column| ≤ modulus)
    (hcongr : ∀ column, cld.size ≤ column.val →
      (∑ row ∈ support, raw row column) = actual column + carry column * modulus)
    (hfloor : ∀ column, cld.size ≤ column.val → ∀ row ∈ support,
      Generated.StrictRecombine.sparseCoeff cld[row.val]!
        (cldSpiralDegree width (column.val - cld.size)) = raw row column / divisor) :
    let coefficients := fun index : Fin matrix.size =>
      if index.val < cld.size then (if index ∈ support then (1 : ZZ) else 0)
      else -carry index
    let residual := (height + divisor - 1) / divisor + 2 * (cld.size : Int)
    dotProduct
      (fun index => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
      (fun index => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ)) ≤
      ((cld.size : ZZ) * scale * scale + (columns : ZZ) * residual * residual : ZZ) := by
  dsimp only
  let coefficients := fun index : Fin matrix.size =>
    if index.val < cld.size then (if index ∈ support then (1 : ZZ) else 0)
    else -carry index
  apply hcanonical.vector_norm_bound coefficients
    ((height + divisor - 1) / divisor + 2 * (cld.size : Int))
  · intro index hi
    by_cases hmem : index ∈ support <;> simp [coefficients, hi, hmem]
  · intro column hc
    rw [hcanonical.data_coordinate coefficients column hc]
    have hsum : (∑ row : Fin matrix.size, if row.val < cld.size then
        coefficients row * Generated.StrictRecombine.sparseCoeff cld[row.val]!
          (cldSpiralDegree width (column.val - cld.size)) else 0) =
        ∑ row ∈ support, raw row column / divisor := by
      calc
        _ = ∑ row : Fin matrix.size, if row ∈ support then raw row column / divisor else 0 := by
          apply Finset.sum_congr rfl
          intro row _
          by_cases hmem : row ∈ support
          · have hr := hsupport row hmem
            simpa only [coefficients, if_pos hr, if_pos hmem, one_mul] using
              hfloor column hc row hmem
          · by_cases hr : row.val < cld.size <;> simp [coefficients, hr, hmem]
        _ = _ := by simp
    rw [hsum]
    have hcoeff : coefficients column = -carry column := by
      simp [coefficients, not_lt_of_ge hc]
    rw [hcoeff]
    simpa only [neg_mul, sub_eq_add_neg] using
      truncated_subset_residual_bound support (fun row => raw row column) cld.size
        modulus divisor height (actual column) (carry column) hcard hm hd hh
        (hactual column hc) (hraw column hc) (hcongr column hc)

end Refinement.StrictRecombine
