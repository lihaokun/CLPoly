import CLPoly.Refinement.Recombine
import CLPoly.Algorithm.VanHoeijBound
import CLPoly.Generated.StrictVanHoeij

set_option autoImplicit false

namespace Refinement.StrictRecombine

open Generated.StrictVanHoeij (retainedPrefixCount)

/-- Product-rule identity for a selected group of monic modular factors.
It holds over ZMod(p^a), including its zero divisors. -/
theorem cld_sum_monic_product {R : Type*} [CommRing R]
    (factors : List (Polynomial R)) (complement : Polynomial R)
    (hmonic : ∀ factor ∈ factors, factor.Monic) :
    (factors.map fun factor =>
      ((complement * factors.prod) /ₘ factor) * Polynomial.derivative factor).sum =
        complement * Polynomial.derivative factors.prod := by
  induction factors generalizing complement with
  | nil => simp
  | cons head tail ih =>
    have hh := hmonic head (by simp)
    have ht : ∀ factor ∈ tail, factor.Monic := fun factor hmem =>
      hmonic factor (List.mem_cons_of_mem head hmem)
    simp only [List.prod_cons, List.map_cons, List.sum_cons]
    have hquotient : (complement * (head * tail.prod)) /ₘ head =
        complement * tail.prod := by
      rw [show complement * (head * tail.prod) = head * (complement * tail.prod) by ring]
      exact Polynomial.mul_divByMonic_cancel_left _ hh
    have htail := ih (complement * head) ht
    simp only [mul_assoc] at htail
    rw [hquotient, htail, Polynomial.derivative_mul]
    ring

/-- Hensel factors need only divide modulo m. No integer divisibility of a
lifted factor is assumed, and the coefficient ring may have zero divisors. -/
theorem generatedModularQuotient_eq_modular_divByMonic
    (termination : Generated.StrictHensel.DivmodTermination)
    (f divisor quotient remainder : SparsePolyZZ) (modulus : Nat)
    (hdivisorSize : 0 < divisor.size)
    (hdivisorDegree : divisor[0]!.1.deg < 2 ^ 63)
    (hdivisorHead : StrictHensel.HeadDominates divisor)
    (hfBound : StrictHensel.DegreesBound f)
    (hinvert : (ZZ.invert 0 divisor[0]!.2 (modulus : Int)).1 = true)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f divisor (modulus : Int) = .ok (quotient, remainder))
    (hremainder : remainder.isEmpty = true)
    (hmonic : (StrictHensel.toPolyMod modulus divisor).Monic) :
    StrictHensel.toPolyMod modulus quotient =
      StrictHensel.toPolyMod modulus f /ₘ StrictHensel.toPolyMod modulus divisor := by
  have hdivision := StrictHensel.__upoly_divmod_mod_raw_ir_refines_of_run
    termination f divisor quotient remainder modulus hdivisorSize
    hdivisorDegree hdivisorHead hfBound hinvert hrun
  have hremainderEmpty : remainder = #[] := Array.isEmpty_iff.mp hremainder
  subst remainder
  have hproduct : StrictHensel.toPolyMod modulus quotient *
      StrictHensel.toPolyMod modulus divisor = StrictHensel.toPolyMod modulus f := by
    simpa using hdivision
  rw [← hproduct, mul_comm, Polynomial.mul_divByMonic_cancel_left _ hmonic]

/-- The physically computed CLD denotes the modular logarithmic derivative,
without requiring the lifted factor to be an integer divisor of f. -/
theorem generatedCldElement_modular
    (termination : Generated.StrictHensel.DivmodTermination)
    (f divisor quotient remainder derivativeRaw derivativeMod product cld : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (hdivisorSize : 0 < divisor.size)
    (hdivisorDegree : divisor[0]!.1.deg < 2 ^ 63)
    (hdivisorHead : StrictHensel.HeadDominates divisor)
    (hfBound : StrictHensel.DegreesBound f)
    (hinvert : (ZZ.invert 0 divisor[0]!.2 (modulus : Int)).1 = true)
    (hdivision : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f divisor (modulus : Int) = .ok (quotient, remainder))
    (hremainder : remainder.isEmpty = true)
    (hderivative : Generated.StrictRecombine.derivativeZZRaw divisor = .ok derivativeRaw)
    (hmod : Generated.StrictRecombine.modCoeffLoop derivativeRaw
      (modulus : Int) 0 #[] = .ok derivativeMod)
    (hproduct : Generated.StrictRecombine.multiplyNormalizeModRaw quotient
      derivativeMod (modulus : Int) = .ok product)
    (hsymmetric : Generated.StrictRecombine.symmetricModRaw product
      (modulus : Int) = .ok cld)
    (hmonic : (StrictHensel.toPolyMod modulus divisor).Monic) :
    StrictHensel.toPolyMod modulus cld =
      (StrictHensel.toPolyMod modulus f /ₘ StrictHensel.toPolyMod modulus divisor) *
        Polynomial.derivative (StrictHensel.toPolyMod modulus divisor) := by
  have hquotient := generatedModularQuotient_eq_modular_divByMonic termination f divisor
    quotient remainder modulus hdivisorSize hdivisorDegree hdivisorHead hfBound hinvert
    hdivision hremainder hmonic
  have hderivativeSemantic := derivativeZZRaw_toPoly divisor derivativeRaw hderivative
  have hderivativeMod := modCoeffLoop_toPolyMod derivativeRaw modulus hmodulus
    0 #[] derivativeMod hmod
  have hd : StrictHensel.toPolyMod modulus derivativeMod =
      Polynomial.derivative (StrictHensel.toPolyMod modulus divisor) := by
    calc
      StrictHensel.toPolyMod modulus derivativeMod =
          StrictHensel.toPolyMod modulus derivativeRaw := by
        simpa [StrictHensel.toPolyMod_eq_termsToPolyMod] using hderivativeMod
      _ = Polynomial.derivative (StrictHensel.toPolyMod modulus divisor) := by
        simp only [StrictHensel.toPolyMod, hderivativeSemantic, Polynomial.derivative_map]
  rw [symmetricModRaw_toPolyMod product cld modulus hmodulus hsymmetric,
    multiplyNormalizeModRaw_correct quotient derivativeMod modulus product hmodulus hproduct,
    hquotient, hd]

theorem fillCldColumnLoop_size
    (cld : Array SparsePolyZZ) (degree column index : Nat)
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (hrun : Generated.StrictVanHoeij.fillCldColumnLoop cld degree column index matrix =
      .ok output) : output.size = matrix.size := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on
      generalizing index matrix with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.fillCldColumnLoop] at hrun
    split at hrun
    next hi =>
      split at hrun
      next hm =>
        split at hrun
        next hc =>
          have htail := ih (cld.size - (index + 1)) (by omega) (index + 1)
            (matrix.set index (matrix[index].set column
              (Generated.StrictRecombine.sparseCoeff cld[index] degree))) hrun rfl
          simpa using htail
        next hc => contradiction
      next hm => contradiction
    next hi => exact congrArg Array.size (Except.ok.inj hrun).symm

theorem fillCldColumnLoop_row
    (cld : Array SparsePolyZZ) (degree column index : Nat)
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (hrun : Generated.StrictVanHoeij.fillCldColumnLoop cld degree column index matrix =
      .ok output) (row : Nat) (hrow : row < matrix.size) :
    output[row]! = if index ≤ row ∧ row < cld.size then
      (matrix[row]!).set! column (Generated.StrictRecombine.sparseCoeff cld[row]! degree)
      else matrix[row]! := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on
      generalizing index matrix with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.fillCldColumnLoop] at hrun
    split at hrun
    next hi =>
      split at hrun
      next hm =>
        split at hrun
        next hc =>
          have htail := ih (cld.size - (index + 1)) (by omega) (index + 1)
            (matrix.set index (matrix[index].set column
              (Generated.StrictRecombine.sparseCoeff cld[index] degree))) hrun
            (by simpa using hrow) rfl
          by_cases heq : row = index
          · subst row
            simpa [hi, getElem!_pos, hm, hc, Array.setIfInBounds_def] using htail
          · have hget :
                (matrix.set index (matrix[index].set column
                  (Generated.StrictRecombine.sparseCoeff cld[index] degree)))[row]! =
                  matrix[row]! := by
              simp [getElem!_pos, hrow, Ne.symm heq]
            rw [hget] at htail
            have hiff : index + 1 ≤ row ↔ index ≤ row := by omega
            simpa only [hiff] using htail
        next hc => contradiction
      next hm => contradiction
    next hi =>
      have hout := Except.ok.inj hrun
      subst output
      have hnot : ¬(index ≤ row ∧ row < cld.size) := by omega
      simp [hnot]

/-- The new row contains the actual modulus, not the old identity entry. -/
theorem appendCldColumn_modulus_row
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (existing degree : Nat) (modulus : ZZ)
    (hrun : Generated.StrictVanHoeij.appendCldColumn matrix cld existing degree modulus =
      .ok output) :
    output.size = matrix.size + 1 ∧
    (output[matrix.size]!)[cld.size + existing]! = modulus := by
  unfold Generated.StrictVanHoeij.appendCldColumn at hrun
  obtain ⟨extended, hext, hextSize⟩ :
      ∃ extended, Generated.StrictRecombine.appendZeroColumn matrix = .ok extended ∧
        extended.size = matrix.size :=
    ⟨_, appendZeroColumn_eq matrix, appendZeroSuffix_size matrix⟩
  rw [hext] at hrun
  simp only [bind, Except.bind] at hrun
  cases hfill : Generated.StrictVanHoeij.fillCldColumnLoop cld degree
      (cld.size + existing) 0 extended with
  | error fault => simp [hfill] at hrun
  | ok filled =>
    have hsize := fillCldColumnLoop_size cld degree (cld.size + existing) 0
      extended filled hfill
    rw [hextSize] at hsize
    simp only [hfill, bind, Except.bind, zeroMatrixRow_size, Nat.lt_succ_self,
      ↓reduceDIte] at hrun
    have hout := Except.ok.inj hrun
    subst output
    constructor
    · simp [hsize]
    · rw [← hsize, arrayPush_getElem!_last]
      simp [getElem!_pos, zeroMatrixRow_size]

/-- Full physical row description of a successful CLD-column append. -/
theorem appendCldColumn_rows
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (existing degree : Nat) (modulus : ZZ)
    (hrun : Generated.StrictVanHoeij.appendCldColumn matrix cld existing degree modulus =
      .ok output) :
    output.size = matrix.size + 1 ∧
    (∀ row, row < matrix.size → output[row]! =
      if row < cld.size then (matrix[row]!.push 0).set! (cld.size + existing)
        (Generated.StrictRecombine.sparseCoeff cld[row]! degree)
      else matrix[row]!.push 0) ∧
    output[matrix.size]! =
      (Generated.StrictRecombine.zeroMatrixRow (cld.size + existing + 1)).set!
        (cld.size + existing) modulus := by
  obtain ⟨extended, hext, hextSize, hextRows⟩ :
      ∃ extended, Generated.StrictRecombine.appendZeroColumn matrix = .ok extended ∧
        extended.size = matrix.size ∧
        ∀ row, row < matrix.size → extended[row]! = matrix[row]!.push 0 :=
    ⟨_, appendZeroColumn_eq matrix, appendZeroSuffix_size matrix,
      appendZeroSuffix_row matrix⟩
  unfold Generated.StrictVanHoeij.appendCldColumn at hrun
  rw [hext] at hrun
  simp only [bind, Except.bind] at hrun
  cases hfill : Generated.StrictVanHoeij.fillCldColumnLoop cld degree
      (cld.size + existing) 0 extended with
  | error fault => simp [hfill] at hrun
  | ok filled =>
    have hsize := fillCldColumnLoop_size cld degree (cld.size + existing) 0
      extended filled hfill
    rw [hextSize] at hsize
    simp only [hfill, zeroMatrixRow_size, Nat.lt_succ_self, ↓reduceDIte] at hrun
    have hout := Except.ok.inj hrun
    subst output
    refine ⟨by simp [hsize], ?_, ?_⟩
    · intro row hrow
      rw [arrayPush_getElem!_lt _ _ row (by omega)]
      have hf := fillCldColumnLoop_row cld degree (cld.size + existing) 0
        extended filled hfill row (by omega)
      simpa only [Nat.zero_le, true_and, hextRows row hrow] using hf
    · rw [← hsize, arrayPush_getElem!_last]
      simp [Array.setIfInBounds_def, zeroMatrixRow_size]

theorem retainedPrefixCount_spec (norms : Array QQ) (bound : QQ)
    (remaining : Nat) (hsize : remaining ≤ norms.size) :
    retainedPrefixCount norms bound remaining ≤ remaining ∧
    ∀ index, index < remaining →
      retainedPrefixCount norms bound remaining ≤ index →
      bound < norms[index]! := by
  induction remaining with
  | zero => simp [retainedPrefixCount]
  | succ remaining ih =>
      have hindex : remaining < norms.size := by omega
      rw [retainedPrefixCount]
      simp only [Nat.succ_pos, ↓reduceDIte,
        Nat.add_sub_cancel, hindex]
      split
      next hlarge =>
        obtain ⟨hle, htail⟩ := ih (by omega)
        refine ⟨by omega, ?_⟩
        intro index hlt hretained
        by_cases hlast : index = remaining
        · subst index
          simpa only [getElem!_pos norms remaining hindex] using hlarge
        · exact htail index (by omega) hretained
      next hsmall => exact ⟨le_rfl, by intro index hlt hge; omega⟩

theorem iotaLoop_spec (count index : Nat) (rows : Array Nat)
    (hsize : rows.size = count) :
    (Generated.StrictVanHoeij.iotaLoop count index rows).size = count ∧
    ∀ position, position < count →
      (Generated.StrictVanHoeij.iotaLoop count index rows)[position]! =
        if position < index then rows[position]! else position := by
  rw [Generated.StrictVanHoeij.iotaLoop]
  split
  next hindex =>
    have hnextSize : (rows.set! index index).size = count := by simp [hsize]
    obtain ⟨hs, hv⟩ := iotaLoop_spec count (index + 1)
      (rows.set! index index) hnextSize
    refine ⟨hs, ?_⟩
    intro position hp
    rw [hv position hp]
    by_cases hbefore : position < index
    · simp [hbefore, show position < index + 1 by omega,
        getElem!_pos, hsize, hp,
        show index ≠ position by omega]
    · by_cases heq : position = index
      · subst position
        simp [hindex, hsize]
      · simp [hbefore, show ¬position < index + 1 by omega]
  next hdone =>
    refine ⟨hsize, ?_⟩
    intro position hp
    simp [show position < index by omega]
termination_by count - index

theorem shortSpanFromNorms_spec (norms : Array QQ) (bound : ZZ) :
    let rows := Generated.StrictVanHoeij.shortSpanFromNorms norms bound
    rows.size ≤ norms.size ∧
    (∀ position, position < rows.size → rows[position]! = position) ∧
    (∀ position, position < norms.size → rows.size ≤ position →
      (bound : QQ) < norms[position]!) := by
  let count := retainedPrefixCount norms (bound : QQ) norms.size
  have hi := iotaLoop_spec count 0 (Array.replicate count 0) (by simp)
  have hc := retainedPrefixCount_spec norms (bound : QQ) norms.size le_rfl
  let rows := Generated.StrictVanHoeij.iotaLoop count 0 (Array.replicate count 0)
  change rows.size ≤ norms.size ∧
    (∀ position, position < rows.size → rows[position]! = position) ∧
    (∀ position, position < norms.size → rows.size ≤ position → _)
  have hs : rows.size = count := hi.1
  rw [hs]
  refine ⟨hc.1, ?_, hc.2⟩
  · intro position hp
    simpa [rows] using hi.2 position hp

/-- The truncation error bound used for every appended CLD coordinate.
`u` is the sum of selected coefficient remainders, `v` the modulus
remainder, and `z` the integral modular carry. -/
theorem truncated_cld_coordinate_bound
    (divisor radius height ceiling actual carry remainderSum modulusRemainder
      residual : Int)
    (hdivisor : 0 < divisor) (hradius : 0 ≤ radius)
    (hheight : height ≤ divisor * ceiling)
    (hactual : -height ≤ actual ∧ actual ≤ height)
    (hcarry : -radius ≤ carry ∧ carry ≤ radius)
    (hsum : 0 ≤ remainderSum ∧ remainderSum ≤ radius * divisor)
    (hmod : 0 ≤ modulusRemainder ∧ modulusRemainder < divisor)
    (hidentity : divisor * residual =
      actual - remainderSum + carry * modulusRemainder) :
    -(ceiling + 2 * radius) ≤ residual ∧ residual ≤ ceiling + 2 * radius := by
  have hupper := mul_le_mul_of_nonneg_right hcarry.2 hmod.1
  have hlower := mul_le_mul_of_nonneg_right hcarry.1 hmod.1
  have hrem := mul_le_mul_of_nonneg_left (le_of_lt hmod.2) hradius
  constructor <;> nlinarith

/-- Safe suffix removal, unlike filtering ordinary basis-row lengths, keeps
every short integral lattice vector in the retained prefix.  This is a
statement about the physical GS state already verified for concrete LLL. -/
theorem short_vector_coordinates_vanish_in_removed_suffix
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (rowCount : Nat) (hrowCount : rowCount ≤ state.matrix.size)
    (coordinates : Fin rowCount → ZZ) (retained : Nat) (bound : QQ)
    (hremoved : ∀ index : Fin rowCount, retained ≤ index.val →
      bound < state.norms[index.val]!)
    (hshort : dotProduct (fun position => (coordinates position : QQ))
        (Matrix.mulVec (gramPrefixMatrixQQ state.matrix rowCount)
          (fun position => (coordinates position : QQ))) ≤ bound) :
    ∀ index : Fin rowCount, retained ≤ index.val → coordinates index = 0 := by
  intro index hindex
  by_contra hnonzero
  have hne : coordinates ≠ 0 := by
    intro hz
    apply hnonzero
    simpa [hz]
  obtain ⟨last, hlast, hafter⟩ := exists_last_nonzero coordinates hne
  have hindexLast : index ≤ last := by
    by_contra hnot
    exact hnonzero (hafter index (lt_of_not_ge hnot))
  have hretained : retained ≤ last.val := le_trans hindex hindexLast
  have hbound := gramSchmidt_norm_le_of_last_nonzero state hvalid rowCount
    hrowCount coordinates last hlast hafter
  exact (not_le_of_gt (hremoved last hretained)) (hbound.trans hshort)

/-- All invariants belong to the same physically executed LLL state. -/
theorem lllReduceWithNorms_invariants
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    ∃ state : Generated.StrictRecombine.LLLState,
      ConcreteLLLExecutionValid state ∧ state.matrix = output.1 ∧
      state.transform = output.2.1 ∧ state.norms = output.2.2.2 ∧
      LLLTransformRel matrix state.matrix state.transform ∧
      LLLTransformUnimodular state.transform := by
  rw [Generated.StrictVanHoeij.lllReduceWithNorms] at hrun
  split at hrun
  next fault hinitialize => contradiction
  next mu norms transform hinitialize =>
    dsimp only at hrun
    split at hrun
    next fault hloop => contradiction
    next state hloop =>
      split at hrun
      next fault hcollect => contradiction
      next rows hcollect =>
        have hvalid := initializeLLL_concrete_valid matrix mu norms transform hinput hinitialize
        have hfinal := concreteLLLMainLoop_preserves_execution_valid
          (Generated.StrictRecombine.LLLState.mk matrix transform mu norms 1)
          state hvalid hloop
        cases Except.ok.inj hrun
        refine ⟨state, hfinal, rfl, rfl, rfl, ?_, ?_⟩
        · exact concreteLLLMainLoop_preserves_transform_rel matrix _ state hvalid
            (initializeLLL_transform_rel matrix mu norms transform hinitialize) hloop
        · exact concreteLLLMainLoop_preserves_transform_unimodular matrix _ state hvalid
            (initializeLLL_transform_rel matrix mu norms transform hinitialize)
            (initializeLLL_transform_unimodular matrix mu norms transform hinput hinitialize)
            hloop

/-- Exported norms are the cache of an actually executed, verified GS state. -/
theorem lllReduceWithNorms_valid
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    ∃ state : Generated.StrictRecombine.LLLState,
      ConcreteLLLExecutionValid state ∧ state.matrix = output.1 ∧
      state.transform = output.2.1 ∧ state.norms = output.2.2.2 := by
  obtain ⟨state, hv, hm, ht, hn, _⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  exact ⟨state, hv, hm, ht, hn⟩

/-- Inverting the physical unimodular transform recovers the original
integer basis coordinates, including the factor-indicator coordinates. -/
theorem integerReducedCoordinates_transform
    (initial matrix transform : Generated.StrictRecombine.LLLMatrix)
    (hrel : LLLTransformRel initial matrix transform)
    (hunimodular : LLLTransformUnimodular transform)
    (coefficients : Fin initial.size → ZZ) :
    Matrix.vecMul (integerReducedCoordinates initial transform coefficients)
      (basisPrefixMatrix transform initial.size initial.size) = coefficients := by
  have hunit := hunimodular.2
  rw [hrel.2.1] at hunit
  unfold integerReducedCoordinates
  rw [Matrix.vecMul_vecMul, Matrix.nonsing_inv_mul _ hunit, Matrix.vecMul_one]

/-- Equal columns of the retained physical transform agree on every short
vector's original integer coordinates. No retained-space oracle is used. -/
theorem retained_transform_equal_original_coordinates
    (initial : Generated.StrictRecombine.LLLMatrix)
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrel : LLLTransformRel initial state.matrix state.transform)
    (hunimodular : LLLTransformUnimodular state.transform)
    (coefficients : Fin initial.size → ZZ) (bound : ZZ)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix initial initial.size initial.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix initial initial.size initial.size) column : ZZ) : QQ)) ≤
      (bound : QQ))
    (left right : Fin initial.size)
    (hequal : ∀ row, row ∈ Generated.StrictVanHoeij.shortSpanFromNorms state.norms bound →
      (state.transform[row]!)[left.val]! = (state.transform[row]!)[right.val]!) :
    coefficients left = coefficients right := by
  let reduced := integerReducedCoordinates initial state.transform coefficients
  let rows := Generated.StrictVanHoeij.shortSpanFromNorms state.norms bound
  have hrowCount : initial.size ≤ state.matrix.size := by rw [hrel.1]
  have hquadratic : dotProduct (fun i => (reduced i : QQ))
      (Matrix.mulVec (gramPrefixMatrixQQ state.matrix initial.size)
        (fun i => (reduced i : QQ))) ≤ (bound : QQ) := by
    rw [gramPrefix_quadratic_eq_vecMul_norm, hrel.1]
    have hreconstruct := integerReducedCoordinates_reconstruct initial state.matrix
      state.transform hrel hunimodular coefficients
    change Matrix.vecMul reduced _ = _ at hreconstruct
    rw [hreconstruct]
    exact hshort
  obtain ⟨_, hentries, hremoved⟩ := shortSpanFromNorms_spec state.norms bound
  have hzero := short_vector_coordinates_vanish_in_removed_suffix state hvalid
    initial.size hrowCount reduced rows.size (bound : QQ)
    (fun i hge => hremoved i.val (by rw [hvalid.norms_size, hrel.1]; exact i.isLt) hge)
    hquadratic
  have hrecover := integerReducedCoordinates_transform initial state.matrix
    state.transform hrel hunimodular coefficients
  change Matrix.vecMul reduced _ = coefficients at hrecover
  rw [← hrecover]
  simp only [Matrix.vecMul, dotProduct]
  apply Finset.sum_congr rfl
  intro i _
  by_cases hi : i.val < rows.size
  · have hentry := hentries i.val hi
    change rows[i.val]! = i.val at hentry
    have hmem : i.val ∈ rows := by
      rw [← hentry, getElem!_pos rows i.val hi]
      exact Array.getElem_mem hi
    change reduced i * (state.transform[i.val]!)[left.val]! =
      reduced i * (state.transform[i.val]!)[right.val]!
    rw [hequal i.val hmem]
  · have hz := hzero i (Nat.le_of_not_gt hi)
    simp only [hz, zero_mul]

/-- Concrete LLL plus the strict retention helper: equality tested on the
returned transform forces equality of a short initial vector's coordinates. -/
theorem lllReduceWithNorms_equal_original_coordinates
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (coefficients : Fin matrix.size → ZZ)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤
      (bound : QQ))
    (left right : Fin matrix.size)
    (hequal : ∀ row, row ∈ Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound →
      (output.2.1[row]!)[left.val]! = (output.2.1[row]!)[right.val]!) :
    coefficients left = coefficients right := by
  obtain ⟨state, hv, _, ht, hn, hrel, hu⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  apply retained_transform_equal_original_coordinates matrix state hv hrel hu
    coefficients bound hshort left right
  simpa only [ht, hn] using hequal

theorem lllReduceWithNorms_candidate_columns_valid
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (factorCount : Nat)
    (hfactorCount : factorCount ≤ matrix.size)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    CandidateColumnsValid output.2.1
      (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) factorCount := by
  obtain ⟨state, hv, _, ht, hn, hrel, _⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  rw [← ht, ← hn]
  obtain ⟨hsize, hentries, _⟩ := shortSpanFromNorms_spec state.norms bound
  intro position hposition
  dsimp only
  rw [hentries position hposition]
  have hp : position < state.transform.size := by
    rw [hrel.2.1]
    rw [hv.norms_size, hrel.1] at hsize
    omega
  refine ⟨hp, ?_⟩
  rw [getElem!_pos state.transform position hp, hrel.2.2.1 position hp]
  exact hfactorCount

/-- The actual candidate groups are constant on the coordinates of every
short initial lattice vector. This is the grouping part of the certificate,
before interpreting those coordinates as true-factor indicators. -/
theorem lllReduceWithNorms_class_members_coordinates
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (factorCount : Nat)
    (hfactorCount : factorCount ≤ matrix.size)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (coefficients : Fin matrix.size → ZZ)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤
      (bound : QQ))
    (candidates : Array (Array Int32))
    (hnonempty : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = false)
    (hextract : Generated.StrictRecombine.extractCandidates
      (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
      output.2.1 factorCount = .ok candidates)
    (classId : Nat) (leftValue rightValue : Int32)
    (hleftMember : CandidateMember candidates classId leftValue)
    (hrightMember : CandidateMember candidates classId rightValue) :
    ∃ (left right : Nat) (hl : left < factorCount) (hr : right < factorCount),
      leftValue = left.toUInt32.toInt32 ∧ rightValue = right.toUInt32.toInt32 ∧
      coefficients ⟨left, lt_of_lt_of_le hl hfactorCount⟩ =
        coefficients ⟨right, lt_of_lt_of_le hr hfactorCount⟩ := by
  let rows := Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound
  obtain ⟨left, right, hl, hr, hlv, hrv, hcompare⟩ :=
    extractCandidates_class_members_equal rows output.2.1 factorCount candidates hnonempty
      (lllReduceWithNorms_candidate_columns_valid matrix hinput bound factorCount
        hfactorCount output hrun) hextract classId leftValue rightValue hleftMember hrightMember
  refine ⟨left, right, hl, hr, hlv, hrv, ?_⟩
  apply lllReduceWithNorms_equal_original_coordinates matrix hinput bound output hrun
    coefficients hshort
  intro row hrow
  obtain ⟨position, hp, heq⟩ := Array.mem_iff_getElem.mp hrow
  have hpoint := candidateColumnsEqual_true_sound output.2.1 rows left right 0 hcompare
    position (Nat.zero_le _) hp
  dsimp only at hpoint
  have hrowEq : rows[position]! = row := by
    simpa only [getElem!_pos rows position hp] using heq
  rw [hrowEq] at hpoint
  exact hpoint.2.2.2

/-- Physical helper output suffices: no caller-supplied suffix bound remains.
Every nonzero coordinate of a short lattice vector occurs in the returned
row-index array. This connects the strict iota execution to safe removal. -/
theorem shortSpanFromNorms_preserves_coordinates
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (coordinates : Fin state.matrix.size → ZZ) (bound : ZZ)
    (hshort : dotProduct (fun position => (coordinates position : QQ))
        (Matrix.mulVec (gramPrefixMatrixQQ state.matrix state.matrix.size)
          (fun position => (coordinates position : QQ))) ≤ (bound : QQ)) :
    ∀ index : Fin state.matrix.size, coordinates index ≠ 0 →
      index.val ∈ Generated.StrictVanHoeij.shortSpanFromNorms state.norms bound := by
  let rows := Generated.StrictVanHoeij.shortSpanFromNorms state.norms bound
  obtain ⟨hsize, hentries, hremoved⟩ := shortSpanFromNorms_spec state.norms bound
  have hzero := short_vector_coordinates_vanish_in_removed_suffix state hvalid
    state.matrix.size le_rfl coordinates rows.size (bound : QQ)
    (fun index hge => hremoved index.val (by simpa [hvalid.norms_size] using index.isLt) hge)
    hshort
  intro index hne
  have hlt : index.val < rows.size := by
    by_contra hnot
    exact hne (hzero index (Nat.le_of_not_gt hnot))
  have hentry : rows[index.val]'hlt = index.val := by
    have h := hentries index.val hlt
    change rows[index.val]! = index.val at h
    simpa only [getElem!_pos rows index.val hlt] using h
  rw [← hentry]
  exact Array.getElem_mem hlt

/-- Composition for the actual exported LLL tuple and the actual selected
row array. The caller only supplies a short integer combination of the
returned basis; no GS-state or removed-row premise is required. -/
theorem lllReduceWithNorms_shortSpan_preserves
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (coordinates : Fin output.1.size → ZZ)
    (hshort : dotProduct (fun position => (coordinates position : QQ))
        (Matrix.mulVec (gramPrefixMatrixQQ output.1 output.1.size)
          (fun position => (coordinates position : QQ))) ≤ (bound : QQ)) :
    ∀ index : Fin output.1.size, coordinates index ≠ 0 →
      index.val ∈ Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound := by
  obtain ⟨state, hvalid, hmatrix, htransform, hnorms⟩ :=
    lllReduceWithNorms_valid matrix hinput bound output hrun
  rcases output with ⟨reduced, transform, rows, norms⟩
  dsimp only at hmatrix htransform hnorms
  subst reduced
  subst norms
  exact shortSpanFromNorms_preserves_coordinates state hvalid coordinates bound hshort

/-- Equality of two factor-coordinate columns propagates to every vector in
the retained rational row span, including each preserved factor indicator. -/
theorem retained_span_equal_columns
    {rows columns : Nat} (basis : Fin rows → Fin columns → QQ)
    (left right : Fin columns)
    (hequal : ∀ row, basis row left = basis row right)
    (value : Fin columns → QQ)
    (hvalue : value ∈ Submodule.span QQ (Set.range basis)) :
    value left = value right := by
  induction hvalue using Submodule.span_induction with
  | mem value hmem =>
      obtain ⟨row, rfl⟩ := hmem
      exact hequal row
  | zero => rfl
  | add x y _ _ hx hy => simpa only [Pi.add_apply, hx, hy]
  | smul a x _ hx => simpa only [Pi.smul_apply, smul_eq_mul, hx]

/-- Preserving every true factor indicator prevents an equal-column block
from merging modular factors belonging to different integer irreducibles. -/
theorem equal_columns_same_true_factor
    {rows columns factors : Nat} (basis : Fin rows → Fin columns → QQ)
    (owner : Fin columns → Fin factors)
    (hpreserved : ∀ factor : Fin factors,
      (fun index => if owner index = factor then (1 : QQ) else 0) ∈
        Submodule.span QQ (Set.range basis))
    (left right : Fin columns)
    (hequal : ∀ row, basis row left = basis row right) :
    owner left = owner right := by
  have h := retained_span_equal_columns basis left right hequal
    (fun index => if owner index = owner left then (1 : QQ) else 0)
    (hpreserved (owner left))
  by_contra hne
  simp [Ne.symm hne] at h

/-- A recovered integer factor's modular support is a union of true factor
supports.  If it is also a nonempty equal-column block, it contains exactly
one true factor.  The hypotheses concern supports, not irreducibility of
the candidate output. -/
theorem complete_block_is_one_true_factor
    {indices factors : Type*} (owner : indices → factors)
    (block : Set indices) (selected : Set factors)
    (hblock : block.Nonempty)
    (hconstant : ∀ left ∈ block, ∀ right ∈ block, owner left = owner right)
    (hunion : ∀ index, index ∈ block ↔ owner index ∈ selected) :
    ∃ factor, block = {index | owner index = factor} := by
  obtain ⟨first, hfirst⟩ := hblock
  refine ⟨owner first, Set.ext fun index => ?_⟩
  constructor
  · intro hindex
    exact hconstant index hindex first hfirst
  · intro hindex
    apply (hunion index).2
    change owner index = owner first at hindex
    rw [hindex]
    exact (hunion first).1 hfirst

end Refinement.StrictRecombine
