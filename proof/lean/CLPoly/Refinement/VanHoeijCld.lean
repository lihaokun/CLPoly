import CLPoly.Refinement.VanHoeijTruncation

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Properties of the actual CLD output needed by truncation and the lattice.
No integer divisibility of the lifted divisor is assumed. -/
structure ModularCldValue (f divisor output : SparsePolyZZ) (modulus : Nat) : Prop where
  canonical : StrictPolynomialMod.SparsePolyZZCanonical output
  semantic : StrictHensel.toPolyMod modulus output =
    (StrictHensel.toPolyMod modulus f /ₘ StrictHensel.toPolyMod modulus divisor) *
      Polynomial.derivative (StrictHensel.toPolyMod modulus divisor)
  bounded : ∀ degree, ((SparsePolyZZ.toPoly output).coeff degree).natAbs * 2 ≤ modulus

theorem generatedCldElement_valid
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
    ModularCldValue f divisor cld modulus := by
  have hproductCanonical := multiplyNormalizeModRaw_canonical quotient derivativeMod product
    (modulus : ZZ) hproduct
  exact ⟨symmetricModRaw_canonical product cld modulus hmodulus hproductCanonical hsymmetric,
    generatedCldElement_modular termination f divisor quotient remainder derivativeRaw
      derivativeMod product cld modulus hmodulus hdivisorSize hdivisorDegree hdivisorHead
      hfBound hinvert hdivision hremainder hderivative hmod hproduct hsymmetric hmonic,
    symmetricModRaw_toPoly_coefficients_bounded product cld modulus hmodulus
      hproductCanonical hsymmetric⟩

theorem cldPolysLoop_modular_valid
    (ops : Generated.StrictRecombine.CldRawOps)
    (f : SparsePolyZZ) (active : Array SparsePolyZZ) (modulus : Nat)
    (hmodulus : 0 < modulus) (hfBound : StrictHensel.DegreesBound f)
    (hsize : ∀ index (hindex : index < active.size),
      0 < active[index].size)
    (hdegree : ∀ index (hindex : index < active.size),
      active[index][0]!.1.deg < 2 ^ 63)
    (hhead : ∀ index (hindex : index < active.size),
      StrictHensel.HeadDominates active[index])
    (hinvert : ∀ index (hindex : index < active.size),
      (ZZ.invert 0 active[index][0]!.2 (modulus : Int)).1 = true)
    (hmonic : ∀ index (hindex : index < active.size),
      (StrictHensel.toPolyMod modulus active[index]).Monic)
    (index : Nat) (result output : Array SparsePolyZZ)
    (hindexLe : index ≤ active.size) (hresultSize : result.size = index)
    (hresult : ∀ position (hposition : position < index),
      ModularCldValue f active[position]! result[position]! modulus)
    (hrun : Generated.StrictRecombine.cldPolysLoop ops f active
      (modulus : Int) index result = .ok output) :
    output.size = active.size ∧
      ∀ position (hposition : position < active.size),
        ModularCldValue f active[position]! output[position]! modulus := by
  induction hmeasure : active.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.cldPolysLoop] at hrun
      split at hrun
      next hindex =>
        cases hdivision : Generated.StrictHensel.__upoly_divmod_mod_raw_ir
            ops.divmodTermination f active[index] (modulus : Int) with
        | error fault => simp [hdivision] at hrun
        | ok division =>
          rcases division with ⟨quotient, remainder⟩
          simp only [hdivision] at hrun
          split at hrun
          next hremainder =>
            cases hderivative : Generated.StrictRecombine.derivativeZZRaw
                active[index] with
            | error fault => simp [hderivative] at hrun
            | ok derivativeRaw =>
              simp only [hderivative] at hrun
              cases hmod : Generated.StrictRecombine.modCoeffLoop derivativeRaw
                  (modulus : Int) 0 #[] with
              | error fault => simp [hmod] at hrun
              | ok derivativeMod =>
                simp only [hmod] at hrun
                cases hproduct :
                    Generated.StrictRecombine.multiplyNormalizeModRaw quotient
                      derivativeMod (modulus : Int) with
                | error fault => simp [hproduct] at hrun
                | ok product =>
                  simp only [hproduct] at hrun
                  cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                      product (modulus : Int) with
                  | error fault => simp [hsymmetric] at hrun
                  | ok cld =>
                    simp only [hsymmetric] at hrun
                    have hcld := generatedCldElement_valid
                      ops.divmodTermination f active[index] quotient remainder
                      derivativeRaw derivativeMod product cld modulus hmodulus
                      (hsize index hindex) (hdegree index hindex)
                      (hhead index hindex) hfBound (hinvert index hindex)
                      hdivision hremainder hderivative hmod hproduct hsymmetric
                      (hmonic index hindex)
                    have hnextResult : ∀ position
                        (hposition : position < index + 1),
                        ModularCldValue f active[position]! (result.push cld)[position]! modulus := by
                      intro position hposition
                      by_cases hbefore : position < index
                      · rw [getElem!_pos (result.push cld) position (by
                            simp [hresultSize]; omega)]
                        rw [Array.getElem_push_lt (by
                          simpa [hresultSize] using hbefore)]
                        rw [← getElem!_pos result position (by
                          simpa [hresultSize] using hbefore)]
                        exact hresult position hbefore
                      · have hpositionEq : position = index := by omega
                        subst position
                        rw [getElem!_pos (result.push cld) index (by
                          simp [hresultSize])]
                        simpa [Array.getElem_push, hresultSize,
                          getElem!_pos active index hindex] using hcld
                    exact ih (active.size - (index + 1)) (by omega)
                      (index + 1) (result.push cld) output (by omega)
                      (by simp [hresultSize]) hnextResult hrun rfl
          next hremainder => contradiction
      next hindex =>
        have hout := Except.ok.inj hrun
        subst output
        refine ⟨by omega, ?_⟩
        intro position hposition
        exact hresult position (by omega)

theorem cldPolys_modular_valid
    (ops : Generated.StrictRecombine.CldRawOps)
    (f : SparsePolyZZ) (active output : Array SparsePolyZZ) (modulus : Nat)
    (hmodulus : 0 < modulus) (hfBound : StrictHensel.DegreesBound f)
    (hsize : ∀ index (hindex : index < active.size), 0 < active[index].size)
    (hdegree : ∀ index (hindex : index < active.size),
      active[index][0]!.1.deg < 2 ^ 63)
    (hhead : ∀ index (hindex : index < active.size),
      StrictHensel.HeadDominates active[index])
    (hinvert : ∀ index (hindex : index < active.size),
      (ZZ.invert 0 active[index][0]!.2 (modulus : Int)).1 = true)
    (hmonic : ∀ index (hindex : index < active.size),
      (StrictHensel.toPolyMod modulus active[index]).Monic)
    (hrun : Generated.StrictRecombine.cldPolys ops f active (modulus : Int) =
      .ok output) :
    output.size = active.size ∧
      ∀ index (hindex : index < active.size),
        ModularCldValue f active[index]! output[index]! modulus := by
  unfold Generated.StrictRecombine.cldPolys at hrun
  exact cldPolysLoop_modular_valid ops f active modulus hmodulus hfBound hsize
    hdegree hhead hinvert hmonic 0 #[] output (by omega) (by simp)
    (by simp) hrun

theorem ModularCldValue.sparseCoeff_abs_le
    {f divisor output : SparsePolyZZ} {modulus : Nat}
    (hvalid : ModularCldValue f divisor output modulus) (degree : Nat) :
    |Generated.StrictRecombine.sparseCoeff output degree| ≤ (modulus : Int) := by
  rw [sparseCoeff_eq_coeff output hvalid.canonical degree]
  have hb := hvalid.bounded degree
  have hle : ((SparsePolyZZ.toPoly output).coeff degree).natAbs ≤ modulus := by omega
  have hcast : (((SparsePolyZZ.toPoly output).coeff degree).natAbs : Int) ≤ modulus := by
    exact_mod_cast hle
  simpa only [Int.natCast_natAbs] using hcast

/-- Sum the actual CLD outputs over a modular support. The complement and
product equation are algebraic inputs, not assumptions about output factors. -/
theorem modularCld_support_sum
    (f : SparsePolyZZ) (active output : Array SparsePolyZZ) (modulus : Nat)
    (hvalid : ∀ index, index < active.size →
      ModularCldValue f active[index]! output[index]! modulus)
    (hmonic : ∀ index, index < active.size →
      (StrictHensel.toPolyMod modulus active[index]!).Monic)
    (indices : List Nat) (hindices : ∀ index ∈ indices, index < active.size)
    (complement : Polynomial (ZMod modulus))
    (hproduct : StrictHensel.toPolyMod modulus f = complement *
      (indices.map fun index => StrictHensel.toPolyMod modulus active[index]!).prod) :
    (indices.map fun index => StrictHensel.toPolyMod modulus output[index]!).sum =
      complement * Polynomial.derivative
        (indices.map fun index => StrictHensel.toPolyMod modulus active[index]!).prod := by
  calc
    _ = (indices.map fun index =>
        (StrictHensel.toPolyMod modulus f /ₘ StrictHensel.toPolyMod modulus active[index]!) *
          Polynomial.derivative (StrictHensel.toPolyMod modulus active[index]!)).sum := by
      congr 1
      apply List.map_congr_left
      intro index hi
      exact (hvalid index (hindices index hi)).semantic
    _ = _ := by
      rw [hproduct]
      have hm : ∀ factor ∈ (indices.map fun index =>
          StrictHensel.toPolyMod modulus active[index]!), factor.Monic := by
        intro factor hf
        obtain ⟨index, hi, rfl⟩ := List.mem_map.mp hf
        exact hmonic index (hindices index hi)
      simpa only [List.map_map, Function.comp_def] using
        cld_sum_monic_product (indices.map fun index =>
          StrictHensel.toPolyMod modulus active[index]!) complement hm

/-- For any true integer divisor with the stated Hensel support, the actual
CLD subset sum equals its integer logarithmic derivative modulo m. The
integer divisor need not be monic or irreducible. -/
theorem modularCld_support_sum_of_factor
    (f : SparsePolyZZ) (active output : Array SparsePolyZZ) (modulus : Nat)
    (hvalid : ∀ index, index < active.size →
      ModularCldValue f active[index]! output[index]! modulus)
    (hmonic : ∀ index, index < active.size →
      (StrictHensel.toPolyMod modulus active[index]!).Monic)
    (indices : List Nat) (hindices : ∀ index ∈ indices, index < active.size)
    (factor complement : Polynomial Int) (scalar : ZMod modulus)
    (hfactor : SparsePolyZZ.toPoly f = complement * factor)
    (hsupport : Polynomial.map (Int.castRingHom (ZMod modulus)) factor = Polynomial.C scalar *
      (indices.map fun index => StrictHensel.toPolyMod modulus active[index]!).prod) :
    (indices.map fun index => StrictHensel.toPolyMod modulus output[index]!).sum =
      Polynomial.map (Int.castRingHom (ZMod modulus)) (complement * Polynomial.derivative factor) := by
  have hproduct : StrictHensel.toPolyMod modulus f =
      (Polynomial.map (Int.castRingHom (ZMod modulus)) complement * Polynomial.C scalar) *
        (indices.map fun index => StrictHensel.toPolyMod modulus active[index]!).prod := by
    simp only [StrictHensel.toPolyMod, hfactor, Polynomial.map_mul, hsupport, mul_assoc]
  rw [modularCld_support_sum f active output modulus hvalid hmonic indices hindices _ hproduct]
  rw [Polynomial.map_mul, ← Polynomial.derivative_map, hsupport,
    Polynomial.derivative_mul, Polynomial.derivative_C, zero_mul, zero_add]
  ring

/-- Extract an integer carry from an equality of actual CLD outputs modulo m.
The sparse coefficient lookup is justified by canonicality, rather than by
assuming arbitrary term arrays have unique degrees. -/
theorem modularCld_coeff_carry
    (output : Array SparsePolyZZ) (modulus : Nat) (indices : List Nat)
    (actual : Polynomial Int)
    (hcanonical : ∀ index ∈ indices,
      StrictPolynomialMod.SparsePolyZZCanonical output[index]!)
    (hsum : (indices.map fun index =>
      StrictHensel.toPolyMod modulus output[index]!).sum =
        Polynomial.map (Int.castRingHom (ZMod modulus)) actual)
    (degree : Nat) :
    ∃ carry : Int,
      (indices.map fun index =>
        Generated.StrictRecombine.sparseCoeff output[index]! degree).sum =
          actual.coeff degree + carry * (modulus : Int) := by
  have hcoeff :
      (((indices.map fun index => Generated.StrictRecombine.sparseCoeff
        output[index]! degree).sum : Int) : ZMod modulus) =
          (actual.coeff degree : ZMod modulus) := by
    have hlookup : ∀ index ∈ indices,
        (Generated.StrictRecombine.sparseCoeff output[index]! degree : ZMod modulus) =
          (StrictHensel.toPolyMod modulus output[index]!).coeff degree := by
      intro index hi
      rw [sparseCoeff_eq_coeff _ (hcanonical index hi)]
      simp [StrictHensel.toPolyMod]
    have hleft :
        (((indices.map fun index => Generated.StrictRecombine.sparseCoeff
          output[index]! degree).sum : Int) : ZMod modulus) =
        ((indices.map fun index => StrictHensel.toPolyMod modulus output[index]!).sum).coeff
          degree := by
      clear hsum hcanonical
      induction indices with
      | nil => simp
      | cons index rest ih =>
        simp only [List.map_cons, List.sum_cons, Int.cast_add, Polynomial.coeff_add]
        rw [hlookup index (by simp)]
        congr 1
        exact ih (by intro i hi; exact hlookup i (by simp [hi]))
    rw [hleft, hsum, Polynomial.coeff_map]
    rfl
  have hdiv := (ZMod.intCast_eq_intCast_iff_dvd_sub
    (actual.coeff degree)
    ((indices.map fun index => Generated.StrictRecombine.sparseCoeff
      output[index]! degree).sum) modulus).mp hcoeff.symm
  obtain ⟨carry, hcarry⟩ := hdiv
  exact ⟨carry, by nlinarith [hcarry]⟩

/-- Compose the physical CLD semantics with a true factor's modular support.
In particular, the carry required by the lattice is derived, not supplied. -/
theorem modularCld_factor_coeff_carry
    (f : SparsePolyZZ) (active output : Array SparsePolyZZ) (modulus : Nat)
    (hvalid : ∀ index, index < active.size →
      ModularCldValue f active[index]! output[index]! modulus)
    (hmonic : ∀ index, index < active.size →
      (StrictHensel.toPolyMod modulus active[index]!).Monic)
    (indices : List Nat) (hindices : ∀ index ∈ indices, index < active.size)
    (factor complement : Polynomial Int) (scalar : ZMod modulus)
    (hfactor : SparsePolyZZ.toPoly f = complement * factor)
    (hsupport : Polynomial.map (Int.castRingHom (ZMod modulus)) factor = Polynomial.C scalar *
      (indices.map fun index => StrictHensel.toPolyMod modulus active[index]!).prod)
    (degree : Nat) :
    ∃ carry : Int,
      (indices.map fun index =>
        Generated.StrictRecombine.sparseCoeff output[index]! degree).sum =
          (complement * Polynomial.derivative factor).coeff degree +
            carry * (modulus : Int) := by
  apply modularCld_coeff_carry output modulus indices
    (complement * Polynomial.derivative factor)
  · intro index hi
    exact (hvalid index (hindices index hi)).canonical
  · exact modularCld_support_sum_of_factor f active output modulus hvalid hmonic
      indices hindices factor complement scalar hfactor hsupport

end Refinement.StrictRecombine
