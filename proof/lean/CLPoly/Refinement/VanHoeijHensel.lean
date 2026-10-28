import CLPoly.Refinement.VanHoeijCld
import CLPoly.Refinement.FactorZZ

set_option autoImplicit false

namespace Refinement.StrictFactorZZ

open CLPoly.Math

theorem legalCombination_indices_nodup
    {upper : Nat} {indices : Array Nat}
    (hlegal : StrictRecombine.LegalCombination upper indices.size indices) :
    indices.toList.Nodup := by
  rw [List.nodup_iff_injective_getElem]
  intro left right heq
  have hl : left.val < indices.size := by simpa using left.isLt
  have hr : right.val < indices.size := by simpa using right.isLt
  have heq' : indices[left.val]! = indices[right.val]! := by
    simpa only [getElem!_pos indices left.val hl, getElem!_pos indices right.val hr,
      Array.getElem_toList] using heq
  apply Fin.ext
  by_contra hne
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · have hgap := hlegal.2.1 left.val right.val hl hr hlt
    omega
  · have hgap := hlegal.2.1 right.val left.val hr hl hgt
    omega

/-- The actual Hensel output supplies a physical modular support for every
integer divisor. Nonmonic factors are handled by cancelling a unit scalar
over the prime-power ring, never by assuming that ring is a field. -/
theorem hensel_divisor_physical_support
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (factor complement : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = factor * complement) :
    ∃ (indices : Array Nat) (exponent : Nat),
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      0 < exponent ∧ output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
      Polynomial.map (Int.castRingHom (ZMod (selection.prime.toNat ^ exponent))) factor =
        Polynomial.C (factor.leadingCoeff : ZMod (selection.prime.toNat ^ exponent)) *
          (indices.toList.map fun index =>
            StrictHensel.toPolyMod (selection.prime.toNat ^ exponent) output.1[index]!).prod := by
  have hleading := hselection.goodPrime.lc_nonzero
  have hmul : ((factor.leadingCoeff : ZMod selection.prime.toNat) *
      (complement.leadingCoeff : ZMod selection.prime.toNat)) ≠ 0 := by
    simpa only [hfactor, Polynomial.leadingCoeff_mul, Int.cast_mul] using hleading
  have hfLeading := left_ne_zero_of_mul hmul
  have hcLeading := right_ne_zero_of_mul hmul
  have hfMapLeading := Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero
    (Int.castRingHom (ZMod selection.prime.toNat)) hfLeading
  have hfMapNe : Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat)) factor ≠ 0 := by
    intro hz
    apply hfLeading
    simpa [hz] using hfMapLeading.symm
  obtain ⟨indices, hlegal, hassociated⟩ :=
    integer_divisor_mod_has_legal_hensel_candidate hcount hfactors hleadingSemantic
      hselection hentry factor ⟨complement, hfactor⟩
  obtain ⟨exponent, he, hm, hscaled⟩ :=
    henselCandidate_scaled_eq_divisor_mod_primePower hcount hfactors hleadingSemantic
      hselection hentry factor complement hfactor hfMapNe hfMapLeading indices
      hlegal hassociated
  refine ⟨indices, exponent, hlegal, he, hm, ?_⟩
  have hcUnit : IsUnit (Polynomial.C
      (complement.leadingCoeff : ZMod (selection.prime.toNat ^ exponent))) :=
    Polynomial.isUnit_C.mpr (intCast_isUnit_primePower_of_ne_zero_prime
      selection.prime.toNat exponent complement.leadingCoeff hcLeading)
  have hselected :
      Polynomial.map (Int.castRingHom (ZMod (selection.prime.toNat ^ exponent)))
        ((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
          SparsePolyZZ.toPoly).prod =
      (indices.toList.map fun index =>
        StrictHensel.toPolyMod (selection.prime.toNat ^ exponent) output.1[index]!).prod := by
    simp only [Polynomial.map_list_prod, List.map_map,
      StrictRecombine.selectSourceIndices, Function.comp_def, StrictHensel.toPolyMod]
    apply congrArg List.prod
    apply List.map_congr_left
    intro index hi
    obtain ⟨position, hp, heq⟩ := List.mem_iff_getElem.mp hi
    have hp' : position < indices.size := by simpa using hp
    have hindex : index < output.1.size := by
      have hb := hlegal.2.2 position hp'
      simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', heq] using hb
    rw [getElem!_pos output.1.toList index (by simpa using hindex),
      Array.getElem_toList hindex, getElem!_pos output.1 index hindex]
  rw [Polynomial.map_mul, Polynomial.map_C, hselected,
    Polynomial.map_mul, Polynomial.map_C] at hscaled
  have hlc : (SparsePolyZZ.toPoly f).leadingCoeff =
      factor.leadingCoeff * complement.leadingCoeff := by
    rw [hfactor, Polynomial.leadingCoeff_mul]
  rw [hlc, map_mul, Polynomial.C_mul] at hscaled
  apply hcUnit.mul_left_cancel
  calc
    _ = _ := hscaled.symm
    _ = _ := by
      change Polynomial.C (factor.leadingCoeff : ZMod (selection.prime.toNat ^ exponent)) *
        Polynomial.C (complement.leadingCoeff : ZMod (selection.prime.toNat ^ exponent)) * _ =
        Polynomial.C (complement.leadingCoeff : ZMod (selection.prime.toNat ^ exponent)) *
          (Polynomial.C (factor.leadingCoeff : ZMod (selection.prime.toNat ^ exponent)) * _)
      ring

/-- Combine the genuine Hensel support with the actual CLD outputs at their
returned modulus. No support equation or carry is an input assumption. -/
theorem hensel_cld_factor_coeff_carry
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (cld : Array SparsePolyZZ)
    (hcld : ∀ index, index < output.1.size →
      StrictRecombine.ModularCldValue f output.1[index]! cld[index]! output.2.toNat)
    (factor complement : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = factor * complement) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      ∀ degree, ∃ carry : Int,
        (indices.toList.map fun index =>
          Generated.StrictRecombine.sparseCoeff cld[index]! degree).sum =
            (complement * Polynomial.derivative factor).coeff degree + carry * output.2 := by
  obtain ⟨indices, exponent, hlegal, he, hm, hsupport⟩ :=
    hensel_divisor_physical_support hcount hfactors hleadingSemantic hselection hentry
      factor complement hfactor
  have hmod : output.2.toNat = selection.prime.toNat ^ exponent := by
    rw [hm]
    exact Int.toNat_natCast _
  have hindices : ∀ index ∈ indices.toList, index < output.1.size := by
    intro index hi
    obtain ⟨position, hp, heq⟩ := List.mem_iff_getElem.mp hi
    have hp' : position < indices.size := by simpa using hp
    have hb := hlegal.2.2 position hp'
    simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', heq] using hb
  have hmonic : ∀ index, index < output.1.size →
      (StrictHensel.toPolyMod output.2.toNat output.1[index]!).Monic := by
    intro index hi
    rw [getElem!_pos output.1 index hi]
    exact hentry.outputToPolyModMonic output.2.toNat index hi
  refine ⟨indices, hlegal, ?_⟩
  intro degree
  have hsupport' : Polynomial.map (Int.castRingHom (ZMod output.2.toNat)) factor =
      Polynomial.C (factor.leadingCoeff : ZMod output.2.toNat) *
        (indices.toList.map fun index => StrictHensel.toPolyMod output.2.toNat output.1[index]!).prod := by
    rw [hmod]
    exact hsupport
  obtain ⟨carry, hc⟩ := StrictRecombine.modularCld_factor_coeff_carry
    f output.1 cld output.2.toNat hcld hmonic indices.toList hindices
    factor complement _ (by simpa only [mul_comm] using hfactor) hsupport' degree
  exact ⟨carry, by simpa only [hmod, hm] using hc⟩

end Refinement.StrictFactorZZ
