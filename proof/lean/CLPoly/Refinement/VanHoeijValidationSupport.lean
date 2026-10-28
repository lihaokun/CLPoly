import CLPoly.Refinement.VanHoeijValidationProduct

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem validatePartition_factor_primitive
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, (SparsePolyZZ.toPoly factor).IsPrimitive := by
  have hp := validatePartition_product f lifted modulus groups before after hcanonical hdegree hrun
  intro factor hf
  apply Polynomial.isPrimitive_of_dvd hprimitive
  rw [← hp]
  exact List.dvd_prod (List.mem_map.mpr ⟨factor, hf, rfl⟩)

theorem partitionProductLoop_projection
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (product output : SparsePolyZZ) (degree finalDegree : Int64)
    (hrun : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus index
      product degree = .ok (output, finalDegree)) :
    Generated.StrictRecombine.trialProductLoop ⟨()⟩ group lifted modulus index product =
      .ok output := by
  induction hmeasure : group.size - index using Nat.strong_induction_on
      generalizing index product output degree finalDegree with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.partitionProductLoop] at hrun
    split at hrun
    next hi =>
      dsimp only at hrun
      split at hrun
      next hnegative => simp at hrun
      next hnonnegative =>
        split at hrun
        next hp =>
          cases hm : Generated.StrictRecombine.multiplyNormalizeModRaw product
              lifted[group[index].toInt64.toNat] modulus with
          | error fault => simp [hm, bind, Except.bind] at hrun
          | ok next =>
            simp only [hm, bind, Except.bind] at hrun
            have ht := ih (group.size - (index + 1)) (by omega)
              (index + 1) next output (degree + get_deg lifted[group[index].toInt64.toNat])
              finalDegree hrun rfl
            have hn : 0 ≤ group[index] := by
              rw [Int32.le_iff_toInt_le]
              rw [Int32.lt_iff_toInt_lt] at hnonnegative
              simpa using hnonnegative
            rw [Generated.StrictRecombine.trialProductLoop, dif_pos hi]
            dsimp only
            rw [dif_pos hn, dif_pos hp, hm]
            exact ht
        next hout => simp at hrun
    next hdone =>
      have he := (Prod.mk.inj (Except.ok.inj hrun)).1
      rw [Generated.StrictRecombine.trialProductLoop, dif_neg hdone, he]

theorem partitionProductLoop_modular
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus : Nat)
    (hm : 0 < modulus) (index : Nat) (product output : SparsePolyZZ)
    (degree finalDegree : Int64)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hrun : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus index
      product degree = .ok (output, finalDegree)) :
    StrictHensel.toPolyMod modulus output = StrictHensel.toPolyMod modulus product *
      SelectedProductMod modulus group lifted index := by
  exact trialProductLoop_refines ⟨()⟩ group lifted modulus hm index product output hvalid
    (partitionProductLoop_projection group lifted modulus index product output degree finalDegree hrun)

/-- Actual modular multiplication, symmetric recovery and primitive division
give the scaled support equation. Invertibility of the content is a separate
obligation; it is not silently cancelled here. -/
theorem partition_recovered_factor_scaled_support
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus base : Nat)
    (hm : 0 < modulus) (hb : 0 < base) (hdivides : base ∣ modulus) (leading content : ZZ)
    (product symmetric factor : SparsePolyZZ) (degree : Int64)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hproduct : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus 0
      #[(⟨0⟩, leading)] 0 = .ok (product, degree))
    (hsymmetric : Generated.StrictRecombine.symmetricModRaw product modulus = .ok symmetric)
    (hprimitive : Generated.StrictRecombine.primitiveRaw symmetric = .ok (content, factor)) :
    Polynomial.C (content : ZMod base) * StrictHensel.toPolyMod base factor =
      Polynomial.C (leading : ZMod base) * SelectedProductMod base group lifted 0 := by
  have hp := trialProductLoop_refines_of_dvd ⟨()⟩ group lifted modulus base hm hb hdivides
    0 #[(⟨0⟩, leading)] product hvalid
    (partitionProductLoop_projection group lifted modulus 0 _ product 0 degree hproduct)
  have hs := symmetricModRaw_toPolyMod_of_dvd product symmetric modulus base hm hb hdivides hsymmetric
  have hprimitiveMap := congrArg (Polynomial.map (Int.castRingHom (ZMod base)))
    (primitiveRaw_toPoly symmetric factor content hprimitive)
  change StrictHensel.toPolyMod base symmetric = _ at hprimitiveMap
  rw [Polynomial.map_mul, Polynomial.map_C] at hprimitiveMap
  change StrictHensel.toPolyMod base symmetric =
    Polynomial.C (content : ZMod base) * StrictHensel.toPolyMod base factor at hprimitiveMap
  rw [← hprimitiveMap, hs, hp]
  simp [StrictHensel.toPolyMod, SparsePolyZZ.toPoly, Polynomial.monomial_zero_left]

/-- A nonzero leading scalar and nonzero selected product force the actual
primitive content to be a unit at the selected prime. -/
theorem scaled_support_associated_at_prime
    (prime : Nat) [Fact (Nat.Prime prime)] (leading content : ZZ)
    (factor : SparsePolyZZ) (selected : Polynomial (ZMod prime))
    (hleading : (leading : ZMod prime) ≠ 0) (hselected : selected ≠ 0)
    (hscaled : Polynomial.C (content : ZMod prime) * StrictHensel.toPolyMod prime factor =
      Polynomial.C (leading : ZMod prime) * selected) :
    IsUnit (content : ZMod prime) ∧ Associated (StrictHensel.toPolyMod prime factor) selected := by
  have hright : Polynomial.C (leading : ZMod prime) * selected ≠ 0 :=
    mul_ne_zero (Polynomial.C_ne_zero.mpr hleading) hselected
  have hcontent : (content : ZMod prime) ≠ 0 := by
    intro hz
    rw [hz, Polynomial.C_0, zero_mul] at hscaled
    exact hright hscaled.symm
  have hc := isUnit_iff_ne_zero.mpr hcontent
  have hl := isUnit_iff_ne_zero.mpr hleading
  refine ⟨hc, ?_⟩
  have ha : Associated (Polynomial.C (content : ZMod prime) * StrictHensel.toPolyMod prime factor)
      (Polynomial.C (leading : ZMod prime) * selected) := hscaled ▸ Associated.refl _
  exact (associated_isUnit_mul_right_iff (Polynomial.isUnit_C.mpr hl)).mp
    ((associated_isUnit_mul_left_iff (Polynomial.isUnit_C.mpr hc)).mp ha)

theorem partition_recovered_factor_associated
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus prime : Nat)
    [Fact (Nat.Prime prime)] (remaining product symmetric factor : SparsePolyZZ)
    (content : ZZ) (degree : Int64)
    (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hirreducible : ∀ index (hi : index < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[index]))
    (hremaining : 0 < remaining.size)
    (hleading : (remaining[0].2 : ZMod prime) ≠ 0)
    (hproduct : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus 0
      #[(⟨0⟩, remaining[0].2)] 0 = .ok (product, degree))
    (hsymmetric : Generated.StrictRecombine.symmetricModRaw product modulus = .ok symmetric)
    (hprimitive : Generated.StrictRecombine.primitiveRaw symmetric = .ok (content, factor)) :
    Associated (StrictHensel.toPolyMod prime factor)
      (SelectedProductMod prime group lifted 0) := by
  exact validationRecoveredFactor_mod_associated_selected ⟨⟨()⟩⟩ group lifted
    modulus prime remaining product symmetric factor content hm hdivides hvalid
    hirreducible hremaining hleading
    (partitionProductLoop_projection group lifted modulus 0 _ product 0 degree hproduct)
    hsymmetric hprimitive

/-- Exact division preserves the good-prime leading coefficient invariant
for both the emitted factor and the unnormalised next remainder. -/
theorem partition_exact_division_leading_nonzero
    (prime : Nat) [Fact (Nat.Prime prime)] (remaining factor quotient : SparsePolyZZ)
    (hleading : ((SparsePolyZZ.toPoly remaining).leadingCoeff : ZMod prime) ≠ 0)
    (hdivide : Generated.StrictRecombine.exactDivmodRaw remaining factor =
      .ok (quotient, #[])) :
    ((SparsePolyZZ.toPoly factor).leadingCoeff : ZMod prime) ≠ 0 ∧
      ((SparsePolyZZ.toPoly quotient).leadingCoeff : ZMod prime) ≠ 0 := by
  have hp := exactDivmodRaw_toPoly remaining factor quotient #[] hdivide
  simp only [SparsePolyZZ.toPoly, Array.toList_empty, List.map_nil, List.sum_nil,
    add_zero] at hp
  change SparsePolyZZ.toPoly remaining =
    SparsePolyZZ.toPoly factor * SparsePolyZZ.toPoly quotient at hp
  have hl := congrArg Polynomial.leadingCoeff hp
  rw [Polynomial.leadingCoeff_mul] at hl
  have hc := congrArg (fun z : Int => (z : ZMod prime)) hl
  dsimp only at hc
  rw [Int.cast_mul] at hc
  exact mul_ne_zero_iff.mp (hc ▸ hleading)

end Refinement.StrictRecombine
