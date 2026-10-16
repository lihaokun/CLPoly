import CLPoly.Refinement.VanHoeijPartitionDegree
import CLPoly.Refinement.VanHoeijPrimitiveRecovery
import CLPoly.Refinement.VanHoeijRecoveryHeight

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- One actual validation iteration recovers the positive primitive factor,
passes all degree/divisibility checks, and recurses with the exact complement.
No successful validation trace is assumed. -/
theorem validatePartitionLoop_step_succeeds
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hm : 0 < modulus) (hd : prime ∣ modulus)
    (lifted : Array SparsePolyZZ) (groups : Array (Array Int32)) (index : Nat)
    (hi : index < groups.size) (remaining : SparsePolyZZ) (previous : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hfits : (SparsePolyZZ.toPoly remaining).natDegree < 2 ^ 63)
    (factor complement : Polynomial Int)
    (hremaining : SparsePolyZZ.toPoly remaining = factor * complement)
    (hprimitive : factor.IsPrimitive)
    (hfactorLeading : (0 : Int) < factor.leadingCoeff)
    (hcomplementLeading : (0 : Int) < complement.leadingCoeff)
    (hdegree : 0 < factor.natDegree)
    (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hvalid : CandidateIndicesValid groups[index] (Array.replicate lifted.size false))
    (hliftedCanonical : ∀ i, i < lifted.size → StrictPolynomialMod.SparsePolyZZCanonical lifted[i]!)
    (hliftedNonempty : ∀ i, i < lifted.size → 0 < lifted[i]!.size)
    (hmonic : ∀ i, i < lifted.size → (SparsePolyZZ.toPoly lifted[i]!).Monic)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) * SelectedProductMod modulus groups[index] lifted 0)
    (hscaled : ∀ degree,
      ((Polynomial.C complement.leadingCoeff * factor).coeff degree).natAbs * 2 < modulus) :
    ∃ recovered quotient,
      SparsePolyZZ.toPoly recovered = factor ∧ SparsePolyZZ.toPoly quotient = complement ∧
      StrictPolynomialMod.SparsePolyZZCanonical recovered ∧
      StrictPolynomialMod.SparsePolyZZCanonical quotient ∧
      Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups index remaining previous =
        Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups (index + 1)
          quotient (previous.push recovered) := by
  have hfactorNe := hprimitive.ne_zero
  have hcomplementNe : complement ≠ 0 :=
    Polynomial.leadingCoeff_ne_zero.mp (ne_of_gt hcomplementLeading)
  have hremainingNe : SparsePolyZZ.toPoly remaining ≠ 0 := by
    rw [hremaining]; exact mul_ne_zero hfactorNe hcomplementNe
  have hr : 0 < remaining.size := by
    by_contra hn
    have he : remaining = #[] := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hn)
    simp [he, SparsePolyZZ.toPoly] at hremainingNe
  have hfactorFits : factor.natDegree < 2 ^ 63 :=
    lt_of_le_of_lt (Polynomial.natDegree_le_of_dvd ⟨complement, hremaining⟩ hremainingNe) hfits
  obtain ⟨product, hproduct, hdegreePositive⟩ := partitionProductLoop_succeeds_with_factor_degree
    prime modulus hm hd groups[index] lifted #[(⟨0⟩, remaining[0].2)] factor hleading hvalid
    hliftedCanonical hliftedNonempty hmonic hsupport hdegree hfactorFits
  have hinitial : StrictPolynomialMod.SparsePolyZZCanonical #[(⟨0⟩, remaining[0].2)] := by
    constructor
    · simp
    · intro term ht
      simp at ht
      subst term
      exact hcanonical.2 remaining[0] (Array.getElem_mem_toList hr)
  have hproductCanonical := trialProductLoop_canonical ⟨()⟩ groups[index] lifted modulus 0
    #[(⟨0⟩, remaining[0].2)] product hinitial
    (partitionProductLoop_projection groups[index] lifted modulus 0 _ product 0 _ hproduct)
  have hmod := partitionProductLoop_modular groups[index] lifted modulus hm 0
    #[(⟨0⟩, remaining[0].2)] product 0 _ hvalid hproduct
  have hinitialMod : StrictHensel.toPolyMod modulus #[(⟨0⟩, remaining[0].2)] =
      Polynomial.C (remaining[0].2 : ZMod modulus) := by
    simp [StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
  rw [hinitialMod, ← sparsePolyZZ_leadingCoeff_eq_head remaining hcanonical hr,
    hremaining, Polynomial.leadingCoeff_mul, Int.cast_mul, Polynomial.C_mul] at hmod
  have hcongruent : (SparsePolyZZ.toPoly product).map (Int.castRingHom (ZMod modulus)) =
      (Polynomial.C complement.leadingCoeff * factor).map (Int.castRingHom (ZMod modulus)) := by
    change StrictHensel.toPolyMod modulus product = _
    rw [hmod, Polynomial.map_mul, Polynomial.map_C, hsupport]
    change Polynomial.C (factor.leadingCoeff : ZMod modulus) *
        Polynomial.C (complement.leadingCoeff : ZMod modulus) *
        SelectedProductMod modulus groups[index] lifted 0 =
      Polynomial.C (complement.leadingCoeff : ZMod modulus) *
        (Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[index] lifted 0)
    ring
  obtain ⟨symmetric, hsymmetric⟩ := symmetricModRaw_complete product modulus (by exact_mod_cast hm)
  have hsymmetricCanonical := symmetricModRaw_canonical product symmetric modulus hm
    hproductCanonical hsymmetric
  have hsymmetricTarget := symmetricModRaw_recovers_strictly_bounded_target product symmetric
    (Polynomial.C complement.leadingCoeff * factor) modulus hm hproductCanonical
    hsymmetric hcongruent hscaled
  obtain ⟨content, recovered, hrecover⟩ := primitiveRaw_complete symmetric hsymmetricCanonical
  obtain ⟨hcontent, hfactor⟩ := primitiveRaw_recovers_positive_scaled_target symmetric recovered
    content complement.leadingCoeff factor hsymmetricCanonical hcomplementLeading hprimitive
    hfactorLeading hsymmetricTarget hrecover
  have hrecCanonical := primitiveRaw_canonical symmetric recovered content hsymmetricCanonical hrecover
  have hrecNonempty : 0 < recovered.size := by
    by_contra hn
    have he : recovered = #[] := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hn)
    simp only [he, SparsePolyZZ.toPoly, Array.toList_empty, List.map_nil, List.sum_nil] at hfactor
    exact hfactorNe hfactor.symm
  have hdivides : SparsePolyZZ.toPoly recovered ∣ SparsePolyZZ.toPoly remaining := by
    rw [hfactor]; exact ⟨complement, hremaining⟩
  obtain ⟨quotient, hdivide⟩ := exactDivmodRaw_complete_of_dvd remaining recovered hcanonical
    hrecCanonical hrecNonempty hdivides
  have hquotientCanonical := exactDivmodRaw_quotient_canonical remaining recovered quotient #[]
    hcanonical.2 hdivide
  have hquotient : SparsePolyZZ.toPoly quotient = complement := by
    have he := exactDivmodRaw_toPoly remaining recovered quotient #[] hdivide
    simp only [SparsePolyZZ.toPoly, Array.toList_empty, List.map_nil, List.sum_nil, add_zero] at he
    change SparsePolyZZ.toPoly remaining = SparsePolyZZ.toPoly recovered * SparsePolyZZ.toPoly quotient at he
    rw [hremaining, hfactor] at he
    exact (mul_left_cancel₀ hfactorNe he).symm
  have hrecEmpty : recovered.isEmpty = false := by simp [Array.isEmpty, Nat.ne_of_gt hrecNonempty]
  have hrecDegree : get_deg recovered = Int64.ofNat factor.natDegree := by
    rw [canonical_get_deg_eq_ofNat recovered hrecCanonical hrecNonempty, hfactor]
  have hdegreeNot : ¬ Int64.ofNat factor.natDegree ≤ 0 := by
    rw [Int64.le_iff_toInt_le, Int64.toInt_ofNat_of_lt hfactorFits]
    change ¬ (factor.natDegree : Int) ≤ 0
    exact not_le_of_gt (by exact_mod_cast hdegree)
  have hremainingDegree : 0 < (SparsePolyZZ.toPoly remaining).natDegree := by
    rw [hremaining, Polynomial.natDegree_mul hfactorNe hcomplementNe]
    omega
  have hremainingNot : ¬ get_deg remaining ≤ 0 := by
    rw [canonical_get_deg_eq_ofNat remaining hcanonical hr, Int64.le_iff_toInt_le,
      Int64.toInt_ofNat_of_lt hfits]
    change ¬ ((SparsePolyZZ.toPoly remaining).natDegree : Int) ≤ 0
    exact not_le_of_gt (by exact_mod_cast hremainingDegree)
  refine ⟨recovered, quotient, hfactor, hquotient, hrecCanonical, hquotientCanonical, ?_⟩
  rw [Generated.StrictVanHoeij.validatePartitionLoop, dif_pos hi, dif_pos hr, if_neg hremainingNot]
  simp only [hproduct, hsymmetric, hrecover, bind, Except.bind]
  simp [hrecEmpty, hrecDegree, hdegreeNot, hdivide, bind, Except.bind]

end Refinement.StrictRecombine
