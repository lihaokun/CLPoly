import CLPoly.Refinement.VanHoeijPartitionRecoveryStep
import CLPoly.Refinement.VanHoeijValidationNonunit

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem canonical_one_passes_validation_checks
    (remaining : SparsePolyZZ) (hc : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hone : SparsePolyZZ.toPoly remaining = 1) :
    remaining.size = 1 ∧ get_deg remaining = 0 ∧ remaining[0]!.2 = 1 := by
  have hn : 0 < remaining.size := by
    by_contra hnot
    have he : remaining = #[] := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hnot)
    simp [he, SparsePolyZZ.toPoly] at hone
  have hheadDegree : remaining[0]!.1.deg = 0 := by
    have h := sparsePolyZZ_natDegree_eq_head remaining hc hn
    rw [hone] at h
    simpa only [Polynomial.natDegree_one, getElem!_pos remaining 0 hn] using h.symm
  have htail : remaining.toList.drop 1 = [] := by
    apply List.eq_nil_iff_forall_not_mem.mpr
    intro term ht
    have h := StrictHensel.canonical_tail_degree_lt_head remaining hn hc term ht
    rw [hheadDegree] at h
    omega
  have hsize : remaining.size = 1 := by
    have h := List.drop_eq_nil_iff.mp htail
    simp only [Array.length_toList] at h
    omega
  refine ⟨hsize, ?_, ?_⟩
  · rw [canonical_get_deg_eq_ofNat remaining hc hn, hone]
    simp
  · have h := sparsePolyZZ_leadingCoeff_eq_head remaining hc hn
    rw [hone] at h
    simpa only [Polynomial.leadingCoeff_one, getElem!_pos remaining 0 hn] using h.symm

theorem positive_leading_list_product (factors : List (Polynomial Int))
    (hpositive : ∀ factor ∈ factors, (0 : Int) < factor.leadingCoeff) :
    (0 : Int) < factors.prod.leadingCoeff := by
  induction factors with
  | nil => simp
  | cons head tail ih =>
    rw [List.prod_cons, Polynomial.leadingCoeff_mul]
    exact mul_pos (hpositive head List.mem_cons_self)
      (ih (fun factor hf => hpositive factor (List.mem_cons_of_mem head hf)))

theorem list_suffix_product_dvd (factors : List (Polynomial Int)) (index : Nat) :
    (factors.drop index).prod ∣ factors.prod := by
  exact ⟨(factors.take index).prod,
    (List.prod_take_mul_prod_drop factors index).symm.trans (mul_comm _ _)⟩

/-- Mathematical targets aligned with the actual physical groups. There is
no irreducibility field and no assumed successful validation execution.
Deriving this entire array from actual grouping remains a separate obligation. -/
structure ExactPartitionTargets (source : Polynomial Int) (modulus : Nat)
    (lifted : Array SparsePolyZZ) (groups : Array (Array Int32))
    (targets : Array (Polynomial Int)) : Prop where
  size_eq : targets.size = groups.size
  product_eq : targets.toList.prod = source
  primitive : ∀ factor ∈ targets.toList, factor.IsPrimitive
  leading_pos : ∀ factor ∈ targets.toList, (0 : Int) < factor.leadingCoeff
  degree_pos : ∀ factor ∈ targets.toList, 0 < factor.natDegree
  support : ∀ index, index < groups.size →
    targets[index]!.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (targets[index]!.leadingCoeff : ZMod modulus) *
        SelectedProductMod modulus groups[index]! lifted 0

theorem validatePartitionLoop_complete_of_targets
    (source : Polynomial Int) (prime modulus : Nat) [Fact (Nat.Prime prime)]
    (hm : 0 < modulus) (hd : prime ∣ modulus)
    (lifted : Array SparsePolyZZ) (groups : Array (Array Int32))
    (targets : Array (Polynomial Int)) (data : ExactPartitionTargets source modulus lifted groups targets)
    (hsourceDegree : 0 < source.natDegree) (hsourceFits : source.natDegree < 2 ^ 63)
    (hleading : (source.leadingCoeff : ZMod prime) ≠ 0)
    (hvalid : ∀ i (hi : i < groups.size), CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hliftedCanonical : ∀ i, i < lifted.size → StrictPolynomialMod.SparsePolyZZCanonical lifted[i]!)
    (hliftedNonempty : ∀ i, i < lifted.size → 0 < lifted[i]!.size)
    (hmonic : ∀ i, i < lifted.size → (SparsePolyZZ.toPoly lifted[i]!).Monic)
    (divisor : Int) (hdivisor : 0 ≤ divisor)
    (hthreshold : separationResultantThreshold source lifted.size divisor < (modulus : Int))
    (index : Nat) (remaining : SparsePolyZZ) (previous : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hremaining : SparsePolyZZ.toPoly remaining = (targets.toList.drop index).prod) :
    ∃ result, Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups index remaining previous =
      .ok (some result) := by
  have hsourceNe : source ≠ 0 := by intro hz; simp [hz] at hleading
  induction hmeasure : groups.size - index using Nat.strong_induction_on
      generalizing index remaining previous with
  | h measure ih =>
    by_cases hi : index < groups.size
    · have ht : index < targets.size := by rw [data.size_eq]; exact hi
      have hmem : targets[index] ∈ targets.toList := Array.getElem_mem_toList ht
      have hfactorDvd : targets[index] ∣ source := by
        rw [← data.product_eq]; exact List.dvd_prod hmem
      have hcomplementDvd : (targets.toList.drop (index + 1)).prod ∣ source := by
        rw [← data.product_eq]; exact list_suffix_product_dvd _ _
      have hremainingDvd : SparsePolyZZ.toPoly remaining ∣ source := by
        rw [hremaining, ← data.product_eq]; exact list_suffix_product_dvd _ _
      have hremainingFits : (SparsePolyZZ.toPoly remaining).natDegree < 2 ^ 63 :=
        lt_of_le_of_lt (Polynomial.natDegree_le_of_dvd hremainingDvd hsourceNe) hsourceFits
      have hsplit : SparsePolyZZ.toPoly remaining = targets[index] *
          (targets.toList.drop (index + 1)).prod := by
        have hdrop : targets.toList.drop index = targets[index] :: targets.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons (l := targets.toList) (i := index) (by simpa using ht)
        rw [hremaining, hdrop, List.prod_cons]
      have hsupport : targets[index].map (Int.castRingHom (ZMod modulus)) =
          Polynomial.C (targets[index].leadingCoeff : ZMod modulus) *
            SelectedProductMod modulus groups[index] lifted 0 := by
        simpa only [getElem!_pos targets index ht, getElem!_pos groups index hi] using data.support index hi
      obtain ⟨recovered, quotient, hfactor, hquotient, _, hquotientCanonical, hstep⟩ :=
        validatePartitionLoop_step_succeeds prime modulus hm hd lifted groups index hi remaining previous
          hcanonical hremainingFits targets[index] (targets.toList.drop (index + 1)).prod hsplit
          (data.primitive _ hmem) (data.leading_pos _ hmem)
          (positive_leading_list_product _ (fun factor hf =>
            data.leading_pos factor (List.mem_of_mem_drop hf))) (data.degree_pos _ hmem)
          (divisor_leading_mod_nonzero prime source targets[index] hleading hfactorDvd)
          (hvalid index hi) hliftedCanonical hliftedNonempty hmonic hsupport
          (separation_threshold_covers_scaled_recovery source targets[index]
            (targets.toList.drop (index + 1)).prod lifted.size divisor hsourceNe hsourceDegree
            hdivisor hfactorDvd hcomplementDvd modulus hthreshold)
      obtain ⟨result, hr⟩ := ih (groups.size - (index + 1)) (by omega) (index + 1)
        quotient (previous.push recovered) hquotientCanonical hquotient rfl
      exact ⟨result, hstep.trans hr⟩
    · have hdrop : targets.toList.drop index = [] :=
        List.drop_eq_nil_iff.mpr (by simpa only [Array.length_toList, data.size_eq] using Nat.le_of_not_gt hi)
      have hone : SparsePolyZZ.toPoly remaining = 1 := by simpa only [hdrop, List.prod_nil] using hremaining
      obtain ⟨hsize, hdegree, hcoeff⟩ := canonical_one_passes_validation_checks remaining hcanonical hone
      have hcoeff' : remaining[0].2 = (1 : ZZ) := by
        simpa only [getElem!_pos remaining 0 (by omega)] using hcoeff
      refine ⟨Generated.StrictRecombine.sortFactorsByDegree previous, ?_⟩
      rw [Generated.StrictVanHoeij.validatePartitionLoop, dif_neg hi]
      simp [hsize, hdegree, hcoeff']

theorem validatePartition_complete_of_targets
    (f : SparsePolyZZ) (prime modulus : Nat) [Fact (Nat.Prime prime)]
    (hm : 0 < modulus) (hd : prime ∣ modulus)
    (lifted : Array SparsePolyZZ) (groups : Array (Array Int32))
    (targets : Array (Polynomial Int))
    (data : ExactPartitionTargets (SparsePolyZZ.toPoly f) modulus lifted groups targets)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hleading : ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod prime) ≠ 0)
    (hvalid : ∀ i (hi : i < groups.size), CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hliftedCanonical : ∀ i, i < lifted.size → StrictPolynomialMod.SparsePolyZZCanonical lifted[i]!)
    (hliftedNonempty : ∀ i, i < lifted.size → 0 < lifted[i]!.size)
    (hmonic : ∀ i, i < lifted.size → (SparsePolyZZ.toPoly lifted[i]!).Monic)
    (divisor : Int) (hdivisor : 0 ≤ divisor)
    (hthreshold : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.size divisor < (modulus : Int))
    (before : Generated.StrictVanHoeij.VanHoeijResult) :
    ∃ after, Generated.StrictVanHoeij.validatePartition f lifted modulus groups before = .ok (true, after) := by
  have hgroups : groups.isEmpty = false := by
    cases he : groups.isEmpty with
    | false => rfl
    | true =>
      have hsize : groups.size = 0 := by rw [Array.isEmpty_iff.mp he]; rfl
      have ht : targets = #[] := Array.size_eq_zero_iff.mp (data.size_eq.trans hsize)
      have hp := data.product_eq
      simp only [ht, Array.toList_empty, List.prod_nil] at hp
      rw [← hp] at hdegree
      simp at hdegree
  obtain ⟨result, hrun⟩ := validatePartitionLoop_complete_of_targets (SparsePolyZZ.toPoly f)
    prime modulus hm hd lifted groups targets data hdegree hfits hleading hvalid
    hliftedCanonical hliftedNonempty hmonic divisor hdivisor hthreshold 0 f #[] hcanonical
    (by simpa only [List.drop_zero] using data.product_eq.symm)
  refine ⟨{ before with factors := result, complete := true }, ?_⟩
  simp [Generated.StrictVanHoeij.validatePartition, hgroups, hrun, bind, Except.bind]

end Refinement.StrictRecombine
