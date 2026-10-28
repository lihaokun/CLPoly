import CLPoly.Refinement.VanHoeijValidationSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine

def ValidatedNonunit (factor : SparsePolyZZ) : Prop :=
  ¬ IsUnit (SparsePolyZZ.toPoly factor)

theorem validatePartitionLoop_nonunit
    (lifted : Array SparsePolyZZ) (modulus prime : Nat) [Fact (Nat.Prime prime)]
    (groups : Array (Array Int32)) (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hirreducible : ∀ i (hi : i < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (index : Nat) (remaining : SparsePolyZZ) (factors result : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hleading : ((SparsePolyZZ.toPoly remaining).leadingCoeff : ZMod prime) ≠ 0)
    (hfactors : ∀ factor ∈ factors.toList, ValidatedNonunit factor)
    (hrun : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups index
      remaining factors = .ok (some result)) :
    ∀ factor ∈ result.toList, ValidatedNonunit factor := by
  induction hmeasure : groups.size - index using Nat.strong_induction_on
      generalizing index remaining factors result with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.validatePartitionLoop] at hrun
    split at hrun
    next hi =>
      split at hrun
      next hr =>
        split at hrun
        next hconstant => simp at hrun
        next hpositive =>
          cases hp : Generated.StrictVanHoeij.partitionProductLoop groups[index] lifted modulus
              0 #[(⟨0⟩, remaining[0].2)] 0 with
          | error fault => simp [hp, bind, Except.bind] at hrun
          | ok pair =>
            obtain ⟨product, degree⟩ := pair
            simp only [hp, bind, Except.bind] at hrun
            cases hs : Generated.StrictRecombine.symmetricModRaw product modulus with
            | error fault => simp [hs] at hrun
            | ok symmetric =>
              simp only [hs] at hrun
              cases ht : Generated.StrictRecombine.primitiveRaw symmetric with
              | error fault => simp [ht] at hrun
              | ok primitive =>
                obtain ⟨content, factor⟩ := primitive
                simp only [ht] at hrun
                split at hrun
                next hbad => simp at hrun
                next hgood =>
                  cases hdiv : Generated.StrictRecombine.exactDivmodRaw remaining factor with
                  | error fault => simp [hdiv] at hrun
                  | ok divided =>
                    obtain ⟨quotient, remainder⟩ := divided
                    simp only [hdiv] at hrun
                    split at hrun
                    next hempty =>
                      have hrzero : remainder = #[] := Array.isEmpty_iff.mp hempty
                      have hdivide : Generated.StrictRecombine.exactDivmodRaw remaining factor =
                          .ok (quotient, #[]) := by simpa only [hrzero] using hdiv
                      have hqcanonical := exactDivmodRaw_quotient_canonical remaining factor
                        quotient remainder hcanonical.2 hdiv
                      have hlead := partition_exact_division_leading_nonzero prime remaining
                        factor quotient hleading hdivide
                      have hhead : (remaining[0].2 : ZMod prime) ≠ 0 := by
                        rw [← sparsePolyZZ_leadingCoeff_eq_head remaining hcanonical hr]
                        exact hleading
                      have hgroup : 0 < groups[index].size := by
                        by_contra hn
                        have hz : groups[index].size = 0 := by omega
                        have hnil : groups[index] = #[] :=
                          Array.isEmpty_iff.mp (by simp [Array.isEmpty, hz])
                        have hpzero := hp
                        simp only [hnil, Generated.StrictVanHoeij.partitionProductLoop,
                          Array.size_empty, Nat.not_lt_zero, ↓reduceDIte] at hpzero
                        have hdeg : (0 : Int64) = degree :=
                          congrArg Prod.snd (Except.ok.inj hpzero)
                        simp [← hdeg] at hgood
                      have hnew : ValidatedNonunit factor :=
                        validationRecoveredFactor_not_isUnit ⟨⟨()⟩⟩ groups[index] lifted
                          modulus prime remaining product symmetric factor content hm hdivides
                          (hvalid index hi) hgroup hirreducible hr hhead
                          (partitionProductLoop_projection groups[index] lifted modulus 0 _
                            product 0 degree hp) hs ht
                      have hfactorsNext : ∀ q ∈ (factors.push factor).toList,
                          ValidatedNonunit q := by
                        intro q hq
                        simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hq
                        rcases hq with hq | hq
                        · exact hfactors q hq
                        · simpa only [hq] using hnew
                      exact ih (groups.size - (index + 1)) (by omega)
                        (index + 1) quotient (factors.push factor) result hqcanonical hlead.2
                        hfactorsNext hrun rfl
                    next hnotempty => simp at hrun
      next hempty => simp at hrun
    next hdone =>
      split at hrun
      next hbad => simp at hrun
      next hgood =>
        have hout := Option.some.inj (Except.ok.inj hrun)
        intro factor hfactor
        rw [← hout] at hfactor
        have hperm := List.mergeSort_perm factors.toList fun left right =>
          left[0]!.1.deg < right[0]!.1.deg
        apply hfactors factor
        exact hperm.mem_iff.mp (by simpa [Generated.StrictRecombine.sortFactorsByDegree] using hfactor)

theorem validatePartition_nonunit
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus prime : Nat)
    [Fact (Nat.Prime prime)] (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hirreducible : ∀ i (hi : i < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hleading : ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod prime) ≠ 0)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, ValidatedNonunit factor := by
  unfold Generated.StrictVanHoeij.validatePartition at hrun
  split at hrun
  · simp at hrun
  · cases hloop : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups 0 f #[] with
    | error fault => simp [hloop, bind, Except.bind] at hrun
    | ok result =>
      cases result with
      | none => simp [hloop, bind, Except.bind] at hrun
      | some factors =>
        simp only [hloop, bind, Except.bind] at hrun
        have heq := (Prod.mk.inj (Except.ok.inj hrun)).2
        rw [← heq]
        exact validatePartitionLoop_nonunit lifted modulus prime groups hm hdivides hvalid
          hirreducible 0 f #[] factors hcanonical hleading (by simp) hloop

theorem primitive_nonunit_natDegree_pos (factor : Polynomial Int)
    (hprimitive : factor.IsPrimitive) (hnonunit : ¬ IsUnit factor) :
    0 < factor.natDegree := by
  by_contra hn
  have hz : factor.natDegree = 0 := Nat.eq_zero_of_not_pos hn
  have he := Polynomial.eq_C_of_natDegree_eq_zero hz
  have hd : Polynomial.C (factor.coeff 0) ∣ factor := ⟨1, by simpa using he⟩
  have hu := Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hprimitive (factor.coeff 0) hd
  apply hnonunit
  rw [he]
  exact Polynomial.isUnit_C.mpr hu

theorem divisor_leading_mod_nonzero (prime : Nat)
    (source divisor : Polynomial Int)
    (hleading : (source.leadingCoeff : ZMod prime) ≠ 0)
    (hdivides : divisor ∣ source) :
    (divisor.leadingCoeff : ZMod prime) ≠ 0 := by
  obtain ⟨quotient, hq⟩ := hdivides
  have hp : (divisor.leadingCoeff : ZMod prime) *
      (quotient.leadingCoeff : ZMod prime) ≠ 0 := by
    simpa only [hq, Polynomial.leadingCoeff_mul, Int.cast_mul] using hleading
  exact left_ne_zero_of_mul hp

end Refinement.StrictRecombine
