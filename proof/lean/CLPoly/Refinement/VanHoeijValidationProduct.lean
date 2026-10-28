import CLPoly.Refinement.VanHoeijValidation

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem validated_final_quotient_one (remaining : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hdegree : (SparsePolyZZ.toPoly remaining).natDegree < 2 ^ 63)
    (hsize : remaining.size = 1) (hget : get_deg remaining = 0)
    (hcoeff : remaining[0]!.2 = 1) : SparsePolyZZ.toPoly remaining = 1 := by
  have hn : 0 < remaining.size := by omega
  have hdeg := sparsePolyZZ_natDegree_eq_head remaining hcanonical hn
  have hempty : remaining.isEmpty = false := by simp [Array.isEmpty, hsize]
  have hcast : get_deg remaining = Int64.ofNat (SparsePolyZZ.toPoly remaining).natDegree := by
    unfold get_deg
    rw [hempty]
    simp only [Bool.false_eq_true, ↓reduceIte, getElem!_pos remaining 0 hn]
    rw [hdeg]
    rfl
  have hzero : (SparsePolyZZ.toPoly remaining).natDegree = 0 := by
    have hz := congrArg Int64.toInt hget
    rw [hcast, Int64.toInt_ofNat_of_lt hdegree] at hz
    simpa using hz
  have hhead : remaining[0]! = (⟨0⟩, (1 : ZZ)) := by
    apply Prod.ext
    · have hd : remaining[0]!.1.deg = 0 := by
        simpa only [getElem!_pos remaining 0 hn] using hdeg.symm.trans hzero
      cases hm : remaining[0]!.1 with
      | mk d =>
        have hd' : d = 0 := by simpa only [hm] using hd
        subst d
        rfl
    · exact hcoeff
  have harray : remaining = #[remaining[0]!] := by
    apply Array.ext (by simpa using hsize)
    intro index hi hj
    have hz : index = 0 := by omega
    subst index
    simp [getElem!_pos remaining 0 hn]
  rw [harray, hhead]
  simp [SparsePolyZZ.toPoly]

theorem partition_sort_product (factors : Array SparsePolyZZ) :
    factorArrayProduct (Generated.StrictRecombine.sortFactorsByDegree factors) =
      factorArrayProduct factors := by
  unfold factorArrayProduct Generated.StrictRecombine.sortFactorsByDegree
  simpa using ((List.mergeSort_perm factors.toList fun left right =>
    left[0]!.1.deg < right[0]!.1.deg).map SparsePolyZZ.toPoly).prod_eq

theorem validatePartitionLoop_product
    (lifted : Array SparsePolyZZ) (modulus : ZZ) (groups : Array (Array Int32))
    (index : Nat) (remaining : SparsePolyZZ) (factors result : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hdegree : (SparsePolyZZ.toPoly remaining).natDegree < 2 ^ 63)
    (hrun : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups index
      remaining factors = .ok (some result)) :
    factorArrayProduct result = factorArrayProduct factors * SparsePolyZZ.toPoly remaining := by
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
                      have hdivision := exactDivmodRaw_toPoly remaining factor quotient remainder hdiv
                      rw [hrzero] at hdivision
                      simp only [SparsePolyZZ.toPoly, Array.toList_empty, List.map_nil,
                        List.sum_nil, add_zero] at hdivision
                      change SparsePolyZZ.toPoly remaining =
                        SparsePolyZZ.toPoly factor * SparsePolyZZ.toPoly quotient at hdivision
                      have hqcanonical := exactDivmodRaw_quotient_canonical remaining factor
                        quotient remainder hcanonical.2 hdiv
                      have hremainNe := sparsePolyZZ_toPoly_ne_zero_of_canonical_nonempty
                        remaining hcanonical hr
                      have hqdvd : SparsePolyZZ.toPoly quotient ∣ SparsePolyZZ.toPoly remaining :=
                        ⟨SparsePolyZZ.toPoly factor, by rw [hdivision, mul_comm]⟩
                      have hqdegree : (SparsePolyZZ.toPoly quotient).natDegree < 2 ^ 63 :=
                        lt_of_le_of_lt (Polynomial.natDegree_le_of_dvd hqdvd hremainNe) hdegree
                      have hnext := ih (groups.size - (index + 1)) (by omega)
                        (index + 1) quotient (factors.push factor) result hqcanonical hqdegree hrun rfl
                      rw [hnext, hdivision]
                      simp only [factorArrayProduct, Array.toList_push, List.map_append,
                        List.map_cons, List.map_nil, List.prod_append, List.prod_cons,
                        List.prod_nil, mul_one, mul_assoc]
                    next hnotempty => simp at hrun
      next hempty => simp at hrun
    next hdone =>
      split at hrun
      next hbad => simp at hrun
      next hgood =>
        have hchecks : remaining.size = 1 ∧ get_deg remaining = 0 ∧ remaining[0]!.2 = 1 := by
          simpa only [Bool.or_eq_true, bne_iff_ne, not_or, not_not, and_assoc] using hgood
        have hone := validated_final_quotient_one remaining hcanonical hdegree
          hchecks.1 hchecks.2.1 hchecks.2.2
        have hout := Option.some.inj (Except.ok.inj hrun)
        rw [← hout, partition_sort_product, hone, mul_one]

theorem validatePartition_product
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    factorArrayProduct after.factors = SparsePolyZZ.toPoly f := by
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
        have hproduct := validatePartitionLoop_product lifted modulus groups 0 f #[] factors
          hcanonical hdegree hloop
        simpa [factorArrayProduct] using hproduct

end Refinement.StrictRecombine
