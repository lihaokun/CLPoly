import CLPoly.Refinement.VanHoeijValidationSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine

noncomputable def HasPartitionSupport (prime : Nat) (lifted : Array SparsePolyZZ)
    (groups : Array (Array Int32)) (factor : SparsePolyZZ) : Prop :=
  ∃ group ∈ groups.toList, Associated (StrictHensel.toPolyMod prime factor)
    (SelectedProductMod prime group lifted 0)

theorem validatePartitionLoop_support
    (lifted : Array SparsePolyZZ) (modulus prime : Nat) [Fact (Nat.Prime prime)]
    (groups : Array (Array Int32)) (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hirreducible : ∀ i (hi : i < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (index : Nat) (remaining : SparsePolyZZ) (factors result : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining)
    (hleading : ((SparsePolyZZ.toPoly remaining).leadingCoeff : ZMod prime) ≠ 0)
    (hfactors : ∀ factor ∈ factors.toList, HasPartitionSupport prime lifted groups factor)
    (hrun : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups index
      remaining factors = .ok (some result)) :
    ∀ factor ∈ result.toList, HasPartitionSupport prime lifted groups factor := by
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
                      have hassoc := partition_recovered_factor_associated groups[index] lifted
                        modulus prime remaining product symmetric factor content degree hm hdivides
                        (hvalid index hi) hirreducible hr hhead hp hs ht
                      have hnew : HasPartitionSupport prime lifted groups factor :=
                        ⟨groups[index], by simp, hassoc⟩
                      have hfactorsNext : ∀ q ∈ (factors.push factor).toList,
                          HasPartitionSupport prime lifted groups q := by
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

theorem validatePartition_support
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
    ∀ factor ∈ after.factors.toList, HasPartitionSupport prime lifted groups factor := by
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
        exact validatePartitionLoop_support lifted modulus prime groups hm hdivides hvalid
          hirreducible 0 f #[] factors hcanonical hleading (by simp) hloop

end Refinement.StrictRecombine
