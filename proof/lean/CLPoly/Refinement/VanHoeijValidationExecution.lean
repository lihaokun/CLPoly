import CLPoly.Refinement.VanHoeijPartitionProductExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeij

theorem validatePartitionLoop_succeeds (lifted : Array SparsePolyZZ) (modulus : Nat)
    (hm : 0 < modulus) (groups : Array (Array Int32))
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (index : Nat) (remaining : SparsePolyZZ) (factors : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical remaining) :
    ∃ output, validatePartitionLoop lifted (modulus : Int) groups index remaining factors = .ok output := by
  have hmInt : (0 : Int) < modulus := by exact_mod_cast hm
  induction hmeasure : groups.size - index using Nat.strong_induction_on
      generalizing index remaining factors with
  | h measure ih =>
    rw [validatePartitionLoop]
    split
    next hi =>
      split
      next hr =>
        split
        · exact ⟨none, rfl⟩
        next hdegree =>
          obtain ⟨product, degree, hp⟩ := partitionProductLoop_succeeds groups[index] lifted
            (modulus : Int) 0 #[(⟨0⟩, remaining[0].2)] 0 (ne_of_gt hmInt) (hvalid index hi)
          have hinitial : StrictPolynomialMod.SparsePolyZZCanonical #[(⟨0⟩, remaining[0].2)] := by
            constructor
            · simp
            · intro term ht
              simp at ht
              subst term
              exact hcanonical.2 remaining[0] (Array.getElem_mem_toList hr)
          have hpc := trialProductLoop_canonical ⟨()⟩ groups[index] lifted (modulus : Int) 0
            #[(⟨0⟩, remaining[0].2)] product hinitial
            (partitionProductLoop_projection groups[index] lifted (modulus : Int) 0 _ product 0 degree hp)
          obtain ⟨symmetric, hs⟩ := symmetricModRaw_complete product (modulus : Int) hmInt
          have hsc := symmetricModRaw_canonical product symmetric modulus hm hpc hs
          obtain ⟨content, factor, hf⟩ := primitiveRaw_complete symmetric hsc
          have hfc := primitiveRaw_canonical symmetric factor content hsc hf
          simp only [hp, hs, hf, bind, Except.bind]
          split
          · exact ⟨none, rfl⟩
          next hcheck =>
            have hfactor : 0 < factor.size := by
              by_contra hn
              have he : factor = #[] := Array.size_eq_zero_iff.mp (by omega)
              simp [he] at hcheck
            obtain ⟨quotient, remainder, hdiv⟩ := exactDivmodRaw_complete remaining factor hcanonical hfc hfactor
            have hqc := exactDivmodRaw_quotient_canonical remaining factor quotient remainder hcanonical.2 hdiv
            simp only [hdiv, bind, Except.bind]
            split
            · exact ih _ (by omega) (index + 1) quotient (factors.push factor) hqc rfl
            · exact ⟨none, rfl⟩
      next hr => exact ⟨none, rfl⟩
    next hi => split <;> exact ⟨_, rfl⟩

theorem validatePartition_succeeds (f : SparsePolyZZ) (lifted : Array SparsePolyZZ)
    (modulus : Nat) (hm : 0 < modulus) (groups : Array (Array Int32))
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (answer : VanHoeijResult) :
    ∃ output, validatePartition f lifted (modulus : Int) groups answer = .ok output := by
  unfold validatePartition
  split
  · exact ⟨_, rfl⟩
  · obtain ⟨output, ho⟩ := validatePartitionLoop_succeeds lifted modulus hm groups hvalid 0 f #[] hcanonical
    simp only [ho, bind, Except.bind]
    cases output <;> exact ⟨_, rfl⟩

end Refinement.StrictRecombine
