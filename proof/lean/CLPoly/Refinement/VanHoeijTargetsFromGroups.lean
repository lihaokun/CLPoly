import CLPoly.Refinement.VanHoeijTargetPermutation

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Construct the entire ordered validation target array from the individual
exact-group conclusions. The grouping run is the actual returned grouping;
there is no assumed target array or successful validation execution. -/
theorem exactPartitionTargets_of_individual_groups
    (source : Polynomial Int) (factors : Array (Polynomial Int))
    (hproduct : factors.toList.prod = source) (hnodup : factors.toList.Nodup)
    (hf : source ≠ 0) (hn : 0 < source.natDegree)
    (hprimitive : ∀ i : Fin factors.size, factors[i.val].IsPrimitive)
    (hpositive : ∀ i : Fin factors.size, (0 : Int) < factors[i.val].leadingCoeff)
    (hdegree : ∀ i : Fin factors.size, 0 < factors[i.val].natDegree)
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hm : 0 < modulus)
    (hdivides : prime ∣ modulus) (lifted : Array SparsePolyZZ)
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (groups : Array (Array Int32)) (hfits : lifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform lifted.size = .ok groups)
    (hprime : ∀ i : Fin lifted.size, Prime (StrictHensel.toPolyMod prime lifted[i.val]))
    (hseparate : ∀ i j : Fin lifted.size,
      Associated (StrictHensel.toPolyMod prime lifted[i.val])
        (StrictHensel.toPolyMod prime lifted[j.val]) → i = j)
    (hatom : ∀ i : Fin lifted.size, StrictHensel.toPolyMod prime lifted[i.val] ∣
      source.map (Int.castRingHom (ZMod prime)))
    (hleading : (source.leadingCoeff : ZMod prime) ≠ 0)
    (divisor : Int) (hd : 0 ≤ divisor)
    (hprecision : separationResultantThreshold source lifted.size divisor < (modulus : Int))
    (hexact : ∀ i : Fin factors.size, ∃ (indices : Array Nat) (classId : Nat),
      LegalCombination lifted.size indices.size indices ∧ classId < groups.size ∧
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          (indices.toList.map fun j => StrictHensel.toPolyMod modulus lifted[j]!).prod ∧
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[classId]! lifted 0 ∧
      ∀ index, index < lifted.size →
        (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList)) :
    ∃ targets, ExactPartitionTargets source modulus lifted groups targets := by
  classical
  choose indices classId hlegal hclass hsupport hgroup hmembers using hexact
  let assigned : Fin factors.size → Fin groups.size := fun i => ⟨classId i, hclass i⟩
  have hl : ∀ i : Fin factors.size, (factors[i.val].leadingCoeff : ZMod prime) ≠ 0 := by
    intro i
    apply divisor_leading_mod_nonzero prime source factors[i.val] hleading
    rw [← hproduct]
    exact List.dvd_prod (Array.getElem_mem_toList i.isLt)
  have hsurj := factorGroupAssignment_surjective source factors hproduct prime modulus
    hdivides lifted transform rows groups hfits hrun hprime hseparate hatom hl assigned
    indices hlegal hsupport hmembers
  have hinj := factorGroupAssignment_injective source factors hproduct hnodup hf hn
    hprimitive hpositive modulus hm lifted groups divisor hd hprecision assigned hgroup
  exact exactPartitionTargets_of_bijective_assignment source factors hproduct hprimitive
    hpositive hdegree modulus lifted groups assigned ⟨hinj, hsurj⟩ hgroup

end Refinement.StrictRecombine
