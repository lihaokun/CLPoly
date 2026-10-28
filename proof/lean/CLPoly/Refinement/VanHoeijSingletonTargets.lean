import CLPoly.Refinement.VanHoeijTargetPermutation
import CLPoly.Refinement.VanHoeijSingletonConditions

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem singletonGroups_selectedProductMod (count : Int32) (groups : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok groups)
    (lifted : Array SparsePolyZZ) (modulus index : Nat) (hi : index < groups.size) :
    SelectedProductMod modulus groups[index]! lifted 0 =
      StrictHensel.toPolyMod modulus lifted[index]! := by
  have hs := singletonGroups_size count groups hrun
  have hc := Int32.toInt_lt count
  have hfits : index < 2 ^ 31 := by omega
  have hdecode : (Int32.ofNat index).toInt64.toNat = index :=
    (nat_toUInt32_toInt32_nonnegative_and_toNat index hfits).2
  rw [getElem!_pos groups index hi, singletonGroups_entries count groups hrun index hi]
  simp [SelectedProductMod, hdecode]

/-- Build all validation targets for the actual singleton-group execution
from per-factor singleton modular supports, using prime-atom coverage. -/
theorem exactPartitionTargets_of_singleton_supports
    (source : Polynomial Int) (factors : Array (Polynomial Int))
    (hproduct : factors.toList.prod = source) (hnodup : factors.toList.Nodup)
    (hf : source ≠ 0) (hn : 0 < source.natDegree)
    (hprimitive : ∀ i : Fin factors.size, factors[i.val].IsPrimitive)
    (hpositive : ∀ i : Fin factors.size, (0 : Int) < factors[i.val].leadingCoeff)
    (hdegree : ∀ i : Fin factors.size, 0 < factors[i.val].natDegree)
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hm : 0 < modulus)
    (hdivides : prime ∣ modulus) (lifted : Array SparsePolyZZ)
    (count : Int32) (groups : Array (Array Int32))
    (hcount : count.toInt.toNat = lifted.size)
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok groups)
    (hprime : ∀ i : Fin lifted.size, Prime (StrictHensel.toPolyMod prime lifted[i.val]))
    (hseparate : ∀ i j : Fin lifted.size,
      Associated (StrictHensel.toPolyMod prime lifted[i.val])
        (StrictHensel.toPolyMod prime lifted[j.val]) → i = j)
    (hatom : ∀ i : Fin lifted.size, StrictHensel.toPolyMod prime lifted[i.val] ∣
      source.map (Int.castRingHom (ZMod prime)))
    (hleading : (source.leadingCoeff : ZMod prime) ≠ 0)
    (divisor : Int) (hd : 0 ≤ divisor)
    (hprecision : separationResultantThreshold source lifted.size divisor < (modulus : Int))
    (hexact : ∀ i : Fin factors.size, ∃ index : Fin lifted.size,
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          StrictHensel.toPolyMod modulus lifted[index.val]!) :
    ∃ targets, ExactPartitionTargets source modulus lifted groups targets := by
  classical
  choose indices hs using hexact
  have hsize : groups.size = lifted.size := (singletonGroups_size count groups hrun).trans hcount
  let assigned : Fin factors.size → Fin groups.size := fun i =>
    ⟨(indices i).val, by rw [hsize]; exact (indices i).isLt⟩
  have hsupport : ∀ i, factors[i.val].map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
        SelectedProductMod modulus groups[(assigned i).val]! lifted 0 := by
    intro i
    rw [singletonGroups_selectedProductMod count groups hrun lifted modulus _ (assigned i).isLt]
    exact hs i
  have hinj := factorGroupAssignment_injective source factors hproduct hnodup hf hn hprimitive
    hpositive modulus hm lifted groups divisor hd hprecision assigned hsupport
  have hsurj : Function.Surjective assigned := by
    intro classId
    let i : Fin lifted.size := ⟨classId.val, by rw [← hsize]; exact classId.isLt⟩
    let atoms := fun j : Fin lifted.size => StrictHensel.toPolyMod prime lifted[j.val]
    have hdiv : atoms i ∣ (factors.toList.map (Polynomial.map (Int.castRingHom (ZMod prime)))).prod := by
      rw [← Polynomial.map_list_prod, hproduct]
      exact hatom i
    obtain ⟨reduced, hmem, hdvd⟩ := prime_dvd_list_product_witness (atoms i) (hprime i) _ hdiv
    obtain ⟨factor, hfactor, he⟩ := List.mem_map.mp hmem
    obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hfactor
    have hp' : position < factors.size := by simpa using hp
    have hv' : factors[position] = factor := by simpa only [Array.getElem_toList hp'] using hv
    let j : Fin factors.size := ⟨position, hp'⟩
    have hfactorDvd : factors[j.val] ∣ source := by
      rw [← hproduct]; exact List.dvd_prod (Array.getElem_mem_toList j.isLt)
    have hl := divisor_leading_mod_nonzero prime source factors[j.val] hleading hfactorDvd
    have hpSupport := selected_support_project prime modulus hdivides
      groups[(assigned j).val]! lifted factors[j.val] (hsupport j)
    rw [singletonGroups_selectedProductMod count groups hrun lifted prime _ (assigned j).isLt] at hpSupport
    have ha : Associated (factors[j.val].map (Int.castRingHom (ZMod prime))) (atoms (indices j)) := by
      rw [hpSupport]
      change Associated (Polynomial.C (factors[j.val].leadingCoeff : ZMod prime) *
        StrictHensel.toPolyMod prime lifted[(indices j).val]!) (atoms (indices j))
      rw [getElem!_pos lifted (indices j).val (indices j).isLt]
      exact (associated_isUnit_mul_left_iff
        (Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hl))).mpr (Associated.refl _)
    have hdvd' : atoms i ∣ factors[j.val].map (Int.castRingHom (ZMod prime)) := by
      simpa only [j, hv', he] using hdvd
    have heq : i = indices j := hseparate i (indices j)
      ((hprime i).associated_of_dvd (hprime (indices j)) (hdvd'.trans ha.dvd))
    refine ⟨j, Fin.ext ?_⟩
    exact (congrArg Fin.val heq).symm
  exact exactPartitionTargets_of_bijective_assignment source factors hproduct hprimitive
    hpositive hdegree modulus lifted groups assigned ⟨hinj, hsurj⟩ hsupport

end Refinement.StrictRecombine
