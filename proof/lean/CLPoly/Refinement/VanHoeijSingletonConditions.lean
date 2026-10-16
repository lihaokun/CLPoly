import CLPoly.Refinement.VanHoeijSingletonConstruction
import CLPoly.Refinement.Recombine

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem singletonGroups_singleton (count : Int32) (output : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok output) :
    ∀ group ∈ output.toList, group.size = 1 := by
  intro group hg
  obtain ⟨i, hi, he⟩ := List.mem_iff_getElem.mp hg
  have hi' : i < output.size := by simpa using hi
  have he' : output[i] = group := by simpa only [Array.getElem_toList] using he
  rw [← he', singletonGroups_entries count output hrun i hi']
  rfl

theorem singletonGroups_indices_valid (count : Int32) (output : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok output) :
    ∀ i (hi : i < output.size),
      CandidateIndicesValid output[i] (Array.replicate count.toInt.toNat false) := by
  intro i hi
  have hs := singletonGroups_size count output hrun
  have hib : i < count.toInt.toNat := by omega
  have hc := Int32.toInt_lt count
  have hfits : i < 2 ^ 31 := by omega
  have hr := nat_toUInt32_toInt32_nonnegative_and_toNat i hfits
  rw [singletonGroups_entries count output hrun i hi]
  intro position hp
  have hz : position = 0 := by simpa using hp
  subst position
  change 0 ≤ Int32.ofNat i ∧ (Int32.ofNat i).toInt64.toNat < _
  change 0 ≤ i.toUInt32.toInt32 ∧ i.toUInt32.toInt32.toInt64.toNat < _
  refine ⟨hr.1, ?_⟩
  simpa only [hr.2, Array.size_replicate] using hib

end Refinement.StrictRecombine
