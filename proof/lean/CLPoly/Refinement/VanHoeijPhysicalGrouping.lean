import CLPoly.Refinement.VanHoeijPhysicalSupport
import CLPoly.Refinement.VanHoeijGroupingExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem candidatePhysicalSupport_uniform
    (count : Nat) (groups : Array (Array Int32)) (classId : Nat)
    (hc : classId < groups.size)
    (hvalid : CandidateIndicesValid groups[classId] (Array.replicate count false))
    (hfits : count ≤ 2 ^ 31) (indices : Array Nat)
    (huniform : ∀ leftValue rightValue,
      CandidateMember groups classId leftValue →
      CandidateMember groups classId rightValue →
      ∃ left right : Nat, left < count ∧ right < count ∧
        leftValue = left.toUInt32.toInt32 ∧ rightValue = right.toUInt32.toInt32 ∧
        (left ∈ indices.toList ↔ right ∈ indices.toList)) :
    ∀ left ∈ candidatePhysicalSupport count groups[classId] hvalid,
      ∀ right ∈ candidatePhysicalSupport count groups[classId] hvalid,
        (left.val ∈ indices.toList ↔ right.val ∈ indices.toList) := by
  intro left hl right hr
  obtain ⟨lp, hlp⟩ := List.mem_ofFn.mp hl
  obtain ⟨rp, hrp⟩ := List.mem_ofFn.mp hr
  have hlm : CandidateMember groups classId groups[classId][lp.val] := by
    refine ⟨hc, lp.val, ?_, ?_⟩
    · simpa only [getElem!_pos groups classId hc] using lp.isLt
    · simp only [getElem!_pos groups classId hc,
        getElem!_pos groups[classId] lp.val lp.isLt]
  have hrm : CandidateMember groups classId groups[classId][rp.val] := by
    refine ⟨hc, rp.val, ?_, ?_⟩
    · simpa only [getElem!_pos groups classId hc] using rp.isLt
    · simp only [getElem!_pos groups classId hc,
        getElem!_pos groups[classId] rp.val rp.isLt]
  obtain ⟨l, r, hlt, hrt, hlv, hrv, heq⟩ := huniform _ _ hlm hrm
  have hleft : groups[classId][lp.val].toInt64.toNat = left.val := congrArg Fin.val hlp
  have hright : groups[classId][rp.val].toInt64.toNat = right.val := congrArg Fin.val hrp
  rw [hlv, (nat_toUInt32_toInt32_nonnegative_and_toNat l
    (lt_of_lt_of_le hlt hfits)).2] at hleft
  rw [hrv, (nat_toUInt32_toInt32_nonnegative_and_toNat r
    (lt_of_lt_of_le hrt hfits)).2] at hright
  simpa only [hleft, hright] using heq

end Refinement.StrictRecombine
