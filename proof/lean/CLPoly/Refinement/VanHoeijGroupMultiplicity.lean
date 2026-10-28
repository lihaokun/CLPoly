import CLPoly.Refinement.VanHoeijExactGrouping

set_option autoImplicit false

namespace Refinement.StrictRecombine

def classSourceSuffix (classes : Array (Option Nat)) (classId start : Nat) : List Nat :=
  ((List.range classes.size).drop start).filter fun index => classes[index]! = some classId

theorem classSourceSuffix_step (classes : Array (Option Nat)) (classId start : Nat)
    (hs : start < classes.size) :
    classSourceSuffix classes classId start =
      if classes[start]! = some classId then start :: classSourceSuffix classes classId (start + 1)
      else classSourceSuffix classes classId (start + 1) := by
  unfold classSourceSuffix
  rw [List.drop_eq_getElem_cons (by simpa using hs : start < (List.range classes.size).length)]
  simp only [List.getElem_range, List.filter_cons, decide_eq_true_eq]

theorem classSourceSuffix_nodup (classes : Array (Option Nat)) (classId start : Nat) :
    (classSourceSuffix classes classId start).Nodup := by
  have hd : ((List.range classes.size).drop start).Nodup := List.nodup_range.drop
  exact hd.filter _

theorem classSourceSuffix_bound (classes : Array (Option Nat)) (classId start index : Nat)
    (hi : index ∈ classSourceSuffix classes classId start) : index < classes.size := by
  exact List.mem_range.mp (List.mem_of_mem_drop (List.mem_filter.mp hi).1)

/-- Exact physical row contents of the actual collection loop, including
the append order and multiplicity of source indices. -/
theorem collectCandidateClasses_row
    (classes : Array (Option Nat)) (start : Nat) (before after : Array (Array Int32))
    (hrun : Generated.StrictRecombine.collectCandidateClasses classes start before = .ok after)
    (classId : Nat) (hc : classId < before.size) :
    after[classId]!.toList = before[classId]!.toList ++
      (classSourceSuffix classes classId start).map (fun index => index.toUInt32.toInt32) := by
  induction hmeasure : classes.size - start using Nat.strong_induction_on
      generalizing start before with
  | h measure ih =>
    rw [Generated.StrictRecombine.collectCandidateClasses] at hrun
    split at hrun
    next hs =>
      split at hrun
      next hn => contradiction
      next current hcurrent =>
        split at hrun
        next hb =>
          let updated := before.set current (before[current].push start.toUInt32.toInt32)
          have ht := ih (classes.size - (start + 1)) (by omega) (start + 1)
            updated hrun (by simpa [updated] using hc) rfl
          rw [classSourceSuffix_step classes classId start hs]
          have hcurrentBang : classes[start]! = some current := by
            simpa only [getElem!_pos classes start hs] using hcurrent
          rw [hcurrentBang]
          by_cases he : classId = current
          · subst classId
            have hrow : updated[current]! = before[current]!.push start.toUInt32.toInt32 := by
              simp [updated, getElem!_pos, hb]
            rw [hrow] at ht
            simpa only [↓reduceIte, Array.toList_push, List.map_cons, List.append_assoc,
              List.singleton_append] using ht
          · have hrow : updated[classId]! = before[classId]! := by
              simp [updated, getElem!_pos, hc, Ne.symm he]
            rw [hrow] at ht
            simpa [Ne.symm he] using ht
        next hb => contradiction
    next hs =>
      have he := Except.ok.inj hrun
      subst after
      have hdrop : (List.range classes.size).drop start = [] :=
        List.drop_eq_nil_iff.mpr (by simpa using Nat.le_of_not_gt hs)
      simp [classSourceSuffix, hdrop]

/-- Decoding the actual collected row gives the exact filtered source list,
so its multiplicities do not depend on a later successful validation. -/
theorem collectCandidateClasses_decoded_row
    (classes : Array (Option Nat)) (classCount : Nat) (groups : Array (Array Int32))
    (hfits : classes.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.collectCandidateClasses classes 0
      (Array.replicate classCount #[]) = .ok groups)
    (classId : Nat) (hc : classId < groups.size) :
    groups[classId]!.toList.map (fun value => value.toInt64.toNat) =
      classSourceSuffix classes classId 0 := by
  have hsize := collectCandidateClasses_size classes 0 (Array.replicate classCount #[]) groups hrun
  have hc' : classId < (Array.replicate classCount (#[] : Array Int32)).size := by omega
  have hrow := collectCandidateClasses_row classes 0 (Array.replicate classCount #[]) groups hrun classId hc'
  have hempty : (Array.replicate classCount (#[] : Array Int32))[classId]! = #[] := by
    rw [getElem!_pos _ classId hc']; simp
  rw [hempty] at hrow
  simp only [Array.toList_empty, List.nil_append] at hrow
  rw [hrow, List.map_map]
  conv_rhs => rw [← List.map_id' (classSourceSuffix classes classId 0)]
  apply List.map_congr_left
  intro index hi
  exact (nat_toUInt32_toInt32_nonnegative_and_toNat index
    (lt_of_lt_of_le (classSourceSuffix_bound classes classId 0 index hi) hfits)).2

theorem extractCandidates_decoded_row_nodup
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32)) (hfits : count ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (classId : Nat) (hc : classId < groups.size) :
    (groups[classId]!.toList.map (fun value => value.toInt64.toNat)).Nodup := by
  rw [Generated.StrictRecombine.extractCandidates] at hrun
  split at hrun
  next hempty =>
    have he := Except.ok.inj hrun
    rw [← he] at hc
    simp at hc
  next hnonempty =>
    cases hp : Generated.StrictRecombine.partitionCandidateColumns transform rows count
        0 0 (Array.replicate count none) with
    | error fault => rw [hp] at hrun; contradiction
    | ok result =>
      obtain ⟨classes, classCount⟩ := result
      rw [hp] at hrun
      have hsize := (partitionCandidateColumns_from_empty_all_assigned transform rows count
        classes classCount hp).1
      rw [collectCandidateClasses_decoded_row classes classCount groups (by omega) hrun classId hc]
      exact classSourceSuffix_nodup classes classId 0

end Refinement.StrictRecombine
