import CLPoly.Refinement.VanHoeijExactGrouping

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Every allocated class retains a physical source-column witness. -/
theorem partitionCandidateColumns_surjective_classes
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count column classCount : Nat) (classes output : Array (Option Nat))
    (outputCount : Nat)
    (hsurj : ∀ id, id < classCount → ∃ source, source < classes.size ∧
      classes[source]! = some id)
    (hrun : Generated.StrictRecombine.partitionCandidateColumns transform rows count
      column classCount classes = .ok (output, outputCount)) :
    ∀ id, id < outputCount → ∃ source, source < output.size ∧ output[source]! = some id := by
  induction hmeasure : count - column using Nat.strong_induction_on
      generalizing column classCount classes output outputCount with
  | h measure ih =>
    rw [Generated.StrictRecombine.partitionCandidateColumns] at hrun
    split at hrun
    next hc =>
      split at hrun
      next hcs =>
        split at hrun
        next existing =>
          exact ih (count - (column + 1)) (by omega) (column + 1) classCount
            classes output outputCount hsurj hrun rfl
        next hnone =>
          let seeded := classes.set column (some classCount)
          have hseed : ∀ id, id < classCount + 1 → ∃ source,
              source < seeded.size ∧ seeded[source]! = some id := by
            intro id hi
            by_cases he : id = classCount
            · subst id
              exact ⟨column, by simpa [seeded] using hcs,
                by simp [seeded, getElem!_pos, hcs]⟩
            · obtain ⟨source, hs, hv⟩ := hsurj id (by omega)
              have hne : source ≠ column := by
                intro h
                subst source
                rw [getElem!_pos classes column hcs] at hv
                simp [hnone] at hv
              refine ⟨source, by simpa [seeded] using hs, ?_⟩
              simpa [seeded, getElem!_pos, hs, Ne.symm hne] using hv
          cases ha : Generated.StrictRecombine.assignCandidateClass transform rows column
              count classCount (column + 1) seeded with
          | error fault => simp [seeded, ha] at hrun
          | ok assigned =>
            simp only [seeded, ha] at hrun
            have hsize := assignCandidateClass_size transform rows column count classCount
              (column + 1) seeded assigned ha
            have hsurj' : ∀ id, id < classCount + 1 → ∃ source,
                source < assigned.size ∧ assigned[source]! = some id := by
              intro id hi
              obtain ⟨source, hs, hv⟩ := hseed id hi
              exact ⟨source, by omega, assignCandidateClass_preserves_some transform rows
                column count classCount (column + 1) seeded assigned ha source id hs hv⟩
            exact ih (count - (column + 1)) (by omega) (column + 1) (classCount + 1)
              assigned output outputCount hsurj' hrun rfl
      next hcs => contradiction
    next hc =>
      have he := Except.ok.inj hrun
      cases he
      exact hsurj

/-- No returned candidate class is empty: each has an actual source column.
The empty retained-row branch returns no classes and satisfies this vacuously. -/
theorem extractCandidates_every_group_has_source
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32))
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (classId : Nat) (hc : classId < groups.size) :
    ∃ source, source < count ∧ CandidateMember groups classId source.toUInt32.toInt32 := by
  by_cases hempty : rows.isEmpty = true
  · simp [Generated.StrictRecombine.extractCandidates, hempty] at hrun
    subst groups
    simp at hc
  · have hn : rows.isEmpty = false := Bool.eq_false_of_not_eq_true hempty
    simp only [Generated.StrictRecombine.extractCandidates, hn,
      Bool.false_eq_true, ↓reduceIte] at hrun
    cases hp : Generated.StrictRecombine.partitionCandidateColumns transform rows count
        0 0 (Array.replicate count none) with
    | error fault => rw [hp] at hrun; contradiction
    | ok result =>
      obtain ⟨classes, classCount⟩ := result
      rw [hp] at hrun
      have hs := partitionCandidateColumns_size transform rows count 0 0
        (Array.replicate count none) classes classCount hp
      have hg := collectCandidateClasses_size classes 0
        (Array.replicate classCount #[]) groups hrun
      simp only [Array.size_replicate] at hs hg
      obtain ⟨source, hsource, hv⟩ := partitionCandidateColumns_surjective_classes transform
        rows count 0 0 (Array.replicate count none) classes classCount
        (by intro id hi; omega) hp classId (by omega)
      exact ⟨source, by omega, collectCandidateClasses_includes_source classes 0
        (Array.replicate classCount #[]) groups hrun source classId (Nat.zero_le _) hsource hv⟩

/-- The same physical column cannot belong to two distinct output classes. -/
theorem extractCandidates_class_unique
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32)) (hfits : count ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (source left right : Nat) (hs : source < count)
    (hl : CandidateMember groups left source.toUInt32.toInt32)
    (hr : CandidateMember groups right source.toUInt32.toInt32) : left = right := by
  by_cases hempty : rows.isEmpty = true
  · simp [Generated.StrictRecombine.extractCandidates, hempty] at hrun
    subst groups
    have := hl.1
    simp at this
  · have hn : rows.isEmpty = false := Bool.eq_false_of_not_eq_true hempty
    simp only [Generated.StrictRecombine.extractCandidates, hn,
      Bool.false_eq_true, ↓reduceIte] at hrun
    cases hp : Generated.StrictRecombine.partitionCandidateColumns transform rows count
        0 0 (Array.replicate count none) with
    | error fault => rw [hp] at hrun; contradiction
    | ok result =>
      obtain ⟨classes, classCount⟩ := result
      rw [hp] at hrun
      have hsize := partitionCandidateColumns_size transform rows count 0 0
        (Array.replicate count none) classes classCount hp
      simp only [Array.size_replicate] at hsize
      have he₁ := (collectCandidateClasses_member_iff_class classes classCount groups
        (by omega) hrun source left (by omega)).mp hl
      have he₂ := (collectCandidateClasses_member_iff_class classes classCount groups
        (by omega) hrun source right (by omega)).mp hr
      exact Option.some.inj (he₁.symm.trans he₂)

end Refinement.StrictRecombine
