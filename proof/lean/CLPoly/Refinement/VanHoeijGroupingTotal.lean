import CLPoly.Refinement.VanHoeijGroupingComplete

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem candidateColumnsEqual_succeeds
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count left right index : Nat) (hvalid : CandidateColumnsValid transform rows count)
    (hl : left < count) (hr : right < count) :
    ∃ equal, Generated.StrictRecombine.candidateColumnsEqual transform rows left right index = .ok equal := by
  induction hmeasure : rows.size - index using Nat.strong_induction_on generalizing index with
  | h measure ih =>
    rw [Generated.StrictRecombine.candidateColumnsEqual]
    by_cases hi : index < rows.size
    · rw [dif_pos hi]
      obtain ⟨hrow, hcols⟩ := hvalid index hi
      have hrow' : rows[index] < transform.size := by
        simpa only [getElem!_pos rows index hi] using hrow
      have hcols' : count ≤ transform[rows[index]].size := by
        simpa only [getElem!_pos rows index hi,
          getElem!_pos transform rows[index] hrow'] using hcols
      rw [dif_pos hrow', dif_pos (lt_of_lt_of_le hl hcols'), dif_pos (lt_of_lt_of_le hr hcols')]
      split
      · exact ih (rows.size - (index + 1)) (by omega) (index + 1) rfl
      · exact ⟨false, rfl⟩
    · rw [dif_neg hi]
      exact ⟨true, rfl⟩

theorem assignCandidateClass_succeeds
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count representative classId start : Nat) (classes : Array (Option Nat))
    (hvalid : CandidateColumnsValid transform rows count) (hrep : representative < count)
    (hsize : count ≤ classes.size) :
    ∃ output, Generated.StrictRecombine.assignCandidateClass transform rows
      representative count classId start classes = .ok output := by
  induction hmeasure : count - start using Nat.strong_induction_on generalizing start classes with
  | h measure ih =>
    rw [Generated.StrictRecombine.assignCandidateClass]
    by_cases hs : start < count
    · rw [dif_pos hs, dif_pos (lt_of_lt_of_le hs hsize)]
      cases hc : classes[start] with
      | some id =>
        exact ih (count - (start + 1)) (by omega) (start + 1) classes hsize rfl
      | none =>
        obtain ⟨equal, he⟩ := candidateColumnsEqual_succeeds transform rows count representative start
          0 hvalid hrep hs
        rw [he]
        cases equal with
        | false => exact ih (count - (start + 1)) (by omega) (start + 1) classes hsize rfl
        | true =>
          exact ih (count - (start + 1)) (by omega) (start + 1)
            (classes.set start (some classId)) (by simpa using hsize) rfl
    · rw [dif_neg hs]
      exact ⟨classes, rfl⟩

theorem partitionCandidateColumns_succeeds
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count column classCount : Nat) (classes : Array (Option Nat))
    (hvalid : CandidateColumnsValid transform rows count) (hsize : classes.size = count)
    (hbound : CandidateClassIdsBound classCount classes) :
    ∃ output outputCount,
      Generated.StrictRecombine.partitionCandidateColumns transform rows count column classCount
        classes = .ok (output, outputCount) ∧ CandidateClassIdsBound outputCount output := by
  induction hmeasure : count - column using Nat.strong_induction_on
      generalizing column classCount classes with
  | h measure ih =>
    rw [Generated.StrictRecombine.partitionCandidateColumns]
    by_cases hc : column < count
    · have hcs : column < classes.size := by omega
      rw [dif_pos hc, dif_pos hcs]
      cases hv : classes[column] with
      | some id =>
        exact ih (count - (column + 1)) (by omega) (column + 1) classCount classes hsize hbound rfl
      | none =>
        let seeded := classes.set column (some classCount)
        have hs : seeded.size = count := by simp [seeded, hsize]
        have hb : CandidateClassIdsBound (classCount + 1) seeded := by
          intro index id hi he
          by_cases hie : index = column
          · subst index
            have he' : some classCount = some id := by
              rw [getElem!_pos seeded column hi] at he
              simpa [seeded] using he
            have hid := Option.some.inj he'
            omega
          · have hi' : index < classes.size := by simpa [seeded] using hi
            have he' : classes[index]! = some id := by
              simpa [seeded, getElem!_pos, hi', Ne.symm hie] using he
            exact Nat.lt_succ_of_lt (hbound index id hi' he')
        obtain ⟨assigned, ha⟩ := assignCandidateClass_succeeds transform rows count column
          classCount (column + 1) seeded hvalid hc (by omega)
        dsimp only
        rw [ha]
        have hsa : assigned.size = count :=
          (assignCandidateClass_size transform rows column count classCount (column + 1)
            seeded assigned ha).trans hs
        exact ih (count - (column + 1)) (by omega) (column + 1) (classCount + 1)
          assigned hsa (assignCandidateClass_ids_bound transform rows column count classCount
            (column + 1) seeded assigned hb ha) rfl
    · rw [dif_neg hc]
      exact ⟨classes, classCount, rfl, hbound⟩

theorem collectCandidateClasses_succeeds
    (classes : Array (Option Nat)) (column : Nat) (result : Array (Array Int32))
    (hbound : CandidateClassIdsBound result.size classes)
    (hassigned : ∀ index, index < classes.size → classes[index]! ≠ none) :
    ∃ output, Generated.StrictRecombine.collectCandidateClasses classes column result = .ok output := by
  induction hmeasure : classes.size - column using Nat.strong_induction_on
      generalizing column result with
  | h measure ih =>
    rw [Generated.StrictRecombine.collectCandidateClasses]
    by_cases hc : column < classes.size
    · rw [dif_pos hc]
      cases hv : classes[column] with
      | none =>
        exact (hassigned column hc (by simpa only [getElem!_pos classes column hc] using hv)).elim
      | some id =>
        have hi : id < result.size := hbound column id hc
          (by simpa only [getElem!_pos classes column hc] using hv)
        dsimp only
        rw [dif_pos hi]
        exact ih (classes.size - (column + 1)) (by omega) (column + 1)
          (result.set id (result[id].push column.toUInt32.toInt32)) (by simpa using hbound) rfl
    · rw [dif_neg hc]
      exact ⟨result, rfl⟩

/-- Array shape alone suffices for the entire generated grouping function to
return normally; no successful-run premise is needed. -/
theorem extractCandidates_succeeds
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat) (count : Nat)
    (hvalid : CandidateColumnsValid transform rows count) :
    ∃ groups, Generated.StrictRecombine.extractCandidates rows transform count = .ok groups := by
  rw [Generated.StrictRecombine.extractCandidates]
  split
  · exact ⟨#[], rfl⟩
  · obtain ⟨classes, classCount, hp, hb⟩ := partitionCandidateColumns_succeeds transform rows
      count 0 0 (Array.replicate count none) hvalid (by simp)
      (by
        intro index id hi he
        rw [getElem!_pos _ index hi] at he
        simp at he)
    rw [hp]
    obtain ⟨hs, ha⟩ := partitionCandidateColumns_from_empty_all_assigned transform rows count
      classes classCount hp
    exact collectCandidateClasses_succeeds classes 0 (Array.replicate classCount #[])
      (by simpa using hb) (fun index hi => ha index (by omega))

/-- The same actual LLL output supplies a successful grouping execution, and
every pair constant on retained rows is put in a common returned group. -/
theorem lllReduceWithNorms_grouping_complete
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (count : Nat) (hcount : count ≤ matrix.size)
    (hnonempty : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = false) :
    ∃ groups,
      Generated.StrictRecombine.extractCandidates
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) output.2.1 count = .ok groups ∧
      ∀ left right : Fin count,
        (∀ row : Fin matrix.size,
          row.val < Generated.StrictVanHoeij.retainedPrefixCount output.2.2.2
            (bound : QQ) output.2.2.2.size →
          (output.2.1[row.val]!)[left.val]! = (output.2.1[row.val]!)[right.val]!) →
        ∃ classId, CandidateMember groups classId left.val.toUInt32.toInt32 ∧
          CandidateMember groups classId right.val.toUInt32.toInt32 := by
  have hvalid := lllReduceWithNorms_candidate_columns_valid matrix hinput bound count hcount output hrun
  obtain ⟨groups, hg⟩ := extractCandidates_succeeds output.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) count hvalid
  refine ⟨groups, hg, ?_⟩
  intro left right hequal
  apply extractCandidates_equal_columns_same_group output.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) count groups
    hvalid hnonempty hg left.val right.val left.isLt right.isLt
  exact lllReduceWithNorms_columns_equal_of_retained matrix hinput bound output hrun
    ⟨left.val, lt_of_lt_of_le left.isLt hcount⟩
    ⟨right.val, lt_of_lt_of_le right.isLt hcount⟩ hequal

end Refinement.StrictRecombine
