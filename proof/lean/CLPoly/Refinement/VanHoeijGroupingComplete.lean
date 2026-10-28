import CLPoly.Refinement.VanHoeijActualSupportSeparation

set_option autoImplicit false

namespace Refinement.StrictRecombine

open scoped Classical

theorem shortSpanFromNorms_size_eq (norms : Array QQ) (bound : ZZ) :
    (Generated.StrictVanHoeij.shortSpanFromNorms norms bound).size =
      Generated.StrictVanHoeij.retainedPrefixCount norms (bound : QQ) norms.size := by
  exact (iotaLoop_spec _ 0 _ (Array.size_replicate ..)).1

/-- Row equality at the retained prefix makes the actual column comparator
return true; its array bounds come from the same LLL execution. -/
theorem lllReduceWithNorms_columns_equal_of_retained
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (left right : Fin matrix.size)
    (hequal : ∀ row : Fin matrix.size,
      row.val < Generated.StrictVanHoeij.retainedPrefixCount output.2.2.2
        (bound : QQ) output.2.2.2.size →
      (output.2.1[row.val]!)[left.val]! = (output.2.1[row.val]!)[right.val]!) :
    Generated.StrictRecombine.candidateColumnsEqual output.2.1
      (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
      left.val right.val 0 = .ok true := by
  have hvalid := lllReduceWithNorms_candidate_columns_valid matrix hinput bound
    matrix.size le_rfl output hrun
  obtain ⟨state, hv, _, _, hn, hrel, _⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  have hnormSize : output.2.2.2.size = matrix.size := by
    rw [← hn, hv.norms_size, hrel.1]
  obtain ⟨hsize, hentries, _⟩ := shortSpanFromNorms_spec output.2.2.2 bound
  apply candidateColumnsEqual_true_complete
  intro position _ hp
  obtain ⟨hrow, hcols⟩ := hvalid position hp
  refine ⟨hrow, lt_of_lt_of_le left.isLt hcols, lt_of_lt_of_le right.isLt hcols, ?_⟩
  dsimp only
  rw [hentries position hp]
  exact hequal ⟨position, by omega⟩ (by
    rw [shortSpanFromNorms_size_eq] at hp
    exact hp)

/-- Exact effect of the literal inner scan on any physical class slot. -/
theorem assignCandidateClass_entry
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (representative count classId start : Nat)
    (classes output : Array (Option Nat))
    (hrun : Generated.StrictRecombine.assignCandidateClass transform rows
      representative count classId start classes = .ok output)
    (column : Nat) (hcolumn : column < classes.size) :
    output[column]! =
      if start ≤ column ∧ column < count ∧ classes[column]! = none ∧
        Generated.StrictRecombine.candidateColumnsEqual transform rows
          representative column 0 = .ok true then some classId else classes[column]! := by
  induction hmeasure : count - start using Nat.strong_induction_on
      generalizing start classes with
  | h measure ih =>
    rw [Generated.StrictRecombine.assignCandidateClass] at hrun
    split at hrun
    next hs =>
      split at hrun
      next hcs =>
        split at hrun
        next existing he =>
          have ht := ih (count - (start + 1)) (by omega) (start + 1)
            classes hrun hcolumn rfl
          by_cases hc : column = start
          · subst column
            simpa [getElem!_pos classes start hcs, he] using ht
          · simpa [show (start + 1 ≤ column) ↔ start ≤ column by omega] using ht
        next hn =>
          cases heq : Generated.StrictRecombine.candidateColumnsEqual transform rows
              representative start 0 with
          | error fault => simp [heq] at hrun
          | ok equal =>
            simp only [heq] at hrun
            cases equal with
            | false =>
              have ht := ih (count - (start + 1)) (by omega) (start + 1)
                classes hrun hcolumn rfl
              by_cases hc : column = start
              · subst column
                simpa [heq] using ht
              · simpa [show (start + 1 ≤ column) ↔ start ≤ column by omega] using ht
            | true =>
              have ht := ih (count - (start + 1)) (by omega) (start + 1)
                (classes.set start (some classId)) hrun (by simpa using hcolumn) rfl
              by_cases hc : column = start
              · subst column
                simpa [getElem!_pos, hcs, hn, hs, heq] using ht
              · have hget : (classes.set start (some classId))[column]! = classes[column]! := by
                  simp [getElem!_pos, hcolumn, Ne.symm hc]
                rw [hget] at ht
                simpa [show (start + 1 ≤ column) ↔ start ≤ column by omega] using ht
      next hcs => contradiction
    next hs =>
      have ho := Except.ok.inj hrun
      subst output
      have hfalse : ¬ (start ≤ column ∧ column < count ∧ classes[column]! = none ∧
          Generated.StrictRecombine.candidateColumnsEqual transform rows
            representative column 0 = .ok true) := by omega
      simp only [hfalse, ↓reduceIte]

/-- Equal columns have identical optional class identifiers, including before
either column has been assigned. -/
def CandidateClassesClosed
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (classes : Array (Option Nat)) : Prop :=
  ∀ left right, left < count → right < count →
    Generated.StrictRecombine.candidateColumnsEqual transform rows left right 0 = .ok true →
    classes[left]! = classes[right]!

/-- One seeded outer iteration assigns precisely the unassigned columns
equivalent to its representative. -/
theorem seed_assignCandidateClass_entry
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count representative classId : Nat) (classes output : Array (Option Nat))
    (hvalid : CandidateColumnsValid transform rows count)
    (hsize : classes.size = count) (hrep : representative < count)
    (hnone : classes[representative]! = none)
    (hprefix : ∀ index, index < representative → classes[index]! ≠ none)
    (hrun : Generated.StrictRecombine.assignCandidateClass transform rows representative
      count classId (representative + 1)
      (classes.set representative (some classId)) = .ok output)
    (column : Nat) (hc : column < count) :
    output[column]! = if classes[column]! = none ∧
        Generated.StrictRecombine.candidateColumnsEqual transform rows
          representative column 0 = .ok true then some classId else classes[column]! := by
  have ht := assignCandidateClass_entry transform rows representative count classId
    (representative + 1) (classes.set representative (some classId)) output hrun
    column (by simpa [hsize] using hc)
  by_cases he : column = representative
  · subst column
    have hrefl := candidateColumnsEqual_refl transform rows count representative hvalid hrep
    rw [hnone]
    simpa [getElem!_pos, hsize, hrep, hrefl] using ht
  · have hget : (classes.set representative (some classId))[column]! = classes[column]! := by
      simp [getElem!_pos, hsize, hc, Ne.symm he]
    rw [hget] at ht
    by_cases hbefore : column < representative
    · simpa [hprefix column hbefore] using ht
    · simpa [show representative + 1 ≤ column by omega, hc] using ht

theorem seed_assignCandidateClass_closed
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count representative classId : Nat) (classes output : Array (Option Nat))
    (hvalid : CandidateColumnsValid transform rows count)
    (hsize : classes.size = count) (hrep : representative < count)
    (hnone : classes[representative]! = none)
    (hprefix : ∀ index, index < representative → classes[index]! ≠ none)
    (hclosed : CandidateClassesClosed transform rows count classes)
    (hrun : Generated.StrictRecombine.assignCandidateClass transform rows representative
      count classId (representative + 1)
      (classes.set representative (some classId)) = .ok output) :
    CandidateClassesClosed transform rows count output := by
  intro left right hl hr heq
  rw [seed_assignCandidateClass_entry transform rows count representative classId classes
    output hvalid hsize hrep hnone hprefix hrun left hl,
    seed_assignCandidateClass_entry transform rows count representative classId classes
    output hvalid hsize hrep hnone hprefix hrun right hr]
  have hvalues := hclosed left right hl hr heq
  have htests :
      (Generated.StrictRecombine.candidateColumnsEqual transform rows representative left 0 = .ok true) ↔
      (Generated.StrictRecombine.candidateColumnsEqual transform rows representative right 0 = .ok true) := by
    constructor
    · intro h
      exact candidateColumnsEqual_trans transform rows representative left right h heq
    · intro h
      exact candidateColumnsEqual_trans transform rows representative right left h
        (candidateColumnsEqual_symm transform rows left right heq)
  simp only [hvalues, htests]

theorem partitionCandidateColumns_preserves_closed
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count column classCount : Nat) (classes output : Array (Option Nat))
    (outputCount : Nat) (hvalid : CandidateColumnsValid transform rows count)
    (hsize : classes.size = count)
    (hprefix : ∀ index, index < column → classes[index]! ≠ none)
    (hclosed : CandidateClassesClosed transform rows count classes)
    (hrun : Generated.StrictRecombine.partitionCandidateColumns transform rows count
      column classCount classes = .ok (output, outputCount)) :
    CandidateClassesClosed transform rows count output := by
  induction hmeasure : count - column using Nat.strong_induction_on
      generalizing column classCount classes with
  | h measure ih =>
    rw [Generated.StrictRecombine.partitionCandidateColumns] at hrun
    split at hrun
    next hc =>
      split at hrun
      next hcs =>
        split at hrun
        next existing he =>
          have hprefix' : ∀ index, index < column + 1 → classes[index]! ≠ none := by
            intro index hi
            by_cases heq : index = column
            · subst index
              simp [getElem!_pos classes column hcs, he]
            · exact hprefix index (by omega)
          exact ih (count - (column + 1)) (by omega) (column + 1) classCount
            classes hsize hprefix' hclosed hrun rfl
        next hn =>
          have hnone : classes[column]! = none := by
            simpa only [getElem!_pos classes column hcs] using hn
          let seeded := classes.set column (some classCount)
          cases ha : Generated.StrictRecombine.assignCandidateClass transform rows column
              count classCount (column + 1) seeded with
          | error fault => simp [seeded, ha] at hrun
          | ok assigned =>
            simp only [seeded, ha] at hrun
            have hclosed' := seed_assignCandidateClass_closed transform rows count column
              classCount classes assigned hvalid hsize hc hnone hprefix hclosed ha
            have hsize' : assigned.size = count := by
              have h := assignCandidateClass_size transform rows column count classCount
                (column + 1) seeded assigned ha
              simpa [seeded, hsize] using h
            have hprefix' : ∀ index, index < column + 1 → assigned[index]! ≠ none := by
              intro index hi
              rw [seed_assignCandidateClass_entry transform rows count column classCount
                classes assigned hvalid hsize hc hnone hprefix ha index (by omega)]
              by_cases heq : index = column
              · subst index
                simp [hnone, candidateColumnsEqual_refl transform rows count column hvalid hc]
              · simp [hprefix index (by omega)]
            exact ih (count - (column + 1)) (by omega) (column + 1) (classCount + 1)
              assigned hsize' hprefix' hclosed' hrun rfl
      next hcs => contradiction
    next hc =>
      have ho := Except.ok.inj hrun
      injection ho with ho _
      subst output
      exact hclosed

/-- Completeness of the actual grouping scan: it never splits equal columns
between distinct physical output groups. -/
theorem extractCandidates_equal_columns_same_group
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32))
    (hvalid : CandidateColumnsValid transform rows count)
    (hnonempty : rows.isEmpty = false)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (left right : Nat) (hl : left < count) (hr : right < count)
    (heq : Generated.StrictRecombine.candidateColumnsEqual transform rows left right 0 = .ok true) :
    ∃ classId, CandidateMember groups classId left.toUInt32.toInt32 ∧
      CandidateMember groups classId right.toUInt32.toInt32 := by
  simp only [Generated.StrictRecombine.extractCandidates, hnonempty,
    Bool.false_eq_true, ↓reduceIte] at hrun
  cases hp : Generated.StrictRecombine.partitionCandidateColumns transform rows count
      0 0 (Array.replicate count none) with
  | error fault => rw [hp] at hrun; contradiction
  | ok result =>
    obtain ⟨classes, classCount⟩ := result
    rw [hp] at hrun
    have hcovered := partitionCandidateColumns_from_empty_all_assigned transform rows count
      classes classCount hp
    have hclosed := partitionCandidateColumns_preserves_closed transform rows count 0 0
      (Array.replicate count none) classes classCount hvalid (by simp)
      (by intro index hi; omega)
      (by intro i j hi hj _; simp [getElem!_pos, hi, hj]) hp
    have hsame := hclosed left right hl hr heq
    cases hc : classes[left]! with
    | none => exact (hcovered.2 left hl hc).elim
    | some classId =>
      refine ⟨classId, ?_, ?_⟩
      · exact collectCandidateClasses_includes_source classes 0 (Array.replicate classCount #[])
          groups hrun left classId (Nat.zero_le _) (by omega) hc
      · exact collectCandidateClasses_includes_source classes 0 (Array.replicate classCount #[])
          groups hrun right classId (Nat.zero_le _) (by omega) (hsame.symm.trans hc)

end Refinement.StrictRecombine
