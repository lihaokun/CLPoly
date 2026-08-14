/-
  No strict recombination refinement theorem is exported yet.

  The former theorem equated the generated C++ result with an arbitrary UFD
  existence witness and contained a placeholder proof.  Both the witness-based
  implementation and the claimed refinement have been removed.  A replacement
  must prove the generated candidate/extraction loops directly.
-/
import CLPoly.Algorithm.Recombine
import CLPoly.Generated.StrictRecombine
import CLPoly.Refinement.Basic
import CLPoly.Refinement.Hensel
import Batteries.Data.Array.Lemmas

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictRecombine

private def positionalCode (base : Nat) : List Nat → Nat
  | [] => 0
  | digit :: rest => digit * base ^ rest.length + positionalCode base rest

private theorem positionalCode_lt_pow (base : Nat) (hbase : 0 < base)
    (digits : List Nat) (hdigits : ∀ digit ∈ digits, digit < base) :
    positionalCode base digits < base ^ digits.length := by
  induction digits with
  | nil => simp [positionalCode, hbase]
  | cons digit rest ih =>
      have hdigit := hdigits digit (by simp)
      have hrest : ∀ value ∈ rest, value < base := by
        intro value hvalue
        exact hdigits value (by simp [hvalue])
      have htail := ih hrest
      simp only [positionalCode, List.length_cons, pow_succ]
      have hpow : 0 < base ^ rest.length := Nat.pow_pos hbase
      nlinarith

private theorem positionalCode_append (base : Nat) (leading suffix : List Nat) :
    positionalCode base (leading ++ suffix) =
      positionalCode base leading * base ^ suffix.length +
        positionalCode base suffix := by
  induction leading with
  | nil => simp [positionalCode]
  | cons digit rest ih =>
      simp only [List.cons_append, positionalCode, List.length_append,
        List.length_cons]
      rw [ih, pow_add]
      ring

private theorem positionalCode_pivot_lt (base : Nat) (hbase : 0 < base)
    (leading leftSuffix rightSuffix : List Nat) (leftPivot rightPivot : Nat)
    (hlength : leftSuffix.length = rightSuffix.length)
    (hpivot : leftPivot < rightPivot)
    (hleftDigits : ∀ digit ∈ leftSuffix, digit < base)
    (hrightDigits : ∀ digit ∈ rightSuffix, digit < base) :
    positionalCode base (leading ++ leftPivot :: leftSuffix) <
      positionalCode base (leading ++ rightPivot :: rightSuffix) := by
  rw [positionalCode_append, positionalCode_append]
  simp only [positionalCode, List.length_cons, hlength]
  have hleftTail := positionalCode_lt_pow base hbase leftSuffix hleftDigits
  have hleftTail' : positionalCode base leftSuffix < base ^ rightSuffix.length := by
    simpa [hlength] using hleftTail
  apply Nat.add_lt_add_left
  calc
    leftPivot * base ^ rightSuffix.length + positionalCode base leftSuffix <
        leftPivot * base ^ rightSuffix.length + base ^ rightSuffix.length :=
      Nat.add_lt_add_left hleftTail' _
    _ = (leftPivot + 1) * base ^ rightSuffix.length := by ring
    _ ≤ rightPivot * base ^ rightSuffix.length :=
      Nat.mul_le_mul_right _ (by omega)
    _ ≤ rightPivot * base ^ rightSuffix.length + positionalCode base rightSuffix :=
      Nat.le_add_right _ _

private theorem array_toList_pivot (values : Array Nat) (pivot : Nat)
    (hpivot : pivot < values.size) :
    values.toList = values.toList.take pivot ++ values[pivot]! ::
      values.toList.drop (pivot + 1) := by
  calc
    values.toList = values.toList.take pivot ++ values.toList.drop pivot :=
      (List.take_append_drop pivot values.toList).symm
    _ = values.toList.take pivot ++ values.toList[pivot] ::
        values.toList.drop (pivot + 1) := by
      rw [List.drop_eq_getElem_cons (by simpa using hpivot)]
    _ = values.toList.take pivot ++ values[pivot]! ::
        values.toList.drop (pivot + 1) := by
      congr 2
      rw [getElem!_pos values pivot hpivot]
      exact Array.getElem_toList hpivot

private theorem positionalCode_array_pivot_lt (base : Nat) (hbase : 0 < base)
    (left right : Array Nat) (pivot : Nat)
    (hpivotLeft : pivot < left.size) (hpivotRight : pivot < right.size)
    (hsize : left.size = right.size)
    (hprefix : left.toList.take pivot = right.toList.take pivot)
    (hpivot : left[pivot]! < right[pivot]!)
    (hleftDigits : ∀ digit ∈ left.toList, digit < base)
    (hrightDigits : ∀ digit ∈ right.toList, digit < base) :
    positionalCode base left.toList < positionalCode base right.toList := by
  rw [array_toList_pivot left pivot hpivotLeft,
    array_toList_pivot right pivot hpivotRight, hprefix]
  apply positionalCode_pivot_lt base hbase
  · simp [hsize]
  · exact hpivot
  · intro digit hdigit
    exact hleftDigits digit (List.mem_of_mem_drop hdigit)
  · intro digit hdigit
    exact hrightDigits digit (List.mem_of_mem_drop hdigit)

theorem initialCombinationLoop_toList (count index : Nat)
    (result : Array Nat) :
    (Generated.StrictRecombine.initialCombinationLoop count index result).toList =
      result.toList ++ List.range' index (count - index) := by
  rw [Generated.StrictRecombine.initialCombinationLoop]
  split
  next hmore =>
    rw [initialCombinationLoop_toList]
    rw [Array.toList_push, List.append_assoc]
    congr 1
    rw [show count - index = (count - (index + 1)) + 1 by omega,
      List.range'_succ]
    simp
  next hdone =>
    have hzero : count - index = 0 := Nat.sub_eq_zero_of_le (by omega)
    simp [hzero]
termination_by count - index
decreasing_by omega

theorem initialCombination_toList (count : Nat) :
    (Generated.StrictRecombine.initialCombination count).toList =
      List.range count := by
  simp [Generated.StrictRecombine.initialCombination,
    initialCombinationLoop_toList, List.range_eq_range']

theorem initialCombination_size (count : Nat) :
    (Generated.StrictRecombine.initialCombination count).size = count := by
  simpa using congrArg List.length (initialCombination_toList count)

theorem resetCombinationSuffix_size (indices : Array Nat)
    (pivot offset : Nat) :
    (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset).size =
      indices.size := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    rw [resetCombinationSuffix_size]
    simp
  next hposition => rfl
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

theorem resetCombinationSuffix_getElem_le (indices : Array Nat)
    (pivot offset index : Nat) (hindex : index < indices.size)
    (hle : index ≤ pivot) :
    (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset)[index]! =
      indices[index]! := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    rw [resetCombinationSuffix_getElem_le
      (indices.set (pivot + 1 + offset) (indices[pivot + offset] + 1)) pivot
      (offset + 1) index (by simpa using hindex) hle]
    have hne : index ≠ pivot + 1 + offset := by omega
    rw [getElem!_pos _ index (by simpa using hindex),
      getElem!_pos indices index hindex]
    rw [Array.getElem_set]
    rw [if_neg (by omega)]
  next hposition => rfl
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

private def ValidCombination (upper count : Nat) (indices : Array Nat) : Prop :=
  indices.size = count ∧
    ∀ index (hindex : index < indices.size),
      indices[index]! ≤ upper - count + index

private theorem resetCombinationSuffix_preserves_bounds (indices : Array Nat)
    (upper count pivot offset : Nat)
    (hbound : ∀ index (hindex : index < indices.size),
      indices[index]! ≤ upper - count + index) :
    ∀ index (hindex : index <
      (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset).size),
      (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset)[index]! ≤
        upper - count + index := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    apply resetCombinationSuffix_preserves_bounds
    intro index hindex
    rw [getElem!_pos _ index hindex, Array.getElem_set]
    split
    next heq =>
      subst index
      have hprevious : pivot + offset < indices.size := by omega
      have := hbound (pivot + offset) hprevious
      rw [getElem!_pos indices (pivot + offset) hprevious] at this
      omega
    next hne =>
      have hi : index < indices.size := by simpa using hindex
      have hb := hbound index hi
      rw [getElem!_pos indices index hi] at hb
      exact hb
  next hposition => exact hbound
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

theorem nextCombination_size (indices : Array Nat) (upper : Nat) :
    (Generated.StrictRecombine.nextCombination indices upper).2.size =
      indices.size := by
  unfold Generated.StrictRecombine.nextCombination
  split
  next hfits =>
    split
    next hpivot => rfl
    next pivot hpivot =>
      split
      next hpivotBounds =>
        dsimp
        rw [resetCombinationSuffix_size]
        simp
      next hpivotBounds => rfl
  next hfits => rfl

theorem nextCombinationPivot_some_lt (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
      inspected = some pivot) :
    pivot < indices.size := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          exact ih (indices.size - (inspected + 1)) (by omega)
            (inspected + 1) hresult rfl
        next havailable =>
          have hpivot := Option.some.inj hresult
          subst pivot
          omega
      next hinspected => contradiction

theorem nextCombinationPivot_some_suffix (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
    inspected = some pivot) :
    ∀ position (hpivot : pivot < position) (hposition : position < indices.size)
      (hunscanned : inspected ≤ indices.size - 1 - position),
      indices[position] = upper - indices.size + position := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          intro position hpivot hposition hunscanned
          by_cases hcurrent : position = indices.size - 1 - inspected
          · subst position
            exact hmaximal
          · exact ih (indices.size - (inspected + 1)) (by omega)
              (inspected + 1) hresult rfl position hpivot hposition (by omega)
        next havailable =>
          intro position hpivot hposition hunscanned
          have hpivotValue := Option.some.inj hresult
          subst pivot
          omega
      next hinspected => contradiction

theorem nextCombinationPivot_some_suffix_zero (indices : Array Nat)
    (upper pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper 0 =
      some pivot) :
    ∀ position (hpivot : pivot < position) (hposition : position < indices.size),
      indices[position] = upper - indices.size + position := by
  intro position hpivot hposition
  exact nextCombinationPivot_some_suffix indices upper 0 pivot hresult
    position hpivot hposition (by omega)

theorem nextCombinationPivot_some_ne (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
      inspected = some pivot) :
    indices[pivot]'(nextCombinationPivot_some_lt indices upper inspected pivot
      hresult) ≠ upper - indices.size + pivot := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          exact ih (indices.size - (inspected + 1)) (by omega)
            (inspected + 1) hresult rfl
        next havailable =>
          have hpivotValue := Option.some.inj hresult
          subst pivot
          exact havailable
      next hinspected => contradiction

theorem nextCombination_true_pivot (indices next : Array Nat) (upper : Nat)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    ∃ pivot, ∃ hpivot : pivot < indices.size,
      next[pivot]! = indices[pivot]! + 1 ∧
        ∀ index, index < pivot → next[index]! = indices[index]! := by
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next hfits =>
    split at hrun
    next hpivotNone => simp at hrun
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds =>
        have hpivotLt := nextCombinationPivot_some_lt indices upper 0 pivot
          hpivotSome
        have hout := Prod.mk.inj hrun
        have hnext := hout.2
        subst next
        refine ⟨pivot, hpivotLt, ?_, ?_⟩
        · rw [resetCombinationSuffix_getElem_le _ pivot 0 pivot
            (by simpa using hpivotLt) (Nat.le_refl _)]
          simp [getElem!_pos, hpivotBounds]
        · intro index hindex
          have hindexBounds : index < indices.size := by omega
          rw [resetCombinationSuffix_getElem_le _ pivot 0 index
            (by simpa using hindexBounds) (Nat.le_of_lt hindex)]
          rw [getElem!_pos _ index (by simpa using hindexBounds),
            getElem!_pos indices index hindexBounds]
          rw [Array.getElem_set]
          rw [if_neg (by omega)]
      next hpivotBounds => simp at hrun
  next hfits => simp at hrun

private theorem nextCombination_valid (upper count : Nat)
    (indices next : Array Nat) (hvalid : ValidCombination upper count indices)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    ValidCombination upper count next := by
  rcases hvalid with ⟨hsize, hbound⟩
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next hfits =>
    split at hrun
    next hpivotNone => simp at hrun
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds =>
        have hout := Prod.mk.inj hrun
        rw [← hout.2]
        constructor
        · simpa [resetCombinationSuffix_size, hsize]
        · apply resetCombinationSuffix_preserves_bounds
          intro index hindex
          rw [getElem!_pos _ index hindex, Array.getElem_set]
          split
          next heq =>
            subst index
            have hpivotLt := nextCombinationPivot_some_lt indices upper 0 pivot
              hpivotSome
            have hne := nextCombinationPivot_some_ne indices upper 0 pivot
              hpivotSome
            have hb := hbound pivot hpivotLt
            rw [getElem!_pos indices pivot hpivotLt] at hb
            rw [hsize] at hne
            omega
          next heq =>
            have hi : index < indices.size := by simpa using hindex
            have hb := hbound index hi
            rw [getElem!_pos indices index hi] at hb
            exact hb
      next hpivotBounds => simp at hrun
  next hfits => simp at hrun

private theorem validCombination_digits (upper count : Nat)
    (indices : Array Nat) (hvalid : ValidCombination upper count indices)
    (hfits : count ≤ upper) :
    ∀ digit ∈ indices.toList, digit < upper + 1 := by
  intro digit hdigit
  obtain ⟨index, hindex, hvalue⟩ := List.mem_iff_getElem.mp hdigit
  rcases hvalid with ⟨hsize, hbound⟩
  have hi : index < indices.size := by simpa using hindex
  have hb := hbound index hi
  rw [getElem!_pos indices index hi] at hb
  have harray := Array.getElem_toList (xs := indices) hi
  rw [hvalue] at harray
  rw [← harray] at hb
  omega

private theorem nextCombination_code_increases (upper count : Nat)
    (indices next : Array Nat) (hvalid : ValidCombination upper count indices)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    positionalCode (upper + 1) indices.toList <
      positionalCode (upper + 1) next.toList := by
  have hfits : count ≤ upper := by
    have hcopy := hrun
    unfold Generated.StrictRecombine.nextCombination at hcopy
    split at hcopy
    next hsizeFits => exact hvalid.1 ▸ hsizeFits
    next hsizeFits => simp at hcopy
  have hnextValid := nextCombination_valid upper count indices next hvalid hrun
  obtain ⟨pivot, hpivot, hpivotValue, hprefixValues⟩ :=
    nextCombination_true_pivot indices next upper hrun
  have hnextSize := nextCombination_size indices upper
  rw [hrun] at hnextSize
  have hnextSize' : next.size = indices.size := by simpa using hnextSize
  have hpivotNext : pivot < next.size := by omega
  apply positionalCode_array_pivot_lt (upper + 1) (by omega) indices next pivot
    hpivot hpivotNext hnextSize'.symm
  · apply List.ext_getElem
    · simp [List.length_take]
      omega
    · intro index hleft hright
      have hi : index < pivot := by
        simp [List.length_take] at hleft
        exact hleft.1
      rw [List.getElem_take, List.getElem_take,
        Array.getElem_toList (xs := indices) (by omega),
        Array.getElem_toList (xs := next) (by omega)]
      have heq := hprefixValues index hi
      rw [getElem!_pos indices index (by omega),
        getElem!_pos next index (by omega)] at heq
      exact heq.symm
  · omega
  · exact validCombination_digits upper count indices hvalid hfits
  · exact validCombination_digits upper count next hnextValid hfits

private theorem initialCombination_valid (upper count : Nat)
    (hfits : count ≤ upper) :
    ValidCombination upper count
      (Generated.StrictRecombine.initialCombination count) := by
  constructor
  · exact initialCombination_size count
  · intro index hindex
    have hindexCount : index < count := by
      simpa [initialCombination_size] using hindex
    have harrayList :
        (Generated.StrictRecombine.initialCombination count)[index]! =
          (Generated.StrictRecombine.initialCombination count).toList[index]! := by
      rw [getElem!_pos _ index hindex,
        getElem!_pos _ index (by simpa using hindex)]
      exact (Array.getElem_toList hindex).symm
    have htoList := congrArg (fun values : List Nat => values[index]!)
      (initialCombination_toList count)
    dsimp at htoList
    have hrange : (List.range count)[index]! = index := by
      rw [getElem!_pos _ index (by simpa using hindexCount)]
      simp
    rw [harrayList, htoList, hrange]
    omega

private def concreteCombinationTermination (upper count : Nat) :
    Generated.StrictRecombine.CombinationTermination upper count where
  valid := ValidCombination upper count
  valid_size := fun _ hvalid => hvalid.1
  rank := fun current =>
    (upper + 1) ^ count - positionalCode (upper + 1) current.toList
  next_valid := fun current next hvalid hrun =>
    nextCombination_valid upper count current next hvalid hrun
  next_decreases := by
    intro current next hvalid hrun
    have hfits : count ≤ upper := by
      have hcopy := hrun
      unfold Generated.StrictRecombine.nextCombination at hcopy
      split at hcopy
      next hsizeFits => exact hvalid.1 ▸ hsizeFits
      next hsizeFits => simp at hcopy
    have hnextValid := nextCombination_valid upper count current next hvalid hrun
    have hcode := nextCombination_code_increases upper count current next hvalid hrun
    apply Nat.sub_lt_sub_left
    · have hdigits := validCombination_digits upper count current hvalid hfits
      simpa [hvalid.1] using
        positionalCode_lt_pow (upper + 1) (by omega) current.toList hdigits
    · exact hcode

theorem removeCombinationLoop_size (candidate : Array Nat)
    (remaining : Nat) (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombinationLoop candidate remaining
      active = .ok output) :
    output.size + remaining = active.size := by
  induction remaining generalizing active output with
  | zero =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      have hout := Except.ok.inj hrun
      subst output
      simp
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      split at hrun
      next hcand =>
        dsimp at hrun
        split at hrun
        next hactive =>
          have htail := ih (active.eraseIdxIfInBounds candidate[remaining])
            output hrun
          simp only [Array.size_eraseIdxIfInBounds, if_pos hactive] at htail
          omega
        next hactive => contradiction
      next hcand => contradiction

theorem removeCombination_strict (candidate : Array Nat)
    (active output : Array SparsePolyZZ) (hnonempty : 0 < candidate.size)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    output.size < active.size := by
  unfold Generated.StrictRecombine.removeCombination at hrun
  have hsize := removeCombinationLoop_size candidate candidate.size active
    output hrun
  omega

/-- The concrete termination package for the source Zassenhaus loops.  Its
combination rank is the complement of the base-`upper+1` positional code;
the outer rank is discharged by the actual successful subset removal. -/
def concreteZassenhausTermination :
    Generated.StrictRecombine.ZassenhausTermination where
  combinations := concreteCombinationTermination
  initial_valid := fun upper count hfits =>
    initialCombination_valid upper count hfits
  removal_decreases := fun active candidate output hnonempty hrun =>
    removeCombination_strict candidate active output hnonempty hrun

theorem sizeReduceAt_preserves_norms_k
    (state output : Generated.StrictRecombine.LLLState) (j : Nat)
    (hrun : Generated.StrictRecombine.sizeReduceAt state j = .ok output) :
    output.norms = state.norms ∧ output.k = state.k := by
  unfold Generated.StrictRecombine.sizeReduceAt at hrun
  split at hrun
  next hk =>
    split at hrun
    next hj =>
      cases hround : Generated.StrictRecombine.roundQQ state.mu[state.k][j] with
      | error fault => simp [hround] at hrun
      | ok q =>
        simp only [hround] at hrun
        split at hrun
        next hzero =>
          have hout := Except.ok.inj hrun
          subst output
          exact ⟨rfl, rfl⟩
        next hnonzero =>
          cases hsubtract : Generated.StrictRecombine.subtractMatrixRows
              state.matrix state.transform state.k j q with
          | error fault => simp [hsubtract] at hrun
          | ok matrices =>
            rcases matrices with ⟨matrix', transform'⟩
            simp only [hsubtract] at hrun
            cases hmu : Generated.StrictRecombine.reduceMuPrefixLoop
                state.mu state.k j q j 0 with
            | error fault => simp [hmu] at hrun
            | ok mu' =>
              simp only [hmu] at hrun
              split at hrun
              next hk' =>
                split at hrun
                next hj' =>
                  have hout := Except.ok.inj hrun
                  subst output
                  exact ⟨rfl, rfl⟩
                next hj' => contradiction
              next hk' => contradiction
    next hj => contradiction
  next hk => contradiction

theorem extraSizeReduceLoop_preserves_norms_k
    (remaining : Nat) (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.extraSizeReduceLoop remaining state =
      .ok output) :
    output.norms = state.norms ∧ output.k = state.k := by
  induction remaining generalizing state output with
  | zero =>
      simp [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      subst output
      exact ⟨rfl, rfl⟩
  | succ remaining ih =>
      rw [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      cases hstep : Generated.StrictRecombine.sizeReduceAt state remaining with
      | error fault => simp [hstep] at hrun
      | ok next =>
        simp only [hstep] at hrun
        have hhead := sizeReduceAt_preserves_norms_k state next remaining hstep
        have htail := ih next output hrun
        exact ⟨htail.1.trans hhead.1, htail.2.trans hhead.2⟩

theorem lllStep_advanced_control
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    output.norms = state.norms ∧ output.k = state.k + 1 := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hreduceControl := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        split at hrun
        next hkNorm =>
          split at hrun
          next hpredNorm =>
            split at hrun
            next hkMu =>
              split at hrun
              next hpredMu =>
                dsimp at hrun
                split at hrun
                next hlovasz =>
                  cases hextra : Generated.StrictRecombine.extraSizeReduceLoop
                      (reduced.k - 1) reduced with
                  | error fault => simp [hextra] at hrun
                  | ok fullyReduced =>
                    simp only [hextra] at hrun
                    have hextraControl :=
                      extraSizeReduceLoop_preserves_norms_k
                        (reduced.k - 1) reduced fullyReduced hextra
                    have hout := Except.ok.inj hrun
                    injection hout with hstate
                    subst output
                    constructor
                    · exact hextraControl.1.trans hreduceControl.1
                    · simp only
                      rw [hextraControl.2, hreduceControl.2]
                next hlovasz =>
                  repeat' first | split at hrun | simp_all
              next hpredMu => contradiction
            next hkMu => contradiction
          next hpredNorm => contradiction
        next hkNorm => contradiction
    · rw [dif_neg hkMatrix] at hrun
      simp [hkMatrix] at hrun
  · rw [dif_neg hkPositive] at hrun
    simp [hkPositive] at hrun

theorem lllStep_swapped_k
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    output.k = Nat.max (state.k - 1) 1 := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun; rfl
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_swapped_norms
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    ∃ reduced,
      Generated.StrictRecombine.sizeReduceAt state (state.k - 1) = .ok reduced ∧
      reduced.norms = state.norms ∧ reduced.k = state.k ∧
      output.norms = Generated.StrictRecombine.normsAfterLovaszSwap
        reduced.norms reduced.k
          ((reduced.mu[reduced.k]!)[reduced.k - 1]!) := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun
        all_goals
          have hkMuState : state.k < reduced.mu.size := by
            rw [← hcontrol.2]
            assumption
          have hpredMuState : state.k - 1 < reduced.mu[state.k].size := by
            simpa only [hcontrol.2] using
              (‹reduced.k - 1 < reduced.mu[reduced.k].size›)
          change Generated.StrictRecombine.normsAfterLovaszSwap state.norms
              state.k ((reduced.mu[state.k]'hkMuState)[state.k - 1]'hpredMuState) =
            Generated.StrictRecombine.normsAfterLovaszSwap state.norms
              state.k ((reduced.mu[state.k]!)[state.k - 1]!)
          congr 1
          rw [getElem!_pos reduced.mu state.k hkMuState]
          rw [getElem!_pos (reduced.mu[state.k]) (state.k - 1) hpredMuState]
    · rw [dif_neg hkMatrix] at hrun
      simp [hkMatrix] at hrun
  · rw [dif_neg hkPositive] at hrun
    simp [hkPositive] at hrun

/-- The classical LLL potential, evaluated on the exact Gram–Schmidt norm
array carried by the generated C++ state. -/
noncomputable def arrayLLLPotential (norms : Array QQ) : QQ :=
  ∏ i : Fin norms.size, norms[i] ^ (norms.size - i.val)

theorem arrayLLLPotential_pos (norms : Array QQ)
    (hpositive : ∀ i (hi : i < norms.size), 0 < norms[i]'hi) :
    0 < arrayLLLPotential norms := by
  unfold arrayLLLPotential
  apply Finset.prod_pos
  intro i _
  exact pow_pos (hpositive i i.isLt) _

/-- View an operational array state as the finite-function state used by the
pure L2 LLL potential theorem.  Matrix and μ reads are total only because this
view is used for the termination algebra; the execution theorem separately
establishes every source-level bound before a read. -/
noncomputable def toPotentialState
    (state : Generated.StrictRecombine.LLLState) :
    LLLState state.norms.size state.matrix.size where
  basis := fun i j => (state.matrix[i.val]!)[j.val]!
  mu := fun i j => (state.mu[i.val]!)[j.val]!
  gs_norm_sq := fun i => state.norms[i]
  k := state.k

theorem toPotentialState_potential
    (state : Generated.StrictRecombine.LLLState) :
    lllPotential' (toPotentialState state) = arrayLLLPotential state.norms := by
  rfl

theorem normsAfterLovaszSwap_size (norms : Array QQ) (k : Nat) (mu : QQ) :
    (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu).size =
      norms.size := by
  unfold Generated.StrictRecombine.normsAfterLovaszSwap
  dsimp
  split <;> simp

theorem normsAfterLovaszSwap_get (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hnew : norms[k] + mu * mu * norms[k - 1] ≠ 0)
    (i : Nat) (hi : i < norms.size) :
    (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[i]! =
      if k - 1 = i then norms[k] + mu * mu * norms[k - 1]
      else if k = i then
        norms[k] * norms[k - 1] /
          (norms[k] + mu * mu * norms[k - 1])
      else norms[i] := by
  unfold Generated.StrictRecombine.normsAfterLovaszSwap
  dsimp
  rw [getElem!_pos norms k hk, getElem!_pos norms (k - 1) hpred]
  rw [if_pos hnew]
  rw [getElem!_pos _ i (by simpa using hi)]
  rw [Array.getElem_setIfInBounds (by simpa using hi)]
  rw [Array.getElem_setIfInBounds hi]

/-- Once the determinant potential is converted to a natural number, this is
the concrete scalar rank used by the generated while-loop: an advancing step
spends the low-order progress component, while a swap must lower the dominant
lattice component. -/
def lllLexRank (latticePotential dimension index : Nat) : Nat :=
  latticePotential * (dimension + 1) + (dimension - index)

theorem lllLexRank_advanced (potential dimension index : Nat)
    (hindex : index < dimension) :
    lllLexRank potential dimension (index + 1) <
      lllLexRank potential dimension index := by
  unfold lllLexRank
  omega

theorem lllLexRank_swap (oldPotential newPotential dimension oldIndex newIndex : Nat)
    (hpotential : newPotential < oldPotential) :
    lllLexRank newPotential dimension newIndex <
      lllLexRank oldPotential dimension oldIndex := by
  unfold lllLexRank
  calc
    newPotential * (dimension + 1) + (dimension - newIndex) ≤
        newPotential * (dimension + 1) + dimension :=
      Nat.add_le_add_left (Nat.sub_le dimension newIndex) _
    _ < newPotential * (dimension + 1) + (dimension + 1) := by omega
    _ = (newPotential + 1) * (dimension + 1) := by
      rw [Nat.add_mul, one_mul]
    _ ≤ oldPotential * (dimension + 1) :=
      Nat.mul_le_mul_right _ (Nat.succ_le_iff.mpr hpotential)
    _ ≤ oldPotential * (dimension + 1) + (dimension - oldIndex) :=
      Nat.le_add_right _ _

/-- Integer Gram matrix of the first `prefix` rows of the current C++ lattice
basis.  Total reads make the definition available on every raw state; the
operational validity invariant later states that the basis is square. -/
noncomputable def gramPrefixMatrix
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin rowCount) Int :=
  fun i j => ∑ c : Fin matrix.size,
    (matrix[i.val]!)[c.val]! * (matrix[j.val]!)[c.val]!

/-- Each prefix Gram determinant is a nonnegative integer.  Their product is
the discrete version of the rational LLL potential
`∏ᵢ Bᵢ^(n-i)`: Gram–Schmidt identifies the determinant of prefix `r`
with `∏ i<r, Bᵢ`. -/
noncomputable def lllDeterminantPotential
    (matrix : Generated.StrictRecombine.LLLMatrix) : Nat :=
  ∏ i : Fin matrix.size, (Matrix.det (gramPrefixMatrix matrix (i.val + 1))).natAbs

noncomputable def concreteLLLRank
    (state : Generated.StrictRecombine.LLLState) : Nat :=
  lllLexRank (lllDeterminantPotential state.matrix) state.matrix.size state.k

noncomputable def prefixNormProduct (norms : Array QQ) (rowCount : Nat) : QQ :=
  ∏ i : Fin rowCount, norms[i.val]!

/-- Concrete semantic invariant carried by the generated LLL execution.  It
states square integer storage, positive Gram–Schmidt norms, and the exact
Gram determinant identity for every prefix.  No result or termination trace
is stored here. -/
structure ConcreteLLLValid
    (state : Generated.StrictRecombine.LLLState) : Prop where
  norms_size : state.norms.size = state.matrix.size
  rows_square : ∀ row (hrow : row < state.matrix.size),
    state.matrix[row].size = state.matrix.size
  norms_positive : ∀ index (hindex : index < state.norms.size),
    0 < state.norms[index]
  gram_prefix : ∀ rowCount, rowCount ≤ state.matrix.size →
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
      prefixNormProduct state.norms rowCount

theorem prefixNormProduct_pos
    (norms : Array QQ) (rowCount : Nat)
    (hrowCount : rowCount ≤ norms.size)
    (hpositive : ∀ index (hindex : index < norms.size), 0 < norms[index]) :
    0 < prefixNormProduct norms rowCount := by
  unfold prefixNormProduct
  apply Finset.prod_pos
  intro i _
  rw [getElem!_pos norms i.val (lt_of_lt_of_le i.isLt hrowCount)]
  exact hpositive i.val (lt_of_lt_of_le i.isLt hrowCount)

theorem prefixNormProduct_succ (norms : Array QQ) (rowCount : Nat) :
    prefixNormProduct norms (rowCount + 1) =
      prefixNormProduct norms rowCount * norms[rowCount]! := by
  unfold prefixNormProduct
  rw [Fin.prod_univ_castSucc]
  rfl

theorem prefixNormProduct_normsAfterLovaszSwap
    (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositive : 0 < k)
    (hnew : norms[k] + mu * mu * norms[k - 1] ≠ 0) :
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k =
      prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) := by
  have hkSplit : k - 1 + 1 = k := Nat.sub_add_cancel (by omega)
  calc
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k =
        prefixNormProduct
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) (k - 1) *
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[k - 1]! := by
      simpa only [hkSplit] using
        (prefixNormProduct_succ
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) (k - 1))
    _ = prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) := by
      congr 1
      · unfold prefixNormProduct
        apply Finset.prod_congr rfl
        intro i _
        rw [normsAfterLovaszSwap_get norms k mu hk hpred hnew
          i.val (lt_trans i.isLt hpred)]
        have hpredNe : k - 1 ≠ i.val := by omega
        have hkNe : k ≠ i.val := by omega
        simp [hpredNe, hkNe, getElem!_pos norms i.val (lt_trans i.isLt hpred)]
      · rw [normsAfterLovaszSwap_get norms k mu hk hpred hnew
          (k - 1) hpred]
        simp

theorem prefixNormProduct_normsAfterLovaszSwap_lt
    (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositiveK : 0 < k)
    (hnormsPositive : ∀ index (hindex : index < norms.size), 0 < norms[index])
    (hfail : norms[k] < ((3 : QQ) / 4 - mu * mu) * norms[k - 1]) :
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k <
      prefixNormProduct norms k := by
  have hkNormPos := hnormsPositive k hk
  have hpredNormPos := hnormsPositive (k - 1) hpred
  have hnewPos : 0 < norms[k] + mu * mu * norms[k - 1] := by
    nlinarith [sq_nonneg mu]
  have hnewLt : norms[k] + mu * mu * norms[k - 1] < norms[k - 1] := by
    nlinarith
  have hprefixPos : 0 < prefixNormProduct norms (k - 1) :=
    prefixNormProduct_pos norms (k - 1) (Nat.le_of_lt hpred) hnormsPositive
  have hkSplit : k - 1 + 1 = k := Nat.sub_add_cancel (by omega)
  rw [prefixNormProduct_normsAfterLovaszSwap norms k mu hk hpred hpositiveK
    (ne_of_gt hnewPos)]
  calc
    prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) <
        prefixNormProduct norms (k - 1) * norms[k - 1] :=
      mul_lt_mul_of_pos_left hnewLt hprefixPos
    _ = prefixNormProduct norms k := by
      simpa only [hkSplit, getElem!_pos norms (k - 1) hpred] using
        (prefixNormProduct_succ norms (k - 1)).symm

theorem lllStep_swapped_prefixNormProduct_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    prefixNormProduct output.norms state.k <
      prefixNormProduct state.norms state.k := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun
        all_goals
          have hreducedPositive : ∀ index (hindex : index < reduced.norms.size),
              0 < reduced.norms[index] := by
            intro index hindex
            have hindexState : index < state.norms.size := by
              simpa only [hcontrol.1] using hindex
            have hpositiveState := hvalid.norms_positive index hindexState
            simpa only [← hcontrol.1] using hpositiveState
          have hkStateNorm : state.k < state.norms.size := by
            simpa only [← hcontrol.1, ← hcontrol.2] using
              (‹reduced.k < reduced.norms.size›)
          have hpredStateNorm : state.k - 1 < state.norms.size := by
            simpa only [← hcontrol.1, ← hcontrol.2] using
              (‹reduced.k - 1 < reduced.norms.size›)
          have hkStateMu : state.k < reduced.mu.size := by
            simpa only [← hcontrol.2] using (‹reduced.k < reduced.mu.size›)
          have hpredStateMu : state.k - 1 < reduced.mu[state.k].size := by
            simpa only [← hcontrol.2] using
              (‹reduced.k - 1 < reduced.mu[reduced.k].size›)
          have hfail : reduced.norms[reduced.k] <
              ((3 : QQ) / 4 -
                reduced.mu[reduced.k][reduced.k - 1] *
                  reduced.mu[reduced.k][reduced.k - 1]) *
                reduced.norms[reduced.k - 1] := by
            have hfailState :
                state.norms[state.k]'hkStateNorm <
                  ((3 : QQ) / 4 -
                    (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu *
                      (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu) *
                    state.norms[state.k - 1]'hpredStateNorm := by assumption
            simpa only [← hcontrol.1, ← hcontrol.2] using hfailState
          have hdecrease := prefixNormProduct_normsAfterLovaszSwap_lt
            reduced.norms reduced.k reduced.mu[reduced.k][reduced.k - 1]
            (by assumption) (by assumption) (by simpa [hcontrol.2] using hkPositive)
            hreducedPositive hfail
          simpa only [hcontrol.1, hcontrol.2] using hdecrease
    · rw [dif_neg hkMatrix] at hrun
      simp at hrun
  · rw [dif_neg hkPositive] at hrun
    simp at hrun

theorem lllStep_swapped_source_index_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    state.k < state.matrix.size := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · exact hkMatrix
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_swapped_matrix_size
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    output.matrix.size = state.matrix.size := by
  rcases lllStep_swapped_norms state output hrun with
    ⟨reduced, hreduce, hnorms, hk, hout⟩
  calc
    output.matrix.size = output.norms.size := houtputValid.norms_size.symm
    _ = reduced.norms.size := by
      rw [hout, normsAfterLovaszSwap_size]
    _ = state.norms.size := congrArg Array.size hnorms
    _ = state.matrix.size := hvalid.norms_size

theorem lllStep_swapped_boundaryGramDet_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    (Matrix.det (gramPrefixMatrix output.matrix state.k)).natAbs <
      (Matrix.det (gramPrefixMatrix state.matrix state.k)).natAbs := by
  have hkMatrix := lllStep_swapped_source_index_lt state output hrun
  have hsize := lllStep_swapped_matrix_size state output hvalid houtputValid hrun
  have hkOutput : state.k ≤ output.matrix.size := by omega
  have hprefix := lllStep_swapped_prefixNormProduct_lt state output hvalid hrun
  have houtGram := houtputValid.gram_prefix state.k hkOutput
  have hinGram := hvalid.gram_prefix state.k (Nat.le_of_lt hkMatrix)
  rw [← houtGram, ← hinGram] at hprefix
  have hdet : Matrix.det (gramPrefixMatrix output.matrix state.k) <
      Matrix.det (gramPrefixMatrix state.matrix state.k) := by
    exact_mod_cast hprefix
  have houtProductPos := prefixNormProduct_pos output.norms state.k
    (by rw [houtputValid.norms_size]; exact hkOutput)
    houtputValid.norms_positive
  have houtDetPos : 0 < Matrix.det
      (gramPrefixMatrix output.matrix state.k) := by
    exact_mod_cast (show (0 : QQ) <
      ((Matrix.det (gramPrefixMatrix output.matrix state.k) : Int) : QQ) by
        rw [houtGram]
        exact houtProductPos)
  exact Int.natAbs_lt_natAbs_of_nonneg_of_lt
    (Int.le_of_lt houtDetPos) hdet

theorem ConcreteLLLValid.gram_prefix_pos
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) (rowCount : Nat)
    (hrowCount : rowCount ≤ state.matrix.size) :
    0 < Matrix.det (gramPrefixMatrix state.matrix rowCount) := by
  have hproduct := prefixNormProduct_pos state.norms rowCount
    (by rw [hvalid.norms_size]; exact hrowCount) hvalid.norms_positive
  have heq := hvalid.gram_prefix rowCount hrowCount
  exact_mod_cast (show (0 : QQ) <
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) by
      rw [heq]
      exact hproduct)

theorem ConcreteLLLValid.gram_prefix_natAbs
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) (rowCount : Nat)
    (hrowCount : rowCount ≤ state.matrix.size) :
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount)).natAbs : Int) =
      Matrix.det (gramPrefixMatrix state.matrix rowCount) :=
  Int.natAbs_of_nonneg
    (Int.le_of_lt (hvalid.gram_prefix_pos rowCount hrowCount))

theorem swapMatrixRows_ok
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output) :
    ∃ (hleft : left < matrix.size) (hright : right < matrix.size),
      output = (matrix.set left matrix[right]).set right matrix[left] (by
        simpa) := by
  unfold Generated.StrictRecombine.swapMatrixRows at hrun
  split at hrun
  next hleft =>
    split at hrun
    next hright =>
      exact ⟨hleft, hright, (Except.ok.inj hrun).symm⟩
    next hright => contradiction
  next hleft => contradiction

theorem swapMatrixRows_size
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output) :
    output.size = matrix.size := by
  rcases swapMatrixRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp

theorem swapMatrixRows_get
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right row : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrow : row < matrix.size) :
    output[row]'(by simpa [swapMatrixRows_size matrix output left right hrun]) =
      if right = row then matrix[left]!
      else if left = row then matrix[right]!
      else matrix[row] := by
  rcases swapMatrixRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp only [getElem!_pos matrix left hleft,
    getElem!_pos matrix right hright, Array.getElem_set]

theorem swapMatrixRows_get_fin
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) (row : Fin rowCount) :
    output[row.val]! = matrix[((Equiv.swap
      (⟨left, hleftCount⟩ : Fin rowCount)
      (⟨right, hrightCount⟩ : Fin rowCount)) row).val]! := by
  have hrow : row.val < matrix.size := lt_of_lt_of_le row.isLt hrowCount
  have houtSize := swapMatrixRows_size matrix output left right hrun
  rw [getElem!_pos output row.val (by simpa [houtSize])]
  rw [swapMatrixRows_get matrix output left right row hrun hrow]
  by_cases hrightRow : row = ⟨right, hrightCount⟩
  · subst row
    simp only [Equiv.swap_apply_right]
    rw [getElem!_pos matrix left (lt_of_lt_of_le hleftCount hrowCount)]
    simp
  · by_cases hleftRow : row = ⟨left, hleftCount⟩
    · subst row
      simp only [Equiv.swap_apply_left]
      rw [getElem!_pos matrix right (lt_of_lt_of_le hrightCount hrowCount)]
      have hrightLeft : right ≠ left := by
        intro heq
        apply hrightRow
        exact Fin.ext heq.symm
      simp [hrightLeft]
    · have hrightVal : right ≠ row.val := by
        intro heq
        apply hrightRow
        exact Fin.ext heq.symm
      have hleftVal : left ≠ row.val := by
        intro heq
        apply hleftRow
        exact Fin.ext heq.symm
      rw [Equiv.swap_apply_of_ne_of_ne hleftRow hrightRow]
      rw [getElem!_pos matrix row.val hrow]
      simp [hrightVal, hleftVal]

theorem gramPrefixMatrix_swap_reindex
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    gramPrefixMatrix output rowCount = Matrix.reindex
      (Equiv.swap (⟨left, hleftCount⟩ : Fin rowCount)
        (⟨right, hrightCount⟩ : Fin rowCount))
      (Equiv.swap (⟨left, hleftCount⟩ : Fin rowCount)
        (⟨right, hrightCount⟩ : Fin rowCount))
      (gramPrefixMatrix matrix rowCount) := by
  funext i j
  unfold gramPrefixMatrix Matrix.reindex
  simp only [Equiv.coe_fn_mk]
  have hsize := swapMatrixRows_size matrix output left right hrun
  rw [hsize]
  apply Finset.sum_congr rfl
  intro c _
  rw [swapMatrixRows_get_fin matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount i]
  rw [swapMatrixRows_get_fin matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount j]
  simp

theorem gramPrefixDet_swap_preserved_of_both
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    Matrix.det (gramPrefixMatrix output rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  rw [gramPrefixMatrix_swap_reindex matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount]
  exact Matrix.det_reindex_self _ _

theorem swapMatrixRows_get_before
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrowCount : rowCount ≤ matrix.size)
    (hbeforeLeft : rowCount ≤ left) (hbeforeRight : rowCount ≤ right)
    (row : Fin rowCount) :
    output[row.val]! = matrix[row.val]! := by
  have hrow : row.val < matrix.size := lt_of_lt_of_le row.isLt hrowCount
  have houtSize := swapMatrixRows_size matrix output left right hrun
  rw [getElem!_pos output row.val (by simpa [houtSize])]
  rw [swapMatrixRows_get matrix output left right row hrun hrow]
  have hleft : left ≠ row.val := by omega
  have hright : right ≠ row.val := by omega
  simp [hleft, hright, getElem!_pos matrix row.val hrow]

theorem gramPrefixDet_swap_preserved_of_before
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrowCount : rowCount ≤ matrix.size)
    (hbeforeLeft : rowCount ≤ left) (hbeforeRight : rowCount ≤ right) :
    Matrix.det (gramPrefixMatrix output rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  congr 1
  funext i j
  unfold gramPrefixMatrix
  have hsize := swapMatrixRows_size matrix output left right hrun
  rw [hsize]
  apply Finset.sum_congr rfl
  intro c _
  rw [swapMatrixRows_get_before matrix output left right rowCount hrun
    hrowCount hbeforeLeft hbeforeRight i]
  rw [swapMatrixRows_get_before matrix output left right rowCount hrun
    hrowCount hbeforeLeft hbeforeRight j]

noncomputable def factorArrayToL2 (factors : Array SparsePolyZZ) :
    List (Polynomial Int) :=
  factors.toList.map SparsePolyZZ.toPoly

def ValidActiveIndices (active : Array Int32)
    (lifted : Array SparsePolyZZ) : Prop :=
  ∀ index (hindex : index < active.size),
    0 ≤ active[index]'hindex ∧
      (active[index]'hindex).toInt64.toNat < lifted.size

private def gatherSuffix (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat) : Array SparsePolyZZ :=
  ((active.toList.drop index).map
    (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray

theorem gatherActiveLoop_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat)
    (result : Array SparsePolyZZ)
    (hvalid : ValidActiveIndices active lifted) :
    Generated.StrictRecombine.gatherActiveLoop active lifted index result =
      .ok (result ++ gatherSuffix active lifted index) := by
  induction hmeasure : active.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.gatherActiveLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        rw [ih (active.size - (index + 1)) (by omega)
          (index + 1) (result.push lifted[active[index].toInt64.toNat]) rfl]
        have hsuffix : active.toList.drop index = active[index] ::
            active.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := active.toList) (i := index) (by simpa using hindex)
        simp [gatherSuffix, hsuffix, Array.push, hentry.2]
      next hindex =>
        have hle : active.size ≤ index := Nat.le_of_not_gt hindex
        simp [gatherSuffix, List.drop_eq_nil_iff.mpr hle]

theorem gatherActive_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (hvalid : ValidActiveIndices active lifted) :
    ∃ output,
      Generated.StrictRecombine.gatherActive active lifted = .ok output ∧
      factorArrayToL2 output = active.toList.map
        (fun sourceIndex => SparsePolyZZ.toPoly
          lifted[sourceIndex.toInt64.toNat]!) := by
  refine ⟨((active.toList.map
      (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray), ?_, ?_⟩
  · simpa [Generated.StrictRecombine.gatherActive, gatherSuffix] using
      gatherActiveLoop_refines active lifted 0 #[] hvalid
  · simp [factorArrayToL2]

theorem gatherActive_size_of_success (active : Array Int32)
    (lifted output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.gatherActive active lifted = .ok output) :
    output.size = active.size := by
  unfold Generated.StrictRecombine.gatherActive at hrun
  have loopSize : ∀ index result,
      Generated.StrictRecombine.gatherActiveLoop active lifted index result =
        .ok output → output.size = result.size + (active.size - index) := by
    intro index result hloop
    induction hmeasure : active.size - index using Nat.strong_induction_on
        generalizing index result with
    | h measure ih =>
        rw [Generated.StrictRecombine.gatherActiveLoop] at hloop
        split at hloop
        next hindex =>
          dsimp at hloop
          split at hloop
          next hnonnegative =>
            split at hloop
            next hsource =>
              have htail := ih (active.size - (index + 1)) (by omega)
                (index + 1) (result.push lifted[active[index].toInt64.toNat])
                hloop rfl
              simp only [Array.size_push] at htail
              omega
            next hsource => simp at hloop
          next hnonnegative => simp at hloop
        next hindex =>
          have hle : active.size ≤ index := Nat.le_of_not_gt hindex
          have hout := Except.ok.inj hloop
          subst output
          have : measure = 0 := by omega
          simp [this]
  simpa using (loopSize 0 #[] hrun)

theorem appendFallbackLoop_refines (fallback : Array SparsePolyZZ)
    (index : Nat) (result : Array SparsePolyZZ) :
    Generated.StrictRecombine.appendFallbackLoop fallback index result =
      .ok (result ++ (fallback.toList.drop index).toArray) := by
  induction hmeasure : fallback.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.appendFallbackLoop]
      split
      next hindex =>
        rw [ih (fallback.size - (index + 1)) (by omega)
          (index + 1) (result.push fallback[index]) rfl]
        have hsuffix : fallback.toList.drop index = fallback[index] ::
            fallback.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := fallback.toList) (i := index) (by simpa using hindex)
        simp [hsuffix, Array.push]
      next hindex =>
        have hle : fallback.size ≤ index := Nat.le_of_not_gt hindex
        simp [List.drop_eq_nil_iff.mpr hle]

theorem appendFallback_refines (fallback result : Array SparsePolyZZ) :
    ∃ output,
      Generated.StrictRecombine.appendFallback fallback result = .ok output ∧
      factorArrayToL2 output =
        factorArrayToL2 result ++ factorArrayToL2 fallback := by
  refine ⟨result ++ fallback, ?_, ?_⟩
  · simpa [Generated.StrictRecombine.appendFallback] using
      appendFallbackLoop_refines fallback 0 result
  · simp [factorArrayToL2]

/-- Pure list meaning of the source reverse-erasure loop. -/
def removeConsumedL2 (active : Array Int32) (consumed : Array Bool) :
    List Int32 :=
  (active.toList.zip consumed.toList).filterMap fun item =>
    if item.2 then none else some item.1

private theorem removeConsumedLoop_refines_prefix
    (consumed : Array Bool) (remaining : Nat) (active : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
        .ok output ∧
      output.size ≤ active.size := by
  induction remaining generalizing active with
  | zero =>
      refine ⟨active, ?_, Nat.le_refl _⟩
      rw [Generated.StrictRecombine.removeConsumedLoop]
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop]
      have hindex : remaining < consumed.size := by omega
      rw [dif_pos hindex]
      split
      next hmarked =>
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive]
        have heraseSize : (active.eraseIdxIfInBounds remaining).size =
            active.size - 1 := by
          simp [hactive]
        have hremaining' : remaining ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rw [heraseSize]
          omega
        rcases ih (active.eraseIdxIfInBounds remaining) (by omega) hremaining' with
          ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize.trans (by rw [heraseSize]; omega)⟩
      next hkept =>
        rcases ih active (by omega) (by omega) with ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize⟩

private theorem removeConsumedLoop_strict_of_marked
    (consumed : Array Bool) (remaining : Nat) (active output : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size)
    (hmarked : ∃ index, ∃ hindex : index < remaining,
      consumed[index] = true)
    (hrun : Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
      .ok output) :
    output.size < active.size := by
  induction remaining generalizing active output with
  | zero =>
      rcases hmarked with ⟨index, hindex, _⟩
      omega
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop] at hrun
      have hindexBound : remaining < consumed.size := by omega
      rw [dif_pos hindexBound] at hrun
      by_cases hlast : consumed[remaining] = true
      · rw [if_pos hlast] at hrun
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive] at hrun
        have htailLe : output.size ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rcases removeConsumedLoop_refines_prefix consumed remaining
              (active.eraseIdxIfInBounds remaining) (by omega) (by
                simp [hactive]; omega) with ⟨tail, htailRun, htailSize⟩
          rw [hrun] at htailRun
          exact Except.ok.inj htailRun ▸ htailSize
        simpa [hactive] using lt_of_le_of_lt htailLe (by
          simp [hactive]
          omega)
      · rw [if_neg hlast] at hrun
        have hmarkedPrefix : ∃ index, ∃ hindex : index < remaining,
            consumed[index] = true := by
          rcases hmarked with ⟨index, hindex, hvalue⟩
          refine ⟨index, ?_, hvalue⟩
          have hne : index ≠ remaining := by
            intro heq
            subst index
            exact hlast hvalue
          omega
        exact ih active output (by omega) (by omega) hmarkedPrefix hrun

theorem removeConsumed_succeeds (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumed active consumed = .ok output ∧
      output.size ≤ active.size := by
  unfold Generated.StrictRecombine.removeConsumed
  rw [dif_pos hsizes]
  exact removeConsumedLoop_refines_prefix consumed active.size active
    (by omega) (Nat.le_refl _)

theorem removeConsumed_strict_of_marked (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size)
    (hmarked : ∃ index, ∃ hindex : index < consumed.size,
      consumed[index] = true)
    (output : Array Int32)
    (hrun : Generated.StrictRecombine.removeConsumed active consumed =
      .ok output) :
    output.size < active.size := by
  unfold Generated.StrictRecombine.removeConsumed at hrun
  rw [dif_pos hsizes] at hrun
  apply removeConsumedLoop_strict_of_marked consumed active.size active output
    (by omega) (Nat.le_refl _) (by simpa [hsizes] using hmarked) hrun

def CandidateIndicesValid (candidate : Array Int32)
    (consumed : Array Bool) : Prop :=
  ∀ index (hindex : index < candidate.size),
    0 ≤ candidate[index] ∧
      candidate[index].toInt64.toNat < consumed.size

theorem candidateAvailableLoop_succeeds (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ available,
      Generated.StrictRecombine.candidateAvailableLoop candidate consumed index =
        .ok available := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index with
  | h measure ih =>
      rw [Generated.StrictRecombine.candidateAvailableLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        split
        next hconsumed => exact ⟨false, rfl⟩
        next hfree =>
          exact ih (candidate.size - (index + 1)) (by omega)
            (index + 1) rfl
      next hindex => exact ⟨true, rfl⟩

theorem markConsumedLoop_succeeds_size (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ output,
      Generated.StrictRecombine.markConsumedLoop candidate index consumed =
        .ok output ∧ output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        let activeNat := candidate[index].toInt64.toNat
        rw [dif_pos hentry.1, dif_pos hentry.2]
        have hvalidSet : CandidateIndicesValid candidate
            (consumed.set activeNat true) := by
          intro next hnext
          have hv := hvalid next hnext
          exact ⟨hv.1, by simpa using hv.2⟩
        rcases ih (candidate.size - (index + 1)) (by omega)
            (consumed.set activeNat true) (index + 1) hvalidSet rfl with
          ⟨output, hrun, hsize⟩
        refine ⟨output, hrun, ?_⟩
        simpa using hsize
      next hindex => exact ⟨consumed, rfl, rfl⟩

theorem markConsumedLoop_size_of_success (candidate : Array Int32)
    (index : Nat) (consumed output : Array Bool)
    (hrun : Generated.StrictRecombine.markConsumedLoop candidate index consumed =
      .ok output) :
    output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed output with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        split at hrun
        next hnonnegative =>
          split at hrun
          next hactive =>
            have htail := ih (candidate.size - (index + 1)) (by omega)
              (index + 1)
              (consumed.set candidate[index].toInt64.toNat true) output
              hrun rfl
            simpa using htail
          next hactive => contradiction
        next hnonnegative => contradiction
      next hindex => exact Except.ok.inj hrun ▸ rfl

theorem validateCandidatesLoop_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    consumed'.size = consumed.size := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have htail := ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega) (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining - candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hrun rfl
                                exact htail.trans
                                  (markConsumedLoop_size_of_success
                                    candidates[candidateIndex] 0 consumed
                                    consumedNext hmark)
                          · simp only [hremainder, if_false] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        exact (congrArg (fun output => output.2.2.size)
          (Except.ok.inj hrun)).symm

theorem validateCandidates_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    consumed.size = activeLifted.size := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact (validateCandidatesLoop_consumed_size ops candidates 0 activeLifted
    modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size hrun).trans
      (by simp)

/-- The concrete gather and validation loops themselves prove that every
successful extraction removes at least one active entry.  No external length
or semantic termination oracle is required. -/
def removalTermination (ops : Generated.StrictRecombine.VanHoeijRawOps) :
    Generated.StrictRecombine.VanHoeijTermination ops := {
  extraction_decreases := by
    intro lifted modulus state activeLifted candidates fStar' result' consumed
      active' hgather hvalidate hfound hremove
    exact removeConsumed_strict_of_marked state.active consumed
      ((validateCandidates_consumed_size ops.validation state.fStar activeLifted
        modulus candidates state.result fStar' result' consumed hvalidate).trans
        (gatherActive_size_of_success state.active lifted activeLifted hgather))
      hfound active' hremove }

/-- Algebraic execution contract for one actual multiply/normalize/mod step.
It relates concrete successful output modulo `m`; it does not supply or choose
a factorization. -/
def TrialProductStepCorrect
    : Prop :=
  ∀ left right (modulus : Nat) output, 0 < modulus →
    Generated.StrictRecombine.multiplyNormalizeModRaw left right
      (modulus : ZZ) = .ok output →
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus left *
        Refinement.StrictHensel.toPolyMod modulus right

noncomputable def SelectedProductMod (modulus : Nat) (candidate : Array Int32)
    (activeLifted : Array SparsePolyZZ) (index : Nat) :
    Polynomial (ZMod modulus) :=
  ((candidate.toList.drop index).map fun activeIndex =>
    Refinement.StrictHensel.toPolyMod modulus
      activeLifted[activeIndex.toInt64.toNat]!).prod

private noncomputable def intTermsToPoly (terms : List (UMonomial × Int)) :
    Polynomial Int :=
  (terms.map fun term => Polynomial.monomial term.1.deg term.2).sum

theorem multiplyRowLoop_toPoly (left : UMonomial × Int)
    (right : SparsePolyZZ) (rightIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyRowLoop
      left right rightIndex terms).toList =
      intTermsToPoly terms.toList +
        Polynomial.monomial left.1.deg left.2 *
          intTermsToPoly (right.toList.drop rightIndex) := by
  induction hmeasure : right.size - rightIndex using Nat.strong_induction_on
      generalizing rightIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyRowLoop]
      split
      next hright =>
        rw [ih (right.size - (rightIndex + 1)) (by omega)
          (rightIndex + 1) (terms.push
            (⟨left.1.deg + right[rightIndex].1.deg⟩,
              left.2 * right[rightIndex].2)) rfl]
        have hsuffix : right.toList.drop rightIndex = right[rightIndex] ::
            right.toList.drop (rightIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := right.toList) (i := rightIndex) (by simpa using hright)
        simp [intTermsToPoly, hsuffix, add_comm, add_left_comm,
          mul_add]
        exact (Polynomial.monomial_mul_monomial _ _ _ _).symm
      next hright =>
        have hle : right.size ≤ rightIndex := Nat.le_of_not_gt hright
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem multiplyTermsLoop_toPoly (left right : SparsePolyZZ)
    (leftIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyTermsLoop
      left right leftIndex terms).toList =
      intTermsToPoly terms.toList +
        intTermsToPoly (left.toList.drop leftIndex) *
          intTermsToPoly right.toList := by
  induction hmeasure : left.size - leftIndex using Nat.strong_induction_on
      generalizing leftIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyTermsLoop]
      split
      next hleft =>
        rw [ih (left.size - (leftIndex + 1)) (by omega)
          (leftIndex + 1)
          (Generated.StrictRecombine.multiplyRowLoop left[leftIndex]
            right 0 terms) rfl]
        rw [multiplyRowLoop_toPoly]
        have hsuffix : left.toList.drop leftIndex = left[leftIndex] ::
            left.toList.drop (leftIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := left.toList) (i := leftIndex) (by simpa using hleft)
        simp [intTermsToPoly, hsuffix]
        ring
      next hleft =>
        have hle : left.size ≤ leftIndex := Nat.le_of_not_gt hleft
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

private theorem intTermsToPoly_modify_add (terms : List (UMonomial × Int))
    (index : Nat) (term : UMonomial × Int) (hindex : index < terms.length)
    (hdegree : terms[index].1.deg = term.1.deg) :
    intTermsToPoly (terms.modify index
      (fun existing => (existing.1, existing.2 + term.2))) =
      intTermsToPoly terms + Polynomial.monomial term.1.deg term.2 := by
  induction terms generalizing index with
  | nil => simp at hindex
  | cons head tail ih =>
      cases index with
      | zero =>
          simp [intTermsToPoly] at hdegree ⊢
          rw [hdegree]
          abel
      | succ index =>
          have htail : index < tail.length := by simpa using hindex
          have hdegree' : tail[index].1.deg = term.1.deg := by
            simpa using hdegree
          simp only [List.modify_succ_cons, intTermsToPoly, List.map_cons,
            List.sum_cons]
          change Polynomial.monomial head.1.deg head.2 +
              intTermsToPoly (tail.modify index
                (fun existing => (existing.1, existing.2 + term.2))) = _
          rw [ih index htail hdegree']
          change _ + (intTermsToPoly tail + _) =
            _ + intTermsToPoly tail + _
          abel

private theorem groupTermsStep_toPoly (acc : SparsePolyZZ)
    (term : UMonomial × Int) :
    intTermsToPoly
      ((match acc.findIdx? (fun t : UMonomial × Int =>
          t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term).toList) =
      intTermsToPoly acc.toList + Polynomial.monomial term.1.deg term.2 := by
  split
  next index hfind =>
    obtain ⟨hindex, hdegree, _⟩ :=
      Array.findIdx?_eq_some_iff_getElem.mp hfind
    simpa using intTermsToPoly_modify_add acc.toList index term
      (by simpa using hindex) (by simpa using hdegree)
  next hfind =>
    simp [intTermsToPoly]

private theorem groupTerms_toPoly_aux (source : List (UMonomial × Int))
    (acc : SparsePolyZZ) :
    intTermsToPoly
      (source.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) acc).toList =
      intTermsToPoly acc.toList + intTermsToPoly source := by
  induction source generalizing acc with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      simp only [List.foldl_cons]
      rw [ih]
      rw [groupTermsStep_toPoly]
      simp [intTermsToPoly]
      abel

private theorem groupTerms_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly
      (terms.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) #[]).toList =
      intTermsToPoly terms.toList := by
  rw [← Array.foldl_toList]
  simpa [intTermsToPoly] using groupTerms_toPoly_aux terms.toList #[]

private theorem filterNonzero_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly (terms.filter (fun term => term.2 ≠ 0)).toList =
      intTermsToPoly terms.toList := by
  rw [Array.toList_filter]
  induction terms.toList with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      rcases head with ⟨monomial, coefficient⟩
      by_cases hzero : coefficient = 0
      · simpa [hzero, intTermsToPoly] using ih
      · simpa [hzero, intTermsToPoly] using ih

private theorem mergeSort_toPoly (terms : List (UMonomial × Int)) :
    intTermsToPoly
        (terms.mergeSort (fun a b => a.fst.deg > b.fst.deg)) =
      intTermsToPoly terms := by
  exact ((List.mergeSort_perm terms _).map
    (fun term => Polynomial.monomial term.1.deg term.2)).sum_eq

theorem normalization_toPoly (terms : SparsePolyZZ) :
    SparsePolyZZ.toPoly (SparsePolyZZ.normalization terms) =
      SparsePolyZZ.toPoly terms := by
  unfold SparsePolyZZ.normalization SparsePolyZZ.toPoly
  rw [List.toList_toArray]
  change intTermsToPoly
      (((terms.foldl (fun acc term =>
          match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
          | some index => acc.modify index
              (fun existing => (existing.1, existing.2 + term.2))
          | none => acc.push term) #[]).filter
        (fun term => term.2 ≠ 0)).toList.mergeSort
          (fun a b => a.fst.deg > b.fst.deg)) = intTermsToPoly terms.toList
  rw [mergeSort_toPoly, filterNonzero_toPoly, groupTerms_toPoly]

theorem multiplyNormalizeRaw_toPoly (left right output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.multiplyNormalizeRaw left right =
      .ok output) :
    SparsePolyZZ.toPoly output =
      SparsePolyZZ.toPoly left * SparsePolyZZ.toPoly right := by
  unfold Generated.StrictRecombine.multiplyNormalizeRaw at hrun
  have houtput := Except.ok.inj hrun
  subst output
  rw [normalization_toPoly]
  change intTermsToPoly
      (Generated.StrictRecombine.multiplyTermsLoop left right 0 #[]).toList = _
  rw [multiplyTermsLoop_toPoly]
  simp [intTermsToPoly, SparsePolyZZ.toPoly]

private theorem emod_cast (modulus : Nat) (hmodulus : 0 < modulus)
    (coefficient : Int) :
    ((coefficient % (modulus : Int) : Int) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
  rw [ZMod.intCast_eq_intCast_iff]
  refine Int.modEq_iff_dvd.mpr ?_
  use coefficient / (modulus : Int)
  have h := Int.mul_ediv_add_emod coefficient (modulus : Int)
  omega

theorem modCoeffLoop_toPolyMod (input : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.modCoeffLoop input (modulus : ZZ)
      index result = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus result +
        Refinement.StrictHensel.termsToPolyMod modulus
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.modCoeffLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hmodulus' =>
          let coefficient := input[index].2 % (modulus : Int)
          by_cases hzero : coefficient = 0
          · change input[index].2 % (modulus : Int) = 0 at hzero
            simp only [hzero, if_true] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1) result output hrun rfl]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            have hcast := emod_cast modulus hmodulus input[index].2
            change (coefficient : ZMod modulus) =
              (input[index].2 : ZMod modulus) at hcast
            dsimp [coefficient] at hcast
            rw [hzero] at hcast
            rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
            rw [← hcast]
            simp
          · change input[index].2 % (modulus : Int) ≠ 0 at hzero
            simp only [hzero, if_false] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1) (result.push (input[index].1, coefficient))
              output hrun rfl]
            rw [Refinement.StrictHensel.toPolyMod_push]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [emod_cast modulus hmodulus input[index].2]
            rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
            abel
        next hmodulus' => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

theorem multiplyNormalizeModRaw_correct : TrialProductStepCorrect := by
  intro left right modulus output hmodulus hrun
  unfold Generated.StrictRecombine.multiplyNormalizeModRaw at hrun
  split at hrun
  next fault hmultiply => contradiction
  next product hmultiply =>
    have hmod := modCoeffLoop_toPolyMod product modulus hmodulus 0 #[]
      output hrun
    have hmul := multiplyNormalizeRaw_toPoly left right product hmultiply
    have hmod' : Refinement.StrictHensel.toPolyMod modulus output =
        Refinement.StrictHensel.toPolyMod modulus product := by
      simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using hmod
    rw [hmod']
    simpa [Refinement.StrictHensel.toPolyMod] using
      congrArg (Polynomial.map (Int.castRingHom (ZMod modulus))) hmul

theorem trialProductLoop_refines
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (product output : SparsePolyZZ)
    (hvalid : CandidateIndicesValid candidate
      (Array.replicate activeLifted.size false))
    (hrun : Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted (modulus : ZZ) index product = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus product *
        SelectedProductMod modulus candidate activeLifted index := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product output with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop] at hrun
      split at hrun
      next hindex =>
        have hentry := hvalid index hindex
        have hactive : candidate[index].toInt64.toNat < activeLifted.size := by
          simpa using hentry.2
        rw [dif_pos hentry.1, dif_pos hactive] at hrun
        split at hrun
        next fault hcall => contradiction
        next product' hcall =>
          let selected := activeLifted[candidate[index].toInt64.toNat]
          have htail := ih (candidate.size - (index + 1)) (by omega)
            (index + 1) product' output hrun rfl
          have hsuffix : candidate.toList.drop index = candidate[index] ::
              candidate.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := candidate.toList) (i := index) (by simpa using hindex)
          have hstep' := multiplyNormalizeModRaw_correct product selected
            modulus product' hmodulus hcall
          rw [htail, hstep']
          have hselected : selected =
              activeLifted[candidate[index].toInt64.toNat]! := by
            simp [selected, hactive]
          rw [hselected]
          simp [SelectedProductMod, hsuffix, mul_assoc]
      next hindex =>
        have hle : candidate.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [SelectedProductMod, List.drop_eq_nil_iff.mpr hle]

private theorem symmetricMod_cast (modulus : Nat) (hmodulus : 0 < modulus)
    (coefficient : Int) :
    (ZZ.symmetricMod coefficient (modulus : ZZ) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
  unfold ZZ.symmetricMod
  have hfmod : ((Int.fmod coefficient (modulus : Int) : Int) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
    rw [Int.fmod_eq_emod_of_nonneg coefficient (by omega)]
    rw [ZMod.intCast_eq_intCast_iff]
    refine Int.modEq_iff_dvd.mpr ?_
    use coefficient / (modulus : Int)
    have h := Int.mul_ediv_add_emod coefficient (modulus : Int)
    omega
  dsimp [ZZ.symmetricMod]
  split
  next hsmall => exact hfmod
  next hlarge =>
    rw [Int.cast_sub, hfmod]
    simp

theorem symmetricModLoop_toPolyMod (input : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.symmetricModLoop input (modulus : ZZ)
      index result = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus result +
        Refinement.StrictHensel.termsToPolyMod modulus
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.symmetricModLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        let coefficient := ZZ.symmetricMod input[index].2 (modulus : ZZ)
        by_cases hzero : coefficient = 0
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) = 0 at hzero
          simp only [hzero, if_true] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) result output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          have hcast := symmetricMod_cast modulus hmodulus input[index].2
          change (ZZ.symmetricMod input[index].2 (modulus : ZZ) : ZMod modulus) =
            (input[index].2 : ZMod modulus) at hcast
          rw [hzero] at hcast
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          rw [← hcast]
          simp
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) ≠ 0 at hzero
          simp only [hzero, if_false] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) (result.push (input[index].1,
              ZZ.symmetricMod input[index].2 (modulus : ZZ))) output
            hrun rfl]
          rw [Refinement.StrictHensel.toPolyMod_push]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [symmetricMod_cast modulus hmodulus input[index].2]
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          abel
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

theorem symmetricModRaw_toPolyMod (input output : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus input := by
  unfold Generated.StrictRecombine.symmetricModRaw at hrun
  rw [dif_pos (by exact_mod_cast hmodulus)] at hrun
  simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using
    symmetricModLoop_toPolyMod input modulus hmodulus 0 #[] output hrun

theorem primitiveDivideLoop_toPoly (input : SparsePolyZZ) (divisor : Int)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.primitiveDivideLoop input divisor index
      result = .ok output) :
    Polynomial.C divisor * SparsePolyZZ.toPoly output =
      Polynomial.C divisor * SparsePolyZZ.toPoly result +
        intTermsToPoly (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.primitiveDivideLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hdivisor =>
          split at hrun
          next hdivides =>
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1)
              (result.push (input[index].1, input[index].2 / divisor))
              output hrun rfl]
            rw [SparsePolyZZ.toPoly]
            simp only [Array.toList_push, List.map_append, List.map_singleton,
              List.sum_append, List.sum_singleton]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [hsuffix]
            simp only [intTermsToPoly, List.map_cons, List.sum_cons]
            rw [mul_add, Polynomial.C_mul_monomial]
            have hcancel : divisor * (input[index].2 / divisor) =
                input[index].2 := by
              exact Int.mul_ediv_cancel_of_dvd hdivides
            rw [hcancel]
            change Polynomial.C divisor * SparsePolyZZ.toPoly result + _ + _ =
              Polynomial.C divisor * SparsePolyZZ.toPoly result + (_ + _)
            abel
          next hdivides => contradiction
        next hdivisor => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem primitiveRaw_toPoly (input primitive : SparsePolyZZ) (content : Int)
    (hrun : Generated.StrictRecombine.primitiveRaw input =
      .ok (content, primitive)) :
    SparsePolyZZ.toPoly input =
      Polynomial.C content * SparsePolyZZ.toPoly primitive := by
  unfold Generated.StrictRecombine.primitiveRaw at hrun
  split at hrun
  next hempty =>
    have hinput : input = #[] := Array.isEmpty_iff.mp hempty
    subst input
    have hout := Except.ok.inj hrun
    cases hout
    simp [SparsePolyZZ.toPoly]
  next hempty =>
    dsimp at hrun
    split at hrun
    next fault hdivide => contradiction
    next primitive' hdivide =>
      have hout : (content, primitive) =
          (if (input[0]'(by
              have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
              omega)).2 < 0 then
            -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
          else (Generated.StrictRecombine.contentLoop input 0 0 : Int), primitive') :=
        (Except.ok.inj hrun).symm
      have hcontent := congrArg Prod.fst hout
      have hprimitive := congrArg Prod.snd hout
      have hsemantic := primitiveDivideLoop_toPoly input
        (if (input[0]'(by
            have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
            omega)).2 < 0 then
          -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
        else (Generated.StrictRecombine.contentLoop input 0 0 : Int)) 0 #[]
        primitive' hdivide
      have hcontent' : content =
          (if (input[0]'(by
              have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
              omega)).2 < 0 then
            -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
          else (Generated.StrictRecombine.contentLoop input 0 0 : Int)) := hcontent
      rw [← hcontent'] at hsemantic
      have hprimitive' : primitive = primitive' := hprimitive
      rw [← hprimitive'] at hsemantic
      simpa [SparsePolyZZ.toPoly, intTermsToPoly] using hsemantic.symm

theorem subtractScaledTermsLoop_toPoly (divisor : SparsePolyZZ)
    (scale : Int) (degreeShift index : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly
        (Generated.StrictRecombine.subtractScaledTermsLoop divisor scale
          degreeShift index terms).toList =
      intTermsToPoly terms.toList -
        Polynomial.monomial degreeShift scale *
          intTermsToPoly (divisor.toList.drop index) := by
  induction hmeasure : divisor.size - index using Nat.strong_induction_on
      generalizing index terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.subtractScaledTermsLoop]
      split
      next hindex =>
        rw [ih (divisor.size - (index + 1)) (by omega)
          (index + 1)
          (terms.push (⟨divisor[index].1.deg + degreeShift⟩,
            -(scale * divisor[index].2))) rfl]
        have hsuffix : divisor.toList.drop index = divisor[index] ::
            divisor.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := divisor.toList) (i := index) (by simpa using hindex)
        simp [intTermsToPoly, hsuffix]
        rw [mul_add, Polynomial.monomial_mul_monomial]
        rw [add_comm divisor[index].1.deg degreeShift]
        abel
      next hindex =>
        have hle : divisor.size ≤ index := Nat.le_of_not_gt hindex
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem subtractScaledNormalize_toPoly (remainder divisor : SparsePolyZZ)
    (scale : Int) (degreeShift : Nat) :
    SparsePolyZZ.toPoly
        (Generated.StrictRecombine.subtractScaledNormalize remainder divisor
          scale degreeShift) =
      SparsePolyZZ.toPoly remainder -
        Polynomial.monomial degreeShift scale * SparsePolyZZ.toPoly divisor := by
  unfold Generated.StrictRecombine.subtractScaledNormalize
  rw [normalization_toPoly]
  change intTermsToPoly
      (Generated.StrictRecombine.subtractScaledTermsLoop divisor scale
        degreeShift 0 remainder).toList = _
  rw [subtractScaledTermsLoop_toPoly]
  simp [intTermsToPoly, SparsePolyZZ.toPoly]

theorem exactDivmodLoop_toPoly (divisor remainder quotient q r : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
      .ok (q, r)) :
    SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly q +
        SparsePolyZZ.toPoly r =
      SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly quotient +
        SparsePolyZZ.toPoly remainder := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient q r with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop] at hrun
      split at hrun
      next hremainder =>
        split at hrun
        next hdivisor =>
          dsimp at hrun
          split at hrun
          next hdegree =>
            split at hrun
            next hnonzero =>
              split at hrun
              next hdivides =>
                split at hrun
                next hdecrease =>
                  let degreeShift := remainder[0].1.deg - divisor[0].1.deg
                  let scale := remainder[0].2 / divisor[0].2
                  let remainder' :=
                    Generated.StrictRecombine.subtractScaledNormalize remainder
                      divisor scale degreeShift
                  let quotient' := quotient.push (⟨degreeShift⟩, scale)
                  have hdec : Generated.StrictRecombine.divisionRank remainder' <
                      measure := by
                    rw [← hmeasure]
                    exact hdecrease
                  have htail := ih
                    (Generated.StrictRecombine.divisionRank remainder')
                    hdec remainder' quotient' q r hrun rfl
                  rw [htail, subtractScaledNormalize_toPoly]
                  simp only [quotient', SparsePolyZZ.toPoly, Array.toList_push,
                    List.map_append, List.map_singleton, List.sum_append,
                    List.sum_singleton]
                  change _ * (_ + Polynomial.monomial degreeShift scale) +
                      (_ - Polynomial.monomial degreeShift scale * _) = _
                  ring
                next hdecrease => contradiction
              next hdivides =>
                have hout := Except.ok.inj hrun
                cases hout
                rfl
            next hnonzero => contradiction
          next hdegree =>
            have hout := Except.ok.inj hrun
            cases hout
            rfl
        next hdivisor => contradiction
      next hremainder =>
        have hout := Except.ok.inj hrun
        cases hout
        rfl

theorem exactDivmodRaw_toPoly (dividend divisor quotient remainder : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.exactDivmodRaw dividend divisor =
      .ok (quotient, remainder)) :
    SparsePolyZZ.toPoly dividend =
      SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly quotient +
        SparsePolyZZ.toPoly remainder := by
  unfold Generated.StrictRecombine.exactDivmodRaw at hrun
  have hsemantic := exactDivmodLoop_toPoly divisor dividend #[] quotient
    remainder hrun
  simpa [SparsePolyZZ.toPoly] using hsemantic.symm

theorem successfulTrialExtraction_toPoly
    (fStar factor quotient quotientPrimitive : SparsePolyZZ)
    (quotientContent : Int)
    (hdivide : Generated.StrictRecombine.exactDivmodRaw fStar factor =
      .ok (quotient, #[]))
    (hprimitive : Generated.StrictRecombine.primitiveRaw quotient =
      .ok (quotientContent, quotientPrimitive)) :
    SparsePolyZZ.toPoly fStar =
      Polynomial.C quotientContent *
        (SparsePolyZZ.toPoly factor * SparsePolyZZ.toPoly quotientPrimitive) := by
  rw [exactDivmodRaw_toPoly fStar factor quotient #[] hdivide]
  rw [primitiveRaw_toPoly quotient quotientPrimitive quotientContent hprimitive]
  simp [SparsePolyZZ.toPoly]
  ring

theorem zassenhausAttempt_extracted_toPoly
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok
        (.extracted factor quotientPrimitive)) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  unfold Generated.StrictRecombine.zassenhausAttempt at hrun
  split at hrun
  next hfstar =>
    dsimp at hrun
    cases hleading : Generated.StrictRecombine.selectedLeadingProductLoop
        candidate activeLifted 0 fStar[0].2 with
    | error fault => simp [hleading] at hrun
    | ok leadingProduct =>
      simp only [hleading] at hrun
      split at hrun
      next hpruned => simp at hrun
      next hleadingAccepted =>
        cases hconstant :
            Generated.StrictRecombine.selectedConstantProductLoop candidate
              activeLifted 0 fStar[0].2 with
        | error fault => simp [hconstant] at hrun
        | ok constantProduct =>
          simp only [hconstant] at hrun
          split at hrun
          next hpruned => simp at hrun
          next hconstantAccepted =>
            cases hconvert : Generated.StrictRecombine.combinationToInt32
                candidate with
            | error fault => simp [hconvert] at hrun
            | ok candidate32 =>
              simp only [hconvert] at hrun
              cases hproduct : Generated.StrictRecombine.trialProductLoop
                  ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, fStar[0].2)] with
              | error fault => simp [hproduct] at hrun
              | ok product =>
                simp only [hproduct] at hrun
                cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                    product modulus with
                | error fault => simp [hsymmetric] at hrun
                | ok symmetric =>
                  simp only [hsymmetric] at hrun
                  cases hprimitive : Generated.StrictRecombine.primitiveRaw
                      symmetric with
                  | error fault => simp [hprimitive] at hrun
                  | ok primitiveResult =>
                    rcases primitiveResult with ⟨symmetricContent,
                      recoveredFactor⟩
                    simp only [hprimitive] at hrun
                    cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                        fStar recoveredFactor with
                    | error fault => simp [hdivmod] at hrun
                    | ok divResult =>
                      rcases divResult with ⟨quotient, remainder⟩
                      simp only [hdivmod] at hrun
                      by_cases hremainder : remainder.isEmpty = true
                      · simp only [hremainder, if_true] at hrun
                        have hremainderEmpty : remainder = #[] :=
                          Array.isEmpty_iff.mp hremainder
                        subst remainder
                        cases hquotientPrimitive :
                            Generated.StrictRecombine.primitiveRaw quotient with
                        | error fault => simp [hquotientPrimitive] at hrun
                        | ok quotientResult =>
                          rcases quotientResult with ⟨quotientContent,
                            recoveredQuotient⟩
                          simp only [hquotientPrimitive] at hrun
                          have hout := Except.ok.inj hrun
                          injection hout with hfactor hquotient
                          subst factor
                          subst quotientPrimitive
                          exact ⟨quotientContent,
                            successfulTrialExtraction_toPoly fStar
                              recoveredFactor quotient recoveredQuotient
                              quotientContent hdivmod hquotientPrimitive⟩
                      · simp only [hremainder, if_false] at hrun
                        simp at hrun
  next hfstar => contradiction

theorem zassenhausAttempt_extracted_unit_scalar
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok
        (.extracted factor quotientPrimitive)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  rcases zassenhausAttempt_extracted_toPoly fStar factor quotientPrimitive
      activeLifted modulus candidate hrun with ⟨scalar, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul] at hcontent
  have hnormalizeUnit : IsUnit (normalize scalar) :=
    IsUnit.of_mul_eq_one _ hcontent.symm
  have hscalarUnit : IsUnit scalar :=
    (associated_normalize scalar).isUnit_iff.mpr hnormalizeUnit
  exact ⟨scalar, hscalarUnit, hproduct⟩

theorem scanZassenhausCombinations_extracted_unit_scalar
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst factor
          subst quotientPrimitive
          exact zassenhausAttempt_extracted_unit_scalar fStar extractedFactor
            extractedQuotient activeLifted modulus start hprimitive hattempt
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next)
              hdecrease next hvalidNext hrun rfl

private noncomputable def factorArrayProduct (factors : Array SparsePolyZZ) :
    Polynomial Int :=
  (factors.toList.map SparsePolyZZ.toPoly).prod

theorem validateCandidatesLoop_product
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            have hremainderEmpty : remainder = #[] :=
                              Array.isEmpty_iff.mp hremainder
                            subst remainder
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                rcases ih (candidates.size - (candidateIndex + 1))
                                    (by omega) (candidateIndex := candidateIndex + 1)
                                    (fStar := quotientPrimitive)
                                    (result := result.push factor)
                                    (consumed := consumedNext)
                                    (remaining := remaining - candidates[candidateIndex].size)
                                    (fStar' := fStar') (result' := result')
                                    (consumed' := consumed') hrun rfl with
                                  ⟨tailScalar, htail⟩
                                refine ⟨quotientContent * tailScalar, ?_⟩
                                rw [successfulTrialExtraction_toPoly fStar factor
                                  quotient quotientPrimitive quotientContent hdivmod
                                  hquotientPrimitive]
                                simp [factorArrayProduct] at htail ⊢
                                calc
                                  _ = (quotientContent : Polynomial Int) *
                                      (quotientPrimitive.toPoly *
                                        ((result.toList.map SparsePolyZZ.toPoly).prod *
                                          factor.toPoly)) := by ring
                                  _ = (quotientContent : Polynomial Int) *
                                      ((tailScalar : Polynomial Int) *
                                        (fStar'.toPoly *
                                          (result'.toList.map SparsePolyZZ.toPoly).prod)) := by
                                    rw [htail]
                                  _ = _ := by ring
                          · simp only [hremainder] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        have hout := Except.ok.inj hrun
        cases hout
        exact ⟨1, by simp⟩

theorem validateCandidates_product
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact validateCandidatesLoop_product ops candidates 0 activeLifted modulus
    fStar fStar' result result' (Array.replicate activeLifted.size false)
    consumed activeLifted.size hrun

/-- A successful concrete validation run cannot hide a non-unit integer
content factor when its incoming accumulated product is primitive.  This is
the content bridge needed to turn `validateCandidates_product` into an
associated-product statement, without assuming that the generated primitive
normalization returned a semantically convenient result. -/
theorem validateCandidates_product_unit_scalar
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  rcases validateCandidates_product ops fStar activeLifted modulus candidates
      result fStar' result' consumed hrun with ⟨scalar, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul] at hcontent
  have hnormalizeUnit : IsUnit (normalize scalar) :=
    IsUnit.of_mul_eq_one _ hcontent.symm
  have hscalarUnit : IsUnit scalar :=
    (associated_normalize scalar).isUnit_iff.mpr hnormalizeUnit
  exact ⟨scalar, hscalarUnit, hproduct⟩

/-- The concrete extraction loop preserves primitivity of the whole live
product.  This is the induction invariant used by the source-shaped van-Hoeij
loop after it removes the consumed lifted factors. -/
theorem validateCandidates_preserves_primitive
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    (SparsePolyZZ.toPoly fStar' * factorArrayProduct result').IsPrimitive := by
  rcases validateCandidates_product_unit_scalar ops fStar activeLifted modulus
      candidates result fStar' result' consumed hprimitive hrun with
    ⟨scalar, hscalarUnit, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul,
    normalize_eq_one.mpr hscalarUnit, one_mul] at hcontent
  exact Polynomial.isPrimitive_iff_content_eq_one.mpr hcontent.symm

end Refinement.StrictRecombine
