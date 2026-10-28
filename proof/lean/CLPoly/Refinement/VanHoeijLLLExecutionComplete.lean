import CLPoly.Refinement.VanHoeijLLLStepTotal
import CLPoly.Refinement.VanHoeijInitialExecution
import CLPoly.Generated.StrictVanHoeij

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine

theorem insertShortRow_succeeds (matrix : Generated.StrictRecombine.LLLMatrix) (row : Nat)
    (sorted : Array Nat) (position : Nat) (hr : row < matrix.size)
    (hs : ShortRowIndicesValid matrix sorted) :
    ∃ output, insertShortRow matrix row sorted position = .ok output := by
  induction hmeasure : sorted.size - position using Nat.strong_induction_on generalizing position with
  | h measure ih =>
    rw [insertShortRow]
    by_cases hp : position < sorted.size
    · have he : sorted[position] < matrix.size := hs _ (Array.getElem_mem_toList hp)
      rw [dif_pos hp, dif_pos hr, dif_pos he,
        dotRows_eq_fin_sum matrix[row] matrix[row] le_rfl]
      dsimp only
      rw [dotRows_eq_fin_sum matrix[sorted[position]] matrix[sorted[position]] le_rfl]
      dsimp only
      split
      · exact ⟨_, rfl⟩
      · exact ih _ (by omega) (position + 1) rfl
    · exact ⟨sorted.push row, by rw [dif_neg hp]⟩

theorem collectShortRows_succeeds (matrix : Generated.StrictRecombine.LLLMatrix) (bound : ZZ)
    (row : Nat) (result : Array Nat) (hresult : ShortRowIndicesValid matrix result) :
    ∃ output, collectShortRows matrix bound row result = .ok output := by
  induction hmeasure : matrix.size - row using Nat.strong_induction_on generalizing row result with
  | h measure ih =>
    rw [collectShortRows]
    by_cases hr : row < matrix.size
    · rw [dif_pos hr, dotRows_eq_fin_sum matrix[row] matrix[row] le_rfl]
      dsimp only
      split
      · obtain ⟨next, hn⟩ := insertShortRow_succeeds matrix row result 0 hr hresult
        rw [hn]
        exact ih _ (by omega) (row + 1) next
          (insertShortRow_indices_valid matrix row result next 0 hr hresult hn) rfl
      · exact ih _ (by omega) (row + 1) result hresult rfl
    · exact ⟨result, by rw [dif_neg hr]⟩

/-- Full execution of the actual repaired LLL interface, including its
legacy short-row collection and the same final cached GS norms. -/
theorem lllReduceWithNorms_succeeds (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (hnonempty : 0 < matrix.size) (bound : ZZ) :
    ∃ output, Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output := by
  obtain ⟨⟨mu, norms, transform⟩, hi⟩ := initializeLLL_succeeds matrix hnonempty hinput.rows_square
  let state := Generated.StrictRecombine.LLLState.mk matrix transform mu norms 1
  have hv : ConcreteLLLExecutionValid state :=
    initializeLLL_concrete_valid matrix mu norms transform hinput hi
  have ht := initializeLLL_transform_rel matrix mu norms transform hi
  obtain ⟨reduced, hr⟩ := lllMainLoop_succeeds matrix state hv ht (by change 0 < (1 : Nat); decide)
  obtain ⟨rows, hc⟩ := collectShortRows_succeeds reduced.matrix bound 0 #[]
    (by intro i hi; simp at hi)
  refine ⟨(reduced.matrix, reduced.transform, rows, reduced.norms), ?_⟩
  unfold Generated.StrictVanHoeij.lllReduceWithNorms
  split
  next fault hbad => rw [hi] at hbad; contradiction
  next mu' norms' transform' hinitialize =>
    have he : (mu', norms', transform') = (mu, norms, transform) :=
      Except.ok.inj (hinitialize.symm.trans hi)
    rcases Prod.mk.inj he with ⟨hmu, he⟩
    rcases Prod.mk.inj he with ⟨hnorms, htransform⟩
    subst mu'; subst norms'; subst transform'
    change (match lllMainLoop concreteLLLTermination state hv with
      | .error fault => Except.error fault
      | .ok next => match collectShortRows next.matrix bound 0 #[] with
        | .error fault => Except.error fault
        | .ok found => Except.ok (next.matrix, next.transform, found, next.norms)) = _
    rw [hr]
    dsimp only
    rw [hc]

end Refinement.StrictRecombine
