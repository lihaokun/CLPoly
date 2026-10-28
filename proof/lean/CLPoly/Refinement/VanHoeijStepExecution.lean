import CLPoly.Refinement.VanHoeijRowExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine

theorem extraSizeReduceLoop_succeeds (initial : Generated.StrictRecombine.LLLMatrix)
    (state : Generated.StrictRecombine.LLLState) (remaining : Nat)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrel : LLLTransformRel initial state.matrix state.transform)
    (hk : state.k < state.matrix.size) (hr : remaining ≤ state.k) :
    ∃ output, extraSizeReduceLoop remaining state = .ok output := by
  induction remaining generalizing state with
  | zero => exact ⟨state, by rw [extraSizeReduceLoop]⟩
  | succ remaining ih =>
    obtain ⟨next, hnext⟩ := sizeReduceAt_succeeds state remaining hvalid hk (by omega)
      (hrel.2.1.trans hrel.1.symm) (by intro i hi; rw [hrel.2.2.1 i hi, hrel.1])
    have hn := sizeReduceAt_preserves_execution_valid state next remaining hvalid (by omega) hnext
    have hnr := sizeReduceAt_preserves_transform_rel initial state next remaining hrel hnext
    have hnk := (sizeReduceAt_preserves_norms_k state next remaining hnext).2
    obtain ⟨output, ho⟩ := ih next hn hnr
      (by rw [hnk, hnr.1, ← hrel.1]; exact hk) (by rw [hnk]; omega)
    exact ⟨output, by rw [extraSizeReduceLoop, hnext]; exact ho⟩

theorem swapMatrixRows_succeeds (matrix : Generated.StrictRecombine.LLLMatrix)
    (left right : Nat) (hl : left < matrix.size) (hr : right < matrix.size) :
    ∃ output, swapMatrixRows matrix left right = .ok output := by
  exact ⟨_, by rw [swapMatrixRows, dif_pos hl, dif_pos hr]⟩

theorem swapQQRows_succeeds_square (matrix : Generated.StrictRecombine.QQMatrix)
    (left right size : Nat) (hl : left < matrix.size) (hr : right < matrix.size)
    (hrows : ∀ i (hi : i < matrix.size), matrix[i].size = size) :
    ∃ output, swapQQRows matrix left right = .ok output ∧ output.size = matrix.size ∧
      ∀ i (hi : i < output.size), output[i].size = size := by
  let output := (matrix.set left matrix[right]).set right matrix[left] (by simpa)
  refine ⟨output, by rw [swapQQRows, dif_pos hl, dif_pos hr], by simp [output], ?_⟩
  intro i hi
  have hi' : i < matrix.size := by simpa [output] using hi
  dsimp only [output]
  rw [Array.getElem_set]
  split
  · exact hrows left hl
  · rw [Array.getElem_set]
    split
    · exact hrows right hr
    · exact hrows i hi'

end Refinement.StrictRecombine
