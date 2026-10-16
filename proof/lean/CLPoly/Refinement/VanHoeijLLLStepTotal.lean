import CLPoly.Refinement.VanHoeijStepExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine

theorem lllStep_succeeds (initial : Generated.StrictRecombine.LLLMatrix)
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrel : LLLTransformRel initial state.matrix state.transform)
    (hpositive : 0 < state.k) (hk : state.k < state.matrix.size) :
    ∃ branch, lllStep state = .ok branch ∧ 0 < branch.state.k := by
  obtain ⟨reduced, hr⟩ := sizeReduceAt_succeeds state (state.k - 1) hvalid hk (by omega)
    (hrel.2.1.trans hrel.1.symm) (by intro i hi; rw [hrel.2.2.1 i hi, hrel.1])
  have hv := sizeReduceAt_preserves_execution_valid state reduced (state.k - 1) hvalid (by omega) hr
  have hrel' := sizeReduceAt_preserves_transform_rel initial state reduced (state.k - 1) hrel hr
  have hkr := (sizeReduceAt_preserves_norms_k state reduced (state.k - 1) hr).2
  have hkn : reduced.k < reduced.matrix.size := by rw [hkr, hrel'.1, ← hrel.1]; exact hk
  have hkp : 0 < reduced.k := by rw [hkr]; exact hpositive
  have hknorm : reduced.k < reduced.norms.size := by rw [hv.norms_size]; exact hkn
  have hpnorm : reduced.k - 1 < reduced.norms.size := by omega
  have hkmu : reduced.k < reduced.mu.size := by rw [hv.mu_size]; exact hkn
  have hpmu : reduced.k - 1 < reduced.mu[reduced.k].size := by
    rw [hv.mu_rows_square reduced.k hkmu]; omega
  rw [lllStep, dif_pos hpositive, dif_pos hk, hr]
  dsimp only
  rw [dif_pos hknorm, dif_pos hpnorm, dif_pos hkmu, dif_pos hpmu]
  split
  next hlovasz =>
    obtain ⟨fully, hf⟩ := extraSizeReduceLoop_succeeds initial reduced (reduced.k - 1)
      hv hrel' hkn (Nat.sub_le _ _)
    refine ⟨.advanced { fully with k := fully.k + 1 }, ?_, by simp [LLLStepResult.state]⟩
    rw [hf]
  next hlovasz =>
    obtain ⟨matrix, hmatrix⟩ := swapMatrixRows_succeeds reduced.matrix reduced.k
      (reduced.k - 1) hkn (by omega)
    have hktu : reduced.k < reduced.transform.size := by rw [hrel'.2.1, ← hrel'.1]; exact hkn
    obtain ⟨transform, htransform⟩ := swapMatrixRows_succeeds reduced.transform reduced.k
      (reduced.k - 1) hktu (by omega)
    obtain ⟨swapped, hswap, hswapSize, hswapRows⟩ := swapQQRows_succeeds_square reduced.mu
      reduced.k (reduced.k - 1) reduced.matrix.size hkmu (by omega) hv.mu_rows_square
    rw [hmatrix, htransform, hswap]
    dsimp only
    have hks : reduced.k < swapped.size := by rw [hswapSize]; exact hkmu
    have hps : reduced.k - 1 < swapped[reduced.k].size := by rw [hswapRows reduced.k hks]; omega
    rw [dif_pos hks, dif_pos hps]
    let old := reduced.mu[reduced.k][reduced.k - 1]
    let newNorm := reduced.norms[reduced.k] + old * old * reduced.norms[reduced.k - 1]
    let newMu := if newNorm ≠ 0 then old * reduced.norms[reduced.k - 1] / newNorm else 0
    let corrected := swapped.set reduced.k (swapped[reduced.k].set (reduced.k - 1) newMu)
    have hshape : GramStorageShape corrected reduced.norms reduced.matrix.size := by
      apply GramStorageShape.setMu (hshape := ⟨hswapSize.trans hv.mu_size, hswapRows, hv.norms_size⟩)
    have hu := updateMuAfterSwapLoop_eq_array corrected reduced.k old newMu (reduced.k + 1)
      (by intro i hi; rw [hshape.mu_rows i hi]; exact hkn)
      (by intro i hi; rw [hshape.mu_rows i hi]; omega)
    refine ⟨.swapped (Generated.StrictRecombine.LLLState.mk matrix transform
      (updateMuAfterSwapArray corrected reduced.k old newMu (reduced.k + 1))
      (normsAfterLovaszSwap reduced.norms reduced.k old) (max (reduced.k - 1) 1)), ?_, ?_⟩
    · change (match updateMuAfterSwapLoop corrected reduced.k old newMu (reduced.k + 1) with
        | .error fault => Except.error fault
        | .ok finalMu => Except.ok (LLLStepResult.swapped (Generated.StrictRecombine.LLLState.mk matrix transform finalMu
            (normsAfterLovaszSwap reduced.norms reduced.k old) (max (reduced.k - 1) 1)))) = _
      rw [hu]
    · change 0 < max (reduced.k - 1) 1
      omega

theorem lllMainLoop_succeeds (initial : Generated.StrictRecombine.LLLMatrix)
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrel : LLLTransformRel initial state.matrix state.transform) (hk : 0 < state.k) :
    ∃ output, lllMainLoop concreteLLLTermination state hvalid = .ok output := by
  induction hmeasure : concreteLLLRank state using Nat.strong_induction_on generalizing state with
  | h measure ih =>
    rw [lllMainLoop]
    by_cases hlt : state.k < state.matrix.size
    · rw [dif_pos hlt]
      obtain ⟨branch, hb, hpositive⟩ := lllStep_succeeds initial state hvalid hrel hk hlt
      split
      next fault hstep => rw [hb] at hstep; contradiction
      next nextBranch hstep =>
        have he : nextBranch = branch := Except.ok.inj (hstep.symm.trans hb)
        subst nextBranch
        have hv := concreteLLLTermination.step_valid state branch hvalid hb
        have hnext := lllStep_preserves_transform_rel initial state branch hvalid hrel hb
        have hd := concreteLLLTermination.step_decreases state branch hvalid hb
        change concreteLLLRank branch.state < concreteLLLRank state at hd
        exact ih _ (by omega) branch.state hv hnext hpositive rfl
    · exact ⟨state, by rw [dif_neg hlt]⟩

end Refinement.StrictRecombine
