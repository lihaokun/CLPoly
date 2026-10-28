import CLPoly.Refinement.Recombine

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine

theorem subtractRowsLoop_succeeds (targetM sourceM targetU sourceU : Array ZZ)
    (coefficient : ZZ) (size index : Nat) (resultM resultU : Array ZZ)
    (htM : size ≤ targetM.size) (hsM : size ≤ sourceM.size)
    (htU : size ≤ targetU.size) (hsU : size ≤ sourceU.size) :
    ∃ output, subtractRowsLoop targetM sourceM targetU sourceU coefficient size index
      resultM resultU = .ok output := by
  induction hmeasure : size - index using Nat.strong_induction_on generalizing index resultM resultU with
  | h measure ih =>
    rw [subtractRowsLoop]
    by_cases hi : index < size
    · rw [dif_pos hi, dif_pos (lt_of_lt_of_le hi htM), dif_pos (lt_of_lt_of_le hi hsM),
        dif_pos (lt_of_lt_of_le hi htU), dif_pos (lt_of_lt_of_le hi hsU)]
      exact ih _ (by omega) (index + 1) _ _ rfl
    · exact ⟨(resultM, resultU), by rw [dif_neg hi]⟩

theorem subtractMatrixRows_succeeds (matrix transform : Generated.StrictRecombine.LLLMatrix)
    (target source : Nat) (coefficient : ZZ)
    (ht : target < matrix.size) (hs : source < matrix.size)
    (hsize : transform.size = matrix.size)
    (hrows : ∀ i (hi : i < matrix.size), matrix[i].size = matrix.size)
    (htransform : ∀ i (hi : i < transform.size), transform[i].size = matrix.size) :
    ∃ output, subtractMatrixRows matrix transform target source coefficient = .ok output := by
  have htU : target < transform.size := by omega
  have hsU : source < transform.size := by omega
  obtain ⟨⟨rowM, rowU⟩, hrun⟩ := subtractRowsLoop_succeeds matrix[target] matrix[source]
    transform[target] transform[source] coefficient matrix.size 0 #[] #[]
    (hrows target ht).ge (hrows source hs).ge (htransform target htU).ge (htransform source hsU).ge
  refine ⟨(matrix.set target rowM, transform.set target rowU), ?_⟩
  rw [subtractMatrixRows, dif_pos ht, dif_pos hs, dif_pos htU, dif_pos hsU, hrun]

theorem sizeReduceAt_succeeds (state : Generated.StrictRecombine.LLLState) (source : Nat)
    (hvalid : ConcreteLLLExecutionValid state) (hk : state.k < state.matrix.size)
    (hs : source < state.k) (hsize : state.transform.size = state.matrix.size)
    (htransform : ∀ i (hi : i < state.transform.size), state.transform[i].size = state.matrix.size) :
    ∃ output, sizeReduceAt state source = .ok output := by
  have hkMu : state.k < state.mu.size := by rw [hvalid.mu_size]; exact hk
  have hsMu : source < state.mu[state.k].size := by rw [hvalid.mu_rows_square state.k hkMu]; omega
  rw [sizeReduceAt, dif_pos hkMu, dif_pos hsMu, roundQQ_eq_round]
  dsimp only
  by_cases hzero : round state.mu[state.k][source] = 0
  · exact ⟨state, by rw [dif_pos hzero]⟩
  · rw [dif_neg hzero]
    obtain ⟨⟨matrix, transform⟩, hsubtract⟩ := subtractMatrixRows_succeeds state.matrix state.transform
      state.k source (round state.mu[state.k][source]) hk (lt_trans hs hk) hsize
      hvalid.rows_square htransform
    rw [hsubtract]
    have hsMuSize : source < state.mu.size := by rw [hvalid.mu_size]; omega
    have hprefix := reduceMuPrefixLoop_eq_array state.mu state.k source
      (round state.mu[state.k][source]) source 0 hkMu hsMuSize hsMu.le
      (by rw [hvalid.mu_rows_square source hsMuSize]; omega) (by omega)
    rw [hprefix]
    let mu := reduceMuPrefixArray state.mu state.k source
      (round state.mu[state.k][source]) source 0
    have hmuSize : mu.size = state.mu.size := reduceMuPrefixArray_size _ _ _ _ _ _
    have hk' : state.k < mu.size := by rw [hmuSize]; exact hkMu
    have hrow : mu[state.k].size = state.mu[state.k].size := by
      have h := reduceMuPrefixArray_row_size state.mu state.k source state.k
        (round state.mu[state.k][source]) source 0 hkMu
      change mu[state.k]!.size = state.mu[state.k]!.size at h
      rw [getElem!_pos mu state.k hk', getElem!_pos state.mu state.k hkMu] at h
      exact h
    have hs' : source < mu[state.k].size := by rw [hrow]; exact hsMu
    refine ⟨Generated.StrictRecombine.LLLState.mk matrix transform
      (mu.set state.k (mu[state.k].set source (mu[state.k][source] -
        (round state.mu[state.k][source] : QQ)))) state.norms state.k, ?_⟩
    change (if h : state.k < mu.size then
      if h' : source < mu[state.k].size then _ else _ else _) = _
    rw [dif_pos hk', dif_pos hs']

end Refinement.StrictRecombine
