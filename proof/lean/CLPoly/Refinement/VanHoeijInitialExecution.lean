import CLPoly.Refinement.Recombine

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine

theorem setInitialDiagonalLoop_succeeds (scale : ZZ) (size index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hsize : matrix.size = size)
    (hrows : ∀ row (hr : row < matrix.size), matrix[row].size = size) :
    ∃ output, setInitialDiagonalLoop scale size index matrix = .ok output := by
  induction hmeasure : size - index using Nat.strong_induction_on generalizing index matrix with
  | h measure ih =>
    rw [setInitialDiagonalLoop]
    by_cases hi : index < size
    · have hr : index < matrix.size := by omega
      have hc : index < matrix[index].size := by rw [hrows index hr]; exact hi
      rw [dif_pos hi, dif_pos hr, dif_pos hc]
      apply ih (size - (index + 1)) (by omega) (index + 1)
        (matrix.set index (matrix[index].set index scale)) (by simpa using hsize) ?_ rfl
      intro row hrow
      have hrow' : row < matrix.size := by simpa using hrow
      rw [Array.getElem_set]
      split
      · simp only [Array.size_set, hrows index hr]
      · exact hrows row hrow'
    · exact ⟨matrix, by rw [dif_neg hi]⟩

theorem makeInitialMatrix_succeeds (size : Nat) (scale : ZZ) :
    ∃ output, makeInitialMatrix size scale = .ok output := by
  apply setInitialDiagonalLoop_succeeds scale size 0 (zeroMatrix size size)
    (zeroMatrix_size size size)
  intro row hr
  have hrow : row < size := by simpa only [zeroMatrix_size] using hr
  have hh := zeroMatrix_get size size row hrow
  rw [getElem!_pos (zeroMatrix size size) row hr] at hh
  rw [hh, zeroMatrixRow_size]

theorem initializeLLL_succeeds (matrix : Generated.StrictRecombine.LLLMatrix) (hnonempty : 0 < matrix.size)
    (hrows : ∀ row (hr : row < matrix.size), matrix[row].size = matrix.size) :
    ∃ output, initializeLLL matrix = .ok output := by
  obtain ⟨transform, ht⟩ := makeInitialMatrix_succeeds matrix.size 1
  have hdot := dotRows_eq_fin_sum matrix[0] matrix[0] le_rfl
  let dot : ZZ := ∑ k : Fin matrix[0].size, matrix[0][k.val] * matrix[0][k.val]
  let mu := zeroQQMatrix matrix.size matrix.size
  let norms := (Array.replicate matrix.size (0 : QQ)).set 0 (dot : QQ)
    (by simpa using hnonempty)
  have hshape : GramStorageShape mu norms matrix.size := by
    refine ⟨zeroQQMatrix_size _ _, ?_, by simp [norms]⟩
    intro row hr
    have hh := zeroQQMatrix_rows matrix.size matrix.size row (by simpa [mu, zeroQQMatrix_size] using hr)
    simpa only [mu, getElem!_pos (zeroQQMatrix matrix.size matrix.size) row hr] using hh
  obtain ⟨output, ho, _⟩ := initializeGramSchmidtLoop_succeeds matrix 1 matrix.size mu norms
    rfl hrows hshape (by omega)
  refine ⟨(output.1, output.2, transform), ?_⟩
  rw [initializeLLL, dif_pos hnonempty, ht, hdot]
  change Except.map _ (initializeGramSchmidtLoop matrix 1 mu norms) = _
  rw [ho]
  rfl

end Refinement.StrictRecombine
