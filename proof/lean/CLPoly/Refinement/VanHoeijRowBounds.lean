import CLPoly.Refinement.VanHoeijReducedBounds

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The actual physical row norm, obtained from the diagonal of the cached
Gram factorization. It is not identified with a GS norm. -/
theorem concreteGramSchmidt_row_square_sum
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state) (row : Fin state.matrix.size) :
    (∑ column : Fin state.matrix.size,
      (((state.matrix[row.val]!)[column.val]! : ZZ) : QQ) ^ 2) =
    ∑ j : Fin state.matrix.size,
      (gsLowerPrefix state state.matrix.size row j) ^ 2 * state.norms[j.val]! := by
  have hg := congrFun (congrFun (hvalid.gram_schmidt state.matrix.size le_rfl) row) row
  rw [Matrix.mul_apply] at hg
  simp only [gsNormDiagonal, Matrix.mul_diagonal, Matrix.transpose_apply] at hg
  simp only [gramPrefixMatrixQQ, gramPrefixMatrix, Matrix.map_apply,
    Int.cast_sum, Int.cast_mul] at hg
  convert hg using 1
  · apply Finset.sum_congr rfl
    intro j _
    ring
  · apply Finset.sum_congr rfl
    intro j _
    ring

/-- A uniform bound for earlier GS norms gives a physical row-norm bound
via the actual size-reduced lower-triangular coefficients. -/
theorem concrete_reduced_row_square_sum_le
    (state : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hsize : ∀ row, row < state.matrix.size → ∀ column, column < row →
      |((state.mu[row]!)[column]!)| ≤ (1 : QQ) / 2)
    (count : Nat) (S : QQ) (hS : 0 ≤ S)
    (hnorm : ∀ index, index < state.matrix.size → index < count → state.norms[index]! ≤ S)
    (row : Fin state.matrix.size) (hrow : row.val < count) :
    (∑ column : Fin state.matrix.size,
      (((state.matrix[row.val]!)[column.val]! : ZZ) : QQ) ^ 2) ≤
      (state.matrix.size : QQ) * S := by
  rw [concreteGramSchmidt_row_square_sum state hvalid row]
  calc
    (∑ j : Fin state.matrix.size,
        (gsLowerPrefix state state.matrix.size row j) ^ 2 * state.norms[j.val]!) ≤
        ∑ _j : Fin state.matrix.size, S := by
      apply Finset.sum_le_sum
      intro j _
      by_cases hj : j.val ≤ row.val
      · have hn := hnorm j.val j.isLt (by omega)
        have hp : 0 ≤ state.norms[j.val]! := by
          have hjs : j.val < state.norms.size := by
            rw [hvalid.norms_size]
            exact j.isLt
          simpa only [getElem!_pos state.norms j.val hjs] using
            le_of_lt (hvalid.norms_positive j.val hjs)
        have hl : (gsLowerPrefix state state.matrix.size row j) ^ 2 ≤ 1 := by
          by_cases hjlt : j.val < row.val
          · simp only [gsLowerPrefix, hjlt, ↓reduceIte]
            obtain ⟨hlo, hhi⟩ := abs_le.mp (hsize row.val row.isLt j.val hjlt)
            nlinarith
          · have heq : row = j := Fin.ext (by omega)
            simp [gsLowerPrefix, heq]
        exact (mul_le_mul_of_nonneg_right hl hp).trans (by simpa using hn)
      · have hjlt : ¬j.val < row.val := by omega
        have hne : row ≠ j := by intro heq; have := congrArg Fin.val heq; omega
        simpa only [gsLowerPrefix, hjlt, hne, ↓reduceIte, zero_pow (by decide : 2 ≠ 0),
          zero_mul] using hS
    _ = _ := by simp

/-- A modulus-independent norm bound for each physical retained row of the
actual LLL output; `bound` is the source's squared short-vector cutoff. -/
theorem lllReduceWithNorms_retained_row_square_sum_le
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (hB : 0 ≤ bound)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (row : Nat) (hrow : row < output.1.size)
    (hkeep : row < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    (∑ column : Fin output.1.size,
      (((output.1[row]!)[column.val]! : ZZ) : QQ) ^ 2) ≤
      (output.1.size : QQ) * (2 : QQ) ^ output.1.size * (bound : QQ) := by
  obtain ⟨state, hv, hm, hn, hsize, _⟩ :=
    lllReduceWithNorms_reduced_state matrix hinput bound output hrun
  let count := Generated.StrictVanHoeij.retainedPrefixCount
    state.norms (bound : QQ) state.norms.size
  have hcount : count ≤ state.matrix.size := by
    have hc := (retainedPrefixCount_spec state.norms (bound : QQ)
      state.norms.size le_rfl).1
    exact hc.trans_eq hv.norms_size
  have hbound : 0 ≤ (bound : QQ) := by exact_mod_cast hB
  have hnorm : ∀ j, j < state.matrix.size → j < count →
      state.norms[j]! ≤ (2 : QQ) ^ state.matrix.size * (bound : QQ) := by
    intro j _ hj
    have hkeepj : j < Generated.StrictVanHoeij.retainedPrefixCount
        output.2.2.2 (bound : QQ) output.2.2.2.size := by
      simpa only [← hn] using hj
    have hgs := lllReduceWithNorms_retained_gs_bound matrix hinput bound output hrun j hkeepj
    rw [← hn] at hgs
    change state.norms[j]! ≤ (2 : QQ) ^ (count - 1 - j) * (bound : QQ) at hgs
    exact hgs.trans (mul_le_mul_of_nonneg_right
      (pow_le_pow_right₀ (by norm_num : (1 : QQ) ≤ 2) (by omega : count - 1 - j ≤ state.matrix.size))
      hbound)
  have hr : row < state.matrix.size := by rwa [hm]
  have hk : row < count := by simpa only [← hn] using hkeep
  have hresult := concrete_reduced_row_square_sum_le state hv hsize count
    ((2 : QQ) ^ state.matrix.size * (bound : QQ)) (by positivity)
    hnorm ⟨row, hr⟩ hk
  rw [← hm]
  simpa only [mul_assoc] using hresult

/-- Every actual retained integer coordinate has an explicit bound in the
dimension and cutoff. No assumption about the factorization output is used. -/
theorem lllReduceWithNorms_retained_coordinate_abs_le
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (hB : 0 ≤ bound)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (row column : Nat) (hrow : row < output.1.size) (hcolumn : column < output.1.size)
    (hkeep : row < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    |(output.1[row]!)[column]!| ≤
      (output.1.size : ZZ) * 2 ^ output.1.size * bound + 1 := by
  have hnorm := lllReduceWithNorms_retained_row_square_sum_le
    matrix hinput bound hB output hrun row hrow hkeep
  have hcoord : (((output.1[row]!)[column]! : ZZ) : QQ) ^ 2 ≤
      (output.1.size : QQ) * (2 : QQ) ^ output.1.size * (bound : QQ) := by
    calc
      _ ≤ ∑ c : Fin output.1.size,
          (((output.1[row]!)[c.val]! : ZZ) : QQ) ^ 2 :=
        Finset.single_le_sum (fun c _ => sq_nonneg _) (Finset.mem_univ ⟨column, hcolumn⟩)
      _ ≤ _ := hnorm
  have hsmall : |(((output.1[row]!)[column]! : ZZ) : QQ)| ≤
      (output.1.size : QQ) * (2 : QQ) ^ output.1.size * (bound : QQ) + 1 := by
    have habs := abs_nonneg (((output.1[row]!)[column]! : ZZ) : QQ)
    have hsquare := sq_nonneg (|(((output.1[row]!)[column]! : ZZ) : QQ)| - 1)
    rw [sub_sq, sq_abs] at hsquare
    nlinarith
  exact_mod_cast hsmall

end Refinement.StrictRecombine
