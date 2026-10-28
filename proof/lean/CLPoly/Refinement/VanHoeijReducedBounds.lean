import CLPoly.Refinement.VanHoeijCertificate

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Exact size reduction and Lovász conditions for the same executed state
whose matrix and GS norms are exported by the repaired helper. -/
theorem lllReduceWithNorms_reduced_state
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    ∃ state : Generated.StrictRecombine.LLLState,
      ConcreteLLLExecutionValid state ∧ state.matrix = output.1 ∧
      state.norms = output.2.2.2 ∧
      (∀ row, row < state.matrix.size → ∀ column, column < row →
        |((state.mu[row]!)[column]!)| ≤ (1 : QQ) / 2) ∧
      (∀ index, 0 < index → index < state.matrix.size →
        ((3 : QQ) / 4 - ((state.mu[index]!)[index - 1]!) ^ 2) *
          state.norms[index - 1]! ≤ state.norms[index]!) := by
  rw [Generated.StrictVanHoeij.lllReduceWithNorms] at hrun
  split at hrun
  next fault hinitialize => contradiction
  next mu norms transform hinitialize =>
    dsimp only at hrun
    split at hrun
    next fault hloop => contradiction
    next state hloop =>
      split at hrun
      next fault hcollect => contradiction
      next rows hcollect =>
        have hvalid := initializeLLL_concrete_valid matrix mu norms transform hinput hinitialize
        let initial := Generated.StrictRecombine.LLLState.mk matrix transform mu norms 1
        have hsizePrefix : SizeReducedPrefix initial := by
          intro row hrow _ column hcolumn
          change row < 1 at hrow
          omega
        have hlovaszPrefix : LovaszPrefix initial := by
          intro index hpositive hindex _
          change index < 1 at hindex
          omega
        have hfinal := concreteLLLMainLoop_preserves_execution_valid
          initial state hvalid hloop
        have hsize := concreteLLLMainLoop_size_reduced_all
          initial state hvalid hsizePrefix hloop
        have hlovasz := concreteLLLMainLoop_lovasz_all
          initial state hvalid hlovaszPrefix hloop
        cases Except.ok.inj hrun
        refine ⟨state, hfinal, rfl, rfl, hsize, ?_⟩
        simpa only [pow_two] using hlovasz

/-- Each preceding GS norm is at most twice the next one in the actual
returned array. This is the quantitative fact needed by eventual separation. -/
theorem lllReduceWithNorms_norms_step_bound
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output) :
    ∀ index, 0 < index → index < output.2.2.2.size →
      output.2.2.2[index - 1]! ≤ 2 * output.2.2.2[index]! := by
  obtain ⟨state, hvalid, _, hnorms, hsize, hlovasz⟩ :=
    lllReduceWithNorms_reduced_state matrix hinput bound output hrun
  rw [← hnorms]
  intro index hpos hindex
  have hm : index < state.matrix.size := by rwa [← hvalid.norms_size]
  have hp : index - 1 < state.norms.size := by omega
  have hmu := hsize index hm (index - 1) (by omega)
  have hl := hlovasz index hpos hm
  have hn : 0 ≤ state.norms[index - 1]! := by
    simpa only [getElem!_pos state.norms (index - 1) hp] using
      le_of_lt (hvalid.norms_positive (index - 1) hp)
  have hsquare : ((state.mu[index]!)[index - 1]!) ^ 2 ≤ (1 : QQ) / 4 := by
    obtain ⟨hlo, hhi⟩ := abs_le.mp hmu
    nlinarith
  have hproduct := mul_le_mul_of_nonneg_right hsquare hn
  nlinarith

/-- A nonempty retained prefix ends precisely at a norm below the cutoff. -/
theorem retainedPrefixCount_last_le
    (norms : Array QQ) (bound : QQ) (remaining : Nat)
    (hsize : remaining ≤ norms.size)
    (hpositive : 0 < Generated.StrictVanHoeij.retainedPrefixCount norms bound remaining) :
    norms[Generated.StrictVanHoeij.retainedPrefixCount norms bound remaining - 1]! ≤ bound := by
  induction remaining with
  | zero => simp [Generated.StrictVanHoeij.retainedPrefixCount] at hpositive
  | succ remaining ih =>
    have hi : remaining < norms.size := by omega
    rw [Generated.StrictVanHoeij.retainedPrefixCount]
    simp only [Nat.succ_pos, ↓reduceDIte, Nat.add_sub_cancel, hi]
    split
    next hlarge =>
      apply ih (by omega)
      rw [Generated.StrictVanHoeij.retainedPrefixCount] at hpositive
      simpa only [Nat.succ_pos, ↓reduceDIte, Nat.add_sub_cancel, hi, hlarge,
        ↓reduceIte] using hpositive
    next hsmall =>
      simpa only [Nat.add_sub_cancel, getElem!_pos norms remaining hi] using
        le_of_not_gt hsmall

/-- Iterate the adjacent GS inequality without assuming monotonicity of the
norms; LLL only guarantees the factor-two comparison. -/
theorem gs_norm_step_bound_iterate
    (norms : Array QQ)
    (hstep : ∀ index, 0 < index → index < norms.size →
      norms[index - 1]! ≤ 2 * norms[index]!)
    (index distance : Nat) (hsize : index + distance < norms.size) :
    norms[index]! ≤ (2 : QQ) ^ distance * norms[index + distance]! := by
  induction distance with
  | zero => simp
  | succ distance ih =>
    have hprevious := ih (by omega)
    have hnext := hstep (index + distance + 1) (by omega) (by omega)
    simp only [Nat.add_sub_cancel] at hnext
    calc
      norms[index]! ≤ (2 : QQ) ^ distance * norms[index + distance]! := hprevious
      _ ≤ (2 : QQ) ^ distance * (2 * norms[index + distance + 1]!) :=
        mul_le_mul_of_nonneg_left hnext (by positivity)
      _ = _ := by rw [pow_succ]; simp only [Nat.add_assoc]; ring

/-- Every actually retained GS norm is bounded in terms of the cutoff and
the retained dimension, independently of the input modulus. -/
theorem lllReduceWithNorms_retained_gs_bound
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (index : Nat)
    (hindex : index < Generated.StrictVanHoeij.retainedPrefixCount
      output.2.2.2 (bound : QQ) output.2.2.2.size) :
    output.2.2.2[index]! ≤ (2 : QQ) ^
      (Generated.StrictVanHoeij.retainedPrefixCount output.2.2.2
        (bound : QQ) output.2.2.2.size - 1 - index) * (bound : QQ) := by
  let norms := output.2.2.2
  let count := Generated.StrictVanHoeij.retainedPrefixCount norms (bound : QQ) norms.size
  have hcount := (retainedPrefixCount_spec norms (bound : QQ) norms.size le_rfl).1
  have hpositive : 0 < count := by change index < count at hindex; omega
  have hlast := retainedPrefixCount_last_le norms (bound : QQ) norms.size le_rfl hpositive
  have hstep := lllReduceWithNorms_norms_step_bound matrix hinput bound output hrun
  have heq : index + (count - 1 - index) = count - 1 := by
    change index < count at hindex
    omega
  have hlastindex : count - 1 < norms.size := by change count ≤ norms.size at hcount; omega
  have hiter := gs_norm_step_bound_iterate norms hstep index (count - 1 - index)
    (by rw [heq]; exact hlastindex)
  rw [heq] at hiter
  exact hiter.trans (mul_le_mul_of_nonneg_left hlast (by positivity))

end Refinement.StrictRecombine
