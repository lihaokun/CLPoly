import CLPoly.Refinement.VanHoeijGroupingExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Validation does not publish a successful prefix when a later group fails. -/
theorem validatePartition_failure_preserves
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (false, after)) : after = before := by
  unfold Generated.StrictVanHoeij.validatePartition at hrun
  split at hrun
  · exact (Prod.mk.inj (Except.ok.inj hrun)).2.symm
  · cases hloop : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups 0 f #[] with
    | error fault => simp [hloop, bind, Except.bind] at hrun
    | ok result =>
      cases result with
      | none =>
        simp only [hloop, bind, Except.bind] at hrun
        exact (Prod.mk.inj (Except.ok.inj hrun)).2.symm
      | some factors => simp [hloop, bind, Except.bind] at hrun

theorem validatePartition_success_frame
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    after.complete = true ∧ after.cldColumns = before.cldColumns ∧
      after.retainedRows = before.retainedRows := by
  unfold Generated.StrictVanHoeij.validatePartition at hrun
  split at hrun
  · simp at hrun
  · cases hloop : Generated.StrictVanHoeij.validatePartitionLoop lifted modulus groups 0 f #[] with
    | error fault => simp [hloop, bind, Except.bind] at hrun
    | ok result =>
      cases result with
      | none => simp [hloop, bind, Except.bind] at hrun
      | some factors =>
        simp only [hloop, bind, Except.bind] at hrun
        have heq := (Prod.mk.inj (Except.ok.inj hrun)).2
        rw [← heq]
        exact ⟨rfl, rfl, rfl⟩

end Refinement.StrictRecombine
