import CLPoly.Refinement.VanHoeijColumnTrace
import CLPoly.Refinement.VanHoeijValidation

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij
open Generated.StrictVanHoeijController

theorem columnRound_failure_complete (execution : ControllerExecution)
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus modulus scale residual : ZZ) (count target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (before after : VanHoeijResult)
    (hrun : columnRound execution f lifted cld originalModulus modulus scale residual
      count target countSize scaleNonzero modulusNonzero before = .ok (false, after)) :
    after.complete = before.complete := by
  obtain ⟨initial, built, hi, hb, ready, reduced, candidates, hl, hc, hv⟩ :=
    columnRound_trace execution f lifted cld originalModulus modulus scale residual
      count target countSize scaleNonzero modulusNonzero before after false hrun
  have he := validatePartition_failure_preserves f lifted originalModulus candidates _ after hv
  have hc := congrArg VanHoeijResult.complete he
  exact hc

theorem columnLoop_complete_origin (execution : ControllerExecution)
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus modulus scale residual : ZZ) (count width target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (targetBound : target.toInt ≤ width.toInt)
    (before after : VanHoeijResult) (hbefore : before.complete = false)
    (hafter : after.complete = true)
    (hrun : columnLoop execution f lifted cld originalModulus modulus scale residual
      count width target countSize scaleNonzero modulusNonzero targetBound before = .ok after) :
    ∃ acceptedTarget acceptedBefore,
      0 < acceptedTarget ∧ acceptedTarget.toInt ≤ width.toInt ∧
      columnRound execution f lifted cld originalModulus modulus scale residual
        count acceptedTarget countSize scaleNonzero modulusNonzero acceptedBefore =
        .ok (true, after) := by
  induction hmeasure : (width.toInt - target.toInt).toNat using Nat.strong_induction_on
      generalizing target before with
  | h measure ih =>
    rw [columnLoop] at hrun
    split at hrun
    next ht =>
      cases hr : columnRound execution f lifted cld originalModulus modulus scale residual
          count target countSize scaleNonzero modulusNonzero before with
      | error fault => simp [hr] at hrun
      | ok result =>
        rcases result with ⟨accepted, next⟩
        cases accepted with
        | true =>
          have he : next = after := by simpa [hr] using hrun
          subst next
          exact ⟨target, before, ht, targetBound, hr⟩
        | false =>
          have hn : next.complete = false :=
            (columnRound_failure_complete execution f lifted cld originalModulus modulus
              scale residual count target countSize scaleNonzero modulusNonzero before next hr).trans hbefore
          simp only [hr, Bool.false_eq_true, ↓reduceIte] at hrun
          split at hrun
          next he =>
            have heq : next = after := Except.ok.inj hrun
            subst next
            simp [hn] at hafter
          next he =>
            have hd := execution.targetDecreases target width ht he targetBound
            rw [hmeasure] at hd
            exact ih _ hd (nextCldTarget target width)
              (execution.targetBound target width ht he targetBound) next hn hrun rfl
    next ht =>
      have he : before = after := Except.ok.inj hrun
      subst before
      simp [hbefore] at hafter

end Refinement.StrictRecombine
