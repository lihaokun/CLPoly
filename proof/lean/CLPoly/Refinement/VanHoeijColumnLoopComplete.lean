import CLPoly.Refinement.VanHoeijColumnLoop

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij
open Generated.StrictVanHoeijController

theorem columnRound_accepted_complete (execution : ControllerExecution)
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus modulus scale residual : ZZ) (count target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (before after : VanHoeijResult)
    (hrun : columnRound execution f lifted cld originalModulus modulus scale residual
      count target countSize scaleNonzero modulusNonzero before = .ok (true, after)) :
    after.complete = true := by
  obtain ⟨initial, built, hi, hb, ready, reduced, candidates, hl, hc, hv⟩ :=
    columnRound_trace execution f lifted cld originalModulus modulus scale residual
      count target countSize scaleNonzero modulusNonzero before after true hrun
  exact (validatePartition_success_frame f lifted originalModulus candidates _ after hv).1

/-- The actual column-growth loop cannot normally exit incomplete if its
full-width round cannot normally reject. This is instantiated with the
sufficient-precision theorem, not retained as an algorithm-output oracle. -/
theorem columnLoop_complete_of_full_width (execution : ControllerExecution)
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus modulus scale residual : ZZ) (count width target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (targetBound : target.toInt ≤ width.toInt) (hpositive : 0 < target)
    (hfull : ∀ before after,
      columnRound execution f lifted cld originalModulus modulus scale residual
        count width countSize scaleNonzero modulusNonzero before = .ok (false, after) → False)
    (before after : VanHoeijResult)
    (hrun : columnLoop execution f lifted cld originalModulus modulus scale residual
      count width target countSize scaleNonzero modulusNonzero targetBound before = .ok after) :
    after.complete = true := by
  induction hmeasure : (width.toInt - target.toInt).toNat using Nat.strong_induction_on
      generalizing target before with
  | h measure ih =>
    rw [columnLoop, dif_pos hpositive] at hrun
    cases hr : columnRound execution f lifted cld originalModulus modulus scale residual
        count target countSize scaleNonzero modulusNonzero before with
    | error fault => simp [hr] at hrun
    | ok result =>
      obtain ⟨accepted, next⟩ := result
      cases accepted with
      | true =>
        have he : next = after := by simpa [hr] using hrun
        subst next
        exact columnRound_accepted_complete execution f lifted cld originalModulus modulus
          scale residual count target countSize scaleNonzero modulusNonzero before after hr
      | false =>
        simp only [hr, Bool.false_eq_true, ↓reduceIte] at hrun
        by_cases he : target = width
        · exact (hfull before next (by simpa only [he] using hr)).elim
        · rw [dif_neg he] at hrun
          have hp := nextCldTarget_progress target width hpositive
            (column_target_lt target width he targetBound)
          have hn : 0 < nextCldTarget target width := by
            apply Int32.lt_iff_toInt_lt.mpr
            have ht := Int32.lt_iff_toInt_lt.mp hpositive
            change (0 : Int) < (nextCldTarget target width).toInt
            change (0 : Int) < target.toInt at ht
            omega
          have hd := execution.targetDecreases target width hpositive he targetBound
          rw [hmeasure] at hd
          exact ih _ hd (nextCldTarget target width)
            (execution.targetBound target width hpositive he targetBound) hn next hrun rfl

end Refinement.StrictRecombine
