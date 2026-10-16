import CLPoly.Generated.StrictVanHoeijRetry

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictHensel
open Generated.StrictVanHoeijEntry Generated.StrictVanHoeijRetry

/-- Every returned decomposition comes from a complete actual van Hoeij
certificate. The controller contains no terminal recombination fallback. -/
def RetrySuccessOrigin (execution : RetryExecution) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) (hp : 2 ≤ p.toNat)
    (output : Array SparsePolyZZ) : Prop :=
  ∃ target current,
    ∃ (hl : __hensel_lift_upoly_raw_ir execution.step execution.tree f factors p target hp = .ok current),
    ∃ after,
    tryRecombine execution.entry f current.1 current.2
      (execution.liftedCountFits f factors p target hp current hl) = .ok after ∧
    after.complete = true ∧ after.factors = output

theorem retryLoop_success_origin (execution : RetryExecution) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) (hp : 2 ≤ p.toNat)
    (full precision : Int32) (output : Array SparsePolyZZ)
    (hrun : retryLoop execution f factors p hp full precision = .ok output) :
    RetrySuccessOrigin execution f factors p hp output := by
  induction hmeasure : (2147483647 - precision.toInt).toNat using Nat.strong_induction_on
      generalizing precision with
  | h measure ih =>
    rw [retryLoop] at hrun
    split at hrun
    next fault hl => contradiction
    next current hl =>
      split at hrun
      next fault hr => contradiction
      next after hr =>
        split at hrun
        next hc => exact ⟨precision, current, hl, after, hr, hc, Except.ok.inj hrun⟩
        next hc =>
          split at hrun
          next fault hn => contradiction
          next next hn =>
            have hg := nextPrecision_increases precision full next
              (execution.precisionNonnegative f factors p precision hp current hl) hn
            have hu := Int32.toInt_lt next
            exact ih _ (by omega) next hrun rfl

theorem lllFactorize_success_origin (execution : RetryExecution) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) (hp : 2 ≤ p.toNat)
    (output : Array SparsePolyZZ)
    (hrun : lllFactorize execution f factors p hp = .ok output) :
    RetrySuccessOrigin execution f factors p hp output := by
  unfold lllFactorize at hrun
  cases hh : Generated.StrictVanHoeijHeuristic.startingPrecision f (Int32.ofNat factors.size) p with
  | error fault => simp only [hh, bind, Except.bind] at hrun; contradiction
  | ok precision =>
    simp only [hh, bind, Except.bind] at hrun
    exact retryLoop_success_origin execution f factors p hp precision.2 precision.1 output hrun

end Refinement.StrictRecombine
