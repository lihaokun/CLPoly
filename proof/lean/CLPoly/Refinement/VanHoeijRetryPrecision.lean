import CLPoly.Refinement.VanHoeijRetryTrace
import CLPoly.Refinement.Hensel

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictVanHoeijRetry

theorem retryLoop_success_nonnegative (execution : RetryExecution) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) (hp : 2 ≤ p.toNat)
    (full precision : Int32) (output : Array SparsePolyZZ)
    (hrun : retryLoop execution f factors p hp full precision = .ok output) :
    0 ≤ precision := by
  rw [retryLoop] at hrun
  split at hrun
  next fault hl => contradiction
  next current hl => exact execution.precisionNonnegative f factors p precision hp current hl

/-- The actual quadratic Hensel lift reaches at least the positive requested
prime power; it may overshoot, so equality is deliberately not asserted. -/
theorem hensel_output_modulus_ge_requested
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {target : Int32} {output : Array SparsePolyZZ × ZZ}
    (hpositive : 0 < target)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p target output) :
    ((p.toNat ^ target.toNatClampNeg : Nat) : Int) ≤ output.2 := by
  have hn : ∀ value, StrictHensel.HenselLiftTargetCorrect f p target value → 0 ≤ value := by
    intro value hv
    cases hv with
    | mignotte leading hz hl =>
      subst target
      exact False.elim (by simpa using hpositive)
    | explicit hp =>
      have hone := Nat.one_le_pow target.toNatClampNeg p.toNat (Fact.out : Nat.Prime p.toNat).pos
      change 0 ≤ ((p.toNat ^ target.toNatClampNeg : Nat) : Int) - 1
      exact_mod_cast (by omega : (0 : Int) ≤ (p.toNat ^ target.toNatClampNeg : Nat) - (1 : Int))
  obtain ⟨value, hv, hm⟩ := hentry.outputModulus_gt_target hn
  cases hv with
  | mignotte leading hz hl =>
    subst target
    exact False.elim (by simpa using hpositive)
  | explicit hp =>
    change ((p.toNat ^ target.toNatClampNeg : Nat) : Int) - 1 < output.2 at hm
    omega

end Refinement.StrictRecombine
