import CLPoly.Refinement.FactorZZ
import CLPoly.Refinement.VanHoeijRetryPrecision

set_option autoImplicit false

namespace Refinement.StrictFactorZZ

/-- Any actual lift requested at or beyond the full heuristic exponent has
the recovery bound. The actual exponent need not equal either heuristic output. -/
theorem heuristic_extended_hensel_liveRecoveryPrecision
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {initial full target : Int32} {output : Array SparsePolyZZ × ZZ}
    (r : Int32)
    (hheuristic : Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir f r p =
      .ok (initial, full))
    (hfullNonnegative : 0 ≤ full)
    (hfull : full.toInt ≤ target.toInt)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p target output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (exponent : Nat)
    (houtput : output.2 = ((p.toNat ^ exponent : Nat) : Int)) :
    LiveRecoveryPrecision (p.toNat ^ exponent) f := by
  by_cases hzero : target = 0
  · have hentryZero : StrictHensel.HenselLiftEntryCorrect termination f factors p 0 output := by
      simpa only [hzero] using hentry
    exact selectionHenselFactors_liveRecoveryPrecision hentryZero hcanonical
      hnonempty hdegree leading hleading exponent houtput
  · have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
      intro hsourceZero
      have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty
      rw [hsourceZero] at hhead
      exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
        (by simpa using hhead.symm)
    obtain ⟨bound, hboundRun, _, _⟩ :=
      StrictRecombine.mignotteBoundRaw_bounds_divisor f hcanonical hnonempty hdegree
        (1 : Polynomial Int) hsourceNe (one_dvd _)
    have hheuristicLarge := heuristic_starting_precision_full_pow_gt
      f r p initial full leading bound hleading hboundRun hheuristic hfullNonnegative
    have htargetNonnegative : ∀ value,
        StrictHensel.HenselLiftTargetCorrect f p target value → 0 ≤ value := by
      intro value htarget
      cases htarget with
      | mignotte sourceLeading htargetZero hsourceLeading =>
          exact False.elim (hzero htargetZero)
      | explicit htargetPositive =>
          have honeNat : 1 ≤ p.toNat ^ target.toNatClampNeg :=
            Nat.one_le_pow _ p.toNat (Fact.out : Nat.Prime p.toNat).pos
          have honeInt : (1 : Int) ≤ ((p.toNat ^ target.toNatClampNeg : Nat) : Int) := by
            exact_mod_cast honeNat
          change 0 ≤ ((p.toNat ^ target.toNatClampNeg : Nat) : Int) - 1
          omega
    obtain ⟨value, htarget, htargetLt⟩ := hentry.outputModulus_gt_target htargetNonnegative
    have hpowerLeOutput : ((p.toNat ^ target.toNatClampNeg : Nat) : Int) ≤ output.2 := by
      cases htarget with
      | mignotte sourceLeading htargetZero hsourceLeading =>
          exact False.elim (hzero htargetZero)
      | explicit htargetPositive =>
          change ((p.toNat ^ target.toNatClampNeg : Nat) : Int) - 1 < output.2 at htargetLt
          omega
    have hpowerLe : p.toNat ^ full.toNatClampNeg ≤ p.toNat ^ exponent := by
      have hmonotone : p.toNat ^ full.toNatClampNeg ≤ p.toNat ^ target.toNatClampNeg :=
        pow_le_pow_right₀ (Nat.Prime.one_le (Fact.out : Nat.Prime p.toNat))
          (Int32.toNatClampNeg_le (Int32.le_iff_toInt_le.mpr hfull))
      have hlast : p.toNat ^ target.toNatClampNeg ≤ p.toNat ^ exponent := by
        rw [houtput] at hpowerLeOutput
        exact_mod_cast hpowerLeOutput
      exact hmonotone.trans hlast
    apply liveRecoveryPrecision_of_generated_mignotte_lt f hcanonical hnonempty
      hdegree leading hleading (p.toNat ^ exponent)
    intro candidateBound hcandidateBoundRun
    have hsame : candidateBound = bound := by
      rw [hboundRun] at hcandidateBoundRun
      exact (Except.ok.inj hcandidateBoundRun).symm
    subst candidateBound
    exact hheuristicLarge.trans_le hpowerLe

end Refinement.StrictFactorZZ
