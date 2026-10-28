import CLPoly.Refinement.VanHoeijRetryExecution
import CLPoly.Refinement.FactorZZ

set_option autoImplicit false

namespace Refinement.StrictRecombine

private theorem sparse_isCoprime_of_prime_eq (p q : UInt64)
    [Fact (Nat.Prime p.toNat)] [Fact (Nat.Prime q.toNat)]
    (hpq : p = q) (left right : SparsePolyZp)
    (hcop : IsCoprime (CLPoly.Math.SparsePolyZp.toPoly q.toNat left)
      (CLPoly.Math.SparsePolyZp.toPoly q.toNat right)) :
    IsCoprime (CLPoly.Math.SparsePolyZp.toPoly p.toNat left)
      (CLPoly.Math.SparsePolyZp.toPoly p.toNat right) := by
  subst q
  exact hcop

/-- Recover the Hensel semantic certificate from the concrete retry callee.
Only recursive execution readiness remains; all source algebra and machine
factor-count checks follow from selection and the actual successful run. -/
theorem concreteRetry_hensel_success
    (this : DenseUPolyZp) (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (target : Int32)
    (hruntime : StrictHensel.HenselLiftRuntimeReadiness this termination mulProvider
      f selection.factors target)
    (output : Array SparsePolyZZ × ZZ)
    (hrun : Generated.StrictHensel.__hensel_lift_upoly_raw_ir
      (StrictHensel.strictHenselRawOps termination)
      (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
      f selection.factors selection.prime target (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output) :
    StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime target output := by
  letI densePrime : Fact (Nat.Prime this._p.toNat) := ⟨by
    simpa only [hprime] using (Fact.out : Nat.Prime selection.prime.toNat)⟩
  have hrun' : Generated.StrictHensel.__hensel_lift_upoly_raw_ir
      (StrictHensel.strictHenselRawOps termination)
      (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
      f selection.factors this._p target (Fact.out : Nat.Prime this._p.toNat).two_le = .ok output := by
    simpa only [hprime] using hrun
  have hcount := StrictHensel.__hensel_lift_upoly_raw_ir_factorCount_of_success
    _ _ _ _ _ _ _ _ hrun'
  have hfits := StrictHensel.__hensel_lift_upoly_raw_ir_factorCountFits_of_success
    _ _ _ _ _ _ _ _ hrun'
  have hexponent := StrictHensel.__hensel_lift_upoly_raw_ir_exponent_of_success
    _ _ _ _ _ _ _ _ hrun'
  have hheadDegree : f[0].1.deg < 2 ^ 63 := by
    rw [sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty] at hdegree
    omega
  have hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff := by
    intro leading hleading
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty, hfront]
  have hinvariant : StrictHensel.HenselLiftEntryInvariant this termination mulProvider
      f selection.factors target := by
    refine
      { sourceLeading := ⟨f[0], Array.getElem?_eq_getElem hnonempty⟩
        exponent := hexponent
        targetNonnegative := ?_
        adjustedCanonical := ?_
        adjustedNonempty := ?_
        factorsMonic := ?_
        adjustedPairwise := ?_
        liftReady := hruntime.liftReady
        normalizeReady := hruntime.normalizeReady }
    · exact StrictFactorZZ.henselLiftTarget_nonnegative hnonempty hheadDegree
    · intro adjusted hadjust
      have hadjust' : StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors selection.prime adjusted := by
        simpa only [hprime] using hadjust
      simpa only [hprime] using hadjust'.canonical hphysical
    · intro adjusted hadjust
      apply StrictFactorZZ.selectionFactors_adjusted_nonempty hphysical hleadingSemantic hselection
      simpa only [hprime] using hadjust
    · intro index hi
      change (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat selection.factors[index]).Monic
      rw [hprime]
      exact StrictFactorZZ.selectionFactors_monic hselection index hi
    · intro adjusted hadjust
      have hadjust' : StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors selection.prime adjusted := by
        simpa only [hprime] using hadjust
      intro i j hi hj hij
      change IsCoprime (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat adjusted[i])
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat adjusted[j])
      exact sparse_isCoprime_of_prime_eq this._p selection.prime hprime _ _
        (StrictFactorZZ.selectionAdjustedFactors_pairwise_isCoprime
          hphysical hleadingSemantic hselection hadjust' i j hi hj hij)
  obtain ⟨actual, hactual, hcorrect⟩ := StrictHensel.__hensel_lift_upoly_raw_ir_refines
    this hcfg termination mulProvider f selection.factors target hinvariant hcount hfits
  have heq : actual = output := Except.ok.inj (hactual.symm.trans hrun')
  subst actual
  simpa only [hprime] using hcorrect

/-- Construct the actual Hensel return from source invariants and the same
recursive execution readiness, without assuming a successful call. -/
theorem concreteRetry_hensel_complete
    (this : DenseUPolyZp) (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (target : Int32)
    (hruntime : StrictHensel.HenselLiftRuntimeReadiness this termination mulProvider
      f selection.factors target)
    (hcount : 2 ≤ selection.factors.size)
    (hfits : selection.factors.size < 2 ^ 31)
    (hpositive : 0 < target) :
    ∃ output, Generated.StrictHensel.__hensel_lift_upoly_raw_ir
      (StrictHensel.strictHenselRawOps termination)
      (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
      f selection.factors selection.prime target (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output ∧
      StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime target output := by
  letI densePrime : Fact (Nat.Prime this._p.toNat) := ⟨by
    simpa only [hprime] using (Fact.out : Nat.Prime selection.prime.toNat)⟩
  have hheadDegree : f[0].1.deg < 2 ^ 63 := by
    rw [sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty] at hdegree
    omega
  have hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff := by
    intro leading hleading
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty, hfront]
  have hinvariant : StrictHensel.HenselLiftEntryInvariant this termination mulProvider
      f selection.factors target := by
    refine
      { sourceLeading := ⟨f[0], Array.getElem?_eq_getElem hnonempty⟩
        exponent := Or.inr hpositive
        targetNonnegative := ?_
        adjustedCanonical := ?_
        adjustedNonempty := ?_
        factorsMonic := ?_
        adjustedPairwise := ?_
        liftReady := hruntime.liftReady
        normalizeReady := hruntime.normalizeReady }
    · exact StrictFactorZZ.henselLiftTarget_nonnegative hnonempty hheadDegree
    · intro adjusted hadjust
      have hadjust' : StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors selection.prime adjusted := by
        simpa only [hprime] using hadjust
      simpa only [hprime] using hadjust'.canonical hphysical
    · intro adjusted hadjust
      apply StrictFactorZZ.selectionFactors_adjusted_nonempty hphysical hleadingSemantic hselection
      simpa only [hprime] using hadjust
    · intro index hi
      change (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat selection.factors[index]).Monic
      rw [hprime]
      exact StrictFactorZZ.selectionFactors_monic hselection index hi
    · intro adjusted hadjust
      have hadjust' : StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors selection.prime adjusted := by
        simpa only [hprime] using hadjust
      intro i j hi hj hij
      change IsCoprime (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat adjusted[i])
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat adjusted[j])
      exact sparse_isCoprime_of_prime_eq this._p selection.prime hprime _ _
        (StrictFactorZZ.selectionAdjustedFactors_pairwise_isCoprime
          hphysical hleadingSemantic hselection hadjust' i j hi hj hij)
  obtain ⟨actual, hactual, hcorrect⟩ := StrictHensel.__hensel_lift_upoly_raw_ir_refines
    this hcfg termination mulProvider f selection.factors target hinvariant hcount hfits
  exact ⟨actual, by simpa only [hprime] using hactual, by simpa only [hprime] using hcorrect⟩

end Refinement.StrictRecombine
