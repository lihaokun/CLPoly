import CLPoly.Refinement.VanHoeijSufficientPrecisionValidation
import CLPoly.Refinement.VanHoeijPrecisionThreshold
import CLPoly.Refinement.VanHoeijColumnTrace

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- A normally returning actual full-width round must accept above the
input-fixed threshold. The canonical basis, cutoff and validation targets
are derived from that round's own execution trace. -/
theorem columnRound_full_width_accepts
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget lifted)
    (hindexFits : lifted.1.size ≤ 2 ^ 31)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (original truncated : Array SparsePolyZZ)
    (data : CldSeparationData f lifted.1 original modulus)
    (height : Int) (count width : Int32)
    (parameters : Generated.StrictVanHoeijParameters.Parameters)
    (hcountValue : count.toInt = (lifted.1.size : Int))
    (hselected : Generated.StrictVanHoeijParameters.selectParameters count height (modulus : Int) = some parameters)
    (htruncate : Generated.StrictVanHoeij.truncateCld original parameters.divisor = .ok truncated)
    (hwidth : width.toInt.toNat = Generated.StrictRecombine.cldSpiralWidth truncated)
    (countSize : count.toInt.toNat = truncated.size)
    (modulusNonzero : parameters.modulus ≠ 0)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hsourceCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hsourcePrimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hsourcePositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hsourceFits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hprecision : vanHoeijInputThreshold (SparsePolyZZ.toPoly f) count height < (modulus : Int))
    (before after : Generated.StrictVanHoeij.VanHoeijResult) (accepted : Bool)
    (hrun : Generated.StrictVanHoeijController.columnRound concreteVanHoeijControllerExecution
      f lifted.1 truncated (modulus : Int) parameters.modulus 32 parameters.residual
      count width countSize (by decide) modulusNonzero before = .ok (accepted, after)) :
    accepted = true := by
  obtain ⟨initial, built, hi, hb, ready, reduced, groups, hl, hg, hv⟩ :=
    columnRound_trace concreteVanHoeijControllerExecution f lifted.1 truncated (modulus : Int)
      parameters.modulus 32 parameters.residual count width countSize (by decide)
      modulusNonzero before after accepted hrun
  let bound := count.toInt * 32 * 32 + (Int32.ofNat built.2).toInt * parameters.residual * parameters.residual
  change Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution built.1 ready bound = .ok reduced at hl
  have hcountNat : count.toInt.toNat = lifted.1.size := by rw [hcountValue]; simp
  have hcountNonneg : 0 ≤ count.toInt := by rw [hcountValue]; positivity
  have hheightNonneg : 0 ≤ height := by rw [hheight]; positivity
  obtain ⟨selected, hs, hresultant, hcoordinates⟩ := selected_parameters_above_input_threshold
    (SparsePolyZZ.toPoly f) count height (modulus : Int) hcountNonneg hheightNonneg hprecision
  have heq : selected = parameters := Option.some.inj (hs.symm.trans hselected)
  subst selected
  rw [hcountNat] at hresultant hcoordinates
  obtain ⟨hh, hd, _, hmod, hres, hm, _⟩ := selectParameters_spec count height (modulus : Int) parameters hselected
  have hinit : Generated.StrictRecombine.makeInitialMatrix truncated.size 32 = .ok initial := by
    simpa only [countSize] using hi
  have hbuild : Generated.StrictVanHoeij.buildCldMatrix initial truncated 0
      (Generated.StrictRecombine.cldSpiralWidth truncated) ((modulus : Int) / parameters.divisor) = .ok built := by
    simpa only [hwidth, hmod] using hb
  have hadded := buildCldMatrix_full_width truncated ((modulus : Int) / parameters.divisor) initial built hbuild
  have haddedInt := buildCldMatrix_added_int32 truncated 0 width parameters.modulus initial built hb
  have hcanonical : CanonicalCldBasis truncated (Generated.StrictRecombine.cldSpiralWidth truncated)
      32 ((modulus : Int) / parameters.divisor) (Generated.StrictRecombine.cldSpiralWidth truncated) built.1 := by
    have hc := buildCldMatrix_canonical truncated 0 (Generated.StrictRecombine.cldSpiralWidth truncated)
      32 ((modulus : Int) / parameters.divisor) initial built
      (makeInitialMatrix_canonical truncated _ 32 _ initial (by decide) hinit) hbuild
    simpa only [Nat.zero_add, hadded] using hc
  have hdimensions := data.rebuild_size_bounds truncated parameters.divisor 32
    (Generated.StrictRecombine.cldSpiralWidth truncated) initial built.1 built.2
    htruncate (by decide) hinit le_rfl hbuild
  have hphysical : bound = (lifted.1.size : Int) * 32 * 32 +
      (built.2 : Int) * parameters.residual * parameters.residual := by
    dsimp only [bound]
    rw [hcountValue, haddedInt]
  obtain ⟨hB, hbound⟩ := selected_cutoff_le_separationCutoff count lifted.1.size
    (SparsePolyZZ.toPoly f).natDegree built.2 height (modulus : Int) parameters
    hcountValue hheightNonneg hselected hdimensions.1
  rw [← hphysical] at hB hbound
  have hcutoff : (lifted.1.size : Int) * 32 * 32 +
      (Generated.StrictRecombine.cldSpiralWidth truncated : Int) *
      ((height + parameters.divisor - 1) / parameters.divisor + 2 * (lifted.1.size : Int)) *
      ((height + parameters.divisor - 1) / parameters.divisor + 2 * (lifted.1.size : Int)) ≤ bound := by
    rw [hphysical, hadded, hres, hcountValue]
  obtain ⟨returned, result, he, hsucceed⟩ := hensel_lll_validation_complete hcount hfactors
    hleading hselection hentry hindexFits modulus hmodulus original truncated data
    32 parameters.divisor height bound built.1 htruncate hcanonical (by decide)
    (by rw [← hmod]; exact hm.le) ready hB hdimensions.2 hbound reduced hl hheight hh hcutoff
    hsourceCanonical hsourcePrimitive hsourcePositive hsourceFits hresultant hcoordinates
    { before with cldColumns := Int32.ofNat built.2, retainedRows := Int32.ofNat (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).size }
  have hgroups : returned = groups := by
    apply Except.ok.inj
    exact he.symm.trans (by simpa only [hcountNat] using hg)
  subst returned
  exact (Prod.mk.inj (Except.ok.inj (hv.symm.trans hsucceed))).1

end Refinement.StrictRecombine
