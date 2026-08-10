import CLPoly.Generated.StrictSquarefreeZp
import CLPoly.Refinement.SquarefreeZp

set_option autoImplicit false

namespace Refinement.StrictSquarefreeGenerated

open CLPoly.Math
open CLPoly.Impl.StrictPolynomialGCDRefinement
open Refinement.StrictSquarefreeZp

/-- Concrete representation conditions for every recursive invocation of the
generated SQF entry. -/
structure EntryInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp) : Prop where
  canonical : SparsePolyZp.Canonical this._p.toNat f
  monic : (SparsePolyZp.toPoly this._p.toNat f).Monic
  nonempty : 0 < f.size
  positive : 0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree
  denseBound : sparseDenseLength f ≤ 2 ^ 63

noncomputable def sourceYunC
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp) :
    Polynomial (ZMod this._p.toNat) :=
  normalize (EuclideanDomain.gcd
    (SparsePolyZp.toPoly this._p.toNat source)
    (Polynomial.derivative (SparsePolyZp.toPoly this._p.toNat source)))

noncomputable def sourceYunTarget
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp) :
    List (Polynomial (ZMod this._p.toNat) × Nat) ×
      Polynomial (ZMod this._p.toNat) :=
  let c := sourceYunC this source
  if hc : c ≠ 0 then
    yunLoop (normalize (SparsePolyZp.toPoly this._p.toNat source /ₘ c))
      c 1 [] hc
  else
    ([], c)

/-- Local representation and machine-word invariant for the generated Yun
state.  Semantic accumulator fields will be added when the loop simulation is
connected; these fields already suffice to execute every raw primitive. -/
structure YunInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) : Prop where
  wCanonical : SparsePolyZp.Canonical this._p.toNat w
  cCanonical : SparsePolyZp.Canonical this._p.toNat c
  wMonic : (SparsePolyZp.toPoly this._p.toNat w).Monic
  cMonic : (SparsePolyZp.toPoly this._p.toNat c).Monic
  wBound : sparseDenseLength w ≤ 2 ^ 63
  cBound : sparseDenseLength c ≤ 2 ^ 63
  multiplicityBudget : multiplicity.toNat +
      Generated.StrictSquarefreeZp.squarefreeMeasure w +
      Generated.StrictSquarefreeZp.squarefreeMeasure c < UInt64.size
  resultCanonical : ∀ item ∈ result.toList,
    SparsePolyZp.Canonical this._p.toNat item.1
  semanticContinuation :
    yunLoop (SparsePolyZp.toPoly this._p.toNat w)
        (SparsePolyZp.toPoly this._p.toNat c) multiplicity.toNat
        (toPolyList result this._p.toNat) cMonic.ne_zero =
      sourceYunTarget this source
  targetDerivativeZero :
    Polynomial.derivative (sourceYunTarget this source).2 = 0
  targetDegreeLeSource :
    (sourceYunTarget this source).2.natDegree ≤
      (SparsePolyZp.toPoly this._p.toNat source).natDegree

/-- Project the sparse result of the checked C++ make-monic entry. -/
def makeMonicRawIR (this : DenseUPolyZp) (f : SparsePolyZp) :
    RawExec SparsePolyZp :=
  match upolyMakeMonicIR this f with
  | .error fault => .error fault
  | .ok (_, monic) => .ok monic

/-- Project the actual sparse output of the checked public raw GCD call.
Missing object output remains a concrete assertion failure. -/
def gcdRawIR
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (left right : SparsePolyZp) : RawExec SparsePolyZp :=
  match yunRawGCDIR this hcfg left right (physical left right) with
  | .error fault => .error fault
  | .ok output =>
    match output.output with
    | none => .error .assertionFailure
    | some result => .ok result

/-- A physical workspace makes the projected raw GCD execute successfully;
the result and its semantics come from that execution. -/
theorem gcdRawIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (left right : SparsePolyZp)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonzero : SparsePolyZp.toPoly this._p.toNat left ≠ 0)
    (hrightNonzero : SparsePolyZp.toPoly this._p.toNat right ≠ 0)
    (hleftBound : sparseDenseLength left ≤ 2 ^ 63)
    (hrightBound : sparseDenseLength right ≤ 2 ^ 63) :
    ∃ result,
      gcdRawIR this hcfg physical left right = .ok result ∧
      SparsePolyZp.Canonical this._p.toNat result ∧
      SparsePolyZp.toPoly this._p.toNat result = normalize
        (EuclideanDomain.gcd (SparsePolyZp.toPoly this._p.toNat left)
          (SparsePolyZp.toPoly this._p.toNat right)) := by
  rcases yunRawGCDWorkspace_succeeds this hcfg left right
      (physical left right) hleftCanonical hrightCanonical hleftNonzero
      hrightNonzero hleftBound hrightBound with
    ⟨output, result, hrun, houtput, hcanonical, hsemantic⟩
  refine ⟨result, ?_, hcanonical, hsemantic⟩
  simp [gcdRawIR, hrun, houtput]

/-- On the monic states used by SQF, the projected generated make-monic call
returns the same concrete sparse array. -/
theorem makeMonicRawIR_eq_of_monic
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat f).Monic) :
    makeMonicRawIR this f = .ok f := by
  have hrun := upolyMakeMonicIR_eq_of_monic this f hcanonical hnonempty
    hmonic
  simp [makeMonicRawIR, hrun]

/-- The generated derivative-zero transition is justified by the outputs of
the actual p-th-root and make-monic executions. -/
theorem derivativeZeroStep
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (f root rootMonic : SparsePolyZp)
    (hinvariant : EntryInvariant this f)
    (hderivative : (derivativeIR this f).isEmpty = true)
    (hroot : extractPthRootIR f = .ok root)
    (hmonic : makeMonicRawIR this root = .ok rootMonic) :
    EntryInvariant this rootMonic ∧
      Generated.StrictSquarefreeZp.squarefreeMeasure rootMonic <
        Generated.StrictSquarefreeZp.squarefreeMeasure f := by
  have hdegreeWord : ∀ term ∈ f.toList, term.1.deg < UInt64.size := by
    intro term hterm
    obtain ⟨index, hindex, htermEq⟩ := List.mem_iff_getElem.mp hterm
    have hi : index < f.size := by simpa using hindex
    have heq : f[index] = term := by
      rw [← Array.getElem_toList hi]
      exact htermEq
    rw [← heq]
    exact lt_trans (Nat.lt_of_lt_of_le
      (sparse_degree_lt_denseLength this._p.toNat f hinvariant.canonical
        index hi) hinvariant.denseBound) (by norm_num [UInt64.size])
  have hderivativeEq : derivativeIR this f = #[] := by
    have hsize : (derivativeIR this f).size = 0 := by
      simpa [Array.isEmpty] using hderivative
    exact Array.size_eq_zero_iff.mp hsize
  have hderivativeSemantic : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat f) = 0 :=
    (derivativeIR_eq_empty_iff this f hcfg hinvariant.canonical
      hdegreeWord).mp hderivativeEq
  rcases sqfDerivativeZeroIR_prepares_recursive_call this hcfg f
      hinvariant.canonical hinvariant.monic hinvariant.nonempty
      hinvariant.positive hinvariant.denseBound hderivativeSemantic with
    ⟨preparedRoot, hpreparedRun, hrootCanonical, hrootSemantic,
      hrootMonic, hpreparedMonic, hmeasure, hrootBound, _⟩
  have hrootEq : root = preparedRoot :=
    Except.ok.inj (hroot.symm.trans hpreparedRun)
  subst root
  have hpreparedProjected : makeMonicRawIR this preparedRoot =
      .ok preparedRoot := by
    simp [makeMonicRawIR, hpreparedMonic]
  have hrootMonicEq : rootMonic = preparedRoot :=
    Except.ok.inj (hmonic.symm.trans hpreparedProjected)
  subst rootMonic
  have hrootPositive :
      0 < (SparsePolyZp.toPoly this._p.toNat preparedRoot).natDegree := by
    have hp := (Fact.out : Nat.Prime this._p.toNat)
    have hexpand := Polynomial.expand_contract this._p.toNat
      hderivativeSemantic hp.ne_zero
    have hdegree :
        (SparsePolyZp.toPoly this._p.toNat f).natDegree =
          (SparsePolyZp.toPoly this._p.toNat preparedRoot).natDegree *
            this._p.toNat := by
      rw [hrootSemantic, ← Polynomial.natDegree_expand
        (R := ZMod this._p.toNat) (p := this._p.toNat), hexpand]
    by_contra hnot
    have hzero :
        (SparsePolyZp.toPoly this._p.toNat preparedRoot).natDegree = 0 :=
      Nat.eq_zero_of_not_pos hnot
    rw [hzero, zero_mul] at hdegree
    have hfPositive := hinvariant.positive
    rw [hdegree] at hfPositive
    omega
  exact ⟨⟨hrootCanonical, hrootMonic,
    sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat preparedRoot
      hrootMonic.ne_zero,
    hrootPositive, hrootBound⟩, by simpa using hmeasure⟩

/-- The nonzero-derivative prefix constructs the initial generated Yun state
from the actual raw GCD and exact-division outputs. -/
theorem yunInitial
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (f derivativeOut c wRaw : SparsePolyZp)
    (hinvariant : EntryInvariant this f)
    (hderivative : (derivativeIR this f).isEmpty = false)
    (hderivativeOut : derivativeOut = derivativeIR this f)
    (hgcd : gcdRawIR this hcfg physical f derivativeOut = .ok c)
    (hdiv : pairVecDivIR this f c = .ok wRaw) :
    YunInvariant this f 1 wRaw.normalization c #[] := by
  subst derivativeOut
  have hdegreeWord : ∀ term ∈ f.toList, term.1.deg < UInt64.size := by
    intro term hterm
    obtain ⟨index, hindex, htermEq⟩ := List.mem_iff_getElem.mp hterm
    have hi : index < f.size := by simpa using hindex
    have heq : f[index] = term := by
      rw [← Array.getElem_toList hi]
      exact htermEq
    rw [← heq]
    exact lt_trans (Nat.lt_of_lt_of_le
      (sparse_degree_lt_denseLength this._p.toNat f hinvariant.canonical
        index hi) hinvariant.denseBound) (by norm_num [UInt64.size])
  have hderivativeSemantic : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat f) ≠ 0 := by
    intro hzero
    have hempty := (derivativeIR_eq_empty_iff this f hcfg
      hinvariant.canonical hdegreeWord).mpr hzero
    have hsize : (derivativeIR this f).size = 0 :=
      Array.size_eq_zero_iff.mpr hempty
    have : (derivativeIR this f).isEmpty = true := by
      simpa [Array.isEmpty] using hsize
    simp [hderivative] at this
  rcases sqfNonzeroDerivativeIR_prepares_yun this hcfg physical f
      hinvariant.canonical hinvariant.monic hinvariant.nonempty
      hinvariant.denseBound hderivativeSemantic with
    ⟨gcdOut, preparedC, preparedW, hrawGCD, houtput, hpreparedDiv,
      hpreparedNorm, hcCanonical, hwCanonical, hcMonic, hwMonic,
      hcBound, hwBound, hbudget, hcSemantic, hwSemantic⟩
  have hprojectedGCD : gcdRawIR this hcfg physical f
      (derivativeIR this f) = .ok preparedC := by
    simp [gcdRawIR, hrawGCD, houtput]
  have hcEq : c = preparedC :=
    Except.ok.inj (hgcd.symm.trans hprojectedGCD)
  subst c
  have hwEq : wRaw = preparedW :=
    Except.ok.inj (hdiv.symm.trans hpreparedDiv)
  subst wRaw
  rw [hpreparedNorm]
  let target := yunLoop
    (SparsePolyZp.toPoly this._p.toNat preparedW)
    (SparsePolyZp.toPoly this._p.toNat preparedC) 1 [] hcMonic.ne_zero
  have htargetDerivative : Polynomial.derivative target.2 = 0 := by
    have hraw := yunLoop_sqf_remainder_derivative_zero
      (SparsePolyZp.toPoly this._p.toNat f) hinvariant.monic.ne_zero
      (by rw [← hcSemantic]; exact hcMonic.ne_zero)
    simpa [target, hcSemantic, hwSemantic] using hraw
  have hcDvdSource : SparsePolyZp.toPoly this._p.toNat preparedC ∣
      SparsePolyZp.toPoly this._p.toNat f := by
    rw [hcSemantic]
    exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left _ _)
  have htargetDegree : target.2.natDegree ≤
      (SparsePolyZp.toPoly this._p.toNat f).natDegree :=
    (yunLoop_c_natDegree_le
      (SparsePolyZp.toPoly this._p.toNat preparedW)
      (SparsePolyZp.toPoly this._p.toNat preparedC) 1 []
      hcMonic.ne_zero).trans
      (Polynomial.natDegree_le_of_dvd hcDvdSource hinvariant.monic.ne_zero)
  have hsourceTarget : sourceYunTarget this f = target := by
    rw [sourceYunTarget, dif_pos (by
      unfold sourceYunC
      rw [← hcSemantic]
      exact hcMonic.ne_zero)]
    simp only [sourceYunC]
    simpa only [← hcSemantic, ← hwSemantic] using (show
      yunLoop (SparsePolyZp.toPoly this._p.toNat preparedW)
        (SparsePolyZp.toPoly this._p.toNat preparedC) 1 []
          hcMonic.ne_zero = target from rfl)
  refine ⟨hwCanonical, hcCanonical, hwMonic, hcMonic, hwBound, hcBound,
    ?_, ?_, ?_, ?_, ?_⟩
  · simpa using hbudget
  · intro item hitem
    simp at hitem
  · rw [hsourceTarget]
    simp [target]
  · rw [hsourceTarget]
    exact htargetDerivative
  · rw [hsourceTarget]
    exact htargetDegree

/-- Common semantic and machine payload of one generated Yun body.  Every
property is derived after identifying the actual raw GCD/division outputs. -/
theorem yunRawStep_prepares
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64))
    (y zRaw cRaw : SparsePolyZp)
    (hinvariant : YunInvariant this source multiplicity w c result)
    (hguard : (!w.isEmpty && get_deg w > 0) = true)
    (hgcd : gcdRawIR this hcfg physical w c = .ok y)
    (hz : pairVecDivIR this w y = .ok zRaw)
    (hc : pairVecDivIR this c y = .ok cRaw) :
    SparsePolyZp.Canonical this._p.toNat y ∧
      SparsePolyZp.Canonical this._p.toNat zRaw ∧
      SparsePolyZp.Canonical this._p.toNat cRaw ∧
      zRaw.normalization = zRaw ∧
      cRaw.normalization = cRaw ∧
      (SparsePolyZp.toPoly this._p.toNat y).Monic ∧
      (SparsePolyZp.toPoly this._p.toNat zRaw).Monic ∧
      (SparsePolyZp.toPoly this._p.toNat cRaw).Monic ∧
      sparseDenseLength y ≤ 2 ^ 63 ∧
      sparseDenseLength cRaw ≤ 2 ^ 63 ∧
      Generated.StrictSquarefreeZp.squarefreeMeasure y +
          Generated.StrictSquarefreeZp.squarefreeMeasure cRaw <
        Generated.StrictSquarefreeZp.squarefreeMeasure w +
          Generated.StrictSquarefreeZp.squarefreeMeasure c ∧
      (multiplicity + 1).toNat = multiplicity.toNat + 1 ∧
      (multiplicity + 1).toNat +
          Generated.StrictSquarefreeZp.squarefreeMeasure y +
          Generated.StrictSquarefreeZp.squarefreeMeasure cRaw <
        UInt64.size ∧
      SparsePolyZp.toPoly this._p.toNat y = normalize
        (EuclideanDomain.gcd (SparsePolyZp.toPoly this._p.toNat w)
          (SparsePolyZp.toPoly this._p.toNat c)) ∧
      SparsePolyZp.toPoly this._p.toNat zRaw =
        SparsePolyZp.toPoly this._p.toNat w /ₘ
          SparsePolyZp.toPoly this._p.toNat y ∧
      SparsePolyZp.toPoly this._p.toNat cRaw =
        SparsePolyZp.toPoly this._p.toNat c /ₘ
          SparsePolyZp.toPoly this._p.toNat y := by
  have hwSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat w
    hinvariant.wMonic.ne_zero
  have hcSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat c
    hinvariant.cMonic.ne_zero
  have hwHeadBound : w[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength w = w[0].1.deg + 1 := by
      simp [sparseDenseLength, hwSize]
    have hwBound := hinvariant.wBound
    rw [hdense] at hwBound
    omega
  have hwPositive := (generated_yun_guard_eq_true_iff this._p.toNat w
    hinvariant.wCanonical hwSize hwHeadBound).mp hguard
  rcases gcdRawIR_refines this hcfg physical w c hinvariant.wCanonical
      hinvariant.cCanonical hinvariant.wMonic.ne_zero
      hinvariant.cMonic.ne_zero hinvariant.wBound hinvariant.cBound with
    ⟨preparedY, hpreparedGCD, hyCanonical, hySemantic⟩
  have hyEq : y = preparedY :=
    Except.ok.inj (hgcd.symm.trans hpreparedGCD)
  subst y
  rcases rawGCDOutput_yun_invariants this._p.toNat w c preparedY
      hinvariant.wMonic.ne_zero hySemantic with
    ⟨hySize, hyMonic, hyDvdW, hyDvdC⟩
  rcases yunPairDivisionsIR_refine this w c preparedY hcfg
      hinvariant.wCanonical hinvariant.cCanonical hyCanonical
      hinvariant.wMonic.ne_zero hySemantic with
    ⟨preparedZ, preparedC, hpreparedZ, hpreparedC, hzCanonical,
      hcCanonical, hzNorm, hcNorm, hzSemantic, hcSemantic⟩
  have hzEq : zRaw = preparedZ :=
    Except.ok.inj (hz.symm.trans hpreparedZ)
  subst zRaw
  have hcEq : cRaw = preparedC :=
    Except.ok.inj (hc.symm.trans hpreparedC)
  subst cRaw
  have hzMonic :
      (SparsePolyZp.toPoly this._p.toNat preparedZ).Monic := by
    rw [hzSemantic]
    exact divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat w)
      (SparsePolyZp.toPoly this._p.toNat preparedY) hinvariant.wMonic
      hyMonic hyDvdW
  have hcMonic :
      (SparsePolyZp.toPoly this._p.toNat preparedC).Monic := by
    rw [hcSemantic]
    exact divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat c)
      (SparsePolyZp.toPoly this._p.toNat preparedY) hinvariant.cMonic
      hyMonic hyDvdC
  have hdecrease := yunNext_generatedMeasure_lt this w c preparedY
    preparedC hinvariant.wCanonical hinvariant.cCanonical hyCanonical
    hcCanonical hwSize hcSize hwPositive hinvariant.cMonic hyMonic hyDvdC
    hcSemantic
  have hnextBounds := yunNext_sparseDenseLength_bounds this w c preparedY
    preparedC hinvariant.wCanonical hinvariant.cCanonical hyCanonical
    hcCanonical hwSize hcSize hinvariant.wMonic hinvariant.cMonic hyMonic
    hyDvdW hyDvdC hcSemantic hinvariant.wBound hinvariant.cBound
  have hnextBudget := yunMultiplicityBudget_step multiplicity
    (Generated.squarefreeMeasure w + Generated.squarefreeMeasure c)
    (Generated.squarefreeMeasure preparedY +
      Generated.squarefreeMeasure preparedC) hdecrease (by
        simpa [Nat.add_assoc] using hinvariant.multiplicityBudget)
  exact ⟨hyCanonical, hzCanonical, hcCanonical, hzNorm, hcNorm, hyMonic,
    hzMonic, hcMonic, hnextBounds.1, hnextBounds.2, by simpa using hdecrease,
    hnextBudget.1, by simpa [Nat.add_assoc] using hnextBudget.2,
    hySemantic, hzSemantic, hcSemantic⟩

theorem yunFactorStep
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64))
    (y zRaw zMonic cRaw : SparsePolyZp)
    (hinvariant : YunInvariant this source multiplicity w c result)
    (hguard : (!w.isEmpty && get_deg w > 0) = true)
    (hgcd : gcdRawIR this hcfg physical w c = .ok y)
    (hz : pairVecDivIR this w y = .ok zRaw)
    (_hzGuard : (!zRaw.normalization.isEmpty &&
      get_deg zRaw.normalization > 0) = true)
    (hmonic : makeMonicRawIR this zRaw.normalization = .ok zMonic)
    (hc : pairVecDivIR this c y = .ok cRaw) :
    YunInvariant this source (multiplicity + 1) y cRaw.normalization
        (result.push (zMonic, multiplicity)) ∧
      Generated.StrictSquarefreeZp.squarefreeMeasure y +
          Generated.StrictSquarefreeZp.squarefreeMeasure cRaw.normalization <
        Generated.StrictSquarefreeZp.squarefreeMeasure w +
          Generated.StrictSquarefreeZp.squarefreeMeasure c := by
  rcases yunRawStep_prepares this hcfg physical source multiplicity w c result y
      zRaw cRaw hinvariant hguard hgcd hz hc with
    ⟨hyCanonical, hzCanonical, hcCanonical, hzNorm, hcNorm, hyMonic,
      hzMonic, hcMonic, hyBound, hcBound, hdecrease, hincrement,
      hbudget, hySemantic, hzSemantic, hcSemantic⟩
  have hzSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat zRaw
    hzMonic.ne_zero
  have hsameMonic : makeMonicRawIR this zRaw.normalization = .ok zRaw := by
    rw [hzNorm]
    exact makeMonicRawIR_eq_of_monic this zRaw hzCanonical hzSize hzMonic
  have hzMonicEq : zMonic = zRaw :=
    Except.ok.inj (hmonic.symm.trans hsameMonic)
  subst zMonic
  have hwSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat w
    hinvariant.wMonic.ne_zero
  have hwHeadBound : w[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength w = w[0].1.deg + 1 := by
      simp [sparseDenseLength, hwSize]
    have hwBound := hinvariant.wBound
    rw [hdense] at hwBound
    omega
  have hwPositive := (generated_yun_guard_eq_true_iff this._p.toNat w
    hinvariant.wCanonical hwSize hwHeadBound).mp hguard
  have hzBound := yunQuotient_sparseDenseLength_bound this w y zRaw
    hinvariant.wCanonical hzCanonical hwSize hzMonic hyMonic hzSemantic
    hinvariant.wBound
  have hzGuard : (!zRaw.isEmpty && get_deg zRaw > 0) = true := by
    simpa [hzNorm] using _hzGuard
  have hnextList := yunNextResult_toPolyList this result zRaw multiplicity
    hzCanonical hzMonic hzBound
  simp only [hzGuard, if_pos] at hnextList
  have hzConcrete : SparsePolyZp.toPoly this._p.toNat zRaw = normalize
      (SparsePolyZp.toPoly this._p.toNat w /ₘ
        SparsePolyZp.toPoly this._p.toNat y) := by
    rw [← hzSemantic, hzMonic.normalize_eq_self]
  have hcConcrete : SparsePolyZp.toPoly this._p.toNat cRaw = normalize
      (SparsePolyZp.toPoly this._p.toNat c /ₘ
        SparsePolyZp.toPoly this._p.toNat y) := by
    rw [← hcSemantic, hcMonic.normalize_eq_self]
  have hloopStep := yunLoop_eq_concrete_step
    (SparsePolyZp.toPoly this._p.toNat w)
    (SparsePolyZp.toPoly this._p.toNat c)
    (SparsePolyZp.toPoly this._p.toNat y)
    (SparsePolyZp.toPoly this._p.toNat zRaw)
    (SparsePolyZp.toPoly this._p.toNat cRaw) multiplicity.toNat
    (toPolyList result this._p.toNat) hinvariant.cMonic.ne_zero
    hcMonic.ne_zero hwPositive hySemantic hzConcrete hcConcrete
  rw [← hnextList, ← hincrement] at hloopStep
  have hcontinuation : yunLoop (SparsePolyZp.toPoly this._p.toNat y)
      (SparsePolyZp.toPoly this._p.toNat cRaw)
      (multiplicity + 1).toNat
      (toPolyList (result.push (zRaw, multiplicity)) this._p.toNat)
      hcMonic.ne_zero =
        sourceYunTarget this source :=
    hloopStep.symm.trans hinvariant.semanticContinuation
  rw [hcNorm]
  refine ⟨⟨hyCanonical, hcCanonical, hyMonic, hcMonic, hyBound, hcBound,
    ?_, ?_, ?_, hinvariant.targetDerivativeZero,
      hinvariant.targetDegreeLeSource⟩, ?_⟩
  · simpa [Nat.add_assoc] using hbudget
  · intro item hitem
    rw [Array.toList_push] at hitem
    rcases List.mem_append.mp hitem with hprevious | hlast
    · exact hinvariant.resultCanonical item hprevious
    · simp at hlast
      subst item
      exact hzCanonical
  · simpa [hcNorm] using hcontinuation
  · simpa using hdecrease

theorem yunNoFactorStep
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64))
    (y zRaw cRaw : SparsePolyZp)
    (hinvariant : YunInvariant this source multiplicity w c result)
    (hguard : (!w.isEmpty && get_deg w > 0) = true)
    (hgcd : gcdRawIR this hcfg physical w c = .ok y)
    (hz : pairVecDivIR this w y = .ok zRaw)
    (_hzGuard : (!zRaw.normalization.isEmpty &&
      get_deg zRaw.normalization > 0) = false)
    (hc : pairVecDivIR this c y = .ok cRaw) :
    YunInvariant this source (multiplicity + 1) y cRaw.normalization result ∧
      Generated.StrictSquarefreeZp.squarefreeMeasure y +
          Generated.StrictSquarefreeZp.squarefreeMeasure cRaw.normalization <
        Generated.StrictSquarefreeZp.squarefreeMeasure w +
          Generated.StrictSquarefreeZp.squarefreeMeasure c := by
  rcases yunRawStep_prepares this hcfg physical source multiplicity w c result y
      zRaw cRaw hinvariant hguard hgcd hz hc with
    ⟨hyCanonical, hzCanonical, hcCanonical, hzNorm, hcNorm, hyMonic,
      hzMonic, hcMonic,
      hyBound, hcBound, hdecrease, hincrement, hbudget, hySemantic,
      hzSemantic, hcSemantic⟩
  have hwSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat w
    hinvariant.wMonic.ne_zero
  have hwHeadBound : w[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength w = w[0].1.deg + 1 := by
      simp [sparseDenseLength, hwSize]
    have hwBound := hinvariant.wBound
    rw [hdense] at hwBound
    omega
  have hwPositive := (generated_yun_guard_eq_true_iff this._p.toNat w
    hinvariant.wCanonical hwSize hwHeadBound).mp hguard
  have hzBound := yunQuotient_sparseDenseLength_bound this w y zRaw
    hinvariant.wCanonical hzCanonical hwSize hzMonic hyMonic hzSemantic
    hinvariant.wBound
  have hzGuard : (!zRaw.isEmpty && get_deg zRaw > 0) = false := by
    simpa [hzNorm] using _hzGuard
  have hnextList := yunNextResult_toPolyList this result zRaw multiplicity
    hzCanonical hzMonic hzBound
  simp only [hzGuard] at hnextList
  have hzConcrete : SparsePolyZp.toPoly this._p.toNat zRaw = normalize
      (SparsePolyZp.toPoly this._p.toNat w /ₘ
        SparsePolyZp.toPoly this._p.toNat y) := by
    rw [← hzSemantic, hzMonic.normalize_eq_self]
  have hcConcrete : SparsePolyZp.toPoly this._p.toNat cRaw = normalize
      (SparsePolyZp.toPoly this._p.toNat c /ₘ
        SparsePolyZp.toPoly this._p.toNat y) := by
    rw [← hcSemantic, hcMonic.normalize_eq_self]
  have hloopStep := yunLoop_eq_concrete_step
    (SparsePolyZp.toPoly this._p.toNat w)
    (SparsePolyZp.toPoly this._p.toNat c)
    (SparsePolyZp.toPoly this._p.toNat y)
    (SparsePolyZp.toPoly this._p.toNat zRaw)
    (SparsePolyZp.toPoly this._p.toNat cRaw) multiplicity.toNat
    (toPolyList result this._p.toNat) hinvariant.cMonic.ne_zero
    hcMonic.ne_zero hwPositive hySemantic hzConcrete hcConcrete
  rw [← hnextList, ← hincrement] at hloopStep
  have hcontinuation : yunLoop (SparsePolyZp.toPoly this._p.toNat y)
      (SparsePolyZp.toPoly this._p.toNat cRaw)
      (multiplicity + 1).toNat (toPolyList result this._p.toNat)
      hcMonic.ne_zero =
        sourceYunTarget this source :=
    hloopStep.symm.trans hinvariant.semanticContinuation
  rw [hcNorm]
  exact ⟨⟨hyCanonical, hcCanonical, hyMonic, hcMonic, hyBound, hcBound,
    by simpa [Nat.add_assoc] using hbudget, hinvariant.resultCanonical,
    by simpa [hcNorm] using hcontinuation,
    hinvariant.targetDerivativeZero, hinvariant.targetDegreeLeSource⟩,
    by simpa using hdecrease⟩

/-- The generated Yun loop can enter the post-loop p-th-root branch only
after its actual C++ loop guard has become false.  The carried continuation
then identifies the concrete remainder with the L2 Yun remainder, so the
recursive SQF call is justified by the real raw root and make-monic runs. -/
theorem remainderRootStep
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (f _derivativeOut : SparsePolyZp)
    (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64))
    (root rootMonic : SparsePolyZp)
    (hentry : EntryInvariant this f)
    (_hderivative : (derivativeIR this f).isEmpty = false)
    (hinvariant : YunInvariant this f multiplicity w c result)
    (hdone : (!w.isEmpty && get_deg w > 0) = false)
    (hcGuard : (!c.isEmpty && get_deg c > 0) = true)
    (hroot : extractPthRootIR c = .ok root)
    (hmonic : makeMonicRawIR this root = .ok rootMonic) :
    EntryInvariant this rootMonic ∧
      Generated.StrictSquarefreeZp.squarefreeMeasure rootMonic <
        Generated.StrictSquarefreeZp.squarefreeMeasure f := by
  have hwSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat w
    hinvariant.wMonic.ne_zero
  have hwHeadBound : w[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength w = w[0].1.deg + 1 := by
      simp [sparseDenseLength, hwSize]
    have hwBound := hinvariant.wBound
    rw [hdense] at hwBound
    omega
  have hwNotPositive : ¬ 0 <
      (SparsePolyZp.toPoly this._p.toNat w).natDegree := by
    intro hwPositive
    have hguardTrue := (generated_yun_guard_eq_true_iff this._p.toNat w
      hinvariant.wCanonical hwSize hwHeadBound).mpr hwPositive
    rw [hdone] at hguardTrue
    contradiction
  have hwDegree : (SparsePolyZp.toPoly this._p.toNat w).natDegree = 0 :=
    Nat.eq_zero_of_not_pos hwNotPositive
  have hcontinuation := hinvariant.semanticContinuation
  rw [yunLoop, dif_pos hwDegree] at hcontinuation
  have hcSemantic : SparsePolyZp.toPoly this._p.toNat c =
      (sourceYunTarget this f).2 := congrArg Prod.snd hcontinuation
  have hcDerivative : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat c) = 0 := by
    rw [hcSemantic]
    exact hinvariant.targetDerivativeZero
  have hcDegree : (SparsePolyZp.toPoly this._p.toNat c).natDegree ≤
      (SparsePolyZp.toPoly this._p.toNat f).natDegree := by
    rw [hcSemantic]
    exact hinvariant.targetDegreeLeSource
  have hcSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat c
    hinvariant.cMonic.ne_zero
  have hcHeadBound : c[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength c = c[0].1.deg + 1 := by
      simp [sparseDenseLength, hcSize]
    have hcBound := hinvariant.cBound
    rw [hdense] at hcBound
    omega
  have hcPositive := (generated_yun_guard_eq_true_iff this._p.toNat c
    hinvariant.cCanonical hcSize hcHeadBound).mp hcGuard
  rcases sqfPostYunIR_prepares_recursive_call this hcfg f c
      hentry.canonical hentry.nonempty hinvariant.cCanonical
      hinvariant.cMonic hcPositive hinvariant.cBound hcDerivative hcDegree with
    ⟨preparedRoot, hpreparedRoot, hrootCanonical, hrootSemantic,
      hrootMonic, hpreparedMonic, hmeasure, hrootBound⟩
  have hrootEq : root = preparedRoot :=
    Except.ok.inj (hroot.symm.trans hpreparedRoot)
  subst root
  have hprojectedMonic : makeMonicRawIR this preparedRoot =
      .ok preparedRoot := by
    simp [makeMonicRawIR, hpreparedMonic]
  have hrootMonicEq : rootMonic = preparedRoot :=
    Except.ok.inj (hmonic.symm.trans hprojectedMonic)
  subst rootMonic
  refine ⟨⟨hrootCanonical, hrootMonic, ?_, ?_, hrootBound⟩, ?_⟩
  · exact sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat preparedRoot
      hrootMonic.ne_zero
  · rw [hrootSemantic]
    have hp := (Fact.out : Nat.Prime this._p.toNat)
    have hexpand : Polynomial.expand (ZMod this._p.toNat) this._p.toNat
        (Polynomial.contract this._p.toNat
          (SparsePolyZp.toPoly this._p.toNat c)) =
          SparsePolyZp.toPoly this._p.toNat c :=
      Polynomial.expand_contract this._p.toNat hcDerivative hp.ne_zero
    have hdegree : (SparsePolyZp.toPoly this._p.toNat c).natDegree =
        (Polynomial.contract this._p.toNat
          (SparsePolyZp.toPoly this._p.toNat c)).natDegree *
            this._p.toNat := by
      conv_lhs => rw [← hexpand]
      rw [Polynomial.natDegree_expand]
    by_contra hnot
    have hzero : (Polynomial.contract this._p.toNat
        (SparsePolyZp.toPoly this._p.toNat c)).natDegree = 0 :=
      Nat.eq_zero_of_not_pos hnot
    rw [hzero, zero_mul] at hdegree
    omega
  · simpa only [Generated.StrictSquarefreeZp.squarefreeMeasure,
      Generated.squarefreeMeasure] using hmeasure

/-- The concrete operation bundle consumed by the generated C++ SQF control
flow.  Every operation is an existing checked L1 implementation; no field
evaluates `sqfZp`, `yunLoop`, or another L2 specification. -/
noncomputable def strictSQFRawOps
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg) :
    Generated.StrictSquarefreeZp.SQFRawOps where
  derivative := derivativeIR this
  extractPthRoot := extractPthRootIR
  makeMonic := makeMonicRawIR this
  gcd := gcdRawIR this hcfg physical
  exactDiv := pairVecDivIR this
  EntryInvariant := EntryInvariant this
  YunInvariant := YunInvariant this
  derivativeZeroStep := derivativeZeroStep this hcfg
  yunInitial := fun f derivativeOut c wRaw hinvariant hderivative
      hderivativeOut hgcd hdiv =>
    yunInitial this hcfg physical f derivativeOut c wRaw hinvariant
      hderivative hderivativeOut hgcd hdiv
  yunFactorStep := yunFactorStep this hcfg physical
  yunNoFactorStep := yunNoFactorStep this hcfg physical
  remainderRootStep := fun f derivativeOut multiplicity w c result root
      rootMonic hentry hderivative hinvariant hdone hcGuard hroot hmonic =>
    remainderRootStep this hcfg f derivativeOut multiplicity w c result root
      rootMonic hentry hderivative hinvariant hdone hcGuard hroot hmonic

theorem sqfCertifyRawExec_ok {α : Type} (run : RawExec α) (output : α)
    (hrun : run = .ok output) :
    Generated.StrictSquarefreeZp.certifyRawExec run =
      .ok ⟨output, hrun⟩ := by
  subst run
  rfl

theorem sqfCertifyBool_eq (value result : Bool)
    (hvalue : value = result) :
    Generated.StrictSquarefreeZp.certifyBool value = ⟨result, hvalue⟩ := by
  subst value
  rfl

/-- The generated well-founded Yun state loop really executes to a terminal
state.  This proof follows every generated raw call and recurses only through
the strict source measure supplied by the transition invariant. -/
theorem generatedYunLoopState_executes
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) :
    ∀ state : Generated.StrictSquarefreeZp.YunRawState
        (strictSQFRawOps this hcfg physical) source,
      ∃ finalState : Generated.StrictSquarefreeZp.YunRawFinalState
          (strictSQFRawOps this hcfg physical) source,
        Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state
            (strictSQFRawOps this hcfg physical) state = .ok finalState := by
  suffices hstrong : ∀ n state,
      Generated.StrictSquarefreeZp.yunStateMeasure state = n →
      ∃ finalState : Generated.StrictSquarefreeZp.YunRawFinalState
          (strictSQFRawOps this hcfg physical) source,
        Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state
            (strictSQFRawOps this hcfg physical) state = .ok finalState by
    intro state
    exact hstrong (Generated.StrictSquarefreeZp.yunStateMeasure state)
      state rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
      intro state hmeasureEq
      rcases state with ⟨multiplicity, w, c, result, hinvariant⟩
      by_cases hguard : (!w.isEmpty && get_deg w > 0) = true
      · rcases gcdRawIR_refines this hcfg physical w c
            hinvariant.wCanonical hinvariant.cCanonical
            hinvariant.wMonic.ne_zero hinvariant.cMonic.ne_zero
            hinvariant.wBound hinvariant.cBound with
          ⟨y, hyRun, hyCanonical, hySemantic⟩
        rcases yunPairDivisionsIR_refine this w c y hcfg
            hinvariant.wCanonical hinvariant.cCanonical hyCanonical
            hinvariant.wMonic.ne_zero hySemantic with
          ⟨zRaw, cRaw, hzRun, hcRun, hzCanonical, hcCanonical,
            hzNorm, hcNorm, _hzSemantic, _hcSemantic⟩
        have hprepared := yunRawStep_prepares this hcfg physical source
          multiplicity w c result y zRaw cRaw hinvariant hguard hyRun hzRun
          hcRun
        rcases hprepared with
          ⟨_hyCanonical, _hzCanonical, _hcCanonical, _hzNorm, _hcNorm,
            _hyMonic, hzMonic, _hcMonic, _hyBound, _hcBound, _hdecrease,
            _hincrement, _hbudget, _hySemantic, _hzSemantic',
            _hcSemantic'⟩
        by_cases hzGuard : (!zRaw.normalization.isEmpty &&
            get_deg zRaw.normalization > 0) = true
        · have hzSize := sparsePolyZp_size_pos_of_toPoly_ne_zero
            this._p.toNat zRaw hzMonic.ne_zero
          have hmakeCore := upolyMakeMonicIR_eq_of_monic this zRaw
            hzCanonical hzSize hzMonic
          have hmake : makeMonicRawIR this zRaw.normalization = .ok zRaw := by
            simp [hzNorm, makeMonicRawIR, hmakeCore]
          have hstep := yunFactorStep this hcfg physical source multiplicity
            w c result y zRaw zRaw cRaw hinvariant hguard hyRun hzRun
            hzGuard hmake hcRun
          let nextState : Generated.StrictSquarefreeZp.YunRawState
              (strictSQFRawOps this hcfg physical) source :=
            ⟨multiplicity + 1, y, cRaw.normalization,
              result.push (zRaw, multiplicity), hstep.1⟩
          rcases ih (Generated.StrictSquarefreeZp.yunStateMeasure nextState)
              (by simpa [nextState,
                Generated.StrictSquarefreeZp.yunStateMeasure] using
                  hstep.2.trans_le (Nat.le_of_eq hmeasureEq))
              nextState rfl with ⟨finalState, hrecursive⟩
          refine ⟨finalState, ?_⟩
          rw [Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state.eq_1]
          rw [sqfCertifyBool_eq _ true hguard]
          simp only [strictSQFRawOps]
          rw [sqfCertifyRawExec_ok _ y hyRun]
          simp only
          rw [sqfCertifyRawExec_ok _ zRaw hzRun]
          simp only
          rw [sqfCertifyBool_eq _ true hzGuard]
          simp only
          rw [sqfCertifyRawExec_ok _ zRaw hmake]
          simp only
          rw [sqfCertifyRawExec_ok _ cRaw hcRun]
          simp only
          simpa [nextState] using hrecursive
        · have hzGuardFalse : (!zRaw.normalization.isEmpty &&
              get_deg zRaw.normalization > 0) = false := by
            exact Bool.eq_false_of_not_eq_true hzGuard
          have hstep := yunNoFactorStep this hcfg physical source multiplicity
            w c result y zRaw cRaw hinvariant hguard hyRun hzRun
            hzGuardFalse hcRun
          let nextState : Generated.StrictSquarefreeZp.YunRawState
              (strictSQFRawOps this hcfg physical) source :=
            ⟨multiplicity + 1, y, cRaw.normalization, result, hstep.1⟩
          rcases ih (Generated.StrictSquarefreeZp.yunStateMeasure nextState)
              (by simpa [nextState,
                Generated.StrictSquarefreeZp.yunStateMeasure] using
                  hstep.2.trans_le (Nat.le_of_eq hmeasureEq))
              nextState rfl with ⟨finalState, hrecursive⟩
          refine ⟨finalState, ?_⟩
          rw [Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state.eq_1]
          rw [sqfCertifyBool_eq _ true hguard]
          simp only [strictSQFRawOps]
          rw [sqfCertifyRawExec_ok _ y hyRun]
          simp only
          rw [sqfCertifyRawExec_ok _ zRaw hzRun]
          simp only
          rw [sqfCertifyBool_eq _ false hzGuardFalse]
          simp only
          rw [sqfCertifyRawExec_ok _ cRaw hcRun]
          simp only
          simpa [nextState] using hrecursive
      · have hguardFalse : (!w.isEmpty && get_deg w > 0) = false := by
          exact Bool.eq_false_of_not_eq_true hguard
        refine ⟨⟨⟨multiplicity, w, c, result, hinvariant⟩,
          hguardFalse⟩, ?_⟩
        rw [Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state.eq_1]
        rw [sqfCertifyBool_eq _ false hguardFalse]

/-- Successful execution of the generated Yun loop returns exactly the
deterministic L2 Yun result associated with its top-level SQF source. -/
theorem generatedYunLoopState_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp)
    (state : Generated.StrictSquarefreeZp.YunRawState
      (strictSQFRawOps this hcfg physical) source) :
    ∃ finalState : Generated.StrictSquarefreeZp.YunRawFinalState
        (strictSQFRawOps this hcfg physical) source,
      Generated.StrictSquarefreeZp._loop___squarefree_Zp_1_raw_ir_state
          (strictSQFRawOps this hcfg physical) state = .ok finalState ∧
      toPolyList finalState.state.result this._p.toNat =
        (sourceYunTarget this source).1 ∧
      SparsePolyZp.toPoly this._p.toNat finalState.state.c =
        (sourceYunTarget this source).2 := by
  rcases generatedYunLoopState_executes this hcfg physical source state with
    ⟨finalState, hrun⟩
  have hwSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat
    finalState.state.w finalState.state.valid.wMonic.ne_zero
  have hwHeadBound : finalState.state.w[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength finalState.state.w =
        finalState.state.w[0].1.deg + 1 := by
      simp [sparseDenseLength, hwSize]
    have hwBound := finalState.state.valid.wBound
    rw [hdense] at hwBound
    omega
  have hwNotPositive : ¬ 0 < (SparsePolyZp.toPoly this._p.toNat
      finalState.state.w).natDegree := by
    intro hwPositive
    have hguardTrue := (generated_yun_guard_eq_true_iff this._p.toNat
      finalState.state.w finalState.state.valid.wCanonical hwSize
      hwHeadBound).mpr hwPositive
    rw [finalState.done] at hguardTrue
    contradiction
  have hwDegree : (SparsePolyZp.toPoly this._p.toNat
      finalState.state.w).natDegree = 0 := Nat.eq_zero_of_not_pos hwNotPositive
  have hsemantic := finalState.state.valid.semanticContinuation
  rw [yunLoop, dif_pos hwDegree] at hsemantic
  refine ⟨finalState, hrun, ?_, ?_⟩
  · exact congrArg Prod.fst hsemantic
  · exact congrArg Prod.snd hsemantic

@[simp] theorem strict_squarefreeMeasure_eq_generated
    (f : SparsePolyZp) :
    Generated.StrictSquarefreeZp.squarefreeMeasure f =
      Generated.squarefreeMeasure f := by
  rfl

/-- The generated result-copy loop is the already verified concrete loop,
not a list-level specification operation. -/
theorem generated_scaleMultiplicityLoop_eq
    (index : Nat) (source output : Array (SparsePolyZp × UInt64))
    (prime : UInt64) :
    Generated.StrictSquarefreeZp.scaleMultiplicityLoop index source output
        prime =
      scaleMultiplicityLoop index source output prime := by
  rw [Generated.StrictSquarefreeZp.scaleMultiplicityLoop,
    scaleMultiplicityLoop]
  by_cases hmore : index < source.size
  · simp only [hmore, dif_pos]
    exact generated_scaleMultiplicityLoop_eq (index + 1) source
      (output.push (source[index].1, source[index].2 * prime)) prime
  · simp [hmore]
termination_by source.size - index
decreasing_by omega

theorem certifyRawExec_ok {α : Type} (run : RawExec α) (output : α)
    (hrun : run = .ok output) :
    Generated.StrictSquarefreeZp.certifyRawExec run =
      .ok ⟨output, hrun⟩ := by
  subst run
  rfl

theorem certifyBool_eq (value result : Bool) (hvalue : value = result) :
    Generated.StrictSquarefreeZp.certifyBool value = ⟨result, hvalue⟩ := by
  subst value
  rfl

/-- The nonzero-derivative L2 branch, expressed through the deterministic
target carried by the generated Yun invariant. -/
theorem sqfZp_eq_sourceYunTarget
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hderivative : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat source) ≠ 0) :
    sqfZp (SparsePolyZp.toPoly this._p.toNat source) =
      if 0 < (sourceYunTarget this source).2.natDegree then
        (sourceYunTarget this source).1 ++
          (sqfZp (Polynomial.contract this._p.toNat
            (sourceYunTarget this source).2)).map
              (fun item => (item.1, item.2 * this._p.toNat))
      else
        (sourceYunTarget this source).1 := by
  have hsourceNonzero : SparsePolyZp.toPoly this._p.toNat source ≠ 0 := by
    intro hzero
    rw [hzero] at hpositive
    simp at hpositive
  have hcNonzero : sourceYunC this source ≠ 0 := by
    unfold sourceYunC
    intro hnormalize
    have hgcdZero := normalize_eq_zero.mp hnormalize
    exact hsourceNonzero
      (zero_dvd_iff.mp (hgcdZero ▸ EuclideanDomain.gcd_dvd_left
        (SparsePolyZp.toPoly this._p.toNat source)
        (Polynomial.derivative (SparsePolyZp.toPoly this._p.toNat source))))
  rw [sqfZp, dif_neg (Nat.ne_of_gt hpositive), dif_neg hderivative]
  unfold sourceYunTarget
  rw [dif_pos hcNonzero]
  unfold sourceYunC
  simp only

set_option maxHeartbeats 800000 in
/-- Semantic composition of the generated top-level derivative-zero branch.
The recursive premise is itself the generated raw state entry. -/
theorem generatedSQFDerivativeZero_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (hentry : EntryInvariant this source)
    (hderivative : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat source) = 0)
    (hrecursive : ∀ root (hrootEntry : EntryInvariant this root),
      Generated.StrictSquarefreeZp.squarefreeMeasure root <
        Generated.StrictSquarefreeZp.squarefreeMeasure source →
      ∃ factors,
        Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
            (strictSQFRawOps this hcfg physical) ⟨root, hrootEntry⟩ =
          .ok factors ∧
        toPolyList factors this._p.toNat =
          sqfZp (SparsePolyZp.toPoly this._p.toNat root) ∧
        ∀ item ∈ factors.toList,
          SparsePolyZp.Canonical this._p.toNat item.1) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
          (strictSQFRawOps this hcfg physical) ⟨source, hentry⟩ = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
  have hdegreeWord : ∀ term ∈ source.toList,
      term.1.deg < UInt64.size := by
    intro term hterm
    obtain ⟨index, hindex, htermEq⟩ := List.mem_iff_getElem.mp hterm
    have hi : index < source.size := by simpa using hindex
    have heq : source[index] = term := by
      rw [← Array.getElem_toList hi]
      exact htermEq
    rw [← heq]
    exact lt_trans (Nat.lt_of_lt_of_le
      (sparse_degree_lt_denseLength this._p.toNat source hentry.canonical
        index hi) hentry.denseBound) (by norm_num [UInt64.size])
  have hderivativeEq : derivativeIR this source = #[] :=
    (derivativeIR_eq_empty_iff this source hcfg hentry.canonical
      hdegreeWord).mpr hderivative
  have hderivativeEmpty : (derivativeIR this source).isEmpty = true := by
    simp [hderivativeEq]
  rcases sqfDerivativeZeroIR_prepares_recursive_call this hcfg source
      hentry.canonical hentry.monic hentry.nonempty hentry.positive
      hentry.denseBound hderivative with
    ⟨root, hrootRun, hrootCanonical, hrootSemantic, hrootMonic,
      hmonicCore, hmeasure, hrootBound, hprime⟩
  have hmonicRun : makeMonicRawIR this root = .ok root := by
    simp [makeMonicRawIR, hmonicCore]
  have hstep := derivativeZeroStep this hcfg source root root hentry
    hderivativeEmpty hrootRun hmonicRun
  rcases hrecursive root hstep.1 (by simpa using hstep.2) with
    ⟨subfactors, hsubRun, hsubSemantic, hsubCanonical⟩
  let factors := Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0
    subfactors #[] source[0]!.2.prime
  have hscaledDegree : (SparsePolyZp.toPoly this._p.toNat root).natDegree *
      this._p.toNat < UInt64.size := by
    have hp := (Fact.out : Nat.Prime this._p.toNat)
    have hexpand := Polynomial.expand_contract this._p.toNat hderivative
      hp.ne_zero
    have heq : (SparsePolyZp.toPoly this._p.toNat root).natDegree *
        this._p.toNat =
        (SparsePolyZp.toPoly this._p.toNat source).natDegree := by
      rw [hrootSemantic, ← Polynomial.natDegree_expand
        (R := ZMod this._p.toNat) (p := this._p.toNat), hexpand]
    rw [heq]
    have hsourceMeasure :=
      sparseDenseLength_eq_squarefreeMeasure_eq_natDegree_succ
        this._p.toNat source hentry.canonical hentry.nonempty
    have hbound := hentry.denseBound
    rw [hsourceMeasure.1, hsourceMeasure.2] at hbound
    exact lt_trans (by omega :
      (SparsePolyZp.toPoly this._p.toNat source).natDegree < 2 ^ 63)
      (by norm_num [UInt64.size])
  refine ⟨factors, ?_, ?_, ?_⟩
  · rw [Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state.eq_1]
    simp only [strictSQFRawOps]
    rw [certifyBool_eq _ true hderivativeEmpty]
    simp only
    rw [certifyRawExec_ok _ root hrootRun]
    simp only
    rw [certifyRawExec_ok _ root hmonicRun]
    simp only
    have hsubRun' :
        Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
            (strictSQFRawOps this hcfg physical) ⟨root, hstep.1⟩ =
          .ok subfactors := by
      simpa only using hsubRun
    have hsubRunExpanded := hsubRun'
    simp only [strictSQFRawOps] at hsubRunExpanded
    rw [hsubRunExpanded]
  · change toPolyList
      (Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0 subfactors #[]
        source[0]!.2.prime) this._p.toNat = _
    rw [generated_scaleMultiplicityLoop_eq]
    have hscaled := scaleMultiplicityLoop_toPolyList_sqfZp
      this._p.toNat subfactors source[0]!.2.prime
      (SparsePolyZp.toPoly this._p.toNat root)
      (congrArg UInt64.toNat hprime) hsubSemantic hscaledDegree
    rw [hscaled]
    have hsourceSqf : sqfZp
        (SparsePolyZp.toPoly this._p.toNat source) =
        (sqfZp (Polynomial.contract this._p.toNat
          (SparsePolyZp.toPoly this._p.toNat source))).map
            (fun item => (item.1, item.2 * this._p.toNat)) := by
      rw [sqfZp, dif_neg (Nat.ne_of_gt hentry.positive),
        dif_pos hderivative]
    rw [hsourceSqf, ← hrootSemantic]
  · intro item hitem
    change item ∈ (Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0
      subfactors #[] source[0]!.2.prime).toList at hitem
    rw [generated_scaleMultiplicityLoop_eq,
      scaleMultiplicityLoop_toList] at hitem
    simp only [Array.toList_empty, List.nil_append, List.mem_map] at hitem
    rcases hitem with ⟨original, horiginal, rfl⟩
    exact hsubCanonical original horiginal

set_option maxHeartbeats 1600000 in
/-- Semantic composition of the generated nonzero-derivative branch.  The
raw prefix and generated Yun loop are executed before the optional recursive
contraction call; the accumulator is the actual generated result array. -/
theorem generatedSQFNonzeroDerivative_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (hentry : EntryInvariant this source)
    (hderivative : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat source) ≠ 0)
    (hrecursive : ∀ root (hrootEntry : EntryInvariant this root),
      Generated.StrictSquarefreeZp.squarefreeMeasure root <
        Generated.StrictSquarefreeZp.squarefreeMeasure source →
      ∃ factors,
        Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
            (strictSQFRawOps this hcfg physical) ⟨root, hrootEntry⟩ =
          .ok factors ∧
        toPolyList factors this._p.toNat =
          sqfZp (SparsePolyZp.toPoly this._p.toNat root) ∧
        ∀ item ∈ factors.toList,
          SparsePolyZp.Canonical this._p.toNat item.1) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
          (strictSQFRawOps this hcfg physical) ⟨source, hentry⟩ = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
  have hdegreeWord : ∀ term ∈ source.toList,
      term.1.deg < UInt64.size := by
    intro term hterm
    obtain ⟨index, hindex, htermEq⟩ := List.mem_iff_getElem.mp hterm
    have hi : index < source.size := by simpa using hindex
    have heq : source[index] = term := by
      rw [← Array.getElem_toList hi]
      exact htermEq
    rw [← heq]
    exact lt_trans (Nat.lt_of_lt_of_le
      (sparse_degree_lt_denseLength this._p.toNat source hentry.canonical
        index hi) hentry.denseBound) (by norm_num [UInt64.size])
  have hderivativeSemantic := derivativeIR_toPoly this source hcfg
    hentry.canonical hdegreeWord
  have hderivativeNotEmpty : (derivativeIR this source).isEmpty = false := by
    have hne : derivativeIR this source ≠ #[] := by
      intro hempty
      apply hderivative
      rw [← hderivativeSemantic]
      simp [hempty]
    simpa [Array.isEmpty_iff, hne]
  rcases sqfNonzeroDerivativeIR_prepares_yun this hcfg physical source
      hentry.canonical hentry.monic hentry.nonempty hentry.denseBound
      hderivative with
    ⟨gcdOut, c, wRaw, hgcdCore, hcOutput, hwRun, hwNorm,
      hcCanonical, hwCanonical, hcMonic, hwMonic, hcBound, hwBound,
      hbudget, hcSemantic, hwSemantic⟩
  have hgcdRun : gcdRawIR this hcfg physical source
      (derivativeIR this source) = .ok c := by
    simp [gcdRawIR, hgcdCore, hcOutput]
  have hyunInitial := yunInitial this hcfg physical source
    (derivativeIR this source) c wRaw hentry hderivativeNotEmpty rfl
    hgcdRun hwRun
  let initialState : Generated.StrictSquarefreeZp.YunRawState
      (strictSQFRawOps this hcfg physical) source :=
    ⟨1, wRaw.normalization, c, #[], hyunInitial⟩
  rcases generatedYunLoopState_refines this hcfg physical source
      initialState with
    ⟨finalState, hyunRun, hresultSemantic, hcRemSemantic⟩
  have hcRemSize := sparsePolyZp_size_pos_of_toPoly_ne_zero this._p.toNat
    finalState.state.c finalState.state.valid.cMonic.ne_zero
  have hcRemHeadBound : finalState.state.c[0].1.deg < 2 ^ 63 := by
    have hdense : sparseDenseLength finalState.state.c =
        finalState.state.c[0].1.deg + 1 := by
      simp [sparseDenseLength, hcRemSize]
    have hbound := finalState.state.valid.cBound
    rw [hdense] at hbound
    omega
  have hcGuardIff := generated_yun_guard_eq_true_iff this._p.toNat
    finalState.state.c finalState.state.valid.cCanonical hcRemSize
    hcRemHeadBound
  have hprime : source[0]!.2.prime.toNat = this._p.toNat := by
    have hsourceNonempty := hentry.nonempty
    have hmem : source[0] ∈ source.toList :=
      Array.getElem_mem_toList hsourceNonempty
    have hp := (hentry.canonical.1 source[0] hmem).1
    simpa [getElem!_pos source 0 hsourceNonempty] using hp
  by_cases hcPositive : 0 <
      (SparsePolyZp.toPoly this._p.toNat finalState.state.c).natDegree
  · have hcGuard : (!finalState.state.c.isEmpty &&
        get_deg finalState.state.c > 0) = true := hcGuardIff.mpr hcPositive
    have hcDerivative : Polynomial.derivative
        (SparsePolyZp.toPoly this._p.toNat finalState.state.c) = 0 := by
      rw [hcRemSemantic]
      exact finalState.state.valid.targetDerivativeZero
    have hcDegree : (SparsePolyZp.toPoly this._p.toNat
        finalState.state.c).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat source).natDegree := by
      rw [hcRemSemantic]
      exact finalState.state.valid.targetDegreeLeSource
    rcases sqfPostYunIR_prepares_recursive_call this hcfg source
        finalState.state.c hentry.canonical hentry.nonempty
        finalState.state.valid.cCanonical finalState.state.valid.cMonic
        hcPositive finalState.state.valid.cBound hcDerivative hcDegree with
      ⟨root, hrootRun, hrootCanonical, hrootSemantic, hrootMonic,
        hmonicCore, hmeasure, hrootBound⟩
    have hmonicRun : makeMonicRawIR this root = .ok root := by
      simp [makeMonicRawIR, hmonicCore]
    have hstep := remainderRootStep this hcfg source
      (derivativeIR this source) finalState.state.multiplicity
      finalState.state.w finalState.state.c finalState.state.result root root
      hentry hderivativeNotEmpty finalState.state.valid finalState.done
      hcGuard hrootRun hmonicRun
    rcases hrecursive root hstep.1 (by simpa using hstep.2) with
      ⟨subfactors, hsubRun, hsubSemantic, hsubCanonical⟩
    let factors := Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0
      subfactors finalState.state.result source[0]!.2.prime
    have hscaledDegree :
        (SparsePolyZp.toPoly this._p.toNat root).natDegree *
          this._p.toNat < UInt64.size := by
      have hp := (Fact.out : Nat.Prime this._p.toNat)
      have hexpand := Polynomial.expand_contract this._p.toNat hcDerivative
        hp.ne_zero
      have heq : (SparsePolyZp.toPoly this._p.toNat root).natDegree *
          this._p.toNat =
          (SparsePolyZp.toPoly this._p.toNat finalState.state.c).natDegree := by
        rw [hrootSemantic, ← Polynomial.natDegree_expand
          (R := ZMod this._p.toNat) (p := this._p.toNat), hexpand]
      rw [heq]
      have hsourceDegreeLt :
          (SparsePolyZp.toPoly this._p.toNat source).natDegree < 2 ^ 63 := by
        rcases sparseDenseLength_eq_squarefreeMeasure_eq_natDegree_succ
          this._p.toNat source hentry.canonical hentry.nonempty with
          ⟨hsourceLength, hsourceMeasure⟩
        have hbound := hentry.denseBound
        rw [hsourceLength, hsourceMeasure] at hbound
        omega
      exact lt_trans (lt_of_le_of_lt hcDegree hsourceDegreeLt)
        (by norm_num [UInt64.size])
    refine ⟨factors, ?_, ?_, ?_⟩
    · rw [Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state.eq_1]
      simp only [strictSQFRawOps]
      rw [certifyBool_eq _ false hderivativeNotEmpty]
      simp only
      rw [certifyRawExec_ok _ c hgcdRun]
      simp only
      rw [certifyRawExec_ok _ wRaw hwRun]
      simp only
      have hyunRunExpanded := hyunRun
      simp only [strictSQFRawOps, initialState] at hyunRunExpanded
      rw [hyunRunExpanded]
      simp only
      rw [certifyBool_eq _ true hcGuard]
      simp only
      rw [certifyRawExec_ok _ root hrootRun]
      simp only
      rw [certifyRawExec_ok _ root hmonicRun]
      simp only
      have hsubRunExpanded := hsubRun
      simp only [strictSQFRawOps] at hsubRunExpanded
      rw [hsubRunExpanded]
    · change toPolyList
        (Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0 subfactors
          finalState.state.result source[0]!.2.prime) this._p.toNat = _
      rw [generated_scaleMultiplicityLoop_eq]
      have hscaled := scaleMultiplicityLoop_toPolyList_sqfZp_append
        this._p.toNat subfactors finalState.state.result
        source[0]!.2.prime (SparsePolyZp.toPoly this._p.toNat root)
        hprime hsubSemantic hscaledDegree
      rw [hscaled, hresultSemantic]
      rw [sqfZp_eq_sourceYunTarget this source hentry.positive hderivative,
        if_pos (by simpa [hcRemSemantic] using hcPositive), hrootSemantic,
        hcRemSemantic]
    · intro item hitem
      change item ∈ (Generated.StrictSquarefreeZp.scaleMultiplicityLoop 0
        subfactors finalState.state.result source[0]!.2.prime).toList at hitem
      rw [generated_scaleMultiplicityLoop_eq,
        scaleMultiplicityLoop_toList] at hitem
      simp only [List.drop_zero, List.mem_append, List.mem_map] at hitem
      rcases hitem with hprevious | ⟨original, horiginal, rfl⟩
      · exact finalState.state.valid.resultCanonical item hprevious
      · exact hsubCanonical original horiginal
  · have hcGuardFalse : (!finalState.state.c.isEmpty &&
        get_deg finalState.state.c > 0) = false :=
      Bool.eq_false_of_not_eq_true (fun htrue =>
        hcPositive (hcGuardIff.mp htrue))
    refine ⟨finalState.state.result, ?_, ?_,
      finalState.state.valid.resultCanonical⟩
    · rw [Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state.eq_1]
      simp only [strictSQFRawOps]
      rw [certifyBool_eq _ false hderivativeNotEmpty]
      simp only
      rw [certifyRawExec_ok _ c hgcdRun]
      simp only
      rw [certifyRawExec_ok _ wRaw hwRun]
      simp only
      have hyunRunExpanded := hyunRun
      simp only [strictSQFRawOps, initialState] at hyunRunExpanded
      rw [hyunRunExpanded]
      simp only
      rw [certifyBool_eq _ false hcGuardFalse]
    · rw [hresultSemantic,
        sqfZp_eq_sourceYunTarget this source hentry.positive hderivative,
        if_neg (by simpa [hcRemSemantic] using hcPositive)]

set_option maxHeartbeats 2000000 in
/-- The complete generated outer state recursion refines `sqfZp`. -/
theorem generatedSQFState_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg) :
    ∀ state : Generated.StrictSquarefreeZp.SQFRawState
        (strictSQFRawOps this hcfg physical),
      ∃ factors,
        Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
            (strictSQFRawOps this hcfg physical) state = .ok factors ∧
        toPolyList factors this._p.toNat =
          sqfZp (SparsePolyZp.toPoly this._p.toNat state.f) ∧
        ∀ item ∈ factors.toList,
          SparsePolyZp.Canonical this._p.toNat item.1 := by
  suffices hstrong : ∀ n state,
      Generated.StrictSquarefreeZp.sqfStateMeasure state = n →
      ∃ factors,
        Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir_state
            (strictSQFRawOps this hcfg physical) state = .ok factors ∧
        toPolyList factors this._p.toNat =
          sqfZp (SparsePolyZp.toPoly this._p.toNat state.f) ∧
        ∀ item ∈ factors.toList,
          SparsePolyZp.Canonical this._p.toNat item.1 by
    intro state
    exact hstrong (Generated.StrictSquarefreeZp.sqfStateMeasure state)
      state rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
      intro state hmeasureEq
      rcases state with ⟨source, hentry⟩
      by_cases hderivative : Polynomial.derivative
          (SparsePolyZp.toPoly this._p.toNat source) = 0
      · apply generatedSQFDerivativeZero_refines this hcfg physical source
          hentry hderivative
        intro root hrootEntry hdecrease
        apply ih (Generated.StrictSquarefreeZp.sqfStateMeasure
          (⟨root, hrootEntry⟩ : Generated.StrictSquarefreeZp.SQFRawState
            (strictSQFRawOps this hcfg physical)))
        · simpa [Generated.StrictSquarefreeZp.sqfStateMeasure] using
            hdecrease.trans_le (Nat.le_of_eq hmeasureEq)
        · rfl
      · apply generatedSQFNonzeroDerivative_refines this hcfg physical source
          hentry hderivative
        intro root hrootEntry hdecrease
        apply ih (Generated.StrictSquarefreeZp.sqfStateMeasure
          (⟨root, hrootEntry⟩ : Generated.StrictSquarefreeZp.SQFRawState
            (strictSQFRawOps this hcfg physical)))
        · simpa [Generated.StrictSquarefreeZp.sqfStateMeasure] using
            hdecrease.trans_le (Nat.le_of_eq hmeasureEq)
        · rfl

/-- Public genuine C++ L1 → Lean L2 refinement contract.  Its name retains
the original C++ `__squarefree_Zp` spelling, including both underscores. -/
theorem __squarefree_Zp_raw_ir_refines_sqfZp
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat source).Monic)
    (hnonempty : 0 < source.size)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hbound : sparseDenseLength source ≤ 2 ^ 63) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
          (strictSQFRawOps this hcfg physical) source
          (fun _ => ⟨hcanonical, hmonic, hnonempty, hpositive, hbound⟩) =
        .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
  let hentry : EntryInvariant this source :=
    ⟨hcanonical, hmonic, hnonempty, hpositive, hbound⟩
  rcases generatedSQFState_refines this hcfg physical
      (⟨source, hentry⟩ : Generated.StrictSquarefreeZp.SQFRawState
        (strictSQFRawOps this hcfg physical)) with
    ⟨factors, hrun, hsemantic, hcanonicalResult⟩
  have hnotEmpty : source.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hnonempty]
  refine ⟨factors, ?_, hsemantic, hcanonicalResult⟩
  simp only [Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir,
    hnotEmpty]
  simpa [hentry] using hrun

end Refinement.StrictSquarefreeGenerated
