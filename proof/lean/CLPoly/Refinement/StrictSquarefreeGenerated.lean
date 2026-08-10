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
  semanticContinuation : ∃
      (targetResult : List (Polynomial (ZMod this._p.toNat) × Nat))
      (targetRemainder : Polynomial (ZMod this._p.toNat)),
    yunLoop (SparsePolyZp.toPoly this._p.toNat w)
        (SparsePolyZp.toPoly this._p.toNat c) multiplicity.toNat
        (toPolyList result this._p.toNat) cMonic.ne_zero =
      (targetResult, targetRemainder) ∧
    Polynomial.derivative targetRemainder = 0 ∧
    targetRemainder.natDegree ≤
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
        index hi) hinvariant.denseBound) (by native_decide)
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
        index hi) hinvariant.denseBound) (by native_decide)
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
  refine ⟨hwCanonical, hcCanonical, hwMonic, hcMonic, hwBound, hcBound,
    ?_, ?_, ⟨target.1, target.2, ?_, htargetDerivative,
      htargetDegree⟩⟩
  · simpa using hbudget
  · intro item hitem
    simp at hitem
  · simp [target]

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
  rcases hinvariant.semanticContinuation with
    ⟨targetResult, targetRemainder, htarget, hderivative, hdegree⟩
  have hcontinuation : yunLoop (SparsePolyZp.toPoly this._p.toNat y)
      (SparsePolyZp.toPoly this._p.toNat cRaw)
      (multiplicity + 1).toNat
      (toPolyList (result.push (zRaw, multiplicity)) this._p.toNat)
      hcMonic.ne_zero =
        (targetResult, targetRemainder) :=
    hloopStep.symm.trans htarget
  rw [hcNorm]
  refine ⟨⟨hyCanonical, hcCanonical, hyMonic, hcMonic, hyBound, hcBound,
    ?_, ?_, ⟨targetResult, targetRemainder, ?_, hderivative,
      hdegree⟩⟩, ?_⟩
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
  rcases hinvariant.semanticContinuation with
    ⟨targetResult, targetRemainder, htarget, hderivative, hdegree⟩
  have hcontinuation : yunLoop (SparsePolyZp.toPoly this._p.toNat y)
      (SparsePolyZp.toPoly this._p.toNat cRaw)
      (multiplicity + 1).toNat (toPolyList result this._p.toNat)
      hcMonic.ne_zero =
        (targetResult, targetRemainder) :=
    hloopStep.symm.trans htarget
  rw [hcNorm]
  exact ⟨⟨hyCanonical, hcCanonical, hyMonic, hcMonic, hyBound, hcBound,
    by simpa [Nat.add_assoc] using hbudget, hinvariant.resultCanonical,
    ⟨targetResult, targetRemainder,
      by simpa [hcNorm] using hcontinuation, hderivative, hdegree⟩⟩,
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
  rcases hinvariant.semanticContinuation with
    ⟨targetResult, targetRemainder, hcontinuation, htargetDerivative,
      htargetDegree⟩
  rw [yunLoop, dif_pos hwDegree] at hcontinuation
  have hcSemantic : SparsePolyZp.toPoly this._p.toNat c =
      targetRemainder := congrArg Prod.snd hcontinuation
  have hcDerivative : Polynomial.derivative
      (SparsePolyZp.toPoly this._p.toNat c) = 0 := by
    rw [hcSemantic]
    exact htargetDerivative
  have hcDegree : (SparsePolyZp.toPoly this._p.toNat c).natDegree ≤
      (SparsePolyZp.toPoly this._p.toNat f).natDegree := by
    rw [hcSemantic]
    exact htargetDegree
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

end Refinement.StrictSquarefreeGenerated
