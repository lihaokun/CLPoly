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
    (multiplicity : UInt64) (w c : SparsePolyZp)
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
