import CLPoly.Impl.StrictHGCDContinuationRefinement
import CLPoly.Generated.StrictHGCDChecked

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic
/-- Complete proof-erased provider value consumed by
exactly, and every computational field of a successful result is retained. -/
def hgcdRecursiveExecValue (run : RawExec HgcdRecursiveResult) :
    RawExec HgcdRecursiveValue :=
  run.map HgcdRecursiveResult.value

/-- Equality of proof-erased recursive executions is enough for equality of
the original executions because `HgcdRecursiveResult.valid` is proof-only. -/
theorem hgcdRecursiveExec_ext_value
    (left right : RawExec HgcdRecursiveResult)
    (hvalue : hgcdRecursiveExecValue left = hgcdRecursiveExecValue right) :
    left = right := by
  cases left with
  | error leftFault =>
      cases right with
      | error rightFault =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          exact congrArg Except.error (Except.error.inj hvalue)
      | ok rightResult =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          contradiction
  | ok leftResult =>
      cases right with
      | error rightFault =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          contradiction
      | ok rightResult =>
          simp only [hgcdRecursiveExecValue, Except.map, Except.ok.injEq]
            at hvalue
          exact congrArg Except.ok
            (HgcdRecursiveResult.ext_value leftResult rightResult hvalue)

/-- Proof erasure of `hgcdRecursiveCallBelowOfCall` is definitionally the
ordinary callback execution at every actual smaller call. -/
theorem hgcdRecursiveCallBelowOfCall_apply (bound : Nat)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (horder : lenB < lenA) (hdecrease : lenA < bound) :
    hgcdRecursiveCallBelowOfCall bound plain M hM computeM A B a b lenA lenB
        W scratch heap horder hdecrease =
      plain M hM computeM A B a b lenA lenB W scratch heap := rfl

/-- The proof-indexed well-founded body is the exact source-shaped body after
erasing only its strict-decrease proofs. -/
theorem hgcdRecursiveBodyBelow_eq_body (this : DenseUPolyZp)
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (firstLength : ¬ lenB < lenA / 2 + 1 → ∀ first,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound below ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first →
      HgcdRecursiveLengthInvariant high.lenA0 first)
    (reconstructionBound : ¬ lenB < lenA / 2 + 1 →
      HgcdFirstReconstructionBoundProvider this a b W
      scratch lenA lenB (HgcdFirstDispatchResult this bound below a b W
        scratch lenA lenB heap))
    (hagrees : ∀ (matrix : HgcdMat) (hMatrix : matrix.Valid)
      (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
      (WNext childScratch : RawPtr UInt64) (childHeap : RawHeap)
      (hchildOrder : lenInputB < lenInputA)
      (hchildDecrease : lenInputA < bound),
      below matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap hchildOrder hchildDecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap) :
    hgcdRecursiveBodyBelow this bound below M hM computeM A B a b lenA lenB
        W scratch heap hbound horder firstLength reconstructionBound =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  apply hgcdRecursiveExec_ext_value
  rw [hgcdRecursiveBodyBelow, hgcdRecursiveBody]
  split
  · rfl
  · simp only [hgcdRecursiveExecValue, hgcdRecursiveDispatchBelow,
      hgcdRecursiveDispatch, hagrees]
    split <;> simp_all
    split <;> simp_all
    split <;> simp_all
    split <;> simp_all



/-- Branch-local admissibility changes no generated computation.  Erasing the
strict-decrease proofs yields exactly the source-shaped recursive body on both
the base and non-base paths. -/
theorem hgcdRecursiveBodyBranchAdmissible_eq_body (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (nonBase : ¬ lenB < lenA / 2 + 1 →
      HgcdRecursiveNonBasePackage this bound below a b W scratch lenA lenB
        heap)
    (recursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound below a b W scratch lenA
        lenB heap (nonBase hbase))
    (hagrees : ∀ (matrix : HgcdMat) (hMatrix : matrix.Valid)
      (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
      (WNext childScratch : RawPtr UInt64) (childHeap : RawHeap)
      (hchildOrder : lenInputB < lenInputA)
      (hchildDecrease : lenInputA < bound),
      below matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap hchildOrder hchildDecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap) :
    hgcdRecursiveBodyBranchAdmissible this bound below M hM computeM A B a b
        lenA lenB W scratch heap hcfg hp hbound horder nonBase
        recursiveRefines =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  unfold hgcdRecursiveBodyBranchAdmissible
  split
  next hbase =>
    exact hgcdRecursiveBodyBelow_eq_body this bound below plain M hM computeM
      A B a b lenA lenB W scratch heap hbound horder
      (fun hnonbase => False.elim (hnonbase hbase))
      (fun hnonbase => False.elim (hnonbase hbase)) hagrees
  next hbase =>
    let package := nonBase hbase
    let providers := hgcdRecursiveFirstCall_providers this bound below a b W
      scratch lenA lenB heap package.inputHighA package.inputHighB
      package.lowPolyA package.lowPolyB hcfg hp horder
      package.workspace (recursiveRefines hbase)
    exact hgcdRecursiveBodyBelow_eq_body this bound below plain M hM computeM
      A B a b lenA lenB W scratch heap hbound horder
      (fun _ => providers.1) (fun _ => providers.2) hagrees

/-- Specialization of branch-local admissibility to the proof-erased view of
one ordinary recursive call.  The resulting execution is exactly the
generated source body using that same call as its recursive callback. -/
theorem hgcdRecursiveBodyBranchAdmissible_belowOfCall_eq_body
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (nonBase : ¬ lenB < lenA / 2 + 1 →
      HgcdRecursiveNonBasePackage this bound
        (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB
        heap)
    (recursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound
        (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB
        heap (nonBase hbase)) :
    hgcdRecursiveBodyBranchAdmissible this bound
        (hgcdRecursiveCallBelowOfCall bound plain) M hM computeM A B a b lenA
        lenB W scratch heap hcfg hp hbound horder nonBase recursiveRefines =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  apply hgcdRecursiveBodyBranchAdmissible_eq_body this bound
    (hgcdRecursiveCallBelowOfCall bound plain) plain M hM computeM A B a b
    lenA lenB W scratch heap hcfg hp hbound horder nonBase recursiveRefines
  intro matrix hMatrix a3 b3 inputA inputB lenInputA lenInputB WNext
    childScratch childHeap hchildOrder hchildDecrease
  rfl

/-- On an ordered source invocation, the checked well-founded call unfolds
to the exact branch-admissible generated body.  Both dynamic checks disappear
using the strict order/decrease proofs already carried by the generated
`Below` dispatch; no fixed-point hypothesis is assumed. -/
theorem hgcdRecursiveCallChecked_eq_branchAdmissible
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (nonBase : ¬ lenB < lenA / 2 + 1 →
      HgcdRecursiveNonBasePackage this bound
        (hgcdRecursiveCallBelowOfCall bound
          (hgcdRecursiveCallChecked this)) a b W scratch lenA lenB heap)
    (recursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound
        (hgcdRecursiveCallBelowOfCall bound
          (hgcdRecursiveCallChecked this)) a b W scratch lenA lenB heap
        (nonBase hbase)) :
    hgcdRecursiveCallChecked this M hM computeM A B a b lenA lenB W scratch
        heap =
      hgcdRecursiveBodyBranchAdmissible this bound
        (hgcdRecursiveCallBelowOfCall bound (hgcdRecursiveCallChecked this))
        M hM computeM A B a b lenA lenB W scratch heap hcfg hp hbound horder
        nonBase recursiveRefines := by
  let guarded : HgcdRecursiveCall :=
    fun childM childHM childCompute childA childB childa childb childLenA
        childLenB childW childScratch childHeap =>
      if hchildOrder : childLenB < childLenA then
        if hchildDecrease : childLenA < lenA then
          hgcdRecursiveCallChecked this childM childHM childCompute childA
            childB childa childb childLenA childLenB childW childScratch
            childHeap
        else
          .error .assertionFailure
      else
        .error .assertionFailure
  have hunfoldChecked :
      hgcdRecursiveCallChecked this M hM computeM A B a b lenA lenB W scratch
          heap =
        hgcdRecursiveBody this guarded M hM computeM A B a b lenA lenB W
          scratch heap := by
    rw [hgcdRecursiveCallChecked]
  have hbody := hgcdRecursiveBodyBranchAdmissible_eq_body this bound
    (hgcdRecursiveCallBelowOfCall bound (hgcdRecursiveCallChecked this))
    guarded M hM computeM A B a b lenA lenB W scratch heap hcfg hp hbound
    horder nonBase recursiveRefines
    (by
      intro childM childHM childa childb inputA inputB childLenA childLenB
        childW childScratch childHeap hchildOrder hchildDecrease
      have hdecrease : childLenA < lenA := by omega
      simp [hgcdRecursiveCallBelowOfCall, guarded, hchildOrder, hdecrease])
  exact hunfoldChecked.trans hbody.symm

/-- Genuine refinement theorem for the concrete checked well-founded HGCD
execution.  It has no abstract callback and no assumed fixed-point unfolding
equation. -/
theorem hgcdRecursiveCallChecked_rawInvariant_wf
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (horder : lenB < lenA)
    (workspace : HgcdRecursiveInvocationWorkspace this
      (hgcdRecursiveCallChecked this) lenA M hM computeM A B a b lenA lenB W
      scratch heap left right hcfg hp rfl horder)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this
      (hgcdRecursiveCallChecked this) hcfg hp lenA)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM computeM A B a b lenA lenB W
      scratch heap = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries computeM
        A B lenA result := by
  let recurse := hgcdRecursiveCallBelowOfCall lenA
    (hgcdRecursiveCallChecked this)
  let firstRecursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this lenA recurse a b W scratch lenA
        lenB heap (workspace.nonBase hbase) := by
    intro hbase
    let package := workspace.nonBase hbase
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    dsimp [HgcdRecursiveNonBaseCallbackRefines]
    intro hchildOrder hchildDecrease
    dsimp [HgcdRecursiveCallbackRefinesAt]
    intro childResult hchildRun
    let childWorkspace := provider high.lenA0 hchildDecrease ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) true ws.a3 ws.b3 high.a0
      high.b0 high.lenB0 ws.next scratch heap package.inputHighA
      package.inputHighB hchildOrder package.workspace.iter.leftRep
      package.workspace.iter.rightRep
    have hchildChecked : hgcdRecursiveCallChecked this ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) true ws.a3 ws.b3 high.a0
        high.b0 high.lenA0 high.lenB0 ws.next scratch heap = .ok childResult := by
      simpa [recurse, hgcdRecursiveCallBelowOfCall] using hchildRun
    exact hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) true ws.a3 ws.b3 high.a0
      high.b0 high.lenA0 high.lenB0 ws.next scratch heap package.inputHighA
      package.inputHighB hchildOrder childWorkspace
      (provider.mono this (hgcdRecursiveCallChecked this) hcfg hp lenA
        high.lenA0 hchildDecrease) childResult hchildChecked
  let secondRecursiveRefines : ∀ (hbase : ¬ lenB < lenA / 2 + 1)
      (bodyResult : HgcdRecursiveResult)
      (hrunBody :
        let package := workspace.nonBase hbase
        let firstCall := package.admissible this lenA recurse a b W scratch
          lenA lenB heap (firstRecursiveRefines hbase)
        let providers := hgcdRecursiveFirstCall_providers this lenA recurse a b
          W scratch lenA lenB heap package.inputHighA package.inputHighB
          package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
          firstCall.recursiveRefines
        hgcdRecursiveBodyBelow this lenA recurse M hM computeM A B a b lenA
          lenB W scratch heap rfl horder (fun _ => providers.1)
            (fun _ => providers.2) = .ok bodyResult),
      HgcdRecursiveContinuationSecondRefines this lenA recurse M hM computeM
        A B a b W scratch lenA lenB heap (workspace.nonBase hbase) rfl horder
        hbase (workspace.continuation hbase (firstRecursiveRefines hbase)
          bodyResult hrunBody) := by
    intro hbase bodyResult hrunBody
    let next := workspace.continuation hbase (firstRecursiveRefines hbase)
      bodyResult hrunBody
    change HgcdRecursiveContinuationSecondRefines this lenA recurse M hM
      computeM A B a b W scratch lenA lenB heap (workspace.nonBase hbase) rfl
      horder hbase next
    cases next with
    | early earlyWorkspace =>
        trivial
    | nonEarly nonEarlyWorkspace =>
        dsimp [HgcdRecursiveContinuationSecondRefines,
          HgcdRecursiveSecondCallbackRefines, HgcdRecursiveCallbackRefinesAt]
        intro highC highD hHighC hHighD childResult hchildRun
        let ws := hgcdRecursiveWorkspace W lenA
        let childWorkspace := provider nonEarlyWorkspace.middle.lenC0
          nonEarlyWorkspace.secondDecrease ws.S
          (hgcdRecursiveWorkspace_S_valid W lenA) true ws.a3 ws.b3
          nonEarlyWorkspace.middle.c0 nonEarlyWorkspace.middle.d0
          nonEarlyWorkspace.middle.lenD0 ws.next scratch
          nonEarlyWorkspace.middle.heap highC highD
          nonEarlyWorkspace.secondOrder hHighC hHighD
        have hchildChecked : hgcdRecursiveCallChecked this ws.S
            (hgcdRecursiveWorkspace_S_valid W lenA) true ws.a3 ws.b3
            nonEarlyWorkspace.middle.c0 nonEarlyWorkspace.middle.d0
            nonEarlyWorkspace.middle.lenC0 nonEarlyWorkspace.middle.lenD0
            ws.next scratch nonEarlyWorkspace.middle.heap = .ok childResult := by
          simpa [recurse, hgcdRecursiveCallBelowOfCall] using hchildRun
        exact hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp ws.S
          (hgcdRecursiveWorkspace_S_valid W lenA) true ws.a3 ws.b3
          nonEarlyWorkspace.middle.c0 nonEarlyWorkspace.middle.d0
          nonEarlyWorkspace.middle.lenC0 nonEarlyWorkspace.middle.lenD0
          ws.next scratch nonEarlyWorkspace.middle.heap highC highD
          nonEarlyWorkspace.secondOrder childWorkspace
          (provider.mono this (hgcdRecursiveCallChecked this) hcfg hp lenA
            nonEarlyWorkspace.middle.lenC0
            nonEarlyWorkspace.secondDecrease) childResult hchildChecked
  have hbody : hgcdRecursiveBodyBranchAdmissible this lenA recurse M hM
      computeM A B a b lenA lenB W scratch heap hcfg hp rfl horder
      workspace.nonBase firstRecursiveRefines = .ok result := by
    rw [← hgcdRecursiveCallChecked_eq_branchAdmissible this lenA M hM
      computeM A B a b lenA lenB W scratch heap hcfg hp rfl horder
      workspace.nonBase firstRecursiveRefines]
    exact hrun
  exact hgcdRecursiveBodyBranchAdmissible_rawInvariant this lenA recurse M hM
    computeM A B a b lenA lenB W scratch heap left right hcfg hp rfl horder
    workspace.positive workspace.base workspace.nonBase workspace.nonBaseInput
    firstRecursiveRefines
    (fun hbase bodyResult hrunBody =>
      workspace.continuation hbase (firstRecursiveRefines hbase) bodyResult
        hrunBody)
    secondRecursiveRefines result hbody
termination_by lenA
decreasing_by
  · exact hchildDecrease
  · exact nonEarlyWorkspace.secondDecrease


end CLPoly.Impl.StrictHGCDRawRefinement
