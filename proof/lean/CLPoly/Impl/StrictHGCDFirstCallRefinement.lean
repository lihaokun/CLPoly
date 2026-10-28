import CLPoly.Impl.StrictHGCDRawRefinement

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
`hgcdRecursiveBodyBelow`.  Its fields come from one successful raw paired
reconstruction plus the real first-child length invariant. -/
theorem hgcdRecursiveFirstReconstruct_invariant_of_execution
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 a b highA highB scratch : RawPtr UInt64)
    (lenA lenB : Nat) (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinputOrder : lenB < lenA)
    (hinvariant : HgcdRecursiveLengthInvariant (lenA - lenA / 2) first)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      a b highA highB scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap)
    (hMatrix : HgcdMatRawDenseRep this first.heap first.matrix entries
      first.valid)
    (hLowA : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) polyLowA)
    (hLowB : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) polyLowB)
    (hHighA : RawDensePolyRep this first.heap highA first.lenA polyHighA)
    (hHighB : RawDensePolyRep this first.heap highB first.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 a b highA highB scratch
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) first.lenA
      first.lenB (lenA / 2) first.matrix first.valid first.sgn first.heap =
        .ok result) :
    HgcdFirstReconstructionInvariant lenA first result := by
  have horder := hgcdRecursiveFirstReconstruct_order_of_invariant this A B
    T0 a b highA highB scratch lenA lenB first result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp hinvariant physical hMatrix hLowA
    hLowB hHighA hHighB hrun
  exact {
    leadingA := horder.1
    positiveA := horder.2.1
    ordered := horder.2.2
    decreases := hgcdRecursiveFirstReconstruct_bound_of_invariant this A B
      T0 a b highA highB scratch lenA lenB first result entries polyLowA
      polyLowB polyHighA polyHighB hcfg hp hinputOrder hinvariant physical
      hMatrix hLowA hLowB hHighA hHighB hrun }

/-- Physical workspace for the reconstruction that follows an actual first
dispatch result. -/
def HgcdRecursiveFirstReconstructWorkspaceProvider (this : DenseUPolyZp)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat)
    (actualFirst : HgcdRecursiveResult → Prop) : Prop :=
  let ws := hgcdRecursiveWorkspace W lenA
  ∀ first, actualFirst first →
    HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0 a b
      ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap

/-- Construct the first-reconstruction bound from the actual child execution,
its induction invariant, and physical frame/workspace facts. -/
theorem hgcdFirstReconstructionBoundProvider_of_actual_dispatch
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (actualFirst : HgcdRecursiveResult → Prop)
    (lowPolyA lowPolyB : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (hLowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA)
    (hLowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB)
    (frame : ∀ first, actualFirst first →
      HgcdRecursiveFirstDispatchFrame a b (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)) heap first)
    (reconstructWorkspace :
      HgcdRecursiveFirstReconstructWorkspaceProvider this a b W scratch lenA
        lenB actualFirst)
    (firstInvariant : ∀ first, actualFirst first →
      ∃ inputHighA inputHighB outputHighA outputHighB entries,
        let ws := hgcdRecursiveWorkspace W lenA
        HgcdRecursiveRawInvariant this inputHighA inputHighB outputHighA
          outputHighB entries true ws.a3 ws.b3 (lenA - lenA / 2) first) :
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      actualFirst := by
  intro first reconstructed hactual hlength hreconstruct
  let ws := hgcdRecursiveWorkspace W lenA
  rcases firstInvariant first hactual with
    ⟨inputHighA, inputHighB, outputHighA, outputHighB, entries, hFirst⟩
  have hframe := frame first hactual
  have hLowAAfter : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA :=
    rawDensePolyRep_of_same_prefix this heap first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA hframe.layout hframe.lowAPrefix hLowA
  have hLowBAfter : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB :=
    rawDensePolyRep_of_same_prefix this heap first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB hframe.layout hframe.lowBPrefix hLowB
  exact hgcdRecursiveFirstReconstruct_invariant_of_execution this ws.a2 ws.b2
    ws.T0 a b ws.a3 ws.b3 scratch lenA lenB first reconstructed entries
    lowPolyA lowPolyB outputHighA outputHighB hcfg hp horder hlength
    (reconstructWorkspace first hactual) (hFirst.matrixSemantics rfl).1
    hLowAAfter hLowBAfter hFirst.aRep hFirst.bRep hreconstruct

/-- Assemble the body reconstruction provider directly from the real first
cutoff dispatch, its recursive induction hypothesis, and physical frames. -/
theorem hgcdFirstReconstructionBoundProvider_of_dispatch
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (iterWorkspace :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      HgcdRecursiveDispatchIterWorkspace this ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 heap
        inputHighA inputHighB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB)
    (frame :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveFirstDispatchFrameProvider this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease a b (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)))
    (reconstructWorkspace : HgcdRecursiveFirstReconstructWorkspaceProvider
      this a b W scratch lenA lenB
        (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
          heap))
    (hLowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA)
    (hLowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB) :
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
        heap) := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  apply hgcdFirstReconstructionBoundProvider_of_actual_dispatch this a b W
    scratch lenA lenB heap
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB heap)
      lowPolyA lowPolyB hcfg hp horder hLowA hLowB
  · intro first hactual
    apply hgcdRecursiveFirstDispatchResult_frame this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap a b (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) frame
      first
    simpa [HgcdFirstDispatchResult, ws, high] using hactual
  · exact reconstructWorkspace
  · intro first hactual
    rcases hgcdRecursiveDispatchResult_rawInvariant this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap inputHighA inputHighB hcfg hp iterWorkspace recursiveRefines first
        (by simpa [HgcdFirstDispatchResult, ws, high] using hactual) with
      ⟨outputHighA, outputHighB, entries, hInvariant⟩
    exact ⟨inputHighA, inputHighB, outputHighA, outputHighB, entries, by
      simpa [high, hgcdRecursiveHighInput] using hInvariant⟩

/-- Admissible physical/representation workspace for the first child call of
one recursive body.  Recursive semantics are deliberately not stored here;
they are supplied only by the well-founded induction hypothesis. -/
structure HgcdRecursiveFirstCallWorkspace (this : DenseUPolyZp)
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat)) : Prop where
  iter :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    HgcdRecursiveDispatchIterWorkspace this ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 heap
      inputHighA inputHighB
  frame :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    ∀ (hchildOrder : high.lenB0 < high.lenA0)
      (hchildDecrease : high.lenA0 < bound),
    HgcdRecursiveFirstDispatchFrameProvider this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap hchildOrder hchildDecrease a b (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2))
  reconstruct : HgcdRecursiveFirstReconstructWorkspaceProvider this a b W
    scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB heap)
  lowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA
  lowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB

/-- First-call data available at one well-founded body invocation.  The
semantic field is precisely the smaller-call induction hypothesis. -/
structure HgcdRecursiveFirstCallAdmissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat)) : Prop where
  workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB
  recursiveRefines :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    ∀ (hchildOrder : high.lenB0 < high.lenA0)
      (hchildDecrease : high.lenA0 < bound),
    HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
      inputHighA inputHighB

/-- Attach only the well-founded child semantic theorem to an already
established physical first-call workspace. -/
theorem HgcdRecursiveFirstCallWorkspace.admissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB) :
    HgcdRecursiveFirstCallAdmissible this bound recurse a b W scratch lenA
      lenB heap inputHighA inputHighB lowPolyA lowPolyB :=
  ⟨workspace, recursiveRefines⟩

/-- A first-call admissible workspace plus the child induction hypothesis
supplies exactly both proof arguments consumed by `hgcdRecursiveBodyBelow`. -/
theorem hgcdRecursiveFirstCall_providers (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB) :
    (∀ first,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first →
      HgcdRecursiveLengthInvariant high.lenA0 first) ∧
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
        heap) := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  constructor
  · intro first
    dsimp only
    intro hchildOrder hchildDecrease hrun
    exact hgcdRecursiveDispatchBelow_lengthInvariant this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap inputHighA inputHighB hcfg hp workspace.iter recursiveRefines first
      hchildOrder hchildDecrease hrun
  · exact hgcdFirstReconstructionBoundProvider_of_dispatch this bound
      recurse a b W scratch lenA lenB heap inputHighA inputHighB lowPolyA
      lowPolyB hcfg hp horder workspace.iter recursiveRefines workspace.frame
      workspace.reconstruct workspace.lowA workspace.lowB

/-- The actual first dispatch simultaneously supplies its recursive semantic
invariant and transports both source low prefixes into the returned heap. -/
theorem hgcdRecursiveFirstCall_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (admissible : HgcdRecursiveFirstCallAdmissible this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (first : HgcdRecursiveResult)
    (hchildOrder :
      (hgcdRecursiveHighInput a b lenA lenB).lenB0 <
        (hgcdRecursiveHighInput a b lenA lenB).lenA0)
    (hchildDecrease :
      (hgcdRecursiveHighInput a b lenA lenB).lenA0 < bound)
    (hrun :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first) :
    ∃ outputHighA outputHighB entries,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      HgcdRecursiveRawInvariant this inputHighA inputHighB outputHighA
          outputHighB entries true ws.a3 ws.b3 high.lenA0 first ∧
        RawDensePolyRep this first.heap a (Nat.min lenA (lenA / 2)) lowPolyA ∧
        RawDensePolyRep this first.heap b (Nat.min lenB (lenA / 2)) lowPolyB := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  rcases hgcdRecursiveDispatchBelow_rawInvariant this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap inputHighA inputHighB hcfg hp hchildOrder hchildDecrease
      admissible.workspace.iter
      (admissible.recursiveRefines hchildOrder hchildDecrease) first hrun with
    ⟨outputHighA, outputHighB, entries, hInvariant⟩
  refine ⟨outputHighA, outputHighB, entries, ?_⟩
  dsimp only
  refine ⟨hInvariant, ?_⟩
  have hframe : HgcdRecursiveFirstDispatchFrame a b
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) heap first := by
    exact admissible.workspace.frame hchildOrder hchildDecrease first hrun
  have hLowAAfter : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA :=
    rawDensePolyRep_of_same_prefix this heap first.heap a
    (Nat.min lenA (lenA / 2)) lowPolyA hframe.layout hframe.lowAPrefix
    admissible.workspace.lowA
  have hLowBAfter : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB :=
    rawDensePolyRep_of_same_prefix this heap first.heap b
    (Nat.min lenB (lenA / 2)) lowPolyB hframe.layout hframe.lowBPrefix
    admissible.workspace.lowB
  exact ⟨hLowAAfter, hLowBAfter⟩


end CLPoly.Impl.StrictHGCDRawRefinement
