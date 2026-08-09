import CLPoly.Generated.StrictGCDHGCD
import CLPoly.Impl.StrictHGCDCheckedRefinement
import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictGCDHGCDRefinement

open Generated.StrictHGCD
open Generated.StrictGCDHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictHGCDRawRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- Physical divrem storage available at every represented state reached by
the source HGCD-GCD loop.  The provider supplies only allocation and aliasing
facts; quotient and remainder semantics still come from the actual raw call. -/
def GcdHgcdLoopDivremWorkspaceProvider (this : DenseUPolyZp)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3) (capacity : Nat) : Prop :=
  ∀ (heap : RawHeap) (lenG lenJ : Nat)
    (left right : Polynomial (ZMod this._p.toNat)),
    lenG ≤ capacity → lenJ ≤ capacity →
    RawDensePolyRep this heap G lenG left →
    RawDensePolyRep this heap J lenJ right →
    EuclidWorkspace heap Q W3 G J R capacity

/-- Physical local-vector storage for the exact `_gcd_euclid` call selected
by the small branch.  It is conditional on the represented operands and does
not state the call's result. -/
def GcdHgcdLoopEuclidWorkspaceProvider (this : DenseUPolyZp)
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (capacity : Nat) : Prop :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat)),
    lenJ ≤ capacity → lenR ≤ capacity →
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    GcdEuclidRawWorkspace heap G J R euclidA euclidB euclidQ euclidR
      euclidW3 capacity

/-- Exact-size physical contract for one dynamically reached loop-head
division.  Unlike the fixed-capacity convenience provider above, this record
remains usable after an HGCD call returns new descriptor lengths. -/
structure GcdHgcdDivremWorkspace (heap : RawHeap)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenG lenJ : Nat) : Prop where
  validQ : heap.ValidU64Slice Q (lenG - (lenJ - 1))
  validR : heap.ValidU64Slice R (Nat.min lenG (lenJ - 1))
  validW3 : heap.ValidWord3Slice W3 lenG
  validGForDivisor : heap.ValidU64Slice G lenJ
  quotientCapacity : lenG - (lenJ - 1) < limbBase
  r_g : R.region ≠ G.region
  w3_g : W3.region ≠ G.region
  w3_j : W3.region ≠ J.region
  q_j : Q.region ≠ J.region
  q_w3 : Q.region ≠ W3.region
  r_w3 : R.region ≠ W3.region
  r_q : R.region ≠ Q.region
  r_j : R.region ≠ J.region
  g_j : G.region ≠ J.region

/-- Dynamic physical safety for every represented loop state.  It carries no
semantic output and is therefore stable under the lengths returned by the
actual checked HGCD execution. -/
def GcdHgcdDynamicDivremWorkspaceProvider (this : DenseUPolyZp)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3) : Prop :=
  ∀ (heap : RawHeap) (lenG lenJ : Nat)
    (left right : Polynomial (ZMod this._p.toNat)),
    RawDensePolyRep this heap G lenG left →
    RawDensePolyRep this heap J lenJ right →
    GcdHgcdDivremWorkspace heap G J Q R W3 lenG lenJ

/-- Dynamic local-vector safety for the small Euclid branch. -/
def GcdHgcdDynamicEuclidWorkspaceProvider (this : DenseUPolyZp)
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) : Prop :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat)),
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    ∃ capacity,
      GcdEuclidRawWorkspace heap G J R euclidA euclidB euclidQ euclidR
          euclidW3 capacity ∧
      lenJ ≤ capacity ∧ lenR ≤ capacity ∧ capacity < limbBase

/-- Concrete recursive safety package for the exact large-branch HGCD call
at one dynamically reached divisor/remainder state. -/
structure GcdHgcdDynamicHgcdWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ) : Type where
  invocation : HgcdRecursiveInvocationWorkspace this
    (hgcdRecursiveCallChecked this) lenJ M hM false G J J R lenJ lenR W
    scratch heap divisor remainder hcfg hp rfl horder
  descendants : HgcdRecursiveInvocationWorkspaceProviderBelow this
    (hgcdRecursiveCallChecked this) hcfg hp lenJ

/-- Dynamic physical provider for every represented large-branch state.  It
contains recursive allocation safety but neither an HGCD result nor a length
claim. -/
def GcdHgcdDynamicHgcdWorkspaceProvider (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64) : Type :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ),
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    GcdHgcdDynamicHgcdWorkspace this hcfg hp M hM G J R W scratch heap lenJ
      lenR divisor remainder horder

/-- The exact generated base arm cannot fault when its physical workspace is
valid.  This is the first total-execution component needed to upgrade checked
HGCD from conditional semantic correctness to total raw refinement. -/
theorem hgcdRecursiveBase_succeeds (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (workspace : HgcdRecursiveBaseCallWorkspace this M hM A B a b lenA lenB
      heap left right) :
    ∃ result, hgcdRecursiveBase M computeM A B a b lenA lenB heap =
      .ok result := by
  cases computeM with
  | false =>
      rcases copyU64_refines_rawDense this heap A a lenA left
          workspace.validA workspace.Aa workspace.leftRep with
        ⟨heap1, hcopyA, hlayout1, hARep⟩
      have hright1 :=
        CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
          this heap heap1 b lenB right hlayout1
        (copyU64_preserves_prefix heap heap1 A a b lenA lenB
          workspace.validA workspace.leftRep.1 workspace.Ab hcopyA)
        workspace.rightRep
      have hB1 : heap1.ValidU64Slice B lenB :=
        (hlayout1 B lenB).mp workspace.validB
      rcases copyU64_refines_rawDense this heap1 B b lenB right hB1
          workspace.Bb hright1 with ⟨heap2, hcopyB, _, _⟩
      exact ⟨⟨heap2, M, lenA, lenB, 1⟩, by
        simp [hgcdRecursiveBase, hcopyA, hcopyB]⟩
  | true =>
      rcases hgcdIterInit_refines this M A B A A 0 a lenA b lenB heap left
          right hM hp workspace.valid0 workspace.valid3 workspace.disjoint03
          workspace.validA workspace.validB workspace.Aa workspace.Bb
          workspace.Ab workspace.BA workspace.row0a workspace.row3a
          workspace.row0b workspace.row3b workspace.aMatrix workspace.bMatrix
          workspace.matrixValid workspace.leftRep workspace.rightRep with
        ⟨initial, hinit, _⟩
      refine ⟨initial.toRecursiveBaseResult, ?_⟩
      rw [hgcdRecursiveBase_true_eq_init, hinit]
      rfl

/-- Dynamic form of the loop-head divrem refinement. -/
theorem gcdHgcdLoop_dynamic_divrem_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenG lenJ : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdDynamicDivremWorkspaceProvider this G J Q R W3)
    (hlenJPos : 0 < lenJ)
    (hleft : RawDensePolyRep this heap G lenG left)
    (hright : RawDensePolyRep this heap J lenJ right)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R G lenG J
          lenJ W3 heap = .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      RawDensePolyRep this heap' J lenJ right ∧
      RawDensePolyRep this heap' R lenR remainder ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd right remainder) ∧
      RawHeap.SameLayout heap heap' ∧ lenR < lenJ := by
  have ws := physical heap lenG lenJ left right hleft hright
  rcases polyDivrem_next_state this Q R G J lenG lenJ W3 heap left right
      hlenJPos hleft hright ws.validQ ws.validR ws.validW3
      ws.quotientCapacity ws.r_g ws.w3_g ws.w3_j ws.q_j ws.q_w3 ws.r_w3
      ws.r_q ws.r_j hcfg with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep, hRRep,
      _, hgcd, hlayout, _, _, hlt⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep,
    hRRep, hgcd, hlayout, hlt⟩

/-- Dynamic form of the exact small raw-Euclid branch refinement. -/
theorem gcdHgcdLoop_dynamic_euclid_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (lenJ lenR : Nat) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdDynamicEuclidWorkspaceProvider this G J R euclidA
      euclidB euclidQ euclidR euclidW3)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenG result,
      strictGcdEuclidRaw this G J R euclidA euclidB euclidQ euclidR
          euclidW3 lenJ lenR heap = .ok ⟨heap', lenG⟩ ∧
      RawDensePolyRep this heap' G lenG result ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize result ∧
      RawHeap.SameLayout heap heap' := by
  rcases physical heap lenJ lenR divisor remainder hdivisor hremainder with
    ⟨capacity, hworkspace, hlenJ, hlenR, hcapacity⟩
  exact strictGcdEuclidRaw_refines this G J R euclidA euclidB euclidQ euclidR
    euclidW3 lenJ lenR capacity heap divisor remainder hworkspace hlenJ hlenR
    hcapacity hdivisor hremainder hcfg

/-- The actual loop-head divrem call produces the next represented pair and
preserves its normalized gcd. -/
theorem gcdHgcdLoop_divrem_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (capacity lenG lenJ : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdLoopDivremWorkspaceProvider this G J Q R W3 capacity)
    (hlenG : lenG ≤ capacity) (hlenJ : lenJ ≤ capacity)
    (hcapacity : capacity < limbBase) (hlenJPos : 0 < lenJ)
    (hleft : RawDensePolyRep this heap G lenG left)
    (hright : RawDensePolyRep this heap J lenJ right)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R G lenG J
          lenJ W3 heap = .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      RawDensePolyRep this heap' J lenJ right ∧
      RawDensePolyRep this heap' R lenR remainder ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd right remainder) ∧
      RawHeap.SameLayout heap heap' ∧ lenR < lenJ := by
  have ws := physical heap lenG lenJ left right hlenG hlenJ hleft hright
  have hqBound : lenG - (lenJ - 1) ≤ capacity :=
    (Nat.sub_le lenG (lenJ - 1)).trans hlenG
  have hrBound : Nat.min lenG (lenJ - 1) ≤ capacity :=
    (Nat.min_le_left lenG (lenJ - 1)).trans hlenG
  rcases polyDivrem_next_state this Q R G J lenG lenJ W3 heap left right
      hlenJPos hleft hright
      (heap.validU64Slice_mono Q capacity _ ws.validQ hqBound)
      (heap.validU64Slice_mono R capacity _ ws.validR hrBound)
      (heap.validU64Slice_mono (RawPtr.reinterpret W3) (3 * capacity)
        (3 * lenG) ws.validW3 (by omega))
      (hqBound.trans_lt hcapacity) ws.regions.a_r.symm ws.regions.w3_a
      ws.regions.w3_b ws.regions.q_b ws.regions.q_w3
      ws.regions.w3_r.symm ws.regions.q_r.symm ws.regions.b_r.symm hcfg with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep, hRRep,
      _, hgcd, hlayout, _, _, hlt⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep,
    hRRep, hgcd, hlayout, hlt⟩

/-- The small source branch is the already-refined exact raw Euclid helper on
the divisor/remainder pair produced by the loop-head division. -/
theorem gcdHgcdLoop_euclid_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (capacity lenJ lenR : Nat)
    (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdLoopEuclidWorkspaceProvider this G J R euclidA
      euclidB euclidQ euclidR euclidW3 capacity)
    (hlenJ : lenJ ≤ capacity) (hlenR : lenR ≤ capacity)
    (hcapacity : capacity < limbBase)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenG result,
      strictGcdEuclidRaw this G J R euclidA euclidB euclidQ euclidR
          euclidW3 lenJ lenR heap = .ok ⟨heap', lenG⟩ ∧
      RawDensePolyRep this heap' G lenG result ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize result ∧
      RawHeap.SameLayout heap heap' := by
  exact strictGcdEuclidRaw_refines this G J R euclidA euclidB euclidQ euclidR
    euclidW3 lenJ lenR capacity heap divisor remainder
    (physical heap lenJ lenR divisor remainder hlenJ hlenR hdivisor
      hremainder)
    hlenJ hlenR hcapacity hdivisor hremainder hcfg

/-- The zero-remainder source branch copies the current divisor into `G`; a
zero-length represented remainder is the zero polynomial, so this divisor is
the normalized gcd of the current pair. -/
theorem gcdHgcdLoop_zero_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R : RawPtr UInt64) (lenJ : Nat) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (hG : heap.ValidU64Slice G lenJ) (hGJ : G.region ≠ J.region)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R 0 remainder) :
    ∃ heap', heap.copyU64 G J lenJ = .ok heap' ∧
      RawDensePolyRep this heap' G lenJ divisor ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize divisor ∧
      RawHeap.SameLayout heap heap' := by
  have hremainderZero : remainder = 0 :=
    slicePolyRep_zero_length heap R this._p.toNat remainder
      hremainder.2.2.1
  subst remainder
  rcases copyU64_refines_rawDense_of_region_ne this heap G J lenJ divisor hG
      hGJ hdivisor with ⟨heap', hcopy, hlayout, hresult⟩
  exact ⟨heap', hcopy, hresult, by rw [EuclideanDomain.gcd_zero_right],
    hlayout⟩

/-- The actual result of a successful large-branch checked HGCD call strictly
decreases `_gcd_hgcd`'s source loop measure.  The bound is obtained from the
concrete well-founded HGCD invariant and the source cutoff, rather than from a
separate execution. -/
theorem hgcdRecursiveCallChecked_lenB_lt_of_large
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (horder : lenB < lenA)
    (hlarge : ¬ lenA < hgcdRecursiveCutoff)
    (workspace : HgcdRecursiveInvocationWorkspace this
      (hgcdRecursiveCallChecked this) lenA M hM false A B a b lenA lenB W
      scratch heap left right hcfg hp rfl horder)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this
      (hgcdRecursiveCallChecked this) hcfg hp lenA)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false A B a b lenA lenB W
      scratch heap = .ok result) :
    result.lenB < lenA := by
  rcases hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp M hM false A B
      a b lenA lenB W scratch heap left right horder workspace provider result
      hrun with ⟨finalA, finalB, entries, hinvariant⟩
  have hhalf : lenA / 2 + 1 < lenA := by
    simp only [hgcdRecursiveCutoff] at hlarge
    omega
  exact hinvariant.stopped.trans hhalf

/-- One successful large source branch simultaneously preserves the raw
polynomial pair's normalized gcd and decreases the exact loop length. -/
theorem hgcdRecursiveCallChecked_gcd_step_of_large
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G J R : RawPtr UInt64) (lenJ lenR : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ)
    (hlarge : ¬ lenJ < hgcdRecursiveCutoff)
    (workspace : HgcdRecursiveInvocationWorkspace this
      (hgcdRecursiveCallChecked this) lenJ M hM false G J J R lenJ lenR W
      scratch heap divisor remainder hcfg hp rfl horder)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this
      (hgcdRecursiveCallChecked this) hcfg hp lenJ)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false G J J R lenJ lenR W
      scratch heap = .ok result) :
    ∃ (finalG finalJ : Polynomial (ZMod this._p.toNat)),
      RawDensePolyRep this result.heap G result.lenA finalG ∧
      RawDensePolyRep this result.heap J result.lenB finalJ ∧
      normalize (EuclideanDomain.gcd divisor remainder) =
        normalize (EuclideanDomain.gcd finalG finalJ) ∧
      result.lenB < lenJ := by
  rcases hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp M hM false G J
      J R lenJ lenR W scratch heap divisor remainder horder workspace provider
      result hrun with ⟨finalG, finalJ, entries, hinvariant⟩
  refine ⟨finalG, finalJ, hinvariant.aRep, hinvariant.bRep,
    hinvariant.gcdPreserved, ?_⟩
  have hhalf : lenJ / 2 + 1 < lenJ := by
    simp only [hgcdRecursiveCutoff] at hlarge
    omega
  exact hinvariant.stopped.trans hhalf

/-- Reachable large-branch termination and semantics are derived from the
dynamic physical provider and the one actual checked HGCD result. -/
theorem gcdHgcdLoop_dynamic_hgcd_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64)
    (physical : GcdHgcdDynamicHgcdWorkspaceProvider this hcfg hp M hM G J R
      W scratch)
    (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ) (hlarge : ¬ lenJ < hgcdRecursiveCutoff)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false G J J R lenJ lenR W
      scratch heap = .ok result) :
    ∃ (finalG finalJ : Polynomial (ZMod this._p.toNat)),
      RawDensePolyRep this result.heap G result.lenA finalG ∧
      RawDensePolyRep this result.heap J result.lenB finalJ ∧
      normalize (EuclideanDomain.gcd divisor remainder) =
        normalize (EuclideanDomain.gcd finalG finalJ) ∧
      result.lenB < lenJ := by
  let node := physical heap lenJ lenR divisor remainder horder hdivisor
    hremainder
  exact hgcdRecursiveCallChecked_gcd_step_of_large this hcfg hp M hM G J R
    lenJ lenR W scratch heap divisor remainder horder hlarge node.invocation
    node.descendants result hrun

end CLPoly.Impl.StrictGCDHGCDRefinement
