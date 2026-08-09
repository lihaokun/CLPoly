import CLPoly.Generated.StrictGCDHGCD
import CLPoly.Impl.StrictHGCDCheckedRefinement
import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictGCDHGCDRefinement

open Generated.StrictHGCD
open Generated.StrictGCDHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictHGCDRawRefinement
open CLPoly.Impl.StrictWordArithmetic

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

end CLPoly.Impl.StrictGCDHGCDRefinement
