import CLPoly.Impl.StrictKarMulRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep

/-- The actual zero-padding loop extends a represented prefix without
changing its L2 polynomial.  Recursion follows the generated loop index. -/
theorem mulZeroPadLoop_refines (bPad : RawPtr UInt64)
    (start count i p : Nat) (heap : RawHeap)
    (poly : Polynomial (ZMod p)) (modulus : UInt64)
    (hi : i ≤ count) (hmodulus : modulus ≠ 0)
    (hPad : heap.ValidU64Slice bPad (start + count))
    (hRep : SlicePolyRep heap bPad (start + i) p poly)
    (hCanonical : CanonicalU64Prefix heap bPad (start + i) modulus) :
    ∃ heap', mulZeroPadLoop bPad start count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' bPad (start + count) p poly ∧
      CanonicalU64Prefix heap' bPad (start + count) modulus := by
  rw [mulZeroPadLoop]
  split
  next hmore =>
    have hThroughNext := heap.validU64Slice_mono bPad (start + count)
      ((start + i) + 1) hPad (by omega)
    rcases heap.writeU64_of_valid bPad (start + count) (start + i) 0 hPad
      (by omega) with ⟨heap1, hwrite⟩
    simp only [hwrite]
    rcases writeZero_extends_slice heap heap1 bPad (start + i) p modulus
      poly hmodulus hThroughNext hRep hCanonical hwrite with
      ⟨hlayout1, hRep1, hCanonical1⟩
    have hPad1 := (hlayout1 bPad (start + count)).mp hPad
    rcases mulZeroPadLoop_refines bPad start count (i + 1) p heap1 poly
      modulus (by omega) hmodulus hPad1 (by simpa [Nat.add_assoc] using hRep1)
      (by simpa [Nat.add_assoc] using hCanonical1) with
      ⟨heap2, hrun, hlayout2, hRep2, hCanonical2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length), hRep2, hCanonical2⟩
  next hdone =>
    have hieq : i = count := by omega
    subst i
    exact ⟨heap, rfl, fun _ _ => Iff.rfl, hRep, hCanonical⟩
termination_by count - i
decreasing_by omega

end CLPoly.Impl.StrictMulRefinement
