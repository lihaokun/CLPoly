/- Total building blocks for the generated Hensel entry point. -/
import CLPoly.Generated.Corpus

set_option autoImplicit false

namespace Generated

/-- Total counterpart of `_loop___hensel_lift_upoly_0_ir`: compute `p^k`. -/
def henselTargetPowSafe (p : UInt64) : Nat → Int
  | 0 => 1
  | k + 1 => henselTargetPowSafe p k * (p.toNat : Int)

/-- The explicit-exponent branch of the generated Hensel target calculation. -/
def henselTargetSafe (p : UInt64) (k : Nat) : Int :=
  henselTargetPowSafe p k - 1

theorem henselTargetPowSafe_eq (p : UInt64) (k : Nat) :
    henselTargetPowSafe p k = (p.toNat : Int) ^ k := by
  induction k with
  | zero => rfl
  | succ k ih => simpa [henselTargetPowSafe, ih] using
      (Int.pow_succ (p.toNat : Int) k).symm

theorem henselTargetSafe_eq (p : UInt64) (k : Nat) :
    henselTargetSafe p k = (p.toNat : Int) ^ k - 1 := by
  rw [henselTargetSafe, henselTargetPowSafe_eq]

/-- Total counterpart of `_loop___hensel_lift_upoly_1_ir`. -/
def henselScaleFirstFactorSafe (factor : SparsePolyZp) (lc : Zp) : SparsePolyZp :=
  factor.map fun term => (term.1, term.2 * lc)

/-- Total counterpart of `_loop___hensel_lift_upoly_3_ir`. -/
def henselScaleFirstLiftSafe (factor : SparsePolyZZ) (lcInv : Int) : SparsePolyZZ :=
  factor.map fun term => (term.1, term.2 * lcInv)

end Generated
