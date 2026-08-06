/-
  Semantic refinement of the fixed-width scalar operations used by the raw
  dense-polynomial implementation.
-/
import CLPoly.Generated.StrictGCD
import Mathlib.Tactic

set_option autoImplicit false

namespace CLPoly.Impl.StrictWordArithmetic

open Generated.StrictGCD

/-- Mathematical base of one C++ `uint64_t` limb. -/
def limbBase : Nat := 2 ^ 64

/-
  Natural-language proof.

  Both operands are unsigned 64-bit words, hence each is below `2^64`; their
  product is therefore below `2^128`.  The generated multiplication first
  zero-extends both operands and multiplies them in `BitVec 128`, so its value
  is the ordinary natural-number product rather than a wrapped product.  The
  returned low word is the remainder modulo `2^64`, and shifting right by 64
  before truncation returns the quotient.  The quotient/remainder identity
  then reconstructs the original product.
-/
theorem umul128_reconstruct (a b : UInt64) :
    let out := dense_upoly_zp__umul128_ir 0 0 a b
    out.2.toNat + limbBase * out.1.toNat = a.toNat * b.toNat := by
  simp only [dense_upoly_zp__umul128_ir, uint128_of_uint64, uint128_lo,
    limbBase, BitVec.toNat_mul, BitVec.toNat_ofNat]
  have ha : a.toNat < 2 ^ 64 := UInt64.toNat_lt a
  have hb : b.toNat < 2 ^ 64 := UInt64.toNat_lt b
  have ha128 : a.toNat < 2 ^ 128 := lt_trans ha (by norm_num)
  have hb128 : b.toNat < 2 ^ 128 := lt_trans hb (by norm_num)
  have hab : a.toNat * b.toNat < 2 ^ 128 := by
    nlinarith [mul_nonneg (show 0 ≤ 2 ^ 64 - a.toNat by omega)
      (show 0 ≤ 2 ^ 64 - b.toNat by omega)]
  have hprod :
      (BitVec.ofNat 128 a.toNat * BitVec.ofNat 128 b.toNat).toNat =
        a.toNat * b.toNat := by
    rw [BitVec.toNat_mul]
    simp only [BitVec.toNat_ofNat]
    rw [Nat.mod_eq_of_lt ha128, Nat.mod_eq_of_lt hb128,
      Nat.mod_eq_of_lt hab]
  have hshift :
      ((BitVec.ofNat 128 a.toNat * BitVec.ofNat 128 b.toNat) >>>
        (64 : UInt128)).toNat =
        a.toNat * b.toNat / 2 ^ 64 := by
    simp [hprod, Nat.shiftRight_eq_div_pow]
  have hhi : a.toNat * b.toNat / 2 ^ 64 < 2 ^ 64 := by
    rw [Nat.div_lt_iff_lt_mul (by positivity)]
    simpa [← pow_add] using hab
  simp only [Nat.mod_eq_of_lt ha128, Nat.mod_eq_of_lt hb128,
    Nat.mod_eq_of_lt hab, UInt64.toNat_ofNat']
  rw [hshift, Nat.mod_eq_of_lt hhi]
  exact Nat.mod_add_div (a.toNat * b.toNat) (2 ^ 64)

end CLPoly.Impl.StrictWordArithmetic
