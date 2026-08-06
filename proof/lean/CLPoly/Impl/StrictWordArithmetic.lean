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

/-- Little-endian natural-number observation of the C++ three-limb word. -/
def word3Value (s : Word3) : Nat :=
  s.lo.toNat + limbBase * s.mid.toNat + limbBase ^ 2 * s.hi.toNat

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

/-
  Natural-language proof.

  Any 128-bit word is below `2^128`.  Truncation to `UInt64` returns its
  remainder modulo `2^64`; shifting right by 64 returns its quotient by that
  base, which is already below `2^64`, so the second truncation is exact.  The
  natural-number quotient/remainder identity reconstructs the original word.
  This lemma is the common carry splitter used twice by `_add_carry3`.
-/
theorem uint128_split_reconstruct (x : UInt128) :
    (uint128_lo x).toNat + limbBase *
      (uint128_lo (x >>> (64 : UInt128))).toNat = x.toNat := by
  have hx : x.toNat < 2 ^ 128 :=
    BitVec.toNat_lt_twoPow_of_le (x := x) (Nat.le_refl 128)
  have hshift : (x >>> (64 : UInt128)).toNat = x.toNat / 2 ^ 64 := by
    simp [Nat.shiftRight_eq_div_pow]
  have hhi : x.toNat / 2 ^ 64 < 2 ^ 64 := by
    rw [Nat.div_lt_iff_lt_mul (by positivity)]
    simpa [← pow_add] using hx
  simp only [uint128_lo, limbBase, UInt64.toNat_ofNat']
  rw [hshift, Nat.mod_eq_of_lt hhi]
  exact Nat.mod_add_div x.toNat (2 ^ 64)

/-
  Natural-language proof.

  Apply `uint128_split_reconstruct` to the two intermediate 128-bit sums.  The
  first identity replaces low-limb addition by a low result plus one base times
  its carry; the second does the same at the middle limb.  After multiplying
  the second identity by the limb base, the carry terms telescope.  The final
  `UInt64` addition wraps the high limb modulo the base, which is precisely a
  wrap modulo `base^3` after weighting that limb by `base^2`.
-/
theorem addCarry3_modEq (s : Word3) (b1 b0 : UInt64) :
    let out := dense_upoly_zp__add_carry3_ir s b1 b0
    Nat.ModEq (limbBase ^ 3) (word3Value out)
      (word3Value s + b0.toNat + limbBase * b1.toNat) := by
  let sum0 : UInt128 := uint128_of_uint64 s.lo + uint128_of_uint64 b0
  let carry0 : UInt64 := uint128_lo (sum0 >>> (64 : UInt128))
  let sum1 : UInt128 := uint128_of_uint64 s.mid + uint128_of_uint64 b1 +
    uint128_of_uint64 carry0
  let carry1 : UInt64 := uint128_lo (sum1 >>> (64 : UInt128))
  have hslo := UInt64.toNat_lt s.lo
  have hsmid := UInt64.toNat_lt s.mid
  have hb0 := UInt64.toNat_lt b0
  have hb1 := UInt64.toNat_lt b1
  have hcarry0 := UInt64.toNat_lt carry0
  have hx0 : s.lo.toNat + b0.toNat < 2 ^ 128 := by omega
  have hx1 : s.mid.toNat + b1.toNat + carry0.toNat < 2 ^ 128 := by omega
  norm_num at hx0 hx1
  have hsum0 : sum0.toNat = s.lo.toNat + b0.toNat := by
    simp [sum0, uint128_of_uint64, BitVec.toNat_add, BitVec.toNat_ofNat,
      Nat.mod_eq_of_lt hx0]
  have hsum1 : sum1.toNat = s.mid.toNat + b1.toNat + carry0.toNat := by
    simp [sum1, uint128_of_uint64, BitVec.toNat_add, BitVec.toNat_ofNat,
      Nat.mod_eq_of_lt hx1]
  have hsplit0 := uint128_split_reconstruct sum0
  have hsplit1 := uint128_split_reconstruct sum1
  have hcarry0eq : carry0.toNat =
      (uint128_lo (sum0 >>> (64 : UInt128))).toNat := rfl
  have hcarry1eq : carry1.toNat =
      (uint128_lo (sum1 >>> (64 : UInt128))).toNat := rfl
  have hhi : (s.hi + carry1).toNat =
      (s.hi.toNat + carry1.toNat) % 2 ^ 64 := UInt64.toNat_add _ _
  simp only [dense_upoly_zp__add_carry3_ir, word3_addCarry_x86,
    word3Value, limbBase]
  change Nat.ModEq ((2 ^ 64) ^ 3)
    ((uint128_lo sum0).toNat + 2 ^ 64 * (uint128_lo sum1).toNat +
      (2 ^ 64) ^ 2 * (s.hi + carry1).toNat)
    (s.lo.toNat + 2 ^ 64 * s.mid.toNat + (2 ^ 64) ^ 2 * s.hi.toNat +
      b0.toNat + 2 ^ 64 * b1.toNat)
  simp only [limbBase] at hsplit0 hsplit1
  norm_num [Nat.ModEq] at hsplit0 hsplit1 hcarry0eq hcarry1eq hhi ⊢
  rw [hsum0] at hsplit0
  rw [hsum1] at hsplit1
  omega

/-
  Natural-language proof.

  The generated numerator places `~pn` in the high limb and an all-one word
  in the low limb.  As a 128-bit bitvector this is exactly the all-one
  128-bit word minus `pn * 2^64`.  Since `pn < 2^64`, that subtraction is
  nonnegative and the shifted product is below `2^128`; converting to `Nat`
  therefore yields `2^128 - 1 - pn*2^64` without modular wrap.
-/
theorem preinvert_numerator_value (pn : UInt64) :
    let num : UInt128 :=
      ((uint128_of_uint64 (~~~pn)) <<< (64 : UInt128)) |||
        uint128_of_uint64 (~~~(0 : UInt64))
    num.toNat = limbBase ^ 2 - 1 - pn.toNat * limbBase := by
  let num : UInt128 :=
    ((uint128_of_uint64 (~~~pn)) <<< (64 : UInt128)) |||
      uint128_of_uint64 (~~~(0 : UInt64))
  change num.toNat = limbBase ^ 2 - 1 - pn.toNat * limbBase
  have hpn := UInt64.toNat_lt pn
  have h64 : (64 : UInt128).toNat = 64 := by decide
  have hhigh :
      (2 ^ 64 - 1 - pn.toNat) <<< 64 < 2 ^ 128 := by
    rw [Nat.shiftLeft_eq]
    calc
      (2 ^ 64 - 1 - pn.toNat) * 2 ^ 64 < (2 ^ 64) * 2 ^ 64 :=
        Nat.mul_lt_mul_of_pos_right (by omega) (by positivity)
      _ = 2 ^ 128 := by rw [← pow_add]
  have hlow : 2 ^ 64 - 1 < 2 ^ 64 := by omega
  have hor := Nat.shiftLeft_add_eq_or_of_lt hlow
    (2 ^ 64 - 1 - pn.toNat)
  dsimp [num, uint128_of_uint64]
  simp only [BitVec.toNat_shiftLeft,
    BitVec.toNat_ofNat, UInt64.toNat_not, UInt64.size]
  norm_num [Nat.shiftLeft_eq, limbBase] at hhigh hor hpn ⊢
  rw [Nat.mod_eq_of_lt hhigh]
  rw [← hor]
  omega

/-
  Natural-language proof.

  By `preinvert_numerator_value`, the generated dividend is
  `N = B^2 - 1 - pn*B`.  Normalization gives `pn ≥ B/2`, hence `N < pn*B`, so
  `N/pn < B` and truncating the quotient to one limb is exact.  Finally
  `B^2 - 1 = N + pn*B`; division by `pn` and `Nat.add_mul_div_left` give
  `(B^2 - 1)/pn = N/pn + B`, which is the FLINT preinverse equation.
-/
theorem preinvert_limb_spec (pn : UInt64)
    (hnorm : limbBase / 2 ≤ pn.toNat) :
    limbBase + (dense_upoly_zp___preinvert_limb_ir pn).toNat =
      (limbBase ^ 2 - 1) / pn.toNat := by
  have hpnlt : pn.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt pn
  have hpnpos : 0 < pn.toNat := by
    have hhalf : 0 < limbBase / 2 := by norm_num [limbBase]
    omega
  let num : UInt128 :=
    ((uint128_of_uint64 (~~~pn)) <<< (64 : UInt128)) |||
      uint128_of_uint64 (~~~(0 : UInt64))
  have hnum : num.toNat = limbBase ^ 2 - 1 - pn.toNat * limbBase :=
    preinvert_numerator_value pn
  have hdecomp :
      limbBase ^ 2 - 1 = num.toNat + pn.toNat * limbBase := by
    rw [hnum]
    norm_num [limbBase] at hpnlt ⊢
    omega
  have hnumlt : num.toNat < pn.toNat * limbBase := by
    rw [hnum]
    norm_num [limbBase] at hnorm hpnlt ⊢
    omega
  have hquotlt : num.toNat / pn.toNat < limbBase := by
    rw [Nat.div_lt_iff_lt_mul hpnpos]
    simpa [Nat.mul_comm] using hnumlt
  change limbBase + (uint128_lo (num / uint128_of_uint64 pn)).toNat = _
  simp only [uint128_lo, uint128_of_uint64, UInt64.toNat_ofNat',
    BitVec.toNat_udiv, BitVec.toNat_ofNat]
  rw [Nat.mod_eq_of_lt (lt_trans hpnlt (by norm_num [limbBase]))]
  change limbBase + (num.toNat / pn.toNat) % limbBase = _
  rw [Nat.mod_eq_of_lt hquotlt]
  calc
    limbBase + num.toNat / pn.toNat = num.toNat / pn.toNat + limbBase :=
      Nat.add_comm _ _
    _ = (num.toNat + pn.toNat * limbBase) / pn.toNat := by
      rw [Nat.add_mul_div_left num.toNat limbBase hpnpos]
    _ = (limbBase ^ 2 - 1) / pn.toNat := by rw [hdecomp]

end CLPoly.Impl.StrictWordArithmetic
