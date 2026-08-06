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

/-- The repeated two-limb quotient-estimate block in `_lll_mod_preinv`. -/
def preinvQuotientPair (u1 u0 pinv : UInt64) : UInt64 × UInt64 :=
  let qm : UInt128 := uint128_of_uint64 u1 * uint128_of_uint64 pinv
  let q1 : UInt64 := uint128_lo (qm >>> (64 : UInt128))
  let q0 : UInt64 := uint128_lo qm
  let q0' : UInt64 := q0 + u0
  let carry : UInt64 := if q0' < u0 then 1 else 0
  let q1' : UInt64 := q1 + (u1 + carry)
  (q1', q0')

/-- One source-level normalized reduction round, retaining the generated
`Int32 → Int64 → UInt64` carry conversion verbatim.  This definition is a
factoring of the C++ machine operations, not a mathematical oracle. -/
def preinvRoundIR (u1 u0 pn pinv : UInt64) : UInt64 :=
  let qm : UInt128 := uint128_of_uint64 u1 * uint128_of_uint64 pinv
  let q1 : UInt64 := uint128_lo (qm >>> (64 : UInt128))
  let q0 : UInt64 := uint128_lo qm
  let q0' : UInt64 := q0 + u0
  let carry : UInt64 :=
    ((if q0' < u0 then (1 : Int32) else 0).toInt64.toUInt64)
  let q1' : UInt64 := q1 + (u1 + carry)
  let r : UInt64 := u0 - ((q1' + 1) * pn)
  let r' : UInt64 := if r > q0' then r + pn else r
  if r' >= pn then r' - pn else r'

/-- The same source round with its quotient/carry prefix shared through
`preinvQuotientPair`; this is the form used by the arithmetic invariant. -/
def preinvReduceNormalized (u1 u0 pn pinv : UInt64) : UInt64 :=
  let q := preinvQuotientPair u1 u0 pinv
  let r : UInt64 := u0 - ((q.1 + 1) * pn)
  let r' : UInt64 := if r > q.2 then r + pn else r
  if r' >= pn then r' - pn else r'

@[simp] theorem uint32_zero_toUInt64 : (0 : UInt32).toUInt64 = 0 := by
  decide

@[simp] theorem uint64_shiftLeft_zero (x : UInt64) :
    x <<< (0 : UInt64) = x := by
  apply UInt64.toNat_inj.mp
  change (x <<< (0 : UInt64)).toNat = x.toNat
  simp

@[simp] theorem uint64_shiftLeft_zero_u32 (x : UInt64) :
    x <<< (0 : UInt32) = x := by
  apply UInt64.toNat_inj.mp
  change (x <<< (0 : UInt64)).toNat = x.toNat
  simp

@[simp] theorem uint64_shiftRight_zero (x : UInt64) :
    x >>> (0 : UInt64) = x := by
  apply UInt64.toNat_inj.mp
  change (x >>> (0 : UInt64)).toNat = x.toNat
  simp

@[simp] theorem uint64_shiftRight_zero_u32 (x : UInt64) :
    x >>> (0 : UInt32) = x := by
  apply UInt64.toNat_inj.mp
  change (x >>> (0 : UInt64)).toNat = x.toNat
  simp

theorem int32Carry_toUInt64 (c : Prop) [Decidable c] :
    ((if c then (1 : Int32) else 0).toInt64.toUInt64) =
      (if c then (1 : UInt64) else 0) := by
  by_cases h : c <;> simp [h] <;> decide

theorem int64_one_toUInt64 : (1 : Int64).toUInt64 = 1 := by
  decide

theorem int64_zero_toUInt64 : (0 : Int64).toUInt64 = 0 := by
  decide

/-
  Natural-language proof.

  The two definitions differ only in how the carry bit is written.  The C++
  translation passes `0` or `1` through `Int32`, `Int64`, and `UInt64`, while
  `preinvQuotientPair` constructs the same UInt64 bit directly.  Splitting on
  the source comparison makes both conversions concrete; every subsequent
  fixed-width operation and correction branch is then identical.
-/
theorem preinvRoundIR_eq_reduceNormalized (u1 u0 pn pinv : UInt64) :
    preinvRoundIR u1 u0 pn pinv =
      preinvReduceNormalized u1 u0 pn pinv := by
  by_cases hcarry :
      uint128_lo (uint128_of_uint64 u1 * uint128_of_uint64 pinv) + u0 < u0
  · simp [preinvRoundIR, preinvReduceNormalized, preinvQuotientPair, hcarry,
      int64_one_toUInt64]
  · simp [preinvRoundIR, preinvReduceNormalized, preinvQuotientPair, hcarry,
      int64_zero_toUInt64]

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

/-
  Natural-language proof.

  Put `x = B^2 - 1` and `m = x / pn`, using `preinvert_limb_spec` to identify
  `m` with `B + pinv`.  Euclidean division gives `m*pn ≤ x`.  Since
  `x < B^2`, this is the strict lower approximation.  Also `m < m+1`; the
  division comparison theorem yields `x < (m+1)*pn`, hence
  `B^2 ≤ (m+1)*pn`.  These are exactly the two bounds consumed by the
  Granlund--Möller reduction step.
-/
theorem preinverse_mul_bounds (pn : UInt64)
    (hnorm : limbBase / 2 ≤ pn.toNat) :
    let pinv := dense_upoly_zp___preinvert_limb_ir pn
    let m := limbBase + pinv.toNat
    m * pn.toNat < limbBase ^ 2 ∧
      limbBase ^ 2 ≤ (m + 1) * pn.toNat := by
  have hspec := preinvert_limb_spec pn hnorm
  have hpnpos : 0 < pn.toNat := by
    have hhalf : 0 < limbBase / 2 := by norm_num [limbBase]
    omega
  let pinv := dense_upoly_zp___preinvert_limb_ir pn
  let m := limbBase + pinv.toNat
  have hm : m = (limbBase ^ 2 - 1) / pn.toNat := by
    simpa [m, pinv] using hspec
  have hlo : m * pn.toNat ≤ limbBase ^ 2 - 1 := by
    rw [hm]
    exact Nat.div_mul_le_self _ _
  have hupper : limbBase ^ 2 - 1 < (m + 1) * pn.toNat := by
    rw [← Nat.div_lt_iff_lt_mul hpnpos, ← hm]
    exact Nat.lt_succ_self m
  change m * pn.toNat < limbBase ^ 2 ∧
    limbBase ^ 2 ≤ (m + 1) * pn.toNat
  have hBpos : 0 < limbBase ^ 2 := by norm_num [limbBase]
  exact ⟨by omega, by omega⟩

/-
  Natural-language proof.

  For nonzero `p`, the count-leading-zeros result `n` is strictly below 64.
  The bit immediately following those `n` zeroes is set, giving
  `2^(63-n) ≤ p`; the definition of `clz` also gives `p < 2^(64-n)`.
  Multiplying both inequalities by `2^n` and using `n ≤ 63` yields
  `2^63 ≤ p*2^n < 2^64`.  The strict upper bound also proves that the C++
  left shift is an exact natural-number multiplication, not a wrapped value.
-/
theorem clz_normalizes_uint64 (p : UInt64) (hp : p ≠ 0) :
    let n := p.toBitVec.clz.toNat
    n < 64 ∧ limbBase / 2 ≤ p.toNat * 2 ^ n ∧
      p.toNat * 2 ^ n < limbBase := by
  let n := p.toBitVec.clz.toNat
  have hpBits : p.toBitVec ≠ (0#64) := by
    intro h
    apply hp
    exact UInt64.toNat_inj.mp (by simpa using congrArg BitVec.toNat h)
  have hn : n < 64 := by
    simpa [n] using (BitVec.clz_lt_iff_ne_zero (x := p.toBitVec)).2 hpBits
  have hlo0 : 2 ^ (64 - 1 - n) ≤ p.toNat := by
    simpa [n] using BitVec.two_pow_sub_clz_le_toNat_of_ne_zero
      (x := p.toBitVec) (by omega : 0 < 64) hpBits
  have hhi0 : p.toNat < 2 ^ (64 - n) := by
    simpa [n] using BitVec.toNat_lt_two_pow_sub_clz (x := p.toBitVec)
  have hlo : 2 ^ 63 ≤ p.toNat * 2 ^ n := by
    have := Nat.mul_le_mul_right (2 ^ n) hlo0
    rw [← pow_add] at this
    have hexp : 64 - 1 - n + n = 63 := by omega
    rw [hexp] at this
    exact this
  have hhi : p.toNat * 2 ^ n < 2 ^ 64 := by
    have := Nat.mul_lt_mul_of_pos_right hhi0 (by positivity : 0 < 2 ^ n)
    rw [← pow_add] at this
    have hexp : 64 - n + n = 64 := by omega
    rw [hexp] at this
    exact this
  exact ⟨hn, by simpa [limbBase] using hlo, by simpa [limbBase] using hhi⟩

/-- Exact machine value assigned to `_norm` by dense C++ precomputation. -/
def denseNorm (p : UInt64) : UInt32 := UInt32.ofNat p.toBitVec.clz.toNat

/-- Exact normalized limb `p << _norm` used by generated preinverse code. -/
def denseNormalizedModulus (p : UInt64) : UInt64 := p <<< denseNorm p

/-
  Natural-language proof.

  `clz_normalizes_uint64` supplies `n<64` and the no-overflow product bound.
  Therefore converting `n` through `UInt32` and then through the generated
  shift instance preserves it exactly.  The machine left shift consequently
  has natural value `p*2^n`, which lies in the normalized half-open interval.
-/
theorem denseNormalizedModulus_spec (p : UInt64) (hp : p ≠ 0) :
    let n := p.toBitVec.clz.toNat
    (denseNorm p).toNat = n ∧
      (denseNormalizedModulus p).toNat = p.toNat * 2 ^ n ∧
      limbBase / 2 ≤ (denseNormalizedModulus p).toNat := by
  let n := p.toBitVec.clz.toNat
  have hnorm := clz_normalizes_uint64 p hp
  have hn : n < 64 := hnorm.1
  have hn32 : n < UInt32.size := lt_trans hn (by norm_num [UInt32.size])
  have hnormNat : (denseNorm p).toNat = n := by
    simp [denseNorm, n, UInt32.toNat_ofNat', Nat.mod_eq_of_lt hn32]
  have hprodlt : p.toNat * 2 ^ n < limbBase := hnorm.2.2
  have hshift : (denseNormalizedModulus p).toNat = p.toNat * 2 ^ n := by
    unfold denseNormalizedModulus
    change (p <<< (denseNorm p).toUInt64).toNat = p.toNat * 2 ^ n
    rw [UInt64.toNat_shiftLeft]
    have hprodlt' : p.toNat * 2 ^ n < UInt64.size := by
      simpa [limbBase, UInt64.size] using hprodlt
    have hnorm64 : (denseNorm p).toUInt64.toNat = n := by
      rw [UInt32.toNat_toUInt64]
      exact hnormNat
    rw [hnorm64, Nat.mod_eq_of_lt hn, Nat.shiftLeft_eq]
    exact Nat.mod_eq_of_lt hprodlt'
  exact ⟨hnormNat, hshift, hshift ▸ hnorm.2.1⟩

/-
  Natural-language proof.

  Split the exact 128-bit product `u1*pinv` into `(q1,q0)`.  Adding `u0` to
  the low limb produces `q0'` and the comparison `q0'<u0` is exactly its carry
  bit.  Adding `u1+carry` to the high limb completes a two-limb addition.
  Reconstructing both limbs cancels the carry; only overflow beyond the high
  limb is discarded, hence equality modulo `B^2`.
-/
theorem preinvQuotientPair_modEq (u1 u0 pinv : UInt64) :
    let out := preinvQuotientPair u1 u0 pinv
    Nat.ModEq (limbBase ^ 2)
      (out.2.toNat + limbBase * out.1.toNat)
      (u1.toNat * pinv.toNat + u0.toNat + limbBase * u1.toNat) := by
  let qm : UInt128 := uint128_of_uint64 u1 * uint128_of_uint64 pinv
  let q1 : UInt64 := uint128_lo (qm >>> (64 : UInt128))
  let q0 : UInt64 := uint128_lo qm
  let q0' : UInt64 := q0 + u0
  let carry : UInt64 := if q0' < u0 then 1 else 0
  let q1' : UInt64 := q1 + (u1 + carry)
  have hprod : q0.toNat + limbBase * q1.toNat =
      u1.toNat * pinv.toNat := by
    simpa [qm, q0, q1, dense_upoly_zp__umul128_ir] using
      umul128_reconstruct u1 pinv
  have hq0 := UInt64.toNat_lt q0
  have hu0 := UInt64.toNat_lt u0
  have hu1 := UInt64.toNat_lt u1
  norm_num [limbBase] at hq0 hu0 hu1
  have hcarry : carry.toNat = (q0.toNat + u0.toNat) / limbBase := by
    by_cases hov : q0' < u0
    · have hwrap : limbBase ≤ q0.toNat + u0.toNat := by
        simp only [q0', UInt64.lt_iff_toNat_lt, UInt64.toNat_add] at hov
        norm_num [limbBase] at hov ⊢
        omega
      have hsumlt : q0.toNat + u0.toNat < 2 * limbBase := by
        norm_num [limbBase]
        omega
      simp [carry, hov, UInt64.toNat_ofNat, limbBase]
      norm_num [limbBase] at hwrap hsumlt ⊢
      omega
    · have hnowrap : q0.toNat + u0.toNat < limbBase := by
        simp only [q0', UInt64.lt_iff_toNat_lt, UInt64.toNat_add] at hov
        norm_num [limbBase] at hov ⊢
        omega
      simp [carry, hov, Nat.div_eq_of_lt hnowrap]
  have hq0' : q0'.toNat = (q0.toNat + u0.toNat) % limbBase := by
    simp [q0', UInt64.toNat_add, limbBase]
  have hq1' : q1'.toNat =
      (q1.toNat + u1.toNat + carry.toNat) % limbBase := by
    simp [q1', UInt64.toNat_add, limbBase, Nat.add_assoc]
  change Nat.ModEq (limbBase ^ 2)
    (q0'.toNat + limbBase * q1'.toNat)
    (u1.toNat * pinv.toNat + u0.toNat + limbBase * u1.toNat)
  simp only [limbBase] at hprod hq0 hu0 hu1 hcarry hq0' hq1' ⊢
  norm_num [Nat.ModEq] at hprod hq0 hu0 hu1 hcarry hq0' hq1' ⊢
  omega

/-
  Natural-language proof.

  The low output is the low product limb plus `u0`, hence the remainder of
  `u1*pinv+u0` modulo `B`.  Its wrap comparison is exactly the quotient (zero
  or one) of that limb addition by `B`.  Adding the product high limb, `u1`,
  and that carry therefore gives the quotient of `u1*pinv+u0` by `B`, plus
  `u1`, with the final UInt64 wrap represented by `% B`.
-/
theorem preinvQuotientPair_components (u1 u0 pinv : UInt64) :
    let out := preinvQuotientPair u1 u0 pinv
    out.2.toNat = (u1.toNat * pinv.toNat + u0.toNat) % limbBase ∧
      out.1.toNat =
        ((u1.toNat * pinv.toNat + u0.toNat) / limbBase + u1.toNat) %
          limbBase := by
  let qm : UInt128 := uint128_of_uint64 u1 * uint128_of_uint64 pinv
  let q1 : UInt64 := uint128_lo (qm >>> (64 : UInt128))
  let q0 : UInt64 := uint128_lo qm
  let q0' : UInt64 := q0 + u0
  let carry : UInt64 := if q0' < u0 then 1 else 0
  let q1' : UInt64 := q1 + (u1 + carry)
  have hprod : q0.toNat + limbBase * q1.toNat =
      u1.toNat * pinv.toNat := by
    simpa [qm, q0, q1, dense_upoly_zp__umul128_ir] using
      umul128_reconstruct u1 pinv
  have hq0 := UInt64.toNat_lt q0
  have hu0 := UInt64.toNat_lt u0
  norm_num [limbBase] at hq0 hu0
  have hcarry : carry.toNat = (q0.toNat + u0.toNat) / limbBase := by
    by_cases hov : q0' < u0
    · have hwrap : limbBase ≤ q0.toNat + u0.toNat := by
        simp only [q0', UInt64.lt_iff_toNat_lt, UInt64.toNat_add] at hov
        norm_num [limbBase] at hov ⊢
        omega
      have hsumlt : q0.toNat + u0.toNat < 2 * limbBase := by
        norm_num [limbBase]
        omega
      simp [carry, hov, UInt64.toNat_ofNat, limbBase]
      norm_num [limbBase] at hwrap hsumlt ⊢
      omega
    · have hnowrap : q0.toNat + u0.toNat < limbBase := by
        simp only [q0', UInt64.lt_iff_toNat_lt, UInt64.toNat_add] at hov
        norm_num [limbBase] at hov ⊢
        omega
      simp [carry, hov, Nat.div_eq_of_lt hnowrap]
  have hlow : q0'.toNat =
      (u1.toNat * pinv.toNat + u0.toNat) % limbBase := by
    simp only [q0', UInt64.toNat_add]
    rw [← hprod]
    simp [limbBase, Nat.add_mod]
  have hquot : (u1.toNat * pinv.toNat + u0.toNat) / limbBase =
      q1.toNat + carry.toNat := by
    rw [← hprod]
    rw [show q0.toNat + limbBase * q1.toNat + u0.toNat =
      (q0.toNat + u0.toNat) + limbBase * q1.toNat by omega]
    rw [Nat.add_mul_div_left _ _ (by norm_num [limbBase])]
    omega
  have hhigh : q1'.toNat =
      ((u1.toNat * pinv.toNat + u0.toNat) / limbBase + u1.toNat) %
        limbBase := by
    have hq1lt : q1.toNat < limbBase := by
      simpa [limbBase] using UInt64.toNat_lt q1
    calc
      q1'.toNat =
          (q1.toNat + (u1.toNat + carry.toNat) % limbBase) % limbBase := by
            simp [q1', UInt64.toNat_add, limbBase]
      _ = (q1.toNat + (u1.toNat + carry.toNat)) % limbBase := by
            rw [Nat.add_mod]
            simp [Nat.mod_eq_of_lt hq1lt]
      _ = ((u1.toNat * pinv.toNat + u0.toNat) / limbBase + u1.toNat) %
          limbBase := by rw [hquot]; congr 1 <;> omega
  change q0'.toNat = _ ∧ q1'.toNat = _
  exact ⟨hlow, hhigh⟩

/-
  Natural-language proof.

  Write `m = B + pinv`.  The quotient estimate is exactly
  `(u1*m+u0)/B`: adding `u1` after division accounts for the `u1*B` part.
  Since `u1 < pn`, multiplying `u1+1 ≤ pn` by `m` and using
  `m*pn < B²` leaves at least `m ≥ B` below `B²`; the low limb `u0 < B`
  fits in that gap.  Thus the numerator is below `B²`, so the estimate is
  strictly below `B` and its UInt64 representation cannot wrap.
-/
theorem preinv_estimate_lt_base (B pn pinv u1 u0 : Nat)
    (hB : 0 < B) (hu0 : u0 < B) (hu1 : u1 < pn)
    (hmul : (B + pinv) * pn < B ^ 2) :
    (u1 * pinv + u0) / B + u1 < B := by
  let m := B + pinv
  have humul : (u1 + 1) * m ≤ pn * m :=
    Nat.mul_le_mul_right m (Nat.succ_le_iff.mpr hu1)
  have hmcomm : pn * m < B ^ 2 := by
    simpa [m, Nat.mul_comm] using hmul
  have hmB : B ≤ m := by simp [m]
  have hnum : u1 * m + u0 < B ^ 2 := by
    nlinarith
  have hrewrite : (u1 * pinv + u0) / B + u1 =
      (u1 * m + u0) / B := by
    rw [show u1 * m + u0 = (u1 * pinv + u0) + B * u1 by
      simp [m]; ring]
    rw [Nat.add_mul_div_left _ _ hB]
  rw [hrewrite]
  exact (Nat.div_lt_iff_lt_mul hB).2 (by simpa [pow_two] using hnum)

/-
  Natural-language proof.

  The component theorem says that the machine high limb is the mathematical
  estimate modulo `B`.  Under the normalized preinverse multiplication bound
  and `u1 < pn`, the preceding theorem proves that estimate is already below
  `B`; consequently the modulo operation is the identity.  This rules out the
  high-limb wrap that would otherwise invalidate the quotient invariant.
-/
theorem preinvQuotientPair_high_noWrap (u1 u0 pn pinv : UInt64)
    (hu1 : u1.toNat < pn.toNat)
    (hmul : (limbBase + pinv.toNat) * pn.toNat < limbBase ^ 2) :
    (preinvQuotientPair u1 u0 pinv).1.toNat =
      (u1.toNat * pinv.toNat + u0.toNat) / limbBase + u1.toNat := by
  have hest := preinv_estimate_lt_base limbBase pn.toNat pinv.toNat
    u1.toNat u0.toNat (by norm_num [limbBase])
    (by simpa [limbBase] using UInt64.toNat_lt u0) hu1 hmul
  have hcomp := (preinvQuotientPair_components u1 u0 pinv).2
  simpa [Nat.mod_eq_of_lt hest] using hcomp

/-
  Natural-language proof.

  The generated C++ multiplication first shifts `a`, forms its exact UInt128
  product with `b`, and splits that product into high and low limbs.  The
  remaining quotient estimate, carry conversion, two correction branches and
  final denormalising shift are exactly `preinvRoundIR`.  Unfolding both
  definitions therefore gives the same fixed-width program; no `% p`, L2
  multiplication, or specification oracle is introduced.
-/
theorem nmod_mul_eq_preinvRoundIR (this : DenseUPolyZp) (a b : UInt64) :
    dense_upoly_zp_nmod_mul_ir this a b =
      let prod : UInt128 :=
        uint128_of_uint64 (a <<< this._norm) * uint128_of_uint64 b
      (preinvRoundIR
        (uint128_lo (prod >>> (64 : UInt128)))
        (uint128_lo prod)
        (this._p <<< this._norm) this._ninv) >>> this._norm := by
  simp [dense_upoly_zp_nmod_mul_ir, preinvRoundIR] <;>
    repeat' split <;> simp_all

end CLPoly.Impl.StrictWordArithmetic
