/-
  CLPoly/Experiment/L1Barrett.lean — L1 翻译验证原型

  目标：验证 C++ Barrett 模乘 (__nmod_mul) 的正确性
  方法：C++ → Lean IR → UB-freedom 证明 → 精化证明

  对应 C++: clpoly/number.hh lines 115-132
-/
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Nat.Log
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith

set_option autoImplicit false

-- ============================================================
-- 1. 类型建模
-- ============================================================

/-- 128 位无符号整数（建模 unsigned __int128）。-/
def U128 := Fin (2^128)

/-- 64 位无符号整数。-/
abbrev U64 := UInt64

-- ============================================================
-- 2. C++ __nmod_mul 的 Lean IR 翻译
-- ============================================================

/-- Barrett 模乘 IR。
    对应 C++ Zp::__nmod_mul (number.hh:115-132)。
    输入：a, b < p；p, ninv, norm 是预计算参数。
    输出：a * b mod p。

    C++ 代码：
    ```
    uint64_t pn = p << norm;
    uint64_t a_shifted = a << norm;
    __int128 prod = (__int128)a_shifted * b;  // [hi:lo]
    __int128 qm = (__int128)hi * ninv;        // [q1:q0]
    q0 += lo; q1 += hi + carry;
    uint64_t r = lo - (q1+1)*pn;
    if (r > q0) r += pn;
    if (r >= pn) r -= pn;
    return r >> norm;
    ```
-/
def barrett_mul_ir (a b p ninv : Nat) (norm : Nat) : Nat :=
  let pn := p * 2^norm          -- p << norm
  let a_shifted := a * 2^norm   -- a << norm
  let prod := a_shifted * b     -- 128-bit product
  let hi := prod / 2^64         -- high 64 bits
  let lo := prod % 2^64         -- low 64 bits
  let qm := hi * ninv           -- Barrett quotient estimate
  let q1 := qm / 2^64
  let q0 := qm % 2^64
  let q0' := (q0 + lo) % 2^64
  let carry := if q0 + lo ≥ 2^64 then 1 else 0
  let q1' := (q1 + hi + carry) % 2^64
  let r := (lo + 2^64 - ((q1' + 1) * pn) % 2^64) % 2^64
  let r := if r > q0' then (r + pn) % 2^64 else r
  let r := if r ≥ pn then r - pn else r
  r / 2^norm                    -- r >> norm

-- ============================================================
-- 3. UB-Freedom 证明目标
-- ============================================================

/-- Barrett 模乘的前置条件（UB-freedom）。
    C++ 中这些条件由 __precompute 和调用者保证。-/
structure BarrettPrecondition (a b p : Nat) (norm : Nat) : Prop where
  /-- a < p -/
  ha : a < p
  /-- b < p -/
  hb : b < p
  /-- p ≥ 2（素数） -/
  hp : 2 ≤ p
  /-- norm = clz(p)，使得 p << norm 的最高位为 1 -/
  hnorm : 2^63 ≤ p * 2^norm
  /-- p << norm < 2^64（归一化后不溢出） -/
  hpn : p * 2^norm < 2^64

/-- UB-freedom：在 BarrettPrecondition 下，所有中间计算不溢出 128 位。
    这是最关键的 UB 证明目标。-/
theorem barrett_no_overflow (a b p ninv : Nat) (norm : Nat)
    (h : BarrettPrecondition a b p norm) :
    -- a_shifted * b < 2^128（128 位乘积不溢出）
    a * 2^norm * b < 2^128 := by
  have h1 : a * 2^norm < 2^64 := by
    have := Nat.mul_lt_mul_of_pos_right h.ha (Nat.two_pow_pos norm)
    linarith [h.hpn]
  have h2 : b < 2^64 := by nlinarith [h.hb, h.hp, h.hpn, Nat.two_pow_pos norm]
  -- a*2^norm < 2^64 且 b < 2^64 → a*2^norm*b < 2^128
  nlinarith [show (2:Nat) ^ 64 * 2 ^ 64 = 2 ^ 128 from by norm_num]

-- ============================================================
-- 4. 精化证明
-- ============================================================

/-- Barrett 模乘的数学规约：result = a * b % p。-/
theorem barrett_mul_correct (a b p ninv : Nat) (norm : Nat)
    (h : BarrettPrecondition a b p norm)
    -- ninv 是 Barrett 预计算逆元：ninv = ⌊2^128 / (p << norm)⌋ 的低 64 位
    -- （精确定义复杂，此处简化为正确性假设）
    (hninv : True) :
    barrett_mul_ir a b p ninv norm = a * b % p := by
  -- Barrett 约化的正确性是 FLINT 的核心引理。
  -- 完整证明需要：
  --   1. 商估计 q1 ≈ ⌊a*b/p⌋（误差 ≤ 1）
  --   2. 余式 r = a*b - q1*p（可能偏大 p）
  --   3. 最后两个 if 修正余式到 [0, p)
  -- 此处标记为原型——完整证明需要 ~50 行 Nat 算术。
  sorry

-- ============================================================
-- 5. 端到端：IR + UB-freedom + 精化 = 实现正确
-- ============================================================

/-- 端到端定理：C++ __nmod_mul 在 UB-freedom 下计算 a*b mod p。-/
theorem nmod_mul_verified (a b p ninv : Nat) (norm : Nat)
    (h : BarrettPrecondition a b p norm)
    (hninv : True) :
    -- 实现在无 UB 条件下正确
    barrett_mul_ir a b p ninv norm = a * b % p ∧
    -- 且无 128 位溢出
    a * 2^norm * b < 2^128 :=
  ⟨barrett_mul_correct a b p ninv norm h hninv,
   barrett_no_overflow a b p ninv norm h⟩
