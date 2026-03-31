-- 手动修正版：验证 nmod_mul 翻译的 Lean 是否可编译

-- UInt128 占位（Phase 2 待实现）
structure UInt128 where
  hi : UInt64
  lo : UInt64
deriving Repr

def uint128_mul (a b : UInt64) : UInt128 :=
  -- 简化：实际需要正确的 128 位乘法
  ⟨0, a * b⟩

def uint128_hi (x : UInt128) : UInt64 := x.hi
def uint128_lo (x : UInt128) : UInt64 := x.lo

-- 翻译后的 nmod_mul（手动修正版）
partial def nmod_mul_ir (a b p ninv : UInt64) (norm : UInt64)
    (hp : p ≠ 0)
    : UInt64 :=
  let pn_1 : UInt64 := p <<< norm
  let a_shifted_1 : UInt64 := a <<< norm
  let prod_1 : UInt128 := uint128_mul a_shifted_1 b
  let hi_1 : UInt64 := uint128_hi prod_1
  let lo_1 : UInt64 := uint128_lo prod_1
  let qm_1 : UInt128 := uint128_mul hi_1 ninv
  let q1_1 : UInt64 := uint128_hi qm_1
  let q0_1 : UInt64 := uint128_lo qm_1
  let q0_2 : UInt64 := q0_1 + lo_1
  let q1_2 : UInt64 := q1_1 + (hi_1 + (if q0_2 < lo_1 then 1 else 0))
  let r_1 : UInt64 := lo_1 - ((q1_2 + 1) * pn_1)
  let r_2 : UInt64 := if r_1 > q0_2 then r_1 + pn_1 else r_1
  let r_3 : UInt64 := if r_2 >= pn_1 then r_2 - pn_1 else r_2
  r_3 >>> norm

-- 测试：可 #eval
#eval nmod_mul_ir 7 11 13 0 0 (by decide)
#eval nmod_mul_ir 100 200 17 0 0 (by decide)
