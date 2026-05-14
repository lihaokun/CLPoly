-- B2B Lean 端函数注册表
--
-- dispatch fn args → Except String Json
-- 找不到函数返回 .error "unknown fn: <name>"

import Lean.Data.Json
import B2B.Types
import CLPoly.Generated.Corpus
import CLPoly.Model

open Lean

namespace B2B

-- args 是 Json.arr，dispatch 已校验
def dispatch (fn : String) (args : Array Json) : Except String Json := do
  match fn with
  | "__make_zp" =>
    let val ← parseInt64 args[0]!
    let p   ← parseUInt64 args[1]!
    return encodeZp (Generated.__make_zp_ir val p)
  -- Zp 算术 wrapper（C++ 端是 operator 重载，B2B 起 ____xx_zp 名）
  | "__add_zp" =>
    let a ← parseZp args[0]!
    let b ← parseZp args[1]!
    return encodeZp (a + b)
  | "__sub_zp" =>
    let a ← parseZp args[0]!
    let b ← parseZp args[1]!
    return encodeZp (a - b)
  | "__mul_zp" =>
    let a ← parseZp args[0]!
    let b ← parseZp args[1]!
    return encodeZp (a * b)
  | "__neg_zp" =>
    let a ← parseZp args[0]!
    return encodeZp (-a)
  | "__inv_zp" =>
    let a ← parseZp args[0]!
    return encodeZp a.inv
  | "__pow_zp" =>
    let a ← parseZp args[0]!
    let i ← parseInt64 args[1]!
    return encodeZp (a ^ i)
  -- ---- ZZ 算术 ----
  | "__gcd_zz" =>
    let a ← parseZZ args[0]!
    let b ← parseZZ args[1]!
    return encodeZZ (Int.gcd a b)
  | "__lcm_zz" =>
    let a ← parseZZ args[0]!
    let b ← parseZZ args[1]!
    return encodeZZ (Int.lcm a b)
  | "__fdiv_q_zz" =>
    let a ← parseZZ args[0]!
    let b ← parseZZ args[1]!
    return encodeZZ (ZZ.fdiv_q 0 a b)
  | "__fdiv_r_zz" =>
    let a ← parseZZ args[0]!
    let b ← parseZZ args[1]!
    return encodeZZ (ZZ.fdiv_r 0 a b)
  | "__fdiv_ui_zz" =>
    let a ← parseZZ args[0]!
    let b ← parseUInt64 args[1]!
    -- Lean ZZ.fdiv_ui 返 ZZ；C++ 返 uint64_t。统一编码为 UInt64
    return encodeUInt64 (ZZ.fdiv_ui a b).toNat.toUInt64
  | "__sizeinbase_zz" =>
    let z    ← parseZZ args[0]!
    let base ← parseInt32 args[1]!
    return encodeUInt64 (ZZ.sizeinbase z base)
  | "__invert_zz" =>
    let a ← parseZZ args[0]!
    let m ← parseZZ args[1]!
    let (ok, inv) := ZZ.invertImpl a m
    return encodeBoolZZ ok inv
  -- ---- SparsePolyZZ ops ----
  | "__cont_zz_poly" =>
    let p ← parseSparsePolyZZ args[0]!
    return encodeZZ (SparsePolyZZ.contImpl p)
  | "__pp_zz_poly" =>
    let p ← parseSparsePolyZZ args[0]!
    return encodeSparsePolyZZ (SparsePolyZZ.ppImpl p)
  | "__derivative_zz_poly" =>
    let p ← parseSparsePolyZZ args[0]!
    return encodeSparsePolyZZ (SparsePolyZZ.derivativeImpl p)
  | "__normalization_zz_poly" =>
    let p ← parseSparsePolyZZ args[0]!
    return encodeSparsePolyZZ (SparsePolyZZ.normalization p)
  | "__polynomial_mod_zz" =>
    let p ← parseSparsePolyZZ args[0]!
    let prime ← parseUInt64 args[1]!
    return encodeSparsePolyZp (polynomial_mod p prime)
  -- ---- SparsePolyZp ops ----
  | "__add_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    return encodeSparsePolyZp (f + g)
  | "__sub_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    return encodeSparsePolyZp (f - g)
  | "__mul_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    return encodeSparsePolyZp (f * g)
  | "__divmod_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    let (q, r) := SparsePolyZp.divmod f g
    return encodePairSPZp q r
  | "__gcd_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    return encodeSparsePolyZp (SparsePolyZp.gcd f g)
  | "__gcd_eea_zp_poly" =>
    let f ← parseSparsePolyZp args[0]!
    let g ← parseSparsePolyZp args[1]!
    let (gcd, s, t) := polynomial_GCD_eea f g #[] #[]
    return encodeTripleSPZp gcd s t
  | "__derivative_zp_poly" =>
    let p ← parseSparsePolyZp args[0]!
    return encodeSparsePolyZp (SparsePolyZp.derivative p)
  | "__normalization_zp_poly" =>
    let p ← parseSparsePolyZp args[0]!
    return encodeSparsePolyZp (SparsePolyZp.normalization p)
  -- ---- A5: MvPoly ops ----
  | "__assign_mv_zz" =>
    let p ← parseMvPolyZZ args[0]!
    let v ← parseVariable args[1]!
    let c ← parseZZ args[2]!
    return encodeMvPolyZZ (assign p v c)
  | "__assign2_mv_zz" =>
    let p ← parseMvPolyZZ args[0]!
    let m ← parseVarMapZZ args[1]!
    return encodeMvPolyZZ (assign2 p m)
  | "__leadcoeff_mv_zz" =>
    let p ← parseMvPolyZZ args[0]!
    let v ← parseVariable args[1]!
    return encodeMvPolyZZ (leadcoeff p v)
  | "__poly_convert_spzp_to_spzz" =>
    let p ← parseSparsePolyZp args[0]!
    return encodeSparsePolyZZ (poly_convert p (#[] : SparsePolyZZ))
  | "__poly_convert3_spzz_to_mvzz" =>
    let p ← parseSparsePolyZZ args[0]!
    let v ← parseVariable args[1]!
    return encodeMvPolyZZ (poly_convert3 p (#[] : MvPolyZZ) v)
  -- ---- A6: Bezout ----
  | "__nat_extgcd" =>
    let a ← parseZZ args[0]!
    let b ← parseZZ args[1]!
    -- Nat.extGcd 接受非负 Nat；ZZ 输入需取 natAbs。但 C++ 端可处理负数（gmp_gcdext 支持）。
    -- 为简化：B2B 测例约定 a,b ≥ 0
    let (g, s, t) := Nat.extGcd a.toNat b.toNat
    return encodeZZTriple (g : Int) s t
  | _ => throw s!"unknown fn: {fn}"

end B2B
