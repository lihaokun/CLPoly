-- B2B Lean 端函数注册表
--
-- dispatch fn args → Except String Json
-- 找不到函数返回 .error "unknown fn: <name>"

import Lean.Data.Json
import B2B.Types
import B2B.StrictRuntime
import CLPoly.Generated.StrictFactorZZ
import CLPoly.Generated.StrictHensel
import CLPoly.Generated.StrictVanHoeij
import CLPoly.Generated.StrictVanHoeijParameters
import CLPoly.Generated.StrictVanHoeijHeuristic
import CLPoly.Refinement.VanHoeijEntryExecution
import CLPoly.Refinement.VanHoeijRetryExecution
import CLPoly.Model

open Lean

namespace B2B

-- args 是 Json.arr，dispatch 已校验
unsafe def dispatch (fn : String) (args : Array Json) : Except String Json := do
  match fn with
  | "__lll_factorize" =>
    let f ← parseSparsePolyZZ args[0]!
    let polys ← args[1]!.getArr?
    let factors ← polys.mapM parseSparsePolyZp
    let prime ← args[2]!.getNat?
    if prime ≥ 2^64 || factors.size < 2 || factors.size ≥ 2^31 || get_deg f ≤ 0 then
      throw "invalid retry test input"
    let p := prime.toUInt64
    if hp : 2 ≤ p.toNat then
      let ctx := StrictRuntime.denseContext p
      let execution := Refinement.StrictRecombine.concreteVanHoeijRetryExecution ctx
        Refinement.StrictHensel.concreteDivmodTermination (StrictRuntime.rawMulProvider ctx)
      match Generated.StrictVanHoeijRetry.lllFactorize execution f factors p hp with
      | .error fault => throw s!"strict van Hoeij retry failed: {repr fault}"
      | .ok result =>
        return Json.arr ((result.map encodeSparsePolyZZ).toList.mergeSort
          fun a b => a.compress < b.compress).toArray
    else throw "prime must exceed one"
  | "__vanhoeij_try_recombine" =>
    let f ← parseSparsePolyZZ args[0]!
    let polys ← args[1]!.getArr?
    let lifted ← polys.mapM parseSparsePolyZZ
    let modulus ← parseZZ args[2]!
    if lifted.size < 2 || modulus ≤ 1 || get_deg f ≤ 0 then
      throw "invalid try-recombine test input"
    if hfits : lifted.size < 2^31 then
      match Generated.StrictVanHoeijEntry.tryRecombine
          (Refinement.StrictRecombine.concreteVanHoeijEntryExecution
            Refinement.StrictHensel.concreteDivmodTermination) f lifted modulus hfits with
      | .error fault => throw s!"strict try-recombine failed: {repr fault}"
      | .ok result =>
        let factors := ((result.factors.map encodeSparsePolyZZ).toList.mergeSort
          fun a b => a.compress < b.compress).toArray
        let expectedComplete ← args[3]!.getBool?
        let expectedJson ← args[4]!.getArr?
        let expected := (expectedJson.toList.mergeSort fun a b => a.compress < b.compress).toArray
        if result.complete != expectedComplete || factors != expected then
          throw "try-recombine differs from independently known factors/completion"
        return Json.arr #[toJson result.complete, Json.arr factors,
          toJson result.cldColumns.toInt, toJson result.retainedRows.toInt]
    else throw "factor count must fit Int32"
  | "__vanhoeij_choose_parameters" =>
    let count ← args[0]!.getNat?
    let height ← parseZZ args[1]!
    let modulus ← parseZZ args[2]!
    if count < 2 || count ≥ 2^31 || height < 0 || modulus ≤ 1 then
      throw "require 2<=count<2^31, height>=0, modulus>1"
    match Generated.StrictVanHoeijParameters.selectParameters count.toUInt32.toInt32 height modulus with
    | none => return Json.arr #[]
    | some output => return Json.arr #[encodeZZ output.divisor, encodeZZ output.modulus,
        encodeZZ output.residual]
  | "__vanhoeij_cld_width" =>
    let input ← args[0]!.getArr?
    let cld ← input.mapM parseSparsePolyZZ
    if cld.any (fun poly => !poly.isEmpty && get_deg poly ≥ 2147483647) then
      throw "degree+1 must fit Int32"
    return toJson (Generated.StrictVanHoeijParameters.cldWidth cld).toInt
  | "__vanhoeij_singleton_groups" =>
    let count ← args[0]!.getNat?
    if count > 4096 then throw "test count must be in 0..4096"
    match Generated.StrictVanHoeij.singletonGroups count.toUInt32.toInt32 with
    | .error fault => throw s!"strict singleton construction failed: {repr fault}"
    | .ok groups => return Json.arr (groups.map fun group => Json.arr (group.map fun i => toJson i.toInt))
  | "__vanhoeij_next_target" =>
    let target ← args[0]!.getNat?
    let width ← args[1]!.getNat?
    if target = 0 || width ≥ 2^31 || target > width then throw "require 0 < target <= width < 2^31"
    return toJson (Generated.StrictVanHoeij.nextCldTarget
      target.toUInt32.toInt32 width.toUInt32.toInt32).toInt
  | "__vanhoeij_clamp_precision" =>
    let bits ← args[0]!.getNat?
    if bits ≥ 2^64 then throw "invalid float bits"
    return toJson (Generated.StrictVanHoeijHeuristic.clampPrecision (Float.ofBits bits.toUInt64)).toInt
  | "__vanhoeij_full_precision" =>
    let prime ← args[0]!.getNat?
    if prime < 2 || prime ≥ 2^64 then throw "invalid prime"
    let target ← parseZZ args[1]!
    match Generated.StrictVanHoeijHeuristic.fullPrecision prime.toUInt64 target with
    | .error fault => throw s!"full precision failed: {repr fault}"
    | .ok exponent => return toJson exponent.toInt
  | "__heuristic_starting_precision" =>
    let f ← parseSparsePolyZZ args[0]!
    let count ← args[1]!.getNat?
    let prime ← args[2]!.getNat?
    if count ≥ 2^31 || prime < 2 || prime ≥ 2^64 then throw "invalid heuristic input"
    match Generated.StrictVanHoeijHeuristic.startingPrecision f (Int32.ofNat count) prime.toUInt64 with
    | .error fault => throw s!"heuristic failed: {repr fault}"
    | .ok result => return Json.arr #[toJson result.1.toInt, toJson result.2.toInt]
  | "__vanhoeij_next_precision" =>
    let precision ← args[0]!.getNat?
    let full ← args[1]!.getNat?
    if precision ≥ 2^31 || full ≥ 2^31 then throw "precision exceeds int range"
    match Generated.StrictVanHoeijRetry.nextPrecision (Int32.ofNat precision) (Int32.ofNat full) with
    | .ok next => return Json.mkObj [("value", toJson next.toInt)]
    | .error .arithmeticOverflow => return Json.mkObj [("overflow", toJson true)]
    | .error fault => throw s!"unexpected precision fault: {repr fault}"
  | "__vanhoeij_validate_partition" =>
    let f ← parseSparsePolyZZ args[0]!
    let liftedJson ← args[1]!.getArr?
    let lifted ← liftedJson.mapM parseSparsePolyZZ
    let modulus ← parseZZ args[2]!
    if modulus ≤ 1 then throw "modulus must exceed one"
    let groupsJson ← args[3]!.getArr?
    let groups ← groupsJson.mapM fun group => do
      let entries ← group.getArr?
      entries.mapM fun entry => do
        let index ← entry.getNat?
        if index ≥ lifted.size || index ≥ 2^31 then throw "invalid physical group index"
        pure index.toUInt32.toInt32
    let complete ← args[4]!.getBool?
    let initial : Generated.StrictVanHoeij.VanHoeijResult :=
      ⟨#[#[(⟨0⟩, 17)]], complete, 7, 9⟩
    match Generated.StrictVanHoeij.validatePartition f lifted modulus groups initial with
    | .error fault => throw s!"strict partition validation failed: {repr fault}"
    | .ok (accepted, answer) =>
      let factors := ((answer.factors.map encodeSparsePolyZZ).toList.mergeSort
        fun a b => a.compress < b.compress).toArray
      return Json.arr #[toJson accepted, toJson answer.complete, Json.arr factors,
        toJson answer.cldColumns.toInt, toJson answer.retainedRows.toInt]
  | "__upoly_norm_l1" =>
    let f ← parseSparsePolyZZ args[0]!
    return encodeZZ (Generated.StrictVanHoeij.normL1 f)
  | "__vanhoeij_cld_bound" =>
    let f ← parseSparsePolyZZ args[0]!
    if get_deg f ≤ 0 || get_deg f > 4096 then throw "test degree must be in 1..4096"
    match Generated.StrictVanHoeij.cldBound f with
    | .error fault => throw s!"strict CLD bound failed: {repr fault}"
    | .ok output => return encodeZZ output
  | "__vanhoeij_truncate_cld" =>
    let input ← args[0]!.getArr?
    let cld ← input.mapM parseSparsePolyZZ
    let divisor ← parseZZ args[1]!
    if divisor ≤ 0 then throw "positive divisor required"
    match Generated.StrictVanHoeij.truncateCld cld divisor with
    | .error fault => throw s!"strict truncation failed: {repr fault}"
    | .ok output => return Json.arr (output.map encodeSparsePolyZZ)
  | "__lll_reduce_with_norms" =>
    let input ← args[0]!.getArr?
    let matrix ← input.mapM fun row => do
      let entries ← row.getArr?
      entries.mapM parseZZ
    let bound ← parseZZ args[1]!
    match Generated.StrictVanHoeij.lllReduceWithNorms
        Refinement.StrictRecombine.concreteLLLExecution matrix
        StrictRuntime.erasedValue bound with
    | .error fault => throw s!"strict LLL execution failed: {repr fault}"
    | .ok (reduced, transform, _, norms) =>
      let encodeMatrix := fun (m : Array (Array ZZ)) =>
        Json.arr (m.map fun row => Json.arr (row.map encodeZZ))
      return Json.arr #[encodeMatrix reduced, encodeMatrix transform,
        Json.arr (norms.map fun norm =>
          Json.arr #[Json.str (toString (QQ.num norm)), Json.str (toString (QQ.den norm))])]
  | "__lll_short_span_from_norms" =>
    let input ← args[0]!.getArr?
    let norms ← input.mapM fun pair => do
      let parts ← pair.getArr?
      let numerator ← parseZZ parts[0]!
      let denominator ← parseZZ parts[1]!
      if denominator = 0 then throw "zero denominator"
      pure ((numerator : QQ) / (denominator : QQ))
    let bound ← parseZZ args[1]!
    return Json.arr ((Generated.StrictVanHoeij.shortSpanFromNorms norms bound).map
      fun row => encodeUInt64 row.toUInt64)
  | "__build_cld_matrix" =>
    let input ← args[0]!.getArr?
    let matrix ← input.mapM fun row => do
      let entries ← row.getArr?
      entries.mapM parseZZ
    let polys ← args[1]!.getArr?
    let cld ← polys.mapM parseSparsePolyZZ
    let current ← parseUInt64 args[2]!
    let target ← parseUInt64 args[3]!
    let modulus ← parseZZ args[4]!
    match Generated.StrictVanHoeij.buildCldMatrix matrix cld current.toNat
        target.toNat modulus with
    | .error fault => throw s!"strict CLD execution failed: {repr fault}"
    | .ok (output, added) =>
      return Json.arr #[Json.arr (output.map fun row => Json.arr (row.map encodeZZ)),
        encodeUInt64 added.toUInt64]
  | "__factor_Zp" =>
    let f ← parseSparsePolyZp args[0]!
    match StrictRuntime.factorZpRuntime f with
    | .ok result => return encodeFactorZpResult result
    | .error fault => throw s!"strict FactorZp execution failed: {repr fault}"
  | "__factor_squarefree_primitive_ZZ" =>
    let f ← parseSparsePolyZZ args[0]!
    let useLargePrime ← parseBool args[1]!
    match StrictRuntime.factorZZRuntime useLargePrime f with
    | .ok result =>
      if args.size > 2 then
        let expectedJson ← args[2]!.getArr?
        let expected := (expectedJson.toList.mergeSort fun a b => a.compress < b.compress).toArray
        let actual := ((result.map encodeSparsePolyZZ).toList.mergeSort
          fun a b => a.compress < b.compress).toArray
        if actual != expected then throw "outer FactorZZ differs from independently known factors"
      return encodeArraySPZZ result
    | .error fault => throw s!"strict FactorZZ execution failed: {repr fault}"
  | "__hensel_factor_count_fits" =>
    let factorCount ← parseUInt64 args[0]!
    return encodeBool
      (Generated.StrictHensel.__hensel_factor_count_fits_ir
        factorCount.toNat)
  | "__make_zp" =>
    let val ← parseInt64 args[0]!
    let p   ← parseUInt64 args[1]!
    return encodeZp (Zp.ofInt val.toInt p)
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
    let output := polynomial_GCD_eea f g #[] #[]
    return encodeTripleSPZp output.1 output.2.1 output.2.2
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
    let output := Nat.extGcd a.toNat b.toNat
    return encodeZZTriple (output.1 : Int) output.2.1 output.2.2
  | _ => throw s!"unknown fn: {fn}"

end B2B
