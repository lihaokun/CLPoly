-- B2B Lean 端类型库：JSON ↔ Lean 类型双向转换
--
-- 协议：每个值在 JSON 中表示为 {"type": "<类型名>", "val": <载荷>}
--   parseX 收 Json 返回 Except String T
--   encodeX 收 T 返回 Json

import Lean.Data.Json
import CLPoly.Model

open Lean

namespace B2B

-- ---- 内部工具 ----

-- 校验 type 字段并返回 val 子节点
def checkType (j : Json) (expected : String) : Except String Json := do
  let ty ← (← j.getObjVal? "type").getStr?
  if ty != expected then
    throw s!"type mismatch: expected {expected}, got {ty}"
  j.getObjVal? "val"

-- 接受整数 JSON 或字符串，返回 Int
-- JsonNumber 内部 = mantissa * 10^-exponent；整数 case exponent=0
def parseAnyInt (j : Json) : Except String Int := do
  match j with
  | .num n =>
    if n.exponent == 0 then return n.mantissa
    else throw s!"parseAnyInt: non-integer JSON number (exponent={n.exponent})"
  | .str s =>
    match s.toInt? with
    | some i => return i
    | none => throw s!"parseAnyInt: not an integer string: {s}"
  | _ => throw s!"parseAnyInt: expected number/string, got {j.compress}"

-- ---- 标量 ----

def parseInt32 (j : Json) : Except String Int32 := do
  let v ← checkType j "Int32"
  let i ← parseAnyInt v
  return i.toInt32

def parseInt64 (j : Json) : Except String Int64 := do
  let v ← checkType j "Int64"
  let i ← parseAnyInt v
  return i.toInt64

def parseUInt64 (j : Json) : Except String UInt64 := do
  let v ← checkType j "UInt64"
  let i ← parseAnyInt v
  if i < 0 then throw s!"UInt64: negative value {i}"
  return i.toNat.toUInt64

def encodeInt32 (v : Int32) : Json := Json.mkObj [("type", "Int32"), ("val", v.toInt)]
def encodeInt64 (v : Int64) : Json := Json.mkObj [("type", "Int64"), ("val", v.toInt)]

-- UInt64 一律用 string 编码避免 IEEE754 精度（>2^53）
def encodeUInt64 (v : UInt64) : Json :=
  Json.mkObj [("type", "UInt64"), ("val", Json.str (toString v.toNat))]

def parseBool (j : Json) : Except String Bool := do
  let v ← checkType j "Bool"
  v.getBool?

def encodeBool (value : Bool) : Json :=
  Json.mkObj [("type", "Bool"), ("val", Json.bool value)]

-- ---- ZZ (= Int) ----

def parseZZ (j : Json) : Except String Int := do
  let v ← checkType j "ZZ"
  parseAnyInt v

def encodeZZ (z : Int) : Json :=
  Json.mkObj [("type", "ZZ"), ("val", Json.str (toString z))]

-- ---- Zp ----

def parseZp (j : Json) : Except String Zp := do
  let v ← checkType j "Zp"
  let arr ← v.getArr?
  if arr.size != 2 then throw s!"Zp val: expected [val, prime], got size {arr.size}"
  let valI ← parseAnyInt arr[0]!
  let primeI ← parseAnyInt arr[1]!
  if primeI < 0 then throw s!"Zp prime: negative {primeI}"
  return Zp.ofInt valI primeI.toNat.toUInt64

def encodeZp (z : Zp) : Json :=
  -- val/prime 都用 UInt64 → 字符串编码（与 C++ 端一致）
  Json.mkObj [
    ("type", "Zp"),
    ("val", Json.arr #[
      Json.str (toString z.val.toNat),
      Json.str (toString z.prime.toNat)
    ])
  ]

-- ---- BoolZZ (用于 ZZ::invert 结果) ----

def encodeBoolZZ (ok : Bool) (z : Int) : Json :=
  Json.mkObj [
    ("type", "BoolZZ"),
    ("val", Json.arr #[Json.bool ok, Json.str (toString z)])
  ]

-- ---- SparsePolyZp ----

def parseSparsePolyZp (j : Json) : Except String SparsePolyZp := do
  let v ← checkType j "SparsePolyZp"
  let arr ← v.getArr?
  let mut result : SparsePolyZp := #[]
  for term in arr do
    let termArr ← term.getArr?
    if termArr.size != 2 then throw s!"SparsePolyZp term: expected [deg, [val,prime]]"
    let degI ← parseAnyInt termArr[0]!
    if degI < 0 then throw s!"SparsePolyZp deg negative: {degI}"
    let zpArr ← termArr[1]!.getArr?
    if zpArr.size != 2 then throw "SparsePolyZp term zp: expected [val, prime]"
    let valI ← parseAnyInt zpArr[0]!
    let primeI ← parseAnyInt zpArr[1]!
    let zp := Zp.ofInt valI primeI.toNat.toUInt64
    result := result.push (UMonomial.mk degI.toNat, zp)
  return result

def encodeSparsePolyZp (p : SparsePolyZp) : Json :=
  let terms := p.map fun (m, c) =>
    Json.arr #[
      Json.num (JsonNumber.fromNat m.deg),
      Json.arr #[
        Json.str (toString c.val.toNat),
        Json.str (toString c.prime.toNat)
      ]
    ]
  Json.mkObj [("type", "SparsePolyZp"), ("val", Json.arr terms)]

-- ---- SparsePolyZZ ----

def parseSparsePolyZZ (j : Json) : Except String SparsePolyZZ := do
  let v ← checkType j "SparsePolyZZ"
  let arr ← v.getArr?
  let mut result : SparsePolyZZ := #[]
  for term in arr do
    let termArr ← term.getArr?
    if termArr.size != 2 then throw s!"SparsePolyZZ term: expected [deg, zz_decimal_str]"
    let degI ← parseAnyInt termArr[0]!
    if degI < 0 then throw s!"SparsePolyZZ deg negative: {degI}"
    let coef ← parseAnyInt termArr[1]!
    result := result.push (UMonomial.mk degI.toNat, coef)
  return result

def encodeSparsePolyZZ (p : SparsePolyZZ) : Json :=
  let terms := p.map fun (m, c) =>
    Json.arr #[
      Json.num (JsonNumber.fromNat m.deg),
      Json.str (toString c)
    ]
  Json.mkObj [("type", "SparsePolyZZ"), ("val", Json.arr terms)]

-- ---- 复合返回类型 ----

def encodePairSPZp (q r : SparsePolyZp) : Json :=
  Json.mkObj [("type", "PairSPZp"),
              ("val", Json.arr #[encodeSparsePolyZp q, encodeSparsePolyZp r])]

def encodeTripleSPZp (g s t : SparsePolyZp) : Json :=
  Json.mkObj [("type", "TripleSPZp"),
              ("val", Json.arr #[encodeSparsePolyZp g,
                                 encodeSparsePolyZp s,
                                 encodeSparsePolyZp t])]

/-- Observable return encoding shared by the C++ and strict Lean
`__factor_Zp` entries. -/
def insertCanonicalJson (item : Json) (items : List Json) : List Json :=
  match items with
  | [] => [item]
  | head :: tail =>
    if item.compress < head.compress then item :: head :: tail
    else head :: insertCanonicalJson item tail

def encodeFactorZpResult
    (result : Zp × Array (SparsePolyZp × UInt64)) : Json :=
  let factors := result.2.toList.map (fun item =>
      Json.arr #[encodeSparsePolyZp item.1, encodeUInt64 item.2])
    |>.foldr insertCanonicalJson [] |>.toArray
  Json.mkObj [
    ("type", "FactorZpResult"),
    ("val", Json.arr #[encodeZp result.1, Json.arr factors])
  ]

/-- Observable return encoding shared by the C++ and strict Lean
`__factor_squarefree_primitive_ZZ` entries. -/
def encodeArraySPZZ (factors : Array SparsePolyZZ) : Json :=
  Json.mkObj [
    ("type", "ArraySPZZ"),
    ("val", Json.arr (factors.map encodeSparsePolyZZ))
  ]

-- ---- A5: Variable ----

def parseVariable (j : Json) : Except String Variable := do
  let v ← checkType j "Variable"
  let i ← parseAnyInt v
  if i < 0 then throw s!"Variable: negative serial {i}"
  return i.toNat.toUInt64

def encodeVariable (v : Variable) : Json :=
  Json.mkObj [("type", "Variable"), ("val", Json.str (toString v.toNat))]

-- ---- A5: Monomial (内部 parse/serialize 不带 type 包装，用于嵌套) ----

def _parseMonomialInner (v : Json) : Except String Monomial := do
  let arr ← v.getArr?
  let mut result : Monomial := #[]
  for entry in arr do
    let entryArr ← entry.getArr?
    if entryArr.size != 2 then throw "Monomial entry: expected [var_serial, exp]"
    let serial ← parseAnyInt entryArr[0]!
    let exp ← parseAnyInt entryArr[1]!
    if serial < 0 then throw s!"Monomial var serial negative: {serial}"
    result := result.push (serial.toNat.toUInt64, exp.toInt64)
  return Monomial.normalize result

def _encodeMonomialInner (m : Monomial) : Json :=
  Json.arr (m.map fun (v, e) =>
    Json.arr #[Json.str (toString v.toNat), Json.num (JsonNumber.fromInt e.toInt)])

def parseMonomial (j : Json) : Except String Monomial := do
  let v ← checkType j "Monomial"
  _parseMonomialInner v

def encodeMonomial (m : Monomial) : Json :=
  Json.mkObj [("type", "Monomial"), ("val", _encodeMonomialInner m)]

-- ---- A5: MvPolyZZ / MvPolyZp ----

def parseMvPolyZZ (j : Json) : Except String MvPolyZZ := do
  let v ← checkType j "MvPolyZZ"
  let arr ← v.getArr?
  let mut result : MvPolyZZ := #[]
  for term in arr do
    let termArr ← term.getArr?
    if termArr.size != 2 then throw "MvPolyZZ term: expected [mono, zz_str]"
    let m ← _parseMonomialInner termArr[0]!
    let coef ← parseAnyInt termArr[1]!
    result := result.push (m, coef)
  return MvPolyZZ.normalization result

def encodeMvPolyZZ (p : MvPolyZZ) : Json :=
  let terms := p.map fun (m, c) =>
    Json.arr #[_encodeMonomialInner m, Json.str (toString c)]
  Json.mkObj [("type", "MvPolyZZ"), ("val", Json.arr terms)]

def parseMvPolyZp (j : Json) : Except String MvPolyZp := do
  let v ← checkType j "MvPolyZp"
  let arr ← v.getArr?
  let mut result : MvPolyZp := #[]
  for term in arr do
    let termArr ← term.getArr?
    if termArr.size != 2 then throw "MvPolyZp term: expected [mono, [val,prime]]"
    let m ← _parseMonomialInner termArr[0]!
    let zpArr ← termArr[1]!.getArr?
    if zpArr.size != 2 then throw "MvPolyZp term zp: expected [val,prime]"
    let valI ← parseAnyInt zpArr[0]!
    let primeI ← parseAnyInt zpArr[1]!
    let zp := Zp.ofInt valI primeI.toNat.toUInt64
    result := result.push (m, zp)
  return MvPolyZp.normalization result

def encodeMvPolyZp (p : MvPolyZp) : Json :=
  let terms := p.map fun (m, c) =>
    Json.arr #[_encodeMonomialInner m,
               Json.arr #[Json.str (toString c.val.toNat),
                          Json.str (toString c.prime.toNat)]]
  Json.mkObj [("type", "MvPolyZp"), ("val", Json.arr terms)]

-- ---- A5: VarMap (Variable → ZZ)，C++ 端是 std::map，Lean 端用 StdMap ----

def parseVarMapZZ (j : Json) : Except String (StdMap Variable Int) := do
  let v ← checkType j "VarMap"
  let arr ← v.getArr?
  let mut result : StdMap Variable Int := StdMap.empty
  for entry in arr do
    let entryArr ← entry.getArr?
    if entryArr.size != 2 then throw "VarMap entry: expected [var_serial, zz_str]"
    let serial ← parseAnyInt entryArr[0]!
    let val ← parseAnyInt entryArr[1]!
    if serial < 0 then throw s!"VarMap var serial negative: {serial}"
    result := result.insert serial.toNat.toUInt64 val
  return result

-- ---- A6: ZZTriple ----

def encodeZZTriple (g s t : Int) : Json :=
  Json.mkObj [("type", "ZZTriple"),
              ("val", Json.arr #[Json.str (toString g),
                                 Json.str (toString s),
                                 Json.str (toString t)])]

end B2B
