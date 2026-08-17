-- B2B 类型库 round-trip 测试（Lean 端）
-- 运行：lake env lean --run B2B/TestTypes.lean

import B2B.Types
import CLPoly.Model

open Lean B2B

-- 简单测试框架
def assertOk {α : Type} [BEq α] [ToString α]
    (label : String) (expected : α) (actual : Except String α) : IO Bool := do
  match actual with
  | .ok v =>
    if v == expected then
      IO.println s!"  PASS: {label}"
      return true
    else
      IO.println s!"  FAIL: {label} — expected {expected}, got {v}"
      return false
  | .error e =>
    IO.println s!"  FAIL: {label} — error: {e}"
    return false

def runTest (name : String) (test : IO Bool) : IO Bool := do
  IO.println s!"--- {name} ---"
  test

def testInt32 : IO Bool := do
  let mut all_ok := true
  for v in [(0 : Int32), 1, -1, 2147483647, -2147483648] do
    let j := encodeInt32 v
    let parsed := parseInt32 j
    let ok ← assertOk s!"Int32 {v}" v parsed
    if !ok then all_ok := false
  return all_ok

def testInt64 : IO Bool := do
  let mut all_ok := true
  for v in [(0 : Int64), 1, -1, 9223372036854775807, -9223372036854775808] do
    let j := encodeInt64 v
    let parsed := parseInt64 j
    let ok ← assertOk s!"Int64 {v}" v parsed
    if !ok then all_ok := false
  return all_ok

def testUInt64 : IO Bool := do
  let mut all_ok := true
  for v in [(0 : UInt64), 1, 13, 18446744073709551615] do
    let j := encodeUInt64 v
    let parsed := parseUInt64 j
    let ok ← assertOk s!"UInt64 {v}" v parsed
    if !ok then all_ok := false
  return all_ok

def testZZ : IO Bool := do
  let mut all_ok := true
  let cases : List Int := [0, 1, -1, 123456789012345, 99999999999999999999999999]
  for v in cases do
    let j := encodeZZ v
    let parsed := parseZZ j
    let ok ← assertOk s!"ZZ {v}" v parsed
    if !ok then all_ok := false
  return all_ok

def testZp : IO Bool := do
  let mut all_ok := true
  let cases : List Zp := [
    Zp.ofInt 3 7,
    Zp.ofInt 0 13,
    Zp.ofInt 6 7
  ]
  for v in cases do
    let j := encodeZp v
    match parseZp j with
    | .ok p =>
      if p.val == v.val ∧ p.prime == v.prime then
        IO.println s!"  PASS: Zp({v.val}, {v.prime})"
      else
        IO.println s!"  FAIL: Zp expected ({v.val},{v.prime}), got ({p.val},{p.prime})"
        all_ok := false
    | .error e =>
      IO.println s!"  FAIL: Zp parse error: {e}"
      all_ok := false
  return all_ok

def testSparsePolyZp : IO Bool := do
  let mut all_ok := true
  let p1 : SparsePolyZp := #[]
  let p2 : SparsePolyZp := #[
    (UMonomial.mk 2, Zp.ofInt 3 7),
    (UMonomial.mk 0, Zp.ofInt 5 7)
  ]
  for v in [p1, p2] do
    let j := encodeSparsePolyZp v
    match parseSparsePolyZp j with
    | .ok p =>
      let mut eq := p.size == v.size
      for i in [:v.size] do
        let (m1, c1) := p[i]!
        let (m2, c2) := v[i]!
        if m1.deg != m2.deg ∨ c1.val != c2.val ∨ c1.prime != c2.prime then
          eq := false
      if eq then
        IO.println s!"  PASS: SparsePolyZp size={v.size}"
      else
        IO.println s!"  FAIL: SparsePolyZp roundtrip mismatch (size={v.size})"
        all_ok := false
    | .error e =>
      IO.println s!"  FAIL: SparsePolyZp parse error: {e}"
      all_ok := false
  return all_ok

def testFactorResultEncodings : IO Bool := do
  let factorZp : SparsePolyZp := #[
    (UMonomial.mk 1, Zp.ofInt 1 5),
    (UMonomial.mk 0, Zp.ofInt 4 5)
  ]
  let encodedZp := encodeFactorZpResult
    (Zp.ofInt 3 5, #[(factorZp, (2 : UInt64))])
  let expectedZp := Json.mkObj [
    ("type", "FactorZpResult"),
    ("val", Json.arr #[encodeZp (Zp.ofInt 3 5),
      Json.arr #[Json.arr #[encodeSparsePolyZp factorZp, encodeUInt64 2]]])
  ]
  let encodedZZ := encodeArraySPZZ
    #[#[(UMonomial.mk 1, (1 : Int)), (UMonomial.mk 0, (-2 : Int))]]
  let expectedZZ := Json.mkObj [
    ("type", "ArraySPZZ"),
    ("val", Json.arr #[encodeSparsePolyZZ
      #[(UMonomial.mk 1, (1 : Int)), (UMonomial.mk 0, (-2 : Int))]])
  ]
  let zpOk := encodedZp == expectedZp
  let zzOk := encodedZZ == expectedZZ
  if zpOk then IO.println "  PASS: FactorZpResult observable encoding"
  else IO.println "  FAIL: FactorZpResult observable encoding"
  if zzOk then IO.println "  PASS: ArraySPZZ observable encoding"
  else IO.println "  FAIL: ArraySPZZ observable encoding"
  return zpOk && zzOk

def testTypeMismatch : IO Bool := do
  let j := encodeInt64 7
  match parseInt32 j with
  | .ok _ =>
    IO.println "  FAIL: type mismatch should raise but didn't"
    return false
  | .error _ =>
    IO.println "  PASS: type mismatch raises"
    return true

def main : IO UInt32 := do
  IO.println "=== B2B Lean types round-trip ==="
  let r1 ← runTest "Int32" testInt32
  let r2 ← runTest "Int64" testInt64
  let r3 ← runTest "UInt64" testUInt64
  let r4 ← runTest "ZZ" testZZ
  let r5 ← runTest "Zp" testZp
  let r6 ← runTest "SparsePolyZp" testSparsePolyZp
  let r7 ← runTest "factor result encodings" testFactorResultEncodings
  let r8 ← runTest "type mismatch" testTypeMismatch
  let all := r1 && r2 && r3 && r4 && r5 && r6 && r7 && r8
  if all then
    IO.println "\n=== All round-trip tests PASSED ==="
    return 0
  else
    IO.println "\n=== FAILED ==="
    return 1
