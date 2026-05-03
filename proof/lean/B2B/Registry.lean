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
  | _ => throw s!"unknown fn: {fn}"

end B2B
