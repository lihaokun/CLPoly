-- B2B Lean driver — stdin/stdout NDJSON dispatcher

import Lean.Data.Json
import B2B.Registry

open Lean B2B

-- 处理一行请求，返回响应 Json
def handleRequest (line : String) : Json :=
  let trimmed := line.trimAscii.toString
  if trimmed.isEmpty then Json.null else
  match Json.parse trimmed with
  | .error e =>
    Json.mkObj [("id", Json.str ""), ("ok", Json.bool false), ("err", Json.str s!"parse: {e}")]
  | .ok req =>
    let id := (req.getObjValD "id").getStr?.toOption.getD ""
    let fnRes := (req.getObjValD "fn").getStr?
    let argsRes := (req.getObjValD "args").getArr?
    match fnRes, argsRes with
    | .ok fn, .ok args =>
      match dispatch fn args with
      | .ok ret =>
        Json.mkObj [("id", Json.str id), ("ok", Json.bool true), ("ret", ret)]
      | .error e =>
        Json.mkObj [("id", Json.str id), ("ok", Json.bool false), ("err", Json.str e)]
    | .error e, _ =>
      Json.mkObj [("id", Json.str id), ("ok", Json.bool false), ("err", Json.str s!"fn: {e}")]
    | _, .error e =>
      Json.mkObj [("id", Json.str id), ("ok", Json.bool false), ("err", Json.str s!"args: {e}")]

partial def loop (stdin : IO.FS.Stream) (stdout : IO.FS.Stream) : IO Unit := do
  let line ← stdin.getLine
  if line.isEmpty then return ()
  let resp := handleRequest line
  if !resp.isNull then
    -- compress 输出单行 NDJSON（toString 会 pretty-print 跨行，破坏协议）
    stdout.putStrLn resp.compress
    stdout.flush
  loop stdin stdout

def main : IO Unit := do
  let stdin ← IO.getStdin
  let stdout ← IO.getStdout
  loop stdin stdout
