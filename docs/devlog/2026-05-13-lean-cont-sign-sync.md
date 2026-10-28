# 2026-05-13: Lean 端 cont 同步 issue #17（删 sign-match）

## 做了什么

同步 PR #18 (master，已合并) 的 cont sign 修复到 Lean 端 `feature/formal-proofs` 分支。删除 `Model.lean` 中 3 处 `contImpl` / `contInt` 的"符号匹配 leading coeff"逻辑，使 Lean cont 与 C++ post-fix 一致（**始终非负**，对齐 Maple/SymPy/FLINT）。

修改：
- `proof/lean/CLPoly/Model.lean:1158` `MvPolyZZ.contImpl`：删 sign-match 行
- `proof/lean/CLPoly/Model.lean:1268` `SparsePolyZZ.contImpl`：删 sign-match 行
- `proof/lean/CLPoly/Model.lean:1626` `MvPolyZZ.contInt`：删 sign-match 行
- `proof/lean/CLPoly/Model.lean:1413+` `#eval` 期望值更新 + 4 个新 case

## 为什么做

PR #18 合并时 PR 描述里承诺：
> Lean Model.lean 端有自己的 `MvPolyZZ.contImpl` / `SparsePolyZZ.contImpl`（在 `feature/formal-proofs` 分支），同样问题。本 PR 后会在 feature/formal-proofs 分支独立同步修。

本日志兑现此承诺。

## 关键决策及其理由

### 决策 1: empty case 保留返 `0`（不对齐 C++ `return 1`）

`MvPolyZZ.contImpl(empty) = 0` 是 Lean 内部 **sentinel 值**，被 wrapper `MvPolyZZ.contMv` 用 `if c = 0 then #[]` 判 empty。若改成 1，contMv 会变成 `#[(#[], 1)]`（常数 1 polynomial），破坏 empty-detection。

承认事实：`contMv` 本身是 incomplete placeholder（Model.lean:1171-1174 注释明确"Phase F-impl-v2 续"），并未完全对齐 C++ multivariate cont。本 fix 仅同步 sign 约定，contMv 完整对齐留给后续 phase。

### 决策 2: 附带修复 `__upoly_primitive_upoly_ir` 的 sign 反转潜伏 bug

审核 doc 时发现 `Corpus.lean:5044` (`__upoly_primitive_upoly_ir`) 是 C++ `polynomial_factorize_univar.hh:714-722` 的翻译，caller **显式做 sign-flip 保 pp 首项 > 0**：
```lean
let c_1 := cont f
if (front.snd < 0) then c_2 = -c_1
```

修前 Lean 行为（trace 输入 `f = [(-x)]`）：
- contImpl 返 `-1`（sign-match），caller 翻成 `+1`，pp = f/+1 = `-x` ❌ pp 首项负

修后行为：
- contImpl 返 `+1`，caller 翻成 `-1`，pp = f/-1 = `+x` ✓ pp 首项正

→ 两层 sign-match 互相抵消的潜伏 bug 被本 fix **附带修复**（非主目标但顺手）。

### 决策 3: 范围限定 — 仅 sign 处理，不动 contMv 完整性 / pp 兜底

- 不动 `MvPolyZZ.contMv` 的 incomplete placeholder（Phase F-impl-v2 范围）
- 不动 `Corpus.lean:6888` `squarefreefactorize_lex_ir` 的 PR #16 sign-flip 兜底（独立机制，与 cont sign 正交，本身 robust）

## 遇到的问题与解决方式

### 问题 1: 编译错误 — Nat → Int 类型推断冲突

初版改写 `(f.foldl ... 0 : Int)`，Lean 把 `0` 推断为 Int 但 lambda 内 `acc : Nat`，导致 Application type mismatch。

解决：保留显式 `let c_nat : Nat := f.foldl ... 0; (c_nat : Int)` 两步式，类型推断明确。

### 问题 2: doc 审核漏 Corpus.lean 调用点

初版 doc §4 仅 audit Model.lean 内部，遗漏 Pass 8 翻译产物 (`Corpus.lean`) 的实际调用方。审核时 grep 整个 `proof/lean/` 发现 3 个 corpus 调用点（其中 1 个含 sign 反转潜伏 bug），doc 升级为 4a-4e 子段。

### 问题 3: cont 行为变化是否影响下游证明

担心改 cont sign 会破坏已有 proof。Grep 确认：除 Model.lean 自身 + Corpus.lean 3 个调用 + Math/MvGCD.lean 占位 spec，无其他 cont/pp 引用 → **无证明依赖**。lake build 验证：3071 jobs 全过。

## 涉及的文件

### Lean 改动
- `proof/lean/CLPoly/Model.lean:1157-1164` `MvPolyZZ.contImpl`（删 sign-match + 注释 + 用 let 解类型推断）
- `proof/lean/CLPoly/Model.lean:1267-1274` `SparsePolyZZ.contImpl`（同上）
- `proof/lean/CLPoly/Model.lean:1625-1632` `MvPolyZZ.contInt`（同上）
- `proof/lean/CLPoly/Model.lean:1413-1416` `#eval` 期望值修正（`-2` → `2`）
- `proof/lean/CLPoly/Model.lean:1418+` 新增 4 个 `#eval` 测例（`cont(-x)`, `cont(-6)`, `pp(-2x²-4)`, 含期望值）

### 新增文档
- `docs/fixes/2026-05-13-lean-cont-sign-sync.md` 修正方案（含半形式化论证 + 调用点审计 + 关键决策）
- `docs/devlog/2026-05-13-lean-cont-sign-sync.md` 本日志

## 度量

- 耗时：~2 小时
  - 切分支 + 调研定位（grep contImpl / Corpus.lean cont 调用）：~20 min
  - 修正方案 doc 第一版：~20 min
  - 审核轮：发现 `__upoly_primitive` 潜伏 bug + 5 处 doc 缺漏 → 修订：~30 min
  - 实施 + 调试（Nat/Int 类型错 → 修 → 全过）：~20 min
  - 验证 #eval 输出 + lake build 回归：~10 min
  - devlog + commit：~10 min
- 迭代：2 轮（type-mismatch 修一次，doc 审核轮一次）
- Lean 改动行数：~30 行（3 函数重写 + 4 个新 #eval）
- 对应 C++ 行数：~6 行（PR #18 的 3 个改动）
- 修正方案文档：~250 行（fix doc）
- 放弃的方案：
  - 改 empty case → 1 对齐 C++（会破坏 contMv 的 empty-detection sentinel，放弃）
  - 在 Lean 端补 pp sign-normalize 兜底（无必要，Corpus.lean 实际调用方已有兜底层或被本 fix 顺手修复）
  - 显式 `(... : Int)` 类型 ascription（编译失败，改用 let 显式 Nat 中间变量）

## 后续

- 不开 PR（feature/formal-proofs 不对 master），直接 push
- 等 Phase F-impl-v2 时一起处理 `MvPolyZZ.contMv` 的 incomplete placeholder（"返回剩余变量多项式"完整语义）
