# 修正方案：Lean 端 cont 同步 issue #17（删 sign-match）

> **状态**：草稿，待用户审核
> 对应 workflow.md §5.1 算法内部错误修正流程
> 上游：[issue #17](https://github.com/lihaokun/CLPoly/issues/17) + [PR #18](https://github.com/lihaokun/CLPoly/pull/18)（C++ 端已合并）
> 分支：`feature/formal-proofs`

---

## 第一部分：背景与定位

### 上游变更
PR #18 把 C++ `cont()` 统一为 **always non-negative**（Maple/SymPy/FLINT 约定）：
- `clpoly/polynomial_gcd.hh:1067` univariate ZZ：`init c = 0` + `gcd(0, c_0) = |c_0|`
- `clpoly/upolynomial.hh:220` template Tc：同样模式
- `clpoly/polynomial_gcd.hh:631+` multivariate：return 前加 sign-normalize

### Lean 端现状
`proof/lean/CLPoly/Model.lean` 有三处 cont 实现，均带 **sign-match leading coeff** 逻辑：

| 位置 | 函数 | sign-match 代码 |
|------|------|----------------|
| L1158-1164 | `MvPolyZZ.contImpl` | `if f[0]!.snd < 0 then -c_int else c_int` |
| L1268-1274 | `SparsePolyZZ.contImpl` | `if f[0]!.snd < 0 then -c_int else c_int` |
| L1626-1632 | `MvPolyZZ.contInt` | `if f[0].snd < 0 then -g else g` |

### 与 C++ 上游的偏离
| 情形 | Lean 现状 | C++ post-fix | Maple/SymPy/FLINT |
|------|-----------|--------------|------------------|
| `cont(-x)` | -1 (sign) | 1 | 1 |
| `cont(-2x²-4)` | -2 (sign) | 2 | 2 |
| `cont(2x+4)` | 2 | 2 | 2 |
| `cont(0)` (empty) | 0 | 1 (univariate) | varies |

**Lean 端 cont 行为偏离 C++ post-fix 与三家工业系统**。本文档同步删除 sign-match，使 Lean 端 cont 始终非负。

### 修复定位升级：附带修复 Corpus.lean 中 `__upoly_primitive` 的 sign 反转 bug

C++ `polynomial_factorize_univar.hh:714-722` 的 `__upoly_primitive`：
```cpp
ZZ c = cont(f);
if (f.front().second < ZZ(0)) c = -c;  // 确保 lc > 0
for (auto& term : f) term.second /= c;
```
caller **显式翻 sign 以保 pp 首项正**。Pass 8 翻译到 `Corpus.lean:5033-5050`：
```lean
let c_1 : ZZ := (cont f)
if (front.snd < 0) then let c_2 := (- c_1)
else c_1
```

**当前 Lean 行为分析**（输入 `f = -x`）：
- `contImpl(-x)` 返 `-1`（sign-match）
- caller 见 leading<0 → `c_2 = -(-1) = +1`
- pp = f / +1 = `-x` → **pp 首项负** ❌（与 C++ 意图相反）

**修复后行为**：
- `contImpl(-x)` 返 `+1`
- caller 翻 → `c_2 = -1`
- pp = f / -1 = `+x` → **pp 首项正** ✓（对齐 C++ post-fix）

→ 本 fix 不仅同步 C++，还**附带修复 `__upoly_primitive` 在 Lean 端的 sign 反转**（两层 sign-match 互相抵消的潜伏 bug）。

---

## 第二部分：修复方案

### 改动 1: `Model.lean:1158-1164` MvPolyZZ.contImpl

**修前**：
```lean
def MvPolyZZ.contImpl (f : MvPolyZZ) : ZZ :=
  if f.isEmpty then 0
  else
    let c_nat := f.foldl (fun (acc : Nat) (term : Monomial × Int) =>
      Nat.gcd acc term.snd.natAbs) 0
    let c_int : Int := c_nat
    if f[0]!.snd < 0 then -c_int else c_int    -- ← 删除
```

**修后**：
```lean
-- issue #17 同步（PR #18 C++ 端）：cont 始终非负
def MvPolyZZ.contImpl (f : MvPolyZZ) : ZZ :=
  if f.isEmpty then 0
  else
    f.foldl (fun (acc : Nat) (term : Monomial × Int) =>
      Nat.gcd acc term.snd.natAbs) 0
```

注释 L1157 `-- MvPolyZZ: cont = signed gcd of all int coeffs; pp = f / cont`
改为：`-- MvPolyZZ: cont = non-negative gcd of all int coeffs (issue #17); pp = f / cont`

### 改动 2: `Model.lean:1268-1274` SparsePolyZZ.contImpl

完全同样的改动模式。注释 L1267 `-- cont = gcd 所有系数（符号匹配 leading coeff）` 改为 `-- cont = gcd 所有系数（始终非负，issue #17 / C++ PR #18 对齐）`。

### 改动 3: `Model.lean:1626-1632` MvPolyZZ.contInt

**修前**：
```lean
def MvPolyZZ.contInt (f : MvPolyZZ) : Int :=
  if h : 0 < f.size then
    let nat_gcd := f.foldl (fun (acc : Nat) (t : Monomial × Int) =>
      Nat.gcd acc t.snd.natAbs) 0
    let g : Int := (nat_gcd : Int)
    if f[0].snd < 0 then -g else g   -- ← 删除
  else 0
```

**修后**：
```lean
-- issue #17 同步：cont 始终非负
def MvPolyZZ.contInt (f : MvPolyZZ) : Int :=
  if 0 < f.size then
    (f.foldl (fun (acc : Nat) (t : Monomial × Int) =>
      Nat.gcd acc t.snd.natAbs) 0 : Int)
  else 0
```

注释 L1625 `-- 整数 cont（所有系数的 gcd 绝对值；首项符号决定整体符号）`
改为：`-- 整数 cont（所有系数 gcd 绝对值，始终非负，issue #17）`

### 改动 4: `Model.lean:1413-1416` #eval 期望值

**修前**：
```lean
-- cont(-2x² - 4) = -2 (sign matches leading coeff)
#eval SparsePolyZZ.contImpl
  (#[(⟨2⟩, (-2 : Int)), (⟨0⟩, (-4 : Int))] : SparsePolyZZ)
-- 期望: -2
```

**修后**：
```lean
-- issue #17 / PR #18：cont 始终非负
#eval SparsePolyZZ.contImpl
  (#[(⟨2⟩, (-2 : Int)), (⟨0⟩, (-4 : Int))] : SparsePolyZZ)
-- 期望: 2
```

---

## 第三部分：empty case 不改（保留返 `0`）

### 偏离决策
C++ univariate `cont(empty) = 1`；Lean `contImpl(empty) = 0`。本 PR **不改 empty case**。

### 理由：Lean 内部 ABI 选择，不涉及 spec 对齐

`contImpl` 返 `0` 是被 wrapper 函数依赖的 **sentinel 值**（empty-detection），不是 spec：

1. **`MvPolyZZ.contMv`** (L1175-1178)：用 `c == 0` 判 empty case 返 `#[]`：
   ```lean
   def MvPolyZZ.contMv (f : MvPolyZZ) : MvPolyZZ :=
     let c : ZZ := MvPolyZZ.contImpl f
     if c = 0 then #[]
     else #[(#[], c)]
   ```
   若改 contImpl(empty) → 1，contMv(empty) 会变成 `#[(#[], 1)]`（常数 1 polynomial）— **打破 empty-detection**。

2. **`SparsePolyZZ.ppImpl`** (L1278-1280) graceful fallback `if c = 0 then f`：空入空出，无副作用。

3. **`SparsePolyZZ.gcd`** (L1367-1386) 显式 `if F.isEmpty then G; else if G.isEmpty then F`：不进入 contImpl 调用。

### 关于 contMv 的"对齐 C++"

需注意：`MvPolyZZ.contMv` 本身是 **incomplete placeholder**（Model.lean:1171-1174 注释明确："对真正多变量输入是不完整占位，Phase F-impl-v2 续"）— **并未完全对齐 C++ multivariate cont** 的"返回剩余变量多项式"语义。本 fix 不涉及 contMv 改动，仅保持现有 sentinel 约定不破。

→ empty case = 0 是 **Lean 内部 ABI 选择**（私有契约），与 C++ univariate `return 1` 的偏离不影响公开语义。

---

## 第四部分：下游影响 audit

### 4a. Model.lean 内部调用点

| 调用点 | 文件:行 | 用 cont 的什么? | sign 变化影响? |
|--------|--------|----------------|---------------|
| `MvPolyZZ.contMv` | Model.lean:1176 | `c = 0` 判 empty + 包成常数 poly | ✅ 透明（仅看零/非零）|
| `MvPolyZZ.ppImpl` | Model.lean:1167 | `f.map (snd / c)` | ⚠️ pp 首项可能变负（issue #17 目标行为）|
| `SparsePolyZZ.ppImpl` | Model.lean:1278 | 同上 | ⚠️ 同上 |
| `SparsePolyZZ.gcd` | Model.lean:1371-1374 | `Int.gcd cF cG` 取绝对值 | ✅ 透明 |
| `MvPolyZZ.polynomialGCD` | Model.lean:1676,1681,1695-1697 | `Nat.gcd c.natAbs cInt.natAbs` | ✅ 透明 |

### 4b. Corpus.lean (Pass 8 翻译) 实际调用点

跨 `find proof/lean/ -name "*.lean" -not -path "*/.lake/*"` 找到 3 个真正 cont/pp 调用：

| 位置 | 函数 | 调用形式 | 行为分析 |
|------|------|---------|---------|
| Corpus.lean:5044 | `__upoly_primitive_upoly_ir` | `cont f` + caller flip | ⚠️ **现有 bug → 本 fix 修复**（见 §1）|
| Corpus.lean:5489 | `_lambda___wang_core_lex_1_ir` | `pp h` + caller sign-normalize loop | ✅ robust（修前修后等价输出 — 兜底层吸收 sign 差异）|
| Corpus.lean:6888 | `squarefreefactorize_lex_ir` | `cont F` + PR #16 sign-flip 层 | ✅ robust（同上）|

**机制**：L5489 / L6888 的 caller 形如：
```lean
let pp_or_quot := ... (cont/pp F)
if (pp_or_quot.front.snd < 0) then 翻 sign else 不翻
```
此 if-flip 把 "pp 首项符号" 强制为 +。无论 cont 返 +|gcd| 还是 -|gcd|，最终都被这层 normalize 兜底 → robust。

L5044 不同：caller flip 后 pass 给 division loop 用作除数。当前两层 sign-match 互相抵消，导致 pp 首项实际为负（**bug**）。本 fix 解开这个抵消。

### 4c. Math/MvGCD.lean 数学 spec

`proof/lean/CLPoly/Math/MvGCD.lean:101-104` 有 `noncomputable def cont_spec`，用 `fromMathlib (toMathlib f)` 占位（非引用 contImpl）。**独立，无依赖**。✅

### 4d. Proof 文件依赖
跨 `find proof/lean/ -name "*.lean" -not -path "*/.lake/*"`：**除 Corpus.lean 3 个调用 + Math/MvGCD.lean 占位 + Model.lean 自身**，无任何其他 cont/pp 引用。

→ **无证明会因 sign 变化而破裂**；唯一行为变化是 L5044 的 bug 修复。

### 4e. pp 行为变化（issue #17 目标）
`pp(-2x²-4)`：
- 修前：`cont = -2, pp = (-2x²-4) / -2 = x²+2`（首项正）
- 修后：`cont = 2, pp = (-2x²-4) / 2 = -x²-2`（首项负，符合 Maple/SymPy/FLINT）

下游若依赖"pp 首项正"：见 §4b — Corpus.lean 3 个调用点已 audit，全部 robust 或被本 fix 顺手修复。

---

## 第五部分：测试计划

### 编译验证
- `lake build` 全过（3071 jobs）
- 若有 corpus 调用方依赖 pp 旧行为，可能产生 `#eval` 输出变化或 proof 失败 — 现场处理

### `#eval` 输出对照

#### 已有 #eval（须更新期望）
- L1405 `cont(2x²+4x+6)` → `2`（修前已正，无变化）
- L1414 `cont(-2x²-4)` → `2`（修前 `-2`，**修复目标**）

#### 应新增 #eval（覆盖关键变化）
- `cont(-x)` → `1`（核心 single-term bug）
- `cont(-6)` (常数 poly) → `6`（与 issue #17 §1 C++ 测例对齐）
- `pp(-2x²-4)` → `-x²-2`（验证 pp 首项变负，对齐 Maple/SymPy/FLINT）
- **`__upoly_primitive_upoly_ir [(-x)]`** → `(-1, +x)`（关键回归 — 验证 pp 首项 > 0，确认 §1 的潜伏 bug 已修）

### C++ 端不动
C++ 端 PR #18 已合并，测试已 cover（`test/test_cont_sign_consistency.cc` 16 用例 / 3 段全过）。Lean 端是独立同步。

---

## 第六部分：执行步骤（待批准）

1. **[等批准]** 用户确认本方案
2. 应用 4 处改动（Model.lean）
3. `lake build` 验证编译
4. `#eval` 输出对照 L1405 / L1414 / 新增几个 case
5. 若 corpus 内有 pp 依赖破裂：分析 → 修兜底（参 C++ PR #16 模式）
6. devlog + commit + 推 feature/formal-proofs

---

## 第七部分：风险与未决问题

1. **Q**: corpus 中 pp 调用方是否依赖"pp 首项正"？
   **A**: ✅ 已 audit（§4b）。Corpus.lean 3 个调用点：2 个有 sign-normalize 兜底层（robust），1 个是潜伏 bug 本 fix 顺手修复。

2. **Q**: 是否需要在 Lean 端补"pp sign-normalize"兜底（参 C++ PR #16 在 sqf 层）？
   **A**: ✅ 不需要。Corpus.lean:6888 `squarefreefactorize_lex_ir` 已是 PR #16 的翻译版（含 if-flip 层），自带兜底。

3. **Q**: 是否同步更新 Lean Spec / 证明（若有）？
   **A**: 当前无 cont 相关 Spec/证明（grep 确认；MvGCD.lean `cont_spec` 是占位）。未来若加 `cont_correct` 类定理，应按"cont ≥ 0"为新 Spec。

4. **Q**: contMv 的 incomplete placeholder（L1171-1174）是否要本 fix 一起处理？
   **A**: 不处理。contMv 是 Phase F-impl-v2 范围（"对真正多变量输入返回剩余变量多项式"），与 cont sign 是正交问题。本 fix 仅同步 sign 约定。

---

## 第八部分：与 C++ PR #18 的关系

PR #18 描述里承诺：
> Lean Model.lean 端有自己的 `MvPolyZZ.contImpl` / `SparsePolyZZ.contImpl`（在 `feature/formal-proofs` 分支），同样问题。本 PR 后会在 feature/formal-proofs 分支独立同步修（不在本 PR 范围）。

本 fix 兑现此承诺。范围：3 处 sign-match 删除 + 注释/`#eval` 更新。
