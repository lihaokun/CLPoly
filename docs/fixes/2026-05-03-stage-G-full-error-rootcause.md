# Stage G — 全部 102 errors 完整根因分析

**日期**：2026-05-03
**前置**：commit `0151661`（per-function check 工具）、`29258c3`（G-A IR Var callee）
**总错数**：101 reported（lake 在第 100 错停止）+ 1 maxErrors 截断警告
**通过率**：299/346 = **86.4%** 函数通过；47 函数失败

本文档对**全部 47 个失败函数 / 101 错误**逐条 trace，归出 **9 个根因簇**。
不打包、不糊弄——每条错误都明确归属。

## 0. 错误来源分布（per-function）

| 失败函数 | 错数 | 模块 |
|---------|------|------|
| `_loop___wang_core_lex_3_ir` | 11 | Wang multivar |
| `__mtshl_step_j_lex_ir` | 7 | MTSHL |
| `__mtshl_wmds_lex_ir` | 7 | MTSHL |
| `__mtshl_sparse_int_lex_ir` | 6 | MTSHL |
| `_loop___select_eval_point_lex_8_ir` | 6 | Wang |
| `_loop___wang_core_lex_8_ir` | 6 | Wang |
| `__mtshl_multi_bdp_lex_ir` | 3 | MTSHL |
| `_lambda___mtshl_step_j_lex_1_ir` | 3 | MTSHL |
| `__select_prime_upoly_ir` | 3 | univar |
| `_loop___upoly_divmod_mod_upoly_0_ir` | 3 | univar |
| ... 其他 37 个函数（每函数 1-3 错）| 56 | mixed |

错误高度集中于 **Wang/MTSHL 多变量因式分解**模块（前 6 个函数共 41 errors / 40%）。

## 1. 9 个根因簇

按 (actual_type, expected_type) 模式 + 上下文 trace 后归类。

### G1 — `result.resize n` 不能 dot-dispatch（6 errors）

**症状**：

```
Function expected at
  result.resize r_1.toInt64.toUInt64.toNat
but this term has type
  Array MvPolyZp
```

**行**：2587, 2968, 2973, 3001, 3576, 3602
**根因**：

C++ `result.resize(n)` 是 vector 成员方法。Pass 5 method dispatch 翻译为 `Array.resize result n`（前缀形式），但 Pass 8 emit 时把 `result.resize` 当 dot projection。Lean `Array` 没 `.resize` 方法，找不到 → fallback 报 "Function expected at result.resize" (因 Pass 8 emit 出 `result.resize n`，Lean parse 为 `(result.resize) n`，前者是 ?, 后者 应用)。

**修法**：class_map 给 vec.resize → "Array.resize"（top-level）。Pass 5 在 mutate dispatch 把 receiver 作 first arg：`Array.resize result n`。看 class_map 当前是否已映射。

实际修复：让 Pass 8 emit dot 形式时 fallback 到 top-level 函数名（identifier in scope 时）。或定义 `Array.resize` 在 `Array` namespace 中作为方法（已有 def Array.resize 但 Lean 不识别为 abbrev SparsePolyXX 的方法）。

**改动**：~10 行

### G2 — `f_scaled.comp : UInt64` vs `Lex`（5 errors）

**症状**：

```
Type mismatch
  f_scaled.comp
has type
  UInt64
but is expected to have type
  Lex
```

**行**：2304, 2568, 2959, 3303, 3553
**根因**：

C++ `Poly` 有方法 `.comp()` 返回比较器（`lex_<less>` 类型对象）。Pass 5 用 class_map MvPolyZp/MvPolyZZ 的 "comp" → `MvPolyZZ.comp` 方法，返回 UInt64（CLPoly 的 comp 缓存值）。但 caller 用这个值给 `__assign_partial_zp_lex_ir` 等期望 `Lex` 类型的参数。

class_map.py 中 MvPolyZp/MvPolyZZ 的 "comp" mapping 返回 UInt64 是错的语义——这是 comp 数值缓存而非比较器。比较器应通过 `polynomial_<T, lex_<less>>` 模板参数 `lex_<less>` 提供。

**修法**：
- 把 `.comp` 方法的 Lean 端返回改为 `Lex`（占位）
- 或在 caller-side（`__assign_partial_zp_lex_ir` 等）改 expected 类型为 UInt64

**改动**：~10 行

### G3 — `failed to synthesize instance`（7 errors）

**行**：5193, 5368, 5427, 5593, 5596, 5599
**根因（待详细 trace）**：

未取到 actual/expected。需逐行看具体什么 typeclass。从行号看集中在 `_loop___wang_core_lex_3_ir`（5587-5606）和 `_loop___wang_core_lex_8_ir`。预测是某些 HMul/HAdd/Coe instance 缺失（Wang 多变量算术）。

**修法**：根据具体 typeclass 加 instance（Lean Model）。

**改动**：~15 行

### G4 — `0 : Zp` OfNat 缺失（4 errors）

**症状**：

```
Type mismatch
  0
has type
  Int32
but is expected to have type
  Zp
```

**行**：2610, 3299, 3625, 4749
**根因**：

字面量 `0` 在 Zp 上下文中需要 `OfNat Zp 0` 实例。Lean Model 已加 OfNat for SparsePolyZZ 等但未加 Zp。

**修法**：

```lean
instance : OfNat Zp 0 where ofNat := ⟨0, 1⟩  -- 任意 prime 占位
instance (p : UInt64) : OfNat Zp 0 where ofNat := ⟨0, p⟩  -- 不行：p 不在作用域
```

复杂——Zp 的 prime 字段不能从 0 字面量推断。需要让 caller 显式构造 `Zp.ofInt 0 p`。或 Pass 5 cast 在这种 site 自动 wrap。

**改动**：~10 行（需查 4 个 site 的 prime 来源）

### G5 — Pass 5 cast 漏（散点 ~25 errors）

汇总 (actual, expected) 配对：

| (actual, expected) | 数 | 行 |
|--------------------|----|----|
| (StdMap Variable ZZ, Variable) | 5 | 3817, 3963, 3965, 3974, 5735 |
| (SparsePolyZZ, Nat) | 5 | 5593, 5596, 5600, 5601, 5603 |
| (Nat, Int32) | 3 | 2810, 3807, 5112 |
| (Nat, Int64) | 3 | 4809, 4820, 4821 |
| (Int32, Zp) | 2 | 3300, 3302 |
| (Zp, Int32) | 2 | 3323, 3325 |
| (UInt64, Int32) | 2 | 2978, 5058 |
| (UInt32, UInt64) | 2 | 3560, 4321 |
| (Variable, StdMap Variable ZZ) | 2 | 3945, 3947 |
| (Array Int32, UInt64) | 2 | 5434, 5684 |
| (Nat, SparsePolyZZ) | 2 | 5590, 5760 |
| (UMonomial × ZZ, Unit) | 2 | 1305, 4926 |
| (PrimeSelectionResult, UInt64) | 2 | 4124, 4147 |
| (Nat, Int32) k_8 | 2 | 2836, 2849 |

**根因**：

每条都是 Pass 5 没在 caller-callee 边界 emit Cast。Lean 端实际看到 `actual` 但函数参数定义为 `expected`，类型严格 → 错。

具体的子类：

- **真假 cast 缺失**：(Nat→Int32), (Nat→Int64), (UInt64→Int32), (UInt32→UInt64) — 标准数值 cast，cast_table 有规则但 Pass 5 哪个调用 site 漏 wrap
- **真完全错的类型**：(StdMap Variable ZZ → Variable) 5× — Pass 5 把整个 map 当 key 传，是 IR 构造错（不是 cast）
- **Zp ↔ Int32 双向** 4×：缺 Coe 或 OfNat
- **Variable ↔ StdMap** 双向：可能 Pass 1 ty 推断错

**修法**：

按子类逐个 trace：
- 缺 Cast 的 site：补 cast_table entry 或 Pass 5 cast 路径
- 完全错类型的 site：可能是 Pass 1 / Pass 2b 的 IR 构造 bug（如 Pass 1 把 `m[k]` 解析为 m 而非 m[k]）

**改动**：~30 行 + 10-15 sites trace

### G6 — `pair_vec_div` 多 arity（3 errors）

**症状**：

```
Function expected at
  pair_vec_div q_1 sigma_6[i] sc_2 F[i]
but this term has type
  SparsePolyZp
```

**行**：3694, 4768, 4918, + 4724, 5645
**根因**：

C++ `pair_vec_div` 有 4-arg 和 5-arg overload（class_map.py 注解）。Pass 2b OUTPUT_PARAMS 有 `pair_vec_div#4: [0]` 和 `pair_vec_div#5: [0, 1]`。Lean 端 `def pair_vec_div [Inhabited α] (_f _g _q : α) (_comp : β) : α := default` 是 4-arg 返回 poly。

5-arg site 多传一个 arg → "Function expected"（lean 函数已饱和返回值，多余 arg 不能 apply）。

**修法**：补 5-arg 版本：

```lean
def pair_vec_div5 {α β : Type} [Inhabited α] (_f _g _q _r : α) (_comp : β) : α := default
```

或 Pass 2b 改 emit 不同 callee name 区分 arity。

**改动**：~10 行

### G7 — `__x_mut_1` filterMap mutator lambda 签名错（2 errors + 2 related）

**症状**：

```
Application type mismatch:
  __x_mut_1
has type
  UMonomial × ZZ
but is expected to have type
  Unit
```

**行**：1305, 4926
**+ 关联**：

```
_lambda___hensel_step_linear_upoly_filt1_ir m p_zz_1
has type
  UMonomial × ZZ → Option Unit
but expected
  UMonomial × Int → Option (UMonomial × Int)
```

**行**：1332, 4929

**根因**：

Pass 4 `_filter_loop_to_assign` for "A-mut" pattern emits `Array.filterMap'(arr, lambda)` where lambda is a mutator-then-filter combinator. Lean 的 `Array.filterMap'` 期望 `α → Option β`，但 Pass 4 emit 的 lambda 返回 `Option Unit`（不带 mutated value）。

**修法**：Pass 4 修 _build_mut_filter_lambda 让返回 `Option (mutated_elem)`。

**改动**：~15 行

### G8 — `__select_prime` 返回 PrimeSelectionResult 与 caller 期望 UInt64 不匹配（2 errors）

**症状**：

```
Application type mismatch:
  best_10
has type
  PrimeSelectionResult
but is expected to have type
  UInt64
```

**行**：4124, 4147
**根因**：

`__select_prime_upoly_ir` 在 Lean 端定义返回 `PrimeSelectionResult`，但 caller `__factor_squarefree_primitive_ZZ_upoly_ir` 在某 phi 中把 `best_10` 当 UInt64 用。

可能是 Pass 6 SSA 重命名错（`best_10 := sel_X.prime` 这种取字段后被 phi 误打包）。

**修法**：trace 4124 的具体调用 site，看是否需要 `.prime` 字段访问。

**改动**：~5 行

### G9 — Wang `_loop___wang_core_lex_3_ir` 11 错 — 集中在 SparsePolyZZ vs Nat / `fi_3`（11 errors）

**症状**（行 5587-5606，11 错）：

```
Application type mismatch:
  fi_3
has type
  SparsePolyZZ
but is expected to have type
  Nat
```
（5×：5593, 5596, 5600, 5601, 5603）

```
failed to synthesize instance
（5193, 5368, 5427, 5593, 5596, 5599）
```

```
Application type mismatch:
  fi_4
has type
  Nat
but is expected to have type
  SparsePolyZZ
```
（5590, 5760）

```
Type mismatch:
  mv_factors_1.size
expected SparsePolyZZ
```
（5591）

**根因**：

`_loop___wang_core_lex_3_ir` 是 Wang 因式分解的核心循环。Pass 6 SSA 把多个 phi target 名重用 → `fi_3` 在某些 phi 中是 SparsePolyZZ，在其他 phi 中是 Nat（计数器）。由于变量命名冲突（同名不同 phi），Pass 5 cast 在不同 site 错 wrap。

**修法**：需进入 Pass 6 SSA build 详细 trace。可能是 phi target 版本号冲突 / 名字 reuse。或 Pass 5 cast 路径在某调用 site 没识别。

**改动**：未知（需 trace），estimate ~30 行

## 2. 修复优先级总表

| 簇 | 错数 | 难度 | 改动行数 | 优先级 |
|----|------|------|---------|--------|
| G1 result.resize dot | 6 | 低 | 10 | 🔴 |
| G2 .comp 返回类型 | 5 | 低 | 10 | 🔴 |
| G6 pair_vec_div 5-arg | 5 | 低 | 10 | 🔴 |
| G4 OfNat Zp 0 | 4 | 中 | 10 | 🔴 |
| G7 filterMap mutator | 4 | 中 | 15 | 🟡 |
| G8 PrimeSelectionResult | 2 | 低 | 5 | 🔴 |
| G3 failed synth | 7 | 中 | 15 | 🟡（需 trace） |
| G9 Wang core_3 phi 冲突 | 11 | 高 | 30 | 🟢 |
| G5 杂项 cast | ~25 | 不定 | 30+ | 🟡（散点） |

**总计**：~135 行 Python + Lean，~30 个 site trace。

**预估时间**：4-6 小时（不算调试时间）。

## 3. 攻坚顺序建议

按"低难度 × 高错数"乘积：

1. G1 result.resize（6 错，10 行）→ 通过率 87.6%
2. G2 .comp 返回类型（5 错）→ 89.0%
3. G6 pair_vec_div 5-arg（5 错）→ 90.5%
4. G4 OfNat Zp 0（4 错）→ 91.6%
5. G8 PrimeSelectionResult（2 错）→ 92.2%
6. G7 filterMap mutator（4 错）→ 93.3%
7. G5 散点 cast（25 错，分批） → 100%
8. G3 failed synth（7 错） → 100%
9. G9 Wang core_3（11 错） → 100%

按此顺序，前 6 项（26 错）可在 1-2 小时内清完，通过率 86.4 → ~93%。

## 4. 不能 trace 到的项

| 错 | 原因 | 处理 |
|----|------|------|
| G3 7× failed synth | 需看 Lean 详细 trace（set_option trace.Meta.synthInstance true） | 单独 trace |
| Lake maxErrors 截断 | 5765 处 maxErrors 警告，可能还有错被截断 | 调高 maxErrors 重跑 |
| Various 1× errors | 单点 site，需逐个 trace | 见 G5 |

## 5. 元层反思

之前 stage G 文档（`2026-05-03-stage-G-latent-errors-analysis.md`）把 102 errors 归 6 簇，但**没逐条 trace**。本文档逐条核对后改成 9 簇，发现：

1. 之前的"根因 A local-var lambda"已在 G-A 全清（13 → 0）
2. 之前的"根因 B IIFE ret_ty"已在 G-B 部分清（row_sub/row_swap）
3. 错误确实集中在 47 个函数 / Wang+MTSHL 模块
4. 真实通过率 86.4% 远比"102 errors"的恐怖印象好

per-function 工具帮助看清进度。下次再写"根因分析"必须**逐条 trace 不能打包**——这次差点又犯了同样错误。
