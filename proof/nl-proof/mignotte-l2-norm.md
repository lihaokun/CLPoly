# Mignotte bound L2 范数版本：M(f) ≤ ‖f‖₂

> 状态：nl-proof v2（调研完成：Parseval 需自证，Jensen/log凹 Mathlib 有）
> 目标：从 M(f) ≤ ‖f‖₁ 加强到 M(f) ≤ ‖f‖₂，精确匹配 C++ `__mignotte_bound`

---

## 0. 当前状态

已证（Recombine.lean）：`‖g‖_∞ ≤ C(n,n/2) · ‖f‖₁`
C++ 用：`B = C(n,n/2) · ‖f‖₂`
需要：`‖g‖_∞ ≤ C(n,n/2) · ‖f‖₂`

差距：Step 3 用了 `M(f) ≤ ‖f‖₁`（`mahlerMeasure_le_sum_norm_coeff`）。
需要改为 `M(f) ≤ ‖f‖₂`（Landau 不等式）。

---

## 1. M(f) ≤ ‖f‖₂ 的证明

### 1.1 定理（Landau 不等式）

```
∀ f : ℂ[X], f ≠ 0 → M(f) ≤ ‖f‖₂
```

其中 ‖f‖₂ = √(Σ|fₖ|²) 是系数的 L2 范数。

### 1.2 证明

**Step 1**：Mahler measure 的对数等于圆上积分（Mathlib `logMahlerMeasure` 定义）

```
log M(f) = (1/2π) ∫₀²π log|f(e^{iθ})| dθ
```

这是 Mathlib 的 `logMahlerMeasure` 定义本身。Mathlib 还有
`logMahlerMeasure_eq_log_leadingCoeff_add_sum_log_roots`
证明了这等于 `log|lc| + Σ log max(1, |αᵢ|)`。
两个定义的等价性由 Jensen 公式保证（Mathlib `JensenFormula.lean`）。

**Step 2**：Parseval 恒等式

```
(1/2π) ∫₀²π |f(e^{iθ})|² dθ = Σₖ |fₖ|² = ‖f‖₂²
```

证明：f(e^{iθ}) = Σ fₖ e^{ikθ}。
|f(e^{iθ})|² = (Σ fₖ e^{ikθ})(Σ f̄ⱼ e^{-ijθ}) = ΣΣ fₖ f̄ⱼ e^{i(k-j)θ}。
积分：(1/2π) ∫ e^{i(k-j)θ} dθ = δ_{kj}。
故 (1/2π) ∫ |f|² dθ = Σ |fₖ|² = ‖f‖₂²。✓

**Lean 路径**：这是有限 Fourier 正交性。Mathlib 可能有，否则自证 ~20 行。

**Step 3**：Jensen 不等式（log 凹性）

```
(1/2π) ∫ log|f(e^{iθ})|² dθ ≤ log((1/2π) ∫ |f(e^{iθ})|² dθ)
```

即 `E[log X] ≤ log E[X]`（Jensen 不等式，log 凹函数）。

取 X = |f(e^{iθ})|²，关于均匀测度 dθ/2π 在 [0, 2π] 上。

**Lean 路径**：Mathlib 应有 Jensen 不等式（`ConcaveOn.integral_le` 或类似）。

**Step 4**：组合

```
2 · log M(f) = (1/2π) ∫ log|f(e^{iθ})|² dθ    [Step 1: 2·logM = ∫ log|f|²]
             ≤ log((1/2π) ∫ |f(e^{iθ})|² dθ)   [Step 3: Jensen]
             = log ‖f‖₂²                         [Step 2: Parseval]
             = 2 · log ‖f‖₂

→ log M(f) ≤ log ‖f‖₂
→ M(f) ≤ ‖f‖₂  ✓
```

**注意**：Step 1 的 `2·log M(f) = ∫ log|f|²` 需要验证。
log|f|² = 2·log|f|。所以 (1/2π)∫ 2·log|f| dθ = 2·(1/2π)∫ log|f| dθ = 2·log M(f)。✓

---

## 2. Mathlib API 调研结果

| 需要 | Mathlib 状态 | 引理名 |
|------|-------------|--------|
| log 凹性 | ✅ | `strictConcaveOn_log_Ioi` |
| Jensen 不等式 | ✅ | `ConcaveOn.le_map_integral` |
| logMahlerMeasure 定义 | ✅ | `Polynomial.logMahlerMeasure` |
| Jensen 公式 | ✅ | `MeromorphicOn.circleAverage_log_norm` |
| M(f) = \|lc\|·∏max(1,\|αᵢ\|) | ✅ | `mahlerMeasure_eq_leadingCoeff_mul_prod_roots` |
| **Parseval 恒等式** | ❌ 不在 Mathlib | Norm.lean 标注 TODO |
| **M(f) ≤ ‖f‖₂** | ❌ 不在 Mathlib | 需自证 |

**阻塞点**：Parseval（`(1/2π)∫|f(e^{iθ})|²dθ = Σ|fₖ|²`）不在 Mathlib。

---

## 3. 选定方案：代数路径（避免积分 Parseval）

### 3.1 思路

不走 `logMahlerMeasure → Jensen 不等式 → Parseval` 的积分路径。
改走**代数路径**：直接证 `M(f)² ≤ ‖f‖₂²`，利用根的代数恒等式。

```
M(f)² = |lc|² · ∏ max(1, |αᵢ|²)        [mahlerMeasure_eq... 的平方]
      ≤ |lc|² · ∏(1 + |αᵢ|²)           [逐项：max(1,x) ≤ 1+x for x ≥ 0]
      = Σ|fₖ|² = ‖f‖₂²                  [代数 Parseval 恒等式]
```

### 3.2 代数 Parseval 恒等式

**引理**：设 f = lc · ∏ᵢ(X - αᵢ) ∈ ℂ[X]，deg f = n。则：

```
Σₖ₌₀ⁿ |fₖ|² = |lc|² · ∏ᵢ₌₁ⁿ (1 + |αᵢ|²)
```

**证明**（对 n 归纳，纯代数，不需要积分）：

**Base n = 0**：f = lc（常数）。LHS = |lc|²。RHS = |lc|² · (空积 = 1)。✓

**Step n → n+1**：设 f = g · (X - α) 其中 g = lc · ∏ᵢ₌₁ⁿ (X - αᵢ)，deg g = n。

IH：Σ|gₖ|² = |lc|² · ∏ᵢ₌₁ⁿ (1 + |αᵢ|²)。

f = g · (X - α)。设 g = Σ gₖ Xᵏ。则：
```
fₖ = gₖ₋₁ - α · gₖ  （其中 g₋₁ = 0, gₙ₊₁ = 0）
```

```
Σ|fₖ|² = Σ|gₖ₋₁ - α · gₖ|²
        = Σ(|gₖ₋₁|² + |α|²|gₖ|² - 2Re(gₖ₋₁ · ᾱ · ḡₖ))
        = Σ|gₖ₋₁|² + |α|² · Σ|gₖ|² - 2Re(ᾱ · Σ gₖ₋₁ḡₖ)
```

**关键引理**：`Σₖ gₖ₋₁ · ḡₖ = 0`。

证明：Σₖ₌₀ⁿ⁺¹ gₖ₋₁ḡₖ = Σₖ₌₁ⁿ⁺¹ gₖ₋₁ḡₖ（k=0 项 = g₋₁ḡ₀ = 0）
= Σⱼ₌₀ⁿ gⱼ · ḡⱼ₊₁（换 j = k-1）。

但 gₙ₊₁ = 0（deg g = n）。所以：
= Σⱼ₌₀ⁿ⁻¹ gⱼ · ḡⱼ₊₁

**这不一定等于 0！** 例如 g = 1 + X：g₀ḡ₁ = 1·1 = 1 ≠ 0。

**修正**：代数 Parseval 的归纳证明不走 "交叉项 = 0" 路径。

**正确的归纳**：

```
Σ|fₖ|² = Σ|gₖ₋₁ - αgₖ|²
        = Σ(|gₖ₋₁|² - αgₖḡₖ₋₁ - ᾱgₖ₋₁ḡₖ + |α|²|gₖ|²)
        = Σ|gₖ₋₁|² + |α|²Σ|gₖ|² - αΣgₖḡₖ₋₁ - ᾱΣgₖ₋₁ḡₖ
```

注意 Σ|gₖ₋₁|² = Σⱼ₌₋₁ⁿ⁻¹ |gⱼ|² = Σ|gⱼ|²（因 g₋₁ = 0）= Σ|gₖ|²。

所以：
```
Σ|fₖ|² = (1 + |α|²) · Σ|gₖ|² - αΣgₖḡₖ₋₁ - ᾱΣgₖ₋₁ḡₖ
```

交叉项 `Σgₖḡₖ₋₁` 和 `Σgₖ₋₁ḡₖ` 是共轭对，差不一定为 0...

**这条路走不通**（交叉项不消）。

### 3.3 修正：直接用积分 Parseval（自证）

代数路径的交叉项不消。回到积分路径。**自证有限 Fourier Parseval**：

**引理 (parseval_poly)**：
```
(1/2π) ∫₀²π |f(e^{iθ})|² dθ = Σₖ |fₖ|²
```

**证明**（初等，不需要 Fourier 分析理论）：

展开 |f(e^{iθ})|²：
```
|f(e^{iθ})|² = (Σⱼ fⱼ e^{ijθ})(Σₖ f̄ₖ e^{-ikθ}) = ΣΣ fⱼf̄ₖ e^{i(j-k)θ}
```

逐项积分：
```
(1/2π) ∫₀²π e^{inθ} dθ = { 1  if n = 0
                           { 0  if n ≠ 0
```

后者证明：n ≠ 0 时 ∫₀²π e^{inθ} dθ = [e^{inθ}/(in)]₀²π = (e^{2πin} - 1)/(in) = 0。

所以：
```
(1/2π) ∫₀²π |f(e^{iθ})|² dθ = ΣΣ fⱼf̄ₖ · δⱼₖ = Σ |fₖ|²  ✓
```

### 3.4 Lean 形式化路径

**Step 1**：证 `∫₀²π e^{inθ} dθ = 0` when `n ≠ 0`。
Lean：`intervalIntegral.integral_exp_mul_complex`（可能在 Mathlib，否则自证 ~15 行）。

**Step 2**：展开 `∫|f(e^{iθ})|²`。
对有限和（多项式的有限项），交换求和与积分。
Lean：`integral_finset_sum`（Mathlib）。

**Step 3**：应用正交性消去交叉项。

**Step 4**：Jensen 不等式。
```
log M(f) = (1/2π) ∫ log|f(e^{iθ})| dθ        [logMahlerMeasure 定义]
         ≤ (1/2) · log((1/2π) ∫ |f(e^{iθ})|² dθ)  [Jensen: E[log X] ≤ log E[X], applied to X = |f|²]
         = (1/2) · log(‖f‖₂²)                   [Parseval]
         = log(‖f‖₂)
```

Lean：`ConcaveOn.le_map_integral`（Mathlib）+ 上面的 Parseval。

**Step 5**：exp 两边得 M(f) ≤ ‖f‖₂。

### 3.5 估计

| 组件 | 行数 |
|------|------|
| ∫ e^{inθ} = 0 (n ≠ 0) | ~20 |
| Parseval 展开 + 正交消去 | ~30 |
| Jensen 不等式应用 | ~20 |
| M(f) ≤ ‖f‖₂ 组合 | ~15 |
| mignotte_bound 修改 | ~15 |
| **总计** | **~100** |

---

## 4. 对 mignotte_bound 的修改

改 Step 3 从 `mahlerMeasure_le_sum_norm_coeff` 到自证的 Landau 不等式：

```lean
-- 自证：Landau 不等式
private lemma landau_inequality (f : ℂ[X]) (hf : f ≠ 0) :
    f.mahlerMeasure ≤ Real.sqrt (∑ i in Finset.range (f.natDegree + 1), ‖f.coeff i‖ ^ 2) := by
  -- Parseval + Jensen + exp
  sorry

-- mignotte_bound 签名不变（仍用 natAbs），但 RHS 改用 L2 范数
-- 或：保持当前签名不变（L1 范数版本已够用），单独声明 L2 版本
```

**注意**：签名中使用 `Nat.sqrt` 处理平方根在 ℕ 中的问题，或者保持在 ℝ 中工作。
C++ 的 `__mignotte_bound` 用 `__isqrt_ceil(norm_sq)` 向上取整。
