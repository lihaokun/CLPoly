# Mignotte bound L2 范数版本：M(f) ≤ ‖f‖₂

> 状态：nl-proof v1
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

## 2. Mathlib API 调研需求

| 需要 | 搜索关键词 | 备注 |
|------|-----------|------|
| Parseval（有限 Fourier） | `integral_exp`, `orthogonal` | 可能需自证 |
| Jensen 不等式 | `ConcaveOn`, `integral_le_log` | Mathlib 应有 |
| log 凹性 | `concaveOn_log`, `StrictConcaveOn` | Mathlib 应有 |
| ∫ log\|f(e^{iθ})\| = logMahlerMeasure | `logMahlerMeasure` 定义 | Mathlib 有 |

---

## 3. 对 mignotte_bound 的修改

当前 Step 3 用 `mahlerMeasure_le_sum_norm_coeff`。改为：

```lean
-- 新 Step 3: M(f) ≤ ‖f‖₂ (Landau 不等式)
have h_landau : (f.map (Int.castRingHom ℂ)).mahlerMeasure ≤
    Real.sqrt ((Finset.range (f.natDegree + 1)).sum
      (fun j => ((f.coeff j).natAbs : ℝ) ^ 2)) := by
  sorry -- Landau 不等式证明（Jensen + Parseval）
```

签名改为：
```lean
theorem mignotte_bound (f g : Polynomial ℤ) (hf : f ≠ 0) (hg : g ∣ f) :
    ∀ i, (g.coeff i).natAbs ≤
      Nat.choose f.natDegree (f.natDegree / 2) *
      Nat.sqrt ((Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs ^ 2))
```

**注**：这需要在 Nat 中处理 √，可能引入向上取整问题。
替代方案：保持 ℝ 中工作，最后用 `Nat.cast` 转换。

---

## 4. 复杂度评估

| 组件 | 行数 | 难度 |
|------|------|------|
| Parseval（有限 Fourier 正交） | ~30 | 中 |
| Jensen 不等式（log 凹 + 积分） | ~20（引用 Mathlib） | 中 |
| Landau 不等式（组合） | ~20 | 低 |
| mignotte_bound 签名修改 + 证明调整 | ~20 | 低 |
| **总计** | **~90** | |

---

## 5. 替代方案

如果 Jensen 不等式在 Lean 中证明困难，可以用替代路径：

**替代**：直接证 M(f)² ≤ ‖f‖₂² 而不经过积分。

M(f)² = |lc|² · ∏ max(1, |αᵢ|²)

需要证：|lc|² · ∏ max(1, |αᵢ|²) ≤ Σ|fₖ|²

这等价于：|lc|² · ∏ max(1, |αᵢ|²) ≤ |lc|² · ∏(1 + |αᵢ|²)

即：∏ max(1, |αᵢ|²) ≤ ∏(1 + |αᵢ|²)

这由 max(1, x) ≤ 1 + x（对 x ≥ 0）逐项成立。✓

**但**需要等式 Σ|fₖ|² = |lc|² · ∏(1 + |αᵢ|²)，这就是 Parseval 的根版本。

证明 Σ|fₖ|² = |lc|² · ∏(1 + |αᵢ|²)：
对于 f = lc · ∏(x - αᵢ)，展开 Σ|fₖ|² 用 Vieta 公式...
这也需要类似 Parseval 的论证。

**结论**：无论哪条路，都需要某种形式的 "Σ|fₖ|² = |lc|² · ∏(1+|αᵢ|²)" 或 Parseval。
Jensen 公式路径最干净（Mathlib 已有基础设施）。
