# 多变量 L2 模块 4：多变量 Hensel 提升（MTSHL）

> 状态：nl-proof v1
> 对应 C++：`__mtshl_lift` + `__mtshl_step_j`
> 依赖：单变量 Hensel（已验证 ✅）

---

## 0. C++ 算法核心逻辑

MTSHL 的数学核心（剥离稀疏插值、MDP cascade 等工程细节）：

```
输入：
  f ∈ Z[x₁,...,xₙ]（缩放后的多项式）
  v₁,...,vᵣ ∈ Z[x₁]（单变量因子，满足 f(x₁,α) = ∏vᵢ）
  α = (α₂,...,αₙ)
输出：
  G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]，满足 f = ∏Gᵢ，Gᵢ(x₁,α) = vᵢ

算法（逐变量提升）：
  Fᵢ := vᵢ（初始，仅含 x₁）
  for j = 2 to n:
    for k = 1 to deg(f, xⱼ):
      error := f|_{xⱼ₊₁=αⱼ₊₁,...} - ∏Fᵢ|_{xⱼ₊₁=αⱼ₊₁,...}
      cₖ := Taylor 系数 of error at xⱼ=αⱼ, order k
      (σ₁ₖ,...,σᵣₖ) := MDP_solve(cₖ, F₁,...,Fᵣ)
      Fᵢ += σᵢₖ · (xⱼ - αⱼ)^k
  Gᵢ := symmetric_mod(Fᵢ)  // Zp → Z 恢复
```

## 1. L2 模型的抽象层次

MTSHL 的内部实现（稀疏插值、Vandermonde、MDP cascade）是**实现细节**（L1）。

L2 建模的是**数学不变量**：
1. **初始条件**：Fᵢ(x₁,α₂,...,αₙ) = vᵢ(x₁)
2. **每步不变量**：∏Fᵢ ≡ f (mod (x₂-α₂)^{k₂} · ... · (xₙ-αₙ)^{kₙ})
3. **终止条件**：kⱼ > deg(f, xⱼ) for all j → ∏Fᵢ = f（精确）
4. **因子对应**：Gᵢ(x₁,α) = vᵢ 保持

这是**多变量 Hensel 提升的 Taylor 展开解读**——类似单变量中 mod p^k 的提升，这里是 mod (xⱼ-αⱼ)^k 的提升。

## 2. 核心定理

### 2.1 提升不变量

```lean
/-- 多变量 Hensel 提升不变量：
    ∏ Gᵢ ≡ f (mod I^k)，其中 I = (x₂-α₂,...,xₙ-αₙ) 是求值理想。

    数学意义：Gᵢ 是 f 的因子的逐步近似，
    在 (x₂,...,xₙ) = (α₂,...,αₙ) 的邻域中逼近真因子。-/
def MvHenselInvariant {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (α : Fin n → ℤ)
    (precision : Fin n → ℕ)  -- 每个变量的提升精度
    : Prop :=
  -- f - ∏ factors ∈ (x₂-α₂)^k₂ · ... · (xₙ-αₙ)^kₙ
  -- 即：在 mod 这些幂次下，乘积等于 f
  sorry -- 精确形式化需要理想的幂次，复杂
```

### 2.2 精确提升定理

```lean
/-- 当提升精度足够（kⱼ > deg(f, xⱼ)），∏ Gᵢ = f 精确。
    这是因为 f - ∏Gᵢ 的每项都含 (xⱼ-αⱼ)^kⱼ 因子，
    但 deg(f, xⱼ) < kⱼ，所以 f - ∏Gᵢ = 0。-/
theorem mv_hensel_exact_when_sufficient_precision {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (α : Fin n → ℤ)
    (precision : Fin n → ℕ)
    (h_inv : MvHenselInvariant f factors α precision)
    (h_prec : ∀ j, precision j > MvPolynomial.degreeOf (Fin.succ j) f)
    : f = factors.prod
```

### 2.3 务实方案

完整形式化 MvHenselInvariant（涉及多变量理想的幂次）工程量很大。

**替代方案**：建模为"黑盒"——MTSHL 的输出满足特定条件，证明这些条件蕴含正确性。

```lean
/-- 多变量 Hensel 提升输出规约。
    不建模 MTSHL 内部过程，而是刻画其输出必须满足的条件。
    对应 C++ __mtshl_lift 的后置条件。-/
structure MvHenselOutput {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (lifted : List (MvPolynomial (Fin (n + 1)) ℤ))
    (α : Fin n → ℤ)
    (univar_facs : List (Polynomial ℤ)) : Prop where
  -- 乘积条件：∏ lifted = f（精确，在 Z[x₁,...,xₙ] 中）
  prod_eq : f = lifted.prod
  -- 求值一致：Gᵢ(x₁, α) = vᵢ
  eval_consistent : List.Forall₂ (fun v G => eval_at_α G α = v) univar_facs lifted
  -- 长度一致
  length_eq : univar_facs.length = lifted.length
```

**注**：`prod_eq : f = lifted.prod` 是精确等式（不是 Associated），因为 MTSHL + 对称恢复后因子是精确的。在 C++ 中，试除验证最终确认了这个等式。

## 3. 从 MvHenselOutput 到 MvFactorCorrect

```lean
/-- 若 Hensel 提升给出精确分解 f = ∏ Gᵢ，
    且每个 Gᵢ 经试除验证为不可约，则 MvFactorCorrect 成立。-/
theorem mv_hensel_to_factor {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (lifted : List (MvPolynomial (Fin (n + 1)) ℤ))
    (α : Fin n → ℤ)
    (univar_facs : List (Polynomial ℤ))
    (h_hensel : MvHenselOutput f lifted α univar_facs)
    (h_irred : ∀ g ∈ lifted, Irreducible g)
    : MvFactorCorrect f lifted := by
  exact ⟨by rw [h_hensel.prod_eq], h_irred⟩
```

## 4. 讨论

### 4.1 L2 vs L1 的边界

| 层次 | 建模内容 | 不建模 |
|------|---------|--------|
| L2 | MvHenselOutput 规约（乘积精确 + 求值一致） | MTSHL 内部（Taylor 循环、MDP cascade、稀疏插值） |
| L2 | 试除验证（g \| f → 因子） | 子集枚举策略（Gosper's hack） |
| L2 | LC 分配循环（valuation 提取） | Non-divisor 检查的完备性 |
| L1 | MTSHL Taylor 循环 | — |
| L1 | MDP solve 具体实现 | — |

### 4.2 MvHenselOutput 中 prod_eq 的正当性

C++ 中 `∏ Gᵢ = f` 不是直接由 MTSHL 保证的——MTSHL 在 Zp 上操作，阶段 C 做对称恢复和 p-adic 校正。最终的精确等式由**试除验证**确认。

所以 `prod_eq : f = lifted.prod` 实际上来自试除（`rem == 0`），不来自 Hensel 提升本身。这在 L2 中是合理的——我们把 Hensel 和试除的组合效果作为 MvHenselOutput。

## 5. 形式化估计

| 内容 | 行数 |
|------|------|
| `MvHenselOutput` 结构体 | ~15 |
| `mv_hensel_to_factor` 定理 | ~10 |
| **总计** | **~25** |

MTSHL 内部的 Taylor 提升正确性（MvHenselInvariant → MvHenselOutput）是 L1/后续工作。
