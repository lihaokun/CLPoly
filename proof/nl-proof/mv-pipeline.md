# 多变量因式分解 Pipeline 框架

> 状态：nl-proof v1
> 目标：多变量 Pipeline 骨架（类似 factor_ZZ_correct）

---

## 0. 对应 C++ 结构

```
__factor_multivar(f):
  1. sqf = squarefreefactorize(f)     -- 无平方分解
  2. for each (gₖ, mₖ) in sqf:
       if univariate(gₖ):
         factors_k = factorize_univar(gₖ)
       else:
         factors_k = wang_core(gₖ)
  3. return ∏ₖ factors_k^mₖ
```

## 1. Pipeline 定理

### 1.1 直接版（UFD）

```lean
/-- 多变量因式分解：任何非零 MvPolynomial 有不可约分解。
    直接由 MvPolynomial σ ℤ 是 UFD 推出。-/
theorem mv_factor_instantiate {σ : Type*} [DecidableEq σ] [Fintype σ]
    (f : MvPolynomial σ ℤ) (hf : f ≠ 0)
    : ∃ result : List (MvPolynomial σ ℤ), MvFactorCorrect f result
```

**证明**：`WfDvdMonoid.exists_factors f hf`，与单变量 `recombine_correct` 完全一致。

**前提**：需要 `MvPolynomial σ ℤ` 是 `WfDvdMonoid`（即 UFD）。
Mathlib 中 `MvPolynomial σ R` 在 `R` 是 UFD + `σ` 是 `Fintype` 时是 UFD。
ℤ 是 UFD ✓，`σ` 要求 `Fintype` ✓。

### 1.2 参数化版（Wang Pipeline）

```lean
/-- 多变量因式分解管线：SQF + Wang 组合。
    假设子过程正确，则组合满足 MvFactorCorrect。

    结构：f → SQF 分解 → 对每个无平方分量调 wang → 合并。-/
theorem mv_factor_via_wang {σ : Type*} [DecidableEq σ] [Fintype σ]
    (f : MvPolynomial σ ℤ) (hf : f ≠ 0)
    -- 假设：任何无平方本原多变量多项式可因式分解
    -- （由 Wang 或单变量因式分解提供）
    (factor_sqfree : ∀ g : MvPolynomial σ ℤ,
        g ≠ 0 → Squarefree g →
        ∃ result, MvFactorCorrect g result)
    : ∃ result : List (MvPolynomial σ ℤ), MvFactorCorrect f result
```

**证明思路**：
1. f ≠ 0 → f 在 UFD 中有有限个不可约因子
2. 设 f = u · p₁^e₁ · ... · pₘ^eₘ（不可约分解，由 UFD）
3. 定义 gₖ = ∏{pᵢ : eᵢ = k} pᵢ（按重数分组，无平方）
4. f = u · ∏ₖ gₖ^k
5. 对每个 gₖ：gₖ 无平方 → factor_sqfree(gₖ) 给出不可约因子
6. 合并 → MvFactorCorrect

**但实际上**这个证明不需要 SQF——直接用 UFD 即可（和 1.1 一样）。参数化版的价值在于展示 factor_sqfree 可被 Wang 实例化。

### 1.3 更有意义的参数化：递归结构

Wang 的递归结构是：n 变量 → 求值得 (n-1) 变量 → ... → 1 变量（单变量）。

```lean
/-- 多变量因式分解：递归减少变量数。
    base case (n=0): 常数，trivial
    step (n+1): 求值 + 单变量分解 + Wang 提升 -/
theorem mv_factor_recursive (n : ℕ)
    (f : MvPolynomial (Fin n) ℤ) (hf : f ≠ 0)
    -- 假设：单变量因式分解正确（已验证 ✅）
    (factor_univar : ∀ g : Polynomial ℤ, g ≠ 0 →
        ∃ result, FactorZZCorrect g result)
    -- 假设：Wang 核心正确（对 n+1 变量，给定 n 变量的递归分解）
    (wang : ∀ m, m < n →
        ∀ g : MvPolynomial (Fin (m + 1)) ℤ, g ≠ 0 →
        ∃ result, MvFactorCorrect g result)
    : ∃ result : List (MvPolynomial (Fin n) ℤ), MvFactorCorrect f result
```

**这个版本最有意义**——它展示了 Wang 的递归下降结构：
- n = 0：常数多项式，trivial
- n = 1：通过 `finSuccEquiv` 转为 `Polynomial ℤ`，用 `factor_univar`
- n ≥ 2：用 Wang（假设 wang 对 < n 变量成立）

## 2. Lean 形式化

### 2.1 直接版

```lean
theorem mv_factor_instantiate {σ : Type*} [DecidableEq σ] [Fintype σ]
    (f : MvPolynomial σ ℤ) (hf : f ≠ 0)
    : ∃ result : List (MvPolynomial σ ℤ), MvFactorCorrect f result := by
  obtain ⟨factors, hirred, hassoc⟩ := WfDvdMonoid.exists_factors f hf
  exact ⟨factors.toList, hassoc.symm, fun g hg => hirred g (Multiset.mem_toList.mp hg)⟩
```

### 2.2 Mathlib 依赖确认

- `WfDvdMonoid (MvPolynomial σ ℤ)`：需要 `[Fintype σ]` + `ℤ` 是 UFD
  - Mathlib: `MvPolynomial.instUniqueFactorizationMonoid` 或 `MvPolynomial.wfDvdMonoid`
  - 需要搜索确认实例名称

### 2.3 工程问题

- `MvPolynomial σ ℤ` 的 `Associated` 和 `Irreducible` 应该继承自代数结构
- `Multiset.toList` + `Multiset.prod_toList` 同单变量情况

## 3. 形式化估计

| 内容 | 行数 |
|------|------|
| `mv_factor_instantiate`（UFD 版） | ~10 |
| Pipeline 文件框架 | ~30 |
| **总计** | **~40** |

后续 Wang L2 完成后，再写 `mv_factor_via_wang` 的完整实例化。
