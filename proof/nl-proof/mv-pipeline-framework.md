# 多变量因式分解参数化框架

> 状态：nl-proof v1
> 目标：mv_factor_correct 参数化定理（类似 factor_ZZ_correct）

---

## 0. C++ 结构回顾

```
__factor_multivar(f):
  1. sqf = squarefreefactorize(f)         // f = unit * ∏ gₖ^mₖ
  2. for each (gₖ, mₖ):
       if univariate(gₖ): factorize(gₖ)  // 递归到单变量
       else: wang_core(gₖ)               // Wang 核心
  3. combine results
```

框架需要抽象掉 Wang 的具体实现，只假设子过程正确。

---

## 1. 需要的子过程规约

### 1.1 多变量 SQF

不需要新定义。CLPoly 的 `squarefreefactorize` 用 GCD 做 content 提取 + Yun SQF，输出 `(gₖ, mₖ)` 列表。
框架只需要：`f = unit * ∏ gₖ^mₖ`，每个 gₖ squarefree。

但形式化 SQF 输出类型比较复杂。**简化方案**：框架不显式接受 SQF，而是假设"任何非零多项式可因式分解"（这是 UFD 给的）。参数化的重点在 Wang。

### 1.2 Wang 核心

Wang 的输入：squarefree primitive f ∈ MvPolynomial
Wang 的输出：不可约因子列表

这和 `MvFactorCorrect` 定义相同。框架说：如果 Wang 对每个 squarefree primitive 分量正确，则整体正确。

### 1.3 递归结构

Wang 递归减少变量数。框架需要体现：
- n = 0：常数，trivial
- n = 1：单变量因式分解（已验证 ✅）
- n ≥ 2：Wang

---

## 2. 框架定理

```lean
/-- 多变量因式分解管线：假设 Wang 对 n+1 变量的 squarefree primitive 多项式正确，
    则任何非零 n+1 变量多项式有不可约分解。

    结构：UFD SQF 分解 → 对每个 squarefree 分量调 wang → 合并。

    注：这里用 UFD 做 SQF（得到 squarefree 分量），用 wang 做每个分量的不可约分解。
    UFD 保证 squarefree 分量存在，wang 保证每个分量可因式分解。-/
theorem mv_factor_correct (n : ℕ)
    (f : MvPolynomial (Fin n) ℤ) (hf : f ≠ 0)
    -- 假设：Wang 对任何非零 squarefree 多项式给出正确分解
    (wang : ∀ g : MvPolynomial (Fin n) ℤ,
        g ≠ 0 → Squarefree g →
        ∃ result, MvFactorCorrect g result)
    : ∃ result : List (MvPolynomial (Fin n) ℤ), MvFactorCorrect f result
```

### 2.1 证明思路

1. f ≠ 0 → f 在 UFD 中有不可约分解 `f ~ ∏ pᵢ^eᵢ`
2. 定义 squarefree 分量：`gₖ = ∏{pᵢ : eᵢ = k} pᵢ`，每个 gₖ squarefree
3. 对每个 gₖ 调 wang → 得到不可约因子
4. 合并 → MvFactorCorrect

**但实际上**，步骤 1 已经给了不可约分解（UFD）。所以这个框架的证明和 `mv_factor_instantiate` 完全一样——直接用 UFD。

### 2.2 更有意义的框架

如果要让框架真正体现 Wang 的价值，需要：

```lean
theorem mv_factor_correct (n : ℕ)
    (f : MvPolynomial (Fin n) ℤ) (hf : f ≠ 0)
    -- Wang 作为函数（不是存在性）
    (wang : MvPolynomial (Fin n) ℤ → List (MvPolynomial (Fin n) ℤ))
    (hwang : ∀ g, g ≠ 0 → Squarefree g → MvFactorCorrect g (wang g))
    : ∃ result, MvFactorCorrect f result
```

这样 wang 是一个具体的函数，hwang 是它的正确性证明。框架说：有了正确的 Wang 函数，就能因式分解任何多项式。

但证明仍然只需 UFD（因为 MvFactorCorrect 是存在性的）。

### 2.3 真正的参数化框架

类比 `factor_ZZ_correct`，应该是：

```lean
theorem mv_factor_correct (n : ℕ)
    (f : MvPolynomial (Fin (n + 1)) ℤ) (hf : f ≠ 0)
    -- 子过程 1：单变量因式分解（已验证）
    (factor_univar : Polynomial ℤ → List (Polynomial ℤ))
    (hfactor_univar : ∀ g, g ≠ 0 → FactorZZCorrect g (factor_univar g))
    -- 子过程 2：Wang 核心
    (wang : MvPolynomial (Fin (n + 1)) ℤ → List (MvPolynomial (Fin (n + 1)) ℤ))
    (hwang : ∀ g, g ≠ 0 → Squarefree g → MvFactorCorrect g (wang g))
    -- 子过程 3：SQF（content 提取 + Yun）
    (sqf : MvPolynomial (Fin (n + 1)) ℤ →
           List (MvPolynomial (Fin (n + 1)) ℤ × ℕ))
    (hsqf : ∀ g, g ≠ 0 →
        Associated g ((sqf g).map (fun pr => pr.1 ^ pr.2)).prod
        ∧ ∀ pr ∈ sqf g, Squarefree pr.1 ∧ pr.1 ≠ 0 ∧ pr.2 ≥ 1)
    : ∃ result, MvFactorCorrect f result
```

**证明**：
1. hsqf 给出 f ~ ∏ gₖ^mₖ，每个 gₖ squarefree
2. 对每个 gₖ：hwang 给出 gₖ 的不可约因子列表 rₖ，gₖ ~ ∏ rₖ
3. 合并：f ~ ∏ₖ (∏ rₖ)^mₖ = ∏ₖ ∏ⱼ rₖⱼ^mₖ
4. 每个 rₖⱼ 不可约 ✓

这需要：`Associated (∏ gₖ^mₖ) f` + `Associated gₖ (∏ rₖⱼ)` → `Associated f (∏ₖⱼ rₖⱼ^mₖ)`。

---

## 3. Lean 形式化

最简版（和单变量 `factor_ZZ_correct` 平行）：

```lean
theorem mv_factor_correct {n : ℕ}
    (f : MvPolynomial (Fin n) ℤ) (hf : f ≠ 0)
    -- SQF 分解
    (sqf : MvPolynomial (Fin n) ℤ → List (MvPolynomial (Fin n) ℤ × ℕ))
    (hsqf : ∀ g, g ≠ 0 →
        Associated g ((sqf g).map (fun pr => pr.1 ^ pr.2)).prod
        ∧ ∀ pr ∈ sqf g, Squarefree pr.1 ∧ pr.1 ≠ 0 ∧ pr.2 ≥ 1)
    -- 每个 squarefree 分量的因式分解（Wang 或单变量）
    (factor_sqfree : MvPolynomial (Fin n) ℤ → List (MvPolynomial (Fin n) ℤ))
    (hfactor : ∀ g, g ≠ 0 → Squarefree g → MvFactorCorrect g (factor_sqfree g))
    : ∃ result : List (MvPolynomial (Fin n) ℤ), MvFactorCorrect f result := by
  -- Step 1: SQF 分解
  obtain ⟨hassoc, hprops⟩ := hsqf f hf
  -- Step 2: 对每个 squarefree 分量调 factor_sqfree
  -- Step 3: 合并
  -- 具体需要：List.bind/flatMap 每个 (gₖ, mₖ) 展开为 (factor_sqfree gₖ) 重复 mₖ 次
  sorry -- 组合证明，~30-50 行
```

---

## 4. 估计

| 内容 | 行数 |
|------|------|
| SQF 相关 Associated 引理 | ~20 |
| mv_factor_correct 框架证明 | ~40 |
| **总计** | **~60** |
