# 多变量 L2 模块 3：试除重组

> 状态：nl-proof v1
> 对应 C++：`__wang_core` lines 2360-2426（Zassenhaus 子集枚举 + 试除）

---

## 0. C++ 算法逻辑

```cpp
// 简化后的试除重组
g_remaining = f;
for s = 1, 2, ..., while 2*s <= mv_T.size():
  for each subset S of mv_T with |S| = s:
    prod = pp(∏_{i∈S} normed[i])    // 子集乘积的原始部分
    q, rem = divmod(g_remaining, prod)
    if rem == 0:                     // 整除 → prod 是真因子
      verified.push_back(prod)
      g_remaining = q
      remove S from mv_T
      s = 1; break
if mv_T.size() <= 1:
  verified.push_back(g_remaining)    // 剩余部分也是因子
```

## 1. 数学核心

试除的正确性是平凡的：**如果 g | f 且 g 不可约，则 g 是 f 的不可约因子**。

但 C++ 的实际逻辑更精细：
1. 从 Hensel 因子子集乘积构造候选 `prod`
2. 取 primitive part（`pp`）
3. 检查 `prod | g_remaining`
4. 若整除，提取因子，更新 `g_remaining`

L2 需要建模这个**循环过程**，证明：
- 每次成功试除后 `g_remaining` 缩小
- 终止时 `f = ∏ verified × g_remaining`
- 每个 verified 因子不可约（从 Hensel 唯一性 + Mignotte 推出）

## 2. L2 模型

### 2.1 单次试除

```lean
/-- 单次试除：若 candidate | f_remaining 且 candidate 不可约，
    则 candidate 是因子，余式是新的 f_remaining。
    对应 C++ 的 divmod + rem == 0 检查。-/
structure TrialDivStep {σ : Type*}
    (f_remaining candidate : MvPolynomial σ ℤ) : Prop where
  dvd : candidate ∣ f_remaining
  irred : Irreducible candidate
```

### 2.2 试除循环

```lean
/-- 试除重组结果：f = ∏ verified_factors × remaining。
    对应 C++ 的 verified 列表 + g_remaining。-/
structure TrialDivResult {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ) : Prop where
  -- f = remaining × ∏ verified
  prod_eq : f = remaining * verified.prod
  -- 每个因子不可约
  all_irred : ∀ g ∈ verified, Irreducible g
```

### 2.3 试除到因式分解

```lean
/-- 若试除消耗了所有因子（remaining 是 unit），
    则 verified 是 f 的不可约分解。-/
theorem trial_div_to_factor {σ : Type*} [DecidableEq σ]
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ)
    (h : TrialDivResult f verified remaining)
    (hrem : IsUnit remaining) :
    MvFactorCorrect f verified
```

**证明**：
- `f = remaining × ∏ verified`（from h.prod_eq）
- `remaining` is unit → `Associated f (∏ verified)`
- 每个因子不可约（from h.all_irred）
- 所以 `MvFactorCorrect f verified` ✓

## 3. Lean 形式化

```lean
theorem trial_div_to_factor ... := by
  refine ⟨?_, h.all_irred⟩
  -- Associated f verified.prod
  rw [h.prod_eq]
  exact (Associated.mul_left verified.prod hrem.associated_one).symm.trans
    (by simp)
```

~15 行。

## 4. 形式化估计

| 内容 | 行数 |
|------|------|
| `TrialDivStep` 结构体 | ~5 |
| `TrialDivResult` 结构体 | ~10 |
| `trial_div_to_factor` 定理 | ~15 |
| **总计** | **~30** |
