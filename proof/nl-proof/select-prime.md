# 素数选取 L2 模型

> 对应 C++: `__select_prime`（polynomial_factorize_univar.hh:1388-1470）

## 算法

```
输入：f ∈ Z[x]
输出：(prime, factors_mod_p, irreducible_flag)

for up to 3 valid primes p:
  1. p ∤ lc(f)                     // 首项系数不为零 mod p
  2. deg(f mod p) = deg(f)          // 度数保持
  3. gcd(f mod p, (f mod p)') = 1   // 无平方 mod p
  4. DDF + EDF → factors_mod_p
  5. 如果 |factors| ≤ 1：返回 (p, [], irreducible=true)
  6. 更新最优（因子数最少的 p）
返回最优
```

## L2 模型

数学保证：**好的素数存在**。

引理：对任何 f ∈ Z[x], deg(f) ≥ 1, lc(f) ≠ 0：
- p ∤ lc(f) 对所有 p > |lc(f)| 成立（有限排除）
- deg(f mod p) = deg(f) 当 p ∤ lc(f)（首项保持）
- f mod p 无平方 ⇔ gcd(f, f') = 1 mod p，对几乎所有 p 成立
  （只有 disc(f) 的有限多个素因子使 f mod p 有重根）

所以好的 p 存在，C++ 枚举一定找到。

```lean
/-- 素数选取规约：好的素数存在。
    对应 C++ __select_prime（lines 1388-1470）。
    条件：p ∤ lc(f)，deg 保持，f mod p 无平方。
    C++ 枚举至多 3 个有效素数，选因子数最少的。
    L2 模型：假设选取成功（存在满足条件的素数）。-/
structure PrimeSelectionResult (f : Polynomial ℤ) where
  p : ℕ
  hp : Nat.Prime p
  factors : List (Polynomial (ZMod p))
  /-- p 不整除首项系数 -/
  lc_nonzero : (f.leadingCoeff : ZMod p) ≠ 0
  /-- f mod p 无平方 -/
  sqfree : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f)
  /-- factors 是 f mod p 的不可约分解 -/
  factored : FactorZpCorrect (Polynomial.map (Int.castRingHom (ZMod p)) f)
              (Polynomial.map (Int.castRingHom (ZMod p)) f).leadingCoeff factors

/-- 好的素数存在：对任何非零非常数 f ∈ Z[x]，存在满足条件的 p。-/
theorem good_prime_exists (f : Polynomial ℤ) (hf : f ≠ 0)
    (hd : 0 < f.natDegree) :
    ∃ r : PrimeSelectionResult f, True := by
  sorry -- 存在性：p > |disc(f)| * |lc(f)| 即可
```

## 审查

1. **数学正确性** ✓：好的素数存在是标准结果（判别式有限素因子）
2. **无跳步** ✓：三个条件各自独立
3. **Lean 可形式化** ✓：`Squarefree`、`leadingCoeff`、`ZMod` 均在 Mathlib
4. **工程问题**：`PrimeSelectionResult` 需要 `[Fact (Nat.Prime p)]` 实例
5. **边界**：f = 0 排除（hf），常数排除（hd > 0）
