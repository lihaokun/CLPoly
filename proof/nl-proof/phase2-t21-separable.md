# T2.1: X^{p^d} - X 在 ZMod p 上是 Separable 的

> 状态：已形式化，0 sorry

---

## 定理陈述

设 p 为素数，d ≥ 1，则 X^(p^d) - X ∈ (ZMod p)[x] 是 separable 的。

```lean
theorem X_pow_sub_X_separable (d : ℕ) (hd : 0 < d) :
    Separable (X ^ (p ^ d) - X : Polynomial (ZMod p))
```

## 定义回顾

Separable f ↔ IsCoprime f f'（f 与其形式导数互素）

## 证明

1. **计算导数**：
   - f = X^(p^d) - X
   - f' = (p^d) · X^(p^d - 1) - 1
   - 在 ZMod p 中，p^d ≡ 0（因为 p ∣ p^d，而 d ≥ 1）
   - 所以 f' = 0 · X^(p^d - 1) - 1 = -1

2. **互素性**：
   - f' = -1 是单位（可逆元）
   - 任何多项式与单位互素：IsCoprime f (-1) 成立
   - 所以 Separable f

## 依赖的 Mathlib 路径

- `Polynomial.derivative_sub`：导数的减法
- `Polynomial.derivative_pow`（或 `derivative_X_pow`）：(X^n)' = n · X^(n-1)
- `Polynomial.derivative_X`：X' = 1
- `CharP.cast_eq_zero`：ZMod p 中 char p 特征 → p 的倍数 = 0
- `Separable` 定义 = `IsCoprime f f'`
- `isCoprime_one_right` 或 `IsCoprime.neg_right`：与单位/可逆元互素

## 预期难度

很短（~20-30 行），核心就是 f' = -1 然后 -1 是单位。
