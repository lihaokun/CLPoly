# FactorZZ.lean 自然语言证明

> 状态：已形式化，0 sorry

---

## factor_ZZ_correct

**定理**：假设 Zp 因式分解、Hensel 提升、因子重组各自正确，则组合结果满足 FactorZZCorrect。

**前置条件**：
- f ∈ Z[x]，f ≠ 0，f 本原
- p 素数，k > 0
- f mod p 无平方（好素数）
- deg(f mod p) = deg(f)（p 不整除 lc(f)）

**证明**：

1. **f mod p ≠ 0**：若 f mod p = 0 则不可能 Squarefree（not_squarefree_zero），矛盾

2. **Zp 因式分解**：对 fp = f mod p 应用 FactorZpCorrect，得到
   - fp = C(lc) * ∏(fᵢ^eᵢ)

3. **构造 Hensel 输入**：
   - facs_p = [C(lc), f₁^e₁, f₂^e₂, ...]
   - 乘积 = C(lc) * ∏(fᵢ^eᵢ) = fp（由 List.prod_cons）

4. **Hensel 提升**：对 facs_p 应用 HenselCorrect，得到
   - f mod p^k = (hensel facs_p).prod

5. **因子重组**：对 Hensel 输出应用 RecombineCorrect
   - RecombineCorrect f result ≡ FactorZZCorrect f result（定义完全相同）
   - 所以直接得到 FactorZZCorrect

**关键观察**：RecombineCorrect 和 FactorZZCorrect 定义完全相同（f = ∏ gᵢ 且每个 gᵢ 不可约），因此不需要任何额外转换。整个证明只是串联三个子过程假设。
