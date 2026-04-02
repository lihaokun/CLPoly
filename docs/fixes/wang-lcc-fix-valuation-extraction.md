# 修正方案：Wang LCC GCD 匹配 → Valuation 提取

> 对应 workflow.md §5.1 修正方案文档
> 调研报告：[`docs/research/wang-lcc-research.md`](../research/wang-lcc-research.md)
> 关联历史：`wang-lc-fix-power-distribution.md` 曾提出 valuation 提取，后被 GCD 匹配取代；本文档恢复 valuation 提取并基于完整调研（GCL §8.7 + SymPy + FLINT）重新论证

---

## 第一部分：复现与定位

### 最小复现用例

`test/test_repro_bivar4fac.cc` 第 87-102 行（确定性复现）：

```cpp
polynomial_ZZ f1 = -two*x*y + y*y - one;     // lc(f1, x) = -2y
polynomial_ZZ f2 = -three*x*y - two*y + one;  // lc(f2, x) = -3y
polynomial_ZZ f3 = x*x + x*y + two*x;         // lc(f3, x) = 1
polynomial_ZZ f4 = three*x*y - three*y - two;  // lc(f4, x) = 3y
polynomial_ZZ f = f1 * f2 * f3 * f4;
auto fac = factorize(f);
// 预期: 5 个非平凡因子 (f3 = x(x+y+2) 可进一步分解)
// 实际: 4 个非平凡因子 (LC 分配失败，某些因子未能拆出)
```

对应的随机测试 `crosscheck_flint_mpoly_factor_random_bivar_4fac` 有 ~1.8% 失败率。

### 出错路径

```
factorize(f)
  → __factor_multivar(f)
    → __wang_core(f)
      → __wang_leading_coeff(f, univar_factors, eval_point, main_var, uni_content)
        → Step 2: factorize(L) → L = 18y³, gamma=18, lc_factors={(y, 3)}
        → Step 3 (line 1400-1445): GCD 匹配循环
          → lj=y, ej=3, lj_val=α
          → gcd(|w[0]|, α³) vs gcd(|w[1]|, α³) vs ...
          → 多个 w[i] 含 α 的幂次 → ambiguous=true → return result (FAIL)
      → 所有求值点失败 → 返回未分解
```

关键失败点：`polynomial_factorize_wang.hh:1400-1445`（Step 3 GCD 匹配循环）。

### 预期行为 vs 实际行为

| | 预期 | 实际 |
|---|------|------|
| Step 3 | `y³` 拆分为 `y¹·y¹·y¹`，分别分配给 u₁, u₂, u₄ | `y³` 整块匹配，多个因子的 GCD 值相等 → ambiguous → FAIL |
| `__wang_leading_coeff` | `success=true`, σ = [y, y, 1, y] | `success=false` |
| `factorize(f)` | 5 个非平凡因子 | 4 个非平凡因子（LC 分配失败导致部分因子合并） |

---

## 第二部分：根因分析

### 症状

`__wang_leading_coeff` 对所有求值点返回 `success=false`，因为 Step 3 的 GCD 匹配在 `eⱼ > 1` 且需要拆分给多个因子时总是报 ambiguous。

### 根因

Step 3（line 1400-1445）将 `(lⱼ, eⱼ)` 作为**原子单元**，试图整体分配给唯一一个因子：

```cpp
// line 1402-1445 (当前代码核心)
for (auto& [lj, ej] : lc_factors)
{
    ZZ lj_pow(1);
    for (uint64_t e = 0; e < ej; ++e)
        lj_pow = lj_pow * abs(lj_val);    // ← 整块: |Eⱼ|^eⱼ

    // 找唯一最佳匹配
    for (size_t i = 0; i < r; ++i)
    {
        ZZ g = gcd(abs(w[i]), lj_pow);    // ← GCD 整块比较
        // ... 找 gcd 最大且唯一的 i
    }
    // ambiguous → FAIL
    sigma[best_i] *= lj^ej;               // ← 整块给一个因子
}
```

当 `eⱼ > 1` 且该 LC 因子需要分配给多个提升因子时（如 `y³ → y·y·y`），不存在"唯一赢家"：多个 `w[i]` 都含有 `lⱼ(α)` 的幂次，GCD 可能相等，代码宣告 ambiguous 并失败。

**这不是可以用 workaround 绕过的问题——整块分配在数学上就无法表达 power splitting。**

### 根因 vs 症状

- **症状**：MTSHL Taylor 循环不收敛（因为 LC 目标 τᵢ 全是常数，校正无效）
- **根因**：上游 LC 分配算法无法拆分幂次 → τᵢ 错误 → 下游必然失败
- 之前的 D_j 上界截断和 LC 校正增强只是在掩盖根因

---

## 第三部分：参考实现对照

### 3.1 GCL §8.7 (pp.377-378)：Wang 原始算法

GCL 描述了 Wang's LC Predetermination 的三个条件：

> **(3)** each factor of the leading coefficient `aᵈ`, when evaluated at `α₂,...,αᵥ`, has a prime number factor which is not contained in the evaluations of the other factors (**identifying primes**)

核心匹配逻辑：

> "By using the identifying primes of condition (3) and checking which of the multivariate factors of the leading coefficient aᵈ(x₂,...,xᵥ), when evaluated, is divisible by each identifying prime, we can attach the correct multivariate leading coefficient to each factor u₁(x₁),...,uₙ(x₁)."

GCL 描述的是**高层启发式**。具体的 power splitting 由"divisibility matching"自然实现——条件 (3) 保证各因子可通过素因子区分，然后用整除检查分配。

### 3.2 SymPy `dmp_zz_wang_lead_coeffs`

```python
for h in H:                              # 外层: 遍历每个单变量因子
    c = dmp_one(v, K)                    # σᵢ = 1
    d = dup_LC(h, K) * cs               # R = lc(hᵢ) * cs

    for i in reversed(range(len(E))):    # 内层: 逆序遍历 LC 不可约因子
        k, e, (t, _) = 0, E[i], T[i]
        while not (d % e):               # ← 提取 Eⱼ-adic valuation
            d, k = d // e, k + 1
        if k != 0:
            c = dmp_mul(c, dmp_pow(t, k, v, K), v, K)   # σᵢ *= lⱼ^k
```

**注意**：SymPy 忽略 `eⱼ`（factorize 返回的幂次），由 valuation 自行确定 `k`。

### 3.3 FLINT `fmpz_mpoly_factor_lcc_wang`

```c
for (j = 0; j < r; j++) {                       // 外层: 每个单变量因子
    fmpz_mul(R, lc(Auf[j]), Auc);               // R = lc(uⱼ) * cs

    for (i = lcAfac->num - 1; i >= 0; i--) {    // 内层: 逆序遍历 LC 因子
        k = fmpz_remove(R, R, Q);               // ← 提取 Eⱼ-adic valuation
        lc_divs[j] *= poly[i]^k;                // σⱼ *= lⱼ^k
    }
}
```

### 3.4 CLPoly vs 参考实现的结构差异

| 方面 | CLPoly（当前，错误） | SymPy / FLINT / GCL |
|------|---------------------|----------------------|
| **外层遍历** | LC 因子 `(lⱼ, eⱼ)` | 单变量因子 `uᵢ` |
| **内层遍历** | 单变量因子 `uᵢ`，找唯一赢家 | LC 因子 `lⱼ`，各取所需 |
| **匹配方式** | `gcd(\|w[i]\|, \|lⱼ(α)\|^eⱼ)` 整块 | `while R % Eⱼ == 0` 逐幂提取 |
| **Power splitting** | 不可能 | 自然支持 |
| **幂次来源** | 使用 `eⱼ` | 由 valuation 自行确定 `k` |
| **遍历顺序** | 正序 | 逆序 |

### 3.5 用复现用例走一遍差异

设求值点 `y = 5`：

```
L = 18y³ → gamma=18, lc_factors={(y, 3)}, lj_val=5
w = [lc(u₀)*cs, lc(u₁)*cs, lc(u₂)*cs, lc(u₃)*cs]
  = [-10, -15, 1, 15]   (示例值，取决于单变量因式分解)
```

**CLPoly（GCD 匹配）**：
```
lj_pow = 5³ = 125
gcd(10, 125) = 5
gcd(15, 125) = 5
gcd(1, 125) = 1
gcd(15, 125) = 5
→ w[0], w[1], w[3] 的 gcd 都是 5 → ambiguous → FAIL
```

**SymPy/FLINT（valuation 提取）**：
```
u₀: R=10, 10%5=0 → R=2, k=1.  2%5≠0 → σ₀ *= y¹
u₁: R=15, 15%5=0 → R=3, k=1.  3%5≠0 → σ₁ *= y¹
u₂: R=1,  1%5≠0  → k=0.
u₃: R=15, 15%5=0 → R=3, k=1.  3%5≠0 → σ₃ *= y¹
→ y³ 拆分为 y¹·y¹·y¹ → SUCCESS
```

---

## 第四部分：修正方案

### 修改位置

`clpoly/polynomial_factorize_wang.hh`，`__wang_leading_coeff` 函数，**Step 3**（line 1400-1445）。

其余代码（Step 1/2/4/5、`__select_eval_point`、MTSHL、Hensel lift）均不改。

> 注：coprime 检查放在 `__wang_leading_coeff` 内部（Step 3 开头），而非 `__select_eval_point`。
> 原因：γ 在 Step 2 才产生，cs 作为参数传入；在此处检查最自然，与 SymPy 架构一致
> （SymPy 的 `non_divisors` 也在 LC 分配函数内部调用，而非选点函数内部）。

### 修正后完整算法

```
__wang_leading_coeff(f, u₁...uᵣ, α, x₁, cs)

输入:
  f ∈ Z[x₁,...,xᵥ]          — 待因式分解的多变量多项式
  u₁,...,uᵣ ∈ Z[x₁]         — f(x₁, α₂,...,αᵥ) 的单变量因子
  α = {x₂→α₂,...,xᵥ→αᵥ}    — 求值点
  x₁                         — 主变量
  cs                          — uni_content = content(f(x₁, α))（单变量像的整数内容）

输出:
  success    — 是否成功
  f_scaled   — δ^(r-1) · f
  v₁,...,vᵣ  — 缩放后的单变量因子: vᵢ = (δ / lc(uᵢ)) · uᵢ
  σ₁,...,σᵣ  — LC 分配: σᵢ ∈ Z[x₂,...,xᵥ]
  τ₁,...,τᵣ  — LC 目标: τᵢ = (δ / σᵢ(α)) · σᵢ

算法:

Step 1: 提取 Leading Coefficient
  L ← lc(f, x₁) ∈ Z[x₂,...,xᵥ]
  δ ← L(α) ∈ Z
  若 δ = 0: FAIL（换点）

Step 2: 递归因式分解 L
  若 L 是常数:
    σ₁ ← L, σ₂ = ... = σᵣ ← 1
    跳到 Step 4
  否则:
    (γ, {(l₁,e₁),...,(lₘ,eₘ)}) ← factorize(L)
    其中 L = γ · l₁^e₁ · ... · lₘ^eₘ, γ ∈ Z, lⱼ ∈ Z[x₂,...,xᵥ] 不可约

  P2 验证:
    检查 |γ · ∏ lⱼ(α)^eⱼ| = |cs · ∏ lc(uᵢ)|
    不等则 FAIL（数据不一致）

Step 3: Valuation 提取（本次修正的核心）
  // 预计算 LC 因子求值值
  对 j = 1,...,m:
    Eⱼ ← |lⱼ(α)|
    若 Eⱼ ≤ 1: FAIL（无法通过整除区分，换点）

  // ★ 互素前置检查（GCL 条件 3 完整版 / SymPy non_divisors / FLINT distilled bases）
  // Eⱼ 两两互素已由 __select_eval_point 保证
  // 此处额外检查 Eⱼ 与 cs·γ 互素，防止 content 污染导致 valuation 过度提取
  对 j = 1,...,m:
    若 gcd(Eⱼ, |cs · γ|) ≠ 1: FAIL（换点）

  // 初始化
  σᵢ ← 1      (i = 1,...,r)
  wᵢ ← |lc(uᵢ) · cs|   (i = 1,...,r)

  // 对每个单变量因子，独立提取各 LC 因子的幂次
  对 i = 1,...,r:
    R ← wᵢ
    对 j = m, m-1,...,1  (逆序):     // ← SymPy/FLINT 均逆序遍历
      kᵢⱼ ← 0
      while R mod Eⱼ = 0:            // 提取 Eⱼ-adic valuation
        R ← R / Eⱼ
        kᵢⱼ ← kᵢⱼ + 1
      σᵢ ← σᵢ · lⱼ^kᵢⱼ
    // R 中的剩余整数因子（来自 cs 和 lc(uᵢ) 中 coprime with all Eⱼ 的部分）
    // 不吸收到 σᵢ — Step 5 的 (δ/σᵢ(α)) 缩放会自动处理

  // ★ 守恒验证：Σᵢ kᵢⱼ = eⱼ
  // 当 Eⱼ 是合数时，while 循环提取的是 Eⱼ-整体 valuation，可能 < 真实 p-adic 分配
  // 此检查以 O(m) 代价捕获所有提取不精确的情形
  对 j = 1,...,m:
    若 Σᵢ kᵢⱼ ≠ eⱼ: FAIL（换点）

  // 吸收 γ
  σ₁ ← σ₁ · γ

Step 4: 缩放
  f_scaled ← δ^(r-1) · f
  对 i = 1,...,r:
    vᵢ ← (δ / lc(uᵢ)) · uᵢ        // 使 lc(vᵢ) = δ

Step 5: 计算 LC 目标
  对 i = 1,...,r:
    τᵢ ← (δ / σᵢ(α)) · σᵢ          // 使 τᵢ(α) = δ = lc(vᵢ)

  return SUCCESS
```

**数学性质**:

- **输入前提 1**: `__select_eval_point` 保证各 `Eⱼ` 两两互素（GCL 条件 3 的因子间互素部分）
- **输入前提 2**: Step 3 互素检查保证 `gcd(Eⱼ, cs·γ) = 1`（GCL 条件 3 的 content 互素部分）
- **σᵢ 只含多项式因子**: σᵢ = (γ if i=1) · ∏ⱼ lⱼ^kᵢⱼ，不含整数余数 R（R 由 Step 5 自动处理）
- **乘积守恒**: 守恒验证通过时 `∏ᵢ σᵢ = γ · ∏ⱼ lⱼ^eⱼ = L`
- **Step 5 整除性**: `δ/σᵢ(α) = γ^{[i≠1]} · ∏ⱼ lⱼ(α)^{eⱼ - kᵢⱼ}` ∈ Z（因 eⱼ - kᵢⱼ ≥ 0）
- **Step 5 保证**: `τᵢ(α) = δ`，与 `lc(vᵢ) = δ` 匹配 → MTSHL 的 LC 校正有正确目标
- **乘积不变式**: `∏ᵢ τᵢ = δ^(r-1) · L = lc(f_scaled, x₁)`
- **Power splitting**: 当 `eⱼ > 1` 时，各 `kᵢⱼ` 由 valuation 自动确定，自然拆分

**正确性证明**（`Σᵢ kᵢⱼ = eⱼ`，限 Eⱼ 为素数时严格成立）:

```
前提:
  (A) gcd(Eⱼ, Eₖ) = 1      对所有 j ≠ k      ← __select_eval_point
  (B) gcd(Eⱼ, cs · γ) = 1   对所有 j            ← Step 3 互素检查
  (C) |δ| = |γ| · ∏ⱼ Eⱼ^eⱼ                     ← L 的因式分解 + 求值
  (D) |δ| = cs · ∏ᵢ |lc(uᵢ)|                    ← P2 验证
  (E) Eⱼ 是素数                                   ← 本证明的额外假设

推导:
  由 (C)(D): cs · ∏ᵢ |lc(uᵢ)| = |γ| · ∏ⱼ Eⱼ^eⱼ

  对固定 j, 取两边的 p-adic valuation（p = Eⱼ，由 (E) 保证 p 是素数）:
    v_p(cs) + Σᵢ v_p(|lc(uᵢ)|) = v_p(|γ|) + eⱼ · v_p(Eⱼ) + Σ_{k≠j} eₖ · v_p(Eₖ)

  由 (B): v_p(cs) = 0, v_p(|γ|) = 0
  由 (A): v_p(Eₖ) = 0 对 k ≠ j
  由 (E): v_p(Eⱼ) = 1

  因此: Σᵢ v_p(|lc(uᵢ)|) = eⱼ

  又 wᵢ = |lc(uᵢ)| · cs, 且 v_p(cs) = 0（由 B），故:
    kᵢⱼ = v_p(wᵢ) = v_p(|lc(uᵢ)|)       （p 是素数，v_p 即 while 循环的提取结果）

  结论: Σᵢ kᵢⱼ = eⱼ  ∎
```

**合数 Eⱼ 的情形**:

当 Eⱼ 是合数（如 Eⱼ = 6 = 2·3）时，`while R mod Eⱼ = 0` 提取的是
`v_{Eⱼ}(R) = min_{p|Eⱼ} ⌊v_p(R) / v_p(Eⱼ)⌋`，不再等于任何单个 p-adic valuation。
若 lc(uᵢ) 中各素因子分布不均（如 v₂ = 2, v₃ = 0），则 `v_6(lc(uᵢ)) = 0`，丢失分配。

此问题 **SymPy 和 FLINT 也存在**（均使用相同的 while 循环 / `fmpz_remove`）。
三个实现都依赖后续的 Hensel lift 失败 → 换点重试来兜底。

**守恒验证** `Σᵢ kᵢⱼ = eⱼ` 在合数情形下**不会**通过，因此直接 FAIL 换点，
不需要等到 Hensel lift 才发现错误。这是对 SymPy/FLINT 的改进。

### 修改内容

将 Step 3 从"外层遍历 LC 因子 → 内层找唯一赢家（GCD 匹配）"改为"外层遍历单变量因子 → 内层逐个提取 LC 因子的 valuation"。

**修改前**（line 1400-1445）：

```cpp
// GCD 匹配 LC 因子 (SymPy dmp_zz_wang_lead_coeffs 风格)
// 对每个 LC 因子 lⱼ^eⱼ, 找 gcd(|w[i]|, |lⱼ(α)|^eⱼ) 最大的唯一 i
for (auto& [lj, ej] : lc_factors)
{
    auto lj_eval = assign(lj, eval_point);
    ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
    if (lj_val == 0 || abs(lj_val) == ZZ(1)) return result;

    ZZ lj_pow(1);
    for (uint64_t e = 0; e < ej; ++e)
        lj_pow = lj_pow * abs(lj_val);

    // 找唯一最佳匹配
    size_t best_i = r;
    ZZ best_g(0);
    bool ambiguous = false;
    for (size_t i = 0; i < r; ++i)
    {
        ZZ g = gcd(abs(w[i]), lj_pow);
        if (g > best_g)
        {
            best_g = g;
            best_i = i;
            ambiguous = false;
        }
        else if (g == best_g && g > ZZ(1))
        {
            ambiguous = true;
        }
    }
    if (best_i >= r || ambiguous || best_g <= ZZ(1))
        return result;

    // σ[best_i] *= lⱼ^eⱼ
    Poly lj_power = lj;
    for (uint64_t e = 1; e < ej; ++e)
    {
        lj_power = lj_power * lj;
        lj_power.normalization();
    }
    sigma[best_i] = sigma[best_i] * lj_power;
    sigma[best_i].normalization();

    // 从 w 中移除已匹配的部分
    w[best_i] /= best_g;
}
```

**修改后**：

```cpp
// Valuation 提取 (GCL §8.7 / SymPy dmp_zz_wang_lead_coeffs / FLINT lcc_wang 风格)
// 预计算各 LC 因子的求值值
std::vector<ZZ> lc_evals(lc_factors.size());
for (size_t j = 0; j < lc_factors.size(); ++j)
{
    auto lj_eval = assign(lc_factors[j].first, eval_point);
    ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
    if (lj_val == 0 || abs(lj_val) <= ZZ(1)) return result;
    lc_evals[j] = abs(lj_val);
}

// 互素前置检查: gcd(Eⱼ, cs·γ) == 1
// 防止 content/γ 污染导致 valuation 过度提取 (SymPy non_divisors / FLINT distilled bases)
ZZ cs_gamma = abs(uni_content * gamma);
for (size_t j = 0; j < lc_evals.size(); ++j)
{
    if (gcd(lc_evals[j], cs_gamma) != ZZ(1))
        return result;  // Eⱼ 与 cs·γ 共享素因子，换点
}

// 对每个单变量因子，逐个提取各 LC 因子的 p-adic valuation
for (size_t i = 0; i < r; ++i)
{
    ZZ R = abs(w[i]);

    for (int j = (int)lc_factors.size() - 1; j >= 0; --j)  // 逆序 (同 SymPy/FLINT)
    {
        ZZ Ej = lc_evals[j];
        int k = 0;
        while (R % Ej == ZZ(0))
        {
            R = R / Ej;
            k++;
        }
        if (k > 0)
        {
            Poly lj_power = lc_factors[j].first;
            for (int e = 1; e < k; ++e)
            {
                lj_power = lj_power * lc_factors[j].first;
                lj_power.normalization();
            }
            sigma[i] = sigma[i] * lj_power;
            sigma[i].normalization();
        }
    }
    // R 中的剩余整数因子不吸收到 sigma — Step 5 的 (δ/σᵢ(α)) 自动处理
}

// 守恒验证: Σᵢ kᵢⱼ = eⱼ（捕获合数 Eⱼ 导致的提取不精确）
for (size_t j = 0; j < lc_factors.size(); ++j)
{
    int total_k = 0;
    for (size_t i = 0; i < r; ++i)
        total_k += k_matrix[i][j];  // 需在提取循环中记录 kᵢⱼ
    if (total_k != (int)lc_factors[j].second)
        return result;  // 幂次不守恒，换点
}
```

> 注：实现时需用 `std::vector<std::vector<int>> k_matrix(r, std::vector<int>(m, 0))` 记录各 kᵢⱼ。

### 为什么这样修

1. **消除根因**：valuation 提取逐幂分配，`y³` 自然拆分为 `y¹·y¹·y¹`——power splitting 不再需要特殊处理
2. **与三个独立来源一致**：GCL 的 identifying primes + divisibility matching、SymPy 的 `while not (d % e)`、FLINT 的 `fmpz_remove` 本质相同
3. **不依赖唯一性假设**：GCD 匹配要求每个 LC 因子有唯一匹配的单变量因子；valuation 提取不需要
4. **互素前置检查消除 content 污染**：`gcd(Eⱼ, cs·γ) == 1` 保证 valuation 精确（参见完整算法的正确性证明）

### 正确性保证

**三层防线**：

| 层 | 检查 | 位置 | 防御目标 |
|----|------|------|----------|
| 1 | `gcd(Eⱼ, Eₖ) == 1` (j≠k) | `__select_eval_point` line 1301-1306 | LC 因子间交叉匹配 |
| 2 | `gcd(Eⱼ, cs·γ) == 1` | Step 3 互素检查（**新增**） | content/γ 污染 valuation |
| 3 | `Σᵢ kᵢⱼ == eⱼ` | Step 3 守恒验证（**新增**） | 合数 Eⱼ 导致的提取不精确 |

- 层 1+2 在 Eⱼ 为素数时**充分保证** `Σᵢ kᵢⱼ = eⱼ`（见正确性证明）
- 层 3 在 Eⱼ 为合数时**兜底**：提取不精确则立即 FAIL，不等 Hensel lift 才报错
- 三层合在一起，对任意 Eⱼ 都能正确检测分配失败

**R ≠ 1 不吸收**：R 代表 wᵢ 中 coprime with all Eⱼ 的整数部分（来自 cs 和 lc(uᵢ) 的整数因子）。
不应乘入 σᵢ——否则 `δ/σᵢ(α)` 不再是整数（审核发现的 bug，已修正）。
Step 5 的 `(δ/σᵢ(α))` 缩放自动将这些整数因子分配到 τᵢ 中。

### 修正后复现用例的执行过程

设求值点 `y = 5`, `cs = 1`：

```
L = 18y³, gamma=18, lc_factors={(y, 3)}
lc_evals = [5]
互素检查: gcd(5, |1·18|) = gcd(5, 18) = 1 ✓
w = [10, 15, 1, 15]   (= |lc(uᵢ) · cs|)

提取:
i=0: R=10.  10%5=0 → R=2, k=1.  2%5≠0.  σ[0] *= y¹.  (R=2, 不吸收)
i=1: R=15.  15%5=0 → R=3, k=1.  3%5≠0.  σ[1] *= y¹.  (R=3, 不吸收)
i=2: R=1.   1%5≠0 → k=0.                                (R=1)
i=3: R=15.  15%5=0 → R=3, k=1.  3%5≠0.  σ[3] *= y¹.  (R=3, 不吸收)

守恒验证: Σk = 1+1+0+1 = 3 = e₁ ✓

γ 吸收: σ[0] *= 18

σ = [18y, y, 1, y]
∏σᵢ = 18y³ = L ✓

Step 5: δ = L(5) = 18·125 = 2250
τ[0] = (2250 / (18·5)) · 18y = 25 · 18y = 450y     ✓ 整数
τ[1] = (2250 / 5)     · y   = 450 · y  = 450y       ✓ 整数
τ[2] = (2250 / 1)     · 1   = 2250                   ✓ 整数
τ[3] = (2250 / 5)     · y   = 450 · y  = 450y       ✓ 整数

∏τᵢ = 450y · 450y · 2250 · 450y = 450³·2250·y³ = 2250³·18·y³ = δ³·L ✓
```

### 不改动的部分

- **Step 1-2**（L 提取、递归因式分解）：正确，保持不变
- **P2 恒等式验证**（line 1383-1398）：正确，保持不变
- **Step 4**（δ^(r-1) 缩放）：正确，保持不变
- **Step 5**（τᵢ 计算）：正确，保持不变
- **`__select_eval_point`**：已有 Eⱼ 两两互素检查，保持不变（Eⱼ vs cs·γ 检查放在 Step 3 内）
- **MTSHL / Hensel lift**：不变
- **`__wang_core` 的主变量轮换、Phase 1/2 结构**：不变
- **eⱼ 降序排序**（line 1375-1376）：对 valuation 提取无意义（GCD 匹配遗留），可删可留，不影响正确性

---

## 附注：与 `wang-lc-fix-power-distribution.md` 的关系

`wang-lc-fix-power-distribution.md` 在早期提出了 per-factor valuation 提取，后被 GCD 匹配方案取代（`wang-fix-gamma-coprimality.md`）。本文档基于完整调研（GCL §8.7 原文 + SymPy + FLINT 源码），确认 valuation 提取是正确方法，GCD 匹配是根因。
