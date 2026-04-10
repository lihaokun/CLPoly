# 调研报告：Wang LCC（Leading Coefficient Correction）算法与 Power Splitting

> 调研日期：2026-03-03
> 目标：修复 CLPoly 多变量因式分解中 LC 分配的 power splitting 缺陷
> 背景：MTSHL Taylor 循环对非 monic 多项式不收敛，根因是 LC 分配不完整
> 一手源码（已直接阅读）：
>   - SymPy `sympy/polys/factortools.py`：`dmp_zz_wang_lead_coeffs`、`dmp_zz_wang_non_divisors`
>   - FLINT `src/fmpz_mpoly_factor/lcc_wang.c`：`fmpz_mpoly_factor_lcc_wang`
>   - FLINT `src/fmpz_mpoly_factor/lcc_kaltofen.c`：Kaltofen 回退方案
>   - FLINT `src/fmpz_mpoly_factor/irred_wang.c`：调度逻辑
>   - CLPoly `clpoly/polynomial_factorize_wang.hh`：`__wang_leading_coeff`（line 1329-1503）
> 参考论文：
>   - Wang 1978, "An Improved Multivariate Polynomial Factoring Algorithm", Math. Comp. 32:1215-1231 — LCC 原始出处（扫描 PDF，未提取全文）
>   - GCL (Geddes/Czapor/Labahn), "Algorithms for Computer Algebra", 1992, §8.7 pp.377-378 — **已读**
>   - Monagan & Tuncer, CASC 2016/2018 — "We use Wang's LCC"
>   - Lee 2013, "Factorization of multivariate polynomials", PhD thesis, TU Kaiserslautern — Algorithm 6.6（GCL Algorithm 6.3 改编）

---

## 1. 问题描述

### 1.1 触发条件

CLPoly 的 `crosscheck_flint_mpoly_factor_random_bivar_4fac` 测试有 ~1.8% 的失败率。
复现用例：

```
f₁ = -2xy + y² - 1     lc(f₁, x) = -2y
f₂ = -3xy - 2y + 1     lc(f₂, x) = -3y
f₃ = x² + xy + 2x      lc(f₃, x) = 1
f₄ = 3xy - 3y - 2       lc(f₄, x) = 3y
f = f₁ · f₂ · f₃ · f₄
```

`lc(f, x) = (-2y)(-3y)(1)(3y) = 18y³`

### 1.2 根因

`__wang_leading_coeff` 使用 GCD 匹配，将 `y³` **整块**分配给一个因子：

```
factorize(18y³) → content=18, factors={(y, 3)}
GCD 匹配: σ[best] = y³, 其余 σ[i] = 1
正确分布: σ₁ = y, σ₂ = y, σ₃ = 1, σ₄ = y
```

导致 τᵢ（LC 目标）全部是常数，MTSHL 的 LC 校正变成 no-op，Taylor 循环无法收敛。

### 1.3 触发频率低的原因

需要同时满足：
1. `lc(f, x₁)` 非常数
2. LC 的不可约因子有重复幂次（如 `y³`）
3. 该幂次需要拆分到多个提升因子

大多数随机多项式的 LC 要么是常数，要么不可约因子互不相同，GCD 匹配足够。

---

## 2. 论文依据

### 2.1 MTSHL 与 monic 假设

**CASC 2016** (p.2):
> "We will assume that a, f, g are monic in x₁ so as not to complicate the MHL algorithm with leading coefficient correction."

**CASC 2018** (p.3):
> "In the paper we assume the input polynomial a is monic in x₁ so as not to complicate the presentation with LCC. We note that what we explain remains true for the non-monic case with slight modifications. Our implementation uses Wang's LCC for the non-monic case."

**结论**：MTSHL 本身没有错。monic 假设由 Wang LCC 在上游保证。问题在 CLPoly 的 LCC 实现。

### 2.2 Taylor 循环上界

论文的循环条件同时有 error≠0 **和**度数上界：

- CASC 2016 Algorithm 1: `for i from 1 to deg(aⱼ, xⱼ) while error ≠ 0 do`
- CASC 2018 Algorithm 1: `for k from 1 while error ≠ 0 and Σᵢ deg(fⱼ,ᵢ, xⱼ) < deg(aⱼ, xⱼ) do`
- CASC 2018 Algorithm 4 (MTSHL-d): `for k = 1, 2, 3, ... while error ≠ 0 and Σᵢ deg(fⱼ,ᵢ, xⱼ) < deg(aⱼ, xⱼ) do`

CLPoly 用 `k ≤ D_j` 作为上界，与论文一致。

### 2.3 LC 校正

CASC 2016 (p.7) 提到 Kaltofen 在循环内做 LC 校正：
> "However in [3] leading coefficient correction is also done in the for loop."

CLPoly 的 in-loop LC 校正有文献依据。

### 2.4 GCL §8.7 (pp.377-378)：Wang's LC Predetermination 原文

GCL 的描述（整理自 OCR）：

> **输入**：`a(x₁,...,xᵥ)` 的单变量像因子 `u₁(x₁),...,uₙ(x₁) ∈ Z[x₁]`
>
> **步骤**：
> 1. 因式分解 LC：`aᵈ(x₂,...,xᵥ) = lc(a, x₁)`（递归调用多变量因式分解）
> 2. 选择求值点 `α₂,...,αᵥ ∈ Z` 满足三个条件：
>    - (1) `aᵈ(α₂,...,αᵥ) ≠ 0`
>    - (2) `a(x₁,α₂,...,αᵥ)` 无重因子
>    - **(3) each factor of the leading coefficient `aᵈ`, when evaluated at `α₂,...,αᵥ`,
>      has a prime number factor which is not contained in the evaluations of the other factors**
> 3. 因式分解每个像因子的 LC：`lc(u₁),...,lc(uₙ)`
> 4. 用条件 (3) 的 **identifying primes** 做整除匹配：检查哪个多变量 LC 因子（求值后）能被
>    该素数整除，从而将正确的多变量 LC 分配给每个因子

原文关键句：
> "By using the identifying primes of condition (3) and checking which of the
> multivariate factors of the leading coefficient aᵈ(x₂,...,xᵥ), when evaluated,
> is divisible by each identifying prime, we can attach the correct multivariate
> leading coefficient to each factor u₁(x₁),...,uₙ(x₁)."

**分析**：

GCL 的描述是**高层启发式**，没有给出 power splitting 的具体伪代码。核心思想是：

1. **条件 (3)** 保证各 LC 因子求值后有**唯一的素因子**（identifying prime）
2. 通过**整除检查**（divisibility）将 LC 因子匹配到单变量因子
3. 条件 (3) 本质上就是 SymPy `non_divisors` / FLINT `distilled bases` 的前置互素检查

GCL 没有描述当条件 (3) 不满足时如何处理（Wang 的原始算法是换求值点重试）。
SymPy/FLINT 的 valuation 提取是这个高层描述的具体实现：用 `while d % E == 0` 提取幂次，
本质上就是用 identifying primes 做整除匹配的精确化版本。

---

## 3. 正确算法：Valuation 提取

### 3.1 核心思想

**SymPy** 和 **FLINT** 都用相同的方法——对每个单变量因子独立提取各不可约 LC 因子的 p-adic valuation：

```
对每个单变量因子 uᵢ (i = 1..r):
    R = |lc(uᵢ)| * cs          // 数值跟踪值
    σᵢ = 1                     // 多项式 LC 分配（累积）

    对每个不可约因子 lⱼ（逆序遍历）:
        Eⱼ = |lⱼ(α)|           // 不可约因子在求值点的值
        k = 0
        while R % Eⱼ == 0:     // 提取 Eⱼ-adic valuation
            R = R / Eⱼ
            k += 1
        σᵢ *= lⱼ^k             // 分配 k 次幂（不是 eⱼ 次幂！）

    // 处理完后 R 应该是 ±1
```

### 3.2 与 CLPoly 当前做法的对比

| 方面 | CLPoly（当前，错误） | SymPy / FLINT（正确） |
|------|---------------------|----------------------|
| **遍历方向** | LC 因子为主：对每个 `lⱼ^eⱼ` 找匹配因子 | **因子为主**：对每个 `uᵢ` 提取各 LC 因子幂次 |
| **匹配方式** | `gcd(\|w[i]\|, \|lⱼ(α)\|^eⱼ)` → 整块分配 | `while R % Eⱼ == 0` → 逐幂提取 |
| **Power splitting** | **不可能**（整块分配给一个因子） | **自然支持**（每个因子独立提取所需幂次） |
| **复用乘数** | 使用 `eⱼ`（factorize 返回的幂次） | **忽略** `eⱼ`，由 valuation 自行确定 |
| **遍历顺序** | 正序 | 逆序（FLINT/SymPy 都是逆序） |

### 3.3 为什么 valuation 提取能自动 power split

以 `lc(f) = 18y³` 为例，evaluation `y = α`：

- `E = |α|`（y 在求值点的值）
- `lc(u₁) * cs` 包含 `α¹` → `k=1` → `σ₁ *= y¹`
- `lc(u₂) * cs` 包含 `α¹` → `k=1` → `σ₂ *= y¹`
- `lc(u₃) * cs` 不含 `α` → `k=0` → `σ₃` 不变
- `lc(u₄) * cs` 包含 `α¹` → `k=1` → `σ₄ *= y¹`

`y³` 被自然拆分为 `y¹ · y¹ · y¹`。

---

## 4. SymPy 实现细节

### 4.1 `dmp_zz_wang_lead_coeffs`

```python
def dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K):
    """Wang/EEZ: Compute correct leading coefficients."""
    C, J, v = [], [0]*len(E), u - 1

    # Phase 1: Valuation 提取（power splitting 核心）
    for h in H:                          # 外层：遍历每个单变量因子
        c = dmp_one(v, K)                # σᵢ = 1
        d = dup_LC(h, K) * cs            # R = lc(hᵢ) * cs

        for i in reversed(range(len(E))):  # 内层：逆序遍历 LC 不可约因子
            k, e, (t, _) = 0, E[i], T[i]  # 注意：multiplicity _ 被忽略！
            while not (d % e):             # 提取 E[i]-adic valuation
                d, k = d // e, k + 1
            if k != 0:
                c, J[i] = dmp_mul(c, dmp_pow(t, k, v, K), v, K), 1

        C.append(c)

    if not all(J):                       # 验证所有 LC 因子都被分配
        raise ExtraneousFactors

    # Phase 2: 调整 content cs
    CC, HH = [], []
    for c, h in zip(C, H):
        d = dmp_eval_tail(c, A, v, K)    # σᵢ(α)
        lc = dup_LC(h, K)                # lc(hᵢ)
        if K.is_one(cs):
            cc = lc // d
        else:
            g = K.gcd(lc, d)
            d, cc = d // g, lc // g
            h, cs = dup_mul_ground(h, d, K), cs // d
        c = dmp_mul_ground(c, cc, v, K)
        CC.append(c)
        HH.append(h)

    if K.is_one(cs):
        return f, HH, CC

    # Phase 3: 剩余 cs 分配给所有因子
    CCC, HHH = [], []
    for c, h in zip(CC, HH):
        CCC.append(dmp_mul_ground(c, cs, v, K))
        HHH.append(dmp_mul_ground(h, cs, 0, K))
    f = dmp_mul_ground(f, cs**(len(H) - 1), u, K)
    return f, HHH, CCC
```

### 4.2 `dmp_zz_wang_non_divisors`（前置互素检查）

```python
def dmp_zz_wang_non_divisors(E, cs, ct, K):
    """确保 E[i] 在去除公共因子后两两互素。"""
    result = [cs * ct]
    for q in E:
        q = abs(q)
        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r
            if K.is_one(q):
                return None       # 求值点无法区分此 LC 因子
        result.append(q)
    return result[1:]
```

此函数确保各 `Eⱼ` 在去除与 `cs * ct` 和已有 `E` 值的公共因子后仍非 1。
若任何 `E[i]` 被完全吸收（变为 1），说明求值点无法区分该 LC 因子，需要换点。

---

## 5. FLINT 实现细节

### 5.1 `fmpz_mpoly_factor_lcc_wang`（Wang 启发式 LCC）

```c
// Phase 1: 评估 LC 因子
for (j = 0; j < lcAfac->num; j++)
    lcAfaceval[j] = poly[j](alpha);    // Eⱼ = lⱼ(α)

// Phase 2: 构建互素 distilled bases
// 确保 d[0], d[1], ..., d[num] 两两互素
fmpz_mul(d + 0, Auc, lcAfac->constant);   // d[0] = cs * γ
for (i = 0; i < lcAfac->num; i++) {
    Q = |lcAfaceval[i]|;
    for (j = i; j >= 0; j--) {             // 去除与已有 d[j] 的公共因子
        R = d[j];
        while (!fmpz_is_one(R)) {
            R = gcd(R, Q);
            Q = Q / R;
            if (fmpz_is_one(Q)) { success = 0; goto cleanup; }
        }
    }
    d[i + 1] = Q;
}

// Phase 3: Valuation 提取（power splitting 核心）
for (j = 0; j < r; j++) {                  // 外层：每个单变量因子
    lc_divs[j] = 1;
    R = lc(Auf[j]) * Auc;                  // R = lc(uⱼ) * cs

    for (i = lcAfac->num - 1; i >= 0; i--) {  // 内层：逆序遍历 LC 因子
        Q = |lcAfaceval[i]|;
        k = fmpz_remove(R, R, Q);          // k = v_Q(R)，R = R / Q^k
        lc_divs[j] *= poly[i]^k;           // σⱼ *= lⱼ^k
    }
}

// Phase 4: 调整剩余整数 content
for (j = 0; j < r; j++) {
    Q = lc(Auf[j]) / gcd(lc(Auf[j]), dtilde[j]);
    lc_divs[j] *= Q;                       // 乘以剩余整数因子
}
```

### 5.2 `fmpz_mpoly_factor_lcc_kaltofen`（Kaltofen 确定性回退）

当 Wang 启发式失败时（通常因为求值点无法区分 LC 因子），FLINT 回退到 Kaltofen 方法：

1. 对每个非主变量 `xₖ`，求值其余变量得到二元多项式 `f(x₁, xₖ)`
2. 因式分解二元多项式，检查各因子的 LC
3. 通过多个二元因式分解交叉确认，**确定性**分配 LC 因子
4. 代价更高（需要多次二元因式分解），但不依赖求值点运气

### 5.3 调度逻辑 (`irred_wang.c`)

```c
if (lcAfac->num > 0) {
    success = 0;
    if (lcAfac_irred)   // 只有 LC 因子确认不可约时才用 Wang 启发式
        success = fmpz_mpoly_factor_lcc_wang(...);
    if (!success)        // 回退到 Kaltofen
        success = fmpz_mpoly_factor_lcc_kaltofen(...);
    if (!success)
        goto next_alpha; // 都失败 → 换求值点
}
```

---

## 6. CLPoly 当前实现分析

### 6.1 `__wang_leading_coeff`（line 1329-1503）

**Step 1-2**（正确）：提取 `L = lc(f, x₁)`，递归因式分解 `L`。

**Step 3**（有 bug）：GCD 匹配，整块分配：

```cpp
// line 1400-1445: GCD 匹配
for (auto& [lj, ej] : lc_factors) {
    ZZ lj_pow(1);
    for (uint64_t e = 0; e < ej; ++e)
        lj_pow = lj_pow * abs(lj_val);   // lj_pow = |Eⱼ|^eⱼ （整块！）

    // 找唯一最佳匹配
    size_t best_i = r;
    ZZ best_g(0);
    for (size_t i = 0; i < r; ++i) {
        ZZ g = gcd(abs(w[i]), lj_pow);
        if (g > best_g) { best_g = g; best_i = i; }
        else if (g == best_g && g > ZZ(1)) { ambiguous = true; }
    }

    sigma[best_i] *= lj^ej;              // 整块分配给一个因子！
    w[best_i] /= best_g;
}
```

**问题**：将 `lⱼ^eⱼ` 作为整体分配给 `gcd` 最大的单个因子，无法拆分幂次。

**Step 4-5**（正确）：缩放和 τᵢ 计算。

### 6.2 缺失功能

1. **Power splitting**：无法将 `y³` 拆分为 `y · y · y` 分配给多个因子
2. **互素前置检查**：没有验证各 `Eⱼ` 两两互素（SymPy 的 `non_divisors` / FLINT 的 distilled bases）
3. **Kaltofen 回退**：Wang 启发式失败时没有备选方案（直接 `return result` 失败）

---

## 7. 修复方案

### 7.1 核心改动：Step 3 替换为 valuation 提取

将 `__wang_leading_coeff` 的 Step 3（line 1400-1445）从 GCD 匹配改为 valuation 提取：

```cpp
// 新: valuation 提取（SymPy/FLINT 风格）
for (size_t i = 0; i < r; ++i) {             // 外层：每个单变量因子
    ZZ R = abs(w[i]);                         // R = |lc(uᵢ)| * cs

    for (int j = (int)lc_factors.size() - 1; j >= 0; --j) {  // 逆序遍历
        auto& [lj, ej] = lc_factors[j];
        ZZ E_j = abs(lc_factor_evals[j]);     // |lⱼ(α)|
        if (E_j <= ZZ(1)) { /* 失败，换点 */ return result; }

        int k = 0;
        while (R % E_j == ZZ(0)) {            // 提取 Eⱼ-adic valuation
            R = R / E_j;
            k++;
        }

        if (k > 0) {
            Poly lj_power = lj;
            for (int e = 1; e < k; ++e) {
                lj_power = lj_power * lj;
                lj_power.normalization();
            }
            sigma[i] = sigma[i] * lj_power;
            sigma[i].normalization();
        }
    }
    // 可选：验证 R == ±1（完全分配）
}
```

### 7.2 新增：互素前置检查

在 valuation 提取之前，验证各 `Eⱼ` 在去除公共因子后两两互素。若不互素，`return result`（失败，换点）。

### 7.3 可选：Kaltofen 回退

FLINT 在 Wang 启发式失败时回退到 Kaltofen 方法。CLPoly 目前换求值点重试即可，
Kaltofen 回退作为后续优化。

---

## 8. 验证计划

1. **复现用例**：`test_repro_bivar4fac` 的确定性测试应得到 5 个非平凡因子（当前得 4 个）
2. **随机测试**：200 次随机二元 4 因子测试全部通过
3. **回归测试**：268 个单元测试 + 258 个 crosscheck 测试全部通过
4. **稳定性测试**：crosscheck 连续 10 次全部通过（当前 master 有 ~1.8% 失败率）

---

## 9. 参考文献

### 一手来源

| 来源 | 内容 | 状态 |
|------|------|------|
| Wang 1978 | LCC 算法原始出处 | 扫描 PDF，未提取全文 |
| GCL §8.7 (pp.377-378) | Wang's LC Predetermination 描述 | **已读** — identifying primes + divisibility matching |
| SymPy `factortools.py` | `dmp_zz_wang_lead_coeffs` 完整实现 | **已读** |
| FLINT `lcc_wang.c` | `fmpz_mpoly_factor_lcc_wang` 完整实现 | **已读** |
| FLINT `lcc_kaltofen.c` | Kaltofen 回退方案 | **已读** |
| Lee 2013 thesis | Algorithm 6.6（GCL 6.3 改编） | 服务器不可达 |

### 论文引用

- Wang, P.S. "An Improved Multivariate Polynomial Factoring Algorithm", Mathematics of Computation, 32:1215–1231, 1978.
- Geddes, K.O., Czapor, S.R., Labahn, G. "Algorithms for Computer Algebra", Kluwer Academic, 1992, Chapter 6.
- Monagan, M., Tuncer, B. "Using Sparse Interpolation in Hensel Lifting", CASC 2016, LNCS 9890:381–400.
- Monagan, M., Tuncer, B. "Factoring Multivariate Polynomials with Many Factors and Huge Coefficients", CASC 2018, LNCS 11077:319–334.
- Lee, M.M.-D. "Factorization of multivariate polynomials", PhD thesis, TU Kaiserslautern, 2013.
- Kaltofen, E. "Sparse Hensel Lifting", EUROCAL '85, LNCS 204:4–17, 1985.
