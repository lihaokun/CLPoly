# CLPoly 因式分解模块：完整设计手册

> **版本**：2026-03-09
> **状态**：活跃开发文档
> **归档**：旧文档已移至 `archive/` 子目录

---

## 目录

1. [模块总览与文件结构](#1-模块总览与文件结构)
2. [单变量因式分解管线 (M1-M4)](#2-单变量因式分解管线)
3. [多变量因式分解管线 (M5-M6)](#3-多变量因式分解管线)
4. [Wang LC 分配算法 — 数学分析](#4-wang-lc-分配算法--数学分析)
5. [MTSHL 稀疏 Hensel 提升](#5-mtshl-稀疏-hensel-提升)
6. [试除重组 (Trial Division Recombination)](#6-试除重组)
7. [历史修复分析与漏洞模式](#7-历史修复分析与漏洞模式)
8. [根本性问题诊断](#8-根本性问题诊断)
9. [已知 Bug 与回归测试](#9-已知-bug-与回归测试)
10. [参考文献](#10-参考文献)

---

## 1. 模块总览与文件结构

### 1.1 文件组织

| 文件 | 行数 | 职责 |
|------|------|------|
| `polynomial_factorize.hh` | ~240 | 公共 API 入口：`factorize(polynomial/upolynomial)` |
| `polynomial_factorize_zp.hh` | ~400 | Zp[x] 原语：SqFree, DDF, EDF, factor_Zp |
| `polynomial_factorize_univar.hh` | ~1640 | Z[x] 单变量：Hensel 提升, Zassenhaus 重组, 素数选择, van Hoeij LLL |
| `polynomial_factorize_wang.hh` | ~2555 | Z[x₁,...,xₙ] 多变量：Wang LC 校正, MTSHL, 试除重组 |

### 1.2 完整调用图

```
factorize(polynomial_<ZZ,lex>)               [polynomial_factorize.hh]
  ├── 单变量 dispatch:
  │   └── __factor_squarefree_primitive_ZZ   [polynomial_factorize_univar.hh]
  │         ├── __select_prime               选 Zp 分解因子最少的素数
  │         │     └── __factor_Zp            [polynomial_factorize_zp.hh]
  │         │           ├── __squarefree_Zp
  │         │           ├── __ddf_Zp
  │         │           └── __edf_Zp
  │         ├── __hensel_lift                二叉树 Hensel 提升
  │         └── __factor_recombine           Zassenhaus 子集枚举
  │               或 __factor_recombine_vanhoeij  LLL 格约化
  │
  └── 多变量 dispatch:
      └── __factor_multivar                  [polynomial_factorize_wang.hh]
            ├── squarefreefactorize          cont 提取 + Yun [polynomial_gcd.hh]
            ├── __extract_monomial_content   纯变量幂提取
            └── __wang_core                  Wang 算法核心
                  ├── __select_eval_point    条件 (a)(b)(b')(d) 验证
                  ├── factorize(f₀)          递归单变量分解
                  ├── __wang_leading_coeff   LC 因子分配 ★★★ 最高风险
                  ├── __mtshl_lift           MTSHL-d 提升 + p-adic 恢复
                  │     ├── __mtshl_step_j   逐变量 Taylor 提升
                  │     │     ├── __mtshl_zp_univar_mdp  单变量 MDP
                  │     │     ├── __mtshl_sparse_int     稀疏插值 MDP
                  │     │     ├── __mtshl_multi_bdp      二变量 BDP 回退
                  │     │     └── __mtshl_wmds           多变量 WMDS 回退
                  │     └── p-adic 提升循环   系数恢复到 Z
                  └── 试除重组               Zassenhaus 子集枚举
```

### 1.3 数据流总览

```
输入 f ∈ Z[x₁,...,xₙ]
  │
  ▼
squarefreefactorize
  ├── cont(f, x₁) → 递归 factorize（变量数 ↓）
  └── pp(f, x₁) = g₁^m₁ · g₂^m₂ · ...
        │
        ▼ 对每个 gₖ（本原无平方）
  __wang_core(gₖ)
  │
  ├── 选主变量 x₁, 选求值点 α
  ├── f₀ = gₖ(x₁, α) → factorize(f₀) = u₁·u₂·...·uᵣ
  │
  ├── __wang_leading_coeff
  │   ├── L = lc(gₖ, x₁), δ = L(α)
  │   ├── factorize(L) → 递归（变量数 ↓）
  │   ├── 分配: σ₁,...,σᵣ 使得 ∏σᵢ = L
  │   ├── 缩放: f_scaled = δ^(r-1)·gₖ, vᵢ = δ·ūᵢ
  │   └── τᵢ = (δ/σᵢ(α))·σᵢ (Hensel LC 目标)
  │
  ├── __mtshl_lift(f_scaled, vᵢ, τᵢ, α, x₁, p)
  │   ├── 逐变量 Zp 提升（MTSHL-d）
  │   └── 对称约化 + p-adic 恢复到 Z
  │   → G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]
  │
  └── 试除重组
      ├── normed[i] = pp(Gᵢ), 正首项
      ├── Zassenhaus 子集枚举: 从 s=1 开始
      │   ∏ normed[subset] | g_remaining ?
      └── 剩余 → 最后一个因子
```

---

## 2. 单变量因式分解管线

### 2.1 M1: Zp[x] 分解

| 步骤 | 函数 | 算法 | 复杂度 |
|------|------|------|--------|
| 无平方 | `__squarefree_Zp` | Yun + Frobenius 根 | O(n²) |
| 按度分组 | `__ddf_Zp` | h ← h^p mod f* | O(n · M(n) · log p) |
| 等度分裂 | `__edf_Zp` | Cantor-Zassenhaus | O(k · M(n) · log p^d) |

### 2.2 M2: Hensel 提升

二叉树结构，二次提升（精度每步翻倍）。Bézout 系数 s·g + t·h ≡ 1 (mod m)。

**Mignotte 界**：B = C(n, ⌊n/2⌋) · ‖f‖₂，提升直到 p^k > 2·|lc(f)|·B。

### 2.3 M3: 因子重组

- **Zassenhaus**：子集枚举，Gosper's hack，r ≤ 64
- **van Hoeij LLL** (P1a)：格约化，多项式时间，r 较大时使用
- 剪枝：LC 整除 + 常数项整除

### 2.4 M4: `__select_prime`

遍历前 5-20 个素数，选模因子数最少的。跳过 lc(f) ≡ 0 (mod p) 和非无平方的。

---

## 3. 多变量因式分解管线

### 3.1 `__factor_multivar` 入口

```cpp
// 输入: f ∈ Z[x₁,...,xₙ]
// 输出: factorization<Poly>

1. sqf ← squarefreefactorize(f)    // cont 提取 + Yun
2. 对每个 (gₖ, mₖ):
   - 常数/单变量 → factorize(gₖ)
   - 多变量:
     a. __extract_monomial_content(gₖ) → 纯变量幂因子
     b. __wang_core(gₖ_reduced)
3. 排序, 验证 (debug assert: ∏fᵢ^eᵢ = f)
```

**关键前置条件**：`squarefreefactorize` 正确提取 cont(f, x₁)。若 cont 提取不完整（例如遗漏 Z[y] 中的因子），则 `__wang_core` 收到的输入不满足"关于 x₁ 本原"的前置条件，后续 LC 分配和 Hensel 提升可能出错。

### 3.2 `__wang_core` 核心循环

```
交错轮换主变量策略:
for 每个主变量 x₁ (deg(g, x₁) > 1):
    for skip = 0, 1, ..., BATCH_SIZE-1:
        α ← __select_eval_point(g, x₁, skip)
        f₀ ← g(x₁, α)
        factors ← factorize(f₀)

        if |factors| ≤ 1:
            x₁ 方向不可约, 标记 dead, break

        lc_result ← __wang_leading_coeff(g, factors, α, x₁)
        if !lc_result.success: continue  // 换点

        mv_factors ← __mtshl_lift(...)
        if empty: continue  // MTSHL 失败, 换点

        verified ← 试除重组(g, mv_factors)
        if |verified| ≥ 2: return verified
```

**终止性保证**：
- 可约多项式至少存在一个好求值点使提升成功
- 不可约多项式的所有主变量最终被标记 dead
- BATCH_SIZE = 200 提供充足的重试空间

---

## 4. Wang LC 分配算法 — 数学分析

### 4.1 问题定义

设 g ∈ Z[x₁,...,xₙ] 本原无平方，L = lc(g, x₁) ∈ Z[x₂,...,xₙ]。
取值后 f₀ = g(x₁, α) 分解为 u₁·u₂·...·uᵣ。

**目标**：找 σ₁,...,σᵣ ∈ Z[x₂,...,xₙ] 使得 ∏σᵢ = L，
且 σᵢ 是 g 的第 i 个不可约因子的 lc(·, x₁)。

### 4.2 当前算法（Valuation 提取）

```
1. factorize(L) → γ · ∏ lⱼ^eⱼ    (lⱼ 不可约)
2. Eⱼ = |lⱼ(α)|
3. Non-divisor 检查: 每个 Eⱼ 剥离与 cs·γ 及之前 E 的共享素因子后 > 1
4. wᵢ = |lc(uᵢ)·cs|
5. 对每个 (lⱼ, eⱼ), 对每个 wᵢ:
   提取 kᵢⱼ = max{k : Eⱼ^k | wᵢ}
   σᵢ ← σᵢ · lⱼ^kᵢⱼ
6. 守恒验证: Σᵢ kᵢⱼ = eⱼ
7. γ 吸收到 σ₀
```

### 4.3 正确性条件（严格数学推导）

**定理**：Valuation 提取算法正确当且仅当以下条件同时满足：

**(C1) 各 Eⱼ 两两 coprime**
: 由 `__select_eval_point` 条件 (d) 保证。若 gcd(Eⱼ, Eₖ) > 1，则 Eⱼ 的幂整除 wᵢ 时无法区分来自 lⱼ 还是 lₖ。

**(C2) 各 Eⱼ 不含 cs·γ 的素因子**
: 由 Non-divisor 检查保证。若 Eⱼ 的某个素因子 q 也出现在 cs·γ 中，则 wᵢ = |lc(uᵢ)·cs| 中 q 的 valuation 包含来自 cs 的"噪声"，导致 kᵢⱼ 过大。

**(C3) 各 Eⱼ ≥ 2**
: 当前代码要求 `|lⱼ(α)| > 1`。若 Eⱼ = 1，则 Eⱼ^k | wᵢ 恒成立，k 无上界。

**(C4) 乘积守恒 Σᵢ kᵢⱼ = eⱼ**
: 验证所有 lⱼ 的幂次被完全分配。若不守恒，说明 (C1)-(C3) 的保证不足（例如合数 Eⱼ 导致幂次提取不精确）。

**C4 失败的根本原因**：当 Eⱼ 是合数且 Eⱼ 的某个素因子也出现在另一个 wᵢ 中时，valuation 提取可能不精确：`Eⱼ^k | wᵢ` 但真正的贡献只有 `lⱼ^(k-1)` 次方。

**示例**：
```
L = (y+2)(y+3)，α: y=1 → E₁=3, E₂=4
u₁ 的 lc = 12 = 3·4，u₂ 的 lc = 1
→ 对 E₁=3: k₁₁ = v₃(12) = 1 ✓
→ 对 E₂=4: k₁₂ = max{k: 4^k | 12} = 1 ✓
→ Σk = 1+1 = 1+1 ✓  (假设 e₁=e₂=1)
```

但若 E₁=6, E₂=3（不 coprime!），条件 (d) 会拒绝。

### 4.4 历史漏洞

| 版本 | Bug | 根本原因 | 修复 |
|------|-----|----------|------|
| `44b0e04` | LC 幂次分配错误 | 初始实现使用 argmax GCD 匹配，不精确 | 改为幂次提取 |
| `5e2a48e` | GCD 匹配三类失败 | argmax 在多个因子同 GCD 时歧义 | 重写为 valuation |
| `9c8bdd6` | Valuation 提取不精确 | 合数 Eⱼ 导致过度提取 | 加守恒验证 |
| `bf52300` | Coprime 条件过严 | 要求 gcd(Eⱼ, cs·γ)=1 拒绝太多点 | 改为 non-divisor 剥离 |
| **B1** | f2*f4 未分裂 | **待诊断** | — |

### 4.5 σᵢ → τᵢ 转换的正确性

```
τᵢ = (δ / σᵢ(α)) · σᵢ

验证:
  τᵢ(α) = (δ/σᵢ(α)) · σᵢ(α) = δ = lc(vᵢ)  ✓
  ∏τᵢ = ∏((δ/σᵢ(α))·σᵢ)
       = (δ^r / ∏σᵢ(α)) · ∏σᵢ
       = (δ^r / δ) · L
       = δ^(r-1) · L = lc(f_scaled, x₁)      ✓
```

**前置条件**：`δ % σᵢ(α) == 0`（精确整除）。当前代码有 assert 检查。若此条件失败，意味着 σᵢ 的分配有误。

---

## 5. MTSHL 稀疏 Hensel 提升

### 5.1 概述

MTSHL (Monagan-Tuncer Sparse Hensel Lifting, CASC 2016/2018) 替代经典 GCL Diophantine 方程求解。核心思想：在 Zp 上逐变量 Taylor 提升，用稀疏插值恢复因子系数。

### 5.2 流程

```
__mtshl_lift(f_scaled, v₁,...,vᵣ, τ₁,...,τᵣ, α, x₁, p):

阶段 A: v₁,...,vᵣ ∈ Z[x₁] → F₁,...,Fᵣ ∈ Zp[x₁]

阶段 B: 逐变量 j=2,...,n:
  __mtshl_step_j(aⱼ, F, lc_tau, xⱼ, αⱼ, x₁, ...)
  │
  ├── LC 校正: lc(Fᵢ, x₁) ← τᵢ|{xⱼ₊₁=αⱼ₊₁,...}
  ├── Taylor 循环 k=1,...,Dⱼ:
  │   ├── error ← aⱼ - ∏Fᵢ
  │   ├── cₖ ← Taylor 系数 of error at xⱼ=αⱼ, order k
  │   ├── MDP 求解:
  │   │   ├── j=2: 单变量 MDP (__mtshl_zp_univar_mdp)
  │   │   └── j≥3: sparse_int → multi_bdp → wmds (级联回退)
  │   └── Fᵢ ← Fᵢ + σₖᵢ · (xⱼ-αⱼ)^k

阶段 C: 系数恢复
  C.1: 对称约化 F → Z
  C.2: 若 p > 2B → 直接返回
  C.3: p-adic 提升循环 (最多 l_max=5 轮)
       误差 = f_scaled - ∏result
       MDP 求解修正项
       result[i] += symmetric_mod(σᵢ, p) · M
```

### 5.3 MDP 求解级联

```
j=2 (单变量):
  __mtshl_zp_univar_mdp:
  δᵢ = sᵢ · c mod vᵢ (偏分式分解)

j≥3 (多变量):
  尝试 1: __mtshl_sparse_int (Vandermonde + θ-array)
  尝试 2: 重试一次 sparse_int (随机性)
  回退:
    aux_vars=1: __mtshl_multi_bdp (二变量 BDP)
    aux_vars≥2: __mtshl_wmds (多变量递归 WMDS)
```

### 5.4 MTSHL 关键限制

1. **Zp 上操作**：MTSHL 在 Zp 上提升，丢失精确 Z 信息。阶段 C 需要额外 p-adic 提升恢复。
2. **非首一问题**：CASC 2018 假设 f 关于 x₁ 首一。非首一时 Taylor 循环不保证收敛（LC 部分无法被 MDP 捕获），依赖后续试除兜底。
3. **素数敏感**：MTSHL 素数 p 不能整除 lc(g, x₁)(α)，否则 Zp 约化丢失信息。
4. **骨架估计**：sparse_int 依赖 Theorem 1（支撑集单调递减），若初始骨架估计过小可能失败。

---

## 6. 试除重组

### 6.1 算法

在 `__wang_core` 中：

```cpp
// mv_factors 来自 __mtshl_lift
// normed[i] = pp(mv_factors[i]), 正首项
// g_remaining 初始为原始输入 g

for s = 1 to |mv_T|/2:
    for 每个大小为 s 的子集 S ⊂ mv_T:
        prod = ∏ normed[S]
        normalize_factor(prod)
        (q, rem) = pair_vec_div(g_remaining, prod)
        if rem == 0:
            verified.push(prod)
            g_remaining = q
            normalize_factor(g_remaining)
            s = 1; break

// 剩余部分作为最后一个因子
if g_remaining 非常数:
    verified.push(pp(g_remaining))

return verified if |verified| ≥ 2
```

### 6.2 正确性分析

**定理**（子集整除 → 不可约）：若 pp(∏_{i∈S} Hᵢ) 整除 g（Z[x₁,...,xₙ] 上精确），且 S 是最小整除子集，则 pp(∏_{i∈S} Hᵢ) 在 Z[x₁,...,xₙ] 上不可约。

**证明**：若 F = pp(∏_{i∈S} Hᵢ) 可约，F = A·B，则 A 和 B 各对应 Hᵢ 的真子集 S_A ⊊ S, S_B ⊊ S（因模求值后 Hᵢ 回到不可约 uᵢ）。故 pp(∏_{i∈S_A} Hᵢ) | g，|S_A| < |S|，与 S 最小矛盾。

**但是**：剩余因子 g_remaining（最后添加的那个）**没有最小性保证**。它是所有未被小子集提取的 Hensel 因子的乘积，可能是**可约的**。

### 6.3 ★ 剩余因子可约性 — 系统性漏洞

当前代码直接将 `g_remaining` 加为因子，**不做进一步递归分解**。这是一个结构性设计缺陷：

```
场景: r=3 Hensel 因子 H₁, H₂, H₃ 对应真因子 F₁, F₂, F₃
  s=1: pp(H₂) | g → 提取 F₂, g_remaining = F₁·F₃
  s=1: pp(H₁) 不整除 g_remaining = F₁·F₃?
       → 取决于 MTSHL 提升精度和 pp 归一化
  s=1: pp(H₃) 不整除 g_remaining?
       → 取决于同上
  → g_remaining = F₁·F₃ 作为单个因子添加 ← BUG
```

**这可能是 B1 bug 的根本原因之一**：MTSHL 提升产生了近似但不精确的 Hensel 因子，使得单个因子的试除失败，但乘积整除成立。

---

## 7. 历史修复分析与漏洞模式

### 7.1 修复时间线

| 提交 | 日期 | 修复内容 | 影响函数 |
|------|------|----------|----------|
| `5c867b1` | Phase 4 | Zassenhaus 双重 lc 乘法 | `__factor_recombine` |
| `44b0e04` | Phase 5 | Wang LC 幂次分配/单项式提取/gamma 互素 | `__wang_leading_coeff` |
| `5e2a48e` | Phase 5 | LC 从 GCD 匹配改为幂次提取 | `__wang_leading_coeff` |
| `e40858a` | Phase 5 | 正确性/鲁棒性审计 | 多处 |
| `e8bfa3f` | Phase 5 | Diophantine 基本情形模逆 | `__multivar_diophantine` |
| `a905765` | Phase 5 | 数学正确性审计 | 多处 |
| `9c8bdd6` | M6 | Valuation 提取不精确 + 守恒验证 | `__wang_leading_coeff` |
| `bf52300` | M6 | Non-divisor 替换 coprime 条件 | `__wang_leading_coeff` |
| `38589c9` | M6b | p-adic + Taylor 早退 + 健壮性 | `__mtshl_lift` |

### 7.2 漏洞模式分类

**模式 A: LC 分配错误** (4 次修复)
- 症状：σᵢ 分配不满足 ∏σᵢ = L 或分配到错误的因子
- 根因：整数 valuation 的组合爆炸 + 求值点质量不足
- 受影响：`__wang_leading_coeff`

**模式 B: Hensel 提升精度不足** (2 次修复)
- 症状：p-adic 恢复后因子系数不精确，试除失败
- 根因：非首一多项式导致 MTSHL Taylor 循环不收敛
- 受影响：`__mtshl_step_j`, `__mtshl_lift`

**模式 C: 试除验证不完整** (1 次修复)
- 症状：可约因子被当作不可约返回
- 根因：剩余因子无递归分解
- 受影响：`__wang_core` 试除重组

**模式 D: 缩放/归一化错误** (1 次修复)
- 症状：首项系数重复乘法
- 根因：lc 在提升和重组中的角色混淆
- 受影响：`__factor_recombine`

### 7.3 热点函数

```
__wang_leading_coeff:  4 次修复  ← 最高风险
__mtshl_lift:          2 次修复
__wang_core 试除:      1 次修复  ← 结构性缺陷
__factor_recombine:    1 次修复
```

---

## 8. 根本性问题诊断

### 8.0 文献与参考实现交叉验证

#### GCL §8.7 / Wang (1978)

GCL p.377 **明确承认**求值点可能导致过度分裂：

> "Note that the polynomial may factor into more factors in an image domain than in the
> original multivariate domain, and hence the issue of combining image factors arises."

Wang 的三个求值点条件（lc 非零、无平方、LC 因子可分辨）**不包含因子数最小化要求**。
算法依赖 Zassenhaus 子集试除来处理多余因子的重组。

#### FLINT `fmpz_mpoly_factor`

FLINT 同样**不要求**求值点因子数最小。其处理方式：

1. **`zassenhaus_prune`**：收集最多 3 个求值点的因子度数模式，用度数集合交集剪枝不可能的子集组合（而非选因子最少的点）
2. **精确 Hensel 提升**：在 Z\_{p^k} 上提升，Mignotte 界保证精度，因子系数在 Z 上精确。因此即使过度分裂，子集乘积也能通过试除正确重组
3. **Wang + Kaltofen 双路径**：LC 分配失败时有 Kaltofen 回退（二变量分解确定 LC）

#### Maple MTSHL (CASC 2016/2018, Tian Chen PhD §4.1)

Maple **从根源上避免**过度分裂，策略是**随机大整数求值点**：

1. **素数**：p = prevprime(2^62 - 1)，63-bit 机器素数
2. **求值点**：α = (α₂,...,αₙ) 从 [1, N̄-1]^{n-1} 中**均匀随机**采样，N̄ ~ 2^62
3. **概率保证**（Cohen bound）：Pr[非 Hilbert 点] < c̄(d) · N̄^{-1/2} · log N̄ ≈ 10^{-8}
4. **失败处理**：换新随机 α 重试（Las Vegas 算法），无需经典 Hensel 回退
5. 对每个提升因子**单独**试除，因此不存在"子集重组"步骤

**关键**：Maple 论文中示例求值点为 `α = [2908, 3830, 2798]`，不是小整数。
Hilbert 不可约性定理的概率保证**仅对大范围随机采样有效**，对小整数枚举无效。

#### ★ 关键差异：CLPoly (MTSHL) vs FLINT (经典 Hensel)

| | 经典 Hensel (FLINT/NTL) | MTSHL (CLPoly) |
|---|---|---|
| 提升域 | Z\_{p^k}（精确到 Mignotte 界） | Zp（单素数） |
| 过度分裂时 | 子集乘积仍有精确 Z 系数 → 试除成功 | 不存在的分解 → Zp 系数无 Z 对应 → 对称约化为大整数 |
| 重组方式 | Zassenhaus 子集枚举 | 单独试除 + 子集枚举 |
| 过度分裂的容忍度 | 高（子集重组精确） | **低**（垃圾系数阻止重组） |

**结论**：CLPoly 的 `__select_eval_point` 沿用了 Wang (1978) 的小整数枚举策略（从 bound=0 开始），
但提升核心已替换为 MTSHL。这两个组件**架构不兼容**：

- Wang 小整数枚举 + 经典 Z\_{p^k} Hensel → **兼容**（过度分裂通过子集试除处理）
- Maple 随机大整数 + MTSHL Zp → **兼容**（Hilbert 定理保证几乎不过度分裂）
- **CLPoly 小整数枚举 + MTSHL Zp → 不兼容**（小整数易触发过度分裂，MTSHL 无法容忍）

**修复方向**：将 `__select_eval_point` 改为 Maple 的随机大整数策略，或增加因子数验证。

### 8.1 核心问题：试除剩余因子不递归分解

`__wang_core` 的试除重组有一个**结构性设计缺陷**：

```cpp
// 当前代码 (__wang_core, line ~2410-2417):
if (!g_remaining.empty() && !is_number(g_remaining))
{
    auto h = g_remaining;
    normalize_factor(h);
    if (!is_number(h))
        verified.push_back({std::move(h), 1});  // ← 不递归分解!
}
```

**问题**：g_remaining 是所有未被小子集提取的 Hensel 因子的乘积。如果某些 Hensel 因子不够精确（MTSHL 提升精度问题或 p-adic 恢复不完整），单个因子的试除可能失败，但它们的乘积恰好整除 g。此时 g_remaining 是可约的。

**类比**：单变量 Zassenhaus (`__factor_recombine`) 同样将剩余添加为最后一个因子，但那里有关键的不同：
- 单变量 Hensel 提升在 Z_{p^k} 上精确（Mignotte 界保证），因子系数完全确定
- 多变量 MTSHL 在 Zp 上操作，系数恢复需要额外 p-adic 提升，精度有限

### 8.2 B1 Bug 分析（已确认）

**输入**：f = (y²-y-2)·(-2x²+xy-y)·(-2x²-2y²+3)·(x²-3y²+y)

**已确认的执行路径**（通过 `/tmp/debug_b1.cc` 诊断程序验证）：

1. squarefreefactorize 正确提取 cont(f, x) = y²-y-2 = (y-2)(y+1) ✓
2. pp(f, x) = f2·f3·f4 传入 __wang_core ✓
3. __wang_core 选主变量 x，选求值点 y=-1（skip=0，第一个候选点）
4. **关键**：f2·f3·f4 在 y=-1 处分解为 **5** 个单变量因子：
   - f2(-1) = -2x²-x+1 = -(2x-1)(x+1) → 2 个因子
   - f3(-1) = -2x²+1 → 1 个因子 (2x²-1)
   - f4(-1) = x²-4 = (x+2)(x-2) → 2 个因子
5. __wang_leading_coeff 成功（lc = 4 是常数）✓
6. __mtshl_lift 返回 5 个 Hensel 因子 ✓（"成功"但有问题）
7. **MTSHL 提升的 4 个线性因子系数全是垃圾值**（~10^18 量级）：
   ```
   H[0] = 4x + 2486369558885913747·y⁵ + ...  ← 垃圾！
   H[1] = 4x - 325948023030939647·y⁵ + ...   ← 垃圾！
   H[4] = 4x² + 4y² - 6                       ← 正确！(= 2·f3)
   ```
   原因：真正只有 3 个多变量因子，不存在 5 因子分解，
   Zp 上的提升方程无 Z 上的解，产生的 Zp 值经对称约化后是大整数。
8. 试除重组：
   - s=1: 只有 normed[4] = 2x²+2y²-3 (f3) 整除 g ✓
   - s=2: 所有 10 对都不整除（因为线性因子的系数是垃圾）✗
   - g_remaining = f2·f4 作为最后一个因子添加
9. **`verified.size() == 2 ≥ 2` → 立即返回，不再尝试其他求值点！**

**验证**：y=1, y=-2, y=2, y=-3 都给出 3 个单变量因子，MTSHL 提升全部正确，试除全部成功。但代码在 y=-1 处就停止了。

### 8.3 根本原因总结

**两个独立问题协同导致了 bug**：

**问题 A: 过早返回**
`__wang_core` 在 `verified.size() >= 2` 时立即返回，不验证剩余因子的不可约性。
一个"部分成功"的求值点阻止了后续更好的求值点被尝试。

**问题 B: 求值点因子数过多**
当求值点使某些真因子进一步分裂（如 f2 在 y=-1 处分裂为两个线性因子），
MTSHL 尝试提升过多因子，导致提升方程无解（Zp 上有解但 Z 上无解），
产生垃圾系数。这本身不是 bug（试除会拒绝），但与问题 A 结合就导致了错误结果。

### 8.4 修复方向

#### 方向 1: 剩余因子递归分解（兜底方案，简单正确）

```cpp
// 修改 __wang_core 的试除重组尾部:
if (!g_remaining.empty() && !is_number(g_remaining))
{
    auto h = g_remaining;
    normalize_factor(h);
    if (!is_number(h) && get_variables(h).size() >= 2)
    {
        // ★ 递归分解剩余因子
        auto sub = __wang_core(h);
        for (auto& [fi, ei] : sub)
            verified.push_back({std::move(fi), ei});
    }
    else if (!is_number(h))
        verified.push_back({std::move(h), 1});
}
```

**优点**：简单、正确、能兜底所有"试除不完整"的情况。
**缺点**：递归调用开销较大；可能掩盖上游问题（如求值点选择不优）。

#### 方向 2: 不完全分裂时继续尝试下一个求值点（推荐方案）

```cpp
// 修改 __wang_core 的返回条件:
// 旧: if (verified.size() >= 2) return verified;
// 新: 检查 g_remaining 是否不可约（即所有 Hensel 因子已被消耗）
if (!mv_T.empty())
{
    // 有 Hensel 因子未被试除消耗 → 此求值点不够好
    // → continue 到下一个 skip（即尝试下一个求值点）
    continue;
}
if (verified.size() >= 2)
    return verified;
```

**优点**：不引入额外递归，让外层循环自然找到更好的求值点。
**缺点**：在极端情况下可能耗尽 BATCH_SIZE 个求值点。

#### 方向 3: 结合方案（建议采用）

1. 优先方向 2：如果有 Hensel 因子未被消耗，continue 到下一个求值点
2. 兜底方向 1：如果方向 2 也失败（BATCH_SIZE 用完），递归分解 g_remaining

```cpp
if (!mv_T.empty() && skip < batch_end - 1)
    continue;  // 还有更多求值点可尝试
// 否则用当前结果，递归分解 g_remaining
```

#### 方向 4: 最小因子数选择（FLINT 单变量策略的多变量推广）

FLINT 单变量分解尝试 3 个素数并选因子数最少的。多变量版：

```
在 __wang_core 中，对每个主变量的前 N 个求值点（如 N=3），
记录因子数，选因子数最少的点作为 MTSHL 提升的起点。
```

**优点**：从根源减少过度分裂的概率。
**缺点**：额外 2 次单变量分解的开销（实际开销小于 MTSHL 提升）。

#### 方向 5: 随机大整数求值点（Maple 策略，治本）

将 `__select_eval_point` 从小整数枚举改为随机大整数采样：

```cpp
// 当前 (Wang 小整数枚举):
for (int bound = 0; ; ++bound)
    for each α ∈ [-bound, bound]^{n-1}: ...

// 改为 (Maple 随机采样):
std::mt19937_64 rng(seed);
uint64_t N_bar = p / 2;  // 采样范围上界
for (int attempt = 0; ; ++attempt)
    for each variable v:
        alpha[v] = ZZ(1 + rng() % (N_bar - 1));
    // 验证条件 (a)(b)(d) + 因子数校验
```

**优点**：从根源解决架构不兼容，与 Maple 对齐，Hilbert 概率保证生效。
**缺点**：求值点系数变大，多变量求值开销增加（但 MTSHL 提升在 Zp 上，不受影响）。
**注意**：需要保留 Wang LC 条件 (d)（LC 因子求值两两互素），Maple 也需要此条件。

#### 安全网分析：即使改用随机求值点，小概率事件能否被捕获？

**不能。** 当前代码对场景 2（部分因子正确+部分垃圾）无防御：

| 检查点 | 能否捕获？ | 原因 |
|--------|-----------|------|
| `__mtshl_lift` 返回空 | 否 | Zp 上总能"成功"提升 |
| 试除 `rem.empty()` | 部分 | 垃圾因子被拒绝，但真因子被提取 |
| `verified.size() >= 2` | **漏洞** | 1 个真因子 + 可约 g_remaining → 返回 |
| debug assert `∏fᵢ = f` | 否 | 可约因子的乘积仍等于 f |

因此方向 5（随机求值点）降低触发概率，但方向 1 或 2 是**必要的安全网**。

#### 建议的修复策略（分层防御）

1. **短期**：方向 2（mv_T 非空时 continue），修复 B1 的直接原因
2. **短期兜底**：方向 1（递归分解 g_remaining），对小概率事件的最后防线
3. **中期**：方向 5（随机大整数求值点），从根源解决架构不兼容 ★推荐
4. **中期备选**：方向 4（最小因子数选择），比方向 5 简单但不彻底
5. **长期**：评估是否需要 Kaltofen 回退（FLINT 双路径策略）

方向 5 + 方向 2 是最佳组合：方向 5 使过度分裂概率降至 ~10^{-8}，方向 2 在极端情况下仍能正确处理。

---

## 9. 已知 Bug 与回归测试

### 9.1 B1: 二变量 4 因子不完全分裂

**测试**：`test/test_factorize_multivar.cc` — `"factorize: B1 bivar 4-factor incomplete (KNOWN BUG)"`

```
f = (y²-y-2)(-2x²+xy-y)(-2x²-2y²+3)(x²-3y²+y)
期望: 5 个不可约因子 [(y-2), (y+1), (2x²-xy+y), (2x²+2y²-3), (x²-3y²+y)]
实际: 4 个因子 [(y-2), (y+1), (2x²+2y²-3), (2x⁴-x³y-6x²y²+3x²y+3xy³-xy²-3y³+y²)]
```

**状态**：KNOWN BUG，已加入回归测试

### 9.2 已修复的历史 Bug

| Bug ID | 描述 | 修复提交 | 回归测试 |
|--------|------|----------|----------|
| — | Zassenhaus 双重 lc 乘法 | `5c867b1` | `test_factorize` |
| — | Wang LC 幂次分配三类失败 | `44b0e04`, `5e2a48e` | `test_factorize_multivar` |
| — | Diophantine 基本情形除法 | `e8bfa3f` | `test_multivar_hensel` |
| — | Valuation 提取不精确 | `9c8bdd6` | `test_factorize_multivar` |
| — | Coprime 条件过严 | `bf52300` | `test_factorize_multivar` |

---

## 10. 参考文献

| 标记 | 文献 |
|------|------|
| **GCL** | Geddes, Czapor & Labahn, *Algorithms for Computer Algebra*, Kluwer 1992, §8 |
| **MCA** | von zur Gathen & Gerhard, *Modern Computer Algebra*, Cambridge 2013, §14-16 |
| **Wang** | Wang, "An Improved Multivariate Polynomial Factoring Algorithm", Math. Comp. 1978 |
| **CASC16** | Monagan & Tuncer, *CASC 2016* — MTSHL-d 算法 |
| **CASC18** | Monagan & Tuncer, *CASC 2018* — MTSHL p-adic 系数恢复 |
| **CZ** | Cantor & Zassenhaus, *Math. Comp.* 1981 — 有限域分解 |
| **Mignotte** | Mignotte, *Math. Comp.* 1974 — 因子系数界 |
| **van Hoeij** | van Hoeij, *J. Number Theory* 2002 — LLL 重组 |

---

## 附录 A: `__wang_leading_coeff` 数学正确性检查清单

每次修改 LC 分配代码时，必须验证以下所有条件：

- [ ] ∏σᵢ = L（精确等式，作为多项式）
- [ ] ∏σᵢ(α) = δ = L(α)
- [ ] 守恒验证 Σᵢ kᵢⱼ = eⱼ 对所有 j
- [ ] δ % σᵢ(α) == 0 对所有 i（τᵢ 的精确整除前置条件）
- [ ] τᵢ(α) = δ 对所有 i
- [ ] ∏τᵢ = δ^(r-1) · L
- [ ] ∏vᵢ = f_scaled(x₁, α)
- [ ] lc(vᵢ) = δ 对所有 i

## 附录 B: 试除重组正确性检查清单

- [ ] normed[i] = pp(mv_factors[i])，正首项
- [ ] pair_vec_div(g_remaining, prod) 在 Z[x₁,...,xₙ] 上精确除法（rem 真正为零）
- [ ] 提取因子后 g_remaining = g / ∏(已提取因子)，无缩放因子残留
- [ ] **剩余因子不可约性**：当前不保证 ← 结构性缺陷
