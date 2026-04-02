# MTSHL 算法正确性分析

> 日期：2026-03-03
> 基于 5 篇 Monagan-Tuncer 原文 PDF 精读 + GCL Algorithm 6.4 对照
> 目的：确认 MTSHL 的真正过程，厘清 LC 校正的归属，严格论证 CLPoly 实现的数学正确性

---

## 一、论文原文精读总结

### 1.1 五篇论文的算法一致性

| 论文 | 核心算法 | monic 假设 | LC 校正 | non-monic 说明 |
|------|---------|-----------|---------|---------------|
| CASC2016 | Alg 1 (MHL) + Alg 4 (SHL) | 明确 monic | 无 | "assuming monic...for simplicity" |
| ICMS2018 | Alg 1 (MHL) | 明确 monic | 无 | "refer to [3] for non-monic"（[3]=GCL） |
| CASC2018 | Alg 1 (MHL) + Alg 4 (MTSHL-d) | 明确 monic | 无 | "remains true for non-monic with slight modifications. Our implementation uses Wang's LCC" |
| MC2019 | （无新算法，综述） | — | — | "integrated with Wang's LCC" |
| JSC2020 | Alg 1 (MHL) + Alg 5 (MTSHL) | 明确 monic | 无 | "assume monic...not to complicate with LCC" |

**结论：所有 5 篇论文的 MTSHL 算法描述都不包含 LC 校正。**

### 1.2 MTSHL 算法（JSC2020 Algorithm 5，最终版）

```
Algorithm 5: j-th step of MTSHL for j > 1

输入:
  aj ∈ Zp[x1,...,xj]
  fj-1, gj-1 ∈ Zp[x1,...,xj-1]
  全部 monic in x1
  aj(xj=αj) = fj-1 · gj-1

输出:
  fj, gj ∈ Zp[x1,...,xj] such that aj = fj·gj, or FAIL

1: if #fj-1 > #gj-1 then interchange fj-1 with gj-1 end if
2: (σ0, τ0) ← (fj-1, gj-1)
3: (fj, gj) ← (fj-1, gj-1)
4: error ← aj - fj·gj; monomial ← 1
5: for i = 1, 2, 3, ... while error ≠ 0 and deg(fj,xj)+deg(gj,xj) < deg(aj,xj) do
6:    monomial ← monomial × (xj - αj)
7:    c ← Taylor coefficient of (xj-αj)^i of error at xj = αj
8:    if c ≠ 0 then
9:        Solve MDP: σi·gj-1 + τi·fj-1 = c for σi, τi
10:       σg ← σi-1  (Strong SHL assumption: Supp(σi) ⊆ Supp(σi-1))
11:       (σi, τi) ← SparseInt(gj-1, fj-1, c, σg)
12:       if FAIL then (σi, τi) ← MDSolver(fj-1, gj-1, c, p) end if
13:       if FAIL then restart with different α end if
14:       (fj, gj) ← (fj + σi×monomial, gj + τi×monomial)
15:       error ← aj - fj·gj
16:   end if
17: end for
18: if error = 0 return (fj, gj) end if
19: return FAIL
```

关键特征：
- **无 LC 校正步骤**
- 循环条件：`while error ≠ 0 and deg(fj,xj)+deg(gj,xj) < deg(aj,xj)`
- MDP 方程：`σi·gj-1 + τi·fj-1 = c`，约束 `degx1(σi) < degx1(fj-1)`
- SparseInt 使用 Strong SHL assumption：`Supp(σi) ⊆ Supp(σi-1)`

### 1.3 CASC2018 多因子版本（Algorithm 4, MTSHL-d）

```
Algorithm 4: j-th step of MTSHL-d for j > 1 (r factors)

1: for i from 1 to r do fj,i ← fj-1,i, σ0,i ← fj-1,i end do
2: error ← aj - ∏ fj,i
3: for k = 1, 2, 3, ... while error ≠ 0 and Σ deg(fj,i, xj) < deg(aj, xj) do
4:    ck ← Taylor coefficient of (xj-αj)^k of error
5:    if ck ≠ 0 then
6:        for i from 1 to r do σk,i ← σk-1,i end do  (Strong SHL assumption)
7:        (σk,1,...,σk,r) ← SparseInt(fj-1,i, ck, σk,i, i=1,...,r)
8:        if FAIL then restart with new α end if
9:        for i from 1 to r do fj,i ← fj,i + σk,i × (xj-αj)^k end do
10:       error ← aj - ∏ fj,i
11:   end if
12: end for
13: if error = 0 then return fj,1,...,fj,r else return FAIL end if
```

与 Algorithm 5 结构完全一致，仅从 2 因子推广到 r 因子。同样无 LC 校正。

### 1.4 循环上界的论文原文

| 论文 | 循环条件 |
|------|---------|
| CASC2016 Alg 1 | `for i from 1 to deg(aj, xj) while error ≠ 0` |
| CASC2016 Alg 4 | `for i from 1 to deg(aj, xj) while error ≠ 0` |
| CASC2018 Alg 1 | `for k from 1 while error ≠ 0 and Σdeg < deg(aj,xj)` |
| CASC2018 Alg 4 | `for k = 1,2,3,... while error ≠ 0 and Σdeg < deg(aj,xj)` |
| JSC2020 Alg 5 | `for i = 1,2,3,... while error ≠ 0 and deg(fj)+deg(gj) < deg(aj,xj)` |

CASC2016 直接使用 `deg(aj, xj)` 作为上界。后续论文改用度数累积条件 `Σ deg(fj,i, xj) < deg(aj, xj)`，本质等价——Taylor 展开的次数不超过 `deg(aj, xj)`。

---

## 二、GCL 的 non-monic Hensel lifting

### 2.1 GCL Algorithm 6.4（多变量 Hensel 提升算法）

论文引用 GCL (Algorithms for Computer Algebra, Geddes-Czapor-Labahn) Chapter 6 处理 non-monic 情况。

Algorithm 6.4 的输入：
```
(1) a(x1,...,xv) ∈ Z[x1,...,xv]，primitive in x1
(2) 求值同态 I = <x2-α2,...,xv-αv>
(3) 素数 p
(4) 精度 l (p^l/2 bounds all coefficients)
(5) u1,...,un ∈ Zp[x1]，pairwise coprime，a ≡ u1·...·un (mod <I,p>)
(6) lcU1,...,lcUn — 正确的多变量 leading coefficients  ← 关键输入
```

Algorithm 6.4 的核心结构（GCL pp.272-273）：
```
for j from 2 to v do
    // ★ LC 校正（Taylor 循环之前，执行一次）
    for m from 1 to n do
        if lcU_m ≠ 1 then
            coef ← evaluate lcU_m at current point mod p^l
            replace lcoeff(U_m, x1) with coef

    // Taylor 循环
    error ← Aj - product(Ui)
    for k from 1 to deg(Aj, xj) while error ≠ 0 do
        c ← Taylor coeff of (xj-αj)^k of error
        if c ≠ 0 then
            ΔU ← MultivariateDiophant(U, c, I, maxdeg, p, l)
            U ← U + ΔU × (xj-αj)^k
            error ← Aj - product(Ui)
    end for
end for
```

### 2.2 关键发现：LC 校正只在 Taylor 循环前执行一次

GCL Algorithm 6.4 的 LC 校正在 Taylor 循环（`for k from 1 ...`）**之前**执行，不是每步都做。

数学原因：MDP 方程 `Σ σk,i · bi = ck` 有约束 `degx1(σk,i) < degx1(fi)`。这意味着 σk,i 不会修改 fi 的 leading coefficient（关于 x1）。因此初始设置正确后，Taylor 循环中的更新 `F[i] += σk,i · (xj-αj)^k` 始终保持 LC 正确。

### 2.3 MTSHL 论文与 GCL 的关系

论文明确说明了分工：

| 职责 | 由谁负责 | 在哪里执行 |
|------|---------|-----------|
| 因式分解 lc(a, x1) | Wang's LCC | MTSHL 之前 |
| 分配 LC 因子到各 ui | Wang's LCC | MTSHL 之前 |
| 计算 τi（LC 目标） | Wang's LCC | MTSHL 之前 |
| 缩放 f, ui 使 lc 匹配 | Wang's LCC | MTSHL 之前 |
| 替换 lcoeff(ui) 为正确值 | GCL Alg 6.4 | 每步 j 的 Taylor 循环前 |
| Taylor 循环 + MDP 求解 | MTSHL | Algorithm 5 |

CASC2018 脚注 1 (line 330-332)：
> "This argument also works for the non-monic case if the leading coefficients of u and w w.r.t. x1 do not vanish at (α2,...,αn) modulo p, conditions which we note are imposed by Wang's LCC."

---

## 三、CLPoly 完整实现与论文的逐步对照

### 3.1 CLPoly 完整流程概览

```
__wang_core (line 2218-2424)
  │
  ├─ 选主变量、MTSHL 素数选择
  ├─ __select_eval_point → 求值点 α，条件 (1)(2)(3)
  ├─ factorize(f(x₁, α)) → 单变量因子 u₁,...,uᵣ
  │
  ├─ __wang_leading_coeff (line 1329-1518)    ← Wang's LCC
  │   ├─ Step 1: L = lc(f, x₁), δ = L(α)
  │   ├─ Step 2: factorize(L) → γ·∏lⱼ^eⱼ
  │   ├─ Step 3: Valuation 提取 → σᵢ = (γ if i=0)·∏lⱼ^kᵢⱼ
  │   ├─ Step 4: f_scaled = δ^(r-1)·f, vᵢ = (δ/lc(uᵢ))·uᵢ
  │   └─ Step 5: τᵢ = (δ/σᵢ(α))·σᵢ
  │
  ├─ __multivar_hensel_lift (line 1980-2126)
  │   ├─ Phase 0: Bézout 链（扩展 GCD）→ sᵢ, denom
  │   ├─ Phase 1: 初始化多变量因子 G[i]
  │   ├─ Phase 2: 确定提升变量顺序 x₂,...,xᵥ
  │   └─ Phase 3: 逐变量循环 → __hensel_lift_one_var
  │       │
  │       └─ __mtshl_step_j (line 763-941)    ← MTSHL
  │           ├─ LC 校正: lcoeff(F[i]) ← τᵢ evaluated
  │           ├─ F_base[i] = F[i]|_{xⱼ=αⱼ}
  │           └─ Taylor 循环 k=1..Dⱼ:
  │               ├─ cₖ = Taylor coeff
  │               ├─ solve MDP (SparseInt / BDP / WMDS)
  │               ├─ F[i] += σₖ[i]·(xⱼ-αⱼ)^k
  │               └─ error = aⱼ - ∏F[i]
  │
  └─ Zassenhaus 子集重组 → 提取真因子
```

### 3.2 逐步差异对照

#### Step A: 求值点选择与单变量因式分解

| 方面 | 论文 / GCL | CLPoly | 一致？ |
|------|-----------|--------|--------|
| 条件 (1): L(α) ≠ 0 | GCL §8.7 | `__select_eval_point` line 1272 | ✅ |
| 条件 (2): a(x₁,α) squarefree | GCL §8.7 | `__select_eval_point` line 1275 | ✅ |
| 条件 (3): Eⱼ pairwise coprime | GCL §8.7 | `__select_eval_point` line 1301-1306 | ✅ |
| 单变量因式分解 | Z[x₁] | `factorize(f₀)` | ✅ |

#### Step B: Wang's LCC（`__wang_leading_coeff`）

| 方面 | GCL §8.7 / SymPy / FLINT | CLPoly (修正后) | 一致？ |
|------|--------------------------|----------------|--------|
| L = lc(f, x₁) | ✓ | line 1350-1351 | ✅ |
| δ = L(α) | ✓ | line 1357 | ✅ |
| factorize(L) | ✓ | line 1368-1370 | ✅ |
| P2 恒等式验证 | SymPy has it | line 1384-1398 | ✅ |
| LC 因子分配 | Valuation 提取 | Valuation 提取 (line 1421-1450) | ✅ (修正后) |
| 分配方式 | `while R % Eⱼ == 0` | `while (R % Ej == ZZ(0))` | ✅ |
| 遍历顺序 | 逆序 (SymPy/FLINT) | 逆序 `j = m-1,...,0` | ✅ |
| R 吸收 | SymPy/FLINT 均不吸收 | 不吸收 (line 1449) | ✅ |
| coprime 检查 | SymPy: non_divisors | line 1412-1419 | ✅ (CLPoly 改进) |
| 守恒验证 | SymPy/FLINT 无 | line 1452-1460 | ✅ (CLPoly 改进) |
| γ 吸收到 σ₀ | SymPy/FLINT 有 | line 1462-1464 | ✅ |
| 缩放 f | δ^(r-1)·f | line 1472-1478 | ✅ |
| 缩放 uᵢ | (δ/lc(uᵢ))·uᵢ | line 1480-1500 | ✅ |
| 计算 τᵢ | (δ/σᵢ(α))·σᵢ | line 1502-1514 | ✅ |

#### Step C: Bézout 链（扩展 GCD）

| 方面 | GCL Alg 6.3 | CLPoly | 一致？ |
|------|------------|--------|--------|
| 计算 sᵢ 使 Σsᵢ·V̂ᵢ = denom | EEAlift | line 1992-2025 | ✅ |
| denom ≠ ±1 时的处理 | 隐含在 mod p^l 中 | 模逆元 (line 1674-1686) | ✅ (等价) |

#### Step D: MTSHL Taylor 循环（`__mtshl_step_j`）

| 方面 | 论文 Alg 5 / GCL Alg 6.4 | CLPoly | 一致？ |
|------|-------------------------|--------|--------|
| LC 校正（循环前） | GCL: 有 | line 800-802 | ✅ |
| 基底因子 | fⱼ₋₁,ᵢ | F_base[i] = F[i]\|_{xⱼ=αⱼ} (line 806-809) | ✅ |
| 循环上界 | deg(aⱼ, xⱼ) | D_j = degree(aj, xj) (line 864) | ✅ |
| Taylor 系数 | coeff of (xⱼ-αⱼ)^k in error | `__taylor_coeff_zp` (line 872) | ✅ |
| MDP 方程 | Σσₖ,ᵢ·bᵢ = cₖ | j=2: univar MDP; j≥3: SparseInt/BDP/WMDS | ✅ |
| SHL assumption | Supp(σₖ) ⊆ Supp(σₖ₋₁) | forms[i] = Supp(σₖ₋₁,ᵢ) (line 924-925) | ✅ |
| 因子更新 | fⱼ,ᵢ += σₖ,ᵢ·(xⱼ-αⱼ)^k | F[i] += σₖ[i]·xk_pow (line 917-920) | ✅ |
| error 更新 | aⱼ - ∏fⱼ,ᵢ | aj - product_F() (line 934) | ✅ |
| LC 校正（循环内） | **无** | **已删除** | ✅ 一致 |
| 早期终止 | error = 0 | error.empty() → break (line 937) | ✅ |

#### Step E: 因子恢复

| 方面 | GCL §8.7 | CLPoly | 一致？ |
|------|---------|--------|--------|
| Symmetric mod | 必须 | `__multivar_hensel_lift` 后处理 | ✅ |
| Trial division | Zassenhaus 子集重组 | line 2316-2408, Gosper 枚举 | ✅ |
| pp 化 + 正规化 | 必须 | line 2362-2380 | ✅ |

### 3.3 原有差异（已消除）

CLPoly 原先在 Taylor 循环内每步都调用 `lc_correct`，论文和 GCL 均不做。
经命题 3 证明这是空操作后，已删除该冗余代码。当前实现与论文/GCL 完全一致。

### 3.4 CLPoly 的改进（相对于 SymPy/FLINT）

1. **coprime 前置检查**：`gcd(Eⱼ, cs·γ) = 1`（line 1412-1419）——SymPy/FLINT 无此检查
2. **守恒验证**：`Σᵢ kᵢⱼ = eⱼ`（line 1452-1460）——捕获合数 Eⱼ 的提取失败
3. **MDP 多级回退**：SparseInt → 重试 → BDP → WMDS（line 896-910）——比论文更健壮

---

## 四、严格数学正确性证明

### 4.1 定理陈述

**定理**：设 `f ∈ Z[x₁,...,xₙ]` 是 squarefree primitive 多项式，`f = f₁·...·fᵣ` 是不可约分解。
若 CLPoly 的 `__wang_core` 返回因子列表 `{g₁,...,gₛ}`，则 `f = ±∏gᵢ^eᵢ`（乘积还原）。

**证明思路**：分步证明每个子步骤保持正确性不变式。

### 4.2 Step B 正确性：Wang's LCC

**命题 1**（LC 分配正确性）：若 `__wang_leading_coeff` 返回 success，则：

```
(P1) ∏ᵢ σᵢ = L = lc(f, x₁)
(P2) τᵢ(α) = δ = L(α)  对所有 i
(P3) ∏ᵢ τᵢ = δ^(r-1) · L = lc(f_scaled, x₁)
(P4) lc(vᵢ, x₁) = δ    对所有 i
```

**证明**：

**(P1)** σᵢ = (γ if i=0) · ∏ⱼ lⱼ^kᵢⱼ（line 1422-1464）。
守恒验证保证 Σᵢ kᵢⱼ = eⱼ（line 1452-1460），因此：
```
∏ᵢ σᵢ = γ · ∏ⱼ lⱼ^(Σᵢ kᵢⱼ) = γ · ∏ⱼ lⱼ^eⱼ = L  ∎
```

**(P2)** τᵢ = (δ/σᵢ(α)) · σᵢ（line 1512）。
求值：τᵢ(α) = (δ/σᵢ(α)) · σᵢ(α) = δ  ∎

**(P3)**
```
∏ᵢ τᵢ = ∏ᵢ (δ/σᵢ(α)) · ∏ᵢ σᵢ
       = (δʳ / ∏ᵢ σᵢ(α)) · L          由 (P1)
       = (δʳ / δ) · L                   因为 ∏ᵢ σᵢ(α) = L(α) = δ
       = δ^(r-1) · L  ∎
```

**(P4)** vᵢ = (δ/lc(uᵢ)) · uᵢ（line 1486-1494），因此 lc(vᵢ) = (δ/lc(uᵢ)) · lc(uᵢ) = δ  ∎

**命题 2**（δ/σᵢ(α) 是整数）：

```
σᵢ(α) = (γ if i=0) · ∏ⱼ Eⱼ^kᵢⱼ

对 i=0: δ/σ₀(α) = (γ·∏ⱼ Eⱼ^eⱼ) / (γ·∏ⱼ Eⱼ^k₀ⱼ) = ∏ⱼ Eⱼ^(eⱼ-k₀ⱼ)
对 i≠0: δ/σᵢ(α) = (γ·∏ⱼ Eⱼ^eⱼ) / (∏ⱼ Eⱼ^kᵢⱼ) = γ · ∏ⱼ Eⱼ^(eⱼ-kᵢⱼ)

由守恒验证 Σᵢ kᵢⱼ = eⱼ，故 eⱼ - kᵢⱼ ≥ 0，商为整数。 ∎
```

### 4.3 Step D 正确性：MTSHL Taylor 循环

**命题 3**（循环内 LC 不变）：若在 Taylor 循环开始前 `lcoeff(F[i], x₁) = τᵢ(α_{j+1},...,αₙ)` 正确，则循环内每步更新 `F[i] += σₖ[i]·(xⱼ-αⱼ)^k` 不改变 `lcoeff(F[i], x₁)`。

**证明**：

MDP 方程 `Σ σₖ,ᵢ · bᵢ = cₖ` 的约束为 `degx₁(σₖ,ᵢ) < degx₁(fⱼ₋₁,ᵢ)`。

设 dᵢ = degx₁(F[i])。更新后：
```
F[i]_new = F[i] + σₖ[i] · (xⱼ - αⱼ)^k
```

由于 `degx₁(σₖ[i]) < dᵢ`，而 `(xⱼ - αⱼ)^k` 不含 x₁，故：
```
degx₁(σₖ[i] · (xⱼ - αⱼ)^k) = degx₁(σₖ[i]) < dᵢ = degx₁(F[i])
```

因此 F[i]_new 的 x₁ 最高次项仍为 F[i] 的最高次项，`lcoeff(F[i], x₁)` 不变。 ∎

**推论**：CLPoly line 931 的循环内 `lc_correct` 在 MDP 约束满足时是空操作。

**命题 4**（Taylor 循环在 Dⱼ 步内收敛）：

设 `aⱼ = ∏ fⱼ,ᵢ` 是真分解。将 fⱼ,ᵢ 做 (xⱼ-αⱼ)-adic 展开：
```
fⱼ,ᵢ = Σ_{k=0}^{dᵢ} σₖ,ᵢ · (xⱼ - αⱼ)^k
```
其中 dᵢ = degxⱼ(fⱼ,ᵢ)。

由 GCL Theorem 6.6，Hensel 提升的唯一性保证：若 fⱼ₋₁,ᵢ 两两互素（mod I），则 Taylor 循环在 k = max(dᵢ) ≤ degxⱼ(aⱼ) = Dⱼ 步内精确恢复所有 σₖ,ᵢ。

CLPoly 使用 `for k=1 to D_j`（line 866），与此上界一致。 ∎

### 4.4 缩放的正确性：monic 等价

**命题 5**（δ 缩放使问题等效 monic）：

缩放后：
```
f_scaled = δ^(r-1) · f
vᵢ = (δ/lc(uᵢ)) · uᵢ
```

由 (P4)，lc(vᵢ, x₁) = δ 对所有 i。

MTSHL 的 Taylor 循环中，LC 校正将 lcoeff(F[i], x₁) 替换为 τᵢ evaluated。
由 (P2)，τᵢ(α) = δ 对所有 i。

因此所有因子的 leading coefficient（关于 x₁，在求值点处）都等于 δ。
除以 δ 后等效 monic，MTSHL 的正确性证明适用。 ∎

### 4.5 端到端正确性

**定理证明**：

1. `__wang_leading_coeff` 返回正确的 f_scaled, vᵢ, τᵢ（命题 1, 2）
2. Bézout 链正确构造 Diophantine 求解器（GCL Alg 6.3）
3. LC 校正在 Taylor 循环前设置正确的 leading coefficients（GCL Alg 6.4）
4. Taylor 循环在 Dⱼ 步内收敛到正确的多变量因子（命题 3, 4, 5）
5. Symmetric mod + trial division 从 mod p 因子恢复 Z 上的因子（GCL §8.7）
6. Zassenhaus 重组处理 spurious 因子合并

每一步的输出满足下一步的输入前提，链式保证端到端正确性。 ∎

---

## 五、结论

### 5.1 MTSHL 的真正过程

MTSHL 是一个 **纯 MDP 求解器**，不包含 LC 校正。它的输入是已经做过 LC 替换的多项式（等效 monic），输出是提升后的因子。

完整的因式分解流程中，MTSHL 的位置：

```
Wang 框架:
  1. 选主变量、求值点、因式分解 a(x₁, α)
  2. Wang's LCC: 分配 LC 因子，计算 τᵢ
  3. 缩放: f_scaled = δ^(r-1)·f, vᵢ = (δ/lc(uᵢ))·uᵢ
  ↓
  4. 对每个变量 j = 2,...,v:
     4a. LC 替换: lcoeff(Fᵢ, x₁) ← τᵢ evaluated (GCL Alg 6.4)
     4b. MTSHL Taylor 循环: 求解 MDP, 累加, 重算 error (论文 Alg 5)
  ↓
  5. Symmetric mod, trial division, 恢复因子
```

### 5.2 论文假设 monic 不是简化——是架构分工

论文说 "assume monic for simplicity" 容易误导。实际上：
- monic 假设不是简化表述，而是**架构分工**的体现
- Wang 框架的上游（LCC + δ 缩放）负责将问题转化为 MTSHL 可处理的形式
- MTSHL 只需要解决转化后的（等效 monic 的）问题

### 5.3 CLPoly 实现的正确性总结

| 检查项 | 结论 |
|--------|------|
| Wang's LCC (valuation 提取) | ✅ 与 GCL/SymPy/FLINT 一致（修正后） |
| δ 缩放 | ✅ 与 GCL 一致 |
| τᵢ 计算 | ✅ 数学上正确（命题 1-2） |
| LC 校正（循环前） | ✅ 与 GCL Alg 6.4 一致 |
| LC 校正（循环内） | ✅ 已删除冗余代码，与论文一致 |
| Taylor 循环上界 | ✅ D_j = deg(aj, xj)，与论文一致（命题 4） |
| MDP 求解 | ✅ SparseInt + 多级回退 |
| Error 更新 | ✅ 全量重算 aj - ∏F[i] |
| 因子恢复 | ✅ Symmetric mod + trial division |

**CLPoly 的 Wang→MTSHL 流程在数学上是正确的，且与论文/GCL 完全一致（原有的循环内冗余 LC 校正已删除）。**

---

## 附录：原文引用

### CASC2016 (Algorithm 4, SHL)
- "aj, fj-1, gj-1 are monic in x1"
- "for i from 1 to deg(aj, xj) while error ≠ 0 do"

### CASC2018 (Algorithm 4, MTSHL-d)
- "aj, fj-1,i are monic in x1"
- "for k = 1, 2, 3, ... while error ≠ 0 and Σ deg(fj,i, xj) < deg(aj, xj) do"
- 脚注: "This argument also works for the non-monic case if the leading coefficients of u and w w.r.t. x1 do not vanish at (α2,...,αn) modulo p, conditions which we note are imposed by Wang's LCC."

### JSC2020 (Algorithm 5, MTSHL)
- "aj, fj-1, gj-1 are monic in x1"
- "for i = 1, 2, 3, ... while error ≠ 0 and deg(fj,xj)+deg(gj,xj) < deg(aj,xj) do"

### GCL Algorithm 6.4 (pp.272-273)
- 输入包含 "lcU, a list of the n correct multivariate leading coefficients"
- LC 替换在 Taylor 循环之前: "U ← updated list U with lcoeff(U_m, x1) replaced by coef"
- Taylor 循环内无 LC 替换
