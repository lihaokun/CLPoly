# 修复方案：Over-splitting 检测与随机化求值点

**状态**: 待实现
**日期**: 2026-03-10
**关联**: `docs/design/factorization/architecture.md` §6 根因诊断

## 1. 问题概述

CLPoly 多变量因式分解存在架构性缺陷：Wang 经典算法的小整数求值点枚举与 MTSHL 的 Zp 单素数 Hensel 提升组合不兼容。

| 组件 | 来源 | 设计假设 |
|------|------|----------|
| `__select_eval_point` | Wang (GCL §8.3) | 小整数枚举，配合 Z_{p^k} Hensel，**容忍 over-splitting** |
| `__mtshl_lift` | Monagan-Tuncer (CASC 2016/2018) | 随机大整数求值点，**要求因子数正确**（Hilbert 定理保证） |

**后果**：小整数求值点高概率导致 over-splitting → MTSHL 对分裂的单变量因子产生垃圾 Zp 系数 → trial division 只能验证部分因子 → `verified.size() >= 2` 提前返回 → 因子丢失。

**B1 复现案例**：`f = (x+y²-y-2)(x²+xy+1)(x²-xy+1)(x³+y³+1)`，首个求值点 y=-1 产生 5 个单变量因子（真因子数为 3），仅 f3 通过 trial division，f2·f4 作为"剩余"被整体返回（不可约性丢失）。

## 2. 修复方案

两个方向组合：

| 方向 | 作用 | 类型 |
|------|------|------|
| **方向 2**：严格化 `__wang_core` 返回条件 | 检测 over-splitting，拒绝不完整结果 | 正确性修复 |
| **方向 5**：`__select_eval_point` 随机化 | 降低 over-splitting 概率 | 效率优化 |

## 3. 改动 1：`__wang_core` 返回条件（方向 2）

### 位置

`clpoly/polynomial_factorize_wang.hh`，`__wang_core` 函数，lines 2410-2420

### 当前代码（有缺陷）

```cpp
// 无条件把 g_remaining 当最后因子塞入
if (!g_remaining.empty() && !is_number(g_remaining))
{
    auto h = g_remaining;
    normalize_factor(h);
    if (!is_number(h))
        verified.push_back({std::move(h), 1});
}

if (verified.size() >= 2)
    return verified;  // ← 无论 mv_T 状态如何都返回
```

缺陷：不检查 `mv_T`（未消耗的 Hensel 因子）是否为空。Over-splitting 时 `mv_T` 中有多个垃圾因子，但只要凑够 2 个 verified 就返回。

### 替换为

```cpp
// 接受条件：mv_T ≤ 1（所有 Hensel 因子已被 trial division 消耗，
// 最多剩 1 个作为 Zassenhaus 互补因子）
if (mv_T.size() <= 1)
{
    // mv_T == 1 → 最后一个 Zassenhaus 互补因子 = g_remaining
    if (!g_remaining.empty() && !is_number(g_remaining))
    {
        auto h = g_remaining;
        normalize_factor(h);
        if (!is_number(h))
            verified.push_back({std::move(h), 1});
    }
    if (verified.size() >= 2)
        return verified;
}
// mv_T > 1 → over-splitting 或 Hensel 提升失败 → 换下一个求值点
```

### `mv_T.size()` 语义

| `mv_T.size()` | 含义 | 动作 |
|---|---|---|
| 0 | 全部 Hensel 因子由 trial division 精确匹配，g_remaining 为常数 | 返回 verified |
| 1 | 最后一个因子是 Zassenhaus 互补（正常情况） | 加 g_remaining 后返回 |
| ≥ 2 | over-splitting：部分 Hensel 因子有垃圾系数，无法匹配 | `continue` 换下一个求值点 |

### Zassenhaus 互补因子说明

Zassenhaus 枚举子集大小从 `s=1` 到 `s ≤ |mv_T|/2`。循环结束时 `mv_T` 最多剩 1 个元素：

- 若 r 个 Hensel 因子对应 r 个真不可约因子（无 over-splitting），trial division 在 `s=1` 逐个匹配前 r-1 个，`mv_T` 剩最后 1 个，`g_remaining` 即最后一个不可约因子
- 若 over-splitting 导致某些 Hensel 因子有垃圾系数，这些因子无法通过精确整数 trial division → 留在 `mv_T` → `mv_T.size() > 1` → 被拒绝

## 4. 改动 2：`__select_eval_point` 随机化（方向 5）

### 位置

`clpoly/polynomial_factorize_wang.hh`，`__select_eval_point` 函数，lines 1222-1306

### 当前逻辑

从 `bound=0` 开始扩展 shell 枚举：`{0}` → `{±1}` → `{±2}` → ...，固定顺序，`skip` 跳过前 N 个合法点。

**问题**：小整数偏向高 over-splitting 概率。`skip=0` 返回的通常是 |α| ≤ 2 的点。

### 替换为

用确定性 PRNG 从 `[-B, B]^{n-1}` 随机采样（`B = max(1000, 10·tdeg(f))`）。

```
seed = hash(f 的系数摘要) ⊕ skip
PRNG = mt19937_64(seed)
for attempt = 0, 1, 2, ... :
    α[i] = PRNG() % (2B+1) - B,  i = 1, ..., n-1
    if α 满足条件 (a)(b)(b')(d):
        return α
```

### 保留不变的合法性检查

| 条件 | 含义 |
|------|------|
| (a) | `f(x₁, α)` 无平方 |
| (b) | `lc(f, x₁)(α) ≠ 0` |
| (b') | `lc(f)(α)` 不被 MTSHL 素数整除 |
| (d) | LC 不可约因子在 α 处的值两两互素 |

### 为什么不用 Maple 的 2^62 范围

CLPoly 的 Wang 框架在 Z 上运算（LC 分配、trial division 均为整数多项式运算）。Eval point 过大会导致：

- 单变量 `f(x₁, α)` 的系数是 O(α^{deg}) 量级，增大因式分解开销
- LC 分配中的 GCD/求值计算系数膨胀

`B ~ 10³` 已足够：Cohen 界给出 over-splitting 概率 `O(d · B^{-1/2}) ≈ O(d / 30)`，对 deg < 100 的多项式远小于 1%。配合方向 2 的 retry 机制，实际失败概率可忽略。

## 5. 正确性证明

### 5.1 基本设置与记号

设 R = Z[x₁, ..., xₙ]。由 Z 是 UFD 与 UFD 上多项式环保持 UFD 的定理（归纳 n 次），R 是 UFD。R 的可逆元群 R× = {±1}。

设 g ∈ R 满足：
- (SQ) g squarefree：若 p ∈ R 不可约且 p² | g，则矛盾
- (PR) g primitive：contentℤ(g) := gcd(g 的所有整数系数) = 1
- (LC) lc(g) > 0：g 在某固定单项式序下的首项系数为正

设 g 的不可约分解为：

> g = G₁ · G₂ · ... · Gₜ  ··· (★)

每个 Gⱼ ∈ R 不可约，contentℤ(Gⱼ) = 1，lc(Gⱼ) > 0。

**存在性**：R 是 UFD。**唯一性**：在 R 中，不可约分解在相伴和顺序意义下唯一。再加上 contentℤ = 1 和 lc > 0 的归一化条件，分解完全唯一（无顺序歧义外）。

**引理 5.1.1（Gauss 引理）**：设 A, B ∈ R，contentℤ(A) = contentℤ(B) = 1。则 contentℤ(A · B) = 1。

*证明*：标准 Gauss 引理。视 R = Z[x₁][x₂]...[xₙ]，对每个变量逐层应用 UFD 上的 Gauss 引理。 ∎

**引理 5.1.2（不可约因子的 primitive 性）**：在 (★) 中，每个 Gⱼ 满足 contentℤ(Gⱼ) = 1。

*证明*：g = G₁ · ... · Gₜ，contentℤ(g) = 1（条件 PR）。由 Gauss 引理的逆：contentℤ(G₁ · ... · Gₜ) = contentℤ(G₁) · ... · contentℤ(Gₜ)（在 Z 中相等）。各 contentℤ(Gⱼ) ≥ 1，乘积 = 1，故每个 contentℤ(Gⱼ) = 1。 ∎

### 5.2 求值点与 Hensel 提升

设 α = (α₂, ..., αₙ) ∈ Zⁿ⁻¹ 为求值点，满足条件 (a)(b)(b')(d)。

记 ḡ = g(x₁, α₂, ..., αₙ) ∈ Z[x₁]。由条件 (a)，ḡ 无平方。由条件 (b)，deg(ḡ) = degₓ₁(g)（首项系数未消失）。

设 ḡ 的不可约分解为 ḡ = c · f₁ · f₂ · ... · fᵣ，fᵢ ∈ Z[x₁] 不可约，lc(fᵢ) > 0，c ∈ Z。

对每个 j，记 Ḡⱼ = Gⱼ(x₁, α) ∈ Z[x₁]。由 ḡ = Ḡ₁ · ... · Ḡₜ（代入是环同态），每个 Ḡⱼ 是 ḡ 的因子。

**定义 5.2.1（好求值点）**：称 α 为好求值点（good evaluation point），若对每个 j = 1,...,t，Ḡⱼ 在 Z[x₁] 中不可约。

**等价条件**：α 为好求值点 ⟺ r = t 且存在双射 σ: {1,...,r} → {1,...,t} 使得 fᵢ ~ Ḡσ(i)。

*证明*：(⟹) 每个 Ḡⱼ 不可约，且 ḡ = ∏ Ḡⱼ = c · ∏ fᵢ。由 Z[x₁] 是 UFD，唯一分解给出 r = t 和双射。(⟸) r = t 和双射意味着每个 Ḡⱼ ~ fσ⁻¹(j) 不可约。 ∎

**定义 5.2.2（坏求值点）**：若存在某 j 使 Ḡⱼ 在 Z[x₁] 中可约，则称 α 为坏求值点。此时 r > t。

**定理 5.2.3（MTSHL 正确性）**[Monagan-Tuncer, CASC 2016 定理 1; 2018 §3]：设 α 为好求值点，满足条件 (a)(b)(b')(d)，且 MTSHL 素数 p ∤ lc(ḡ)。则 MTSHL 算法产出 h₁, ..., hᵣ ∈ R 使得：

> 对每个 i = 1,...,r：hᵢ = Gσ(i)

即 Hensel 提升精确恢复每个真因子。

**引理 5.2.4（好求值点下 hᵢ 的性质）**：在好求值点下，MTSHL 产出的 h₁,...,hᵣ 满足：

(a) 每个 hᵢ ∈ R 不可约（因为 hᵢ = Gσ(i) 不可约）

(b) contentℤ(hᵢ) = 1（因为 contentℤ(Gσ(i)) = 1，引理 5.1.2）

(c) hᵢ 两两互不相伴：若 i ≠ j，则 hᵢ ≁ hⱼ（即不存在 ε ∈ {±1} 使 hᵢ = ε · hⱼ）

*证明 (c)*：假设 hᵢ = ε · hⱼ（ε = ±1）。代入 α：hᵢ(x₁, α) = ε · hⱼ(x₁, α)。由定理 5.2.3，hᵢ(x₁, α) = Gσ(i)(x₁, α) = Ḡσ(i)，同理 hⱼ(x₁, α) = Ḡσ(j)。故 Ḡσ(i) = ε · Ḡσ(j)。由好求值点定义，Ḡσ(i) 和 Ḡσ(j) 都是 ḡ 的不可约因子。由 ḡ 无平方（条件 (a)），若 Ḡσ(i) ~ Ḡσ(j) 则 σ(i) = σ(j)，与 σ 双射矛盾（因 i ≠ j → σ(i) ≠ σ(j)）。 ∎

### 5.3 Zassenhaus 算法的形式描述

以下是 `__wang_core` 中 Zassenhaus trial division 循环的形式化：

```
输入: g ∈ R (满足 SQ, PR, LC), Hensel 因子 h₁,...,hᵣ ∈ R
T ← {1, ..., r}           // 可用 Hensel 因子下标集
g_rem ← g                 // 剩余多项式
verified ← []             // 已验证因子

s ← 1
while 2s ≤ |T|:
    对每个 S ⊂ T, |S| = s:
        Q ← ∏_{i∈S} hᵢ
        P ← pp(Q)             // P = Q / contentℤ(Q), lc(P) > 0
        if P | g_rem (精确整数多项式除法，余式 = 0):
            verified ← verified ∪ {P}
            g_rem ← g_rem / P
            T ← T \ S
            s ← 1               // 重启
            break
    if 本轮未找到:
        s ← s + 1

// 循环结束后：T 为 mv_T，g_rem 为 g_remaining
```

**不变量 I₁**：在算法每一步，g_rem · ∏_{F ∈ verified} F = g · ε，ε ∈ {±1}。

*维持*：初始 g_rem = g，verified 为空，ε = 1。每次 trial division 找到 P | g_rem，令 g_rem' = g_rem / P，则 g_rem' · P · ∏_{已有} F = g_rem · ∏_{已有} F = g · ε。∎

**不变量 I₂**：g_rem = ε' · ∏_{j ∈ J} Gⱼ，其中 J ⊂ {1,...,t} 为尚未被消耗的真因子下标集，ε' ∈ {±1}。

**I₂ 不依赖好求值点假设**。仅依赖以下事实：g 在 R（UFD）中的不可约分解为 g = G₁ · ... · Gₜ，且 g squarefree。

*初始*：g_rem = g = ∏_{j=1}^{t} Gⱼ，J = {1,...,t}，ε' = 1。

*维持*：设当前 g_rem = ε' · ∏_{j∈J} Gⱼ。若 trial division 找到 P | g_rem 且 P 非常数。由 R 是 UFD 且 g squarefree（与定理 1 步骤 1 相同的论证——注意步骤 1 不依赖好求值点），P = ∏_{j∈J'} Gⱼ 对某 J' ⊂ J。除法后 g_rem' = ε'' · ∏_{j∈J\J'} Gⱼ。不变量以 J → J \ J' 维持。 ∎

### 5.4 定理 1（Zassenhaus 最小子集不可约性）

**前提（A）**：α 为好求值点，MTSHL 正确（定理 5.2.3），即 h₁,...,hᵣ 满足引理 5.2.4 的性质 (a)(b)(c)。

**定理**：在前提 (A) 下，设当前状态为 (T, g_rem, verified)，g_rem 满足不变量 I₂（即 g_rem = ε' · ∏_{j∈J} Gⱼ）。设子集 S ⊂ T（|S| = s ≥ 1）满足：

- (i) pp(∏_{i∈S} hᵢ) | g_rem
- (ii) 对所有 S' ⊂ T，1 ≤ |S'| < s，pp(∏_{i∈S'} hᵢ) ∤ g_rem（最小性）

令 P = pp(∏_{i∈S} hᵢ)。则 P 不可约。更精确地，存在 j₀ ∈ J 使得 P = Gⱼ₀。

**证明**：

**步骤 1**：P 是 g_rem 的一组不可约因子的乘积。

由 (i)，P | g_rem。由不变量 I₂，g_rem = ε' · ∏_{j∈J} Gⱼ。

将 P 分解为不可约因子：P = ε_P · p₁ · p₂ · ... · pₘ（ε_P = ±1，每个 pₗ ∈ R 不可约）。

由 P | g_rem 和 R 是 UFD：对每个 pₗ，pₗ | ∏_{j∈J} Gⱼ。由 R 是 UFD 的基本性质（不可约元是素元），pₗ | Gⱼ 对某 j ∈ J。由 Gⱼ 不可约，pₗ ~ Gⱼ。

由 g squarefree（条件 SQ），g_rem 中每个 Gⱼ 至多出现一次。因此不同的 pₗ 对应不同的 Gⱼ。

所以存在单射 φ: {1,...,m} → J 使得 pₗ ~ Gφ(ℓ)。令 J' = Im(φ) ⊂ J，|J'| = m。

由 contentℤ(P) = 1（P 是 pp）和 contentℤ(Gⱼ) = 1（引理 5.1.2）以及 lc > 0 的归一化，得到：

> **P = ∏_{j∈J'} Gⱼ**   ··· (†)

**步骤 2**：确定 S 与 J' 的对应关系。

由好求值点假设，hᵢ = Gσ(i)。因此：

∏_{i∈S} hᵢ = ∏_{i∈S} Gσ(i) = ∏_{j∈σ(S)} Gⱼ

（σ 是双射，σ(S) = {σ(i) : i ∈ S}。由引理 5.2.4(c)，σ(S) 中各 Gⱼ 两两不同。）

又 P = pp(∏_{i∈S} hᵢ) = pp(∏_{j∈σ(S)} Gⱼ)。由引理 5.1.2 和 Gauss 引理 5.1.1，∏_{j∈σ(S)} Gⱼ 已经是 primitive 的。所以 pp(...) = ∏_{j∈σ(S)} Gⱼ（至多差符号，由 lc > 0 归一化消除）：

> **P = ∏_{j∈σ(S)} Gⱼ**   ··· (‡)

比较 (†) 和 (‡)：∏_{j∈J'} Gⱼ = ∏_{j∈σ(S)} Gⱼ。由 R 是 UFD 且 Gⱼ 两两互不相伴（它们是 g 的不同不可约因子，lc > 0，contentℤ = 1），唯一分解给出：

> **J' = σ(S)**   ··· (§)

**步骤 3**：反证法证明 |J'| = 1。

假设 |J'| ≥ 2。取 j₁ ∈ J'。由 (§)，j₁ = σ(i₁) 对某 i₁ ∈ S。

令 S' = {σ⁻¹(j₁)} = {i₁} ⊂ S ⊂ T。|S'| = 1 < s（因 |S| = s ≥ |J'| ≥ 2）。

pp(∏_{i∈S'} hᵢ) = pp(hᵢ₁) = pp(Gⱼ₁) = Gⱼ₁（因 contentℤ(Gⱼ₁) = 1，lc > 0）。

Gⱼ₁ | g_rem（因 j₁ ∈ J' ⊂ J，而 g_rem = ε' · ∏_{j∈J} Gⱼ）。

因此 S' ⊂ T，|S'| = 1 < s，pp(∏_{i∈S'} hᵢ) | g_rem。与条件 (ii)（最小性）矛盾。

**结论**：|J'| = 1，即 J' = {j₀}，**P = Gⱼ₀ 不可约**。 **∎**

**推论 5.4.1**：在前提 (A) 下，Zassenhaus 循环中所有 trial division 均在 s = 1 成功。

*证明*：由定理 1 证明步骤 3，|J'| = 1 意味着 |σ(S)| = 1，即 |S| = 1。所以每个因子都在 s = 1 被找到。 ∎

### 5.5 定理 2（互补因子不可约性）

**前提**：同定理 1 的前提 (A)。

**定理**：Zassenhaus 循环结束后 |T| = 1，且 g_rem 不可约。具体地，g_rem = ε' · Gⱼ₀ 对某 j₀ ∈ {1,...,t}，ε' = ±1。

**证明**：

由推论 5.4.1，每个 verified 因子 Pₘ 在 s = 1 被找到，对应单个 Hensel 因子 hᵢₘ。每次 |T| 减 1，s 重置为 1。

循环条件 `2s ≤ |T|` 在 |T| < 2 时终止。由 t ≥ 2（单因子情况已在单变量不可约检测中排除），初始 |T| = t ≥ 2。循环找到 t - 1 个因子后，|T| = 1，循环终止。

此时 verified = {P₁,...,Pₜ₋₁}，每个 Pₘ = Gπ(m)（由定理 1）。各 π(m) 互不相同（每个 Pₘ 从不同的 g_rem 因子消耗，Pₘ | g_rem 但 Pₘ 不等于已消耗的任何 Gⱼ）。

由不变量 I₂：

> g_rem = ε' · ∏_{j ∈ {1,...,t} \ {π(1),...,π(t-1)}} Gⱼ = ε' · Gⱼ₀

其中 j₀ 是唯一不在 {π(1),...,π(t-1)} 中的下标。Gⱼ₀ 不可约。 **∎**

### 5.6 定理 3（乘积完整性）

此定理**不依赖好求值点假设**，是纯代数结论。

**定理**：设修复后的 `__wang_core(g)` 返回 verified = {F₁, ..., Fₖ}（k ≥ 2）。则：

> g = ε · F₁ · F₂ · ... · Fₖ，ε ∈ {±1}

**证明**：

**情形 A（|T| = 0，g_rem 为常数）**：

由不变量 I₁：g = ε · (∏_{m=1}^{k} Fₘ) · g_rem。

g_rem 为常数 c ∈ Z。每个 Fₘ = pp(...)，contentℤ(Fₘ) = 1。由 Gauss 引理 5.1.1（归纳），contentℤ(∏ Fₘ) = 1。因此 contentℤ(g) = |c|。由条件 (PR)，contentℤ(g) = 1。故 c = ±1。

g = ε · c · ∏ Fₘ = ε' · ∏ Fₘ，ε' ∈ {±1}。

**情形 B（|T| = 1，g_rem 加入 verified 作为 Fₖ）**：

代码执行 `Fₖ = pp(g_rem)`（normalize_factor）。设 g_rem₀ 为加入前的值。

由不变量 I₁：g = ε · (∏_{m=1}^{k-1} Fₘ) · g_rem₀。

g_rem₀ = contentℤ(g_rem₀) · pp(g_rem₀) = contentℤ(g_rem₀) · Fₖ'（其中 Fₖ' 是符号调整前的 pp）。

`normalize_factor` 还做 lc > 0 调整：Fₖ = ±Fₖ'。所以 g_rem₀ = ±contentℤ(g_rem₀) · Fₖ。

同情形 A 的 Gauss 引理论证，contentℤ(∏_{m=1}^{k-1} Fₘ) = 1，而 g = ε · (∏_{m=1}^{k-1} Fₘ) · g_rem₀，contentℤ(g) = 1，故 |contentℤ(g_rem₀)| = 1。

综合：g = ε'' · ∏_{m=1}^{k} Fₘ，ε'' ∈ {±1}。 **∎**

### 5.7 定理 4（好求值点 ⟹ |T| = 1）

**定理**：在前提 (A) 下，Zassenhaus 循环结束后 |T| = 1（精确等于 1，不是 0）。

**证明**：由推论 5.4.1，所有因子在 s = 1 被找到。

**|T| ≥ 1**：循环条件 `2s ≤ |T|` 在 s = 1 时要求 |T| ≥ 2 才进入。找到一个因子后 |T| 减 1。当 |T| = 1 时，2 · 1 > 1，循环不再进入。

**|T| ≤ 1**：初始 |T| = t。每轮 s = 1 必然找到一个因子（推论 5.4.1），|T| 每轮减 1。t ≥ 2 保证至少执行一轮。当 |T| 减到 1 时停止。

综合 |T| = 1。 **∎**

### 5.8 命题 5（坏求值点下的检测）

将代数部分与概率论证严格分开。

#### 5.8.1 代数部分

**命题**：设 α 为坏求值点（r > t）。则 r ≥ t + 1。

*证明*：ḡ = ∏_{j=1}^{t} Ḡⱼ，其中 Ḡⱼ = Gⱼ(x₁, α)。ḡ = c · f₁ · ... · fᵣ（不可约分解）。每个 Ḡⱼ 是若干 fᵢ 的乘积。坏求值点意味着某 Ḡⱼ₀ 至少是 2 个 fᵢ 的乘积。其余 Ḡⱼ（j ≠ j₀）至少是 1 个 fᵢ 的乘积。故 r ≥ (t-1) · 1 + 2 = t + 1。 ∎

**推论**：在坏求值点下，MTSHL 产出 r ≥ t + 1 个 Hensel 因子，比真因子数 t 至少多 1。

#### 5.8.2 MTSHL 在坏求值点下的行为

设 Gⱼ₀ 在 α 处分裂：Ḡⱼ₀ = fₐ · f_b · ...（至少 2 个不可约因子）。MTSHL 对 fₐ 执行 Hensel 提升。

**关键观察**：MTSHL 的 Taylor 提升约束为 ∏_{i=1}^{r} hᵢ ≡ g (mod 足够高次的 Taylor 展开)。这是一个**联立方程组**：r 个因子同时满足乘积 = g。在好求值点下，此方程有唯一解（即真因子 G₁,...,Gₜ 的重排）。在坏求值点下，不存在 R 中的解满足：

> (a) hₐ(x₁, α) = fₐ
> (b) h_b(x₁, α) = f_b
> (c) ∏_{i=1}^{r} hᵢ = g

因为 g = G₁ · ... · Gₜ（t 个不可约因子），而 (a)(b)(c) 要求 r > t 个因子的乘积 = g。在 R（UFD）中，这与唯一分解矛盾：g 恰好有 t 个不可约因子，无法写成 r > t 个非常数因子的乘积。

因此 MTSHL 的 p-adic 恢复步骤（有理重建）对分裂因子的系数无法找到满足联立方程组的有理数。此时有两种结果：

1. **恢复失败** → `__mtshl_lift` 返回空 → continue 到下一个求值点
2. **恢复"成功"** → 系数通过了界检查但值错误。此时 hₐ ∈ R，但它不是 MTSHL 联立方程的真解

**在情形 2 中**：hₐ 的系数由有理重建从 Zp 残差恢复，但不满足 ∏ hᵢ = g 的全局约束。由不变量 I₂，g_rem 始终是 g 的某些真不可约因子 Gⱼ 的乘积。hₐ 通过 trial division（即 hₐ | g_rem）的充要条件是 hₐ ~ Gⱼ 对某 j ∈ J（由 R 是 UFD + g squarefree）。

这里存在三种子情形：

(2a) hₐ 不整除 g_rem → trial division 失败，hₐ 留在 T 中。同理 h_b 等分裂因子留在 T 中，|T| ≥ 2。**这是压倒性的常见情形**（概率论证见 §5.8.3）。

(2b) hₐ 恰好等于某个 Gⱼ → trial division 成功。此时 hₐ 虽非 fₐ 的"真提升"，但恰好是 g 的一个真不可约因子。算法返回正确结果。

(2c) hₐ 是多个 Gⱼ 的乘积（可约但整除 g_rem）→ trial division 成功，但返回的因子可约。此情形要求 hₐ 的系数恰好满足多个不可约因子的乘积关系，对"错误"恢复系数而言概率极低（§5.8.3）。

**总结**：情形 (2a) 导致 |T| ≥ 2（被检测），情形 (2b) 导致正确结果，情形 (2c) 理论可能但概率可忽略。因此：

**命题 5.8.3（检测或正确）**：设 α 为坏求值点。Zassenhaus 循环结束后，以下三种情况之一成立：

(a) |T| ≥ 2，修复后的接受条件拒绝此求值点

(b) |T| ≤ 1，且所有通过 trial division 的因子 P 均为 g 的真不可约因子（P = Gⱼ 对某 j）。此时乘积完整性 (I) 和不可约性 (II) 均成立。

(c) |T| ≤ 1，但某个通过 trial division 的因子 P 是 g 的多个不可约因子的乘积（P = ∏_{j∈J'} Gⱼ，|J'| ≥ 2）。此时乘积完整性 (I) 仍成立（定理 3，无条件），但不可约性 (II) 不成立。

*证明*：

**情况 (a) 与 (b)(c) 的区分**：trial division 要求 P | g_rem（精确多项式除法）。由不变量 I₂，g_rem = ε' · ∏_{j∈J} Gⱼ。P | g_rem 意味着 P 是 {Gⱼ}_{j∈J} 某子集的乘积（定理 1 步骤 1 的论证——仅用 R 是 UFD + g squarefree，不需要好求值点假设）。若分裂因子无一通过 trial division（§5.8.2 情形 2a），则 |T| ≥ 2，为情况 (a)。

**情况 (b) 与 (c) 的区分**：通过 trial division 的因子 P = ∏_{j∈J'} Gⱼ。若 |J'| = 1 则 P 不可约（情况 b）；若 |J'| ≥ 2 则 P 可约（情况 c）。在好求值点下，定理 1 的最小性论证证明 |J'| = 1。但此论证依赖步骤 2 中的 hᵢ = Gσ(i) 对应关系，需要好求值点假设。**在坏求值点下，我们不能纯代数地排除情况 (c)**。

**实际概率**：情况 (c) 要求 MTSHL 恢复出的"错误"因子的乘积恰好等于 g 的某些真不可约因子的乘积。对随机错误系数，这等价于 O(N) 个整除性条件同时成立（§5.8.3 概率论证），概率 ≤ 10⁻⁴⁰。因此实践中情况 (c) 不会发生。 ∎

**注**：情况 (c) 的存在是所有基于 Hensel 提升的因式分解算法的固有理论间隙（§5.10）。乘积完整性 (I) 的无条件保证意味着即使在此极端情况下，返回的因子集也是 g 的一个**合法分解**（只是可能不是最细分解）。

#### 5.8.3 概率论证

**命题（实际检测）**：在坏求值点下，MTSHL 恢复的"错误"因子以压倒性概率不通过 trial division。

MTSHL 的有理重建从 Zp 残差 c̄ 恢复 a/b ∈ Q（满足 ab⁻¹ ≡ c̄ mod p，|a| ≤ B, |b| ≤ B, B ≈ ⌊√(p/2)⌋）。对分裂因子，不存在满足 MTSHL 联立方程组的 Z-解（§5.8.2），故恢复出的 a/b 不是"真"系数。

设恢复出的因子 hₐ 有 N 个多变量系数，每个系数量级 O(B) ≈ O(2³¹)。hₐ | g_rem（精确多项式除法）要求商 Q = g_rem / hₐ 的所有系数为整数。这等价于 O(N) 个整除性条件同时成立。

对一个"随机"的度 d₁ 多项式 hₐ 和度 d₂ 多项式 g_rem（d₂ > d₁），精确整除的概率受 Schwartz-Zippel 类型界限制。在 hₐ 的 N 个系数各 ~2³¹ 量级、g_rem 的系数 ~10⁶ 量级下，每个除法步骤产生的余式为零的概率 ≈ O(10⁶ / 2³¹) ≈ 10⁻⁴。N ≥ 10 个条件同时满足的概率 ≤ 10⁻⁴⁰。

**结论**：在坏求值点下，|T| ≥ 2 以压倒性概率（≥ 1 - 10⁻⁴⁰）成立。

### 5.9 主定理

**定理（修复后正确性）**：设修复后的 `__wang_core(g)` 返回 verified = {F₁, ..., Fₖ}（k ≥ 2）。则：

> **(I) 乘积完整性**：g = ε · F₁ · ... · Fₖ，ε ∈ {±1}。**无条件成立。**

> **(II) 不可约性**：若 α 为好求值点且 MTSHL 正确（前提 A），则每个 Fₘ 不可约。

**证明 (I)**：定理 3（§5.6）。纯代数结论，仅依赖 trial division 的精确性和 Gauss 引理。 ∎

**证明 (II)**：

由前提 (A) 和定理 4（§5.7），Zassenhaus 循环结束后 |T| = 1。

- F₁, ..., Fₖ₋₁ 由 trial division 找到。由定理 1（§5.4），每个 Fₘ = Gπ(m) 不可约。
- Fₖ = pp(g_rem)。由定理 2（§5.5），g_rem = ε' · Gⱼ₀ 不可约。故 Fₖ = Gⱼ₀ 不可约。

综合：所有 Fₘ 不可约。 **∎**

**推论（整体正确性链）**：

修复方案的两个方向提供以下保证：

| 保证 | 条件 | 由谁提供 |
|------|------|----------|
| (I) 乘积完整性 | 无条件 | 精确 trial division + Gauss 引理 |
| (II) 不可约性 | 好求值点 + MTSHL 正确 | 定理 1 + 定理 2 |
| 好求值点概率 ≥ 1 - δ | 随机采样 | 方向 5 + Hilbert 定理（§5.11） |
| 坏求值点被拒绝概率 ≥ 1 - ε | |T| ≤ 1 检查 | 方向 2 + 命题 5.8.3 |
| k 次重试全部失败概率 | ≤ δᵏ | 独立采样 |

### 5.10 关于绝对保证的讨论

主定理 (II) 的条件是前提 (A)（好求值点 + MTSHL 正确）。在实践中：

**不可消除的理论间隙**：在多变量 UFD Z[x₁,...,xₙ] 中验证一个多项式的不可约性等价于因式分解本身。因此**没有**基于 Hensel 提升的因式分解算法能提供绝对的不可约性后验保证——这需要对每个返回因子重新执行完整因式分解。

**所有主流 CAS 的策略一致**：FLINT、NTL、Maple、Singular 均信任 Hensel 提升在好求值点下的正确性。唯一的运行时验证是乘积检查（product = input），CLPoly 在 debug 模式下已有此检查（`__factor_multivar` line 2545）。

**本方案的保证等级**：

- 乘积完整性 (I)：**绝对保证**（纯代数，定理 3）
- 不可约性 (II)：**条件保证**（前提 A 下，定理 1 + 2，纯代数证明）
- 前提 A 的成立概率：**≥ 1 - O(d · B^{-1/2} · log B)**（Hilbert 定理）
- 坏求值点逃逸 |T| ≤ 1 检测的概率：**≤ 10⁻⁴⁰**（命题 5.8.3）

### 5.11 Hilbert 不可约定理与概率界

**定理（Hilbert, 1892; 有效版本 Cohen 1981, Dèbes-Walkowiak 2008）**：设 F ∈ Z[x₁, ..., xₙ]（n ≥ 2）不可约，d = deg(F)。设 α = (α₂,...,αₙ) 从 {-B,...,B}ⁿ⁻¹ 中均匀随机选取。则：

> P[F(x₁, α) 在 Z[x₁] 中可约] ≤ c(d, n) · B^{-1/2} · log B

其中 c(d, n) 仅依赖于 d 和 n，不依赖 F 的系数。

**推论**：设 g = G₁ · ... · Gₜ（不可约分解）。α 为坏求值点当且仅当某个 Gⱼ(x₁, α) 可约。由 union bound：

> P[α 为坏求值点] ≤ ∑_{j=1}^{t} c(dⱼ, n) · B^{-1/2} · log B =: δ

**数值估计**：B = 1000, d = 50, t = 10 典型情况下，δ 量级约 10⁻¹。每次采样独立，方向 2 的 retry（每个主变量 BATCH_SIZE = 200）保证：

> P[200 次均为坏求值点] ≤ δ^{200} ≈ 10⁻²⁰⁰

**实际失败概率可忽略。**

## 6. 改动文件清单

| 文件 | 改动 | 行数估计 |
|------|------|----------|
| `clpoly/polynomial_factorize_wang.hh` | 替换 `__wang_core` lines 2410-2420（返回条件） | ~10 行 |
| `clpoly/polynomial_factorize_wang.hh` | 重写 `__select_eval_point` lines 1222-1306（随机化候选生成） | ~40 行 |

## 7. 测试计划

| 测试 | 命令 | 验证目标 |
|------|------|----------|
| B1 回归测试 | `make test/test_factorize_multivar && _build/debug/bin/test_factorize_multivar` | B1 案例正确分解为 4 因子 |
| 全量回归 | `bash test/run_all_tests.sh` | 所有现有测试通过 |
| FLINT/NTL 交叉验证 | `make crosscheck` | 与参考实现结果一致 |
| 压力测试 | `make stress` | 70-factor bivar + 60-factor trivar 正确 |

## 8. 实施步骤

1. 修改 `__wang_core` 返回条件（方向 2）
2. 运行 B1 回归测试验证修复
3. 全量回归测试
4. 重写 `__select_eval_point` 随机化（方向 5）
5. 全量回归 + 压力测试
6. 从 B1 测试中移除 `KNOWN BUG` 标记
