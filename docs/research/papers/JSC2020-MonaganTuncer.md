# 论文摘要：The Complexity of Sparse Hensel Lifting and Sparse Polynomial Factorization

> Monagan & Tuncer, J. Symbolic Computation 99:189–230 (2020)
> PDF: `JSC2020-MonaganTuncer.pdf`
> 角色：**MTSHL 复杂度分析**，修正 Zippel 假设，给出平均代价 Theorem 19，理论验收标准

---

## 主题

- 研究多变量多项式在逐步求值时项数的期望值和方差（Section 4）
- 修正 Zippel 1979 稀疏插值复杂度分析中的错误假设（Section 4）
- 给出 MTSHL（Algorithm 5）的平均复杂度（Sections 5–7）
- 实验数据与 Wang/Magma/Singular 对比（Section 8）

---

## 逐步求值的项数变化（Section 4）

**实验观察**（12 变量，10000 项多项式，αi 非零随机）：

```
j   : 12    11    10     9     8     7     6     5     4    3    2   1
#aj : 10000 9996  9954  9802  9207  7550  4837  2478  978  304  68  12
```

关键现象：前半（j=12..7）项数几乎不变；后半（j=6..1）快速减少。

**Taylor 系数稀疏性**（同一多项式的 σi 序列）：

```
i    :  0     1     2     3     4    5    6    7    8   9  10  11
#σi  : 9996  5526  2988  1504  760  343  158   60   28   8   3   1
```

σi 的项数快速指数衰减 → MTSHL 代价随 i 增大而剧减。

**Theorem（Section 4）**：对 Tp 项、n 变量多项式，在 xn = α（随机）后：

```
E[#an-1] = Tp · (1 - 1/p)^{Tp-1} ≈ Tp / e  （对大 Tp，p 远大于 Tp）
```

当 Tp << p（稀疏多项式）时，E[#an-1] ≈ Tp（几乎无减少），与实验一致（前半几乎不变）。

---

## WMDS 的指数性（Wang MDP 求解器）

Algorithm 2（WMDS）的 Euclidean 算法调用次数：

```
M(x1) = 1
M(xj) ≤ (dj - 1) · M(xj-1)
M(xn-1) ≤ ∏_{i=2}^{n-1} (di - 1)  →  指数于 n
```

当所有 αj ≠ 0 时，Taylor 系数稠密（ci ≠ 0），导致完全指数爆炸。

---

## Algorithm 5 (MTSHL)：完整算法描述

JSC 2020 将整个 MTSHL 作为 Algorithm 5 正式定义：

```
Algorithm 5 MTSHL（j-th step，2 因子）：

1: fj ← fj-1; gj ← gj-1
2: 用 Supp(fj-1) 作为骨架 (form_σ) 初始化
3: for k = 1, 2, ... while error ≠ 0:
4:   ck ← Taylor coeff of (xj - αj)^k in (aj - fj·gj)
5:   if ck ≠ 0:
6:     [σk, τk] ← SparseInt(fj-1, gj-1, ck, form_σ, form_τ)
        → if FAIL: SparseInt with new random points
        → if still FAIL: WMDS（Wang 回退）
7:     form_σ ← Supp(σk)  ← 更新骨架（下次复用）
8:     fj += σk·(xj - αj)^k; gj += τk·(xj - αj)^k
9:     error ← aj - fj·gj
```

---

## 平均复杂度（Theorem 19）

**主定理**（JSC 2020 核心结果）：

设 f, g ∈ Zp[x1,...,xn] 各有 Tg 项，次数上界 d，理想 type 2（αi 随机非零）。MTSHL 一次 j-step 的期望代价为：

```
O(n²/d · Tg³)   （evaluation dominates）
```

**与 Wang 的对比**：

| 算法 | ideal type 1 (αi=0) | ideal type 2 (αi≠0) |
|------|---------------------|---------------------|
| Wang WMDS | O(Tg^n)（指数！） | O(Tg^n)（指数！） |
| MTSHL | 优势小（αi=0 时稀疏保持稀疏） | **O(n²/d·Tg³)（多项式！）** |

---

## Ideal Type 区分（§重要概念）

| 类型 | 定义 | 特征 |
|------|------|------|
| **Type 1** | αi = 0（零理想） | 稀疏多项式保持稀疏；Taylor 展开少有非零项；Wang 快；SHL 不使用 |
| **Type 2** | αi ≠ 0（随机非零理想） | 即便稀疏输入，Taylor 展开也稠密；Wang 指数；MTSHL 多项式时间 |

**实际影响**：在因式分解算法中，为避免 lc(f) 为零或 f 有重根，必须选择非零 αi，**几乎所有实际情形都是 Type 2**。

---

## 性能数据（Section 8）

**Table 6/7**（与 Wang, Magma, Singular 对比，ideal type 2）：

| 用例 | Wang | Magma | Singular | MTSHL |
|------|------|-------|----------|-------|
| n=7, Tg=100, d=35 | 800s | 382s | >1000s | 1.58s |
| n=5, Tg=100, d=35 | 88.1s | 45.2s | - | 1.16s |

---

## Zippel 复杂度分析的修正

JSC 2020 指出 Zippel 1979 的复杂度分析有一个错误假设：Zippel 假设 `#ai << #a` 当 i 约为 n/2 时，但实验和理论均表明 `#ai ≈ #a` 直到 i 约为 n/2，之后才快速减少。JSC 2020 给出了正确的期望值公式。

---

## 多因子限制（JSC 2020 §7）

对 r > 2 因子（MTSHL-d，CASC 2018），JSC 2020 指出：
- 主要代价仍在求值
- 期望代价：O(n²/d · r · Tg³) — 与 r 线性
- MTSHL-d 比归约到 2 因子快约 (r-1) 倍（与 CASC 2018 数据一致）

---

## 对 CLPoly 的意义

- **Theorem 19** = Phase 2 Zippel/MTSHL 实现后的验收标准参考（多项式 vs 指数）
- **ideal type 1/2 区分** = `__select_eval_point` 的评估点选取策略（可以考虑两类分别处理）
- **σi 项数衰减** 解释了为何 MTSHL 在高精度时仍有优势（σk 的 Supp 越来越小）
- **Zippel 错误假设修正**：CLPoly Phase 2 实现时，项数估算不能用 Zippel 的原始假设
- **JSC 2020 是 MTSHL 的最终权威描述**（覆盖复杂度证明 + 完整算法），实现时应以此为主要参考
