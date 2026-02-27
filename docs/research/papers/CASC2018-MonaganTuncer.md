# 论文摘要：Factoring Multivariate Polynomials with Many Factors and Huge Coefficients

> Monagan & Tuncer, CASC 2018, LNCS 11077:319–334
> PDF: `CASC2018-MonaganTuncer.pdf`
> 角色：**MTSHL 多因子扩展 + 大系数 p-adic 提升**，Algorithm 3 (SparseInt) + Algorithm 4 (MTSHL-d) 是 Maple 2019 默认算法核心

---

## 主题

将 CASC 2016 的 2 因子 MTSHL 扩展到以下两种重要情形：
1. **r > 2 因子**：直接多因子 MDP 求解（MTSHL-d）
2. **大整数系数**：p-adic 提升优化（机器素数 p + sparse lift to Z_{p^l}）

---

## 问题重申（多因子 MDP）

**r 因子 MDP**（Algorithm 1, line 8）：给定 c ∈ Zp[x1,...,xj-1]，求 σk,i 使得：

```
σ1,i·b1 + σ2,i·b2 + ... + σr,i·br = c
其中 bk = ∏_{l≠k} fj-1,l
```

**旧方法**（GCL §6）：将多因子 MDP 归约为 r-1 个 2 因子 MDP。代价：(r-1) 次双因子 MDP 求解。

**MTSHL-d 方法**（Algorithm 4，本文新贡献）：直接同时求解所有 r 个 σk,i，代价降低至多 r-1 倍。

---

## Algorithm 3 (SparseInt)：MTSHL 的 MDP 求解器

稀疏插值求解 MDP σu + τw = c mod I^{dj+1}（2因子情形）：

```
SparseInt(u, w, c, p, α, form_σ, form_τ):
  1. 使用 form_σ（= Supp(σ_{i-1,j})）构建线性方程组
  2. 选随机点 β ∈ Zp，求值：u(β), w(β), c(β) → 解出 σ(β)
  3. 取足够多 β（|Supp(form_σ)| 个点）
  4. 解 Vandermonde 方程组 → 恢复 σ 系数
  5. τ = (c - σ·u) / w（验证整除）
  6. if FAIL: 回退到 Wang's WMDS（Algorithm 2）
```

**强 SHL 假设**：Supp(σk,i) ⊆ Supp(σ_{k-1,i})（Lemma 1，概率 > 1 - T·d/p）

---

## Algorithm 4 (MTSHL-d)：直接多因子 SHL

核心改进：对 r 个因子同时求解 MDP，而非逐对归约：

```
Algorithm 4 MTSHL-d（j-th step，multi-factor）：
for k = 1, 2, ... while error ≠ 0:
  ck ← Taylor coeff of (xj - αj)^k in error
  if ck ≠ 0:
    form_σk,i ← Supp(σ_{k-1,i})   ← 强 SHL 假设复用前一步骨架
    σk,1,...,σk,r ← SparseInt-multi(b1,...,br, ck, form)
      若 FAIL → multi-BDP（稠密多因子双变量求解器）
    更新每个 fj,i += σk,i·(xj - αj)^k
    更新 error
```

**关键区别（vs 旧方法）**：
- 旧方法：σk = (σk,1,...,σk,r-1)，通过 r-1 次 2 因子 WMDS 求解（串行）
- MTSHL-d：同时处理所有 r 个 σk,i（并行友好，节省 r-1 倍乘法）

---

## multi-BDP：多因子双变量 Diophantine 求解器

BDP 的 r 因子推广：

```
multi-BDP(f1,...,fr, c, p):
  在 Zp[x1,x2] 中求解 σ1b1 + ... + σrbr = c
  其中 bi = ∏_{l≠i} fl

  基本情形（j=1）：扩展欧几里得算法（r 因子版本）
  递归：同 BDP 但处理 r 个系数
  代价：O(r·d1·d2²)（vs 2-BDP: O(d1·d2²)）
```

---

## Algorithm 5：大系数 p-adic 提升

**问题**：系数界 B 很大时（||f|| > p_machine），经典做法选 p^l > 2B，导致 Z_{p^l} 下大数运算。

**新方法**（基于稀疏 p-adic 结构）：

```
Theorem 1: 若 p 随机选，f = Σ fk·p^k（p-adic 展开），则
  Pr[Supp(fk) ⊆ Supp(fk-1) for all k] > 1 - t·l / (π(2^{m+1}) - π(2^m) - l)
其中 t = 项数，l = p-adic 位数，m = 机器 word 大小

对 31-bit 机器素数，这个概率极高（接近 1）。
```

**Algorithm 5 流程**：

```
1. 在 Zp 中完整运行 MTSHL-d（p = 63-bit 机器素数）
   → 得 f1,...,fr mod p（机器精度）
2. 对每个 p-adic 位 k = 1, 2, ..., l:
   dk = (a - Π fi) / p^k  （差商）
   用 SparseInt 解 MDP：Σ σk,i·bi = dk（在 Zp 中）
   fi ← fi + σk,i·p^k
3. 最终 f1,...,fr 满足 Π fi = a（over Z）
```

**收益**：99%+ 的运算在 machine prime p 下进行，仅最后 l 步做模 p 运算（不做大数运算）。

---

## 性能数据（Table 3–5）

**多因子 (r>2) 对比**：

| 用例 | 旧 MTSHL（归约） | MTSHL-d（直接） | 加速 |
|------|----------------|----------------|------|
| 5var, 50term, r=4 | 1.23s | 0.41s | 3x |
| 5var, 50term, r=6 | 2.31s | 0.48s | ~5x |
| 5var, 50term, r=8 | 4.18s | 0.61s | ~7x |

加速比接近 r-1，验证了理论预测。

**大系数 p-adic 对比**：
- 当 ||fi|| 很大时（如 Swinnerton-Dyer 多项式），p-adic 方法比 Z_{p^l} 方法快 3–10x

---

## 与其他论文的关系

| 论文 | 贡献 | 关系 |
|------|------|------|
| CASC 2016 | 2 因子 SHL（Algorithm 4 旧版） | 本文 Algorithm 3 = 2 因子 SparseInt；Algorithm 4 = 多因子扩展 |
| ICMS 2018 | 并行 HenselLift1 | 与本文正交（本文重点是多因子 MDP，ICMS 重点是并行求值） |
| MC 2019 | Maple 2019 全貌 | MTSHL-d（本文 Algorithm 4）是 Maple 2019 的默认算法 |
| JSC 2020 | 复杂度证明 | 本文 Algorithm 4 的复杂度在 JSC 2020 中分析 |

---

## 对 CLPoly 的意义

- **Algorithm 3 (SparseInt)** = P2-Sparse-B 的核心实现蓝图（多变量 Diophantine 稀疏求解）
- **Algorithm 4 (MTSHL-d)** = P2-Sparse-B 的实际目标（直接多因子，非归约到 2 因子）
- **multi-BDP** = SparseInt 的回退求解器（2 因子 BDP 的 r 因子扩展）
- **Algorithm 5 (p-adic)** = Phase 1（模 Bézout）的精细化版本；CLPoly Phase 1 已覆盖基础需求
- **Theorem 1** 为 p-adic 支撑假设提供概率保证（与 Phase 1 的理论基础相同，仅数量级不同）
