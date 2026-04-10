# 论文摘要：Using Sparse Interpolation in Hensel Lifting

> Monagan & Tuncer, CASC 2016, LNCS 9890:381–400
> PDF: `CASC2016-MonaganTuncer.pdf`
> 角色：**MTSHL 算法起源**，2-因子情形完整实现，Maple 2019 基础

---

## 主题

在 Wang 的多变量 Hensel 提升（MHL）框架内，用**稀疏插值**替代 Wang 的递归 MDP 求解器，将 MDP 求解从指数复杂度降为多项式复杂度。

---

## 问题定义（MDP）

多变量 Diophantine 问题：给定 u, w, c ∈ Zp[x1,...,xj]，求 σ, τ满足

```
σu + τw = c  mod I_j^{dj+1}
```

Wang 的 MDP 解法是递归 ideal-adic 展开，当 αk ≠ 0 时关于 n 指数复杂。

---

## 核心观察（Lemma 1）

设 f ∈ Zp[x1,...,xn]，α 随机选自 Zp，令 f = Σ bi(xn-α)^i，则：

```
Pr[Supp(bj+1) ⊄ Supp(bj)] ≤ |Supp(bj+1)| · (dn-j)/(p-dn+j+1)
```

当 p 足够大时，Supp(bj+1) ⊆ Supp(bj) 以高概率成立。这说明每步 MDP 解的支撑是前一步的子集，可以复用前一步的骨架。

---

## 算法结构

### Algorithm 4：j-th step of MTSHL（2-因子）

```
for i = 1, 2, ... while error ≠ 0:
  c ← i-th Taylor coefficient of error at xj = αj
  if c ≠ 0:
    σg ← skeleton of τj,i-1          ← 强 SHL 假设：复用前一步支撑
    用稀疏插值（§3.2 + BDP）求解 MDP σj0·τji + τj0·σji = c
    if FAIL → BSDiophant（稠密回退）
    if BSDiophant FAIL → 重新选 ideal
    更新 (σj, τj) 和 error
```

### Algorithm 2：BSDiophant（递归 MDP 求解器）

```
BSDiophant(u, w, c, p):
  if n=2: 调用 BDP（密集双变量求解，O(d³)，C 实现）
  否则：
    取随机 β₁，递归 BSDiophant(u|_{xj=β1}, ...)
    以 σ₁ 的骨架为形式 σf
    repeat:
      取 βk，用稀疏插值（带形式 σf）求 σk
      Newton 插值更新 σ
    until σ 稳定 且 w | (c - σu)
    τ ← (c - σu)/w
```

### §3.2 双变量投影（bivariate projection）

关键改进：将稀疏插值投影到 Zp[x1,x2]（而非 Zp[x1]），通过 BDP 求解双变量 MDP，再稀疏插值恢复 x3,...,xj 的系数。评估点数从 t 降至 s < t（t/s² 复杂度降低）。

---

## Evaluation 优化（§3.3）

对 t 个几何点 (α₃ʲ,...,αₙʲ), j=1,...,t 的批量求值：

```
预计算：θi = ∏k αk^{exponent of xk in term i}
初始：c(0) = 各项系数
迭代：c(j)_i = c(j-1)_i · θi  （每步 s 次乘法）
```

代价从 s(n-2)t 降至 **st**，节省因子 (n-2)。Roman Pearce 用 C 实现。

---

## 误差增量更新（§4.3）

Taylor 系数更新公式（3x 快于导数公式）：

```
c^(i+1)_j = (c^(i)_j - U^(i)) / (xj - αj)
其中 U^(i) = σj^(i-1)·τji + τj^(i-1)·σji
```

---

## Ideal Type 限制

| Ideal 类型 | αi 取值 | 结果 |
|-----------|---------|------|
| Type 1 | αi = 0 | 稀疏多项式保持稀疏，MDP 少，SHL 不使用；Wang 更快 |
| Type 2 | αi ≠ 0（随机） | Taylor 展开稠密，MDP 多；MTSHL 大幅胜出 |

---

## 性能数据（Table 4，ideal type 2）

| n/d/t | Wang | tNBS（本文 SHL） | 加速 |
|--------|------|---------------|------|
| 3/35/100 | 2.87s | 0.32s | 9x |
| 5/35/100 | 88.1s | 1.16s | 76x |
| 7/35/100 | 800s | 1.58s | **506x** |

---

## 对 CLPoly 的意义

- Algorithm 4 = P2-Sparse-B 的直接实现蓝图
- BDP（密集双变量求解器）= 性能关键，需 C 实现
- Evaluation θ-array 优化 = 必要工程细节
- Ideal type 1 限制 = 分流策略的依据
