# 论文摘要：Sparse Multivariate Hensel Lifting: A High-Performance Design and Implementation

> Monagan & Tuncer, ICMS 2018, LNCS 10931:378–387
> PDF: `ICMS2018-MonaganTuncer.pdf`
> 角色：**MTSHL 高性能实现架构**，双变量 Hensel 提升单元 + 并行求值框架，Maple 2019 C 实现基础

---

## 主题

将 MTSHL（CASC 2016 的 2 因子稀疏 Hensel 提升）重新组织为高性能并行实现。核心思路：将多变量 Hensel 提升分解为两个独立可并行化的核心计算——**多变量多项式求值**和**双变量 Hensel 提升**（HenselLift1）。

---

## 增量误差更新（§2，Bernardin 改进）

在 MHL 的每一步 i，不重新计算完整乘积 fj×gj，而是利用已有的 σk, τk 系数：

```
ci = a^(i)(αj)/i! - Σ_{k=1}^{i-1} σk · τ_{i-k}
```

其中 a^(i) 是 aj 关于 xj 的 i 阶导数。这将误差更新从 O(d1dj²) 降至 O(d1dj + id1²)。

---

## HenselLift1（Algorithm 2）：双变量 Hensel 提升

**输入**：p, αj ∈ Zp，a ∈ Zp[x1,xj]，f0, g0 ∈ Zp[x1]（首一），满足 a(x1,αj) = f0g0

**输出**：fj, gj ∈ Zp[x1,xj] 使得 a = fj·gj

```
Algorithm 2 HenselLift1:
1: Solve sg0 + tf0 = 1 in Zp[x1] via Extended Euclidean Alg.      O(d1²)
2: for i = 1, 2, ... while df + dg < da do
3:   a ← ∂a/∂xj                                                    O(d1dj)
4:   ci ← a(x1, αj)/i! - Σ_{k=1}^{i-1} σk(x1)τ_{i-k}(x1)        O(d1dj) + O(id1²)
5:   [σi ← (ci·s) rem f0;  τi ← (ci - σi·g0) quo f0]             O(d1·deg(f0))
     【优化】：先求值、后差值 x1：
       σil ← σ_{i-1}(l), cil ← Σσkl×τ(i-k)l, 再插值             O(d1²) + O(d1dj)
6: end for
总复杂度：O(d1²dj + d1dj²)
```

**关键优化（Step 5）**：将 σk 在 x1 = 0,...,d1 处的值缓存，下一步迭代直接复用；Horner 法求值 + Newton 插值，将 Σ 的代价从 O(id1²) 降为 O(d1²) + O(d1dj)。

---

## Algorithm 3：多变量 → 双变量归约

将提升 xj（在 Zp[x1,...,xj] 中）归约为 s 个独立的双变量 Hensel 提升：

```
【输入】aj ∈ Zp[x1,...,xj], fj-1, gj-1 ∈ Zp[x1,...,xj-1]

步骤：
1. 在 s 个随机点 (β₁k,...,β^k_{j-1}) 处对 x3,...,xj-1 求值
   → 得 s 个双变量像 aj(x1,x2,βk,xj), fj-1(x1,x2,βk), gj-1(x1,x2,βk)
2. 对每个双变量像，调用 HenselLift1（并行！）
   → 得 fj(x1,x2,βk,xj), gj(x1,x2,βk,xj)
3. 再对 x2 在 deg(aj,x2)+1 个点求值 → 调用更多 HenselLift1
4. 对 x3,...,xj-1 做稀疏插值，对 x2 做稠密插值
   → 恢复 fj, gj ∈ Zp[x1,...,xj]
```

**分层结构（Figure 1，同构图）**：
```
aj(x1,...,xj)
    ↓ evaluate x3,...,xj-1 at s points
aj(x1,x2,βk,xj) → [HenselLift1] → fj(x1,x2,βk,xj)
    ↓ evaluate x2 at deg+1 points
aj(x1,γl,βk,xj) → [HenselLift1] → fj(x1,γl,βk,xj)
    ↑ dense interpolate x2
    ↑ sparse interpolate x3,...,xj-1
fj(x1,...,xj) ←
```

**求值次数 s**：与因子的项数 T 相关，实践中 s << deg(aj,xj)。

---

## 数据结构（§2.1）

```
Zp[x1]（单变量）：
  dense array of long long int: [c0, c1, ..., cd]
  8B/项，连续内存，SIMD 友好

Zp[x1,...,xn]（多变量）：
  pair-of-arrays: A[t] + X[t]
  A[i] = 系数（64-bit int）
  X[i] = 单项式（64-bit packed 整数）
    e.g., x1^i x2^j x3^k → 2^42·i + 2^21·j + k
  每个 monomial 占 8B，避免指针跳跃
```

**对比 CLPoly**：CLPoly 用 `vector<pair<umonomial, Zp>>`，24B/项 + 堆分配，Monagan 用 C 原始数组，8B/项，无动态分配。

---

## 并行化（Cilk C）

每个双变量提升（HenselLift1 调用）是独立任务，分配到 Cilk 工作线程：
- 单次任务粒度 ≥ 10³ 次 Zp 乘法（满足 Cilk 启动开销阈值）
- **内存预分配**：每个 HenselLift1 实例的所有多项式空间在调用前分配好；HenselLift1 内部不分配内存（无 malloc）
- N 个评估点 → N 个并行 HenselLift1

**评估并行化**：将多项式分成 N = 核心数 个块，使用 `cilk_spawn` 并行求值。

---

## 性能数据

对 2 因子、6–15 变量、100–8000 项多项式：
- 单核：MTSHL ≫ Wang（稀疏情形，αi ≠ 0）
- 12 核：额外 10–12x 并行加速（HenselLift1 独立并行）
- 求值是最大瓶颈（observation from implementation）：几乎所有时间在 step 8（evaluate aj at points Yk）

---

## 对 CLPoly 的意义

- **HenselLift1** = 核心计算单元，替换 CLPoly 的 `__hensel_lift_one_var` 中的双变量部分
- **数据结构**：P1c（DDF/EDF 稠密化）可参考此处的密集 `long[]` 设计
- **Algorithm 3 双变量归约** = P2-Sparse-D 的具体化（bivariate projection 方案）
- **内存预分配**：CLPoly 性能优化的重要参考（避免稀疏分配开销）
- 本文的实现是 MTSHL 被 Maple 采用的直接原因（Baris Tuncer 2018 年通过 MITACS 将其集成到 Maple 2019）
