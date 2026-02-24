# 调研报告：MTSHL 多变量因式分解算法全貌

> 调研日期：2026-02-24
> 前置报告：`multivar-factorization-research.md`（2026-02-23，根因分析 + Phase 1/2 路线）
> 本报告：5 篇 Monagan-Tuncer 论文的综合整理，为 Phase 2 实现提供完整算法细节
> 参考架构：`docs/design/multivar-opt/architecture.md`

---

## 1. 论文系列总览

| 论文 | 关键贡献 | 对 CLPoly 的意义 |
|------|----------|----------------|
| **CASC 2016** (LNCS 9890:381–400) | MTSHL 核心思想；2 因子 SHL；BDP；θ-array 求值优化；Lemma 1 | P2-Sparse-B 直接实现蓝图 |
| **ICMS 2018** (LNCS 10931:378–387) | 高性能实现架构；HenselLift1；Algorithm 3 双变量归约；Cilk C 并行 | P2-Sparse-D 数据结构 + 求值优化 |
| **CASC 2018** (LNCS 11077:319–334) | 多因子扩展（MTSHL-d）；SparseInt；multi-BDP；p-adic 提升 | P2-Sparse-B 多因子版本；Phase 1 精细化 |
| **MC 2019** (Extended Abstract) | Maple 2019 集成报告；两组基准（稀疏随机 + 循环矩阵） | 实现目标 + 性能基准 |
| **JSC 2020** (JSC 99:189–230) | 完整复杂度分析；Theorem 19；修正 Zippel 假设；ideal type 1/2 | 验收标准 + 复杂度参照 |

PDF 存档：`docs/research/papers/`
摘要 MD：`docs/research/papers/{CASC2016,ICMS2018,CASC2018,MC2019,JSC2020}-MonaganTuncer.md`

---

## 2. 算法族谱与核心思想

### 2.1 Wang MHL（当前 CLPoly 实现）

```
Wang MHL（GCL §6，CLPoly __multivar_hensel_lift）：
  for j = 2 to n:  ← 逐变量提升 xj
    for k = 1, 2, ...:  ← Taylor 阶
      ck = Taylor coeff of (xj - αj)^k in error
      if ck ≠ 0:
        [σk, τk] ← WMDS(fj-1, gj-1, ck)  ← Wang MDP 求解器（递归，GCL Alg）
        fj += σk·(xj-αj)^k
        gj += τk·(xj-αj)^k
```

**WMDS 代价**：M(xj) ≤ ∏_{i=2}^{j-1} (di-1) → **指数于 j**（当所有 αi ≠ 0 时）

**CLPoly 额外问题**：Bézout 链在 Z[x] 上构建（非 Zp[x]），bezout_s[0] 次数 O(r²)、系数 O(r log r!) bit，导致大整数运算爆炸（详见 `multivar-factorization-research.md §2.2`）

### 2.2 MTSHL 核心观察（Lemma 1）

**关键引理**（CASC 2016 Lemma 1，所有后续论文的基础）：

设 f ∈ Zp[x1,...,xn]，α 随机选自 Zp，令 f = Σ bi(xn-α)^i。则：

```
Pr[Supp(bj+1) ⊄ Supp(bj)] ≤ |Supp(bj+1)| · (dn-j) / (p - dn+j+1)
```

**强 SHL 假设**：当 p 足够大时，Supp(σk,i) ⊆ Supp(σ_{k-1,i}) 以高概率成立。

**实际验证**（JSC 2020，Tg=10000 项多项式的 σi 序列）：

```
i    :  0     1     2     3     4    5    6    7    8   9  10  11
#σi  : 9996  5526  2988  1504  760  343  158   60   28   8   3   1
```

→ σi 项数指数衰减，MTSHL 代价随 i 增大而剧减。

### 2.3 MTSHL 用 Supp(σ_{k-1}) 作为骨架求解 MDP

```
SparseInt（CASC 2018 Algorithm 3 / JSC 2020 §2）：
  form = Supp(σ_{k-1,j})       ← 已知骨架（t 个单项式）
  1. 在 t 个随机点 β_l 处求值：solve σ(β_l)·fj-1(β_l) + τ(β_l)·gj-1(β_l) = ck(β_l)
     → t 个标量 σ(β_l)（Euclidean Alg in Zp[x1]）
  2. 构造 t×t Vandermonde 方程组（关于单项式系数）
  3. 解 Vandermonde 系统 → 恢复稀疏 σk            O(t²)
  4. τk = (ck - σk·fj-1) / gj-1                   验证整除
  5. if FAIL: 回退到 BDP/WMDS
```

**代价**：O(t²) per step（t = 因子项数），vs Wang O(d^{j-1}) per step。

---

## 3. 算法变体详解

### 3.1 CASC 2016 Algorithm 4（2 因子 SHL，v1）

```
for i = 1, 2, ...:
  c ← i-th Taylor coeff of error
  if c ≠ 0:
    σg ← skeleton of τ_{j,i-1}       ← 复用前一步 Supp
    稀疏插值（§3.2 + BDP）求 σji, τji
    if FAIL: BSDiophant（稠密回退）
    更新 (σj, τj) 和 error
```

**§3.2 双变量投影**：evaluate x3,...,xj-1 → 在 Zp[x1,x2] 上用 BDP 求解 → 稀疏插值恢复高维系数。

**BSDiophant（Alg 2）**：递归 MDP 求解器，基本情形 = BDP（密集双变量求解，O(d³)，C 实现）。

**θ-array 求值优化**（§3.3）：预计算 θi = ∏k αk^{exp(xk, term_i)}，每步求值仅需 t 次乘法（vs 朴素 t(n-2) 次）。

### 3.2 ICMS 2018 Algorithm 2 (HenselLift1) + Algorithm 3

**HenselLift1**（核心计算单元）：在 Zp[x1,xj] 中做双变量 Hensel 提升

```
代价：O(d1²dj + d1dj²)
优化：σkl 值缓存复用 → Σ 的代价从 O(id1²) → O(d1²) + O(d1dj)
```

**Algorithm 3**（多变量 → 双变量归约）：

```
提升 xj:
  for k = 1..s:
    evaluate x3,...,xj-1 at (β3k,...,βj-1k) → bivariate image
    for l = 1..deg(aj,x2):
      evaluate x2 at γl → univariate image
      HenselLift1(univariate images)  ← 并行执行
    dense interpolate x2
  sparse interpolate x3,...,xj-1
```

**数据结构**：

```
Zp[x1]：dense long[] (8B/项)
Zp[x1,...,xn]：pair (A[], X[])，X[i] = 2^42·ex1 + 2^21·ex2 + ex3（64-bit packed）
→ 无动态内存分配，SIMD 友好
```

### 3.3 CASC 2018 Algorithm 4 (MTSHL-d)：直接多因子

对 r > 2 因子，同时求解所有 σk,i（而非归约到 r-1 个 2 因子问题）：

```
for k = 1, 2, ...:
  ck ← Taylor coeff of (xj - αj)^k in error
  σk,1,...,σk,r ← SparseInt-multi(b1,...,br, ck, forms)
    bi = ∏_{l≠i} fj-1,l   ← 不显式构造，只需计算求值点处的 bi(β)
    if FAIL: multi-BDP（稠密 r 因子双变量求解，O(r·d1·d2²)）
  for i = 1..r: fj,i += σk,i·(xj-αj)^k
```

**加速比**：O(r-1) vs 旧归约法（CASC 2016）。实测：r=8 时 ~7x 加速。

### 3.4 CASC 2018 Algorithm 5：大系数 p-adic 提升

```
1. 用机器素数 p（63-bit）完整运行 MTSHL-d → 得 f1,...,fr mod p
2. 对 k = 1, 2, ..., l（l = ⌈log_p(2B)⌉）:
   dk = (a - Π fi) / p^k
   用 SparseInt 解 Σ σk,i·bi = dk（在 Zp 中）  ← Supp(σk,i) ⊆ Supp(σk-1,i) w.h.p.
   fi += σk,i·p^k
```

**Theorem 1（p-adic 支撑假设）**：

```
Pr[Supp(fk) ⊆ Supp(fk-1) for all k] > 1 - t·l / (π(2^{m+1}) - π(2^m) - l)
```

对 31-bit 素数，概率极高（实验中 SparseInt 从未在此处 FAIL）。

---

## 4. 复杂度对比

### 4.1 Theorem 19（JSC 2020）

设 a = f·g，f,g 各有 Tg 项，次数 d，n 变量，大素数 p，失败概率 ≤ (n-2)d²Tf/(2(p-2d+1)) + 1/(p-1)：

```
E[MTSHL cost] = O(n²/d · Tg³  +  n²·Tg²  +  n²·d²·Tg  +  nd·Tg²  +  nd⁴)
                   ↑                ↑           ↑           ↑         ↑
                求值（主项）    Vandermonde   评估支撑     插值       固定开销
```

**稀疏主导情形**（Tg << d²）：代价 ~ O(n²/d · Tg³)

**与 Wang 对比**：

| 算法 | 复杂度 | ideal type 1 (αi=0) | ideal type 2 (αi≠0) |
|------|--------|---------------------|---------------------|
| Wang WMDS | O(d₂·d₃·...·dj per step) | 快（稀疏保持稀疏） | **指数爆炸** |
| MTSHL | O(n²/d·Tg³) | 较慢（评估仍需 t 个点） | **多项式** |

### 4.2 Ideal Type 区分

| 类型 | αi 取值 | 实际场景 | 策略 |
|------|---------|---------|------|
| **Type 1** | αi = 0 | 主变量 lc = 1 时可能；稀疏输入有时允许 | Wang 或全零理想专用路径 |
| **Type 2** | αi ≠ 0（随机） | **大多数实际情形**（lc 约束 + 不可分条件迫使非零 α） | MTSHL |

**关键结论**：对工业级 CAS，type 2 是常态，MTSHL 的优势在实践中普遍成立。

### 4.3 MC 2019 基准数据

**稀疏 2 因子（100 项，degree 15）**：

| n | Wang | MTSHL | 加速 |
|---|------|-------|------|
| 5 | 4.87s (89%MDP) | 0.51s | 9.6x |
| 8 | 35.0s (95%MDP) | 0.72s | 49x |
| 12 | 170s (99%MDP) | 1.78s | 95x |
| 14 | 604s (99%MDP) | 2.37s | **255x** |

**稠密循环矩阵行列式 det(Cn)**（验证 MTSHL 对稠密情形也优）：

| n | MTSHL | Wang | 结论 |
|---|-------|------|------|
| 8-10 | 0.14-3.0s | 0.10-1.0s | Wang 略快 |
| **11** | **1.33s** | **12.4s** | MTSHL 开始胜出 |
| 13 | 10.2s | 212s | MTSHL **20x** |
| 14 | 666s | 1364s | MTSHL **2x** |

**结论修正**（vs 架构文档之前的假设）：即使因子稠密，MTSHL 对 n≥11 也更快。原因：Wang MDP 即使 Taylor 系数稠密，其指数性在 n 上仍然体现。

---

## 5. CLPoly 实施路线（更新版）

### 5.1 Phase 1：模 Bézout（P2-Dense-A）— 当前优先项

**目标**：消除 Z[x] Bézout 链系数爆炸（主要瓶颈）

**具体改动**：将 Bézout 链从 Z[x] 移至 Zp[x]（CLPoly 已有 `pa` 参数路径，改为默认）

**理论支撑**：CASC 2018 Algorithm 5 的精简版本（仅 p-adic Step 0：在 Zp 中工作）

**验收标准**：bivar-70 从 65s → <1s；crosscheck 全部通过

### 5.2 Phase 2A：Zippel 稀疏路径（P2-Sparse-D）— 最小交付目标

**目标**：稀疏多项式走完全不同路径（绕开 Hensel 提升和 Bézout）

**参考**：FLINT `factor_zippel.c`（非 MTSHL 路径，纯 Zippel 稀疏插值）

**流程**：
```
1. 在 Zp 上随机求值点 → 单变量因式分解 → 得 r 个模因子
2. 对每个模因子：逐变量 Zippel 稀疏插值恢复高维系数
3. 多个 p 的 CRT + 有理数重建
4. 试除验证
（完全无 Hensel 提升、无 Bézout 链、无 Zassenhaus 重组）
```

**验收标准**：bivar-70 接近 FLINT（<20ms）；crosscheck 通过

### 5.3 Phase 2B：MTSHL-d（P2-Sparse-B）— 可选增强

**目标**：替换 `__multivar_diophantine` 内部的稠密 Bézout 路径为 SparseInt

**参考**：CASC 2018 Algorithm 3+4，CASC 2016 Algorithm 4，JSC 2020 Algorithm 5

**核心组件**：
- `SparseInt`：形式骨架 + Vandermonde 插值 + BDP 回退
- `multi-BDP`：r 因子双变量稠密求解器（C 实现，O(r·d1·d2²)）
- θ-array 批量求值（CASC 2016 §3.3）

**验收标准**：稠密大-r 用例（cyclic matrix n≥11）进一步加速；P2-Sparse-A 完成后按需实施

---

## 6. 各阶段性能预测（更新版）

| 优化阶段 | bivar-70 | trivar-60 | 14var 100term |
|---------|----------|-----------|--------------|
| **当前** | ~65s | ~6s | N/A |
| **Phase 1（模 Bézout）** | <1s | <0.5s | 不变 |
| **Phase 2A（Zippel 路径）** | ~5ms | ~5ms | ~2-3s |
| **Phase 2B（MTSHL-d）** | <1ms | <1ms | ~0.7s |
| **Maple 2019（参考）** | <10ms | <5ms | ~2.4s (n=14) |
| **FLINT（参考）** | ~9ms | ~1.3ms | - |

---

## 7. 实现注意事项

### 7.1 评估点选择（Ideal Type 影响）

- 为满足 lc(f)(α) ≠ 0 和 f(x1,α) 无重根，通常被迫选非零 α（Ideal Type 2）
- MTSHL 在 Type 2 下有压倒性优势；Type 1 (αi=0) 时 Wang 可能更快
- **CLPoly `__select_eval_point`**：当前优先 αi=0（利用稀疏性），这是 Type 1 策略；改为 MTSHL 后不需要强制零理想

### 7.2 BDP 实现（性能关键）

BDP = dense bivariate Diophantine solver in Zp[x1,x2]：
- CASC 2016 指出 BDP 的 C 实现是性能关键
- 代价：O(d1³·d2) —— 需要密集数组实现（非稀疏 pair-vector）
- CLPoly Phase 2B 实现时需要针对 BDP 的特殊优化

### 7.3 σi 项数衰减利用

JSC 2020 的 σi 衰减分析可用于提前终止：
- 当 #σi 降至 1 时，后续所有 σk 都是常数 → 可以提前退出内层循环
- 当 #σi 稳定不再减小时，切换为稠密求解

### 7.4 多因子 Zassenhaus 重组问题

MTSHL-d（CASC 2018）虽然解决了 MDP 问题，但在 r 个 Zp 因子对应 s>1 个真因子时仍需 Zassenhaus 重组。Maple 2019 通过**双变量基底**保证 s=1：
- 双变量基底：评估到 Zp[x1,x2]（保留 x2），而非 Zp[x1]
- 双变量 Zp[x1,x2] 因式分解更能保留因子结构，好求值点下 s=1 以高概率成立
- CLPoly Phase 1/2A 暂不需要双变量基底；Phase 2B（完整 MTSHL）才需要

---

## 8. 参考文献完整列表

| 论文 | 引用 | 对应摘要 |
|------|------|---------|
| Monagan & Tuncer (2016). "Using Sparse Interpolation in Hensel Lifting." *CASC 2016*, LNCS 9890:381–400. | MTSHL 核心思想 + 2 因子 | `CASC2016-MonaganTuncer.md` |
| Monagan & Tuncer (2018a). "Sparse Multivariate Hensel Lifting: A High-Performance Design and Implementation." *ICMS 2018*, LNCS 10931:378–387. | 高性能实现架构 + 并行 | `ICMS2018-MonaganTuncer.md` |
| Monagan & Tuncer (2018b). "Factoring Multivariate Polynomials with Many Factors and Huge Coefficients." *CASC 2018*, LNCS 11077:319–334. | 多因子 + 大系数 | `CASC2018-MonaganTuncer.md` |
| Monagan & Tuncer (2019). "Polynomial Factorization in Maple 2019." *MC 2019* (Extended Abstract). | Maple 2019 集成 + 基准 | `MC2019-MonaganTuncer.md` |
| Monagan & Tuncer (2020). "The Complexity of Sparse Hensel Lifting and Sparse Polynomial Factorization." *J. Symbolic Comput.* 99:189–230. | 完整复杂度分析 | `JSC2020-MonaganTuncer.md` |
| Wang (1975, 1978). Wang EEZ 算法原始文献 | CLPoly 当前实现基础 | - |
| Zippel (1979). "Probabilistic Algorithms for Sparse Polynomials." *EUROSAM '79*:216–226. | Zippel 稀疏插值（FLINT 主路径） | - |
| Geddes, Czapor, Labahn. *Algorithms for Computer Algebra*, Ch.6. | Wang MHL 完整描述 | - |

---

## 9. 结论

CLPoly 多变量因式分解与 FLINT/Maple 的差距有清晰的理论根因和明确的解决路径：

1. **Phase 1（已立项）**：模 Bézout 链 → 消除 Z[x] 系数爆炸 → bivar-70 从 65s → <1s
2. **Phase 2A（最小目标）**：Zippel 稀疏路径 → bivar-70 接近 FLINT（<20ms）
3. **Phase 2B（长期目标）**：MTSHL-d → 与 Maple 2019 对标，n 变量时保持竞争力

MTSHL 系列论文（5 篇，2016–2020）已提供完整的算法描述、复杂度证明和性能数据。实现路线清晰，风险主要在工程复杂度（Phase 2B 约 500–800 行新代码），而非算法设计。
