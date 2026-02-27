# 调研报告：多变量因式分解性能优化

> 调研日期：2026-02-23
> 目标：定位 CLPoly 多变量因式分解与 FLINT 7000x 性能差距的根因，确定优化方向
> 参考实现：FLINT `fmpz_mpoly_factor.c`（`factor_zippel.c` / `factor_zassenhaus.c`）；Maple 2019 MTSHL
> 当前代码：`clpoly/polynomial_factorize_wang.hh`（`__wang_core`、`__multivar_hensel_lift`、`__multivar_diophantine`）
>
> **⚠ 后续报告**：`multivar-mtshl-research.md`（2026-02-24）整合了 5 篇 Monagan-Tuncer 论文全貌，包含完整算法描述、复杂度证明和更新后的实施路线。本报告作为根因分析基础仍有效，但算法细节请以新报告为准。

---

## 1. 性能差距实测数据

来源：`benchmarks/2026-02-23-043839.txt`（release build `-O3 -DNDEBUG`）

### 1.1 小规模用例（CLPoly 与 FLINT 接近）

| 用例 | CLPoly | FLINT | 比率 |
|------|--------|-------|------|
| factor x²-y² | 0.038ms | 0.023ms | 1.67x |
| factor bivar deg3*deg3 | 0.293ms | 0.160ms | 1.83x |
| factor bivar deg5*deg5 | 0.834ms | 0.325ms | 2.56x |
| factor trivar known | 0.130ms | 0.109ms | 1.20x |
| factor SymPy f_1 (3var) | 0.645ms | 0.471ms | 1.37x |

**结论**：小规模（低因子数）多变量因式分解，CLPoly 与 FLINT 差距仅 1.2x–2.6x，属数据结构开销（稀疏 pair-vector vs 密集数组）。

### 1.2 大因子数用例（7000x 差距）

| 用例 | CLPoly | FLINT（估计） | 比率 |
|------|--------|-------------|------|
| bivar 70 线性因子（deg70） | **65,263ms** | ~9ms | **~7000x** |
| trivar 60 线性因子（deg60） | **5,891ms** | ~1.3ms | **~4500x** |

`bivar 70 factors = ∏_{i=1}^{35}(x+iy)(x-iy)`：度数 70，**仅 36 项**（极稀疏），但 70 个线性因子。

**核心问题**：因子数（r）大时，CLPoly 指数级退化；FLINT 保持多项式时间。

---

## 2. 根因分析

### 2.1 CLPoly 当前管线（Wang EEZ + Zassenhaus）

```
factorize(f ∈ Z[x₁,...,xₙ])
  → squarefreefactorize
  → __wang_core(g)
        ├─ __select_eval_point      选随机求值点 α（α₂,...,αₙ）
        ├─ factorize(g|_{x₂=α₂,...})  单变量因式分解，得 r 个模因子 v₁,...,vᵣ
        ├─ __wang_leading_coeff     LC 分配
        ├─ __multivar_hensel_lift   多变量 Hensel 提升 v₁,...,vᵣ → G₁,...,Gᵣ
        │       └─ __hensel_lift_one_var × (n-1) 个提升变量
        │               └─ __multivar_diophantine × deg(f, xₖ) 次
        └─ Zassenhaus 子集枚举重组  O(2^r) 最坏情形
```

### 2.2 根因 A：Bézout 链系数爆炸（主要瓶颈）

**位置**：`__multivar_hensel_lift`，lines 835–865：

```cpp
g_acc = v[0];
for (i = 1 to r-1):
    GCD(g_acc, v[i], alpha, beta);   // alpha: 常数，beta: 度 i-1
    for (j = 0 to i-1):
        bezout_s[j] = bezout_s[j] * beta;  // 🔴 累积乘法
    bezout_s[i] = denom * alpha;
    g_acc = g_acc * v[i];
```

**问题**：对 r 个因子构建 Bézout 链（∑ sᵢ·V̂ᵢ = denom，V̂ᵢ = ∏_{j≠i} vⱼ），在 Z[x] 上（非 Z_p[x]）有灾难性的系数增长：

| 量 | r=2 | r=10 | r=70 |
|----|-----|------|------|
| bezout_s[0] 的次数 | 0 | 0+1+...+8=36 | 0+1+...+68=**2346** |
| denom 的位数（估算） | O(1) | O(10 log 10)≈50 bit | O(70 × log 70!)≈**6000+ bit** |

> **denom 的大小**：对于线性因子 vᵢ = x+kᵢ（kᵢ ∈ {1,...,35}），Bézout 链的 denom 本质上是 Vandermonde 行列式 V = ∏_{i<j}(kⱼ-kᵢ)。
> log₂(V) = ∑_{d=1}^{69} (70-d)·log₂(d) ≈ **6000 bits**（约 1800 位十进制数）。

**直接代价**：`__multivar_diophantine` 基本情形（line 499）中：

```cpp
auto si_h = bezout_s[i] * h_upoly;  // degree-2346, 6000-bit coefficients × polynomial
upoly_prem(rem_upoly, si_h, v_factors[i], main_var);  // pseudo-division: 2346 bignum ops
```

- 每次 Diophantine 调用：70 × (2346 bignum 乘法 × 6000-bit 整数) = **~1亿次机器操作**
- bivar-70 共需约 2 次调用（j=1..2 Taylor 阶）
- 总计：约 **2 亿次大整数操作** → 数十秒完全合理

**对比**：FLINT 在 Z_p[x] 上工作，系数有界（< p），无系数爆炸。

### 2.3 根因 B：Zassenhaus O(2^r) 重组（次要瓶颈）

**位置**：`__wang_core`，lines 1175–1222，枚举所有 C(r, s) 个子集。

对 bivar-70（每个 Hensel 因子恰好对应一个真因子），s=1 即全部找到，**Zassenhaus 退化为 O(r) 而非 O(2^r)**。所以对此特定用例，Zassenhaus 不是主因。

但对一般多因子多变量多项式：若真因子是模因子的复杂组合（s > 1 时 r 较大），Zassenhaus 会成为主要瓶颈。

### 2.4 根因 C：算法选择不匹配（架构级问题）

bivar-70 = ∏(x+iy)(x-iy) 是**极稀疏**多项式（36 项 / 71×71=5041 个单项式格点 = **0.7% 密度**）。

- CLPoly：无论稀疏/稠密，一律走 Wang EEZ（稠密 Hensel + Zassenhaus）
- FLINT：检测密度，**稀疏时走 Zippel 插值**，完全绕开 Hensel 提升

Zippel 算法对稀疏输入：复杂度 O(r·T²·n·d)，其中 T=最大因子项数。对 bivar-70（每个因子 2 项：T=2），复杂度 O(70×4×1×1) = O(280) 个求值/插值操作 → 毫秒级。

---

## 3. 参考实现分析

### 3.1 FLINT `fmpz_mpoly_factor`

FLINT 的多变量因式分解调度策略（`fmpz_mpoly_factor.c`）：

```
fmpz_mpoly_factor(f):
  若 f 是稀疏（density < threshold）：
    → factor_zippel(f)       ← Zippel 稀疏插值，在 F_p 上工作
  否则：
    → factor_zassenhaus(f)   ← Wang + 模 Hensel + Zassenhaus
```

**factor_zippel 核心思路**（`factor_zippel.c`）：
1. 选随机素数 p，在 F_p[x₁,...,xₙ] 中因式分解
2. 通过多组求值点构造 Vandermonde 插值系统
3. 对每个因子的每个系数做稀疏插值（Prony/Ben-Or-Tiwari）
4. 在 Z 上用 CRT + 有理数重建恢复真因子
5. 全程无 Hensel 提升、无 Bézout 链、无 Zassenhaus 重组

**factor_zassenhaus 核心改进**（vs CLPoly）：
- Bézout 链在 F_{p^a}[x] 上构建（非 Z[x]），系数有界
- Zassenhaus 重组结合 van Hoeij 思路（当 r 大时切换）

### 3.2 Maple 2019 双框架

Maple 2019 的多变量因式分解（Monagan-Tuncer, MC 2019 Extended Abstract）使用 MTSHL 算法：

```
Maple factor(f ∈ Z[x₁,...,xₙ]):

  1. 选择求值点 α₂,...,αₙ（随机，满足 squarefree + coprime 条件）
  2. 单变量因式分解：a₁ = f(x₁,α₂,...,αₙ) → f₁·g₁ ∈ Z[x₁]
  3. 逐变量 Hensel 提升（MTSHL-d）：
     - 每步用 SparseInt（稀疏插值 + Vandermonde）求解 MDP
     - 骨架收缩（Theorem 1）保证求解高效
  4. p-adic 提升恢复 Z 系数
  5. Trial division 验证
```

> **勘误**：此前文档声称 Maple 使用"双变量基底 Zp[x₁,x₂]"（保留一个额外变量做双变量因式分解），
> 经核实 Monagan & Tuncer 全部 5 篇论文（CASC 2016/2018, ICMS 2018, MC 2019, JSC 2020），
> **均使用单变量基底 Zp[x₁]**。MC 2019 原文明确写道：
> "MHL first chooses integers α₂,...,αₙ ... and factors the univariate image a₁ in Z[x₁]."
> ICMS 2018 的"双变量归约"是指 MDP 求解策略（Algorithm 3 将多变量 MDP 归约为
> 双变量 HenselLift1 调用），而非使用双变量因式分解作为基底。该错误推断已修正。

**MTSHL 复杂度**（Monagan-Tuncer 2020, JSC 99:189-230）：
- 稀疏路径（T 项）：O(r·n·d·T²) — 多项式时间
- 当前 CLPoly（单变量基底 + Z[x] Bézout）：Bézout 链 O(r²) 次数、O(r log r) bit 系数 → 实际指数退化

### 3.3 NTL `ZZXFactoring`

NTL 仅支持单变量因式分解，无多变量实现（此对比不适用）。

---

## 4. 算法选型对比

| 算法 | 原理 | 实现难度 | 适用场景 | 预期收益 |
|------|------|----------|----------|----------|
| **Zippel 稀疏插值** | F_p 上求值+Vandermonde 插值 | 高（需要稀疏插值基础设施） | 稀疏多项式（T << d^n） | bivar-70: ~7000x → <2x |
| **模 Diophantine 修复** | 将 Bézout 链移到 F_{p^a} | 中（局部修改 `__multivar_diophantine`） | 所有情形 | bivar-70: ~7000x → ~100x（Zassenhaus s=1 仍可接受） |
| **MTSHL** | 稀疏 Diophantine（Zippel in-Hensel） | 高（重写 Diophantine 为 SparseInt） | 稀疏多项式，保留 Hensel 框架 | 类似 Zippel，实现更复杂 |

**关键结论：优先级排序**

1. **系数爆炸**（根因 A）必须修复，否则大 r 情形不可用——Phase 1 解决
2. **算法选择**（根因 C）需要加 density-based 分流——Phase 2 解决
3. **Zassenhaus O(2^r)**（根因 B）对大 r 是潜在瓶颈，但 Phase 2 Zippel 路径对稀疏输入会完全绕开（非替换），稠密大-r 场景暂不处理

---

## 5. 推荐优化路线

### Phase 1（优先，中等难度）：模 Bézout 链改造

**目标**：修复系数爆炸，Zassenhaus 重组保持不变。

**具体改动**：

```
P2a: __multivar_hensel_lift 改为模提升
  - Bézout 链在 F_{p^a}[x] 上构建（类似单变量 __hensel_lift 的 p-adic 路线）
  - __multivar_diophantine 改用模运算（已有 pa 参数支持！）
  - Zassenhaus 重组保持不变（bivar-70 中 s=1，不构成瓶颈）
```

**说明**：van Hoeij LLL 的 CLD 方法是单变量专有的（CLD = f'/f 只在单变量情形下定义），
没有文献提出将其直接推广到 Z[x₁,...,xₙ]。Maple/MTSHL 使用单变量基底 Zp[x₁]，
随机求值点使 s=1 以高概率成立，无需 LLL 重组。此阶段 Zassenhaus 保持不变；
若未来出现 s>1 爆炸场景，可在稀疏路径中一并解决。

**预期收益**：
- bivar-70: 65s → <1s（系数爆炸消除 → 约 1000x 改善；bivar-70 中 Zassenhaus s=1 足矣）
- 稠密通用情形（r 小）：同比例加速，Zassenhaus O(2^r) 对小 r 无影响

**已有基础**：
- `pa` 参数已在 `__multivar_diophantine` 中支持（见 lines 518-529），但当前仅在 `denom ≠ ±1` 时激活，改造目标是令其成为默认路径

### Phase 2（长期，高难度）：Zippel 稀疏插值路径

**目标**：对稀疏多项式走完全不同的算法路径（绕开 Hensel 提升）。

**具体改动**：
```
新增: __factor_zippel(f, p)   ← 全新模块，约 500-800 行
  - 稀疏 Ben-Or-Tiwari/Prony 插值
  - 在 F_p 上工作，CRT 重建
  - 需要稀疏插值基础设施（目前 CLPoly 没有）

修改: __wang_core 入口：
  density = num_terms / (deg_bound_x × deg_bound_y × ... )
  if density < THRESHOLD: __factor_zippel(f)
  else: 现有 Wang EEZ
```

**预期收益**：
- bivar-70: 65s → ~几ms（完全对标 FLINT）
- 稀疏多变量通用情形：O(r·T²·n·d) 多项式时间

**实施难度**：高——需要新建稀疏插值子系统，与现有代码几乎完全独立。

---

## 6. 各阶段性能预测

| 优化 | bivar-70 | trivar-60 | bivar deg5*deg5 |
|------|----------|-----------|-----------------|
| 当前 | 65s | 5.9s | 0.83ms |
| Phase 1（模 Bézout） | <1s | <0.5s | ~0.4ms |
| Phase 2（Zippel 路径） | ~5ms | ~5ms | ~0.4ms |
| FLINT（参考） | ~9ms | ~1.3ms | ~0.33ms |

> **注**：Phase 1 的估算基于：消除 6000-bit bignum 运算后，Bézout 链在 F_p 上约 64-bit 运算，减少约 1000x；加之只有约 2 次 Diophantine 调用。Phase 2 的估算基于 FLINT Zippel 实测数据类比。

---

## 7. 推荐实施顺序与关键决策

### 决策 1：Phase 1 是否值得做？

**是**。Phase 1（模 Bézout 链）改动局部（主要在 `__multivar_diophantine` 和 Bézout 链构建），风险低，收益巨大：bivar-70 从 65s 降至 <1s。相比 Phase 2（完全重写），性价比极高。

### 决策 2：~~双变量基底是否要做？~~ [已废弃]

> **勘误**：经核实论文原文，Monagan & Tuncer 全部论文均使用**单变量基底 Zp[x₁]**，
> 不存在"双变量基底"算法。此决策基于错误推断，已废弃。
> CLPoly 的单变量基底与 Maple MTSHL 完全一致。

### 决策 3：Zippel 是否现在做？

**暂缓**。实现难度高（~500-800 行新代码 + 稀疏插值子系统），且 Phase 1 已能解决大部分瓶颈。Zippel 的优势在于超稀疏情形（T << d^n），Phase 1 已将 bivar-70 降至可接受范围。可在 Phase 1 稳定后，通过 profiling 确认剩余差距后再决定是否推进。

### 决策 4：MTSHL 与 Phase 1 的关系

MTSHL 本质上是 Hensel 框架 + SparseInt Diophantine（稀疏插值求解 MDP）的组合。
CLPoly M5 已实现完整 MTSHL-d（CASC 2018 全部 5 个算法），M6 补齐 p-adic 提升和 Zp 64-bit。
与 Maple MTSHL 算法完全一致（均使用单变量基底 Zp[x₁]）。

---

## 8. 关键参考文献

| 文献 | 角色 |
|------|------|
| Wang (1978). *EUROSAM '78*:137–151 | Wang EEZ 算法（CLPoly 当前实现基础） |
| Zippel (1979). *EUROSAM '79*:216–226 | 稀疏插值算法（FLINT 主路径） |
| Zippel (1993). *Effective Polynomial Computation* §7 | Zippel 算法完整描述 |
| Monagan & Tuncer (2016). *CASC 2016*, LNCS 9890:171–186 | MTSHL 核心算法 |
| Monagan & Tuncer (2018). *ICMS 2018*, LNCS 10931:378–387 | MTSHL 实现（Maple 采用版本） |
| Monagan & Tuncer (2019). *ISSAC 2019* | Maple 2019 算法全貌 |
| Monagan & Tuncer (2020). *J. Symbolic Comput.* 99:189–230 | MTSHL 复杂度证明 |
| GCL §6.3–§6.7 | 多变量 Hensel 提升与 Diophantine（当前实现参考） |

---

## 9. 结论

CLPoly 多变量因式分解 7000x 差距的**主要根因**是：

1. **Bézout 链系数爆炸**（根因 A，占 ~99% 的时间）：在 Z[x] 上为 r 个因子构建 Bézout 链，bezout_s[0] 的次数达 O(r²)，系数位数达 O(r log r)，导致 Diophantine 基本情形每次调用执行亿次大整数操作。**修复**：将 Bézout 链和 Diophantine 移至 F_{p^a}（已有 `pa` 参数骨架）。

2. **算法选择不匹配**（根因 C）：对稀疏输入（bivar-70 密度 0.7%），CLPoly 仍走稠密 Wang EEZ，FLINT 走 Zippel 稀疏插值。**修复**：Phase 2 实现 Zippel + 密度分流。

**下一步**：架构文档已完成（`docs/design/multivar-opt/architecture.md`），按 workflow 进入细化 → 实施阶段，从 Phase 1 M1（模 Bézout 链）开始。
