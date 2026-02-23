# 多变量因式分解优化：架构文档

> 阶段：架构（workflow.md §2.2）
> 基础调研：`docs/research/multivar-factorization-research.md`
> 目标：复现 Maple 2019 双路径架构（Dense Wang + Sparse MTSHL）
> 目标文件：`clpoly/polynomial_factorize_wang.hh`

---

## 总体架构

```
factorize(f ∈ Z[x₁,...,xₙ])
  → squarefreefactorize
  → __wang_core(g)
        ├─ 【Phase 2 新】density = num_terms(g) / prod(deg(g,xᵢ)+1)
        ├─ density < SPARSE_THRESHOLD
        │     → 【Phase 2】__factor_zippel(g)         ← 全新稀疏路径
        └─ density ≥ SPARSE_THRESHOLD
              → 【已有+Phase 1 改造】Wang EEZ 稠密路径
                    ├─ __select_eval_point
                    ├─ factorize(univar)（单变量基底，保持不变）
                    ├─ __wang_leading_coeff
                    ├─ 【Phase 1】__multivar_hensel_lift（模 Bézout 版）
                    └─ Zassenhaus 重组（保持不变）
```

**实施顺序**：
- **Phase 1（优先）**：改造稠密路径——模 Bézout 链修复系数爆炸，Zassenhaus 重组保持不变
- **Phase 2（后续）**：实现稀疏路径——Zippel 稀疏插值 + 密度分流

---

## Phase 1：稠密路径改造

> 对应调研 §5 "Phase 1（优先，中等难度）"
> 改动范围：仅 P2-Dense-A（模 Bézout 链），Zassenhaus 重组保持不变
> 预期收益：bivar-70 从 65s → <1s（系数爆炸消除 ~1000x）

### 关于评估基底（Evaluation Base）

CLPoly 当前使用**单变量基底**：`assign(g, eval)` 将所有非主变量代入，调用 `factorize(upolynomial)`。
Maple 使用**双变量基底**：保留一个额外变量，内层做双变量 Z_p[x₁,x₂] 因式分解。

双变量基底带来 Zassenhaus s=1 的保证（因子数不分裂），是 Maple MTSHL 的基础。
但 **Phase 1（模 Bézout）在单变量基底下即可工作**，无需改变基底。
双变量基底为**未来工作**（仅在实现完整 Maple MTSHL 稠密路径时才需要）。

### P2-Dense-A：模 Bézout 链（修复系数爆炸）

**问题**：`__multivar_hensel_lift` 中 Bézout 链在 Z[x] 上构建（lines 836–865），bezout_s[0] 次数达 O(r²)，系数达 6000+ bit，导致每次 Diophantine 执行亿次大整数运算。

**方案**：将 Bézout 链移至 F_p[x]，在 Z_{p^a} 上做模算术。

```
改前（Z[x]）：
  g_acc, bezout_s[i] 都是 ZZ 系数多项式
  denom 可达 ~2^6000（Vandermonde 行列式量级）

改后（F_{p^a}[x]）：
  1. 选素数 p（与所有 lc(vᵢ) 互素，与 denom 互素）
  2. 选精度 a 使 p^a > 2·B（B = 系数界）
  3. 在 F_p[x] 上构建 Bézout 链（系数 ∈ {0,...,p-1}，无爆炸）
  4. __multivar_diophantine 的基本情形全程用模运算
```

**代码改动位置**：
- `__multivar_hensel_lift`（line 836–865）：Bézout 链构建改为 F_p[x] 版本
- `__multivar_diophantine` 基本情形（line 490–544）：已有 `pa` 参数路径，改为强制走此路径
- 新增辅助：`__bezout_chain_modular(factors, p)`

**接口**（保持外部签名不变）：

```cpp
// 改动前后签名相同，内部改为模运算
std::vector<Poly>
__multivar_hensel_lift(
    const Poly& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::vector<Poly>& lc_targets,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

---

### Phase 1 模块依赖图

```
P2-Dense-A（唯一 Phase 1 改动）
────────────────────────────────────────────
__bezout_chain_modular  [新：F_p[x] 上构建 Bézout 链]
        ↓
__multivar_hensel_lift  [改：Bézout 链改为模版本]
        ↓
__multivar_diophantine  [改：强制走 pa 模路径]
（已有 pa 路径，改为默认启用）

Zassenhaus 重组：保持不变（bivar-70 中 s=1，不构成瓶颈）
```

### Phase 1 预期性能

| 用例 | 当前 | Phase 1 后 | 说明 |
|------|------|-----------|------|
| bivar 70 factors | 65s | <1s | 系数爆炸消除（~1000x），Zassenhaus s=1 |
| trivar 60 factors | 5.9s | <0.5s | 类似 |
| bivar deg5*deg5 | 0.83ms | ~0.5ms | r 小，改善有限 |

---

## Phase 2：稀疏路径（MTSHL）

> 对应调研 §5 "Phase 2（长期，高难度）"
> 依赖：Phase 1 完成且稳定
> 预期收益：bivar-70 <1s → ~5ms（对标 FLINT ~9ms）

### P2-Sparse-A：Zippel 稀疏插值基础设施

多变量稀疏多项式的基本构件，MTSHL 的核心依赖。

```
给定：f ∈ Z_p[x₁,...,xₙ] 的若干点值 (β_j, f(β_j))
输出：f 的稀疏表示 Σ cₐ x^a（仅非零项）

算法（Zippel 1979 / Ben-Or-Tiwari 1988）：
  1. 固定 x₂,...,xₙ = r₂,...,rₙ（随机），对 x₁ 做单变量插值 → 得候选项
  2. 逐变量归纳：固定更少变量，确认每个候选项的各变量次数
  3. 线性方程组求解（Vandermonde 求逆）确定系数
```

**新增模块**：

```cpp
// 稀疏插值：从点值重建稀疏多项式
// samples[j] = (eval_point_j, value_j)
// deg_bound[i] = 第 i 个变量的次数上界
// T_bound = 非零项数上界（可从输入 f 的项数估计）
// 返回：稀疏多项式（pair-vector 格式）
Poly
__zippel_interpolate(
    const std::vector<std::pair<std::map<variable, ZZ>, ZZ>>& samples,
    const std::map<variable, int>& deg_bound,
    int T_bound,
    uint32_t p);
```

---

### P2-Sparse-B：MTSHL Diophantine 求解器

替换 `__multivar_diophantine` 的稠密 Bézout 路径，用稀疏插值求解 Diophantine 方程。

```
当前（稠密 Bézout）：
  解 h = Σ δᵢ·Ĝᵢ → 全稠密 pseudo-division → O(r·d^k) per call

MTSHL 稀疏路径：
  1. Probe eval：在随机点 β 处求值 h(β), Ĝᵢ(β)
     解小线性方程组 → 得 δᵢ(β)（r 个标量）
  2. 取足够多 β 点（T_δ 个，T_δ = δᵢ 的项数上界）
  3. Zippel 插值 → 重建稀疏 δᵢ
  4. 验证：h ≈ Σ δᵢ·Ĝᵢ（mod p）
  → 代价 O(T_δ·n·d) per call（vs 稠密 O(r·d^k)）
```

**新增函数**：

```cpp
// MTSHL 稀疏 Diophantine 求解
// 同 __multivar_diophantine 接口，内部用稀疏插值
std::vector<Poly>
__multivar_diophantine_sparse(
    const Poly& h,
    const std::vector<Poly>& G,
    const variable& main_var,
    const std::vector<std::pair<variable, ZZ>>& eval_vars,
    uint32_t p, const ZZ& pa,
    int T_bound);
```

---

### P2-Sparse-C：密度估计 + 路径分流

**新增函数**：

```cpp
// 计算多项式密度 = 实际项数 / 全稠密项数上界
// density ∈ (0, 1]
double __poly_density(const Poly& f);

// 分流阈值（待 profiling 调整，初始参考 FLINT：0.5）
constexpr double SPARSE_THRESHOLD = 0.5;
```

**`__wang_core` 改动**（新增入口判断）：

```cpp
// 在现有 __wang_core 顶部插入：
double density = __poly_density(g);
if (density < SPARSE_THRESHOLD)
    return __factor_zippel_core(g);
// else: 现有稠密路径（已经 Phase 1 改造过）
```

---

### P2-Sparse-D：Zippel 主控（`__factor_zippel_core`）

稀疏因式分解的顶层控制流，与 Wang 路径并列：

```cpp
// 稀疏路径顶层，对应 FLINT factor_zippel.c 的逻辑
// 1. 在多个 F_p 下工作（选好的 p）
// 2. 对每个 p：随机求值点 → 单变量分解 → 稀疏插值重建因子
// 3. 多个 p 的结果 CRT 重建 → 有理数重构
// 4. 试除验证
std::vector<std::pair<Poly, uint64_t>>
__factor_zippel_core(const Poly& f);
```

---

### Phase 2 模块依赖图

```
P2-Sparse-A（Zippel 插值基础设施）
        ↓
P2-Sparse-B（MTSHL Diophantine）   P2-Sparse-D（Zippel 主控）
        ↓                                   ↓
        └──────────────┬────────────────────┘
                       ↓
              P2-Sparse-C（密度分流）
                       ↓
               __wang_core（分流入口）

注：P2-Sparse-B 和 P2-Sparse-D 可并行开发（不互相依赖）
    P2-Sparse-C 最后集成
```

---

## 与现有代码的关系

| 现有函数 | Phase 1 改动 | Phase 2 改动 |
|---------|------------|------------|
| `__multivar_diophantine` | 强制走 `pa` 模路径 | 新增 `_sparse` 变体 |
| `__multivar_hensel_lift` | Bézout 链改为模版本 | 接口不变，内部可选稀疏 Diophantine |
| `__hensel_lift_one_var` | 不变 | 不变（接口稳定） |
| `__wang_core` | 不变（Zassenhaus 保持） | 新增 Zippel 分流入口 |
| `__vanhoeij_recombine` | 不变 | 不变 |
| `__wang_leading_coeff` | 不变 | 不变 |
| `__select_eval_point` | 不变 | 不变 |

**新增文件**：无（所有改动在 `polynomial_factorize_wang.hh` 内）

---

## 实施里程碑

### Phase 1 里程碑

```
M1: 模 Bézout 链（P2-Dense-A）
    验收：bivar-70 从 65s → <1s；trivar-60 从 5.9s → <0.5s；crosscheck 通过
    测试：make crosscheck + make stress + make bench-all
```

### Phase 2 里程碑

```
M3: Zippel 基础设施（P2-Sparse-A）
    验收：单元测试通过（给定点值能正确重建稀疏多项式）

M4: MTSHL Diophantine（P2-Sparse-B）
    验收：crosscheck 通过，对稀疏用例加速可观

M5: Zippel 主控 + 密度分流（P2-Sparse-C/D）
    验收：bivar-70 接近 FLINT（<20ms）；稠密用例不退化
    测试：make bench-all（对比 2026-02-23 基线）
```

---

## 参考

| 文献 | 对应模块 |
|------|---------|
| GCL §6.3–§6.7 | 当前 Wang 实现基础；模 Bézout 的理论依据 |
| Zippel (1979/1993) | P2-Sparse-A/D |
| Monagan-Tuncer (2016, CASC) | P2-Sparse-B（MTSHL Diophantine）|
| Monagan-Tuncer (2018, ICMS) | P2-Sparse-B/D（实现细节）|
| Monagan-Tuncer (2019, ISSAC) | 双路径架构总览 |
| Monagan-Tuncer (2020, JSC) | 复杂度证明（验收标准参考）|
| FLINT `fmpz_mpoly_factor/` | P2-Sparse-D（开源参考实现）|
