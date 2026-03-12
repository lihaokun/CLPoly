# CLPoly TODO

## 性能优化

### P1: 单变量因式分解优化

#### P1a: van Hoeij LLL 重组（**已完成**）

**实现结果**（2026-02-22 benchmark）：

| 用例 | 实现前 | 实现后 | 提升 |
|------|--------|--------|------|
| 单变量 70 因子 | 80 ms | 77 ms | 无回退 ✓ |
| Wilkinson W(15) | 13.9 ms | 4.4 ms | 3.2x |
| x^24-1 (cyclotomic) | 24.6 ms | 9.2 ms | 2.7x |

架构文档：`docs/design/vanhoeij/architecture.md`
调研报告：`docs/research/vanhoeij-lll-research.md`

实现任务（全部在 `clpoly/polynomial_factorize_univar.hh`）：

- [x] **M1** `__cld_polys`：CLD 多项式 C_i = (f/h_i)·h_i' mod m
- [x] **M2** `__build_cld_matrix`：内螺旋列喂入格矩阵（初版无列过滤）
- [x] **M3** `__lll_reduce`：整数 LLL（Cohen §2.6，QQ Gram-Schmidt）
- [x] **M4** `__extract_candidates`：U_short 列等价类分组
- [x] **M5** `__vanhoeij_recombine`：主控循环（批量提取 + 对角 LLL 预通道）
- [x] `__lll_factorize`：入口函数（二次 Hensel 提升 + van Hoeij LLL，全 r 路由）

**关键设计**：
- 初始 J_target=0（对角 LLL）免费执行 s=1 单因子检验，避免全分裂多项式触发昂贵 CLD 计算
- 批量提取：单次 LLL 调用中处理所有候选，将 r 次重复 LLL 降为 O(1) 次

#### P1b: 线性 Hensel 提升（调研完成，实现保留为基础设施）

调研与设计（Monagan 2019）已完成：`docs/research/linear-hensel-research.md`、`docs/design/vanhoeij/detailed-design-p1b.md`

- [x] 线性 Hensel 基础设施：`__hensel_step_linear`、`__hensel_lift_linear_recursive`（保留供未来真正线性交织优化）
- [x] `__heuristic_starting_precision`：FLINT 启发式精度估算

**结论**：van Hoeij LLL 需在满精度（Mignotte 界）下运行才能高效收敛；低精度 LLL 导致 O(n^4) 规约退化。P1a 的二次 Hensel + van Hoeij LLL 已达到实用性能，P1b 完整实现（精度自适应 + 早期因子检测）可作为后续优化项。

参考：van Hoeij 2002, Belabas 2009, Hart 2011（P1a）；Monagan 2019（P1b）

**P1b v2：Monagan 2022 矩阵优化**（依赖 P1b 完整实现）

- **Monagan 2022 (ISSAC)**: *"Linear Hensel Lifting for Z_p[x,y] for n Factors with Cubic Cost"* — 将 2019 论文从 2 因子扩展到 n 因子，通过对所有 D 步纠正序列做多项式矩阵乘法，总代价从 O(r·n²·D) 降至 O(n·D²)
- 设计文档：`docs/research/linear-hensel-research.md §6.3`，架构文档 `docs/design/vanhoeij/architecture.md §10`
- **实用价值限制**：早期终止时 D 很小，矩阵乘法优势有限；仅在 D 较大（低精度阶段）时才显著

- [ ] P1b v2：实现 Monagan 2022 多因子矩阵线性 Hensel 提升

#### P1c: Zp DDF/EDF 稠密化（`__select_prime` 瓶颈）

**Profiling 结果**（2026-02-23，`-O2 -DCLPOLY_PROFILE`）：

| 用例 | select_prime | hensel_lift | recombine |
|------|-------------|-------------|-----------|
| ~deg15 (3 fac) | **59.8%** | 31.4% | 8.8% |
| ~deg29 (5 fac) | **53.2%** | 38.5% | 8.4% |
| W(20) | 38.4% | **52.2%** | 9.4% |

`select_prime` 已调为 `max_tries = 3`（与 FLINT 一致）。瓶颈在于 DDF/EDF 使用 CLPoly 通用稀疏 pair-vector 表示 `vector<pair<umonomial, Zp>>`，而 Zp 多项式在小素数下几乎全密，稀疏结构没有任何收益，反而有大量开销。

**根本原因**：
- CLPoly `upolynomial_<Zp>`：`vector<pair<umonomial, Zp>>`，24B/项，内存跳跃，非 SIMD 友好
- FLINT `nmod_poly_t`：`ulong[]` 密集数组，8B/项，顺序访问，Barrett 模约减（~7 cycles vs ~25 cycles）

**候选优化方案（三选一）**：

- **Option A（推荐）：DDF/EDF 内部使用密集 `vector<uint32_t>`**
  - 在 `__ddf_Zp` / `__edf_Zp` / `__upoly_powmod` 入口处将 `upolynomial_<Zp>` 转为密集数组
  - 所有 powmod / GCD 在密集数组上运算，结果再转回
  - 预计收益：3–5x 加速 `select_prime`，整体因式分解加速 ~40–50%
  - 参考：FLINT `nmod_poly_t`，GCL §8.2

- **Option B（最快实现）：在 `__select_prime` 里调用 FLINT `nmod_poly_factor_squarefree`**
  - 直接用已链接的 FLINT 做 DDF+EDF，只需格式转换（pair-vector ↔ nmod_poly_t）
  - 一周内可实现，效果接近理论上限
  - 缺点：引入深层 FLINT 依赖，不利于自主实现目标

- **Option C：BSGS-DDF**（Shoup 1995）
  - O(√d) 次 powmod 代替 O(d) 次，对高次分组有效
  - 对主要用例（d=1 线性因子）无帮助，优先级低

**相关文件**：`clpoly/polynomial_factorize_zp.hh`（`__ddf_Zp`, `__edf_Zp`），`clpoly/polynomial_factorize_univar.hh`（`__select_prime`）

- [ ] 选定方案并实现 Zp DDF/EDF 稠密化

#### P1d: Hensel lift 优化（第二大瓶颈）

W(20) 等高次多项式中 `hensel_lift` 占 52%。这与 **P1b（线性 Hensel 提升）** 直接相关——P1b 的完整实现正是针对此瓶颈的解决方案。

**已有调研**（Monagan 2019 论文，详见 `docs/research/linear-hensel-research.md`）：
- **Monagan 2019 (ISSAC)**: *"Linear Hensel Lifting for F_p[x,y] and Z[x] with Cubic Cost"* — P1b 核心算法，Bézout 系数固定 mod p，每步代价从 O(n²d) → O(n²)（d = 精度 bit 数）
- **Monagan & Tuncer 2019 (CASC)**: *"Polynomial Factorization in Maple 2019"* — Maple 2019 的线性 Hensel + van Hoeij LLL 交织架构，早期终止减少提升步数

**当前状态**：P1b 已实现启发式起始精度（`__heuristic_starting_precision`）和 Phase 2 fallback，但完整的线性 Hensel + LLL 交织循环（`__linear_hensel_lift_with_lll`）尚未实现。详细设计见 `docs/design/vanhoeij/detailed-design-p1b.md`。

- [ ] 实现 P1b 完整线性 Hensel 交织循环（M0–M3，约 247 行，见 linear-hensel-research.md §5）

### P2: 多变量因式分解算法选择与稀疏路径

多变量因式分解当前使用 Wang EEZ + Zassenhaus 重组，多因子时性能极差。

| 用例 | CLPoly | FLINT | 比率 |
|------|--------|-------|------|
| 双变量 70 因子 | 71 s | 10 ms | 7200x |
| 三变量 60 因子 | 6.7 s | 1.3 ms | 5100x |

当前 CLPoly 无论密度一律走 Wang + Zassenhaus。需要：1) 实现稀疏算法路径；2) 根据密度/变量数自动分流。

**稀疏算法候选：**

- **Zippel 稀疏插值**（FLINT 默认）：在多个随机点求值 → 单变量因式分解 → Vandermonde 插值重建多变量因子。完全绕开重组问题，复杂度 O(t³)，t = 项数。参考：Zippel 1979/1993, FLINT `irred_zippel.c`
- **MTSHL 稀疏 Hensel lift**（Maple 2019+）：计算大量双变量像 → 双变量 Hensel lift → 稀疏插值重建多变量因子。避免多变量 Diophantine 方程，高变量数时可能更优。参考：Monagan-Tuncer 2019 (ISSAC)

**算法分流策略**（参考 FLINT/Maple）：稠密 → Wang EEZ，稀疏 → Zippel/MTSHL，回退 → Zassenhaus。

- [ ] 评估 Zippel vs MTSHL 的适用场景
- [ ] 实现选定的稀疏算法作为 Wang EEZ 的备选路径
- [ ] 实现密度计算 + 根据密度/变量数自动选择算法路径

详见：`docs/research/polynomial-systems-comparison.md` §4, `docs/design/multivariate-factorization-design.md`

### P3: 多项式乘法与除法优化

#### P3-mul: NTT 乘法（基础设施）

- [ ] FFT / NTT 大次数单变量乘法
- [ ] 多变量乘法堆排序优化

#### P3-div: Newton 快速除法（依赖 P3-mul）

当前 `_poly_divrem` 是 schoolbook O(n²)，是 HGCD 无法达到理论复杂度的瓶颈。

**原理**：Newton 迭代除法将 divrem 转化为模逆 + 两次乘法，代价 O(M(n))。
- Q = A · rev(B)^{-1} mod x^{deg(Q)+1}（模逆用 Newton 迭代，底层为 NTT 乘法）
- R = A - Q·B

**依赖**：P3-mul（NTT 乘法）。用 Karatsuba 做 Newton 除法可得 O(n^1.585)，用 NTT 可得 O(n log n)。

**预期效果**：

| 组合 | divrem | HGCD 总复杂度 | 预估 α |
|------|--------|--------------|--------|
| 当前（schoolbook divrem） | O(n²) | O(n²) | 2.0 |
| Newton + Karatsuba | O(n^1.585) | O(n^1.585 · log n) | ~1.7 |
| Newton + NTT | O(n log n) | O(n log² n) | ~1.2 |

参考：FLINT `nmod_poly_div_newton`，GCL §3.1，Modern Computer Algebra §9.1

- [ ] 实现 `_poly_inv_series`（多项式模逆，Newton 迭代）
- [ ] 实现 `_poly_divrem_newton`（基于模逆的快速除法）
- [ ] 在 `_hgcd_recursive` 和 `_gcd_hgcd` 中按次数阈值分派 schoolbook vs Newton

#### P3b: HGCD（已完成）

**实现结果**（2026-03-09 benchmark）：

| 用例 | Euclid (改前) | HGCD (现在) | FLINT | 改前 ratio | 现在 ratio |
|------|--------------|-------------|-------|-----------|-----------|
| deg500+common250 | 0.57ms | 0.43ms | 0.40ms | 1.44x | 1.01x |
| deg1000+common500 | 2.27ms | 1.53ms | 0.85ms | 2.67x | 1.31x |
| deg3000+common1500 | 20.0ms | 13.5ms | 4.70ms | 4.28x | 2.36x |
| deg5000+common2500 | 55.5ms | 37.8ms | 10.5ms | 5.28x | 3.27x |
| deg10000+common5000 | 221ms | 148ms | 32.2ms | 6.87x | 4.61x |

HGCD 给出常数因子加速（~1.5x），但有效指数 α ≈ 2.0（与 Euclid 相同），未达到理论 α ≈ 1.6。根因：`_poly_divrem` schoolbook O(n²) 主导了 HGCD 的 O(M(n) log n) 乘法贡献。需 P3-div（Newton 快速除法）解决。

- [x] 三层架构：`gcd` → `_gcd_hgcd` → `_hgcd_recursive` / `_hgcd_iter`
- [x] 矩阵指针稳定化修复（`_hgcd_iter` 后两阶段暂存+拷回）
- [x] 单元测试：`test/test_dense_upoly_hgcd.cc`（10 用例，含 2000 组随机扫描）

### P4: GCD 性能

当前多变量 GCD 比 FLINT 慢 3-9x。

- [ ] 调研 FLINT 的模 GCD 实现
- [ ] 考虑稀疏模 GCD（SPMOD）或 Zippel GCD

## Bug 修复

### B1: 多变量因式分解 over-splitting（**已修复**）

**根因**：Wang 小整数 eval point + MTSHL Zp 提升的架构不兼容。小整数 eval point 高概率导致 over-splitting → MTSHL 产出垃圾 → Zassenhaus trial division 不完全 → `verified.size() >= 2` 无条件返回。

**修复**（方向 2 + 方向 5）：
- 方向 2：`__wang_core` 返回前检查 `mv_T.size() <= 1`，拒绝 over-splitting 的 eval point
- 方向 5：`__select_eval_point` 从 `[-B, B]^{n-1}` 随机采样（暂用 B=50，待 B2 修复后可增大）
- 设计文档：`docs/design/factorization/fix-oversplit.md`

**测试**：B1 回归用例改为正常断言，全量 + 压力测试通过。

### B2: 单变量因式分解大系数不完全分解（**未修复**）

**发现于**：B1 修复过程中，方向 5 使用 B=1000 随机 eval point 暴露。`crosscheck_flint_mpoly_factor_random_bivar` trial 8 和 `crosscheck_flint_mpoly_factor_random_trivar_3fac` trial 2 失败。

**复现**：
```
f = 96000*x^4 + 38464001*x^3 - 61401631999*x^2 - 24616960640000*x - 24555520640000
# 真因子: (x+1)(x-800)(x+800)(96000*x+38368001)
factorize(f) → content=1, 2 factors: (x+1), (96000x^3 + ...)  # 错误：cubic 未继续分解
```

**根因**：`__lll_factorize` (polynomial_factorize_univar.hh:1488) Phase 2 触发条件过严。
- `__select_prime` 正确找到 p=7，r=4 个 Zp 因子
- Phase 1：启发式精度 a_h=6（模 7^6 ≈ 2^17），van Hoeij 返回 **2 因子**（精度不足）
- Phase 2 触发条件 `result.size() == 1 && r > 1`：Phase 1 返回 2 ≠ 1，**未触发**
- Phase 2（Mignotte 精度 2^90）能正确返回 4 因子

**修复方案**：将 Phase 2 触发条件改为 `result.size() < r`（或更精确的条件）。需要注意 x^4+1 等"mod p 多分裂但 Z 上不可约"的情况——对这些多项式 Phase 2 应直接返回单因子而非无限重试。

**B1 修复的临时规避**：方向 5 的 B 暂设为 50（避免大系数触发 B2），待 B2 修复后可增大到 100-1000。

- [ ] 修改 `__lll_factorize` Phase 2 触发条件
- [ ] 补充大系数单变量因式分解回归测试

### 方向 5 进一步测试（待 B2 修复后）

方向 5（`__select_eval_point` 随机化）当前使用 B=50 作为临时安全值。B2 修复后：

- [ ] 将 B 恢复到 `max(100, 2*tdeg_sum)` 或更大
- [ ] 重新运行 crosscheck（`make crosscheck`）验证无回归
- [ ] 运行压力测试（`make stress`）
- [ ] 运行随机复现测试（`test_repro_bivar4fac 100`、`test_repro_trivar3fac 100`）
- [ ] 考虑对 B 做参数化，支持用户配置或根据系数大小自适应

### H1: GCDHEU `(uint64_t)(-v)` 有符号溢出（UB）

`clpoly/polynomial_gcd.cc:57`：当 `v == INT64_MIN` 时，`-v` 在 `int64_t` 域溢出，是 C++ 未定义行为。

**快速修复**：用 `clpoly_zz_to_mpz` 替代手动 small/large 分支，避免有符号取负。

**长期方案**：给 ZZ 补上位操作（`mul_2exp`/`tdiv_q_2exp`/`tdiv_r_2exp`），GCDHEU 直接用 ZZ，参考 FLINT fmpz 接口。

- [ ] 快速修复：GCDHEU 系数加法改用 `clpoly_zz_to_mpz`
- [ ] 长期：ZZ 补充位操作接口（`mul_2exp`, `tdiv_q_2exp`, `tdiv_r_2exp`）

## 文档完善

- [ ] API 文档：为核心头文件添加 Doxygen 注释
- [ ] 算法文档：补充各模块采用的算法描述和参考文献
- [ ] 示例代码：在 `examples/` 目录添加典型用例
- [ ] 性能指南：记录各操作的复杂度和适用场景
