# CLPoly 性能现状报告

> 日期：2026-03-04
> 分支：`feature/gcd-optimization`
> 基准数据：`benchmarks/2026-03-04-020307.txt`
> 对比基线：`benchmarks/2026-02-23-043839.txt`（优化前）、`benchmarks/2026-02-28-032628.txt`（中期）

---

## 1. 已完成优化总结

### 1.1 单变量因式分解

| 优化项 | 阶段 | 效果 |
|--------|------|------|
| van Hoeij LLL 重组 (P1a) | 已完成 | recombine 仅占 1-15%，无显著加速 |
| 启发式起始精度 (P1b) | 已完成 | uni-70: 21ms → 13ms (38% 加速) |

### 1.2 GCD

| 优化项 | 阶段 | 效果 |
|--------|------|------|
| M0.1 ZZ 接口扩展 | 已完成 | `is_small()`, `get_val()`, `get_mpz_v()` |
| M0.2 素数生成组件 | 已完成 | `next_prime_64()`, `prev_prime_64()` |
| M0.3 boost::math::prime 迁移 | 已完成 | 完全移除 boost 素数依赖 |
| M0.4 Zp 大素数安全修复 | 已完成 | 支持全范围 uint64_t 素数 |
| M0.5 MTSHL 素数选择修复 | 已完成 | 从 2^64-59 向后选取 |
| P0 大素数 GCD | 已完成 | 单/多变量均从 2^64-59 开始 |

### 1.3 多变量因式分解 Bug 修复

| 修复项 | 效果 |
|--------|------|
| Wang LCC valuation extraction | GCD 整块匹配 → per-factor valuation 提取 |
| Wang LCC non-divisor | coprime 检查 → GCL §8.7 non-divisor 算法，修复无限循环 |
| MTSHL 内循环 LC 校正 | 删除冗余 no-op |

---

## 2. 当前性能数据

### 2.1 单变量 GCD — CLPoly vs NTL

| 用例 | CLPoly | NTL | ratio | 结论 |
|------|--------|-----|-------|------|
| gcd deg50+common25 | 0.031ms | 0.083ms | **0.37x** | ✅ 已反超 |
| gcd deg200+common100 | 0.270ms | 0.569ms | **0.47x** | ✅ 已反超 |

**优化前（02-23）**: gcd deg200+common100 = 3.849ms → **现在 0.270ms（14.3x 加速）**

### 2.2 多变量 GCD — CLPoly vs FLINT

| 用例 | CLPoly | FLINT | ratio | 结论 |
|------|--------|-------|-------|------|
| gcd deg8+common4 (2var) | 0.385ms | 0.064ms | 6.06x | ❌ 仍有差距 |
| gcd deg15+common8 (2var) | 1.585ms | 0.159ms | 9.95x | ❌ 仍有差距 |

**P0 对多变量 GCD 基本无效**——瓶颈不在 CRT 素数个数，而在递归求值 + Lagrange 插值。

### 2.3 单变量因式分解 — CLPoly vs FLINT

| 用例 | CLPoly | FLINT | ratio |
|------|--------|-------|-------|
| ~deg15 (3 fac) | 0.631ms | 0.290ms | **2.17x** |
| ~deg21 (4 fac) | 3.434ms | 0.317ms | **10.83x** |
| ~deg29 (5 fac) | 10.672ms | 0.881ms | **12.11x** |
| Wilkinson W(10) | 0.288ms | 0.128ms | **2.25x** |
| Wilkinson W(15) | 0.675ms | 0.313ms | **2.16x** |
| Wilkinson W(20) | 1.074ms | 0.417ms | **2.57x** |
| Wilkinson W(25) | 1.579ms | 0.681ms | **2.32x** |
| uni 70 factors (stress) | 13.512ms | 3.423ms | **3.95x** |

注：~deg21/~deg29 为随机生成多项式，波动大（历史范围 1.2-14.7ms），Wilkinson 系列更稳定。

### 2.4 单变量因式分解 — CLPoly vs NTL

| 用例 | CLPoly | NTL | ratio |
|------|--------|-----|-------|
| ~deg15 (3 fac) | 0.520ms | 0.462ms | 1.12x |
| ~deg21 (4 fac) | 1.126ms | 0.795ms | 1.42x |
| ~deg29 (5 fac) | 5.704ms | 1.310ms | **4.35x** |
| Wilkinson W(10) | 0.287ms | 0.428ms | **0.67x** ✅ |
| Wilkinson W(15) | 0.672ms | 0.549ms | 1.22x |
| Wilkinson W(20) | 1.079ms | 0.896ms | 1.20x |
| x^15-1 | 0.313ms | 0.351ms | **0.89x** ✅ |
| Swinnerton-Dyer S3 | 0.348ms | 0.313ms | 1.11x |

经典多项式（Wilkinson/cyclotomic/SD）基本持平或反超 NTL；随机多因子多项式有 1.4-4.4x 差距。

### 2.5 多变量因式分解 — CLPoly vs FLINT

| 用例 | CLPoly | FLINT | ratio |
|------|--------|-------|-------|
| bivar deg3*deg3 | 0.395ms | 0.128ms | 3.07x |
| bivar deg5*deg5 | 0.958ms | 0.355ms | 2.70x |
| trivar known | 0.150ms | 0.044ms | 3.38x |
| SymPy f_1 (3var) | 0.490ms | 0.418ms | 1.17x |
| bivar 70 fac (stress) | 345.8ms | 9.5ms | **36.4x** |
| trivar 60 fac (stress) | 184.1ms | 1.1ms | **166.3x** |
| 10var disjoint | 1.601ms | 0.377ms | 4.25x |

### 2.6 其他算子 — CLPoly vs FLINT/NTL

| 算子 | CLPoly vs FLINT | CLPoly vs NTL | 结论 |
|------|----------------|---------------|------|
| 多变量 add | 1.1-2.1x | — | 接近 |
| 多变量 mul | 1.4-3.7x | — | 小差距 |
| 单变量 mul deg500 | — | 6.8x | ❌ 差距大 |
| 单变量 eval deg500 | — | 3.4x | 中等差距 |
| 单变量 add/deriv | — | 0.2-0.5x | ✅ 已反超 |

---

## 3. 优化前后对比（关键里程碑）

| 用例 | 02-23 (基线) | 03-04 (当前) | 加速 |
|------|-------------|-------------|------|
| uni gcd deg200+common100 | 3.849ms | 0.270ms | **14.3x** |
| uni gcd deg80+common40 | 0.594ms | 0.124ms | **4.8x** |
| bivar 70 fac (stress) | 65263ms | 345.8ms | **189x** |
| trivar 60 fac (stress) | 5891ms | 184.1ms | **32x** |
| uni 70 fac (stress) | 12.4ms | 13.5ms | ≈1x |

---

## 4. 差距分析与瓶颈定位

### 4.1 单变量因式分解（vs FLINT 2-4x）

已知瓶颈（profiling 数据）：
- **`__select_prime`**: 占 77-94% 的总时间。CLPoly 尝试 30-50 个素数 × 每个素数做完整 Zp DDF+EDF。FLINT 的 `_fmpz_poly_factor_zassenhaus` 使用更高效的素数筛选策略。
- 单变量 Zp GCD（朴素欧几里德）是 `__select_prime` 内的核心操作。

GCD 优化路线图中的 P1 (GCDHEU) 和 P2 (CRT 优化) 可能间接帮助（减少 Zp 操作开销），但直接优化 `__select_prime` 本身（减少素数尝试次数、更快的无平方判定）可能更有效。

### 4.2 多变量 GCD（vs FLINT 6-10x）

P0 大素数对多变量几乎无效。瓶颈在递归求值 + Lagrange 插值结构：
- CLPoly 用 Brown/Zippel 稠密插值（需 deg+1 个求值点）
- FLINT 有多种策略选择（Brown/Zippel/BFSS/Hensel），根据稀疏度自动切换
- 可能的优化：Zippel 稀疏插值、更高效的求值点选择

### 4.3 多变量因式分解（stress 36-166x）

大量因子的 stress test 差距极大，主要原因：
- CLPoly 的 Wang/EEZ 在因子数多时需大量求值点重试
- FLINT 使用 Zassenhaus + multivariate Hensel，且底层算术更快
- 单个小多项式（bivar deg5*deg5 = 2.7x）差距小，说明算法框架本身合理，瓶颈在底层效率

### 4.4 单变量乘法（vs NTL 6.8x at deg500）

CLPoly 使用朴素 O(n²) 乘法，NTL 在高次使用 Karatsuba/FFT。这影响因式分解中所有涉及多项式乘法的步骤。

---

## 5. 待实施优化（按优先级）

### GCD 路线图（architecture.md 已规划）

| 优先级 | 模块 | 预期效果 | 适用范围 |
|--------|------|---------|---------|
| 1 | P1 GCDHEU | 小系数 2-5x | 单变量 |
| 2 | P2 CRT 优化 | 1.2-1.5x | 单/多变量 |
| 3 | P3 HGCD | 高次 (deg>1725) | 单变量 Zp |

### 因式分解（未规划）

| 优先级 | 方向 | 预期效果 | 说明 |
|--------|------|---------|------|
| 高 | 优化 `__select_prime` | 2-4x | 减少素数尝试次数、更快的 Zp 操作 |
| 中 | 单变量快速乘法 | 间接受益 | Karatsuba/FFT，影响所有高次操作 |
| 低 | 多变量 GCD 策略选择 | 多变量 2-5x | Zippel 稀疏插值 / 策略自动选择 |
