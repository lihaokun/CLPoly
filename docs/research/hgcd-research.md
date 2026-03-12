# P3 调研：HGCD (Half-GCD) 算法用于 Zp 多项式 GCD

> 日期：2026-03-05（初版），2026-03-08（更新）
> 分支：`feature/gcd-optimization-p1`
> 上下文：P3a（dense_upoly_zp）和 P3-mul（Karatsuba + lazy divrem）已完成。P3-HGCD 是下一步优化方向。

---

## 1. 问题背景

### 1.1 当前瓶颈

CLPoly 的单变量 ZZ GCD 使用模 GCD 框架（CRT 循环），每轮调用 Zp 多项式 Euclid GCD：

```
polynomial_GCD(F, G)  // ZZ[x]
  while CRT 未稳定:
    __polynomial_GCD(Pout_mod, f_p, g_p, lc_gcd_p, Pout_d)  // Zp[x], 稠密 Euclid
```

P3a + P3-mul 完成后，Zp Euclid GCD 已接近 FLINT 的 Euclid base case 性能。但 FLINT 在 deg >= 340 时切换到 HGCD（O(M(n) log n)），而 CLPoly 仍使用 O(n^2) Euclid，高次差距随 degree 增长。

### 1.2 性能数据（2026-03-08，P3a + P3-mul 完成后）

**CLPoly Euclid vs FLINT HGCD（FLINT 在 deg >= 340 自动使用 HGCD）**：

| 用例 | CLPoly | FLINT | ratio |
|------|--------|-------|-------|
| gcd Zp deg200+common100 | 0.066ms | 0.053ms | 1.25x |
| gcd Zp deg500+common250 | 0.571ms | 0.397ms | 1.44x |
| gcd Zp deg1000+common500 | 2.27ms | 0.85ms | 2.67x |
| gcd Zp deg3000+common1500 | 20.0ms | 4.70ms | 4.28x |
| gcd Zp deg5000+common2500 | 55.5ms | 10.5ms | 5.28x |
| gcd Zp deg10000+common5000 | 221ms | 32.2ms | 6.87x |

差距从 deg500 的 1.44x 增长到 deg10000 的 6.87x —— 这正是 O(n^2) vs O(M(n) log n) 的渐近差异。

---

## 2. HGCD 算法概述

### 2.1 核心思想

Half-GCD（半 GCD）是一种分治 GCD 算法：

1. 只处理输入多项式的**高半部分**（前 n/2 项系数），递归计算变换矩阵
2. 用矩阵乘法将变换应用到**完整多项式**
3. 递归合并

类比：Euclid GCD 之于 HGCD，如同选择排序之于归并排序。

### 2.2 复杂度

| 算法 | 时间复杂度 | 条件 |
|------|-----------|------|
| Euclid GCD | O(n^2) | -- |
| HGCD | O(M(n) log n) | M(n) = 多项式乘法复杂度 |

其中 M(n) 取决于乘法算法：

| 乘法算法 | M(n) | HGCD 复杂度 |
|---------|------|------------|
| Schoolbook | O(n^2) | **O(n^2 log n)** -- 比 Euclid 更差！ |
| Karatsuba | O(n^1.585) | O(n^1.585 log n) |
| NTT/FFT | O(n log n) | O(n log^2 n) |

**关键结论**：HGCD 必须搭配亚二次乘法才有收益。CLPoly 已有 Karatsuba（P3-mul），满足此前提。

### 2.3 算法伪代码（递归 Half-GCD）

```
HGCD(A, B):
  n = deg(A), m = n/2
  if deg(B) < m: return identity matrix

  // 上半部分递归
  A_top = A >> m,  B_top = B >> m     // 取高 n-m 项
  R = HGCD(A_top, B_top)              // 递归
  (A', B') = R * (A, B)               // 矩阵乘法应用变换

  if deg(B') < m: return R

  // 一步 Euclid
  (q, r) = divmod(A', B')
  Q = [[0,1],[1,-q]] * R

  // 下半部分递归
  k = 2m - deg(B')
  A''_top = B' >> k,  B''_top = r >> k
  S = HGCD(A''_top, B''_top)

  return S * Q
```

---

## 3. FLINT 3.0.1 HGCD 实现分析

### 3.1 GCD 分派（`src/nmod_poly/gcd.c`）

FLINT 3.0.1 的 GCD 分派非常简洁：

```c
slong _nmod_poly_gcd(mp_ptr G, mp_srcptr A, slong lenA,
                     mp_srcptr B, slong lenB, nmod_t mod) {
    slong cutoff = NMOD_BITS(mod) <= 8
                 ? NMOD_POLY_SMALL_GCD_CUTOFF   // 200
                 : NMOD_POLY_GCD_CUTOFF;        // 340
    if (lenB < cutoff)
        return _nmod_poly_gcd_euclidean(G, A, lenA, B, lenB, mod);
    else
        return _nmod_poly_gcd_hgcd(G, A, lenA, B, lenB, mod);
}
```

**Cutoff 常量**（`src/nmod_poly.h`）：

| 常量 | 值 | 说明 |
|------|---|------|
| `NMOD_POLY_GCD_CUTOFF` | **340** | 适用于所有 >8-bit 素数（包括 64-bit） |
| `NMOD_POLY_HGCD_CUTOFF` | **100** | HGCD 递归内部 base case |
| `NMOD_POLY_SMALL_GCD_CUTOFF` | 200 | 仅 <=8-bit 素数 |

> **重要纠正**：旧版调研文档引用的 "1725 for 64-bit primes" 来自更早版本的 FLINT tuning 表。
> FLINT 3.0.1 已简化为统一 cutoff：对所有 >8-bit 素数使用 340。
> 这意味着 HGCD 在 deg >= 340 即有收益，而非之前认为的 deg >= 1725。

### 3.2 三层架构

FLINT 的 HGCD 实现分三层：

#### Layer 1: GCD 入口（`src/gr_poly/gcd_hgcd.c`）

```c
_gr_poly_gcd_hgcd(G, lenG, A, lenA, B, lenB, inner_cutoff=100, cutoff=340, ctx):
  divrem(Q, R, A, B)           // 先做一次 divrem 使 deg(A) >= deg(B)
  if R == 0: G = B; return

  hgcd(G, J, B, R, inner_cutoff=100)  // 第一次 HGCD

  while J != 0:
    divrem(Q, R, G, J)
    if R == 0: G = J; break
    if len(J) < cutoff:        // 340
      gcd_euclidean(G, J, R)   // fallback Euclid
      break
    hgcd(G, J, J, R, inner_cutoff=100)  // 继续 HGCD
```

**关键设计**：HGCD 不直接产出 GCD，而是快速缩减多项式度数。交替使用 HGCD + divrem，直到度数够小时 fallback Euclid。

#### Layer 2: HGCD 递归（`src/gr_poly/hgcd.c :: _gr_poly_hgcd_recursive`）

核心递归算法：

1. 设 `m = len(a) / 2`
2. 若 `len(b) < m + 1`，base case：返回单位矩阵
3. 取 a, b 的高位部分（右移 m 位得到 `a0`, `b0`）
4. 若 `len(a0) < cutoff(100)`：调用迭代 base case；否则递归
5. 得到矩阵 R 和缩减后的 `a3`, `b3`
6. 重构完整 `a2`, `b2`：将 R 应用到低位部分，再与高位合并
7. 若 `len(b2) < m + 1`：完成，返回 R
8. 否则：做一次 divrem → 取新的高位 → 再次递归得矩阵 S
9. 合并 `M = R * S`（2x2 多项式矩阵乘法）

**工作区大小**：每层 `6 * lena + 10 * (lena + 1) / 2` 元素，总计 `22 * lena + 16 * ceil(log2(lena) + 1)`。

#### Layer 3: 迭代 base case（`src/gr_poly/hgcd.c :: _gr_poly_hgcd_recursive_iter`）

当 `len < 100` 时使用：

```
_hgcd_recursive_iter(M, A, B, a, b):
  m = len(a) / 2
  M = identity
  A = a, B = b
  while len(B) >= m + 1:
    divrem(Q, T, A, B)
    swap(A, B); swap(B, T)       // 标准 Euclid 步
    // 更新矩阵：M' = [[0,1],[1,Q]] * M
    M[2], M[3] = M[3] + Q*M[2], M[2]    // (简化表达)
    M[0], M[1] = M[1] + Q*M[0], M[0]
```

本质是带矩阵追踪的 Euclid 迭代。

### 3.3 2x2 矩阵乘法（`__mat_mul`）

两种策略自适应：

```c
__mat_mul(C, A, B):
  min = min(lenA[0..3], lenB[0..3])
  if min < 20:
    __mat_mul_classical(C, A, B)   // 8 muls + 4 adds
  else:
    __mat_mul_strassen(C, A, B)    // 7 muls + 15 adds/subs
```

- **Classical**: `C[i][j] = sum_k A[i][k] * B[k][j]`，8 次多项式乘法 + 4 次加法
- **Strassen**: 7 次多项式乘法 + 15 次加减法，当矩阵元素度数 >= 20 时更优

### 3.4 关键子操作依赖

| 子操作 | FLINT 实现 | CLPoly 现状 |
|--------|-----------|-------------|
| 多项式乘法 `_gr_poly_mul` | schoolbook → Karatsuba → KS4 → NTT | Karatsuba (threshold=16) |
| 多项式 divrem `_gr_poly_divrem` | schoolbook | lazy 3-word 稠密 divrem |
| 多项式加法 `_gr_poly_add` | 逐系数加 | **需新增** |
| 多项式减法 `_gr_poly_sub` | 逐系数减 | **需新增** |
| 2x2 矩阵乘法 | classical / Strassen | **需新增** |
| 工作区管理 | 指针切片大块内存 | **需设计** |

---

## 4. NTL 参考

### 4.1 GCD 分派

**文件**：`src/ZZ_pX.cpp`, `src/lzz_pX.cpp`

NTL 通过 PrimeCnt 分级（素数大小分类）：
- GCD crossover: cnt=0: ~1400, cnt=1: ~800, cnt=2: ~400, cnt=3: ~600, cnt=4: ~1200

### 4.2 HGCD 实现

- 基于 Thull & Yap (1990) 的 iterative 变体
- 矩阵累乘在内层循环，批量应用
- `PlainHalfGCD` 处理 base case（小度数回退 Euclid）

### 4.3 三级乘法

| 度数范围 | 算法 | 复杂度 |
|---------|------|--------|
| < 16 | Schoolbook | O(n^2) |
| 16 ~ 150-500 | Karatsuba | O(n^1.585) |
| > 150-500 | Multi-prime FFT | O(n log n) |

---

## 5. 论文参考

### 5.1 核心论文

| 论文 | 贡献 | 备注 |
|------|------|------|
| **Schonhage (1971)** "Schnelle Berechnung von Kettenbruchentwicklungen" | 首次提出 O(M(n) log n) GCD | 原始论文，德文 |
| **Moenck (1973)** "Fast computation of GCDs", STOC | 将 Schonhage 的思想推广到多项式 | 标准引用 |
| **Thull & Yap (1990)** "A unified approach to HGCD algorithms" | 统一框架，NTL 的基础 | NTL 直接引用 |
| **Moller (2008)** "On Schonhage's algorithm and subquadratic integer GCD computation" | 现代改进，实用优化 | 工程实现参考 |
| **von zur Gathen & Gerhard** "Modern Computer Algebra" Ch.8-11 | 教科书级阐述 | FLINT 的主要理论来源 |

### 5.2 乘法前置论文

| 论文 | 贡献 |
|------|------|
| **Karatsuba & Ofman (1962)** | O(n^1.585) 乘法 |
| **Schonhage & Strassen (1971)** | FFT 乘法 O(n log n log log n) |
| **Harvey & van der Hoeven (2021)** | O(n log n) 整数乘法（理论最优） |

### 5.3 综述/教材

| 文献 | 说明 |
|------|------|
| **GCL** (Geddes, Czapor, Labahn) S7.4 | 模 GCD 框架（CLPoly 当前基础） |
| **MCA** (von zur Gathen & Gerhard) S8-11 | HGCD 完整理论 + 复杂度分析 |
| **Brent & Zimmermann** "Modern Computer Arithmetic" | 整数 HGCD，可类推到多项式 |

---

## 6. 演进路线与当前状态

### 6.1 路线总览

```
P3a: dense_upoly_zp + schoolbook mul + dense Euclid GCD   [已完成]
  |  消除稀疏表示开销，实测 2.15-3.1x
P3-mul: Karatsuba (threshold=16) + lazy 3-word divrem     [已完成]
  |  mul deg100: 3.5x → 1.07x vs FLINT
  |  divrem deg1000/500: 6.1x → 0.97x vs FLINT
P3-HGCD: Half-GCD (cutoff ~340)                           [当前]
  |  依赖 Karatsuba，预计 deg1000+ 显著加速
P3-KS: KS4 乘法 → GMP mpn_mul (crossover ~43)              [未来]
  |  进一步加速 HGCD + powmod
```

### 6.2 各阶段状态

| 阶段 | 任务 | 实测收益 | 状态 |
|------|------|---------|------|
| **P3a** | `dense_upoly_zp` 稠密类型 + schoolbook | 2.15-3.1x（消除稀疏开销） | **已完成** |
| **P3-mul** | Karatsuba + lazy 3-word divrem | mul 3.3x, divrem 6.3x, gcd 2.8x | **已完成** |
| **P3-HGCD** | Half-GCD | 预计 deg1000+ 2-5x | **调研中** |
| **P3-KS** | KS4 乘法 → GMP mpn_mul | -- | 未开始 |

### 6.3 HGCD 收益预估

| 度数 | 当前 ratio vs FLINT | HGCD 后预估 ratio | 说明 |
|------|--------------------|--------------------|------|
| deg200 | 1.25x | ~1.25x | 低于 cutoff，不进 HGCD |
| deg500 | 1.44x | ~1.2x | 刚过 cutoff，小幅改善 |
| deg1000 | 2.67x | ~1.5x | HGCD 显著缩减 Euclid 步数 |
| deg3000 | 4.28x | ~1.5-2x | HGCD 主力区间 |
| deg5000 | 5.28x | ~1.5-2x | 同上 |
| deg10000 | 6.87x | ~2x | 剩余差距来自 NTT（FLINT mul 用 KS4/NTT） |

注：HGCD 后的剩余差距主要来自乘法层（CLPoly Karatsuba O(n^1.585) vs FLINT KS4→GMP O(n^~1.4)），需 P3-KS 进一步消除。

---

## 7. CLPoly HGCD 实现要点

### 7.1 需要新增的操作

在 `dense_upoly_zp` 上：

| 操作 | 说明 | 用于 |
|------|------|------|
| `add(C, A, B)` | 逐系数模加 | 矩阵合并 |
| `sub(C, A, B)` | 逐系数模减 | 矩阵合并 |
| `rem(R, A, B)` | 仅取余式（省略商） | HGCD 间隔步 |

### 7.2 HGCD 特有组件

| 组件 | 说明 |
|------|------|
| 2x2 多项式矩阵类型 | 4 个 `dense_upoly_zp`（或 raw `uint64_t*`） |
| 矩阵乘法 | classical (min < 20) / Strassen (min >= 20) |
| 工作区管理 | 预分配连续内存块，指针切片 |
| `hgcd_recursive` | 核心递归算法 |
| `hgcd_iter` | 迭代 base case (len < 100) |

### 7.3 设计决策

1. **Cutoff 参数**：初始沿用 FLINT 的 inner=100, outer=340，后续可 benchmark 微调
2. **内存策略**：预分配 `std::vector<uint64_t>` 大块工作区 + 指针切片（C++ 惯例，避免裸 malloc）
3. **矩阵乘法**：初始仅实现 classical（8 mul + 4 add），Strassen 视 benchmark 需要再加
4. **接口**：`dense_upoly_zp::gcd` 内部自动分派：`deg < 340` → Euclid, `deg >= 340` → HGCD

### 7.4 复杂度估算

- 新增代码量：~400-500 行（hgcd 递归 + 迭代 base case + 矩阵操作 + add/sub）
- 集成改动：`dense_upoly_zp::gcd` 内部加分派，~10 行
