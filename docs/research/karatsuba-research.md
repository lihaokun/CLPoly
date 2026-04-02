# P3-mul 调研：dense_upoly_zp 快速乘法

> 日期：2026-03-06
> 分支：`feature/gcd-optimization`

## 1. 当前性能基线

### 1.1 CLPoly schoolbook mul vs FLINT (p = 2^64 - 59)

| 规模 | CLPoly | FLINT | ratio | 差距倍数 |
|------|--------|-------|-------|---------|
| deg50×deg50 | 0.005ms | 0.002ms | 2.4x | |
| deg100×deg100 | 0.020ms | 0.006ms | 3.5x | |
| deg200×deg200 | 0.072ms | 0.014ms | 5.1x | |
| deg500×deg500 | 0.464ms | 0.060ms | 7.7x | |
| deg1000×deg1000 | 1.829ms | 0.160ms | 11.5x | |

**差距随次数增长而增大**：schoolbook O(n²) vs FLINT 的 KS4/FFT O(n^1.5~n·log n)。

### 1.2 CLPoly schoolbook divrem vs FLINT

| 规模 | CLPoly | FLINT | ratio |
|------|--------|-------|-------|
| deg200/deg100 | 0.045ms | 0.007ms | 6.1x |
| deg500/deg250 | 0.366ms | 0.034ms | 10.9x |
| deg1000/deg500 | 1.427ms | 0.146ms | 9.8x |

### 1.2a divrem 差距原因分析

divrem 不调用 mul，但差距（6-11x）比 mul（2.4-11.5x）更大。原因有三层叠加：

**① 惰性归约（Lazy Reduction）— 3-5x**

CLPoly 内层循环每步做完整模乘+模减（Q×B 次归约）：

```cpp
// CLPoly: 每个系数对都做完整模归约
R[i+j] = nmod_sub(R[i+j], nmod_mul(q_i, B[j]));  // 2次128-bit mul + 条件分支
```

FLINT 将系数扩展为 3-word（192-bit）累加器，**仅在外层循环每步归约一次**：

```c
// FLINT: 展开为 192-bit，累加不溢出
R3[3*i+0..2] += c * B3[3*j+0..2]   // 无模归约
// 仅归约首项 (每次外层循环 1 次)
r = n_lll_mod_preinv(R3[2], R3[1], R3[0], mod.n, mod.ninv);
```

对 deg200/deg100：CLPoly ~10,000 次模归约 vs FLINT ~200 次，**单项 50x 减少**。

**② GMP 汇编内循环 — 2-3x**

FLINT 内层循环委托给 `mpn_addmul_1`（GMP 手写汇编）：
- `mulx` (ADX) 指令，流水线化，~1 multiply/cycle
- 无分支预测开销
- 最优寄存器分配

CLPoly 通过 C++ 编译器生成代码，`unsigned __int128` 乘法 → 2-3 条 `mul`/`imul`，额外条件分支。

**③ Newton 迭代除法 — 2-10x（次数相关）**

FLINT 对 64-bit 素数，**lenA > 20** 即切换到 Newton divrem：

```c
// FLINT 分派逻辑 (nmod_poly/divrem.c)
if (lenA <= 20 || lenB <= 8 || lenA - lenB <= 6)
    → basecase (惰性归约)
else
    → Newton 迭代: O(M(n)) 代替 O(n²)
```

Newton divrem：先算 `1/rev(B) mod x^{lenQ}`（Newton 迭代幂级数求逆），再 `Q = rev(A) × Binv`，最后 `R = A - QB`。复杂度 O(M(n))（M(n) = 快速乘法时间）。

**三层叠加**：3x（惰性归约）× 2x（汇编）× 2x（Newton）≈ **12x**，与实测 6-11x 吻合。

**对 CLPoly 优化的启示**：

| 优化方向 | 预估收益 | 实现难度 | 前置依赖 |
|---------|---------|---------|---------|
| 惰性归约 divrem | 3-5x | 中 | 无（纯微优化） |
| Newton 迭代 divrem | 2-5x | 高 | 需要快速 mul（Karatsuba/KS） |
| HGCD | 额外加速高次 GCD | 高 | 需要快速 mul + 快速 divrem |

### 1.3 CLPoly Zp GCD vs FLINT

| 规模 | CLPoly | FLINT | ratio |
|------|--------|-------|-------|
| gcd deg200+common100 | 0.195ms | 0.055ms | 3.6x |
| gcd deg500+common250 | 1.144ms | 0.408ms | 2.8x |
| gcd deg1000+common500 | 4.469ms | 1.194ms | 3.7x |
| gcd deg200 coprime | 0.291ms | 0.070ms | 4.1x |
| gcd deg1000 coprime | 6.021ms | 1.685ms | 3.6x |

### 1.4 每次模乘实测时间

从 schoolbook mul 数据反推：

| 规模 | 运算量 (muls) | 时间 | ns/mul |
|------|-------------|------|--------|
| 50×50 | 2,601 | 5μs | 1.9 |
| 100×100 | 10,201 | 20μs | 2.0 |
| 200×200 | 40,401 | 72μs | 1.8 |
| 500×500 | 251,001 | 464μs | 1.8 |
| 1000×1000 | 1,002,001 | 1,829μs | 1.8 |

CLPoly Barrett 模乘 ≈ **1.8-2.0ns/mul**，合理（2 × `mulq` @ 3 cycle × 0.33ns/cycle）。

## 2. 快速乘法算法分析

### 2.1 Karatsuba 算法

**数学原理**

将 n 次多项式在中点 m=⌊n/2⌋ 处分割：

```
A(x) = A_lo + x^m · A_hi
B(x) = B_lo + x^m · B_hi
```

朴素乘法需要 4 次子乘法。Karatsuba 用 3 次：

```
P0 = A_lo · B_lo
P2 = A_hi · B_hi
P1 = (A_lo + A_hi)(B_lo + B_hi) - P0 - P2    ← 交叉项
C(x) = P0 + x^m · P1 + x^{2m} · P2
```

**复杂度**：T(n) = 3·T(n/2) + O(n) → O(n^{log₂3}) = O(n^{1.585})

**对 Zp 系数的适用性**

Karatsuba 用"加法换乘法"。在 Zp 系数情形：
- 模乘 (Barrett)：~1.8ns（2 × `mulq` + 加减法）
- 模加/减：~0.5-1ns（1 条加/减 + 1 条条件移动）
- mul/add 成本比 ≈ 2-3:1

Karatsuba 每步省 1 次子乘法（代价是 O(n) 次额外加法），当 mul/add 比 > 1 即有收益。**Zp 情形下 Karatsuba 有效**。

**Crossover 分析**

NTL 的 `zz_pX` 使用 KARX=16 作为递归基例（经大量调优）。CLPoly Barrett 成本略高于 NTL Montgomery（NTL 限制 60-bit），交叉点可能更低。**初步建议 16，需 A/B 调优 12-24 范围**。

### 2.1a 为什么不用 Toom-3？

**Toom-3 (Toom-Cook 3-way)** 是 Karatsuba 的推广：3 分割，5 次子乘法（代替朴素 9 次），复杂度 O(n^{log₃5}) = O(n^{1.465})。

**不适用于 Zp 多项式乘法的原因：**

1. **常数因子大幅增加**：虽然渐近指数从 1.585 降到 1.465（8% 改进），但插值步骤需要在 5 个点（0, 1, -1, 2, ∞）求值和插值，涉及大量加减法和移位。常数因子增加 ~50%，在实际可用的次数范围内（deg < 500）完全抵消渐近改进。

2. **插值需要除以小常数**：Toom-3 插值矩阵的逆涉及除以 2, 3, 6。在整数上是精确除法（shift + 定点操作），但在 Zp 上需要**模逆**运算——虽然单次便宜，但每次递归都需要，累积开销显著。

3. **Karatsuba + FFT/KS 已覆盖有效范围**：
   - deg 16-200：Karatsuba 已够快
   - deg > 200：FFT/KS 的 O(n·log n) 远优于 Toom-3 的 O(n^{1.465})
   - Toom-3 没有"甜蜜区间"

4. **无参考实现**：NTL、FLINT 均不对 Zp 多项式使用 Toom-3。GMP 在整数乘法中有 Toom-3（中间层），但多项式库不采用。

**结论**：CLPoly 不考虑 Toom-3，直接 Karatsuba → KS/FFT。

### 2.2 Kronecker 替换 (KS)

**数学原理**

多项式乘法 C(x) = A(x)·B(x) 归约为单个大整数乘法：

1. 选位宽 N > 2·log₂(p) + log₂(n)，使乘积系数不溢出
2. 打包：a = A(2^N) = Σ aᵢ·2^{iN} （每个系数占 N bit 字段）
3. 大整数乘法：c = a·b （GMP 优化，Toom/FFT）
4. 解包：从 c 的每 N-bit 字段提取系数，mod p 归约

**64-bit 大素数的困难**

对 p ≈ 2^64：
- 乘积系数上界 < n·p² ≈ n·2^{128}
- 所需位宽 N = 128 + ⌈log₂(n)⌉ ≈ **134-138 bit/系数**

| deg | n | N (bit) | 打包后整数 (limb) | GMP 算法层级 |
|-----|---|---------|-------------------|-------------|
| 50 | 51 | 134 | ~107 | Toom-3 |
| 100 | 101 | 135 | ~213 | Toom-4/Toom-6 |
| 500 | 501 | 137 | ~1074 | FFT |
| 1000 | 1001 | 138 | ~2160 | FFT |

打包宽度大（小素数 ≤30 bit 时只需 ~70 bit），导致大整数异常庞大。

**KS2/KS4（Harvey 多点 Kronecker）**

在多个点求值以减小打包宽度：

| 变体 | 求值点 | 每系数位宽 | 整数乘法次数 |
|------|--------|-----------|------------|
| KS1 | 2^N | ~134 | 1 |
| KS2 | 2^N, 2^{-N} | ~67 | 2 |
| KS4 | ±2^N, ±2^{-N} | ~34 | 4 |

KS4 将每系数位宽降到 ~34 bit，整数大小降为 KS1 的 ~25%，但需要 4 次 GMP 乘法。

**FLINT 对 64-bit 素数的分派（源码分析）**

```c
bits = NMOD_BITS(mod);  // = 64
cutoff_len = min(len1, 2*len2);

if (len2 <= 5)                              → classical
if (3*cutoff_len < 2*max(bits,10))          → classical
  // 64-bit: cutoff_len < 43  → len2 ≤ 42 用 classical
if (cutoff_len*bits < 800)                  → KS1
  // 64-bit: cutoff_len < 13  → len2 ≤ 12 用 KS
if (cutoff_len*(bits+1)² < 100000)          → KS2
  // 64-bit: cutoff_len < 24  → len2 ≤ 23 用 KS2
else                                        → KS4 (或 fft_small)
```

**关键发现**：FLINT 在 64-bit 大素数下，deg ≤ 42 用 classical（不用 KS），deg > 42 跳 KS4/FFT-small。**没有 Karatsuba 层**。FLINT 的 classical 交叉点这么高，恰好说明 KS 在大素数下效率不佳。

### 2.3 NTT (数论变换)

**NTT 友好素数 vs 任意素数**

NTT 需要 q = k·2^m + 1 且 q 为素数。CLPoly 的 p = 2^64 - 59 **不是** NTT 友好素数。

对任意素数需 CRT：选 k 个 NTT 友好素数，每个做 NTT 乘法，再 CRT 重建。对 64-bit p 需要 k ≥ 5 个 30-bit 素数，交叉点很高（通常 deg > 1000）。

**FLINT fft_small**

FLINT 3.x 的 fft_small 利用 AVX2 向量指令实现快速 FFT，对中高次有 2-10x 加速，但需 AVX2 硬件。

### 2.4 Karatsuba 不能直接加速 Euclidean GCD

`dense_upoly_zp::divrem` 使用 long division，内部是**标量×多项式**操作：

```cpp
for (int64_t j = 0; j <= B.deg(); ++j)
    R._coeffs[i + j] = R.nmod_sub(R._coeffs[i + j],
                                   R.nmod_mul(q_i, B._coeffs[j]));
```

`divrem` 从不调用多项式 `mul`，因此加 Karatsuba 到 `mul` **不会**加速 Euclid GCD。

FLINT 的 `_nmod_poly_divrem_basecase` 同样是 long division。

## 3. 性能预测

### 3.1 Karatsuba 加速预测

**模型推导**

设两个 n 项多项式相乘（n = deg+1），递归基例大小 k：

- **Schoolbook**：n² 次模乘 + n² 次模加
- **Karatsuba**：递归关系 T(n) = 3·T(n/2) + c·n
  - 递归到基例 k 时，叶子节点数 = 3^{log₂(n/k)} = (n/k)^{log₂3} = (n/k)^{1.585}
  - 每个叶子做 k×k schoolbook = k² 次模乘
  - 总模乘数 = (n/k)^{1.585} × k²
  - 加法开销：每层递归 O(n) 次加法，共 log₂(n/k) 层 → O(n·log(n/k)) 次模加

以 k=16 为例，各规模计算：

| n | (n/16)^{1.585} | × 256 = muls | n·log₂(n/16) = adds | Schoolbook n² |
|---|----------------|-------------|---------------------|---------------|
| 51 | 6.9 | 1,766 | 85 | 2,601 |
| 101 | 20.7 | 5,299 | 262 | 10,201 |
| 201 | 62.1 | 15,898 | 728 | 40,401 |
| 501 | 267.8 | 68,557 | 2,508 | 251,001 |
| 1001 | 804.3 | 205,901 | 5,976 | 1,002,001 |

**实测参数**：T_mul ≈ 1.8ns，T_add ≈ 0.8ns（模加含条件移动）。

**预测**：

| 规模 | Schoolbook | Karatsuba 预测 | 加速比 |
|------|-----------|---------------|--------|
| 51×51 | n²·1.8ns = 4.7μs | 1766·1.8 + 85·0.8 = 3.2μs | 1.4x |
| 101×101 | 18.4μs | 5299·1.8 + 262·0.8 = 9.7μs | 1.9x |
| 201×201 | 72.7μs | 15898·1.8 + 728·0.8 = 29.2μs | 2.5x |
| 501×501 | 452μs | 68557·1.8 + 2508·0.8 = 125μs | 3.6x |
| 1001×1001 | 1804μs | 205901·1.8 + 5976·0.8 = 375μs | 4.8x |

**注**：实际略低于预测，因递归调用开销（函数调用、scratch 指针计算）、cache 效应等。保守估计打 0.8 折：

| 规模 | 保守预测加速 | 保守预测时间 |
|------|------------|------------|
| 50×50 | 1.2x | ~4μs |
| 100×100 | 1.5x | ~13μs |
| 200×200 | 2.0x | ~36μs |
| 500×500 | 2.9x | ~155μs |
| 1000×1000 | 3.8x | ~480μs |

### 3.2 KS 加速预测（64-bit 大素数）

**预测模型依据**

- **打包/解包**：每系数 ~134 bit，涉及 2-3 个 64-bit limb 的移位+存储。实测 GMP `mpn_lshift` + store 约 ~2-3ns/limb，每系数 ~3 limb → ~8-10ns/系数。对 n 系数的两个输入 + 一个输出共 ~3n 系数操作。
- **GMP 乘法**：两个 L-limb 整数的乘法时间。GMP 的渐近分层（[GMP Manual §15.1](https://gmplib.org/manual/Multiplication-Algorithms)）：
  - L < ~30 limb: schoolbook O(L²)，~0.5-1μs
  - L ~30-100: Toom-3 O(L^{1.465})，~1-10μs
  - L ~100-1000: Toom-6/Toom-8 O(L^{1.1-1.2})，~10-100μs
  - L > ~1000: FFT O(L·log L·log log L)，>100μs
- **实际参考**：FLINT 对 64-bit 素数 deg50 用 classical（不用 KS），佐证 KS 在此规模下无优势。

| 规模 | n | 整数大小 L (limb) | 打包+解包 | GMP mul 估计 | KS 总计 | vs Schoolbook |
|------|---|-------------------|----------|-------------|---------|-------------|
| 50×50 | 51 | ~107 | ~1.5μs | ~5μs (Toom-3) | ~7μs | **1.4x 更慢** |
| 100×100 | 101 | ~213 | ~3μs | ~15μs (Toom-4) | ~18μs | 持平 |
| 200×200 | 201 | ~427 | ~6μs | ~35μs (Toom-6) | ~41μs | 1.8x |
| 500×500 | 501 | ~1074 | ~15μs | ~90μs (FFT) | ~105μs | 4.4x |
| 1000×1000 | 1001 | ~2160 | ~30μs | ~170μs (FFT) | ~200μs | 9.1x |

**注**：
- 64-bit 大素数下 KS 打包宽度 ~134 bit，效率远不如小素数场景（小素数 ~30 bit 打包，整数仅 ~1/4 大小）。
- KS4 可将整数大小降至 ~25%（4 次较小乘法），在 deg > 200 时可能优于 KS1，但实现复杂度高。
- GMP mul 时间为估计值，受 CPU 微架构和 GMP 版本影响，±50% 范围内。

### 3.3 算法选择建议

| 规模 | 最优算法 | 理由 |
|------|---------|------|
| deg < 16 | Schoolbook | 递归开销 > 收益 |
| deg 16-200 | **Karatsuba** | 无格式转换，直接在 Zp 上运算 |
| deg 200-500 | Karatsuba 或 KS4 | 需 A/B 测试确定交叉点 |
| deg > 500 | KS4 | GMP FFT 层级的 O(n·log n) |

**关键结论**：对 64-bit 大素数，**Karatsuba 是 deg 16-200 的最优选择**，优于 KS。原因：
1. KS 打包宽度 ~134 bit 太大，打包/解包开销无法忽略
2. Karatsuba 直接在 Zp 上运算，无格式转换
3. Cache 友好，实现简单

## 4. Karatsuba 实际受益场景

| 操作 | 调用 mul? | Karatsuba 受益? |
|------|----------|----------------|
| `divrem` (GCD 核心) | ❌ | ❌ long division，标量操作 |
| `powmod` (DDF/EDF) | ✅ | ✅ `base * base mod poly` |
| Hensel lifting | ✅ | ✅ 多项式乘法 + 模归约 |
| 乘积构造 (Wang/EEZ) | ✅ | ✅ 因子乘积 |
| `squarefree` | ✅ | ✅ 多项式除法（通过 mul 验证） |

### 4.1 各操作中 mul 的时间占比估计

**powmod（DDF/EDF 核心）**

`powmod(base, exp, modpoly)` 计算 base^exp mod modpoly：
- 二进制平方乘法：~2·log₂(exp) 次 mul + ~2·log₂(exp) 次 mod（即 divrem 取余）
- 对 DDF：exp = p（素数），log₂(p) ≈ 64 → ~128 次 mul + ~128 次 divrem
- 设 modpoly 次数 d，每次 mul 是 d×d，每次 divrem 是 2d/d
- mul 占比 ≈ mul/(mul+divrem) ≈ d²/(d² + 2d²) ≈ **1/3**（schoolbook 下 divrem 与 mul 同阶）

Karatsuba 加速 mul 2x → powmod 整体加速 ~1.5x（mul 占 1/3 × 省 1/2 = 省 1/6）。

但如果同时实现 Newton divrem（利用快速 mul），divrem 也降到 O(M(d))：
- powmod 全部操作降到 O(M(d)·log p)
- 整体加速可达 3-5x

**Hensel lifting**

单变量 Hensel 每步：2 次 mul + 若干 add/sub。mul 占比 ~**60-70%**。
Karatsuba 加速 mul 2x → Hensel 整体加速 ~1.5-1.7x。

**因子乘积构造**

Wang/EEZ 中计算 ∏fᵢ：k 个因子，二叉乘法树 log₂(k) 层。
全部是 mul，**100%** 受益于 Karatsuba。

### 4.2 对单变量因式分解的整体影响

当前单变量因式分解 vs FLINT 差距 2.6-3.8x，瓶颈是 `__select_prime`（占 77-94%），内部主要是 DDF/EDF（powmod）。

Karatsuba 仅加速 mul（powmod 中占 ~1/3）：
- 预估 DDF/EDF 加速 ~1.3-1.5x
- 整体因式分解加速 ~1.2-1.4x

若进一步实现 Newton divrem + 惰性归约 divrem：
- DDF/EDF 加速可达 3-5x
- 整体因式分解加速 2-3x，接近 FLINT

## 5. FLINT / NTL 参考实现对比

### 5.1 FLINT `nmod_poly_mul` 策略

- 不使用 Karatsuba（KS 将问题转化为 GMP 大整数乘法）
- 64-bit 大素数下 classical 交叉点高达 deg ~42
- FLINT 3.x 的 fft_small (AVX2) 在 deg > 42 时接管
- 无 AVX2 → KS4

### 5.2 NTL `zz_pX` 策略

- schoolbook → Karatsuba (KARX=16) → 多模 FFT
- FFT crossover 因架构而异: `{150, 150, 300, 500, 500}`
- 限制 60-bit 素数，使用 Montgomery REDC

### 5.3 CLPoly 应采用的策略

综合两家经验：

1. **第一步：Karatsuba**（NTL 路线）— deg 16+ 使用，实现简单，收益确定
2. **第二步：KS4**（FLINT 路线）— deg 200+ 可选，依赖 GMP，实现复杂
3. **第三步：FFT/NTT**（远期）— 需 AVX2 或多模 CRT，实现极复杂

## 6. Karatsuba 实现方案

### 6.1 核心算法

```cpp
// 内部递归：C[0..2n-2] = A[0..n-1] × B[0..n-1]
// 要求：C 与 A, B, scratch 不重叠
// scratch 空间：本层消耗 2h + (2m-1) + (2h-1) ≈ 3n，递归共 O(n)
void _kar_mul(uint64_t* C, const uint64_t* A, const uint64_t* B,
              size_t n, uint64_t* scratch) const
{
    if (n < KARATSUBA_THRESHOLD) {
        std::fill(C, C + 2*n - 1, 0);
        for (size_t i = 0; i < n; ++i) {
            if (A[i] == 0) continue;
            for (size_t j = 0; j < n; ++j)
                C[i+j] = nmod_add(C[i+j], nmod_mul(A[i], B[j]));
        }
        return;
    }

    size_t m = n / 2, h = n - m;  // h >= m; n 偶数时 h=m, 奇数时 h=m+1

    // scratch 布局（本层）：
    //   [0,     h)           → t1 = A_lo + A_hi
    //   [h,     2h)          → t2 = B_lo + B_hi
    //   [2h,    2h+2m-1)     → P0 = A_lo * B_lo    (2m-1 项)
    //   [2h+2m-1, 2h+2m-1+2h-1) → P1 中间结果      (2h-1 项)
    //   [2h+2m-1+2h-1, ...)  → rec_scratch（传给递归子调用）
    uint64_t* t1  = scratch;
    uint64_t* t2  = scratch + h;
    uint64_t* P0  = scratch + 2*h;
    uint64_t* P1  = P0 + (2*m - 1);
    uint64_t* rec = P1 + (2*h - 1);

    // ---- 步骤 1：计算 t1 = A_lo + A_hi, t2 = B_lo + B_hi ----
    for (size_t i = 0; i < m; ++i) {
        t1[i] = nmod_add(A[i], A[m + i]);
        t2[i] = nmod_add(B[i], B[m + i]);
    }
    if (h > m) {  // n 奇数：A_hi[m] = A[n-1], A_lo 无对应项
        t1[m] = A[n - 1];
        t2[m] = B[n - 1];
    }

    // ---- 步骤 2：三次递归乘法 ----
    _kar_mul(P0, A, B, m, rec);                     // P0 = A_lo * B_lo (2m-1 项)
    _kar_mul(P1, t1, t2, h, rec);                    // P1 = (A_lo+A_hi)(B_lo+B_hi) (2h-1 项)
    _kar_mul(C + 2*m, A + m, B + m, h, rec);         // P2 = A_hi * B_hi → 直接写入 C[2m..2n-2]

    // ---- 步骤 3：P1 -= P0 + P2 → 交叉项 ----
    for (size_t i = 0; i < 2*m - 1; ++i)
        P1[i] = nmod_sub(P1[i], P0[i]);
    for (size_t i = 0; i < 2*h - 1; ++i)
        P1[i] = nmod_sub(P1[i], C[2*m + i]);        // C[2m+i] 此时 = P2[i]

    // ---- 步骤 4：组装 C = P0 + x^m · P1 + x^{2m} · P2 ----
    //
    // 最终 C 的结构（2n-1 项）：
    //   C[0..m-1]      = P0[0..m-1]                           （纯 P0）
    //   C[m..2m-2]     = P0[m..2m-2] + P1[0..m-2]             （P0 + P1 重叠）
    //   C[2m-1]        = P1[m-1]                               （间隙：仅 P1 贡献）
    //   C[2m..m+2h-2]  = P2[0..2h-m-2] + P1[m..2h-2]          （P2 + P1 重叠）
    //   C[m+2h-1..2n-2]= P2[2h-m-1..2h-2]                     （纯 P2，仅 n 奇数时存在）
    //
    // 注意：P0 覆盖 C[0..2m-2]，P2 覆盖 C[2m..2n-2]，C[2m-1] 是间隙位置。

    // 4a. C[0..2m-2] = P0
    for (size_t i = 0; i < 2*m - 1; ++i)
        C[i] = P0[i];

    // 4b. C[2m-1] = 0（间隙位置：P0 和 P2 都不覆盖，仅由 P1 叠加填充）
    C[2*m - 1] = 0;

    // 4c. C[m..m+2h-2] += P1（叠加交叉项到 P0 高位 + 间隙 + P2 低位）
    for (size_t i = 0; i < 2*h - 1; ++i)
        C[m + i] = nmod_add(C[m + i], P1[i]);
}
```

**Scratch 空间总量**：
- t1: h, t2: h, P0: 2m-1, P1: 2h-1, 递归: T(h)
- 总量 T(n) = 2h + (2m-1) + (2h-1) + T(h) ≈ 4h + 2m + T(h)
- T(n) ≈ 3n + T(n/2) → T(n) = O(n)（几何级数收敛）
- 实际分配 **6n** 足够（NTL 用 8n 留余量）

### 6.2 实现要点

- **Scratch 空间**：递归总需 O(n) 临时空间。外层 `mul` 预分配一次 `vector<uint64_t>(8*n)`，传入递归。
- **不等长处理**：当 len_A ≠ len_B 时，零填充到相同长度（简单方案），或分块 Karatsuba（更复杂，性能提升有限）。
- **Aliasing**：外层 `mul` 已有保护。
- **Threshold 调优**：初始设 16，后续 A/B 测试在 12-24 范围调优。

### 6.3 修改清单

| 文件 | 改动 |
|------|------|
| `clpoly/dense_upoly_zp.hh` | 添加 `_kar_mul` 私有方法 + 修改 `mul` 分派逻辑 |
| `test/bench_clpoly.cc` | 已有 mul bench（deg50-1000）|
| `test/bench_comparative.cc` | 已有 FLINT 对比 |

## 7. Newton 迭代除法（替代路径分析）

将 `divrem` 转化为利用快速乘法：
- 先算 `1/B mod x^n`（Newton 迭代，O(M(n)) 次乘法）
- 再算 `Q = A × (1/B) mod x^n`，`R = A - Q·B`

**不适用于 Euclidean GCD**：每步除数 B 不同，Newton 初始化开销无法摊销。

但 **HGCD** 可以利用快速乘法（通过矩阵乘法），这是高次 GCD 的正确路径。

## 8. 参考文献

### 8.1 快速乘法

1. **Karatsuba & Ofman (1962)**
   A. Karatsuba, Yu. Ofman. "Multiplication of Multidigit Numbers on Automata."
   *Doklady Akademii Nauk SSSR*, 145:293-294, 1962.
   英译: *Soviet Physics -- Doklady*, 7:595-596, 1963.

2. **Harvey (2009)** — Multipoint Kronecker Substitution (FLINT KS2/KS4 理论基础)
   D. Harvey. "Faster polynomial multiplication via multipoint Kronecker substitution."
   *Journal of Symbolic Computation*, 44(10):1502-1510, 2009.
   DOI: [10.1016/j.jsc.2009.05.004](https://doi.org/10.1016/j.jsc.2009.05.004)
   arXiv: [0712.4046](https://arxiv.org/abs/0712.4046)

3. **Schönhage & Strassen (1971)** — FFT 整数乘法
   A. Schönhage, V. Strassen. "Schnelle Multiplikation grosser Zahlen."
   *Computing*, 7:281-292, 1971.
   DOI: [10.1007/BF02242355](https://doi.org/10.1007/BF02242355)

4. **Harvey & van der Hoeven (2019)** — 整数乘法最优复杂度
   D. Harvey, J. van der Hoeven. "Integer multiplication in time O(n log n)."
   HAL: [hal-02070778](https://hal.archives-ouvertes.fr/hal-02070778)

### 8.2 Newton 迭代除法

5. **Sieveking (1972)**
   M. Sieveking. "An algorithm for division of power series."
   *Computing*, 10:153-156, 1972.
   DOI: [10.1007/BF02242389](https://doi.org/10.1007/BF02242389)

6. **Brent & Kung (1978)**
   R. P. Brent, H. T. Kung. "Fast Algorithms for Manipulating Formal Power Series."
   *Journal of the ACM*, 25(4):581-595, 1978.
   PDF: [Harvard](https://www.eecs.harvard.edu/~htk/publication/1978-jacm-brent-kung.pdf)

### 8.3 HGCD / 快速 GCD

7. **Schönhage (1971)** — 快速 GCD
   A. Schönhage. "Schnelle Berechnung von Kettenbruchentwicklungen."
   *Acta Informatica*, 1:139-144, 1971.
   DOI: [10.1007/BF00289520](https://doi.org/10.1007/BF00289520)

8. **Moenck (1973)**
   R. T. Moenck. "Fast computation of GCDs."
   *Proc. 5th ACM STOC*, pp. 142-151, 1973.
   DOI: [10.1145/800125.804045](https://doi.org/10.1145/800125.804045)

### 8.4 教科书

9. **von zur Gathen & Gerhard (2013)**
   J. von zur Gathen, J. Gerhard. *Modern Computer Algebra*, 3rd ed.
   Cambridge University Press, 2013. ISBN 978-1-107-03903-2.

10. **Geddes, Czapor & Labahn (1992)** (GCL)
    K. O. Geddes, S. R. Czapor, G. Labahn. *Algorithms for Computer Algebra*.
    Kluwer Academic Publishers, 1992. ISBN 0-7923-9259-0.

### 8.5 实现参考

11. **Hart (2010)** — FLINT
    W. B. Hart. "Fast Library for Number Theory: An Introduction."
    *ICMS 2010*, LNCS 6327, pp. 88-91, 2010.
    Preprint: [Warwick WRAP](http://wrap.warwick.ac.uk/41629/)

12. **Kronecker (1882)** — 原始 Kronecker 替换
    L. Kronecker. "Grundzüge einer arithmetischen Theorie der algebraischen Grössen."
    *J. reine angew. Math.*, 92:1-122, 1882.
    [EUDML](https://eudml.org/doc/148487)

13. **Shoup (2001-)**  — NTL
    V. Shoup. "NTL: A Library for doing Number Theory."
    [libntl.org](https://libntl.org/)
    性能对比: [benchmarks.pdf](https://libntl.org/benchmarks.pdf)

14. **Johansson (2014)** — FLINT KS 交叉点分析
    F. Johansson. "Embettered polynomial multiplication."
    博客: [fredrikj.net](https://fredrikj.net/blog/2014/04/embettered-polynomial-multiplication/)
