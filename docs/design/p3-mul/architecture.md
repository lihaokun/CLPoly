# P3-mul 架构文档：dense_upoly_zp 快速算术

> 状态：待确认
> 调研依据：`docs/research/karatsuba-research.md`

---

## 1. 核心问题

P3a 引入 `dense_upoly_zp` 后，Zp GCD 已获 1.6-3x 加速（稠密表示 vs 稀疏）。但与 FLINT 仍有显著差距：

| 操作 | CLPoly | FLINT | ratio | 差距来源 |
|------|--------|-------|-------|---------|
| mul deg100×100 | 0.020ms | 0.006ms | 3.5x | schoolbook O(n²) vs KS4/FFT |
| mul deg500×500 | 0.464ms | 0.060ms | 7.7x | 同上，高次差距更大 |
| divrem deg200/100 | 0.045ms | 0.007ms | 6.1x | 逐元素模归约 + 无 Newton |
| divrem deg500/250 | 0.366ms | 0.034ms | 10.9x | 同上 |
| gcd deg200+common100 | 0.195ms | 0.055ms | 3.6x | divrem 效率 |
| gcd deg1000 coprime | 6.021ms | 1.685ms | 3.6x | 同上 |

调研（§1.2a）识别出 divrem 差距的三层原因：

1. **逐元素模归约**（3-5x）：CLPoly 内循环每步做完整 Barrett 模乘+模减；FLINT 用 3-word 累加器延迟归约，逐系数 `umul_ppmm` + `add_sssaaaaaa` 累加
2. **底层指令效率**（~1.5x）：`umul_ppmm` 编译为单条 `mulq`，`add_sssaaaaaa` 用 inline asm 保证 `adcq` 链（纯 C 无法利用进位标志），流水线友好
3. **算法层级**（2-10x）：FLINT 对 64-bit 素数 lenA > 20 即用 Newton divrem O(M(n))

## 2. 优化路线图

两个模块共享 M2.0 辅助函数（`word3`/`_umul128`/`_add_carry3`/`_lll_mod_preinv`），在此之上独立实现：

```
M1: Karatsuba mul                M2: 惰性归约 divrem
  修改 dense_upoly_zp::mul         修改 dense_upoly_zp::divrem
  deg < 16: schoolbook             3-word 累加器 + _umul128
  deg >= 16: Karatsuba             仅首项做模归约
  预期: mul 加速 1.5-4x            预期: divrem 加速 2-3x
       ↓                                ↓
  受益: powmod, Hensel,            受益: Euclid GCD（直接加速）
       因子乘积                          DDF/EDF powmod 中的 mod
       ↓                                ↓
  ┌────┴────────────────────────────────┘
  ↓
  后续: Newton divrem (依赖 M1)
  后续: HGCD (依赖 M1 + Newton divrem)
```

**M1 和 M2 共享底层辅助函数**：`_classical_mul`（M1 基例）和 `divrem`（M2）均使用 3-word 惰性累加（`word3`/`_umul128`/`_add_carry3`/`_lll_mod_preinv`）。辅助函数作为 M2.0 先行实现，M1 和 M2 在此之上独立开发。

## 3. 模块设计

### M1: Karatsuba 乘法

#### M1 功能规约

```
模块名称：Karatsuba mul

功能描述：为 dense_upoly_zp::mul 添加 Karatsuba 分治乘法，
         deg >= THRESHOLD 时自动切换，替代 schoolbook O(n²)。

前置条件（Requires）：
  - A._p == B._p（同一素数域）
  - A, B 非空（空输入由外层 mul 处理）

后置条件（Ensures）：
  - C(x) = A(x) · B(x) mod _p
  - C 无 leading zero（__strip 保证）
  - A, B 不被修改（const 引用）

不变式（Invariants）：
  - _classical_mul 基例使用 192-bit 惰性累加器，中间值可超过 uint64_t，
    每输出系数归约到 [0, _p)（通过 _lll_mod_preinv）
  - _kar_mul 递归层使用逐元素 nmod_add/nmod_sub，值始终在 [0, _p) 内
  - scratch 空间仅在 mul 调用期间使用，不持久化

副作用：无（纯函数，输出通过引用参数）
```

#### M1 内部结构

**新增私有成员**：

```cpp
static constexpr size_t KARATSUBA_THRESHOLD = 16;  // 初始值，需 A/B 调优

// Karatsuba 核心递归
// C[0..2n-2] = A[0..n-1] × B[0..n-1]，scratch 为预分配临时空间
void _kar_mul(uint64_t* C, const uint64_t* A, const uint64_t* B,
              size_t n, uint64_t* scratch) const;

// schoolbook 基例（惰性累加，点积形式，复用 M2.0 辅助函数）
void _classical_mul(uint64_t* C, const uint64_t* A, size_t len_a,
                    const uint64_t* B, size_t len_b) const;
```

**修改 `mul` 分派逻辑**：

```cpp
static void mul(dense_upoly_zp& C,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    // aliasing 保护（不变）
    // 空输入检查（不变）

    size_t n = std::max(A.size(), B.size());
    if (n < KARATSUBA_THRESHOLD) {
        // schoolbook（现有逻辑）
        C._classical_mul(...);
    } else {
        // Karatsuba
        // 1. 零填充 A, B 到相同长度 n
        // 2. 分配 scratch: vector<uint64_t>(6 * n)
        // 3. 调用 _kar_mul(C._coeffs.data(), ...)
        // 4. __strip()
    }
}
```

#### M1 Karatsuba 算法

输入：A[0..n-1], B[0..n-1]（升幂系数）
输出：C[0..2n-2]（乘积系数）

```
_kar_mul(C, A, B, n, scratch):
  if n < THRESHOLD:
    _classical_mul(C, A, n, B, n)
    return

  m = n / 2, h = n - m    // h = m 或 m+1

  // scratch 布局：t1[h] | t2[h] | P0[2m-1] | P1[2h-1] | rec[...]
  t1 = A_lo + A_hi         // h 项，逐元素 nmod_add
  t2 = B_lo + B_hi         // h 项

  P0 = A_lo × B_lo         // 递归，2m-1 项 → scratch
  P1 = t1 × t2             // 递归，2h-1 项 → scratch
  P2 = A_hi × B_hi         // 递归，2h-1 项 → C[2m..2n-2]

  P1 -= P0                  // 逐元素 nmod_sub
  P1 -= P2                  // 逐元素 nmod_sub → P1 = 交叉项

  C[0..2m-2] = P0
  C[2m-1] = 0               // 间隙位置
  C[m..m+2h-2] += P1        // 叠加交叉项
```

**复杂度**：O(n^{log₂3}) = O(n^{1.585}) 次模乘 + O(n·log n) 次模加

**Scratch 空间**：T(n) = 2h + (2m-1) + (2h-1) + T(h) ≈ 3n + T(n/2) = O(n)。外层分配 6n 足够。

#### M1 不等长处理

当 len_A ≠ len_B 时，短的一方零填充到 max(len_A, len_B)。

零填充开销 = O(n) memset，远小于乘法本身。若 len_A >> len_B（极端不对称），schoolbook 可能更优——但当前不特殊处理，后续按需添加分块乘法。

---

### M2: 惰性归约 divrem

#### M2 功能规约

```
模块名称：惰性归约 divrem

功能描述：优化 dense_upoly_zp::divrem 的内层循环，
         将逐元素模归约改为延迟归约（3-word 累加器），
         内循环使用 _umul128 + _add_carry3 逐系数累加。

前置条件（Requires）：
  - B 非空（除数不为零）
  - A._p == B._p
  - deg(A) >= deg(B)（否则由外层处理：Q=0, R=A）

后置条件（Ensures）：
  - A(x) = Q(x)·B(x) + R(x) mod _p
  - deg(R) < deg(B) 或 R = 0
  - Q, R 无 leading zero

不变式（Invariants）：
  - 3-word 累加器中间值可能超过 _p，但最终归约后 ∈ [0, _p)
  - 每次外层迭代仅归约首项（1 次 _lll_mod_preinv）

副作用：无
```

#### M2 核心思路

**当前实现**（逐元素归约）：

```
外层循环 i = q_len-1 downto 0:
    q_i = nmod_mul(R[i+deg(B)], inv_lc)       // 1 次模乘
    内层循环 j = 0 to deg(B):
        R[i+j] = nmod_sub(R[i+j], nmod_mul(q_i, B[j]))  // 每项: 1 模乘 + 1 模减
```

每次内层循环：1 次 Barrett 模乘（2 × 128-bit mul）+ 1 次模减 + 条件分支。
总模归约次数 = q_len × (deg(B)+1)。

**惰性归约实现**：

```
// 每个系数用独立的 3-word (lo, mid, hi) 累加器，系数间不传播进位
// R3[i] = {lo, mid, hi} 表示 192-bit 未归约值

// 1. 初始化 3-word 表示
for each i: R3[i] = {R[i], 0, 0}

外层循环 i = q_len-1 downto 0:
    // 仅归约首项（1 次 3-word → 1-word 模归约）
    r = _lll_mod_preinv(R3[i+deg(B)].hi, R3[...].mid, R3[...].lo)
    q_i = nmod_mul(r, inv_lc)

    // negmod：将减法转为加法，保证累加器始终非负
    c = (q_i == 0) ? 0 : p - q_i

    // 内层循环：逐系数 3-word 累加，系数间无进位传播
    for j = 0 to deg(B):
        _umul128(p1, p0, c, B[j])                   // 128-bit 乘积
        _add_carry3(R3[i+j], p1, p0)               // 3-word 加法

// 2. 最终归约：仅余式部分
for i = 0 to deg(B)-1:
    R[i] = _lll_mod_preinv(R3[i].hi, R3[i].mid, R3[i].lo)
```

**两个关键设计要点**：

1. **为什么不用 `mpn_addmul_1` 跨系数调用**：`mpn_addmul_1` 将数组视为连续大整数，进位会从系数 j 的高位字泄漏到系数 j+1 的低位字，导致最终逐系数归约时读到被污染的值。必须逐系数独立累加。

2. **为什么用 negmod + 加法而非直接减法**：`sub_dddmmmsss` 会产生负值（unsigned wrap around 到接近 2^192 的大数）。由于 2^192 mod p ≠ 0，`_lll_mod_preinv` 对 wrapped 值的归约结果与真实负值的 mod p 不一致。FLINT 的做法：先取 `c = p - q_i`（negmod），再用 `add_sssaaaaaa` 加上 `c * B[j]`，保证累加器始终非负，`_lll_mod_preinv` 归约正确。

**关键优化点**：

| 指标 | 当前 | 惰性归约 |
|------|------|---------|
| 模归约次数 | q_len × (deg(B)+1) | q_len + deg(B) |
| 内层乘法 | Barrett（2×128-bit mul + 归约逻辑）| `_umul128`（单条 `mulq`，无归约）|
| 内层分支 | 每元素 2 个条件分支（模加减） | 无（3-word 加减无条件分支） |
| 内存用量 | 1x（原地 uint64_t）| 3x（3-word 扩展）|

#### M2 3-word 归约函数

```cpp
// 将 3-word 值 (hi, mid, lo) 归约为 r ∈ [0, p)
// 等价于 (hi * 2^128 + mid * 2^64 + lo) mod p
static uint64_t _lll_mod_preinv(uint64_t hi, uint64_t mid, uint64_t lo,
                                uint64_t p, uint64_t pinv, uint32_t norm)
{
    // 两次 2-word Barrett 归约：
    // 第一步：(hi, mid) → r1 ∈ [0, p)
    // 第二步：(r1, lo) → r2 ∈ [0, p)
    // 参考 FLINT n_lll_mod_preinv
}
```

#### M2 编译器内建函数与 inline asm

内循环依赖以下辅助函数（M2.0 实现，M1/M2 共用）：

```cpp
// _umul128: 128-bit 乘法（__int128，编译器自动生成单条 mulq）
unsigned __int128 prod = (unsigned __int128)a * b;
uint64_t p0 = (uint64_t)prod, p1 = (uint64_t)(prod >> 64);

// _add_carry3: 3-word 加法（平台分派 inline asm）
//   x86_64:  addq + adcq + adcq（3 条指令）
//   AArch64: adds + adcs + adc （3 条指令）
//   其他:    __int128 fallback  （~16 条指令）
```

**为什么 `_add_carry3` 必须用 inline asm**：C 语言无法表达 CPU 进位标志，编译器无法将溢出检测（`s.lo < old_lo` 等比较链）优化为 `adcq`。实测纯 C 比较链生成 ~18 条指令，在 deg199 内循环中比 inline asm 慢 68%。

无需新增外部依赖。不直接使用 FLINT 宏，而是自行实现等价的平台分派 inline asm + fallback。

注：不使用 3-word 减法（FLINT `sub_dddmmmsss`），因为 unsigned 减法产生的 wrapped 负值在 `_lll_mod_preinv` 归约时不正确（2^192 mod p ≠ 0）。改用 negmod + `_add_carry3` 保证累加器始终非负。

#### M2 溢出安全性分析

每个系数的 3-word 累加器独立，需要保证不溢出 192 bit。

由于使用 negmod + `_add_carry3`，累加器始终非负且单调递增。

初始值：R[i] ∈ [0, p) < 2^64。
单次加法上界：`c = p - q_i < p < 2^64`，`B[j] < p < 2^64` → 乘积 `c × B[j]` < 2^{128}。
q_len 次累加：总值 < p + q_len × p^2 < q_len × 2^{128}。
192-bit 容量 = 2^{192}，溢出条件：q_len > 2^{64} ≈ 1.8 × 10^{19}。

**不可能溢出**。实际 q_len ≤ deg(A) ≤ 10^6 量级，远小于 2^{64}。

注：
- 累加器始终非负（negmod + 加法），`_lll_mod_preinv` 归约正确
- 逐系数独立累加（非 `mpn_addmul_1` 跨系数），无进位泄漏

---

## 4. 接口规约

### 4.1 公开接口（不变）

```
接口：调用方 → dense_upoly_zp

mul(C, A, B):
  输入：A, B ∈ dense_upoly_zp，同一素数域
  输出：C = A × B mod p
  协议：调用方保证 A._p == B._p；C 可与 A/B 别名（内部处理）

divrem(Q, R, A, B):
  输入：A, B ∈ dense_upoly_zp，B 非空，同一素数域
  输出：A = Q × B + R，deg(R) < deg(B)
  协议：调用方保证 B 非空；Q/R 可与 A/B 别名（内部处理）

gcd(G, A, B):
  输入：A, B ∈ dense_upoly_zp，同一素数域
  输出：G = gcd(A, B)（monic 归一化由调用方负责）
  说明：内部调用 divrem，自动受益于 M2 优化
```

**所有公开接口签名和语义不变**。M1/M2 均为内部实现优化，调用方零感知。

### 4.2 内部接口（新增）

```
接口：mul → _kar_mul

输入数据：
  - C: uint64_t* 输出数组，长度 2n-1（调用方预分配）
  - A, B: const uint64_t* 输入数组，长度 n（零填充后等长）
  - n: size_t 输入长度
  - scratch: uint64_t* 临时空间，长度 >= 6n（调用方预分配）

协议约定：
  - 调用方保证 C 与 A, B, scratch 不重叠
  - 调用方保证 scratch 长度足够
  - 被调用方保证 C[0..2n-2] = A × B mod p
  - 被调用方不修改 A, B
  - scratch 内容在调用后未定义
```

```
接口：divrem 内循环 → _umul128 / _add_carry3

_umul128(hi, lo, a, b):
  输入：a, b ∈ uint64_t
  输出：(hi, lo) = a × b，128-bit 无符号乘积
  协议：unsigned __int128 乘法，编译为单条 mulq（x86-64）

_add_carry3(s, b1, b0):
  输入：s ∈ word3（192-bit 累加器），(b1, b0) ∈ 2-word 128-bit 值
  输出：s += (b1 << 64) | b0
  协议：结果始终非负（因使用 negmod + 加法策略）；
        最终由 _lll_mod_preinv 归约到 [0, p)
  说明：x86_64/AArch64 使用 inline asm（3 条指令），其他平台 __int128 fallback；
        输出可与输入别名
```

---

## 5. 关键设计决策

### 决策 1：Karatsuba 而非 KS

**选择**：M1 使用 Karatsuba 递归乘法，不使用 Kronecker 替换。

**理由**（详见调研 §2.1 vs §2.2）：
- 64-bit 大素数下 KS 打包宽度 ~134 bit/系数，打包/解包开销大
- FLINT 在 64-bit 下 classical 交叉点高达 42，佐证 KS 效率不佳
- Karatsuba 直接在 Zp 上运算，无格式转换，cache 友好
- NTL 对 zz_pX 使用 Karatsuba（KARX=16），实践验证有效

### 决策 2：Threshold 初始值 16

**选择**：`KARATSUBA_THRESHOLD = 16`，后续 A/B 调优。

**理由**：
- NTL 经大量调优确定 KARX=16
- CLPoly Barrett 模乘成本与 NTL Montgomery 相当（~1.8ns vs ~1.5ns）
- 调优范围 12-24，用 bench_clpoly 的 mul 基准 A/B 测试

### 决策 3：Scratch 空间外层一次性分配

**选择**：在 `mul` 入口分配 `vector<uint64_t>(6 * n)`，传指针给递归。

**理由**：
- 避免递归中反复 malloc/free
- 6n 足够（T(n) ≈ 3n + T(n/2) → 收敛到 6n）
- `vector` 管理生命周期，无泄漏风险

### 决策 4：逐系数 negmod + `_add_carry3` 累加

**选择**：M2 内循环先取 `c = p - q_i`（negmod），再逐系数 `_umul128` + `_add_carry3` 累加。

**理由**：
- **不跨系数**：`mpn_addmul_1` 将数组视为连续大整数，进位会从系数 j 泄漏到系数 j+1，导致最终逐系数归约错误
- **不用减法**：`sub_dddmmmsss` 产生 unsigned wrapped 负值，由于 2^192 mod p ≠ 0，`_lll_mod_preinv` 对 wrapped 值归约不正确。negmod + 加法保证累加器始终非负
- FLINT `_nmod_poly_divrem_basecase` 同样使用 negmod + `add_sssaaaaaa`
- `_umul128`（`__int128`）编译为单条 `mulq`；`_add_carry3` 使用 inline asm 保证 `addq` + `adcq` 链，已充分利用硬件流水线
- 主要加速来自延迟模归约（从 q×B 次减到 q+B 次），而非汇编级差异

### 决策 5：M1 和 M2 独立实现

**选择**：两个模块在 M2.0 共享辅助函数之上独立实现，分别修改 `mul` 和 `divrem`。

**理由**：
- M2.0 辅助函数（`word3`/`_umul128`/`_add_carry3`/`_lll_mod_preinv`）先行实现，M1 和 M2 共用
- `divrem` 内循环是标量×多项式操作，不调用 `mul`，M2 不依赖 M1
- M1 加速 `mul`，受益路径是 powmod/Hensel/因子乘积
- M2 加速 `divrem`，受益路径是 Euclid GCD
- M2.0 就绪后，M1 和 M2 可并行开发，各自 A/B 测试

### 决策 6：不实现 Toom-3

**选择**：跳过 Toom-3，直接 Karatsuba。

**理由**（详见调研 §2.1a）：
- 渐近改进仅 8%（n^{1.585} → n^{1.465}），常数因子增加 ~50%
- 插值需除以 2/3/6，Zp 上需模逆
- Karatsuba + KS/FFT 已覆盖有效范围，Toom-3 无甜蜜区间
- NTL、FLINT 均不对 Zp 多项式使用 Toom-3

---

## 6. 改动文件清单

| 文件 | M1 改动 | M2 改动 |
|------|--------|--------|
| `clpoly/dense_upoly_zp.hh` | 新增 `_kar_mul`, `_classical_mul`；修改 `mul` 分派 | 重写 `divrem` 内循环（3-word 累加器） |
| （共享 M2.0）| `_classical_mul` 和 `divrem` 共用 | 新增 `word3`, `_umul128`, `_add_carry3`, `_lll_mod_preinv` |
| `test/bench_clpoly.cc` | 已有 mul bench | 已有 divrem bench |
| `test/bench_comparative.cc` | 已有 FLINT 对比 | 已有 FLINT 对比 |

**不改动的文件**：
- `polynomial_gcd.hh`：`__polynomial_GCD` 调用 `gcd`/`divrem`，自动受益
- `polynomial_factorize_zp.hh`：后续迁移到 dense 时才改
- `upolynomial.hh`：稀疏类型不变

---

## 7. 预期收益

### M1: Karatsuba mul（含惰性 schoolbook 基例）

注：`_classical_mul` 已改为 3-word 惰性累加（点积形式），与 FLINT `_NMOD_VEC_DOT3` 同等效率。
预测分两层：schoolbook 基例加速 ~2-3x（实测 2.15-3.1x）+ Karatsuba 分治加速。

| 规模 | 当前 | 预测 | 总加速 | FLINT | vs FLINT |
|------|------|------|--------|-------|---------|
| deg50×50 | 0.005ms | ~0.003ms | 1.7x | 0.002ms | 1.5x |
| deg100×100 | 0.020ms | ~0.008ms | 2.5x | 0.006ms | 1.3x |
| deg200×200 | 0.072ms | ~0.022ms | 3.3x | 0.014ms | 1.6x |
| deg500×500 | 0.464ms | ~0.090ms | 5.1x | 0.060ms | 1.5x |
| deg1000×1000 | 1.829ms | ~0.280ms | 6.5x | 0.160ms | 1.8x |

低次（<16）纯受益于惰性 schoolbook；高次叠加 Karatsuba 分治效果更显著。

**间接受益**：powmod 加速 ~1.5-2x → 单变量因式分解加速 ~1.3-1.5x

### M2: 惰性归约 divrem

| 规模 | 当前 | 预测 | 加速 |
|------|------|------|------|
| divrem deg200/100 | 0.045ms | ~0.015-0.023ms | 2-3x |
| divrem deg500/250 | 0.366ms | ~0.120-0.180ms | 2-3x |
| divrem deg1000/500 | 1.427ms | ~0.480-0.710ms | 2-3x |

加速来源：模归约次数从 q_len×(deg(B)+1) 降为 q_len+deg(B)，内循环去除条件分支。`_umul128` 单条 `mulq` vs Barrett 的 2×128-bit mul + 归约逻辑。

**间接受益**：Euclid GCD 加速 ~1.5-2x → 单变量/多变量 ZZ GCD 加速 ~1.3-1.5x

### 两者叠加

| 场景 | 当前 vs FLINT | M1+M2 后预测 | 目标 |
|------|-------------|-------------|------|
| Zp GCD deg200 | 3.6x | ~1.2-1.8x | 接近 FLINT |
| Zp GCD deg1000 | 3.6x | ~1.2-1.8x | 同上 |
| mul deg500 | 7.7x | ~1.5x | 剩余差距来自 KS4/FFT |
| divrem deg500 | 10.9x | ~3.5-5.5x | 剩余差距来自 Newton divrem |
| 单变量因式分解 | 2.6-3.8x | ~1.5-2.5x | 显著缩小 |

**剩余差距**主要来自：
1. mul: FLINT 用 KS4/FFT-small（亚二次），CLPoly Karatsuba 仍是 O(n^{1.585})
2. divrem: FLINT deg>20 用 Newton O(M(n))，CLPoly 惰性归约仍是 O(n²)
3. 标量模乘: FLINT 小素数用 Montgomery REDC，略快于 Barrett

这些差距留待后续阶段（Newton divrem → HGCD）解决。
