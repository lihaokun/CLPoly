# Zp 类 64-bit 扩展调研报告

> 日期：2026-02-27
> 关联：M6a（`docs/design/multivar-opt/post-m5-optimizations.md` 差异 3）
> 动机：MTSHL-d 需要 63-bit 机器素数（CASC 2018 Algorithm 5），当前 Zp 仅支持 31-bit

---

## 1. CLPoly 现有 Zp 实现

**位置**：`clpoly/number.hh` L85-262

### 1.1 数据成员

```cpp
uint32_t _i;     // 值 ∈ [0, p)
uint32_t _p;     // 素数 p（0 = 未初始化哨兵）
uint64_t _ninv;  // Barrett 常数：UINT64_MAX / p
```

对象大小：24 bytes（含 padding）。

### 1.2 乘法（Barrett 归约）

```cpp
// 预计算（构造时）：
_ninv = UINT64_MAX / p;

// 归约：
uint32_t __barrett_reduce(uint64_t product) const {
    uint64_t q = (unsigned __int128)product * _ninv >> 64;
    uint64_t r = product - q * _p;
    return (uint32_t)(r >= _p ? r - _p : r);
}
```

工作原理：`q ≈ product × (2⁶⁴/p) / 2⁶⁴ = product/p`。误差至多 1，一次修正。
前提：`product < p² < 2⁶²`（p < 2³¹），乘积 `a*b` 可放入 uint64。
已使用 `__int128`：`(unsigned __int128)product * _ninv`。

### 1.3 逆元

```cpp
// number.hh L72
inline uint64_t inv_prime(uint64_t _i, uint32_t _p) {
    uint64_t a=_p, b=_i, c;
    uint64_t s1=0, s2=1, s3;
    while (c=(a%b)) {
        s3 = (s1 + _p - (s2 * (a/b)) % _p) % _p;
        a=b; b=c; s1=s2; s2=s3;
    }
    return s2;
}
```

扩展 Euclidean 算法。中间乘法 `s2 * (a/b)` 在 p < 2³¹ 时不溢出 uint64。

### 1.4 约束

`assert(p >= 2 && p < (1u << 31))` — 硬编码 31-bit 上限。

### 1.5 使用点统计

| 模块 | 文件 | 素数来源 | 素数范围 |
|------|------|---------|---------|
| 单变量因式分解 | `factorize_univar.hh` | `boost::math::prime` 表 | < 10⁵ |
| 多变量 GCD | `polynomial_gcd.hh/.cc` | `boost::math::prime` 表 | < 10⁵ |
| Zp 因式分解 (DDF/EDF) | `factorize_zp.hh` | 来自 `__select_prime` | < 10⁵ |
| MTSHL | `factorize_wang.hh` | 固定 `mtshl_p` | 当前 ~2×10⁹，目标 ~9.2×10¹⁸ |

**关键发现**：只有 MTSHL 需要大素数。其余模块全部使用小素数（< 10⁵）。

---

## 2. 外部参考实现

### 2.1 NTL `zz_p`

**存储**：单个 `long` 值 + 全局模数（`zz_p::init(p)` 设置）。

**乘法**：预计算逆元归约（`PrepMulMod`），内部分三档：
- 0-23 bit：AVX 浮点
- 24-31 bit：64-bit 整数乘法
- 32+ bit：标准整数运算

**模数上限**：60-bit（保守设计，避免 `long` 溢出）。

**模数管理**：`zz_pContext` 保存/恢复全局模数，支持多模数切换。

**参考价值**：NTL 的 60-bit 限制是人为保守设计，非技术必要。

### 2.2 FLINT `nmod`

**存储**：

```c
typedef struct {
    mp_limb_t n;      // 模数
    mp_limb_t ninv;   // 预计算逆元
    flint_bitcnt_t norm;  // 归一化位移量
} nmod_t;
```

**核心思想**：归一化 Barrett 归约。

```
初始化：
  norm = clz(n)              // 前导零数
  pn   = n << norm           // 左移使最高位为 1
  ninv = preinvert(pn)       // floor(2^128 / pn) - 2^64

  由于 pn ∈ [2^63, 2^64)，2^128/pn ∈ [2^64, 2^65)，
  所以 ninv = 2^128/pn - 2^64 ∈ [0, 2^64)，恰好放入 uint64。

乘法 NMOD_MUL_PRENORM(res, a, b, mod)：
  a_shifted = a << norm            // a < p → a_shifted < pn ≤ 2^64
  hi:lo = a_shifted * b            // __int128 乘 (umul_ppmm)
  q1:q0 = hi * ninv               // __int128 乘
  q1:q0 += hi:lo                  // 128-bit 加
  r = lo - (q1+1) * pn            // 近似余数（可能偏差 ≤ 2）
  if r > q0: r += pn              // 修正 1（下溢）
  if r >= pn: r -= pn             // 修正 2（超范围）
  res = r >> norm                  // 反归一化
```

**模数上限**：全 64-bit，无人为限制。

**性能**：2 次 `__int128` 乘法 + 若干 64-bit 运算 ≈ 15 cycles。

**逆元预计算** `n_preinvert_limb_prenorm(n)`：
```c
// n 必须最高位为 1
udiv_qrnnd(ninv, dummy, ~n, ~(mp_limb_t)0, n);
// 即 ninv = (2^128 - n*2^64 - 1) / n
```

**参考价值**：FLINT 是最成熟的开源 64-bit 模运算实现，CLPoly 应直接采用此方案。

### 2.3 简单方案 `(unsigned __int128)a * b % p`

编译为 x86-64 `mulq` + `divq`。
- `mulq`：~3 cycles（64×64→128 乘法）
- `divq`：~35-90 cycles（128/64 除法，CPU 和数值相关）

总计约 40-50 cycles。无需预计算，代码最简。

---

## 3. 方案对比

### 3.1 性能

| 方案 | 小素数 (p < 2³¹) | 大素数 (p ~ 2⁶³) | 预计算开销 |
|------|-----------------|-----------------|-----------|
| 当前 CLPoly Barrett | ~10 cycles | 不支持 | `UINT64_MAX / p` |
| FLINT 归一化 Barrett | ~15 cycles | ~15 cycles | `clz + preinvert` |
| `__int128 %`（divq） | ~40 cycles | ~40 cycles | 无 |

### 3.2 综合评价

| 方案 | 优势 | 劣势 |
|------|------|------|
| A: 统一 FLINT 归一化 Barrett | 全素数范围最优性能；代码路径唯一 | 小素数比现状慢 ~1.5x；对象 +8 bytes；`preinvert` 实现需仔细 |
| B: 双路径（小=现有，大=divq） | 小素数零退化；代码简单 | 大素数慢 ~2.7x；分支逻辑 |
| C: 全部 divq | 代码最简 | 小素数慢 ~4x（DDF/EDF 性能退化明显） |

---

## 4. 推荐方案

**方案 A：统一 FLINT 归一化 Barrett。**

理由：

1. **代码路径唯一**：无分支，维护简单，不会出现"小素数走路径 X、大素数走路径 Y"的分裂
2. **大素数最优**：MTSHL 的 63-bit 素数也享受 ~15 cycles 乘法，与 FLINT 对齐
3. **小素数退化可接受**：10 → 15 cycles（1.5x），DDF/EDF 主要瓶颈是 `powmod` 的算法复杂度 O(d²·log p) 而非单次乘法
4. **与 FLINT 完全对齐**：FLINT 是 CLPoly 的性能基准，使用相同底层策略便于对比分析
5. **可渐进优化**：若后续 profiling 显示 15 cycles 是瓶颈，可为小素数加特化分支，不影响接口

潜在风险：
- 对象从 24 bytes → 32 bytes（增加 `_norm` 字段），多项式中每个 Zp 系数增大 8 bytes
- 需仔细实现 `preinvert`（参考 FLINT `n_preinvert_limb_prenorm`）
