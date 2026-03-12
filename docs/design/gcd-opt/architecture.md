# GCD 优化架构文档

> 状态：已确认（2026-02-28）
> 调研依据：`docs/research/gcd-optimization-research.md`

---

## 1. 核心问题

CLPoly GCD 在 Benchmark 中与 FLINT 存在 4-200x 性能差距。调研报告（§5）识别出的核心原因：

| 瓶颈 | 当前状况 | FLINT 做法 | 差距根源 |
|------|----------|-----------|---------|
| 素数太小 | 从 prime(0)=2 开始递增 | 从 2^63 开始 | 每素数贡献 ~1 bit vs ~63 bits |
| 所需素数过多 | ~200 个素数 | ~4 个素数 | 循环次数差 50x |
| content 重复计算 | 多变量入口 `cont(F)` 计算两次（hh:297-299） | 复用已有变量 | 2x content 开销 |
| 素数表上限 | boost 表仅 9999 个（p<104729） | 无上限（n_nextprime） | 大系数输入会耗尽 |

## 2. 优化路线图

按优先级排列，每个模块可独立实现和测试：

```
M0.1: ZZ 接口扩展 ─┐
M0.2: 素数生成组件 ─┤
M0.3: boost 迁移   ─┤─→ P0: 大素数 + content 修复 ── GCD 核心优化（预期消除 90% 差距）
M0.4: Zp 修复 ─────┤                                  ↓
                    ↓                          P1: GCDHEU ── 启发式快速路径（小系数提速 2-5x）
M0.5: MTSHL 素数 ──┘                                  ↓
   (依赖 M0.2+M0.4)                           P2: CRT 优化 ── 增量 content + 稳定检测
                                                       ↓
                                               P3: 稠密 Zp 快速路径 ── 长期（deg>200 才有意义）
```

## 3. 模块设计

### M0: 基础设施

#### M0.1: ZZ 接口扩展

在 `clpoly/number/ZZ.hh` 中新增以下 public 接口：

```cpp
// 只读访问器
int64_t   get_val() const;     // 返回 _val（小模式下是值，大模式下未定义）
mpz_srcptr get_mpz() const;    // 返回 _mpz（小模式下返回 nullptr）
bool      is_small() const;    // 等价于 _mpz == nullptr

// 构造函数
ZZ(mpz_srcptr z);              // 从 mpz 构造，自动 demote_if_small
```

**设计决策**：

- `get_mpz()` 返回 `nullptr` 可作为模式判定条件，调用方无需额外 `is_small()` 检查
- `ZZ(mpz_srcptr)` 构造函数调用 `_fits_si` 判定后决定走 `_val` 还是 `_mpz` 分支，与 `ZZ(const char*)` 构造函数模式一致
- `get_val()` 仅在小模式下有意义；大模式下行为未定义（由调用方负责先检查 `is_small()`）

#### M0.2: 素数生成组件

在 `clpoly/number/ZZ.hh` 中（或独立头文件 `clpoly/number/prime.hh`）新增：

```
函数签名                              实现方式
────────────────────────────────────────────────────────
uint64_t next_prime_64(uint64_t n)    mpz_nextprime + overflow 检查
uint64_t prev_prime_64(uint64_t n)    mpz_prevprime（GMP 6.3+）
ZZ       next_prime(const ZZ& n)      mpz_nextprime + demote_if_small
```

**`next_prime_64` 规约**：
- 前置条件：`n >= 2`
- 后置条件：返回最小素数 `p > n`
- 异常：若 `p` 超出 `uint64_t` 范围（`mpz_sizeinbase(z,2) > 64`），抛出 `std::overflow_error`
- 实现：直接调用 GMP C API，**不复用** ZZ 版本（避免性能损失）

**`prev_prime_64` 规约**：
- 前置条件：`n >= 3`
- 后置条件：返回最大素数 `p < n`
- 实现：调用 `mpz_prevprime`（GMP ≥ 6.3.0，已确认安装版本满足）
- 用途：MTSHL 素数选择（从 `2^64 - 59` 向前找不整除 lc 的素数）

**`next_prime(ZZ)` 规约**：
- 前置条件：`n >= 2`
- 后置条件：返回最小素数 `p > n`
- 实现：通过 `get_mpz()` / `get_val()` 获取 GMP 值，调用 `mpz_nextprime`，构造 ZZ 返回

**三个函数完全独立**，各自直接调用 GMP，不互相依赖。

#### M0.3: boost::math::prime 完全迁移

移除 `boost/math/special_functions/prime.hpp` 依赖。迁移方式：
**不再使用查表**，全部改为 `next_prime_64` 顺序调用。

**迁移清单**：

| # | 文件 | 行号 | 当前代码 | 迁移方式 |
|---|------|------|---------|---------|
| 1 | `polynomial_gcd.hh` | 10 | `#include <boost/...prime.hpp>` | 删除 |
| 2 | `polynomial_factorize_univar.hh` | 17 | `#include <boost/...prime.hpp>` | 删除 |
| 3 | `polynomial_gcd.cc` | 43-44 | `p_index=0; prime=boost::math::prime(p_index)` | `prime = INITIAL_PRIME; // 2^63` |
| 4 | `polynomial_gcd.cc` | 59-60 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 5 | `polynomial_gcd.cc` | 75-76 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 6 | `polynomial_gcd.cc` | 194-195 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 7 | `polynomial_gcd.hh` | 312-318 | 用 `p_index` 查表找 > degree 的素数 | `prime = INITIAL_PRIME; // 2^63` |
| 8 | `polynomial_gcd.hh` | 333-334 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 9 | `polynomial_gcd.hh` | 349-350 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 10 | `polynomial_gcd.hh` | 463-464 | `++p_index; prime=boost::math::prime(p_index)` | `prime = next_prime_64(prime)` |
| 11 | `polynomial_factorize_univar.hh` | 1397 | `boost::math::prime((unsigned)idx)` | `next_prime_64(p)` |
| 12 | `polynomial_factorize_wang.hh` | 2024 | `boost::math::prime(idx)` | `next_prime_64(p)` |
**附加清理**：
- 删除所有 `p_index` 变量和 `>= 9999` 的边界检查（`next_prime_64` 无表限制）
- 因式分解中的 `PRIME_TABLE_SIZE` 常量可移除

**风险**：无。`mpz_nextprime` 是 GMP 核心函数，生产级可靠。Miller-Rabin ~1-10μs/调用，但 GCD/因式分解通常只需 3-10 个素数，总开销 < 100μs，远小于单次 DDF 或 CRT 的 ms 级开销。

#### M0.4: Zp 类大素数安全修复

当前 `Zp` 类（`number.hh`）在 `p ≥ 2^63` 时存在多处缺陷。本轮 GCD 优化（P0）使用 `p ≈ 2^62`，在安全范围内；但 **M0.5 MTSHL 素数修复要求 `p = 2^64 - 59`，必须先完成 M0.4**。

**缺陷清单**（`clpoly/number.hh`）：

| # | 行号 | 缺陷 | 触发条件 | 修复方案 |
|---|------|------|---------|---------|
| 1 | 145-146 | `(int64_t)p` 溢出 | `p > INT64_MAX` | FLINT `nmod_set_si` 方案：unsigned 取绝对值 → unsigned 取模 → 符号补偿 |
| 2 | 164-165 | `(int64_t)this->_p` 溢出 | 同上 | 同上 |
| 3 | 193-194 | `op1._i + op2._i` 溢出 | `p > 2^63` | FLINT `nmod_add` 方案：用 `p - a` 与 `b` 比较代替直接加法 |
| 4 | 200-201 | `this->_i + op2._i` 溢出 | 同上 | 同上（`+=` 版本） |
| 5 | 109 | `assert(p < (1ULL << 63))` | — | 修复 1-4 后放宽为 `assert(p >= 2)` |
| 6 | 136 | `assert(p < (1ULL << 63))` | — | 同上 |

**减法无问题**：`p - op2._i + op1._i` 在 `op1._i < op2._i` 分支中结果恒 `< p < 2^64`，无溢出。

**修复方案参考**（FLINT `nmod.h`）：

```c
// FLINT nmod_set_si — 有符号数 → 模值（全程 unsigned）
ulong nmod_set_si(slong x, nmod_t mod) {
    ulong res = (x >= 0) ? (ulong)x : -(ulong)x;  // 取绝对值
    NMOD_RED(res, res, mod);                         // unsigned Barrett 取模
    return (res == 0 || x > 0) ? res : mod.n - res;  // 负数补偿
}

// FLINT nmod_add — 无溢出加法
ulong nmod_add(ulong a, ulong b, nmod_t mod) {
    const ulong neg = mod.n - a;   // p - a，不溢出
    if (neg > b)
        return a + b;              // a + b < p，安全
    else
        return b - neg;            // = a + b - p
}
```

**对应 CLPoly 修复**：

```cpp
// 构造函数 + operator=（替换 int64_t 取模）
uint64_t abs_i = (i >= 0) ? (uint64_t)i : -(uint64_t)i;
uint64_t r = abs_i % p;
_i = (r == 0 || i > 0) ? r : p - r;

// operator+（替换直接加法）
uint64_t neg = op1._p - op1._i;
if (neg > op2._i)
    op1._i = op1._i + op2._i;
else
    op1._i = op2._i - neg;
```

#### M0.5: MTSHL 素数选择修复（正确性）

**依赖**：M0.2（`prev_prime_64`）+ M0.4（Zp 支持 `p > 2^63`）

**问题**：`polynomial_factorize_wang.hh:2272` 硬编码 `mtshl_p = 2^63 - 25`，无 bad-prime 检查。若 `mtshl_p | lc(scaled_factors[i])`，`__mtshl_lift` 内 `Zp(coeff, p)` 将 leading term 静默丢弃（line 1053），导致因子次数下降、提升结果错误。当 `lc(f)` 为常数且为 `mtshl_p` 的倍数时，所有求值点均失败，算法无法恢复。

**修复流程**：
```
// 在 eval point 循环外，选 p 一次
p = 2^64 - 59                                  // 最大 64-bit 素数（已验证）
while p | content(lc(f, x1)):                   // 检查 p 是否整除 lc 所有系数
    p = prev_prime_64(p)

for skip = ...:
    __select_eval_point(g, x1, skip, p)
        ├── lc(f)(α) ≠ 0 over Z                // 已有
        └── lc(f)(α) % p ≠ 0                   // 新增，复用已有 delta，极少触发
    factorize(f(x₁,α))
    __wang_leading_coeff → scaled_factors
    __mtshl_lift(..., p)
```

**数学依据**（MTSHL 论文）：
- Theorem 19（JSC 2020）：失败概率 ∝ 1/p，p 越大越好，仅要求 p > 2d
- CASC 2018 Algorithm 5 中 "63-bit 机器素数" 是 Maple signed int64 工程限制，非数学要求

---

### P0: 大素数 + content 修复

**改动范围**：`polynomial_gcd.cc`（单变量）、`polynomial_gcd.hh`（多变量模板）

#### P0.1: 素数起始位置

```
当前：prime = boost::math::prime(0)        // prime = 2
改后：prime = 2^64 - 59                    // 最大 64-bit 素数，向后推进
```

**理论依据**：Brown §4.4 eq.43 推荐 `p ∈ [½(β+1), β]`。M0.4 修复 Zp 后 `β = 2^64 - 1`，起始点选 `2^64 - 59`（最大 64-bit 素数）。与 MTSHL 统一从最大素数向后选取。

**效果估算**：

| 指标 | 当前 | 改后 |
|------|------|------|
| 每素数贡献 bits | ~1-17（依次增长） | ~64 |
| Mignotte bound 400 bits 需要素数数 | ~200 | ~7 |
| CRT 循环次数 | ~200 | ~7 |

#### P0.2: content 重复计算修复

多变量 GCD 入口（`polynomial_gcd.hh:297-299`）存在 content 重复计算 bug：

```cpp
polynomial_<ZZ,...> F_cont = cont(F);   // 第一次计算
polynomial_<ZZ,...> G_cont = cont(G);   // 第一次计算
polynomial_<ZZ,...> cont_gcd = polynomial_GCD(cont(F), cont(G));  // 重新计算！
```

第 299 行应改为 `polynomial_GCD(F_cont, G_cont)`，复用已有结果。

**注**：CRT 循环中的 `cont()` 调用（`polynomial_gcd.cc:169`、`polynomial_gcd.hh:437`）已经位于 `if (tmp_Pout_ == Pout_)` 稳定检测内部，仅在终止时执行一次，无需优化。

#### P0.3: 移除 p_index 和 9999 限制

`p_index` 变量及 `>= 9999` 分支全部删除。改用 `prime = next_prime_64(prime)` 推进。循环终止条件改为 CRT 乘积超过 Mignotte bound 时自然终止（而非素数表耗尽）。

#### P0 接口规约

```
接口不变：
  单变量：polynomial_GCD(upolynomial_<ZZ>, upolynomial_<ZZ>) → upolynomial_<ZZ>
  多变量：polynomial_GCD(polynomial_<ZZ,lex_<var_order>>, ...) → polynomial_<ZZ,lex_<var_order>>

内部行为变化：
  - 素数从 2^64 - 59 开始，向后推进（M0.4 修复后）
  - CRT 循环次数从 ~200 降至 ~7
  - 多变量入口 content 重复计算已修复
```

---

### P1: GCDHEU（启发式 GCD）

**调研依据**：Parisse 2002 + FLINT `fmpz_poly/gcd_heuristic.c`

在小系数单变量 GCD 前插入启发式快速路径：将多项式 bit-pack 为大整数，用 GMP `mpn_gcd_full` 求整数 GCD，然后 unpack 验证。

**适用条件**：单变量 Z[x]，系数 bits + 6 < word-pack 阈值（FLINT 用 `pack_bits < 32` 时按 bit 打包，`≥ 32` 时按 limb 打包）。

**模块位置**：`polynomial_gcd.cc` 中 `polynomial_GCD` 入口处，P0 循环之前。

**接口**：
```
int gcdheu(upolynomial_<ZZ>& result, const upolynomial_<ZZ>& F, const upolynomial_<ZZ>& G);
// 返回 1 = 成功（result 为 GCD），0 = 失败（回退到模算法）
```

**优先级说明**：P1 独立于 P0，可并行开发。但对大系数/高次多项式无效（FLINT 在大于阈值时直接跳过），因此 P0 是更通用的优化。

---

### P2: CRT 优化 — ❌ 废弃

**结论（2026-03-05）**：经 A/B 测试验证，P2 的三项优化（P2a 稳定检测合并、P2b content early exit + lc_gcd 快速路径、P2c trial division 预检查）均无可测量的性能差异。根因：CRT 循环瓶颈在 Zp 多项式 Euclid GCD（O(n²)），bookkeeping 开销（content、比较、trial division）占比 <1%，优化无意义。详见 `p2-crt-optimization.md`。

---

### P3: 稠密 Zp 快速路径（长期）

**调研结论**：HGCD 在 deg > 1725（FLINT 阈值）时才有意义。当前 benchmark 均 deg ≤ 200，故降级为长期项。

若未来需要：
- 实现 `_nmod_poly_gcd_euclidean` 风格的原生 Zp[x] GCD（替代当前通过 `Zp` 类的间接实现）
- 高次时接入 HGCD（half-GCD）亚二次算法

---

## 4. 关键设计决策

### 决策 1：GCD 和 MTSHL 统一使用 `2^64 - 59` 起始

FLINT 用 `UWORD(1) << 63`，但 FLINT 全程 unsigned。CLPoly 的 `Zp` 类在 M0.4 修复前存在 `(int64_t)p` 溢出和加法溢出缺陷，限制 `p < 2^63`。

- **M0.4 修复后**：`Zp` 支持全范围 `uint64_t` 素数。GCD 和 MTSHL 统一从 `2^64 - 59`（最大 64-bit 素数）开始，用 `prev_prime_64` 向后推进。每素数贡献 ~64 bits，优于 FLINT 的 ~63 bits

### 决策 2：`next_prime_64` 放在 ZZ.hh 还是独立头文件

推荐放在 `clpoly/number/ZZ.hh`，理由：
- 函数依赖 GMP（ZZ 已 include gmp.h）
- ZZ 版本 `next_prime` 直接操作 ZZ 内部状态
- 避免新增头文件

### 决策 3：先做 P0 再做 P1

P0（大素数）消除的是系统性瓶颈（~200 次循环 → ~7 次），预期提速 20-30x。
P1（GCDHEU）是特定场景优化（小系数），提速 2-5x 但不适用大系数。
实现顺序：M0 → P0 → P1。

### 决策 4：保持 GCD 接口不变

所有优化均为内部实现改动，`polynomial_GCD` 的签名和语义不变。调用方无需修改。

---

## 5. 预期收益

| 模块 | 预期提速 | 依赖 | 风险 |
|------|---------|------|------|
| M0 | 无直接提速（修复 Zp 后可用更大素数） | 无 | 低（接口扩展 + bug fix） |
| P0 | 20-30x（大系数）| M0 | 低（FLINT 已验证路线） |
| P1 | 2-5x（小系数）| M0 | 低（FLINT gcdheu 成熟） |
| P2 | ~~1.2-1.5x~~ 废弃（实测无效） | P0 | — |
| P3 | 高次多项式 | P0 | 中（实现复杂度高） |

目标：P0 完成后，CLPoly GCD 与 FLINT 的差距从 4-200x 缩小到 2-5x。
