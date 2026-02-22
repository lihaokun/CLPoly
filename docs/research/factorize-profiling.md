# 单变量因式分解管线 Profiling 调研报告

> 日期：2026-02-22
> 目标：定位 CLPoly 与 FLINT/NTL 8-16x 性能差距的根因，评估 P1a/P1b/Monagan 2022 的实际影响范围，确定下一步优化方向。
> 参考源码：FLINT `src/fmpz_poly_factor/factor_zassenhaus.c`、`factor_van_hoeij.c`；NTL `src/ZZXFactoring.cpp`

---

## 1. 现有基础设施与实测数据

### 1.1 单变量因式分解管线（CLPoly 当前实现）

```
factorize(f)
  └─ __factor_squarefree_primitive_ZZ(f)
       ├─ __select_prime(f)          §8.1  — 选择最优模素数 + Zp 因式分解
       ├─ __hensel_lift(f, …)        §6    — 二次 Hensel 提升到 Mignotte 精度
       └─ __lll_factorize(f, …)     §8.3  — van Hoeij LLL 重组（P1a）
            └─ __vanhoeij_recombine  §8.3  — LLL 主控循环
```

### 1.2 实测 Profiling 数据（release build，中位数 7 次）

工具：`test/profile_factorize.cc`，分别独立计时各步骤。

| 用例 | total | sel_prime | % | hensel | % | recombine | % | r | a_mig |
|------|-------|-----------|---|--------|---|-----------|---|---|-------|
| Wilkinson W(10) | 2.36ms | 2.22ms | **94%** | 0.09ms | 4% | 0.02ms | 1% | 10 | 10 |
| Wilkinson W(15) | 4.50ms | 4.00ms | **89%** | 0.26ms | 6% | 0.04ms | 1% | 15 | 14 |
| Wilkinson W(20) | 7.54ms | 6.33ms | **84%** | 0.77ms | 10% | 0.07ms | 1% | 20 | 19 |
| x^15-1 | 2.61ms | 2.31ms | **88%** | 0.05ms | 2% | 0.08ms | 3% | 5 | 15 |
| x^24-1 | 9.15ms | 7.65ms | **84%** | 0.19ms | 2% | 1.37ms | 15% | 13 | 7 |
| ~deg15 3 fac | 3.19ms | 2.51ms | **79%** | 0.12ms | 4% | 0.07ms | 2% | 5 | 5 |
| ~deg25 5 fac | 10.02ms | 7.76ms | **77%** | 0.29ms | 3% | 1.18ms | 12% | 8 | 8 |

**结论：`__select_prime` 占总时间 77-94%，是绝对瓶颈。**

### 1.3 与 FLINT/NTL 的差距（来源：benchmarks/2026-02-21.txt）

| 用例 | CLPoly | FLINT | NTL | CLPoly/FLINT | CLPoly/NTL | NTL/FLINT |
|------|--------|-------|-----|--------------|-----------|-----------|
| W(10) | 2.567ms | 0.166ms | 0.210ms | **15.5x** | **12.2x** | 1.3x |
| W(15) | 5.263ms | 0.633ms | 0.892ms | **8.3x** | **5.9x** | 1.4x |
| ~deg15 3 fac | 2.439ms | 0.153ms | 0.150ms | **15.9x** | **16.3x** | 1.0x |

---

## 2. P1a / P1b / Monagan 2022 的实际影响范围

### 2.1 P1a：van Hoeij LLL 重组（已实现）

- **优化范围**：`recombine` 步骤
- **实测占比**：1-15%
- **理论上限**：即使 recombine 降为 0，总时间最多减少 15%
- **实际观测**：与 Feb-21 Zassenhaus 基线（W(15)=5.023ms）相比，P1a 实现后 W(15)=4.50ms，**在误差范围内无显著差别**
- **评价**：P1a 消除了 Zassenhaus 2^r 爆炸（uni-70 从理论上限 2^70 降到实际 77ms），但在典型用例（r ≤ 20）上 Zassenhaus 本身就快，LLL 重组不是瓶颈

### 2.2 P1b：线性 Hensel 提升（基础设施已实现，未完整集成）

- **优化范围**：`hensel` 步骤
- **实测占比**：2-10%
- **理论上限**：即使 hensel 降为 0，总时间最多减少 10%
- **评价**：对整体性能影响极为有限。线性 Hensel 的真正价值在于超大精度（a_mig >> 100）下减少中间数膨胀，但当前测试用例 a_mig ≤ 19，不在线性 Hensel 的适用场景内

### 2.3 Monagan 2022 矩阵优化

- **优化范围**：线性 Hensel 提升中的 Bézout 更新步骤（O(r·n²·D) → O(n·D²)）
- **实测占比**：属于 `hensel` 步骤（2-10%）的一部分
- **评价**：影响范围是 hensel 步骤内的子步骤，实际贡献极小。在 `__select_prime` 不解决的情况下，实现 Monagan 2022 无法带来可观测的整体提升

---

## 3. `__select_prime` 瓶颈根因分析

### 3.1 CLPoly 当前实现（polynomial_factorize_univar.hh §8.1）

```
__select_prime(f):
  max_tries = 30（若 best_count > deg/2 则扩展到 50）
  for each prime p in 素数表（顺序遍历）:
    (1) lc(f) mod p ≠ 0 检查
    (2) f mod p 的次数保持检查
    (3) GCD(fp, fp') == 1 无平方检查     → 通过则 ++tried
    (4) DDF(fp) — 按次数分组
    (5) EDF(fp) — 等次随机分裂，得到 r_p 个因子
    记录 best = argmin(r_p)
  return best（最小-r 策略）
```

**关键参数**：~~固定尝试 30-50 个通过无平方检查的素数~~ → **[已改为 max_tries=3，2026-02-23]**，每次做完整 DDF+EDF。

### 3.2 最小-r 策略的历史背景

最小-r 策略来自 Zassenhaus 算法时代：重组需枚举 2^r 个子集，r 越小越关键。van Hoeij LLL 对任意 r 都是多项式复杂度，实测 recombine 只占 1-15%，最小化 r 对 LLL 路径收益极小。

### 3.3 参考实现对比（源码核实）

#### FLINT（`src/fmpz_poly_factor/factor_zassenhaus.c`）

```c
// 外层循环：固定尝试 3 次
for (i = 0; i < 3; i++)
{
    // 内层循环：找下一个满足无平方条件的素数
    for ( ; ; p = n_nextprime(p, 0))
    {
        // lc(f) mod p ≠ 0，f mod p 次数不变，无平方检查
        if (suitable) { DDF+EDF; break; }
    }
    // 记录模因子数最少的结果
    if (temp_fac->num <= r) { ... best = this; }
    p = n_nextprime(p, 0);
}
// r > cutoff(=8) 且 use_van_hoeij=1 时调用 van Hoeij
if (r > cutoff && use_van_hoeij)
    fmpz_poly_factor_van_hoeij(final_fac, fac, f, exp, p);
```

**要点**：
- **仅尝试 3 个素数**（无论 r 多大）
- 仍保留最小-r 策略（但只在 3 个中选，而非 30 个）
- `factor_van_hoeij.c` 本身**不做素数选择**，prime 由调用方传入
- cutoff=8：r ≤ 8 用 Zassenhaus（快），r > 8 用 van Hoeij（处理大 r）

#### NTL（`src/ZZXFactoring.cpp`）

```cpp
NTL_CHEAP_THREAD_LOCAL long ZZXFac_InitNumPrimes = 7;   // 初始试验次数
NTL_CHEAP_THREAD_LOCAL long ZZXFac_MaxNumPrimes = 50;   // 最大试验次数

// SmallPrimeFactorization() 中：
long minr = n + 1;
for (; NumPrimes < ZZXFac_InitNumPrimes;) {
    // DDF+EDF，更新 minr
    if (r < minr) { minr = r; best = this; }
}
```

**要点**：
- **初始试 7 个**，如有必要最多试 50 个，选最小 r
- 策略与 CLPoly 相同（最小-r），但初始次数（7）远少于 CLPoly（30）

### 3.4 分层根因：两个独立的性能差距

对比三方数据：

| 对比 | 差距 | 主因 |
|------|------|------|
| CLPoly vs NTL | **5.9-16x** | Zp 算术：数据结构 + 标量模运算算法（见下） |
| NTL vs FLINT | **1.0-1.4x** | 试验次数（7-50 vs 3）+ FLINT `nmod_poly` 的额外优化 |
| CLPoly vs FLINT | **8-16x** | 上述两者的乘积 |

**关键推论**：NTL 也用最小-r 策略，试 7-50 次；NTL 仍比 CLPoly 快 6-16x。**主导因素是 Zp 算术速度**。

#### Zp 算术速度差距：详细根因分析

**根因 A：标量模乘算法**

CLPoly（`number.hh`）：
```cpp
// Zp operator*:
op1._i *= op2._i;   // 64-bit multiply
op1._i %= op1._p;   // x86 div 指令 ≈ 20-40 cycles
// Zp operator+:
op1._i += op2._i;
op1._i %= op1._p;   // 加法后也用 div！≈ 20-40 cycles
```

FLINT（`nmod.h`，`NMOD_MUL_PRENORM` 宏）：
```c
// 预计算: ninv = -1/n mod 2^64 (Barrett 常数)
umul_ppmm(p_hi, p_lo, a, b);      // 128-bit 乘法，≈ 3 cycles
umul_ppmm(q1, q0, ninv, p_hi);    // 商估计
add_ssaaaa(q1, q0, ...);           // 128-bit 加法
r = p_lo - (q1 + 1) * n;          // 条件修正，无 div 指令
// 总计 ≈ 6-8 cycles（无除法）
```

每次标量模乘：CLPoly ≈ 25 cycles，FLINT ≈ 7 cycles，**差距约 3-4x**。

NTL `zz_pX` 使用类似的 Montgomery/Barrett 归约，差距与 FLINT 相近。

**根因 B：数据结构布局**

CLPoly `upolynomial_<Zp>`：
```
std::vector<std::pair<umonomial, Zp>>
  每个元素 = { uint64_t deg, uint64_t _i, uint32_t _p, uint32_t padding }
           = 24 bytes/系数，且每个元素都存储 prime（冗余）
```

FLINT `nmod_poly`：
```c
struct { ulong* coeffs; slong len; nmod_t mod; }
  // coeffs = 紧密 uint64_t 数组，8 bytes/系数
  // mod（含 Barrett 常数）只存一次
```

对度为 n 的 Zp 多项式：
- CLPoly：24n bytes，非紧密，prime 冗余 n 次
- FLINT：8n bytes，紧密连续，cache 友好
- **内存占用 3x，缓存效率更差**

**根因 C：多项式乘法与 mod 例程**

CLPoly `__upoly_mod`（`polynomial_factorize_zp.hh §4.1.2`）：
```cpp
pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
```
使用通用 `pair_vec_div`——为整系数精确除法设计的例程，做了大量指针跳转和 degree 比较，**不是为密集 Zp 多项式优化的**。

CLPoly `__upoly_powmod`（`§4.1.7`）：
```cpp
ZZ e = exp;
while (e > 0) {
    if (e % 2 != 0) { result = result * b; result = __upoly_mod(result, modpoly); }
    e = e / 2;   // ZZ 整数除法，非位移
    b = b * b; b = __upoly_mod(b, modpoly);
}
```
- 二进制快速幂，正确，但 `e = e / 2` 经过 `ZZ::operator/`（非 `>> 1`）
- 每次 `__upoly_mod` 都分配新的 `vector`（堆分配/释放）
- 每次 `result * b` 生成临时对象

FLINT `nmod_poly_powmod_ui_binexp`（标准库实现）：
- 直接操作 `ulong` 数组，原地更新，无堆分配
- 多项式乘法内联 Barrett 归约

DDF 每次迭代的主要操作是 `powmod(h, p, f_star)`，对度-n 多项式：
- `log₂(p)` 次平方 × 每次平方 ≈ O(n²) 次标量乘法
- 加上 `log₂(p)` 次多项式 mod（另一个 O(n²)）
- n=15, p≈50: `log₂(50) ≈ 6`，每次 powmod ≈ 12 × 225 = 2700 次标量乘法
- 加上 GCD（Euclidean，≈ 3-4 次 poly div）：约 1000 次标量乘法
- **每次 DDF ≈ 3700 次标量乘法**

CLPoly per DDF：3700 × 25 cycles = 92,500 cycles ≈ 0.031ms @3GHz
加上数据结构、堆分配、`pair_vec_div` 开销，实测 0.133ms/trial（合理）

FLINT per DDF（估算）：3700 × 7 cycles + cache 优势 ≈ 26,000 cycles ≈ 0.009ms
FLINT 3 trials × 0.009ms ≈ 0.027ms，加 Hensel+recombine ≈ 0.633ms（FLINT 实测，说明 Hensel 提升本身也更快）

**根因 D：DDF 算法（大度数时）**

- CLPoly：标准 DDF（逐次计算 x^(p^d) mod f，O(n²) per degree）
- FLINT：提供 `nmod_poly_factor_kaltofen_shoup`（Kaltofen-Shoup 1998，baby-step-giant-step DDF），复杂度 O(√n) 次 poly 乘法，对大度数（n ≥ 64）有显著优势
- 对当前测试用例（n ≤ 24）影响有限，但对更大多项式会拉大差距

---

## 4. 推荐优化方向与优先级

### P2a：减少素数试验次数（极低难度，中等收益）

> ✅ **[已实现 2026-02-23]** `max_tries` 已从 30-50 改为 3，移除了 `> deg/2` 时扩展到 50 的逻辑。与 FLINT 对齐。

**方案**：参照 FLINT，将 `max_tries` 从 30-50 降到 3-5。仍保留最小-r（在少量试验中选最优）。

**预期收益**（实现前预估）：
- `sel_prime` 时间减少约 **85-90%**（30次→3次 DDF+EDF）
- 总时间减少约 **65-75%**
- CLPoly vs FLINT 差距：从 8-16x 缩小到约 **2-4x**（剩余差距来自 Zp 算术，即根因 A/B/C）

**实施难度**：极低——`max_tries = 30` 改为 `max_tries = 3`，约 2 行改动。

**风险**：r 可能略大于全局最优，但 recombine 只占 1-15%，FLINT 实践验证 3 次足够。

> **[后续 profiling 待做]**：P2a 实施后应重新 profiling 确认剩余差距，以决定是否推进 P2b。（当前 §1.2 profiling 数据对应 max_tries=30，P2a 后数据将显著不同。）

### P2b：原生 Zp 多项式算术（中等难度，大收益）

**方案**：实现密集 `uint64_t` 系数 Zp 多项式类型（参照 FLINT `nmod_poly`），用于 `__select_prime` 内的 DDF、EDF、GCD 运算：
1. `nmod_upoly`：`std::vector<uint64_t>` 系数 + Barrett 常数（`ninv`, `norm`）存一次
2. 标量模乘改用 Barrett 归约（消除 `div` 指令）
3. 多项式乘法改为密集数组操作（原地，无堆分配）
4. `powmod` 用位移代替 `ZZ::operator/`，减少临时对象

**预期收益**：
- 根因 A（标量乘法）：3-4x
- 根因 B（数据结构/缓存）：1.5-2x
- 根因 C（例程优化）：1.5-2x
- 综合估计：**5-10x per DDF 调用**，接近 NTL 水平
- 与 P2a 叠加：总差距可望缩小到 FLINT 的 **1.5-2x** 以内

**实施难度**：中等，需新建 `nmod_upoly` 类型并修改 DDF/EDF/GCD 实现。

**建议时机**：P2a 完成并 profiling 后视剩余差距决定。

### P1a/P1b/Monagan 2022：维持现状

- **P1a（LLL 重组）**：已实现，消除 2^r 爆炸，维持现状
- **P1b（线性 Hensel）**：基础设施已有，完整集成 ROI 极低（影响 <10% 的步骤），暂缓
- **Monagan 2022**：影响范围更小（hensel 子步骤），暂缓

---

## 5. 结论

| 根因 | 具体表现 | 对应差距 | 修复方案 | 难度 |
|------|----------|----------|----------|------|
| 试验次数过多 | ~~30-50 次~~ → **已修复（max_tries=3）** | **~10x** 浪费 | ✅ P2a：max_tries=3 | 极低 |
| 标量模乘使用 `div` 指令 | `%` ≈ 25 cycles vs Barrett ≈ 7 cycles | **3-4x** per 乘法 | P2b A：Barrett 归约 | 中 |
| 稀疏 pair-vector 数据结构 | 24 bytes/coeff vs 8 bytes/coeff | **1.5-2x** cache 效率 | P2b B：密集 `uint64_t[]` | 中 |
| 通用 `pair_vec_div` 做 poly mod | 非 Zp 专用，含冗余开销 | **~1.5x** per 操作 | P2b C：专用 poly mod | 中 |
| 合计 | | CLPoly vs FLINT **8-16x** | P2a + P2b | — |

**行动计划**：
1. **立即**：实施 P2a（2 行改动），再次 profiling 确认效果
2. **评估后**：若 vs FLINT 仍有 >3x 差距，实施 P2b（原生 Zp 算术），预期接近 NTL 水平
3. **长期**：P1b / Monagan 2022 仅在超大精度用例（a_mig > 100）下有实际意义
