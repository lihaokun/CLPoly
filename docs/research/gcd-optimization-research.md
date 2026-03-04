# 调研报告：GCD 算法优化

> 调研日期：2026-02-28
> 目标：定位 CLPoly 多项式 GCD 与 FLINT/NTL 2-10x 性能差距的根因，调研参考系统算法路线，确定优化方向与优先级。
> 一手源码（已直接阅读）：
>   - FLINT `src/fmpz_poly/gcd.c`（选择逻辑）、`gcd_heuristic.c`（GCDHEU）、`gcd_modular.c`（模算法）
>   - FLINT `src/nmod_poly/gcd.c`（Zp GCD 含 HGCD 阈值表）、`hgcd.c`（HGCD 实现）
>   - FLINT `src/fmpz_mpoly_factor/gcd_algo.c`（多变量 GCD 选择逻辑）
>   - NTL `src/ZZX1.cpp:2838-2919`（ZZX::GCD 实现）
> 已读论文（全文）：
>   - Brown 1971, "On Euclid's Algorithm and the Computation of Polynomial Greatest Common Divisors", J. ACM 18(4), pp.478-504 — 模 GCD 算法原始论文，复杂度分析
>   - Zippel 1979, "Probabilistic Algorithms for Sparse Polynomials", EUROSAM'79, pp.216-226 — 稀疏模 GCD 算法原始论文
>   - Parisse 2002, "A correct proof of the heuristic GCD algorithm" (arxiv:cs/0206032v1, 5pp) — GCDHEU 正确性证明 + 界
>   - Monagan 2022, "Speeding up polynomial GCD, a crucial operation in Maple" (Maple Trans. 2(1), 18pp) — Hu/Monagan 算法完整描述 + benchmark
> 参考文献（二手/摘要）：Char-Geddes-Gonnet 1989（付费墙）

---

## 1. 现有基础设施与实测数据

### 1.1 当前 GCD 管线架构

**单变量 ZZ (`polynomial_gcd.cc`)**：
```
polynomial_GCD(F, G)  — upolynomial_<ZZ>
  ├─ 边界处理: 空/单项式
  ├─ 提取 content: cont(F), cont(G), cont_gcd = gcd(f_cont, g_cont)
  ├─ 除以 content → 本原化
  ├─ lc_gcd = gcd(lc(F), lc(G))
  └─ 模 GCD + CRT 主循环:
       ├─ 跳过整除 lc 的素数
       ├─ polynomial_mod → Zp 多项式
       ├─ __polynomial_GCD(Zp) — 朴素欧几里德
       ├─ 度数跟踪 (Pout_d): 度降 → 重置; 度不变 → CRT 合并
       ├─ 对称模约简
       └─ 稳定 → cont → 试除 F,G → 返回
```

**多变量 ZZ (`polynomial_gcd.hh:246-468`)**：
```
polynomial_GCD(F, G)  — polynomial_<ZZ, lex>
  ├─ 边界处理: 空/单项式/变量不同
  ├─ 提取 content: F_cont, G_cont, cont_gcd = polynomial_GCD(cont(F), cont(G))
  ├─ 本原化 + lc_gcd
  ├─ 素数选择: p_index = max(deg) / ln(max(deg)), prime > max(deg)
  └─ 模 GCD + CRT 主循环:
       ├─ polynomial_mod → Zp
       ├─ __polynomial_GCD(Zp, multivar) — 递归求值+Lagrange 插值
       ├─ CRT 合并 + 对称模约简
       └─ 稳定 → cont → 试除 → 返回
```

**多变量 Zp (`polynomial_gcd.hh:833-1074`)**：
```
__polynomial_GCD(Pout, F, G, Lc_gcd, deg)  — polynomial_<Zp, lex>
  ├─ 若为同一变量（单变量）→ 朴素欧几里德
  └─ 若为多变量:
       ├─ 选择自由变量 v（最小 lex 排序的非首变量）
       ├─ v_d = max(deg_v(F), deg_v(G)) + 1（需要的求值点数）
       ├─ Fisher-Yates 随机抽样求值点（O(prime) 空间）
       ├─ 对每个求值点 α: assign(F,v,α) → 递归 __polynomial_GCD
       ├─ 度数跟踪 + 收集求值结果
       └─ 收够 v_d 个点 → Lagrange 插值重建
```

### 1.2 实测 Benchmark 数据（release build, 2026-02-28）

**多变量 GCD（2 变量, polynomial_ZZ）**：

| 用例 | 时间 | 内存 |
|------|------|------|
| gcd deg8+common4 | 0.428ms | 18.0KB |
| gcd deg15+common8 | 0.979ms | 29.5KB |
| gcd deg25+common12 | 7.154ms | 44.5KB |

**单变量 GCD（upolynomial_ZZ）**：

| 用例 | 时间 | 内存 |
|------|------|------|
| gcd deg30+common15 | 0.086ms | 13.7KB |
| gcd deg80+common40 | 0.655ms | 40.7KB |
| gcd deg200+common100 | 4.745ms | 33.9KB |

### 1.3 与 FLINT/NTL 的差距

**多变量 GCD（vs FLINT `fmpz_mpoly_gcd`）**：

| 用例 | CLPoly | FLINT | 倍数 | CLPoly 内存 | FLINT 内存 |
|------|--------|-------|------|------------|-----------|
| deg8+common4 | 0.335ms | 0.068ms | **4.92x** | 12.8KB | 248.8KB |
| deg15+common8 | 1.477ms | 0.172ms | **8.58x** | 19.5KB | 269.9KB |

**单变量 GCD（vs NTL `GCD`）**：

| 用例 | CLPoly | NTL | 倍数 | CLPoly 内存 | NTL 内存 |
|------|--------|-----|------|------------|---------|
| deg50+common25 | 0.190ms | 0.090ms | **2.11x** | 19.3KB | 169.6KB |
| deg200+common100 | 5.072ms | 0.532ms | **9.53x** | 34.3KB | 166.2KB |

**趋势**：差距随度数增长显著扩大（单变量从 2x → 10x；多变量从 5x → 9x），说明瓶颈在算法复杂度而非常数开销。

**注**：CLPoly 内存占用远低于 FLINT/NTL（0.05-0.2x），说明 CLPoly 的稀疏表示在内存上有优势，但在 GCD 计算中反而因数据结构开销拖慢速度。

---

## 2. 根因分析

### 2.1 根因 A：单变量 Zp GCD 使用朴素欧几里德（无 HGCD）

**位置**：`polynomial_gcd.hh:1093-1122`

```cpp
// 朴素欧几里德循环 — O(n^2) 域运算
Pout=G;
pair_vec_div(Pout_2.data(),Pout_1.data(),F.data(),Pout.data(),Pout.comp());
while(!Pout_1.empty())
{
    std::swap(Pout.data(),Pout_.data());
    std::swap(Pout.data(),Pout_1.data());
    pair_vec_div(Pout_2.data(),Pout_1.data(),Pout_.data(),Pout.data(),Pout.comp());
}
```

CLPoly 的 Zp GCD 使用朴素欧几里德算法，每步调用 `pair_vec_div` 做多项式除法。复杂度为 **O(n^2)** 域运算（n = 多项式度数）。

**FLINT HGCD 阈值**（源码：`src/nmod_poly/gcd.c` 中 `nmod_poly_gcd_hgcd_cutoff_tab[64]`）：

FLINT 的 HGCD 阈值**不是固定值**，而是按模数 bit-size 查表。`_nmod_poly_gcd` 的选择逻辑：

```c
// FLINT src/nmod_poly/gcd.c:62-74
slong _nmod_poly_gcd(..., nmod_t mod) {
    if (lenB < nmod_poly_gcd_hgcd_cutoff(mod))  // 按 bit-size 查表
        // 欧几里德（或 redc_fast 变体）
    else
        return _nmod_poly_gcd_hgcd(G, A, lenA, B, lenB, mod);
}
```

实际阈值范围 470-2199，随模数 bit-size 变化（64 个条目，索引 = NMOD_BITS - 1）：
- 1-bit 模数：470
- 32-bit 模数：1170
- 63-bit 模数：1289
- **64-bit 模数（模 GCD 实际使用）：1725**

**关键发现**：FLINT 的模算法 `gcd_modular.c` 使用 64-bit 素数（`p = 2^63` 起始），对应 HGCD 阈值 **1725**。因此在 deg200 时，FLINT **同样使用欧几里德算法**，HGCD 并非当前 benchmark 用例（deg ≤ 200）的性能差距来源。

HGCD 仅在 deg ≥ 1725 时对 64-bit 素数生效。对于更小的素数（如 32-bit），阈值更低（~1170），但模算法不使用这些小素数。

此外，`pair_vec_div` 在稀疏 pair vector 上操作（`basic.hh:614-695`），对密集 Zp 多项式（每项都非零）引入不必要的单项式比较和堆操作。密集 Zp 数组可用简单指针减法替代。

**影响**：对当前 benchmark 用例（deg ≤ 200），HGCD **不是差距来源**。差距主要来自 pair_vec_div 开销（约 **1.5-2x**）和其他根因。HGCD 在 deg ≥ 2000 时才有显著价值。

### 2.2 根因 B：单变量 ZZ GCD 从 prime=2 开始

**位置**：`polynomial_gcd.cc:43-44`

```cpp
std::uint32_t p_index=0;
std::uint64_t prime=boost::math::prime(p_index);  // prime = 2
```

单变量 ZZ GCD 从最小素数 2 开始。小素数贡献的比特数少（prime=2 仅贡献 1 bit），导致 CRT 乘积增长缓慢，需要更多素数迭代才能达到系数界。

**FLINT 的做法**（源码：`src/fmpz_poly/gcd_modular.c:117-118`）：

```c
// FLINT gcd_modular.c
pbits = FLINT_BITS - 1;
p = (UWORD(1)<<pbits);   // p = 2^63 ≈ 9.2×10^18
// 然后: p = n_nextprime(p, 0);  每个素数贡献 ~63 bits
```

FLINT 从 **2^63** 开始选素数，每个素数贡献约 63 bits。系数界使用 Yap 的公式：
`bound = (n0+3)*max(nb1,nb2) + (n0+1)`（源码 line 139），其中 nb 是 2-范数比特数。

**NTL 的做法**（源码：`src/ZZX1.cpp:2879`）：

```cpp
zz_p::FFTInit(i);  // 选择第 i 个 FFT 素数（形如 k·2^m + 1，约 60 bits）
```

NTL 使用 FFT 素数（支持 NTT 加速乘法），同样是 ~60 bits 的大素数。

对比**CLPoly 多变量版本**（`polynomial_gcd.hh:310-314`）从较大素数开始：

```cpp
std::uint32_t p_index=tmp_x/std::log(tmp_x);  // deg200 → p_index≈37, prime≈157
```

对于 deg200、系数界 B ≈ 200 bits 的多项式：
- **CLPoly 单变量** 从 prime=2：需约 200 个素数（前 200 个素数的 log₂ 和 ≈ 200 bits）
- **CLPoly 多变量** 从 prime≈157：需约 27 个素数（每个 ~7.3 bits）
- **FLINT** 从 prime≈2^63：仅需 **4 个素数**（每个 ~63 bits，4×63=252 > 200）
- **NTL** 从 FFT 素数（~60 bits）：仅需 **4 个素数**

CLPoly 单变量需要约 **50x** 更多的 Zp GCD 迭代（200 vs 4）。考虑到提前终止启发式，实际差距约 **5-10x**。这是 deg200 单变量 9.53x 差距的**首要根因**。

**Brown 1971 的理论依据**：Brown 在 §4.4 eq.43 中明确推荐素数选自区间 `α = ½(β+1) < p ≤ β`，其中 β 是机器字长最大整数（64-bit 机器即 2^63）。Brown 的 eq.94 证明：所需素数个数 `n̄ ≤ 4l + 2`，其中 l 是系数最大长度（以 α 为底的对数）。对于 64-bit 素数，l ≈ 1（系数长度相对于 2^63 极小），因此仅需 ~6 个素数。CLPoly 从 prime=2 开始严重违背了 Brown 的设计意图。

**影响**：约 **3-8x**（是最大的单一差距来源）。

### 2.3 根因 C：多变量无 Zippel/Brown 分算法

**位置**：`polynomial_gcd.hh:880-1072`

CLPoly 的多变量 Zp GCD 使用单一策略：对自由变量做随机求值 → 递归 GCD → Lagrange 插值。本质是 **密集 Brown 算法**（Brown 1971, §4.5 Algorithm P），但实现有以下低效点：

**(a) Fisher-Yates 全域扫描（lines 907-918）**：

```cpp
std::vector<int> v_bool(prime, 1);  // O(prime) 空间
std::random_device rd;               // 硬件熵，昂贵
std::mt19937 gen(rd());
```

为每次 GCD 调用分配大小为 `prime` 的布尔向量，并使用硬件随机设备。对于 prime ≈ 200，这不是大问题，但对于大素数则引入不必要的开销。

**(b) Lagrange 插值 O(v_d^2) 空间 + O(v_d^2 · #mono) 时间（lines 1000-1056）**：

对于 deg8 的双变量多项式，v_d ≈ 9，需要 9 个递归 GCD 图像 + 81 元素的 Lagrange 基矩阵。这是 Brown 密集策略的固有代价：对 v 个变量、每变量度数 d 的多项式，需要 (d+1)^v 个独立求值（Brown 1971, §1 "an exponential worst case behavior since they need as many as (d+1)^v independent evaluations"）。

**(c) 无稀疏插值**：Zippel 1979 算法对 t 项的 GCD 仅需 O(t) 求值点（而非 O(d+1)），代价渐近为 O(t³)。Zippel 的关键洞察是：若某幂次在首次单变量 GCD 中系数为零，则以高概率该幂次在完整多项式中恒为零——因此只需插值非零项。这将稀疏多项式的 GCD 从指数复杂度降至多项式复杂度。

**(d) 无早期终止**：总是收集 v_d 个点后一次性插值，即使 GCD 度数比预期低也不提前返回。

**Brown vs Zippel 复杂度对比**（源自原始论文）：

| | Brown（密集） | Zippel（稀疏） |
|---|---|---|
| 求值点数（每变量） | d+1 | t（已知非零项数） |
| 总求值数（v 变量） | (d+1)^v | O(dvt) |
| 渐近代价 | O(l²(d+1)^v + l(d+1)^(v+1))（Brown 1971, eq.95） | O(t³)（Zippel 1979, §3.2），其中 t ≪ (d+1)^v |
| 适用场景 | 密集、变量少 | 稀疏、变量多 |

Zippel 1979 的实验数据（10 变量、deg≤3、9 项 GCD）显示：Modular（Brown）在 7 变量时耗时 2409s，Sparse Modular 仅 4.19s（575x 加速），且 Brown/EZ 在 6+ 变量时内存溢出。

**FLINT 做法**：自动选择 Brown（密集双变量快）、Zippel（稀疏多变量快）、BMA/Zippel2（极稀疏最快），根据估计代价动态切换。

**影响**：多变量 deg15 的 8.58x 差距中，此项贡献约 **2-5x**。

### 2.4 根因 D：`cont(F)` 重复计算

**位置**：`polynomial_gcd.hh:297-299`

```cpp
polynomial_<ZZ,lex_<var_order>> F_cont=cont(F);    // 计算 cont(F)
polynomial_<ZZ,lex_<var_order>> G_cont=cont(G);    // 计算 cont(G)
polynomial_<ZZ,lex_<var_order>> cont_gcd=polynomial_GCD(cont(F),cont(G));  // 重新计算 cont(F) 和 cont(G)！
```

第 299 行调用 `polynomial_GCD(cont(F), cont(G))` 重新计算了 `cont(F)` 和 `cont(G)`，而非使用已有的 `F_cont` 和 `G_cont`。这使 content 提取的代价翻倍。

**修复**：改为 `polynomial_GCD(F_cont, G_cont)`。

**影响**：约 **1.2-1.5x**（取决于 content 计算在总时间中的占比）。

### 2.5 根因 E：无 GCDHEU（启发式 GCD）快速路径

CLPoly 完全缺失启发式 GCD 算法。

**FLINT 的 GCDHEU 实现**（源码：`src/fmpz_poly/gcd_heuristic.c`）：

FLINT 的 GCDHEU **不是**传统的"在大整数 ξ 处求值"方法。它使用 **bit-packing + 整数 GCD**：

```
输入: 本原多项式 A, B（已除以 content）
1. 计算 pack_bits = max(min(bits₁, bits₂) + 6, max(bits₁, bits₂) + 1)
   // +6 是启发式选择；理论上 +3 即可（参考 arxiv:cs/0206032v1）
2. _fmpz_poly_bit_pack(array1, A, len1, pack_bits)  → 密集打包为大整数
   _fmpz_poly_bit_pack(array2, B, len2, pack_bits)  → 密集打包为大整数
3. limbsg = flint_mpn_gcd_full(arrayg, array1, array2)  → GMP 大整数 GCD
4. _fmpz_poly_bit_unpack(G, glen, arrayg, pack_bits)    → 解包为多项式
5. 除以 content(G)
6. 验证: flint_mpn_divides(q, array1, arrayg) 且 flint_mpn_divides(q, array2, arrayg)
   如果整数商精确 → 解包商并检查 bits_G + bits_Q + log(len) < pack_bits
   否则 → 乘法验证 multiplies_out(A, Q, G)
7. 两次验证都通过 → 返回 G；否则返回失败
```

这本质上是在 ξ = 2^pack_bits 处求值（bit-packing 等价于 2 的幂基求值），利用位操作实现 O(1) 打包/解包，比任意基展开高效得多。

**FLINT 的选择条件**（源码：`src/fmpz_poly/gcd.c:61`）：

```c
if (b1 + b2 < 2 * FLINT_BITS)  // 即 bitlen(A) + bitlen(B) < 128
    if (_fmpz_poly_gcd_heuristic(res, poly1, len1, poly2, len2))
        return;  // 成功则直接返回
// 失败则回退到模算法
```

对于我们的 benchmark 用例（系数 ∈ [-20, 20]，bitlen ≈ 5），5+5 = 10 ≪ 128，**FLINT 总是优先尝试 GCDHEU**。单次大整数 GCD 代价远低于模算法的多次素数迭代。

**GCDHEU 正确性界**（论文：Parisse 2002, arxiv:cs/0206032v1, Theorem 1）：

设 P, Q 为整系数多项式，z 为整数满足 |z| ≥ 2·min(|P|, |Q|) + 2（其中 |P| 是 P 的最大系数绝对值）。若 gcd(P(z), Q(z)) 的 z-adic 对称重建的本原部分 G 整除 P 和 Q，则 G = gcd(P, Q)。该定理对高斯整数同样成立。

对于 FLINT 的 bit-packing 方案（ξ = 2^pack_bits），pack_bits = min(bits₁, bits₂) + 6 保证了 ξ 远超上述界（因为 ξ ≥ 2^(bits+6) ≫ 2·max|coeff| + 2），所以 bit-packing 实质上选择了一个远超理论需求的求值点。

**影响**：对小/中系数用例（deg8+common4、deg30+common15 等），预期 **2-5x** 提速。

### 2.6 根因 F：CRT 合并循环冗余计算

**位置**：`polynomial_gcd.cc:106-161`、`polynomial_gcd.hh:374-429`

CRT 更新公式：`coeff_new = coeff_old + (coeff_mod - coeff_old) * inv * Pout_prime`

其中 `inv * Pout_prime` 对所有系数不变，但在循环中每个系数都重新计算 `tmp_inv.number() * Pout_prime`。

```cpp
// .cc:119 — 每个系数都计算 Pout_ptr->second * tmp_inv.number() * Pout_prime
tmp_Pout_.push_back({Pout_ptr->first,
    Pout_ptr->second - Pout_ptr->second * tmp_inv.number() * Pout_prime});
```

此外，`Pout_prime / 2`（对称约简界）也应预计算。

**影响**：约 **1.2-1.3x**。

### 2.7 根因 G：Zp 对象 32 字节/元素（缓存低效）

**位置**：`number.hh` Zp 类定义

每个 `Zp` 值存储：`_i`(8B) + `_p`(8B) + `_ninv`(8B) + `_norm`(4B) ≈ 32 字节。FLINT 的 `nmod` 仅存值 (8B)，模上下文单独传递。

200 项 Zp 多项式：CLPoly 6400B vs FLINT 1600B（系数部分），4x 内存膨胀直接影响 L1/L2 缓存命中率。

**影响**：约 **1.2-1.5x**（带宽/缓存效应在大多项式上更明显）。

### 2.8 根因 H：`pair_vec_div` 每次堆分配

**位置**：`basic.hh:614-619`

```cpp
VHC<...> **heap = new VHC<...>*[v2_.size()-1];
VHC<...> *node  = new VHC<...>[v2_.size()-1];
VHC<...> **lin  = new VHC<...>*[v2_.size()-1];
// ... 使用 ...
delete [] lin; delete [] node; delete [] heap;
```

欧几里德 GCD 循环中每步除法做 3 次堆分配 + 3 次释放。deg200 的 GCD 约需 O(n) 步，共 ~600 次 new/delete。

**影响**：约 **1.1-1.3x**。

### 2.9 根因 I：多变量 `assign` 无幂次预计算

**位置**：`polynomial.hh:720-723`

```cpp
for (auto& j:i.first)
    if (j.first==v)
        z*=pow(c, j.second);  // 每项独立计算 c^k
```

对每个单项式独立调用 `pow(c, j.second)` 计算幂次。若预计算 `c_powers[k] = c^k`（k = 0..d），每项仅需一次查表乘法。

**影响**：约 **1.1-1.2x**。

---

## 3. 参考系统算法调研

### 3.1 FLINT

#### 3.1.1 单变量 ZZ GCD (`fmpz_poly_gcd`)

**算法选择逻辑**（源码：`src/fmpz_poly/gcd.c:16-69`）：

```c
void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, slong len1,
                    const fmpz * poly2, slong len2) {
    // 先剥离公因子 x^k
    if (len1 < 6)
        _fmpz_poly_gcd_subresultant(res, poly1, len1, poly2, len2);
    else {
        b1 = FLINT_ABS(_fmpz_vec_max_bits(poly1, len1));
        b2 = FLINT_ABS(_fmpz_vec_max_bits(poly2, len2));
        if (b1 + b2 < 2 * FLINT_BITS)  // < 128 bits
            if (_fmpz_poly_gcd_heuristic(res, ...)) return;
        _fmpz_poly_gcd_modular(res, ...);
    }
}
```

三级策略：

| 条件 | 算法 | 源码文件 |
|------|------|---------|
| deg < 6 | Subresultant（精确，无模运算） | `gcd_subresultant.c` |
| bitlen(A)+bitlen(B) < 128 | **GCDHEU**（bit-packing + mpn_gcd） | `gcd_heuristic.c` |
| 默认 | **模算法 + CRT** | `gcd_modular.c` |

**模算法关键细节**（源码：`src/fmpz_poly/gcd_modular.c`）——本质上是 Brown 1971 §4.3 Algorithm M 的优化实现：

1. **素数起始**：`p = 2^63`（line 117-118），然后 `n_nextprime(p, 0)`。每个素数 ~63 bits。这与 Brown §4.4 eq.43 推荐的 `p ∈ [½(β+1), β]` 一致
2. **系数界**：Yap 公式 `(n0+3)*max(nb1,nb2) + (n0+1)`（line 139），其中 nb 是 2-范数比特数。对小多项式还精确计算 `_fmpz_vec_dot` 得到准确 2-范数。对应 Brown §4.3 Step (5) 的 `μ̄ = 2ḡ·max|φ|`
3. **启发式小界**：求值 A(-1)、B(-1) 的 GCD 作为提前终止条件（line 88-113）
4. **CRT 优化**：`_fmpz_poly_CRT_ui` 是针对单精度模数优化的增量 CRT（line 234）。对应 Brown §4.8 CRA
5. **每素数 Zp GCD**：调用 `_nmod_poly_gcd`，对 64-bit 素数在 deg < 1725 时使用**欧几里德**（非 HGCD）。对应 Brown §4.7 Algorithm U

**HGCD 与 Zp GCD 实际阈值**（源码：`src/nmod_poly/gcd.c:37-44`）：

每素数的 Zp GCD 算法选择取决于 `nmod_poly_gcd_hgcd_cutoff_tab[NMOD_BITS(mod)-1]`：
- 对 64-bit 素数（模 GCD 使用）：阈值 = **1725**（欧几里德 → HGCD）
- 在 HGCD 内部，基础情形阈值 `nmod_poly_hgcd_iter_recursive_cutoff_tab` 在 68-244 范围
- 低于 HGCD 阈值时，若满足 `NMOD_POLY_GCD_EUCLIDEAN_USE_REDC_FAST` 则使用 Montgomery REDC 优化的欧几里德

结论：对当前 benchmark 的 deg ≤ 200，**FLINT 同样使用欧几里德算法**做 Zp GCD。性能差距的首要原因是素数大小（4 次 vs 200 次迭代），其次是 GCDHEU 快速路径。

#### 3.1.2 多变量 ZZ GCD (`fmpz_mpoly_gcd`)

**源码**：`src/fmpz_mpoly_factor/gcd_algo.c`（Daniel Schultz, 2018-2021）

**早期退出链**（`_fmpz_mpoly_gcd_algo_small`, lines 1426-1761）：

```
1. 单项式 A 或 B → _do_monomial_gcd（O(1)）
2. ess(A) == ess(B) → _try_monomial_cofactors（检测余因子为单项式）
3. A,B 无共同变量 → _do_trivial
4. 仅一个共同变量 → _do_univar（退化为单变量 GCD）
5. 某变量仅在 A 或仅在 B 中 → _try_missing_var
6. _set_estimates: 在大素数下随机求值，估计 GCD 度数/项数/密度
7. 所有度数界=0 → _do_trivial
8. 度数界匹配 A 或 B → _try_divides（整除性检查）
9. 算法选择（见下）
```

每个早期退出都避免了昂贵的完整 GCD 计算。CLPoly 完全缺失步骤 1-3、5、7-8。

**算法选择逻辑**（lines 1676-1720）：

```c
mpoly_gcd_info_measure_brown(I, ...);
mpoly_gcd_info_measure_bma(I, ...);

if (I->mvars < 3) {
    _try_brown(G, Abar, Bbar, A, B, I, ctx);  // 双变量优先 Brown
} else if (zippel2_time < brown_time &&
           (density < 1 || zippel2_time < 0.01*brown_time)) {
    _try_bma(...) || _try_brown(...);  // 稀疏优先 BMA
} else {
    _try_brown(...) || _try_bma(...);  // 密集优先 Brown
}
// 最终兜底: _try_zippel(...)
```

| 算法 | 适用场景 | 入口函数 |
|------|---------|---------|
| **Brown** | 密集双变量（vars < 3）或密集多变量 | `_try_brown` (line 1352) |
| **BMA/Zippel2** | 稀疏多变量（zippel2_time < brown_time） | `_try_bma` (line 1093) |
| **Zippel** | 兜底/显式指定 | `_try_zippel` (line 967) |
| **Hensel** | 显式指定 | `_try_hensel` (line 1233) |
| **PRS** | 显式指定 | `_try_prs` (line 831) |

**FLINT 多变量 Zp GCD 本质是 Brown 1971 §4.5 Algorithm P 的优化实现**。Brown 证明 Algorithm P 的复杂度为 `P(v,d) ≲ (d+1)^(v+1)`（§5.6 eq.91），而经典算法 Algorithm C 的复杂度为 `C(v,l,d) ≲ l²(d+1)^(4v)·2^(2v²)·3^v`（§5.5 eq.80）。Brown §5.8 关键结论：**M（模算法）严格被 C（经典算法）的第一次伪除法支配**（"the maximum computing time for the modular algorithm is strictly dominated by the maximum computing time for the first pseudo-division in the classical algorithm"），证明 `C_/M+ = (d+1)^(v-2)`，对 v ≥ 2 成立。

**求值优化——LUT 幂次预计算**（`fmpz_mpoly_evals`, lines 30-180）：

```c
// 预计算 alpha[j]^k 的查找表（直接 LUT 或二幂 LUT）
if (use_direct_LUT) {
    LUTvalue[j][k] = alpha[j]^k  // k = 0..Amax_exp[j]
} else {
    // 二幂 LUT: alpha[j]^(2^i), 然后组合
}
// 对每个单项式: meval = ∏ LUTvalue[j][varexp] — 查表乘法
```

这正是 CLPoly 在 §2.9 中缺失的优化。FLINT 根据最大指数自动选择直接 LUT 或二幂 LUT。

#### 3.1.3 Zp 算术

**源码**：`src/nmod_poly/gcd.c`、`src/nmod.h`

FLINT 的 `nmod` 类型仅存值（8 字节 `ulong`），模上下文 `nmod_t` 作为参数传递：

```c
typedef struct { ulong n; ulong ninv; flint_bitcnt_t norm; } nmod_t;
// n = 模数, ninv = 预计算逆元（Barrett 约简）, norm = 规范化移位
```

对于欧几里德 GCD 有两种变体（`_nmod_poly_gcd`, line 62-74）：
- **标准欧几里德**：`_nmod_poly_gcd_euclidean` — 用 `_nmod_poly_rem` 做多项式取余
- **REDC-fast 欧几里德**：`_nmod_poly_gcd_euclidean_redc_fast` — 先做一次标准取余平衡度数，再转为 Montgomery REDC 表示做后续取余（条件：模数奇数且 ≤ 62 bits）

REDC 变体避免了每次乘法后的 Barrett 除法，改用更快的 Montgomery 约简。源码 `gcd.c:313-356` 显示：先做初始余式（标准表示），然后调用 `gr_ctx_init_nmod_redc_fast` 切换到 REDC，最后 `nmod_redc_fast_normalise` 转回标准表示。

### 3.2 NTL

**源码**：`src/ZZX1.cpp:2838-2919`（`GCD` 函数）

NTL 使用**模算法 + CRT**（无 GCDHEU、无 HGCD）：

```cpp
// NTL src/ZZX1.cpp:2878-2913
for (i = 0; ;i++) {
    zz_p::FFTInit(i);          // 选择第 i 个 FFT 素数
    long p = zz_p::modulus();  // 大素数，~60 bits，形如 k·2^m + 1
    if (divide(LeadCoeff(f1), p) || divide(LeadCoeff(f2), p)) continue;  // 跳过坏素数
    conv(F1, f1); conv(F2, f2);  // Z → Zp
    GCD(G, F1, F2);              // Zp GCD（欧几里德或 HGCD，按 deg 选择）
    mul(G, G, LD);               // 乘以 LC 的 GCD
    if (deg(G) == 0) { set(res); break; }  // 互素
    if (FirstTime || deg(G) < deg(g)) {
        conv(prod, p); BalCopy(g, G);  // 新的度数界，重置 CRT
    } else if (deg(G) > deg(g)) continue;  // 坏素数
    else if (!CRT(g, prod, G)) {          // CRT 合并，返回 false 表示稳定
        PrimitivePart(res, g);
        if (divide(f1, res) && divide(f2, res)) break;  // 试除验证
    }
}
```

关键设计：
- **FFT 素数**：`zz_p::FFTInit(i)` 选择形如 k·2^m + 1 的大素数（~60 bits），支持 NTT 加速多项式乘法。每个素数贡献 ~60 bits 到 CRT
- **CRT 稳定检测**：`CRT(g, prod, G)` 返回 false 当系数未变化 → 触发试除验证。比固定界更高效
- **Zp GCD**：NTL 的 `GCD(G, F1, F2)` 内部也使用 HGCD（阈值更低），但对 deg200 可能仍用欧几里德

NTL 优于 CLPoly 的主要原因：**大 FFT 素数**（4 次 vs 200 次迭代）+ 基于 NTT 的快速多项式乘法 + 优化的 CRT 实现。

### 3.3 Maple / Hu-Monagan 算法

**论文**：Monagan 2022, "Speeding up polynomial GCD, a crucial operation in Maple", Maple Trans. 2(1)

Maple 的 GCD 算法演进：

| 时期 | 算法 | 参考 |
|------|------|------|
| 1980s-2005 | Wang EEZ-GCD | Wang 1980 |
| 2005-2020 | **Zippel 稀疏模 GCD** | Zippel 1979, de Kleine+ 2005 |
| 2021+ | **Hu/Monagan**（新默认） | Hu-Monagan 2021, Monagan 2022 |

**算法核心流程**（论文 §3, Algorithms GCD → MGCD → PGCD）：

```
Algorithm GCD(A, B):
  1. ca, cb = cont(A, x₀), cont(B, x₀)
  2. A, B = A/ca, B/cb  （本原化）
  3. repeat: P = MGCD(A, B) until P|A and P|B
  4. return P × gcd(ca, cb)

Algorithm PGCD(A, B):  // 第一个素数，确定支撑
  1. 选择 Kronecker 替换 Kr: Z[x₀,...,xₙ] → Z[x, y]，ri > deg(A, xi)
  2. 选择光滑素数 p > ∏ri（需要 p-1 无大素因子，用于离散对数）
  3. 选择 Zp* 的生成元 α
  4. 计算 gj = gcd(Kr(A)(x, αʲ), Kr(B)(x, αʲ)) for j = 1,2,...
  5. 同时尝试插值 Kr(H) = LC(Ā,x)·Kr(G) 和 Kr(C) = LC(G,x)·Kr(Ā)
     在 j ∈ {4,6,8,10,...} 时调用 Modified Ben-Or/Tiwari
  6. 输出支撑较小者（H 或 C）

Algorithm MGCD(A, B):  // 后续素数，CRT 恢复系数
  1. (cof, H₁, p₁) = PGCD(A, B)
  2. for k = 2,3,...: 选新素数 pk
     如果 cof=false: Hk = SGCD(Ak, Bk, H₁)  // Zippel 插值 G
     如果 cof=true:  Hk = SCOF(Ak, Bk, H₁)  // Zippel 插值 Ā
  3. CRT 合并 → 对称模约简 → 稳定后除以 content → 返回
```

**关键技术细节**：

1. **Kronecker 替换**：Kr(f) = f(x, y, y^r₁, y^(r₁r₂), ...)，将 n+1 变量映射到 2 变量。需要素数 p > ∏ri（对 8 变量 deg10：p > 11⁸ ≈ 2.14×10⁸，可用 62-bit 硬件运算）

2. **Modified Ben-Or/Tiwari 插值**（论文 §2.1）：
   - 选光滑素数 p，取 Zp* 生成元 α
   - 从求值序列 f(αʲ) 用 Berlekamp-Massey 算法计算最小多项式 λ(z)
   - 对 λ(z) 做 Zp 上因式分解得到根 mi = α^(ei)
   - 解离散对数得指数 ei（用 Pohlig-Hellman，需 p-1 光滑）
   - 解移位 Vandermonde 系统得系数 ai

3. **同时插值 G 和 Ā**（论文 §3.3，核心创新）：
   - PGCD 同时计算 gj（GCD 图像）和 dj = aj/gj（余因子图像）
   - 对两者同时尝试 BT 插值，**哪个先成功就用哪个**
   - 代价取决于 t = min(max(#hi), max(#ci))，即 G 与 Ā 中项数较少者
   - 当 #Ā ≪ #G 时，插值 Ā 比 G 快得多

4. **求值优化**（论文 §4.2）：
   - 预计算单项式求值 mi = Mi(β₁,...,βn)
   - 利用几何序列：Mi(β₁ʲ,...,βnʲ) = miʲ
   - 迭代更新 Ci ← Ci × mi，避免重复幂次计算
   - 实现为 C 子程序 `getsupport64s`

**论文 Benchmark（Table 1，8 变量稀疏多项式）**：

| #G | #Ā | Maple Zippel | Hu/Monagan GCD2 | 加速比 |
|----|----|-------------|-----------------|--------|
| 10¹ | 10³ | 1.227s | 0.067s | **18x** |
| 10³ | 10¹ | 26.97s | 0.094s | **287x** |
| 10⁴ | 10¹ | 1041s | 0.454s | **2293x** |
| 10⁵ | 10¹ | >24h | 4.462s | **>19000x** |
| 10⁶ | 10¹ | N/A | 53.64s | — |

关键观察：当 GCD 大而余因子小（#G ≫ #Ā）时，同时插值策略获得巨大加速，因为只需 O(#Ā) 个图像而非 O(#G)。

**与 FLINT BMA/Zippel2 的关系**：FLINT 的 `_try_bma` 实现了类似的 BMA 支撑检测思路，但没有 Kronecker 替换和同时插值余因子的优化。Hu/Monagan 是目前已知最完整的稀疏 GCD 实现。

---

## 4. 推荐优化方向与优先级

### P0：修复单变量素数起始位置 + cont 重复计算（即时修复）

**目标**：消除单变量 ZZ GCD 的冗余 CRT 迭代

**方案 a — 大素数起始**：将 `polynomial_gcd.cc:43-44` 的素数起始策略改为参考 FLINT（从大素数开始）：

```cpp
// 当前（.cc:43-44）
std::uint32_t p_index=0;  // prime = 2，每素数贡献 ~1 bit
// 方案 1: 与多变量版一致
std::uint32_t tmp_x = std::max(get_deg(F), get_deg(G));
if (tmp_x < 2) tmp_x = 2;
std::uint32_t p_index = tmp_x / std::log(tmp_x);  // prime ≈ 157 for deg200
// 方案 2: 参考 FLINT，从尽可能大的素数开始（需改用 64-bit 素数生成）
// p = nextprime(2^63)  // 每素数贡献 ~63 bits，4 次即可
```

方案 1 简单可行（复用多变量逻辑），将迭代从 ~200 次减到 ~27 次。
方案 2 与 FLINT 一致（~4 次迭代），但需要 64-bit 素数生成器（`boost::math::prime` 仅覆盖到第 10000 个素数 ≈ 104729）。

**方案 b — 修复 cont 重复计算**（`polynomial_gcd.hh:299`）：

```cpp
// 当前（重复计算 content）
polynomial_<ZZ,...> cont_gcd = polynomial_GCD(cont(F), cont(G));
// 改为（使用已有变量）
polynomial_<ZZ,...> cont_gcd = polynomial_GCD(F_cont, G_cont);
```

- **预期收益**：**3-8x**（单变量，大素数方案）；1.2x（多变量 cont 修复）
- **难度**：极低（几行代码）
- **风险**：无（仅优化，不改算法正确性）
- **依据**：FLINT `gcd_modular.c` line 117-118，NTL `ZZX1.cpp` line 2879

### P1：实现 GCDHEU（启发式 GCD）

**目标**：为小/中系数多项式提供快速路径

**方案（推荐 FLINT bit-packing 方案，参考 §2.5）**：

```
在 polynomial_GCD 入口处，若 max_bits(A) + max_bits(B) < 128:
1. pack_bits = max(min(bits₁, bits₂) + 6, max(bits₁, bits₂) + 1)
2. 将 A, B 的系数 bit-pack 为大整数 a_int, b_int（ξ = 2^pack_bits）
3. g_int = mpz_gcd(a_int, b_int)  — GMP 大整数 GCD
4. 从 g_int 以 pack_bits 为单位解包恢复多项式 G
5. 除以 content(G)
6. 验证: mpz_divisible(a_int, g_int) && mpz_divisible(b_int, g_int)
   可选: 解包商并检查 bits_G + bits_Q + log(len) < pack_bits
7. 通过 → 返回 G·d（乘以 content GCD）；失败 → 回退到模算法
```

替代方案（经典 GCDHEU）：取 ξ = 2·max|coeff|·(deg+1)，计算 gcd(A(ξ), B(ξ))，base-ξ 展开。数学等价但打包/解包效率低于 bit-packing。

适用范围：`max_bits(A) + max_bits(B) < 128`（参考 FLINT `gcd.c:61` 阈值）。

- **预期收益**：2-5x（小系数，如 benchmark 中 coeff ∈ [-20, 20] 的用例）
- **难度**：中（需要 bit-pack/unpack + mpz_gcd + 验证，约 100-150 行代码）
- **风险**：低（GCDHEU 失败时回退到已有模算法）
- **依据**：FLINT `gcd_heuristic.c`（bit-packing）、arxiv:cs/0206032v1（理论界）

### P2：CRT 循环优化 + 大素数

**目标**：减少 CRT 常数开销

**方案**：
1. 预计算 `inv_times_prime = tmp_inv.number() * Pout_prime` 和 `half_prime = Pout_prime / 2`
2. 使用更大的素数起点（参考 §2.2 分析）
3. 将 `pair_vec_div` 的堆分配改为栈缓冲或预分配工作区

- **预期收益**：1.3-1.5x
- **难度**：低
- **风险**：无

### P3：密集单变量 Zp GCD 路径

**目标**：绕过稀疏数据结构的开销

**方案**：当 Zp 多项式为密集（非零项 / (deg+1) > 0.5）时，将系数提取到密集 `uint64_t` 数组，用简单指针算术做欧几里德除法，再转换回稀疏表示。

- **预期收益**：1.5-2x（消除 `pair_vec_div` 堆操作 + 单项式比较开销）
- **难度**：中（需要密集/稀疏转换函数，~150 行代码）
- **风险**：低
- **备注**：这是实现 HGCD 的基础设施前提

### P4：实现 HGCD（半 GCD）— 长期

**目标**：将单变量 Zp GCD 从 O(n^2) 降到 O(M(n) log n)

**方案**：实现递归 Half-GCD 算法（参考 GCL §11.1 或 von zur Gathen-Gerhard §11）。

**源码参考**：FLINT `src/nmod_poly/hgcd.c` 委托给 `_gr_poly_hgcd`（generic ring 实现）。HGCD 基础情形阈值 `hgcd_iter_recursive_cutoff` 在 68-244 范围。GCD 切换到 HGCD 的阈值：
- `nmod_poly_gcd_hgcd_cutoff_tab`：470-2199（按 bit-size）
- 对 64-bit 素数：**1725**

**重要修正**：FLINT 在 deg < 1725（64-bit 素数）时使用欧几里德，**不使用 HGCD**。因此 HGCD 对当前 benchmark（deg ≤ 200）**没有直接收益**。HGCD 的价值在于 deg ≥ 2000 的大多项式。

依赖 P3（需要密集 Zp 数组基础设施）。

- **预期收益**：deg ≤ 200 时**无收益**；deg ≥ 2000 时 **3-10x**
- **难度**：高（HGCD 实现复杂，递归矩阵运算，约 400-600 行代码）
- **风险**：中（算法正确性需仔细验证）
- **优先级调整**：根据源码分析，从"中期"降级为"长期"，P0（大素数）+ P1（GCDHEU）足以解决 deg ≤ 200 的差距

### P5：Zippel 稀疏插值（多变量）

**目标**：将多变量 GCD 从密集插值升级为稀疏插值

**方案**：实现 Zippel 1979 Algorithm S（稀疏模 GCD）。以下是论文 §3 的精确描述：

```
Algorithm S（稀疏模 GCD，Zippel 1979 §3）:
输入: F₁, F₂ ∈ Z[x₁,...,xᵥ]
输出: G = gcd(F₁, F₂)

1. 选择主变量 x₁（度数最大者）
2. 选素数 p，随机赋值 b₂,...,bᵥ ∈ Zp
3. 计算骨架(skeleton)：
   g₁(x₁) = gcd(F₁(x₁,b₂,...,bᵥ), F₂(x₁,b₂,...,bᵥ)) ∈ Zp[x₁]
   记录非零项的指数集 S = {e : [x₁^e]g₁ ≠ 0}
4. 逐变量恢复（k = 2,...,v）:
   对每个 e ∈ S:
     已知 G 中 x₁^e 的系数是 Cₑ(x₂,...,xₖ₋₁,xₖ,...,xᵥ) ∈ Z[x₂,...,xₖ]
     从前一步已确定 Cₑ 对 x₂,...,xₖ₋₁ 的结构
     需要确定 Cₑ 中 xₖ 的幂次
   选择 tₖ = (Cₑ 中 xₖ 的已知项数) 个新求值点 αₖ,₁,...,αₖ,tₖ ∈ Zp
   对每个求值点 αₖ,ⱼ:
     特化 xₖ = αₖ,ⱼ（保留 x₁,...,xₖ₋₁ 和随机的 bₖ₊₁,...,bᵥ）
     递归计算 GCD 图像
   用 Vandermonde 系统求解 Cₑ 的 xₖ 系数
5. CRT 合并多个素数的结果 → 对称模约简 → 试除验证
```

**关键假设**（Zippel 1979, Theorem 1）：若某幂次 x₁^e 在第一次单变量 GCD 中系数为零（即 e ∉ S），则以概率 ≥ 1 - v²td/B（B 是 Zp 大小）该幂次在完整多项式中恒为零。这使得插值点数从 O(d+1) 降至 O(t)。

**错误概率界**（Zippel 1979, Theorem 1）：坏求值点数 N_v(B) ≤ B^v - (B-D)^v，其中 D = max(deg)。对大素数 p ≫ D，概率 ≈ vD/p，可通过重试消除。

长期可考虑 Hu/Monagan (BMA + Kronecker) 作为终极方案。

- **预期收益**：2-10x（稀疏多变量）；密集情况无改善
- **难度**：高（归一化问题处理复杂，约 500-800 行代码）
- **风险**：中（需要大量测试验证稀疏情况的正确性）

---

## 5. 结论

### 差距汇总

| 根因 | 影响范围 | 差距贡献 | 修复方案 | 难度 | 优先级 | 源码依据 |
|------|---------|---------|---------|------|--------|---------|
| 单变量从 prime=2 开始 | 单变量 | **3-8x** | 大素数起始 | 极低 | **P0** | FLINT `gcd_modular.c:117` |
| `cont(F)` 重复计算 | 多变量 | 1.2-1.5x | 用已有变量 | 极低 | **P0** | — |
| 无 GCDHEU | 全局（小系数） | 2-5x | bit-packing + mpn_gcd | 中 | **P1** | FLINT `gcd_heuristic.c` |
| CRT 循环冗余计算 | 全局 | 1.2-1.3x | 预计算不变量 | 低 | **P2** | — |
| 堆分配开销 | 全局 | 1.1-1.3x | 栈缓冲/预分配 | 低 | **P2** | — |
| 稀疏数据结构做密集除法 | 单变量 Zp | 1.5-2x | 密集 Zp 路径 | 中 | **P3** | FLINT `nmod_poly` 密集数组 |
| 无 HGCD | 单变量（**deg ≥ 2000**） | 3-10x | 实现 HGCD | 高 | **P4** | FLINT 阈值 1725 |
| 无 Zippel | 多变量（稀疏） | 2-10x | 实现 Zippel | 高 | **P5** | FLINT `gcd_algo.c` |
| Zp 32B/元素 | 全局 | 1.2-1.5x | Zp64 分离上下文 | 中 | (已有设计) | FLINT `nmod_t` 8B |
| 无幂次预计算 | 多变量求值 | 1.1-1.2x | LUT 幂次表 | 低 | 附带 | FLINT `fmpz_mpoly_evals` LUT |
| 无早期退出（多变量） | 多变量 | 1.5-3x | 单项式/整除检查 | 中 | 附带P5 | FLINT 8 级早期退出链 |

### 复合效应估算（基于源码分析修正）

**关键修正**：此前认为 HGCD 贡献 3-5x（§2.1 旧版），现已确认 FLINT 在 deg200 也用欧几里德。**素数大小是首要差距来源**。

**Brown 1971 §5.8 的理论支撑**：Brown 证明模算法 M 的复杂度 `M(v,l,d) ≲ l²(d+1)^v + l(d+1)^(v+1)` 严格被经典算法 C 的第一次伪除法支配（对 v ≥ 2）。差距因子 `C_/M+ = (d+1)^(v-2)`，即变量数越多，模算法优势越大。对于 CLPoly 的 benchmark（2 变量），`(d+1)^0 = 1`，模算法仅在常数上优于经典算法；但对 3+ 变量，优势呈指数增长。这说明 CLPoly 目前的模 GCD 框架是正确的算法路线，性能差距主要来自实现细节（素数大小、GCDHEU、早期退出等）。

**单变量 deg200**（当前 5.072ms vs NTL 0.532ms = 9.53x）：
- P0 大素数（5x: 200→4 次迭代） × P1 GCDHEU 快速路径（跳过模算法 ≈ 2x） × P2 CRT 优化（1.3x）
- P0+P1+P2 ≈ **~10x** → 可将 5.072ms 降至 ~0.5ms（接近 NTL 的 0.532ms）
- P3（密集 Zp 路径）可进一步降低每次 Zp GCD 的常数

**多变量 deg15**（当前 1.477ms vs FLINT 0.172ms = 8.58x）：
- P0/cont 修复（1.3x） × P1 GCDHEU（2x） × P2（1.2x） × P5/早期退出（3x）
- P0+P1+P2+P5 ≈ **~9x** → 可将 1.477ms 降至 ~0.16ms（接近 FLINT 的 0.172ms）

### 行动计划

1. **立即（P0）**：大素数起始 + cont 重复计算 — **影响最大的单一修复**，几行代码，0 风险
2. **短期（P1+P2）**：实现 GCDHEU（参考 FLINT bit-packing 方案） + CRT 优化 — 预期最大收益/工作量比
3. **中期（P3）**：密集 Zp 路径 — 消除 pair_vec_div 常数开销
4. **长期（P4+P5）**：HGCD（仅 deg ≥ 2000 有价值）+ Zippel 稀疏插值 — 大规模/多变量竞争力
