# CLPoly 多项式因式分解设计方案

> 目标：为 CLPoly 添加完整的多项式因式分解支持，覆盖 Z[x] 单变量和 Z[x₁,...,xₙ] 多变量情况。
> 设计参考了 Maple、Mathematica、FLINT、NTL、Singular/Factory 五个主流代数系统。

---

## 1. 背景与动机

### 1.1 当前能力

CLPoly 已具备因式分解所需的大部分基础设施：

| 组件 | 状态 | 位置 |
|---|---|---|
| 无平方分解 `squarefreefactorize` | ✅ | `polynomial_gcd.hh` |
| 多变量 GCD（模 + Hensel） | ✅ | `polynomial_gcd.hh` |
| 有限域 `Zp` 类 | ✅ | `number.hh` |
| `polynomial_mod(f, p)` 模约化 | ✅ | `polynomial.hh` |
| 单变量 Zₚ 上的 Euclidean GCD | ✅ | `polynomial_gcd.hh` |
| 结式 / 子结式 / 判别式 | ✅ | `resultant.hh` |
| 伪余式 / 伪商 | ✅ | `polynomial.hh` |

### 1.2 缺失部分

完整因式分解还需要四个核心模块：

1. **有限域 Zₚ 上的不可约分解** — Cantor-Zassenhaus 算法
2. **Hensel 提升** — 从 mod p 因子恢复到 mod p^k
3. **因子重组** — 从模因子重建 Z 上真因子
4. **集成与多变量扩展** — 端到端流程 + Wang 算法

---

## 2. 主流代数库设计对比

### 2.1 通用管线

所有主流系统的因式分解都遵循同一宏观流程：

```
输入 f ∈ Z[x₁, ..., xₙ]
  │
  ▼
(1) 内容提取: f = cont(f) · pp(f)
  │
  ▼
(2) 无平方分解: pp(f) = ∏ gᵢ^i
  │
  ▼
(3) 对每个无平方因子 gᵢ:
  │  ├─ [多变量] 求值 x₂=α₂,...,xₙ=αₙ 得到单变量像
  │  ├─ 选素数 p，在 Zₚ 上因式分解
  │  ├─ Hensel 提升: mod p → mod p^k
  │  ├─ 因子重组: 组合提升因子为真因子
  │  └─ [多变量] 多变量 Hensel 提升恢复各变量
  │
  ▼
输出: cont(f) · ∏ (不可约因子ⱼ)^(指数ⱼ)
```

### 2.2 各系统设计选择

| | Maple | Mathematica | FLINT | NTL | Singular |
|---|---|---|---|---|---|
| **Zₚ 分解** | Berlekamp/CZ | 可能 CZ | CZ + Berlekamp(p=2) + Kaltofen-Shoup(大度数) | CZ (NewDDF+EDF) | Berlekamp, CZ |
| **Hensel 提升** | 二次 | 可能二次 | 二次 | 二次 (MultiLift) | 线性+二次 |
| **单变量重组** | van Hoeij | 可能 van Hoeij | Zassenhaus + van Hoeij 备选 | van Hoeij(默认) | Zassenhaus |
| **多变量策略** | MTSHL | 可能 Wang/EEZ | fmpz_mpoly_factor | 不支持 | Wang + 稀疏 Hensel |
| **因子表示** | 内部 DAG | `{cont,{fᵢ,eᵢ}}` | C 结构体 | `vec_pair<ZZX,long>` | CFFList |

### 2.3 关键设计决策分析

**Cantor-Zassenhaus vs Berlekamp：**
- Berlekamp 基于矩阵核空间计算，时间 O(n² p + n³)，适合小 p
- Cantor-Zassenhaus 基于概率 GCD，时间 O(n² log n log p)，适合大 p
- **选择 CZ**：CLPoly 使用的素数通常较大，CZ 更通用且实现更简洁

**Zassenhaus vs van Hoeij 重组：**
- Zassenhaus：枚举 2^r 个子集，最坏指数级，但实现简单
- van Hoeij：LLL 格约化，多项式时间，但需要 LLL 实现
- **选择分阶段**：先实现 Zassenhaus（覆盖绝大多数实际情况），后加 van Hoeij

**多变量策略：**
- Wang 算法：经典、确定性，被广泛实现
- MTSHL（Maple 2019）：稀疏情况更优，但更复杂
- **选择 Wang**：成熟可靠，与 CLPoly 已有的 Hensel 基础设施契合

---

## 3. 模块设计

### 3.1 模块总览

```
                      ┌─────────────────┐
                      │  M5: 多变量分解   │  factor(polynomial_ZZ)  [多变量]
                      │  factor_multivar │
                      └────────┬────────┘
                               │
                ┌──────────────┼──────────────┐
                ▼              ▼              ▼
       ┌──────────────┐  ┌─────────┐  ┌──────────────┐
       │ M4: Z[x] 分解 │  │ 求值/代入│  │ 多变量 Hensel │
       │ factor_ZZ     │  │ (已有)   │  │  提升 (M5b)  │
       └──┬────┬────┬──┘  └─────────┘  └──────────────┘
          │    │    │
          ▼    ▼    ▼
      ┌─────┐┌──┐┌───┐
      │ M1  ││M2││M3 │
      │Zₚ  ││提││重  │
      │分解 ││升││组  │
      └──┬──┘└──┘└───┘
         │
      ┌──┼──┐
      ▼  ▼  ▼
    DDF EDF SqFree
   (M1a)(M1b)(已有)
```

### 3.2 文件组织

沿用 CLPoly 现有的头文件组织方式（header-only），新增：

```
clpoly/
  polynomial_factorize.hh    — M1 + M2 + M3 + M4 集成
```

单变量算法集中在一个头文件内，通过命名空间或 `__` 前缀区分内部函数（与
`polynomial_gcd.hh` 中 `__polynomial_GCD` 的风格一致）。多变量扩展（M5）
视复杂度决定是否单独成文件。

---

## 4. 辅助函数层：模运算基础设施

在进入核心模块之前，需要若干辅助函数。它们服务于多个模块，统一放在
`polynomial_factorize.hh` 的前段。

### 4.1 单变量 Zₚ 多项式运算

已有基础设施中，`upolynomial_<Zp>` 的算术 (+, -, *, /)、`polynomial_mod`、
`__polynomial_GCD`（Euclidean GCD）都已可用。还需补充以下操作：

#### 4.1.1 `__upoly_make_monic` — 首一化

```cpp
// 将 f 除以其首项系数，使其首一
// 前置: f 非零
// 后置: f 首一 (lc = 1), 返回原 lc
// 用途: M1 (DDF/EDF 要求首一输入)
inline Zp __upoly_make_monic(upolynomial_<Zp>& f)
{
    Zp lc = f.front().second;
    Zp lc_inv = lc.inv();
    for (auto& term : f)
        term.second *= lc_inv;
    return lc;
}
```

#### 4.1.2 `__upoly_mod` — 多项式取模 (f mod g)

```cpp
// 计算 f mod g in Zₚ[x]（即 f 除以 g 的余式）
// 前置: g 非零
// 后置: 返回 r 满足 deg(r) < deg(g), 存在 q 使 f = q·g + r
// 用途: M1 (DDF 中 h^p mod f* 需要模约化), M2 (Hensel 中 Bézout 系数维护)
// 实现: 复用已有 pair_vec_div（它同时产生商和余式）
inline upolynomial_<Zp> __upoly_mod(
    const upolynomial_<Zp>& f,
    const upolynomial_<Zp>& g)
{
    upolynomial_<Zp> q, r;
    pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
    return r;
}
```

#### 4.1.3 `__upoly_divmod` — 同时返回商和余式

```cpp
// 计算 f = q·g + r in Zₚ[x]
// 前置: g 非零
// 后置: f = q·g + r, deg(r) < deg(g)
// 用途: M2 (Hensel 提升中 divmod(s·e, h))
inline void __upoly_divmod(
    upolynomial_<Zp>& q,
    upolynomial_<Zp>& r,
    const upolynomial_<Zp>& f,
    const upolynomial_<Zp>& g)
{
    pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
}
```

#### 4.1.4 `__upoly_divmod_mod` — ZZ 系数的模除法

```cpp
// 在 Z_m[x] 中计算 f = q·g + r，即所有系数运算都 mod m
// 前置: g 非零, lc(g) 与 m 互素 (即 lc(g) mod m 可逆)
//       m > 0
// 后置: f ≡ q·g + r (mod m), deg(r) < deg(g)
//       q, r 的系数 ∈ [0, m)
// 用途: M2 (Hensel 提升 __hensel_step 中 divmod(s·e, h) in Z_m[x])
//       Hensel 提升到高精度时 m = p^{2^j} 远超 uint32_t，
//       无法使用 Zp 类，因此需要直接在 ZZ 系数上做模运算
// 实现: 与 pair_vec_div 类似的长除法，但每步系数用 ZZ 模逆
//       lc_inv ← lc(g) 的模逆 (通过扩展 GCD 或 ZZ::invert 计算)
//       每次消去首项后对所有系数 mod m
inline void __upoly_divmod_mod(
    upolynomial_<ZZ>& q,
    upolynomial_<ZZ>& r,
    const upolynomial_<ZZ>& f,
    const upolynomial_<ZZ>& g,
    const ZZ& m)
```

#### 4.1.5 `__upoly_gcd_Zp` — Zₚ 上简单 GCD 包装

```cpp
// 简化的 Zₚ GCD 接口: 计算 gcd(a, b) in Zₚ[x]，返回首一结果
// 前置: a 或 b 非零
// 后置: 返回 g = gcd(a, b)，首一化
// 说明: 内部调用已有 __polynomial_GCD(Pout, G, F, Lc_gcd, deg)，
//       自动设置 Lc_gcd = Zp(1, p) 和 deg = min(deg(a), deg(b))
// 用途: M1 (DDF/EDF 中频繁需要 gcd 计算)
inline upolynomial_<Zp> __upoly_gcd_Zp(
    const upolynomial_<Zp>& a,
    const upolynomial_<Zp>& b)
{
    if (a.empty()) return b;
    if (b.empty()) return a;
    upolynomial_<Zp> g;
    Zp lc_gcd(1, a.front().second.prime());
    int64_t deg = std::min(get_deg(a), get_deg(b));
    __polynomial_GCD(g, a, b, lc_gcd, deg);
    __upoly_make_monic(g);
    return g;
}
```

#### 4.1.6 `__upoly_gcd_extended` — 扩展 GCD (Bézout 系数)

```cpp
// 扩展欧几里得算法: 求 s, t 使得 s·a + t·b = gcd(a, b) in Zₚ[x]
// 前置: a, b 非零, gcd(a,b) = 1 (即 a, b 互素)
// 后置: s·a + t·b = 1, deg(s) < deg(b), deg(t) < deg(a)
// 用途: M2 (Hensel 提升需要初始 Bézout 系数)
inline void __upoly_gcd_extended(
    upolynomial_<Zp>& s,
    upolynomial_<Zp>& t,
    const upolynomial_<Zp>& a,
    const upolynomial_<Zp>& b)
```

**算法：** 标准扩展 Euclidean，在 `__polynomial_GCD` 基础上增加系数追踪：

```
r₀ ← a,    r₁ ← b
s₀ ← 1,    s₁ ← 0
t₀ ← 0,    t₁ ← 1
while r₁ ≠ 0:
    q, r₂ ← divmod(r₀, r₁)
    s₂ ← s₀ - q·s₁
    t₂ ← t₀ - q·t₁
    (r₀,s₀,t₀) ← (r₁,s₁,t₁)
    (r₁,s₁,t₁) ← (r₂,s₂,t₂)
// 此时 r₀ = gcd = 1 (up to scalar)
// 归一化使 r₀ = 1
c ← r₀.front().second.inv()
s ← c · s₀
t ← c · t₀
```

#### 4.1.7 `__upoly_powmod` — 模幂运算

```cpp
// 计算 base^exp mod modpoly in Zₚ[x]
// 前置: modpoly 非零
// 后置: 返回 base^exp mod modpoly, 度数 < deg(modpoly)
// 用途: M1 (DDF 中计算 x^{p^d} mod f; EDF 中计算 r^{(p^d-1)/2} mod f)
// 注意: exp 可能非常大 (p^d)，必须用 ZZ 类型
inline upolynomial_<Zp> __upoly_powmod(
    const upolynomial_<Zp>& base,
    const ZZ& exp,
    const upolynomial_<Zp>& modpoly)
```

**算法：** 平方-乘法（binary exponentiation），每步做 mod modpoly：

```
result ← 1
b ← base mod modpoly
e ← exp
while e > 0:
    if e 是奇数:
        result ← (result · b) mod modpoly
    b ← (b · b) mod modpoly
    e ← e >> 1
return result
```

**复杂度：** O(log(exp) · M(n))，其中 M(n) = O(n²) 为度数 n 多项式乘法。

#### 4.1.8 `__upoly_random` — 生成随机 Zₚ 多项式

```cpp
// 生成度数 < max_deg 的随机 Zₚ[x] 多项式
// 前置: p > 0, max_deg > 0
// 后置: 返回随机多项式, 0 ≤ deg < max_deg, 系数在 [0, p) 中均匀分布
// 用途: M1 (EDF/Cantor-Zassenhaus 中的随机元素选取)
inline upolynomial_<Zp> __upoly_random(
    int64_t max_deg,
    uint32_t p,
    std::mt19937& rng)
```

**说明：** 使用 `std::mt19937`（与 `polynomial_gcd.hh` 中已有的 `std::random`
用法一致），`std::uniform_int_distribution<uint64_t>(0, p-1)` 生成系数。

### 4.2 Z 上的对称模约化

#### 4.2.1 `__symmetric_mod` — 标量对称模

```cpp
// 将 a 约化到 (-m/2, m/2] 范围
// 前置: m > 0
// 后置: 返回 r ≡ a (mod m), -m/2 < r ≤ m/2
// 用途: M3 (Zassenhaus 重组中将提升系数映射回 Z)
inline ZZ __symmetric_mod(const ZZ& a, const ZZ& m)
{
    ZZ r;
    ZZ::fdiv_r(r, a, m);    // r ∈ [0, m)   (必须用 fdiv_r, 不能用 operator%
                             //               后者对负数 a 返回负余数)
    if (r > m / 2)
        r -= m;
    return r;
}
```

#### 4.2.2 `__upoly_symmetric_mod` — 多项式对称模

```cpp
// 对多项式每个系数做对称模约化
// 前置: m > 0
// 后置: 返回的多项式每个系数 ∈ (-m/2, m/2], 且 ≡ 原系数 (mod m)
// 用途: M3 (Zassenhaus 重组)
inline upolynomial_<ZZ> __upoly_symmetric_mod(
    const upolynomial_<ZZ>& f,
    const ZZ& m)
```

### 4.3 范数计算

```cpp
// L1 范数 (系数绝对值之和)
// 用途: M3 (Zassenhaus 范数剪枝: ‖g‖₁ · ‖h‖₁ ≤ ‖f‖₁)
inline ZZ __upoly_norm_l1(const upolynomial_<ZZ>& f)

// L2 范数的平方 (系数平方和), 避免开方
// 用途: M2 (Mignotte 界计算)
inline ZZ __upoly_norm_l2_sq(const upolynomial_<ZZ>& f)

// L∞ 范数 (系数绝对值最大值)
// 用途: M2 (简化 Mignotte 界)
inline ZZ __upoly_norm_linf(const upolynomial_<ZZ>& f)
```

### 4.4 Mignotte 界计算

```cpp
// 计算 f 的因子系数绝对值上界 (Mignotte bound)
// 前置: f 非零
// 后置: 返回 B 使得 f 的任意因子 g 的系数绝对值 ≤ B
//       具体: B = C(n, ⌊n/2⌋) · ‖f‖₂
// 用途: M2 (确定 Hensel 提升精度 k: 需 p^k > 2B)
//       M3 (范数剪枝参考值)
inline ZZ __mignotte_bound(const upolynomial_<ZZ>& f)
```

**实现：** 使用简化界 B = C(n, ⌊n/2⌋) · ‖f‖₂，其中 C(n,k) 用 ZZ 精确计算。
完整的 Mignotte 界为 C(n-1, ⌊(n-1)/2⌋) · ‖f‖₂ + C(n-1, ⌊(n-1)/2⌋-1) · |lc(f)|，
此处使用的简化版可能略大，但不影响正确性（只是多提升若干步）。

### 4.5 本原化

```cpp
// 提取内容并本原化 (对 upolynomial_<ZZ>)
// 前置: f 非零
// 后置: 返回 (content, primitive_part), content · primitive_part = f
//       primitive_part 首项系数 > 0
// 备注: cont(upolynomial_<ZZ>) 已在 upolynomial.hh 中实现
// 用途: M3, M4
inline std::pair<ZZ, upolynomial_<ZZ>> __upoly_primitive(upolynomial_<ZZ> f)
{
    ZZ c = cont(f);
    if (f.front().second < 0) c = -c;  // 确保 lc > 0
    for (auto& term : f)
        term.second /= c;
    return {c, std::move(f)};
}
```

---

## 5. 模块 M1：有限域 Zₚ 上的单变量分解

### 5.1 职责

给定 f ∈ Zₚ[x]，返回其完整的不可约因子分解。

### 5.2 函数总览

```
__factor_Zp                    M1 顶层入口
  ├── __upoly_make_monic       首一化 (§4.1.1)
  ├── __squarefree_Zp          无平方分解
  │     ├── __upoly_gcd_Zp     GCD 包装 (§4.1.5)
  │     ├── derivative         求导 (已有)
  │     └── __extract_pth_root p 次根提取
  ├── __ddf_Zp                 按度数分组
  │     ├── __upoly_powmod     模幂 (§4.1.7)
  │     └── __upoly_gcd_Zp     GCD 包装 (§4.1.5)
  └── __edf_Zp                 等度分裂 (Cantor-Zassenhaus)
        ├── __upoly_random     随机多项式 (§4.1.8)
        ├── __upoly_powmod     模幂 (§4.1.7)
        └── __upoly_gcd_Zp     GCD 包装 (§4.1.5)
```

### 5.3 `__squarefree_Zp` — Zₚ 上的无平方分解

```cpp
// Zₚ 上无平方分解
// 前置: f 首一非零
// 后置: ∏ fᵢ^eᵢ = f, 每个 fᵢ 首一无平方且两两互素, eᵢ > 0
// 注意: 与特征 0 (squarefreefactorize) 不同，特征 p 下 f' 可能为零
//       需要处理 f(x) = g(x^p) 的情况
std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
__squarefree_Zp(const upolynomial_<Zp>& f)
```

**算法（详细）：**

```
输入: 首一 f ∈ Zₚ[x]
输出: [(s₁,e₁), (s₂,e₂), ...]

1.  f' ← derivative(f)

2.  if f' = 0:
        // f = g(x^p), 所有指数被 p 整除
        g ← __extract_pth_root(f)
        对 __squarefree_Zp(g) 的每个 (sᵢ, eᵢ):
            输出 (sᵢ, eᵢ · p)
        return

3.  c ← __upoly_gcd_Zp(f, f')
    w ← f / c
    i ← 1
    result ← []

4.  while w ≠ 1:
        y ← __upoly_gcd_Zp(w, c)
        z ← w / y          // z 是 exactly multiplicity i 的因子之积
        if z ≠ 1:
            __upoly_make_monic(z)
            result.push_back({z, i})
        w ← y
        c ← c / y
        i ← i + 1

5.  if c ≠ 1:
        // c = c₀(x^p)
        g ← __extract_pth_root(c)
        对 __squarefree_Zp(g) 的每个 (sⱼ, eⱼ):
            result.push_back({sⱼ, eⱼ · p})

6.  return result
```

#### 5.3.1 `__extract_pth_root` — p 次根提取

```cpp
// 提取 f(x) = g(x^p) 中的 g
// 前置: f 的所有项的度数均为 p 的倍数
// 后置: 返回 g 使得 g(x^p) = f(x)
// 实现: 将每个 (x^{kp}, aₖ) 映射为 (x^k, aₖ)
//       Zₚ 上 Frobenius 逆是恒等映射，系数不变
inline upolynomial_<Zp> __extract_pth_root(
    const upolynomial_<Zp>& f)
{
    uint32_t p = f.front().second.prime();
    upolynomial_<Zp> g;
    g.reserve(f.size());
    for (auto& term : f)
    {
        assert(term.first.deg() % p == 0);
        g.push_back({umonomial(term.first.deg() / p), term.second});
    }
    return g;
}
```

### 5.4 `__ddf_Zp` — 按度数分组 (Distinct-Degree Factorization)

```cpp
// DDF: 将 f 的不可约因子按度数分组
// 前置: f ∈ Zₚ[x], 首一, 无平方
// 后置: 返回 [(g₁,1), (g₂,2), ...] 其中 gd = ∏{不可约 h: deg(h)=d, h|f} h
//       ∏ gd = f
//       每个 gd 首一
//       跳过 gd = 1 的项 (即列表中不含平凡因子)
std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
__ddf_Zp(const upolynomial_<Zp>& f)
```

**算法（详细）：**

```
输入: 首一无平方 f ∈ Zₚ[x], p = f 的系数素数
输出: [(g₁,1), (g₂,2), ...]

1.  h ← x (即 upolynomial_<Zp> 为 {(1, Zp(1,p))})
    f* ← f
    result ← []

2.  for d = 1, 2, ... :
        if deg(f*) < 2·d:
            break                      // deg(f*) < 2d 意味着 f* 本身不可约
        h ← __upoly_powmod(h, ZZ(p), f*)   // h = x^{p^d} mod f*
        // 构造 (h - x) mod f*
        h_minus_x ← h
        h_minus_x 的 x^1 项系数减 1     // h - x
        gd ← __upoly_gcd_Zp(h_minus_x, f*)

        if gd ≠ 1 (即 deg(gd) > 0):
            __upoly_make_monic(gd)
            result.push_back({gd, d})
            f* ← f* / gd
            h ← h mod f*              // 缩减 h 的度数

3.  if deg(f*) > 0:
        __upoly_make_monic(f*)
        result.push_back({f*, (uint64_t)get_deg(f*)})

4.  return result
```

**关键优化：**
- 步骤 2 中 h 是累积的：每步 h ← h^p mod f*，即 h 依次等于
  x^p, x^{p²}, x^{p³}, ...。无需每次从头计算。
- 当找到因子 gd 后，f* 缩小，h 也做 mod f* 以保持度数小。
- 循环在 deg(f*) < 2d 时提前退出：此时 f* 本身必然不可约。

### 5.5 `__edf_Zp` — 等度分裂 (Equal-Degree Factorization, Cantor-Zassenhaus)

```cpp
// EDF: 将 f 分裂为全部度数为 d 的不可约因子
// 前置: f ∈ Zₚ[x], 首一, 无平方, f 的所有不可约因子度数均为 d
//       deg(f) = d·k (k 为因子个数)
// 后置: 返回恰好 k = deg(f)/d 个首一不可约因子, ∏ = f
//       因子顺序不确定（随机算法）
// 注意: 概率算法，期望迭代次数 O(k) (每次分裂成功概率 ≥ 1 - 1/p^d ≈ 1/2)
void __edf_Zp(
    std::vector<upolynomial_<Zp>>& result,
    const upolynomial_<Zp>& f,
    uint64_t d,
    std::mt19937& rng)
```

**算法（详细）：**

```
输入: f, d, rng
输出: result (追加到 result 末尾)

1.  if deg(f) = d:
        result.push_back(f)
        return

2.  p ← f 的系数素数
    n ← deg(f)

3.  repeat:
        r ← __upoly_random(n, p, rng)    // 度数 < n 的随机多项式
        if deg(r) = 0:
            continue

        if p = 2:
            // 特征 2: 计算 trace map T(r) = r + r² + r⁴ + ... + r^{2^{d-1}} mod f
            g ← r
            for i = 1 to d-1:
                g ← __upoly_mod(g * g + r, f)   // g = g² + r mod f
            g ← __upoly_gcd_Zp(g, f)
        else:
            // 奇特征: 计算 g = gcd(r^{(p^d-1)/2} - 1, f)
            exp ← (ZZ::pow(p, d) - 1) / 2
            g_pow ← __upoly_powmod(r, exp, f)
            // g_pow - 1
            g_pow 的常数项减 1
            g ← __upoly_gcd_Zp(g_pow, f)

        if 0 < deg(g) < deg(f):
            // 成功分裂
            h ← f / g
            __upoly_make_monic(g)
            __upoly_make_monic(h)
            __edf_Zp(result, g, d, rng)
            __edf_Zp(result, h, d, rng)
            return
    // repeat 直到成功（期望 O(1) 次）
```

### 5.6 `__factor_Zp` — Zₚ 上完整分解（顶层入口）

```cpp
// Zₚ 上完整不可约分解
// 前置: f 非零, f 的系数类型为 Zp (所有系数共享同一素数 p)
// 后置: 返回 (lc, factors) 满足:
//       lc · ∏ factors[i].first ^ factors[i].second = f
//       每个 factors[i].first 首一不可约
//       factors[i].second > 0
//       factors 按 (degree, then lex) 排列
std::pair<Zp, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
__factor_Zp(upolynomial_<Zp> f)
```

**算法：**

```
1.  if deg(f) = 0:
        return (f.front().second, {})        // 常数

2.  lc ← __upoly_make_monic(f)              // 提取并首一化

3.  sqf ← __squarefree_Zp(f)                // 无平方分解

4.  result ← []
    rng ← std::mt19937(seed)                // 固定或随机 seed

5.  for (sⱼ, eⱼ) in sqf:
        ddf ← __ddf_Zp(sⱼ)                  // 按度数分组
        for (gk, dk) in ddf:
            __edf_Zp(factors_k, gk, dk, rng) // 等度分裂
            for hᵢ in factors_k:
                result.push_back({hᵢ, eⱼ})

6.  排序 result (按 degree 升序，同度数按字典序)

7.  return (lc, result)
```

**参数传递说明：** f 按值传入（与 `polynomial_GCD` 的风格一致），内部会修改。

**优化提示：** 当 `__factor_Zp` 被 `__select_prime` 调用时，`__select_prime` 已
确认 `f mod p` 无平方，此时步骤 3 的 `__squarefree_Zp` 会直接返回 `{(f, 1)}`，
开销很小。如需进一步优化，可增加 `bool already_squarefree = false` 参数跳过此步。

---

## 6. 模块 M2：Hensel 提升

### 6.1 职责

将 Zₚ 上的因式分解结果提升到 Z_{p^k}，精度足以恢复 Z 上的系数。

### 6.2 函数总览

```
__hensel_lift                  M2 顶层入口 (返回 pair<factors, modulus>)
  ├── __mignotte_bound         系数界 (§4.4)
  ├── __hensel_tree_build      构建二叉树 + 初始 Bézout 系数
  │     └── __upoly_gcd_extended   扩展 GCD (§4.1.6)
  └── __hensel_step            单步二次提升
        ├── __upoly_divmod_mod ZZ 系数模除法 (§4.1.4)
        └── __symmetric_mod    对称模 (§4.2.1)
```

### 6.3 Hensel 树节点

```cpp
// Hensel 提升二叉树的内部节点
// 存储一对因子 (g, h) 及其 Bézout 系数 (s, t): s·g + t·h ≡ 1 (mod m)
struct __hensel_node {
    upolynomial_<ZZ> g;     // 左子树因子之积
    upolynomial_<ZZ> h;     // 右子树因子之积
    upolynomial_<ZZ> s;     // Bézout 系数: s·g + t·h ≡ 1
    upolynomial_<ZZ> t;     // Bézout 系数
    int left;               // 左子节点索引 (-1 = 叶子)
    int right;              // 右子节点索引 (-1 = 叶子)
    int leaf_start;         // 对应的叶子范围起始
    int leaf_end;           // 对应的叶子范围结束
};
```

### 6.4 `__hensel_tree_build` — 构建初始提升树

```cpp
// 构建 Hensel 提升的二叉树结构
// 前置: factors 是 f mod p 的首一不可约因子列表，两两互素, |factors| ≥ 2
// 后置: 返回二叉树节点数组 (长度 |factors|-1)
//       每个节点的 g, h 是 mod p 下的因子子积
//       s·g + t·h ≡ 1 (mod p) (由扩展 GCD 计算)
// 树结构: nodes[0] 是根 (g·h ≡ f mod p)
//         叶子对应 factors 中的单个因子
std::vector<__hensel_node> __hensel_tree_build(
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t p)
```

**树的组织：**

对 r 个因子 f₁,...,fᵣ，构建 r-1 个内部节点。采用平衡二叉树：
- 节点 0（根）：g = f₁·...·f_{⌊r/2⌋}, h = f_{⌊r/2⌋+1}·...·fᵣ
- 递归分割直到叶子

每个内部节点的 Bézout 系数 s, t 通过 `__upoly_gcd_extended(g, h)` 计算
（在 mod p 下，类型为 `upolynomial_<Zp>`），然后逐项将系数从 Zp 转为 ZZ
（即 `ZZ(term.second.number())`，值在 `[0, p)` 范围），存入 `__hensel_node`
的 `upolynomial_<ZZ>` 字段供后续提升使用。g, h 同理做 Zp→ZZ 转换。

### 6.5 `__hensel_step` — 单步二次 Hensel 提升

```cpp
// 二次 Hensel 提升: 将 f ≡ g·h (mod m) 提升到 f ≡ g*·h* (mod m²)
// 同时更新 Bézout 系数: s·g + t·h ≡ 1 (mod m) → s*·g* + t*·h* ≡ 1 (mod m²)
//
// 前置: node.g, node.h 满足 f ≡ g·h (mod m)
//       node.s, node.t 满足 s·g + t·h ≡ 1 (mod m)
//       deg(s) < deg(h)
//       m² 不溢出 ZZ (ZZ 无上限)
// 后置: node 中的 g, h, s, t 被更新到 mod m²
//       f ≡ node.g · node.h (mod m²)
//       node.s · node.g + node.t · node.h ≡ 1 (mod m²)
void __hensel_step(
    __hensel_node& node,
    const upolynomial_<ZZ>& f,
    const ZZ& m)            // 当前模数 (提升后变为 m²)
```

**算法（详细）：**

```
输入: node = {g, h, s, t}, f, m (当前 f ≡ g·h mod m)
输出: node 更新到 mod m²

// --- 第一部分: 提升因子 ---

1. e ← f - g·h                        // Z[x] 上精确计算
   对 e 的每个系数: eᵢ ← eᵢ / m       // 精确整除 (由不变量保证)
   对 e 的每个系数: eᵢ ← eᵢ mod m     // 映射到 [0, m)

2. // 解 σ·g + τ·h ≡ e (mod m), 其中 σ = s·e mod h, τ 相应
   se ← s · e                          // Z[x] 乘法
   __upoly_divmod_mod(q, r, se, h, m)  // r = s·e mod h (mod m)

3. g ← g + m · (t·e + q·g)            // 对每个系数 mod m²
   h ← h + m · r                      // 对每个系数 mod m²

// --- 第二部分: 提升 Bézout 系数 ---

4. e' ← 1 - s·g - t·h                 // Z[x] 上精确计算
   对 e' 的每个系数: e'ᵢ ← e'ᵢ / m    // 精确整除
   对 e' 的每个系数: e'ᵢ ← e'ᵢ mod m

5. __upoly_divmod_mod(q', r', s · e', h, m)
   s ← s + m · r'                     // mod m²
   t ← t + m · (t·e' + q'·g)         // mod m²
```

### 6.6 `__hensel_lift` — 多因子 Hensel 提升（M2 入口）

```cpp
// 多因子 Hensel 提升: 从 mod p 提升到 mod p^k
// 前置: f ∈ Z[x] 本原, lc(f) > 0
//       factors 是 f mod p 的首一不可约因子列表
//       lc(f) mod p ≠ 0
//       ∏ factors ≡ (f / lc(f)) (mod p)
//       |factors| ≥ 2
// 后置: 返回 (lifted_factors, modulus) 其中:
//       lifted_factors 是 Z[x] 上的因子列表 (系数 ∈ [0, modulus))
//       modulus = p^{2^j} ≥ p^k (二次提升的实际精度)
//       ∏ Fᵢ ≡ f (mod modulus)
//       每个 Fᵢ 的首项系数 ≡ lc(f) 的对应分配 (mod modulus)
//       k 自动由 __mignotte_bound(f) 确定: p^k > 2·|lc(f)|·B
// 说明: 返回的因子系数在 [0, modulus) 范围，调用方需做 __upoly_symmetric_mod
//       才能得到 Z 上的真实系数
std::pair<std::vector<upolynomial_<ZZ>>, ZZ>
__hensel_lift(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t p)
```

**算法：**

```
1.  B ← __mignotte_bound(f)
    target ← 2 · |lc(f)| · B

2.  // 处理首项系数: 将 lc(f) 分配到 factors[0] 的全部系数
    // 这保证 factors[0]·factors[1]·...·factors[r-1] ≡ f (mod p)
    // (原始的首一因子满足 lc(f)·f₁·...·fᵣ ≡ f mod p，
    //  乘 lc(f) 到 f₁ 后积变为 f 本身)
    factors_adj ← factors 的副本
    factors_adj[0] 的每个系数乘以 lc(f) mod p
    // 注: 是乘以全部系数，不是仅设置首项系数

3.  nodes ← __hensel_tree_build(factors_adj, p)
    // 此时 nodes[0].g · nodes[0].h ≡ f (mod p)

4.  m ← p
    while m ≤ target:
        __hensel_lift_recursive(nodes, 0, f, m)  // 自顶向下递归提升
        m ← m²

    // __hensel_lift_recursive 的递归逻辑:
    //   def __hensel_lift_recursive(nodes, idx, target, m):
    //       __hensel_step(nodes[idx], target, m)
    //       if nodes[idx].left ≠ -1:
    //           __hensel_lift_recursive(nodes, left, nodes[idx].g, m)
    //       if nodes[idx].right ≠ -1:
    //           __hensel_lift_recursive(nodes, right, nodes[idx].h, m)
    //
    // 关键: 根节点以 f 为目标，子节点以父节点更新后的 g/h 为目标。
    // 不能用扁平循环 (所有节点共用 f)，否则 >2 因子时结果错误。

5.  从 nodes 的叶子中提取最终因子 (mod m)

6.  return (因子列表, m)              // m 是最终模数，调用方需要它
```

---

## 7. 模块 M3：因子重组

### 7.1 职责

给定 f ∈ Z[x] 及其 mod p^k 的因子列表（精度足够），恢复 Z[x] 上的
真正不可约因子。

### 7.2 核心难点

Hensel 提升保持因子数量：如果 f mod p 有 r 个不可约因子，提升后仍有
r 个。但 Z 上的真因子可能需要组合多个模因子：

例：f = x⁴ + 1 在 Z[x] 上不可约，但 mod 2: f ≡ (x+1)⁴，有 4 个模因子。

### 7.3 函数总览

```
__factor_recombine                  M3 入口 (Zassenhaus)
  ├── __upoly_symmetric_mod         对称模约化 (§4.2.2)
  ├── __upoly_norm_l1               L1 范数 (§4.3)
  ├── __upoly_primitive             本原化 (§4.5)
  └── __subset_product_mod          子集乘积 mod m
```

### 7.4 `__subset_product_mod` — 子集乘积

```cpp
// 计算 factors 中下标在 subset 中的因子之积, 系数 mod m
// 前置: subset ⊆ {0, ..., factors.size()-1}, m > 0
// 后置: 返回 lc_f · ∏_{i∈subset} factors[i] (mod m)
//       系数经过 __symmetric_mod 约化到 (-m/2, m/2]
inline upolynomial_<ZZ> __subset_product_mod(
    const std::vector<upolynomial_<ZZ>>& factors,
    const std::vector<size_t>& subset,
    const ZZ& lc_f,
    const ZZ& m)
```

### 7.5 `__factor_recombine` — Zassenhaus 因子重组

```cpp
// Zassenhaus 因子重组
// 前置: f 本原, lc(f) > 0, 无平方
//       lifted 是 f mod m 的因子分解 (由 __hensel_lift 产出)
//       m 足够大: m > 2 · |lc(f)| · __mignotte_bound(f)
// 后置: 返回 f 在 Z[x] 上的不可约因子列表
//       ∏(返回因子) = f (精确)
//       每个因子本原, lc > 0
//       因子按 (degree, then lex) 排列
std::vector<upolynomial_<ZZ>>
__factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ& m)
```

**算法（详细）：**

```
输入: f, lifted = [F₁,...,Fᵣ], m
输出: f 的不可约因子列表

1.  T ← {0, 1, ..., r-1}               // 可用因子下标集
    f* ← f                              // 剩余待分解多项式
    result ← []

2.  for s = 1 to ⌊|T|/2⌋:

3.      for S ⊆ T, |S| = s:            // 枚举所有大小为 s 的子集
            // === 快速剪枝 ===

            // 剪枝 1: 首项系数检查
            lc_prod ← lc(f*) · ∏_{i∈S} lc(Fᵢ) mod m
            lc_prod ← __symmetric_mod(lc_prod, m)
            if lc_prod 不整除 lc(f*)²:
                continue

            // 剪枝 2: 常数项检查
            c_prod ← lc(f*) · ∏_{i∈S} Fᵢ(0) mod m
            c_prod ← __symmetric_mod(c_prod, m)
            if c_prod 不整除 (lc(f*) · f*(0)):
                continue

            // === 完整验证 ===

            g ← __subset_product_mod(lifted, S, lc(f*), m)
            // g = lc(f*) · ∏_{i∈S} Fᵢ mod m, 系数已对称约化

            // 试除: g | lc(f*)·f* in Z[x]?
            lc_f_star ← lc(f*)
            f_scaled ← lc_f_star · f*
            q, r ← divmod(f_scaled, g) in Z[x]
            if r = 0:
                // 找到真因子 g
                pp_g ← __upoly_primitive(g).second
                result.push_back(pp_g)
                f* ← __upoly_primitive(q).second   // f* ← pp(h) 等价
                T ← T_comp
                goto step_2_restart                 // 重新从 s=1 开始

4.  // 剩余的 f* 本身不可约
    if deg(f*) > 0:
        result.push_back(f*)

5.  排序 result

6.  return result
```

**子集枚举实现：** 使用位掩码枚举所有大小为 s 的子集：
```cpp
// Gosper's hack: 枚举 n 位中恰好 s 个 1 的所有位掩码
// 限制: r ≤ 64 (uint64_t 位宽)
// 对于 r > 64 的极端情况 (如 Swinnerton-Dyer 多项式)，
// 应切换到 van Hoeij 算法 (Phase 6)
uint64_t mask = (1ULL << s) - 1;
while (mask < (1ULL << r)):
    // 处理 mask
    uint64_t c = mask & (-mask);
    uint64_t r_ = mask + c;
    mask = (((r_ ^ mask) >> 2) / c) | r_;
```

### 7.6 剪枝优化方向

当前实现仅使用剪枝 1（首项系数整除）和剪枝 2（常数项整除）。

> **注：** 原设计中的"剪枝 3: 范数检查"（`||g||₁·||h||₁ > ||lc(f*)·f*||₁`）
> 已移除。L1 范数满足次可乘性 `||g·h||₁ ≤ ||g||₁·||h||₁`，对真因子也几乎
> 总有严格不等式（因乘法展开时正负消去），导致有效分解被错误剪掉。

以下是主流代数库的做法，按性价比排序，可作为后续增强参考：

**1. Newton 迹检查（PARI/GP, NTL）— 推荐优先实现**

从提升因子的前 1-2 个系数预计算 Newton 迹（trace）：
- 1 阶迹：`t₁(S) = Σ_{i∈S} coeff_{n-1}(Fᵢ)` （对应根之和）
- 2 阶迹：`t₂(S) = Σ_{i∈S} coeff_{n-2}(Fᵢ)` （对应根平方和相关量）

对真因子 g（deg = k），有 `|t₁| ≤ k·||g||∞`，可用 Mignotte bound 给出阈值。
仅需 ulong 算术，O(|S|) 计算量，能淘汰绝大多数假子集。

**2. 试除时逐系数 bound 检查（PARI/GP）**

在执行多项式长除法的同时，每算出一个商系数就与 Mignotte bound 比较，
超限立即中止。无需额外计算，只需修改试除函数增加 bound 参数。
PARI 的 `ZX_divides_i(f, g, bound)` 即是此做法。

**3. 多素数度数可行性（FLINT, NTL）**

对 f 分别 mod p₁, p₂, ... 做 Zp 分解，收集各素数下因子度数的并集。
真因子的度数必须是所有素数下因子度数子集和的交集中的元素。
用位掩码 DP 实现，setup 开销较大但对 r 大时效果显著。
- FLINT 使用 3 个素数
- NTL 动态添加最多 50 个素数

**4. Beauzamy bound 替代 Mignotte bound（PARI/GP）**

PARI 同时计算 Mignotte bound 和 Beauzamy bound 并取较小者：
```
B_beauzamy = |lc(f)| · √(3^(3/2+n) · [f]₂² / (4nπ))
```
其中 `[f]₂² = Σ aᵢ²/C(n,i)` 为加权 L2 范数。实际中往往比 Mignotte 更紧。

**5. Meet-in-the-middle 查表（NTL）**

对 r 中最后 t 个因子（`ZZXFac_MaxPrune=10`），预建大小 ~t·2^(t-1) 的
查表。枚举子集时查表判断 trace 值是否可行。运行时间缩减约 2^t 倍。
仅在 r 较大时有收益。

**各库剪枝策略总结：**

| 库 | 度数 | 常数项 | Newton 迹 | 系数 bound | 查表 |
|---|---|---|---|---|---|
| Singular | ✓ | | | | |
| FLINT | ✓(多素数) | | | | |
| SymPy | | ✓ | | | |
| PARI/GP | | ✓ | ✓(1+2 阶) | ✓(试除时) | |
| NTL | ✓(多素数) | ✓ | ✓(1+2 阶) | ✓ | ✓ |

### 7.7 van Hoeij 重组（远期增强）

当模因子数 r 较大时（例如 r > 15-20），Zassenhaus 可能过慢。
van Hoeij 算法通过格约化在多项式时间内完成重组：

**核心思想：**

1. 构造 r×(r+m) 格矩阵 M：
   - 前 r 列：单位矩阵 I_r（缩放为 N·I_r，N 是选定的缩放参数）
   - 后 m 列：从提升因子的第 j 个系数中提取的信息
2. 对 M 做 LLL 约化
3. 约化基中的短向量，其前 r 个分量中非零位置指示了哪些模因子组合成同一个
   真因子

**依赖：** LLL 格基约化算法的实现。

**实现策略：** 作为 Phase 3+ 的增强。如果实际使用中遇到 r 很大的情况再添加。
可以考虑调用外部 LLL 库（如 fplll）或自行实现。

---

## 8. 模块 M4：单变量 Z[x] 因式分解（集成）

### 8.1 职责

端到端的单变量整系数多项式因式分解，整合 M1-M3 及已有的无平方分解。

### 8.2 函数总览

```
factorize(polynomial_<ZZ,comp>)      用户入口 (任意 ordering)
  └── factorize(polynomial_<ZZ,lex>) lex 特化版本
        ├── cont / pp                内容提取 (已有)
        ├── squarefreefactorize      无平方分解 (已有)
        └── __factor_squarefree_ZZ   对单个无平方因子分解
              ├── __select_prime     素数选择
              │     ├── polynomial_mod    模约化 (已有)
              │     ├── __squarefree_Zp   Zp 无平方检测
              │     └── __factor_Zp       M1 入口
              ├── __hensel_lift      M2 入口
              └── __factor_recombine M3 入口
```

### 8.3 `__select_prime` — 素数选择

```cpp
// 选择最优素数用于因式分解
// 前置: f 本原, 无平方, lc(f) > 0, deg(f) ≥ 2
// 后置: 返回 (best_p, best_factors)
//       best_p 是素数, lc(f) mod best_p ≠ 0, f mod best_p 无平方
//       best_factors 是 f mod best_p 的不可约因子列表 (不含 lc)
//       best_factors.size() 在所有候选素数中最小
//       如果某个素数下只有 1 个因子, 则 best_factors.size() = 1
//       (表示 f 不可约, 调用方据此直接返回)
//
// 常量:
//   INIT_NUM_PRIMES = 5   初始尝试素数个数
//   MAX_NUM_PRIMES = 20   最大尝试素数个数
std::pair<uint32_t, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
__select_prime(const upolynomial_<ZZ>& f)
```

**算法：**

```
1.  lc_f ← lc(f)
    best_p ← 0
    best_count ← ∞
    best_factors ← {}

2.  i ← 0                              // 尝试的素数计数
    prime_idx ← 0                      // boost::math::prime 的索引

3.  while i < INIT_NUM_PRIMES (or i < MAX_NUM_PRIMES if best_count > deg(f)/2):
        p ← boost::math::prime(prime_idx++)

        // 筛选坏素数
        if lc_f mod p = 0: continue

        f_p ← polynomial_mod(f, p)      // 已有: upolynomial.hh
        f_p ← f_p / lc(f_p)             // 首一化

        // 检查 f mod p 是否无平方
        f_p_deriv ← derivative(f_p)
        g ← gcd(f_p, f_p_deriv)         // __polynomial_GCD for Zp
        if deg(g) > 0: continue          // f mod p 不无平方，跳过

        // 分解
        factors_p ← __factor_Zp(f_p)

        // 不可约快速路径
        if factors_p.second.size() = 1:
            return (p, factors_p.second)

        // 记录最优
        if factors_p.second.size() < best_count:
            best_p ← p
            best_count ← factors_p.second.size()
            best_factors ← factors_p.second

        i++

4.  return (best_p, best_factors)
```

### 8.4 `__factor_squarefree_ZZ` — 分解单个无平方因子

```cpp
// 分解一个无平方本原多项式
// 前置: f 本原, 无平方, lc(f) > 0, deg(f) ≥ 2
// 后置: 返回 f 的不可约因子列表
//       ∏(因子) = f
//       每个因子本原, lc > 0
std::vector<upolynomial_<ZZ>>
__factor_squarefree_ZZ(const upolynomial_<ZZ>& f)
```

**算法：**

```
1.  (p, mod_factors) ← __select_prime(f)

2.  if mod_factors.size() = 1:
        return {f}                       // f 不可约

3.  // 提取首一因子 (去掉 __factor_Zp 返回的 lc)
    factors_monic ← mod_factors 中每个 pair 的 .first

4.  (lifted, modulus) ← __hensel_lift(f, factors_monic, p)   // M2

5.  result ← __factor_recombine(f, lifted, modulus)          // M3

6.  return result
```

### 8.5 `factorize` — 用户入口（lex 特化）

```cpp
// 单变量因式分解 (lex ordering 特化, 核心实现)
// 前置: f 非零
// 后置: 返回 factorization 满足 content · ∏ fᵢ^eᵢ = f
//       每个 fᵢ 不可约, 本原, lc > 0
//       eᵢ > 0
//       因子按 (degree asc, then lex) 排列
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
factorize(const polynomial_<ZZ, lex_<var_order>>& f)
```

**算法：**

```
1.  if f.empty():
        return {ZZ(0), {}}

2.  if is_number(f):
        return {f.front().second, {}}    // 常数

3.  vars ← get_variables(f)
    if vars.size() > 1:
        return __factor_multivar(f)      // M5 (远期)

4.  // --- 单变量路径 ---
    // 转换为 upolynomial_<ZZ>
    upolynomial_<ZZ> uf;
    poly_convert(f, uf);

5.  c ← cont(uf)                        // 已有: upolynomial.hh
    if uf.front().second < 0: c = -c
    uf ← uf / c                         // 本原化, lc > 0

6.  // 方案 A (§8.7): 转为 polynomial_<ZZ,lex> 调用已有 squarefreefactorize
    polynomial_<ZZ, lex> f_lex;
    poly_convert_back(uf, f_lex, var)     // upolynomial → polynomial (带回变量)
    sqf_lex ← squarefreefactorize(f_lex)  // 已有，成熟且有测试覆盖
    // 将各因子转回 upolynomial_<ZZ>
    sqf ← [(poly_convert(gᵢ_lex) → upolynomial, eᵢ) for (gᵢ_lex, eᵢ) in sqf_lex]

7.  result ← factorization<...>{c, {}}

8.  for (gᵢ, eᵢ) in sqf:
        if deg(gᵢ) ≤ 0:
            result.content *= gᵢ.front().second ^ eᵢ
            continue
        if deg(gᵢ) = 1:
            result.factors.push_back({convert_to_poly(gᵢ), eᵢ})
            continue

        irr_factors ← __factor_squarefree_ZZ(gᵢ)
        for h in irr_factors:
            result.factors.push_back({convert_to_poly(h), eᵢ})

9.  排序 result.factors

10. return result
```

### 8.6 `factorize` — 用户入口（通用 ordering 包装）

```cpp
// 通用 ordering 版本: 转换到 lex, 计算, 转换回来
// 与 polynomial_GCD 的 generic wrapper 风格一致
template<class comp>
factorization<polynomial_<ZZ, comp>>
factorize(const polynomial_<ZZ, comp>& f)
{
    polynomial_<ZZ, lex> f_lex;
    poly_convert(f, f_lex);
    auto result_lex = factorize(f_lex);

    factorization<polynomial_<ZZ, comp>> result;
    result.content = result_lex.content;
    result.factors.reserve(result_lex.factors.size());
    polynomial_<ZZ, comp> tmp(f.comp_ptr());
    for (auto& [fi, ei] : result_lex.factors)
    {
        poly_convert(fi, tmp);
        result.factors.push_back({std::move(tmp), ei});
        tmp = polynomial_<ZZ, comp>(f.comp_ptr());
    }
    return result;
}
```

### 8.7 QQ[x] 因式分解入口

```cpp
// QQ[x] 因式分解: 约化为 ZZ[x] 分解
// 前置: f ∈ QQ[x], 非零
// 后置: 返回 factorization 满足 content · ∏ fᵢ^eᵢ = f
//       content ∈ QQ, 每个 fᵢ ∈ QQ[x] 首一不可约
// 算法: (1) 提取有理内容并乘以 LCD 得 g ∈ ZZ[x]
//       (2) 分解 g
//       (3) 将因子转为 QQ[x] 首一多项式，吸收 lc 到 content
template<class comp>
factorization<polynomial_<QQ, comp>>
factorize(const polynomial_<QQ, comp>& f)
```

**算法：**

```
1.  若 f = 0: return {QQ(0), {}}
2.  lcd ← f 所有系数分母的 LCM
3.  g ← lcd · f ∈ ZZ[x]
4.  result_zz ← factorize(g)              // ZZ[x] 分解
5.  // 构造 QQ 结果: content = result_zz.content / lcd
    // 每个 ZZ 因子 → QQ 首一因子 (除以 lc), lc 吸收到 content
6.  return QQ 格式的 factorization
```

### 8.8 无平方分解的 upolynomial 适配

已有的 `squarefreefactorize` 工作在 `polynomial_<ZZ, lex>` 上。对于单变量
流程需要在 `upolynomial_<ZZ>` 上操作。两个方案：

- **方案 A（推荐）：** 在 `factorize` 内部先用 `poly_convert` 转为
  `polynomial_<ZZ, lex>` 做无平方分解，再转回 `upolynomial_<ZZ>` 做后续
- **方案 B：** 新写 `__squarefree_ZZ_upoly`

选择方案 A，因为 `squarefreefactorize` 已经成熟且有测试覆盖，无需重复实现。

---

## 9. 模块 M5：多变量因式分解（远期）

### 9.1 概述

多变量因式分解通过将问题约化为单变量来解决。核心思想：

1. 取值 x₂=α₂, ..., xₙ=αₙ 得到单变量像
2. 对像进行单变量分解（M4）
3. 通过多变量 Hensel 提升恢复各变量的贡献

### 9.2 函数总览

```
__factor_multivar(polynomial_<ZZ,lex>)     M5 入口
  ├── cont(f, x₁)                         内容提取 (已有)
  ├── factorize(cont)                      递归分解内容
  ├── __select_eval_point                  选取值点
  │     ├── assign(f, v, c)               代入求值 (已有)
  │     └── is_squarefree                 无平方检测 (已有)
  ├── factorize(f₀)                        单变量分解 (M4)
  ├── __wang_leading_coeff                 首项系数分配
  └── __multivar_hensel_lift               多变量 Hensel 提升
        └── __upoly_gcd_extended           扩展 GCD (§4.1.6)
```

### 9.3 `__select_eval_point` — 选取值点

```cpp
// 选择多变量 → 单变量的取值点
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2, x₁ 是主变量
// 后置: 返回 α = {x₂→α₂, ..., xₙ→αₙ} 使得:
//       (a) f(x₁, α₂, ..., αₙ) 无平方
//       (b) lc(f, x₁) |_{xᵢ=αᵢ} ≠ 0
//       (c) deg(f(x₁,...)) = deg(f(x₁, α₂, ...)) (首项不消失)
// 策略: 从小整数开始尝试 (0, ±1, ±2, ...)
std::map<variable, ZZ>
__select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& main_var)
```

### 9.4 `__multivar_hensel_lift` — 多变量 Hensel 提升

```cpp
// Wang 多变量 Hensel 提升
// 将 f(x₁,α₂,...,αₙ) 的因子提升为 f(x₁,...,xₙ) 的因子
// 前置: univar_factors 是 f(x₁,α₂,...,αₙ) 的因子
//       lc_info 包含首项系数的分配信息
// 后置: 返回 f(x₁,...,xₙ) 在 Z[x₁,...,xₙ] 上的候选因子
//       (需要调用方做试除验证)
// 算法: 逐个变量提升 — 先恢复 x₂, 再 x₃, ..., 最后 xₙ
//       每个变量通过对 (xᵢ - αᵢ) 的幂展开做线性 Hensel 提升
std::vector<polynomial_<ZZ, lex_<var_order>>>
__multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var)
```

### 9.5 `__factor_multivar` — 多变量分解入口

```cpp
// 多变量因式分解 (Wang 算法)
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2
// 后置: 返回 factorization
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f)
```

**算法：**

```
1.  x₁ ← 选主变量 (度数最高的, 用 get_variables)
    c ← cont(f, x₁)              // ∈ Z[x₂,...,xₙ]
    f ← pp(f, x₁)

2.  // 递归分解内容
    if !is_number(c):
        cont_factors ← factorize(c)

3.  // 选取值点
    eval ← __select_eval_point(f, x₁)

4.  // 单变量分解
    f₀ ← assign(f, eval)          // f₀ ∈ Z[x₁]
    uni_result ← factorize(f₀)

5.  if uni_result.factors.size() ≤ 1:
        // f 关于 x₁ 不可约 (但可能关于其他变量可分解)
        return ...

6.  // 首项系数校正
    __wang_leading_coeff(f, uni_result.factors, eval, x₁)

7.  // 多变量 Hensel 提升
    mv_factors ← __multivar_hensel_lift(f, uni_result, eval, x₁)

8.  // 试除验证
    verified ← []
    f* ← f
    for g in mv_factors:
        if g | f*:
            verified.push(g)
            f* ← f* / g
    if deg(f*) > 0:
        verified.push(f*)

9.  // 合并内容因子
    return combine(cont_factors, verified)
```

### 9.6 Wang 算法接口

```cpp
// 首项系数校正: 将 lc(f, x₁) 的因式分解信息分配到单变量因子
// 这是 Wang 算法的关键步骤
// 前置: lc(f, x₁) 的因子已知, 单变量因子列表, 取值点
// 后置: 修改 univar_factors 的首项系数以包含多变量信息
void __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var)
```

---

## 10. 公共数据结构

### 10.1 因式分解结果类型

```cpp
// 因式分解结果
// 表示 f = content · ∏ factors[i].first ^ factors[i].second
// 遵循 CLPoly 的 pair vector 风格
template<typename Poly>
struct factorization {
    typename Poly::coeff_type content;                   // 内容 (整数/有理数)
    std::vector<std::pair<Poly, uint64_t>> factors;      // (不可约因子, 正整数重数)
    // 使用 uint64_t 与 squarefreefactorize 的返回类型一致

    // 规范化保证:
    // - 每个因子本原（primitive），首项系数为正
    // - 因子按 (degree asc, then polynomial < operator) 排列
    // - 重数 > 0
    // - content 可以是负数（吸收符号）
    // - 若 f = 0, 则 content = 0, factors 为空

    // 便利方法
    bool is_irreducible() const { return factors.size() == 1 && factors[0].second == 1; }
    Poly expand() const;                                 // 重建原多项式 (用于验证)
};
```

### 10.2 与已有类型的兼容

当前 `squarefreefactorize` 返回 `vector<pair<polynomial, uint64_t>>`。
新的 `factorize` 使用 `factorization<Poly>` 结构体包装，增加 content 字段。

两个接口独立：`squarefreefactorize` 保持不变（无平方分解 ≠ 不可约分解），
`factorize` 是新增功能。重数类型统一使用 `uint64_t`，与已有风格一致。

---

## 11. 实现路线图

| 阶段 | 内容 | 新增函数 | 依赖 |
|---|---|---|---|
| **Phase 1** | M1: Zₚ 上分解 | `__squarefree_Zp`, `__extract_pth_root`, `__ddf_Zp`, `__edf_Zp`, `__factor_Zp` + 辅助函数 (§4.1: `__upoly_make_monic`, `__upoly_mod`, `__upoly_divmod`, `__upoly_gcd_Zp`, `__upoly_powmod`, `__upoly_random`) | Zp 类, upolynomial |
| **Phase 2** | M2: Hensel 提升 | `__hensel_step`, `__hensel_tree_build`, `__hensel_lift` + `__mignotte_bound`, `__upoly_gcd_extended`, `__upoly_divmod_mod` | M1 |
| **Phase 3** | M3: Zassenhaus 重组 | `__factor_recombine`, `__subset_product_mod` + `__symmetric_mod`, `__upoly_symmetric_mod`, 范数函数 | M2 |
| **Phase 4** | M4: 单变量集成 | `__select_prime`, `__factor_squarefree_ZZ`, `factorize` (用户入口) + `factorization<>` 结构体 | M1+M2+M3 |
| **Phase 5** | M5: 多变量 Wang | `__select_eval_point`, `__wang_leading_coeff`, `__multivar_hensel_lift`, `__factor_multivar` | M4 |
| **Phase 6** | 增强：van Hoeij | `__factor_recombine_van_hoeij` + LLL 实现 | M3 替换 |

**Phase 1-4 是最小可用版本**，完成后可处理所有单变量 Z[x] 的因式分解。

### 11.1 每阶段验证策略

- 使用 Mathematica (`/home/ker/.local/bin/math`) 生成测试数据
- 扩展 `test/generate_testdata.wl`，新增因式分解测试用例
- 关键测试类别：
  - 常数、线性、二次多项式
  - 已知不可约多项式（如分圆多项式 Φₙ(x)）
  - 完全分裂多项式（如 ∏(x - aᵢ)）
  - Swinnerton-Dyer 多项式（模因子数极多，测试重组）
  - 含重因子的多项式
  - 系数很大的多项式

### 11.2 每阶段可独立测试的单元

| 阶段 | 可独立测试的函数 | 验证方法 |
|---|---|---|
| Phase 1 | `__squarefree_Zp` | 与 Mathematica `FactorSquareFree[f, Modulus->p]` 对比 |
| Phase 1 | `__ddf_Zp` | 检查 ∏gd = f, 且每个 gd 的因子度数均为 d |
| Phase 1 | `__edf_Zp` | 检查返回的因子均不可约且 ∏ = f |
| Phase 1 | `__factor_Zp` | 与 Mathematica `Factor[f, Modulus->p]` 对比 |
| Phase 2 | `__upoly_gcd_extended` | 验证 s·a + t·b = 1 |
| Phase 2 | `__hensel_step` | 验证 f ≡ g·h (mod m²) |
| Phase 2 | `__hensel_lift` | 验证 ∏Fᵢ ≡ f (mod p^k) |
| Phase 3 | `__factor_recombine` | 验证 ∏(结果) = f in Z[x] |
| Phase 4 | `factorize` | 与 Mathematica `FactorList[f]` 完整对比 |

---

## 12. 设计决策记录

| 决策 | 选择 | 理由 |
|---|---|---|
| Zₚ 分解算法 | Cantor-Zassenhaus | 适合大 p，实现简洁，性能好 |
| Hensel 提升方式 | 二次提升（精度每步翻倍） | 标准选择，收敛快 |
| 多因子提升结构 | 二叉树 | 平衡了提升次数和因子合并开销 |
| 初始重组算法 | Zassenhaus | 实现简单，实际中 r 通常很小 |
| 多变量策略 | Wang | 经典成熟，与现有 Hensel 基础设施兼容 |
| 文件组织 | 单文件 `polynomial_factorize.hh` | 与 `polynomial_gcd.hh` 风格一致 |
| 结果类型 | `factorization<Poly>` 模板结构体 | 通用，类型安全，可扩展到 QQ 等 |
| 素数选择 | 多素数竞选，选模因子数最少的 | NTL/FLINT 验证过的策略 |
| 重数类型 | `uint64_t` | 与已有 `squarefreefactorize` 一致 |
| 内部函数命名 | `__` 前缀 | 与 `__polynomial_GCD` 风格一致 |
| 参数传递 | 值传入可修改的，const ref 传入只读的 | 与 `polynomial_GCD(F, G)` 风格一致 |
| 单变量表示 | `upolynomial_<Zp>` / `upolynomial_<ZZ>` | 避免多变量开销，算法更高效 |
| 排序 ordering | lex 特化 + generic wrapper | 与 `polynomial_GCD` / `squarefreefactorize` 一致 |

---

## 附录 A: 完整函数签名索引

### A.1 辅助函数 (§4)

```cpp
// §4.1.1 首一化
inline Zp __upoly_make_monic(upolynomial_<Zp>& f);

// §4.1.2 多项式取模
inline upolynomial_<Zp> __upoly_mod(
    const upolynomial_<Zp>& f, const upolynomial_<Zp>& g);

// §4.1.3 商和余式
inline void __upoly_divmod(
    upolynomial_<Zp>& q, upolynomial_<Zp>& r,
    const upolynomial_<Zp>& f, const upolynomial_<Zp>& g);

// §4.1.4 ZZ 系数模除法
inline void __upoly_divmod_mod(
    upolynomial_<ZZ>& q, upolynomial_<ZZ>& r,
    const upolynomial_<ZZ>& f, const upolynomial_<ZZ>& g,
    const ZZ& m);

// §4.1.5 Zₚ 上简单 GCD 包装
inline upolynomial_<Zp> __upoly_gcd_Zp(
    const upolynomial_<Zp>& a, const upolynomial_<Zp>& b);

// §4.1.6 扩展 GCD
inline void __upoly_gcd_extended(
    upolynomial_<Zp>& s, upolynomial_<Zp>& t,
    const upolynomial_<Zp>& a, const upolynomial_<Zp>& b);

// §4.1.7 模幂
inline upolynomial_<Zp> __upoly_powmod(
    const upolynomial_<Zp>& base, const ZZ& exp,
    const upolynomial_<Zp>& modpoly);

// §4.1.8 随机多项式
inline upolynomial_<Zp> __upoly_random(
    int64_t max_deg, uint32_t p, std::mt19937& rng);

// §4.2 对称模
inline ZZ __symmetric_mod(const ZZ& a, const ZZ& m);
inline upolynomial_<ZZ> __upoly_symmetric_mod(
    const upolynomial_<ZZ>& f, const ZZ& m);

// §4.3 范数
inline ZZ __upoly_norm_l1(const upolynomial_<ZZ>& f);
inline ZZ __upoly_norm_l2_sq(const upolynomial_<ZZ>& f);
inline ZZ __upoly_norm_linf(const upolynomial_<ZZ>& f);

// §4.4 Mignotte 界
inline ZZ __mignotte_bound(const upolynomial_<ZZ>& f);

// §4.5 本原化
inline std::pair<ZZ, upolynomial_<ZZ>> __upoly_primitive(upolynomial_<ZZ> f);
```

### A.2 M1: Zₚ 分解 (§5)

```cpp
// 顶层入口
std::pair<Zp, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
__factor_Zp(upolynomial_<Zp> f);

// 无平方分解
std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
__squarefree_Zp(const upolynomial_<Zp>& f);

// p 次根提取
inline upolynomial_<Zp> __extract_pth_root(const upolynomial_<Zp>& f);

// DDF
std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
__ddf_Zp(const upolynomial_<Zp>& f);

// EDF
void __edf_Zp(
    std::vector<upolynomial_<Zp>>& result,
    const upolynomial_<Zp>& f, uint64_t d, std::mt19937& rng);
```

### A.3 M2: Hensel 提升 (§6)

```cpp
// Hensel 树节点
struct __hensel_node {
    upolynomial_<ZZ> g, h, s, t;
};

// 构建提升树
std::vector<__hensel_node> __hensel_tree_build(
    const std::vector<upolynomial_<Zp>>& factors, uint32_t p);

// 单步二次提升
void __hensel_step(
    __hensel_node& node, const upolynomial_<ZZ>& f, const ZZ& m);

// 多因子提升入口
std::pair<std::vector<upolynomial_<ZZ>>, ZZ> __hensel_lift(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<Zp>>& factors, uint32_t p);
```

### A.4 M3: 因子重组 (§7)

```cpp
// 子集乘积
inline upolynomial_<ZZ> __subset_product_mod(
    const std::vector<upolynomial_<ZZ>>& factors,
    const std::vector<size_t>& subset, const ZZ& lc_f, const ZZ& m);

// Zassenhaus 重组
std::vector<upolynomial_<ZZ>> __factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted, const ZZ& m);
```

### A.5 M4: 单变量集成 (§8)

```cpp
// 素数选择
std::pair<uint32_t, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
__select_prime(const upolynomial_<ZZ>& f);

// 无平方因子分解
std::vector<upolynomial_<ZZ>>
__factor_squarefree_ZZ(const upolynomial_<ZZ>& f);

// 用户入口 (lex 特化)
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
factorize(const polynomial_<ZZ, lex_<var_order>>& f);

// 用户入口 (generic ordering)
template<class comp>
factorization<polynomial_<ZZ, comp>>
factorize(const polynomial_<ZZ, comp>& f);

// 用户入口 (QQ, §8.7)
template<class comp>
factorization<polynomial_<QQ, comp>>
factorize(const polynomial_<QQ, comp>& f);
```

### A.6 M5: 多变量分解 (§9, 远期)

```cpp
template<class var_order>
std::map<variable, ZZ> __select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f, const variable& main_var);

template<class var_order>
void __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>> __multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f);
```

---

## 附录 B: 参考文献

- **Cantor-Zassenhaus**: Cantor & Zassenhaus, "A New Algorithm for Factoring Polynomials Over Finite Fields", Math. Comp. 1981
- **Hensel Lifting**: Zassenhaus, "On Hensel Factorization I", J. Number Theory 1969
- **van Hoeij**: van Hoeij, "Factoring Polynomials and the Knapsack Problem", J. Number Theory 2002
- **Wang**: Wang, "An Improved Multivariate Polynomial Factoring Algorithm", Math. Comp. 1978
- **Mignotte Bound**: Mignotte, "An Inequality About Factors of Polynomials", Math. Comp. 1974
- **Kaltofen-Shoup**: Kaltofen & Shoup, "Subquadratic-Time Factoring of Polynomials over Finite Fields", Math. Comp. 1998
- **NTL**: Shoup, NTL: A Library for doing Number Theory, https://libntl.org
- **FLINT**: Hart et al., FLINT: Fast Library for Number Theory, https://flintlib.org
- **GCL**: Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992 (§8, §15, §16)
- **MCA**: von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013 (§14-§16)

## 附录 C: 各系统因子表示对比

| 系统 | 内容(content) | 因子列表 | 规范化 |
|---|---|---|---|
| FLINT | `fmpz_t c` 在结构体内 | `fmpz_poly_struct *p` + `slong *exp` | 因子本原，lc > 0 |
| NTL | 单独 `ZZ& c` 参数 | `vec_pair_ZZX_long` | 因子首一 (monic) |
| Mathematica | `FactorList` 第一项 | `{{f₁,e₁},...}` | 因子首项正 |
| Singular | `CFFList` 含 unit | `CFFactor` 列表 | 取决于域 |
| **CLPoly (本设计)** | `factorization::content` | `vector<pair<Poly,uint64_t>>` | 因子本原，lc > 0 |

注意：NTL 选择首一化因子（content 吸收所有首项系数），FLINT 选择本原化因子
（content 只含整数部分）。本设计选择 FLINT 风格（本原 + lc > 0），因为与
CLPoly 已有 `squarefreefactorize` 的风格一致。

## 附录 D: 已有函数依赖清单

以下是本设计直接依赖的已有函数及其位置：

| 函数 | 文件 | 用途 |
|---|---|---|
| `polynomial_mod(f, p)` | `polynomial.hh:677`, `upolynomial.hh:104` | ZZ→Zp 系数约化 |
| `__polynomial_GCD(Pout, G, F, Lc, deg)` | `polynomial_gcd.hh:836` | Zₚ 上 Euclidean GCD |
| `polynomial_GCD(F, G)` | `polynomial_gcd.hh:246` | Z 上 GCD |
| `squarefreefactorize(F)` | `polynomial_gcd.hh:101` | Z 上无平方分解 |
| `cont(F)` | `polynomial_gcd.hh:468`, `upolynomial.hh:160` | 内容提取 |
| `derivative(f)` | `upolynomial.hh:141` | 求导 |
| `pair_vec_div(q, r, f, g, comp)` | `basic.hh:568,698` | 多项式除法 (商+余) |
| `poly_convert(in, out)` | `upolynomial.hh:117` | polynomial ↔ upolynomial 转换 |
| `get_variables(f)` | `polynomial_.hh:219` | 获取变量列表 |
| `assign(f, v, c)` | `polynomial.hh:707` | 变量代入 |
| `is_number(f)` | `upolynomial.hh:155` | 常数检测 |
| `is_squarefree(f)` | `polynomial_gcd.hh:77` | 无平方检测 |
| `pow(Zp, n)` | `number.hh:216` | Zₚ 快速幂 |
| `inv_prime(i, p)` | `number.hh:72` | 模逆 |
| `boost::math::prime(n)` | `<boost/math/...>` | 第 n 个素数 |
