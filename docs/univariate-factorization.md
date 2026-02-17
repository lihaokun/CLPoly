# CLPoly 单变量多项式因式分解

> 本文档描述 CLPoly 中已实现的单变量多项式因式分解功能。
> 所有代码位于 `clpoly/polynomial_factorize.hh`。

---

## 1. 背景与总览

### 1.1 已有基础设施

CLPoly 在实现因式分解之前已具备以下基础设施：

| 组件 | 位置 |
|---|---|
| 无平方分解 `squarefreefactorize` | `polynomial_gcd.hh` |
| 多变量 GCD（模 + Hensel） | `polynomial_gcd.hh` |
| 有限域 `Zp` 类 | `number.hh` |
| `polynomial_mod(f, p)` 模约化 | `polynomial.hh`, `upolynomial.hh` |
| 单变量 Zₚ 上的 Euclidean GCD | `polynomial_gcd.hh` |
| 结式 / 子结式 / 判别式 | `resultant.hh` |
| 伪余式 / 伪商 | `polynomial.hh` |

### 1.2 因式分解管线

单变量因式分解的完整流程：

```
输入 f ∈ Z[x]
  │
  ▼
(1) 内容提取: f = cont(f) · pp(f)
  │
  ▼
(2) 无平方分解: pp(f) = ∏ gᵢ^i
  │
  ▼
(3) 对每个无平方因子 gᵢ:
  │  ├─ 选素数 p，在 Zₚ 上因式分解     (M1)
  │  ├─ Hensel 提升: mod p → mod p^k   (M2)
  │  └─ 因子重组: 组合提升因子为真因子    (M3)
  │
  ▼
输出: cont(f) · ∏ (不可约因子ⱼ)^(指数ⱼ)
```

### 1.3 设计选择

| 决策 | 选择 | 理由 |
|---|---|---|
| Zₚ 分解算法 | Cantor-Zassenhaus | 适合大 p，实现简洁 |
| Hensel 提升方式 | 二次提升（精度每步翻倍） | 标准选择，收敛快 |
| 多因子提升结构 | 二叉树 | 平衡了提升次数和因子合并开销 |
| 初始重组算法 | Zassenhaus | 实现简单，实际中因子数通常很小 |
| 文件组织 | 单文件 `polynomial_factorize.hh` | 与 `polynomial_gcd.hh` 风格一致 |
| 结果类型 | `factorization<Poly>` 模板结构体 | 通用，类型安全 |
| 素数选择 | 多素数竞选，选模因子数最少的 | NTL/FLINT 验证过的策略 |

### 1.4 模块结构

```
                ┌──────────────┐
                │ M4: Z[x] 分解 │  factorize(polynomial / upolynomial)
                │ factor_ZZ     │
                └──┬────┬────┬──┘
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
```

---

## 2. 公共数据结构

### 2.1 因式分解结果类型

> 代码: `polynomial_factorize.hh:23-27`

```cpp
template<class Poly>
struct factorization {
    typename Poly::coeff_type content;
    std::vector<std::pair<Poly, uint64_t>> factors;
};
```

表示 `f = content · ∏ factors[i].first ^ factors[i].second`。

规范化保证：
- 每个因子本原（primitive），首项系数为正
- 因子按 degree 升序排列
- 重数 > 0
- content 可以是负数（吸收符号）
- 若 f = 0，则 content = 0，factors 为空

### 2.2 与已有类型的兼容

`squarefreefactorize` 返回 `vector<pair<polynomial, uint64_t>>`。
`factorize` 使用 `factorization<Poly>` 结构体包装，增加 content 字段。
两个接口独立，重数类型统一使用 `uint64_t`。

---

## 3. 辅助函数层

### 3.1 单变量 Zₚ 多项式运算

#### `__upoly_make_monic` — 首一化

> 代码: `polynomial_factorize.hh:37-46`

将 f 除以其首项系数，使其首一。返回原 lc。

用途: M1 (DDF/EDF 要求首一输入)

#### `__upoly_mod` — 多项式取模

> 代码: `polynomial_factorize.hh:49-56`

计算 f mod g in Zₚ[x]。复用已有 `pair_vec_div`。

用途: M1 (DDF 中 h^p mod f*), M2 (Hensel 中 Bézout 系数维护)

#### `__upoly_divmod` — 商和余式

> 代码: `polynomial_factorize.hh:59-66`

计算 f = q·g + r in Zₚ[x]。

用途: M2 (Hensel 提升中 divmod(s·e, h))

#### `__upoly_divmod_mod` — ZZ 系数的模除法

> 代码: `polynomial_factorize.hh:581-676`

在 Z_m[x] 中计算 f = q·g + r，所有系数运算都 mod m。
Hensel 提升到高精度时 m = p^{2^j} 远超 uint32_t，无法使用 Zp 类，
因此需要直接在 ZZ 系数上做模运算。实现为长除法，每步系数用 ZZ 模逆。

辅助函数 `__upoly_mod_coeff`（行 563-578）将 ZZ 多项式各系数 mod m。

用途: M2 (Hensel 提升 `__hensel_step`)

#### `__upoly_gcd_Zp` — Zₚ 上 GCD 包装

> 代码: `polynomial_factorize.hh:69-81`

简化的 Zₚ GCD 接口：计算 gcd(a, b) in Zₚ[x]，返回首一结果。
内部调用已有 `__polynomial_GCD`。

用途: M1 (DDF/EDF 中频繁需要 gcd 计算)

#### `__upoly_gcd_extended` — 扩展 GCD

> 代码: `polynomial_factorize.hh:84-125`

扩展欧几里得算法：求 s, t 使得 s·a + t·b = 1 in Zₚ[x]。

算法：标准扩展 Euclidean，追踪系数：

```
r₀ ← a,    r₁ ← b
s₀ ← 1,    s₁ ← 0
t₀ ← 0,    t₁ ← 1
while r₁ ≠ 0:
    q, r₂ ← divmod(r₀, r₁)
    s₂ ← s₀ - q·s₁;  t₂ ← t₀ - q·t₁
    (r₀,s₀,t₀) ← (r₁,s₁,t₁)
    (r₁,s₁,t₁) ← (r₂,s₂,t₂)
// 归一化使 r₀ = 1
c ← r₀.lc⁻¹;  s ← c·s₀;  t ← c·t₀
```

用途: M2 (Hensel 提升需要初始 Bézout 系数)

#### `__upoly_powmod` — 模幂运算

> 代码: `polynomial_factorize.hh:128-154`

计算 base^exp mod modpoly in Zₚ[x]。使用平方-乘法（binary exponentiation）。
exp 使用 ZZ 类型（可能非常大，如 p^d）。

复杂度：O(log(exp) · M(n))，M(n) = O(n²)。

用途: M1 (DDF 中计算 x^{p^d} mod f; EDF 中计算 r^{(p^d-1)/2} mod f)

#### `__upoly_random` — 随机多项式

> 代码: `polynomial_factorize.hh:157-171`

生成度数 < max_deg 的随机 Zₚ[x] 多项式。使用 `std::mt19937`。

用途: M1 (EDF/Cantor-Zassenhaus 中的随机元素选取)

### 3.2 Z 上的对称模约化

#### `__symmetric_mod` — 标量对称模

> 代码: `polynomial_factorize.hh:466-475`

将 a 约化到 (-m/2, m/2] 范围。使用 `ZZ::fdiv_r`（保证非负余数）。

用途: M3 (Zassenhaus 重组中将提升系数映射回 Z)

#### `__upoly_symmetric_mod` — 多项式对称模

> 代码: `polynomial_factorize.hh:478-491`

对多项式每个系数做对称模约化。

用途: M3 (Zassenhaus 重组)

### 3.3 范数计算

> 代码: `polynomial_factorize.hh:498-504` (L2²), `polynomial_factorize.hh:966-976` (L1)

```cpp
inline ZZ __upoly_norm_l2_sq(const upolynomial_<ZZ>& f);  // 系数平方和
inline ZZ __upoly_norm_l1(const upolynomial_<ZZ>& f);     // 系数绝对值之和
```

用途: M2 (Mignotte 界), M3 (Zassenhaus 剪枝参考)

### 3.4 Mignotte 界

> 代码: `polynomial_factorize.hh:548-556`

```cpp
inline ZZ __mignotte_bound(const upolynomial_<ZZ>& f);
```

计算 f 的因子系数绝对值上界：B = C(n, ⌊n/2⌋) · ‖f‖₂。

辅助函数：
- `__binomial(n, k)`（行 511-523）：ZZ 精确二项式系数
- `__isqrt_ceil(n)`（行 526-545）：整数平方根上界（Newton 法）

用途: M2 (确定 Hensel 提升精度 k: 需 p^k > 2·|lc(f)|·B)

### 3.5 本原化

> 代码: `polynomial_factorize.hh:979-987`

```cpp
inline std::pair<ZZ, upolynomial_<ZZ>> __upoly_primitive(upolynomial_<ZZ> f);
```

提取内容并本原化，确保 lc > 0。

用途: M3, M4

### 3.6 类型转换辅助

#### `__upoly_Zp_to_ZZ` — Zp→ZZ 多项式转换

> 代码: `polynomial_factorize.hh:682-689`

将 `upolynomial_<Zp>` 的每个系数从 Zp 转为 ZZ。

用途: M2 (Hensel 树构建时初始因子转换)

#### `__upoly_to_poly` — upolynomial→polynomial 转换

> 代码: `polynomial_factorize.hh:1190-1199`

委托 `poly_convert` 将 `upolynomial_<ZZ>` 转为 `polynomial_<ZZ,lex>`。

用途: M4 (输出结果时将内部表示转回用户类型)

#### `__make_zp` — Zp 构造辅助

> 代码: `polynomial_factorize.hh:30`

构造 `Zp(int64_t, uint32_t)` 避免重载歧义。

### 3.7 其他辅助

#### `__upoly_subtract_x` / `__upoly_subtract_one`

> 代码: `polynomial_factorize.hh:252-312`

从多项式中减去 x 或减去 1 的辅助函数。

用途: M1 (DDF 中构造 h-x, EDF 中构造 g^exp - 1)

#### `__upoly_mul_mod` — ZZ 多项式乘法并 mod m

> 代码: `polynomial_factorize.hh:708-716`

用途: M2 (Hensel 树构建)

#### `__upoly_const_term` — 获取常数项

> 代码: `polynomial_factorize.hh:1008-1013`

用途: M3 (常数项剪枝)

---

## 4. M1：有限域 Zₚ 上的单变量分解

给定 f ∈ Zₚ[x]，返回其完整的不可约因子分解。

### 4.1 函数调用关系

```
__factor_Zp                    M1 顶层入口
  ├── __upoly_make_monic       首一化
  ├── __squarefree_Zp          无平方分解
  │     ├── __upoly_gcd_Zp     GCD 包装
  │     ├── derivative         求导 (已有)
  │     └── __extract_pth_root p 次根提取
  ├── __ddf_Zp                 按度数分组
  │     ├── __upoly_powmod     模幂
  │     └── __upoly_gcd_Zp     GCD 包装
  └── __edf_Zp                 等度分裂 (Cantor-Zassenhaus)
        ├── __upoly_random     随机多项式
        ├── __upoly_powmod     模幂
        └── __upoly_gcd_Zp     GCD 包装
```

### 4.2 `__extract_pth_root` — p 次根提取

> 代码: `polynomial_factorize.hh:178-189`

提取 f(x) = g(x^p) 中的 g：将每个 (x^{kp}, aₖ) 映射为 (x^k, aₖ)。
Zₚ 上 Frobenius 逆是恒等映射，系数不变。

### 4.3 `__squarefree_Zp` — Zₚ 上无平方分解

> 代码: `polynomial_factorize.hh:192-249`

与特征 0 (squarefreefactorize) 不同，特征 p 下 f' 可能为零，
需要处理 f(x) = g(x^p) 的情况。

算法：

```
1. f' ← derivative(f)
2. if f' = 0:  // f = g(x^p)
       g ← __extract_pth_root(f)
       递归 __squarefree_Zp(g)，各重数乘以 p
3. c ← gcd(f, f'), w ← f/c, i ← 1
4. while w ≠ 1:
       y ← gcd(w, c), z ← w/y
       if z ≠ 1: result.push(z, i)
       w ← y, c ← c/y, i++
5. if c ≠ 1:
       g ← __extract_pth_root(c)
       递归 __squarefree_Zp(g)，各重数乘以 p
```

### 4.4 `__ddf_Zp` — 按度数分组 (Distinct-Degree Factorization)

> 代码: `polynomial_factorize.hh:315-360`

将 f 的不可约因子按度数分组：返回 [(g₁,1), (g₂,2), ...] 其中
gd = ∏{不可约 h: deg(h)=d, h|f} h。

算法：

```
h ← x, f* ← f
for d = 1, 2, ... :
    if deg(f*) < 2d: break           // f* 本身不可约
    h ← h^p mod f*                    // h = x^{p^d} mod f*
    gd ← gcd(h - x, f*)
    if gd ≠ 1:
        result.push(gd, d)
        f* ← f*/gd, h ← h mod f*    // 缩减
if deg(f*) > 0:
    result.push(f*, deg(f*))
```

关键优化：h 是累积的（每步 h ← h^p mod f*）；找到因子后 f* 缩小。

### 4.5 `__edf_Zp` — 等度分裂 (Equal-Degree Factorization, Cantor-Zassenhaus)

> 代码: `polynomial_factorize.hh:363-422`

将 f 分裂为全部度数为 d 的不可约因子。概率算法，期望迭代 O(k) 次。

算法：

```
if deg(f) = d: return f              // 已是不可约
repeat:
    r ← random(deg < n)
    if p = 2:
        // 特征 2: trace map T(r) = r + r² + r⁴ + ... + r^{2^{d-1}}
        g ← r; for i=1..d-1: g ← g²+r mod f
        g ← gcd(g, f)
    else:
        // 奇特征: g = gcd(r^{(p^d-1)/2} - 1, f)
        g_pow ← r^{(p^d-1)/2} mod f
        g ← gcd(g_pow - 1, f)
    if 0 < deg(g) < deg(f):
        h ← f/g
        递归 edf(g, d), edf(h, d)
        return
```

### 4.6 `__factor_Zp` — Zₚ 上完整分解

> 代码: `polynomial_factorize.hh:425-459`

顶层入口，串联 squarefree → DDF → EDF。

```
1. lc ← make_monic(f)
2. sqf ← __squarefree_Zp(f)
3. for (sⱼ, eⱼ) in sqf:
       ddf ← __ddf_Zp(sⱼ)
       for (gk, dk) in ddf:
           edf(result, gk, dk, rng)  // 每个因子标记重数 eⱼ
4. 排序，返回 (lc, result)
```

注意: 当被 `__select_prime` 调用时，f mod p 已确认无平方，
squarefree 步骤会直接返回 {(f, 1)}。

---

## 5. M2：Hensel 提升

将 Zₚ 上的因式分解结果提升到 Z_{p^k}，精度足以恢复 Z 上的系数。

### 5.1 函数调用关系

```
__hensel_lift                  M2 顶层入口
  ├── __mignotte_bound         系数界
  ├── __hensel_tree_build      构建二叉树 + 初始 Bézout 系数
  │     └── __upoly_gcd_extended   扩展 GCD
  └── __hensel_step            单步二次提升
        ├── __upoly_divmod_mod ZZ 系数模除法
        └── __upoly_mod_coeff  系数 mod m
```

### 5.2 `__hensel_node` — Hensel 树节点

> 代码: `polynomial_factorize.hh:696-705`

```cpp
struct __hensel_node {
    upolynomial_<ZZ> g;     // 左子树因子之积
    upolynomial_<ZZ> h;     // 右子树因子之积
    upolynomial_<ZZ> s;     // Bézout 系数: s·g + t·h ≡ 1 (mod m)
    upolynomial_<ZZ> t;     // Bézout 系数
    int left;               // 左子节点索引 (-1 = 叶子)
    int right;              // 右子节点索引 (-1 = 叶子)
    int leaf_start;         // 叶子范围起始
    int leaf_end;           // 叶子范围结束
};
```

### 5.3 `__hensel_tree_build` — 构建初始提升树

> 代码: `polynomial_factorize.hh:719-786`

对 r 个因子构建 r-1 个内部节点的平衡二叉树。
每个节点的 Bézout 系数通过 `__upoly_gcd_extended` 在 mod p 下计算，
然后通过 `__upoly_Zp_to_ZZ` 将系数从 Zp 转为 ZZ。

递归实现 `__hensel_tree_build_recursive` 自动分割因子列表。

### 5.4 `__hensel_step` — 单步二次 Hensel 提升

> 代码: `polynomial_factorize.hh:789-882`

将 f ≡ g·h (mod m) 提升到 f ≡ g\*·h\* (mod m²)，
同时更新 Bézout 系数。

算法分两部分：

**第一部分：提升因子**

```
1. e ← (f - g·h) / m, then e mod m
2. se ← s·e, divmod(se, h, m) → (q, r)
3. g ← g + m·(t·e + q·g) mod m²
   h ← h + m·r mod m²
```

**第二部分：提升 Bézout 系数**

```
4. e' ← (1 - s·g - t·h) / m, then e' mod m
5. divmod(s·e', h, m) → (q', r')
   s ← s + m·r' mod m²
   t ← t + m·(t·e' + q'·g) mod m²
```

### 5.5 辅助函数

#### `__hensel_extract_factors` — 从树中提取叶子因子

> 代码: `polynomial_factorize.hh:885-900`

递归遍历 Hensel 树，收集所有叶子节点的 g 和 h。

#### `__hensel_lift_recursive` — 自顶向下递归提升

> 代码: `polynomial_factorize.hh:905-917`

关键：根节点以 f 为目标，子节点以父节点更新后的 g/h 为目标。
不能用扁平循环（所有节点共用 f），否则 >2 因子时结果错误。

### 5.6 `__hensel_lift` — 多因子 Hensel 提升（M2 入口）

> 代码: `polynomial_factorize.hh:920-959`

算法：

```
1. B ← __mignotte_bound(f), target ← 2·|lc(f)|·B
2. 处理首项系数: factors[0] 的每个系数乘以 lc(f) mod p
3. nodes ← __hensel_tree_build(factors_adj, p)
4. m ← p; while m ≤ target: 递归提升, m ← m²
5. 提取叶子因子, return (因子, m)
```

返回的因子系数在 [0, m) 范围，调用方需做 `__upoly_symmetric_mod` 映射回 Z。

---

## 6. M3：因子重组 (Zassenhaus)

给定 f ∈ Z[x] 及其 mod p^k 的因子列表，恢复 Z[x] 上的真正不可约因子。

### 6.1 核心难点

Hensel 提升保持因子数量。但 Z 上的真因子可能需要组合多个模因子：
例 f = x⁴ + 1 在 Z[x] 上不可约，但 mod 2: f ≡ (x+1)⁴。

### 6.2 函数调用关系

```
__factor_recombine                  M3 入口 (Zassenhaus)
  ├── __upoly_symmetric_mod         对称模约化
  ├── __upoly_primitive             本原化
  ├── __upoly_const_term            常数项
  └── __subset_product_mod          子集乘积 mod m
```

### 6.3 `__subset_product_mod` — 子集乘积

> 代码: `polynomial_factorize.hh:990-1005`

计算 lc_f · ∏_{i∈subset} factors[i] (mod m)，经过对称约化。

### 6.4 `__factor_recombine` — Zassenhaus 因子重组

> 代码: `polynomial_factorize.hh:1016-1184`

算法：

```
1. T ← {0,...,r-1} (位掩码), f* ← f
2. for s = 1 to ⌊|T|/2⌋:
3.   for S ⊆ T, |S| = s:    // Gosper's hack 枚举
       // 剪枝 1: lc 整除检查
       // 剪枝 2: 常数项整除检查
       // 完整验证: 计算子集积，试除 lc(f*)·f*
       if r = 0:
           pp(g) → result, f* ← pp(q), T ← T\S
           重新从 s=1 开始
4. 剩余 f* 不可约, 加入 result
5. 排序
```

子集枚举使用 Gosper's hack（限制 r ≤ 64）。

### 6.5 剪枝策略

当前实现使用两种剪枝：
1. **首项系数检查**：lc_prod 是否整除 lc(f*)²
2. **常数项检查**：c_prod 是否整除 lc(f*)·f*(0)

> **注：** L1 范数剪枝已移除——L1 范数次可乘性导致有效分解被错误剪掉。

后续可考虑的增强方向（按优先级）：

| 优化 | 来源 | 说明 |
|---|---|---|
| Newton 迹检查 | PARI/GP, NTL | 从前 1-2 个系数预计算迹，O(|S|) 计算量 |
| 试除时逐系数 bound | PARI/GP | 每算出一个商系数就与 bound 比较 |
| 多素数度数可行性 | FLINT, NTL | 收集多个素数下因子度数的交集 |
| Beauzamy bound | PARI/GP | 比 Mignotte 更紧的系数界 |
| van Hoeij (LLL) | van Hoeij | r 很大时的多项式时间重组（需 LLL 实现） |

---

## 7. M4：单变量 Z[x] 因式分解（集成）

端到端的单变量整系数多项式因式分解，整合 M1-M3 及已有的无平方分解。

### 7.1 函数调用关系

```
factorize(polynomial_<ZZ,comp>)              用户入口 (任意 ordering)
  └── factorize(polynomial_<ZZ,lex>)         lex 特化版本
        ├── cont / __upoly_primitive         内容提取
        ├── squarefreefactorize              无平方分解 (已有)
        └── __factor_squarefree_primitive_ZZ 对单个无平方本原因子分解
              ├── __select_prime             素数选择
              │     ├── polynomial_mod       模约化 (已有)
              │     ├── __upoly_gcd_Zp       无平方检测
              │     ├── __ddf_Zp             DDF
              │     └── __edf_Zp             EDF
              ├── __hensel_lift              M2 入口
              └── __factor_recombine         M3 入口

factorize(upolynomial_<ZZ>)                  upolynomial 直接入口
```

### 7.2 `__select_prime` — 素数选择

> 代码: `polynomial_factorize.hh:1205-1287`

选择使模因子数最少的素数。

```cpp
struct __prime_selection_result {
    uint32_t prime;
    std::vector<upolynomial_<Zp>> factors;
    bool irreducible;
};
```

算法：

```
1. 遍历前 5 个素数（如果因子数 > deg/2 则扩展到 20 个）
2. 跳过 lc(f) mod p = 0 或度数下降的素数
3. 检查 f mod p 是否无平方: gcd(fp, fp') = 1
4. 已知无平方后，直接 DDF + EDF（跳过 __factor_Zp 中冗余的 squarefree）
5. 如果只有 1 个因子 → 立即返回 irreducible
6. 选因子数最少的素数
```

### 7.3 `__factor_squarefree_primitive_ZZ` — 分解无平方本原因子

> 代码: `polynomial_factorize.hh:1294-1310`

```
1. sel ← __select_prime(f)
2. if sel.irreducible: return {f}
3. (lifted, modulus) ← __hensel_lift(f, sel.factors, sel.prime)
4. return __factor_recombine(f, lifted, modulus)
```

### 7.4 `factorize` — 用户入口（lex 特化）

> 代码: `polynomial_factorize.hh:1316-1407`

```
1. 零/常数检查
2. 多变量检查（当前抛异常，待 M5 实现后 dispatch）
3. 转 upolynomial，提取内容，本原化
4. 转回 polynomial 做 squarefreefactorize（方案 A: 复用已有成熟实现）
5. 对每个 deg≥2 的无平方因子调用 __factor_squarefree_primitive_ZZ
6. 排序输出
```

### 7.5 `factorize` — 通用 ordering 包装

> 代码: `polynomial_factorize.hh:1413-1433`

转换到 lex → 调用 lex 特化版 → 转换回原 ordering。
与 `polynomial_GCD` 的 generic wrapper 风格一致。

### 7.6 `factorize` — QQ[x] 入口

> 代码: `polynomial_factorize.hh:1439-1496`

```
1. poly_convert(QQ→ZZ) 乘以 LCD
2. 分解 ZZ 多项式
3. 各因子除以 lc 得首一 QQ 因子, lc^mult 吸收到 content
4. content = content_zz / lcd · ∏ lc(fi)^ei
```

### 7.7 `factorize` — upolynomial 直接入口

> 代码: `polynomial_factorize.hh:1502-1571`

允许用户直接对 `upolynomial_<ZZ>` 做因式分解。内部使用临时变量
`__x("x")` 转为 `polynomial_<ZZ,lex>` 调用 `squarefreefactorize`，
对 deg≥2 因子调用 `__factor_squarefree_primitive_ZZ`。

---

## 8. 完整函数签名索引

### 辅助函数

```cpp
inline Zp __make_zp(int64_t val, uint32_t p);                                    // L30

inline Zp __upoly_make_monic(upolynomial_<Zp>& f);                              // L37
inline upolynomial_<Zp> __upoly_mod(const upolynomial_<Zp>& f,
    const upolynomial_<Zp>& g);                                                  // L49
inline void __upoly_divmod(upolynomial_<Zp>& q, upolynomial_<Zp>& r,
    const upolynomial_<Zp>& f, const upolynomial_<Zp>& g);                      // L59
inline upolynomial_<Zp> __upoly_gcd_Zp(const upolynomial_<Zp>& a,
    const upolynomial_<Zp>& b);                                                  // L69
inline void __upoly_gcd_extended(upolynomial_<Zp>& s, upolynomial_<Zp>& t,
    const upolynomial_<Zp>& a, const upolynomial_<Zp>& b);                      // L84
inline upolynomial_<Zp> __upoly_powmod(const upolynomial_<Zp>& base,
    const ZZ& exp, const upolynomial_<Zp>& modpoly);                            // L128
inline upolynomial_<Zp> __upoly_random(int64_t max_deg,
    uint32_t p, std::mt19937& rng);                                              // L157

inline ZZ __symmetric_mod(const ZZ& a, const ZZ& m);                            // L466
inline upolynomial_<ZZ> __upoly_symmetric_mod(const upolynomial_<ZZ>& f,
    const ZZ& m);                                                                // L478
inline ZZ __upoly_norm_l2_sq(const upolynomial_<ZZ>& f);                        // L498
inline ZZ __upoly_norm_l1(const upolynomial_<ZZ>& f);                           // L966
inline ZZ __mignotte_bound(const upolynomial_<ZZ>& f);                          // L548
inline std::pair<ZZ, upolynomial_<ZZ>> __upoly_primitive(upolynomial_<ZZ> f);    // L979

inline void __upoly_mod_coeff(upolynomial_<ZZ>& f, const ZZ& m);                // L563
inline void __upoly_divmod_mod(upolynomial_<ZZ>& q, upolynomial_<ZZ>& r,
    const upolynomial_<ZZ>& f, const upolynomial_<ZZ>& g, const ZZ& m);         // L581
inline upolynomial_<ZZ> __upoly_Zp_to_ZZ(const upolynomial_<Zp>& f);            // L682
inline upolynomial_<ZZ> __upoly_mul_mod(const upolynomial_<ZZ>& a,
    const upolynomial_<ZZ>& b, const ZZ& m);                                    // L708
inline ZZ __upoly_const_term(const upolynomial_<ZZ>& f);                        // L1008
```

### M1: Zₚ 分解

```cpp
inline upolynomial_<Zp> __extract_pth_root(const upolynomial_<Zp>& f);          // L178
inline std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
    __squarefree_Zp(const upolynomial_<Zp>& f);                                 // L192
inline std::vector<std::pair<upolynomial_<Zp>, uint64_t>>
    __ddf_Zp(const upolynomial_<Zp>& f);                                        // L315
inline void __edf_Zp(std::vector<upolynomial_<Zp>>& result,
    const upolynomial_<Zp>& f, uint64_t d, std::mt19937& rng);                  // L363
inline std::pair<Zp, std::vector<std::pair<upolynomial_<Zp>, uint64_t>>>
    __factor_Zp(upolynomial_<Zp> f);                                             // L425
```

### M2: Hensel 提升

```cpp
struct __hensel_node { ... };                                                     // L696
inline void __hensel_tree_build_recursive(...);                                   // L719
inline std::vector<__hensel_node> __hensel_tree_build(
    const std::vector<upolynomial_<Zp>>& factors, uint32_t p);                  // L777
inline void __hensel_step(__hensel_node& node,
    const upolynomial_<ZZ>& f, const ZZ& m);                                    // L789
inline void __hensel_extract_factors(...);                                        // L885
inline void __hensel_lift_recursive(...);                                         // L905
inline std::pair<std::vector<upolynomial_<ZZ>>, ZZ> __hensel_lift(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<Zp>>& factors, uint32_t p);                  // L920
```

### M3: 因子重组

```cpp
inline upolynomial_<ZZ> __subset_product_mod(
    const std::vector<upolynomial_<ZZ>>& factors,
    const std::vector<size_t>& subset, const ZZ& lc_f, const ZZ& m);            // L990
inline std::vector<upolynomial_<ZZ>> __factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted, const ZZ& m);                  // L1016
```

### M4: 单变量集成

```cpp
struct __prime_selection_result { ... };                                           // L1205
inline __prime_selection_result __select_prime(const upolynomial_<ZZ>& f);       // L1211
inline std::vector<upolynomial_<ZZ>>
    __factor_squarefree_primitive_ZZ(const upolynomial_<ZZ>& f);                 // L1294

template<class var_order>
factorization<polynomial_<ZZ,lex_<var_order>>>
    factorize(const polynomial_<ZZ,lex_<var_order>>& f);                         // L1316
template<class comp>
factorization<polynomial_<ZZ,comp>>
    factorize(const polynomial_<ZZ,comp>& f);                                    // L1413
template<class comp>
factorization<polynomial_<QQ,comp>>
    factorize(const polynomial_<QQ,comp>& f);                                    // L1439
factorization<upolynomial_<ZZ>>
    factorize(const upolynomial_<ZZ>& F);                                        // L1502
```

---

## 附录 A: 已有函数依赖清单

| 函数 | 文件 | 用途 |
|---|---|---|
| `polynomial_mod(f, p)` | `polynomial.hh`, `upolynomial.hh` | ZZ→Zp 系数约化 |
| `__polynomial_GCD(Pout, G, F, Lc, deg)` | `polynomial_gcd.hh` | Zₚ 上 Euclidean GCD |
| `squarefreefactorize(F)` | `polynomial_gcd.hh` | Z 上无平方分解 |
| `cont(F)` | `polynomial_gcd.hh`, `upolynomial.hh` | 内容提取 |
| `derivative(f)` | `upolynomial.hh` | 求导 |
| `pair_vec_div(q, r, f, g, comp)` | `basic.hh` | 多项式除法 |
| `poly_convert(in, out)` | `upolynomial.hh` | polynomial ↔ upolynomial 转换 |
| `get_variables(f)` | `polynomial_.hh` | 获取变量列表 |
| `is_number(f)` | `upolynomial.hh` | 常数检测 |
| `boost::math::prime(n)` | `<boost/math/...>` | 第 n 个素数 |

## 附录 B: 参考文献

- **Cantor-Zassenhaus**: Cantor & Zassenhaus, "A New Algorithm for Factoring Polynomials Over Finite Fields", Math. Comp. 1981
- **Hensel Lifting**: Zassenhaus, "On Hensel Factorization I", J. Number Theory 1969
- **van Hoeij**: van Hoeij, "Factoring Polynomials and the Knapsack Problem", J. Number Theory 2002
- **Mignotte Bound**: Mignotte, "An Inequality About Factors of Polynomials", Math. Comp. 1974
- **GCL**: Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992 (§8, §15, §16)
- **MCA**: von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013 (§14-§16)

## 附录 C: 各系统因子表示对比

| 系统 | 内容(content) | 因子列表 | 规范化 |
|---|---|---|---|
| FLINT | `fmpz_t c` 在结构体内 | `fmpz_poly_struct *p` + `slong *exp` | 因子本原，lc > 0 |
| NTL | 单独 `ZZ& c` 参数 | `vec_pair_ZZX_long` | 因子首一 (monic) |
| Mathematica | `FactorList` 第一项 | `{{f₁,e₁},...}` | 因子首项正 |
| Singular | `CFFList` 含 unit | `CFFactor` 列表 | 取决于域 |
| **CLPoly** | `factorization::content` | `vector<pair<Poly,uint64_t>>` | 因子本原，lc > 0 |
