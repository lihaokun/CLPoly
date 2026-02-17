# CLPoly 多变量因式分解设计方案

> **状态：设计阶段，尚未实现。**
>
> 本文档描述 CLPoly 多变量因式分解（M5）的设计方案，基于 Wang 算法。
> 单变量因式分解（M1-M4）已实现，详见 [univariate-factorization.md](univariate-factorization.md)。

---

## 1. 概述

多变量因式分解通过将问题约化为单变量来解决。核心思想：

1. 取值 x₂=α₂, ..., xₙ=αₙ 得到单变量像
2. 对像进行单变量分解（M4，已实现）
3. 通过多变量 Hensel 提升恢复各变量的贡献

### 1.1 设计选择

| 决策 | 选择 | 理由 |
|---|---|---|
| 多变量策略 | Wang 算法 | 经典成熟，与 CLPoly 已有的 Hensel 基础设施兼容 |
| 主变量选择 | lex 首变量 (非最高度数) | 匹配 `cont()` 实现，避免变量重排复杂度 |
| Wang lc 校正返回值 | `__wang_lc_result` 结构体 | 需返回成功标志 + 缩放后的 f + lc 分配列表 |
| Hensel 提升方式 | 逐变量线性提升 | 标准 Wang 方式，每变量提升到 `deg(f,xₖ)` 阶 |
| 失败处理 | 换求值点重试 (最多 10 次) | 简单可靠；远期可加 EEZ-Wang |

### 1.2 主流系统参考

| | Maple | FLINT | Singular |
|---|---|---|---|
| **多变量策略** | MTSHL | fmpz_mpoly_factor | Wang + 稀疏 Hensel |

---

## 2. 函数总览

```
__factor_multivar(polynomial_<ZZ,lex>)     M5 入口
  ├── cont(f, x₁)                         内容提取 (已有)
  ├── factorize(cont)                      递归分解内容 (M4, 已实现)
  ├── __select_eval_point                  选取值点
  │     ├── assign(f, v, c)               代入求值 (已有)
  │     └── is_squarefree                 无平方检测 (已有)
  ├── factorize(f₀)                        单变量分解 (M4, 已实现)
  ├── __wang_leading_coeff                 首项系数分配
  └── __multivar_hensel_lift               多变量 Hensel 提升
        └── __upoly_gcd_extended           扩展 GCD (已实现)
```

---

## 3. 需要新增的辅助函数

### 3.1 `pp` — 多变量本原部分

```cpp
// 多变量本原部分: f / cont(f)
// 前置: f ∈ Z[x₁,...,xₙ], lex 排序
// 后置: 返回 f 除以其关于首变量的内容后的本原部分
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
pp(const polynomial_<ZZ, lex_<var_order>>& f);
```

实现方式：`pp(f) = f / cont(f)`，使用多项式除法 `pair_vec_div`。

> **注：** `cont()` 已存在（`polynomial_gcd.hh:468`），只需补充对应的 `pp()`。

### 3.2 `__poly_coeff_l1_norm` — 多变量系数 L1 范数

```cpp
// 多变量多项式系数绝对值之和
// 用于试除验证时的快速剪枝
template<class var_order>
ZZ __poly_coeff_l1_norm(const polynomial_<ZZ, lex_<var_order>>& f);
```

---

## 4. `__select_eval_point` — 选取值点

```cpp
// 选择多变量 → 单变量的取值点
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2, x₁ 是主变量
// 后置: 返回 α = {x₂→α₂, ..., xₙ→αₙ} 使得:
//       (a) f(x₁, α₂, ..., αₙ) 无平方
//       (b) lc(f, x₁) |_{xᵢ=αᵢ} ≠ 0
//       (c) deg(f(x₁,...)) = deg(f(x₁, α₂, ...)) (首项不消失)
// 策略: 从小整数开始尝试 (0, ±1, ±2, ...)
template<class var_order>
std::map<variable, ZZ>
__select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& main_var);
```

---

## 5. Wang 算法首项系数校正

### 5.1 问题背景

设 f ∈ Z[x₁,...,xₙ] 本原（关于 x₁），其首项系数 `lc(f, x₁) ∈ Z[x₂,...,xₙ]`
可能是非平凡的多变量多项式。设取值后 `f₀ = f(x₁, α₂,...,αₙ)` 在 Z[x₁] 上分解为
`f₀ = c · g₁ · g₂ · ··· · gᵣ`。Hensel 提升需要知道每个因子 gᵢ 对应的多变量首项系数，
否则无法正确恢复多变量因子。

### 5.2 算法

```
__wang_leading_coeff(f, univar_factors, eval_point, main_var):

1.  L ← lc(f, x₁)                          // ∈ Z[x₂,...,xₙ]
    if is_number(L):
        return                               // 首项系数为常数，无需校正

2.  // 递归分解首项系数
    lc_factors ← factorize(L)               // 递归调用多变量 factorize
    // lc_factors = {δ, [(l₁,e₁), (l₂,e₂), ...]}
    // 其中每个 lⱼ ∈ Z[x₂,...,xₙ] 是不可约的

3.  // 在求值点处计算各 lc 因子的值
    for each (lⱼ, eⱼ):
        vⱼ ← assign(lⱼ, eval_point) ^ eⱼ   // lⱼ(α₂,...,αₙ)^eⱼ ∈ Z

4.  // 将 lc 因子分配到单变量因子
    // 原则: lc(gᵢ) = 某些 lⱼ^eⱼ 在 eval_point 处的乘积
    for i = 1 to r:
        σᵢ ← 1                              // 分配给 gᵢ 的 lc 多项式 ∈ Z[x₂,...,xₙ]
    for each (lⱼ, eⱼ):
        // 找到 lc(gᵢ) 能被 vⱼ 整除的 gᵢ
        found ← false
        for i = 1 to r:
            if lc(gᵢ) mod vⱼ == 0:
                σᵢ ← σᵢ · lⱼ^eⱼ
                lc(gᵢ) ← lc(gᵢ) / vⱼ
                found ← true
                break
        if !found:
            // 分配失败 → 换求值点重试
            return FAIL

5.  // 缩放: 令 δ ← lc_factors.content (整数部分)
    // 将 f 乘以 δ^(r-1) 使提升可行
    // 修改各 univar_factor 的首项系数为 σᵢ(α₂,...,αₙ)
    for i = 1 to r:
        gᵢ ← (σᵢ(α₂,...,αₙ) / lc(gᵢ)) · gᵢ
```

### 5.3 函数签名

```cpp
// 首项系数校正结果
template<class var_order>
struct __wang_lc_result {
    bool success;
    polynomial_<ZZ, lex_<var_order>> f_scaled;  // 缩放后的 f
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments; // 每个因子的 lc
};

// 首项系数校正
// 前置: f 关于 x₁ 本原, univar_factors 是 f(x₁,α...) 的因子
// 后置: 返回 true 并填充 lc_assignments（每个因子对应的多变量 lc）
//       返回 false 表示分配失败（需换求值点）
// 修改 univar_factors 使其首项系数与 lc 分配一致
// 修改 f_scaled 为缩放后的 f（乘以适当的 δ^(r-1)）
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

---

## 6. 多变量 Hensel 提升

### 6.1 算法框架

逐变量提升：将 `f(x₁, α₂,...,αₙ)` 的因子逐步恢复为 `f(x₁, x₂, α₃,...,αₙ)` 的因子，
再恢复为 `f(x₁, x₂, x₃, α₄,...,αₙ)` 的因子，以此类推。

对每个变量 `xₖ`（k = 2, 3, ..., n），做如下提升：

```
设 f* 为当前待分解的多项式（可能已被 lc 缩放）
设 g₁,...,gᵣ 为当前因子（满足 g₁·...·gᵣ ≡ f* mod (xₖ-αₖ)^j）
设 dₖ = deg(f, xₖ) 为变量 xₖ 的度数上界

for j = 1 to dₖ:
    // 计算误差
    e ← f* - g₁ · g₂ · ··· · gᵣ
    // 提取 (xₖ - αₖ)^j 的系数
    eⱼ ← coeff(e, (xₖ - αₖ)^j)

    if eⱼ == 0:
        continue

    // 将误差分配到各因子（使用 Bézout 系数）
    // 需要: s₁·ĝ₁ + s₂·ĝ₂ + ... + sᵣ·ĝᵣ ≡ 1 mod I
    // 其中 ĝᵢ = ∏_{j≠i} gⱼ mod (xₖ-αₖ)
    for i = 1 to r:
        δᵢ ← sᵢ · eⱼ mod gᵢ|_{xₖ=αₖ}   // 在 Z[x₁] 上做模运算
        gᵢ ← gᵢ + δᵢ · (xₖ - αₖ)^j
```

### 6.2 终止条件

每个变量 `xₖ` 的提升精度为 `deg(f, xₖ)`，因为 f 的任何因子在 `xₖ` 上的度数
不超过 `deg(f, xₖ)`。

### 6.3 函数签名

```cpp
// 多变量 Hensel 提升
// 前置: univar_factors 是 f 在 eval_point 处的因子
//       lc_assignments 是 __wang_leading_coeff 计算的 lc 分配
// 后置: 返回 f 在 Z[x₁,...,xₙ] 上的候选因子（需试除验证）
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>>
__multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_assignments,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

---

## 7. `__factor_multivar` — 多变量分解入口

```cpp
// 多变量因式分解 (Wang 算法)
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2
// 后置: 返回 factorization
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f);
```

### 7.1 完整算法

```
__factor_multivar(f):

1.  // 选主变量: 使用 lex 序的首变量
    // (避免变量重排; cont() 依赖 lex 首变量)
    x₁ ← get_first_var(f)
    c ← cont(f)                   // cont 关于 x₁ (lex 首变量)
    f ← f / c                     // 本原部分 (用多项式除法)

2.  // 递归分解内容
    cont_factors ← {}
    if !is_number(c):
        cont_factors ← factorize(c)    // 递归

3.  // 选取值点
    eval ← __select_eval_point(f, x₁)

4.  // 单变量分解
    f₀ ← assign(f, eval)               // f₀ ∈ Z[x₁]
    uni_result ← factorize(f₀)

5.  if uni_result.factors.size() ≤ 1:
        // 单变量像不可约 ⇒ f 不可约
        // (f 本原 + 求值点满足条件 ⇒ 理论保证)
        return {1, [(f, 1)]} ∪ cont_factors

6.  // 首项系数校正 (Wang 核心步骤)
    lc_result ← __wang_leading_coeff(f, uni_result.factors, eval, x₁)
    if !lc_result.success:
        // 换求值点重试 (最多 MAX_RETRY 次)
        goto 3

7.  // 多变量 Hensel 提升
    mv_factors ← __multivar_hensel_lift(
        lc_result.f_scaled, uni_result.factors,
        lc_result.lc_assignments, eval, x₁)

8.  // 试除验证 + 去缩放
    verified ← []
    f* ← lc_result.f_scaled
    for g in mv_factors:
        g_prim ← pp(g, x₁)            // 去缩放: 取本原部分
        if g_prim | f:                  // 注意: 对原始 f 试除，非 f_scaled
            verified.push(g_prim)
            f ← f / g_prim
    if deg(f) > 0:
        verified.push(f)

9.  // 合并内容因子
    return combine(cont_factors, verified)
```

---

## 8. 失败处理与重试策略

Wang 算法可能在以下情况失败：

| 失败点 | 原因 | 处理 |
|---|---|---|
| `__select_eval_point` | 小整数范围内无合法点 | 扩大搜索范围 (\|αᵢ\| > 100) |
| `__wang_leading_coeff` | lc 因子分配不唯一 | 换求值点重试 |
| 试除全部失败 | 提升精度不足或数值问题 | 换求值点重试 |

最大重试次数建议 `MAX_RETRY = 10`。超过后抛出异常。

> **注：** 远期可考虑实现 EEZ-Wang (Extended Zassenhaus for Wang) 变体，
> 在 lc 分配困难时使用"延迟 lc 分配"策略，但初始实现不需要。

---

## 9. `factorize` 入口集成

当前 `factorize`（`polynomial_factorize.hh:1340-1341`）在 `vars.size() > 1` 时
抛异常。M5 完成后，应修改为：

```cpp
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
factorize(const polynomial_<ZZ, lex_<var_order>>& F)
{
    ...
    auto vars = get_variables(F);
    if (vars.size() > 1)
        return __factor_multivar(F);     // ← 新增: dispatch 到多变量
    ...
    // 原有单变量逻辑
}
```

QQ[x₁,...,xₙ] 入口无需修改——其内部先转换为 ZZ 多项式再调用
`factorize(ZZ)`，转换逻辑对多变量同样有效。

### 9.1 与 `squarefreefactorize` 的兼容性

现有 `squarefreefactorize`（`polynomial_gcd.hh:101`）对 lex 首变量求导，
递归处理内容。对多变量多项式，它能正确检测**关于首变量**的重因子。

**已知限制：** 如果 f 关于非首变量有重因子，`squarefreefactorize` 可能将其视为无平方。
但这对 Wang 算法不构成问题——Wang 的前置条件只要求**单变量像** `f(x₁, α...)` 无平方，
这在 `__select_eval_point` 中已保证。

`__factor_multivar` 仍应在开头调用 `squarefreefactorize`，因为它能检测
首变量方向的重因子并降低后续提升的规模。

---

## 10. 实现路线图

| 阶段 | 内容 | 新增函数 | 依赖 |
|---|---|---|---|
| **Phase 5** | M5: 多变量 Wang | `pp`, `__select_eval_point`, `__wang_leading_coeff` (含 `__wang_lc_result`), `__multivar_hensel_lift`, `__factor_multivar`, `factorize` 多变量 dispatch | M4 (已实现) |
| **Phase 6** | 增强：van Hoeij 重组 | `__factor_recombine_van_hoeij` + LLL 实现 | M3 替换 |

### 10.1 测试计划

| 可独立测试的函数 | 验证方法 |
|---|---|
| `pp(f)` | 验证 `cont(f) · pp(f) == f` |
| `__select_eval_point` | 验证返回点满足三个条件 (无平方, lc 非零, 度数不降) |
| `__wang_leading_coeff` | 构造已知分解的多项式，验证 lc 分配正确 |
| `__multivar_hensel_lift` | 对二变量多项式，验证提升结果在试除后正确 |
| `__factor_multivar` | 与 Mathematica `Factor[f]` 对比 |
| `factorize` (多变量入口) | 与 Mathematica `FactorList[f]` 对比 |

### 10.2 测试用例

```
// 简单二变量
f = x² - y²                         // = (x-y)(x+y)
f = x² + 2xy + y²                   // = (x+y)²
f = x³ + y³                         // = (x+y)(x²-xy+y²)

// lc 非平凡
f = (2y+1)x² + 3x + y              // lc(f,x) = 2y+1

// 三变量
f = x²-y²-z²+1                     // 检查不可约
f = (x+y+z)(x-y+z)                 // 简单可分解

// 稀疏
f = x¹⁰ + y¹⁰ - 1                  // 高度数稀疏

// 含内容
f = (y+1)(x²-y²)                   // cont = y+1

// 含重因子
f = (x+y)²(x-y)                    // 通过 squarefreefactorize 先拆分
```

---

## 附录 A: 完整函数签名索引

```cpp
// §3.1 多变量本原部分
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
pp(const polynomial_<ZZ, lex_<var_order>>& f);

// §3.2 多变量系数 L1 范数
template<class var_order>
ZZ __poly_coeff_l1_norm(const polynomial_<ZZ, lex_<var_order>>& f);

// §4 选取值点
template<class var_order>
std::map<variable, ZZ> __select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f, const variable& main_var);

// §5 首项系数校正结果
template<class var_order>
struct __wang_lc_result {
    bool success;
    polynomial_<ZZ, lex_<var_order>> f_scaled;
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments;
};

// §5 首项系数校正
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

// §6 多变量 Hensel 提升
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>> __multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_assignments,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

// §7 多变量分解入口
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f);
```

## 附录 B: 参考文献

- **Wang**: Wang, "An Improved Multivariate Polynomial Factoring Algorithm", Math. Comp. 1978
- **Kaltofen-Shoup**: Kaltofen & Shoup, "Subquadratic-Time Factoring of Polynomials over Finite Fields", Math. Comp. 1998
- **GCL**: Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992 (§16)
- **MCA**: von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013 (§16)

## 附录 C: 已有函数依赖清单

以下是 M5 直接依赖的已有函数：

| 函数 | 文件 | 用途 |
|---|---|---|
| `factorize(polynomial_<ZZ>)` | `polynomial_factorize.hh` | 单变量分解（M4，已实现） |
| `squarefreefactorize(F)` | `polynomial_gcd.hh` | 无平方分解 |
| `cont(F)` | `polynomial_gcd.hh` | 内容提取 |
| `assign(f, v, c)` | `polynomial.hh` | 变量代入 |
| `is_squarefree(f)` | `polynomial_gcd.hh` | 无平方检测 |
| `get_variables(f)` | `polynomial_.hh` | 获取变量列表 |
| `pair_vec_div(q, r, f, g, comp)` | `basic.hh` | 多项式除法 |
| `__upoly_gcd_extended(s, t, a, b)` | `polynomial_factorize.hh` | 扩展 GCD（已实现） |
