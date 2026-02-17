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

### 1.2 范围与正确性保证

**支持范围：**
- `Z[x₁,...,xₙ]` 和 `Q[x₁,...,xₙ]`（QQ 通过 LCD 转 ZZ 后调用）
- 变量数 n ≥ 2（n = 1 由 M4 处理，通过 `factorize` 入口自动 dispatch）
- 不支持有限域 `Zp[x₁,...,xₙ]`，不支持代数扩张

**正确性不变量：**

对非零输入 `f ∈ Z[x₁,...,xₙ]`，输出满足：
```
f = content · ∏ fᵢ^eᵢ
```
其中每个 `fᵢ` 在 `Z[x₁,...,xₙ]` 上**不可约**，本原（`cont(fᵢ) = 1`），`lc(fᵢ) > 0`。
因子按 `(degree, 字典序)` 升序排列，重数 `eᵢ > 0`。

**性能预期（初始实现）：**
- 目标：总度数 ≤ 50，变量数 ≤ 6 可在秒级完成
- 瓶颈在多变量 Hensel 提升（逐变量线性）和 LC 分配
- 不追求大规模性能，留给远期 MTSHL

**已知限制：**
- LC 分配采用贪心策略，特定 lc 结构可能需多次换求值点
- 无 EEZ-Wang 后备
- 主变量固定为 lex 首变量，不做度数最优选择

### 1.3 主流系统参考

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

### 2.1 数据流与缩放不变量

Wang 算法中 `f` 经历多次变换，各模块操作的对象不同，必须严格区分：

```
f_input                          用户输入
  │
  ├─ c = cont(f, x₁)            内容 ∈ Z[x₂,...,xₙ]
  ▼
f_prim = pp(f, x₁)              本原部分，后续所有操作基于此
  │
  ├─ f₀ = f_prim(x₁,α)         单变量像 ∈ Z[x₁]
  ├─ u₁,...,uᵣ = factorize(f₀)  单变量因子（首一）
  ▼
__wang_leading_coeff
  ├─ 输入:  f_prim, u₁,...,uᵣ, α
  ├─ 输出:  f_scaled = δ^(r-1) · f_prim       ← 缩放后的多项式
  │         σ₁,...,σᵣ ∈ Z[x₂,...,xₙ]          ← 各因子的 lc 分配
  │         v₁,...,vᵣ ∈ Z[x₁]                  ← 修改后的单变量因子
  │           (vᵢ 首项系数 = σᵢ(α), 满足 ∏vᵢ = f_scaled(x₁,α))
  ▼
__multivar_hensel_lift
  ├─ 输入:  f_scaled, v₁,...,vᵣ, σ₁,...,σᵣ, α
  ├─ 不变量: ∏ Gᵢ ≡ f_scaled (mod 提升理想)
  ├─ 输出:  G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]          ← 候选因子
  ▼
试除验证
  ├─ 对 f_prim（非 f_scaled！）做试除
  ├─ pp(Gᵢ) | f_prim ?
  └─ 输出: 不可约因子列表
```

**关键不变量：**
- `f_scaled = δ^(r-1) · f_prim`，其中 `δ = lc(f_prim, x₁)(α)`
- Hensel 提升在 `f_scaled` 上进行，因为缩放保证各因子首项系数为多变量多项式
- 试除在 `f_prim` 上进行：`pp(Gᵢ)` 消去了缩放因子 `δ`，结果整除 `f_prim`
- 这与单变量 M3 的 `pp(g) | f*` 试除逻辑一致（见 univariate §6.4）

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

### 4.1 求值点条件

选择 `α = {x₂→α₂, ..., xₙ→αₙ}` 满足以下全部条件：

| 条件 | 要求 | 保证 |
|---|---|---|
| (a) 无平方 | `f(x₁, α)` 在 Z[x₁] 中无平方 | Hensel 提升需要因子两两互素 |
| (b) lc 非零 | `lc(f, x₁)(α) ≠ 0` | 首项不消失，deg 不降 |
| (c) 度数守恒 | `deg(f(x₁,α)) = deg(f, x₁)` | 由 (b) 蕴含 |
| (d) lc 因子可分辨 | `lc(f, x₁)` 的不可约因子在 α 处的值**两两互素** | LC 分配的贪心算法能正确工作 |

> **条件 (d) 的重要性：** 若两个不同的 lc 因子 lⱼ、lₖ 在 α 处的值有公因子，
> 则 §5.2 的贪心分配可能将 lⱼ 错误分配给 uᵢ，导致分配失败。
> 实现时，先对 `lc(f, x₁)` 做因式分解（可递归调用 `factorize`），
> 再检查各因子求值后两两 `gcd = 1`。

### 4.2 搜索策略

从小整数开始尝试 `αᵢ ∈ {0, 1, -1, 2, -2, ...}`。
每个候选点检查条件 (a)-(d)，找到第一个满足的即返回。

**保证终止：** Z 是无限域，不满足条件的 α 构成一个代数集合（零测集），
因此有限步内必然找到合法点。实践中几乎总在 |αᵢ| ≤ 5 内找到。

### 4.3 函数签名

```cpp
// 选择多变量 → 单变量的取值点
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2, x₁ 是主变量, f 关于 x₁ 本原
// 后置: 返回 α 满足条件 (a)-(d)
// 注: 保证终止（Z 无限域），实践中搜索范围很小
template<class var_order>
std::map<variable, ZZ>
__select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& main_var);
```

---

## 5. Wang 算法首项系数校正

### 5.1 问题背景

设 `f ∈ Z[x₁,...,xₙ]` 本原（关于 x₁），其首项系数 `L = lc(f, x₁) ∈ Z[x₂,...,xₙ]`
可能是非平凡的多变量多项式。设取值后 `f₀ = f(x₁, α)` 在 `Z[x₁]` 上分解为
`f₀ = c · u₁ · u₂ · ··· · uᵣ`（uᵢ 首一）。

**为什么需要 LC 校正：** Hensel 提升要求各因子的首项系数（在 x₁ 上）已知。
如果 L 是非平凡多项式，则 f 的各多变量因子的 lc 是 L 的某种"分拆"。
不提前分配 lc，提升过程无法确定每步应将多少首项系数"还给"各因子。

### 5.2 算法

**前置条件：** 求值点 α 满足条件 (d)：L 的不可约因子在 α 处两两互素。

```
__wang_leading_coeff(f, u₁,...,uᵣ, α, x₁):

1.  L ← lc(f, x₁)                          // ∈ Z[x₂,...,xₙ]
    δ ← L(α)                               // ∈ Z (非零，由条件 (b) 保证)
    if is_number(L):
        // 首项系数为常数，无需校正，也无需缩放
        σᵢ ← L  for all i (仅 σ₁=L, 其余 σᵢ=1 亦可)
        return SUCCESS, f_scaled = f, σ₁,...,σᵣ

2.  // 递归分解首项系数
    lc_fac ← factorize(L)                  // 递归调用 factorize
    // lc_fac = {γ, [(l₁,e₁), (l₂,e₂), ...]}
    // γ ∈ Z 是整数内容，每个 lⱼ ∈ Z[x₂,...,xₙ] 不可约本原

3.  // 在求值点处计算各 lc 因子的值
    for each (lⱼ, eⱼ):
        vⱼ ← lⱼ(α) ^ eⱼ                   // ∈ Z

4.  // 初始化: 每个单变量因子分配到的 lc 多项式
    lcᵢ ← lc(uᵢ) · (δ / ∏ lc(uⱼ))^...    // 先记录各因子的数值 lc
    // 实际做法: 逐个分配 lⱼ^eⱼ

5.  // 将 lc 因子分配到单变量因子
    for i = 1 to r:
        σᵢ ← 1                              // 累积分配给 uᵢ 的 lc ∈ Z[x₂,...,xₙ]
        wᵢ ← lc(uᵢ)                         // 待分配的数值部分

    // 按 eⱼ 从大到小排序（高次幂优先分配，减少歧义）
    sort (lⱼ, eⱼ) by eⱼ descending

    for each (lⱼ, eⱼ):
        vⱼ ← lⱼ(α) ^ eⱼ                   // 该因子求值后的值
        // 找唯一的 uᵢ 使得 vⱼ | wᵢ
        candidates ← {i : vⱼ | wᵢ}
        if |candidates| = 0:
            return FAIL                      // 分配失败
        if |candidates| > 1:
            return FAIL                      // 不唯一（条件 (d) 应已排除此情况，
                                             //   但仍需防御性检查）
        i ← candidates 中唯一的元素
        σᵢ ← σᵢ · lⱼ^eⱼ
        wᵢ ← wᵢ / vⱼ

    // 分配完成后，wᵢ 应为 ±1 · γ 的某种分配
    // 将剩余整数部分 (γ 和各 wᵢ) 乘入 σ 和 δ
    // 简化处理: 将 γ 乘入 σ₁

6.  // 缩放
    f_scaled ← δ^(r-1) · f                  // 乘 lc(f,x₁)(α)^(r-1)
    for i = 1 to r:
        vᵢ ← (σᵢ(α) / lc(uᵢ)) · uᵢ       // 替换 lc: 现在 lc(vᵢ) = σᵢ(α)
    // 验证: ∏ vᵢ = f_scaled(x₁, α)
    return SUCCESS, f_scaled, σ₁,...,σᵣ, v₁,...,vᵣ
```

> **与单变量 M3 的类比：** 单变量 Hensel 提升也将 `lc(f)` 乘到 `factors[0]`。
> 多变量版本的缩放 `δ^(r-1)` 起类似作用，但更系统化——每个因子有明确的 lc 分配。
> 试除时同样用 `pp()` 消去缩放因子（见 §2.1 数据流）。

### 5.3 正确性条件

| 条件 | 保证 |
|---|---|
| `L(α) ≠ 0` | 由 eval_point 条件 (b) |
| 各 lⱼ(α) 两两互素 | 由 eval_point 条件 (d)，保证唯一分配 |
| `∏ σᵢ(α) = δ` | 分配完整性，所有 lc 因子被完全消耗 |
| `∏ vᵢ = f_scaled(x₁, α)` | 缩放正确性，Hensel 提升的前置条件 |

### 5.4 函数签名

```cpp
// 首项系数校正结果
template<class var_order>
struct __wang_lc_result {
    bool success;
    polynomial_<ZZ, lex_<var_order>> f_scaled;          // δ^(r-1) · f
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments;  // σ₁,...,σᵣ
    std::vector<upolynomial_<ZZ>> scaled_factors;       // v₁,...,vᵣ (修改后的单变量因子)
};

// 首项系数校正
// 前置: f 关于 x₁ 本原, univar_factors = [u₁,...,uᵣ] 是 f(x₁,α) 的首一因子
//       eval_point 满足条件 (a)-(d)
// 后置: success=true 时，f_scaled/lc_assignments/scaled_factors 已填充
//       success=false 时，需换求值点重试
// 不修改输入参数（所有输出通过返回值）
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

> **接口设计说明：** 对比之前版本，`univar_factors` 改为 `const &`，
> 修改后的因子通过 `scaled_factors` 返回。避免副作用 + 返回值混用的问题。

---

## 6. 多变量 Hensel 提升

### 6.1 总体结构

逐变量提升：将 `f_scaled(x₁, α₂,...,αₙ)` 的因子逐步恢复为
`f_scaled(x₁, x₂, α₃,...,αₙ)` 的因子，再恢复为
`f_scaled(x₁, x₂, x₃, α₄,...,αₙ)` 的因子，以此类推。

```
for k = 2 to n:
    f_curr ← f_scaled(x₁,...,xₖ, αₖ₊₁,...,αₙ)    // 部分求值
    dₖ ← deg(f_curr, xₖ)
    __hensel_lift_one_var(f_curr, G₁,...,Gᵣ, ĝ₁,...,ĝᵣ, s₁,...,sᵣ, σ₁,...,σᵣ, xₖ, αₖ, dₖ)
    // 每步后 Gᵢ 从 Z[x₁,...,xₖ₋₁] 扩展为 Z[x₁,...,xₖ]
```

### 6.2 Bézout 系数：计算一次，全局复用

**关键：** Bézout 系数基于**完全求值后的单变量因子**计算，只算一次。

设 `ĝᵢ = uᵢ(x₁)` 为原始单变量因子（求值点 α 处），在 Z[x₁] 中两两互素。
Bézout 系数 `s₁,...,sᵣ ∈ Z[x₁]` 满足偏分式恒等式：

```
Σᵢ sᵢ · (∏_{j≠i} ĝⱼ) = 1    in Z[x₁]
```

其中 `deg(sᵢ) < deg(ĝᵢ)`。此恒等式精确成立（非模某个理想），
因为 ĝᵢ 两两互素。

**计算方法（逐对 XGCD 链，与 Singular 一致）：**

```
g_acc ← ĝ₁
s[1] ← 1
for i = 2 to r:
    (α, β) ← XGCD(g_acc, ĝᵢ)      // α·g_acc + β·ĝᵢ = 1 in Z[x₁]
    for j = 1 to i-1:
        s[j] ← s[j] · α mod ĝⱼ     // 已有系数乘以 α
    s[i] ← β mod ĝᵢ
    g_acc ← g_acc · ĝᵢ
```

> **注：** 当前 `__upoly_gcd_extended` 实现在 Zₚ[x] 上。
> 多变量 Hensel 中需要在 **Z[x₁]** 上做 XGCD（系数是整数，非 Zp）。
> 实现时需要补充 ZZ 版本的 XGCD，或先对 lc 做模逆后在 Q[x₁] 上操作。

这些 sᵢ 和 ĝᵢ 在所有变量的提升过程中保持不变——它们只依赖求值点处的单变量因子。

### 6.3 单变量提升步（核心）

对变量 `xₖ` 的提升。此时：
- `Gᵢ ∈ Z[x₁,...,xₖ₋₁]`（已经过前面变量的提升）
- 目标：将 `Gᵢ` 扩展为 `Z[x₁,...,xₖ]` 上的多项式
- `f_curr = f_scaled(x₁,...,xₖ, αₖ₊₁,...,αₙ) ∈ Z[x₁,...,xₖ]`

**算法：**

```
__hensel_lift_one_var(f_curr, G₁,...,Gᵣ, ĝ₁,...,ĝᵣ, s₁,...,sᵣ, σ₁,...,σᵣ, xₖ, αₖ, dₖ):

    for j = 1 to dₖ:
        // 步骤 A: 计算误差
        e ← f_curr - G₁ · G₂ · ··· · Gᵣ

        if e = 0: break                     // 提前终止

        // 步骤 B: 提取 (xₖ-αₖ)^j 的 Taylor 系数
        //   e 在每步是 (xₖ-αₖ)^j 的倍数（归纳不变量）
        //   eⱼ = [e / (xₖ-αₖ)^j] |_{xₖ=αₖ}  ∈ Z[x₁,...,xₖ₋₁]
        //   实现: 对 xₖ 做 j 次精确除以 (xₖ-αₖ)，再代入 xₖ=αₖ
        eⱼ ← __taylor_coeff(e, xₖ, αₖ, j)

        if eⱼ = 0: continue

        // 步骤 C: 解多变量丢番图方程
        //   求 δ₁,...,δᵣ 使得 Σ δᵢ·(∏_{j≠i} ĝⱼ) = eⱼ
        //   解: δᵢ = (sᵢ · eⱼ) rem ĝᵢ
        //
        //   注意: eⱼ ∈ Z[x₁,...,xₖ₋₁] 是多变量的，
        //   但 sᵢ, ĝᵢ ∈ Z[x₁] 是单变量的。
        //   "rem ĝᵢ" 指的是 以 x₁ 为主变量 做多项式取余，
        //   Z[x₂,...,xₖ₋₁] 部分作为系数环不参与除法。
        for i = 1 to r:
            δᵢ ← (sᵢ · eⱼ) rem ĝᵢ         // in Z[x₁,...,xₖ₋₁], mod ĝᵢ(x₁)

        // 步骤 D: 更新因子
        for i = 1 to r:
            Gᵢ ← Gᵢ + δᵢ · (xₖ - αₖ)^j

        // 步骤 E: LC 校正（Wang 核心创新）
        //   提升可能扰动 Gᵢ 的首项系数。强制恢复:
        for i = 1 to r:
            lc_target ← σᵢ(x₂,...,xₖ, αₖ₊₁,...,αₙ)  // 部分求值的 lc 分配
            replace lc(Gᵢ, x₁) with lc_target

    // 不变量: 提升结束后 G₁·...·Gᵣ = f_curr (精确等式，非模)
```

**步骤 B 归纳不变量证明：**
- j=1 时：`e = f_curr - ∏Gᵢ`。由于 `Gᵢ|_{xₖ=αₖ}` 是提升前的因子，
  `∏Gᵢ|_{xₖ=αₖ} = f_curr|_{xₖ=αₖ}`，所以 `(xₖ-αₖ) | e`。
- j→j+1：步骤 D 添加 `δᵢ·(xₖ-αₖ)^j` 项，恰好消去 `e` 中 `(xₖ-αₖ)^j` 的贡献，
  使新误差 `e' = f_curr - ∏Gᵢ'` 被 `(xₖ-αₖ)^{j+1}` 整除。

### 6.4 Taylor 系数提取

需要新增辅助函数：

```cpp
// 计算 f 在 xₖ = αₖ 处的第 j 阶 Taylor 系数
// 即: 将 f 写成 Σ cⱼ·(xₖ-αₖ)^j，返回 cⱼ
// 实现: 反复精确除以 (xₖ-αₖ) 再代入
// 前置: (xₖ-αₖ)^j | f
template<class var_order>
polynomial_<ZZ, lex_<var_order>> __taylor_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& xk, const ZZ& alpha_k, int j);
```

实现方式：
```
__taylor_coeff(f, xₖ, αₖ, j):
    g ← f
    for t = 1 to j:
        g ← g / (xₖ - αₖ)          // 精确多项式除法 (pair_vec_div)
    return g |_{xₖ = αₖ}            // assign(g, xₖ, αₖ)
```

> **优化（可选）：** 不必每步从头计算 `e = f_curr - ∏Gᵢ`。
> 可以用增量更新：`e_new = e_old - (∏Gᵢ_new - ∏Gᵢ_old)`。
> 初始实现使用直接计算，性能优化留给后续。

### 6.5 步骤 C 的多变量模运算详解

步骤 C 中 `(sᵢ · eⱼ) rem ĝᵢ` 的含义需要精确说明：

- `sᵢ ∈ Z[x₁]`（单变量，由 §6.2 Bézout 计算得到）
- `eⱼ ∈ Z[x₁,...,xₖ₋₁]`（多变量）
- `ĝᵢ ∈ Z[x₁]`（单变量）

将 `eⱼ` 视为 `(Z[x₂,...,xₖ₋₁])[x₁]` 中的多项式（x₁ 为主变量，
其他变量的多项式作为系数），然后与 `sᵢ` 相乘后对 `ĝᵢ(x₁)` 取余。

这等价于"逐系数"操作：对 eⱼ 的每个关于 x₂,...,xₖ₋₁ 的单项式，
分别与 sᵢ 相乘后 mod ĝᵢ。

> **实现注意：** 这里需要的是 Z[x₁] 上的精确多项式除法（系数可能很大），
> 而非 Zp[x₁] 上的。需要确保 `pair_vec_div` 能正确处理
> 除数为单变量、被除数为多变量的情况（它可以，见 §附录 C）。

### 6.6 终止条件

每个变量 `xₖ` 的提升精度为 `deg(f_scaled, xₖ)`，因为 f 的任何因子在 `xₖ`
上的度数不超过 `deg(f, xₖ)`。

若某步 `e = 0`，表示提升已精确完成，可提前终止当前变量的提升。

### 6.7 函数签名

```cpp
// 多变量 Hensel 提升
// 前置: f_scaled = δ^(r-1)·f_prim
//       scaled_factors = v₁,...,vᵣ (lc 已校正的单变量因子)
//       lc_assignments = σ₁,...,σᵣ (各因子的多变量 lc)
//       ∏ vᵢ = f_scaled(x₁, α)  (eval_point 处乘积等于 f_scaled)
// 后置: 返回 G₁,...,Gᵣ ∈ Z[x₁,...,xₙ] 使得 ∏ Gᵢ = f_scaled
//       lc(Gᵢ, x₁) = σᵢ
//       Gᵢ(x₁,...,α) = vᵢ
//       调用方需对 pp(Gᵢ) 做试除验证
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>>
__multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_assignments,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

### 6.8 需要新增的辅助函数

| 函数 | 用途 | 实现依赖 |
|---|---|---|
| `__taylor_coeff(f, xₖ, αₖ, j)` | 提取 Taylor 系数 | `pair_vec_div` + `assign` |
| `__upoly_gcd_extended_ZZ(s, t, a, b)` | Z[x₁] 上扩展 GCD | 已有 Zp 版本，需 ZZ 适配 |
| `__poly_mod_univar(f, g, x₁)` | 多变量 f 对单变量 g 关于 x₁ 取模 | `pair_vec_div` |

---

## 7. `__factor_multivar` — 多变量分解入口

```cpp
// 多变量因式分解 (Wang 算法)
// 前置: f ∈ Z[x₁,...,xₙ], n ≥ 2, f 非零非常数
// 后置: 返回 factorization<Poly>，满足正确性不变量（§1.2）
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f);
```

### 7.1 完整算法

```
__factor_multivar(f_input):

1.  // 选主变量 + 内容提取
    x₁ ← get_variables(f_input).front().first    // lex 首变量
    c ← cont(f_input, x₁)                        // ∈ Z[x₂,...,xₙ]
    f_prim ← f_input / c                          // 本原部分 (pair_vec_div, 精确)

2.  // 递归分解内容
    //   递归终止: c 的变量数 < f 的变量数（cont 消去了 x₁）
    //   若 c ∈ Z（常数），直接作为 content
    cont_factors ← {}
    if !is_number(c):
        cont_factors ← factorize(c)               // 递归调用（变量数递减）

    retry_count ← 0

3.  // 选取值点
    eval ← __select_eval_point(f_prim, x₁)        // 满足条件 (a)-(d)

4.  // 单变量分解
    f₀ ← assign(f_prim, eval)                     // f₀ ∈ Z[x₁]
    uni_fac ← factorize(f₀)                       // 单变量 M4
    u₁,...,uᵣ ← uni_fac.factors (首一化)

5.  if r ≤ 1:
        // 单变量像不可约 ⇒ f_prim 不可约
        // 理由: f_prim 本原 + eval 满足条件 (a)(b)(c)
        //   若 f_prim = g·h 非平凡分解，则 f₀ = g(x₁,α)·h(x₁,α) 也非平凡
        return {content(c), [(f_prim, 1)] ∪ cont_factors}

6.  // 首项系数校正 (§5)
    lc_result ← __wang_leading_coeff(f_prim, [u₁,...,uᵣ], eval, x₁)
    if !lc_result.success:
        retry_count++
        if retry_count ≥ MAX_RETRY (=10):
            throw "Wang LC distribution failed"
        goto 3                                     // 换求值点

7.  // 多变量 Hensel 提升 (§6)
    mv_factors ← __multivar_hensel_lift(
        lc_result.f_scaled,                        // δ^(r-1) · f_prim
        lc_result.scaled_factors,                  // v₁,...,vᵣ (lc 已校正)
        lc_result.lc_assignments,                  // σ₁,...,σᵣ
        eval, x₁)

8.  // 试除验证 + 去缩放
    //   不变量: pp(Gᵢ) | f_prim（因为 Gᵢ 包含缩放因子 δ 的贡献）
    verified ← []
    f_remaining ← f_prim
    for G in mv_factors:
        g ← pp(G, x₁)                             // 取本原部分，消去缩放
        q, r ← divmod(f_remaining, g)              // 对 f_prim (非 f_scaled!) 试除
        if r = 0:
            verified.push(g)
            f_remaining ← q
    if deg(f_remaining) > 0:
        verified.push(pp(f_remaining, x₁))         // 剩余部分也是因子

    if verified is empty:
        // 提升失败（数值问题），换求值点重试
        retry_count++
        if retry_count ≥ MAX_RETRY: throw
        goto 3

9.  // 合并内容因子 + 排序
    all_factors ← cont_factors ∪ {(g, 1) : g ∈ verified}
    content ← (c 的整数部分) · (f_remaining 若为常数)
    return {content, all_factors}  // 按 (degree, 字典序) 排序
```

### 7.2 递归终止性

`__factor_multivar` 的递归调用路径：

1. **步骤 2**：`factorize(c)` 其中 `c = cont(f, x₁) ∈ Z[x₂,...,xₙ]`，
   变量数严格递减（c 不含 x₁）。
2. **步骤 6**：`factorize(L)` 其中 `L = lc(f, x₁) ∈ Z[x₂,...,xₙ]`，
   变量数同样严格递减。

因此递归以变量数为序数函数，基础情况是 n = 1（由 M4 处理）或 n = 0（常数）。

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
| **Phase 5** | M5: 多变量 Wang | `pp`, `__taylor_coeff`, `__upoly_gcd_extended_ZZ`, `__poly_mod_univar`, `__select_eval_point`, `__wang_leading_coeff` (含 `__wang_lc_result`), `__multivar_hensel_lift` (含 `__hensel_lift_one_var`), `__factor_multivar`, `factorize` 多变量 dispatch | M4 (已实现) |
| **Phase 6** | 增强：van Hoeij 重组 | `__factor_recombine_van_hoeij` + LLL 实现 | M3 替换 |
| **Phase 7** | 增强：Zippel 后备 | 稀疏插值模块 + Zippel 算法 | Phase 5 后备 |
| **Phase 8** | 终极：MTSHL | 二变量 Hensel 提升 + 稀疏插值驱动的多变量分解 | 替换 Phase 5 |

### 10.1 测试计划

| 可独立测试的函数 | 验证方法 |
|---|---|
| `pp(f)` | 验证 `cont(f) · pp(f) == f` |
| `__taylor_coeff(f, xₖ, αₖ, j)` | 验证 `Σ coeff_j · (xₖ-αₖ)^j == f` |
| `__upoly_gcd_extended_ZZ` | 验证 `s·a + t·b = gcd`，与 Zp 版本结果一致 |
| `__poly_mod_univar(f, g, x₁)` | 验证 `f = q·g + r`，`deg(r,x₁) < deg(g)` |
| `__select_eval_point` | 验证条件 (a)-(d)：无平方、lc 非零、度数守恒、lc 因子可分辨 |
| `__wang_leading_coeff` | 构造已知分解的多项式，验证 `∏σᵢ(α) = δ`，`∏vᵢ = f_scaled(x₁,α)` |
| `__multivar_hensel_lift` | 二变量多项式提升后 `∏Gᵢ = f_scaled`，`pp(Gᵢ) \| f_prim` |
| `__factor_multivar` | 与 Mathematica `Factor[f]` 对比 |
| `factorize` (多变量入口) | 与 Mathematica `FactorList[f]` 对比；`verify_factorization` 重组检查 |

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

## 11. 远期目标：MTSHL（Maple 路线）

### 11.1 动机

Wang 算法有两个根本性瓶颈：

1. **LC 分配问题**：将 lc(f, x₁) 的因子正确分配给各模因子，分配可能失败或不唯一，
   需要换求值点重试。Singular 为此实现了 4 级级联启发式，FLINT 维护了 Wang + Kaltofen
   两条路径。

2. **多变量丢番图问题 (MDP)**：每步 Hensel 提升需要解
   Σ sᵢ·δᵢ ≡ e mod ∏gⱼ，对稀疏多项式可能退化为指数级复杂度。

Maple 在 2019 年引入的 MTSHL（Monagan-Tuncer Sparse Hensel Lifting）算法
从根本上绕过了这两个问题。

### 11.2 核心思想

MTSHL 用**稀疏插值**替代经典 MDP：

```
经典 Wang:
  单变量像 → 逐变量 Hensel 提升 (每步解 MDP) → 试除验证

MTSHL:
  多个二变量像 → 二变量 Hensel 提升 (BHL) → 稀疏插值恢复因子 → 试除验证
```

关键优势：
- **不需要 LC 预分配**——稀疏插值自然恢复每个因子的首项系数
- **不需要解 MDP**——用 Vandermonde 系统替代
- **复杂度取决于项数而非变量数**——对稀疏多项式（实际中绝大多数情况）优势巨大

### 11.3 所需基础设施

| 组件 | 状态 | 说明 |
|---|---|---|
| 稀疏插值 (Ben-Or/Tiwari 或 Zippel) | 未实现 | MTSHL 的核心依赖 |
| 二变量 Hensel 提升 (BHL) | 未实现 | 从 `f(x₁, α₂+t·x₂)` 提升为 `f(x₁, x₂)` 的因子 |
| Hilbert 点选取 | 未实现 | 概率框架，保证失败概率可控 |
| 非首一二变量提升 | 未实现 | 处理 lc(f, x₁) 非常数的情况 |

### 11.4 演进路线

```
Phase 5 (Wang)  ──→  Phase 7 (Zippel 后备)  ──→  Phase 8 (MTSHL)
     │                      │                          │
     │                      ▼                          ▼
     │               稀疏插值模块               完整替换 Wang
     │               (可独立使用)              (保留 Wang 作后备)
     ▼
  多变量分解可用
```

Phase 7（Zippel 后备）是过渡步骤：它引入稀疏插值模块，先作为 Wang 失败时的后备
路径（FLINT 的做法：Wang → Zippel → Zassenhaus 三级级联），同时为 Phase 8 的
MTSHL 积累基础设施。

### 11.5 参考文献

- Monagan & Tuncer, "Using Sparse Interpolation in Hensel Lifting", CASC 2016
- Monagan & Tuncer, "Polynomial Factorization in Maple 2019", MACIS 2019
- Tian Chen, "Sparse Hensel Lifting Algorithms for Multivariate Polynomial
  Factorization", SFU Master's Thesis, 2019
  (https://www.cecm.sfu.ca/CAG/theses/tian.pdf)

---

## 附录 A: 完整函数签名索引

### 辅助函数（§3, §6）

```cpp
// §3.1 多变量本原部分
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
pp(const polynomial_<ZZ, lex_<var_order>>& f);

// §3.2 多变量系数 L1 范数
template<class var_order>
ZZ __poly_coeff_l1_norm(const polynomial_<ZZ, lex_<var_order>>& f);

// §6.4 Taylor 系数提取
template<class var_order>
polynomial_<ZZ, lex_<var_order>> __taylor_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& xk, const ZZ& alpha_k, int j);

// §6.8 Z[x] 上扩展 GCD
inline void __upoly_gcd_extended_ZZ(
    upolynomial_<ZZ>& s, upolynomial_<ZZ>& t,
    const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b);

// §6.8 多变量 f 对单变量 g 关于 x₁ 取模
template<class var_order>
polynomial_<ZZ, lex_<var_order>> __poly_mod_univar(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const upolynomial_<ZZ>& g, const variable& main_var);
```

### 核心模块（§4-§7）

```cpp
// §4 选取值点
template<class var_order>
std::map<variable, ZZ> __select_eval_point(
    const polynomial_<ZZ, lex_<var_order>>& f, const variable& main_var);

// §5 首项系数校正结果
template<class var_order>
struct __wang_lc_result {
    bool success;
    polynomial_<ZZ, lex_<var_order>> f_scaled;                      // δ^(r-1)·f
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments;   // σ₁,...,σᵣ
    std::vector<upolynomial_<ZZ>> scaled_factors;                   // v₁,...,vᵣ
};

// §5 首项系数校正 (所有输入 const, 输出通过返回值)
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

// §6 多变量 Hensel 提升
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>> __multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_assignments,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

// §7 多变量分解入口
template<class var_order>
factorization<polynomial_<ZZ, lex_<var_order>>>
__factor_multivar(const polynomial_<ZZ, lex_<var_order>>& f);
```

---

## 附录 B: 参考文献

- **Wang**: Wang, "An Improved Multivariate Polynomial Factoring Algorithm", Math. Comp. 1978
- **Kaltofen-Shoup**: Kaltofen & Shoup, "Subquadratic-Time Factoring of Polynomials over Finite Fields", Math. Comp. 1998
- **GCL**: Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992 (§16)
- **MCA**: von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013 (§16)
- **MTSHL**: Monagan & Tuncer, "Using Sparse Interpolation in Hensel Lifting", CASC 2016
- **FLINT**: Hart et al., FLINT: Fast Library for Number Theory, https://flintlib.org
- **Singular/Factory**: Singular Team, Factory Library, https://www.singular.uni-kl.de

## 附录 C: 已有函数依赖清单

以下是 M5 直接依赖的已有函数：

| 函数 | 文件 | 用途 | 备注 |
|---|---|---|---|
| `factorize(polynomial_<ZZ>)` | `polynomial_factorize.hh` | 单变量分解 | M4，已实现 |
| `squarefreefactorize(F)` | `polynomial_gcd.hh` | 无平方分解 | 支持多变量 |
| `cont(F)` | `polynomial_gcd.hh` | 内容提取 (关于 lex 首变量) | 返回 `polynomial_<ZZ,lex>` |
| `assign(f, v, c)` | `polynomial.hh` | 单变量代入 | |
| `assign(f, map)` | `polynomial.hh` | 多变量批量代入 | 接受 `map<variable, Tc>` |
| `is_squarefree(f)` | `polynomial_gcd.hh` | 无平方检测 | 支持多变量 |
| `get_variables(f)` | `polynomial_.hh` | 获取变量列表 | 返回 `list<pair<variable, int64_t>>` |
| `pair_vec_div(q, r, f, g, comp)` | `basic.hh` | 多项式除法 | 支持多变量精确除法 |
| `__upoly_gcd_extended(s, t, a, b)` | `polynomial_factorize.hh` | Zₚ 上扩展 GCD | 需补充 ZZ 版本 |
| `is_number(f)` | `upolynomial.hh` | 常数检测 | |
| `poly_convert(in, out)` | `upolynomial.hh` | polynomial ↔ upolynomial | |
