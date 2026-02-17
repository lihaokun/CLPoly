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
  │
  ├── squarefreefactorize(f)               无平方分解 (已有, 内部自动提取 content)
  │     返回 [(g₁,m₁), (g₂,m₂), ...]
  │     对每个 gₖ：
  │       不含 x₁ → factorize(gₖ) 递归（变量数递减）
  │       含 x₁, 单变量 → factorize(gₖ) (M4)
  │       含 x₁, 多变量 → __wang_core(gₖ)   ← gₖ 已本原且无平方
  │
  └── __wang_core(g)                       Wang 算法核心（输入须本原无平方）
        ├── __select_eval_point            选取值点
        │     ├── assign(f, v, c)         代入求值 (已有)
        │     ├── is_squarefree           无平方检测 (已有)
        │     └── factorize(lc)           lc 因子分解 (条件 d, 已有)
        ├── factorize(f₀)                  单变量分解 (M4, 已实现)
        ├── __wang_leading_coeff           首项系数分配
        ├── __multivar_hensel_lift         多变量 Hensel 提升
        │     ├── __upoly_gcd_extended     Z[x₁] pseudo-XGCD (ZZ 重载，新增)
        │     ├── __taylor_coeff           Taylor 系数提取 (新增)
        │     └── __poly_prem_univar      多变量 prem 单变量 (新增)
        └── 试除验证: pp(Gᵢ) | g            pp (新增, polynomial_gcd.hh)
```

### 2.1 数据流与缩放不变量

Wang 算法中 `f` 经历多次变换，各模块操作的对象不同，必须严格区分：

```
f_input                              用户输入
  │
  ▼
squarefreefactorize(f_input)         无平方分解（内部自动 cont 提取 + Yun）
  │
  ├─ 返回 [(g₁,m₁), (g₂,m₂), ...]
  │   每个 gₖ 是不可约因子的无平方积
  │
  │   分类处理每个 gₖ：
  │     常数 / 单变量       → factorize(gₖ) (M4)，重数 ×mₖ
  │     多变量但不含 x₁     → factorize(gₖ) 递归（变量数递减），重数 ×mₖ
  │     多变量且含 x₁       → __wang_core(gₖ)，重数 ×mₖ
  │                           ↓
  ▼ ════════════════════════════════════════
  __wang_core(g)               g 已本原、无平方、多变量
  │
  ├─ α = __select_eval_point   选求值点，满足 (a)-(d)
  ├─ f₀ = g(x₁, α)            单变量像 ∈ Z[x₁]
  ├─ u₁,...,uᵣ = factorize(f₀) 单变量因子（本原, lc>0, 非首一）
  │   若 r ≤ 1 → g 不可约，直接返回
  ▼
  __wang_leading_coeff
  ├─ 输入:  g, u₁,...,uᵣ, α
  ├─ 输出:  f_scaled = δ^(r-1) · g              ← 缩放后的多项式
  │         σ₁,...,σᵣ ∈ Z[x₂,...,xₙ]            ← 各因子的 lc 分配 (∏σᵢ = L)
  │         τ₁,...,τᵣ ∈ Z[x₂,...,xₙ]            ← 提升用 lc 目标 (∏τᵢ = δ^(r-1)·L)
  │           τᵢ = (δ / σᵢ(α)) · σᵢ,  满足 τᵢ(α) = δ
  │         v₁,...,vᵣ ∈ Z[x₁]                    ← lc 统一为 δ 的单变量因子
  │           (vᵢ = δ · ūᵢ, 满足 ∏vᵢ = f_scaled(x₁,α))
  │   失败 → 换求值点重试
  ▼
  __multivar_hensel_lift
  ├─ 输入:  f_scaled, v₁,...,vᵣ, τ₁,...,τᵣ, α
  ├─ Bézout: sᵢ ∈ Z[x₁], Σ sᵢ·V̂ᵢ = denom     （基于 vᵢ 计算一次，V̂ᵢ = ∏_{j≠i} vⱼ）
  ├─ 逐变量 x₂,...,xₙ 线性提升（步骤 A-E）
  ├─ 输出:  G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]            ← 候选因子
  ▼
  试除验证
  ├─ 对 g（非 f_scaled!）做试除: pp(Gᵢ) | g ?
  ├─ 自洽性检查，失败 → 换求值点重试
  └─ 输出: g 的不可约因子列表
```

**关键不变量：**
- `f_scaled = δ^(r-1) · g`，其中 `δ = lc(g, x₁)(α)`
- `lc(f_scaled, x₁) = δ^(r-1) · L`，其中 `L = lc(g, x₁)`
- Hensel 提升在 `f_scaled` 上进行，各因子的 lc 目标为 `τᵢ = (δ/σᵢ(α))·σᵢ`
  - `∏τᵢ = δ^(r-1) · L = lc(f_scaled, x₁)` ← 与 f_scaled 一致
  - `τᵢ(α) = δ` ← 与 vᵢ 的 lc 一致
  - 注意 `∏σᵢ = L ≠ lc(f_scaled, x₁)`，因此 σᵢ 不能直接用作 Hensel lc 目标
- 试除在 `g`（原始本原无平方多项式）上进行：`pp(Gᵢ)` 消去缩放因子 `δ`
- 这与单变量 M3 的试除逻辑一致（见 univariate §6.4）

---

## 3. 需要新增的已有模块补充函数

### 3.1 `pp` — 多变量本原部分（`polynomial_gcd.hh`，`cont()` 伴侣）

```cpp
// 多变量本原部分: f / cont(f)
// 前置: f ∈ Z[x₁,...,xₙ], lex 排序, f 非零
// 后置: 返回 f 除以其关于首变量的内容后的本原部分
// 位置: polynomial_gcd.hh，紧随 cont() 定义之后
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
pp(const polynomial_<ZZ, lex_<var_order>>& f);
```

实现方式：内部调用 `cont(f)` 后精确除，避免调用方重复计算内容。

> **注：** `cont()` 已存在（`polynomial_gcd.hh:468`）。`pp()` 作为其伴侣函数，
> 所有需要本原部分的场景（试除验证等）直接调用 `pp()`，
> 不必手动 `cont()` + `pair_vec_div` 两步。

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

### 4.2 算法

```
__select_eval_point(f, x₁):
    vars ← get_variables(f) \ {x₁}      // n-1 个需要赋值的变量
    L ← lc(f, x₁)                       // ∈ Z[x₂,...,xₙ]

    // 预处理: 若 L 非常数，分解 L 用于条件 (d)
    lc_irr_factors ← []
    if !is_number(L):
        lc_fac ← factorize(L)           // 递归调用（变量数更少）
        lc_irr_factors ← [lⱼ for (lⱼ, eⱼ) in lc_fac.factors]

    // 枚举候选点（按 L∞ 范数递增）
    for bound = 0, 1, 2, ...:
        for each α ∈ {v → a : v ∈ vars, a ∈ [-bound..bound]}:
            // 跳过已检查过的更小范数的点
            if max(|αᵢ|) < bound: continue

            // 条件 (b): lc 非零
            δ ← L(α)                    // assign(L, α) 若 L 是多项式，否则 δ = L
            if δ == 0: continue

            // 条件 (a): 无平方
            f₀ ← assign(f, α)           // f₀ ∈ Z[x₁]
            if !is_squarefree(f₀): continue

            // 条件 (d): lc 因子求值两两互素（仅当 L 非常数）
            if lc_irr_factors 非空:
                vals ← [|lⱼ(α)| for lⱼ in lc_irr_factors]
                if 存在 i ≠ j 使得 gcd(vals[i], vals[j]) > 1:
                    continue

            return α

    // 不可达（Z 无限域保证终止）
```

> **实现注意：** 对 n-1 个变量枚举所有 `[-bound, bound]` 组合是指数级的。
> 实践中绝大多数情况在 bound ≤ 2 内找到合法点。
> 可优化为：先尝试全零，再逐个变量非零，最后才尝试多变量非零组合。

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
`f₀ = c · u₁ · u₂ · ··· · uᵣ`（uᵢ 本原，lc > 0，一般非首一）。

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
        // 首项系数为常数，无需 lc 因子分配，直接进入步骤 4 缩放
        σᵢ ← 1 for all i; σ₁ ← L
        // 跳过步骤 2-3，进入步骤 4

2.  // 递归分解首项系数
    lc_fac ← factorize(L)                  // 递归调用 factorize
    // lc_fac = {γ, [(l₁,e₁), (l₂,e₂), ...]}
    // γ ∈ Z 是整数内容，每个 lⱼ ∈ Z[x₂,...,xₙ] 不可约本原

3.  // 分配 lc 因子到单变量因子
    for i = 1 to r:
        σᵢ ← 1                              // 累积: 分配给 uᵢ 的 lc 多项式 ∈ Z[x₂,...,xₙ]
        wᵢ ← lc(uᵢ)                         // 数值跟踪: uᵢ 的 lc 中尚未被解释的部分

    // 按 eⱼ 从大到小排序（高次幂优先，减少歧义）
    sort (lⱼ, eⱼ) by eⱼ descending

    for each (lⱼ, eⱼ):
        zⱼ ← lⱼ(α) ^ eⱼ                   // 该 lc 因子求值后的整数值
        // 找唯一的 uᵢ 使得 zⱼ | wᵢ
        candidates ← {i : zⱼ | wᵢ}
        if |candidates| ≠ 1:
            return FAIL                      // 无法唯一分配 → 换求值点
        i ← candidates 中唯一的元素
        σᵢ ← σᵢ · lⱼ^eⱼ
        wᵢ ← wᵢ / zⱼ

    // 吸收整数内容 γ 到 σ₁
    σ₁ ← γ · σ₁
    // 现在: ∏ σᵢ = γ · ∏ lⱼ^eⱼ = L (精确等式，作为多项式)
    //       ∏ σᵢ(α) = L(α) = δ

    // wᵢ 的剩余值: ∏wᵢ = ∏lc(uᵢ) / ∏(分配的 vⱼ) = ∏lc(uᵢ) / (δ/γ)
    // 这些剩余量不影响算法——σᵢ 经步骤 5 转化为 τᵢ 用于 Hensel LC 校正，vᵢ 用统一缩放

4.  // 缩放: 统一的 δ^(r-1) 方案
    f_scaled ← δ^(r-1) · f
    // 将单变量因子首一化后乘以 δ
    for i = 1 to r:
        ūᵢ ← uᵢ / lc(uᵢ)                  // 首一化
        vᵢ ← δ · ūᵢ                        // 所有因子 lc = δ (整数!)
    // 验证:
    //   由 f₀ = c₀·∏uᵢ 和 lc(f₀) = c₀·∏lc(uᵢ) = δ，得 ∏ūᵢ = ∏(uᵢ/lc(uᵢ)) = f₀/δ
    //   ∴ ∏ vᵢ = δ^r · ∏ ūᵢ = δ^r · (f₀/δ) = δ^(r-1) · f₀ = f_scaled(x₁, α) ✓

5.  // 构造 Hensel 提升用 lc 目标 τᵢ
    //
    // 问题：∏σᵢ = L，但 lc(f_scaled, x₁) = δ^(r-1)·L ≠ L
    // 若直接用 σᵢ 做 LC 校正，则 ∏lc(Gᵢ) ≠ lc(f_scaled)，误差 e 在 x₁ 最高次项
    // 不会消去，导致 deg(eⱼ, x₁) = deg(∏vᵢ)，MDP 无解。
    //
    // 解法：将 δ^(r-1) 因子分摊到各 σᵢ 中，构造 τᵢ = (δ / σᵢ(α)) · σᵢ
    for i = 1 to r:
        τᵢ ← (δ / σᵢ(α)) · σᵢ             // 精确整除（δ = ∏σⱼ(α)，σᵢ(α) | δ）
    // 验证:
    //   τᵢ(α) = (δ/σᵢ(α)) · σᵢ(α) = δ                        ← 与 lc(vᵢ) = δ 一致 ✓
    //   ∏τᵢ = ∏((δ/σᵢ(α))·σᵢ) = (∏(δ/σᵢ(α)))·(∏σᵢ)
    //        = (δ^r / ∏σᵢ(α)) · L = (δ^r / δ) · L = δ^(r-1)·L
    //        = lc(f_scaled, x₁)                                  ← 与 f_scaled 一致 ✓

    return SUCCESS, f_scaled, σ₁,...,σᵣ, τ₁,...,τᵣ, v₁,...,vᵣ
```

> **与单变量 M3 的类比：** 单变量 Hensel 提升也将 `lc(f)` 乘到 `factors[0]`。
> 多变量版本的缩放 `δ^(r-1)` 起类似作用，但更系统化——每个因子有明确的 lc 分配。
> 试除时同样用 `pp()` 消去缩放因子（见 §2.1 数据流）。

### 5.3 正确性条件

| 条件 | 保证 |
|---|---|
| `L(α) ≠ 0` | 由 eval_point 条件 (b) |
| 各 lⱼ(α) 两两互素 | 由 eval_point 条件 (d)，保证唯一分配 |
| `∏ σᵢ = L` | 分配完整性，所有 lc 因子被完全消耗 |
| `∏ σᵢ(α) = δ` | 由 `∏σᵢ = L` 和 `L(α) = δ` 推出 |
| `τᵢ(α) = δ` | 与 `lc(vᵢ) = δ` 一致，Hensel 初始条件 |
| `∏ τᵢ = δ^(r-1) · L` | 与 `lc(f_scaled, x₁)` 一致，确保提升后 lc 乘积正确 |
| `∏ vᵢ = f_scaled(x₁, α)` | 缩放正确性，Hensel 提升的前置条件 |

### 5.4 函数签名

```cpp
// 首项系数校正结果
template<class var_order>
struct __wang_lc_result {
    bool success;
    polynomial_<ZZ, lex_<var_order>> f_scaled;          // δ^(r-1) · f
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_assignments;  // σ₁,...,σᵣ (原始 lc 分配)
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_targets;      // τ₁,...,τᵣ (提升用 lc 目标)
    std::vector<upolynomial_<ZZ>> scaled_factors;       // v₁,...,vᵣ (修改后的单变量因子)
};

// 首项系数校正
// 前置: f 关于 x₁ 本原, univar_factors = [u₁,...,uᵣ] 是 f(x₁,α) 的本原因子（lc>0，非首一）
//       eval_point 满足条件 (a)-(d)
// 后置: success=true 时，f_scaled/lc_assignments/lc_targets/scaled_factors 已填充
//       success=false 时，需换求值点重试
// 不修改输入参数（所有输出通过返回值）
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

> **接口设计说明：** `lc_assignments`（σᵢ）保留用于调试和验证，
> 实际 Hensel 提升使用 `lc_targets`（τᵢ）。两者的关系：`τᵢ = (δ/σᵢ(α))·σᵢ`。

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
    __hensel_lift_one_var(f_curr, G₁,...,Gᵣ, v₁,...,vᵣ, s₁,...,sᵣ, τ₁,...,τᵣ, xₖ, αₖ, dₖ)
    // 每步后 Gᵢ 从 Z[x₁,...,xₖ₋₁] 扩展为 Z[x₁,...,xₖ]
```

### 6.2 Bézout 系数：计算一次，全局复用

**关键：** Bézout 系数基于缩放后的单变量因子 `vᵢ` 计算，只算一次。

设 `vᵢ = δ · ūᵢ = δ · uᵢ/lc(uᵢ) ∈ Z[x₁]`（§5.2 step 4 产出），
所有 `lc(vᵢ) = δ`，且 vᵢ 两两互素（因 ūᵢ 两两互素）。

> **为什么用 vᵢ 而非 ūᵢ？** 当 uᵢ 非首一（lc(uᵢ) > 1）时，
> ūᵢ = uᵢ/lc(uᵢ) ∈ Q[x₁] \ Z[x₁]。例如 u₁ = 2x₁+1 则 ū₁ = x₁+½。
> 而 vᵢ = δ·ūᵢ 保证 vᵢ ∈ Z[x₁]，可直接在 Z[x₁] 上做 XGCD。

Bézout 系数 `s₁,...,sᵣ ∈ Z[x₁]` 满足偏分式恒等式：

```
Σᵢ sᵢ · (∏_{j≠i} vⱼ) = denom    in Z[x₁]
```

其中 `deg(sᵢ) < deg(vᵢ)`。此恒等式精确成立（非模某个理想），
因为 vᵢ 两两互素。

**计算方法（逐对 XGCD 链，与 Singular 一致）：**

```
g_acc ← v₁
s[1] ← 1
denom ← 1
for i = 2 to r:
    (α, β, c_k) ← XGCD_ZZ(g_acc, vᵢ)  // α·g_acc + β·vᵢ = c_k ∈ Z
    for j = 1 to i-1:
        s[j] ← s[j] · β                   // 不做 mod 约化（见下方说明）
    s[i] ← denom · α                      // denom_old · α
    denom ← denom · c_k                    // 公分母累积
    g_acc ← g_acc · vᵢ
// 最终: Σ s[i]·V̂_i = denom   (V̂ᵢ = ∏_{j≠i} vⱼ)
```

> **为什么不做 mod 约化？** 对首一多项式，`rem` 在 Z[x₁] 上精确。但 vᵢ 非首一
> （lc = δ），Z[x₁] 上 `rem vⱼ` 需要 prem，每次 prem 引入 δ^k 乘子，
> 需要同步乘到所有其他 s[m] 上，导致 denom 爆炸式增长。
>
> 去掉 mod 后，deg(sᵢ) 可能超过 deg(vᵢ)，但不影响正确性——
> MDP 步骤的 prem 会处理度数约化。代价是 sᵢ 系数稍大，
> 对实际中 r ≤ 5 的情况影响很小。
>
> 若需优化（r 较大），可在 Q[x₁] 上对首一化的 ūᵢ = vᵢ/δ 做链计算
> （此时 mod 精确），再统一通分清分母。

> **XGCD 实现决策：在 Z[x₁] 上计算，显式跟踪公分母。**
>
> Z[x₁] 上 vᵢ 非首一（lc = δ），Euclidean 算法会产生中间分数。
> 使用 pseudo-XGCD（类似 subresultant PRS）保持整系数：
>
> `__upoly_gcd_extended(s, t, c, a, b)` (ZZ 重载) 满足 `s·a + t·b = c`，
> 其中 `s, t ∈ Z[x₁]`，`c ∈ Z \ {0}`。
>
> 对 r 个因子的 Bézout 构造累积公分母 `denom = ∏ cₖ`，
> 使得 `Σ sᵢ·V̂ᵢ = denom`（V̂ᵢ = ∏_{j≠i} vⱼ）。
>
> **MDP 中的使用：** `δᵢ = prem(sᵢ · eⱼ, vᵢ) / (δ^k · denom)`，
> 其中 prem 是伪余式，`δ^k` 是 prem 引入的 lc(vᵢ)^k 因子。
> 由偏分式分解唯一性，除法为精确整除。
>
> **优势：** 全程使用 `mpz_class`，避免 `mpq_class` 每次运算的 GCD 规约开销。
> 实现方式：复用已有 `__upoly_gcd_extended` 的 Euclidean 框架，
> 在每步做 pseudo-division（乘以 lc 后除）以保持整系数。

这些 sᵢ 和 vᵢ 在所有变量的提升过程中保持不变——它们只依赖求值点处的单变量因子。

### 6.3 单变量提升步（核心）

对变量 `xₖ` 的提升。此时：
- `Gᵢ ∈ Z[x₁,...,xₖ₋₁]`（已经过前面变量的提升）
- 目标：将 `Gᵢ` 扩展为 `Z[x₁,...,xₖ]` 上的多项式
- `f_curr = f_scaled(x₁,...,xₖ, αₖ₊₁,...,αₙ) ∈ Z[x₁,...,xₖ]`

**算法：**

```
__hensel_lift_one_var(f_curr, G₁,...,Gᵣ, v₁,...,vᵣ, s₁,...,sᵣ, τ₁,...,τᵣ, xₖ, αₖ, dₖ):

    // 步骤 E（预步骤）: LC 校正——在提升开始前强制替换各因子的 lc
    //
    //   关键: LC 校正必须在误差计算 **之前** 执行，否则 e 的 x₁ 最高次
    //   系数不为零（因为 ∏lc(Gᵢ) ≠ lc(f_curr)），导致 deg(eⱼ,x₁) ≥ deg(∏vᵢ)，
    //   MDP 无解。
    //
    //   使用 τᵢ（而非 σᵢ）作为 lc 目标，因为 ∏τᵢ = δ^(r-1)·L = lc(f_scaled)，
    //   保证提升后 ∏lc(Gᵢ) = lc(f_curr)。
    for i = 1 to r:
        d ← deg(Gᵢ, x₁)
        lc_target ← assign(τᵢ, {xₖ₊₁→αₖ₊₁,...,xₙ→αₙ})  // τᵢ 的部分求值
        lc_current ← coeff(Gᵢ, x₁, d)                     // 当前首项系数
        if lc_target ≠ lc_current:
            Gᵢ ← Gᵢ + (lc_target - lc_current) · x₁^d
        //
        // 实现: lex 序下 x₁ 最高的项在 pair_vec 前端,
        // 找到所有 deg(x₁)=d 的项, 替换为 lc_target·x₁^d 的各项。
        // lc_target ∈ Z[x₂,...,xₖ], lc_current ∈ Z[x₂,...,xₖ],
        // 两者之差乘以 x₁^d 后做 pair_vec 加法即可。

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

        // 步骤 C: 解多变量丢番图方程 (MDP)
        //   求 δ₁,...,δᵣ 使得 Σ δᵢ·V̂ᵢ = eⱼ  (V̂ᵢ = ∏_{j≠i} vⱼ)
        //
        //   因为 vᵢ 非首一（lc = δ），对 vᵢ 做 rem 需要 pseudo-remainder：
        //   prem(f, vᵢ) 满足 δ^k · f = q·vᵢ + prem(f, vᵢ)，k = deg(f) - deg(vᵢ) + 1
        //
        //   解法: δᵢ = prem(sᵢ · eⱼ, vᵢ) / (δ^k · denom)
        //   其中 δ^k 来自 prem 引入的 lc 系数幂，denom 来自 Bézout 公分母。
        //   由偏分式唯一性，此除法为精确整除。
        //
        //   注意: eⱼ ∈ Z[x₁,...,xₖ₋₁] 是多变量的，
        //   但 sᵢ, vᵢ ∈ Z[x₁] 是单变量的。
        //   "prem vᵢ" 指的是 以 x₁ 为主变量 做伪除法取余，
        //   Z[x₂,...,xₖ₋₁] 部分作为系数环不参与除法。
        for i = 1 to r:
            kᵢ ← max(deg(sᵢ·eⱼ, x₁) - deg(vᵢ) + 1, 0)
            δᵢ ← prem(sᵢ · eⱼ, vᵢ) / (δ^kᵢ · denom)  // 精确整除，kᵢ 依赖具体 i

        // 步骤 D: 更新因子
        for i = 1 to r:
            Gᵢ ← Gᵢ + δᵢ · (xₖ - αₖ)^j

        // 步骤 E（再执行一次 LC 校正，修正步骤 D 引入的 lc 偏差）
        for i = 1 to r:
            d ← deg(Gᵢ, x₁)
            lc_target ← assign(τᵢ, {xₖ₊₁→αₖ₊₁,...,xₙ→αₙ})
            lc_current ← coeff(Gᵢ, x₁, d)
            if lc_target ≠ lc_current:
                Gᵢ ← Gᵢ + (lc_target - lc_current) · x₁^d

    // 不变量: 提升结束后 G₁·...·Gᵣ = f_curr (精确等式，非模)
```

> **LC 校正位置说明（关键设计决策）：**
>
> LC 校正（步骤 E）出现两次：j-loop 之前（预步骤）和每步 D 之后。
> - **预步骤**：在首次计算误差之前就将 lc(Gᵢ) 设为正确值，
>   确保 e 的 x₁ 最高次项被正确消去，使 deg(eⱼ, x₁) < deg(∏vᵢ)。
> - **步骤 D 之后**：更新 δᵢ·(xₖ-αₖ)^j 可能扰动 lc(Gᵢ)（若 δᵢ 含 x₁^d 项），
>   需再次校正。
>
> 这与 "Maths From Nothing" 的参考实现一致：
> "target leading coefficients are applied at the beginning of each outer loop iteration,
> before computing residuals"。

**步骤 B 归纳不变量证明：**
- j=1 时：步骤 E 预处理使 `∏lc(Gᵢ) = lc(f_curr)`。
  `e = f_curr - ∏Gᵢ`，由 `Gᵢ|_{xₖ=αₖ} = (提升前因子)` 且 lc 已校正，
  `∏Gᵢ|_{xₖ=αₖ} = f_curr|_{xₖ=αₖ}`，所以 `(xₖ-αₖ) | e`。
  又因 lc 已对齐，`deg(eⱼ, x₁) < deg(∏vᵢ)`，MDP 有解。
- j→j+1：步骤 D+E 消去 `e` 中 `(xₖ-αₖ)^j` 的贡献并保持 lc 正确，
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

步骤 C 中 `δᵢ = prem(sᵢ · eⱼ, vᵢ) / (δ^k · denom)` 的含义需要精确说明：

- `sᵢ ∈ Z[x₁]`（单变量，由 §6.2 Bézout 计算得到，满足 `Σ sᵢ·V̂ᵢ = denom`）
- `eⱼ ∈ Z[x₁,...,xₖ₋₁]`（多变量）
- `vᵢ ∈ Z[x₁]`（单变量，lc = δ，**非首一**）
- `denom ∈ Z`（Bézout 公分母）
- `δ = lc(g, x₁)(α) ∈ Z`（缩放常数）

将 `eⱼ` 视为 `(Z[x₂,...,xₖ₋₁])[x₁]` 中的多项式（x₁ 为主变量，
其他变量的多项式作为系数），然后与 `sᵢ` 相乘后对 `vᵢ(x₁)` 做伪除法取余。

因为 vᵢ 非首一（lc = δ），普通 rem 在 Z 上不精确。使用 pseudo-remainder：
```
prem(f, vᵢ): δ^k · f = q · vᵢ + prem(f, vᵢ),  k = max(deg(f,x₁) - deg(vᵢ) + 1, 0)
```
prem 在 Z[x₁,...,xₖ₋₁] 上精确（无需分数）。
最后除以 `δ^k · denom` 也是精确的（由偏分式唯一性保证 δᵢ 为整系数）。

这等价于"逐系数"操作：对 eⱼ 的每个关于 x₂,...,xₖ₋₁ 的单项式，
分别与 sᵢ 相乘后 prem vᵢ，再整体除以 `δ^k · denom`。

> **实现注意：**
> - 使用 `prem`（CLPoly 已有）而非 `pair_vec_div`，因为 vᵢ 非首一。
> - 也可以先将 eⱼ 乘以 δ^k "预乘法"，然后用 `pair_vec_div`（lex 首变量长除法），
>   效果等价。`pair_vec_div` 能正确处理除数为单变量、被除数为多变量的情况，
>   **前提是 x₁ 在 lex 序中是最高变量**。
> - CLPoly 的 `cont()` 和 `squarefreefactorize` 也依赖此 lex 首变量约定，设计一致。
>
> **优化提示：** 实践中可以在 Bézout 构造阶段将 denom 和 δ 的幂合并跟踪，
> 避免 MDP 中每次重新计算 δ^k。

### 6.6 终止条件

每个变量 `xₖ` 的提升精度为 `deg(f_scaled, xₖ)`，因为 f 的任何因子在 `xₖ`
上的度数不超过 `deg(f, xₖ)`。

若某步 `e = 0`，表示提升已精确完成，可提前终止当前变量的提升。

### 6.7 函数签名

```cpp
// 单变量提升步（内部函数）
// 前置: Gᵢ ∈ Z[x₁,...,xₖ₋₁], f_curr ∈ Z[x₁,...,xₖ]
//       ∏ Gᵢ = f_curr |_{xₖ=αₖ} (仅在 LC 校正后成立)
//       vᵢ ∈ Z[x₁] 是缩放后的单变量因子（lc=δ, 非首一，不随提升变化）
//       sᵢ ∈ Z[x₁] 是 Bézout 系数，满足 Σ sᵢ·V̂ᵢ = bezout_denom
//       τᵢ ∈ Z[x₂,...,xₙ] 是 LC 目标（∏τᵢ = δ^(r-1)·L，不随提升变化）
// 后置: Gᵢ 扩展为 Z[x₁,...,xₖ]，∏ Gᵢ = f_curr
// 修改: Gᵢ 原地更新
template<class var_order>
void __hensel_lift_one_var(
    const polynomial_<ZZ, lex_<var_order>>& f_curr,
    std::vector<polynomial_<ZZ, lex_<var_order>>>& G,        // 因子（原地更新）
    const std::vector<upolynomial_<ZZ>>& bezout_s,           // Bézout sᵢ ∈ Z[x₁]
    const ZZ& bezout_denom,                                   // Bézout 公分母
    const std::vector<upolynomial_<ZZ>>& v_factors,           // vᵢ ∈ Z[x₁]（lc=δ, 非首一）
    const ZZ& delta,                                          // δ = lc(g, x₁)(α)
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_tau,    // τᵢ (LC 目标)
    const variable& xk, const ZZ& alpha_k, int dk);

// 多变量 Hensel 提升（外层入口）
// 前置: f_scaled = δ^(r-1)·g
//       scaled_factors = v₁,...,vᵣ (所有 lc = δ 的单变量因子)
//       lc_targets = τ₁,...,τᵣ (各因子的 lc 目标, ∏τᵢ = δ^(r-1)·L)
//       ∏ vᵢ = f_scaled(x₁, α)
// 后置: 返回 G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]
//       lc(Gᵢ, x₁) = τᵢ, 调用方需对 pp(Gᵢ) 做试除验证
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>>
__multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_targets,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var);
```

### 6.8 需要新增的辅助函数

| 函数 | 用途 | 实现依赖 |
|---|---|---|
| `__taylor_coeff(f, xₖ, αₖ, j)` | 提取 Taylor 系数 | `pair_vec_div` + `assign` |
| `__upoly_gcd_extended(s, t, c, a, b)` | Z[x₁] pseudo-XGCD: s·a+t·b=c (ZZ 重载) | polynomial_gcd.hh；Zp 版本同时迁入 |
| `__poly_prem_univar(f, g, x₁)` | 多变量 f 对单变量 g 关于 x₁ 做伪余 | `prem` 的薄封装（lex 首变量；见 §6.5 说明）|

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

算法分两层：外层做无平方分解 + 分派，内层是 Wang 核心。

```
__factor_multivar(f_input):

1.  // 无平方分解（内部自动提取 content + Yun）
    sqf ← squarefreefactorize(f_input)           // [(g₁,m₁), (g₂,m₂), ...]
    // 每个 gₖ 无平方；来自 content 的因子不含 x₁，来自 Yun 的因子关于 x₁ 本原

    all_factors ← []
    content ← 1

    for (gₖ, mₖ) in sqf:
        if is_number(gₖ):
            content ← content · gₖ^mₖ
        else if get_variables(gₖ).size() == 1 or !contains(gₖ, x₁):
            // 单变量 或 不含 x₁ 的多变量 → 递归 factorize（变量数递减）
            sub ← factorize(gₖ)
            content ← content · sub.content^mₖ
            for (fᵢ, eᵢ) in sub.factors:
                all_factors.push(fᵢ, eᵢ · mₖ)
        else:
            // 多变量 + 含 x₁ + 无平方 + 本原 → Wang 核心
            wang_factors ← __wang_core(gₖ)
            for (fᵢ, eᵢ) in wang_factors:
                all_factors.push(fᵢ, eᵢ · mₖ)  // eᵢ = 1（gₖ 无平方）

    return {content, all_factors}
```

```
__wang_core(g):
    // 前置: g ∈ Z[x₁,...,xₙ], n ≥ 2, g 本原无平方
    x₁ ← get_variables(g).front().first

    retry_count ← 0

1.  // 选取值点
    eval ← __select_eval_point(g, x₁)             // 满足条件 (a)-(d)

2.  // 单变量分解
    f₀ ← assign(g, eval)                          // f₀ ∈ Z[x₁]
    uni_fac ← factorize(f₀)                       // 单变量 M4
    u₁,...,uᵣ ← uni_fac.factors                   // 本原, lc>0（非首一）

3.  if r ≤ 1:
        return [(g, 1)]                            // g 不可约

4.  // 首项系数校正 (§5)
    lc_result ← __wang_leading_coeff(g, [u₁,...,uᵣ], eval, x₁)
    if !lc_result.success:
        retry_count++
        if retry_count ≥ MAX_RETRY (=10): throw "Wang LC distribution failed"
        goto 1

5.  // 多变量 Hensel 提升 (§6)
    mv_factors ← __multivar_hensel_lift(
        lc_result.f_scaled,                        // δ^(r-1) · g
        lc_result.scaled_factors,                  // v₁,...,vᵣ
        lc_result.lc_targets,                      // τ₁,...,τᵣ
        eval, x₁)

6.  // 试除验证 + 去缩放
    verified ← []
    g_remaining ← g
    for G in mv_factors:
        h ← pp(G)                                    // 本原化，消去缩放因子 δ
        q, r ← divmod(g_remaining, h)              // 对 g（非 f_scaled!）试除
        if r = 0:
            verified.push(h)
            g_remaining ← q

    if deg(g_remaining, x₁) > 0:
        verified.push(pp(g_remaining))

    // 自洽性检查
    if ∏ verified ≠ g / (整数常数):
        retry_count++
        if retry_count ≥ MAX_RETRY: throw "Wang factorization failed"
        goto 1

    return [(h, 1) for h in verified]
```

### 7.2 递归终止性

`__factor_multivar` 的递归调用路径：

1. **步骤 1**：`squarefreefactorize` 内部调用 `factorize(cont)`，变量数严格递减。
2. **步骤 1**：不含 x₁ 的因子调用 `factorize(gₖ)`，变量数严格递减。
3. **__wang_core 步骤 4**：`factorize(L)` 其中 `L = lc(g, x₁) ∈ Z[x₂,...,xₙ]`，变量数严格递减。

所有递归都以变量数严格递减，基础情况是 n = 1（由 M4 处理）或 n = 0（常数）。

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

现有 `squarefreefactorize`（`polynomial_gcd.hh:101`）能正确处理多变量多项式：

1. 先提取 `cont(F, x₁)` 并递归调用自身（变量数递减，处理不含 x₁ 的重因子）
2. 对本原部分做 Yun 算法：`gcd(F_, derivative(F_, x₁))`

这两步组合可以检测**所有方向**的重因子——不含 x₁ 的在步骤 1 处理，
含 x₁ 的在步骤 2 处理（因为 g² | F_ 且 g 含 x₁ ⇒ g | derivative(F_, x₁)）。

`__factor_multivar` 以 `squarefreefactorize` 作为第一步（§7.1），
其返回的因子自动分为 content 部分（不含 x₁）和本原无平方部分（含 x₁）。
后者直接传入 `__wang_core`，满足 Wang 的无平方前置条件。

---

## 10. 实现路线图

| 阶段 | 内容 | 新增函数 | 依赖 |
|---|---|---|---|
| **Phase 5** | M5: 多变量 Wang | `pp` (polynomial_gcd.hh), `__upoly_gcd_extended` (ZZ 重载, polynomial_gcd.hh), `__taylor_coeff`, `__poly_prem_univar`, `__select_eval_point`, `__wang_leading_coeff` (含 `__wang_lc_result` + τᵢ), `__multivar_hensel_lift` (含 `__hensel_lift_one_var`), `__factor_multivar`, `factorize` 多变量 dispatch | M4 (已实现) |
| **Phase 6** | 增强：van Hoeij 重组 | `__factor_recombine_van_hoeij` + LLL 实现 | M3 替换 |
| **Phase 7** | 增强：Zippel 后备 | 稀疏插值模块 + Zippel 算法 | Phase 5 后备 |
| **Phase 8** | 终极：MTSHL | 二变量 Hensel 提升 + 稀疏插值驱动的多变量分解 | 替换 Phase 5 |

### 10.1 测试计划

| 可独立测试的函数 | 验证方法 |
|---|---|
| `pp(f)` | 验证 `cont(f) · pp(f) == f` |
| `__taylor_coeff(f, xₖ, αₖ, j)` | 验证 `Σ coeff_j · (xₖ-αₖ)^j == f` |
| `__upoly_gcd_extended` (ZZ 重载) | 验证 `s·a + t·b = c`（c ∈ Z），s, t ∈ Z[x₁]，`Σ sᵢ·V̂ᵢ = denom` |
| `__poly_prem_univar(f, g, x₁)` | 验证 `δ^k·f = q·g + r`，`deg(r,x₁) < deg(g)` |
| `__select_eval_point` | 验证条件 (a)-(d)：无平方、lc 非零、度数守恒、lc 因子可分辨 |
| `__wang_leading_coeff` | 构造已知分解的多项式，验证 `∏σᵢ = L`，`∏τᵢ = δ^(r-1)·L`，`τᵢ(α) = δ`，`∏vᵢ = f_scaled(x₁,α)` |
| `__multivar_hensel_lift` | 二变量多项式提升后 `∏Gᵢ = f_scaled`，`pp(Gᵢ) \| g` |
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
// §3.1 多变量本原部分 (polynomial_gcd.hh, cont() 伴侣)
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
pp(const polynomial_<ZZ, lex_<var_order>>& f);

// §6.4 Taylor 系数提取
template<class var_order>
polynomial_<ZZ, lex_<var_order>> __taylor_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& xk, const ZZ& alpha_k, int j);

// §6.8 Z[x₁] pseudo-XGCD: s·a + t·b = c (c ∈ Z)
// 位置: polynomial_gcd.hh（与 Zp 版本统一归入 GCD 模块）
// Zp 版本: __upoly_gcd_extended(s, t, a, b)       → s·a + t·b = 1  (从 polynomial_factorize.hh 迁入)
// ZZ 版本: __upoly_gcd_extended(s, t, c, a, b)    → s·a + t·b = c, c ∈ Z  (新增)
inline void __upoly_gcd_extended(
    upolynomial_<ZZ>& s, upolynomial_<ZZ>& t, ZZ& c,
    const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b);

// §6.8 多变量 f 对单变量 g 关于 x₁ 做伪余（prem 的薄封装）
template<class var_order>
polynomial_<ZZ, lex_<var_order>> __poly_prem_univar(
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
    std::vector<polynomial_<ZZ, lex_<var_order>>> lc_targets;       // τ₁,...,τᵣ
    std::vector<upolynomial_<ZZ>> scaled_factors;                   // v₁,...,vᵣ
};

// §5 首项系数校正 (所有输入 const, 输出通过返回值)
template<class var_order>
__wang_lc_result<var_order> __wang_leading_coeff(
    const polynomial_<ZZ, lex_<var_order>>& f,
    const std::vector<upolynomial_<ZZ>>& univar_factors,
    const std::map<variable, ZZ>& eval_point, const variable& main_var);

// §6 单变量 Hensel 提升步 (内层循环)
template<class var_order>
void __hensel_lift_one_var(
    const polynomial_<ZZ, lex_<var_order>>& f_curr,
    std::vector<polynomial_<ZZ, lex_<var_order>>>& G,
    const std::vector<upolynomial_<ZZ>>& bezout_s,
    const ZZ& bezout_denom,
    const std::vector<upolynomial_<ZZ>>& v_factors,          // vᵢ (lc=δ, 非首一)
    const ZZ& delta,                                          // δ
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_tau,    // τᵢ
    const variable& xk, const ZZ& alpha_k, int dk);

// §6 多变量 Hensel 提升 (外层入口)
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>> __multivar_hensel_lift(
    const polynomial_<ZZ, lex_<var_order>>& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::vector<polynomial_<ZZ, lex_<var_order>>>& lc_targets,    // τᵢ
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
| `prem(f, g, x)` | `polynomial_gcd.hh` | 伪余式 | MDP 中 vᵢ 非首一时使用 |
| `__upoly_gcd_extended(s, t, a, b)` | `polynomial_gcd.hh` | Zₚ 上扩展 GCD | 从 polynomial_factorize.hh 迁入；ZZ 重载 `(s, t, c, a, b)` 新增 |
| `is_number(f)` | `upolynomial.hh` | 常数检测 | |
| `poly_convert(in, out [, var])` | `upolynomial.hh` | polynomial ↔ upolynomial 双向转换 | upolynomial→polynomial 需传 `variable` 参数 |
