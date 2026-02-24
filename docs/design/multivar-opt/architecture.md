# MTSHL-d 多变量 Hensel 提升：架构文档

> 阶段：架构（workflow.md §2.2）
> 调研基础：`docs/research/multivar-mtshl-research.md`（5 篇 Monagan-Tuncer 论文综合）
> 目标文件：`clpoly/polynomial_factorize_wang.hh`

---

## 1. 核心流程

### 1.1 上下文：MTSHL-d 在 Wang 框架中的位置

```
__wang_core(f)
    ├─ __select_eval_point          [保留不变]
    ├─ factorize(f|_{x2=α2,...})    [保留不变，得到单变量因子 f0[i]]
    ├─ __wang_leading_coeff         [保留不变，LC 校正后得到 scaled_factors[i]]
    ├─ 【替换】__mtshl_lift         [替换 __multivar_hensel_lift]
    └─ Zassenhaus 重组              [保留不变]
```

### 1.2 __mtshl_lift 端到端数据流

```
输入：
  f_scaled          ∈ Z[x1,...,xn]           经 LC 校正的目标多项式
  scaled_factors[i] ∈ Z[x1]，i=1..r          单变量因子（LC 校正后），满足 ∏scaled_factors[i] = f_scaled|_{x2=α2,...}
  eval_point        = {x2=α2,...,xn=αn}     主求值点（ideal alphas，固定不变）
  main_var          = x1
  p                 ∈ 机器素数               满足 p > 2·Mignotte_bound(f_scaled)

阶段 A：初始化（模素数）
  F[i] ← polynomial_mod(scaled_factors[i], p)     ∈ Zp[x1]

  注：forms[i] 由 __mtshl_step_j 在每步入口自行初始化（见模块 B 规约：__mtshl_step_j），
      __mtshl_lift 不负责管理 forms。

阶段 B：逐变量提升，j = 2,...,n
  aj ← f_scaled|_{x_{j+1}=α_{j+1},...,xn=αn}（mod p）  ∈ Zp[x1,...,xj]

  __mtshl_step_j(aj, F, αj, xj, p)
  → 对 F[i] 做 in-place 更新
  → 循环后不变量：∏F[i] = aj in Zp[x1,...,xj]

阶段 C：系数恢复
  ∏F[i] = f_scaled (mod p) in Zp[x1,...,xn]
  G[i] ← symmetric_mod(F[i], p)    ∈ Z[x1,...,xn]
  （单机器素数 p > 2·Mignotte_bound 保证对称约化唯一性）

输出：{G[i]}，满足 ∏G[i] = f_scaled in Z[x1,...,xn]
```

**关于机器素数 p**：选取满足 `p > 2·Mignotte_bound(f_scaled)` 的单个机器素数（uint32_t），
全程使用 64-bit 整数运算，最终一次对称约化即可恢复 Z 系数。
若需要更大的系数界，CASC 2018 提供 p-adic 后处理方案，但属于扩展功能，不在本阶段（M1–M5）实现范围内。

### 1.3 __mtshl_step_j 内部流程

```
入口：
  forms[i] ← Supp(F[i])           （JSC 2020 Alg.5 line 2：每步入口从当前因子支撑初始化）

error ← aj - ∏F[i]

for k = 1, 2, ... while error ≠ 0:    （JSC 2020 Alg.5；CASC 2016 Algorithm 4；上界 deg(aj,xj) 由 Hensel 提升终止定理保证）
    ck ← Taylor 系数：[(xj - αj)^k] in error   ∈ Zp[x1,...,x_{j-1}]
    if ck == 0: continue

    【j=2 路径】若 j==2：F[i] ∈ Zp[x1]，ck ∈ Zp[x1]，为单变量多因子 Diophantine。
               θ-array/Vandermonde 框架在此退化（无辅助变量，θ_t ≡ 空乘积 = 1，
               Vandermonde v[l] = Σ c_t·1^l = Σ c_t 为常数，矩阵奇异，无法恢复各项系数）。
               必须直接在 Zp[x1] 中求解多因子 Diophantine（CASC 2018 §3 multi-BDP 基础情形）：
               - r=2：单次 EEA（ICMS 2018 HenselLift1 Algorithm 2 第 1 步）
               - r>2：r-1 次 EEA 逐对归约（先将 r 因子问题归约为 2 因子，再递归求解），
                       代价 O(r·d1²)，与 multi-BDP 在 j=2 时的基础情形一致
               success ← 多因子 Zp[x1] Diophantine（复用现有单变量 Diophantine 代码，支持 r≥2）

    【j≥3 主路径】若 j≥3：
        【主路径】success ← __sparse_int(F, ck, forms, p, σk)
        【重试】  if !success:
                      重新选取 sparse_betas，再调用一次 __sparse_int
                      （JSC 2020 Alg.5：失败多因随机点碰撞，换点即可；避免不必要地升级到指数回退）
        【回退 j=3】if !success and j==3:
            success ← __multi_bdp(F, ck, p, σk)
        【回退 j>3】if !success and j>3:
            success ← __wmds(F, ck, p, σk)         （Wang 递归 MDP，现有实现复用/改写）
    if !success: return false

    for i: F[i] += σk[i] · (xj - αj)^k
    forms[i] ← Supp(σk[i])        （JSC 2020 Alg.5 line 7：更新骨架；j=2 时此行无实质作用，因 EEA 路径不使用 forms）
    error ← aj - ∏F[i]            （全量重算；Taylor 系数 ck 的提取可用 ICMS 2018 增量公式加速，见细化阶段）

return true
```

---

## 2. 模块划分与功能规约

本节采用 workflow.md §3.1 规约格式。

类型别名（全文统一）：
- `Poly`   = `polynomial_<ZZ, lex_<var_order>>`
- `PolyZp` = `polynomial_<Zp, lex_<var_order>>`（prime 嵌入 Zp 类型）
- `Mono`   = `basic_monomial<lex_<var_order>>`
- `UPZp`   = `upolynomial_<Zp>`

---

### 模块 A：`__mtshl_lift`

**功能描述**：MTSHL-d 顶层驱动。将 Z[x1] 中的单变量因子 `{scaled_factors[i]}` 逐变量提升为 n 变量整系数因子 `{G[i]}`。

**前置条件（Requires）**：
- `f_scaled` ∈ Z[x1,...,xn]，squarefree，primitive
- `scaled_factors[i]` ∈ Z[x1]（数据结构：`upolynomial_<ZZ>`），由 `__wang_leading_coeff` 产出，满足：
  - `∏scaled_factors[i] = f_scaled(x1, α2,...,αn)` in Z[x1]
  - `lc(scaled_factors[i], x1)` 为 `lc(G[i], x1)` 在 `(α2,...,αn)` 处的求值（LC 已正确分配）
- `p` 为机器素数（uint32_t），`p ∤ lc(f_scaled) · ∏lc(scaled_factors[i])`（确保 LC mod p 非零）
- `p > 2 · Mignotte_bound(f_scaled)`（由 `__wang_core` 在选 p 时保证）
- `polynomial_mod(scaled_factors[i], p)` 两两互素（由 `__select_eval_point` 保证）

**后置条件（Ensures）**：
- 若正常返回：`{G[i]}` 满足 `∏G[i] = f_scaled` in Z[x1,...,xn]，`deg(G[i], x1) = deg(scaled_factors[i])`
- 若返回空：提升失败，`__wang_core` 须重新选求值点

**副作用**：无（纯函数，内部中间量不外泄）。

---

### 模块 B：`__mtshl_step_j`

**功能描述**：MTSHL-d 第 j 步。将 `Zp[x1,...,x_{j-1}]` 中的因子提升到 `Zp[x1,...,xj]`，使其乘积与 `aj` 一致。

**前置条件（Requires）**：
- `aj` ∈ Zp[x1,...,xj]，等于 `f_scaled` 在 `{x_{j+1}=α_{j+1},...}` 处的求值
- `F[i]` ∈ Zp[x1,...,x_{j-1}]
- `∏F[i] = aj|_{xj=αj}` in Zp[x1,...,x_{j-1}]（由上一步或初始化保证）
- `F[i] mod p` 两两互素

**后置条件（Ensures）**：
- 若返回 `true`：`∏F[i] = aj` in Zp[x1,...,xj]
- 若返回 `false`：`__wang_core` 须重新选求值点

**不变式（Invariants）**：
- 循环开始前（k=0）及每次 Taylor 迭代后：`∏F[i] ≡ aj (mod (xj - αj)^{k+1})` in Zp[x1,...,xj]
  （k=0 初值由 Requires 中 `∏F[i] = aj|_{xj=αj}` 保证，即 `∏F[i] ≡ aj (mod (xj - αj)^1)`）

**副作用**：in-place 更新 `F[i]`。

**注**：`forms[i]` 为模块内部状态，在入口处初始化为 `Supp(F[i])`，不作为参数暴露给调用方。

**注（j=2 路径）**：当 j=2 时，F[i] ∈ Zp[x1]，MDP 为单变量 Diophantine，直接用 EEA 求解（不调用 `__sparse_int`）；`forms[i]` 在 j=2 时不使用。

---

### 模块 C：`__sparse_int`

**功能描述**：用稀疏插值求解多因子 MDP：在 Zp[x1,...,x_{j-1}] 中求 `{σi}` 使 `Σ σi·bi = c`，`bi = ∏_{l≠i}F[l]`。

每次调用时在 Zp 中**内部随机生成** `sparse_betas = (β2,...,β_{j-1})`（βk ≠ 0，各变量独立，每次 MDP 独立选取），用于 θ-array 求值与 Vandermonde 恢复，与 ideal_alphas 无关。

**前置条件（Requires）**：
- `F[i]` ∈ Zp[x1,...,x_{j-1}]，j ≥ 3（j=2 时 `__mtshl_step_j` 直接用 EEA，不调用本模块），两两互素 mod p
- `c` ∈ Zp[x1,...,x_{j-1}]，`deg(c, x1) < Σ_i deg(F[i], x1)`
- `forms[i]` = 期望 `Supp(σi)`：k=1 时为 `Supp(F[i])`（步入口初始化），k>1 时为 `Supp(σ_{k-1,i})`（前一 Taylor 迭代更新）

**后置条件（Ensures）**：
- 若返回 `true`：`Σ result[i]·bi = c` in Zp[x1,...,x_{j-1}]，`Supp(result[i]) ⊆ forms[i]`，且 `deg(result[i], x1) < deg(F[i], x1)`
- 若返回 `false`（两种情形，均不保证 result 内容）：
  1. 真实 `Supp(σi) ⊄ forms[i]`（骨架假设不成立）
  2. Vandermonde 矩阵奇异（存在 θ_t 值碰撞，sparse_betas 随机选取恰好不当）

**副作用**：无（`result` 为输出参数）。

---

### 模块 D1：`__multi_bdp`

**功能描述**：稠密双变量多因子 Diophantine 求解（`__sparse_int` 在 j=3 时的回退）。
在 Zp[x1,x2] 中求 `{σi}` 使 `Σ σi·bi = c`，`bi = ∏_{l≠i}F[l]`。

**前置条件（Requires）**：
- `F[i]` ∈ Zp[x1,x2]（双变量，j=3 时因子恰好在此空间）
- `c` ∈ Zp[x1,x2]，`deg(c, x1) < Σ_i deg(F[i], x1)`
- `F[i]` 两两互素 mod p

**后置条件（Ensures）**：
- 若返回 `true`：`Σ result[i]·bi = c` in Zp[x1,x2]，且 `deg(result[i], x1) < deg(F[i], x1)`
- 若返回 `false`：极罕见（某个 x1 求值点处因子共根），上层重新选点

**副作用**：无。

**注**：仅在 j=3 时调用。j>3 时 SparseInt 失败的回退为模块 D2。

---

### 模块 D2：`__wmds`（Wang 递归 MDP，j>3 时的回退）

**功能描述**：Wang 多维 Diophantine 求解器（WMDS）。处理 j>3 时 SparseInt 失败的回退，对任意维度均可求解。

**前置条件（Requires）**：
- `F[i]` ∈ Zp[x1,...,x_{j-1}]，j > 3（j=3 时用 `__multi_bdp`）
- `c` ∈ Zp[x1,...,x_{j-1}]，`deg(c, x1) < Σ_i deg(F[i], x1)`
- `F[i]` 两两互素 mod p

**后置条件（Ensures）**：
- 若返回 `true`：`Σ result[i]·bi = c` in Zp[x1,...,x_{j-1}]，且 `deg(result[i], x1) < deg(F[i], x1)`
- 若返回 `false`：退化情形，上层重新选点

**副作用**：无。

**注**：WMDS 本质上是 `__multivar_diophantine` 的 Zp 版本，可在改写时从现有代码演化。

---

### 模块 E：`__theta_array_eval`

**功能描述**：利用 θ-array 技术，将 `f` ∈ Zp[x1,...,x_{j-1}] 在 s 个几何点处批量求值。
第 l 号求值点（l=1,...,s）为 `(x2=β2^l, x3=β3^l,...,x_{j-1}=β_{j-1}^l)`，其中 β2,...,β_{j-1} 为各变量的独立随机基。
输出 s 个单变量像 `f(x1, β2^l,...,β_{j-1}^l)`，l=1,...,s。

**θ-array 原理**（CASC 2016 §3.3）：对 f 中每个单项 `c·x1^d1·x2^d2·...·x_{j-1}^{d_{j-1}}`，预计算
`θ_term = ∏_{k=2}^{j-1} βk^{d_k}`，则第 l 号求值点处该项贡献为 `c·x1^d1·θ_term^l`。
迭代更新：`eval_l = eval_{l-1} · θ_term`，代价 O(s·t)（t = |Supp(f)|），相比朴素求值节省因子 (j-2)。

这里的 `β2,...,β_{j-1}` 是 `__sparse_int` 内部为本次 MDP 新选的随机基（**sparse_betas**），与 eval_point 中的 `α2,...,α_{j-1}`（**ideal_alphas**）无关。

**前置条件（Requires）**：
- `f` ∈ Zp[x1,...,x_{j-1}]，j ≥ 3（j=2 时 θ_t ≡ 1 导致 Vandermonde 奇异，本模块不适用；j=2 由 `__mtshl_step_j` 直接用 EEA 处理）
- `sparse_betas = (β2,...,β_{j-1})`，`βk ≠ 0` in Zp，j-2 个独立非零随机基（由 `__sparse_int` 生成后作为参数传入，与 ideal_alphas 无关）
- `s = max_i(|forms[i]|)`，`s ≥ 1`（由 `__sparse_int` 根据各因子骨架大小计算后作为参数传入）

**后置条件（Ensures）**：
- `images[l] = f(x1, β2^l,...,β_{j-1}^l)` ∈ Zp[x1]，l=1,...,s

**副作用**：无。

---

### 模块 F：`__vandermonde_solve`

**功能描述**：从稀疏多项式在 s 个几何点处的求值恢复其系数。
已知 `v[l] = Σ_t c_t · θ_t^l`（l=1,...,s），恢复 `{c_t}`，其中 `θ_t = ∏_k βk^{e_{t,k}}`。

**前置条件（Requires）**：
- `s = |terms|`（插值点数等于待恢复项数）
- `θ_t` 两两不同（由 `βk` 随机选取以高概率保证，Vandermonde 矩阵非奇异）
- `θ_t ≠ 0`

**后置条件（Ensures）**：
- 返回 `{c_t}` 满足 `Σ_t c_t · θ_t^l = v[l]`，l=1,...,s

**副作用**：无。

---

## 3. 模块间接口规约

本节采用 workflow.md §3.2 规约格式。

---

### 接口 0：`__wang_core` → `__mtshl_lift`

**输入数据**：
- `f_scaled`: Poly，LC 校正后的目标多项式
- `scaled_factors[]`: `upolynomial_<ZZ>[]`，单变量因子（由 `__wang_leading_coeff` 保证已有正确 LC）
- `eval_point`: `map<variable, ZZ>`，主求值点 {x2=α2,...}
- `main_var`: variable，= x1
- `p`: uint32_t，机器素数

**注**：`lc_targets` 不传入 `__mtshl_lift`。设计依据：

多因子 MDP 次数约束（CASC 2018 §3）要求 `deg(σk,i, x1) < deg(F[i], x1)`。推论链：
1. `deg(σk,i·(xj−αj)^k, x1) = deg(σk,i, x1) < deg(F[i], x1)`（xj−αj 不含 x1，不增加 x1 次数）
2. 因此 `σk,i·(xj−αj)^k` 只贡献于 F[i] 的**低于最高次**的 x1 项
3. F[i] 的 x1 最高次项系数（= `lc(F[i], x1)`）不受任何低次修正影响，保持不变

故 LC 在整个提升过程中自动保持，无需显式校正。

（各 MTSHL 论文未独立陈述该性质，但由次数约束逻辑导出。若 M4 验收 `∏G[i] = f_scaled` 未通过，优先排查 MDP 次数约束是否在 multi-BDP、WMDS 所有回退情形中均被保证。）

**输出数据**：
- `vector<Poly>`：Z 系数多变量因子 {G[i]}；空表示失败

**协议约定**：
- 调用方（`__wang_core`）责任：
  - 保证 `p > 2 · Mignotte_bound(f_scaled)`
  - 保证 `polynomial_mod(scaled_factors[i], p)` 两两互素
  - `p ∤ lc(f_scaled) · ∏lc(scaled_factors[i])`（确保所有因子 LC 在 mod p 后非零）
- 被调用方（`__mtshl_lift`）责任：
  - 正常返回时满足 `∏G[i] = f_scaled` in Z[x1,...,xn]
  - 返回空时不修改任何外部状态

---

### 接口 A → B：`__mtshl_lift` → `__mtshl_step_j`

**输入数据**：
- `aj`: PolyZp，`f_scaled` 在 `{x_{j+1}=α_{j+1},...}` 处的 Zp 求值
- `F[]`: PolyZp[]，当前 j-1 变量因子（in-place 更新）
- `αj`: Zp，当前变量的 ideal alpha
- `xj`: variable，当前提升变量
- `p`: uint32_t

**输出数据**：
- `bool`：true = 成功；false = 求值点不适用
- `F[]`（修改）：提升后的 xj-变量因子

**协议约定**：
- 调用方（`__mtshl_lift`）责任：
  - 进入时保证 `∏F[i] = aj|_{xj=αj}` in Zp[x1,...,x_{j-1}]
  - `forms` 不由调用方传入，由 `__mtshl_step_j` 内部管理
- 被调用方（`__mtshl_step_j`）责任：
  - 返回 true 时保证 `∏F[i] = aj` in Zp[x1,...,xj]
  - 在入口处初始化 `forms[i] ← Supp(F[i])`，在每次 Taylor 迭代后更新 `forms[i] ← Supp(σk[i])`

---

### 接口 B → C：`__mtshl_step_j` → `__sparse_int`

**输入数据**：
- `F[]`: PolyZp[]，当前 j-1 变量因子（只读）
- `ck`: PolyZp，当前 Taylor 系数（MDP 右端项）
- `forms[]`: `vector<Mono>[]`，当前期望各 σi 的支撑（只读）
- `p`: uint32_t

**输出数据**：
- `bool`：true = 成功；false = 支撑假设不成立
- `result[]`: PolyZp[]，解 {σi}

**协议约定**：
- 调用方（`__mtshl_step_j`）责任：
  - `forms[i]` 在入口时由 `Supp(F[i])` 初始化，之后由 `Supp(σk[i])` 更新
  - `deg(ck, x1) < Σ_i deg(F[i], x1)`
- 被调用方（`__sparse_int`）责任：
  - **内部随机生成** `sparse_betas = (β2,...,β_{j-1})`（不由调用方传入，每次调用独立选取，与 ideal_alphas 无关）
  - 返回 true 时 `Σresult[i]·bi = ck`，`Supp(result[i]) ⊆ forms[i]`，且 `deg(result[i], x1) < deg(F[i], x1)`
  - 返回 false 时不保证 result 内容

---

### 接口 B → D1：`__mtshl_step_j` → `__multi_bdp`（j=3 时的回退）

**触发条件**：`__sparse_int` 首次失败后重试仍失败，且当前 j=3

**输入数据**：
- `F[]`: PolyZp[]，双变量因子（j=3 时 F[i] ∈ Zp[x1,x2]，直接可用）
- `ck`: PolyZp，MDP 右端项 ∈ Zp[x1,x2]
- `p`: uint32_t

**输出数据**：
- `bool`，`result[]`: PolyZp[]

**协议约定**：
- 调用方（`__mtshl_step_j`）责任：
  - 仅在 j=3 时调用（`__sparse_int` 失败且当前 j=3）
  - 保证 `F[i]` 两两互素 mod p
  - 保证 `deg(ck, x1) < Σ_i deg(F[i], x1)`
- 被调用方（`__multi_bdp`）责任：返回 true 时 `Σresult[i]·bi = ck` in Zp[x1,x2]，且 `deg(result[i], x1) < deg(F[i], x1)`

---

### 接口 B → D2：`__mtshl_step_j` → `__wmds`（j>3 时的回退）

**触发条件**：`__sparse_int` 首次失败后重试仍失败，且当前 j>3

**输入数据**：
- `F[]`: PolyZp[]，j-1 变量因子
- `ck`: PolyZp，MDP 右端项 ∈ Zp[x1,...,x_{j-1}]
- `p`: uint32_t

**输出数据**：
- `bool`，`result[]`: PolyZp[]

**协议约定**：
- 调用方（`__mtshl_step_j`）责任：
  - 仅在 j>3 时调用（`__sparse_int` 失败且维度高于 3）
  - 保证 `F[i]` 两两互素 mod p
  - 保证 `deg(ck, x1) < Σ_i deg(F[i], x1)`
- 被调用方（`__wmds`）责任：返回 true 时 `Σresult[i]·bi = ck` in Zp[x1,...,x_{j-1}]，且 `deg(result[i], x1) < deg(F[i], x1)`

---

### 接口 C → E：`__sparse_int` → `__theta_array_eval`

**输入数据**：
- `f`: PolyZp，待求值多项式，`f` ∈ Zp[x1,...,x_{j-1}]（含 x1 及辅助变量 x2,...,x_{j-1}）
- `sparse_betas`: Zp[]，`(β2,...,β_{j-1})`，由 `__sparse_int` 内部生成的本次 MDP 随机基，非零
- `s`: int，求值点数 = `max_i(|forms[i]|)`

**输出数据**：
- `images[]`: UPZp[]，s 个单变量像

**协议约定**：
- 调用方（`__sparse_int`）责任：`sparse_betas = (β2,...,β_{j-1})` 为本次内部新选随机值（各变量独立基），`s = max_i(|forms[i]|)`；仅在 j≥3 时调用（j=2 路径由 `__mtshl_step_j` 直接用 EEA 处理，不经过 `__sparse_int`）
- 被调用方（`__theta_array_eval`）责任：`images[l] = f(x1, β2^l,...,β_{j-1}^l)` ∈ Zp[x1]，l=1,...,s

---

### 接口 C → F：`__sparse_int` → `__vandermonde_solve`

**输入数据**：
- `values[]`: Zp[]，s 个求值点处的系数值 `v[l]`
- `thetas[]`: Zp[]，各项的几何公比 `θ_t = ∏_k βk^{e_{t,k}}`
- `s`: int

**输出数据**：
- `coeffs[]`: Zp[]，各项系数 `{c_t}`

**协议约定**：
- 调用方（`__sparse_int`）责任：`thetas` 两两不同（由 `sparse_betas` 随机选取以高概率保证）
- 被调用方（`__vandermonde_solve`）责任：返回满足 Vandermonde 方程的唯一解

---

## 4. 关键设计决策

### 决策 1：全程使用单个机器素数 p（不用 p-adic 提升）

**结论**：选取满足 `p > 2·Mignotte_bound(f_scaled)` 的单个素数 p，全程在 Zp 工作，最终一次对称约化恢复 Z 系数。

**理由**：标准 Hensel 提升理论（Mignotte 界）保证：若 `p > 2·Mignotte_bound(f_scaled)`，则 Z 系数因子可由 Zp 系数唯一对称约化恢复，无需多素数 CRT 或 p-adic 迭代。（注：CASC 2018 Algorithm 5 是处理 Mignotte 界超出机器字长时的 p-adic 提升方案，与此处的单素数充分性无关；p-adic 扩展不在 M1–M5 范围内。）

**实现约束**：CASC 2018 推荐 63-bit 素数（uint64_t），使 Zp 乘法可用 `__int128` 完成。CLPoly 当前 `Zp` 类为 uint32_t 实现，乘法在 int64_t 内可行（p < 2^32，p² < 2^64）。

两者权衡：
- uint32_t：实现简单，但 Mignotte 界若超过约 2×10⁹ 则无法找到合适素数，需回落 p-adic（超出 M1–M5 范围）
- uint64_t：覆盖更大的 Mignotte 界，但需修改 `Zp` 类，属独立工程任务

**M1–M5 决策**：沿用 uint32_t `Zp`；若 `p > 2·Mignotte_bound` 无解，报告失败并由 `__wang_core` 重选（罕见场景推迟处理）。

**影响**：`__wang_core` 选 p 时需验证 `p > 2·Mignotte_bound`；若当前 uint32_t 素数集中找不到满足条件的 p，返回失败。

---

### 决策 2：ideal_alphas 与 sparse_betas 双重求值点体系

**结论**：维护两套独立的求值点：
- **ideal_alphas** = eval_point 中的 {α2,...,αn}：仅用于 **Taylor 系数提取**（将 f_scaled 逐步代入，固定，贯穿整个提升过程）
- **sparse_betas** = 每次 `__sparse_int` 调用时在 Zp 中新选的随机基 β2,...,β_{j-1}：用于 **θ-array 求值和 Vandermonde 恢复**（每次 MDP 独立选取）

**理由**：CASC 2016 §3.3 明确了 θ-array 中各 βk ∈ Zp 独立随机选取（非零）；ICMS 2018 Algorithm 3 明确在"s 个随机点 (β₁k,...)"处求值，区别于 ideal αk；CASC 2018 Algorithm 3 (SparseInt) 在函数内部"选随机点 β ∈ Zp"，与 ideal α 参数分开。混用两套求值点会导致：① ideal_alphas 与因子的 LC/根有关（由 `__select_eval_point` 选定），若复用为 θ-array 基，θ_t 值会因多项式结构而产生规律性碰撞，使 Vandermonde 矩阵奇异；② 稀疏插值的随机性假设失效，插值结果与真实 σi 不一致。θ-array 的 θ_term = ∏k βk^{d_k} 中的 βk 是 sparse_betas，不是 ideal_alphas。

---

### 决策 3：forms[] 为 `__mtshl_step_j` 的内部状态

**结论**：`forms[i]`（σ 的支撑骨架）在 `__mtshl_step_j` 入口初始化为 `Supp(F[i])`，在每次 Taylor 迭代 k 后更新为 `Supp(σk[i])`；不作为参数在 `__mtshl_lift` 与 `__mtshl_step_j` 之间传递。

**理由**：JSC 2020 Algorithm 5 第 2、7 行明确了这一生命周期：骨架在每个 j 步入口重新初始化，在步内随迭代收缩。将其暴露为参数会增加接口复杂性且无收益。

---

### 决策 4：回退策略分级（j=3 用 multi-BDP，j>3 用 WMDS）

**结论**：
- 所有 j：`__sparse_int` 失败时，先用**新 sparse_betas 重试一次**（JSC 2020 Alg.5）——SparseInt 失败主因是随机点碰撞，换点即可解决，无需立即升级到指数回退。
- j=3（MDP 在 Zp[x1,x2]）：重试仍失败时回退到 `__multi_bdp`（CASC 2018 §3 multi-BDP，代价 O(r·d1·d2²)）
- j>3（MDP 在 Zp[x1,...,x_{j-1}]）：重试仍失败时回退到 `__wmds`（Wang 递归 MDP 的 Zp 版）

**理由**：CASC 2018 §3 明确定义了 j=3 的 multi-BDP 回退；对于 j>3，论文未指定回退策略，但以下三点推论确定 WMDS 是正确选择：
1. multi-BDP 是纯双变量（Zp[x1,x2]）求解器，其设计本身不支持 j>3 的高维 MDP；
2. JSC 2020 Algorithm 5 在提及 WMDS 时未加维度限制，表明其作为通用求解器的适用性；
3. WMDS 本质上是递归 MDP 求解器（`__multivar_diophantine` 的 Zp 版），维度 j-1 对算法正确性无约束，只影响复杂度（O(T^{j-1})）。

"j>3 用 WMDS"是本架构的确定设计决策（非推测）。若 M4 验收时发现 WMDS 正确性有问题，应立即报告为 bug 并分析根因，而非在此修改架构决策。WMDS 可从现有 `__multivar_diophantine` 演化为 Zp 版本。

---

### 决策 5：不做 M0（跳过模 Bézout 快速里程碑）

**结论**：直接实现 MTSHL-d，不做过渡性的模 Bézout 修改。

**理由**：M0 的 Diophantine 解的系数界（pa 选取）正确性无把握；M0 代码在 MTSHL-d 完成后全部删除，投入产出比低。

---

## 5. 与现有代码的边界

### 保留（复用）

| 函数 | 位置 | 说明 |
|------|------|------|
| `__wang_core` | L1063 | 外层控制，需做以下修改：①将 `__multivar_hensel_lift` 调用替换为 `__mtshl_lift`；②在选 p 时新增 Mignotte 界计算，确保 `p > 2·Mignotte_bound(f_scaled)`；③调用 `__mtshl_lift` 时不传 `lc_targets`（已不需要） |
| `__select_eval_point` | L60 | 保证 `f0[i] mod p` 两两互素 |
| `__wang_leading_coeff` | L188 | LC 校正 |
| `__factor_multivar` | L1253 | 入口 + 无平方分解 |
| `__taylor_coeff` | — | Taylor 系数提取，`__mtshl_step_j` 复用 |
| `__symmetric_mod` | — | 对称约化，阶段 C 复用 |
| Zassenhaus 重组 | L1175 | 因子重组，MTSHL 不负责 |

### 删除

| 函数/结构 | 删除原因 |
|---------|---------|
| `__multivar_hensel_lift`（L824–971） | 被 `__mtshl_lift` 完全替换 |
| `__hensel_lift_one_var`（L698–818） | 被 `__mtshl_step_j` 替换 |
| `__multivar_diophantine`（L474–664） | 被 `__sparse_int` + `__wmds` 替换 |
| Bézout 链（bezout_s, denom）（L836–865） | MTSHL 全程在 Zp，无 Bézout 链 |
| `__hensel_lc_correct` | Wang Z 提升特有；MTSHL-d MDP 次数约束天然保持 LC，无需显式校正 |

### 新增（全部在 `polynomial_factorize_wang.hh`）

| 函数 | 说明 | 代码量估计 |
|------|------|-----------|
| `__vandermonde_solve` | Vandermonde 系数恢复 | ~50 行 |
| `__theta_array_eval` | θ-array 批量求值 | ~50 行 |
| `__multi_bdp` | 双变量多因子 Diophantine（j=3 回退） | ~80 行 |
| `__wmds` | Wang 递归 MDP 的 Zp 版（j>3 回退） | ~80 行（从现有 `__multivar_diophantine` 演化） |
| `__sparse_int` | 稀疏插值 MDP 求解器 | ~150 行 |
| `__mtshl_step_j` | 单变量提升步骤 | ~80 行 |
| `__mtshl_lift` | 顶层提升驱动 | ~80 行 |

---

## 6. 实施里程碑

```
【M1】__vandermonde_solve + __theta_array_eval（基础设施，~100 行）
    验收：单元测试——给定几何点求值能正确恢复稀疏多项式系数

【M2】__multi_bdp（双变量多因子 Diophantine，~80 行）
    验收：单元测试——Σσi·bi = c 验证（r=2,3,4，多组随机输入）

【M3】__sparse_int（稀疏 MDP 求解器，~150 行）
    验收：单元测试——与 __multi_bdp 结果交叉验证（j=3，稀疏输入时结果一致）

【M4】__mtshl_step_j + __mtshl_lift + __wmds（主提升循环，~240 行）
    验收：test_multivar_hensel 全部通过

【M5】接入 __wang_core（替换调用，~20 行）
    验收：test_factorize_multivar 全部通过；test_stress_factorize 通过
    性能：bivar-70 < 20ms；bench-all 无退化

【M6】性能调优 + 基准存档
    验收：make bench-save；与 FLINT 对比 ≤ 5x
```

---

## 7. 验收标准

| 用例 | 当前 | 目标 | FLINT 参考 |
|------|------|------|-----------|
| bivar 70 factors | ~65s | < 20ms | ~9ms |
| trivar 60 factors | ~6s | < 10ms | ~1.3ms |
| 14var 100term（稀疏） | N/A | < 3s | ~2.4s（Maple） |
| bivar deg5*deg5 | 0.83ms | < 2ms | ~0.33ms |
| test_factorize_multivar | PASS | PASS | — |
| crosscheck（FLINT/NTL） | PASS | PASS | — |

---

## 参考文献

| 文献 | 对应模块 | 关键贡献 |
|------|---------|---------|
| Monagan-Tuncer CASC 2016 | `__sparse_int`（r=2）；`__theta_array_eval` | 2-因子 SparseInt；θ-array 批量求值 |
| Monagan-Tuncer ICMS 2018 | `__mtshl_step_j`（j=2 EEA 路径） | HenselLift1；误差增量更新；sparse_betas 独立于 ideal_alphas |
| Monagan-Tuncer CASC 2018 | `__sparse_int`（r>2）；`__mtshl_lift`；`__multi_bdp` | Algorithm 3 多因子 SparseInt；MTSHL-d；multi-BDP；j=2 多因子 Diophantine |
| Monagan-Tuncer JSC 2020 | 验收标准；`__wmds` 回退策略 | Theorem 19；forms[] 生命周期；WMDS 回退层次；循环终止条件 |
