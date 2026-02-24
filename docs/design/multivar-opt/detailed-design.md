# MTSHL-d 细化设计文档

> 阶段：细化（workflow.md §2.3）
> 前置文档：`architecture.md`（模块划分、接口规约、设计决策）
> 目标文件：`clpoly/polynomial_factorize_wang.hh`
> 实施里程碑：M1 → M5（顺序推进，每里程碑验收后进入下一个）

---

## 0. 类型别名与公共约定

所有新函数放在 `polynomial_factorize_wang.hh`，与现有函数同文件、同 `namespace clpoly`。

```cpp
// 文件内局部别名（每个模板函数内部 using 定义）
template<class var_order>
using PolyZp = polynomial_<Zp, lex_<var_order>>;

using UPZp = upolynomial_<Zp>;

template<class var_order>
using Mono = basic_monomial<lex_<var_order>>;

// forms[i] 类型：第 i 个因子的期望支撑骨架
template<class var_order>
using Forms = std::vector<std::vector<Mono<var_order>>>;
```

**命名约定**：
- 新函数统一用 `__mtshl_` 前缀（主流程）或 `__si_` 前缀（SparseInt 内部辅助）
- `aux_vars`：当前 j 步中辅助变量列表（x2,...,x_{j-1}）
- `sparse_betas`：SparseInt 内部随机选取的几何基底，`std::vector<Zp>` 与 `aux_vars` 等长
- `alpha_j_zp`：ideal alpha 的 Zp 投影，`Zp(alpha_j % p, p)`

---

## M1：`__si_vandermonde_solve` + `__si_theta_array_eval`

里程碑验收：单元测试——给定几何点求值能正确恢复稀疏多项式系数。

---

### M1-A：`__si_vandermonde_solve`

**功能**：求解转置 Vandermonde 系统 `v[l] = Σ_t c_t · θ_t^l`（l=1,...,s），恢复系数 `{c_t}`。

```cpp
// 求解 s×s Vandermonde 系统
// 输入：values[l-1] = v[l]（l=1..s），thetas[t-1] = θ_t（t=1..s）
// 输出：coeffs[t-1] = c_t；返回 false 表示矩阵奇异（θ_t 存在碰撞）
bool __si_vandermonde_solve(
    const std::vector<Zp>& values,   // v[1..s]，长度 s
    const std::vector<Zp>& thetas,   // θ_1,...,θ_s，长度 s
    std::vector<Zp>& coeffs);        // 输出 c_1,...,c_s，长度 s
```

**算法**：

矩阵变换：令 `V[l,t] = θ_t^l`（1-indexed）。将系数转换 `d_t = c_t · θ_t`，则：
```
v[l] = Σ_t d_t · θ_t^{l-1}
```
成为标准 Vandermonde 问题（矩阵 `W[l,t] = θ_t^{l-1}`，从 l=1 开始）。

标准 Vandermonde 解法（O(s²) Lagrange 插值）：
```
1. 检查 θ_t 两两不同；若有碰撞 → return false
2. 构造差积 ω_t = ∏_{j≠t} (θ_t - θ_j)；若 ω_t = 0 → return false
3. d_t = Σ_l v[l] · (Σ_{l'≠l} ∏_{l''≠l,l'} (·)) / ω_t  —— 实际用 Gaussian 消元更简
```

实际实现用 **带主元 Gaussian 消元**（O(s²)，s 小时更稳定）：
```
1. 构造增广矩阵 A[l,t] = θ_t^{l-1}（l=0..s-1，t=0..s-1），右端 b[l] = v[l+1]
   （注意：这里用标准 Vandermonde W，即从 θ^0 开始）
2. 前向消元 + 回代，在 Zp 中（利用 Zp::inv()）
3. 得到 d_t
4. c_t = d_t / θ_t（若 θ_t = 0 则矩阵奇异，但 Requires 保证 βk ≠ 0 → θ_t ≠ 0）
```

**错误处理**：
- θ_t 有碰撞 → return false（调用方会重新选 sparse_betas）
- 消元时主元为 0 → return false

---

### M1-B：`__si_theta_array_eval`

**功能**：利用 θ-array 技术，将 `f ∈ Zp[x1,...,x_{j-1}]` 在 s 个几何点处批量求值，
得到 s 个 `Zp[x1]` 像（每个点 (x2=β2^l,...,x_{j-1}=β_{j-1}^l)，l=1,...,s）。

```cpp
// θ-array 批量求值
// f ∈ Zp[x1,...,x_{j-1}]（x1 为主变量，aux_vars = [x2,...,x_{j-1}]）
// sparse_betas[k] = βk+2（即 β2,...,β_{j-1}，长度 j-2）
// s = 求值点数；输出 images[l-1] = f(x1, β2^l,...,β_{j-1}^l)，l=1..s
template<class var_order>
void __si_theta_array_eval(
    const PolyZp<var_order>& f,
    const std::vector<variable>& aux_vars,   // [x2,...,x_{j-1}]，长度 j-2
    const std::vector<Zp>& sparse_betas,     // β2,...,β_{j-1}，长度 j-2
    int s,
    std::vector<UPZp>& images);              // 输出，长度 s
```

**算法**（CASC 2016 §3.3 θ-array）：

预处理（每个非 x1 单项的 θ 值，只需算一次）：
```
for 每个 term m = c·x1^{e1}·x2^{e2}·...·x_{j-1}^{e_{j-1}} in f:
    theta_m = ∏_{k=0}^{j-3} sparse_betas[k] ^ (aux_vars[k]对应的指数)
              即 β2^{e2} · β3^{e3} · ... · β_{j-1}^{e_{j-1}}
    running_m = 1_zp  （θ_m^0，l=0 的初始值）
```

批量求值（逐 l 迭代，利用运行乘积）：
```
images.resize(s)
for l = 1 to s:
    images[l-1] ← 空 UPZp
    for 每个 term m（e1 为 x1 次数，c 为系数）:
        running_m *= theta_m           （变为 θ_m^l）
        images[l-1].add_coeff(e1, c * running_m)
    images[l-1].normalization()         （合并相同 x1 次数的项）
```

**代价**：O(s · |Supp(f)|)（每步一次乘法/项），相比朴素求值（需要指数运算）节省因子 (j-2)。

**注意**：`add_coeff(e1, val)` 意为在 images[l-1] 中将 `x1^{e1}` 项的系数增加 `val`。
实现时先按 e1 分组，最后 normalization。

---

## M2：`__mtshl_multi_bdp`

里程碑验收：单元测试——Σσi·bi = c 验证（r=2,3,4，多组随机输入）。

**功能**：双变量多因子 Diophantine 求解器，在 `Zp[x1,x2]` 中求解
`Σ σi·bi = c`，`bi = ∏_{l≠i}F[l]`，`deg(σi,x1) < deg(F[i],x1)`。

```cpp
// 多因子双变量 MDP（j=3 时的稠密回退）
// F[i], c ∈ Zp[x1,x2]；alpha2 = ideal α2（ZZ → Zp 转换在内部做）
// 返回 false 极罕见（某个求值点处因子共根）
template<class var_order>
bool __mtshl_multi_bdp(
    const std::vector<PolyZp<var_order>>& F,   // 当前因子 F[i] ∈ Zp[x1,x2]
    const PolyZp<var_order>& c,                 // RHS ∈ Zp[x1,x2]
    const variable& x1,
    const variable& x2,
    const Zp& alpha2,                           // ideal α2 mod p
    std::vector<PolyZp<var_order>>& result);    // 输出 σ[i] ∈ Zp[x1,x2]
```

**算法**（Taylor 展开 + 逐阶单变量 MDP）：

```
// 基础情形：评估 x2 = alpha2，在 Zp[x1] 中解多因子 Diophantine
F0[i] = F[i]|_{x2=alpha2}   ∈ Zp[x1]    （调用 assign(F[i], x2, alpha2)）
c0    = c|_{x2=alpha2}       ∈ Zp[x1]

// 解基础 r-因子单变量 Diophantine：Σ σ0[i]·b0[i] = c0
// b0[i] = ∏_{l≠i} F0[l]
success ← __mtshl_zp_univar_mdp(F0, c0, x1, sigma0)
if !success: return false

// 初始化结果：σ[i] 从常数项（在 x2 方向）开始
for i: result[i] = sigma0[i]   （作为 Zp[x1,x2] 常数 w.r.t. x2 的多项式）

// Taylor 提升：逐阶修正 (x2-alpha2)^k
error = c - Σ result[i] · bi（bi = ∏_{l≠i}F[l] 完整双变量版本）
for k = 1, 2, ..., deg(c, x2):
    ck = [(x2-alpha2)^k] error    ∈ Zp[x1]    （调用 __taylor_coeff_zp(error, x2, alpha2, k)）
    if ck == 0: continue

    // 修正项：Σ δk[i]·b0[i] = ck - (来自前序修正与 F[i] 高阶项的交叉项)
    // 简化：delta_rhs = ck（精确到 O((x2-alpha2)^{k+1}) 误差无影响）
    success ← __mtshl_zp_univar_mdp(F0, ck, x1, delta_k)
    if !success: return false

    for i: result[i] += delta_k[i] · (x2-alpha2)^k
    error = c - Σ result[i] · bi    （全量重算）

return true
```

**依赖**：`__mtshl_zp_univar_mdp`（内部辅助，见下）。

---

### 辅助：`__mtshl_zp_univar_mdp`

**功能**：在 `Zp[x1]` 中解 r 因子多项式 Diophantine：`Σ σi·∏_{l≠i}F[l] = c`。

```cpp
// Zp[x1] 中 r 因子 Diophantine（r-1 次 EEA 逐对归约）
// F[i], c ∈ Zp[x1]（单变量，UPZp）；输出 sigma[i]
// 失败（极罕见：F[i] 在 Zp 中不互素）返回 false
bool __mtshl_zp_univar_mdp(
    const std::vector<UPZp>& F,
    const UPZp& c,
    std::vector<UPZp>& sigma);
```

**算法**（GCL §6.3 推广到 r 因子）：

r=2 时：直接 EEA：`σ1·F1 + σ2·F0 = c` → 解 `σ0·F1 + σ1·F0 = c`，用 `polynomial_GCD` 扩展版。

r>2 时：逐对归约（O(r·d²)）：
```
// 构造 Bézout 链：s[i]·∏_{l≠i}F[l] 满足 Σ s[i]·∏_{l≠i}F[l] = 1 (mod all pairs coprime)
// 实际：sigma[i] = s[i]·c mod F[i]（逐步归约）
// 步骤（从 CASC 2016 Algorithm 2 单变量基础情形）：
1. g0 = F[0]
2. for i = 1..r-1:
     EEA: alpha·g_{i-1} + beta·F[i] = gcd  →  gcd ∈ Zp（常数，因为 F[i] 两两互素）
     s[0..i-1] *= beta / gcd
     s[i] = alpha / gcd
     g_i = g_{i-1} · F[i]
3. for i in 0..r-1:
     sigma[i] = (s[i] · c) rem F[i]    （保证 deg < deg(F[i])）
```

注：`polynomial_GCD` 已有单变量 ZZ 版本（`polynomial_factorize_univar.hh`）；
需用 Zp 版本或调用现有 `gcd(UPZp, UPZp)` （待确认是否已有 Zp EEA）。

**M2 新增辅助**：`__si_zp_eea`（若现有 EEA 不支持 UPZp，新增约 20 行）：
```cpp
// Zp[x1] 扩展 Euclidean：返回 gcd，输出 s·a + t·b = gcd
Zp __si_zp_eea(const UPZp& a, const UPZp& b, UPZp& s, UPZp& t);
```

---

## M3：`__mtshl_sparse_int`

里程碑验收：单元测试——与 `__mtshl_multi_bdp` 结果交叉验证（j=3，稀疏输入时结果一致）。

**功能**：稀疏插值 MDP 求解器，在 `Zp[x1,...,x_{j-1}]` 中求解
`Σ σi·bi = c`，用 θ-array + Vandermonde 恢复。

```cpp
// 稀疏插值 MDP 求解器（j≥3）
// F[i] ∈ Zp[x1,...,x_{j-1}]；forms[i] = 期望 Supp(σi)
// aux_vars = [x2,...,x_{j-1}]，sparse_betas 由本函数内部生成
// 返回 false = 支撑假设不成立或 Vandermonde 奇异
template<class var_order>
bool __mtshl_sparse_int(
    const std::vector<PolyZp<var_order>>& F,
    const PolyZp<var_order>& c,
    const std::vector<std::vector<Mono<var_order>>>& forms,   // forms[i]
    const variable& x1,
    const std::vector<variable>& aux_vars,   // [x2,...,x_{j-1}]
    uint32_t p,
    std::vector<PolyZp<var_order>>& result); // 输出 σ[i]
```

**算法**（CASC 2018 Algorithm 3 + CASC 2016 §3.3 θ-array，r 因子版）：

```
s = max_i(|forms[i]|)        （求值点数 = 最大骨架大小）
if s == 0: 所有 forms 为空 → σ[i] = 0 → 直接返回 true

// Step 0：内部随机选取 sparse_betas（每次调用独立）
sparse_betas = [random_nonzero_zp(p) for _ in aux_vars]

// Step 1：θ-array 批量求值
// 对 c 和所有 F[i] 在 s 个几何点处求值
images_c[l] = __si_theta_array_eval(c, aux_vars, sparse_betas, s)      // s 个 UPZp
images_F[i][l] = __si_theta_array_eval(F[i], aux_vars, sparse_betas, s) // r×s 个 UPZp

// Step 2：逐点单变量 MDP
// 在每个求值点 l (l=1..s) 处解 Zp[x1] r-因子 Diophantine
sigma_vals[i][l] ∈ UPZp，l=1..s
for l = 1 to s:
    F_at_l[i] = images_F[i][l-1]   // F[i](x1, β2^l,...) ∈ Zp[x1]
    c_at_l    = images_c[l-1]
    success ← __mtshl_zp_univar_mdp(F_at_l, c_at_l, sigma_vals_at_l)
    if !success: return false
    sigma_vals[i][l] = sigma_vals_at_l[i]   for each i

// Step 3：Vandermonde 恢复
// 对每个因子 i，逐 x1 次数 d 恢复 σi 中所有含 x1^d 的项的系数
result[i] ← 空 PolyZp

for i = 0..r-1:
    按 x1 次数分组 forms[i] → groups[d] = {terms in forms[i] with x1-degree d}
    for 每个 d 出现在 groups[i]:
        terms_d = groups[d]        // 该 x1 次数下的所有期望项
        t = |terms_d|
        // 预计算每个 term 的 θ_t = ∏_k βk^{e_{t,k}}（aux_vars 部分）
        thetas = [compute_theta(term, sparse_betas) for term in terms_d]
        // 从 s 个点值中取 t 个（用前 t 个）
        // v[l] = 点 l 处 sigma_vals[i][l] 中 x1^d 的系数
        values = [coeff(sigma_vals[i][l], d) for l in 1..t]
        // 解 Vandermonde
        coeffs = []
        ok ← __si_vandermonde_solve(values, thetas[0..t-1], coeffs)
        if !ok: return false
        // 将恢复的系数写入 result[i]
        for k, term in enumerate(terms_d):
            if coeffs[k] != 0:
                result[i].push_back({term, coeffs[k]})
    result[i].normalization()

// Step 4：验证（CASC 2018 Algorithm 3 step 5）
// 计算 Σ result[i]·bi，检查是否等于 c
product = compute_sum_sigma_bi(result, F)   // Σ result[i]·∏_{l≠i}F[l]
if product != c: return false               // 骨架假设失败

return true
```

**关键实现细节**：

1. `compute_theta(term, sparse_betas)`：对于 `basic_monomial` 中 aux_vars 对应的指数，计算 `∏_k βk^{e_k}`（O(j) 次 Zp 乘法）。
2. `coeff(uply, d)`：在 `UPZp` 中查找 x1^d 的系数（O(|f|) 遍历）。
3 验证步骤：计算 `Σ result[i]·∏_{l≠i}F[l]` 然后和 c 比较（PolyZp 相等判断）。

---

## M4：`__mtshl_step_j` + `__mtshl_lift` + `__mtshl_wmds`

里程碑验收：`test_multivar_hensel` 全部通过。

---

### M4-A：`__mtshl_taylor_coeff_zp`（辅助）

`__taylor_coeff` 仅支持 ZZ。新增 Zp 版本（算法完全相同，系数域不同）：

```cpp
// Zp 多项式的第 j 阶 Taylor 系数（复用现有 ZZ 版逻辑，系数域改为 Zp）
template<class var_order>
PolyZp<var_order> __mtshl_taylor_coeff_zp(
    const PolyZp<var_order>& f,
    const variable& xk,
    const Zp& alpha_k,
    int j);
```

**实现**：与 `__taylor_coeff`（L24–52）逐字复制，将 `ZZ` → `Zp`，`-alpha_k` → `Zp(0,p) - alpha_k`。

---

### M4-B：`__mtshl_wmds`（j>3 时的稠密回退 MDP）

**功能**：Wang 递归 MDP 求解器的 Zp 版本，处理 j>3 时 SparseInt 失败的回退。
将 `Zp[x1,...,x_{j-1}]` 中 r-因子 MDP 递归归约到 `Zp[x1,...,x_{j-2}]`。

```cpp
// Wang 递归多变量 MDP（Zp 版，j>3 时回退）
// F[i] ∈ Zp[x1,...,x_{j-1}]；aux_vars = [x2,...,x_{j-1}]（长度 j-2）
// ideal_alphas_zp: x2,...,x_{j-1} 对应的 ideal alpha (Zp 值)，用于递归评估
template<class var_order>
bool __mtshl_wmds(
    const std::vector<PolyZp<var_order>>& F,
    const PolyZp<var_order>& c,
    const variable& x1,
    const std::vector<variable>& aux_vars,       // [x2,...,x_{j-1}]
    const std::vector<Zp>& ideal_alphas_zp,      // α2,...,α_{j-1} mod p
    std::vector<PolyZp<var_order>>& result);
```

**算法**（从现有 `__multivar_diophantine` 演化，改为 Zp）：

```
if aux_vars 只含 x2（即 j=3）:
    // 单变量基础情形：调用 __mtshl_multi_bdp（j=3 专用，见 M2）
    // 但此分支不应被调用（j=3 时直接走 multi_bdp，j>3 才走 wmds）
    internal error

// 递归情形：归约到一维更少的 MDP
xj_prev = aux_vars.back()          // x_{j-1}
alpha_prev = ideal_alphas_zp.back()
aux_vars_sub = aux_vars[0..end-2]  // x2,...,x_{j-2}
alphas_sub   = ideal_alphas_zp[0..end-2]

// 基础求值：把 x_{j-1} 代入 ideal alpha，得到更低维的 F 和 c
F_base[i] = assign(F[i], xj_prev, alpha_prev)   ∈ Zp[x1,...,x_{j-2}]
c_base    = assign(c,    xj_prev, alpha_prev)    ∈ Zp[x1,...,x_{j-2}]

// 递归求解低维 MDP
success ← __mtshl_wmds(F_base, c_base, x1, aux_vars_sub, alphas_sub, result_base)
if !success: return false

// 初始化 result[i] = result_base[i]（在 Zp[x1,...,x_{j-1}] 中 x_{j-1} 独立）
for i: result[i] = result_base[i]

// Taylor 提升：逐阶修正 (x_{j-1} - alpha_prev)^k
error = c - Σ result[i] · bi   (bi = ∏_{l≠i} F[l]，在 Zp[x1,...,x_{j-1}] 中)
for k = 1 to deg(c, xj_prev):
    ck = __mtshl_taylor_coeff_zp(error, xj_prev, alpha_prev, k)  ∈ Zp[x1,...,x_{j-2}]
    if ck == 0: continue
    // 递归修正
    success ← __mtshl_wmds(F_base, ck, x1, aux_vars_sub, alphas_sub, delta_k)
    if !success: return false
    for i: result[i] += delta_k[i] · (xj_prev - alpha_prev)^k
    error = c - Σ result[i] · bi

return true
```

**注**：基础情形（`|aux_vars| == 1`，即 j=3 时的基础求值层）直接调用 `__mtshl_zp_univar_mdp`（单变量）。

---

### M4-C：`__mtshl_step_j`

**功能**：MTSHL-d 第 j 步。将 `Zp[x1,...,x_{j-1}]` 中的因子 in-place 提升到 `Zp[x1,...,xj]`。

```cpp
// MTSHL-d 第 j 步（j=2..n）
// F[i]：in-place 从 Zp[x1,...,x_{j-1}] 更新为 Zp[x1,...,xj]
// aux_vars = [x2,...,x_{j-1}]（当 j=2 时为空）
// ideal_alphas_zp[k] = αk+2 mod p（与 aux_vars 等长）
// 返回 false = 求值点不适用，__mtshl_lift 须换求值点重试
template<class var_order>
bool __mtshl_step_j(
    const PolyZp<var_order>& aj,               // f_scaled 代入 {x_{j+1},...} 后的 Zp 像
    std::vector<PolyZp<var_order>>& F,         // in-place 更新
    const variable& xj,
    const Zp& alpha_j,                         // ideal αj mod p
    const variable& x1,
    const std::vector<variable>& aux_vars,     // [x2,...,x_{j-1}]（j=2 时为空）
    const std::vector<Zp>& ideal_alphas_zp,    // α2,...,α_{j-1} mod p（j=2 时为空）
    uint32_t p);
```

**算法**（§1.3 伪代码的 C++ 直译）：

```
r = F.size()

// 入口：初始化 forms[i] = Supp(F[i])
forms[i] = {m for (m,c) in F[i]}   for each i

// 计算初始误差
error = aj - product(F)   // PolyZp 减法

for k = 1; error 不为零; ++k:

    // 提取第 k 个 Taylor 系数
    ck = __mtshl_taylor_coeff_zp(error, xj, alpha_j, k)  ∈ Zp[x1,...,x_{j-1}]
    if ck.empty(): continue   // ck == 0

    sigma_k = []    // 将被填充

    // ————————————————————————————
    //  j = 2 路径：直接单变量 MDP
    // ————————————————————————————
    if j == 2:
        // F[i] ∈ Zp[x1]，ck ∈ Zp[x1]
        success ← __mtshl_zp_univar_mdp(F_as_upoly, ck_as_upoly, sigma_k_upoly)
        // F_as_upoly[i]：从 PolyZp 转为 UPZp（x1 是主变量，F[i] 不含其他变量）
        if !success: return false
        // 将 sigma_k_upoly 转回 PolyZp
        ...

    // ————————————————————————————
    //  j ≥ 3 路径：SparseInt + 回退
    // ————————————————————————————
    else:
        // 主路径：SparseInt
        success ← __mtshl_sparse_int(F, ck, forms, x1, aux_vars, p, sigma_k)

        // 重试一次（换 sparse_betas，失败多因随机碰撞）
        if !success:
            success ← __mtshl_sparse_int(F, ck, forms, x1, aux_vars, p, sigma_k)

        // 回退
        if !success:
            if j == 3:
                success ← __mtshl_multi_bdp(F, ck, x1, aux_vars[0], ideal_alphas_zp[0], sigma_k)
            else:   // j > 3
                success ← __mtshl_wmds(F, ck, x1, aux_vars, ideal_alphas_zp, sigma_k)

        if !success: return false

    // 更新因子和 forms
    for i:
        F[i] += sigma_k[i] * (xj - alpha_j)^k   // PolyZp 多项式更新
        forms[i] = Supp(sigma_k[i])              // 更新骨架

    // 全量重算误差
    error = aj - product(F)

return true
```

**辅助操作**：
- `product(F)`：计算 `∏F[i]`（PolyZp 乘法，r-1 次）
- `make_monomial(xj, k, alpha_j, comp_ptr)`：构造 `(xj - alpha_j)^k` 的 PolyZp 表示

**j=2 路径的类型转换**：
- `PolyZp → UPZp`：F[i] 只含 x1，遍历每项 `(m, c)` 取 `x1` 的指数，构造 `UPZp`
- `UPZp → PolyZp`：反向，每项 `(x1^d, c)` 转为 `({x1→d}, c)` 的 PolyZp

---

### M4-D：`__mtshl_lift`

**功能**：MTSHL-d 顶层驱动，将单变量因子逐变量提升到 n 变量 Z 系数因子。

```cpp
// MTSHL-d 顶层提升
// 返回空 vector 表示提升失败（__wang_core 须重新选求值点）
template<class var_order>
std::vector<polynomial_<ZZ, lex_<var_order>>>
__mtshl_lift(
    const polynomial_<ZZ, lex_<var_order>>& f_scaled,
    const std::vector<upolynomial_<ZZ>>& scaled_factors,
    const std::map<variable, ZZ>& eval_point,
    const variable& main_var,
    uint32_t p);
```

**算法**：

```
using Poly   = polynomial_<ZZ, lex_<var_order>>;
using PolyZp = polynomial_<Zp, lex_<var_order>>;

r = scaled_factors.size()
comp_ptr = f_scaled.comp_ptr()

// 构造 lex_ 变量顺序（从 eval_point 中获取 x2,...,xn）
vars_in_order = [main_var] + [v for v in eval_point keys，按 lex 顺序]
aux_var_list = vars_in_order[1..]   // x2,...,xn

// 阶段 A：初始化，将单变量因子转为 Zp[x1]
F[i] = polynomial_mod_zp(scaled_factors[i], p, comp_ptr)   // UPZp → PolyZp（仅含 x1）

// ideal_alphas_zp：将 eval_point 中的 ZZ 转为 Zp，与 aux_var_list 对应
ideal_alphas_zp[k] = Zp(eval_point[aux_var_list[k]] % p, p)

// 阶段 B：逐变量提升 j = 2,...,n
for j = 2 to n:
    xj = aux_var_list[j-2]                           // j=2 → x2，j=3 → x3，...
    alpha_j = ideal_alphas_zp[j-2]                   // αj mod p

    // aj = f_scaled 在 {x_{j+1}=α_{j+1},...,xn=αn} 处的 Zp 像
    aj = assign_partial_zp(f_scaled, aux_var_list[j-1..], ideal_alphas_zp[j-1..], p)

    aux_sub = aux_var_list[0..j-3]                  // x2,...,x_{j-1}
    alphas_sub = ideal_alphas_zp[0..j-3]

    success ← __mtshl_step_j(aj, F, xj, alpha_j, main_var, aux_sub, alphas_sub, p)
    if !success: return {}                           // 提升失败

// 阶段 C：系数恢复（对称约化）
result[i] = symmetric_mod_poly(F[i], p)    // PolyZp → Poly（每个系数用 __symmetric_mod）

return result
```

**辅助**：
- `polynomial_mod_zp(upoly, p, comp_ptr)` — 将 `upolynomial_<ZZ>` 转为仅含 `x1` 的 `PolyZp`
- `assign_partial_zp(f, vars, alphas, p)` — 将 `f ∈ Z[x1,...,xn]` 代入部分变量，返回 `PolyZp`
- `symmetric_mod_poly(poly_zp, p)` — 逐系数调用 `__symmetric_mod`，返回 `Poly`

---

## M5：接入 `__wang_core`

里程碑验收：test_factorize_multivar 全部通过；bench 无退化。

**改动位置**：`__wang_core`（L1062–1248），约 20 行：

**改动 1**：在调用 `__wang_leading_coeff` 之前，计算 Mignotte 界并选取素数 p：
```cpp
// 新增：计算 Mignotte_bound(g) 并选 p > 2 * bound
uint32_t mtshl_p = __select_mtshl_prime(g);
if (mtshl_p == 0) continue;   // 找不到合适素数（罕见，换求值点）
```

**改动 2**：替换 `__multivar_hensel_lift` 调用（L1136–1138）：
```cpp
// 旧代码（删除）：
// auto mv_factors = __multivar_hensel_lift(
//     lc_result.f_scaled, lc_result.scaled_factors,
//     lc_result.lc_targets, eval, x1);

// 新代码：
auto mv_factors = __mtshl_lift(
    lc_result.f_scaled, lc_result.scaled_factors,
    eval, x1, mtshl_p);
```

**改动 3**：新增辅助 `__select_mtshl_prime`（约 20 行）：
```cpp
// 计算 Mignotte_bound(f) 并返回满足 p > 2*bound 的最小机器素数
// 若 uint32_t 范围内无合适素数则返回 0
template<class var_order>
uint32_t __select_mtshl_prime(
    const polynomial_<ZZ, lex_<var_order>>& f);
```

**Mignotte 界计算**（标准公式）：
```
B(f) = 2^n · ||f||_2 · ∏ max(1, ||fi||_2)  （n = deg(f)）
     ≈ 2^deg · ||f||_inf · r   （简化估计，保守版）
```

实际实现参考：`Mignotte_bound = (n+1)^{1/2} · 2^n · ||f||_inf`（GCL §6.1）。
选 p：从预定义素数表（31-bit 素数，如 `p = 2147483647` 开始往下找）中选满足条件的第一个。

---

## 错误处理策略汇总

| 函数 | 返回值 | 失败含义 | 调用方处理 |
|------|--------|----------|-----------|
| `__si_vandermonde_solve` | `bool` | θ_t 碰撞/奇异 | `__mtshl_sparse_int` 换 sparse_betas 重试 |
| `__si_theta_array_eval` | void | 不可失败 | — |
| `__mtshl_zp_univar_mdp` | `bool` | F[i] 不互素 (mod p) | 上层返回 false |
| `__mtshl_multi_bdp` | `bool` | 极罕见共根 | `__mtshl_step_j` 返回 false |
| `__mtshl_sparse_int` | `bool` | 骨架假设失败/奇异 | `__mtshl_step_j` 重试→回退 |
| `__mtshl_wmds` | `bool` | 递归 MDP 失败 | `__mtshl_step_j` 返回 false |
| `__mtshl_step_j` | `bool` | 任意 MDP 失败 | `__mtshl_lift` 返回 {} |
| `__mtshl_lift` | `vector<Poly>` | 空=失败 | `__wang_core` 换求值点继续 |

---

## 复用点汇总

| 复用对象 | 位置 | 在哪里复用 |
|---------|------|-----------|
| `__taylor_coeff` | wang.hh L24 | — (新增 `__mtshl_taylor_coeff_zp` 为 Zp 副本) |
| `pair_vec_div` | polynomial.hh | `__mtshl_taylor_coeff_zp` 内部（精确除以 (xk-αk)） |
| `assign(poly, var, val)` | polynomial.hh | `__mtshl_multi_bdp`（基础求值）、`__mtshl_wmds` |
| `polynomial_<>::normalization()` | basic_polynomial.hh | 所有 PolyZp 更新后 |
| `__symmetric_mod` | factorize_univar.hh L95 | `__mtshl_lift` 阶段 C |
| `polynomial_mod(upoly, p)` | upolynomial.hh L104 | 改写为 PolyZp 版本 |
| `Zp::inv()` | number.hh | Vandermonde 消元、EEA |

---

## 测试规划

**M1 测试**（`test_mtshl_m1`）：
- `__si_vandermonde_solve`：给定 θ = [2,3,5]（Zp=97），v[l] = Σ c_t·θ_t^l，恢复 c
- `__si_theta_array_eval`：已知 f = 3x1²x2 + 2x1x3 + x2²（j=4），β2=2,β3=3，s=3，验证 images[l]

**M2 测试**（`test_mtshl_m2`）：
- `__mtshl_zp_univar_mdp`：r=2,3,4，随机因子，验证 Σσi·bi = c
- `__mtshl_multi_bdp`：j=3，稀疏和稠密输入，验证 Σσi·bi = c

**M3 测试**（`test_mtshl_m3`）：
- `__mtshl_sparse_int`：与 `__mtshl_multi_bdp` 交叉验证（j=3，相同输入，结果应一致）

**M4 测试**（`test_multivar_hensel`，替换现有单元测试）：
- `__mtshl_lift`：给定 f_scaled 和 scaled_factors，验证 ∏G[i] = f_scaled

**M5 测试**：
- 全量 `test_factorize_multivar`：覆盖 bivar/trivar/多变量情形
- `make bench-all`：无性能退化
