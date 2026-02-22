# P1a van Hoeij LLL 重组：细化设计

> 阶段：细化（§2.3）
> 依赖：`docs/design/vanhoeij/architecture.md`（架构文档，不得变更模块划分）
> 实现目标文件：`clpoly/polynomial_factorize_univar.hh`

---

## 0. 文件级变更摘要

| 变更 | 说明 |
|------|------|
| 新增函数 `__cld_polys` | M1，CLD 多项式计算 |
| 新增函数 `__build_cld_matrix` | M2，格矩阵列喂入 |
| 新增函数 `__lll_reduce` | M3，整数 LLL + 幺模变换 |
| 新增函数 `__extract_candidates` | M4，候选因子子集提取 |
| 新增函数 `__vanhoeij_recombine` | M5，主控循环 |
| 重命名现有函数 | `__factor_recombine` → `__zassenhaus_recombine` |
| 新增路由函数 `__factor_recombine` | 替换原签名，内部路由 Zassenhaus / van Hoeij |

所有新代码位于 `polynomial_factorize_univar.hh`，在现有 `§7.5 Zassenhaus 因子重组` 之后插入。无新文件，无新头文件依赖（QQ 通过现有 `#include <clpoly/polynomial.hh>` 已可用）。

---

## 1. 辅助常量与类型

```cpp
// 编译期常量，可通过 benchmark 调整
static constexpr int ZASSENHAUS_THRESHOLD = 10;

// 格矩阵：行主序，M[i] 是第 i 个基向量（整数，行向量）
using LLLMatrix = std::vector<std::vector<ZZ>>;
```

---

## 2. M1 `__cld_polys`

### 函数签名

```cpp
// 计算 CLD 多项式：C_i = (f_star / h_i) · h_i'  mod m，对称约化
// 参数：
//   f_star        — 当前待分解多项式（ZZ 系数，本原，lc > 0）
//   active_factors — Hensel 提升因子（mod m，对称约化，首一）
//   m             — Hensel 模数
// 返回：r = active_factors.size() 个 CLD 多项式（mod m，对称约化）
std::vector<upolynomial_<ZZ>>
__cld_polys(
    const upolynomial_<ZZ>&              f_star,
    const std::vector<upolynomial_<ZZ>>& active_factors,
    const ZZ&                            m);
```

### 实现步骤

```cpp
for each h_i in active_factors:
    // 步骤 1：精确除法 q_i = f_star / h_i  in Z_m[x]
    //   复用：__upoly_divmod_mod(q, r, f_star, h_i, m)
    //   断言：r 为空（h_i | f_star mod m，Hensel 提升保证）
    upolynomial_<ZZ> q_i, r_i;
    __upoly_divmod_mod(q_i, r_i, f_star, h_i, m);
    assert(r_i.empty());  // DEBUG 断言

    // 步骤 2：形式导数 h_i' = derivative(h_i)，再对系数 mod m
    //   复用：derivative(h_i)（upolynomial.hh:176）
    //   复用：__upoly_mod_coeff(h_prime, m)
    auto h_prime = derivative(h_i);
    __upoly_mod_coeff(h_prime, m);

    // 步骤 3：C_i = q_i · h_i'  mod m，然后对称约化
    //   复用：__upoly_mul_mod(q_i, h_prime, m)
    //   复用：__upoly_symmetric_mod(C_i, m)
    auto C_i = __upoly_mul_mod(q_i, h_prime, m);
    C_i = __upoly_symmetric_mod(C_i, m);
    result.push_back(std::move(C_i));
```

### 复用点

| 调用 | 来源 |
|------|------|
| `__upoly_divmod_mod(q, r, f, g, m)` | 同文件 §4.1.4（ZZ 系数 Z_m 除法） |
| `derivative(f)` | `upolynomial.hh:176` |
| `__upoly_mod_coeff(f, m)` | 同文件 §4.1.4 |
| `__upoly_mul_mod(a, b, m)` | 同文件 §6（乘并 mod） |
| `__upoly_symmetric_mod(f, m)` | 同文件 §4.2.2 |

### 错误处理

- `assert(r_i.empty())`：若 Hensel 提升有误，r_i 非空，DEBUG 下触发断言（Release 下不检查，上层调用方保证 precondition）。

---

## 3. M2 `__build_cld_matrix`

### 函数签名

```cpp
// 按内螺旋顺序将 CLD 系数列喂入格矩阵 M，扩展 M 直到加入 J_target 新列。
// 参数：
//   M        — 格矩阵，原地扩展（当前 (r+J_cur)×(r+J_cur)；扩展后 (r+J_cur+J_new)×(r+J_cur+J_new)）
//   cld      — r 个 CLD 多项式（mod m，对称约化），r 从 cld.size() 推导
//   J_cur    — 已喂入列数（用于定位螺旋起始位置；非零时矩阵已包含 r+J_cur 行/列）
//   J_target — 本次目标新增列数（> 0）
//   m        — Hensel 模数（用于列过滤）
// 返回：J_new（实际新增列数，≤ J_target；若全部被过滤则为 0）
// 注：过滤被跳过的列不计入 J_new，但仍消耗螺旋位置（k 仍递增）。
//     J_cur 仅计已实际加入的列数，不计被过滤位置。下次调用以 J_cur+J_new 开始螺旋。
int __build_cld_matrix(
    LLLMatrix&                           M,
    const std::vector<upolynomial_<ZZ>>& cld,
    int                                  J_cur,
    int                                  J_target,
    const ZZ&                            m);
```

### 矩阵存储结构

M 是行主序方阵，每行是一个格基向量（整数坐标）。

初始状态（由 `__vanhoeij_recombine` 构造，不由 M2 负责）：

```
M = [ 2^U_exp * I_r ]   （r×r，仅 M[i][i] = 2^U_exp，其余为 0）
```

每次加入第 j 列（j = J_cur, J_cur+1, ...，当前 j = J_cur + offset）：
1. **扩展现有行**：对每个现有行 M[i]，`M[i].push_back(ZZ(0))`
2. **新增数据行**：创建长度 `r + j + 1` 的新行：
   - 前 r 项：`upoly_coeff(cld[k], col_idx)`（k = 0..r-1）
   - 接下来 j 项：均为 ZZ(0)
   - 末项：ZZ(1)（对角线 identity 块）

**注**：`upolynomial_<ZZ>` 不支持 `operator[]` 按次数索引，需如下辅助函数（M2 内部 lambda）：
```cpp
auto upoly_coeff = [](const upolynomial_<ZZ>& p, int deg) -> ZZ {
    for (auto& term : p)
        if ((int)term.first.deg() == deg) return term.second;
    return ZZ(0);
};
```

### 内螺旋列选取

CLD 多项式次数为 `deg(f) - 1`，系数索引 0..N-1（N = deg(f)）。
螺旋序列（0-indexed 位置 k → 系数索引）：

```
k = 0: idx = 0
k = 1: idx = N-1
k = 2: idx = 1
k = 3: idx = N-2
...
k 偶数: idx = k / 2
k 奇数: idx = N - 1 - (k-1) / 2
```

M2 从螺旋位置 `k = J_cur` 开始遍历，尝试加入至多 `J_target` 个有效列（跳过被过滤列），直到：
- 加入了 `J_target` 个列，或
- 螺旋位置耗尽（`k >= N`，必须显式 `break`）

```cpp
int r = (int)cld.size();  // 从 cld.size() 推导，不作为参数传入

// N = 螺旋总位置数 = deg(f_star)。
// C_i 次数 ≤ deg(f)-1，故从所有 cld[i] 最高次 + 1 推导：
int N = 0;
for (const auto& c : cld)
    if (!c.empty())
        N = std::max(N, (int)c.front().first.deg() + 1);
// 若所有 CLD 均为零（理论上不可能），N = 0，循环不执行，返回 0。

int J_new = 0;
for (int k = J_cur; J_new < J_target && k < N; ++k) {
    int col_idx = (k % 2 == 0) ? (k / 2) : (N - 1 - (k - 1) / 2);
    // 列过滤（见下节）...
    // 若通过：扩展矩阵，++J_new
}
return J_new;
```

### 列过滤条件（来自 FLINT `CLD_mat.c`）

对螺旋位置 k（系数索引 col_idx）：

```
bound = max_{i=0..r-1} |cld[i][col_idx]| · sqrt(N + 1)
阈值  = m / 2^(1.5 * r)

若 bound > 阈值：跳过该列（过滤）
```

> **实现注**：`sqrt(N+1)` 和 `m / 2^(1.5*r)` 均用浮点估算（`double`），仅用于过滤决策，不影响正确性。过滤被跳过的列不计入 J_new，但仍消耗螺旋位置（即 `k` 仍递增）。

> **初版简化**：可以先实现无过滤版本（所有列均接受），正确性不受影响，后续优化时加过滤条件。

### 返回值语义

`J_new` = 实际加入的列数（≤ J_target）。若 `J_new == 0`（全被过滤或螺旋耗尽），调用方（M5）倍增 J_target 后重试；若 J_target > J_max 则 fallback Zassenhaus（见 §6 M5 节）。

---

## 4. M3 `__lll_reduce`

### 函数签名

```cpp
// 对整数矩阵 M 做 LLL 基规约（Cohen §2.6 算法 2.6.3，δ = 3/4），
// 记录幺模变换矩阵 U（M_new = U · M_old），返回短向量行下标。
// 参数：
//   M  — 格基矩阵（原地规约，n×n）
//   U  — 输出：幺模变换矩阵（n×n，本函数内初始化为 I_n）
//   B  — 短向量阈值：返回所有 ‖M[i]‖² ≤ B 的行
// 返回：short_rows（满足范数条件的行下标，按范数²升序）
std::vector<int>
__lll_reduce(
    LLLMatrix& M,
    LLLMatrix& U,
    const ZZ&  B);
```

### 内部数据结构

```cpp
int n = M.size();
// Gram-Schmidt 系数和平方范数（QQ 精确算术）
std::vector<std::vector<QQ>> mu(n, std::vector<QQ>(n, QQ(0)));
std::vector<QQ>              B_gs(n, QQ(0));

// U 初始化为单位矩阵
U.assign(n, std::vector<ZZ>(n, ZZ(0)));
for (int i = 0; i < n; ++i) U[i][i] = ZZ(1);
```

### Cohen §2.6 算法步骤

```
// 辅助：QQ 四舍五入到最近 ZZ（round to nearest, ties round up）
// 等价于 floor(q + 1/2) = floor((2*num + den) / (2*den))
// ZZ::fdiv_q 签名：static void fdiv_q(ZZ& result, const ZZ& a, const ZZ& b)
auto round_qq = [](const QQ& q) -> ZZ {
    ZZ a = q.get_num() * ZZ(2) + q.get_den();  // 2*num + den（den > 0 保证）
    ZZ b = q.get_den() * ZZ(2);                  // 2*den > 0
    ZZ result;
    ZZ::fdiv_q(result, a, b);                    // floor division，向负无穷取整
    return result;
};

// 辅助：行 a 内积（ZZ 精确）
auto dot = [](const std::vector<ZZ>& a, const std::vector<ZZ>& b) -> ZZ { ... };

// 辅助：行操作（同步作用于 M 和 U）：M[i] -= c * M[j], U[i] -= c * U[j]
auto row_sub = [&](int i, int j, const ZZ& c) {
    for (int k = 0; k < n; ++k) {
        M[i][k] -= c * M[j][k];
        U[i][k] -= c * U[j][k];
    }
};
auto row_swap = [&](int i, int j) {
    std::swap(M[i], M[j]);
    std::swap(U[i], U[j]);
};

// 初始化 Gram-Schmidt
B_gs[0] = QQ(dot(M[0], M[0]), ZZ(1));
for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
        QQ num(dot(M[i], M[j]), ZZ(1));
        // 减去已投影部分
        for (int l = 0; l < j; ++l)
            num -= mu[i][l] * mu[j][l] * B_gs[l];
        mu[i][j] = (B_gs[j] == QQ(0)) ? QQ(0) : num / B_gs[j];
    }
    B_gs[i] = QQ(dot(M[i], M[i]), ZZ(1));
    for (int j = 0; j < i; ++j)
        B_gs[i] -= mu[i][j] * mu[i][j] * B_gs[j];
}

// LLL 主循环（Cohen 算法 2.6.3）
int k = 1;
while (k < n) {
    // 大小规约：M[k] -= round(mu[k][k-1]) * M[k-1]
    ZZ q = round_qq(mu[k][k-1]);
    if (q != ZZ(0)) {
        row_sub(k, k-1, q);
        // 更新 mu[k][0..k-2]
        for (int j = 0; j < k-1; ++j)
            mu[k][j] -= QQ(q, ZZ(1)) * mu[k-1][j];
        mu[k][k-1] -= QQ(q, ZZ(1));
    }

    // Lovász 条件：B_gs[k] >= (3/4 - mu[k][k-1]^2) * B_gs[k-1]
    QQ lhs = B_gs[k];
    QQ rhs = (QQ(3, 4) - mu[k][k-1] * mu[k][k-1]) * B_gs[k-1];

    if (lhs >= rhs) {
        // 额外大小规约（对 j = k-2 down to 0）
        for (int j = k-2; j >= 0; --j) {
            ZZ q2 = round_qq(mu[k][j]);
            if (q2 != ZZ(0)) {
                row_sub(k, j, q2);
                for (int l = 0; l < j; ++l)
                    mu[k][l] -= QQ(q2, ZZ(1)) * mu[j][l];
                mu[k][j] -= QQ(q2, ZZ(1));
            }
        }
        ++k;
    } else {
        // 交换 M[k] 和 M[k-1]，更新 Gram-Schmidt
        QQ mu_old = mu[k][k-1];
        QQ B_new  = B_gs[k] + mu_old * mu_old * B_gs[k-1];
        // mu_new：新的 μ_{k,k-1}（swap 前先计算，swap 后需要显式写入）
        QQ mu_new = (B_new != QQ(0)) ? mu_old * B_gs[k-1] / B_new : QQ(0);
        if (B_new != QQ(0)) {
            B_gs[k]   = B_gs[k] * B_gs[k-1] / B_new;
            B_gs[k-1] = B_new;
        }
        row_swap(k, k-1);
        std::swap(mu[k], mu[k-1]);
        // swap 后：mu[k] = 原 mu[k-1]（其 [k-1] 位置 = 0）
        // 必须显式写入新的 μ_{k,k-1}，否则该位置是 0
        mu[k][k-1] = mu_new;
        // 修正 mu[j][k] 和 mu[j][k-1] (j > k)
        for (int j = k+1; j < n; ++j) {
            QQ t = mu[j][k];
            mu[j][k]   = mu[j][k-1] - mu_old * t;
            mu[j][k-1] = t + mu_new * mu[j][k];  // 用 mu_new，不用 mu[k][k-1]
        }
        k = std::max(k-1, 1);
    }
}

// 收集短向量
std::vector<int> short_rows;
for (int i = 0; i < n; ++i) {
    ZZ norm_sq = dot(M[i], M[i]);
    if (norm_sq <= B)
        short_rows.push_back(i);
}
std::sort(short_rows.begin(), short_rows.end(),
    [&](int a, int b) { return dot(M[a],M[a]) < dot(M[b],M[b]); });
// 注：dot 是本函数内的局部 lambda，必须用 [&] 或 [&M, &dot] 捕获；[&M] 不够。
return short_rows;
```

> **实现注**：`dot` 和 `QQ(3,4)` 使用已有的 ZZ/QQ 算术。QQ 构造 `QQ(p, q)` 需 `q > 0`，CLPoly QQ 构造函数应自动规范化。`QQ::operator>=` 和 `-` 需可用（参考现有 QQ 实现）。

### 复用点

| 调用 | 来源 |
|------|------|
| `ZZ::fdiv_q` | `clpoly/number/ZZ.hh`（整数除法向下取整） |
| `QQ` 算术（+, -, *, /, >=） | `clpoly/number/QQ.hh` |

---

## 5. M4 `__extract_candidates`

### 函数签名

```cpp
// 从 LLL 幺模变换矩阵 U 的短向量行中，提取候选因子对应的 active 因子子集。
// 参数：
//   short_rows — M3 返回的短向量行下标
//   U          — M3 输出的幺模变换矩阵（n×n）
//   r          — 当前 active 因子数（使用 U 的前 r 列）
// 返回：候选子集列表，每项为 active 因子的下标集合（可能为空）
std::vector<std::vector<int>>
__extract_candidates(
    const std::vector<int>& short_rows,
    const LLLMatrix&        U,
    int                     r);
```

### 实现步骤

```cpp
if (short_rows.empty()) return {};

// 构造 U_short：short_rows 中每行的前 r 列
// U_short[k][j] = U[short_rows[k]][j]，j < r
int s = short_rows.size();
std::vector<std::vector<ZZ>> U_short(s, std::vector<ZZ>(r));
for (int k = 0; k < s; ++k)
    for (int j = 0; j < r; ++j)
        U_short[k][j] = U[short_rows[k]][j];

// 列相等性分组：part[j] = 该列的等价类编号
// 两列相等 ⟺ 对所有 k：U_short[k][j1] == U_short[k][j2]
std::vector<int> part(r, -1);
int num_classes = 0;
for (int j = 0; j < r; ++j) {
    if (part[j] != -1) continue;
    part[j] = num_classes;
    for (int j2 = j+1; j2 < r; ++j2) {
        if (part[j2] != -1) continue;
        bool equal = true;
        for (int k = 0; k < s; ++k)
            if (U_short[k][j] != U_short[k][j2]) { equal = false; break; }
        if (equal)
            part[j2] = num_classes;
    }
    ++num_classes;
}

// 按等价类编号收集候选子集
std::vector<std::vector<int>> candidates(num_classes);
for (int j = 0; j < r; ++j)
    candidates[part[j]].push_back(j);

return candidates;
```

> **实现注**：`fmpz_mat_col_partition` 的语义（FLINT 文档）与上述逻辑完全对应。若等价类全为单元素（每个因子独立），则各候选子集大小为 1，退化为逐因子检验，仍然正确。

---

## 6. M5 `__vanhoeij_recombine`

### 函数签名

```cpp
// van Hoeij LLL 因子重组主控循环。
// 前置条件：f 本原，deg(f) >= 2，lc(f) > 0；
//           lifted 是 f 的 Hensel 提升因子（mod m），|lifted| >= 2；
//           m > 2·lc(f)·B_Mig(f)。
std::vector<upolynomial_<ZZ>>
__vanhoeij_recombine(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ&                            m);
```

### 控制流与局部变量

```cpp
int r = (int)lifted.size();
int N = (int)get_deg(f);

// === 初始化参数 ===
int U_exp = (int)ZZ(r > 20 ? r : 20).sizeinbase(2);
ZZ  B     = ZZ(r + 1) * (ZZ(1) << (2 * U_exp));
int J_max = (N + 1) / 2;
int J0    = (3 * r > N + 1) ? 30 : 10;
J0        = std::min(J0, J_max);

// === 活跃集合（下标到 lifted 数组的映射）===
std::vector<int> active(r);                    // active[k] = lifted 中的原始下标
std::iota(active.begin(), active.end(), 0);

upolynomial_<ZZ>              f_star = f;
std::vector<upolynomial_<ZZ>> result;

// === 构造初始格矩阵 ===
auto make_initial_M = [](int rr, int U_exp_) -> LLLMatrix {
    LLLMatrix M(rr, std::vector<ZZ>(rr, ZZ(0)));
    ZZ scale = ZZ(1) << U_exp_;
    for (int i = 0; i < rr; ++i)
        M[i][i] = scale;
    return M;
};

LLLMatrix M = make_initial_M(r, U_exp);
int J_cur    = 0;
int J_target = J0;
```

### 主循环

```cpp
while ((int)active.size() > 1) {

    // [M1] 计算当前 active 因子的 CLD 多项式
    std::vector<upolynomial_<ZZ>> active_lifted;
    for (int k : active)
        active_lifted.push_back(lifted[k]);
    auto cld = __cld_polys(f_star, active_lifted, m);

    // [M2] 喂入列（r 由 M2 内部从 cld.size() 推导）
    int J_new = __build_cld_matrix(M, cld, J_cur, J_target, m);
    J_cur += J_new;

    // [M3] LLL 规约
    LLLMatrix U;
    auto short_rows = __lll_reduce(M, U, B);

    // [M4] 提取候选子集
    auto candidates = __extract_candidates(short_rows, U, (int)active.size());

    // === 候选子集验证 ===
    bool found_any = false;
    for (auto& cand : candidates) {
        if (cand.empty() || (int)cand.size() >= (int)active.size())
            continue;

        // 构造候选因子，与 §7.4 __subset_product_mod（line 564）逻辑完全一致：
        //   prod = lc_f × ∏ lifted[k]，每步 × 后 normalization + mod_coeff，最终 symmetric_mod
        // operator*(upolynomial_<ZZ>, upolynomial_<ZZ>) 通过 basic_polynomial 模板提供，已验证存在。
        ZZ lc_fstar = f_star.front().second;
        upolynomial_<ZZ> g_trial;
        g_trial.push_back({umonomial(0), lc_fstar});
        for (int k : cand) {
            g_trial = g_trial * active_lifted[k];
            g_trial.normalization();
            __upoly_mod_coeff(g_trial, m);
        }
        g_trial = __upoly_symmetric_mod(g_trial, m);

        auto [c_g, pp_g] = __upoly_primitive(std::move(g_trial));

        // 试除
        upolynomial_<ZZ> q_trial, r_trial;
        pair_vec_div(q_trial.data(), r_trial.data(),
                     f_star.data(), pp_g.data(), f_star.comp());

        if (!r_trial.empty()) continue;

        // 成功：记录因子
        result.push_back(std::move(pp_g));
        auto [c_q, pp_q] = __upoly_primitive(std::move(q_trial));
        f_star = std::move(pp_q);

        // 缩减 active：按逆序删除（从大下标到小下标，避免 erase 后下标失效）
        std::vector<int> cand_desc = cand;
        std::sort(cand_desc.rbegin(), cand_desc.rend());
        for (int j : cand_desc)
            active.erase(active.begin() + j);

        // 重建格矩阵参数
        int r_new   = (int)active.size();
        int U_exp_n = (int)ZZ(r_new > 20 ? r_new : 20).sizeinbase(2);
        B           = ZZ(r_new + 1) * (ZZ(1) << (2 * U_exp_n));
        M           = make_initial_M(r_new, U_exp_n);
        J_cur       = 0;
        J_target    = J0;
        found_any   = true;
        break;  // 重新开始循环
    }

    if (found_any) continue;

    // 本轮无新因子：若 J_new == 0（列全被过滤），跳过 LLL 是 future 优化；
    // 当前实现倍增 J_target 后下一轮重试（正确但多做了一次无效 LLL）。
    J_target *= 2;
    if (J_target > J_max) {
        // Fallback：对剩余 f_star 用 Zassenhaus
        // 注：active_lifted 在本轮循环开头已构造完毕，直接复用
        auto zass = __zassenhaus_recombine(f_star, active_lifted, m);
        for (auto& g : zass)
            result.push_back(std::move(g));
        return result;
    }
}

// active 剩余 <= 1：若还有 f_star 未分解
if (!f_star.empty() && get_deg(f_star) > 0)
    result.push_back(std::move(f_star));

std::sort(result.begin(), result.end(),
    [](const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b) {
        return get_deg(a) < get_deg(b);
    });
return result;
```

> **候选验证说明**：上述构造与 Zassenhaus 验证逻辑完全相同（乘积 → 本原化 → `pair_vec_div`），可以从现有 `__factor_recombine` 中直接复用验证代码，提取为小型 lambda 或内联。

---

## 7. 路由层 `__factor_recombine`

### 重命名现有函数

将现有 `§7.5 __factor_recombine` 的函数体不变，仅改名为 `__zassenhaus_recombine`：

```cpp
// 原 §7.5（仅改名，内容不变）
inline std::vector<upolynomial_<ZZ>>
__zassenhaus_recombine(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ&                            m);
```

### 新路由函数

```cpp
// §7.6 因子重组路由（van Hoeij 或 Zassenhaus）
inline std::vector<upolynomial_<ZZ>>
__factor_recombine(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ&                            m)
{
    if ((int)lifted.size() <= ZASSENHAUS_THRESHOLD)
        return __zassenhaus_recombine(f, lifted, m);
    else
        return __vanhoeij_recombine(f, lifted, m);
}
```

`__factor_squarefree_primitive_ZZ` 调用 `__factor_recombine` 的代码不变，路由透明。

---

## 8. 错误处理策略

| 情形 | 处理方式 |
|------|---------|
| M1：`r_i` 非空（精确除法失败） | `assert(r_i.empty())`（DEBUG），precondition 违反，不处理 |
| M2：`J_new == 0`（全列过滤） | 调用方（M5）按"本轮无进展"处理，倍增 J_target；实践中 Hensel 精度充足时极罕见 |
| M3：矩阵奇异（B_gs[k] == 0） | 跳过对应 mu 更新（除零保护，不影响正确性：奇异说明重复因子，已在上层排除） |
| M5：`J_target > J_max` | 对剩余因子调用 `__zassenhaus_recombine`（正确性保证，性能退化但有界） |
| M4：`short_rows` 为空 | 返回空 `candidates`，M5 倍增 J_target（正常情形） |

---

## 9. 复用点总览

| 现有函数 | 来源 | M1-M5 使用场景 |
|---------|------|----------------|
| `__upoly_divmod_mod` | 同文件 §4.1.4 | M1：精确除法 q_i = f/h_i mod m |
| `derivative(f)` | `upolynomial.hh:176` | M1：形式导数 |
| `__upoly_mod_coeff` | 同文件 §4.1.4 | M1、M5：系数取模 |
| `__upoly_mul_mod` | 同文件 §6 | M1：C_i 乘法 |
| `__upoly_symmetric_mod` | 同文件 §4.2.2 | M1、M5：对称约化 |
| `abs(const ZZ&)` | `ZZ.hh:904`（自由函数，friend of ZZ） | M2 列过滤中计算系数绝对值（CLD 系数已对称约化，直接取 abs） |
| `__upoly_primitive` | 同文件 §4.5 | M5：候选因子本原化 |
| `pair_vec_div` | `basic.hh:568` | M5：候选因子试除 |
| `get_deg` | `upolynomial.hh` | M2/M5：次数查询 |
| `ZZ::sizeinbase(2)` | `ZZ.hh` | M5：U_exp 计算（bit 长度） |
| `upolynomial_<ZZ>::normalization()` | `upolynomial.hh` | M5：规范化 |
| `operator*(upolynomial_, upolynomial_)` | `polynomial.hh`（via `basic_polynomial`） | M5：g_trial 乘法（同 `__subset_product_mod` §7.4 第574行） |
| `ZZ::fdiv_q(ZZ&, const ZZ&, const ZZ&)` | `ZZ.hh:633`（静态方法，输出参数形式） | M3：round_qq 中的 floor 除法 |
