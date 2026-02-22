# P1b 线性 Hensel 提升：细化设计文档

> 阶段：细化（workflow.md §2.3）
> 基础文档：`docs/design/vanhoeij/architecture.md` Part B
> 调研报告：`docs/research/linear-hensel-research.md`
> 目标文件：`clpoly/polynomial_factorize_univar.hh`
> 方案：**n 因子直接线性步**（Monagan 2022 路线，无二叉树）

---

## 0. 实现顺序

按依赖从底向上：

```
M4 __heuristic_starting_precision   ← 无依赖，纯计算
M0 __linear_bezout_chain            ← 依赖 Zp 多项式 GCD（现有）
M1 __hensel_step_linear_nfactor     ← 依赖 __upoly_divmod_mod 等（现有）
M3 __linear_hensel_lift_with_lll    ← 调用 M0/M1/M4 + P1a 全部模块
接口修改                             ← 改 __factor_squarefree_primitive_ZZ
```

每个模块完成后独立测试再进入下一个（workflow.md §4.1）。

---

## 1. M4：`__heuristic_starting_precision`

### 函数签名

```cpp
inline int __heuristic_starting_precision(
    const upolynomial_<ZZ>& f,
    int                     r,
    uint32_t                p);
```

### 参数语义

| 参数 | 含义 |
|------|------|
| `f` | 待分解多项式（用于取 `f.size()` = 项数） |
| `r` | 模因子数 |
| `p` | 素数（Hensel 底） |

### 算法步骤

```cpp
// 1. 用 double 计算 FLINT 启发式 a_h
double logp = std::log((double)p);
int    min_b = (int)ZZ(p).sizeinbase(2);          // ≈ log₂(p)
int    N     = (int)f.size() - 1;                  // 项数 - 1 ≈ deg(f)
double a_h_d = std::ceil(
    (2.5 * r + min_b) * std::log(2.0) / logp
    + std::log((double)(N + 1)) / (2.0 * logp));
int a_h = std::max(1, (int)a_h_d);

// 2. 计算 Mignotte 精度 a_mig（循环逼近）
ZZ B_mig = __mignotte_bound(f);
ZZ lc_f  = f.front().second;
if (lc_f < ZZ(0)) lc_f = -lc_f;
ZZ target = ZZ(2) * lc_f * B_mig;
int a_mig = 0;
ZZ  pa(1);
while (pa <= target) { pa *= ZZ(p); ++a_mig; }

// 3. 取较小值（保证不超出 Mignotte 精度）
return std::min(a_mig, a_h);
```

### 复用点

- `__mignotte_bound(f)`（line 122）
- `ZZ::sizeinbase(2)`（GMP `mpz_sizeinbase`，bit 宽度）

---

## 2. M0：`__linear_bezout_chain`

### 函数签名

```cpp
inline std::vector<upolynomial_<ZZ>>
__linear_bezout_chain(
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t                             p);
```

### 参数语义

| 参数 | 含义 |
|------|------|
| `factors` | lc_f 调整后的 r 个 Zp 因子（首一或 lc_f·首一） |
| `p` | 素数 |

### 返回值

`vector<upolynomial_<ZZ>>`：r 个 Bézout 系数 s_i，系数以 ZZ 存储（值在 `[0, p)`），满足：

```
Σᵢ s[i] · ∏_{j≠i} factors[j] ≡ 1  (mod p)
deg s[i] < deg factors[i]
```

### 算法步骤

```cpp
int r = (int)factors.size();

// 步骤 1：计算全积 P = ∏ factors[i]  (mod p)
upolynomial_<Zp> P = factors[0];
for (int i = 1; i < r; ++i)
    P = P * factors[i];   // Zp 多项式乘法，结果自动 mod p

// 步骤 2：对每个 i 计算 s[i]
std::vector<upolynomial_<ZZ>> result;
result.reserve(r);
for (int i = 0; i < r; ++i) {
    // 2a. 计算余积 P_i = P / factors[i]  (mod p 精确除法)
    upolynomial_<Zp> P_i = poly_exact_div_zp(P, factors[i]);
    //    poly_exact_div_zp：__upoly_divmod_mod 的 Zp 版本，余式应为空

    // 2b. 扩展 GCD：找 a_i 使 a_i * P_i ≡ 1 (mod factors[i], mod p)
    //     使用现有 polynomial_GCD(P_i, factors[i], s_i, t_i)
    //     其返回值满足 s_i * P_i + t_i * factors[i] = gcd = 1（因互素）
    upolynomial_<Zp> s_i_zp, t_i_zp;
    polynomial_GCD(P_i, factors[i], s_i_zp, t_i_zp);
    // s_i_zp 已满足 deg s_i_zp < deg factors[i]（扩展 GCD 的标准输出）

    // 2c. 转为 ZZ（系数 ∈ [0, p)）
    result.push_back(__upoly_Zp_to_ZZ(s_i_zp));
}
return result;
```

**关键说明**：
- `polynomial_GCD(g, h, s, t)` 是 CLPoly 现有 Zp 扩展 GCD（`__hensel_tree_build_recursive` line 314 已使用）。
- 返回值满足 `s·g + t·h = gcd`，当 gcd = 1（因子互素）时即 `s·P_i + t·factors[i] = 1`，故 `s_i_zp = s`。
- 无需再做 `s_i mod factors[i]`——扩展 GCD 已保证 `deg s_i < deg h`（即 `deg factors[i]`）。

### 复用点

- `polynomial_GCD`（已有，Zp 扩展 GCD，`__hensel_tree_build_recursive` line 314 已用）
- `__upoly_Zp_to_ZZ`（已有）
- Zp 多项式乘法（`*` 运算符，已有）

### 错误处理

- 若 gcd ≠ 1（因子不互素）：`polynomial_GCD` 返回非平凡 gcd，此时 `s_i_zp` 不满足需求，后续提升将产生错误结果。由 `__select_prime` 保证不发生，调用方无需额外检查。

---

## 3. M1：`__hensel_step_linear_nfactor`

### 函数签名

```cpp
inline void __hensel_step_linear_nfactor(
    std::vector<upolynomial_<ZZ>>&       h,
    const std::vector<upolynomial_<ZZ>>& s,
    const upolynomial_<ZZ>&              f,
    uint32_t                             p,
    const ZZ&                            p_a);
```

### 参数语义

| 参数 | 含义 |
|------|------|
| `h` | r 个提升因子（原地更新），系数 ∈ `[0, p^a)` |
| `s` | r 个 Bézout 系数（只读），系数 ∈ `[0, p)` |
| `f` | 目标多项式（ZZ，完整） |
| `p` | 素数 |
| `p_a` | 当前模数 p^a |

### 算法步骤

```cpp
ZZ p_zz(p);
ZZ p_a1 = p_a * p_zz;    // p^(a+1)

// ── 步骤 1：计算 ∏ h[i]  (mod p^(a+1))，然后算误差 e ────────────
// 为避免全精度大整数乘法，只需 mod p^(a+1) 精度
upolynomial_<ZZ> prod = h[0];
__upoly_mod_coeff(prod, p_a1);
for (int i = 1; i < (int)h.size(); ++i) {
    prod = prod * h[i];
    __upoly_mod_coeff(prod, p_a1);   // 每步截断，控制系数增长
}

// e = (f - prod) / p_a  mod p
upolynomial_<ZZ> err = f - prod;
for (auto& term : err.data()) {
    ZZ::fdiv_q(term.second, term.second, p_a);   // 精确整除（Hensel 不变式）
    ZZ::fdiv_r(term.second, term.second, p_zz);  // mod p → [0, p)
}
// 去零项
{ auto it = err.data().begin(), out = it;
  for (; it != err.data().end(); ++it)
      if (it->second) { if (out != it) *out = std::move(*it); ++out; }
  err.data().erase(out, err.data().end()); }

if (err.empty()) return;    // f ≡ ∏ h[i]  (mod p^(a+1))，无需修正

// ── 步骤 2：对每个 i 独立计算 σ[i] 并更新 h[i] ──────────────────
for (int i = 0; i < (int)h.size(); ++i) {
    // σ[i] = s[i] * e  mod h[i]  (mod p)
    upolynomial_<ZZ> se = s[i] * err;
    // h[i] 系数 ∈ [0, p^a)，divmod 内部会先 mod p 处理被除式和除式首项
    upolynomial_<ZZ> q, sigma;
    __upoly_divmod_mod(q, sigma, se, h[i], p_zz);

    // h'[i] = h[i] + p_a * σ[i]  (mod p^(a+1))
    for (auto& term : sigma.data())
        term.second *= p_a;
    h[i] = h[i] + sigma;
    __upoly_mod_coeff(h[i], p_a1);
}
```

### 关键实现说明

**产品截断（步骤 1）**：每次中间乘积后立即 mod p^(a+1)，使系数不超过 p^(a+1)。这是正确的，因为：
- 我们只需要 `prod mod p^(a+1)` 来计算误差 `(f - prod) / p_a mod p`
- `(f - prod) mod p^(a+1)` 的每个系数都被 p^a 整除（Hensel 不变式），mod p^(a+1) 后再除以 p^a 才得到 e

**`__upoly_divmod_mod(q, sigma, se, h[i], p_zz)` 中 `h[i]` 的首项**：
- `h[i]` 的系数 ∈ `[0, p^a)`，首项系数 mod p = 原始 Zp 因子的首项系数（= 1 对于首一因子，或 = lc_f·1 对于 factors_adj[0]）
- `__upoly_divmod_mod` 内部计算 `ZZ::invert(lc_inv, h[i].front().second, p_zz)`，由 `__select_prime` 保证可逆

**s[i] 系数小**：s[i] 系数 ∈ `[0, p)`，err 系数 ∈ `[0, p)`，故 `se = s[i] * err` 的系数 ∈ `[0, p²)`，不溢出 int64_t（p < 2^30，p² < 2^60 < 2^63）。

### 复用点

- `__upoly_mod_coeff`（line 137）：非对称系数取模
- `__upoly_divmod_mod`（line 155）：模 p 带余除法（已处理大系数除式）
- `ZZ::fdiv_q`, `ZZ::fdiv_r`（ZZ.hh）

---

## 4. M3：`__linear_hensel_lift_with_lll`

### 函数签名

```cpp
inline std::vector<upolynomial_<ZZ>>
__linear_hensel_lift_with_lll(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t                             p);
```

### 算法步骤（分段）

#### 段 A：初始化

```cpp
assert(factors.size() >= 2 && !f.empty());
int r = (int)factors.size();

// A1. 首项系数处理（与现有 __hensel_lift line 511-515 相同）
ZZ lc_f = f.front().second;
if (lc_f < ZZ(0)) lc_f = -lc_f;
std::vector<upolynomial_<Zp>> factors_adj = factors;
Zp lc_mod_p(f.front().second, p);
for (auto& term : factors_adj[0])
    term.second *= lc_mod_p;
factors_adj[0].normalization();

// A2. 预计算 Bézout 链（一次性，固定 mod p）
auto s = __linear_bezout_chain(factors_adj, p);    // [M0]

// A3. 初始提升基：Zp → ZZ
std::vector<upolynomial_<ZZ>> h;
h.reserve(r);
for (auto& fi : factors_adj)
    h.push_back(__upoly_Zp_to_ZZ(fi));

// A4. 精度参数
int a0    = __heuristic_starting_precision(f, r, p);  // [M4]
int a_mig = 0;
{ ZZ pa(1), tgt = ZZ(2) * lc_f * __mignotte_bound(f);
  while (pa <= tgt) { pa *= ZZ(p); ++a_mig; } }

// A5. van Hoeij LLL 状态（与 __vanhoeij_recombine 相同）
int N     = (int)get_deg(f);
int U_exp = (int)ZZ(r > 20 ? r : 20).sizeinbase(2);
ZZ  B     = ZZ(r + 1) * (ZZ(1) << (2 * U_exp));
int J_max = (N + 1) / 2;
int J0    = std::min((3 * r > N + 1) ? 30 : 10, J_max);

std::vector<int>              active(r);
std::iota(active.begin(), active.end(), 0);
upolynomial_<ZZ>              f_star = f;
std::vector<upolynomial_<ZZ>> result;

auto make_initial_M = [](int rr, int u_exp) -> LLLMatrix {
    LLLMatrix M(rr, std::vector<ZZ>(rr, ZZ(0)));
    ZZ scale = ZZ(1) << u_exp;
    for (int i = 0; i < rr; ++i) M[i][i] = scale;
    return M;
};

LLLMatrix M      = make_initial_M(r, U_exp);
int       J_cur  = 0, J_target = J0;
```

#### 段 B：初始线性提升（a0 步）

```cpp
ZZ p_a(p);      // 开始时 p^1（第一步的 "当前模数 p^a" 是 p^1，提升后到 p^2）
int a_cur = 1;

for (int step = 1; step <= a0; ++step) {
    __hensel_step_linear_nfactor(h, s, f, p, p_a);   // [M1]
    p_a   *= ZZ(p);    // p_a → p^(step+1)
    a_cur  = step + 1;
}
// 段 B 结束：h[i] 系数 mod p^a0，p_a = p^a0，a_cur = a0
```

**注**：第一次调用 M1 时 `p_a = p^1`，步骤完成后 `h[i]` 升至 mod p^2，然后 `p_a *= p` 变为 p^2，为下一步做准备。

#### 段 C：主交织循环

```cpp
while ((int)active.size() > 1) {

    // C1. 收集当前活跃因子（对称约化供 P1a 模块使用）
    std::vector<upolynomial_<ZZ>> active_h;
    for (int k : active) {
        auto hi = __upoly_symmetric_mod(h[k], p_a);
        active_h.push_back(std::move(hi));
    }

    // C2. P1a 模块调用（与 __vanhoeij_recombine 完全相同）
    auto cld        = __cld_polys(f_star, active_h, p_a);
    int  J_new      = __build_cld_matrix(M, cld, J_cur, J_target, p_a);
    J_cur          += J_new;
    LLLMatrix U;
    auto short_rows = __lll_reduce(M, U, B);
    auto candidates = __extract_candidates(short_rows, U, (int)active.size());

    // C3. 因子验证（与 __vanhoeij_recombine line 1086-1133 完全相同）
    bool found_any = false;
    for (auto& cand : candidates) {
        if (cand.empty() || (int)cand.size() >= (int)active.size()) continue;
        ZZ lc_fstar = f_star.front().second;
        upolynomial_<ZZ> g_trial;
        g_trial.push_back({umonomial(0), lc_fstar});
        for (int k : cand) {
            g_trial = g_trial * active_h[k];
            g_trial.normalization();
            __upoly_mod_coeff(g_trial, p_a);
        }
        g_trial = __upoly_symmetric_mod(g_trial, p_a);
        auto [c_g, pp_g] = __upoly_primitive(std::move(g_trial));
        upolynomial_<ZZ> q_trial, r_trial;
        pair_vec_div(q_trial.data(), r_trial.data(),
                     f_star.data(), pp_g.data(), f_star.comp());
        if (!r_trial.empty()) continue;
        // 成功
        result.push_back(std::move(pp_g));
        auto [c_q, pp_q] = __upoly_primitive(std::move(q_trial));
        f_star = std::move(pp_q);
        std::vector<int> cand_desc = cand;
        std::sort(cand_desc.rbegin(), cand_desc.rend());
        for (int j : cand_desc) active.erase(active.begin() + j);
        int r_new   = (int)active.size();
        int U_exp_n = (int)ZZ(r_new > 20 ? r_new : 20).sizeinbase(2);
        B           = ZZ(r_new + 1) * (ZZ(1) << (2 * U_exp_n));
        M           = make_initial_M(r_new, U_exp_n);
        J_cur       = 0; J_target = J0;
        found_any   = true;
        break;
    }
    if (found_any) continue;

    // C4. 本轮无进展：先倍增 J_target
    if (J_new > 0) {
        J_target *= 2;
        if (J_target > J_max) {
            // 安全网：fallback Zassenhaus
            std::vector<upolynomial_<ZZ>> active_h_sym;
            for (int k : active)
                active_h_sym.push_back(__upoly_symmetric_mod(h[k], p_a));
            auto zass = __zassenhaus_recombine(f_star, active_h_sym, p_a);
            for (auto& g : zass) result.push_back(std::move(g));
            goto sort_and_return;
        }
        continue;
    }

    // C5. 列已耗尽：增量提升 k 步
    {
        static constexpr int LIFT_STEP_K = 5;
        for (int step = 0; step < LIFT_STEP_K && a_cur <= a_mig; ++step) {
            __hensel_step_linear_nfactor(h, s, f, p, p_a);
            p_a   *= ZZ(p);
            a_cur += 1;
        }
        if (a_cur > a_mig) {
            // 安全网：fallback Zassenhaus
            std::vector<upolynomial_<ZZ>> active_h_sym;
            for (int k : active)
                active_h_sym.push_back(__upoly_symmetric_mod(h[k], p_a));
            auto zass = __zassenhaus_recombine(f_star, active_h_sym, p_a);
            for (auto& g : zass) result.push_back(std::move(g));
            goto sort_and_return;
        }
        // 重置格矩阵（精度改变，CLD 需重新计算）
        J_cur    = 0;
        J_target = J0;
        M        = make_initial_M((int)active.size(), U_exp);
    }
}

// D. 最后一个因子
if (!f_star.empty() && get_deg(f_star) > 0)
    result.push_back(std::move(f_star));

sort_and_return:
std::sort(result.begin(), result.end(),
    [](const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b) {
        return get_deg(a) < get_deg(b); });
return result;
```

**关于 `goto`**：两个安全网出口处使用 `goto sort_and_return` 以避免代码重复。若风格偏好不用 goto，可改用辅助函数或 lambda 封装排序+返回逻辑。

### 调用关系

```
__linear_hensel_lift_with_lll (M3)
  ├── __linear_bezout_chain (M0)
  ├── __heuristic_starting_precision (M4)
  ├── __upoly_Zp_to_ZZ                    ← 现有
  ├── __hensel_step_linear_nfactor (M1)   ← 段 B 初始提升 + 段 C5 增量提升
  ├── __upoly_symmetric_mod               ← 现有
  ├── __cld_polys         (P1a M1)
  ├── __build_cld_matrix  (P1a M2)
  ├── __lll_reduce        (P1a M3)
  ├── __extract_candidates (P1a M4)
  ├── __upoly_primitive                   ← 现有
  ├── pair_vec_div                        ← 现有
  └── __zassenhaus_recombine              ← 现有（安全网）
```

---

## 5. 接口修改：`__factor_squarefree_primitive_ZZ`

### 修改位置

`polynomial_factorize_univar.hh` line 1302-1305。

### 修改前

```cpp
auto [lifted, modulus] = __hensel_lift(f, sel.factors, sel.prime);
return __factor_recombine(f, lifted, modulus);
```

### 修改后

```cpp
return __linear_hensel_lift_with_lll(f, sel.factors, sel.prime);
```

### 保留原函数

`__hensel_lift`, `__factor_recombine`, `__hensel_tree_build`, `__hensel_step` 均保留不删除（调试对比和 fallback 路径备用）。

---

## 6. 复用点汇总

| 复用函数 | 所在位置 | P1b 使用方式 |
|---------|---------|------------|
| `polynomial_GCD`（Zp 扩展 GCD）| 已有 | M0：计算 Bézout 系数 sᵢ |
| `__upoly_Zp_to_ZZ` | 已有 | M0 输出转换；M3 段 A 初始化 h |
| `__mignotte_bound` | line 122 | M4 + M3 段 A：精度计算 |
| `__upoly_mod_coeff` | line 137 | M1：积截断、h 更新 |
| `__upoly_divmod_mod` | line 155 | M1：mod-p 带余除法（计算 σ） |
| `ZZ::fdiv_q`, `ZZ::fdiv_r` | ZZ.hh | M1 步骤 1：精确除再 mod p |
| `__upoly_symmetric_mod` | 已有 | M3 段 C1：对称约化供 CLD |
| `__upoly_primitive` | 已有 | M3 段 C3：本原化 |
| `pair_vec_div` | 已有 | M3 段 C3：试除 |
| `__cld_polys` | P1a M1 | M3 段 C2 |
| `__build_cld_matrix` | P1a M2 | M3 段 C2 |
| `__lll_reduce` | P1a M3 | M3 段 C2 |
| `__extract_candidates` | P1a M4 | M3 段 C2 |
| `__zassenhaus_recombine` | 已有 | M3 安全网 |
| `ZZ(1) << u_exp` | ZZ.hh | M3 格矩阵初始化 |

**不再依赖**：`__hensel_tree_build`, `__hensel_node`, `__hensel_step`, `__hensel_lift_recursive`, `__hensel_extract_factors`（二叉树结构全部从 P1b 路径中移除）。

---

## 7. 测试策略

### T1. M0 单元测试

```
输入：factors = [(x-1) mod 5, (x+1) mod 5, (x-2) mod 5]，p = 5
手动验证：Σ s[i] * ∏_{j≠i} factors[j] ≡ 1 (mod 5)
代码检查：打印每个 s[i]，手动计算乘积之和
```

### T2. M1 单元测试

```
选取小例子：f = (x-1)(x+1)(x-2) = x³-2x²-x+2，p = 5
初始因子 mod 5：h = [(x-1),(x+1),(x-2)]，s = M0 的输出
手动计算一步线性提升，验证 ∏ h'[i] ≡ f (mod 25)
```

### T3. M3 正确性测试（与 P1a 结果比对）

```
策略：对同一组多项式，分别运行
  - 旧路径：__hensel_lift + __factor_recombine
  - 新路径：__linear_hensel_lift_with_lll
比较因子集合是否相同（多重集合）
测试集：test_factorize.cc 全部 256 个用例
命令：make test/test_factorize
```

### T4. 性能目标

```
目标：(x-1)...(x-70) < 20 ms（P1a+P1b 合计，对比 P1a 前的 83 ms）
命令：make bench-clpoly
```

### T5. 压力与回归测试

```
make stress    ← release 模式全部通过
make crosscheck ← FLINT/NTL 交叉验证
```

---

## 8. 错误处理策略

| 情形 | 处理方式 |
|------|---------|
| `factors.size() <= 1` 或 `irreducible` | 由 `__factor_squarefree_primitive_ZZ` 在调用前拦截（不变） |
| M0 中因子不互素（gcd ≠ 1）| `__select_prime` 保证不发生；如发生，后续 s[i] 不满足 Diophantine，提升失败，安全网 fallback |
| M1 中 `fdiv_q` 非精确整除 | 仅通过 assert 发现（Hensel 不变式由调用方保证） |
| M1 中 `__upoly_divmod_mod` 首项不可逆 | `__select_prime` 保证素数不整除各因子首项；如触发断言，属于 prime 选取 bug |
| a_cur > a_mig | fallback Zassenhaus（正确性保底，理论上不触发） |
| J_target > J_max | fallback Zassenhaus（正确性保底） |

---

## 9. TODO：Monagan 2022 矩阵优化（v2）

当前 M1 对每步 a → a+1 独立计算 σ[i]，总代价为 O(a_stop · r · n²/p)。

Monagan 2022 通过将全部 D 步的修正序列 {σ[i]^(1), ..., σ[i]^(D)} 组织为**多项式矩阵乘法**，把总 Bézout 更新代价从 O(r · n² · D) 降至 O(r · n · D²/r) = O(n · D²)（"立方代价"）。

此优化与早期终止结合后（a_stop ≪ D_Mig），收益有限。建议在 P1b v1 性能验证后，若 uni-70 仍有 >2x 差距再考虑实施。
