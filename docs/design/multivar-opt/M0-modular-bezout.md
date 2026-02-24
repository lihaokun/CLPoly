# M0 细化设计：模 Bézout 链（消除系数爆炸）

> 对应里程碑：M0（可选快速里程碑）
> 文件：`clpoly/polynomial_factorize_wang.hh`
> 改动范围：约 80 行删除 + 35 行新增，其余不变（所有辅助函数均已存在，无需新增）
> 验收：bivar-70 < 1s；test_factorize_multivar 全部通过；test_multivar_hensel 全部通过

---

## 1. 根因回顾

`__multivar_hensel_lift`（L836-865）在 **Z[x]** 上构建 Bézout 链，对 r=70 个线性因子：

| 量 | 值 | 原因 |
|---|---|---|
| `bezout_s[0]` 次数 | ~2346 | 每步累乘 beta，degree 累加 |
| `bezout_s[0]` 系数位数 | ~6000 bit | denom = Vandermonde det(k1,...,k35) |
| 每次 Diophantine 调用 | ~1亿次大整数运算 | `si_h = bezout_s[0] * h_upoly`（L499） |

**根本原因**：Bézout 链在 **Z[x]** 上运算。改到 **Fp[x]** 后，系数恒在 `[0, p)` 范围内，自然消除爆炸。

---

## 2. 核心思路

**改前**：
```
bezout chain in Z[x]: Σ si·V̂i = denom  (denom ~2^6000)
Diophantine base: si_h = bezout_s[i] * h  (巨大整数运算)
                  rem = upoly_prem(si_h, vi)  (伪除法)
                  δi = rem / (denom · lc(vi)^k)  (大整数除法)
```

**改后**：
```
bezout chain in Fp[x]: Σ si_p·V̂i_p ≡ 1 (mod p)  (denom = 1 in Fp!)
Diophantine base: h_p = reduce(h, p)
                  si_h_p = si_p · h_p  mod p   (小整数乘法)
                  rem_p = si_h_p rem vi_p  mod p  (Fp 精确除法，无伪除法！)
                  δi = symmetric_mod(rem_p, pa)   (对称约化到 Z)
```

**正确性**：设 `pa = p^a > 2·B`（B 为 δi 的系数界），则 `symmetric_mod(rem_p, pa)` 唯一恢复 Z 上的真实解。

---

## 3. 代码改动清单

### 3.1 新增：`__bezout_chain_modular`（约 55 行）

位置：`__multivar_hensel_lift` 函数之前（L820 前）。

```cpp
// Fp[x] Bézout 链: 计算 si_p ∈ Fp[x] 使得 Σ si_p·V̂i_p ≡ 1 (mod p)
// 同时计算 pa = p^a > 2·B（B = 系数界）
// 输出 bezout_s_p：系数以 ZZ 存储但值域在 [0, p)
// 输出 pa：模数；输出 p_out：素数（供 Diophantine 使用）
static void
__bezout_chain_modular(
    const std::vector<upolynomial_<ZZ>>& factors,   // 输入：Z[x] 因子
    uint32_t& p_out,                                 // 输出：所选素数
    std::vector<upolynomial_<ZZ>>& bezout_s_p,      // 输出：Fp[x] Bézout 系数（值域 [0,p)）
    ZZ& pa)                                          // 输出：模数 p^a > 2B
{
    size_t r = factors.size();

    // ── Step 1: 选素数 p（与所有 lc(vi) 互素）────────────────────────────
    uint32_t p = 0;
    for (unsigned idx = 0; ; ++idx)
    {
        uint32_t candidate = boost::math::prime(idx);
        bool bad = false;
        for (auto& vi : factors)
        {
            ZZ mod_lc;
            ZZ::fdiv_r(mod_lc, vi.front().second, ZZ(candidate));
            if (mod_lc == ZZ(0)) { bad = true; break; }
        }
        if (!bad) { p = candidate; break; }
    }
    p_out = p;

    // ── Step 2: 将 factors 映射到 upolynomial_<Zp> ───────────────────────
    // 使用已有 polynomial_mod(f, p)（upolynomial.hh L104）
    std::vector<upolynomial_<Zp>> vp(r);
    for (size_t i = 0; i < r; ++i) vp[i] = polynomial_mod(factors[i], p);

    // ── Step 3: 在 upolynomial_<Zp> 上迭代构建 Bézout 链 ────────────────
    // 使用已有 polynomial_GCD(F, G, s, t)（polynomial_gcd.hh L719）
    // 该函数返回首一 GCD，并满足 s·F + t·G = gcd（已自动归一化）
    // 若 gcd = 1，则 s·g_acc + t·vp[i] = 1，无需手动缩放
    Zp one_p(1, p);
    std::vector<upolynomial_<Zp>> bezout_s_zp(r);
    upolynomial_<Zp> g_acc_p = vp[0];
    bezout_s_zp[0] = {{umonomial(0), one_p}};   // s_zp[0] = 1

    for (size_t i = 1; i < r; ++i)
    {
        upolynomial_<Zp> s, t;
        // polynomial_GCD 返回首一 gcd；因子两两互素 mod p → gcd = 1
        // → s·g_acc_p + t·vp[i] = 1（polynomial_GCD 已自动归一化，无需再乘逆）
        polynomial_GCD(g_acc_p, vp[i], s, t);

        // 更新已有系数: bezout_s_zp[j] *= t（j < i）
        for (size_t j = 0; j < i; ++j)
            bezout_s_zp[j] = bezout_s_zp[j] * t;
        bezout_s_zp[i] = std::move(s);

        // 积累: g_acc_p *= vp[i]（upolynomial_<Zp> 的 operator* 已定义）
        g_acc_p = g_acc_p * vp[i];
    }

    // ── Step 4: 将 Zp 结果转回 upolynomial_<ZZ>（值在 [0, p)）────────────
    // 使用新增的 poly_convert(upolynomial_<Zp> → upolynomial_<ZZ>)（upolynomial.hh）
    bezout_s_p.resize(r);
    for (size_t i = 0; i < r; ++i)
        poly_convert(bezout_s_zp[i], bezout_s_p[i]);

    // ── Step 5: 计算 pa = p^a > 2B ────────────────────────────────────────
    // B = ∏(||vi||₁ + 1)：Diophantine 解的系数界
    ZZ B(1);
    for (auto& vi : factors)
    {
        ZZ vi_norm(0);
        for (auto& [m, c] : vi) vi_norm += abs(c);
        B *= vi_norm + ZZ(1);
    }
    ZZ two_B = ZZ(2) * B;
    pa = ZZ(1);
    while (pa <= two_B) pa *= ZZ(p);
}
```

**依赖函数**（均已存在，无需新增）：

| 函数 | 来源 | 说明 |
|------|------|------|
| `polynomial_GCD(F, G, s, t)` for `upolynomial_<Zp>` | `polynomial_gcd.hh` L719 | Fp[x] 扩展 GCD，返回首一 gcd，已自动归一化 |
| `operator*` on `upolynomial_<Zp>` | `polynomial_gcd.hh` L737 | Fp[x] 乘法，已用于 GCD 内部 |
| `Zp(int64_t, uint32_t)` | `number.hh` | 构造 Fp 元素 |

> **注**：`polynomial_GCD` for `Zp` 的 prime 嵌入在 `Zp` 类型中（`F.front().second.prime()`），
> 无需额外传递 `p` 参数；转换时使用 `c.value()` 取回整数值。

---

### 3.2 修改：`__multivar_hensel_lift`（L836-909）

**删除**（约 74 行）：整个 Z[x] Bézout 链构建 + pa 计算 + denom 逻辑（L836-909）。

**替换为**（约 8 行）：

```cpp
// ————————————————————————————————————————
// §6.2 Fp[x] Bézout 链（模版本，消除系数爆炸）
// ————————————————————————————————————————
uint32_t p_lift;
std::vector<upolynomial_<ZZ>> bezout_s_p;
ZZ pa;
__bezout_chain_modular(scaled_factors, p_lift, bezout_s_p, pa);
// 注: pa > 0 始终成立; 不再需要 denom 和 Z[x] bezout_s
```

**`__hensel_lift_one_var` 调用处改动**（L965-967）：

```cpp
// 改前:
__hensel_lift_one_var(
    f_curr, G, bezout_s, denom, pa, scaled_factors,
    delta, lc_tau_curr, main_var, xk, alpha_k, dk, prev_eval);

// 改后:
__hensel_lift_one_var(
    f_curr, G, bezout_s_p, p_lift, pa, scaled_factors,
    delta, lc_tau_curr, main_var, xk, alpha_k, dk, prev_eval);
```

---

### 3.3 修改：`__hensel_lift_one_var` 签名（L698-709）

```cpp
// 改前:
void __hensel_lift_one_var(
    ...,
    const std::vector<upolynomial_<ZZ>>& bezout_s,
    const ZZ& bezout_denom,
    const ZZ& pa,
    ...);

// 改后:
void __hensel_lift_one_var(
    ...,
    const std::vector<upolynomial_<ZZ>>& bezout_s_p,  // Fp 版本
    uint32_t p,                                         // 新增
    const ZZ& pa,
    ...);
```

内部仅做参数透传：`__multivar_diophantine` 调用处把新参数传下去（L792-794）。

---

### 3.4 修改：`__multivar_diophantine` 签名 + 基本情形（L474-544）

**签名改动**：

```cpp
// 改前:
__multivar_diophantine(h, G, bezout_s, bezout_denom, v_factors,
                       main_var, eval_vars, pa, depth)

// 改后:
__multivar_diophantine(h, G, bezout_s_p, uint32_t p, v_factors,
                       main_var, eval_vars, pa, depth)
// bezout_denom 参数删除（Fp 下 denom = 1，不再需要）
```

**基本情形重写**（L490-544，约 55 行 → 约 25 行）：

```cpp
if (depth >= eval_vars.size())
{
    std::vector<Poly> result(r, Poly(comp_ptr));
    upolynomial_<ZZ> h_upoly;
    poly_convert(h, h_upoly);

    // ── 改前（Z[x] 伪除法，60 行大整数运算）──
    // auto si_h = bezout_s[i] * h_upoly;         ← 天文数字
    // upoly_prem(rem, si_h, v_factors[i], ...);   ← 海量计算
    // ZZ divisor = bezout_denom * lc_pow;         ← 6000-bit 除法
    // ZZ::invert(divisor_inv, divisor, pa);

    // ── 改后（Fp[x] 精确除法）──────────────────
    for (size_t i = 0; i < r; ++i)
    {
        // 1. 将 h 和 vi 转为 upolynomial_<Zp>
        //    使用已有 polynomial_mod(f, p)（upolynomial.hh L104）
        auto h_p  = polynomial_mod(h_upoly, p);
        auto vi_p = polynomial_mod(v_factors[i], p);

        // 2. si_p · h_p in Fp[x]（bezout_s_p[i] 系数在 [0,p)，先转 Zp）
        auto si_p   = polynomial_mod(bezout_s_p[i], p);
        auto si_h_p = si_p * h_p;                // upolynomial_<Zp> operator*

        // 3. (si_p · h_p) rem vi_p in Fp[x]
        //    使用已有 __upoly_mod(f, g)（polynomial_factorize_zp.hh L39）
        auto rem_p = __upoly_mod(si_h_p, vi_p);

        // 4. Zp → ZZ，再对称约化（pa > 2·||δi||_∞ 保证唯一性）
        upolynomial_<ZZ> rem_zz;
        poly_convert(rem_p, rem_zz);         // 新增 poly_convert（upolynomial.hh）
        upolynomial_<ZZ> delta_i_upoly;
        for (auto& [m, c] : rem_zz)
        {
            ZZ coeff = __symmetric_mod(c, pa);
            if (coeff != ZZ(0))
                delta_i_upoly.push_back({m, coeff});
        }

        poly_convert(delta_i_upoly, result[i], main_var);
    }
    return result;
}
```

**递归情形**（L547 以后）：仅需把签名中的 `bezout_denom` 替换为 `p`，其余不变。

---

### 3.5 辅助函数

**已有，直接使用**：

| 函数 | 来源 | 用途 |
|------|------|------|
| `polynomial_mod(f, p)` for `upolynomial_<ZZ>` | `upolynomial.hh` L104 | ZZ→Zp 转换 |
| `polynomial_GCD(F, G, s, t)` for `upolynomial_<Zp>` | `polynomial_gcd.hh` L719 | Fp[x] 扩展 GCD |
| `operator*` on `upolynomial_<Zp>` | 已有（L737 中已使用） | Fp[x] 乘法 |
| `__upoly_mod(f, g)` | `polynomial_factorize_zp.hh` L39 | Fp[x] 取余 |
| `__symmetric_mod(c, pa)` | `polynomial_factorize_wang.hh` | 对称约化到 Z |

**需新增 1 个**（约 8 行）：

```cpp
// upolynomial.hh —— 补全 poly_convert 体系的缺口
// 现有 polynomial_convert.hh L101 有 poly_convert(polynomial_<Zp,comp> → polynomial_<T2,comp>)
// upolynomial_ 对应重载缺失，__upoly_Zp_to_ZZ 定义在 polynomial_factorize_univar.hh
// 无法在 polynomial_factorize_wang.hh 中访问，需补全
inline void poly_convert(const upolynomial_<Zp>& p_in, upolynomial_<ZZ>& p_out)
{
    p_out.clear();
    p_out.reserve(p_in.size());
    for (auto& term : p_in)
        p_out.push_back({term.first, ZZ(static_cast<int64_t>(term.second.number()))});
}
```

> 位置：`upolynomial.hh`，紧接 `polynomial_mod` for `upolynomial_`（L115 之后）。
> 同时补对称方向：`poly_convert(const upolynomial_<ZZ>&, upolynomial_<Zp>&, uint32_t p)` 可复用 `polynomial_mod`，但 M0 仅需 Zp→ZZ 方向。

---

## 4. 正确性论证

**命题**：改后基本情形输出与改前相同（mod 当前 pa）。

**证明**：
1. Fp[x] Bézout 链满足 `Σ si_p·V̂i_p ≡ 1 (mod p)`（由扩展 GCD 保证，两两互素保证 gcd=1）
2. `δi_p := (si_p · h_p) rem vi_p` 满足 `Σ δi_p · V̂i_p ≡ h_p (mod p)`（Bézout 恒等式 + 同余）
3. 真实整数解 `δi` 满足 `δi ≡ δi_p (mod p)`（由解的唯一性）
4. 设 `pa = p^a > 2·||δi||_∞`，则 `symmetric_mod(δi_p, pa) = δi`（唯一表示）

**B 的选取**：当前使用 `B = ∏(||vi||_1 + 1)`，这是 Diophantine 解的一个安全上界（参考 GCL Lemma 6.17）。

---

## 5. 性能预期

**bivar-70（r=70 个线性因子）**：

| 步骤 | 改前 | 改后 | 说明 |
|------|------|------|------|
| Bézout 链构建 | ~0.1s（Z[x] GCD 累积） | < 1ms（Fp[x] GCD） | 系数 1 bit vs 0 bit |
| 每次 Diophantine | ~1 亿大整数 ops | ~100 万小整数 ops | bezout_s[0] 6000-bit→64-bit |
| 总计（2 次 Diophantine） | ~65s | **< 0.5s** | 约 100-200x 加速 |

---

## 6. 实施步骤

```
Step 1: 实现 __bezout_chain_modular
        → 单元测试：验证 Σ si_p·V̂i_p ≡ 1 (mod p) 对各 r 值
          (约 1 小时)

Step 3: 修改 __multivar_diophantine 基本情形
        → 运行 test_multivar_hensel（全 10 个测试）
          (约 2 小时)

Step 4: 修改 __multivar_hensel_lift + __hensel_lift_one_var 签名
        → 运行 test_factorize_multivar（全 ~50 个测试）
          (约 1 小时)

Step 5: 性能验证
        → make stress（bivar-70 目标 < 1s）
        → make bench-save（存档基线）
```

**总估时**：1 天

---

## 7. 测试矩阵

| 测试 | 命令 | 验收标准 |
|------|------|---------|
| Fp[x] GCD 单元测试 | `_build/debug/bin/test_multivar_hensel` | 10/10 PASS |
| 端到端因式分解 | `_build/debug/bin/test_factorize_multivar` | ~50/50 PASS |
| crosscheck | `make crosscheck` | PASS |
| 压力测试 | `make stress` | bivar-70 < 1s |
| 性能基准 | `make bench-save` | 存档 bivar-70 时间 |

---

## 8. 风险

| 风险 | 概率 | 缓解 |
|------|------|------|
| `polynomial_GCD<Zp>` 接口不兼容 | 极低 | 已验证 L719 接口，直接使用 `polynomial_GCD(F, G, s, t)` |
| p 被所有 lc(vi) 整除（极小概率） | 极低 | 循环选下一个素数（已有模式）|
| pa 界不够紧 | 低 | 放大 B 公式（当前已保守） |
| D5/LC 回归测试失败 | 低 | 首先运行 test_multivar_hensel 的 lc 测试 |
