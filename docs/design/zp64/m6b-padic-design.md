# M6b 细化设计：p-adic 提升 + Taylor 循环对齐

> 阶段：细化（workflow.md §2.3）
> 前置：M6a 已完成（Zp 64-bit，mtshl_p = 63-bit 素数）
> 目标文件：`clpoly/polynomial_factorize_wang.hh`
> 论文：CASC 2018 Algorithm 5 + Algorithm 4 Line 3

---

## 论文 vs 当前实现对比

### Algorithm 5（p-adic 提升）

```
论文流程：
1. 选 p 为 63-bit 机器素数
2. 在 Zp 中运行 MTSHL-d（Algorithm 4）→ F_zp[i] mod p
3. F_zz[i] = symmetric_mod(F_zp[i], p)
4. 计算 B = Mignotte 界, M = p
5. while M ≤ 2B:
     error = f_scaled - ∏ F_zz[i]          ← Z 精确计算
     if error = 0: break
     dk = error / M                          ← Z 精确整除（p-adic 保证）
     dk_zp = dk mod p                        ← 投影到 Zp
     SparseInt 解: Σ σi · bi = dk_zp         ← bi = ∏_{j≠i} F_base[j]
     F_zz[i] += symmetric_mod(σi, p) · M
     M ← M · p
6. return F_zz[i]

当前实现：
1. ✓ mtshl_p = 63-bit 素数（M6a 完成）
2. ✓ Zp 中运行 MTSHL-d（__mtshl_lift 阶段 A+B）
3. ✓ symmetric_mod（阶段 C）
4-5. ✗ 无 p-adic 迭代——单次 symmetric_mod 后直接返回
6. ✗ 若系数 |c| > p/2，对称约化错误 → trial division 拒绝 → 重试求值点
```

### Algorithm 4 Line 3（Taylor 循环终止）

```
论文：for k = 1, 2, ... while error ≠ 0 AND ∑deg(fi,xj) < deg(aj,xj) do
当前：for (int k = 1; k <= D_j; ++k)（无 error=0 提前退出）
```

---

## MDP 求解器复用分析

p-adic 的 MDP 与 MTSHL step_j 中的 MDP 结构相同：

```
给定: dk_zp ∈ Zp[x1,...,xn], F_base[i] ∈ Zp[x1]（单变量基础因子）
求解: Σ σi · (∏_{j≠i} F_base[j]) = dk_zp
```

`__mtshl_sparse_int` 已能处理此问题：
- 输入 `F_base`（单变量 Zp[x1]）+ `ck`（多变量 Zp[x1,...,xn]）+ `forms`（骨架）
- 通过 θ-array 求值将多变量 ck 归约到单变量
- 逐点调用 `__mtshl_zp_univar_mdp` 求解
- Vandermonde 恢复多变量系数

p-adic 调用时：
- `F_base` = scaled_factors mod p（与 MTSHL step j=2 的基础因子相同）
- `ck` = dk_zp（全变量误差 mod p）
- `forms` = F[i] 的当前骨架（Theorem 1：逐步缩小）
- `aux_vars` = [x2,...,xn]（全部辅助变量）

**结论：可以直接复用 `__mtshl_sparse_int` + 回退链，无需新增 MDP 求解器。**

---

## Step 1：Taylor 循环提前退出

**文件**：`wang.hh` `__mtshl_step_j` 函数

在 Taylor 循环末尾（误差重算之后）加：

```cpp
    // 全量重算误差
    error = aj - product_F();
    error.normalization();
+
+   // 论文 CASC 2018 Alg.4 Line 3: error = 0 时提前退出
+   if (error.empty()) break;
}
```

---

## Step 2：p-adic 提升

### 2a. 多变量 Mignotte 界

需要新增 `__mtshl_coeff_bound`：估算 f_scaled 因子系数的上界。

对多变量多项式 f ∈ Z[x1,...,xn]，因子系数的界：
- 简化使用：B = ‖f_scaled‖∞（最大系数绝对值）× 2^deg(f_scaled, x1)
- 论文精确使用 Landau-Mignotte 界，但对 l = ⌈log_p(2B)⌉ 的影响很小

```cpp
template<class var_order>
ZZ __mtshl_coeff_bound(const polynomial_<ZZ, lex_<var_order>>& f)
{
    ZZ max_coeff(0);
    for (const auto& term : f)
    {
        ZZ a = abs(term.second);
        if (a > max_coeff) max_coeff = a;
    }
    // B = 2^(deg+1) · max_coeff（保守上界）
    int64_t d = get_first_deg(f);  // x1 的次数
    ZZ bound = max_coeff;
    for (int64_t i = 0; i <= d; i++)
        bound *= ZZ(2);
    return bound;
}
```

### 2b. `__mtshl_lift` 阶段 C 改造

当前阶段 C（L1070-1075）：
```cpp
// 阶段 C: 系数恢复（对称约化）
std::vector<Poly> result(r, Poly(comp_ptr));
for (int i = 0; i < r; i++)
    result[i] = __symmetric_mod_poly(F[i], p);

return result;
```

改为：

```cpp
// 阶段 C: 系数恢复 — 初始对称约化 + p-adic 提升

// C.1: 初始对称约化
std::vector<Poly> result(r, Poly(comp_ptr));
for (int i = 0; i < r; i++)
    result[i] = __symmetric_mod_poly(F[i], p);

// C.2: 检查是否需要 p-adic 提升
ZZ B = __mtshl_coeff_bound(f_scaled);
ZZ M((int64_t)p);   // 当前已约化到 mod p

if (M > ZZ(2) * B)
    return result;   // 单次约化足够

// C.3: 保存基础因子和骨架（p-adic MDP 复用）
// F_base_zp = scaled_factors mod p → Zp[x1]（单变量，与 step j=2 相同）
std::vector<PolyZp> F_base_zp(r, PolyZp(comp_ptr));
for (int i = 0; i < r; i++)
{
    for (const auto& term : scaled_factors[i])
    {
        Zp coeff(term.second, p);
        if (coeff.number() == 0) continue;
        basic_monomial<lex_<var_order>> m(comp_ptr);
        if (term.first.deg() > 0)
            m.push_back({main_var, term.first.deg()});
        m.normalization();
        F_base_zp[i].push_back({m, coeff});
    }
    F_base_zp[i].normalization();
}

// forms = 当前因子骨架（从 MTSHL-d 结果 F[i] 提取）
std::vector<std::vector<basic_monomial<lex_<var_order>>>> forms(r);
for (int i = 0; i < r; i++)
    for (const auto& term : F[i])
        forms[i].push_back(term.first);

// C.4: p-adic 提升循环（CASC 2018 Algorithm 5）
constexpr int l_max = 5;  // 63-bit 素数 l=5 → 可恢复 ~10^95

for (int l = 1; l < l_max && M <= ZZ(2) * B; l++)
{
    // C.4.1: 计算 Z 上精确误差
    Poly product = result[0];
    for (int i = 1; i < r; i++)
    {
        product = product * result[i];
        product.normalization();
    }
    Poly error_zz = f_scaled - product;
    error_zz.normalization();
    if (error_zz.empty()) break;  // 已精确恢复

    // C.4.2: dk = error / M, dk_zp = dk mod p
    PolyZp dk_zp(comp_ptr);
    for (const auto& term : error_zz)
    {
        ZZ q;
        ZZ::fdiv_q(q, term.second, M);  // 精确整除
        Zp coeff_zp(q, p);               // ZZ → Zp（用 fdiv_ui）
        if (coeff_zp.number() != 0)
            dk_zp.push_back({term.first, coeff_zp});
    }
    dk_zp.normalization();
    if (dk_zp.empty()) break;

    // C.4.3: MDP 求解（复用现有求解链）
    std::vector<PolyZp> sigma_k;

    bool success = __mtshl_sparse_int(
        F_base_zp, dk_zp, forms, main_var, aux_var_list, p, sigma_k);
    if (!success)
        success = __mtshl_sparse_int(
            F_base_zp, dk_zp, forms, main_var, aux_var_list, p, sigma_k);
    if (!success)
    {
        if (aux_var_list.size() == 1)
            success = __mtshl_multi_bdp(
                F_base_zp, dk_zp, main_var, aux_var_list[0],
                ideal_alphas_zp[0], sigma_k);
        else
            success = __mtshl_wmds(
                F_base_zp, dk_zp, main_var, aux_var_list,
                ideal_alphas_zp, sigma_k);
    }
    if (!success) break;  // MDP 失败 → 放弃 p-adic，走 trial division

    // C.4.4: 修正 result[i] += symmetric_mod(σi, p) · M
    for (int i = 0; i < r; i++)
    {
        if (i >= (int)sigma_k.size() || sigma_k[i].empty()) continue;
        Poly correction = __symmetric_mod_poly(sigma_k[i], p);
        for (auto& term : correction)
            term.second *= M;
        result[i] = result[i] + correction;
        result[i].normalization();

        // 更新骨架（Theorem 1 强支撑：Supp(σk) ⊆ Supp(σk-1)）
        forms[i].clear();
        for (const auto& term : sigma_k[i])
            forms[i].push_back(term.first);
    }

    M *= ZZ((int64_t)p);
}

return result;
```

### 2c. `__wang_core` 注释更新

```cpp
// 旧（L2161-2164）：
// MTSHL prime: 使用最大 31-bit 素数
// 不做 Mignotte 界检查——依靠 trial division 验证因子正确性。
// 若 p 不足以恢复真实系数，对称约化给出错误因子，
// trial division 会拒绝并让 __wang_core 重试其他求值点。

// 新：
// MTSHL prime: 63-bit 机器素数（CASC 2018 Algorithm 5）
// 系数恢复: symmetric_mod + p-adic 提升（|c| > p/2 时自动迭代）
// 极端情况由 trial division 兜底。
```

---

## Step 3：文档同步

更新 `post-m5-optimizations.md`：
- M6a 状态 → 已完成
- M6b 状态 → 已完成
- 差异 1 状态 → 已完成

---

## 改动汇总

| Step | 文件 | 改动 |
|------|------|------|
| 1 | `wang.hh` `__mtshl_step_j` | Taylor 循环加 `if (error.empty()) break;` |
| 2a | `wang.hh` 新函数 | `__mtshl_coeff_bound` 多变量系数界 |
| 2b | `wang.hh` `__mtshl_lift` 阶段 C | p-adic 提升循环（~50 行新代码） |
| 2c | `wang.hh` `__wang_core` | 注释更新 |
| 3 | `post-m5-optimizations.md` | M6 状态同步 |

验收：
1. `test_number` 通过（334 tests）
2. `test_factorize` + `test_factorize_zp` 通过
3. `make bench-clpoly` 性能无退化
4. 可选：构造大系数测试用例验证 p-adic 路径
