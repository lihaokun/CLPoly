# 修正方案：Wang LC 分配幂次拆分

> 对应 workflow.md §5.1 修正方案文档

---

## 第一部分：复现与定位

### 最小复现用例

```cpp
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

int main() {
    variable x("x"), y("y"), z("z");
    PolyZZ f;
    poly_convert(
        polynomial_ZZ(x*y + z*z + ZZ(1)) *
        polynomial_ZZ(ZZ(2)*x*x - z*z + ZZ(2)*y) *
        polynomial_ZZ(-x*y - ZZ(2)*y*y - ZZ(1)), f);
    auto fac = factorize(f);
    // 预期: fac.factors.size() >= 3
    // 实际: fac.factors.size() == 1 (返回原多项式未分解)
}
```

### 出错路径

```
factorize(f)
  → __factor_multivar(f)
    → squarefreefactorize(f) → 1 个无平方分量
    → __wang_core(g)
      → 遍历主变量 x, y, z
        → Phase 1: 扫描 50 个求值点，收集候选
        → Phase 2: 对每个候选调用 __wang_leading_coeff()
          → 第 319-341 行：匹配循环 → 全部返回 success=false
      → 所有主变量均失败 → 返回 {{g, 1}} (未分解)
```

关键失败点：`polynomial_factorize.hh:319-341`，`__wang_leading_coeff` 的 LC 因子匹配循环。

### 预期行为 vs 实际行为

| | 预期 | 实际 |
|---|------|------|
| LC 匹配 | `y²` 拆分为 `y¹·y¹`，分配给 u₁ 和 u₃ | 试图把 `y²` 整体分配给一个因子，歧义 → FAIL |
| `__wang_core` | 返回 3 个因子 | 返回 1 个因子（原多项式） |

---

## 第二部分：根因分析

### 症状

`__wang_leading_coeff` 对所有求值点、所有主变量返回 `success=false`。

### 根因

第 319-341 行的匹配循环将 `(lⱼ, eⱼ)` 作为**原子单元**，试图整体分配给唯一一个单变量因子：

```cpp
// 第 320-340 行 (当前代码)
size_t best_idx = r;
int best_count = 0;
for (size_t i = 0; i < r; ++i)
{
    int cnt = 0;
    ZZ tmp = w[i];
    while (tmp % lj_val == ZZ(0))
    { tmp /= lj_val; ++cnt; }

    if (cnt > best_count)
    { best_count = cnt; best_idx = i; }
    else if (cnt == best_count && cnt > 0)
    { best_idx = r; }   // ← 歧义即失败
}
if (best_idx >= r) return result;  // FAIL
// ...
sigma[best_idx] = sigma[best_idx] * lj_power;  // lⱼ^eⱼ 整体给 best_idx
```

问题：当同一个不可约因子 `lⱼ` 的 `eⱼ` 次幂需要**拆分**给多个因子时（如 `y²` 拆为 `y·y`），不存在唯一赢家，代码宣告歧义并失败。

这不是 workaround 能解决的——原子分配在数学上就不对。

---

## 第三部分：参考实现对照

### 对照对象：SymPy `dmp_zz_wang_lead_coeffs`

SymPy 源码 (`sympy/polys/factortools.py`) 的匹配循环：

```python
for i, h in enumerate(H):          # 外层: 遍历每个单变量因子
    c = dmp_one(v, K)               # σᵢ = 1
    d = dup_LC(h, K) * cs           # wᵢ = lc(uᵢ) × content

    for j in reversed(range(len(E))):   # 内层: 遍历每个不可约 LC 因子
        k, e, (t, _) = 0, E[j], T[j]
        while not (d % e):          # 从 wᵢ 中提取 lⱼ(α) 的最大幂次
            d, k = d // e, k + 1
        if k != 0:
            c = dmp_mul(c, dmp_pow(t, k, v, K), v, K)  # σᵢ *= lⱼ^kᵢ
```

### 关键差异

**循环结构不同**：

| | CLPoly（当前） | SymPy / FLINT / Wang 原文 |
|---|------|------|
| 外层 | 遍历 LC 因子 `(lⱼ, eⱼ)` | 遍历单变量因子 `uᵢ` |
| 内层 | 遍历 `uᵢ`，找唯一赢家 | 遍历 LC 因子，各取所需 |
| 分配 | `lⱼ^eⱼ` 整体给一个 `uᵢ` | 每个 `uᵢ` 取 `lⱼ^kᵢ`，`Σkᵢ = eⱼ` |
| 歧义处理 | 宣告失败 | 不需要唯一性，无歧义问题 |

**用复现用例走一遍差异**（求值点 y=3, z=1，content=-1）：

```
w = [lc(u₀)×cs, lc(u₁)×cs, lc(u₂)×cs] = [-3, -2, -3]
lⱼ = y, eⱼ = 2, lⱼ(α) = 3
```

CLPoly（当前）：
```
对 (y, 2)：
  i=0: cnt=1, i=2: cnt=1 → 歧义 → FAIL ✗
```

SymPy 方式：
```
对 u₀ (w₀=-3): 提取 3 的幂次 → k₀=1, w₀ → -1.  σ₀ *= y¹
对 u₁ (w₁=-2): 提取 3 的幂次 → k₁=0.             σ₁ 不变
对 u₂ (w₂=-3): 提取 3 的幂次 → k₂=1, w₂ → -1.  σ₂ *= y¹
验证: Σkᵢ = 1+0+1 = 2 = eⱼ ✓ SUCCESS!
```

### 对照对象：FLINT `fmpz_mpoly_factor_lcc_wang`

FLINT 的逻辑与 SymPy 一致（外层遍历因子，内层遍历 LC）：

```c
for (j = 0; j < r; j++) {               // 外层: 遍历每个单变量因子
    fmpz_mul(R, lc(Auf[j]), Auc);        // R = lc(uⱼ) × content
    for (i = lcAfac->num - 1; i >= 0; i--) {  // 内层: 遍历每个 LC 因子
        k = fmpz_remove(R, R, Q);        // 提取 lⱼ(α) 的最大幂次
        lc_divs[j] *= lcAfac[i]^k;       // σⱼ *= lⱼ^k
    }
}
```

### content 污染问题

当 `content(f₀)` 与 `lⱼ(α)` 共享素因子时（如 content=-2, lⱼ(α)=-2），逐因子提取会**过度提取**（总幂次 > eⱼ）。

SymPy 的处理：`dmp_zz_wang_non_divisors` 检查 `cs × γ` 与各 `lⱼ(α)` 的互素性，拒绝有冲突的求值点。

FLINT 的处理：有 GCD 预处理步骤剥离共享因子；Wang LC 失败则切换 Kaltofen 方法。

**我们的处理**：已有换点重试机制（50 个点），只需在提取后验证 `Σkᵢ = eⱼ`，不等则拒绝该点。50 个点中总有 content 与 lⱼ(α) 互素的点。

---

## 第四部分：修正方案

### 修改位置

`clpoly/polynomial_factorize.hh`，`__wang_leading_coeff` 函数，第 309-355 行（LC 因子匹配循环）。

### 修改内容

将"外层遍历 LC 因子、内层找唯一赢家"改为"外层遍历 LC 因子、内层遍历所有 uᵢ 各取所需"。

**修改前**（第 309-355 行）：

```cpp
for (auto& [lj, ej] : lc_factors)
{
    auto lj_eval = assign(lj, eval_point);
    ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
    if (lj_val == 0) return result;

    // 计算 zⱼ = lj_val^eⱼ（实际未使用，仅用于原子分配逻辑）
    ZZ zj(1);
    for (uint64_t e = 0; e < ej; ++e) zj *= lj_val;

    // 找唯一赢家
    size_t best_idx = r;
    int best_count = 0;
    for (size_t i = 0; i < r; ++i)
    {
        int cnt = 0;
        ZZ tmp = w[i];
        while (tmp % lj_val == ZZ(0)) { tmp /= lj_val; ++cnt; }
        if (cnt > best_count) { best_count = cnt; best_idx = i; }
        else if (cnt == best_count && cnt > 0) { best_idx = r; }
    }
    if (best_idx >= r) return result;  // FAIL

    // lⱼ^eⱼ 整体给 best_idx
    Poly lj_power = lj;
    for (uint64_t e = 1; e < ej; ++e) { lj_power = lj_power * lj; ... }
    sigma[best_idx] = sigma[best_idx] * lj_power;
    for (int c = 0; c < best_count; ++c) w[best_idx] /= lj_val;
}
```

**修改后**：

```cpp
for (auto& [lj, ej] : lc_factors)
{
    auto lj_eval = assign(lj, eval_point);
    ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
    if (lj_val == 0 || abs(lj_val) == ZZ(1)) return result;
    // |lⱼ(α)| ≤ 1 时无法从数值区分分配，拒绝该点 (FLINT 同理: fmpz_cmp_ui(Q,2)<0)

    // 对每个 uᵢ 提取 lⱼ(α) 的最大幂次
    std::vector<int> k(r, 0);
    int total_k = 0;
    for (size_t i = 0; i < r; ++i)
    {
        while (w[i] % lj_val == ZZ(0))
        {
            w[i] /= lj_val;
            ++k[i];
        }
        total_k += k[i];
    }

    // 验证总幂次守恒（不等说明 content 污染，拒绝该求值点）
    if (total_k != (int)ej)
        return result;

    // σᵢ *= lⱼ^kᵢ
    for (size_t i = 0; i < r; ++i)
    {
        if (k[i] == 0) continue;
        Poly lj_power = lj;
        for (int e = 1; e < k[i]; ++e)
        {
            lj_power = lj_power * lj;
            lj_power.normalization();
        }
        sigma[i] = sigma[i] * lj_power;
        sigma[i].normalization();
    }
}
```

### 修正后复现用例的执行过程

求值点 y=3, z=1, content=-1:

```
w = [-3, -2, -3], lⱼ = y, eⱼ = 2, lⱼ(α) = 3

i=0: w[0]=-3, -3%3=0 → w[0]=-1, k[0]=1. -1%3≠0 → 停止.
i=1: w[1]=-2, -2%3≠0 → k[1]=0.
i=2: w[2]=-3, -3%3=0 → w[2]=-1, k[2]=1. -1%3≠0 → 停止.

total_k = 1+0+1 = 2 = eⱼ ✓

σ[0] *= y¹ → σ[0] = y
σ[1] 不变  → σ[1] = 1
σ[2] *= y¹ → σ[2] = y

γ·σ[0]·σ[1]·σ[2] = (-2)·y·1·y = -2y² = L ✓
```

求值点 y=-2, z=-1, content=-2（content 污染的坏点）:

```
w = [-2, -4, -4], lⱼ(α) = -2

i=0: k[0]=1, w[0]=-1
i=1: k[1]=2, w[1]=-1
i=2: k[2]=2, w[2]=-1

total_k = 1+2+2 = 5 ≠ eⱼ = 2 → FAIL, 换点 ✓（正确拒绝）
```

### 正确性保证

`Σkᵢ = eⱼ` 验证是必要条件，不是充分条件。当 `content(f₀)` 与 `lⱼ(α)` 共享素因子时，幂次可能被错误分配（如本应给 u₀ 的幂次被 u₁ "偷走"）。这种情况由后续的 Hensel 提升 + 试除验证兜底——错误的 LC 分配会导致提升结果无法整除原多项式，从而触发换点重试。条件 (d)（不同 `lⱼ(α)` 两两互素）保证不同不可约因子之间不会交叉污染。

### 不改动的部分

- `w[i]` 初始化仍乘以 `uni_content`（与 SymPy/FLINT 一致）
- P2 恒等式验证（第 292-307 行）保留
- `__wang_core` 的主变量轮换、Phase 1/2 结构、因子重组均保留（这些是独立的有价值改进）
- `__select_eval_point` 不变
