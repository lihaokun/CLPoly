# 修正方案：LC 分配 gamma 互素性问题

> **状态：已完成。** 最终采用 GCD 匹配替代幂次提取 + skip 参数 + 交错主变量轮换。
> 本文档保留为历史记录，权威描述见 [multivariate-factorization-design.md](../design/multivariate-factorization-design.md)。
>
> 对应 workflow.md §5.1 修正方案文档

---

## 第一部分：复现与定位

### 最小复现用例

```cpp
variable x("x"), y("y"), z("z");

// Fail 1: L(x 主变量) = 6z
auto f1 = polynomial_ZZ(ZZ(3)*x*z + ZZ(2)*z*z - z + ZZ(2));
auto f2 = polynomial_ZZ(ZZ(2)*x*x - ZZ(2)*x*y - ZZ(3)*y + ZZ(1));
auto f = make_lex(f1 * f2);
// factorize(f) 返回 1 因子（未分解）

// Fail 2: L(x 主变量) = 6y
auto g1 = polynomial_ZZ(ZZ(3)*x*y - y*y - ZZ(3)*y*z + ZZ(2));
auto g2 = polynomial_ZZ(ZZ(2)*x*x + ZZ(2)*z*z + x + ZZ(3));
auto g = make_lex(g1 * g2);
// factorize(g) 返回 1 因子（未分解）
```

### 出错路径（以 Fail 2, x 主变量为例）

```
__wang_core(f), x 主变量:
  L = lc(f, x) = 6y
  factorize(L): gamma = 6, lc_factors = [(y, 1)]

  Phase 1: 扫描 50 个求值点
    每个点: lj_val = y(α), 单变量分解得 2 因子, cs = content(f₀)
    → 全部送入 Phase 2

  Phase 2: 对每个候选调用 __wang_leading_coeff
    对每个点:
      lj_val = y(α) ∈ {±2, ±3, ±4, ...}
      提取 lj_val 的幂次 → total_k ≠ 1 = ej → 拒绝该点
      → 50 个候选全部 lc_ok=0

  y 主变量: deg(f,y) = 2, f2 不含 y → 特化后总是不可约 → continue
  z 主变量: L = -6y, 同样 50 个点全部 lc_ok=0

  所有主变量失败 → 返回 {{f, 1}}
```

---

## 第二部分：根因分析

### 症状

`total_k ≠ ej` 验证对所有 50 个求值点都失败。

### 根因：gamma 与 lj_val 共享素因子

P2 恒等式（数据一致性）：

```
gamma · ∏ lj_val^ej = cs · ∏ lc(u_i)
```

当前代码 `w[i] = lc(u_i) * cs`，从 w[i] 中提取 lj_val 的最大幂次。

**问题**：gamma = 6 = 2·3，lj = y。对于 y(α) = -2：

```
w[0] = lc(u_0) * cs,  w[1] = lc(u_1) * cs
prod(w) = cs² · prod(lc(u_i)) = cs · gamma · lj_val^ej = cs · 6 · (-2) = -12cs

提取 lj_val = -2:
  prod(w) 中 (-2) 的幂次 = v₂(cs) + v₂(6) + ej·v₂(-2)
                          = v₂(cs) + 1 + 1 = v₂(cs) + 2

个别提取: total_k = v₂(cs) + 2
期望: ej = 1
差值: v₂(cs) + 1 ≥ 1 → total_k > ej → 拒绝
```

无论 cs 是什么，只要 `gcd(gamma, lj_val) > 1`，gamma 中的共享素因子就会"泄漏"到 w 值中，导致 total_k 超过 ej。

### 为什么 50 个点不够

L = 6y, gamma = 6 = 2·3。需要 y(α) 满足：
- `|y(α)| > 1`（避免无穷提取）
- `gcd(6, y(α)) = 1`（避免 gamma 污染）

即 y(α) 不能被 2 或 3 整除，且绝对值 > 1。最小的满足条件的值：**±5, ±7, ±11, ...**

点选择器按范围 2, 3, 4, 5, 6, ... 扫描。在范围 [-4, 4] 中，没有满足条件的 y 值。范围 [-5, 5] 才出现 y=±5。50 个 skip 值可能刚好不够到达范围 5 的有效点。

### SymPy 如何处理

SymPy 在 `dmp_zz_wang` 中调用 `dmp_zz_wang_non_divisors`，核心算法：

```python
def dmp_zz_wang_non_divisors(E, cs, ct, K):
    """E = [(lj_val, ej), ...], cs = content(f₀), ct = gamma"""
    result = []
    for (e, _) in E:
        q = cs * ct           # 初始化 q = cs * gamma
        # 从 q 中剥离所有已处理的 lj_val 的因子
        for (r, _) in result:
            while K.rem(q, r) == 0:
                q = K.quo(q, r)
        # 从 q 中剥离当前 lj_val 的因子
        while K.rem(q, e) == 0:
            q = K.quo(q, e)
        if K.abs(e) <= K.abs(q):
            result.append((e, _))
        else:
            return None  # |lj_val| > |剩余 q| → 拒绝该求值点
    return result
```

**实质**：对每个 LC 因子的求值 `lj_val`，从 `cs·gamma` 中剥离与 `lj_val` 共享的素因子。如果剥离后 `|lj_val| > |q|`，说明 `lj_val` 的因子不够独立（被 gamma 或 cs "吞没"），拒绝该点。

**求值点增长**：SymPy 用 `mod` 参数控制范围（从 1 开始，每次失败 `mod += 2`），范围扩大直到找到好点。

### FLINT 如何处理

FLINT 在 `fmpz_mpoly_factor_lcc_wang` 中实现**相同的 non-divisors 算法**：

```c
// fmpz_mpoly_factor_lcc_wang (factor/fmpz_mpoly_factor_lcc_wang.c)
for (j = 0; j < r; j++)
{
    fmpz_set(Q, lcAevals + j);       // Q = lc(u_j)
    fmpz_abs(Q, Q);
    // 剥离与 gamma·cs 共享的因子
    fmpz_mul(t, Acontent, ct);        // t = cs * gamma
    if (fmpz_cmp_ui(Q, 2) < 0)       // |Q| < 2 → 拒绝
        goto next_alpha;
    for (k = 0; k < j; k++)
    {
        // 剥离与已处理 lc(u_k) 共享的因子
        fmpz_set(Q2, lcAevals + k);
        while (fmpz_divisible(Q, Q2))
            fmpz_divexact(Q, Q, Q2);
    }
    // 检查剩余值
    if (fmpz_cmp_ui(Q, 2) < 0)
        goto next_alpha;
}
```

**求值点增长**：`alpha_modulus` 从 1 开始，每次失败 +1，遍历所有组合。此外 FLINT 有 Kaltofen/Zassenhaus 兜底算法。

### CLPoly 缺失的检查

CLPoly 的 `__select_eval_point` 检查了 LC 因子求值之间的两两互素性（条件 d），但**未检查与 gamma 的互素性**。当 `gcd(gamma, lj_val) > 1` 时，gamma 的素因子泄漏到 w 值中，导致提取过度。

SymPy 和 FLINT 的 non-divisors 算法本质是：**确保每个 `lj_val` 在剥离 `gamma·cs` 的共享因子后仍有足够大的独立部分**，这恰好解决了 gamma 污染问题。

---

## 原型实验结果

使用 `/tmp/proto_gamma.cc` 对两个失败用例扫描 300 个 skip 值：

| 用例 | 主变量 | 首次成功 skip | 求值点 | 原因 |
|------|--------|-------------|--------|------|
| Fail1 | x | 54 | y=-5, z=-5 | `lj_val=z(α)=-5`, `gcd(6,5)=1` |
| Fail1 | z | 无 (300点) | — | LC 因子 `2x²-2xy-3y+1` 是多变量，情况更复杂 |
| Fail2 | x | 54 | y=-5, z=-5 | `lj_val=y(α)=-5`, `gcd(6,5)=1` |
| Fail2 | z | 54 | x=-5, y=-5 | `lj_val=-y(α)=5`, `gcd(6,5)=1` |

**关键发现**：

1. 当前 `SCAN_SIZE=50`，首次好点在 skip=54 → **差 4 个点**，移除硬上限后自然覆盖
2. 好点集中在 `|lj_val|=5`（第一个与 6=2·3 互素且 >1 的整数）
3. Fail1 z-main 失败是因为 LC 因子本身是多变量多项式（`2x²-2xy-3y+1`），non-divisors 检查对此不适用，但 x-main 成功就够了
4. 当前 factorize 返回 1 因子，确认 bug 存在

---

## 第三部分：修正方案

### 修改位置

`clpoly/polynomial_factorize.hh`，`__wang_leading_coeff` 函数。

### 修改内容

两级修改：non-divisors 预过滤 + 增大扫描范围。

#### 修改 1：non-divisors 预过滤（`__wang_leading_coeff`）

在提取循环开始前，增加 SymPy/FLINT 的 non-divisors 检查。对每个 LC 因子求值 `lj_val`，从 `cs·gamma` 中剥离共享因子，若剩余值太小则拒绝。

```cpp
// 在提取循环之前增加:

// non-divisors 预过滤 (SymPy dmp_zz_wang_non_divisors)
// 确保每个 lj_val 在剥离 gamma·cs 共享因子后仍有足够独立性
{
    std::vector<ZZ> checked_vals;
    for (auto& [lj, ej] : lc_factors)
    {
        auto lj_eval = assign(lj, eval_point);
        ZZ lj_val = is_number(lj_eval) ? lj_eval.front().second : ZZ(0);
        if (lj_val == 0 || abs(lj_val) == ZZ(1)) return result;

        ZZ q = uni_content * gamma;
        // 剥离已处理的 lj_val 的因子
        for (auto& prev : checked_vals)
            while (q % prev == ZZ(0))
                q /= prev;
        // 剥离当前 lj_val 的因子
        while (q % lj_val == ZZ(0))
            q /= lj_val;
        // 若 |lj_val| > |q|，说明独立部分不足 → 拒绝
        if (abs(lj_val) > abs(q))
            return result;

        checked_vals.push_back(lj_val);
    }
}
```

**这是 early reject 优化**——被拒绝的点本来也会被 `total_k != ej` 拒绝，新检查只是提前退出。

#### 修改 2：交错轮换主变量 + 移除单变量硬上限（`__wang_core`）

当前 `SCAN_SIZE=50` 是硬上限，用完即放弃该主变量。两个问题：
1. 50 个点可能不够到达好点（如 gamma=6 需 skip≥54）
2. 在"坏"主变量上浪费全部配额后才换，效率低

改为**交错轮换**：每个主变量分配 `BATCH_SIZE=200` 个求值点，用完后换下一个主变量，一轮全部用完后扩展配额继续（最多 `MAX_ROUNDS=3` 轮）。

```cpp
// 修改前:
const int SCAN_SIZE = 50;
for (var_idx = 0; var_idx < n_vars; ++var_idx)  // 顺序尝试每个主变量
    for (skip = 0; skip < SCAN_SIZE; ++skip)     // 50 点硬上限
        ...

// 修改后: 交错轮换
const int BATCH_SIZE = 200;
const int MAX_ROUNDS = 3;
// 每个主变量维护自己的 skip 位置和不可约状态
for (round = 0; round < MAX_ROUNDS; ++round)
    for (vi = 0; vi < n_vars; ++vi)           // 轮换每个主变量
        if (!var_dead[vi])
            for (skip = var_skip[vi]; skip < var_skip[vi] + BATCH_SIZE; ++skip)
                // 立即尝试完整流程
```

**示例**：三变量 (x, y, z)，x-main 在 skip=54 才有好点：

```
Round 0: x 试 0-199 → 在 skip=54 成功 → return
```

若 x-main 200 点全部失败：

```
Round 0: x 试 0-199 (失败), y 试 0-199 → 可能在 y-main 成功
Round 1: x 试 200-399, y 试 200-399, ...
```

**终止保证**：
- 不可约检测（3 连续单因子 → 标记 `var_dead`）
- 所有主变量 dead → 退出
- `MAX_ROUNDS` 轮次上限（共 `n_vars × BATCH_SIZE × MAX_ROUNDS` 点）
- 可约多项式理论保证好点无穷多

### 修正后复现用例的执行过程（原型实测）

```
Fail 2: f = (3xy-y²-3yz+2)(2x²+2z²+x+3), x 主变量

L = 6y, gamma = 6, lj = y, ej = 1

skip=0..53: y ∈ {±2, ±3, ±4}
  → gcd(6, y(α)) > 1 → non-divisors 快速拒绝（不再被硬上限截断）

skip=54: y=-5, z=-5
  → lj_val = -5, gcd(6, 5) = 1 → non-divisors 通过
  → 提取: total_k = 1 = ej ✓
  → Hensel lift 得 2 因子 → 试除验证 2/2 ✓ → SUCCESS
```

### 正确性保证

1. **non-divisors 不改变语义**：被拒绝的点本来也会被 `total_k != ej` 拒绝，只是提前退出节省计算。仅在 `|γ| > 1` 时启用，避免 γ=1 时退化（所有点被错误拒绝）。
2. **交错轮换不影响正确性**：每个主变量的 skip 独立递增，不会遗漏任何求值点。
3. **不可约多项式不会死循环**：连续 3 个单因子 → 标记该主变量 dead；所有主变量 dead 或 `MAX_ROUNDS` 轮结束 → 返回 `{{f, 1}}`。
4. **Hensel + 试除验证兜底**：即使 LC 分配在某个点上不正确，试除验证会捕获，继续下一个点。
