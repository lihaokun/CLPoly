# 修正方案：non-divisors 预过滤回归 + 因式分解模块全面审查

> 对应 workflow.md §5.1 修正方案文档

---

## 第一部分：复现与定位

### 压力测试发现 8 个失败用例

gamma coprimality 修复后压力测试 (N=1300)：
- 三变量 2因子: 0/500 (0%) — **已修复**
- 二变量 2因子: 1/198 (0.5%)
- 二变量 3因子: 4/193 (2.1%)
- 三变量 3因子: 2/193 (1.0%)
- 非平凡lc: 1/89 (1.1%)

### 代表性失败追踪

#### bivar2_82: L=2y, gamma=2

| skip | cs  | lj_val | gcd(γ,lj) | gcd(γcs,lj) | lc_ok |
|------|-----|--------|-----------|-------------|-------|
| 2    | -1  | -3     | 1         | 1           | 0     |
| 3    | 1   | 3      | 1         | 1           | 0     |
| 10   | -1  | -7     | 1         | 1           | 0     |
| 11   | 1   | 7      | 1         | 1           | 0     |

**`gcd(γ,lj)=1` 且 `gcd(γcs,lj)=1` 的点仍然 lc_ok=0。**

#### nontrivlc_198: L=6y+10, gamma=2, LC因子 (3y+5)

| skip | cs | lj_val | gcd(γ,lj) | gcd(γcs,lj) | lc_ok |
|------|-----|--------|-----------|-------------|-------|
| 0    | 1   | 5      | 1         | 1           | 0     |
| 3    | 1   | 11     | 1         | 1           | 0     |
| 7    | 1   | 17     | 1         | 1           | 0     |

y-main: L=20x²-25, gamma=5, LC因子 (4x²-5)。lj_val 是多变量多项式的求值，值很大（11, 31, 59, ...），全部 lc_ok=0。

#### bivar3_15: 混合成功/失败

x-main: gamma=-4, LC因子 (y+2)。**11/20 成功，9/20 失败。**
- 成功点: skip=1 (lj=5,cs=-2), skip=4 (lj=-3,cs=2), ...
- 失败点: skip=9 (lj=9,cs=-2), skip=13 (lj=11,cs=-2), skip=16 (lj=-9,cs=2), ...

这说明 **bivar3_15 本身可以分解**，只是某些好点被 non-divisors 过滤掉了。

---

## 第二部分：根因分析

### 问题 1: non-divisors 预过滤过度拒绝

SymPy 的 `dmp_zz_wang_non_divisors` 条件：`|lj_val| <= |d_stripped|`
其中 `d_stripped = cs·γ` 去掉与 lj_val 的公因子。

对 bivar2_82 skip=2 (gamma=2, cs=-1, lj_val=-3)：
```
d = cs * gamma = -2
strip lj_val=-3 from d=-2: -2 % -3 ≠ 0 → d stays -2
check: |−3| <= |−2| → 3 <= 2 → false → REJECT
```

**但 gcd(γ, lj_val)=1，提取不会过度。** 验证：
```
w[i] = lc(u_i) * cs = lc(u_i) * (-1)
prod(w) = (-1)^r * prod(lc(u_i))
identity: γ * lj_val = cs * prod(lc(u_i)) → 2*(-3) = -1 * prod(lc) → prod(lc) = 6
提取 lj_val=-3: 能从某个 lc(u_i) 中提取恰好 1 次 → total_k=1=ej ✓
```

**根因**：SymPy non-divisors 算法为 SymPy 自身的 LC 分配算法设计，与 CLPoly 的"逐因子幂次提取"算法语义不同。SymPy 的条件是：确保 `lj_val` 的大小不超过 `cs·γ` 的大小（用于后续的"除法分配"）。CLPoly 不需要这个条件，CLPoly 只需要确保 **`lj_val` 的素因子不与 `γ` 共享**（否则 γ 的因子泄漏到所有 w[i] 中导致过度提取）。

### 问题 2: 多变量 LC 因子的"无穷"lj_val

nontrivlc_198 y-main: LC因子 `4x²-5`，lj_val = 4α²-5。
- x=2: lj_val=11, x=3: lj_val=31, x=4: lj_val=59, ...
- 这些值随 α 二次增长，但 `gcd(γ,lj_val)` 大多为 1
- **所有 20 个点全部 lc_ok=0**，即使 non-divisors 不拒绝

这说明还有 **non-divisors 之外的根因**。`total_k != ej` 在 gcd 正常的情况下也会发生。

### 问题 3: 提取算法本身的缺陷

当 LC 因子是多变量多项式时（如 `3x+4`, `4x²-5`），lj_val 的绝对值可以很大。
提取循环 `while (w[i] % lj_val == 0)` 在 lj_val 很大时几乎不会触发，因为 w[i] = lc(u_i) * cs 的绝对值通常远小于 lj_val。

结果: total_k = 0 ≠ 1 = ej → 拒绝。

**这不是 γ 污染问题，而是"提取无法工作"问题**：当 lj_val 太大，w[i] 中根本没有 lj_val 的因子可提取。

### 根因总结

当前 `__wang_leading_coeff` 的三类失败模式：

| 失败模式 | 原因 | 是否已修复 |
|----------|------|-----------|
| A: γ 污染 | `gcd(γ, lj_val) > 1` 导致过度提取 | ✅ GCD 匹配天然免疫 |
| B: non-divisors 过度拒绝 | SymPy 条件与 CLPoly 提取算法不匹配 | ✅ 已移除，改用 GCD 匹配 |
| C: 大 lj_val 无法提取 | 多变量 LC 因子求值太大，w[i] 中无因子 | ✅ GCD 匹配对大值也能正确工作 |

---

## 第三部分：全面审查

### 3.1 与 SymPy/FLINT 的 LC 分配算法对比

**SymPy** (`dmp_zz_wang_lead_coeffs`):
- 不使用幂次提取
- 使用 `GCD(lj_val, w[i])` 进行匹配分配
- non-divisors 检查保证 GCD 匹配的唯一性

**FLINT** (`fmpz_mpoly_factor_lcc_wang`):
- 也不使用幂次提取
- 使用 `fmpz_remove(w[i], lj_val)` (类似 GCD)
- 有更复杂的 content 预处理

**CLPoly** (当前):
- 使用 `while (w[i] % lj_val == 0)` 提取最大幂次
- 检查 `total_k == ej`
- 这个方法在 lj_val 很大或有公因子时容易失败

**关键差异**: SymPy/FLINT 的匹配是基于 GCD/remove，对大值也能工作。CLPoly 的纯整除提取在大值时脆弱。

### 3.2 CLPoly 历次 bug 的共同根因

| 轮次 | Bug | 根因 |
|------|-----|------|
| D5 | LC 原子捆绑匹配 | 匹配算法不支持 content≠1 |
| D6 | 单项式 content | Wang 前缺少预处理 |
| D7 | gamma 互素性 | 搜索范围太小 (SCAN_SIZE=50) |
| D8 | non-divisors 回归 | 错误移植 SymPy 算法 |
| D8 | 大 lj_val | 提取算法对多变量 LC 因子本身不可靠 |

**共同根因**: CLPoly 的 LC 分配算法（逐因子幂次提取）与 SymPy/FLINT 的算法不同，
在多种边界条件下不可靠。逐个修补不如**按 SymPy 的 GCD 匹配重写**。

### 3.3 修正方案

#### 方案 A: 最小修复（只修 non-divisors 回归）

将 non-divisors 预过滤从 SymPy 条件改为仅检查 γ 互素性：
```cpp
// 替换当前的 non-divisors 块:
// 仅检查 gcd(γ, lj_val) > 1 → 拒绝
for (auto& [lj2, ej2] : lc_factors)
{
    auto lj2_eval = assign(lj2, eval_point);
    ZZ lj2_val = is_number(lj2_eval) ? lj2_eval.front().second : ZZ(0);
    if (lj2_val == 0 || abs(lj2_val) == ZZ(1)) return result;
    if (gcd(abs(gamma), abs(lj2_val)) != ZZ(1))
        return result;
}
```

**优点**: 修复 non-divisors 回归，不引入新问题。
**缺点**: 不解决"大 lj_val 无法提取"问题。

#### 方案 B: 重写 LC 分配（按 SymPy GCD 匹配）

将提取循环从 `while (w[i] % lj_val == 0)` 改为 SymPy 风格的 GCD 匹配：

```cpp
// SymPy 风格: 找唯一的 gcd(w[i], lj_val^ej) 最大的 i
for (auto& [lj, ej] : lc_factors)
{
    ZZ lj_val = ...; // 求值
    ZZ lj_pow = pow(abs(lj_val), ej);  // lj_val^ej

    size_t best_i = r;  // 无效标记
    ZZ best_g(0);
    for (size_t i = 0; i < r; ++i)
    {
        ZZ g = gcd(abs(w[i]), lj_pow);
        if (g > best_g) { best_g = g; best_i = i; }
        else if (g == best_g && g > ZZ(1)) { best_i = r; }  // 歧义
    }
    if (best_i >= r) return result;  // 无唯一候选

    // 将 lj_pow 分配给 sigma[best_i]
    sigma[best_i] *= lj^ej;
    w[best_i] /= best_g;
}
```

**优点**: 与 SymPy/FLINT 语义一致，对大 lj_val 也能正确匹配。
**缺点**: 改动较大，需要验证正确性。

#### 实施结果: 方案 B（GCD 匹配）

最终实施了方案 B，因为方案 A（仅修 non-divisors）仍存在盲区：
- bivar3_11: gamma=-4, lc_factors=[(y,1),(y-1,1)]，对任意整数 y，
  一定有 `gcd(4, y)>1` 或 `gcd(4, y-1)>1`，gamma 互素性检查永远失败。

方案 B 直接消除了三类失败模式的根因：
- **移除** non-divisors 预过滤和 gamma 互素性检查
- **移除** 幂次提取循环 (`while (w[i] % lj_val == 0)`)
- **替换为** GCD 匹配: `gcd(|w[i]|, |lj_val|^ej)` 最大的唯一 i

压力测试 (N=1300, 5 轮): **0 failures / 5 runs**。
