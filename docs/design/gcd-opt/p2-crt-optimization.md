# P2: CRT 循环优化 — 详细设计

> 日期：2026-03-05
> 分支：`feature/gcd-optimization-p1`
> 依赖：P0（大素数）已完成

## 1. 问题分析

CLPoly 的模 GCD CRT 循环（`polynomial_gcd.cc:193-331`，`polynomial_gcd.hh:323-450`）存在两处常数因子浪费：

### 1.1 稳定检测：全项比较 O(n)

```cpp
// polynomial_gcd.cc:301
if (tmp_Pout_==Pout_)   // 比较所有 n 项系数
```

每轮 CRT 合并后都做全量 `==` 比较。典型 deg200 多项式有 ~200 项，每项是 ZZ 比较。

**参考实现**：
- FLINT：跟踪 `curr_bits = max_bits(coefficients)`，仅比较两个整数。`new_bits == curr_bits` 才做进一步检查。
- NTL：CRT 函数内部维护 `modified` 标志，在合并过程中检测 `h != 0`（CRT 增量是否为零）。若所有增量为零，返回 `modified = 0` 表示稳定。

### 1.2 content 计算：无 early exit

```cpp
// upolynomial.hh:220-230
template <class Tc>
Tc cont(const upolynomial_<Tc> & F)
{
    auto I = (ptr++).second;
    for (;ptr!=F.end();++ptr)
        I = gcd(I, ptr->second);   // 即使 I==1 也继续循环
    return I;
}
```

content 为 1 时（primitive 多项式，CRT 重构的常见情况），应在 `I == 1` 时立即退出。

**参考实现**：
- NTL `content()` (ZZX1.cpp)：`if (IsOne(res)) break;` — 通常 1-2 项即退出。

## 2. 优化方案

### P2a: CRT 内联稳定检测

**思路**：在 CRT 合并过程中直接检测"是否所有系数都未变化"，而非事后全量比较。

**当前代码**（`polynomial_gcd.cc:242-297`）：CRT 合并构建 `tmp_Pout_`，然后 line 301 做 `tmp_Pout_ == Pout_`。

**优化**：在 CRT 合并循环内维护 `bool modified = false` 标志。当某系数的 CRT 增量 `h ≠ 0` 时设 `modified = true`。合并完成后直接判断 `modified`，无需额外比较。

这与 NTL 的方法一致。关键洞察：CRT 合并公式为

```
new_coeff = old_coeff + (image_coeff - old_coeff) * inv(old_modulus, new_prime) * old_modulus
```

当 `image_coeff ≡ old_coeff (mod new_prime)` 时，增量为零，系数不变。检测增量是否为零是 CRT 合并的自然副产物。

**改动位置**：`polynomial_gcd.cc:242-301`（单变量），`polynomial_gcd.hh:356-416`（多变量）

**具体实现**（单变量，`polynomial_gcd.cc`）：

```cpp
// 在 CRT 合并循环中加入 modified 标志
bool modified = false;
// ... 现有合并逻辑 ...

// 在每个系数合并处，用差值整除检测稳定性
// 当 Pout_ptr->first == Pm_ptr->first 时：
ZZ diff = Pm_ptr->second.number() - Pout_ptr->second;
if (diff % prime != 0) modified = true;
tmp_Pout_.push_back({std::move(Pm_ptr->first),
    Pout_ptr->second + diff * tmp_inv.number() * Pout_prime});

// 当只有 Pout 有项（Pm 中缺失，隐式 image = 0）时：
if (Pout_ptr->second % prime != 0) modified = true;
// ... 原有 CRT 计算 ...

// 当只有 Pm 有项（Pout 中缺失）时：
modified = true;  // 新项 = 一定改变

// 合并完成后：
Pout_prime *= prime;
// balanced reduction ...

if (!modified)   // 替代 tmp_Pout_ == Pout_
{
    // 进入 content + division test
}
```

**数学正确性**：

稳定条件是 `old ≡ image (mod prime)`，等价于 `prime | (image - old)`。

- `(image - old) % prime == 0` 检测此条件。检测 `== 0` 不受 truncated/floor division 符号差异影响（`a % b == 0 ⟺ b | a` 对两种语义均成立）。
- 若对所有系数 i 均有 `prime | (image_i - old_i)`，则 CRT 增量 `delta_i = (image_i - old_i) * inv * M` 被 `prime * M` 整除，balanced reduction 后 `new_i == old_i`。因此 `!modified` 等价于 `tmp_Pout_ == Pout_`。
- 反之，若 `prime ∤ (image_i - old_i)`，则 `delta_i` 不被 `prime * M` 整除，`new_i ≠ old_i`。

**注意**：不能用 `delta == 0`（ZZ 层面）替代 `diff % prime == 0`。反例：`old = -3, prime = 7, image = 4`，`diff = 7 ≠ 0` 但 `7 % 7 == 0`（系数实际未变）。

### P2b: content early exit + lc_gcd == ±1 快速路径

**改动 1**：`upolynomial.hh` 的 `cont()` 函数加 early exit。

```cpp
template <class Tc>
Tc cont(const upolynomial_<Tc> & F)
{
    if (F.empty())
        return 1;
    auto ptr = F.begin();
    auto I = (ptr++).second;
    for (; ptr != F.end(); ++ptr)
    {
        I = gcd(I, ptr->second);
        if (I == 1 || I == -1) return I;   // ← early exit
    }
    return I;
}
```

**数学正确性**：`gcd(c₀, ..., cₙ) | gcd(c₀, ..., cₖ)`。若 `gcd(c₀, ..., cₖ) = ±1`，则全体 GCD 也是 ±1。

**影响范围**：所有调用 `cont()` 的地方（GCD、因式分解等）都自动受益。改动安全，无副作用。

**改动 2**：`lc_gcd == ±1` 时跳过 content 计算。

当 `gcd(lc(F), lc(G)) == ±1` 时，CRT 重构的多项式已通过 LC 归一化，结果是 primitive 的。

在 `polynomial_gcd.cc` 的稳定检测分支中：

```cpp
if (!modified)
{
    if (Pout_d > 0)
    {
        if (lc_gcd == 1 || lc_gcd == -1)
        {
            // LC GCD 为 ±1 → 结果已是 primitive，跳过 content
            // 直接做 trial division
            pair_vec_div(tmp.data(), R.data(), F.data(), Pout_.data(), F.comp());
            if (R.empty())
            {
                pair_vec_div(tmp.data(), R.data(), G.data(), Pout_.data(), F.comp());
                if (R.empty())
                {
                    tmp = Pout_ * upolynomial_<ZZ>({{0, cont_gcd}});
                    return tmp;
                }
            }
        }
        else
        {
            // 现有逻辑：cont() + division
            auto cont_ = cont(Pout_);
            // ...
        }
    }
}
```

**数学正确性**：FLINT 的 `g_pm1` 快速路径基于相同原理。当 `gcd(lc(A), lc(B)) = 1` 时，CRT 归一化后 `lc(C_p) = gcd(lc(A), lc(B)) = 1` 对所有好素数成立。CRT 重构的多项式首项系数为 1（monic），content 必为 1。

**注意**：对多变量情况（`polynomial_gcd.hh`），`lc_gcd` 是多项式而非整数。`lc_gcd == ±1` 的快速路径仅适用于单变量。多变量可用 `is_number(lc_gcd)` 判断。

### P2c: trial division 预检查

**思路**：在完整多项式除法之前，用廉价的整数整除检查快速拒绝错误候选。

**参考实现**：NTL `PlainDivide` (ZZX1.cpp) 在做多项式除法前依次检查：
1. `lc(candidate) | lc(F)` — 首项系数整除
2. `ConstTerm(candidate) | ConstTerm(F)` — 常数项整除
3. content 整除检查

若任一检查失败，立即返回"不整除"，跳过 O(n²) 的完整除法。

**改动位置**：`polynomial_gcd.cc` trial division 之前（line 310），`polynomial_gcd.hh` 对应位置

**具体实现**（单变量）：

```cpp
// 在 pair_vec_div 之前加预检查
// candidate = tmp_Pout_（已除 content），被除数 = F（primitive）
// 检查 lc(candidate) | lc(F)
if (F.begin()->second % tmp_Pout_.begin()->second != 0)
    goto division_failed;
// 检查常数项整除（如果都有常数项）
if (!tmp_Pout_.empty() && !F.empty()) {
    ZZ cand_const = tmp_Pout_.back().first.deg() == 0 ? tmp_Pout_.back().second : ZZ(0);
    ZZ f_const = F.back().first.deg() == 0 ? F.back().second : ZZ(0);
    if (cand_const != 0 && f_const % cand_const != 0)
        goto division_failed;
}
pair_vec_div(tmp.data(), R.data(), F.data(), tmp_Pout_.data(), F.comp());
if (!R.empty()) goto division_failed;
// ... 同理检查 G ...
```

**数学正确性**：若 `H | F`，则 `lc(H) | lc(F)` 且 `H(0) | F(0)`。逆否命题：若首项系数或常数项不整除，则 H 不整除 F。

**成本**：两次 ZZ 取模（O(1)），远低于 `pair_vec_div` 的 O(n²)。

**收益场景**：当 CRT 结果假稳定（modulus 不够大）时，预检查可在 O(1) 内拒绝，避免无效的完整除法。对 64-bit 大素数，假稳定较少发生，但当发生时节省显著。

## 3. 改动文件清单

| 文件 | 改动 | 范围 |
|------|------|------|
| `clpoly/upolynomial.hh` | P2b: `cont()` 加 early exit | 1 行 |
| `clpoly/polynomial_gcd.cc` | P2a: CRT modified 标志 + P2b: lc_gcd 快速路径 + P2c: 预检查 | ~40 行 |
| `clpoly/polynomial_gcd.hh` | P2a: 多变量 CRT modified 标志 + P2c: 预检查 | ~25 行 |

## 4. 预期效果

| 优化 | 场景 | 预期收益 |
|------|------|---------|
| P2a CRT 内联检测 | 高次多项式（deg200+） | 省去 O(n) 全项比较，但 CRT 本身是 O(n)，实际节省 ~30% CRT 开销 |
| P2b content early exit | primitive 多项式（常见） | content 从 O(n) 降到 O(1)，但仅在稳定时触发一次 |
| P2b lc_gcd 快速路径 | lc(F), lc(G) 互素（常见） | 完全跳过 content 计算 |
| P2c trial division 预检查 | 假稳定时 | 从 O(n²) 降到 O(1) 拒绝，但假稳定本身罕见 |

总体预期：单变量 GCD 加速 **1.1-1.3x**（CRT 迭代次数不变，每次迭代常数因子减小）。对高次多项式效果更明显。

## 5. 测试方案

1. 全量单元测试：`bash test/run_all_tests.sh`
2. 全量 crosscheck：`make crosscheck`
3. 性能基准：`make bench-clpoly` 对比优化前后
4. 重点关注：`gcd deg200+common100`、`gcd deg500+common250`（高次用例）

## 6. 完成后对比

| 技术 | FLINT | NTL | CLPoly 现状 | CLPoly P2 后 |
|------|-------|-----|------------|-------------|
| 稳定检测 | `new_bits == curr_bits`（单整数比较） | CRT 返回 `modified` 标志（零成本） | `tmp_Pout_ == Pout_`（O(n) 全项比较） | CRT 内联 `modified` 标志（同 NTL） |
| content 计算 | 仅在稳定触发时计算 | early exit：`IsOne(res)` | 全量 `cont()`，无 early exit | early exit：`I == ±1`（同 NTL） |
| lc\_gcd == ±1 快速路径 | `g_pm1` → 跳过 content | 无 | 无 | `lc_gcd == ±1` → 跳过 content（同 FLINT） |
| trial division 预检查 | 无（用 bit bound 兜底） | LC 整除 + 常数项整除 | 直接做完整除法 | LC 整除 + 常数项整除（同 NTL） |

## 7. 一手资料出处

- **FLINT** `src/fmpz_poly/gcd_modular.c`：`new_bits == curr_bits` 稳定检测，`g_pm1` 快速路径
- **NTL** `src/ZZX1.cpp`：`CRT()` 返回 `modified` 标志，`content()` 的 `IsOne(res)` early exit，`PlainDivide` 的 LC/常数项预检查
- **GCL** Algorithm 7.1 (p.307)：模 GCD 算法框架
