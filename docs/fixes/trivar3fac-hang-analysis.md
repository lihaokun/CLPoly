# trivar3fac 卡住问题分析

> 日期：2026-03-04
> 复现用例：trial 387 of test_repro_trivar3fac (500 trials, release)

## 1. 复现多项式

```
f1 = -3xy - x + 2
f2 = 3xy + 3yz - y
f3 = 2xy - 3yz - 3
f  = f1 * f2 * f3
```

## 2. 手动演算

### 2.1 主变量与 LC 分析

lex 序 x > y > z，主变量 = x，每个因子对 x 均为 1 次：

```
lc_x(f1) = -3y - 1
lc_x(f2) = 3y
lc_x(f3) = 2y

L = lc(f, x) = (-3y-1)(3y)(2y) = -6y²(3y+1)
```

因此：
- γ = -6 = -2·3
- l₁ = y，e₁ = 2
- l₂ = 3y+1，e₂ = 1

### 2.2 `__select_eval_point` 枚举

非主变量 = {y, z}。按 bound=0,1,2,... 枚举 (α_y, α_z)。

条件 (d')：|lⱼ(α)| ≥ 2 + 两两互素。

- l₁ = y → E₁ = |α_y|，要求 |α_y| ≥ 2
- l₂ = 3y+1 → E₂ = |3α_y + 1|

| α_y | E₁ | E₂ | coprime(E₁,E₂) | `__select_eval_point` 结果 |
|-----|----|----|-----------------|--------------------------|
| 0   | 0  | 1  | —               | 拒绝 (E₁=0, E₂≤1) |
| ±1  | 1  | 4/2| —               | 拒绝 (E₁≤1) |
| 2   | 2  | 7  | gcd(2,7)=1 ✓   | **通过** |
| -2  | 2  | 5  | gcd(2,5)=1 ✓   | **通过** |
| 3   | 3  | 10 | gcd(3,10)=1 ✓  | **通过** |
| -3  | 3  | 8  | gcd(3,8)=1 ✓   | **通过** |
| 4   | 4  | 13 | gcd(4,13)=1 ✓  | **通过** |
| 5   | 5  | 16 | gcd(5,16)=1 ✓  | **通过** |

大部分 |α_y| ≥ 2 的点都能通过 `__select_eval_point`。

### 2.3 `__wang_leading_coeff` coprime 检查

Valuation 提取的前置条件：`gcd(Eⱼ, |cs·γ|) = 1`。

γ = -6，|γ| = 6 = 2·3。因此 gcd(E₁, |cs·γ|) = gcd(|α_y|, 6·|cs|)。

| α_y | E₁ | gcd(E₁, 6) | `__wang_leading_coeff` 结果 |
|-----|----|-----------|-----------------------------|
| 2   | 2  | 2         | **失败** (coprime 不满足) |
| -2  | 2  | 2         | **失败** |
| 3   | 3  | 3         | **失败** |
| -3  | 3  | 3         | **失败** |
| 4   | 4  | 2         | **失败** |
| -4  | 4  | 2         | **失败** |
| 5   | 5  | 1         | **可能通过** (若 5 ∤ cs) |
| -5  | 5  | 1         | **可能通过** |
| 6   | 6  | 6         | **失败** |
| 7   | 7  | 1         | **可能通过** |

### 2.4 合法点密度

只有 |α_y| 与 6 互素的点才能通过：|α_y| ∈ {5, 7, 11, 13, 17, 19, 23, 25, ...}。

在 |α_y| ∈ [2, N] 范围内，不与 6 互素的比例 = 1 - φ(6)/6 = 1 - 2/6 = 2/3。即 2/3 的 α_y 值一定失败。

更关键的是：小值 2,3,4 全部失败。`__select_eval_point` 优先返回小值，所以前面数百个合法点都会被 `__wang_leading_coeff` 拒绝。

### 2.5 浪费的计算

每次 `__wang_leading_coeff` 失败（返回 !success），`__wang_core` 要完整执行：
1. `__select_eval_point` → 枚举求值点（快）
2. `assign(g, eval)` → 求值（快）
3. `factorize(f₀)` → **单变量因式分解**（慢！deg=3 不算慢，但累积量大）
4. `__wang_leading_coeff` → coprime 检查失败 → continue

每个主变量分配 BATCH_SIZE=200 个 skip。当 l₁=y 时，前 ~200 个合法点中大部分 α_y ∈ {2,3,4,6,8,9,...}，只有极少数 α_y ∈ {5,7,11,...} 能通过 coprime。

如果 200 个 skip 中没有一个 α_y 值与 6 互素，则该主变量 batch 用完 → 换主变量。

### 2.6 主变量旋转也不能解决

换主变量 y 后：
```
lc_y(f) = ...
```
lc(f, y) 可能包含类似的小素数因子，导致同样的 coprime 问题。

三个主变量轮流尝试，每个 200 次，一轮 600 次无效尝试后又回到第一个主变量（skip 从 200 开始），**无限循环**。

## 3. 根因

**`__select_eval_point` 和 `__wang_leading_coeff` 的 coprime 条件不一致。**

- `__select_eval_point` 检查：E_j 两两互素 + E_j ≥ 2
- `__wang_leading_coeff` 额外检查：gcd(E_j, |cs·γ|) = 1

当 LC 含纯变量因子（如 l₁ = y）且 γ 含小素数因子（如 6 = 2·3）时，大量求值点通过 `__select_eval_point` 但被 `__wang_leading_coeff` 拒绝，浪费大量计算。

更严重的是：这不是性能问题，而是**活性问题**（liveness）——200 个点的 batch 中可能没有任何一个能通过 coprime 检查，导致永远无法提升成功。

## 4. SymPy/FLINT 的处理方式

### 4.1 SymPy: 累积素因子剥离（non-divisor 算法）

SymPy 的 `dmp_zz_wang_non_divisors(E, cs, ct)` 算法：

```python
def dmp_zz_wang_non_divisors(E, cs, ct, K):
    result = [cs * ct]          # ct = γ
    for q in E:
        q = abs(q)
        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r      # 从 q 中剥离与 r 的公因子
            if K.is_one(q):
                return None     # q 被完全吸收 → 拒绝
        result.append(q)
    return result[1:]           # 返回各 Eⱼ 的独立残余
```

**关键区别**：SymPy 不要求 `gcd(Eⱼ, cs·γ) = 1`，只要求 Eⱼ 在剥离与 `cs·γ` 和之前 E 值的共享素因子后，残余 > 1。

### 4.2 FLINT: 相同算法 + Kaltofen 回退

FLINT 的 `fmpz_mpoly_factor_lcc_wang` 使用完全相同的累积剥离：

```c
fmpz_mul(d + 0, Auc, lcAfac->constant);   // d[0] = cs * γ
for (i = 0; i < lcAfac->num; i++) {
    fmpz_abs(Q, lcAfaceval + i);
    for (j = i; j >= 0; j--) {
        fmpz_set(R, d + j);
        while (!fmpz_is_one(R)) {
            fmpz_gcd(R, R, Q);
            fmpz_divexact(Q, Q, R);
            if (fmpz_is_one(Q)) { success = 0; goto cleanup; }
        }
    }
    fmpz_set(d + i + 1, Q);   // 残余
}
```

此外 FLINT 有 **Kaltofen 回退**：当 Wang LCC 失败时，调用 `lcc_kaltofen`（使用双变量因式分解），允许最多 4 次 Kaltofen 失败才换 alpha。

### 4.3 手算对比：α_y = 7

```
E₁ = |7| = 7,  E₂ = |3·7+1| = 22 = 2·11
cs·γ = cs·6 (假设 cs = 1, 则 cs·γ = 6)
```

**SymPy non-divisor**:
```
result = [6]

E₁ = 7:
  r = 6: gcd(6, 7) = 1 → break
  q = 7 > 1 ✓
  result = [6, 7]

E₂ = 22:
  r = 7: gcd(7, 22) = 1 → break
  r = 6: gcd(6, 22) = 2, q = 22/2 = 11; gcd(2, 11) = 1 → break
  q = 11 > 1 ✓
  result = [6, 7, 11]

→ return [7, 11]   ← 接受！残余 7 和 11 是独立素因子
```

**CLPoly coprime 检查**:
```
gcd(E₁=7, 6) = 1  ✓
gcd(E₂=22, 6) = 2 ≠ 1  → 拒绝！
```

**SymPy/FLINT 接受 α_y=7，CLPoly 拒绝。**

### 4.4 手算对比：α_y = 10

```
E₁ = 10 = 2·5,  E₂ = 31
```

**SymPy non-divisor**:
```
result = [6]

E₁ = 10:
  r = 6: gcd(6, 10) = 2, q = 10/2 = 5; gcd(2, 5) = 1 → break
  q = 5 > 1 ✓
  result = [6, 5]

E₂ = 31:
  r = 5: gcd(5, 31) = 1 → break
  r = 6: gcd(6, 31) = 1 → break
  q = 31 > 1 ✓

→ return [5, 31]   ← 接受！
```

**CLPoly**: `gcd(10, 6) = 2` → **拒绝！**

### 4.5 手算对比：α_y = 5（两边都拒绝）

```
E₁ = 5,  E₂ = |3·5+1| = 16 = 2⁴
```

**SymPy non-divisor**:
```
result = [6]

E₁ = 5:
  r = 6: gcd(6, 5) = 1 → break
  q = 5 > 1 ✓
  result = [6, 5]

E₂ = 16:
  r = 5: gcd(5, 16) = 1 → break
  r = 6: gcd(6, 16) = 2, q = 8; gcd(2, 8) = 2, q = 4; gcd(2, 4) = 2, q = 2; gcd(2, 2) = 2, q = 1
  q == 1 → return None   ← 拒绝！（E₂=2⁴ 没有独立于 6 的素因子）
```

**CLPoly**: `gcd(16, 6) = 2` → 拒绝

两边一致。E₂ = 16 = 2⁴ 的所有素因子都是 2，被 γ=6 完全吸收。

### 4.6 合法点密度对比

| α_y | E₁ | E₂ | SymPy | CLPoly |
|-----|----|----|-------|--------|
| 2   | 2  | 7  | 拒绝 (E₁=2 被 6 吸收) | 拒绝 |
| 3   | 3  | 10 | 拒绝 (E₁=3 被 6 吸收) | 拒绝 |
| 4   | 4  | 13 | 拒绝 (E₁=4=2² 被 6 吸收) | 拒绝 |
| 5   | 5  | 16 | 拒绝 (E₂=2⁴ 被 6 吸收) | 拒绝 |
| 6   | 6  | 19 | 拒绝 (E₁=6 被 6 吸收) | 拒绝 |
| **7**   | **7**  | **22=2·11** | **接受** (残余 7, 11) | **拒绝** |
| 8   | 8  | 25 | 拒绝 (E₁=2³ 被 6 吸收) | 拒绝 |
| 9   | 9  | 28 | 拒绝 (E₁=3² 被 6 吸收) | 拒绝 |
| **10**  | **10=2·5** | **31** | **接受** (残余 5, 31) | **拒绝** |
| **11**  | **11** | **34=2·17** | **接受** (残余 11, 17) | **拒绝** |
| 12  | 12 | 37 | 拒绝 (E₁=2²·3 被 6 吸收) | 拒绝 |
| **13**  | **13** | **40=2³·5** | **接受** (残余 13, 5) | **拒绝** |

**SymPy 从 α_y=7 开始就能找到合法点，CLPoly 要等到 α_y 与 6 互素 AND E₂ 与 6 互素——极其稀疏。**

实际上 CLPoly 需要同时满足：
- |α_y| 与 6 互素 → α_y ∈ {5, 7, 11, 13, ...}
- |3α_y+1| 与 6 互素 → 3α_y+1 不含因子 2 或 3

检查：
- α_y=5: 3·5+1=16=2⁴ → gcd(16,6)=2 ✗
- α_y=7: 3·7+1=22=2·11 → gcd(22,6)=2 ✗
- α_y=11: 3·11+1=34=2·17 → gcd(34,6)=2 ✗
- α_y=13: 3·13+1=40=2³·5 → gcd(40,6)=2 ✗

**注意到 3α_y+1 对任何奇数 α_y 都是偶数！** 因为 3·odd+1 = odd+1 = even。

所以当 α_y 是奇数时，E₂ 一定含因子 2，而 6 也含因子 2 → CLPoly 永远拒绝。

当 α_y 是偶数时，E₁ = |α_y| 含因子 2，而 6 也含因子 2 → CLPoly 永远拒绝。

**结论：对于 l₁=y, l₂=3y+1, γ=-6 的情况，CLPoly 的 coprime 条件不可能同时对两个 LC 因子满足。CLPoly 会无限循环。** 这是一个严格的不可满足条件（impossible condition），不是概率稀疏问题。

### 4.7 数学证明：CLPoly 条件不可满足

设 α_y = a（整数，|a| ≥ 2）。

- E₁ = |a|
- E₂ = |3a + 1|

CLPoly 要求：gcd(E₁, 6) = 1 且 gcd(E₂, 6) = 1。

**情况 1**：a 是偶数 → E₁ 含因子 2 → gcd(E₁, 6) ≥ 2 → 失败。
**情况 2**：a 是奇数 → 3a 是奇数 → 3a+1 是偶数 → E₂ 含因子 2 → gcd(E₂, 6) ≥ 2 → 失败。

**两种情况都失败，条件不可满足。QED.**

## 5. 根因修正

### 5.1 问题本质

CLPoly 的 `gcd(Eⱼ, cs·γ) = 1` 条件过于严格。正确的条件（Wang 原始算法 / SymPy / FLINT）是累积素因子剥离后残余 > 1，而不是完全互素。

### 5.2 修正方案

将 `__wang_leading_coeff` 中的 coprime 前置检查替换为 SymPy/FLINT 的 non-divisor 算法：

**当前代码**（过于严格）:
```cpp
ZZ cs_gamma = abs(uni_content * gamma);
for (size_t j = 0; j < m; ++j)
    if (gcd(lc_evals[j], cs_gamma) != ZZ(1))
        return result;  // 拒绝
```

**修正为** non-divisor 算法：
```cpp
// d[0] = |cs * γ|
std::vector<ZZ> d = { abs(uni_content * gamma) };
for (size_t j = 0; j < m; ++j) {
    ZZ q = lc_evals[j];  // 已经是 abs
    for (int k = (int)j; k >= 0; --k) {
        ZZ r = d[k];
        while (r != ZZ(1)) {
            r = gcd(r, q);
            q = q / r;
            if (q == ZZ(1)) return result;  // 被完全吸收 → 拒绝
        }
    }
    d.push_back(q);  // 残余
}
// 通过 non-divisor 检查，继续 valuation 提取
```

### 5.3 Valuation 提取的适配

non-divisor 算法通过后，valuation 提取使用 `d[j+1]`（剥离后的残余）而不是原始 `lc_evals[j]`。因为残余保证了各 Eⱼ 之间没有共享素因子，valuation 提取的唯一性得到保证。

（待实现确认）
