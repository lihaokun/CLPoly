# 调研报告：P1b 线性 Hensel 提升

> 调研日期：2026-02-22
> 目标：实现增量 Hensel 提升 + 早期因子检测（Monagan 2019）
> 状态：基于 FLINT 源码分析、Monagan 2019 论文摘要、已有 CLPoly 代码阅读完成

---

## 1. 背景与动机

### 1.1 当前单变量因式分解管线

```
factorize(f ∈ Z[x])
  → squarefreefactorize(f)
  → 对每个无平方因子 g:
      → __select_prime(g)                §8.3  选 30-50 个素数，取模因子数最少者
      → __hensel_lift(g, factors_Zp, p)  §6.6  二次 Hensel 提升 → (lifted_ZZ[], m)
      → __factor_recombine(g, lifted, m) §7.5/8.1  van Hoeij LLL 重组（P1a）
```

其中 `__hensel_lift` 当前实现为**二叉树二次提升**：

```cpp
// polynomial_factorize_univar.hh:520-526
ZZ m(p);
while (m <= target)           // target = 2·lc(f)·B_Mig(f)
{
    __hensel_lift_recursive(nodes, 0, f, m);
    m = m * m;                // 每轮精度平方: p → p² → p⁴ → ...
}
```

### 1.2 Mignotte 界的过提升问题

**Mignotte 界**（`__mignotte_bound`，polynomial_factorize_univar.hh:122）：

```
B = C(n, ⌊n/2⌋) · ‖f‖₂
```

提升目标为 `m > 2·lc(f)·B`，即 `p^k > 2·lc(f)·B`，故 `k = ⌈log_p(2·lc(f)·B)⌉`。

**uni-70 案例分析**：`f = (x-1)(x-2)···(x-70)` ∈ Z[x]，deg = 70。

- `‖f‖₂ ≈ 70!^(1/2)`，对数约 `70·log₂(70)/2 ≈ 228 bit`
- `C(70, 35) ≈ 2^69`（中间二项式系数，约 69 bit）
- `B ≈ 2^(228+69) = 2^297`
- `lc(f) = 1`，故 target ≈ `2^298`，约 300 bit
- 取 p = 2 时，需 `k = 298` 步；取 p ≈ 2^30（GMP 多精度）时，需 k ≈ 10 步

实际 CLPoly 选用 p ≈ 几千的素数，k ≈ 6 轮（每轮精度平方：p → p² → p⁴ → p⁸ → p¹⁶ → p³² > target）。每轮的二次提升在已达到 400+ bit 精度的算术下进行，代价高昂。

**关键问题**：对于 uni-70 这种 *均匀稀疏因子*（每个因子 x-i 次数仅为 1），van Hoeij LLL 实际上在精度很低时（约 60–100 bit）就能收敛。二次提升到 400 bit 是不必要的，导致**大量精度浪费**。

### 1.3 benchmark 数据

| 测试用例 | CLPoly | FLINT | 比值 |
|---------|--------|-------|------|
| f = (x-1)(x-2)···(x-70) | 83 ms | 4.5 ms | ~18x |

P1a（van Hoeij LLL）实现后，重组阶段已基本消除，瓶颈将转移到 Hensel 提升阶段。P1b 的目标是将精度从"全 Mignotte 界"降低到"LLL 收敛所需的最小精度"。

---

## 2. 算法：线性 Hensel 单步

### 2.1 基本设置

设 `f ∈ Z[x]`，`g, h ∈ Z[x]` 满足：

```
f ≡ g·h  (mod p^a)
s·g + t·h ≡ 1  (mod p)           ← Bézout 系数仅在 mod p 下保持
deg(s) < deg(h), deg(t) < deg(g)  ← 次数约束（分而治之后的简化形式）
```

目标：计算 `g', h'` 使得：

```
f ≡ g'·h'  (mod p^(a+1))
g' ≡ g  (mod p^a)
h' ≡ h  (mod p^a)
```

### 2.2 线性提升步骤（两因子版本）

**步骤 1**：计算误差项并约化

```
e = (f - g·h)  /  p^a   ∈ Z/p[x]    （精确整除，因为 f ≡ g·h mod p^a）
```

注意 e 在 Z/p 上（不是 Z/p^a），这是线性提升与二次提升的关键区别：
- **线性提升**：e ∈ Z/p（当前精度一层的误差）
- **二次提升**：e ∈ Z/p^a（当前精度全层的误差）

**步骤 2**：用 Bézout 系数求修正量

```
设 σ, τ ∈ Z/p[x] 使得  σ·g + τ·h ≡ e  (mod p)
```

具体计算：

```
σ = s·e mod h    （在 Z/p[x] 中，做带余除法，取余数，确保 deg σ < deg h）
τ = t·e + q·g    其中 q 是 s·e 除以 h 的商
```

验证：`σ·g + τ·h = (s·e - q·h)·g + (t·e + q·g)·h = s·e·g + t·e·h = e·(sg+th) = e·1 = e  (mod p)` ✓

**步骤 3**：更新因子

```
g' = g + p^a · τ    (mod p^(a+1))
h' = h + p^a · σ    (mod p^(a+1))
```

验证：

```
g'·h' = (g + p^a·τ)·(h + p^a·σ)
      = g·h + p^a·(g·σ + h·τ) + p^(2a)·τ·σ
      ≡ g·h + p^a·e   (mod p^(a+1))   （因为 2a ≥ a+1 当 a ≥ 1，故 p^(2a)·σ·τ ≡ 0 (mod p^(a+1)) 消失）
      ≡ f              (mod p^(a+1))   ✓
```

**注意**：当 a=1 时 2a = 2 = a+1，故 `p^(2a) = p^(a+1) ≡ 0 (mod p^(a+1))` 平凡成立；当 a ≥ 2 时 2a > a+1，同样成立。故对所有 a ≥ 1 均正确。

### 2.3 Bézout 系数的处理（线性 vs 二次的核心区别）

| 特性 | 线性提升 | 二次提升 |
|------|--------|--------|
| Bézout 系数 s, t 精度 | 保持 **mod p**（固定不变） | 提升至 mod p^a（每步更新） |
| e 的精度 | Z/p 上（1 层） | Z/p^a 上（全层） |
| 每步代价（多项式乘法） | O(n) × O(log p)  | O(n) × O(a·log p) |
| Bézout 更新代价 | 无（不需要更新） | O(n²) 每步 |
| 总步数 | d/log p 步 | log d 步 |

**关键洞察**：在线性提升中，s 和 t **始终保持 mod p**，不需要随精度增长而更新。这意味着：
1. 每步只做 mod p 下的多项式运算（系数是单精度整数）
2. 修正量 σ, τ 是 mod p 下的小多项式，乘以 p^a 再加入 g, h
3. g, h 的系数在 Z/p^(a+1) 下（需要大整数），但修正量计算本身用小整数

**Monagan 2019 的核心贡献**：证明线性提升中 Bézout 系数不需要更新，且通过适当的算法设计（矩阵方法、移位多项式技巧）将总复杂度从 O(n²d²)（经典线性）降至 O(nd²)（立方），比二次提升的 O(nd²·log d) 在中等规模更优。

### 2.4 多因子情形（r > 2 因子）

当有 r ≥ 2 个因子时，CLPoly 当前使用**二叉树**结构（`__hensel_node`）将 r 因子的问题归约为 log₂r 层的两因子问题。

线性提升可以同样应用于每个内部节点，但 Monagan 2019 和 2022 提出了直接处理 n 个因子的方法，利用以下 Diophantine 方程系统：

```
Σᵢ σᵢ · ∏_{j≠i} hⱼ ≡ e   (mod p)
deg σᵢ < deg hᵢ
```

其解通过 `σᵢ = sᵢ·e mod hᵢ` 得到（Bézout 链 `{sᵢ}` 在 mod p 下预计算固定）。

**CLPoly P1b 方案（架构评估后更新）**：**采用 n 因子直接方法**，不保留二叉树结构。评估结论：

| 维度 | 二叉树方案 | n 因子直接方案 |
|------|-----------|--------------|
| 代码结构 | 树递归（log r 层） | 平铺循环（更简单） |
| 预计算 | `__hensel_tree_build`（已有）| Bézout 链（新写，约 30 行）|
| 每步误差 | 每树节点各算一次 | 全局仅算一次 |
| Bézout 系数数量 | 2 × (r-1) 个树节点对 | r 个直接系数 |
| 与 Monagan 2022 兼容 | 不兼容（树 vs 矩阵） | 直接兼容（v2 矩阵优化路线）|

n 因子方案不仅实现更简洁，且为后续 Monagan 2022 矩阵优化（v2）铺路，因此选择此方案。详见 `docs/design/vanhoeij/architecture.md` Part B §7 M0-M1。

---

## 3. 完整算法：线性提升 + van Hoeij 交织（Maple 2019）

### 3.1 算法框架

```
算法：LinearHenselWithEarlyDetection(f, h₁,...,hᵣ mod p, p)
前置：f ∈ Z[x]，hᵢ ∈ Z/p[x]，f ≡ ∏hᵢ (mod p)，gcd(hᵢ, hⱼ) = 1 (mod p) ∀i≠j

输入：f, 模因子列表 [h₁,...,hᵣ] in Z/p[x], 素数 p
输出：f 的不可约因子列表 ∈ Z[x]

1. 计算初始精度 a₀（见 §3.2）
2. 线性提升 a 步（a₀ 步），得到 g₁,...,gᵣ mod p^a
3. 进入交织循环：
   While 未完成：
       3a. 构建 CLD 矩阵（van Hoeij P1a 部分）
       3b. 运行 LLL，检测短向量
       3c. 若 LLL 找到因子：提取、验证、从 active 移除
           若所有因子找到：return 结果（早期终止）
       3d. 若 LLL 未收敛：继续提升 k 步（增量）
           精度 a → a + k
           用当前精度的 gᵢ mod p^a 再提升 k 步得 gᵢ mod p^(a+k)
       3e. 若 a > Mignotte 界所需精度且仍未收敛：返回 Zassenhaus（极罕见）
```

### 3.2 初始精度 a₀ 的启发式选取

FLINT 使用以下启发式（`_heuristic_van_hoeij_starting_precision`，factor_van_hoeij.c）：

```c
a_heuristic = (2.5*r + min_b) * log(2) + log(f->length) / 2.0 / log(p)
```

其中：
- `r` = 模因子数
- `min_b` = r 个因子中最小绝对系数的 bit 数（约等于 log₂(p)）
- `f->length` = f 的项数（多项式长度）

**最终初始精度**：

```c
a = min(fmpz_clog_ui(B, p),             // Mignotte 精度 ⌈log_p B⌉
        a_heuristic)                     // 启发式精度
```

即取 Mignotte 精度和启发式精度中的**较小值**作为初始值，保守地比 LLL 收敛所需精度稍大，避免过多不必要的初始提升。

**CLPoly 策略**：初版可以从较保守的 `a₀ = ⌈log_p(p^10)⌉ = 10`（即提升到 10 步），然后每次 +k 步直到 LLL 收敛。也可以直接使用 FLINT 启发式（在 P1b 详细设计阶段确定）。

### 3.3 精度倍增安全网（FLINT 外层循环）

FLINT 的完整外层结构（从 `factor_van_hoeij.c` 提取的逻辑）：

```c
// 初始提升
prev_exp = _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, fac, a);

// 主循环：直到找到完整分解
while (!fmpz_poly_factor_van_hoeij_check_if_solved(...))
{
    // 内层：逐步喂入 CLD 列，运行 LLL
    // （内螺旋列喂入，与 P1a 相同）

    // 若内层 LLL 仍未收敛：精度倍增（quadratic scaling back）
    prev_exp = _fmpz_poly_hensel_continue_lift(lifted_fac, link, v, w,
                                               f, prev_exp, a, 2*a, fp);
    a = 2*a;    // a 倍增（a *= 2，而非 a += k）
    fmpz_pow_ui(P, fp, a);    // P = p^a 更新
}
```

**关键观察**：FLINT 在 van Hoeij 中使用的 Hensel 提升外层仍然是**精度倍增**（a *= 2），不是单步线性（a += 1）。这与 Monagan 2019 的立方代价 LHL 是不同的两个概念：
- FLINT `factor_van_hoeij.c` 的 "双精度安全网"：每次外层循环精度翻倍（从当前 a 到 2a），依赖 `_fmpz_poly_hensel_continue_lift` 的高效实现
- Monagan 2019 的"线性提升"：内部提升算法本身每步仅增加一阶（从 a 到 a+1），Bézout 系数保持 mod p

因此 P1b 的设计需要区分两个层次：
1. **提升策略层**：每次触发提升时增加多少精度（FLINT 用倍增，Monagan 2019 用 +1 或 +k）
2. **单步提升算法**：如何高效完成单次精度增量（线性 vs 二次）

### 3.4 早期终止条件

早期终止触发当 LLL 在当前精度下找到所有因子时：

```
__vanhoeij_recombine 中：
    若 active.size() == 1：
        最后一个因子 = f_star（剩余商）
        return result   ← 早期终止，未达完整 Mignotte 界
```

对于 uni-70（70 个线性因子），每个因子 (x-i) 在 LLL 中对应独立的短向量，理论上只需约 10–15 步线性提升（精度约 10·log₂p ≈ 150 bit）即可收敛，而不是完整的 300 bit Mignotte 精度。

---

## 4. FLINT 实现分析

### 4.1 `factor_van_hoeij.c` 结构解析

**函数 `_heuristic_van_hoeij_starting_precision`**（约 15 行）：
- 输入：`f`（多项式），`r`（模因子数），`p`（素数）
- 计算：`a = (2.5*r + min_b) * log(2) + log(f->length) / 2.0 / log(p)`
- 作用：估计 LLL 收敛所需最小精度，避免从 Mignotte 界开始（省去大量无效提升）

**函数 `fmpz_poly_factor_van_hoeij`**（主入口，约 150 行）：

```c
// 初始化
a = fmpz_clog_ui(B, p);                                     // Mignotte 精度
a = FLINT_MIN(a, _heuristic_van_hoeij_starting_precision(f, r, p));  // 取启发式最小值

// 初始提升（到精度 a）
prev_exp = _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, fac, a);

// 外层主循环
while (!check_if_solved(...))
{
    // 内层：CLD 矩阵构建 + 逐列喂入 + LLL
    // （与 P1a 算法完全相同，见 vanhoeij-lll-research.md §2.4）

    // 若内层未解决：精度倍增
    prev_exp = _fmpz_poly_hensel_continue_lift(lifted_fac, link, v, w,
                                               f, prev_exp, a, 2*a, fp);
    a = 2*a;
    fmpz_pow_ui(P, fp, a);
}
```

**`check_if_solved` 函数**：检查是否所有 r 个模因子已被归入某个确认的真因子子集，即 active 集为空或剩余 f_star 为常数。

### 4.2 `CLD_mat.c` 结构解析

函数 `_fmpz_poly_factor_CLD_mat`（约 120 行），构造 (r+1)×(lo_n + hi_n) 数据矩阵：

```c
// 1. 预计算 CLD 界（末行）
for (i = 0; i < k; i++)
{
    fmpz_poly_CLD_bound(entry(r, i),       f, i);           // 低次端第 i 列界
    fmpz_poly_CLD_bound(entry(r, 2k-i-1), f, f->length-i-2); // 高次端第 i 列界
}

// 2. 低次端列过滤（lo_n）
bound = fmpz_bits(P) - bit_r - bit_r/2;     // bit_r = ⌈1.5·log₂r⌉
for (lo_n = 0; lo_n < k; lo_n++)
{
    fmpz_mul_ui(t, entry(r, lo_n), sqrt(f->length));
    if (fmpz_bits(t) > bound) break;         // 超出阈值则停止
}

// 3. 高次端列过滤（hi_n）类似

return lo_n + hi_n;  // 可用列总数
```

**列过滤本质**：阈值 `bound = log₂(p^a) - 1.5·log₂r` 对应于：

```
界 × √N ≤ p^a / 2^(1.5·r)
```

即列的 CLD 系数界必须在当前 Hensel 精度与格短向量阈值的"信噪比"可接受范围内（与 P1a 文档 §2.2 方案 B 一致）。

### 4.3 精度需求的变化：P1a 的影响

P1a（van Hoeij LLL）实现后，控制流变为：
```
__factor_squarefree_primitive_ZZ:
    → __hensel_lift(f, factors_Zp, p)  ← 提升到 Mignotte 精度
    → __factor_recombine(f, lifted, m) ← P1a：LLL 重组
```

P1b 修改为：
```
    → __linear_hensel_lift_with_lll(f, factors_Zp, p)
          内部交织：线性提升 + LLL 早期检测
          早期终止后不再调用外部 __factor_recombine
```

或等价地，将控制流移入 `__factor_squarefree_primitive_ZZ`，使其管理提升-检测的交织循环。

---

## 5. CLPoly 实现方案

### 5.1 改动范围

**主要改动**：`polynomial_factorize_univar.hh`

| 函数 | 改动类型 | 说明 |
|------|---------|------|
| `__hensel_step` | **废弃/改名** → `__hensel_step_quadratic` | 保留二次步骤备用 |
| `__hensel_step_linear` | **新增** | 线性单步（§5.2） |
| `__hensel_lift` | **保留，新增重载** | 原二次提升保留；新增线性变体 |
| `__linear_hensel_lift_with_lll` | **新增** | P1b 主函数（§5.3） |
| `__factor_squarefree_primitive_ZZ` | **修改** | 调用新 P1b 主函数 |

**不改动**：
- P1a 的所有模块（M1-M5，`__vanhoeij_recombine` 等）——P1b 通过接口调用
- `__zassenhaus_recombine`——保留作为 fallback
- 多变量因式分解部分（`polynomial_factorize_wang.hh` 等）

### 5.2 M1：`__hensel_step_linear`（线性单步）

**函数规约**

```
函数名称：__hensel_step_linear

功能描述：执行一步线性 Hensel 提升，将 g·h ≡ f (mod p^a) 升至 (mod p^(a+1))。
          Bézout 系数 s, t 保持 mod p（不随 a 更新）。

前置条件（Requires）：
  - f ∈ Z[x]，g, h ∈ Z[x]（系数 mod p^a，对称约化）
  - s·g + t·h ≡ 1  (mod p)，deg s < deg h，deg t < deg g
  - f ≡ g·h  (mod p^a)
  - p_a = p^a（作为 ZZ）

后置条件（Ensures）：
  - g 原地更新为 g'，h 原地更新为 h'
  - g'·h' ≡ f  (mod p^(a+1))
  - s, t 不变（仍为 mod p 的 Bézout 系数）

副作用：原地修改 node.g, node.h（不修改 node.s, node.t）
```

**函数签名**

```cpp
// 参数：
//   node  — 包含 g, h（mod p^a），s, t（固定 mod p，不更新）
//   f     — 原多项式（ZZ 系数）
//   p     — 素数（小整数）
//   p_a   — p^a（当前精度，用于计算 e = (f - g*h)/p^a）
void __hensel_step_linear(
    __hensel_node&           node,
    const upolynomial_<ZZ>&  f,
    uint32_t                 p,          // 素数（mod p 运算用）
    const ZZ&                p_a);       // p^a（当前模数）
```

**实现步骤（伪代码）**

```cpp
// 步骤 1：计算误差项 e = (f - g*h) / p^a  (mod p)
upolynomial_<ZZ> gh = node.g * node.h;
upolynomial_<ZZ> err_big = f - gh;           // f - g*h (在 Z 中)
// 每个系数精确整除 p^a，然后取 mod p
for (auto& term : err_big)
    term.second = (term.second / p_a) % p;   // 整数：先除 p^a，再取 mod p
// 去零项
upolynomial_<ZZ> e = normalize_zero(err_big);

// 步骤 2：计算 σ = s*e mod h  (mod p)
// 先将 e, s, h 转换到 Z/p 表示，做带余除法
upolynomial_<Zp> e_p  = ZZ_to_Zp(e, p);
upolynomial_<Zp> s_p  = ZZ_to_Zp(node.s, p);  // s 本身已是 mod p
upolynomial_<Zp> h_p  = ZZ_to_Zp(node.h, p);  // h mod p
upolynomial_<Zp> se_p = s_p * e_p;
upolynomial_<Zp> q_p, sigma_p;
divmod(q_p, sigma_p, se_p, h_p);               // sigma = se mod h (mod p)

// 步骤 3：计算 τ = t*e + q*g  (mod p)
upolynomial_<Zp> t_p   = ZZ_to_Zp(node.t, p);
upolynomial_<Zp> g_p   = ZZ_to_Zp(node.g, p);
upolynomial_<Zp> te_p  = t_p * e_p;
upolynomial_<Zp> qg_p  = q_p * g_p;
upolynomial_<Zp> tau_p = te_p + qg_p;         // 结果已 mod p（Zp 运算）
// 不需要 mod h（tau 对应 g 这边，次数约束由 deg t < deg g 保证）

// 步骤 4：更新 g, h
// g' = g + p^a * tau  (mod p^(a+1))
upolynomial_<ZZ> tau_ZZ = Zp_to_ZZ(tau_p);
for (auto& term : tau_ZZ)
    term.second *= p_a;                        // 乘以 p^a
node.g = node.g + tau_ZZ;
ZZ p_a1 = p_a * ZZ(p);                        // p^(a+1)
__upoly_mod_coeff(node.g, p_a1);              // 系数 mod p^(a+1)，保对称

// h' = h + p^a * sigma  (mod p^(a+1))
upolynomial_<ZZ> sigma_ZZ = Zp_to_ZZ(sigma_p);
for (auto& term : sigma_ZZ)
    term.second *= p_a;
node.h = node.h + sigma_ZZ;
__upoly_mod_coeff(node.h, p_a1);
```

**注**：`ZZ_to_Zp` 将 ZZ 系数多项式转换到 Z/p（类似现有 `__upoly_Zp_to_ZZ` 的逆操作），可通过 `Zp(coeff % p, p)` 实现。`Zp_to_ZZ` 则取 `static_cast<int64_t>(coeff.number())`（与现有 `__upoly_Zp_to_ZZ` 相同）。

**代价**：每步主要代价为两次 mod p 多项式乘法（`s_p * e_p` 和 `t_p * e_p`）+ 一次带余除法（`se_p / h_p`）。所有运算在 Z/p 上（单精度整数），代价 O(n²) 在 Z/p 中（对比二次提升：O(n²) 在 Z/p^a 中，多精度算术代价更高）。

### 5.3 M3：`__linear_hensel_lift_with_lll`（主控函数）

**函数规约**

```
函数名称：__linear_hensel_lift_with_lll

功能描述：对 f ∈ Z[x] 的模因子列表，通过线性 Hensel 提升（精度逐步增加）
          与 van Hoeij LLL（P1a）交织，实现早期因子检测。
          一旦 LLL 找到所有因子立即返回，无需提升到完整 Mignotte 界。

前置条件（Requires）：
  - f 本原，lc(f) > 0，deg(f) ≥ 2
  - factors 是 f mod p 的模因子（互素，首一）
  - |factors| ≥ 2

后置条件（Ensures）：
  - 返回 f 的不可约因子列表（ZZ 系数，本原，lc > 0）
  - 等价于先做完整提升再调用 __factor_recombine，但通常更早完成

副作用：无（不修改输入）
```

**函数签名**

```cpp
std::vector<upolynomial_<ZZ>>
__linear_hensel_lift_with_lll(
    const upolynomial_<ZZ>&        f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t                       p);
```

**控制流伪代码**

```
输入：f, [h₁,...,hᵣ] mod p, p

1. 预处理：
   - lc_f = lc(f)；对 h₁ 乘以 lc_f（处理首项系数分配，与现有 __hensel_lift 相同）
   - 构建二叉树（__hensel_tree_build）
   - 计算初始精度 a₀（FLINT 启发式或固定值 a₀ = 10）
   - 计算 Mignotte 精度 a_mig = ⌈log_p(2·lc_f·B_Mig)⌉

2. 初始线性提升：将所有叶子因子从 mod p 提升到 mod p^a₀
   for step = 1 to a₀:
       __hensel_step_linear_recursive(nodes, f, p, current_p^a)
       current_p^a *= p

3. 初始化 van Hoeij LLL 状态（P1a 的 M5 起始状态）：
   r = 叶子数（模因子数）
   U_exp = bit_length(max(r, 20))
   B = (r+1) · 2^(2·U_exp)
   M = 2^U_exp · I_r
   J_cur = 0
   J_target = J₀（= 30 if 3r > N+1 else 10）
   active = {0,...,r-1}
   f_star = f
   result = []

4. 主交织循环：
   while active.size() > 1:
       4a. [M1 P1a] cld = __cld_polys(f_star, active_factors mod p^a, p^a)
       4b. [M2 P1a] J_new = __build_cld_matrix(M, cld, J_cur, J_target, p^a)
       4c. [M3 P1a] short_rows = __lll_reduce(M, U, B)
       4d. [M4 P1a] candidates = __extract_candidates(short_rows, U, r')
       4e. 验证并提取因子（与 P1a __vanhoeij_recombine 相同逻辑）
           对每个候选子集 S：
               g_trial = lc(f_star) · ∏_{i∈S} hᵢ mod p^a，对称约化，本原化
               若 g_trial | f_star（试除成功）：
                   result.push_back(g_trial)
                   f_star /= g_trial
                   active -= S
                   重置 M, J_cur, J_target
                   break → 继续主循环（早期因子找到）
       4f. 若本轮无进展：
               若 J_target 可增大（< J_max）：
                   J_target *= 2
               否则（列已耗尽）：
                   继续提升 k 步（增量）：
                   for step = 1 to k:
                       __hensel_step_linear_recursive(nodes, f, p, current_p^a)
                       current_p^a *= p
                       a += 1
                   若 a > a_mig：  // 超出 Mignotte 界（不应发生，安全网）
                       对剩余 active 因子回退 __zassenhaus_recombine
                       return result + zassenhaus_result

5. 若 active.size() == 1：
   result.push_back(f_star)   // 最后一个因子

return result
```

### 5.4 M2：`__hensel_step_linear_recursive`（树状线性提升）

类比现有 `__hensel_lift_recursive`（二次），对整棵提升树做一步线性提升：

```cpp
void __hensel_step_linear_recursive(
    std::vector<__hensel_node>& nodes,
    int idx,
    const upolynomial_<ZZ>& target,
    const ZZ& p_a,       // 当前模数 p^a
    uint32_t p)           // 素数（mod p 用）
{
    // 对当前节点做一步线性提升（target ≡ g·h mod p^a → mod p^(a+1)）
    __hensel_step_linear(nodes[idx], target, p, p_a);

    // 递归左右子树（以新的 g, h 作为子节点的 target）
    if (nodes[idx].left != -1)
        __hensel_step_linear_recursive(nodes, nodes[idx].left,
                                       nodes[idx].g, p_a, p);
    if (nodes[idx].right != -1)
        __hensel_step_linear_recursive(nodes, nodes[idx].right,
                                       nodes[idx].h, p_a, p);
}
```

### 5.5 接口改动：`__factor_squarefree_primitive_ZZ`

```cpp
// 原代码（P1a 后）：
auto [lifted, modulus] = __hensel_lift(f, sel.factors, sel.prime);  // 二次提升到全精度
return __factor_recombine(f, lifted, modulus);                       // LLL 重组

// P1b 修改后：
return __linear_hensel_lift_with_lll(f, sel.factors, sel.prime);    // 线性提升+早期检测
// （__linear_hensel_lift_with_lll 内部完成提升+LLL，直接返回因子列表）
```

原 `__hensel_lift` 和 `__factor_recombine` 的接口不变，作为 fallback 保留（当 P1b 内部安全网触发时可调用）。

### 5.6 模块划分（M0-M4，n 因子直接方案，与 P1a 风格一致）

> 架构评估后更新：采用 n 因子直接方法，原 M1（2 因子单步）和 M2（树递归）合并为 M0（Bézout 链）+ M1（n 因子步）。

| 模块 | 函数名 | 行数估计 | 说明 |
|------|--------|--------|------|
| **M0** | `__linear_bezout_chain` | ~30 行 | Bézout 链预计算；一次性，固定 mod p |
| **M1** | `__hensel_step_linear_nfactor` | ~50 行 | n 因子直接线性步；无树递归 |
| **M3** | `__linear_hensel_lift_with_lll` | ~150 行 | 主控循环（提升+LLL交织+早期终止） |
| **M4** | `__heuristic_starting_precision` | ~15 行 | 初始精度 a₀ 计算（FLINT 启发式） |
| 接口 | `__factor_squarefree_primitive_ZZ` 修改 | ~2 行 | 调用改向 M3 |

**复用 P1a 模块**（M3 内部直接调用）：
- `__cld_polys`（P1a M1）
- `__build_cld_matrix`（P1a M2）
- `__lll_reduce`（P1a M3）
- `__extract_candidates`（P1a M4）

**[TODO v2] Monagan 2022 矩阵优化**：将 M1 的逐步修正序列组织为矩阵乘法，总代价从 O(a_stop · r · n²) 降至 O(r · n · a_stop²)。与早期终止叠加后收益有限，留作后续。

### 5.7 测试策略

1. **回归测试**：所有现有 `test_factorize` 用例（256 个）全部通过
   - 线性提升结果与二次提升结果应相同（因子集合一致）
   - 用 `make test/test_factorize` 验证

2. **性能对比**：benchmark uni-70
   - 目标：`(x-1)(x-2)···(x-70)` < 10 ms（P1a+P1b 合计）
   - 当前 P1a 后预计约 30-40 ms（提升阶段仍是 83ms 的主要部分）
   - P1b 后目标 < 10 ms

3. **早期终止验证**：
   - 打印提升步数（DEBUG 模式）：验证 uni-70 确实在 < 20 步时停止（而非 Mignotte 界 ~60 步）
   - 压力测试 `make stress` 全部通过（release 模式）

4. **正确性单元测试**（新增）：
   - `test_hensel_step_linear`：对小例子验证单步线性提升的正确性
   - 比较 `g * h mod p^(a+1)` 与 `f mod p^(a+1)` 是否相等

---

## 6. 复杂度分析

### 6.1 当前二次提升复杂度

设 `n = deg(f)`，`D = log₂(m_target) ≈ Mignotte bit 数`，`r = 模因子数`。

**二次提升（当前 `__hensel_lift`）**：

- 提升步数：`k = ⌈log₂ D⌉`（每步精度平方）
- 第 j 步精度：`d_j = 2^j · log₂ p` bit
- 第 j 步每个节点代价：`O(n · d_j)` bit 操作（多精度多项式乘法）
- 树中节点数：`O(r)` 个
- 总代价：`O(r · n · Σ d_j) = O(r · n · D)` bit 操作

对 `n = D = 300`（uni-70 近似）：`O(r · 300²) = O(r · 90000)` bit 操作。

**Bézout 更新代价**（二次提升每步都需要）：

- 第 j 步：O(n² · d_j) bit 操作（Bézout 更新 = 多项式乘法×3 + 除法×2 in Z/p^(2d_j)）
- 总计：O(n² · D) bit 操作（占总代价的主要部分）

### 6.2 线性提升复杂度

**线性提升（P1b `__hensel_step_linear`）**：

- 提升步数：`a_stop`（早期终止精度，实际约 10–20 步对 uni-70）
- 每步代价：O(n²) in Z/p（单精度，小整数）
  - `s_p * e_p`：O(n²) mod p 操作
  - `t_p * e_p`：O(n²) mod p 操作
  - 带余除法：O(n²) mod p 操作
  - 修正量乘以 p^a：O(n) 次大整数乘法（当前精度 a·log p bit）
- 每步总代价：O(n² + n·a·log p) bit 操作
- **无 Bézout 更新代价**（s, t 固定不变）

**总代价（停在 a_stop）**：`O(a_stop · (n² + n·a_stop·log p))` bit 操作。

| 参数 | uni-70 实际值 |
|------|-------------|
| n = deg(f) | 70 |
| D = Mignotte 精度 (bit) | ≈ 300 |
| p (典型取值) | ≈ 10³ 量级 |
| a_mig = D / log₂p | ≈ 30 步 |
| a_stop (早期终止) | ≈ 10-15 步（LLL 在低精度收敛） |

**对比**（uni-70）：

| 阶段 | 二次提升（当前） | 线性提升+早期终止（P1b）| 加速比 |
|------|--------------|----------------------|-------|
| 提升步数 | 6（二次步数） | 10–15（线性步数） | 步数更多 |
| 每步代价 | O(n²·D) bit（大整数） | O(n²) bit（小整数 mod p） | 每步约 300/10 ≈ 30x 更快 |
| Bézout 更新 | 有（每步 O(n²·D)） | 无 | 完全消除 |
| 精度上限 | 300 bit（Mignotte） | 150 bit（早期终止） | 2x 精度减少 |
| **总代价估计** | ~6 × O(n²·D) | ~12 × O(n²) | 约 **10-20x** |

### 6.3 Monagan 2019 的复杂度结果

**经典线性提升（GCL §6.5）**：每步 Bézout 更新代价 O(n²) in Z/p^a，总计 O(n²·D) = O(n²d²)（其中 d = D/log p）

**Monagan 2019（立方代价 LHL）**：通过矩阵方法将全部 D 步线性提升的 **总 Bézout 更新** 代价降至 O(n·D²) = O(nd²)——不是 O(n²d²)。关键技术：
- 将 r 步的修正序列组织成多项式矩阵乘法
- 利用"移位多项式"（shifted polynomial）在 Zp 上高效压缩计算
- 使得 n 个因子的立方代价为 O(d·n·d) = O(nd²) bit 操作（d = 精度 bit 数）

**对 CLPoly 初版（P1b）的含义**：CLPoly 初版不实现 Monagan 的矩阵优化，采用逐步单步线性提升。其总代价为 O(a_stop · n²)（经典线性），比 Monagan 2019 多一个 n 因子——但由于早期终止大幅减少 a_stop，实际性能已有显著改善。Monagan 的矩阵优化作为 P1b 的后续增强项（v2）。

---

## 7. 实现注意事项

### 7.1 初始 Bézout 系数的获取

线性提升要求 s, t 满足 `s·g + t·h ≡ 1 (mod p)`（不是 mod p^a）。

**现有 `__hensel_tree_build_recursive` 已计算 s, t**（通过 `polynomial_GCD(g_zp, h_zp, s_zp, t_zp)` 扩展 GCD），且结果是 Z/p 下的 Bézout 系数。

**区别**：
- 当前二次提升：调用 `__hensel_step` 时，s, t 已经是 mod p（初始），但 `__hensel_step` 内部会 **将 s, t 升级到 mod m²**（Bézout 提升部分，polynomial_factorize_univar.hh:414–455）
- P1b 线性提升：`__hensel_step_linear` **不更新 s, t**，它们始终保持 mod p 的初始值

因此，`__hensel_node` 结构不需要改变，但 `__hensel_step_linear` 中需要跳过 `node.s`、`node.t` 的更新步骤。

### 7.2 ZZ_to_Zp 转换

在 `__hensel_step_linear` 中，需要将 g, h（mod p^a，大整数系数）转换到 Z/p（单精度）：

```cpp
// 将 ZZ 多项式的系数 mod p，得到 Zp 多项式
upolynomial_<Zp> ZZ_to_Zp_poly(const upolynomial_<ZZ>& f, uint32_t p) {
    Zp zero(0, p);
    upolynomial_<Zp> result;
    for (auto& term : f) {
        ZZ c = term.second;
        // 对称约化的系数：先对 p 取正余数
        ZZ::fdiv_r(c, c, ZZ(p));  // c in [0, p)
        Zp c_p(static_cast<uint64_t>(c._val >= 0 ? c._val : c._val + p), p);
        if (c_p != zero)
            result.push_back({term.first, c_p});
    }
    return result;
}
```

注：`g mod p` 的含义是取 g 的系数模 p（非对称模，取 [0,p)），因为 Bézout 等式是在 F_p 上的。

### 7.3 树状提升与叶子因子提取

线性提升完成后，提取叶子因子的方式与现有 `__hensel_extract_factors` 相同（不变）。

不同之处：提升结束时 `m = p^(a_stop)`（早期终止精度），而非 `m = p^(a_mig)`（Mignotte 精度）。`__linear_hensel_lift_with_lll` 内部跟踪当前精度 `p_a`，用于 CLD 矩阵构建和列过滤。

### 7.4 对称约化策略

当前 `__upoly_mod_coeff` 对系数取 [0, m) 区间（非对称），`__upoly_symmetric_mod` 才做对称约化到 (-m/2, m/2]。

在 `__hensel_step_linear` 中：
- g', h' 的系数应对 `p^(a+1)` 做非对称模（[0, p^(a+1))），与现有 `__hensel_step` 的 `__upoly_mod_coeff(node.g, m2)` 一致
- 但在传递给 CLD 计算时，应先做对称约化（`__upoly_symmetric_mod`），与 P1a 的期望格式一致

### 7.5 增量提升的步长 k

主循环中"继续提升 k 步"的 k 值选取：

| 策略 | 说明 | 推荐 |
|------|------|------|
| k = 1 | 每次提升 1 步，最细粒度 | 最灵活，但循环开销大 |
| k = 5 | 每次提升 5 步，然后检测 | 初版推荐 |
| k = a_current / 2 | 自适应（倍增倍减） | 最接近 FLINT 的外层循环 |

**初版推荐 k = 5**：每提升 5 步检测一次，平衡提升代价与检测频率。

---

## 8. P1a / P1b 接口总结

**P1a 完成后**（当前状态）：`__factor_recombine` → van Hoeij LLL 重组，接口不变

**P1b 完成后**（目标状态）：

```
__factor_squarefree_primitive_ZZ:
    [P1b] __linear_hensel_lift_with_lll(f, sel.factors, sel.prime)
            内部调用：
            [P1b M1] __hensel_step_linear（线性单步）
            [P1a M1] __cld_polys（CLD 多项式）
            [P1a M2] __build_cld_matrix（格矩阵构建）
            [P1a M3] __lll_reduce（LLL 规约）
            [P1a M4] __extract_candidates（候选提取）
            [P1a M5逻辑] 因子验证 + 早期终止
```

P1a 的 `__vanhoeij_recombine`（M5）被 P1b 的 `__linear_hensel_lift_with_lll` 吸收，不再独立调用。P1b 完成后，顶层 `__factor_recombine` 的调用（在 `__factor_squarefree_primitive_ZZ` 中）改为调用 `__linear_hensel_lift_with_lll`。

**P1a 与 P1b 互不干扰**（初版独立实现）：P1b 实现前，P1a 可独立运行（先完整提升到 Mignotte 精度，再 LLL）。P1b 实现后，两者通过 `__linear_hensel_lift_with_lll` 的交织架构合并。

---

## 9. 参考文献

| 文献 | 内容 | 与 P1b 的关系 |
|------|------|-------------|
| **Monagan, M. (2019). Linear Hensel Lifting for F_p[x,y] and Z[x] with Cubic Cost. ISSAC 2019, pp. 299–306.** | **P1b 直接参考**：线性提升立方代价算法，Bézout 保持 mod p，矩阵优化 | 核心算法文献 |
| **Monagan, M. & Tuncer, B. (2019/2020). Polynomial Factorization in Maple 2019. Proc. CASC 2019.** | Maple 2019 算法架构：线性提升 + 早期检测 + van Hoeij LLL 交织 | P1b 架构参考 |
| Monagan, M. (2022). Linear Hensel Lifting for Z_p[x,y] for n Factors with Cubic Cost. ISSAC 2022. | n 因子直接线性提升（P1b v2 参考） | 后续增强 |
| FLINT `src/fmpz_poly_factor/factor_van_hoeij.c` | 外层精度倍增循环、初始精度启发式、LLL 集成 | 实现参考 |
| FLINT `src/fmpz_poly/hensel_lift_without_inverse.c` | 线性提升单步公式（`G' = g + p·correction` 模式） | 实现参考 |
| FLINT `src/fmpz_poly_factor/CLD_mat.c` | CLD 矩阵构建与列过滤 | P1a/P1b 共用 |
| CLPoly `clpoly/polynomial_factorize_univar.hh` | 现有 `__hensel_step`（二次）、`__hensel_lift`、`__hensel_tree_build` | 基础设施 |
| docs/design/vanhoeij/architecture.md | P1a 架构（M1-M5 定义） | 接口约定 |
| docs/research/vanhoeij-lll-research.md §5.3-5.4 | P1b 定位分析 | 背景 |

---

## 10. 结论与下一步

### 结论

P1b 线性 Hensel 提升的核心价值来自**两处叠加效应**：

1. **每步代价降低**：线性步中 Bézout 系数固定 mod p，每步仅需 Z/p 上的多项式运算（单精度），不需要 mod p^a 的大整数 Bézout 更新
2. **早期终止**：交织 LLL 后，许多情形无需提升到完整 Mignotte 精度即可完成分解

对 uni-70 案例，理论加速约 10-20x（P1a 基础上），预期将 P1a 后的瓶颈（提升阶段，约 40 ms）降至 2-5 ms 量级。

### 实施顺序（n 因子直接方案）

1. 确认 P1a 完整通过（test_factorize 256 例，bench_comparative 对比）
2. 实现 P1b M4（`__heuristic_starting_precision`）并验证
3. 实现 P1b M0（`__linear_bezout_chain`）并单元测试（验证 Σ sᵢ·∏_{j≠i}hⱼ ≡ 1）
4. 实现 P1b M1（`__hensel_step_linear_nfactor`）并单元测试（验证单步提升正确性）
5. 实现 P1b M3（`__linear_hensel_lift_with_lll`），先不含早期终止
6. 验证 M3 正确性（与原 P1a 结果相同，256 例）
7. 加入早期终止逻辑（交织 LLL）
8. 性能测试：uni-70 及 stress 测试

### 实现范围估计（n 因子直接方案）

| 模块 | 行数 | 难度 |
|------|------|------|
| `__heuristic_starting_precision` (M4) | ~15 行 | 低 |
| `__linear_bezout_chain` (M0) | ~30 行 | 低 |
| `__hensel_step_linear_nfactor` (M1) | ~50 行 | 低 |
| `__linear_hensel_lift_with_lll` (M3) | ~150 行 | 中 |
| 接口修改 | ~2 行 | 低 |
| 总计 | **~247 行** | 中 |

预期改善（P1a + P1b 合计）：uni-70 80 ms → < 10 ms，与 FLINT 4.5 ms 差距缩小至 <3x。
