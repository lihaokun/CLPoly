# CLPoly Hensel 提升算法设计与分析

> 本文档分析 CLPoly 当前的 Hensel 提升实现，讨论其优缺点，
> 并与其他库的实现进行对比，为后续优化和多变量扩展提供参考。

---

## 1. 当前实现概述

### 1.1 位置与入口

所有 Hensel 提升代码位于 `clpoly/polynomial_factorize.hh`，核心函数：

| 函数 | 行号 | 职责 |
|---|---|---|
| `__hensel_node` | L696-705 | 二叉树节点结构体 |
| `__hensel_tree_build` | L777-785 | 构建初始提升树（mod p） |
| `__hensel_step` | L789-882 | 单步二次提升：mod m → mod m² |
| `__hensel_lift_recursive` | L905-917 | 自顶向下递归提升 |
| `__hensel_lift` | L920-959 | M2 入口：确定精度 + 提升 + 提取因子 |

### 1.2 算法结构

采用**二叉树 + 二次提升**方案：

```
r 个因子 f₁, ..., fᵣ (mod p)
        │
        ▼
构建平衡二叉树 (r-1 个内部节点)
每个节点: g = 左子积, h = 右子积, s·g + t·h ≡ 1
        │
        ▼
循环: m ← p; while m ≤ target: 递归提升, m ← m²
        │
        ▼
提取叶子因子 (mod m)
```

### 1.3 数据结构

```cpp
struct __hensel_node {
    upolynomial_<ZZ> g;          // 左子树因子之积
    upolynomial_<ZZ> h;          // 右子树因子之积
    upolynomial_<ZZ> s;          // Bézout 系数: s·g + t·h ≡ 1 (mod m)
    upolynomial_<ZZ> t;          // Bézout 系数
    int left;                    // 左子节点索引 (-1 = 叶子)
    int right;                   // 右子节点索引 (-1 = 叶子)
    int leaf_start, leaf_end;    // 对应的叶子范围
};
```

### 1.4 二次提升步骤 (`__hensel_step`)

给定 f ≡ g·h (mod m)，s·g + t·h ≡ 1 (mod m)，提升到 mod m²：

**第一部分：提升因子**

```
1. e ← (f - g·h) / m  (精确整除)
   e ← e mod m
2. se ← s·e
   (q, r) ← divmod(se, h) in Z_m[x]
3. τ ← t·e + q·g  (mod m)
   g ← g + m·τ     (mod m²)
   h ← h + m·r     (mod m²)
```

**第二部分：提升 Bézout 系数**

```
4. e' ← (1 - s·g - t·h) / m  (精确整除)
   e' ← e' mod m
5. (q', r') ← divmod(s·e', h) in Z_m[x]
   s ← s + m·r'               (mod m²)
   t ← t + m·(t·e' + q'·g)   (mod m²)
```

### 1.5 提升精度确定

使用 Mignotte 界：

```
B = C(n, ⌊n/2⌋) · ‖f‖₂
target = 2 · |lc(f)| · B
提升直到 m > target
```

### 1.6 首项系数处理

提升前将 lc(f) 分配到 factors[0] 的所有系数：

```cpp
Zp lc_mod_p(f.front().second, p);
for (auto& term : factors_adj[0])
    term.second *= lc_mod_p;
```

---

## 2. 性能分析

### 2.1 复杂度

设 f 度数为 n，r 个因子，提升 k 步（二次提升，实际步数 ≈ log₂(k)）：

| 操作 | 每步代价 | 总代价 |
|---|---|---|
| `f - g·h`（多项式乘法） | O(n²) | O(n² log k) |
| `s·e`（多项式乘法） | O(n²) | O(n² log k) |
| `divmod(se, h)` in Z_m[x] | O(n²) | O(n² log k) |
| Bézout 更新 | O(n²) | O(n² log k) |

系数运算使用 GMP 大整数，系数增长到 O(k log p) 位后，
乘法代价为 O(M(k log p))，其中 M(b) 为 b 位乘法代价。

### 2.2 当前实现的问题

#### 问题 1：`__hensel_step` 中大量临时对象

```cpp
// 当前: 每步创建多个临时多项式
upolynomial_<ZZ> gh = node.g * node.h;      // 临时 1
upolynomial_<ZZ> e = f - gh;                 // 临时 2
upolynomial_<ZZ> se = node.s * e;            // 临时 3
...
upolynomial_<ZZ> sg = node.s * node.g;       // 临时 6
upolynomial_<ZZ> th = node.t * node.h;       // 临时 7
```

每步 `__hensel_step` 创建约 **12 个临时多项式对象**，涉及大量内存分配。
对比 FLINT 的原地操作方式，这是主要性能差距来源。

#### 问题 2：`__upoly_divmod_mod` 的 O(n²) 实现

```cpp
// 当前: 长除法中每步重建 new_r
upolynomial_<ZZ> new_r;
new_r.reserve(r.size());
// ... 归并式逐项计算 ...
r.data() = std::move(new_r.data());
```

每次消去一项都创建新的 `new_r` 向量。可以改为原地减法。

#### 问题 3：缺少提前终止

当前实现固定提升到 `m > target`。FLINT 和 Singular 都有
**自适应提升界**和**提前因子检测**：

- FLINT: `_hlift_quartic` 在每步后检查是否已有因子可提取
- Singular: `earlyFactorDetect()` + `liftBoundAdaption()`

#### 问题 4：`e / m` 使用 `fdiv_q` 逐系数除

```cpp
for (auto& term : e.data())
{
    ZZ::fdiv_q(term.second, term.second, m);
    ZZ::fdiv_r(term.second, term.second, m);
}
```

这里对每个系数做两次 GMP 除法。可以合并为 `fdiv_qr` 一次完成。

---

## 3. 其他库的实现对比

### 3.1 FLINT

```
结构:   扁平数组 (不用二叉树)
提升:   二次提升
特点:
  - 按因子数选择算法: r=2 用 _hlift_quartic2, r<20 用 _hlift_quartic, r≥20 用 _hlift_quintic
  - 使用偏分式结构 fmpz_mpoly_pfrac_t 管理 Bézout 系数
  - 原地操作，极少临时分配
  - 多变量提升: fmpz_mpoly_hlift 逐变量线性提升
```

### 3.2 NTL

```
结构:   多因子同时提升 (MultiLift)
提升:   二次提升
特点:
  - MulLift1 处理 2 因子的基础情况
  - MultiLift 递归分组处理 r>2 因子
  - 使用 zz_pX 做 Zp 运算，ZZX 做 ZZ 运算
  - 有 IsStableDiv 检查提升是否稳定
```

### 3.3 Singular/Factory

```
结构:   逐变量提升 (henselLift23, henselLift)
提升:   线性提升为主
特点:
  - henselStep12 处理单步度数增量
  - 逐对 extgcd 累积 Bézout 系数
  - 自适应提升界 liftBoundAdaption
  - 提前因子检测 earlyFactorDetect
```

### 3.4 对比总结

| | CLPoly | FLINT | NTL | Singular |
|---|---|---|---|---|
| **树结构** | 平衡二叉树 | 扁平数组 | 递归分组 | 逐变量 |
| **提升方式** | 二次 | 二次 | 二次 | 线性/二次 |
| **Bézout 管理** | 树节点 s,t | 偏分式结构 | MultiLift 递归 | 逐对 XGCD |
| **内存管理** | 大量临时对象 | 原地操作 | 中等 | 中等 |
| **提前终止** | 无 | 有 | 有 | 有 |
| **r>2 优化** | 统一二叉树 | 按 r 选算法 | 统一递归 | 统一 |
| **多变量支持** | 仅单变量 | 逐变量线性 | 无 | 逐变量线性 |

---

## 4. 改进方向

### 4.1 近期优化（不改变算法结构）

#### 4.1.1 减少临时对象

将 `__hensel_step` 中的临时多项式改为函数级缓冲区，在多步提升间复用：

```cpp
struct __hensel_workspace {
    upolynomial_<ZZ> e, se, tau, ep, sep;  // 可复用的缓冲区
};
```

#### 4.1.2 合并系数操作

将 `fdiv_q` + `fdiv_r` 合并为单次 `fdiv_qr`：

```cpp
// 之前:
ZZ::fdiv_q(term.second, term.second, m);
ZZ::fdiv_r(term.second, term.second, m);

// 之后:
ZZ q_val;
ZZ::fdiv_qr(q_val, term.second, term.second, m);
term.second = q_val;  // 或根据后续使用选择
```

#### 4.1.3 `__upoly_divmod_mod` 原地化

避免每步长除法重建 `new_r`，改为原地减法 + 尾部清理。

### 4.2 中期改进

#### 4.2.1 提前因子检测

在提升循环中，每隔若干步尝试对称约化 + 试除，提前提取已确定的因子：

```
while m ≤ target:
    __hensel_lift_recursive(...)
    m ← m²
    if m > √target:   // 精度已过半
        尝试提取 degree = 1 的因子 (线性因子通常最先收敛)
```

#### 4.2.2 自适应提升界

对每个因子单独跟踪其系数是否已稳定（连续两步不变），
稳定的因子可以提前提取并从树中移除。

### 4.3 远期：多变量 Hensel 提升

单变量提升扩展到多变量时，需要支持逐变量线性提升（详见
`multivariate-factorization-design.md` §6）。

核心区别：

| | 单变量提升 | 多变量提升 |
|---|---|---|
| **提升变量** | 精度 p → p^k | 变量 xₖ - αₖ 的阶 |
| **提升方式** | 二次（精度每步翻倍） | 线性（每步增加一阶） |
| **系数域** | Z | Z[x₁, ..., xₖ₋₁] |
| **Bézout 系数** | upolynomial_<ZZ> | polynomial_<ZZ, lex> |
| **误差计算** | (f - g·h) / m | coeff(f - ∏gᵢ, (xₖ-αₖ)^j) |

需要新增的函数：

```cpp
// 多变量线性 Hensel 提升（逐变量）
// 将 f(x₁, αₖ) 的因子提升为 f(x₁, xₖ) 的因子
template<class var_order>
void __multivar_hensel_step(
    std::vector<polynomial_<ZZ, lex_<var_order>>>& factors,
    const polynomial_<ZZ, lex_<var_order>>& f,
    const variable& lift_var,
    const ZZ& alpha,
    int64_t max_degree);
```

这需要将单变量 Hensel 的核心逻辑（Bézout 系数管理、误差分解）
抽象为可处理不同"提升维度"的通用框架。

### 4.4 终极：稀疏 Hensel 提升 (MTSHL)

MTSHL 用稀疏插值替代经典 MDP，需要：

1. 二变量 Hensel 提升 (BHL) 作为子程序
2. 稀疏插值模块（Ben-Or/Tiwari 或 Zippel）
3. Vandermonde 求解器

详见 `multivariate-factorization-design.md` §11。

---

## 5. 总结

当前实现的 Hensel 提升**算法正确**，所有测试通过。主要改进空间在于：

1. **性能**：减少临时对象、原地操作、合并 GMP 调用（近期，不改结构）
2. **提前终止**：自适应提升界 + 提前因子检测（中期）
3. **多变量扩展**：逐变量线性提升（Phase 5 需要）
4. **MTSHL**：稀疏插值驱动的提升（Phase 8 终极目标）

---

## 附录 A: 参考文献

- **Zassenhaus**: "On Hensel Factorization I", J. Number Theory 1969
- **GCL**: Geddes, Czapor & Labahn, "Algorithms for Computer Algebra", Kluwer 1992 (§6.3, §15.5)
- **MCA**: von zur Gathen & Gerhard, "Modern Computer Algebra", Cambridge 2013 (§15)
- **NTL**: Shoup, "A New Polynomial Factorization Algorithm and its Implementation", 1995
- **FLINT**: Hart et al., "FLINT: Fast Library for Number Theory", https://flintlib.org
- **MTSHL**: Monagan & Tuncer, "Using Sparse Interpolation in Hensel Lifting", CASC 2016
