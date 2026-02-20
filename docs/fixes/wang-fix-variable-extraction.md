# 修正方案：变量幂次提取预处理 + 不可约启发式修正

> **状态：已完成。** `__extract_monomial_content` 已实现，不可约启发式改为单次求值证明。
> 本文档保留为历史记录，权威描述见 [multivariate-factorization-design.md](../design/multivariate-factorization-design.md)。
>
> 对应 workflow.md §5.1 修正方案文档

---

## 第一部分：复现与定位

### Bug 1：变量幂次提取缺失

**最小复现用例**

```cpp
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;
PolyZZ make_lex(const polynomial_ZZ& p) { PolyZZ r; poly_convert(p, r); return r; }

int main() {
    variable x("x"), y("y"), z("z");

    auto g = make_lex(
        polynomial_ZZ(-ZZ(2)*x*y + ZZ(2)*y*z - ZZ(2)*y + ZZ(3)*z) *
        polynomial_ZZ(ZZ(2)*x*x - ZZ(2)*x*y - ZZ(3)*x*z + x));

    // f2 = x(2x - 2y - 3z + 1)，所以 g = x · f1 · (2x - 2y - 3z + 1)
    // g 的每一项都包含 x（x 整除 g），但 factorize 未能分解
    auto fac = factorize(g);
    // 预期: fac.factors.size() >= 3 (含 x 因子)
    // 实际: fac.factors.size() == 1 (返回原多项式未分解)
}
```

**出错路径**

```
factorize(g)
  → __factor_multivar(g)
    → squarefreefactorize(g) → 1 个无平方分量 gk = g
    → get_variables(gk) → {x, y, z}，gk_vars.size() > 1
    → __wang_core(gk)            ← 没有提取 x 的公共幂次
      → 遍历主变量 x, y, z
        → x 主变量: L = 6z(2y-3z)(-2y+2z-2), 复杂 LC, 50 个点全部 lc_ok=0
        → y 主变量: deg(g,y)=2, 特化后 g|_{x=a,z=b} 总是不可约 → irred
        → z 主变量: L = 3 (常数), Hensel 成功但试除失败
      → 所有主变量均失败 → 返回 {{g, 1}}
```

**关键问题**：g 的每一项都含 x（`min_term(deg_x) = 1`），即 `g = x · h(x,y,z)`，但 `squarefreefactorize` 和 `cont()` 都不提取这个 x 因子。`cont()` 计算的是首变量的多项式系数的 GCD，不是各变量的最小幂次。

### Bug 2：不可约启发式误判（已修复，补文档）

**复现用例**：同上。y 主变量下，因子 f2 = x(2x-2y-3z+1) 不含 y，特化 y=α 后 f2 变为常数，导致单变量分解总给出 1 个因子。连续 3 个求值点不可约 → 启发式误判为不可约 → `return {{g, 1}}`。

**已修复**（第 2492 行 `return {{g,1}}` → `continue`），但未经正式流程。此文档补充记录。

---

## 第二部分：根因分析

### Bug 1 根因：缺失单项式内容提取

多变量因式分解的标准流程：

```
f → 提取整数 content → 提取单项式 content → 无平方分解 → Wang/Hensel
```

CLPoly 当前流程：

```
f → squarefreefactorize (内含整数 content) → Wang/Hensel
         ↑
     缺失：单项式 content 提取
```

**单项式 content (monomial content)** 定义：

```
monomial_content(f) = x₁^a₁ · x₂^a₂ · ... · xₙ^aₙ
其中 aᵢ = min_{所有项 t} deg_{xᵢ}(t)
```

若任一 aᵢ > 0，则 xᵢ^aᵢ 是 f 的因子。提取后 `f' = f / monomial_content(f)` 满足"每个变量至少有一项中幂次为 0"。

**为什么 `cont()` 不做这件事**：

CLPoly 的 `cont(f)` (`polynomial_gcd.hh:592-653`) 计算的是 f 视为首变量 x₁ 的多项式时，各项系数（本身是关于 x₂,...,xₙ 的多项式）的 GCD。这是**多项式 content**，不是**单项式 content**。

例：`f = x²y + x²z` → `cont(f) = y + z`（x² 的系数），`pp(f) = x²`。
但：`f = x²y + xy² = xy(x+y)` → `cont(f) = y`（x² 系数为 y，x 系数为 y²，GCD=y），`pp(f) = x² + xy = x(x+y)`。此时 x 因子未被提取。

### Bug 2 根因：不可约启发式作用域过广

3-连续不可约启发式的假设是"如果 3 个随机点都给出不可约，f 很可能关于该主变量不可约"。但当某个因子不含该主变量时，该因子特化后必为常数，单变量结果**永远**只有 1 个因子——启发式 100% 误触发。

修复方案（已实施）：`return {{g,1}}` → `continue`，只跳过该主变量而非宣告整体不可约。

---

## 第三部分：参考实现对照

### SymPy：`dmp_terms_gcd`

**文件**：`sympy/polys/densebasic.py`

```python
def dmp_terms_gcd(f, u, K):
    """Remove GCD of terms from f in K[X]."""
    if dmp_ground_TC(f, u, K) or dmp_zero_p(f, u):
        return (0,) * (u + 1), f    # 有常数项 → 无公共变量幂

    F = dmp_to_dict(f, u, K)
    G = monomial_min(*list(F.keys()))  # 各变量取最小幂次
    # ... 减去 G，返回 (G, f_reduced) ...
```

**调用位置**：`dmp_factor_list` 的**第一步**（在 Wang 算法之前）：

```python
def dmp_factor_list(f, u, K0):
    J, f = dmp_terms_gcd(f, u, K0)      # ← 第一步
    cont, f = dmp_ground_primitive(f, u, K0)
    # ... Wang/Hensel ...
    # 最后将 J 中的变量幂次作为因子追加回结果
```

### FLINT：`mpoly_remove_var_powers`

**文件**：`src/mpoly/remove_var_powers.c`

```c
void mpoly_remove_var_powers(fmpz *var_powers, ulong *Aexps, ...)
{
    // 1. mpoly_min_fields_fmpz: 扫描所有项，各变量取最小幂次
    // 2. 若最小幂次不全为 0，从每项减去最小幂次（原地修改）
    // 3. 各变量最小幂次写入 var_powers 数组
}
```

**调用位置**：`fmpz_mpoly_factor_content` 的**最开头**（在 Wang/Zassenhaus 之前）：

```c
int fmpz_mpoly_factor_content(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, ...) {
    // 1. 提取整数 content
    // 2. mpoly_remove_var_powers(...)   ← 提取变量幂次
    //    对每个 v: if (var_powers[v] > 0) → 追加因子 (x_v, var_powers[v])
    // 3. 递归 content 分离...
}
```

### 对照总结

| | CLPoly（当前） | SymPy | FLINT |
|---|---|---|---|
| 单项式 content 提取 | **缺失** | `dmp_terms_gcd`（第一步） | `mpoly_remove_var_powers`（第一步） |
| 整数 content 提取 | `squarefreefactorize` 内部 | `dmp_ground_primitive` | `_fmpz_vec_content` |
| 提取时机 | - | Wang 之前 | Wang 之前 |

### 用复现用例走一遍

```
g = x · (-2xy + 2yz - 2y + 3z) · (2x - 2y - 3z + 1)

展开后 g 的所有项: x³(系数...) + x²(...) + x(...)
每项都含 x，min_deg_x = 1

SymPy: dmp_terms_gcd → J = (1, 0, 0)
  g' = g / x → (-2xy + 2yz - 2y + 3z) · (2x - 2y - 3z + 1)
  因子 (x, 1) 记入结果
  对 g' 做 Wang → 成功分解为 2 个因子

FLINT: mpoly_remove_var_powers → var_powers = [1, 0, 0]
  同上

CLPoly (当前):
  squarefreefactorize → gk = g (整个多项式)
  __wang_core(g) → 所有主变量失败 → 返回原多项式
```

---

## 第四部分：修正方案

### 修改位置

`clpoly/polynomial_factorize.hh`，`__factor_multivar` 函数，第 2627 行 `else` 分支内（调用 `__wang_core` 之前）。

### 修改内容

在 `__wang_core(gk_pos)` 调用前，增加单项式 content 提取步骤。

**新增辅助函数** `__extract_monomial_content`：

```cpp
// 提取多项式的单项式内容：对每个变量计算所有项中的最小幂次，
// 将公共变量幂次作为因子提取，返回除去公共幂次后的多项式。
//
// 例：f = x²yz + x³y²z = x²yz(1 + xy) 中
//   min_deg(x) = 2, min_deg(y) = 1, min_deg(z) = 1
//   提取 x²yz 后 f' = 1 + xy
//   var_factors = [(x, 2), (y, 1), (z, 1)]
template<class var_order>
polynomial_<ZZ, lex_<var_order>>
__extract_monomial_content(
    const polynomial_<ZZ, lex_<var_order>>& f,
    std::vector<std::pair<variable, int64_t>>& var_factors)
{
    using Poly = polynomial_<ZZ, lex_<var_order>>;
    var_factors.clear();

    if (f.empty()) return f;

    // Step 1: 收集所有变量及其在各项中的最小幂次
    // 从第一项的变量集合开始，逐项取 min
    std::map<variable, int64_t> min_deg;
    bool first_term = true;

    for (const auto& [mono, coeff] : f)
    {
        if (first_term)
        {
            for (const auto& [var, deg] : mono)
                min_deg[var] = deg;
            first_term = false;
        }
        else
        {
            // 对已知变量取 min
            std::set<variable> present;
            for (const auto& [var, deg] : mono)
            {
                present.insert(var);
                auto it = min_deg.find(var);
                if (it != min_deg.end())
                    it->second = std::min(it->second, deg);
                // 新变量不加入 min_deg（说明之前有项不含它，min=0）
            }
            // 不在本项中的变量 → min 降为 0
            auto it = min_deg.begin();
            while (it != min_deg.end())
            {
                if (present.find(it->first) == present.end())
                    it = min_deg.erase(it);
                else
                    ++it;
            }
        }
    }

    // Step 2: 移除 min_deg 为 0 的变量
    for (auto it = min_deg.begin(); it != min_deg.end(); )
    {
        if (it->second == 0)
            it = min_deg.erase(it);
        else
            ++it;
    }

    if (min_deg.empty()) return f;  // 无公共变量幂

    // Step 3: 构造提取后的多项式（每项的单项式减去 min_deg）
    Poly result(f.comp_ptr());
    for (const auto& [mono, coeff] : f)
    {
        typename Poly::monomial_type new_mono;
        for (const auto& [var, deg] : mono)
        {
            auto it = min_deg.find(var);
            int64_t subtract = (it != min_deg.end()) ? it->second : 0;
            if (deg > subtract)
                new_mono.push_back({var, deg - subtract});
        }
        result.data().push_back({std::move(new_mono), coeff});
    }
    result.normalization();

    // Step 4: 记录提取的变量因子
    for (const auto& [var, deg] : min_deg)
        var_factors.push_back({var, deg});

    return result;
}
```

**修改 `__factor_multivar`**（第 2627-2647 行）：

```cpp
// 修改前
else
{
    // 确保传给 __wang_core 的多项式 lc > 0
    Poly gk_pos = gk;
    bool negated = false;
    if (!gk_pos.empty() && gk_pos.front().second < 0)
    {
        for (auto& term : gk_pos.data())
            term.second = -term.second;
        negated = true;
    }
    auto wang_factors = __wang_core(gk_pos);
    for (auto& [fi, ei] : wang_factors)
        result.factors.push_back({fi, ei * mk});
    if (negated)
    {
        if (mk % 2 == 1)
            result.content = -result.content;
    }
}

// 修改后
else
{
    // 预处理：提取单项式 content（公共变量幂次）
    std::vector<std::pair<variable, int64_t>> var_factors;
    Poly gk_reduced = __extract_monomial_content(gk, var_factors);

    // 将提取的变量幂次作为因子加入结果
    for (auto& [var, vdeg] : var_factors)
    {
        Poly var_poly(gk.comp_ptr());
        typename Poly::monomial_type mono;
        mono.push_back({var, 1});
        var_poly.data().push_back({std::move(mono), ZZ(1)});
        var_poly.normalization();
        result.factors.push_back({std::move(var_poly), vdeg * mk});
    }

    // 检查提取后是否还需要多变量分解
    auto reduced_vars = get_variables(gk_reduced);
    if (reduced_vars.size() <= 1)
    {
        // 降为单变量，用单变量分解
        auto sub = factorize(gk_reduced);
        for (int64_t e = 0; e < mk; ++e)
            result.content *= sub.content;
        for (auto& [fi, ei] : sub.factors)
            result.factors.push_back({fi, ei * mk});
    }
    else
    {
        // 确保传给 __wang_core 的多项式 lc > 0
        Poly gk_pos = gk_reduced;
        bool negated = false;
        if (!gk_pos.empty() && gk_pos.front().second < 0)
        {
            for (auto& term : gk_pos.data())
                term.second = -term.second;
            negated = true;
        }
        auto wang_factors = __wang_core(gk_pos);
        for (auto& [fi, ei] : wang_factors)
            result.factors.push_back({fi, ei * mk});
        if (negated)
        {
            if (mk % 2 == 1)
                result.content = -result.content;
        }
    }
}
```

### 修正后复现用例的执行过程

```
g = x · (-2xy + 2yz - 2y + 3z) · (2x - 2y - 3z + 1)

Step 1: squarefreefactorize(g) → [(gk=g, mk=1)]

Step 2: __extract_monomial_content(gk)
  遍历所有项，计算各变量最小幂次:
    min_deg(x) = 1, min_deg(y) = 0, min_deg(z) = 0
  var_factors = [(x, 1)]
  gk_reduced = g / x = (-2xy + 2yz - 2y + 3z) · (2x - 2y - 3z + 1)
  结果追加因子 (x, 1)

Step 3: get_variables(gk_reduced) → {x, y, z}，仍然多变量

Step 4: __wang_core(gk_reduced)
  gk_reduced 不含公共变量幂，LC 结构更简单
  → 成功分解为 2 个因子

最终结果: x · (-2xy + 2yz - 2y + 3z) · (2x - 2y - 3z + 1) ✓
```

### 提取时机：squarefree 之后（与 SymPy/FLINT 不同）

SymPy/FLINT 在 squarefree **之前**提取单项式 content，我们在 squarefree **之后**对每个无平方分量分别提取。

**正确性**：若 gk 是无平方的，则任一变量的 min_deg 至多为 1（否则该变量的平方整除 gk，与无平方矛盾）。因此两种顺序均正确——提取的变量幂次最终都会被正确记录。

**选择后提取的理由**：避免在提取后重复做 squarefree 分解。代价仅为对每个无平方分量遍历一次各项取 min，开销可忽略。

### 边界情况：gk 为纯单项式

当 gk 是纯单项式（如 `xy`）时，提取后 `gk_reduced = 1`（常数），`get_variables` 返回空集（size=0 ≤ 1），走 `factorize(gk_reduced)` 分支。`factorize` 对常数多项式返回 content=常数值、无因子，行为正确。

### 不改动的部分

- `cont()` / `pp()` 函数不变（它们的语义是多项式 content，有其他用途）
- `__wang_core` 内部逻辑不变
- `__select_eval_point` 不变
- 不可约启发式 `continue` 修正保留（Bug 2 已修复，仅补文档）

---

## 第五部分：主变量轮换的必要性分析

### 背景

在 LC 原子分配 bug 修复（per-factor 幂次提取）之前，添加了主变量轮换机制（遍历所有变量作为主变量）。现在 LC 分配已修正，需要分析轮换是否仍然必要。

### 其他库的做法

| 系统 | 轮换主变量？ | 选择策略 |
|------|------------|---------|
| Wang 原文 (1978) | 否 | 固定 x₁ |
| SymPy | 否 | 固定第一个变量（DMP 序），失败后只换求值点 |
| FLINT | 否 | 固定第一个变量（强制 LEX 序），失败后换随机 α 值 |
| Singular/Factory | **是** | LC 大小 + 度数排序 → 搜索无平方变量 → `swapvar`，最终兜底域扩展 |
| CLPoly（当前） | **是** | 自然变量序，无排序优化 |

**关键发现**：SymPy 和 FLINT 都**不做**主变量轮换。Singular 是主流库中唯一做类似设计的，但它的策略更精细（按 LC 大小和度数排序，搜索无平方变量）。

文献推荐的标准启发式（Lucks 1986）："选度数最低的变量作为主变量，按度数升序提升其余变量。"但同时承认"很难事先判断哪个变量最优"。

### 分析

主变量轮换解决的核心问题是：**不同主变量下 LC 结构不同，某些主变量的 LC 更简单、更容易匹配**。

即使 per-factor 提取修正了匹配算法，主变量选择仍然影响：

1. **LC 复杂度**：主变量 x₁ 下 LC = lc(f, x₁)，可能是高次多变量多项式，其不可约分解可能很复杂。换一个主变量可能得到更简单（甚至常数）的 LC。

2. **|lⱼ(α)| ≤ 1 拒绝率**：如果 LC 的不可约因子在求值点处值为 ±1 或 0，会被拒绝。不同主变量 → 不同 LC → 不同拒绝率。

3. **Content 污染率**：content(f₀) 与 lⱼ(α) 共享素因子的概率取决于 LC 和主变量的选择。

4. **Hensel 提升效率**：LC 越简单 → 缩放因子越小 → Hensel 提升收敛越快 → 试除成功率越高。

### 结论

**主变量轮换保留**，理由：

1. 它提供的价值独立于 LC 匹配算法——即使匹配完美，不同主变量的 LC 复杂度差异仍然影响整体成功率。
2. CLPoly 不像 FLINT 那样有 Kaltofen 备选算法兜底，也不像 Singular 那样有域扩展兜底。轮换是我们唯一的"第二次机会"。
3. 实现已就绪，代码简洁，运行时开销仅在首选主变量失败时才产生。

**可优化方向**（非本次修复范围）：
- 当前按 `get_variables()` 的自然序遍历（按变量名排序），可优化为按 `degree(g, xᵢ)` 升序排列（参照 Lucks 1986 启发式）。
- 可参考 Singular 的多标准排序：主按 LC 复杂度、次按度数。

### 不可约启发式的调整

主变量轮换 + `continue`（而非 `return`）是正确的组合：
- 某个主变量下不可约启发式触发 → 仅表示该主变量不适合 → `continue` 换下一个主变量
- 所有主变量都试完后才宣告不可约 → `return {{g, 1}}`（第 2594 行）
- 这与 Singular 的 `newMainVariableSearch` 的逻辑一致：某变量不可用时搜索下一个
