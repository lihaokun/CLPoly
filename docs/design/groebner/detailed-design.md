# Gröbner 基细化设计文档

> 前置：[架构文档](architecture.md)
> 本文档逐模块列出需实现的函数、签名、算法、调用关系和复用点。

---

## M1: 单项式 LCM

**所在文件：** `clpoly/monomial.hh`（`is_divexact` 之后）

### 函数

```cpp
template<class compare>
basic_monomial<compare> lcm(
    const basic_monomial<compare>& m1,
    const basic_monomial<compare>& m2)
```

**功能：** 计算 lcm(m1, m2)，逐变量取指数最大值。

**算法：** 双指针合并（参考 `pair_vec_add`，basic.hh:166），将加法替换为 max：

```
result.clear()
i ← 0, j ← 0, deg ← 0
while i < |m1| and j < |m2|:
    if m1[i].var == m2[j].var:
        e ← max(m1[i].exp, m2[j].exp)
        result.push_back({var, e})
        deg += e
        i++, j++
    else if comp(m1[i].var, m2[j].var):
        result.push_back(m1[i])
        deg += m1[i].exp
        i++
    else:
        result.push_back(m2[j])
        deg += m2[j].exp
        j++
拷贝 m1 或 m2 的剩余部分，累加 deg
设置 result.__deg = deg
设置 result.__comp_ptr（从 m1 或 m2 获取）
return result
```

**复用：**
- 参考 `gcd(basic_monomial, basic_monomial)`（polynomial_gcd.hh:20）的结构，取 min→取 max
- 参考 `pair_vec_add`（basic.hh:166）的双指针合并模式

**复杂度：** O(v)，v = 变量个数

---

## 辅助函数: __make_monic

**所在文件：** `clpoly/groebner.hh`

```cpp
template<class Tc, class comp>
inline void __make_monic(polynomial_<Tc, comp>& f)
```

**功能：** 将多项式首一化（逐系数除以 LC）。

**说明：** CLPoly 没有 `polynomial /= Tc` 运算符，因此无法直接写 `f /= LC(f)`。
此辅助函数遍历 `f.data()`，将每个系数除以 `f.front().second`。
若 `f` 为零多项式或已首一（`LC = 1`）则直接返回。

**算法：**
```
if f.empty(): return
lc ← f.front().second
if lc == 1: return
for (auto& term : f.data()):
    term.second = term.second / lc
```

**复用：** `Tc::operator/`（QQ.hh:139）

---

## M2: S-多项式

**所在文件：** `clpoly/groebner.hh`

### 函数 1: 内部版本

```cpp
template<class Tc, class comp>
inline polynomial_<Tc, comp> __s_polynomial(
    const polynomial_<Tc, comp>& f,
    const polynomial_<Tc, comp>& g)
```

**功能：** 计算 S(f, g) = (L/LT(f))·f - (L/LT(g))·g，其中 L = lcm(LM(f), LM(g))。

**算法：**

```
LM_f ← f.front().first          // basic_polynomial::front()
LC_f ← f.front().second
LM_g ← g.front().first
LC_g ← g.front().second

L ← lcm(LM_f, LM_g)            // M1

m_f ← L / LM_f                  // basic_monomial::operator/ (basic_monomial.hh:199)
m_g ← L / LM_g

// 逐项乘（通过 data() 获取可写 vector 引用）
// (1/LC_f) · m_f · f
poly_f ← f                            // 拷贝
for (auto& term : poly_f.data()):     // data() 返回 vector<pair<Tm,Tc>>&
    term.first = term.first * m_f     // 单项式乘法
    term.second = term.second / LC_f  // 系数除法

// (1/LC_g) · m_g · g
poly_g ← g                            // 拷贝
for (auto& term : poly_g.data()):
    term.first = term.first * m_g
    term.second = term.second / LC_g

S ← poly_f - poly_g             // pair_vec_sub (basic.hh:244)
return S
```

**实现说明：**
- 逐项乘同一单项式保持项的相对顺序不变（单项式序的平移不变性：
  a < b ⟹ a·m < b·m），因此乘后无需重新排序，`pair_vec_sub` 的已排序前提成立。
- 逐项乘比构造单项多项式再调 `pair_vec_multiplies` 更高效，
  避免 O(n log n) 的卷积开销，只需 O(n) 逐项乘。
- 单项式乘法: `basic_monomial::operator*`（basic_monomial.hh:183）
- 系数除法: `Tc::operator/`（QQ.hh:139）
- 多项式减法: `pair_vec_sub`（basic.hh:244）

**调用关系：** M1(lcm) → M2(__s_polynomial)

### 函数 2: 公开版本

```cpp
template<class Tc, class comp>
inline polynomial_<Tc, comp> s_polynomial(
    const polynomial_<Tc, comp>& f,
    const polynomial_<Tc, comp>& g)
```

**功能：** 公开 API，直接调用 `__s_polynomial`。

---

## M3: Normal Form（多除数约化）

**所在文件：** `clpoly/groebner.hh`

### 函数 1: 基本版

```cpp
template<class Tc, class comp>
inline polynomial_<Tc, comp> __normal_form(
    const polynomial_<Tc, comp>& f,
    const std::vector<polynomial_<Tc, comp>>& G)
```

**功能：** 计算 NF(f, G)。供 M6 Ideal 类和 M5 互约使用。

**算法：**

```
h ← f
r ← 0（空多项式，设置 comp_ptr）

while h ≠ 0:                      // h.empty() 判零
    divided ← false
    for i = 0 to |G|-1:
        basic_monomial<comp> quot_m
        if is_divexact(quot_m, h.front().first, G[i].front().first):
            // LM(G[i]) | LM(h)
            c ← h.front().second / G[i].front().second    // Tc::operator/
            // h ← h - c · quot_m · G[i]
            // 构造 c·quot_m·G[i] 然后用 pair_vec_sub
            tmp ← G[i] 逐项乘 quot_m 再乘 c
            h ← h - tmp                                    // pair_vec_sub
            divided ← true
            break
    if not divided:
        // 首项不可约化，移入余式：
        r.data().push_back(h.front())          // 直接追加到 r 尾部（有序）
        h.data().erase(h.data().begin())       // 删去 h 的首项

return r
```

**实现说明：**
- `is_divexact`（monomial.hh:68）同时返回商单项式 `quot_m`
- "不可约首项移入余式"直接用 `data().push_back` + `data().erase(begin())` 操作，
  简单直接。`erase(begin())` 对 vector 是 O(n)，但此分支执行频率低于约化分支
- r 的各项按单项式序单调递减追加（因为 h 的首项在循环中严格递减），
  因此 r 天然有序，无需排序

**复用：**
- `is_divexact`（monomial.hh:68）
- `basic_monomial::operator*`（basic_monomial.hh:183）
- `pair_vec_sub`（basic.hh:244）
- `basic_polynomial::front()`（basic_polynomial.hh:140）
- `basic_polynomial::empty()`（basic_polynomial.hh:153）

### 函数 2: Sugar 跟踪版

```cpp
template<class Tc, class comp>
inline polynomial_<Tc, comp> __normal_form_with_sugar(
    const polynomial_<Tc, comp>& f,
    const std::vector<polynomial_<Tc, comp>>& G,
    const std::vector<int64_t>& sugar_vec,
    int64_t& f_sugar)
```

**功能：** 同 `__normal_form`，但在约化过程中更新 sugar 度数。供 M5 主循环使用。

**与基本版的差异：** 在约化步骤中增加 sugar 跟踪：

```
// 约化前先保存 LM(h) 的度数
lm_h_deg ← h.front().first.deg()

// 执行约化
h ← h - c · quot_m · G[i]

// 用保存的度数更新 sugar
f_sugar ← max(f_sugar, sugar_vec[i] + lm_h_deg - G[i].front().first.deg())
```

**关键：** `lm_h_deg` 必须在约化之前取值，因为约化后 h 的首项已改变。

### 函数 3: 公开版本

```cpp
template<class Tc, class comp>
inline polynomial_<Tc, comp> normal_form(
    const polynomial_<Tc, comp>& f,
    const std::vector<polynomial_<Tc, comp>>& G)
```

**功能：** 公开 API，直接调用 `__normal_form`。

---

## M4: 临界对管理

**所在文件：** `clpoly/groebner.hh`

### 结构体

```cpp
struct __critical_pair {
    size_t i, j;           // G 中多项式索引, i < j
    int64_t sugar;         // Sugar 度数
    int64_t lcm_deg;       // lcm(LM(G[i]), LM(G[j])) 的总度数
};
```

**比较器（用于选择最小对）：**

```cpp
struct __cp_compare {
    bool operator()(const __critical_pair& a, const __critical_pair& b) const {
        if (a.sugar != b.sugar) return a.sugar < b.sugar;    // sugar 小优先
        return a.lcm_deg < b.lcm_deg;                       // lcm_deg 小优先
    }
};
```

**数据结构：** 使用 `std::vector<__critical_pair>` 存储所有对，选择操作（I5a）时用 `std::min_element` + `__cp_compare` 线性扫描取最小值并删除。不用优先队列，因为 `__update` 需要遍历和过滤旧对。

### 函数: Update（I5b）

```cpp
template<class Tc, class comp>
inline void __update(
    std::vector<polynomial_<Tc, comp>>& G,
    std::vector<__critical_pair>& pairs,
    const polynomial_<Tc, comp>& h,
    std::vector<int64_t>& sugar_vec,
    int64_t h_sugar)
```

**功能：** 将新元素 h 加入基 G，用 Gebauer-Möller 三准则过滤临界对。

**算法：**

```
new_idx ← |G|

// 1. 生成新对候选
new_pairs ← []
for k = 0 to |G|-1:
    L ← lcm(G[k].front().first, h.front().first)           // M1
    s ← max(sugar_vec[k] + L.deg() - G[k].front().first.deg(),
            h_sugar + L.deg() - h.front().first.deg())
    new_pairs.push({k, new_idx, s, L.deg()})

// 2. LCM 准则: 用 h 淘汰旧对
// 遍历 pairs，移除满足以下条件的对 (i,j):
//   LM(h) | lcm(LM(G[i]), LM(G[j]))
//   且 lcm(LM(G[i]), LM(h)) ≠ lcm(LM(G[i]), LM(G[j]))
//   且 lcm(LM(h), LM(G[j])) ≠ lcm(LM(G[i]), LM(G[j]))
// 注: is_divexact(op, a, b) 检查 b | a（即 op = a/b）
filtered_pairs ← []
for (i,j,s,d) in pairs:
    L_ij ← lcm(G[i].front().first, G[j].front().first)
    L_ih ← lcm(G[i].front().first, h.front().first)
    L_hj ← lcm(h.front().first, G[j].front().first)
    basic_monomial<comp> quot_m
    if is_divexact(quot_m, L_ij, h.front().first)    // LM(h) | L_ij
       and L_ih != L_ij and L_hj != L_ij:
        continue                 // 淘汰此对
    filtered_pairs.push({i,j,s,d})
pairs ← filtered_pairs

// 3. 在新对中自身淘汰（保留 lcm 最小的）
// 若 lcm(LM(G[k1]), LM(h)) 整除 lcm(LM(G[k2]), LM(h))，淘汰 (k2, new_idx)
// 实现: 按 lcm_deg 排序后，用已接受集合过滤
对 new_pairs 按 lcm_deg 升序排序
minimal_new ← []
for p in new_pairs:
    L_p ← lcm(G[p.i].front().first, h.front().first)
    被淘汰 ← false
    for q in minimal_new:
        L_q ← lcm(G[q.i].front().first, h.front().first)
        if is_divexact(_, L_p, L_q):      // L_q | L_p
            被淘汰 ← true
            break
    if not 被淘汰:
        minimal_new.push(p)
new_pairs ← minimal_new

// 4. 乘积准则: gcd(LM(G[k]), LM(h)) = 1 的对丢弃
// 判定: gcd(monomial, monomial) 为空单项式（deg = 0）
//       用 polynomial_gcd.hh:20 的 gcd(basic_monomial, basic_monomial)
final_new ← []
for p in new_pairs:
    g ← gcd(G[p.i].front().first, h.front().first)   // polynomial_gcd.hh:20
    if g.deg() != 0:                                   // 非互素
        final_new.push(p)
    // 互素 → 丢弃（乘积准则）
new_pairs ← final_new

// 5. 将 h 加入 G，新对加入 pairs
G.push_back(h)
sugar_vec.push_back(h_sugar)
pairs.insert(pairs.end(), new_pairs.begin(), new_pairs.end())
```

**注意：** 步骤 4 有一个微妙之处——如果步骤 3 后某些对被淘汰，但它们是步骤 3 中淘汰其他对的"依据"，那么乘积准则应在步骤 3 之后独立应用。上述顺序是正确的（先最小化，再乘积准则）。

**调用关系：**
- M1(lcm)
- `is_divexact`（monomial.hh:68）
- `gcd(basic_monomial)`（polynomial_gcd.hh:20）

**数据结构选择：** 由于 __update 需要遍历和过滤旧对（步骤 2），优先队列不适合（无法遍历）。改用 `std::vector<__critical_pair>` 存储所有对，选择操作（I5a）时线性扫描取最小值并删除。这牺牲了选择的 O(log n) 复杂度，但使过滤操作自然。对于首次实现，这是合理的折衷。

---

## M5: Buchberger 主循环

**所在文件：** `clpoly/groebner.hh`

### 函数 1: 互约

```cpp
template<class Tc, class comp>
inline void __interreduce(std::vector<polynomial_<Tc, comp>>& G)
```

**功能：** 将 Gröbner 基化简为约化 Gröbner 基。

**算法：**

```
// 步骤 1: 删除冗余（LM(gi) | LM(gj) 则删 gj）
i ← 0
while i < |G|:
    redundant ← false
    for j = 0 to |G|-1, j ≠ i:
        basic_monomial<comp> tmp
        if is_divexact(tmp, G[i].front().first, G[j].front().first):
            // LM(G[j]) | LM(G[i])，G[i] 冗余
            redundant ← true
            break
    if redundant:
        G.erase(G.begin() + i)
    else:
        i++

// 步骤 2: 完全约化每个元素
for i = 0 to |G|-1:
    others ← G[0..i-1] ∪ G[i+1..|G|-1]     // G 除去第 i 个
    G[i] ← __normal_form(G[i], others)       // M3 基本版
    // 说明: 步骤 1 已删除冗余元素（LM 被其他元素整除的），
    // 因此 G[i] 的首项不会被 others 约化掉，NF 结果非零。
    __make_monic(G[i])                        // 首一化
```

**调用关系：** M3(__normal_form), `is_divexact`, `__make_monic`

### 函数 2: Buchberger 主循环

```cpp
template<class Tc, class comp>
inline std::vector<polynomial_<Tc, comp>> __buchberger(
    std::vector<polynomial_<Tc, comp>> F)
```

**功能：** Buchberger 算法核心，返回约化 Gröbner 基。

**算法：**

```
// 1. 预处理：过滤零多项式，首一化
G ← []
for f in F:
    if not f.empty():
        __make_monic(f)             // 首一化
        G.push_back(f)
if G.empty(): return G

// 2. 初始化 sugar 度数
sugar_vec ← []
for g in G:
    sugar_vec.push_back(degree(g))   // polynomial.hh:40

// 3. 初始化临界对（逐个加入，复用 __update 的 Gebauer-Möller 过滤）
pairs ← []
// 从 G[0] 开始，逐个加入 G[1], G[2], ...
// 用 __update 保证三准则在初始化阶段也被完整应用
init_G ← [G[0]]
init_sugar ← [sugar_vec[0]]
for k = 1 to |G|-1:
    __update(init_G, pairs, G[k], init_sugar, sugar_vec[k])
G ← init_G
sugar_vec ← init_sugar

// 4. 主循环
while pairs 非空:
    // I5a: 选择 (sugar, lcm_deg) 最小的对
    min_idx ← argmin over pairs by (sugar, lcm_deg)
    cp ← pairs[min_idx]
    pairs.erase(pairs.begin() + min_idx)

    // S-多项式
    h ← __s_polynomial(G[cp.i], G[cp.j])             // M2

    // Sugar 度数
    h_sugar ← cp.sugar

    // Normal form with sugar
    h ← __normal_form_with_sugar(h, G, sugar_vec, h_sugar)    // M3

    if not h.empty():
        __make_monic(h)             // 首一化
        __update(G, pairs, h, sugar_vec, h_sugar)              // M4

// 5. 互约
__interreduce(G)

// 6. 返回
return G
```

**调用关系：**
- M2(__s_polynomial)
- M3(__normal_form_with_sugar, __normal_form via __interreduce)
- M4(__update)（内部调用 M1(lcm) 和 gcd(basic_monomial)）
- `degree()`（polynomial.hh:40）

---

## M6: API 层

**所在文件：** `clpoly/groebner.hh`

### 函数 1: QQ 入口

```cpp
template<class comp>
inline std::vector<polynomial_<QQ, comp>> groebner_basis(
    const std::vector<polynomial_<QQ, comp>>& generators)
```

**功能：** QQ 上的 Gröbner 基，直接调用 `__buchberger`。

**实现：**
```cpp
return __buchberger(generators);
```

### 函数 2: ZZ 入口

```cpp
template<class comp>
inline std::vector<polynomial_<ZZ, comp>> groebner_basis(
    const std::vector<polynomial_<ZZ, comp>>& F_zz)
```

**功能：** ZZ 上的 Gröbner 基，ZZ→QQ→ZZ 工作流。

**算法：**

```
// 1. ZZ → QQ
F_qq ← []
for f in F_zz:
    polynomial_<QQ, comp> fq
    poly_convert(f, fq)              // polynomial_convert.hh:63
    F_qq.push_back(fq)

// 2. 在 QQ 上计算
G_qq ← __buchberger(F_qq)

// 3. QQ → ZZ（清分母 + 本原化）
G_zz ← []
for g in G_qq:
    polynomial_<ZZ, comp> gz
    poly_convert(g, gz)              // polynomial_convert.hh:118, 自动 LCD 通分

    // 本原化：计算所有系数的 GCD
    ZZ c = abs(gz.front().second)    // 从 LC 开始
    for (auto& term : gz):
        c = gcd(c, abs(term.second))
        if c == 1: break
    if gz.front().second < 0:
        c = -c                       // 确保 lc > 0
    // 逐系数除以 c
    for (auto& term : gz):
        term.second /= c

    G_zz.push_back(gz)

// 4. 返回
return G_zz
```

**复用：**
- `poly_convert(ZZ→QQ)`（polynomial_convert.hh:64）
- `poly_convert(QQ→ZZ)`（polynomial_convert.hh:128）
- `gcd(ZZ, ZZ)`（mpz_class 的 gcd）

**说明：** 多变量 `cont()`（polynomial_gcd.hh:586）返回 `polynomial_<ZZ, lex_>`，
是按主变量的多项式 content，不适用于此处。这里需要的是整数 content（所有系数的 GCD），
与 `cont(upolynomial_<ZZ>)`（polynomial_gcd.hh:1069）逻辑相同，但对象不同。
直接在 ZZ 入口中内联实现，无需新增辅助函数。

### 函数 3: Ideal 类

```cpp
template<class Tc, class comp = grlex>
class Ideal {
    std::vector<polynomial_<Tc, comp>> generators_;
    mutable std::vector<polynomial_<Tc, comp>> gb_cache_;
    mutable bool gb_computed_ = false;

public:
    Ideal(std::vector<polynomial_<Tc, comp>> gens)
        : generators_(std::move(gens)) {}

    Ideal(std::initializer_list<polynomial_<Tc, comp>> gens)
        : generators_(gens) {}

    const std::vector<polynomial_<Tc, comp>>& generators() const {
        return generators_;
    }

    const std::vector<polynomial_<Tc, comp>>& groebner_basis() const {
        if (!gb_computed_) {
            gb_cache_ = clpoly::groebner_basis(generators_);
            gb_computed_ = true;
        }
        return gb_cache_;
    }

    bool contains(const polynomial_<Tc, comp>& f) const {
        auto nf = normal_form(f, groebner_basis());     // M3 公开版
        return nf.empty();
    }

    polynomial_<Tc, comp> reduce(const polynomial_<Tc, comp>& f) const {
        return normal_form(f, groebner_basis());         // M3 公开版
    }
};
```

---

## 附录: groebner.hh 完整结构

```
/**
 * @file groebner.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief Gröbner basis computation
 */
#ifndef CLPOLY_GROEBNER_HH
#define CLPOLY_GROEBNER_HH

#include <clpoly/polynomial.hh>
#include <clpoly/polynomial_gcd.hh>
#include <clpoly/polynomial_convert.hh>
#include <vector>
#include <algorithm>
#include <cstdint>

namespace clpoly {

    // ===== M4: 临界对结构 =====
    struct __critical_pair { ... };
    struct __cp_compare { ... };

    // ===== 辅助 =====
    template<...> inline void __make_monic(...);

    // ===== M2: S-多项式 =====
    template<...> inline polynomial_<Tc,comp> __s_polynomial(...);

    // ===== M3: Normal Form =====
    template<...> inline polynomial_<Tc,comp> __normal_form(...);
    template<...> inline polynomial_<Tc,comp> __normal_form_with_sugar(...);

    // ===== M4: 临界对管理 =====
    template<...> inline void __update(...);

    // ===== M5: Buchberger =====
    template<...> inline void __interreduce(...);
    template<...> inline std::vector<polynomial_<Tc,comp>> __buchberger(...);

    // ===== M6: 公开 API =====
    template<...> inline polynomial_<Tc,comp> s_polynomial(...);
    template<...> inline polynomial_<Tc,comp> normal_form(...);
    template<...> inline std::vector<polynomial_<QQ,comp>> groebner_basis(...);
    template<...> inline std::vector<polynomial_<ZZ,comp>> groebner_basis(...);
    template<...> class Ideal { ... };

} // namespace clpoly
#endif
```

---

## 附录: 修改文件清单

### clpoly/monomial.hh

在 `is_divexact` 之后新增 `lcm` 函数（M1）。

### clpoly/clpoly.hh

在最后一个 `#include` 之后新增：
```cpp
#include <clpoly/groebner.hh>
```

### test/run_all_tests.sh

新增 `test/test_groebner`。

### .gitignore

新增 `test/test_groebner`。

---

## 附录: 已有函数复用索引

| 已有函数 | 位置 | 被哪些新函数调用 |
|---------|------|----------------|
| `is_divexact(op, m1, m2)` | monomial.hh:68 | M3(__normal_form), M4(__update), M5(__interreduce) |
| `gcd(m1, m2)` (monomial) | polynomial_gcd.hh:20 | M4(__update) |
| `basic_monomial::operator*` | basic_monomial.hh:183 | M2(__s_polynomial), M3(__normal_form) |
| `basic_monomial::operator/` | basic_monomial.hh:199 | M2(__s_polynomial) |
| `basic_monomial::deg()` | basic_monomial.hh:37 | M4(__update), M5(__buchberger) |
| `basic_polynomial::front()` | basic_polynomial.hh:140 | M2, M3, M4, M5 (获取 LT/LM/LC) |
| `basic_polynomial::empty()` | basic_polynomial.hh:153 | M3, M5 (零多项式判定) |
| `pair_vec_sub` | basic.hh:244 | M2(__s_polynomial), M3(__normal_form) |
| `poly_convert(ZZ→QQ)` | polynomial_convert.hh:64 | M6(ZZ 入口) |
| `poly_convert(QQ→ZZ)` | polynomial_convert.hh:128 | M6(ZZ 入口) |
| `degree(f)` | polynomial.hh:40 | M5(__buchberger 初始化 sugar) |
| `Tc::operator/` | QQ.hh:139 | M2, M3 (系数除法) |
| `Tc::operator/` | QQ.hh:139 | __make_monic, M2, M3 (系数除法/首一化) |
