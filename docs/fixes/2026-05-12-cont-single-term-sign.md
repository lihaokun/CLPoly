# 修正方案：cont 单项情形 sign 不一致

> **状态**：草稿，待用户审核。
> 对应 workflow.md §5.1 算法内部错误修正流程
> 对应 GitHub issue: [#17](https://github.com/lihaokun/CLPoly/issues/17)
> 对应分支：`fix/cont-single-term-sign-issue-17`（从 master 起，独立 PR 对 master）

---

## 第一部分：现象与定位

### 现象

CLPoly 的 `cont()` 对单项与多项多项式返回不一致的 sign：

| 输入 | CLPoly `cont()` | Maple `content()` | SymPy `Poly.content()` | FLINT `fmpz_poly_content` |
|------|----------------|------------------|----------------------|----------------------|
| `x` | 1 | 1 | 1 | 1 |
| `-x` | **-1** ✗ | 1 | 1 | 1 |
| `2*x` | 2 | 2 | 2 | 2 |
| `-2*x` | **-2** ✗ | 2 | 2 | 2 |
| `-x+1` | 1 | 1 | 1 | 1 |
| `-2*x+2` | 2 | 2 | 2 | 2 |

**单项时 CLPoly 与三家不一致**。

### 实测验证

```python
# SymPy
>>> from sympy import symbols, Poly, ZZ
>>> x = symbols('x')
>>> [Poly(f, x, domain=ZZ).content() for f in [x, -x, 2*x, -2*x, x-1, -x+1, -2*x+2]]
[1, 1, 2, 2, 1, 1, 2]
```

```cpp
// CLPoly C++ (实测，univariate ZZ)
cont(x)        = 1     ✓
cont(-x)       = -1    ✗（应 1）
cont(2*x)      = 2     ✓
cont(-2*x)     = -2    ✗（应 2）
cont(-x+1)     = 1     ✓
cont(-2*x+2)   = 2     ✓
cont(常数 -6)  = -6    ✗（应 6 — 常数也是单项）
```

### 根因（C++ 源码）

`clpoly/polynomial_gcd.hh:1067`：

```cpp
inline ZZ cont(const upolynomial_<ZZ> &G) {
    if (G.empty()) return 1;
    auto ptr = G.begin();
    ZZ c = (ptr++)->second;          // ← c = leading coef（带符号）
    for (; ptr != G.end(); ++ptr) {  // ← 单项时此循环不执行
        c = gcd(c, ptr->second);     // 多项时 gcd 折正
    }
    return c;                        // 单项 → 直接返带符号 c
}
```

**单项 polynomial 时**：
- `c = first coef`（initial value 含 sign）
- 循环不执行
- 返回 c，sign 保留

**多项 polynomial 时**：
- `c = first coef`，然后 `c = gcd(c, ...)` 在循环中折正（GMP `mpz_gcd` 总返非负）
- 返回正 c

`clpoly/upolynomial.hh:220` 的 template 版本有同样逻辑问题（虽然该 template 可能不被实例化，因为有 ZZ 特化覆盖）。

`clpoly/polynomial_gcd.hh:581` 多变量 `cont()` 递归调用 polynomial_GCD，间接也受影响。

### 与 issue #14 的关系

issue #14（squarefreebasis unit leak）的根因之一就是 cont sign 不一致 → factor 非 canonical → 与 `polynomial_GCD` 归正首项后产生 associate 字面不等。PR #16 在 `squarefreefactorize` 层加 sign-normalization 兜底修了表象，但 cont 本身的不一致没修。

---

## 第二部分：参考实现对比

### Maple

`content(p, x)` 返回 GCD of coefficients of `p` w.r.t. `x`。**总是非负**。

### SymPy

`Poly.content()` 计算 ground domain (ZZ) 上系数的 GCD。**总是非负**。

实测见 §1。

### FLINT

`fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly)` 计算 polynomial 系数的 GCD。**总是非负**（基于 `fmpz_gcd` 取绝对值）。

文档原文：
> Sets res to the content of the polynomial poly. The content is the greatest common divisor of the absolute values of the coefficients of poly. **The content is always non-negative.**

### 共同约定

**所有工业系统：`content/cont` 始终返回非负值**（首项 sign 由 `sqf_list/sqrfree` 的 unit 部分单独处理）。

---

## 第三部分：修复方案

### 改动 1: `clpoly/polynomial_gcd.hh:1067`（ZZ 特化版）

```cpp
inline ZZ cont(const upolynomial_<ZZ> &G) {
    ZZ c = 0;  // ← 初始 c = 0，gcd(0, x) = |x|（mpz_gcd 取绝对值）
    for (auto &t : G) {
        c = gcd(c, t.second);
    }
    if (c == 0) return 1;  // 空 polynomial 约定返回 1
    return c;
}
```

#### 正确性半形式化证明

**Spec**（数学定义）：
```
cont(F) =
  if F 表示零多项式（无项 或 所有项 c_i=0）:
    return 1   -- CLPoly 约定（让 pp(f)=f/cont(f) 对零多项式良定义；与 Maple
              -- "cont(0)=0" 不同，CLPoly 历来约定为 1 — 本 PR 保留）
  else:
    return gcd(|c_0|, |c_1|, ..., |c_{n-1}|)   -- 系数绝对值的 GCD
```
其中 `gcd` 满足 GMP 约定：
- `gcd(a, b) ≥ 0` 总成立
- `gcd(0, x) = |x|`
- `gcd(a, b) = gcd(|a|, |b|)`（取绝对值）

**实现**（提议）：
```
c := 0
for t in G:
    c := gcd(c, t.second)
if c == 0:
    return 1
return c
```

**待证**：实现 = Spec。

**证明**：

*Lemma 1（loop invariant）*：迭代 k 次后，`c = gcd(|c_0|, |c_1|, ..., |c_{k-1}|)`。

证明（对 k 归纳）：
- **Base** (k=0)：未迭代，c=0。Spec: 空集 gcd 约定为 0。✓
- **Step**：假设迭代 k-1 后 `c_old = gcd(|c_0|, ..., |c_{k-2}|)` (k ≥ 1)。
  - 第 k 次迭代：`c := gcd(c_old, c_{k-1})`
  - 由 GMP 约定：`gcd(a, b) = gcd(|a|, |b|)`
  - 故 `c = gcd(|c_old|, |c_{k-1}|) = gcd(c_old, |c_{k-1}|)`（因 c_old ≥ 0 by IH）
  - `= gcd(gcd(|c_0|, ..., |c_{k-2}|), |c_{k-1}|) = gcd(|c_0|, ..., |c_{k-1}|)`（gcd 结合性）✓ □

*Theorem（主结论）*：实现 = Spec。

证明：分两 case。
- **Case 1: F empty**：循环不执行，c=0。`if c==0 return 1`，与 Spec 一致。✓
- **Case 2: F 非空**，n 个项 c_0, ..., c_{n-1}：
  - 由 Lemma 1（k=n）：循环结束时 `c = gcd(|c_0|, ..., |c_{n-1}|)`。
  - 此值非负（gcd 非负）。
  - 子 case 2a: c == 0。意味着所有 `|c_i| = 0` 即所有 c_i = 0。
    - 在 SparsePoly 表示中，coef=0 的项应被 normalization 剔除（参 `SparsePolyZZ.normalization`）。
    - 故此 case 仅在 F 含被规范化遗留的 0-coef 项时发生（异常态）。返回 1 是 graceful fallback。
  - 子 case 2b: c > 0。直接返回 c。与 Spec 完全一致。✓

故实现 = Spec。 ∎

#### 与原实现的等价性（多项情形）

证明修复**不破坏**多项情形的现有行为：

*Claim*：对 F 含多个项（`F.size() ≥ 2`），原实现与新实现返回值相等。

证明：
- **原实现** 第一项后：`c_orig = c_0`，可能为负。第二项后：`c_orig = gcd(c_0, c_1) = gcd(|c_0|, |c_1|) = |c_0_or_minus|`（GMP 约定取绝对值）。后续 k 步：`c_orig = gcd(|c_0|, ..., |c_{k}|)`。终局 `c_orig = gcd(|c_0|, ..., |c_{n-1}|)`。
- **新实现** 由 Lemma 1：终局 `c_new = gcd(|c_0|, ..., |c_{n-1}|)`。
- `c_orig = c_new`。✓ □

故修复对多项情形无副作用，仅改单项情形（原实现返回 c_0 带 sign，新实现返回 |c_0| = gcd(0, c_0)）。

#### 边缘情形枚举

| F | 原实现 | 新实现 | Spec |
|---|--------|--------|------|
| 空 | 1（特殊分支）| 1（c=0 兜底）| 1 |
| 单项 `[c_0]`，c_0 > 0 | c_0 | gcd(0, c_0) = c_0 | c_0 |
| 单项 `[c_0]`，c_0 < 0 | **c_0（负）✗** | gcd(0, c_0) = -c_0 | -c_0 = |c_0| |
| 单项 `[0]`（异常态）| 0 | 1（c=0 兜底）| 1（约定）|
| 多项 | gcd(|c_i|) | 同（by Claim） | 同 |

修复点确认：行为改变涵盖两类（均为单项情形）：
1. **单项 c_0 < 0**：返 `|c_0|` 而非 `c_0`（主要修复目标，含常数 `-N` 情形）
2. **单项 c_0 = 0**（异常态）：返 1 而非 0（normalization 应剔零项，此情形不应正常出现；新行为更接近 Spec）

多项情形和正系数单项情形行为完全保持。

### 改动 2: `clpoly/upolynomial.hh:220`（template 版本）

同样的修复模式：
```cpp
template <class Tc>
Tc cont(const upolynomial_<Tc> & F) {
    Tc c = 0;
    for (auto &t : F) {
        c = gcd(c, t.second);
    }
    if (c == 0) return 1;
    return c;
}
```

原代码 L225 `auto I=(ptr++).second;` 的 `.second` 看起来应该是 `->second`（iterator 解引用）。但由于 grep 确认本 template 未被实例化（无 `cont(upolynomial_<Zp>)` 等非-ZZ 调用），这处潜在编译错误并未触发。新实现用 range-for（`for (auto &t : F)`）规避此问题。

### 改动 3: 多变量 `cont`（polynomial_gcd.hh:581）

**本 PR 同步修**。和单变量是同一个 bug（cont sign 不遵循 Maple/SymPy/FLINT 的"always positive"约定），分 PR 没有真实风险隔离意义。

#### 实测 (cont multivar 行为)

- `cont(-x mvpoly) = -1` ← single-group **同样 sign bug**
- `cont(-2xy+3x mvpoly) = 1`（实测 polynomial_GCD 路径有时归正到整数 cont）
- `cont(-2x+1 mvpoly) = 1`（multi-group 时 `polynomial_GCD(cont, tmp)` 自动折正）

#### 根因（控制流分析）

```cpp
for (auto &i : F_) {
    if (匹配当前组) {
        tmp.push_back({..., i.second});  // 直接 push，coef 带原 sign
    } else {
        if (tmp_deg == deg) {
            cont = std::move(tmp);       // ① first group：直接 move，无 normalize
        } else {
            cont = polynomial_GCD(cont, tmp);  // ② 后续 group：polynomial_GCD 自动归正
        }
        ...
    }
}
// 收尾同样两条路径
if (tmp_deg == deg) cont = std::move(tmp);     // ③ 同 ①
else cont = polynomial_GCD(cont, tmp);          // ④ 同 ②
return cont;
```

- 路径 ②④（polynomial_GCD）：单项分支 `c = gcd(c0, ...)` 走 GMP `mpz_gcd` 取绝对值返非负 → 返回首项正；多项 Brown's 模算法路径理论上也归正首项（标准约定）✓
- 路径 ①③（std::move）：直接传 tmp 不归正 ✗（single-group 输入 only 走这条）

**注**：path ②④ "polynomial_GCD 自归正" 假设依赖标准 Brown 算法约定，本 PR 未深入证明该 invariant。若 polynomial_GCD 在某 corner case 未归正，本 PR 的末尾 normalize 兜底仍能修复（normalize 对已归正情形 no-op）。

#### 修复

在 return 前加 3 行 sign normalize（cover 所有路径，对已归正情形是 no-op）：

```cpp
// ... 收尾后:
if (tmp_deg == deg) cont = std::move(tmp);
else cont = polynomial_GCD(cont, tmp);

// === 新增 sign 归一化 ===
if (!cont.empty() && cont.front().second < ZZ(0)) {
    for (auto& term : cont) term.second = -term.second;
}

return cont;
```

#### 正确性半形式化证明（多变量补充）

**Spec**：多变量 `cont(F)` 关于主变量 v 应返回剩余变量上的"canonical GCD 多项式"（约定：首项系数 > 0，对齐 Maple/SymPy/FLINT）。

**实现 invariant**：
- 在 return 前的 normalize 步骤后，`cont.empty() || cont.front().second > 0`。

证明：
- 若 cont.empty()：直接返回 ✓
- 若 cont.front().second > 0：if 条件 false，不进入 normalize，保持 ✓
- 若 cont.front().second < 0：进入 normalize，所有 coef 取反 → 新 front coef > 0 ✓
- cont.front().second == 0：normalization 不变（不应发生，因为 normalization 应剔零 coef）

**等价性**（修后不破坏多项情形）：
- 多项情形（path ②④）：polynomial_GCD 已归正，cont.front().second > 0 入口 → normalize 是 no-op。
- 单项情形（path ①③）：tmp 可能带负 leading；normalize 后归正。✓

**结合 univariate cont 修复**（改动 1）：两者共同保证 cont 永远返回 canonical form（univariate 返回非负整数；multivariate 返回首项 > 0 多项式）。

### 不需要改

- `squarefreefactorize` 的 sign-normalization（PR #16）—— 保留，仍是 robust 兜底
- `polynomial_GCD` 的逻辑 —— 已自归正首项
- 各 test 文件 —— 见 §4 调用点扫描决定是否更新 expected 值

---

## 第四部分：调用点 sign-敏感性扫描

实测 grep `clpoly/*.hh test/*.cc` 找到所有 `cont(` 调用：

| 位置 | 上下文 | 对 cont sign 敏感? |
|------|--------|-------------------|
| `polynomial_factorize_univar.hh:717` | `ZZ c = cont(f); ... f / c ...` | **可能** —— 需 verify |
| `polynomial_gcd.hh:88` | sqf 内部（PR #16 sign-normalize 层） | ❌（兜底层会处理）|
| `polynomial_gcd.hh:111` | sqf 内部（PR #16） | ❌ |
| `polynomial_gcd.hh:299/303` | polynomial_GCD 主变量同步 | **可能** —— 需 verify |
| `polynomial_gcd.hh:306-308` | polynomial_GCD content 提取 | ❌（接着 F=F/F_cont 都重新算）|
| `polynomial_gcd.hh:429` | polynomial_GCD CRT 内 | ❌（同上）|
| `polynomial_gcd.hh:700` | pp 函数（pp = f / cont(f)） | **可能** —— 需 verify |
| `polynomial_gcd.hh:1067` | univariate ZZ cont 定义 | 修复目标 1 |
| `polynomial_gcd.hh:581` | **多变量 cont 定义** | **修复目标 3**（新加 sign normalize）|
| `upolynomial.hh:220` | template cont 定义 | 修复目标 2 |
| `charset.hh:59` (sqrfree) | 用 sqf 输出，间接受影响 | ❌ |
| 测试文件 | 各种 assertion | **需 verify 测试断言是否依赖 sign** |

### 详细 verify

**`pp(f) = f / cont(f)`**（polynomial_gcd.hh:700）:
- 修前：单项 `-2x`, cont=-2, pp = -2x/-2 = x。结果 x（正）。
- 修后：cont=2, pp = -2x/2 = -x。结果 -x（负）。
- **行为改变** —— pp 返回带 sign 的 primitive part vs 正首项 primitive part。
- 这是 issue #14 根因层！PR #16 在 sqf 内做 sign-normalize 处理这种 case。
- 调用 pp 的下游若依赖"pp 总返正首项"会破。需 grep 验证。

**`F = cont(F)`**（polynomial_gcd.hh:299/303）:
- 这是 polynomial_GCD 内 reduce variable count 的步骤
- F 变成它的 cont（多变量）
- cont 本身的 sign 由内部 polynomial_GCD 递归决定
- 修复后多变量 cont 调用方式不变，仍正确

**测试断言**：
- `test_polynomial_gcd.cc:177,185`、`test_upolynomial.cc:93`、`test_multivar_helpers.cc:30,43,51`、`test_factorize_trace.cc:42`
- 需各跑一遍看修后 expected 值是否变化

---

## 第五部分：测试计划

### 主回归用例

新增 `test/test_cont_sign_consistency.cc`：

| 用例 | 输入 | 期望（修后）|
|------|------|----------|
| `cont_pos_single` | `x` | 1 |
| `cont_neg_single_x` | `-x` | **1**（修前 -1）|
| `cont_neg_single_2x` | `-2*x` | **2**（修前 -2）|
| `cont_neg_multi_minus_x_plus_1` | `-x+1` | 1 |
| `cont_neg_multi_minus_2x_plus_2` | `-2*x+2` | 2 |
| `cont_constant_negative` | `-6`（const poly）| **6**（修前 -6 — 常数 poly 是单项，sign bug 适用）|
| `cont_zero_poly` | `0` | 1（约定）|

### pp 联动测试

| 用例 | 修前 `pp(f)` | 修后 `pp(f)` |
|------|------------|------------|
| `pp(-x)` | x | **-x** |
| `pp(-2x+2)` | x-1 | **-x+1**（首项负） |

### 多变量 cont 回归用例

新增到同一测试文件：

| 用例 | 输入 | 修前 | 修后 |
|------|------|------|------|
| `mv_cont_neg_single` | `-x` (mvpoly) | -1 | **1**（const poly）|
| `mv_cont_neg_term_x` | `-2*x*y` (single x-deg group) | -2y | **2y**（首项 y 系数 > 0）|
| `mv_cont_mixed_groups` | `-2*x+1` | 1 | 1（多组本已正确）|
| `mv_cont_x_in_groups` | `-x*y + x + 1` | 1 | 1（同上）|

### 现有测试回归

- `test_polynomial_gcd.cc`：可能需 update 期望值
- `test_upolynomial.cc:93`：同上
- `test_multivar_helpers.cc`：同上
- `test_factorize_trace.cc`：仅 trace，应不影响

### sqf / squarefreebasis 回归

- `test_squarefreebasis_bugfix.cc`（PR #16 加的）：应继续通过，因为 sign-normalize 层修后仍处理 pp 可能返负首项的情况
- `test_factorize*`：因式分解链路应不变（cont sign 改变不影响因子分解结果）

### 全量回归

`bash test/run_all_tests.sh` 全绿。

---

## 第六部分：执行步骤（待用户批准）

1. **[等批准]** 用户确认修复方案
2. 实施 §3 三处改动
3. 新增 `test/test_cont_sign_consistency.cc`
4. 跑新测试 + 全量回归
5. **关键**：检查 `pp` 的调用方在 "pp 可能返负首项" 后是否仍正确（重点：polynomial_gcd.hh 里的 pp 用例 + 测试）
6. 如发现下游依赖 pp 正首项，加 sign-normalize 兜底（在 pp 后）
7. devlog + commit + PR 关联 issue #17

---

## 第七部分：风险与未决问题

1. **Q**: `pp` 在 cont 修后会返回带 sign 的 primitive part（首项可能负）。下游会破吗？
   **A**: PR #16 在 sqf 已加 sign-normalize，应能处理。但其他 pp 调用方需 grep verify。

2. **Q**: `polynomial_GCD` 内部多处用 cont/pp，sign 改变会导致正确性问题吗？
   **A**: polynomial_GCD 主要用 cont/pp 做 content extraction（数学上 sign 不重要，因为最终又乘回 cont_gcd），不应受影响。但需测试 verify。

3. **Q**: 同步修 Lean Model.lean 端？
   **A**: 独立 PR 后，在 feature/formal-proofs 分支同步修 Lean 端（删除 contImpl 的 sign 处理）。本 PR 仅 C++。
