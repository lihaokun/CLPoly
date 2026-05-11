# 修正方案：squarefreebasis unit leak

> **状态**：草稿，待用户审核。
> 对应 workflow.md §5.1 算法内部错误修正流程
> 对应 GitHub issue: [#14](https://github.com/lihaokun/CLPoly/issues/14)
> 对应分支：`fix/squarefreebasis-issue-14`（从 master 起，独立 PR 对 master）

---

## 第一部分：复现与定位

### 现象

issue #14：`squarefreebasis` 在输出 basis 中混入常数 `-1`，且数量随输入规模增长。

```cpp
variable x("x");
auto [b1, _] = squarefreebasis({-x+1, (x-1)*x});
// 实际：[-1, x-1, x]
// 期望：[x-1, x] 或其 associates 等价集合，不应有 unit 元素
```

### 定位

`clpoly/polynomial_gcd.hh:158-211` 的 `squarefreebasis` 主循环 L176-195：

```cpp
for (size_t j=0;j<lst.size();++j) {
    auto F1 = polynomial_GCD(lst[j], tmp);
    if (!is_number(F1)) {
        tmp = tmp / F1;
        if (F1 != lst[j]) {              // ← 字面比较，漏判 associates
            lst[j] = lst[j] / F1;        // ← unit 残留进 basis
            ...
        } else { ... }
        ...
    }
}
```

### 实测追踪（用 issue 例 1 `{-x+1, (x-1)·x}`）

通过 `/tmp/probe_sqf.cc` 实测确认：

```
sqf(-x+1) = [(1, 1), (-x+1, 1)]    ← factor 首项负，非 canonical
polynomial_GCD(-x+1, x²-x) = x-1   ← gcd 内部归到正首项
(-x+1) / (x-1) = -1                ← associates 比例为 unit
```

**bug 路径**：

| step | 状态 |
|------|------|
| F[0]=-x+1 | sqf 跳过 cont=1 项 → `lst = [-x+1]`（非 canonical） |
| F[1]=(x-1)·x | sqf 跳过 cont=1 项 → tmp=x²-x |
| inner j=0 | F1 = gcd(-x+1, x²-x) = **x-1**（gcd 归正） |
| | F1 (x-1) **!=** lst[0] (-x+1) → 走"拆分 basis"分支 |
| | lst[0] = lst[0]/F1 = -1 ← **unit leak** |
| 最终 | `lst = [-1, x-1, x]` |

### 根因链

1. `polynomial_GCD` 内部规范化首项系数为正
2. `squarefreefactorize` 的 factor 不规范化（仅 cont 取 gcd 绝对值）
3. → 同一个数学多项式（associates）在 sqf 和 polynomial_GCD 之间的字面表示不一致
4. → sqfb 的 `F1 != lst[j]` 字面比较误判为"应拆分"
5. → 拆分时 `lst[j]/F1` 是 unit 残留进 basis

---

## 第二部分：参考实现对比

按 §5.1 必须对照参考实现：

### SymPy 实测（`sqf_list`）

```python
sqf_list(x)              = (1,  [(x, 1)])
sqf_list(-x)             = (-1, [(x, 1)])
sqf_list(1-x)            = (-1, [(x-1, 1)])           ← cont=-1，factor canonical
sqf_list(-2*x)           = (-2, [(x, 1)])
sqf_list(-6*(x-1)*(x-2)) = (-6, [(x²-3x+2, 1)])      ← cont=-6，factor canonical
sqf_list((x-1)**2*x)     = (1,  [(x, 1), (x-1, 2)])
```

**SymPy 约定**：cont 完全吸收 unit/sign，factors 永远首项正。

### Maple `sqrfree`

返回 `[lc, [[f1,e1],[f2,e2],...]]`，`lc` 是首项 unit（含 sign），`fi` 是 monic（首项=1）。

**Maple 约定**：等价于 SymPy，sign 全提到 lc，factors 永远 canonical。

### FLINT `fmpz_poly_factor_squarefree`

输出结构 `fmpz_poly_factor_t {c, p[]}`：
- `c`：内容部分（带 sign）
- `p[]`：每个 factor 经 `fmpz_poly_factor_canonicalise` 归一为 cont=1 + 首项正

**FLINT 约定**：同 SymPy/Maple。

### 共同模式

**三家完全一致**：cont/lc 吸收所有 unit/sign，factor list 内部多项式永远 canonical（首项正 over ZZ，或 monic over field）。

### CLPoly 当前 vs 业界

CLPoly 实测（`/tmp/probe_sqf.cc`）：

| 输入 | CLPoly 当前 | SymPy/Maple/FLINT |
|------|------------|-------------------|
| `x` | `[(1,1), (x,1)]` | 一致 ✓ |
| `-x` | `[(-1,1), (x,1)]` | 一致 ✓ |
| `x-1` | `[(1,1), (x-1,1)]` | 一致 ✓ |
| **`-x+1`** | **`[(1,1), (-x+1,1)]`** | **`[(-1,1), (x-1,1)]`** ✗ 不一致 |
| `2x`, `-2x` | 一致 ✓ | — |
| **`-6(x-1)(x-2)`** | **`[(6,1), (-x²+3x-2,1)]`** | **`[(-6,1), (x²-3x+2,1)]`** ✗ 不一致 |
| 重复因子各种 | 一致 ✓ | — |

**差异定位**：仅在「F 是多项 + 首项系数为负」的情形，CLPoly 没把首项 sign 提到 cont。其余情形与业界完全一致。

---

## 第三部分：修复方案（Plan B：上游规范化）

让 `squarefreefactorize` 在算出 `f_cont` 和 `F_` 之后多做一步：把 `F_` 的首项 sign 提到 `f_cont`，使 `F_` 永远首项正。

```cpp
// clpoly/polynomial_gcd.hh:109-111 修改：
auto f_cont = cont(F);
auto F_ = F / f_cont;
// === Plan B: 把 F_ 的首项 sign 提到 f_cont ===
// 对齐 SymPy/Maple/FLINT 约定：factor 永远 canonical（首项正）
if (!F_.empty() && F_.front().second < ZZ(0)) {
    F_ = -F_;           // polynomial unary -
    f_cont = -f_cont;   // cont 也跟着取负，保持 cont * F_ = F 不变
}
auto lst = squarefreefactorize(f_cont);  // 递归 sqf cont 自身保持 sign 一致
```

### 为什么这一处改动够了

1. `F_` 初始首项正 ⇒ 所有由 `F_` 出发的 GCD/quotient 都自动首项正（`polynomial_GCD` 规范化）
2. 循环里 `lst.push_back({F_3/F_, ...})` 等表达式：分子分母均首项正 ⇒ 商首项正
3. 所有 push 进 lst 的 factor 都 canonical
4. → sqfb 拿到的 lst[j] 永远 canonical
5. → `F1 != lst[j]` 字面比较准确，不再有 associate 误判
6. → unit 不会进 basis

### 不需要改的代码

- `squarefreebasis` 本身：不动
- `cont` 的 sign convention：不动（仍是 ZZ.gcd 的 sign 行为）
- 任何调用 `squarefreefactorize` 的代码：不动（实测 7 处调用点全部对 factor sign 不敏感）

### 调用点 sign-敏感性扫描（实测）

| 文件:行 | 用法 | 敏感? |
|---------|------|------|
| `charset.hh:53,67` | `G = G * factor` 乘积重建 | ❌（乘积语义保持等价）|
| `polynomial_factorize.hh:69` | `is_number` 检查 + 递归 factor | ❌ |
| `polynomial_factorize_wang.hh:2450` | 同上多变量入口 | ❌ |
| `polynomial_factorize_univar.hh:1596` | 已有独立 sign-normalize 层（L1618-1621） | ❌ |
| `polynomial_gcd.hh:111` | sqf 内部递归 cont | ❌（递归 sign 自动跟随）|
| `polynomial_gcd.hh:169` | **squarefreebasis bug 触发点** | 修后自然正确 |
| `test/test_polynomial_gcd.cc:132,147` | `divides` + `is_squarefree` 断言 | ❌（不依赖 sign）|
| `test/bench_clpoly.cc:482-485` | 仅 benchmark，丢 result | ❌ |

**结论**：所有调用方不需修改。

---

## 第四部分：测试计划

### 主回归用例（issue #14 复现）

新增 `test/test_squarefreebasis_bugfix.cc`：

| 用例名 | 输入 | 期望（必须无 `is_number` 元素）|
|--------|------|------------------------------|
| `issue14_ex1` | `{-x+1, (x-1)·x}` | basis 中无 unit |
| `issue14_ex2` | 4 polys：`{-x+1, (x-1)·x, -x+2, (x-2)·x}` | 同上 |
| `issue14_ex3` | 6 polys：含 `-x+3, (x-3)·x` | 同上 |

### sqf canonical 形式断言

| 用例名 | 输入 | 期望 |
|--------|------|------|
| `sqf_neg_lead_canonical` | `-x+1` | `sqf` 返回的非常数 factor 首项系数 > 0 |
| `sqf_neg_multivar` | `-(x-1)·(x-2)` | 同上 |
| `sqf_product_invariant` | 任意 F | `∏ factor^mult = F` 仍成立 |

### Sanity（不破坏正向用例）

| 用例名 | 输入 | 期望 |
|--------|------|------|
| `sanity_pos_lead` | `{x+1, (x+1)·x}` | basis = `[x+1, x]`，与修复前同 |
| `sanity_repeated_factor` | `{(x-1)², x·(x-1)}` | multiplicity 正确 |
| `sanity_coprime` | `{x, x+1}` | basis = `[x, x+1]` |

### 全量回归

- `bash test/run_all_tests.sh`：所有现有测试不破坏
- 重点关注：`test_polynomial_gcd.cc`（直接测 sqf）、`test_factorize.cc`、`test_factorize_multivar.cc`（间接通过 factorize 链路）

---

## 第五部分：执行步骤

1. **[等批准]** 用户确认本方案
2. 在 `fix/squarefreebasis-issue-14` 分支：
   1. 改 `clpoly/polynomial_gcd.hh:109-111` 加 sign 归一化（约 4 行）
   2. 新增 `test/test_squarefreebasis_bugfix.cc`（约 80 行，含表中所有用例）
   3. 编译：`make test/test_squarefreebasis_bugfix`
   4. 跑新测试：期望 6/6 PASS
   5. 全量回归：`bash test/run_all_tests.sh`，期望全绿
3. commit + push 分支
4. 开 PR 对 master，关联 issue #14
5. 写 devlog `docs/devlog/2026-05-XX-squarefreebasis-unit-leak-fix.md`

---

## 第六部分：风险与未决问题

1. **Q**: `cont(-x) = -1`（单项情形 cont 已带 sign）此修复是否冲突？
   **A**: 不冲突。修复后 `F_ = F/cont(F)` 在单项情形得到 `1`，`!is_number(F_)` 为 false，跳过 sign-normalize 步。`is_number(F)` 检查在 L107 已经 return，不会进入 L109。

2. **Q**: `polynomial_<ZZ>` 是否有 unary `-` 算子？
   **A**: 需要在实施时验证。若无，备选：`F_ = F_ * ZZ(-1)` 或对每个 coef 取反的循环。

3. **Q**: Plan B 改完后，对原本"首项正"的输入完全无影响？
   **A**: 是。修复仅触发在 `F_.front().second < 0` 时；首项正的 F 完全走原路径。

4. **Q**: cont 的 sign 在 sqf 输出"第一项 (cont, 1)"中的语义如何？
   **A**: 保持现有：第一项是 cont 多项式（递归 sqf 后可能拆成多个整数 factor）。修复后 cont 多项式可能整体取负（吸收了 sign），但形式不变。
