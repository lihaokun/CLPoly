# CLPoly Lean 4 形式化验证

## 项目定位

CLPoly 因式分解模块的 Lean 4 机器检查证明。验证目标是 C++ 实现的数学正确性，不是重写算法。

依赖：Mathlib4 (stable)，Lean 4.28.0+

## 构建指令

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd proof/lean
lake build              # 编译全部
lake env lean FILE.lean # 单文件检查
```

## 验证架构：三层自顶向下

```
L3  数学基石    CLPoly/Math/         Mathlib 定理组合，纯数学定义与性质
L2  算法模型    CLPoly/Algorithm/    CLPoly 算法逻辑，抽象掉实现细节，证明满足 L3 Spec
L1  实现模型    CLPoly/Impl/         1:1 对应 C++（uint64 语义、数组越界、move），精化 L2
```

**工作顺序**：Spec（接口规约）→ Pipeline（顶层骨架）→ Math（L3）→ Algorithm（L2）→ Impl（L1）

详见 `proof/docs/implementation-roadmap.md`

## 证明工作流

对每个需要证明的定理或 sorry，必须先写**自然语言证明草稿**（关键步骤、依赖引理、Mathlib 路径），确认后再形式化为 Lean 代码。禁止跳过草稿直接写 tactic proof。

## 开发日志规范

每完成一个 Phase 或一组定理后，在 `docs/devlog/` 下写开发日志。除常规内容（做了什么、为什么、关键决策）外，**必须包含 `## 度量` 段**：

```markdown
## 度量
- 耗时：~X 小时（含调研/草稿/形式化/调试）
- 迭代：Y 轮编译-修复循环（每轮 = 一次 lake env lean 到下一次修改）
- Lean 新增/修改行数：Z 行
- 对应 C++ 行数：W 行（如适用，指被验证的 C++ 代码量）
- 放弃的方案：（列出尝试过但未采用的方案及原因，若无则写"无"）
```

此数据用于未来撰写论文时的量化分析（验证成本、放大系数、方法论评估）。

## 编码约定

### Lean 风格
- `autoImplicit = false`，所有变量显式声明
- 使用 `Nat` 而非 `ℕ`（除非在 Mathlib import 作用域内）
- `noncomputable` 按需添加（`EuclideanDomain`、`GCDMonoid` 等必须标注）
- 素数前置：`variable (p : ℕ) [Fact (Nat.Prime p)]`

### 文件组织
- `CLPoly/Spec.lean` — 子过程接口规约（Prop 定义）
- `CLPoly/Pipeline/` — 顶层正确性定理
- `CLPoly/Math/` — L3 纯数学定理
- `CLPoly/Algorithm/` — L2 算法模型（抽象算法逻辑）
- `CLPoly/Impl/` — L1 实现模型（1:1 对应 C++）
- `CLPoly/Experiment/` — Phase 0 实验（保留参考）

### 命名
- 规约 Prop：`SquarefreeDecomp`、`DDFCorrect`、`EDFCorrect`
- 数学定理：`irreducible_dvd_X_pow_sub_X`、`gcd_X_pow_sub_X_factors`
- L2 算法函数：`ddfLoop`、`squarefree_Zp`（对应 CLPoly 算法）
- L1 实现函数：`ddfLoop_impl`、`subtract_x_impl`（1:1 对应 C++ 函数名 + 控制流）

### L2 算法模型原则
算法模型捕获 CLPoly 的算法逻辑（DDF 循环、EDF 随机分裂、Hensel 提升），但抽象掉 C++ 实现细节（整数表示、数组布局、内存管理）。操作对象是 Mathlib 的 `Polynomial (ZMod p)` 等数学类型。

### L1 实现模型：1:1 对应原则
实现模型的控制流必须与 C++ 一一对应，不允许重写为"更优雅的 Lean 风格"。操作对象是 `U64`、`Vec` 等 C++ 语义模型。精化证明证明 L1 行为与 L2 算法一致。

| C++ 构造 | L2 算法模型 | L1 实现模型 |
|---------|-----------|-----------|
| `for` 循环 | 数学归纳 / 递归 | 顶层递归 + `termination_by`（1:1 控制流） |
| `while(true)` | 存在性论证 | `partial def` 或 fuel 参数 |
| `polynomial_GCD(a, b)` | `GCDMonoid.gcd a b` | 调用 Euclid GCD 的 L1 模型 |
| `(int64_t)(p - 1)` | `p - 1`（自然数） | `cast_u64_to_i64 (p - 1)`（显式溢出语义） |
| `vec[i]` | `list.get i` | `Vec.get i (proof : i < size)` |
| `std::move(v)` | 无对应（纯函数式） | `Ownership.move_from v (proof : ¬moved)` |

## Mathlib API 速查

| 需要 | 正确路径 |
|------|---------|
| ZMod p 是域 | `instance : Fact (Nat.Prime p)` 前置 |
| gcd | `GCDMonoid.gcd` 或 `EuclideanDomain.gcd` |
| 域上多项式除法 | `divByMonic` / `modByMonic`（`Ring` 足够） |
| 恒等式 | `modByMonic_add_div : p %ₘ q + q * (p /ₘ q) = p` |
| squarefree ↔ nodup | `UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors` |
| 有限域根 | `roots_X_pow_card_sub_X` |
| 域嵌入 | `nonempty_algHom_iff_finrank_dvd` |
| 层间投影 | `ZMod.castHom` (需 `dvd` 证明) |

## 已知陷阱

- `Polynomial.gcd` 不存在，用 `GCDMonoid.gcd`
- `EuclideanDomain` / `GCDMonoid` 是 noncomputable，不影响证明但阻碍 `#eval`
- `autoImplicit = false` 下 `ℕ` 需要 import Mathlib 才可用
- `let rec` 内的 `termination_by` 语法受限，改用顶层递归函数
- `modByMonic_add_div` 的顺序是 `f %ₘ g + g * (f /ₘ g) = f`，注意加法顺序

## Lean 形式化经验教训（Phase 3 总结）

### normalize 和 Associated 的系统性摩擦

`EuclideanDomain.gcd` 返回的不是 monic 多项式。`normalize` 使其 monic，但引入 `Associated`（非精确等式）。**所有涉及 gcd 的 dvd/乘积证明都需要 `normalize_dvd_iff`、`normalize_associated`、`Associated.trans` 等桥接**。这是 SQF 证明最大的时间消耗源。

**规则**：
- 精确除法用 `modByMonic_add_div`（需要 `Monic` 前提）
- `normalize a ∣ b ↔ a ∣ b`：用 `normalize_dvd_iff`
- `Associated (normalize a) a`：用 `normalize_associated`（注意方向！`.dvd` 给 `normalize a ∣ a`）
- **绝不对含 `sqfZp`/`yunLoop` 递归调用的目标用 `rw [h]`**——会无限展开或替换内部的 `f`。用 isolated `have` + `rw` 代替

### `let (a, b) := e` 不产生 definitional equality

`let (a, b) := e` 编译为 `match`，`a` 和 `e.1` 不是 definitionally equal。**在递归函数中，始终用 `.1`/`.2` 显式投影**：
```lean
-- 错误：let (result, c_rem) := yunLoop w c 1 [] hc
-- 正确：
let output := yunLoop w c 1 [] hc
let result := output.1
let c_rem := output.2
```

### `decreasing_by` 中的 bound 应放入函数体

若终止性需要中间结果的 bound（如 `yunLoop` 输出的 `c_rem.natDegree ≤ c.natDegree`），在函数体中 `have hbound := ...` 比在 `decreasing_by` 中重新推导更可靠。`decreasing_by` 的 context 受 let/match 影响，变量名不可预测。

### functional induction (`f.induct`) vs `Nat.strongRecOn`

- `f.induct`：自动匹配函数分支，但 inaccessible name 数量不可预测（let 绑定也算），需要 `rename_i` 逐个试
- `Nat.strongRecOn`：完全掌控命名，但需要手动 `rw [f]` + `split` 展开函数。**推荐用于不变量较多的递归证明**

### `rename_i` 计数经验

| 函数结构 | inaccessible 数量 |
|---------|-----------------|
| 4 个参数 + 1 个 `if` guard | 5 |
| 4 个参数 + 1 guard + 2 let + 1 inner guard | 8-9 |
| 加 `hc : Prop` 参数 | +1 |
| 有 `let f_new :=` 在 split 内 | +1（仅该分支） |

### 非线性 omega 的处理

`omega` 不处理 `a * b` 形式（非线性）。`a < a * p`（`a ≥ 1, p ≥ 2`）需要：
```lean
have key := Nat.mul_lt_mul_of_pos_left (show 1 < p from hp.out.one_lt) ha_pos
simp only [Nat.mul_one] at key
```

### 严格按 nl-proof 翻译，遇到问题回到 nl-proof 修正

**这是最重要的教训。** Lean 翻译卡住时：
1. **第一反应**：检查 nl-proof 该步是否正确（不是绕路/换策略/加 sorry）
2. 发现错误 → **停止编码** → 修正 nl-proof → 审核确认 → 再翻译
3. nl-proof 中的"同某个已有证明" → **必须验证递归结构完全一致**，不能假设
4. 每个 sorry 解决后再开下一个。不要并行 sorry（sorry 越多越慢）

**反面案例**：SQF 的 `derivative(c_rem) = 0` 花了 8 小时，因为 nl-proof §3.2.1 的"无穷递升"论证有数学错误（步进量 = 0，不是 +1），但我没有回去检查 nl-proof，而是反复尝试绕路，导致代码膨胀 + sorry 增殖。

**正面案例**：修正 nl-proof 后，严格按 steps 2a-2i 逐步翻译，6 轮编译-修复即完成全部 4 sorry（~2 小时）。**nl-proof 正确 + 严格翻译 = 高效形式化**。

### 优先攻克关键瓶颈

nl-proof 完成后，**先识别 Lean 形式化的关键瓶颈**（通常是最深的数学引理 + 最复杂的 Associated 链），从那里开始，而不是自底向上。DDF 成功（0 sorry）因为结构简单；SQF 的 `derivative(c_rem) = 0` 和 Associated 拼接应该最先攻克。

### `set`-bound 变量与 `rw` 的冲突

`set crem := (yunLoop ...).2` 后，`crem` 在 context 中是定义式别名。**`rw [h]` 找不到 `crem` 的展开形式**（Lean 只匹配语法，不展开 `set` 定义）。

**解法**：
- 用 `Associated.of_eq h` 建立 `Associated crem (...)` 桥接，而非 `rw`
- 将需要跨 goal 使用的 `Associated` 提到 `refine ⟨?_, ...⟩` **之前**（否则只在一个 goal 内可见）
- `set`-bound 变量出现在引理参数中时，提取为独立引理（clean 参数签名）

### `Nat.find` 需要 `DecidablePred`

多项式整除性 `q ^ k ∣ f` 在 Lean 中默认非 decidable。使用 `Nat.find` 时需显式提供：
```lean
set v := @Nat.find (fun n => ¬(q ^ (n + 1) ∣ crem)) (Classical.decPred _) hfin
```
`Nat.find_spec` 和 `Nat.find_min` 同理需要 `@` + `Classical.decPred _`。

### Lean 4 Mathlib List API 变化

- `List.mem_cons_self` 不是函数，是 `a ∈ a :: l` 的证明。**用 `by simp` 替代**
- `List.mem_map_of_mem` 同理。**用 `List.mem_map.mpr ⟨x, h, rfl⟩` 替代**
- `List.dvd_prod` 接受 `h : a ∈ l` 而非 `List.mem_map_of_mem` 的结果

### `derivative_pow` 后的简化

`Polynomial.derivative_pow` 展开为 `C ↑n * f ^ (n-1) * derivative f`。之后 `mul_zero` / `zero_mul` 可能因乘积结结构不匹配而失败。**用 `ring` 替代逐步 `rw`**：
```lean
-- 错误：rw [derivative_pow, hq_zero, mul_zero, mul_zero]  -- 可能模式不匹配
-- 正确：rw [derivative_pow, hq_zero]; ring
```

### `IsCoprime.pow_left` 参数名

Mathlib 中参数名是 `m`，不是 `n`：
```lean
-- 错误：hcop.pow_left (n := v + 1)
-- 正确：hcop.pow_left (m := v + 1)
```
查 API 时注意命名参数（尤其是 `IsCoprime` 系列）。

### 非线性 omega 的补充模式

除了 `a < a * p` 外，`(n+1) * deg(q) ≤ n` 也需要辅助：
```lean
have h3 : crem.natDegree + 1 ≤ (crem.natDegree + 1) * q.natDegree :=
  Nat.le_mul_of_pos_right _ (Irreducible.natDegree_pos hq_irr)
omega
```

### Mathlib API 补充速查

| 需要 | 正确路径 |
|------|---------|
| `normalize a ∣ b ↔ a ∣ b` | `normalize_dvd_iff` |
| `a ∣ normalize b ↔ a ∣ b` | `dvd_normalize_iff` |
| `Associated (normalize a) a` | `normalize_associated` |
| `degree(normalize f) = degree(f)` | `degree_eq_degree_of_associated (normalize_associated _)` |
| `natDegree(normalize f) = natDegree(f)` | 自定义 `natDegree_normalize_eq`（Mathlib 无） |
| `expand p f = f^p` in F_p[X] | `map_frobenius_expand` + `frobenius (ZMod p) p = id`（`ZMod.pow_card`） |
| `expand_contract` (f' = 0 → f = expand p (contract p f)) | `@expand_contract _ _ p _ _ f hderiv hp.ne_zero` |
| `Squarefree (a*b) → IsRelPrime a b` | `squarefree_mul_iff` + `.1` |
| `IsRelPrime → IsCoprime` | `IsRelPrime.isCoprime`（需 PID/Bezout） |
| `IsCoprime a b, c ∣ b → IsCoprime a c` | `IsCoprime.of_isCoprime_of_dvd_right` |
| 非零非 unit 多项式有不可约因子 | `WfDvdMonoid.exists_irreducible_factor` |
| `d^k ∣ f → d^{k-1} ∣ f'` | 自定义 `pow_dvd_derivative_of_pow_succ_dvd`（`derivative_pow` + `dvd_add`） |
| `a = b → Associated a b` | `Associated.of_eq` |
| `IsCoprime a b, c ∣ b → IsCoprime a c` | `IsCoprime.of_isCoprime_of_dvd_right` |
| `IsCoprime a b → IsCoprime (a^n) b` | `IsCoprime.pow_left (m := n)` |
| `∃ min n with P n`（非 decidable P） | `@Nat.find _ (Classical.decPred _) hfin` |
| `f / gcd(f, f')` squarefree | 自定义 `squarefree_div_gcd_derivative`（无穷递升 + 度数矛盾） |
| `q^v ∣ f, q' ≠ 0, p∤v → q^v ∤ f'` | 自定义 `not_pow_dvd_derivative_of_separable` |

## 参考文档

- `proof/docs/implementation-roadmap.md` — 实施路线（自顶向下）
- `proof/docs/research-report.md` — Mathlib API 调研
- `proof/docs/mathlib-gap-report.md` — 缺口报告
- `proof/docs/phase0-experiment-plan.md` — Phase 0 实验设计
- `docs/design/factorization/formal-proof-ddf-edf.md` — DDF/EDF 半形式化证明
- `docs/design/factorization/formal-proof-univar-factorization.md` — Z[x] 半形式化证明
