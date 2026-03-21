# SQF derivative(c_rem) = 0：形式化受阻的教训

> 日期：2026-03-21 ~ 2026-03-22
> 分支：`feature/formal-proofs`

---

## 做了什么

围绕 `derivative(c_rem) = 0` 这一个定理，花费了大量时间但未完成。

**已完成（0 sorry）的新引理**：
- `yunLoop_preserves_pow_dvd`：q ∤ w, q^k | c → q^k | c_rem（yunLoop 归纳）
- `not_pow_dvd_derivative_of_separable`：q^v | f (exact), q 可分, p∤v → q^v ∤ f'
- `yunLoop_acc_subset`：yunLoop 只追加 acc
- `yunLoop_c_ne_zero`：yunLoop 保持 c ≠ 0
- `derivative_of_yun_remainder_eq_zero` 的 Steps (a)-(f) + Case 1

**未完成（sorry）**：
- Case 2 的核心 dvd 链（~15 行）
- `derivative(c_rem) = 0` 在 sqf_correct 中的调用
- Case C.1 的下游 3 sorry

## 关键发现：nl-proof §3.2.1 有数学错误

在尝试形式化 Case 2 的"无穷递升"论证时，发现 nl-proof 的步进不成立：

```
nl-proof 声称：q^{n+1} | crem → q^n | crem' → q^n | gcd → q^{n+2} | crem
实际上只能得到：q^{n+1} | crem → q^n | crem' → q^n | gcd
  → q^n * q | gcd * (crem/gcd) = crem → q^{n+1} | crem（回到起点，没有递增！）
```

**根因**：`q | (crem/gcd)` 只提供 1 次 q，`q^n | gcd` 提供 n 次，总共 n+1 次 = 起始值。无法递增到 n+2。

**正确路径**应该是 nl-proof §3.2.2 的 emultiplicity 方法（直接比较 v_q(c₀) 和 v_q(f')），不是无穷递升。

---

## 核心教训

### 教训 1（最重要）：严格按 nl-proof 翻译，遇到问题回到 nl-proof 修正

**问题**：我在 Lean 中遇到困难时，不是回到 nl-proof 检查数学正确性，而是：
- 尝试绕过（换路径、换策略）
- 写大量注释解释为什么某个方法应该行但实际不行
- 在 sorry 和代码之间反复横跳

**正确做法**：
1. nl-proof 必须是**可直接翻译的蓝图**，每一步都对应一个 Lean tactic/lemma
2. 如果 Lean 翻译卡住，**第一反应是检查 nl-proof 该步是否正确**
3. 发现 nl-proof 错误 → 停止编码 → 修正 nl-proof → 审核 → 再翻译
4. 绝不在 nl-proof 错误的基础上"创造性地绕过"

### 教训 2：nl-proof 的"一步到位"假设必须展开

nl-proof §3.2.1 写了"无穷递升，同 squarefree_div_gcd_derivative"。但实际上两个定理的递升结构不同：
- `squarefree_div_gcd_derivative`：`d^2 | w → d | f' → d | c → d^3 | f → ...`（每步 +1 通过 `f = w * c` 的两个因子）
- `derivative(c_rem) = 0` Case 2：`q^{n+1} | crem → q^n | crem'`（只有一个因子，无法通过 `crem = w' * gcd` 递增）

**规则**：nl-proof 中每次引用"同某个已有证明"时，必须验证**递归结构完全一致**。不能假设。

### 教训 3：sorry 越多越慢

从 0 sorry → 4 sorry → 5 sorry → 6 sorry 的过程中，每个新 sorry 都让后续推理更困难（不知道哪些假设是可靠的）。

**规则**：一个 sorry 解决后再开下一个。不要并行 sorry。

### 教训 4：set-bound 变量是独立引理的信号

每次遇到 `set crem := (yunLoop ...).2` 导致类型不匹配，都说明这块证明需要提取为独立引理。花时间在 inline proof 中调试 set-bound 是浪费。

## 后续计划

1. 修正 nl-proof §3.2.1 的 Case 2：回到 §3.2.2 的 emultiplicity 路径
2. 严格审核修正后的 nl-proof
3. 按修正后的 nl-proof 逐步翻译
4. 如果 emultiplicity API 仍然困难，考虑重新设计证明路径

## 度量

- 耗时：~8 小时（围绕 derivative(c_rem)=0 一个定理）
- 迭代：~50 轮编译-修复循环
- 净新增行数：~100 行已验证引理 + ~大量被删除的注释/sorry 代码
- 发现的 nl-proof 错误：1 个（§3.2.1 无穷递升步进错误）
- 放弃的方案：
  - 无穷递升法（步进不成立）
  - 纯 dvd 避免 emultiplicity（无法推出 q | w₀）
  - 独立引理 derivative_yunLoop_remainder_eq_zero（位置 + 参数匹配问题）
  - inline 证明在 Case C.1 内部（set-bound 类型不匹配）
