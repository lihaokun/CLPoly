# divmodAux 终止性证明完成（decreasing_by admit 清零）

日期：2026-07-26

## 做了什么

填掉 `CLPoly/Model.lean` 中 `divmodAux`（域上稀疏多项式长除法主循环）的 `decreasing_by admit`，
使其成为**真正的良基递归**（0 admit，`#print axioms` 无 `sorryAx`）。连带 `divmod` 恢复可计算、`#eval` 通过。

## 为什么做

`feature/l1-l2-refinement` 分支前 3 个提交搭建了三座不变量塔（Sorted / ReducedB / NonZeroB）及其
在各多项式运算下的保持引理，目的正是为长除法主循环的终止性提供「首项真抵消 ⇒ 度数严格下降」所需的语义前提。
本次把这些弹药组装成完整的终止性证明。

## 关键决策及理由

1. **度量选择 `if r.isEmpty then 0 else r[0]!.fst.deg+1`**：`+1` 编码使「r' 变空」与「r' 首项度 < r 首项度」
   两种情形都严格递减，覆盖 constant÷constant（dr=dg=0）边界。
2. **不变量进签名而非 decreasing_by 现推**：`NonZeroB r`（保证首项系数 ≠ 0、真抵消）与 `ReducedB r pm`
   （同素数、约化，保证模乘结论）必须归纳维持，故作为 `divmodAux` 参数，递归调用处用
   `ReducedB_subImpl` / `NonZeroB_subImpl` / `*_single_mul` 证明对 `r' = r - term*g` 保持。
3. **核心引理 `divmod_step_drop` 独立于 divmodAux**：把「r' 全元素度数 < dr」抽为独立定理，
   decreasing_by 只需一次 `generalize` + `split` + `omega` 收尾，避免 WF 上下文里 let/match 命名不可预测的陷阱。
4. **lambda 匹配靠顶层 `def scaleLambda`**：Model.lean 不 import Mathlib（无 `set` 战术），
   scaleByMonomial 内部 filterMap 的 lambda 提取为顶层 def，靠 defeq 在 `show` 处对接，避免手写 `have prod`/match 形式。

## 证明结构（nl-proof: docs/design/divmodAux-termination-nlproof.md）

- 引理 B `Zp_mul3_lc`：`(lc·gh).val=1 ⇒ ((a·lc)·gh).val = a.val`（三层 `Nat.mul_mod`+`mul_assoc`+`mod_mod`）。
- 引理 C `Zp_add_neg_cancel`：`c` 与 `a` 同 val 同 prime ⇒ `(a+(-c)).val = 0`（UInt64 减法不回绕 + 分情形）。
- 引理 D `mergeAdd_cancel_lead`：首项等度、系数和为 0 ⇒ mergeAdd 结果全 `< d`（走 else-else 抵消分支 + `mergeAdd_lt_all`）。
- 引理 A `scaleByMonomial_head_drop`：`term*g` 首项 = (dr, coeff·lc(g))，尾全 `< dr`（filterMap 头提取 + `sortedListB_filterMapShift`）。
- 组装 `divmod_step_drop`：r/g cons 分解 → 首项系数 coeff·lc(g) = r 首项系数（引理 B）→ 相减为 0（引理 C）
  → mergeAdd 抵消（引理 D）→ ALL-DROP。

## 遇到的问题与解决

- **`set` 不可用**（Mathlib 战术，Model 不 import）→ 改用顶层 `def scaleLambda` + 核心 `let`-free 写法。
- **`▸` 无法改写 opaque def**（`ReducedB`/`Sorted`/`NonZeroB`）→ 先 `have h' : reducedListB pm r.toList = true := hr_red` 展开。
- **`List.mem_cons_self` 不是函数**（本版本）→ 用 `by simp` 证首项 ∈。
- **decreasing_by goal 中 r' 被展开为完整表达式** → `generalize ... at hdrop ⊢` 抽象后 `split`。
- **olean 陈旧导致 #eval/#print axioms 误报 sorryAx** → `lake env lean` 只类型检查不更新 olean，须 `lake build CLPoly.Model` 刷新。

## 验证

- `lake build` 全 3464 jobs 通过，0 error。
- `#print axioms SparsePolyZp.divmod` = `[propext, Classical.choice, Quot.sound]`（**无 sorryAx**）。
- `#eval divmod (x²-1) (x-1)` = `(x+1, 0)` ✓；`#eval gcd (x²-1) (x-1)` = `x+4 = x-1` over F_5 ✓。

## 度量

- 耗时：~3 小时（调研三塔基础设施 + nl-proof + 形式化 + 调试）。
- 迭代：~12 轮编译-修复（多数在 scratch 文件快速验证单引理，再合入 Model 全量编译）。
- Lean 新增/修改行数：+255 / -23（Model.lean）。含 6 个辅助引理 + 4 个核心引理 + divmodAux 签名重构 + divmod 守卫扩展。
- 对应 C++ 行数：`upolynomial_<Zp>::divmod` 长除法主循环 ~30 行（被验证的算法：域上首一化长除法终止性）。
- 放弃的方案：
  - `set` 命名 lambda（Mathlib 战术不可用）；
  - 在 decreasing_by 内直接推 bound（WF 上下文命名不可预测，改为独立引理 `divmod_step_drop`）；
  - 显式手写 filterMap lambda（`have prod`/match 形式难精确复制，改用顶层 def）。

## 涉及文件

- `proof/lean/CLPoly/Model.lean`（主要变更）
- `proof/lean/docs/design/divmodAux-termination-nlproof.md`（nl-proof 草稿，新建）
- `proof/lean/docs/devlog/2026-07-26-divmodAux-termination-complete.md`（本文）

## 剩余债务（不属本次范围）

- `Model.lean:2065/2074`：`UInt64_toInt64_toNatClampNeg_le_toNat` / `_eq_toNat_of_lt` 两个 BitVec 转换引理仍 `admit`
  （SquarefreeZp 度数界用，独立子问题）。
