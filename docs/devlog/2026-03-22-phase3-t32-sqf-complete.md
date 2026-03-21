# Phase 3 T3.2 SQF 完成：4 sorry → 0 sorry

> 日期：2026-03-22
> 分支：`feature/formal-proofs`

---

## 做了什么

完成 SQF（Squarefree Decomposition）算法模型的 Lean 4 形式化证明，从 4 sorry 降至 **0 sorry**。

**修复的核心问题**：`derivative_of_yun_remainder_eq_zero` Case 2 的"无穷递升"证明有数学错误（步进量 = 0），替换为正确的 v-based 纯 dvd 证明。

**具体改动**：

1. **新增 `yunLoop_crem_dvd_c` 引理**（~50 行）
   - 证明 `(yunLoop w c i acc hc).2 ∣ c`（yunLoop 残余整除初始 c）
   - yunLoop 归纳：每步 c' = normalize(c /ₘ gcd)，c' | c by dvd_trans

2. **重写 `derivative_of_yun_remainder_eq_zero` Case 2**（~80 行替换 ~20 行）
   - 旧：无穷递升 `q^{n+1} | crem → q^{n+2} | crem`（步进 = 0，类型错误）
   - 新：nl-proof §3.2.1 steps 2a-2i 的精确翻译
   - 关键链：找 max v → q^v | c₀（backward via crem | c₀）→ q^{v+1} ∤ c₀（forward contrapositive）→ p∤v, q'≠0（derivative formula）→ q^v ∤ f'（not_pow_dvd_derivative_of_separable）→ 矛盾 c₀ | f'

3. **填补 `sqf_correct` 3 个下游 sorry**
   - Associated f (yun_prod * pth_prod)：通过 `hcrem_assoc_gp.trans hpth.1` 链接
   - 跨组互素（yun vs pth-root）：pr₂.1 | pth_prod ~ crem + IsCoprime pr₁.1 crem → IsCoprime

4. **清理 nl-proof**：删除废弃的 §3.2.2 emultiplicity 方案和 §3.2.1 后半段探索性笔记

## 为什么做

SQF 4 sorry 全部阻塞在 `derivative(c_rem) = 0` 这一个定理。前一轮（~8 小时）发现了 nl-proof §3.2.1 的数学错误但未修正代码。本轮按修正后的 nl-proof 严格翻译，一次通过。

## 关键决策

1. **纯 dvd 方案 vs emultiplicity 方案**：选择纯 dvd。避免引入 ℕ∞ 算术和 `emultiplicity_mul` 等 API。用 `Nat.find` + `Classical.decPred` 获取 max power v。

2. **新增 `hcrem_dvd_c₀` 参数**：`derivative_of_yun_remainder_eq_zero` 需要 backward 方向（q^v | crem → q^v | c₀），通过新引理 `yunLoop_crem_dvd_c` 提供。

3. **`Associated.of_eq` 处理 set-bound 变量**：`crem` 是 `set`-bound，`rw` 找不到 pattern。用 `Associated.of_eq hcrem_eq_pow` 建立 `Associated crem (contract p crem ^ p)`，绕开 `rw` 限制。

## 遇到的问题

| 问题 | 解法 |
|------|------|
| `Nat.find` 需要 `DecidablePred` | `@Nat.find ... (Classical.decPred _) hfin` |
| `List.mem_cons_self` 不是函数（Lean 4 API 变化） | 改用 `by simp` |
| `List.mem_map_of_mem` 同上 | 改用 `List.mem_map.mpr ⟨pr, h, rfl⟩` |
| `derivative_pow` 后 `mul_zero` 模式匹配失败 | 改用 `ring` 替代逐步 `rw` |
| `IsCoprime.pow_left` 参数名 `n` 不存在 | 改为 `m` |
| set-bound `crem` 阻碍 `rw [hcrem_eq_pow]` | 提前用 `Associated.of_eq` 建立桥接 |

## 涉及的文件

- `proof/lean/CLPoly/Algorithm/SquarefreeZp.lean` — 核心改动（+240 行/-20 行）
- `proof/lean/CLPoly.lean` — 注释掉 E4 import（与 Algorithm.DDF 冲突）
- `proof/nl-proof/phase3-t32-sqf.md` — 清理废弃方案（-190 行）

## 度量

- 耗时：~2 小时（含 6 轮编译-修复循环）
- 迭代：6 轮编译-修复
- Lean 新增/修改行数：~240 行净增
- 放弃的方案：无（严格按修正后的 nl-proof 翻译，一次成功）

## 最终状态

| 文件 | 行数 | sorry |
|------|------|-------|
| DDF.lean | 435 | 0 |
| SquarefreeZp.lean | 1501 | 0 |
| **合计** | **1936** | **0** |

`lake build` 全量通过（1910 jobs）。
