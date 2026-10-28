# Phase 1 完成：框架组合证明 0 sorry

> 日期：2026-03-17
> 分支：`feature/formal-proofs`

---

## 做了什么

将 `FactorZp.lean` 和 `FactorZZ.lean` 中所有 sorry 替换为完整证明，Phase 1 框架组合验证达到 **0 sorry**。

### FactorZp.lean — 5 个辅助引理 + key 证明

| 引理 | 证明方法 |
|------|---------|
| `poly_unit_eq_C` | `coeff_coe_units_zero_ne_zero` + `eq_C_of_natDegree_eq_zero` + `natDegree_coe_units` |
| `eq_of_associated_monic` | 单位是常数 → 比较 `leadingCoeff` → `c = 1` → `mul_one` |
| `monic_list_prod` | 列表归纳 + `Monic.mul` |
| `list_prod_pow` | 列表归纳 + `mul_pow` |
| `list_prod_flatMap` | 列表归纳 + `List.prod_append` |
| `key`（主定理 Step 4） | `hfactors_def` 展开 let → `simp_rw [List.map_flatMap, List.map_map]` → `list_prod_flatMap` → `List.map_congr_left` 逐项 + `conv_rhs` 定向重写 + `list_prod_pow` |

### FactorZZ.lean — 已在上一轮完成（0 sorry）

## 关键技术点

1. **`let` 绑定在 tactic 中不可 `rw` 展开**：解法是 `have hfactors_def : factors = ... := rfl` 提供一个可 `rw` 的等式。

2. **`rw` 替换所有匹配项**：`rw [(component_ok ge hge).1]` 会同时替换 LHS 和 RHS 中的 `ge.1`。解法是 `conv_rhs => rw [...]` 限定重写范围。

3. **首一 Associated → 相等**：这是 DDF→EDF 组合的核心。`Associated a b` 给出 `a * u = b`；单位 `u` 是常数 `C c`；比较 `leadingCoeff` 得 `c = 1`；`mul_one` 得 `a = b`。

## Phase 1 验收状态

| 验收标准 | 状态 |
|---------|------|
| Spec.lean 全部 Prop 定义编译通过 | ✓ |
| FactorZp.lean 主定理 0 sorry | ✓ |
| FactorZZ.lean 主定理 0 sorry | ✓ |
| 子过程接口规约已锁定 | ✓ |
| **框架可组合性验证**（新增标准） | ✓ |

## 度量
- 耗时：~5 小时（辅助引理 ~2h、key 证明 ~2h、FactorZZ sorry 消除 ~1h）
- 迭代：~6 轮编译-修复循环（`conv_rhs` 定向重写技巧发现前尝试了 3 种失败方案）
- Lean 新增/修改行数：~45 行净增（FactorZp sorry 替换为证明，FactorZZ 同）
- 对应 C++ 行数：N/A（组合证明无直接对应）
- 放弃的方案：(1) 直接 `simp` 展开 flatMap 乘积——失败，simp 无法处理嵌套 map；(2) `rw` 不加 `conv` 限定——替换了不想替换的位置

## 涉及文件

| 文件 | 操作 |
|------|------|
| `proof/lean/CLPoly/Pipeline/FactorZp.lean` | 填充全部 sorry → 0 sorry |
| `proof/lean/CLPoly/Pipeline/FactorZZ.lean` | 无变更（已在上一轮完成） |
