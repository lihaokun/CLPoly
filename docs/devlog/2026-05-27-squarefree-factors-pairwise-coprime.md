# squarefree_factors_pairwise_coprime 补全 + Hensel 提升类型设计修正

## 日期
2026-05-27

## 做了什么

### 1. 补全 squarefree_factors_pairwise_coprime 引理
在 `FactorZZ.lean` 中完成了 `squarefree_factors_pairwise_coprime` 引理的全套证明，
该引理证明 Squarefree Zp 因式分解的各因子两两互素（IsCoprime）。

新增辅助引理：
- `monic_irreducible_coprime`：Z/pZ 上的首一不可约多项式若不等则互素
- `IsCoprime_of_isUnit_left / IsCoprime_of_isUnit_right`：单位与任意元素互素
- `list_dvd_prod_of_mem / list_prod_mul_of_two_pos`：列表中两不同位置元素整除列表乘积

同时将 `factor_ZZ_correct` 的 `hhensel` 参数签名从
`(facs_p.Pairwise (fun a b => IsCoprime a b)) → HenselCorrect f k facs_p (hensel facs_p)`
改为显式传递 `hcop` 参数。

### 2. Hensel 提升类型设计修正
`FactorZZInstantiate.lean` 中解决了 `hensel_lift` 与 `hensel_lift_correct`
之间的 `Classical.choice` 非确定性核心矛盾。

#### 关键决策
将提升逻辑分拆为两层：
- `hensel_lift_with_proof`：带证明参数 `hne hprod hcop` 的构造器，内含 `hensel_multifactor` 调用
- `hensel_lift`：无证明参数的包装函数（`by classical; if ... then ... else ...`），
  与 `hensel_lift_correct` 通过调用 `hensel_lift_with_proof_correct` 完成证明

这种分拆确保 `hensel_multifactor` 在构造与证明中使用**同一组证明参数**，
从而 `Classical.choice` 返回一致的结果。

#### 其他修正
- 将 `by_cases` 改为 `by classical; if ... else ...` 以支持 `Decidable` 合成
- 多个 `calc` block 语法问题（Lean 4 `calc` 的限制，改用 `refine` + `calc` 或 `simp`）
- `canonical_lift_map` 参数命名冲突（使用 `(m := p) (_hm := ...) (p := x)`）
- `ext` 与 `funext` 的选择（多项式 `ext` 按系数展开，函数需要 `funext`）
- 多个 `List` 引理查找（`List.length_eq_zero` 在 Batteries 中不存在，改用 `simpa`）

### 3. 构建状态
`lake build CLPoly` 通过。仅剩 4 个 refinement `sorry`（非阻塞项）：
`__binomial_ir_refines`, `__isqrt_ceil_ir_refines`, `__squarefree_Zp_ir_refines`,
`__ddf_Zp_ir_refines`。

## 涉及文件
- `proof/lean/CLPoly/Pipeline/FactorZZ.lean`
- `proof/lean/CLPoly/Pipeline/FactorZZInstantiate.lean`
