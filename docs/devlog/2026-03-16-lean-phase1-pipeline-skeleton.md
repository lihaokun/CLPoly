# Lean 4 形式化验证：Phase 1 完成 — 顶层骨架 + 接口锁定

> 日期：2026-03-16
> 分支：`feature/formal-proofs`
> 前置：T1.1 Spec.lean 完成（`b0a3797`）

---

## 1. 完成的任务

### T1.1 补充：Spec.lean 编译修正

- `Polynomial.natDeg` → `Polynomial.natDegree`：Mathlib 4 中完整名称是 `natDegree`

### T1.2 Zp[x] 顶层骨架

**文件**：`CLPoly/Pipeline/FactorZp.lean`

```lean
theorem factor_Zp_correct
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    (sqf : ...) (hsqf : ∀ g, g ≠ 0 → SquarefreeDecomp g (sqf g))
    (ddf : ...) (hddf : ∀ g, Monic g → Squarefree g → DDFCorrect g (ddf g))
    (edf : ...) (hedf : ∀ g d, Monic g → Squarefree g → 等度条件 → EDFCorrect g d (edf g d))
    : ∃ (lc : ZMod p) (factors : ...), FactorZpCorrect f lc factors
```

结构：假设 SQF/DDF/EDF 各自正确 → 组合结果满足 `FactorZpCorrect`。证明体 `sorry`。

### T1.3 Z[x] 顶层骨架

**文件**：`CLPoly/Pipeline/FactorZZ.lean`

```lean
theorem factor_ZZ_correct
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    {p : ℕ} [Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    (hgood : Squarefree (map f mod p)) (hdeg : 度数保持)
    (factor_zp : ...) (hfzp : ∀ g, g ≠ 0 → FactorZpCorrect ...)
    (hensel : ...) (hhensel : 乘积条件 → HenselCorrect ...)
    (recombine : ...) (hrecombine : 乘积条件 → RecombineCorrect ...)
    : ∃ result, FactorZZCorrect f result
```

设计要点：
- `p` 和 `k` 作为隐式参数（算法选择的素数和提升指数）
- 好素数条件：`Squarefree (f mod p)` + 度数保持
- Hensel/Recombine 假设带前置条件（`map f = facs.prod`），非无条件全称量化

### T1.4 接口一致性审查

逐条检查 Spec 规约之间的衔接：

| 衔接 | 结果 |
|------|------|
| SQF → DDF（Monic + Squarefree） | ✓ 匹配 |
| DDF → EDF（度数条件） | ✓ 匹配 |
| DDF → EDF（Monic） | **缺口**：DDFCorrect 未声明 |
| DDF → EDF（Squarefree） | ✓ 可推导（UFD 中无平方的因子也无平方） |
| FactorZp → Hensel | ✓ 可连接 |
| Hensel → Recombine | ✓ 匹配 |

**修正**：DDFCorrect 添加第 5 条 `∀ pr ∈ result, Monic pr.1`。

## 2. Phase 1 验收

| 验收标准 | 状态 |
|---------|------|
| Spec.lean 所有规约 Prop 编译通过 | ✓ |
| FactorZp.lean 顶层定理陈述编译通过 | ✓ |
| FactorZZ.lean 顶层定理陈述编译通过 | ✓ |
| `lake build` 全量通过（1908 jobs） | ✓ |
| 接口一致性审查 + 修正 | ✓ |

## 3. 未完成：框架组合证明

Phase 1 的两个 `sorry` 是**框架组合证明**——假设子过程满足规约，证明组合达成顶层目标。这不依赖算法实现，纯粹是规约间的逻辑衔接。

应在进入 Phase 2（数学基石）前完成，以验证接口设计的完备性。已知需要的辅助引理：
- 无平方多项式的因子也无平方（UFD 性质）
- `Associated` 的传递性（三层乘积还原组合）

## 度量
- 耗时：~6 小时（T1.2 ~2h、T1.3 ~1.5h、T1.4 接口审查 ~1.5h、Spec 修正 ~1h）
- 迭代：~5 轮编译-修复循环（DDFCorrect Monic 缺口发现后需修改 Spec + 重新对齐 Pipeline）
- Lean 新增行数：89 行（FactorZp 37 + FactorZZ 48 + Spec 修正 4）
- 对应 C++ 行数：~200 行（polynomial_factorize_zp.hh 的管线主函数 + polynomial_factorize_univar.hh 的 factor_univar）
- 放弃的方案：最初 FactorZZ 将 `p, k` 作为显式参数，后改为隐式（算法内部选择，顶层无需暴露）

## 4. 涉及文件

| 文件 | 操作 |
|------|------|
| `proof/lean/CLPoly/Spec.lean` | 修改：natDeg→natDegree，DDFCorrect 添加 Monic |
| `proof/lean/CLPoly/Pipeline/FactorZp.lean` | 新建 |
| `proof/lean/CLPoly/Pipeline/FactorZZ.lean` | 新建 |
| `proof/lean/CLPoly.lean` | 修改：添加 Pipeline 导入 |
