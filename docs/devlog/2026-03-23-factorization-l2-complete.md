# 2026-03-23 因式分解 L2 算法验证完成

## 做了什么

完成 CLPoly 单变量因式分解的 L2 算法模型形式化验证，全部 0 sorry。

本次会话新增/修改：
1. **hensel_step_with_degree H6**：Hensel 单步保持 leadingCoeff（~6 行改动）
2. **symmetric_recovery + factor_recovery**：对称模恢复引理，验证 C++ `__zassenhaus_recombine` 核心数学步骤（~40 行新增）
3. **hensel_unique 泛化**：从 `Monic B` 放宽为 `B1.leadingCoeff = B2.leadingCoeff ∧ ¬(p ∣ lc)`（~20 行改动）
4. **hensel_multifactor**：多因子 Hensel 提升，对 List 归纳，使用 hensel_two_factor 递归分裂（~55 行新增）

## 为什么做

L2 层目标是验证 C++ 因式分解算法的数学正确性。之前存在多个 L2-C++ 差距（Mignotte 界、Hensel 构造、lc-baking、度数保持、因子恢复），需要逐个修正。本次会话完成了所有剩余差距的修正。

## 关键决策

1. **hensel_unique 泛化而非特化**：原版要求 Monic（lc=1），但 lc-baking 场景下 B 可能有 lc = ℓ/ℓⱼ。泛化为"相同 lc + p 不整除 lc"，证明结构不变，适用范围更广。

2. **factor_recovery 以 hmod 为假设**：不在 factor_recovery 内部构造对比对（需要泛化 hensel_unique 的完整应用链），而是将 `map_m(A) = map_m(C(c)*g)` 作为外部假设。模块化更好，核心数学（对称恢复）独立可用。

3. **hensel_multifactor 使用 hensel_two_factor 递归**：多因子提升通过 head/tail 分裂 + 递归实现，避免建模 C++ 的 binary tree 结构。数学等价但形式化更简单。

## 最终状态

### 文件统计（全部 0 sorry）

| 文件 | 行数 | 内容 |
|------|------|------|
| Spec.lean | 133 | 子过程接口规约 |
| SquarefreeZp.lean | 1501 | SQF 算法模型 + 正确性 |
| DDF.lean | 435 | DDF 算法模型 + 正确性 |
| EDF.lean | 257 | EDF 算法模型 + 正确性 |
| Hensel.lean | 751 | Hensel 提升（单步、迭代、多因子）|
| Recombine.lean | 635 | Mignotte bound、Hensel 唯一性、因子恢复 |
| FactorZp.lean | 232 | Zp 因式分解 pipeline |
| FactorZpInstantiate.lean | 72 | Zp 端到端实例化 |
| FactorZZ.lean | 71 | ZZ 因式分解 pipeline |
| FactorZZInstantiate.lean | 23 | ZZ 端到端实例化 |
| **总计** | **4110** | |

### 已验证的定理

| 定理 | 对应 C++ | 说明 |
|------|---------|------|
| sqfZp_correct | __polynomial_sqf_yun | Yun SQF 分解 |
| ddf_correct | __ddf_shoup | Shoup DDF |
| edf_correct | __edf_cantor_zassenhaus | Cantor-Zassenhaus EDF |
| hensel_step / hensel_step_with_degree | __hensel_step | 2-factor Hensel 单步 (H1-H6) |
| hensel_two_factor | hensel_lift 循环 | 2-factor mod p → mod p^k |
| hensel_multifactor | hensel_lift 多因子 | r-factor mod p → mod p^k |
| mignotte_bound_l2 | __mignotte_bound | M(f) ≤ ‖f‖₂ (Parseval + Jensen) |
| hensel_unique | Hensel 唯一性 (implicit) | 泛化版：相同 lc + p ∤ lc |
| symmetric_recovery | __symmetric_mod | 小系数 + 同余 → 精确相等 |
| factor_recovery | __zassenhaus_recombine 核心 | Hensel 因子恢复真因子 |
| recombine_correct | __zassenhaus_recombine | Z[x] 不可约分解 |
| factor_Zp_correct | __polynomial_factorize_Zp | Zp 完整因式分解 |
| factor_ZZ_correct | polynomial_factorize_univar | ZZ 完整因式分解 |

### L2-C++ 差距表（全部关闭）

| 差距 | 状态 |
|------|------|
| Mignotte 界 | ✅ Landau 不等式 |
| Hensel 构造 | ✅ canonical_lift + modByMonic |
| Hensel 度数保持 | ✅ H4 + H5 |
| Hensel lc 保持 | ✅ H6 + 多因子归纳 |
| Hensel 唯一性 | ✅ 泛化版 |
| 因子恢复 | ✅ symmetric_recovery |
| 多因子 Hensel | ✅ hensel_multifactor |
| Recombination | ✅ UFD + factor_recovery |

### 待做（非因式分解）

- GCD L2 验证（Euclidean / HGCD）
- 多项式算术 L2 验证（pair_vec_div / multiplies）
- L1 实现模型（1:1 对应 C++ 控制流）

## 度量

- 耗时：~3 小时（含 nl-proof 编写、审核、形式化、调试）
- 迭代：~12 轮编译-修复循环
- Lean 新增/修改行数：~120 行
- 对应 C++ 行数：~600 行（polynomial_factorize_univar.hh 全部因式分解流程）
- 放弃的方案：无

## 涉及文件

- `proof/lean/CLPoly/Algorithm/Hensel.lean` — H6 + hensel_multifactor
- `proof/lean/CLPoly/Algorithm/Recombine.lean` — symmetric_recovery + factor_recovery + hensel_unique 泛化
- `proof/lean/CLPoly/Spec.lean` — 未改动
- `proof/CLAUDE.md` — 差距表更新
- `proof/nl-proof/hensel-lc-baking.md` — lc 保持 nl-proof
- `proof/nl-proof/factor-recovery.md` — 因子恢复 nl-proof
- `proof/nl-proof/hensel-unique-generalize.md` — 唯一性泛化 nl-proof
