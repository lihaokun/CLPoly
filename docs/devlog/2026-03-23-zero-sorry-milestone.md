# 里程碑：全项目 0 sorry

> 日期：2026-03-23
> 分支：`feature/formal-proofs`

---

## 做了什么

完成全部剩余 sorry，达成 **3257 行 Lean，0 sorry，lake build 3026 jobs 全通过**。

今日完成：
1. **Mignotte bound**（mignotte_bound，0 sorry）
   - Mathlib Mahler measure API：`norm_coeff_le_choose_mul_mahlerMeasure` + `mahlerMeasure_mul` + `mahlerMeasure_le_sum_norm_coeff`
   - `one_le_mahlerMeasure_of_ne_zero`（NumberTheory.MahlerMeasure）
   - L1 版本：‖g‖_∞ ≤ C(n,n/2) · ‖f‖₁

2. **Hensel 唯一性**（hensel_unique，~130 行，0 sorry）
   - 对 k 归纳，每步：
     - `exists_C_mul_of_map_eq_zero` 提取误差 E, G
     - 乘积展开 → `p | (E·B₂ + A₂·G)` 的每个系数
     - `IsCoprime` → `Ā | Ē`, `B̄ | Ḡ`
     - B monic + 度数约束 → `G_bar = 0`, `E_bar = 0`
     - lift: `map_p(Q) = 0 → p^{k+1} | C(p^k)·Q`

3. **recombine_correct**（0 sorry）
   - Z[x] UFD 构造不可约分解

## 完整验证链

```
F_p[x] 因式分解：
  SQF (1501 行) → DDF (435 行) → EDF (257 行)
  → Pipeline 组合 (232 行) → 端到端实例化 (72 行)

Z[x] 因式分解：
  Hensel 提升 (366 行) → Mignotte bound (80 行)
  → Hensel 唯一性 (130 行) → Recombination (323 行)
  → Pipeline ZZ (71 行)
```

## 最终统计

| 模块 | 行数 | sorry |
|------|------|-------|
| SquarefreeZp.lean | 1501 | 0 |
| DDF.lean | 435 | 0 |
| EDF.lean | 257 | 0 |
| Hensel.lean | 366 | 0 |
| Recombine.lean | 323 | 0 |
| FactorZp.lean | 232 | 0 |
| FactorZZ.lean | 71 | 0 |
| FactorZpInstantiate.lean | 72 | 0 |
| **合计** | **3257** | **0** |

## 关键教训（本次 Phase 4）

1. **先对齐 C++，不要简化**：UFD 存在性和弱界都不能验证 C++ 算法
2. **先查 Mathlib 再声明 sorry**：Mahler measure + Jensen 公式已完整形式化
3. **精确条件很重要**：Hensel 唯一性需要 B 首一 in Z_m[x]（不仅 mod p）
4. **nl-proof 审核 5 项标准**：数学正确 + 无跳步 + Lean 可形式化 + 工程问题 + 边界情况
5. **调研要彻底**：每个 Mathlib API 找到精确名称和签名后再写代码

## 度量

- 耗时：~20 小时（Phase 3-4 总计，含 nl-proof 审核 + 形式化 + 调试）
- Lean 总行数：3257 行
- nl-proof 总轮次：SQF 4 版 + EDF 2 版 + Hensel 4 版 + Recombine 5 版 + Mignotte 2 版
- Mathlib 依赖：ZMod, Polynomial, MahlerMeasure, JensenFormula, UFD, IsCoprime
