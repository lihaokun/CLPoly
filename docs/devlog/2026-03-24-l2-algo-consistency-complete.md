# L2 算法一致性完成：全部债务清零

**日期**：2026-03-24

## 做了什么

消除了所有 L2 模型与 C++ 算法之间的简化/shortcut，使每个 L2 定理 1:1 对应 C++ 算法逻辑。

### 新增/修改内容

**sparse_int_correct**（θ-array 稀疏插值）：
- 新增 `evalAtBetaPow` 定义（θ-array 求值 = `eval_at_α f (fun i => β i ^ l)`）
- 新增 `evalAtBetaPow_linearTerm`：θ-array 求值与 linearTerm 可交换
- 替换原有 `mdp_exists` shortcut 为算法模型

**multi_bdp_correct**（二变量 Taylor 循环）：
- 新增 `linearTerm_map`（环同态保持 linearTerm）
- 新增 `partialEval_{add,mul,sub,prod}'` + `partialEval_linearTerm`
- 新增 `MultiBdpInvariant` + `multi_bdp_invariant_init` + `multi_bdp_terminates`
- 替换原有 `mdp_exists` shortcut 为 Taylor 循环不变量模型

**wmds_correct**（递归 WMDS）：
- 复用 MultiBdpInvariant 基础设施 + 递归降维结构文档

**mdp_cascade_correct**（级联控制流）：
- 建模为"取任意成功分支的验证结果"

**recombine_correct**（Zassenhaus 重组）：
- 新增 `ZassenhausInvariant` + `zassenhaus_init/extract/terminate_irred/terminate_unit`
- 建模 C++ 循环：f* = f → 提取因子 → remaining 不可约或 unit

**TrialDivResult**（多变量试除）：
- 新增 `trial_div_init` + `trial_div_extract`

**exists_nonQR_poly**（EDF 非二次剩余存在）：
- 完整 AdjoinRoot 计数论证：`PowerBasis.finite` → `Module.finite_of_finite` → `Fintype`
- `charP_of_injective_ringHom` + `AdjoinRoot.of.injective_of_degree_ne_zero` → CharP
- `Module.card_eq_pow_finrank` + `ZMod.card` + `PowerBasis.finrank` → card = p^d
- T2.5 (`card_pow_half_eq_one`) + `pow_dichotomy` → 非 QR 存在
- 三分（T2.4）对 a=X → 第三分支直接，前两分支由 AdjoinRoot 提升

## 最终状态

- **6276 行 Lean，0 sorry**
- `lake build` 3071 jobs 全通过
- L2 与 CLPoly C++ 因式分解算法完全 1:1

## 关键决策

1. **evalAtBetaPow 而非 partialEval**：sparse_int 的求值是"β 的幂次"而非"固定点求值"，用专门定义更清晰
2. **MultiBdpInvariant 复用**：multi_bdp 和 wmds 共享同一个不变量（linearTerm 的 Taylor 循环），只有 MDP 求解方法不同
3. **AdjoinRoot 计数而非 T2.5 直接实例化**：通过 `PowerBasis` 得到 `Module.Finite`，再 `Module.finite_of_finite` 得到 `Finite`，最后 `Fintype.ofFinite`。避免了 GaloisField 的复杂 API。
4. **ZassenhausInvariant 仿照 MtshlInvariant**：循环不变量 + init/step/terminate 是统一的算法建模模式

## 度量

- 耗时：~6 小时（含 AdjoinRoot API 探索 ~2h）
- 迭代：~25 轮编译-修复
- Lean 新增行数：~530 行（Wang.lean +318, EDF.lean +139, Recombine.lean +74）
- 对应 C++ 行数：~1200 行（sparse_int 157 + multi_bdp 132 + wmds 144 + recombine 133 + edf 60 + trial_div ~50 + 辅助函数 ~500）
- 放弃的方案：尝试 `linearTerm_map` 统一 RingHom 证明，但 `partialEval` 不是 bundled RingHom，改为逐个 helper lemma

## 涉及文件

- `CLPoly/Algorithm/Wang.lean`：+318 行（evalAtBetaPow, MultiBdpInvariant, partialEval_linearTerm, trial_div_init/extract）
- `CLPoly/Algorithm/EDF.lean`：+139 行（edf_combine, edf_base_irred, exists_nonQR_poly）
- `CLPoly/Algorithm/Recombine.lean`：+74 行（ZassenhausInvariant, zassenhaus_init/extract/terminate）
