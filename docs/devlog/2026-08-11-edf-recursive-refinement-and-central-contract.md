# EDF 递归精化与集中公开契约

## 日期

2026-08-11

## 做了什么

- 证明 cpp2lean 生成的 `__edf_Zp_raw_ir_state` 完整递归执行满足
  `EDFCorrect`，同时精确保留输入累加器前缀和 RNG 状态转移。
- 证明候选因子、精确除法及两次首一化调用确实成功，并证明两个递归
  子问题的乘积重新组成原多项式。
- 增加入口证明 `strictEDFEntryIR_refines_edf`。
- 由 Pass 9 在 `CLPoly/Refinement/Generated.lean` 集中生成公开定理
  `__edf_Zp_raw_ir_refines_edf`。
- 停止生成分散且带 `sorry` 的 `GeneratedSkeletons/*` 文件。

## 为什么做

最终精化定理此前散落在算法证明文件中，而且旧式 skeleton 与已证明契约
并存，难以判断哪个才是可信入口。现在所有已经闭合的公开 L1→L2 契约均
由同一生成器集中放在 `Refinement/Generated.lean`；具体证明仍按算法保留在
`Refinement/EDF.lean`，便于维护而不污染公开索引。

## 关键决策

- EDF 的 L2 契约使用 `EDFCorrect`，因为随机算法的具体因子顺序依赖 RNG，
  不应伪造为与某个确定性因子列表逐项相等。
- 递归严格使用 `edfMeasure` 的良基递归；重试终止条件携带真实有限执行
  trace，没有 fuel、超时回退或 L2 见证替换。
- 公开定理保留原 C++ 函数名中的双下划线，以便机械追踪。

## 遇到的问题与解决方式

- `Int64` 的顺序不能直接套用通用 `Preorder` 反证；通过
  `Int64.le_iff_toInt_le` 和 `Int64.lt_iff_toInt_lt` 转到整数后由 `omega`
  消除矛盾。
- 嵌套 `certifyRawExec` 匹配需要在每次重写后进行 iota 化简，才能继续
  展开下一次真实调用。

## 验证

- `lake build CLPoly.Refinement.EDF CLPoly.Refinement.Generated`：成功。
- Pass 9 重跑只生成 `Refinement/Generated.lean`。
- EDF 严格链路扫描无 `sorry`、`admit`、fuel、fallback 或 oracle。
- 三层定理的公理审计均仅包含 `propext`、`Classical.choice`、`Quot.sound`。

## 涉及文件

- `proof/lean/CLPoly/Refinement/EDF.lean`
- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass9_refine_gen.py`
