# 生成严格 SQF raw 控制流

日期：2026-08-10

## 做了什么

- 新增 `build_strict_squarefree.py`，生成 `Generated.StrictSquarefreeZp`。
- 生成 C++ `__squarefree_Zp` 的完整 raw 控制流，包括：
  - 导数为零时的 p 次根、make-monic、递归和重数缩放；
  - 非零导数时的 GCD、精确除法、Yun 循环；
  - Yun 循环后的剩余因子递归和结果拼接。
- 将 SQF 顶层递归和 Yun 循环分别包装为携带不变式的依赖状态，并使用源
  measure 做良基递归。

## 为什么做

此前最终 SQF 定理的执行侧是 refinement 层手写的 `strictSquarefreeZpIR`，
没有直接落在 cpp2lean 生成层。新生成模块为后续 raw 操作实例化和最终
`__squarefree_Zp` 精化提供真正的生成 L1 入口。

## 关键决策

- 生成定义中不引用 `sqfZp` 或其他 L2 算法。
- GCD、精确除法、p 次根和 make-monic 都保留为 `RawExec` 操作；不变式接口
  只能证明实际返回值的安全性和 measure 下降，不能选择返回值。
- 删除旧严格包装中的“measure 不下降则返回空结果”行为。生成入口在有效
  不变式下直接证明递减并执行源递归，因此不会改变 C++ 的有效输入语义。
- 所有递归均为 `termination_by` 良基递归，不引入 fuel。

## 验证

- `lake build CLPoly.Generated.StrictSquarefreeZp`：通过。
- 生成器拒绝 `sorry`、`admit`、`axiom`、`partial def` 和 `sqfZp`。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_squarefree.py`
- `proof/lean/CLPoly/Generated/StrictSquarefreeZp.lean`
- `docs/devlog/2026-08-10-generated-strict-sqf-raw-control-flow.md`
