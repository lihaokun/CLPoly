# SQF 生成 raw 边界与导数零递归桥

日期：2026-08-10

## 做了什么

- 新增生成 SQF 入口所需的具体入口不变式和 Yun 局部不变式。
- 将生成层的 GCD 操作绑定到真实 public raw GCD，并保留对象输出缺失时的
  assertion failure。
- 将生成层 make-monic 操作绑定到实际 checked C++ 控制流。
- 证明生成 result-copy 循环与既有具体数组循环一致。
- 完成导数为零分支的 raw 执行桥：实际 derivative、p 次根和 make-monic
  返回值产生 canonical/monic/非空/正次数/机器界以及严格 measure 下降。

## 为什么做

生成的 `__squarefree_Zp_raw_ir` 只有在这些具体接口被 raw 实现实例化后，才是
可用于最终 L1→L2 定理的真实执行入口。

## 关键决策

- 新入口不变式要求正次数。这与 C++ SQF 的真实递归有效域一致，并排除旧包装
  对常数输入使用“递减失败则返回空结果”的非源代码行为。
- 所有输出相等性都通过 `RawExec.ok` 的单射从实际执行方程得到，不预测或选择
  GCD、p 次根或 make-monic 的结果。

## 验证

- `lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean`：通过。

## 涉及文件

- `proof/lean/CLPoly/Refinement/StrictSquarefreeGenerated.lean`
- `docs/devlog/2026-08-10-sqf-generated-raw-boundary.md`
