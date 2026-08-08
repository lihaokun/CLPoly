# HGCD 首次重构严格下降

## 结论

首次递归 HGCD 返回后，真实的成对重构结果满足
`reconstructed.lenB < lenA`，因此可直接作为后续良基递归边的下降证据，
不需要 fuel、运行时下降保护或 L2 回退。

## 证明来源

- C++ 调用入口已有的严格次序 `lenB < lenA`；
- 首次 HGCD 返回值的停止界与 `lenA > 0`；
- 四个真实矩阵描述符中的第 0、2 项长度界；
- `_mul` 与 `_poly_add` 原始堆执行导出的精确乘积 `- 1` 长度界；
- C++ 使用的两个低半部长度 `min lenA m`、`min lenB m`。

`HgcdRecursiveLengthInvariant` 新增的正长度字段在真实 base/iterator
执行分支中分别由输入正长度和欧几里得循环不变量导出。首次重构定理显式消费
源入口的 `lenB < lenA`，没有把它藏进 workspace、规格 oracle 或额外公理。

## 验证

- `lake build CLPoly.Impl.StrictHGCDRawRefinement`
- `python3 proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `#print axioms`：严格下降算术定理仅依赖 `propext`、`Quot.sound`；
  原始堆精化定理仅再依赖 `Classical.choice`，无 `sorryAx`。

## 后续 leading-A 接口

原始表示长度桥同时推广为：只要整个低部的规范化长度严格小于
`shift + highLength`，移位高部的最高项就不可能被低部抵消，输出长度精确等于
`shift + highLength`。旧的 `lowLength ≤ shift` 版本保留为推论。这个推广
直接对应最终重构所需的次数分离，不引入“结果非零”之类的规格前提。

## 最终 A 首项保持

对照 C++ 第二次重构公式，新增并从真实迭代执行证明统一矩阵系数界：
四个 `_mat_one`/`_mat_row_update` 描述符长度均不超过
`inputLength - inputLength / 2`。每一步证明只使用源循环守卫、真实 divrem
商长度和两次真实行更新的长度结果。

递归长度契约现在同时携带 `aboveHalf` 和该系数界。于是最终重构的
`S[3] * lowA`、`S[1] * lowB` 均严格低于 `X^k * secondA` 的最高项，
由规范化 raw 表示推出 `result.lenA = k + second.lenA`，特别得到
`0 < result.lenA`。这不是对输出非零的假设，而是源矩阵更新与物理重构的推论。
代入源码的 `k = 2*m-lenB2+1` 后还得到
`outerLength / 2 < result.lenA`，因此该不变量可继续传给递归父调用。
