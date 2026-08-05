# SQF Branch B：除法与 GCD 语义桥

日期：2026-08-05

## 做了什么

- 为实际的依赖式 `SparsePolyZp.divmodAux` 证明 `f = q * g + r` 的多项式语义恒等式。
- 证明除法余式保持排序、非零系数和模素数约化不变量，并在非空时严格降低首项次数。
- 为实际的良基 `SparsePolyZp.gcdAux` 证明表示不变量保持及其结果与 `EuclideanDomain.gcd` 相伴。
- 将证明延伸到 `polynomial_GCD`：证明 `makeMonic` 乘入的逆首项系数是单位，因此最终结果仍与数学 GCD 相伴。
- 证明规范稀疏表示的数组首项等于 `toPoly` 的数学首项，进而把非零 `polynomial_GCD` 加强为与 `normalize gcd` 严格相等。
- 证明整除情形下实际 `divmod` 的商严格等于 mathlib 的 `divByMonic` 商。
- 证明 `SparsePolyZp.normalization`（过滤零系数）不改变 `toPoly`。

## 为什么做

SQF 的导数非零分支调用真实的 `polynomial_GCD` 和 `pair_vec_div`，不能用影子算法或抽象规格替代。Yun 循环精化需要先获得真实除法恒等式、余式不变量和 GCD 相伴关系，才能对每一轮的 `y`、`z`、`c` 建立语义对应。

## 关键决策

- 使用统一的 `CanonicalRep` 汇总排序、非零系数和 `AllReduced`，作为除法与 Euclid 递归的不变量。
- GCD 采用 `Associated` 而非强行证明具体代表相等；`makeMonic` 的语义通过非零常数多项式是单位来闭合。
- 所有定理直接针对生产定义 `divmodAux`、`divmod`、`gcdAux` 和 `polynomial_GCD`。

## 验证

- `lake env lean CLPoly/Refinement/SquarefreeZp.lean`：通过；仅保留原有 Branch B 的一个 `admit`。
- `git diff --check`：通过。

## 度量

- Lean 新增：约 560 行。
- 编译修复：多轮单文件验证。
- 已闭合层次：原始有限域逆元、除法恒等式、余式规范形、Euclid GCD、首一化 GCD、精确首一商。
- 下一层：C++ partial Yun 循环与 L2 `yunLoop` 的逐轮对应。
