# Hensel 误差系数循环证书

## 日期

2026-08-11

## 做了什么

- 在严格生成语义中显式增加除法后的未约化数组 `divideCoeffs`。
- 证明 `divideThenReduceCoeffs_toPolyMod`：C++ 的除法、模约化和零项压缩
  输出，在 `ZMod m` 多项式下等于未约化误差商。
- 证明 `scaleCoeffs_divideCoeffs`：当每个输入系数可被 `m` 整除时，商的
  每项重新乘 `m` 后逐项恢复原始数组。
- 汇总为 `divideThenReduceCoeffs_certificate`，同时给出精确整数重构与模
  `m` 解码性质。

## 为什么做

`__hensel_step` 的两次核心误差修正都依赖 `(difference / m) mod m`。
如果只以 L2 存在性选一个误差商，就会跳过 C++ 实际的系数循环。这组引理
直接证明生成循环产生的具体数组就是所需误差商的模表示。

## 关键决策

- 保留 GMP floor division 对应的 `Int.fdiv`，没有替换为任意整除见证。
- 零项压缩只通过其在 `ZMod m` 下不改变多项式来证明，而不伪称压缩后的
  整数数组仍等于未约化商。
- 精确重构定理显式要求源代码注释所述的逐系数整除前提。

## 验证

- 严格 Hensel 生成文件同步。
- `CLPoly.Refinement.Hensel` 编译成功。
- 新增定理不含 `sorry`、oracle、fallback 或 fuel。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
