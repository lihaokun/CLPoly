# Hensel modular divmod trace conservation

## 做了什么

- 证明非空稀疏整数多项式的 `get_deg.toNatClampNeg` 精确等于首项自然次数。
- 证明 active 分支中的 signed `get_deg r - get_deg g` 精确等于自然次数差，且无机器字回绕。
- 证明 quotient coefficient 为零时，具体 `eraseIdxIfInBounds 0` 不改变模多项式语义。
- 定义 trace 的机器次数安全谓词，并对 `DivmodTrace` 完成商余守恒归纳。

## 为什么做

生成的 `__upoly_divmod_mod_raw_ir` 通过精确有限 `DivmodTrace` 执行 C++ long division。要证明真正的精化，必须证明两种递归分支都保持 `remainder + quotient * divisor`，同时解释源程序的 `Int64` 次数运算。

## 关键决策

- `DivmodTraceDegreeSafe` 只记录机器次数范围事实，不包含或预测最终商余式。
- `vanished` 分支使用实际 `ZZ.invert` 和实际 coefficient 计算证明被删除系数模 `m` 为零。
- `subtract` 分支直接复用已证明的 `divmodRemainder_eq_sub`，并跟踪具体 quotient array push。
- 下一步仍需从 canonical/sorted 输入不变式推出整个 trace 的 degree-safe 证书；当前不会把该证书冒充最终定理。

## 遇到的问题与解决方式

- `Int64` subtraction 是带回绕的机器运算；通过 `Int64.toInt_sub` 和 signed 区间界证明 active 分支不回绕。
- 零 coefficient 并不意味着整数首项为零；先用具体 inverse 得到其 `ZMod` cast 为零，再证明删除保持语义。

## 涉及文件

- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `docs/devlog/2026-08-11-hensel-divmod-trace-conservation.md`
