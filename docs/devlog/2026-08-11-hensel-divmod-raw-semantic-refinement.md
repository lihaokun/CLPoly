# Hensel modular divmod raw semantic refinement

## 做了什么

- 定义稀疏整数数组的全项次数界 `DegreesBound` 和首项最大性 `HeadDominates`。
- 对 `pushNonzero`、六分支 `divmodMergeLoop`、`divmodRemainder`、首项删除和 coefficient reduction 证明次数界保持。
- 从普通输入表示条件自动构造整条 `DivmodTraceDegreeSafe`。
- 建立最终 raw-entry 定理 `__upoly_divmod_mod_raw_ir_refines`。

## 为什么做

此前商余守恒定理仍显式接收 trace 的机器次数安全证书。真正的 C++ L1 到 Lean L2 精化必须从输入表示不变式推出该证书，不能把中间执行安全性留给调用方假设。

## 关键决策

- `DegreesBound` 检查数组全部项，而不是只检查当前首项，因此 erase 和下一轮 remainder 可直接继承。
- `HeadDominates g` 足以证明所有移位 divisor 次数不超过当前 remainder 首项；无需借用 L2 多项式除法。
- `divmodMergeLoop_degreesBound` 再次按生成函数自身的良基归纳覆盖六个具体分支。
- raw-entry 最终等式直接针对生成 trace 返回的具体数组：`r + q * g = f`（模 `m`）。

## 遇到的问题与解决方式

- 负乘积 coefficient 分支的双重整数 remainder 会被 Lean 规范化；同步规范化递归归纳假设后保持完全相同的生成执行。
- `Array.getElem!` 的成员性通过边界证明转换到 `toList`，保证全项次数界可应用。

## 涉及文件

- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `docs/devlog/2026-08-11-hensel-divmod-raw-semantic-refinement.md`
