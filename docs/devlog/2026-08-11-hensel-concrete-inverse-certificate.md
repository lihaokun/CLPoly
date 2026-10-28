# Hensel concrete modular inverse certificate

## 做了什么

- 为 `ZZ.invert` 的实际 `invertImpl` 实现证明 success 分支的模逆元正确性。
- 证明 Hensel modular divmod 中实际计算的 `divmodCoefficient` 会消去除数首项系数。

## 为什么做

外层 `DivmodTrace` 的商余守恒依赖首项消去。如果直接假设 inverse 正确，就会在 C++ L1 与 Lean L2 之间留下规格 oracle。本步骤展开实际扩展欧几里得实现，并从已有 Bézout 等式推出 `ZMod` 逆元性质。

## 关键决策

- success bit、inverse 数值都取自生成代码实际调用的 `ZZ.invert`。
- 分别排除无效模数、零剩余和 gcd 非一的失败分支。
- success 分支通过 `Nat.extGcd_bezout` 转入 `ZMod`，不新增公理或抽象 inverse provider。

## 遇到的问题与解决方式

- `invertImpl` 会先做两次 floor remainder 规范化；利用整数模运算化简到同一个 residue。
- Bézout 使用 `natAbs`，通过 residue 非负性恢复整数值后再 cast 到 `ZMod`。

## 涉及文件

- `proof/lean/CLPoly/Math/Bigint.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `docs/devlog/2026-08-11-hensel-concrete-inverse-certificate.md`
