# Hensel modular divmod merge refinement

## 做了什么

- 将生成语义中的 `pushNonzero` 公开为可证明的具体 L1 helper。
- 定义数组后缀的 `ZMod` 多项式解码 `termsToPolyMod`。
- 对生成的良基递归 `divmodMergeLoop` 的六个控制流分支完成直接语义证明。
- 将内层不变式提升到 `divmodRemainder`，并证明首项消去条件下它等于完整多项式减法。

## 为什么做

`__upoly_divmod_mod_raw_ir` 是 `__hensel_step` 两次模多项式除法的真实执行边界。只有直接证明其数组合并循环，外层 `DivmodTrace` 才能建立真正的 C++ L1 到 Lean L2 商余语义；不能用规格 oracle 或 L2 除法结果替代。

## 关键决策

- 证明严格跟踪 `rIndex`、`gIndex` 对应的 `Array.toList.drop` 后缀。
- 使用生成函数自身的 `.induct` 原理，因此归纳顺序与 C++ lowering 的良基递归完全一致。
- 分别处理复制余项、补负乘积项、两种次数比较、同次数相减和结束分支。
- 不引入 fuel，也不引入抽象商余式输出作为执行结果。

## 遇到的问题与解决方式

- `Array.getElem!` 与列表后缀首项的表示不同：建立受边界证明保护的精确索引桥。
- GMP floor remainder 在 `ZMod` 中需要消去：分别证明正项、负乘积项和差项的 cast 引理。
- 同次数分支的系数合并通过 `Polynomial.monomial` 的加法同态与环正规化闭合。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `docs/devlog/2026-08-11-hensel-divmod-merge-refinement.md`
