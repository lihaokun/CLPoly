# Strict dense `__strip` 良基翻译

## 日期

2026-08-07

## 做了什么

- 将 `dense_upoly_zp::empty` 和 `dense_upoly_zp::__strip` 加入严格 GCD 生成闭包。
- 直接从 `clpoly/dense_upoly_zp.hh` 的 Clang AST 生成 Lean 控制流。
- 将 `__strip` 的 `while` 生成为带 `termination_by` 的总函数。
- 严格生成门继续拒绝 `sorry`、`partial def` 和未解析调用。

## 为什么做

`divrem` 在每次商和余数计算完成后调用 `__strip`。若跳过该函数或用规范层操作代替，就不能证明后续 GCD 是真实 C++ L1 执行到 Lean L2 的精化。

## 自然语言证明草稿

循环条件为系数数组非空且末项为零。条件成立时，C++ 循环体恰好执行一次 `pop_back`；因此数组长度由 `n` 变为 `n - 1`。非空条件给出 `n ≠ 0`，所以 `n - 1 < n`。以数组长度为自然数度量，递归调用严格下降，故循环良基终止。Lean 生成代码不增加任何运行时默认分支。

## 关键决策

- 度量使用当前 `this._coeffs.size`，而不是 fuel；它直接对应每轮 `pop_back` 的状态变化。
- 将 C++ 循环 guard 命名为 `h_strip`，只用于证明递减，不改变可执行控制流。
- 不调用 L2 多项式正规化函数；结果仍由生成的 C++ 数组操作产生。

## 问题与解决

初次生成后 Lean 需要显式证明 `size - 1 < size`。通过 `Bool.and_eq_true` 提取原始 guard 的非空分量，再用 `Array.size_eq_zero_iff` 排除长度为零，最后用 `Array.size_pop` 和 `omega` 关闭义务。

## 涉及文件

- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/cpp2lean_v2/tests/build_strict_gcd.py`
- `proof/lean/CLPoly/Generated/StrictGCD.lean`

## 度量

- 耗时：约 0.5 小时
- 迭代：4 轮生成、编译与递减证明修复
- Lean 新增/修改行数：严格生成文件增加约 18 行
- 对应 C++ 行数：`empty` 1 行，`__strip` 4 行
- 放弃的方案：`partial def` 与 fuel 均未采用，因为用户要求良基递归且严格精化链禁止非总执行替代。
