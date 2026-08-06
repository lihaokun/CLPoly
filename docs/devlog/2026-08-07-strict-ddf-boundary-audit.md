# Strict DDF 边界审计

## 结论

`Generated.StrictDDF` 先前仍包含一条不合格路径：生成的
`__upoly_mod_ir` 经 `pair_vec_div5` 和 `HasPolyDivmod SparsePolyZp` 直接调用
手写 L2 `SparsePolyZp.divmod`。这不是 C++ L1 到 Lean L2 的精化。

本次从 strict 生成根中移除了 `__upoly_mod`、`__upoly_powmod`，并进一步
移除了尚依赖 `Array.get!`、`Array.set!`、`front!` 的 monic/subtract-x 操作。
当前 strict DDF 只保留无上述退化路径的 `__make_zp_ir`。

生成检查器现在显式拒绝多项式除法 typeclass dispatch、手写 L2 divmod，
以及带默认越界语义的数组访问。被移除的操作只有在 RawHeap 除法语义和
数组边界不变量闭合后才能重新进入 strict 边界。

## 验证

- `build_strict_ddf.py --check` 通过；
- `lake build CLPoly.Generated.StrictDDF CLPoly.Refinement.DDF` 通过；
- 没有使用 fuel、oracle、默认结果或 L2 算法替代 C++ 执行。
