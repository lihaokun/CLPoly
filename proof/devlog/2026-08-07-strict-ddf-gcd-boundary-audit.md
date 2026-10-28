# Strict DDF GCD 边界审计

## 结论

`Generated.StrictDDF._loop___ddf_Zp_0_ir` 曾直接调用根命名空间中的
`polynomial_GCD`。该名字由 `class_map.py` 的 `direct` 映射解析到
`CLPoly.Model.polynomial_GCD`，而不是 cpp2lean 生成的 C++ L1 定义。因此，
此前的 DDF 主循环虽然自身采用良基递归，仍会在 GCD 调用处退回 L2 实现，
不能作为严格 C++ L1→L2 精化的证明对象。

当前修复先从 StrictDDF 生成闭包中移除 `__ddf_Zp`，并加入生成期拒绝规则：
严格产物中只要出现 `polynomial_GCD` 就失败。已经闭合的 make-monic、powmod
和 subtract-X 生成控制流不受影响。DDF 主入口只会在 C++ GCD 依赖也进入
生成闭包后恢复。

## C++ AST 验证

Clang 可以取得两个所需的具体实例：

- `polynomial_GCD(const upolynomial_<Zp>&, const upolynomial_<Zp>&)`；
- `__polynomial_GCD(upolynomial_<Zp>&, ..., int64_t)`。

外层 `polynomial_GCD` 已能通过 cpp2lean Pass 1–8；其非空分支调用内层
`__polynomial_GCD`。内层也能通过控制流各 Pass，但当前类型映射缺少
`dense_upoly_zp`，生成结果会在稀疏/稠密转换、`gcd`、模乘/模逆、缩放和
稠密转稀疏处产生占位。因此下一阶段必须补齐 dense 类型及其实际 C++ 方法
闭包，不能把这些调用映射回 `SparsePolyZp.gcd`。

## 度量

- StrictDDF 可复现检查：通过；
- StrictDDF 中 `polynomial_GCD`：0；
- StrictDDF 中 `__ddf_Zp_ir`：0（在 GCD 闭包完成前隔离）；
- StrictDDF 中 `partial def` / `sorry` / `fuel` / `Safe`：0；
- `CLPoly.Refinement.DDF` 构建：通过；
- `CLPoly.Pipeline.L1` 构建：通过；
- 已识别的 inner-GCD 翻译缺口：3 个 dense constructor、6 个 dense method、
  1 个 `bool→int64` cast。
