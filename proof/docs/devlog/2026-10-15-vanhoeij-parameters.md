# van Hoeij 参数选择及除数循环

日期：2026-10-15。

## 变更

- [参数模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_parameters.py) 与 [生成模型](../../lean/CLPoly/Generated/StrictVanHoeijParameters.lean)：保留 m>H 检查、D 从 1 按 2 倍增、UInt64 虚拟位数及模数位数计算、同一除数下的模数和残差计算，以及两条信息量检查。另加入实际 Int32 宽度扫描。
- [参数证明](../../lean/CLPoly/Refinement/VanHoeijParameters.lean)：倍增循环正性及退出界、最终除数正性，以及实际成功返回的全部参数等式和精度不等式。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录倍增循环变式、max 保持正性与成功分支推理。

## 证明与实现决定

倍增循环以 `(height-divisor).toNat` 严格递减，正性参数被擦除，不能
选择输出。规格证明使用实际函数的归纳原理，不另建替代循环。
轻量生成模型中的 ZZ 类型别名需要在算术证明中显式展开；只改证明
归一化，没有加入断言、修改源算法或扩大递归限制。

宽度扫描保留 Int64 到 Int32 的转换和加一；其无溢出条件与自然数
宽度对应仍待证明。目标平台 size_t 按 UInt64 建模。模型仍为手写
控制流模板，不能称为 AST 翻译或 C++ 等价性证明。

## 验证与剩余义务

模型及规格定向构建通过（601 jobs），日志为
`/tmp/clpoly-vh-parameters-build.log`。本轮没有修改生产代码，也没有
新增参数块 B2B 验收；需在后续把相同源块接入可直接对照的测试入口。
后续已完成源块提取与对照接入，见
[参数生产 B2B 日志](2026-10-15-vanhoeij-parameters-b2b.md)。
倍增函数及三条规格定理公理审计通过，仅依赖 `propext`、`Quot.sound`，
无新增公理或 sorryAx，日志为 `/tmp/clpoly-vh-parameters-audit.log`。
模板 `--check` 与 `git diff --check` 通过。
完整入口、精度重试、充分精度最终分离、参考对照及性能验收仍待完成。
