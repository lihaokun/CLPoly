# HGCD 递归基础分支统一不变量

日期：2026-08-09

## 做了什么

- 保留 `_hgcd_recursive` 矩阵基础分支真实 `_mat_one` 初始化已经建立的恒等变换、带符号行列式和矩阵长度界。
- 将 `computeM = true` 的真实矩阵初始化及 A/B 拷贝封装为 `HgcdRecursiveRawInvariant`。
- 将 `computeM = false` 的真实 A/B 拷贝封装为同一不变量；矩阵语义与长度义务仅在矩阵实际被请求时成立。
- 增加按真实 `computeM` 分支统一入口的基础情形定理，供完整递归主体直接调用。
- 将第一递归调用后的真实 early-return 拷贝出口封装为统一递归不变量；矩阵描述符长度通过实际四项矩阵拷贝保持。
- 证明带符号行列式控制的两套 C++ 低半重构公式互为实际 HGCD 矩阵的逆，并在加回移位高半后恢复完整多项式 transform。
- 把两个统一基础分支定理加入严格 HGCD 源码门禁。

## 为什么做

良基递归总精化需要所有返回分支提供同一种语义归纳结果。此前基础分支只导出了部分 raw 表示，丢弃了矩阵初始化已经证明的语义和长度证据，无法直接作为递归主体的基础情形。

## 关键决策

- 矩阵启用分支的所有语义均来自真实 `hgcdRecursiveBase` 执行以及已有的 `hgcdIterInit_refines`，没有另行计算规格结果。
- 矩阵禁用分支不声称输入矩阵具有任何代数语义，只保留源程序确实执行的两次拷贝及 GCD 不变性。
- 长度归纳信息来自真实 `_mat_one` 描述符长度和拷贝后的返回长度。

## 验证

- `lake build CLPoly.Impl.StrictHGCDRawRefinement`
- 严格 HGCD 源码门禁（在提交前运行）
- 新增定理公理审计（在提交前运行）

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-hgcd-recursive-base-invariant.md`
