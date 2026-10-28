# van Hoeij 完整分区验证模型

日期：2026-10-15。目标 active，未合并。

## 生产接口与模型

把生产验证 lambda 原样提取为 __vanhoeij_validate_partition，两个
生产调用点直接使用它。保留原分支、更新时机及排序，不写测试副本。
严格模板新增分组模乘积与 Int64 次数累加、完整分区循环、最终余商
检查及答案写回。真实调用整数试除和本原化，不包含结果 oracle。

答案结构同时跟踪 factors、complete、cldColumns、retainedRows。
失败不发布成功前缀，也不修改旧答案；成功仅替换 factors 并置 complete。
新增状态保持与成功更新证明，后续还需证明乘积/支持与不可约性。

## B2B 与回归

106 个独立构造的完整验证用例通过：空分区、空组、缺组、重复组、
后续组试除失败、最终余商非 1、零/常数输入、低精度失败、非首一
因子及随机合法分区。每例分别测试初始 complete=false/true，答案
带非空 sentinel 因子及统计字段，检查失败时保持原状态。

测试允许粗分组通过验证器，不错误要求验证器独立证明不可约性。
同次数因子的 std::sort tie 顺序未规定，B2B 两端只在序列化边界将
因子按一致键排序以比较多重集合；其他字段逐值比较。

生产代码变动后重新编译并运行原生重组测试：382/382。
模板 --check、git diff --check 通过。本次修复七个 helper 的统一 B2B
复跑 476/476 通过。完整证明目标构建通过（3560 jobs，含全部已有
依赖），新增状态证明模块 48 行。总耗时未精确记录，不估造。
两项状态定理公理审计仅含 propext、Quot.sound，无 sorry 或新公理。

## 文件

- [生产接口](../../../clpoly/polynomial_factorize_univar.hh)
- [严格模板](../../cpp2lean_v2/tests/build_strict_vanhoeij.py)
- [验证器证明](../../lean/CLPoly/Refinement/VanHoeijValidation.lean)
- [B2B 向量](../../b2b/vectors/__vanhoeij_validate_partition.json)
