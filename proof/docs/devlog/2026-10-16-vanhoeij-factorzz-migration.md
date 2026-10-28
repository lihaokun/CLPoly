# 外层 squarefree FactorZZ 及公开证明 API 迁移

日期：2026-10-16。

## 内容

新增 [外层模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_factor_zz.py)、
[严格模型](../../lean/CLPoly/Generated/StrictVanHoeijFactorZZ.lean) 和
[具体细化](../../lean/CLPoly/Refinement/VanHoeijFactorZZ.lean)。先在
[自然语言草稿](../../nl-proof/vanhoeij-certificate-repair.md) 展开分支证明，
再形式化实际 selectPrime、早退及新 retry 的组装。

成功选择提供同一 selection 的素性和语义/物理条件；不可约标志来源
保证早退确有至多一个模因子。其余分支实例化按选中素数索引的稠密
运算和工作区，再应用新的 retry 正确性。没有输出不可约性回调，
没有输出计数、次数或恢复精度假设。

更新 [注册表](../../cpp2lean_v2/class_map.py)、
[公开包装生成器](../../cpp2lean_v2/passes/pass9_refine_gen.py)、
[生成包装](../../lean/CLPoly/Refinement/Generated.lean)、
[Pipeline](../../lean/CLPoly/Pipeline/L1.lean)。最终审计移除旧控制器
项，保留仍用到的旧 Hensel/恢复/Zassenhaus 数学引理。

[StrictRuntime](../../lean/B2B/StrictRuntime.lean) 执行新外层模型。
测试中擦除 Prop 的设施与正式证明隔离。新增
[向量生成器](../../b2b/generate_vanhoeij_factorzz_vectors.py) 和十八项
外层向量；两端都独立断言已知因子。

## 验证

- 外层模型、具体细化及首次公开 Pipeline 构建通过。
- 证明字段重命名后最终构建再次通过，3615 jobs。
- [新链条审计](../../lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)
  21 项仅有标准 Lean/Mathlib 公理；最终公开 theorem 审计也无新公理。
- 严格边界、公开注册、模板/向量重现检查通过。
- profiling 模式原生回归 378/378；50 项 B2B 全部通过（外层 23、重试 27）。

## 边界

此为条件式实际执行 refinement，仍要求输入规范、本原、次数范围和
底层递归 Hensel 执行就绪。控制流模型由手写模板生成；并未声称通用
AST 翻译等价。纯 van Hoeij 充分精度分离尚未完成。

## 度量

- 耗时：本步骤墙钟主要包含大素数解释执行；未单独计时，不补造精确数字。
- 迭代：模型、细化、公开包装首次构建通过；字段重命名后重新构建一次。
- Lean 新增行数：两份新外层模块合计 144 行，另有公开接线与审计修改。
- 对应 C++：外层约 15 行核心控制流；另删除旧 helper 和旧断言。
- 放弃的方案：弱化严格边界检查；改为消除证明字段与旧回调的命名冲突。
