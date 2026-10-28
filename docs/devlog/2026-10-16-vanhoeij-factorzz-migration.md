# 外层 squarefree FactorZZ 迁移

日期：2026-10-16。

## 变更与原因

外层实际素数选择、单模因子早退和新的 van Hoeij 重试现已组装到同一
实际返回数组的乘积还原及不可约性。公开 refinement 和 Pipeline 切换
到此实现，避免用旧数量 safety net 控制器的定理作为新代码的证据。

删除生产头文件中已无调用者的数量 helper、四项原生断言及旧 B2B
接口。新增九个已知不可约分解族，各覆盖小素数和大素数模式；两端
各自检查已知因子，再保留原输出顺序作 C++/Lean 对照。

## 验证与边界

- 新原生回归 378/378；减少的四项是已删除 helper 的断言。
- 公开 Pipeline、最终公理审计、新链条审计、B2B 驱动构建成功：3615 jobs。
- 新链条 21 项公理审计只依赖 `propext`、`Classical.choice`、`Quot.sound`。
- 严格边界检查、公开 refinement 注册检查和新模板/向量重现检查通过。
- 外层与重试合计 50 项 B2B 全部通过：旧外层 5 项、新外层 18 项、重试 27 项。

证明仍明确保留源输入约束、物理算术实例和底层递归 Hensel 执行就绪
条件。手写控制流模板不是通用 AST 翻译。有限重试后回退也不是纯
van Hoeij 充分精度必然分离定理；参考对照及性能仍待完成。随后修正
实根头文件定义顺序，完整库构建及全部 27 个原生测试文件均已通过，
详见[原生回归日志](2026-10-16-vanhoeij-full-native-regression.md)。

## 关键决策

素性参数仅作擦除证明，不增加生产代码没有的运行分支。严格边界检查
拒绝了与旧输出回调同名的证明字段；将字段改名为
`zassenhausTermination`，保留检查规则，随后重新构建通过。

涉及：生产头文件、原生重组测试、两个 B2B registry、StrictRuntime、
新外层模板/模型/refinement、公开生成器注册和 Pipeline、两份公理审计。
详细路径及证明度量见[证明日志](../../proof/docs/devlog/2026-10-16-vanhoeij-factorzz-migration.md)。
