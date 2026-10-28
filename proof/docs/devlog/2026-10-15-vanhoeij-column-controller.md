# van Hoeij 列循环控制器与实际执行来源

日期：2026-10-15。

## 变更

- [独立生成模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_controller.py) 与 [控制器模型](../../lean/CLPoly/Generated/StrictVanHoeijController.lean)：实际逐轮初始化、构格、LLL、范数保留、候选提取、完整验证；验证成功返回，失败则到全宽停止或更新 target。
- [具体执行实例](../../lean/CLPoly/Refinement/VanHoeijControllerExecution.lean)：由构格合法性和机器整数递增证明实例化全部擦除的证明参数。没有可选择候选或结果的回调。
- [单轮来源](../../lean/CLPoly/Refinement/VanHoeijColumnTrace.lean)：追踪同一轮实际构格、LLL、提取和验证运行等式。
- [循环来源](../../lean/CLPoly/Refinement/VanHoeijColumnLoop.lean)：失败验证保持 complete；从未完成开始的成功循环结果必须来自某轮完整验证接受。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录循环变式、证明参数来源和成功来源归纳。

## 决策与范围

控制器独立于稳定辅助模型，避免每个新增控制步骤重建全部不可约性链。
保留 added、retained.size 的 Int32 转换和统计字段更新。构格仍用现有
自然数索引辅助模型，其与源机器整数对应所需范围条件不能省略。
尚未建模全入口参数选择、CLD 宽度扫描、三次精度重试及 fallback。
本模型是显式手写模板，不是 AST 翻译或 C++ 等价性证明；新控制器尚无
B2B 对照，旧辅助函数 750/750 不能充当新循环验收。

单轮/循环来源定理已由后续
[列循环不可约性组合](2026-10-15-vanhoeij-column-irreducibility.md)
连接到既有不可约性入口；全入口控制器仍未完成。
有限列循环终止不能代替充分 Hensel 精度下最终分离的证明。

## 实现问题

带初始化证明的依赖 match 不能靠外部 case 等式直接 simp 消去；证明改为
先展开实际函数，再逐层 split，保留该分支自带的执行等式。局部 bound
定义在匹配 LLL 结果前显式展开，避免简化器把相同算式视作不同模式。
验证器失败保留状态后，先对完整等式取 complete 投影，再化简统计更新。

## 验证记录

控制器和具体执行实例构建通过（3538 jobs）；单轮/循环来源定理定向构建
通过（3565 jobs），日志为 `/tmp/clpoly-vh-column-loop-build.log`。
模板 `--check` 与 `git diff --check` 已通过；本轮未改生产代码。
具体执行实例、单轮来源、失败保持及循环成功来源的公理审计通过，均只
依赖 `propext`、`Classical.choice`、`Quot.sound`，无新增公理或 sorryAx。
审计日志为 `/tmp/clpoly-vh-column-controller-audit.log`。
剩余目标保持，包括充分精度终止、参考实现对照、性能验收和最终合并。
