# 完整尝试入口的成功来源及不可约性

日期：2026-10-15。

## 变更

- [完成来源](../../lean/CLPoly/Refinement/VanHoeijEntryTrace.lean)：展开实际入口运行，取得成功高度、单元素组构造和实际验证等式。最终完成必来自单元素成功验证，或从未完成状态开始的实际列循环；参数不足返回及所有 error 分支不能伪造完成。
- [入口不可约性](../../lean/CLPoly/Refinement/VanHoeijEntryIrreducibility.lean)：两种实际来源分别接入既有单元素/列循环不可约性定理。参数选择成功提供 H<m、D>0、模数与残差公式，构格合法性来自具体执行，均不再作为入口调用者的独立假设。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录分支追踪与参数条件消去。

## 范围与剩余义务

在真实选素数/Hensel 正确性、规范表示、本原输入以及明确次数/计数
范围条件下，实际 try-recombine 返回 complete=true 推出乘积等于输入
且所有输出不可约。不假设任何返回因子不可约，不使用旧控制器证明。

该结论是实际尝试入口的条件正确性，尚不是全部目标完成：lifted.size
的 Int32 范围、提升因子的机器次数范围仍明确保留，需继续从更上游
入口导出；全入口 C++/Lean B2B、外层精度重试、充分精度最终分离、
参考实现对照和性能验收仍待完成。手写模板不是 AST 等价性证明。

## 构建

两条新定理及依赖构建通过（3592 jobs），日志为
`/tmp/clpoly-vh-entry-irred-build.log`。展开源入口时先消去 count 的局部
let，再按外层分支展开，避免 split 先选到内部依赖 match。此调整只
影响证明策略，没有修改模型或生产代码。
两条新定理公理审计通过，仅依赖 `propext`、`Classical.choice`、
`Quot.sound`，没有新增公理或 sorryAx；日志为
`/tmp/clpoly-vh-entry-irred-audit.log`。`git diff --check` 通过。
后续完整尝试入口 51/51 生产 B2B 已通过，见
[入口 B2B 日志](2026-10-16-vanhoeij-entry-b2b.md)。
