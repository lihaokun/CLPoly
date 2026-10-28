# van Hoeij 格构造研究

- 日期：2026-10-15。
- 从干净提交 `428a27fe2c22dcbeee100d1100e310d8a08f1e77` 创建 `research/vanhoeij-lattice-certificate`。
- 静态确认数量 safety net 会重算正常的模因子组合结果；格构造忽略模数，短行筛选缺少完备性依据。
- 读取 FLINT 三个实际源码入口，对照加入模约束、列分组与提升循环；确认旧调研不能作为修复依据。
- 核实旧算法层定理接收输出不可约性前提；未重新运行整个 PR 的证明。
- 决策：保留当前正确性保障，提出有真因子向量保留不变量及完成证书的修正方案，待确认后实现。
- 环境：系统 Git/clang++ 因 Command Line Tools 缺失不可用；通过 Nix 用户 Git 成功创建分支。未运行 C++/Lean 或性能基准。
- 文件：[修正方案](../fixes/vanhoeij-fix-lattice-certificate.md)；本日志。未修改生产实现或证明。
