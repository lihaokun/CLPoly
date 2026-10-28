# 完整原生构建与回归恢复

日期：2026-10-16。

## 原因与修复

van Hoeij 原生验收此前只能用临时缩窄头文件的因式分解测试副本，
因为 [realroot.hh](../../clpoly/realroot.hh) 的模板在 `uroot` 仍是不完整
类型时构造它，Clang 拒绝编译。将该模板原样移动到类定义之后，删除
因此多余的前置声明，并清理移动区域的行尾空格。没有改变函数逻辑。

## 验证

- 使用已有 Nix clang 19.1.7、cctools、GMP 6.3.0、Boost 1.87.0。
- 标准 `make -j2` 成功，debug/release 静态库与共享库均构建成功。
- 仓库原样 [run_all_tests.sh](../../test/run_all_tests.sh) 成功：
  **27 个测试文件通过，0 失败**。
- 其中原样因式分解测试 256/256，实根隔离 32/32，普通模式重组 376/376。
  单独 profiling 模式重组为 378/378，多出的两项是路径计数断言。
- 只通过环境提供编译器、include/lib 路径及现有归档工具；未更改系统
  工具链、测试断言或测试脚本。

日志：`/tmp/clpoly-vh-full-native-build.log`、
`/tmp/clpoly-vh-full-native-tests.log`。

这消除了先前完整原生构建的验收障碍；纯 van Hoeij 分离终止定理、
固定参考实现的动态对照及性能验收仍是合并前义务。
