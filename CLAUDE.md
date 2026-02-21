# CLPoly 项目规范

## 项目定位

CLPoly 是一个 C++ 多变量多项式计算库，采用泛型模板设计，支持不同单项式序和系数域（Z, Q, Z/pZ）。

核心能力：多项式算术、GCD、因式分解（单变量）、结式/子结式、特征列、实根隔离。

对比对象：FLINT、NTL、Singular（开源库，详见 README.md 对比表）；Maple、Mathematica（商业系统，作为正确性参考和算法调研来源）。

目标：成为一个通用、高效、完善的 CAS 代数库。

## 构建指令

```bash
make                        # 构建全部库 (lib/ 下 debug/release × .a/.so)
make test/test_name         # 构建单个测试 → _build/debug/bin/test_name
_build/debug/bin/test_name  # 运行单个测试
bash test/run_all_tests.sh  # 构建并运行全部测试
make bench-clpoly           # 构建并运行性能基准 (release)
make clean                  # 清除所有构建产物
```

## 工作流程

遵循 [通用开发工作流程规范](docs/workflow.md)，涵盖设计三阶段、实现流程、修复迭代的完整规范。
