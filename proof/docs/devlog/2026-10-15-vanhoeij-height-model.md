# van Hoeij 高度 helper 严格模型

日期：2026-10-15。目标仍 active，未合并。

## 本轮变化

严格模板新增 L1 逐项循环及高度入口，保留生产代码的系数复制、
条件取负、累加顺序、次数向 Int32 的转换、断言及整数左移。
未将机器次数计算替换为无界自然数公式。该模板仍为显式模型，
不是 AST 翻译等价证明。

证明 L1 循环任意起点的后缀绝对值和不变量，以及高度成功运行的
实际输出式。后续须利用规范稀疏表示和机器次数范围，把数组求和及
转换次数接到上一轮使用的数学高度公式；当前不声称该最后一步完成。

## B2B

注册并调用实际生产 __upoly_norm_l1 和 __vanhoeij_cld_bound。
独立 Python 整数运算生成期望值：L1 60 例、高度 58 例，共 118/118
通过。覆盖空/常数多项式（仅 L1）、正负系数、大整数、稀疏输入和
次数 1、24、32、512、4096。B2B 高度入口限制测试次数 1..4096，
避免测试构造产生无界大输出；该测试限制不属于算法模型前提。

生成模板 --check 及 git diff --check 通过。未改生产 C++；新增
测试重新编译完整 C++ B2B driver。L1 证明单文件检查通过。
完整证明目标构建通过（3555 jobs，含已有依赖）；三个定理的公理
审计仅含 propext、Classical.choice、Quot.sound，无 sorry 或新公理。
新增证明模块 47 行。总耗时未精确记录。

## 文件

- [严格模板](../../cpp2lean_v2/tests/build_strict_vanhoeij.py)
- [高度证明](../../lean/CLPoly/Refinement/VanHoeijHeight.lean)
- [L1 向量](../../b2b/vectors/__upoly_norm_l1.json)
- [高度向量](../../b2b/vectors/__vanhoeij_cld_bound.json)
