# 参数块与宽度扫描的生产入口 B2B

日期：2026-10-15。

## 变更与目的

- [生产代码](../../../clpoly/polynomial_factorize_univar.hh)：提取实际参数选择块及 CLD 宽度扫描，原生产入口调用这两个辅助函数。保持原分支顺序、ZZ 运算和机器整数转换；参数失败不发布部分输出。
- [C++ 注册表](../../b2b/registry/cpp_registry.cc) 与 [Lean 注册表](../../lean/B2B/Registry.lean)：直接比较新生产辅助函数与已有严格模型。参数失败编码为空数组；成功比较 D、m'、T 的精确整数值。C++ 哨兵状态另检查失败不写输出。
- [向量生成器](../../b2b/generate_vanhoeij_parameters_vectors.py)：闭式最小二次幂/整数位数公式作为独立预期值，不复写源倍增循环。包含 1264 个参数用例及 108 个宽度用例。
- [参数向量](../../b2b/vectors/__vanhoeij_choose_parameters.json) 与 [宽度向量](../../b2b/vectors/__vanhoeij_cld_width.json)：覆盖阈值两侧、虚拟位数边界、大整数、r=24/32、Int32 最大计数、空 CLD、稀疏多项式及最大安全宽度。

## 验证与发现

首次新旧十一组 B2B 合计 2122/2122 通过。检查框架发现 JSON null 作为
expected 时跳过独立预期比较；因此参数拒绝值改用 []，以覆盖拒绝分支
的独立预期校验。修订后的新两组重跑 1372/1372 通过，旧九组源/模型
未再变动。最终日志为 `/tmp/clpoly-vh-parameters-b2b-final.log`；首轮
十一组合并日志为 `/tmp/clpoly-vh-parameters-b2b.log`。

原生重组回归已重新编译运行，382/382 通过，日志为
`/tmp/clpoly-vh-parameters-native.log`。首个编译命令的 GMP include 路径
拼写错误已修正，成功结果来自修正后新编译的二进制。
Lean 注册表/驱动构建通过（3552 jobs）。向量生成器及严格模板的
`--check`、`git diff --check` 均通过。

## 范围

本轮 B2B 比较实际生产调用的辅助函数，不是测试中另抄一份参数块；
但它仍不是全入口控制器验收或 C++/Lean 等价性证明。宽度扫描的机器
范围对应证明、全入口、精度重试、充分精度最终分离、参考对照和性能
验收仍未完成。未合并。
