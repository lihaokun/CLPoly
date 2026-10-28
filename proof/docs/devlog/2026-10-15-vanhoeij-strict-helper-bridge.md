# van Hoeij 严格 helper 模型与安全保留桥接

日期：2026-10-15。用户已确认修复，本轮继续执行已批准方案。

## 变更与动机

旧 StrictRecombine 中的矩阵与顶层控制器仍描述修复前算法。新增独立
StrictVanHoeij 模块作为迁移目标，避免在旧定理仍依赖旧布局时将其编译成功
误报为新实现已经验证。生成脚本是显式模板，不是 Clang AST 翻译器；模板
一致性检查也不是 C++/Lean 等价证明。

新模型覆盖：

- 在已有因子行写入 CLD 列，并追加对角元为实际 modulus 的新行。
- 螺旋顺序和分批追加的 bounded loop。
- descending while 保留前缀及分配后 iota 写入。
- LLL 导出实际最终 GS 缓存，保留原有 short-row 收集调用。

新增机器检查证明：

- 成功填列保持行数；成功追加使行数增一，新增行指定坐标等于 modulus。
- iota 输出恰为连续前缀；所有遗漏下标的 GS 范数严格大于界。
- 导出的范数来自实际 LLL 执行状态，满足已有具体 GS 不变式。
- 将该执行状态与安全移除连接：对返回矩阵的任一界内整数线性组合，
  每个非零坐标均出现在实际 helper 返回的行号中。无需调用者假设 GS
  正确性或被删后缀的范数性质。

## 验证

- 新增 190 个 B2B 用例全部通过：106 个前缀用例、60 个矩阵追加用例、
  24 个 LLL 缓存用例。前两类附独立预期值；LLL 比较最终矩阵、变换矩阵
  和精确有理范数。忽略调用者不用的旧 short-row 返回值，因 std::sort
  对等范数行的次序不作保证。
- LLL 用例采用非零对角的三角矩阵加整数行操作，满足满秩前提；B2B
  运行时沿用现有 erasedValue 技术，证明模块不导入 B2B。
- 单文件 Lean 检查退出 0；相关目标 lake build 成功（3553 jobs，含已有依赖）；生成模板 --check 通过。
- 六个新增关键定理的 `#print axioms` 均仅包含 `propext`、`Classical.choice`、`Quot.sound`；无 sorry 或新不可约性假设。
- 曾在 C++ 重编译完成前启动 runner，导致 24 个新注册入口报告 unknown
  function；等待构建完成后重新运行 190/190 通过。这是旧进程问题，不是
  通过改动测试绕过失败。

## 尚未完成

严格顶层尚未切换到新模块。完整 CLD 块布局、满秩、真因子短向量与初始
格的连接、完整分组验证、三次提升控制器和端到端不可约性仍须继续迁移。
当前不能合并或宣称整个新 C++ 已形式验证。纯 van Hoeij 最终分离定理与
原 r=24 性能回归也未解决。

## 度量

- 耗时：本轮未记录精确起止时间，不估造。
- 迭代：多轮单文件/相关目标编译修复，未独立记录准确轮数；另有 B2B 自动构建。
- Lean 行数：新严格模块 97 行；证书模块由 150 行增长到 324 行。
- 对应 C++：本轮未改生产算法，覆盖矩阵构造、前缀选择和 LLL 缓存导出。
- 放弃方案：未直接替换旧矩阵函数，避免旧布局定理与新实现混用；后续
  必须显式迁移调用方，独立模块不是最终闭环。

## 文件

- [生成模板](../../cpp2lean_v2/tests/build_strict_vanhoeij.py)
- [严格模型](../../lean/CLPoly/Generated/StrictVanHoeij.lean)
- [证书证明](../../lean/CLPoly/Refinement/VanHoeijCertificate.lean)
- [自然语言草稿](../../nl-proof/vanhoeij-certificate-repair.md)
- [C++ 注册表](../../b2b/registry/cpp_registry.cc)
- [Lean 注册表](../../lean/B2B/Registry.lean)
- [B2B 用例](../../b2b/vectors/)
