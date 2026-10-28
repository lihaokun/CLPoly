# van Hoeij 原生证书路径

日期：2026-10-15。用户已确认修复方案。当前是研究分支的阶段性实现，尚不可合并。

## 变更与原因

- 修正 CLD 格方向和模数行，使用 `[sI | C; 0 | mI]`。
- 引入实际 CLD 系数界、同步截断及误差界，依据 GS 范数删除整个超界后缀，保留全部相关短向量。
- 完整验证整个列分区后才返回成功；不再把结果数量少于模因子数量视为失败。
- 主控制器最多进行三次精度尝试，最终未完成时仅一次全精度 Zassenhaus；计时覆盖回退，增加实际路径计数。
- 新增已知不可约因子、低精度、非首一和大子集重组测试。

## 决策与限制

- 保留每轮规范重建，避免坐标来源混淆。
- 复用 LLL 已计算的精确 GS 范数，另外独立重算用于回归对照。
- r+16 虚拟精度、两列起步；降低虚拟精度和单列起步的实测更慢，已撤回。
- 三次尝试是工程回退预算，不冒充纯 van Hoeij 的终止性证明。
- 旧数量 helper 暂留给旧 B2B 接口，生产控制器已不调用。生成模型与端到端证明尚未同步，因此不合回原分支。

## 验证

重组 382/382；既有因式分解测试在仅缩窄包含头文件的临时副本上 256/256。新增数学引理构建与标准公理审计通过。性能及环境细节见[测试报告](../test-reports/vanhoeij-certificate-native.md)。

系统工具链缺失通过现有 Nix 安装解决。原测试总头文件引入已有 realroot 不完整类型错误，未修改该无关模块。

## 涉及文件

- [原生实现](../../clpoly/polynomial_factorize_univar.hh)
- [回归](../../test/test_recombine.cc)、[基准](../../test/bench_vanhoeij_certificate.cc)
- [方案](../fixes/vanhoeij-fix-lattice-certificate.md)、[测试报告](../test-reports/vanhoeij-certificate-native.md)
- [架构](../design/vanhoeij/architecture.md)、[细化设计](../design/vanhoeij/detailed-design.md)、[旧调研更正](../research/vanhoeij-lll-research.md)
- [自然语言证明](../../proof/nl-proof/vanhoeij-certificate-repair.md)、[CLD 界](../../proof/lean/CLPoly/Algorithm/VanHoeijBound.lean)、[证书引理](../../proof/lean/CLPoly/Refinement/VanHoeijCertificate.lean)

## 度量

- 耗时：本轮未记录精确墙钟起止，不估填。
- 迭代：多轮原生编译、精度/列数实验及 Lean 编译；关键失败和撤回原因已在上文记录。
- Lean 新增：279 行（两个证明文件）；严格生成控制器尚未修改。
- 对应 C++：核心矩阵、证书和控制器约 440 行 diff（包含删除的旧控制器）。
- 放弃的方案：普通行范数过滤、按数量判完整、较低虚拟精度、单列起步。
