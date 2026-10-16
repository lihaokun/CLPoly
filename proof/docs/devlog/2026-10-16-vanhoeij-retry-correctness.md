# 新重试入口正确性组装及输出次数前提消去

日期：2026-10-16。

## 证明结构

先更新[自然语言草稿](../../nl-proof/vanhoeij-certificate-repair.md)，再实现：

1. [输出次数](../../lean/CLPoly/Refinement/VanHoeijHenselDegree.lean)：实际首一
   Hensel 因子模 p 后次数不变；模乘积相伴给出每个输出模因子整除 f mod p；
   后者非零，整除次数界与映射次数不增合成输出次数≤输入次数。物理首项
   和数学次数对应后，消去[尝试入口](../../lean/CLPoly/Refinement/VanHoeijEntryIrreducibility.lean)
   独立的输出首项 Int64 范围前提。
2. [Hensel 连接](../../lean/CLPoly/Refinement/VanHoeijRetryHensel.lean)：固定真实
   算术与树构造，从成功运行得到计数与指数检查；源表示和选择正确性建立
   调整规范性、非空性、首一性、互素性。与已有递归执行就绪条件构成入口
   不变量，实际细化及运行确定性给出同一返回数组的语义证书。
3. [重试入口](../../lean/CLPoly/Refinement/VanHoeijRetryCorrectness.lean)：加强
   来源定理区分完成和回退。完成支从真实 Hensel、单次尝试证书推出乘积
   及不可约性；回退支由起始精度非负和启发式 min 推出 full 非负，再以
   请求指数≥full 得到恢复界，连接实际 Zassenhaus 正确性。

新入口定理不假设输出不可约、输出次数、输出计数或回退恢复精度。
保留的条件是源规范性、非空、本原性、输入 degree<2^31、选择正确性、
稠密算术配置，以及仅对实际成功提升要求的递归 Hensel 执行就绪条件。
不能据此声称无条件总正确性、纯 van Hoeij 最终分离或普遍 C++/Lean 等价。

## 验证

- 次数界：3547 jobs 构建通过。
- 新重试正确性：3602 jobs 构建通过，日志
  `/tmp/clpoly-vh-retry-correctness-build.log`。
- 扩展[公理审计](../../lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)：3603 jobs
  构建通过，18 项仅依赖 `propext`、`Classical.choice`、`Quot.sound` 的
  子集，无新增公理或 sorryAx。日志 `/tmp/clpoly-vh-retry-correctness-audit.log`。
- 本轮未改运行函数；完整重试 27/27 B2B 的已通过证据继续适用。

## 工程问题

匿名 letI 默认名称遮蔽稠密上下文参数 this，改为命名 densePrime。
素数域实例依赖具体素数，直接改写结构投影等式不产生合法 motive；
把运输独立为两个素数参数上的互素引理，再代入投影，避免对大型执行
上下文进行依赖消去。主重试定理首次到达编译检查即通过。

## 度量

- 耗时：约 0.3 小时（推导、类型转换修正及构建）。
- 迭代：次数模块 1 次；Hensel 连接/组装 6 次相关构建；公理审计单独构建。
- Lean 新增：270 行（三个模块），另调整尝试定理和审计清单。
- 对应 C++：生产代码修改 0 行；验证已有新重试及其 Hensel 调用。
- 放弃方案：在含素数域实例的大上下文直接 rw/cases 结构投影等式，
  改为独立等式运输引理；未改变算法或增加语义假设。
