# 显式充分精度阈值与实际分组完备性

日期：2026-10-16。

## 变更与理由

在既有实际 Hensel/CLD 分离输入证书上，先补充
[自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §11–13，
再形式化统一界、结式分离和实际分组循环。新增七个模块、33 项定理。

- [VanHoeijUniformBounds](../../lean/CLPoly/Refinement/VanHoeijUniformBounds.lean)：
  从实际完整列构造得到 `N≤r+n`，从实际参数选择得到 cutoff 的统一界
  `B0=r·32²+n·(32+2r)²`，再把实际 LLL 行坐标界统一为
  `K0=(r+n)·2^(r+n)·B0+1`，残差高度为 `D·(2+3(r+n))·K0`。
- [VanHoeijSeparationHeight](../../lean/CLPoly/Refinement/VanHoeijSeparationHeight.lean)：
  使用已有 Mignotte 定理控制真实因子系数，用源 L1 范数控制导数，
  对真实因子 F 与残差差值 `H−e_i f'` 给出统一系数界 B；结式阈值
  是 `(2n)! B^(2n)`。这里 B 仍含 D，没有把随精度增长的 D 偷换成常数。
- [VanHoeijRowSeparation](../../lean/CLPoly/Refinement/VanHoeijRowSeparation.lean)：
  模数超过结式阈值时，局部公共因子迫使整数结式为零；F 的本原性与
  Q 上不可约性推出 `F | H−e_i f'`。支持内另一局部因子结合导数互素性
  给出 `m | e_i−e_j`；再由 `m>2K0` 得 `e_i=e_j`。
- [VanHoeijActualSupportSeparation](../../lean/CLPoly/Refinement/VanHoeijActualSupportSeparation.lean)：
  支持由同一次实际 Hensel 执行导出，并接到同一次 LLL 返回的 transform。
  输入中的不可约性描述被分析的真实整数因子，绝不假设算法输出不可约。
- [VanHoeijGroupingComplete](../../lean/CLPoly/Refinement/VanHoeijGroupingComplete.lean)：
  先证明连续保留行索引使实际列比较器返回 true，再证明内循环对每个
  class 槽的精确更新公式。外循环保持“相等列有相同可选标号”的闭包
  不变量，故相等列不会被拆进不同输出组。包含播种、跳过已赋值项、
  数组更新、实际 collect 成员覆盖，未替换成抽象等价类算法。
- [VanHoeijGroupingTotal](../../lean/CLPoly/Refinement/VanHoeijGroupingTotal.lean)：
  依次证明比较器、内扫描、外扫描、收集和完整 extractCandidates 的
  无故障返回；所需条件只是实际数组形状。外扫描同时证明 classId 上界。
  与真实 LLL 的形状证明连接后，不再要求调用方假设分组成功运行。
- [VanHoeijActualSupportGrouping](../../lean/CLPoly/Refinement/VanHoeijActualSupportGrouping.lean)：
  合并上述链条，从真实 Hensel 支持、充分精度和真实 LLL 执行，产出
  实际 groups 数组与运行等式，并证明支持内任意两成员属于同一输出组。

## 验证与边界

七个模块均构建通过；最终实际支持分组目标为 3580 jobs，含依赖缓存。
先行 83 项分离公理审计已机器核对仅含标准公理；本轮合并审计
97/97 通过，3581 jobs。已对最终日志机器解析，全部仅依赖
`propext`、`Classical.choice`、`Quot.sound`。日志为
`/tmp/clpoly-vh-support-grouping-audit.log`；先行审计日志为
`/tmp/clpoly-vh-actual-support-audit.log`。

strict refinement boundary、verified refinements、`git diff --check` 通过。
本轮没有修改生产 C++、模板、B2B 注册或已有测试，因此未重复原生回归。

**仍未完成的义务**：

1. 组合定理的“保留行非空”仍是显式前提，需要由非零真实因子短向量
   的保留性消去；支持内任意两成员同组还需要和已有反向性质组合，
   得到完整真实因子支持与输出组对应。
2. 证明实际对称恢复、primitive、精确除法和完整 partition 验证在该
   分组上必定成功，不能只引用“若成功则正确”的旧定理。
3. 当前 D 仍可能随 m 增长；当前三次重试后回退也未改变。因此尚未
   证明生产控制器能达到上述阈值，有限回退不代替纯 van Hoeij 终止。
4. 生成对应关系最终审计和性能验收仍待完成；研究分支尚未合并。

## 度量

- 耗时：未单独记录净证明时间；最终组合模块单次编译约 14–16 秒，
  分组闭包模块 24 秒，无故障模块 31–34 秒，不能相加当作完整开发耗时。
- 迭代：本阶段共 15 次目标构建调用，含两次合并审计；其中分离模块
  6 次、83 项审计 1 次、分组闭包 4 次、无故障分组 2 次、组合 1 次、
  最终合并审计 1 次。另有一次在最终审计结束前读取日志，未得到结果，
  不作为审计通过或失败的证据。
- Lean：新增七个模块共 993 行、33 项定理和一个闭包谓词；审计新增
  33 项输出及两条 import。
- 对应 C++：本轮无原生逻辑变更；新覆盖已有列比较、内部赋值扫描、
  外部分组扫描和 class 收集的执行完备性。
- 问题与解决：依赖 ZMod 类型用重写整个目标搬运；Int 绝对值先显式
  转为 natAbs；数组 `getElem!` 使用明确索引界；match 先做受控化简再
  重写递归运行，避免把类型搬运或循环展开交给不受控 simp。
- 放弃的方案：不把输出组直接定义为数学等价类；不将运行成功或输出
  不可约性加入 oracle；不宣称当前随模数增长的截断参数已满足可达阈值。
