# van Hoeij 零列宽分支的实际单元素验证完成性

日期：2026-10-16。

## 问题与证明路径

正目标列循环完成性不覆盖截断 CLD 扫描宽度为零的入口。实际入口
先验证单元素分组，失败后才计算截断；若宽度为零，列循环不会运行。
因此必须证明充分精度下这两个条件不能同时出现。

先编写并审核[自然语言证明 §29–31](../../nl-proof/vanhoeij-eventual-separation.md)，
随后实现四个模块：

- [VanHoeijZeroWidth](../../lean/CLPoly/Refinement/VanHoeijZeroWidth.lean)：
  零宽度使每个原始 CLD 系数绝对值小于 D；单位坐标向量的完整残差
  等于单个原始 CLD。将它送入已有结式分离定理，证明同一真实
  不可约因子中的两个局部模因子索引必须相等。
- [VanHoeijZeroWidthSupport](../../lean/CLPoly/Refinement/VanHoeijZeroWidthSupport.lean)：
  从同一次实际 Hensel 提升构造完整支持；排除空支持，并用实际合法
  索引数组的严格递增性排除两项以上支持，得到单元素完整模支持。
- [VanHoeijSingletonTargets](../../lean/CLPoly/Refinement/VanHoeijSingletonTargets.lean)：
  识别实际 singletonGroups 中每一组的解码乘积。以模支持唯一性证
  单射、模 p 素因子的乘积整除与互异性证满射，重排输入导出的因子
  族，构造与所有实际组对应的 ExactPartitionTargets。
- [VanHoeijZeroWidthValidation](../../lean/CLPoly/Refinement/VanHoeijZeroWidthValidation.lean)：
  从输入本原、正首项与好素数条件构造真实因子族，实例化上述支持
  与目标构造，连接实际恢复验证，得到实际返回 `(true, after)`。

此分支没有引入实际未执行的 LLL，也没有通过虚拟 extractCandidates
替代入口的 singletonGroups。数学因子族只用于证明，实际返回仍由
生产恢复验证对应的严格模型计算。

## 验证与边界

四个模块目标构建全部通过，最后组合目标为 3617 jobs。合并公理
审计 165/165 通过，3631 jobs，逐项解析确认仅依赖 propext、
Classical.choice、Quot.sound。新增模块无 sorry 或新公理声明，
补丁格式与新增文件行尾空白检查通过。严格边界与 verified-refinements
检查通过。本阶段仅新增证明与记录，没有重复运行已通过且未受修改
影响的原生或 B2B 测试。

最终组合定理已产生验证成功执行，但仍接收同一次 CLD/truncate 的
执行事实。下一步需要从真实 tryRecombine 入口 trace 自动取得这些
事实，联合正目标列循环定理，证明充分精度入口正常返回即完成。
内部 RawExec 错误不可达、外层精度指数到达门槛、机器溢出边界及
最终对应审计与性能验收仍须完成，尚未合并。

## 度量

- 耗时：未单独记录可靠活动工时，不把工具与构建等待误报为编码时间。
- 迭代：6 次模块目标构建；零宽度核心和 Hensel 支持各 2 次，目标
  构造及组合验证各 1 次。另 1 次合并公理审计。
- Lean：新增 4 模块、374 行、7 项定理；合并审计追加 7 项。
- C++：未修改；本阶段覆盖既有零宽度分支的数学与实际验证条件。
- 问题及解决：`positivity` 不接受一般 `1≤...`，先证非负项再线性
  推导；lambda 先化简再分支；不存在的数组/列表读取快捷引理改用
  已验证的有界读取与 toList 映射引理。
- 放弃方案：没有构造源程序未执行的零宽度 LLL，直接利用 CLD 残差
  和实际单元素分组，使证明保持实际控制流的边界。
