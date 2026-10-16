# van Hoeij 正目标列循环的充分精度完成性

日期：2026-10-16。

## 改动与理由

此前已证明充分精度下实际分组及完整恢复验证成功，但尚未连接实际
列循环。本次先写[自然语言证明 §27–28](../../nl-proof/vanhoeij-eventual-separation.md)，
再实现三个模块：

- [VanHoeijFullWidthRound](../../lean/CLPoly/Refinement/VanHoeijFullWidthRound.lean)：
  从实际全宽轮次 trace 提取初始格、追加格、LLL 和分组；实例化已证
  统一阈值及恢复验证，得到该轮正常返回时 accepted=true。
- [VanHoeijColumnLoopComplete](../../lean/CLPoly/Refinement/VanHoeijColumnLoopComplete.lean)：
  接受轮次返回 complete=true；按实际剩余宽度下降归纳，证明正目标
  列循环在全宽不可拒绝条件下正常返回必完成。
- [VanHoeijColumnLoopSufficient](../../lean/CLPoly/Refinement/VanHoeijColumnLoopSufficient.lean)：
  以实际充分精度定理消去通用循环中的全宽接受条件，不留下成功 oracle。

新增定理均分析实际正常返回，不把执行完备性作为已证结论。
尚须处理零目标入口、内部 RawExec 错误不可达及外层精度增长可达性。

## 验证与问题

三个模块的目标构建通过；合并公理审计 158/158，仅 propext、
Classical.choice、Quot.sound，3625 jobs。全宽证明的首次编译遇到
record update 换行解析和 ZZ 别名算术处理问题，分别通过明确布局、
直接使用已有不等式修复。Int32 added 转换精确性由实际目标界给出，
不再添加无关的源次数 2^31 假设。

上一阶段运行结果已收齐：五组 B2B 1738/1738；新参数固定 FLINT
对照 18/18，证书完成且零回退。有后台负载的计时仅作诊断，见
[截断上界记录](2026-10-16-vanhoeij-bounded-divisor.md)。本阶段没有
修改 C++，没有无理由重复已通过的原生回归。

## 度量

- 耗时：未单独记录可靠活动工时，不将后台 B2B 等待当作编码时间。
- 迭代：3 次模块目标构建，首次失败后修复；另 1 次合并公理审计。
- Lean：新增 3 模块、236 行、4 项定理；审计追加 4 项。
- C++：本阶段没有修改，连接已有实际 columnRound/columnLoop 模型。
- 放弃方案：无。正目标和正常返回的条件明确保留，未扩大结论范围。
