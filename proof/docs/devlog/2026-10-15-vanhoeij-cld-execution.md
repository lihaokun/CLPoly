# van Hoeij 实际 CLD 循环与支持求和

日期：2026-10-15。持续目标仍 active，未合并。

## 本轮证明

新增 ModularCldValue 同时记录规范稀疏表示、模环中的对数导数语义、
对称代表元系数界。沿实际 cldPolysLoop 的运行等式归纳，证明输出长度
与全部元素均满足该性质。除法仍是实际模除法；不要求提升因子在整数
多项式环中整除原输入。这修正了旧循环证明不适用于一般 Hensel 因子的
前提，未修改或依赖旧错误适用范围的结论。

支持列表上的实际输出之和等于补因子乘所选模因子乘积的导数。
若整数因子 g 的模支持为常数倍的该乘积，由 f=h*g 推出该和等于
h*g' 的模像，允许非首一整数因子。再利用规范表示的 sparseCoeff
对应语义系数及 ZMod 整数投影判据，推导整数进位
sum(C_i)=coeff(h*g')+carry*m。进位由整除见证给出，不作为额外假设。

## 适用边界与后续

当前支持等式仍是明确的输入前提，尚未在新模块接入真实 Hensel 执行。
已定位 FactorZZ 中 henselCandidate_scaled_eq_divisor_mod_primePower
及候选存在性定理，可用于后续组合；不能把本轮结论称为端到端证明。
列表到无重复物理支持集合、实际截断格短向量、完整候选验证及顶层
控制器仍需连接。旧顶层严格模型尚未迁移，纯 van Hoeij 分离界未完成。

## 验证与度量

- 完整模块构建成功（3537 jobs，含已有依赖）。
- 首次系数求和归纳保留了不适用递归尾部的总和等式，Lean 拒绝；
  清除该无关前提后重新检查通过，未弱化目标。
- 本轮只新增证明及文档，未改生产 C++，未重复原生测试。
- 总耗时未精确记录；不估造。新增 Lean 模块 307 行。
- 五个关键定理的公理审计均仅依赖 propext、Classical.choice、
  Quot.sound，无 sorry 或新公理。git diff --check 通过。

## 文件

- [CLD 执行证明](../../lean/CLPoly/Refinement/VanHoeijCld.lean)
- [证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)
- [既有 Hensel 支持证明](../../lean/CLPoly/Refinement/FactorZZ.lean)
