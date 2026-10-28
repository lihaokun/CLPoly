# van Hoeij 纯证书重试与可表示精度下的成功终止

日期：2026-10-16。

## 问题与实现

三次失败后回退只保证控制器结束，不能证明纯 van Hoeij 必完成。
在入口充分精度成功已证明之后，本阶段按
[自然语言证明 §43–47](../../nl-proof/vanhoeij-eventual-separation.md)
删除生产 `__lll_factorize` 的三次预算与终端 Zassenhaus 分支。

未完成时调用新的 `__vanhoeij_next_precision`：通常取
max(full,max(1,2*precision))，超过 INT_MAX/2 时先尝试 INT_MAX，
仅在 INT_MAX 仍未获证时抛溢出异常。乘法只在可表示分支执行；
正常 Hensel 调用保证请求非负。profile 用 first 区分初次和后续调用。

严格生成模型执行相同分支，并以实际 (INT_MAX−precision).toNat
严格下降构造递归。RetryExecution 不再携带 Zassenhaus 操作；新增
precisionNonnegative 是同一次真实 Hensel 成功的已有定理，不是
可以选择输出或完成状态的回调。成功来源定理只剩完整 van Hoeij
证书，公共 FactorZZCorrect 已迁移，未保留回退正确性作为替代。

## 成功存在性及边界

- [VanHoeijRetryHensel](../../lean/CLPoly/Refinement/VanHoeijRetryHensel.lean)
  从源与选择条件、同一请求的既有 Hensel runtime readiness，构造
  实际 Hensel 返回及语义证书，不再假设调用已经成功。
- [VanHoeijRetryPrecision](../../lean/CLPoly/Refinement/VanHoeijRetryPrecision.lean)
  证明实际提升模数至少为正请求的 p^a，允许二次提升超过请求。
- [VanHoeijRetryTotal](../../lean/CLPoly/Refinement/VanHoeijRetryTotal.lean)
  在输入阈值 T<p^INT_MAX 的明确机器可行条件下，构造纯重试成功
  返回。若最大指数仍未完成，入口充分精度完成定理给出矛盾；否则
  正常增大指数并应用同一真实 rank 归纳。无输出不可约性假设、
  固定重试次数或终端回退。
- 公共 lllFactorize 的存在性包装进一步连接 FactorZZCorrect。
  目前仍明确给出实际启发式前缀返回及正起始指数；它们不是重组
  成功假设，但前缀的机器范围对应检查仍需单独完成。

核对发现原启发式完整指数的严格模板按 UInt32 范围检查，C++ 实际
使用 signed int 自增。下一阶段必须修正该范围及溢出语义，并核对
浮点启发式到 int 的转换，不能把当前包装条件隐去或宣告普遍等价。
手写模板的控制流核对和 B2B 仍不等于普遍 AST 翻译等价证明。

## 验证

- 新重试审计 25/25，仅 propext、Classical.choice、Quot.sound；3667 jobs。
- 公共 Pipeline.L1、B2B.Registry/Driver 组合构建 3610 jobs 通过。
- 新构建 profile 原生重组测试 390/390，包含 7 项精度边界检查。
- 原完整原生脚本 27 个测试文件全部通过。
- C++ B2B driver 已强制重建。新增 200 个独立整数 oracle 精度向量，
  覆盖半上限、饱和、最大指数溢出与固定种子的随机值；与已有 27 个
  重试、18 个公开入口向量合并运行，245/245 全部通过。此记录对应
  纯重试修改、启发式前缀范围修复之前的构建；前缀修改后须重跑。
- 严格证明边界、verified-refinements、模板/向量 freshness 和补丁
  格式检查通过。最终安静环境性能验收尚未进行，尚未合并。

## 度量

- 耗时：未单独记录可靠活动工时。
- 迭代：7 次显式 Lean 目标/组合构建，含生成模型、公共正确性、
  流水线、成功存在性和最终审计；B2B 自身另做 freshness 构建。
- Lean：新增 RetryTotal 161 行；RetryHensel 增加实际返回构造；
  重写 RetryTrace/Arithmetic/Precision，删除过渡版回退来源证明。
  这五个模块当前共 502 行，另同步生成模型、实例、审计与公共正确性。
- C++：修改外层控制器并新增 checked 精度 helper，新增 7 项原生断言；
  B2B 新增 200 条边界/随机向量及两端 dispatch。
- 放弃方案：不在半上限直接报错，先饱和到 INT_MAX 以覆盖可表示
  请求；不使用尝试次数作为纯算法成功证明；不保留回退来源分支。

涉及文件另包括严格重试生成脚本、两端 B2B 注册、精度向量生成脚本、
原生重组测试，以及修正方案和详细设计。完整证据日志存于
`/tmp/clpoly-vh-certificate-retry-*`。
