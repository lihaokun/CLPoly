# Pass 1-6 全量烟测

- 目标 HIRs：**67**（factorize 展开 3 实例）
- OK: **67** / FAIL: **0**
- 平均 phi 节点数: **11.3** / 函数
- 平均 BasicBlock 数: **25.7** / 函数
- 总 phi: 758 / 总 BB: 1723

## Silent-bug 残余指标

- B2 phi.sources ver=0 非-param: **0**
- B6 aux_dropped (HIR.aux_lambdas != MIR.aux_defs): **0**
- B7 __sideeff_*=__write__(...) (root 非 def): **0**
- P0-1 ranged-for 漏写回原 container: **0**
- P0(lambda) 含 ref 的 lifted lambda 调用点未 destructure: **0**
- 其他 __sideeff_* (多为 discarded return): 0
