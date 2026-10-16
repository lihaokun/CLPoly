# van Hoeij 修复最终验收

日期：2026-10-16。研究分支 `research/vanhoeij-lattice-certificate`，目标分支 `feature/l1-l2-refinement-codex`。

状态：修复及技术验收完成，可以合回原分支。最终 B2B 3107/3107，进程退出码 0；代码、证明及原生回归通过；三版本性能记录包含旧版超时、轻输入退化和主机负载限制，不表示所有输入均加速。

## 需求与证据

| 验收项 | 证据 | 状态 |
|---|---|---|
| 真因子数小于模因子数时保留完整证书 | 四个模因子恢复两个整数不可约因子；主路径不比较两者数量 | 原生通过 |
| 整数不可约但模上可分解 | `x^4+1` 的非平凡 CLD 证书 | 原生通过 |
| 非首一、精度不足、长行短组合及多列 | 专项重组回归及严格模型对应向量 | 原生及最终 B2B 通过 |
| 完整输出而非仅乘积验证 | 精确已知不可约因子匹配；整个分区除法且余项为 1；端到端不可约性定理 | 通过 |
| 纯 van Hoeij 成功终止 | 输入固定的充分精度阈值、实际入口成功、实际 Hensel 与提高精度循环；保留机器可表示性和既有 callee 契约 | Lean 通过 |
| 源码与执行模型对应 | 人工严格模板审查、生成一致性、调用隔离、机器范围证明、实际两端执行 | 静态检查及最终 B2B 通过 |
| 三版本 r=24/32 分阶段性能 | 同一 p=97、输入族、编译配置；每项计划三次；分别记录 Hensel、LLL、内外 Zassenhaus；FLINT 3.1.2 对照 | 已完成；15 组合完成、3 旧版组合超时；主机负载限制见报告 |

## 已完成检查

- 全量 `lake build`：3612 jobs，成功。
- 纯重试审计：34/34；公共流水线审计：340/340；分离及机器下标审计：212/212。逐名匹配审计声明，只使用 `propext`、`Classical.choice`、`Quot.sound`；其中公共审计两项不依赖任何公理。
- 原生回归：27 个测试文件通过、0 失败。最新 profiling 重组测试：406/406。
- 七个 van Hoeij 严格生成器的 `--check`、已验证细化注册检查、严格细化边界检查、B2B 物理运行时隔离检查通过。
- 三版本 18 个组合均已启动并取得终态，47 个已完成样本通过准确因子及 FLINT 对照。修复版 18/18 次通过，实际 LLL 调用非零、内外 Zassenhaus 均零。三项旧版超时保留截断记录；不足三次不填中位数。

性能详见[三版本分阶段报告](vanhoeij-three-version-performance.md)。测量期间本任务 B2B 已挂起，但系统服务和其他应用仍有明显负载；小幅时间差不能当作稳定加速保证。轻输入 r=24 观测中位数仍有约 11%–12% 退化，没有隐藏此结果，也没有将性能验收解释成所有输入必须更快。

最终 B2B 共 18 文件、3107 用例，全部通过且进程退出码 0。兼容入口 profiling 计数修正不改变不带 `CLPOLY_PROFILE` 的源码路径，故不使本批 B2B 或完整原生回归失效。

## 归档证据

[证据目录](vanhoeij-final-evidence/)包含原生、B2B、全量 Lean 构建、三组公理审计及三版本逐项计时原始日志。
[SHA-256 清单](vanhoeij-final-evidence/sha256.json)覆盖 32 个证据文件；
[基准上下文](vanhoeij-final-evidence/run-context.json)关联修复源码提交 `2673cbe3571d7cf0b5af21a4766a6bea1381ad86`，
[源码快照](vanhoeij-final-evidence/manifest.json)给出生产头文件与基准副本哈希。
归档前核对了所有组合终态、修复版每次的实际 LLL/回退调用、最终 B2B 汇总及 B2B 恢复记录。
五份较大的构建/审计/原生日志使用 `.log.gz` 无损保存，以保留工具输出的原始字节。
`static-checks.log` 记录了逐份压缩还原检查和未压缩内容的 SHA-256；其余日志可直接阅读。

## 范围说明

本修复不将旧模型成功、输出不可约性假设或固定次数后的 Zassenhaus 回退当作纯算法的证明。固定精度兼容入口仍有独立回退，纯主路径不调用它。完整前提和人工翻译边界见[源码对应审查](vanhoeij-source-correspondence.md)。

最终声明核对：

- [实际入口不可约性](../../proof/lean/CLPoly/Refinement/VanHoeijEntryIrreducibility.lean)中的 `hensel_tryRecombine_complete_irreducible` 从实际执行及完成标志推导完整乘积和各因子不可约性。
- [实际入口完整执行](../../proof/lean/CLPoly/Refinement/VanHoeijEntryTotal.lean)构造正常返回，并在输入阈值之上构造完整不可约分解；没有要求调用者预先提供成功结果。
- [纯重试总正确性](../../proof/lean/CLPoly/Refinement/VanHoeijRetryTotal.lean)中的 `concrete_lllFactorize_succeeds` 自行构造检查式启发式前缀，结论为实际执行存在且满足 `FactorZZCorrect`。保留的 Hensel 就绪契约涉及提升与归一化递归不变量；它不是 van Hoeij 输出不可约性假设。
- [公共入口](../../proof/lean/CLPoly/Refinement/VanHoeijFactorZZ.lean)与[流水线](../../proof/lean/CLPoly/Pipeline/L1.lean)使用新的纯重试路径。公共入口的成功运行正确性与纯重试的成功存在性分别审计，未以旧 StrictFactorZZ 控制器替代。

用户原始性能输入未提供；固定合成输入的对比不能直接声称复现了用户报告的约 50% 退化。历史诊断结果与最终三版本结果须分别标记，详见[参考报告](vanhoeij-flint-reference.md)。
