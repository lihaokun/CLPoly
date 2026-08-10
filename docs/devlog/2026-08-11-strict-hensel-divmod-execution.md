# 严格 Hensel 模长除执行

## 日期

2026-08-11

## 做了什么

- 生成 C++ `__upoly_divmod_mod` 的真实 raw 语义。
- 将内层有序余式合并循环翻译为 `divmodMergeLoop`，下降量为两个输入后缀
  的剩余长度之和。
- 将外层长除循环翻译为精确分支 trace `DivmodTrace`，覆盖结束、首项模零
  删除、消去首项三个源分支。
- 增加 `DivmodTermination`，只提供有限执行 trace，不携带商余式结果。
- 将 Hensel 单步中的两次抽象 `divmodMod` 调用替换为严格生成入口。
- 证明 `__upoly_divmod_mod_raw_ir_terminates`：非空除数且真实 GMP 逆元调用
  成功时，生成执行确实返回具体商余式。

## 为什么做

模长除是 `__hensel_step` 中最后一个复杂内部循环。旧 Corpus 版本为
`partial def`，不能提供真实终止执行证据，也无法区分 assert 失败。本步骤
使 assert、GMP 逆元结果、每次首项消去和最终商余式均来自源执行。

## 关键决策

- 内层循环使用良基递归，不使用 fuel。
- 外层终止证据采用与 EDF 相同的有限精确 trace 思路；trace 的构造器只
  保存真实分支条件及下一状态，不保存或选择最终输出，因此不是结果 oracle。
- 空除数和逆元失败均返回 `RawFault.assertionFailure`，不回退到 L2 除法。

## 验证

- `build_strict_hensel.py --check` 验证生成同步。
- `lake build CLPoly.Generated.StrictHensel CLPoly.Refinement.Hensel` 成功。
- 严格链路无 `sorry`、oracle、fallback 或 fuel。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
