# Hensel 模系数真实精化

## 日期

2026-08-11

## 做了什么

- 将 C++ `__upoly_mod_coeff` 从抽象的 `HenselStepRawOps` 调用替换为严格
  生成函数 `__upoly_mod_coeff_raw_ir`。
- 生成函数逐项执行 GMP floor remainder 并剔除零系数，完全对应 C++ 循环。
- 证明公开内部桥 `__upoly_mod_coeff_raw_ir_refines`：真实 raw 调用必定成功，
  且输出解码到 `ZMod m[x]` 后与输入相同。
- `__hensel_step_raw_ir` 中六处模系数调用现在全部直接使用这一真实实现。

## 为什么做

`HenselStepRawOps` 最初保留了尚未证明的底层调用边界。本步骤移除其中一个
边界，使 Hensel 单步更接近完全具体的 C++ 执行，同时避免把模约化结果当成
规格 oracle 注入。

## 关键决策

- 使用 `ZZ.fdiv_r`，保持 GMP floor remainder 的源语义。
- 通过逐项 `ZMod` 同余证明零项压缩无语义影响，而不是假设输出数组等于
  输入数组。
- 失败行为仍然显式：本函数自身是总循环，因此精确返回 `.ok output`。

## 验证

- 严格生成文件 `--check` 同步。
- `lake build CLPoly.Generated.StrictHensel CLPoly.Refinement.Hensel` 成功。
- 严格 Hensel 文件无 `sorry`、oracle、fallback 或 fuel。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
