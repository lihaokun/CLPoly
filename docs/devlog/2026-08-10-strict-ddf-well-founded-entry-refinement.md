# Strict DDF 良基递归入口精化闭合

日期：2026-08-10

## 做了什么

- 将生成的 DDF 递归状态打包为携带完整不变式的 `DDFRawState`，并以
  `ddfRawStateMeasure` 实现良基递归，没有使用 fuel。
- 让生成循环逐次执行真实的 raw powmod、subtract-X、GCD、make-monic、
  exact-div 和 mod 调用，并通过成功执行等式驱动分支证明。
- 完成生成循环到 L2 `ddfLoop` 的全程语义模拟，返回终止状态、不变式和
  终止 guard 证据。
- 完成生成入口 `__ddf_Zp_raw_ir` 到 L2 `ddf` 的顶层精化，包括入口非空
  检查、初始状态、尾部 make-monic、结果 push 和 degree 机器字解码。
- 由生成器新增 `StrictDDFRefinement.lean`，自动声明公开定理
  `__ddf_Zp_raw_ir_refines_ddf`，原样保留 C++ 函数名中的双下划线。

## 为什么做

此前 DDF 已有各 raw 子调用的局部语义，但缺少生成良基循环和公开入口的
完整执行桥。没有这两层，无法声称 C++ 生成 L1 真正精化到 Lean L2。

## 关键决策

- `certifyRawExec` 只给实际 `RawExec` 成功值附带运行等式；证明在运行时擦除，
  不提供结果、不调用规格，也不改变错误行为。
- 最终状态继续携带 `DDFLoopInvariant`，使入口尾声可以从真实循环结果证明
  canonical、monic、degree bound，而不是重新假设这些性质。
- 长证明保留在 refinement 层；生成器产出精确命名的公开契约，并直接使用该
  已验证证明。这样契约名称与 C++ AST 符号稳定对应，同时避免生成 L1 与证明
  层形成循环导入。

## 遇到的问题与解决方式

- Lean 4.29 无法直接重写依赖 match 中的 raw 执行等式。将实际执行结果及其
  等式包装为 subtype，并在每层 match 后做定义归约，保持原控制流不变。
- 顶层 prime 重写依赖初始不变式证明。使用依赖安全的 `simp` 重写，并将循环
  初始状态显式命名后对齐生成入口。

## 验证

- `lake build CLPoly.Generated.StrictDDFRefinement`：通过。
- `python3 proof/cpp2lean_v2/tests/build_strict_ddf.py --check`：两个生成文件均可复现，
  且生成器占位检查通过。
- 禁用项扫描：生成 L1 中无 `sorry`、`admit`、`axiom`、`partial def`、fuel、
  规格 oracle 或 L2 算法回退。
- 公理审计：循环、内部入口证明和自动生成公开定理均只依赖
  `propext`、`Classical.choice`、`Quot.sound`。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_ddf.py`
- `proof/lean/CLPoly/Generated/StrictDDF.lean`
- `proof/lean/CLPoly/Generated/StrictDDFRefinement.lean`
- `proof/lean/CLPoly/Refinement/DDF.lean`
- `docs/devlog/2026-08-10-strict-ddf-well-founded-entry-refinement.md`
