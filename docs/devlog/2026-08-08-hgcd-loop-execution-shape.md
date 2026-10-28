# HGCD 良基循环执行分解

日期：2026-08-08

## 做了什么

- 证明 `hgcdIterLoop_stop`：C++ `_hgcd_iter` 的循环条件为假时，生成执行原样返回当前状态。
- 证明 `hgcdIterLoop_step_shape`：每次成功的非终止迭代都来自一次真实 `_poly_divrem`、两次真实 `_mat_row_update`、源码中的指针轮换和对更短余式的递归调用。
- 将两个执行分解定理加入 HGCD 源码门禁。

## 为什么做

SQF 所依赖的 GCD/HGCD 必须从生成 C++ L1 执行逐步精化到 Lean L2。循环语义不能由规格调用、fuel 或假定的结果替代，因此需要先把良基循环的实际执行路径完整暴露给后续不变式证明。

## 关键决策

- 保留 `hgcdIterLoop` 中绑定 divrem 方程的依赖匹配；该方程正是余式长度严格下降的终止证据。
- 在证明层用 `split at hrun` 分解依赖匹配，不改写运行定义，也不增加任何额外分支。
- 定理结论保留两次 row update 的实际执行等式和递归尾调用等式，供后续逐项应用 divrem 与矩阵行更新语义定理。

## 遇到的问题与解决方式

普通 `match` 会让 Lean 的终止性目标失去实际 divrem 方程；直接用外部等式化简依赖 `match` 又无法消除内部证明绑定。最终直接分裂已展开执行中的依赖匹配，同时保留其返回值和方程，完成了无语义改动的执行分解。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-08-hgcd-loop-execution-shape.md`
