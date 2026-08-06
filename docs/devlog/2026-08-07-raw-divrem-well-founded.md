# Raw `_poly_divrem` 良基执行语义

## 日期

2026-08-07

## 做了什么

- 按 C++ `_poly_divrem` 的四个循环建立 `RawExec` L1：W3 初始化、内层乘加、倒序商循环、余式归约。
- 所有原始指针访问均经过 `RawHeap.read/write`，错误传播为 `RawFault`。
- 增加 raw limb slice 有效性谓词，以及“有效 slice 内读写必成功”的基础桥接引理。
- 证明单-limb 与三-limb写入保持所有 raw slice 的分配区长度不变量。
- 完成 W3 初始化循环的无故障定理：合法 A/W3 布局下必返回 `.ok`，且两项布局不变量保持。
- 完成内层乘加循环的无故障定理：由 B 长度和 `i+d` 工作区界推出每次 B/W3 访问合法，并保持布局不变量。
- 完成倒序商循环的无故障定理：逐项保持 Q/B/W3 三个布局不变量，并在非零商项分支调用已证的乘加循环。
- 完成余式归约循环的无故障定理：W3/R 边界由 `d≤lenW3` 推出，并保持两项布局。
- 完成 `normaliseU64` 与 `copyU64` 的无故障桥；复制前后 `SameLayout`。
- 完成整个 `_poly_divrem` 的终止/内存安全定理：按 C++ 规定的 A/B/Q/R/W3 容量，两个控制分支均必返回 `.ok`，且输出长度不超过缓冲区容量。
- 增加无默认值的 `readU64s` 与 `SliceRep` raw→safe 内容关系；证明合法 slice 存在唯一安全系数数组且数组长度精确等于 C++ length。
- 将安全系数数组按小端下标解释为现有 L2 使用的 `Polynomial (ZMod p)`，并证明每个合法 raw slice 唯一表示一个该类型的多项式。
- 增加严格准入检查器，将专用 lowering 绑定到稳定化后的完整 Clang AST 哈希，并检查原始指针签名、四循环形状、Lean artifact 哈希和禁用构造。

## 为什么做

旧 lowering 将 C++ 指针当作 Lean Array，越界时会产生默认值。新的执行语义保留地址、别名、工作区写入和失败路径，是证明 C++ L1 行为精化 L2 除法的必要基础。

## 自然语言证明草稿

初始化循环在 `i < lenA` 时读取 `A[i]` 并写入 `W3[i] = (A[i],0,0)`，度量 `lenA-i` 下降。倒序循环以 `ii` 为剩余商项数，处理 `i=ii-1`，写 `Q[i]`；非零商项进入 `j=0..d` 的乘加循环，其度量为 `d+1-j`。最后余式循环处理 `i=0..d-1`，度量 `d-i`。每次访问失败立即传播，因而不会得到伪结果。后续用 A/B/Q/R/W3 的 slice 长度、非别名和模数不变量证明每个访问均在界内，并归纳得到整个执行为 `.ok`。

## 关键决策

- 使用 Nat 良基/结构递归，不使用 fuel。
- 保留 `lenB=0` 的 C++ assert 失败路径为 `RawFault.assertionFailure`。
- 当前是 AST 绑定的专用 raw-state lowering；在通用 RawHeap 状态提升 pass 完成前，不把它描述为通用后端能力，也不宣称 raw→L2 精化已经完成。

## 遇到的问题与解决方式

- Clang AST 的声明地址每次进程不同：哈希前移除 `id`、`*Id` 和 `referencedMemberDecl`，同时保留其余完整结构。
- `word3 *` 访问跨度为三个 limb：复用已有 `readWord3/writeWord3`，保持 reinterpret 后的地址语义。

## 涉及文件

- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Generated/StrictDivrem.lean`
- `proof/lean/CLPoly/Impl/StrictDivremRefinement.lean`
- `proof/lean/CLPoly/Impl/RawPolynomialRep.lean`
- `proof/cpp2lean_v2/tests/check_strict_divrem.py`

## 度量

- 耗时：约 1 小时
- 迭代：6 轮 Lean/AST/checker 编译检查
- Lean 新增/修改行数：约 650 行
- 对应 C++ 行数：约 50 行 `_poly_divrem` 与 `_poly_normalise/memcpy` 依赖
- 放弃的方案：继续使用 Array `get!`/`set!`；原因是越界默认值不等于 C++ 语义
