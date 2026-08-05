# 2026-08-06: SQF L1→L2 精化闭合

## 做了什么

- 完成 `__squarefree_Zp_ir_refines`的两个主分支：导数为零的 Frobenius/p 次根分支，以及导数非零的 Yun 分解分支。
- 将稀疏表示的 GCD、精确除法、导数、`extract_pth_root`、`make_monic` 和两个输出循环连接到 `sqfZp`的数学语义。
- 在递归边界上传递 `Sorted` / `WellFormed` / `AllReduced` / 非零系数以及 UInt64/Int64 次数界。
- 将生成层外循环改为以 `squarefreeMeasure` 为度量的总函数；递归前增加防御式严格递减守卫，并在精化定理的合法输入上证明守卫恒真。原来两个 `decreasing_by sorry` 已清除。

## 自然语言证明摘要

1. 对输入多项式的次数做强归纳。稀疏表示的规范性保证首项次数等于多项式次数。
2. 导数为零时，所有非零项的次数均被 `p` 整除。`extract_pth_root` 精化为 `contract p`，其次数严格下降；归纳结果的重数再乘 `p`，与 `sqfZp` 的 Frobenius 分支一致。
3. 导数非零时，稀疏 GCD 精化为归一化多项式 GCD，精确除法得到 Yun 循环初值。循环不变式证明已产生的因子列表正确，剩余多项式为 monic、非零且导数为零。
4. 若剩余次数为正，对其 `contract p` 后递归；`natDegree(contract) * p = natDegree(residual)` 和 `p > 1` 给出严格递减。若剩余为常数，直接返回 Yun 已产生的列表。
5. `make_monic` 只乘以非零单位，而 `sqfZp` 对这种标量乘法不变，因此实现层归一化不改变分解语义。

## 验证边界

- SQF 链路（`Algorithm/SquarefreeZp.lean`、`Generated/Corpus.lean` 中 SQF 定义、`Refinement/SquarefreeZp.lean`、`Pipeline/L1.lean`）不再使用 `sorry`/`admit`。
- `lake build` 用于检查全库兼容性。仓库中 EDF、DDF、Hensel、Recombine 和部分 ZZ 算术的 L1 精化仍有独立 `sorry`；它们不在本次 SQF 闭合的声明范围内。

## 度量

- 耗时：~1 个持续开发会话（含调研、自然语言证明、形式化和调试）
- 迭代：本次收尾 12 轮编译-修复循环；前置工作按分段提交保留
- Lean 新增/修改行数：本次收尾约 354 行（不含本日志）
- 对应 C++ 行数：未单独统计；对应生成的 `__squarefree_Zp_ir` 及其 3 个循环
- 放弃的方案：直接以 `get_deg` 证明生成函数对所有数组终止；非规范输入上次数不具备所需单调性，改用防御式 `squarefreeMeasure` 守卫。
