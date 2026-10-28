# 严格 Hensel 入口首项系数边界

日期：2026-08-11

## 做了什么

- 严格翻译 `__hensel_lift` 构树前对 factor zero 的首项系数 baking：读取 `lc(f)`、模 `p` 转换、逐系数相乘、normalization 和数组写回。
- 严格翻译提取后对 result zero 的首一归一化：短路判断、`front()`、`ZZ::invert` 断言、逐系数缩放、模约化和数组写回。
- 为两个边界分别建立 raw→safe 执行不变量和精确语义 trace。

## 为什么做

这两个代码块是完整 C++ `__hensel_lift` 入口不可跳过的前处理与后处理。只证明树提升循环和叶子提取不足以构成完整入口精化。

## 关键决策

- 数组越界和模逆失败继续表现为 raw fault；证明通过对应的 C++ 前置条件排除错误，不改变程序定义。
- 最终归一化结果必须来自 `ZZ.invert`、`scaleCoeffs` 和 `__upoly_mod_coeff_raw_ir` 的真实执行等式，语义关系不携带预测输出。
- 不使用 fuel、`partial def`、oracle 或 L2 fallback。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
