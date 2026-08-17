# `divrem` 的 `word3` L1 表示

## 日期

2026-08-07

## 做了什么

- 为 C++ `dense_upoly_zp::word3` 增加逐字段一致的 Lean L1 结构 `Word3`。
- 将 Clang 类型 `word3` 映射到该结构，并增加解析回归断言。
- 用 `CLPoly.Model`、现有 `StrictGCD` 构建和严格生成可复现检查验证改动。

## 为什么做

C++ `divrem` 用三个 `uint64_t` limb 保存延迟归约累加器。缺少该类型会使生成代码的工作数组元素类型退化为 `sorry`，无法进入真正的严格精化闭包。

## 关键决策

Lean 字段严格保持 C++ 的 `lo`、`mid`、`hi` 顺序和 `UInt64` 宽度，不把它提前抽象成数学整数；后续精化引理再证明三 limb 的数值解释。

## 问题与解决

无。地址别名和三字算术操作属于后续独立步骤，本步骤不为它们提供替代实现。

## 涉及文件

- `proof/cpp2lean_v2/passes/pass1_parse.py`
- `proof/cpp2lean_v2/tests/test_pass1_parse.py`
- `proof/lean/CLPoly/Model.lean`

## 度量

- 耗时：约 0.2 小时
- 迭代：2 轮探测与构建
- Lean 新增/修改行数：8 行
- 对应 C++ 行数：1 行结构定义
- 放弃的方案：未采用 `Int` 或 `UInt128` 单值替代，因为它们会丢失 C++ 三 limb 状态布局。
