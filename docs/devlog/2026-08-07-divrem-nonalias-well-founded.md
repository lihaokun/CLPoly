# `divrem` 非别名路径良基闭合

## 日期

2026-08-07

## 做了什么

- 从 C++ `dense_upoly_zp::divrem` AST 生成明确命名的 `dense_upoly_zp_divrem_nonalias_ir`。
- 生成器精确验证并移除四项地址比较组成的 alias-protection 分支；不匹配源码形状时立即失败。
- 将四个 C++ 循环全部改为 `def + termination_by`。
- 对 inclusive `j <= d` 循环线程擦除前提 `d.toInt < INT64_MAX`，排除机器加一回绕。
- 严格产物通过可复现检查、Lean 构建与零占位审计。

## 为什么做

GCD 的 Euclid 分支调用 `divrem(q, r, a, b)`，四者是不同局部对象，因此执行的是非别名路径。直接删除地址判断而仍称一般 `divrem` 会掩盖 C++ 语义；显式专化并在 GCD 调用点验证声明身份，才能保持真实精化边界。

## 自然语言证明草稿

- `R3` 初始化循环：自然数度量 `r_len-i`。
- 内层乘加循环：度量 `(d-j+1).toNat`。guard 给出 `j≤d`；表示前提给出 `d<INT64_MAX`，因此 `j+1` 不回绕且度量减一。
- 商系数降序循环：度量 `(i+1).toNat`；guard `i≥0` 保证机器减一与整数减一一致。
- 余式归约循环：度量 `(d-i).toNat`；guard `i<d` 保证机器加一不溢出。
- alias guard 恰为 `&Q==&A ∨ &Q==&B ∨ &R==&A ∨ &R==&B`。非别名专化假设四项均假，故执行原函数的后续直线路径。

## 关键决策

- 产物命名带 `nonalias`，不对一般别名输入作未证明声明。
- `d<INT64_MAX` 是证明参数，不参与运行时分支；后续由稠密数组表示不变量提供。
- 不使用 fuel、`partial def`、默认值或 L2 多项式除法。

## 问题与解决

循环体的递归边最初位于生成的局部 continuation 中，遮蔽了外层 guard 证据。代码生成器按 continuation 的结构形状将其内联回 guard 分支，使递减证明直接使用当前 `i≥0`。SSA/BB 编号变化不会影响该结构检查。

## 涉及文件

- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/cpp2lean_v2/tests/build_strict_gcd.py`
- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Generated/StrictGCD.lean`

## 度量

- 耗时：约 0.8 小时
- 迭代：6 轮生成与 Lean 编译
- Lean 新增/修改行数：生成产物约 120 行，模型终止引理约 23 行
- 对应 C++ 行数：`divrem` 约 65 行
- 放弃的方案：未建立全局引用堆模型来覆盖一般 alias 调用，因为当前严格 GCD 唯一调用点可由局部声明身份证明非别名；未将 `d=INT64_MAX` 定义成默认结果，因为源 C++ 在该状态会回绕。
