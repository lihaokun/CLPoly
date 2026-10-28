# `divrem` 状态传播与循环类型修复

## 日期

2026-08-07

## 做了什么

- 将 `word3 {lo, mid, hi}` 初始化生成为 `Word3.mk`，不再误生成为数组。
- 修正 `std::vector<word3>(n)` 为 `n` 个值初始化元素，而非错误的数组类型标注。
- 注册 `_umul128`、`_add_carry3` 和 `divrem` 的引用输出参数，使 limb 和工作数组更新进入 SSA 状态。
- 从版本化 SSA 定义恢复 loop phi 类型，避免三个同名 C++ 局部变量 `i` 相互污染。
- 添加升序 `Int64` 循环的机器加一差值递减引理，并为三个可由 guard 独立证明终止的 `divrem` 循环准备良基度量。

## 为什么做

此前探测代码虽然调用三字辅助函数，却丢弃了输出引用的修改；同时首个 `size_t i` 被后续 `int64_t i` 覆盖成错误类型。这两者都会产生可编译但不对应 C++ 状态变化的伪精化，必须在接入严格产物前消除。

## 自然语言证明草稿

- 初始化循环以 `r_len-i` 为度量，每轮 `i++` 且 guard 给出 `i<r_len`。
- 降序商循环以 `(i+1).toNat` 为度量，`i≥0` 时机器减一不溢出。
- 升序余式循环以 `(d-i).toNat` 为度量；`i<d` 且 `d` 可表示为 `Int64`，所以 `i+1` 不溢出并使差值减一。
- 内层 `j≤d` 在 `d=Int64.max` 时可能回绕，因此不能只靠 guard 总化；后续必须显式使用稠密表示大小不变量。

## 关键决策

- 引用输出通过返回 tuple 与 SSA 写回表达，不把 mutating helper 当成无副作用函数。
- phi 类型以 `(name, version)` 对应的真实定义为权威，裸变量名表只作为回退。
- 不为 `j≤d` 的边界状态加入 fuel、默认返回或截断执行。

## 问题与解决

Pass 6 的临时类型表按裸名称索引，后声明的同名 `i` 覆盖了早先 `size_t i`。SSA 版本本身没有冲突，因此 Pass 7 现在从每个 phi 的版本化 incoming definition 恢复唯一一致类型。

## 涉及文件

- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/constructor_map.py`
- `proof/cpp2lean_v2/passes/pass1_parse.py`
- `proof/cpp2lean_v2/passes/pass7_loop_lower.py`
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Generated/StrictGCD.lean`

## 度量

- 耗时：约 0.7 小时
- 迭代：5 轮 AST/HIR/MIR 探测与构建
- Lean 新增/修改行数：约 25 行（含递减引理）
- 对应 C++ 行数：`divrem` 约 65 行中的初始化及三个循环骨架
- 放弃的方案：未把 `d=Int64.max` 的回绕分支定义为默认结果；该状态必须由表示不变量排除。
