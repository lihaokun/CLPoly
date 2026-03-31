# UB 证明目标目录

> 翻译器（ub_collector.py）检测的全部 UB 类型。
> 每种 UB 有 C++ 标准依据、检测规则、Lean 翻译。
> **不在此目录中的 = 不检测。**

---

## 检测的 UB

### UB-1: 整数除以零

| 项目 | 内容 |
|------|------|
| C++ 标准 | [expr.mul]/4: "If the second operand of / or % is zero, the behavior is undefined" |
| 触发 | `a / b`, `a % b`（signed 和 unsigned 均适用） |
| 检测规则 | `BinOp("/", lhs, rhs)` 或 `BinOp("%", lhs, rhs)` → 对 `rhs` 生成目标 |
| Lean 证明目标 | `rhs ≠ 0` |
| 示例 | `uint64_t q = a / p;` → `require hp : p ≠ 0` |

### UB-2: 数组越界访问

| 项目 | 内容 |
|------|------|
| C++ 标准 | [expr.sub]/1 + [dcl.ref]: 指针/引用超出对象范围为 UB |
| 触发 | `arr[i]`（包括 `vector::operator[]`，无边界检查） |
| 检测规则 | `ArrayAccess(arr, idx)` → 对 `idx` 和 `arr` 生成目标 |
| Lean 证明目标 | `idx < arr.size` |
| 示例 | `poly[i].coeff` → `require hi : i < poly.size` |

### UB-3: 空容器 front/back

| 项目 | 内容 |
|------|------|
| C++ 标准 | [sequence.reqmts]/2: "The expression a.front() is well-defined only if a is not empty" |
| 触发 | `v.front()`, `v.back()` 当容器为空时 |
| 检测规则 | `FieldAccess(obj, "front!")` 或 `FieldAccess(obj, "back!")` → 对 `obj` 生成目标 |
| Lean 证明目标 | `¬ obj.isEmpty` |
| 示例 | `f.front().second.prime()` → `require hf : ¬ f.isEmpty` |

### UB-4: 移位量越界

| 项目 | 内容 |
|------|------|
| C++ 标准 | [expr.shift]/1: "If the value of the right operand is negative or is greater than or equal to the width of the promoted left operand, the behavior is undefined" |
| 触发 | `a << n`, `a >> n`（signed 和 unsigned 均适用） |
| 检测规则 | `BinOp("<<", lhs, rhs)` 或 `BinOp(">>", lhs, rhs)` → 对 `rhs` 生成目标 |
| Lean 证明目标 | `rhs < 64`（假设操作数为 64 位；32 位时为 `rhs < 32`） |
| 示例 | `p << norm` → `require h_shift : norm < 64` |
| 注意 | 当前统一按 64 位处理。如果操作数是 32 位，需要检查 Clang AST 中的类型信息调整位宽。 |

### UB-6: 有符号整数溢出

| 项目 | 内容 |
|------|------|
| C++ 标准 | [expr]/4: "If during the evaluation of an expression, the result is not mathematically defined or not in the range of representable values for its type, the behavior is undefined" |
| 触发 | `a + b`, `a - b`, `a * b` 当操作数为 signed int 且结果超出 `[INT_MIN, INT_MAX]` |
| 检测规则 | `BinOp("+"/"-"/"*", lhs, rhs)` 且操作数含 `INT64` 类型 Cast → 生成目标 |
| Lean 证明目标 | `Int.noOverflow (lhs op rhs)` |
| 示例 | `int i = a + b;`（signed）→ `require h : Int.noOverflow (a + b)` |
| 注意 | CLPoly 因式分解代码极少用 signed 算术，预期此类目标数量为 0 或极少 |

### UB-7: unsigned → signed 类型转换溢出

| 项目 | 内容 |
|------|------|
| C++ 标准 | [conv.integral]/3: C++17 前为 UB，C++17 起为实现定义（通常为 mod 2^N 截断）。保守处理。 |
| 触发 | `(int64_t)x` 当 `x` 是 `uint64_t` 且值 > `INT64_MAX` |
| 检测规则 | `Cast(expr, INT64, UINT64)` 或 `Cast(expr, INT64, UINT128)` |
| Lean 证明目标 | `expr.val ≤ Int64.max` |
| 示例 | `int64_t y = (int64_t)big_uint;` → `require h : big_uint.val ≤ Int64.max` |
| 背景 | CLPoly DDF/EDF 的 64 位溢出 bug 就是此类问题（p^d 接近 2^64 时） |

### UB-8: assert 失败

| 项目 | 内容 |
|------|------|
| C++ 标准 | assert 宏在条件为 false 时调用 abort()。非标准意义的 UB，但我们要求证明条件成立。 |
| 触发 | `assert(cond)` |
| 检测规则 | 在 SSA 变换阶段（ssa_transform.py）检测，提升为函数签名的 `require` 参数。不在 ub_collector 中重复。 |
| Lean 证明目标 | `cond`（由 require 参数表示） |
| 示例 | `assert(p != 0)` → `(h_p_ne : p ≠ 0)` |

---

## 不检测的（及原因）

### 不检测-1: 无符号整数溢出/下溢

| 项目 | 内容 |
|------|------|
| C++ 标准 | [basic.fundamental]/4: "Unsigned integers shall obey the laws of arithmetic modulo 2^N" |
| 原因 | **不是 UB。** C++ 标准明确规定 unsigned 运算是 mod 2^N。Lean `UInt64` 有相同语义。 |
| 涉及运算 | `uint64_t` 的 `+`, `-`, `*` |

### ~~不检测-2: 有符号整数溢出~~ → 已改为检测（UB-6）

**已补充为 UB-6。**

### 不检测-3: 空指针解引用

| 项目 | 内容 |
|------|------|
| C++ 标准 | [expr.unary.op]/1 |
| 原因 | CLPoly 因式分解代码**不使用裸指针**（参数为值传递或引用）。翻译器已将 `.data()` 等指针返回方法映射为 no-op。 |

### 不检测-4: 悬垂引用/迭代器失效

| 项目 | 内容 |
|------|------|
| 原因 | 纯函数式翻译中不存在此问题。每次修改数组都生成新 SSA 版本。 |

### 不检测-5: 数据竞争

| 项目 | 内容 |
|------|------|
| 原因 | CLPoly 因式分解代码**无多线程**。 |

### 不检测-6: throw（显式异常）

| 项目 | 内容 |
|------|------|
| 处理方式 | 翻译为 `Except.error`，不生成证明目标。正确性定理条件化为 `.ok`。 |

---

## 覆盖性论证

CLPoly C++ 子集（blueprint §5a.1）的全部 UB 源：

| UB 源 | 检测？ | 编号 |
|-------|-------|------|
| 整数除以零 | ✅ | UB-1 |
| 数组越界 | ✅ | UB-2 |
| 空容器 front/back | ✅ | UB-3 |
| 移位越界 | ✅ | UB-4 |
| assert 失败 | ✅ | UB-5（SSA 阶段） |
| 有符号溢出 | ✅ | UB-6 |
| unsigned→signed 转换溢出 | ✅ | UB-7 |
| unsigned 溢出 | 不需要 | 不检测-1 |
| signed 溢出 | 暂不检测 | 不检测-2 |
| 空指针 | 不适用 | 不检测-3 |
| 悬垂引用 | 不适用 | 不检测-4 |
| 数据竞争 | 不适用 | 不检测-5 |
| throw | Except 处理 | 不检测-6 |

**完整性声明**：对于 CLPoly 因式分解代码子集（无裸指针、无多线程、无 signed 密集运算），上述 4 种检测 + 6 种不检测覆盖了 C++ 标准中的全部 UB 来源。

---

## 统计（polynomial_factorize_zp.hh）

自动提取结果（13 个函数）：

| UB 类型 | 数量 | 典型来源 |
|--------|------|---------|
| ARRAY_OOB | 8 | range-for 循环的 `__coll[__idx]!` |
| DIV_BY_ZERO | 1 | `extract_pth_root` 的 `deg / p` |
| SHIFT_OOB | 0 | （nmod_mul 不在此文件中） |
| **总计** | **9** | |
