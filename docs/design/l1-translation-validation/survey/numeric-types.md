# 基础数值类型使用扫描

## 全局频次

| 类型 | 总出现 | 作为变量 | 作为参数 | 作为返回类型 |
|---|---|---|---|---|
| `int32_t` | 3017 | 172 | 30 | 1 |
| `bool` | 1052 | 20 | 1 | 6 |
| `uint64_t` | 471 | 23 | 21 | 0 |
| `int64_t` | 352 | 25 | 5 | 0 |
| `size_t` | 322 | 33 | 0 | 0 |
| `char` | 69 | 0 | 0 | 0 |
| `double` | 30 | 2 | 0 | 0 |
| `uint32_t` | 25 | 0 | 0 | 0 |

## 运算符结果类型（基础数值类型）

| Operator | 结果类型 | 次数 |
|---|---|---|
| `BinaryOperator::<` | `bool` | 146 |
| `UnaryOperator::++` | `int32_t` | 106 |
| `UnaryOperator::!` | `bool` | 69 |
| `BinaryOperator::-` | `int32_t` | 50 |
| `BinaryOperator::==` | `bool` | 40 |
| `BinaryOperator::>` | `bool` | 37 |
| `BinaryOperator::&&` | `bool` | 34 |
| `BinaryOperator::+` | `int32_t` | 31 |
| `BinaryOperator::=` | `bool` | 29 |
| `BinaryOperator::>=` | `bool` | 24 |
| `UnaryOperator::++` | `size_t` | 22 |
| `BinaryOperator::<=` | `bool` | 19 |
| `BinaryOperator::=` | `int32_t` | 18 |
| `BinaryOperator::!=` | `bool` | 18 |
| `UnaryOperator::-` | `int32_t` | 12 |
| `UnaryOperator::--` | `int32_t` | 11 |
| `BinaryOperator::||` | `bool` | 9 |
| `BinaryOperator::=` | `int64_t` | 9 |
| `UnaryOperator::++` | `uint64_t` | 8 |
| `BinaryOperator::-` | `int64_t` | 6 |
| `BinaryOperator::/` | `int32_t` | 6 |
| `BinaryOperator::-` | `uint64_t` | 6 |
| `BinaryOperator::*` | `int32_t` | 6 |
| `UnaryOperator::++` | `int64_t` | 4 |
| `BinaryOperator::*` | `uint64_t` | 4 |
| `BinaryOperator::=` | `uint64_t` | 4 |
| `BinaryOperator::+` | `int64_t` | 3 |
| `BinaryOperator::*` | `double` | 3 |
| `BinaryOperator::%` | `int32_t` | 2 |
| `BinaryOperator::+` | `double` | 2 |
| `BinaryOperator::/` | `double` | 2 |
| `BinaryOperator::+` | `size_t` | 2 |
| `CompoundAssignOperator::*=` | `int32_t` | 2 |
| `CompoundAssignOperator::+=` | `int32_t` | 2 |
| `BinaryOperator::%` | `uint64_t` | 1 |
| `BinaryOperator::/` | `uint64_t` | 1 |
| `BinaryOperator::/` | `size_t` | 1 |
| `BinaryOperator::/` | `int64_t` | 1 |
| `CompoundAssignOperator::/=` | `int32_t` | 1 |
| `BinaryOperator::=` | `size_t` | 1 |
| `UnaryOperator::--` | `int64_t` | 1 |
| `CompoundAssignOperator::-=` | `int32_t` | 1 |

## 按类型聚合的运算符

### `bool` (合计 425)

- `BinaryOperator::<`: 146
- `UnaryOperator::!`: 69
- `BinaryOperator::==`: 40
- `BinaryOperator::>`: 37
- `BinaryOperator::&&`: 34
- `BinaryOperator::=`: 29
- `BinaryOperator::>=`: 24
- `BinaryOperator::<=`: 19
- `BinaryOperator::!=`: 18
- `BinaryOperator::||`: 9

### `int32_t` (合计 248)

- `UnaryOperator::++`: 106
- `BinaryOperator::-`: 50
- `BinaryOperator::+`: 31
- `BinaryOperator::=`: 18
- `UnaryOperator::-`: 12
- `UnaryOperator::--`: 11
- `BinaryOperator::/`: 6
- `BinaryOperator::*`: 6
- `BinaryOperator::%`: 2
- `CompoundAssignOperator::*=`: 2
- `CompoundAssignOperator::+=`: 2
- `CompoundAssignOperator::/=`: 1
- `CompoundAssignOperator::-=`: 1

### `size_t` (合计 26)

- `UnaryOperator::++`: 22
- `BinaryOperator::+`: 2
- `BinaryOperator::/`: 1
- `BinaryOperator::=`: 1

### `int64_t` (合计 24)

- `BinaryOperator::=`: 9
- `BinaryOperator::-`: 6
- `UnaryOperator::++`: 4
- `BinaryOperator::+`: 3
- `BinaryOperator::/`: 1
- `UnaryOperator::--`: 1

### `uint64_t` (合计 24)

- `UnaryOperator::++`: 8
- `BinaryOperator::-`: 6
- `BinaryOperator::*`: 4
- `BinaryOperator::=`: 4
- `BinaryOperator::%`: 1
- `BinaryOperator::/`: 1

### `double` (合计 7)

- `BinaryOperator::*`: 3
- `BinaryOperator::+`: 2
- `BinaryOperator::/`: 2
