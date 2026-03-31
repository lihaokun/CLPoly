# 背靠背测试设计

## 1. 目标

验证翻译器的**翻译忠实性**：对于相同的输入，C++ 函数和翻译后的 Lean 函数产出相同的输出。

**不验证**：算法正确性（由 L2 精化证明保证）。

## 2. 测试范围

### 2.1 被测函数

`polynomial_factorize_zp.hh` 中的 13 个函数。

### 2.2 外部依赖处理

依赖的外部函数（`derivative`、`polynomial_GCD` 等）用**测试桩**替代：

| 外部函数 | 测试桩策略 | 理由 |
|---------|----------|------|
| `derivative` | 朴素求导（度数×系数，度数-1） | 简单正确 |
| `polynomial_GCD` | 返回第二个参数（假 GCD） | 只需类型匹配 + 可执行 |
| `pair_vec_div` | no-op（不修改输出） | 除法是外部原语 |
| `normalize` | identity（返回自身） | 简化 |
| `get_deg` | 取首项度数 | 简单正确 |
| `comp` | 返回 0 | 比较函数不影响逻辑 |

**关键**：C++ 端和 Lean 端必须使用**相同的桩**。否则比较无意义。

实现方式：
- C++ 端：编译一个特殊的测试程序，用桩替代真实函数
- Lean 端：在生成的 .lean 文件中用 `def` 替代 `opaque`

### 2.3 不测试的

- 依赖外部函数结果才能产生有意义输出的函数（如 `__squarefree_Zp` 的 while 循环依赖 `polynomial_GCD` 的正确性）
- 这些函数的背靠背测试只能验证**控制流翻译**是否正确（分支、循环结构），不能验证**最终结果**

## 3. 测试流程

```
1. 生成测试向量
   ├── 手动：从 CLPoly test suite 提取关键用例
   └── 随机：随机生成小多项式 + 小素数

2. C++ 执行
   ├── 编译 test_back2back.cpp（含桩函数）
   ├── 对每个函数 + 每个测试向量：打印输入→输出
   └── 输出为 JSON

3. Lean 执行
   ├── 生成 eval_back2back.lean（含桩函数 + #eval 调用）
   ├── lake env lean eval_back2back.lean
   └── 解析输出

4. 比较
   ├── 逐字段比较 C++ 和 Lean 输出
   ├── 一致 → PASS
   └── 不一致 → FAIL + 打印差异
```

## 4. 测试向量

### 4.1 格式

```json
{
  "function": "__make_zp",
  "inputs": {"val": 7, "p": 13},
  "expected_output": {"val": 7, "prime": 13}
}
```

### 4.2 生成策略

| 函数 | 测试向量策略 |
|------|------------|
| `__make_zp` | val ∈ {0, 1, p-1, p, 2p}，p ∈ {3, 5, 13, 17} |
| `__upoly_make_monic` | 小多项式（1-5 项），lc ∈ {1, 2, p-1} |
| `__extract_pth_root` | 度数全是 p 的倍数的小多项式 |
| `__upoly_subtract_x` | 小多项式 ± 含/不含 degree-1 项 |
| `__upoly_subtract_one` | 小多项式 ± 含/不含常数项 |

### 4.3 require 参数

Lean 函数有 `require` 参数（UB 证明目标）。对具体测试值，用 `by decide` / `by omega` / `by native_decide` 生成证明：

```lean
#eval __make_zp_ir 7 13  -- 无 require，直接调用
#eval __extract_pth_root_ir test_poly (by decide) (by decide)  -- 有 require
```

## 5. 比较规则

| 类型 | 比较方式 |
|------|---------|
| `UInt64` | 精确相等 |
| `Zp` | `val` 和 `prime` 分别比较 |
| `UMonomial` | `deg` 比较 |
| `SparsePolyZp` | 逐项比较（按度数降序） |
| `Array` | 逐元素比较 |
| `Bool` | 精确相等 |

## 6. 可行性

### 6.1 C++ 端

编写 `test_back2back.cpp`：
- include CLPoly 头文件
- 用 `#define` 或链接时替换外部函数为桩
- 对每个测试向量调用被测函数
- 以 JSON 输出结果

### 6.2 Lean 端

生成 `eval_back2back.lean`：
- include 翻译后的 IR 文件（桩版本，非 opaque）
- 对每个测试向量生成 `#eval` 调用
- require 参数用 `by decide` / `by omega`

### 6.3 限制

- `partial def` 可以 `#eval`，但如果函数实际不终止（死循环），`#eval` 会挂起。测试时加超时。
- 桩函数可能导致被测函数进入非预期路径（如 `polynomial_GCD` 返回错误结果导致 `__squarefree_Zp` 的 while 循环行为异常）。这不是翻译错误，是桩的限制。
- `Zp` 构造函数 `Zp.mk val prime` 在 Lean 中不检查 `val < prime`。测试向量应保证合法。

## 7. 预期产出

```
back2back test report:
  __make_zp: 8/8 PASS
  __extract_pth_root: 4/4 PASS
  __upoly_subtract_x: 6/6 PASS
  __upoly_subtract_one: 6/6 PASS
  __upoly_make_monic: 5/5 PASS
  total: 29/29 PASS
```

如果有 FAIL，说明翻译器有 bug——C++ 和 Lean 的行为不一致。

## 8. 实施步骤

1. 编写 C++ 测试桩 + 测试驱动程序
2. 生成 Lean 测试桩版 IR 文件
3. 编写测试向量（JSON）
4. 编写比较脚本 `back2back.py`
5. 运行 + 报告
