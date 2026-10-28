# 背靠背测试设计

## 1. 目标

验证翻译器的**翻译忠实性**：对于相同的输入，C++ 函数和翻译后的 Lean 函数产出相同的输出。

**不验证**：算法正确性（由 L2 精化证明保证）。

## 2. 测试范围

### 2.1 被测函数

`polynomial_factorize_zp.hh` 中的**全部 13 个函数**。

### 2.2 外部依赖处理

依赖的外部函数（`derivative`、`polynomial_GCD` 等）是**可信原语**——C++ 端使用真实 CLPoly 实现，Lean 端提供**正确的可执行实现**（非 opaque）。

| 外部函数 | C++ 端 | Lean 端 | 正确性保证 |
|---------|--------|---------|----------|
| `derivative` | 真实 CLPoly | Lean 朴素求导实现 | 数学定义直接实现 |
| `polynomial_GCD` | 真实 CLPoly | Lean 欧几里得 GCD | Mathlib `EuclideanDomain.gcd` |
| `pair_vec_div` | 真实 CLPoly | Lean 多项式除法 | 数学定义直接实现 |
| `normalize` | 真实 CLPoly | 除以首项系数 | 数学定义 |
| `get_deg` | 真实 CLPoly | 取首项度数 | 一行实现 |
| `comp` | 真实 CLPoly | 返回比较函数 | CLPoly 内部约定 |
| `inv` | 真实 CLPoly | 模逆（扩展欧几里得） | 数学定义 |
| `number` | 真实 CLPoly | `Zp.mk (a % p) p` | 数学定义 |

**关键原则**：
- C++ 端跑**真实代码**（可信基）
- Lean 端跑**翻译代码 + 正确的原语实现**（可信基）
- 两端的原语都是正确的 → 如果结果一致 → 翻译忠实
- 不使用假桩（arbitrary stub），所有实现基于数学定义

### 2.3 Lean 端原语实现来源

| 来源 | 适用 |
|------|------|
| Mathlib | `EuclideanDomain.gcd`（需要 SparsePolyZp → Polynomial 桥接） |
| 手写正确实现 | derivative, normalize, get_deg（简单，几行） |
| CLPoly 对应的数学定义 | pair_vec_div（多项式长除法） |

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
| `__squarefree_Zp` | 小无平方多项式（deg ≤ 10，p ∈ {3, 5, 7}） |
| `__ddf_Zp` | 无平方多项式（deg ≤ 15） |
| `__edf_Zp` | 等度无平方多项式（deg = k*d，d ∈ {1,2,3}） |
| `__factor_Zp` | 小多项式（deg ≤ 10）完整因式分解 |

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
- include 翻译后的 IR 文件（原语用正确实现替代 `opaque`）
- 对每个测试向量生成 `#eval` 调用
- require 参数用 `by decide` / `by omega` / `by native_decide`

### 6.3 限制

- `partial def` 可以 `#eval`，但如果函数实际不终止（死循环），`#eval` 会挂起。测试时加超时。
- Lean 端的原语实现是**独立正确实现**（非 CLPoly 翻译），与 C++ 端的实现**算法可能不同但结果相同**。如果两端原语实现有差异导致结果不同，需要排查是原语差异还是翻译错误。
- `Zp` 构造函数 `Zp.mk val prime` 在 Lean 中不检查 `val < prime`。测试向量应保证合法。

## 7. 预期产出

```
back2back test report:
  __make_zp:            8/8 PASS
  __upoly_make_monic:   5/5 PASS
  __extract_pth_root:   4/4 PASS
  __squarefree_Zp:      3/3 PASS
  __upoly_subtract_x:   6/6 PASS
  __upoly_subtract_one: 6/6 PASS
  __ddf_Zp:             4/4 PASS
  __edf_Zp:             3/3 PASS
  __factor_Zp:          3/3 PASS
  ...
  total: 50+/50+ PASS
```

FAIL 的含义：
- **翻译器 bug**：C++ 和 Lean 对同一输入产出不同结果
- **原语实现差异**：Lean 端原语与 CLPoly 实现行为不完全一致（需排查）
- **不是算法 bug**：算法正确性由 L2 证明保证

## 8. 实施步骤

1. 编写 C++ 测试桩 + 测试驱动程序
2. 生成 Lean 测试桩版 IR 文件
3. 编写测试向量（JSON）
4. 编写比较脚本 `back2back.py`
5. 运行 + 报告
