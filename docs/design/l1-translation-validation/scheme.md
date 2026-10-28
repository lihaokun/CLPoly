# L1 翻译验证方案：C++ → Lean IR → 精化证明

> 状态：方案设计 + 原型实验
> 目标：验证 CLPoly C++ 实现正确实现 L2 算法模型

---

## 1. 架构

```
C++ 源码                    Lean 形式化
─────────                   ──────────
polynomial_factorize_*.hh
    │                       L2 算法模型（已完成，0 sorry）
    │ [Clang AST 提取]            │
    ↓                             │
简化 IR                           │
    │ [翻译器]                    │
    ↓                             ↓
Lean IR def                 L2 def
    │                             │
    └──── 精化证明 ────────────────┘
          toPoly(impl_result) = l2_result
    │
    └──── UB-freedom 证明
          ∀ UB 点：前置条件成立
```

## 2. C++ 子集建模

CLPoly 使用的 C++ 受限子集：

### 2.1 类型

| C++ 类型 | Lean 模型 | UB 点 |
|----------|----------|-------|
| `uint64_t` | `Fin (2^64)` 或 `UInt64` | 溢出（乘法） |
| `int64_t` | `Int64` | 符号溢出 |
| `Zp` (Barrett) | `{ val : UInt64, p : UInt64, ninv : UInt64 }` | Barrett 溢出条件 |
| `ZZ` (dual-field) | `Sum Int64 (Ref MPZ)` | null 指针（大整数模式） |
| `vector<T>` | `Array T` | 越界访问 |
| `pair<A,B>` | `A × B` | 无 |
| `upolynomial_<Zp>` | `Array (Nat × UInt64)` | 不变量（降序、无零项） |

### 2.2 操作

| C++ 操作 | Lean IR | UB 证明目标 |
|----------|---------|------------|
| `a + b` (uint64) | `UInt64.add a b` | 无（wrapping） |
| `a * b` (uint64) | `UInt64.mul a b` | 无（wrapping） |
| `nmod_mul(a, b, p, ninv)` | `barrett_mul a b p ninv` | `a < p ∧ b < p ∧ p ≤ 2^63` |
| `vec[i]` | `arr.get ⟨i, proof⟩` | `i < arr.size` |
| `vec.push_back(x)` | `arr.push x` | 无（自动扩容） |
| `vec.size()` | `arr.size` | 无 |
| `a / b` (int) | `Int64.div a b` | `b ≠ 0` |

### 2.3 控制流

| C++ 构造 | Lean IR |
|----------|---------|
| `for (int i=0; i<n; i++)` | `Fin.foldl n f init` 或 `Nat.fold` |
| `while (cond)` | `partial def` + fuel 或 `termination_by` |
| `if (cond)` | `if cond then ... else ...` |
| `return` | 函数返回值 |
| 赋值 `x = e` | `let x := e` 或 `StateM` |

## 3. UB-Freedom 作为证明目标

### 3.1 原则

不建模 UB 行为。假设程序是 well-defined 的，然后在每个潜在 UB 点插入前置条件作为证明目标。

### 3.2 Barrett 模乘的 UB 分析

C++ `Zp::__nmod_mul`（number.hh ~line 200）：
```cpp
static uint64_t __nmod_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t ninv, uint32_t norm) {
    // Barrett reduction: compute a*b mod p without 128-bit multiply
    uint64_t hi, lo;
    __asm__("mulq %3" : "=d"(hi), "=a"(lo) : "a"(a), "rm"(b));  // [hi:lo] = a*b
    // ... Barrett reduction steps ...
}
```

UB 点：
1. `a*b` 的 128 位乘积本身无 UB（用 asm 获取 hi:lo）
2. Barrett 约化中的中间步骤可能溢出 — 需要 `p` 归一化条件

证明目标：`a < p ∧ b < p → result = a * b % p ∧ result < p`

### 3.3 多项式除法的 UB 分析

C++ `pair_vec_div`（basic.hh）：
```cpp
// f = q * g + r
pair_vec_div(q, r, f, g, comp);
```

UB 点：
1. `g` 非空（除数 ≠ 0）→ 证明目标：`!g.empty()`
2. `lc(g)` 可逆（首项系数非零）→ 证明目标：`g.front().second != 0`

## 4. 精化证明结构

### 4.1 表示函数

```lean
/-- 将 IR 多项式转换为 L2 数学多项式。-/
def toPoly : Array (Nat × UInt64) → Polynomial (ZMod p)
/-- 将 IR 整数转换为数学整数。-/
def toZZ : Int64OrMPZ → ℤ
```

### 4.2 精化定理

```lean
/-- __hensel_step 的 IR 实现精化 L2 模型。-/
theorem hensel_step_impl_refines
    (g h s t : Array (Nat × UInt64)) (m : UInt64)
    -- UB-freedom 前置条件
    (h_bounds : ∀ i, (g.get i).2.val < m.val)
    (h_m_pos : 0 < m.val) :
    -- 精化：impl 结果 map 到 L2 = L2 的 hensel_step
    toPoly (hensel_step_impl g h s t m) =
    hensel_step (toPoly g) (toPoly h) (toPoly s) (toPoly t) (m.val)
```

## 5. 原型实验

### 5.1 目标函数

选择 `__nmod_mul`（Barrett 模乘）：
- 最小（~15 行 C++）
- 有明确 UB 关注（溢出）
- 数学规约清晰（a * b mod p）
- 被全部因式分解代码使用

### 5.2 实验步骤

1. 手动将 `__nmod_mul` 翻译为 Lean IR
2. 定义 UB-freedom 证明目标
3. 证明 UB-freedom（Barrett 条件 → 无溢出）
4. 证明精化（IR 结果 = a * b mod p）

### 5.3 成功标准

- UB-freedom 证明：0 sorry
- 精化证明：0 sorry
- 证明可以自动/半自动生成
