# derivative_toPoly_eq 证明完成

## 日期
2026-05-30

## 做了什么
完成了 `derivative_toPoly_eq` 引理的证明：`SparsePolyZp.toPoly p (SparsePolyZp.derivative f) = Polynomial.derivative (SparsePolyZp.toPoly p f)`。

## 为什么做
这是 `__squarefree_Zp_ir_refines` 正确性证明所需的核心引理之一。在上一版本中，该引理有两处 `sorry` 待补。

## 关键决策与理由

### 1. 添加 `h_no_overflow` 条件
`SparsePolyZp.derivative` 在 UInt64 空间中做 `c.val * m.deg % p`，但 `c.val * m.deg` 可能溢出（UInt64 乘法 wrap around）。为避免这层复杂性，添加了预条件：
```
h_no_overflow : ∀ (m, c) ∈ f.toList, c.val.toNat * m.deg < 2 ^ 64
```
这个条件在 C++ 的实际计算中通常满足（degree 较小，coeff < p ≤ 2^63），但作为泛型引理必须显式要求。

### 2. 编写 UInt64→ZMod 辅助引理
```lean4
UInt64_mul_mod_toZMod_eq (a b m : UInt64) (p' : ℕ)
    (hm : m.toNat = p') (h_no_ov : a.toNat * b.toNat < 2 ^ 64) :
    ((a * b % m).toNat : ZMod p') = (a.toNat : ZMod p') * (b.toNat : ZMod p')
```
链式使用 `UInt64.toNat_mod`、`UInt64.toNat_mul`、`ZMod.natCast_mod` 和 `Nat.mod_eq_of_lt` 完成转换。

### 3. 结构
- 按 `f.toList` 归纳
- 逐项处理：deg=0 跳过、deg>0 且新系数=0 跳过（both sides zero）、deg>0 且新系数≠0 产生对应单项式
- UInt64 端的关键障碍：`(m.deg.toUInt64).toNat` 是三段式（Nat→UInt64→Nat）会截断 mod 2^64；利用 `h_no_overflow` 保证 `m.deg < 2^64` 从而 `(m.deg.toUInt64).toNat = m.deg`

## 涉及的文件
- `CLPoly/Refinement/SquarefreeZp.lean`：`derivative_toPoly_eq` 引理 + `UInt64_mul_mod_toZMod_eq` 辅助引理
