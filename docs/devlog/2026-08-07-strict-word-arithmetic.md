# 严格位宽算术精化

## 日期

2026-08-07

## 做了什么

- 为生成的 C++ `_umul128` 建立首个机器字语义定理。
- 证明返回的高、低两个 `uint64_t` limb 精确重构输入的 128 位乘积。
- 证明直接展开 `BitVec 128` 的乘法、右移与截断，没有使用抽象乘法 oracle。
- 证明任意 `UInt128` 的高低 limb 商余重构定理。
- 证明生成的 `_add_carry3` 将三-limb 值增加 `b0 + 2^64*b1`，并精确按 `2^192` 回绕。
- 证明生成的 `__preinvert_limb` numerator 的精确 128 位自然数值。
- 在规范化除数最高位为一的条件下，证明生成预逆值满足 FLINT 公式 `B + pinv = (B^2-1)/pn`。
- 从该公式推导 Granlund--Möller 约简所需的严格双边界：`m*pn < B^2 ≤ (m+1)*pn`，其中 `m=B+pinv`。

## 为什么做

RawHeap `_poly_divrem` 的循环虽然已经终止且无越界，但其多项式不变量还依赖乘加和三-limb 取模的数学语义。`_umul128` 是该链的最底层乘法步骤。

## 关键决策

- 以 `2^64` 为 limb 基数，定理陈述直接观察生成函数的返回值。
- 利用两个 64 位输入的乘积小于 `2^128` 排除机器乘法回绕，再用商余等式重构。
- 将 `_add_carry3` 的两次 carry 分别视为 128 位中间和除以 `2^64` 的商；两条商余等式加权后消去 carry。
- 不引入 fuel、默认返回值或未证明的算术接口。

## 遇到的问题与解决方式

- 生成代码的移位量本身是 `UInt128`，不是 `Nat`；证明中保留该精确类型后再化简为自然数除法。
- `UInt64.ofNat` 和 `BitVec.ofNat` 都带模语义；分别用输入和乘积的位宽上界消去模运算。
- 预逆 numerator 的高低 limb 拼接使用 `Nat.shiftLeft_add_eq_or_of_lt`，显式证明低 64 位与高移位部分不重叠。
- 预逆 quotient 截断前先由 `pn ≥ B/2` 证明 quotient `< B`，排除低 limb 投影回绕。
- 双边界直接使用自然数 Euclidean division 的 `div_mul_le_self` 与 `div_lt_iff_lt_mul`，没有引用外部正确性定理。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictWordArithmetic.lean`

## 度量

- 耗时：约 20 分钟
- 迭代：6 轮 Lean 编译-修复
- Lean 新增/修改行数：约 265 行
- 对应 C++ 行数：约 15 行 `_umul128`、`_add_carry3` 与 `__preinvert_limb`
- 放弃的方案：直接位爆破 128 位乘法；通用量化命题应使用算术界与商余定理
