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
- 证明非零 `p` 的 `clz` 移位满足 `norm<64`，且 `p*2^norm ∈ [2^63,2^64)`。
- 将该自然数结论连接到实际 `UInt32 _norm` 与 `UInt64 (p << _norm)`，证明 C++ 左移没有回绕。
- 抽取 `_lll_mod_preinv` 两次使用的双-limb 近似商机器块，证明低 limb 比较精确表示 carry。
- 证明生成的 `(q1,q0)` 按 `B^2` 同余于 `u1*pinv + u0 + B*u1`，完整覆盖乘法、低位加法和高位带进位加法。
- 将生成的 `nmod_mul_ir` 精确分解为 UInt128 乘积拆分、原样预逆
  quotient/carry/correction round 和最终反归一化右移；该结构定理直接按
  生成分支证明，没有改写成 `% p`。
- 证明生成 round 中经 `Int32 → Int64 → UInt64` 的 carry 与此前已证明
  同余性质的 `preinvQuotientPair` 完全相同，从而把真实 IR 接到后续数值
  不变量所用的共享归约定义。
- 加强 quotient-pair 语义为逐 limb 精确公式：低 limb 是
  `(u1*pinv+u0) % B`，高 limb 是
  `((u1*pinv+u0)/B+u1) % B`。该定理显式证明低位加法 carry 与嵌套
  UInt64 加法的两层回绕。
- 从 `(B+pinv)*pn < B²` 和 `u1<pn` 推导商估计严格小于 `B`，并据此
  消去高 quotient limb 的 `% B`；生成机器商在有效输入域内不会回绕。
- 证明预逆估计倍数处于 `N-B < q*pn ≤ N+pn` 的 correction 窗口。
- 证明任意规范化 `pn`（`B≤2*pn`）下，生成 round 的最终条件减法输出
  严格小于 `pn`；该结论已从共享定义转回保留原始整数 cast 的真实 IR。
- 将 correction 控制流证明为独立的自然数定理：在估计倍数高于输入时，
  `B-(qd-N)` 经 add-back 回绕为 `d-(qd-N)`；估计倍数低于输入时，
  差值或差值加 `d` 经最后一次条件减法成为同一个 `N%d`。这精确保留
  C++ 的“先比较并可能加回、再至多减一次”顺序。
- 剩余局部义务已隔离为真实 `r>q0` 比较的两个 detector 性质，不能以
  `B` 在模 `d` 下为零来偷换（一般并不为零）。
- 证明负误差 detector：balance identity
  `B*k+A=d*(B-q0)` 强制 `q0<B-k`，因此机器比较必定识别减法回绕。
- 证明正误差 add-back 安全性：由
  `B*delta+d*(B-q0)=A` 与 `A+d*q0<B²` 推出 `delta+d<B`，排除加回
  模数时的额外 UInt64 回绕。
- 从 `X=u1*m+u0=B*t+q0` 推导统一且无自然数截断的 balance identity，
  并按 `q*d≤N` 与 `N<q*d` 分解为 detector 使用的正、负方程。
- 证明预逆 deficit `e=B²-m*d` 满足 `0<e≤d`；再用
  `u1<d`、`u0<B` 与正 balance 方程直接推出比较触发时
  `delta+d<B`，不再把该性质作为外部假设。

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
- `clz` 规范化使用 BitVec 已证明的首个置位上下界；机器整数转换均逐层证明保持自然数值。
- 双-limb 商估计没有抽象成数学除法；定理直接展开与生成代码一致的 UInt128/UInt64 操作并保留高位回绕。
- 生成器将最终右移复制到 correction 的各个 continuation 分支；证明通过
  穷尽这些真实分支对齐共享的 `preinvRoundIR`，而不是依赖大规模定义等价
  展开。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictWordArithmetic.lean`

## 度量

- 耗时：约 20 分钟
- 迭代：6 轮 Lean 编译-修复
- Lean 新增/修改行数：约 430 行
- 对应 C++ 行数：约 15 行 `_umul128`、`_add_carry3` 与 `__preinvert_limb`
- 放弃的方案：直接位爆破 128 位乘法；通用量化命题应使用算术界与商余定理
