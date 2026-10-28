# DDF L1→L2 精化自然语言证明

日期：2026-08-06

## 1. 目标

对规范的稀疏 `Zp[x]` 输入 `f`，证明 C++ 翻译函数的输出与
L2 DDF 算法一致：

```lean
toPolyList (ddfZpSafe p f) p = ddf (SparsePolyZp.toPoly p f)
```

`ddfZpSafe` 保留 `__ddf_Zp_ir` 的控制流，但把内循环写成有良基度量的
总函数。生成语料库中原函数位于一个大型 `mutual partial` 块，Lean 不允许只将
DDF 一项改为 total，因此不对整个生成文件做无关重构。

最终再用控制流等式（对满足前置条件的输入）将 safe 结果与
`Generated.__ddf_Zp_ir` 对齐；若 Lean 无法对 `partial` 包装生成可用的归纳原理，
Pipeline 直接使用 safe 总函数，与 SQF 的常数保护包装同层级。

## 2. 需要的前置条件

现有骨架只有 `WellFormed` / `AllReduced` / `2*p ≤ 2^64`，不足以证明实现语义。
完整接口还需：

- `Sorted f`：`get_deg` 与首项才能代表 `natDegree`；
- `p*p ≤ 2^64`：稀疏多项式乘法/GCD/divmod 需要；
- 次数 `< 2^63` 与 `< 2^64`：`Int32/Int64/UInt64` 转换需要；
- `p * degree < 2^64` 及系数乘次数界：powmod/算术循环需要；
- 非零系数：稀疏表示规范性和首项语义需要。

Pipeline 的 `toSparsePolyZp` 可以从一个多项式的 `natDegree` 统一推导这些界，
但 `ddf_l1_correct` 需显式接受相应的机器整数安全假设。

## 3. 状态对应

循环状态 `(d, f★, h, result)` 与 L2
`ddfLoop H F d acc` 对应：

1. `toPoly f★ = F`；
2. `toPoly h = H`；
3. `toPolyList result = acc`；
4. `d.toNat = d_math`；
5. `f★` 和 `h` 保持 canonical representation 及所有机器整数界。

初始状态中 `h = X`、`f★ = f`、`d = 1`、`result = []`，上述五条直接成立。

## 4. 单步精化

在 `deg(F) ≥ 2*d` 时：

1. `__upoly_powmod_ir h p f★` 精化为 `(H ^ p) %ₘ F`。
   对二进制幂循环归纳，不变式为
   `result * base^e = H^p (mod F)`；每次乘法后用 divmod 余式桥接。
2. `__upoly_subtract_x_ir h' p` 精化为 `H' - X`。
   对降幂项列遍历：次数 1 项的系数减 1；若不存在则在正确位置插入
   `-1·X`；零系数项被删除。遍历不变式同时给出排序性。
3. `polynomial_GCD (h'-X) f★` 经 SQF 已完成的 GCD 桥接得
   `normalize(gcd (H'-X) F)`。
4. `get_deg gd > 0` 与 `0 < natDegree GD` 等价。

若 `GD` 次数为正：

- `make_monic` 不改变 `GD`，因为 GCD 输出已 monic；
- push `(gd,d)` 与 L2 `acc ++ [(GD,d)]` 一致；
- 精确除法给出 `Fnew = F /ₘ GD`，`normalization` 不改变 monic 商；
- `h mod Fnew` 精化为 `H' %ₘ Fnew`；
- 由 `GD ∣ F` 且 `deg GD > 0`，`deg Fnew < deg F`。

若 `GD` 次数为零，不修改 `F` 和 `acc`，只令 `d := d+1`。

## 5. 终止性

度量为：

```text
natDegree(F) + 1 - 2*d
```

- split 分支：`Fnew` 次数严格下降，`d` 增 1，度量严格下降；
- no-split 分支：`F` 不变、`d` 增 1，度量至少减 2。

对非规范数组和发生 UInt64 wraparound 的状态，safe 循环在递归前检查严格递减，
不成立时安全返回当前状态。在上述合法不变式下，证明该守卫恒真。

## 6. 终止分支

当 `deg(F) < 2*d`：

- `deg(F)=0` 时返回 `acc`；
- `deg(F)>0` 时追加 `(makeMonic F, deg F)`。不变式保证 `F` 已 monic，
  `get_deg` 的 UInt64 标签精确等于 `natDegree F`。

这与 `ddfLoop` 的终止分支完全一致。

## 7. 审核

1. **数学正确性**：分支、GCD、商、余式与 L2 定义逐项一致。
2. **无跳步**：powmod 和 subtract-x 各需独立循环不变式；不将它们当作未证黑盒。
3. **Lean 可形式化**：GCD/divmod/make-monic 桥接已在 SQF 精化中编译验证；
   需将可复用结果从 private 调整为公开定理。
4. **工程问题**：原 DDF 位于大型 partial mutual 块，不对整块做无关总化；
   采用精化层 safe total loop。
5. **边界情况**：空多项式、常数、`GD=1`、最后一个非常数剩余项、
   `d+1` 的 UInt64 wraparound 均由前置条件或防御守卫处理。

## 8. 实现顺序

1. 公开 SQF 中已证的 canonical/GCD/divmod 桥接定理。
2. 证 `subtract_x` 的列表语义和 canonical 保持。
3. 证 powmod 内循环及顶层语义。
4. 定义 DDF safe total loop，证明单步对应及度量递减。
5. 组装顶层精化、Pipeline 边界、公理审计和全构建。
