# EDF L1→L2 精化证明设计

> 状态：实现前审核
> 对应 C++：`clpoly/polynomial_factorize_zp.hh` 中 `__edf_Zp`

## 1. 现有骨架不能成立

`Generated.__edf_Zp_ir` 使用 RNG 反复生成候选多项式，找到非平凡 gcd 后递归分裂。
L2 `edf f d splits` 则显式消费 `splits`。旧骨架声称

```text
Generated EDF 的输出 = edf f d []
```

但当 `degree f > d` 时，`edf f d []` 按定义立即返回 `[f]`，而生成 EDF 可能成功
拆成多个因子。因此该命题一般为假，不能用填 tactic 的方式完成。

## 2. 终止性边界

C++ 包含 `while (true)`。当 RNG 不再产生有效候选时，该循环没有决定性终止
证明。当前 Lean RNG 是 xorshift64；种子 `0` 永远返回 `0` 并保持在 `0`，所以
生成的 partial EDF 在该输入上确实可发散。因而不声称对任意 RNG 的无条件
终止精化。

## 3. 可证明的总化接口

实现分为两层：

1. `Generated.edfZpTrySafe fuel` 是可执行的有界 L1 控制流：它总化了随机
   多项式生成、特征 2 trace、奇特征 powmod/GCD、非平凡分裂和两侧递归。
   每次尝试和递归都消耗 fuel；耗尽返回 `none`。失败尝试会保留已推进的
   RNG，不会重复同一候选。
2. `Refinement.edfZpSafe` 运行该 L1 路径。若成功结果通过 `EDFCorrect`
   认证则直接采用；若 fuel 耗尽或结果未通过认证，则返回 L2 已证明
   存在的首一不可约分解。
3. 常数多项式在 L1 路径中返回空列表，与生成代码的基本分支一致。

对应定理不再是错误的“输出等于 `edf ... []`”，而是：

```text
Monic g → Squarefree g → EqualDegree g d
  → EDFCorrect g d (edfZpSafe g sparse_g d fuel rng)
```

## 4. Pipeline 闭环

`Pipeline/L1.lean` 的 `edf_l1` 以 `toSparsePolyZp g`、固定种子 42 和
`8 * (degree g + 1)` fuel 调用 `Refinement.edfZpSafe`，而不再直接在
Pipeline 内展开 `edf_correct_unconditional.choose`。生成路径成功时使用其输出，
其余情况走显式的认证回退。

## 5. 完成标准

- `Refinement/EDF.lean` 中无 `sorry` / `admit` / `native_decide`。
- 错误的 `__edf_Zp_ir_refines ... []` 命题被删除。
- Pipeline 实际调用 `edfZpSafe`。
- EDF 相关定理的 `#print axioms` 只含项目允许的基础公理。
- `lake build` 全量通过，并记录 devlog 和提交。
