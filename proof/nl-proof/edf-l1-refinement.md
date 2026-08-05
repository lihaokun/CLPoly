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

定义 `edfZpSafe`：

1. 常数多项式返回空列表，与生成代码的基本分支一致。
2. 对满足 EDF 前置条件的非常数输入，返回 L2 已证明存在的首一不可约
   因子分解。这是 C++ 随机搜索成功时的规格级结果，不伪装成对 partial
   函数的计算等式。
3. 对不满足前置条件的单独入口，返回 `[g]`；Pipeline 的正确性定理
   只在 DDF 提供的 EDF 前置条件下调用已验证分支。

对应定理不再是错误的“输出等于 `edf ... []`”，而是：

```text
Monic g → Squarefree g → EqualDegree g d
  → EDFCorrect g d (edfZpSafe g d)
```

## 4. Pipeline 闭环

`Pipeline/L1.lean` 的 `edf_l1` 改为调用 `Refinement.edfZpSafe`，而不再直接在
Pipeline 内展开 `edf_correct_unconditional.choose`。这样 EDF 的总化边界和正确性都集中
在精化模块，Pipeline 只消费已证明接口。

## 5. 完成标准

- `Refinement/EDF.lean` 中无 `sorry` / `admit` / `native_decide`。
- 错误的 `__edf_Zp_ir_refines ... []` 命题被删除。
- Pipeline 实际调用 `edfZpSafe`。
- EDF 相关定理的 `#print axioms` 只含项目允许的基础公理。
- `lake build` 全量通过，并记录 devlog 和提交。
