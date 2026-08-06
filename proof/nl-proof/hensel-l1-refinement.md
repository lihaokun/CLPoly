# Hensel L1→L2 精化证明设计

> 状态：实现前审核
> 对应 C++/Lean 生成入口：`Generated.__hensel_lift_upoly_ir`

## 1. 旧骨架的两个错误

旧 `hensel_l2` 定义成了 `List (Polynomial ℤ)` 的空列表，然后要求生成
函数的整系数输出等于它。这既没有表达 Hensel 规约，也在非空输入上为假。

生成入口返回：

```text
Array SparsePolyZZ × ZZ
```

即整系数提升因子和实际模数 `m`。L2 `HenselCorrect f k facs_p facs_pk`
则要求 `facs_pk : List (Polynomial (ZMod (p^k)))`。正确的 bridge 是将生成的
整系数因子 map 到 `ZMod (p^k)`，不是与某个整系数列表做字面等式。

## 2. 实际模数与 `p^k`

`a_target = k` 时，目标值是 `p^k - 1`。提升循环从 `m=p` 开始，每次执行

```text
m := m * m
```

直到 `m > p^k - 1`。所以返回的模数形如 `p^(2^r)`，其指数不小于 `k`，
但通常不等于 `k`。需证明：

```text
p^k ∣ m
```

然后用 `ZMod.castHom` 将模 `m` 的乘积不变式投影到模 `p^k`。

## 3. 实现分层

1. 总化 target 幂循环，证明 `a_target=k` 得到 `p^k-1`。
2. 总化因子首项的 leading-coefficient 调整循环。
3. 总化 Hensel tree 构建、单步提升、树递归和叶子提取，所有数组索引
   前置条件显式化。
4. 建立节点不变式：左右乘积等于目标（模当前 `m`），Bézout 关系保持，
   叶子顺序与输入因子一致。
5. 顶层证明输出数量与输入因子数量一致，且 map 到 `ZMod (p^k)` 后乘积
   还原 `f`。
6. Pipeline 使用已验证的总化入口；对生成 partial 函数不声称无前置条件
   的计算等式。

## 4. 顶层定理形状

最终定理应直接表达规约：

```text
facs_p ≠ []
→ map_p f = prod facs_p
→ Pairwise IsCoprime facs_p
→ HenselCorrect f k facs_p (henselZpSafe ...)
```

而不是错误的：

```text
generated_integer_factors = []
```

## 5. 完成标准

- `Generated/HenselSafe.lean` 覆盖所有可达 Hensel 循环和递归。
- `Refinement/Hensel.lean` 无 `sorry` / `admit` / `native_decide`。
- Pipeline 实际调用 Hensel safe 入口。
- 目标定理和端到端定理通过公理审计。
- `lake build` 全量通过，开发日志与提交完整。
