# __binomial_ir_refines 证明完成

Date: 2026-05-26

## 做了什么

**Corpus.lean (生成代码改造):**
- 拆分巨型 `mutual` 块：在 binomial 函数前后加 `end`/`mutual` 截断
- 前置 `_loop___polynomial_to_zp_lex_0_ir` 和 `__polynomial_to_zp_lex_ir` 到第一个 mutual 块（因为被 `__assign_partial_zp_lex_ir` 依赖）
- `_loop___binomial_0_ir`: `partial def` → `def`，用 `Nat` 结构递归（`go : Nat → ...`）替代无法终止检查的 Int64 递归
- `__binomial_ir`: `partial def` → `def`（非递归，直接可改）

**ZZArith.lean (组合数学引理):**
- `my_mul_add`: Nat 乘法分配律 `a*(b+c) = a*b + a*c`
- `myChoose_one`: `C(n, 1) = n`
- `myChoose_eq_zero_of_lt`: `n < j → C(n, j) = 0`
- `myChoose_mul_eq`: `C(n, k+1)*(k+1) = C(n, k)*(n-k)`——核心恒等式，双重归纳证明
- `__binomial_ir_refines`: 完整证明，拆分为 Int64 条件分支 + 组合恒等式 + `calc` 链

## 关键决策

- 用 `Nat` 结构递归代替 Int64 上的 `termination_by`——因为 `omega` 不支持 Int64
- 组合恒等式用双重归纳（外 `k` 内 `n`）而非单层归纳
- 需要 `by_cases k+1 ≤ n` 分支处理 `C(n, k+1)=0` 时的 Nat 代数等式失效问题

## 涉及文件

- `proof/lean/CLPoly/Generated/Corpus.lean`
- `proof/lean/CLPoly/Refinement/ZZArith.lean`
