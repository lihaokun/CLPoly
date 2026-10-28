# 修正方案：L1 SQF 精化的 UInt64 溢出边界（乙类 admit）

状态：待用户确认（接口/架构层面变更）
涉及模块：`Pipeline/L1.lean`、`Refinement/SquarefreeZp.lean`、`Pipeline/FactorZp.lean`、
`Pipeline/FactorZZ.lean`、`Pipeline/FactorZpInstantiate.lean`、`Pipeline/FactorZZInstantiate.lean`

## 第一部分：复现与定位

三处相关 admit（“乙类”）：

1. `Pipeline/L1.lean:182` `h_no_overflow : ∀ x ∈ (toSparsePolyZp f).toList, x.2.val.toNat * x.1.deg < 2^64`
2. `Pipeline/L1.lean:188` `h_deg_bound : ∀ x ∈ (toSparsePolyZp f).toList, x.1.deg < 2^64`
3. `Refinement/SquarefreeZp.lean:2354`（`h_toPolyList_specific` 内）
   `(e * UInt64.ofNat p).toNat = e.toNat * p`，需 `e.toNat * p < 2^64`

调用链：
`factor_ZZ_cpp_correct` → `factor_ZZ_correct`（L2 骨架，`hfzp : ∀ g, g≠0 → FactorZpCorrect …`）
→ `factor_zp_l1_func` → `factor_Zp_l1` → `factor_Zp_correct`（L2 骨架，
`hsqf : ∀ g, g≠0 → SquarefreeDecomp g (sqf g)`）→ `sqfZp_l1_correct` →
`__squarefree_Zp_ir_refines`（其签名已把 (1)(2) 当**输入假设**）。

预期行为：C++ 翻译代码在“输入规模适配硬件”时数学正确。
实际：(1)(2) 在 `sqfZp_l1_correct` 处无法从 `toSparsePolyZp f` 无条件兑现；(3) 无法从
现有 in-scope 假设兑现。

## 第二部分：根因分析

**根因（非症状）**：C++ 用 `UInt64` 存指数/度数，对度数 `n` 与系数 `v`（`v < p`）作
`v * n`、指数 `e * p` 等乘法。这些在 `n`、`e` 足够大时**真的会溢出**——这是 C++ 实现的
硬件前置条件，不是证明技巧问题。数学上（`Polynomial (ZMod p)`）指数/度数无 2^64 上界，
所以“无条件正确”本身为假。正确的形式化陈述必须携带前置条件。

两个**相互独立**的子问题：

- **(甲) L1 边界兑现 (1)(2)**：可从一个顶层度数前置条件推出。
  `x.1.deg = n ≤ f.natDegree`，`x.2.val.toNat = ZMod.val (f.coeff n) < p`。
  取前置 `f.natDegree * p < 2^64` 即可：
  - (2) `n ≤ natDegree ≤ natDegree*p < 2^64`（`p ≥ 1`）。
  - (1) `x.2.val.toNat * n ≤ (p-1)*natDegree ≤ natDegree*p < 2^64`。
- **(乙) 内部输出指数界 (3)**：`e` 是 `__squarefree_Zp_ir_safe p g_2` 的**输出乘数**，
  被 ×p。需 `e * p < 2^64`。这**不能**由 (1)(2) 得出——(1)(2) 约束的是输入的
  coeff×deg 与 deg，(3) 约束的是**输出乘数**。需要一条新的结构不变量：
  “`__squarefree_Zp_ir_safe` 的输出乘数 ≤ natDegree(输入)”，再配合度数前置条件。
  该不变量需随强归纳一起传播（每层 p 次根把乘数 ×p，故实际需 `natDegree * p^depth` 级别的
  界，或改为对“输出乘数”单独立不变量）。

## 第三部分：参考实现对照

参考 FLINT `nmod_poly_factor` / Maple `Squarefree`：工业实现同样假设指数/度数落在机器字长内
（FLINT 用 `slong` 存指数，隐含 `< 2^63`）。即“输入适配硬件”是 CAS 的通行前提，不是我们独有的
简化。因此**加前置条件**与参考实现的隐含假设一致，非 workaround。

## 第四部分：修正方案

### 方案 A（推荐，对应用户已选“顶层加度数前置条件”）

1. **`__squarefree_Zp_ir_refines`（已有 (1)(2) 假设）**：新增一条输出指数不变量或前置
   `h_exp_bound`，用于消解 (3)。具体：在强归纳中额外维护
   “结果中每个乘数 e 满足 e * p < 2^64”，base/step 由度数前置推出。
2. **`sqfZp_l1_correct`**：签名新增 `hdeg : f.natDegree * p < 2^64`（或更强的 `* p^?`），
   用它兑现 (1)(2)（如第二部分 (甲)）。
3. **接口弱化（关键，触及 0-sorry L2 核心）**：
   - `factor_Zp_correct`：`hsqf : ∀ g, g≠0 → …` ⇒ `hsqf : SquarefreeDecomp f (sqf f)`
     （体内仅用 `hsqf f hf`，弱化安全）。
   - `factor_ZZ_correct`：`hfzp : ∀ g, g≠0 → …` ⇒ `hfzp : FactorZpCorrect fp …`
     （体内仅用 `hfzp fp hfp_ne`）。
   - 修 3 处 `factor_Zp_correct` 调用（`FactorZpInstantiate` ×2、`L1`）与
     `factor_ZZ_correct` 调用（`FactorZZInstantiate` ×2、`L1`），改为在 f/fp 处提供。
   - `factor_ZZ_cpp_correct` / `factor_zp_l1_func_correct` / `factor_Zp_l1`：加 `hdeg` 前置。

   影响面：5 文件。弱化假设对定理**真值**恒安全（仅调用点需在 f/fp 处提供），风险在
   编译期由 build 兜底；但因触及 0-sorry L2 核心，需先确认。

### 方案 B（不改接口，条件回退包装）

`sqfZp_l1 f := if f.natDegree * p < 2^64 then <C++ 路径> else sqfZp f`，
使 `sqfZp_l1_correct` 无条件成立，`factor_Zp_correct` 接口不动。
代价：`factor_ZZ_cpp_correct` 验证的是“C++＋L2 混合函数”，稀释“验证纯 C++”目标。
与用户“加前置条件”意向不符，仅列作对照。

### 建议的执行时序

(3) 的输出指数不变量与 Branch B（`2371`，导数≠0 全链）**同处 `__squarefree_Zp_ir_refines`
证明体内**，且 Branch B 完成后该定理的归纳骨架才完整、更易一并加输出指数不变量。
故建议：**先做 Branch B（自包含，不动 L2 核心）→ 再统一做方案 A 的接口弱化 + (甲)(乙) 兑现**。
