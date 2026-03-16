# CLPoly Lean 4 形式化验证

## 项目定位

CLPoly 因式分解模块的 Lean 4 机器检查证明。验证目标是 C++ 实现的数学正确性，不是重写算法。

依赖：Mathlib4 (stable)，Lean 4.28.0+

## 构建指令

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd proof/lean
lake build              # 编译全部
lake env lean FILE.lean # 单文件检查
```

## 验证架构：三层自顶向下

```
L3  数学基石    CLPoly/Math/         Mathlib 定理组合，纯数学定义与性质
L2  算法模型    CLPoly/Algorithm/    CLPoly 算法逻辑，抽象掉实现细节，证明满足 L3 Spec
L1  实现模型    CLPoly/Impl/         1:1 对应 C++（uint64 语义、数组越界、move），精化 L2
```

**工作顺序**：Spec（接口规约）→ Pipeline（顶层骨架）→ Math（L3）→ Algorithm（L2）→ Impl（L1）

详见 `proof/docs/implementation-roadmap.md`

## 证明工作流

对每个需要证明的定理或 sorry，必须先写**自然语言证明草稿**（关键步骤、依赖引理、Mathlib 路径），确认后再形式化为 Lean 代码。禁止跳过草稿直接写 tactic proof。

## 编码约定

### Lean 风格
- `autoImplicit = false`，所有变量显式声明
- 使用 `Nat` 而非 `ℕ`（除非在 Mathlib import 作用域内）
- `noncomputable` 按需添加（`EuclideanDomain`、`GCDMonoid` 等必须标注）
- 素数前置：`variable (p : ℕ) [Fact (Nat.Prime p)]`

### 文件组织
- `CLPoly/Spec.lean` — 子过程接口规约（Prop 定义）
- `CLPoly/Pipeline/` — 顶层正确性定理
- `CLPoly/Math/` — L3 纯数学定理
- `CLPoly/Algorithm/` — L2 算法模型（抽象算法逻辑）
- `CLPoly/Impl/` — L1 实现模型（1:1 对应 C++）
- `CLPoly/Experiment/` — Phase 0 实验（保留参考）

### 命名
- 规约 Prop：`SquarefreeDecomp`、`DDFCorrect`、`EDFCorrect`
- 数学定理：`irreducible_dvd_X_pow_sub_X`、`gcd_X_pow_sub_X_factors`
- L2 算法函数：`ddfLoop`、`squarefree_Zp`（对应 CLPoly 算法）
- L1 实现函数：`ddfLoop_impl`、`subtract_x_impl`（1:1 对应 C++ 函数名 + 控制流）

### L2 算法模型原则
算法模型捕获 CLPoly 的算法逻辑（DDF 循环、EDF 随机分裂、Hensel 提升），但抽象掉 C++ 实现细节（整数表示、数组布局、内存管理）。操作对象是 Mathlib 的 `Polynomial (ZMod p)` 等数学类型。

### L1 实现模型：1:1 对应原则
实现模型的控制流必须与 C++ 一一对应，不允许重写为"更优雅的 Lean 风格"。操作对象是 `U64`、`Vec` 等 C++ 语义模型。精化证明证明 L1 行为与 L2 算法一致。

| C++ 构造 | L2 算法模型 | L1 实现模型 |
|---------|-----------|-----------|
| `for` 循环 | 数学归纳 / 递归 | 顶层递归 + `termination_by`（1:1 控制流） |
| `while(true)` | 存在性论证 | `partial def` 或 fuel 参数 |
| `polynomial_GCD(a, b)` | `GCDMonoid.gcd a b` | 调用 Euclid GCD 的 L1 模型 |
| `(int64_t)(p - 1)` | `p - 1`（自然数） | `cast_u64_to_i64 (p - 1)`（显式溢出语义） |
| `vec[i]` | `list.get i` | `Vec.get i (proof : i < size)` |
| `std::move(v)` | 无对应（纯函数式） | `Ownership.move_from v (proof : ¬moved)` |

## Mathlib API 速查

| 需要 | 正确路径 |
|------|---------|
| ZMod p 是域 | `instance : Fact (Nat.Prime p)` 前置 |
| gcd | `GCDMonoid.gcd` 或 `EuclideanDomain.gcd` |
| 域上多项式除法 | `divByMonic` / `modByMonic`（`Ring` 足够） |
| 恒等式 | `modByMonic_add_div : p %ₘ q + q * (p /ₘ q) = p` |
| squarefree ↔ nodup | `UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors` |
| 有限域根 | `roots_X_pow_card_sub_X` |
| 域嵌入 | `nonempty_algHom_iff_finrank_dvd` |
| 层间投影 | `ZMod.castHom` (需 `dvd` 证明) |

## 已知陷阱

- `Polynomial.gcd` 不存在，用 `GCDMonoid.gcd`
- `EuclideanDomain` / `GCDMonoid` 是 noncomputable，不影响证明但阻碍 `#eval`
- `autoImplicit = false` 下 `ℕ` 需要 import Mathlib 才可用
- `let rec` 内的 `termination_by` 语法受限，改用顶层递归函数
- `modByMonic_add_div` 的顺序是 `f %ₘ g + g * (f /ₘ g) = f`，注意加法顺序

## 参考文档

- `proof/docs/implementation-roadmap.md` — 实施路线（自顶向下）
- `proof/docs/research-report.md` — Mathlib API 调研
- `proof/docs/mathlib-gap-report.md` — 缺口报告
- `proof/docs/phase0-experiment-plan.md` — Phase 0 实验设计
- `docs/design/factorization/formal-proof-ddf-edf.md` — DDF/EDF 半形式化证明
- `docs/design/factorization/formal-proof-univar-factorization.md` — Z[x] 半形式化证明
