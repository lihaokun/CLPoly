# CLPoly Lean 4 形式化验证

## 项目定位

CLPoly 因式分解模块的 Lean 4 机器检查证明。验证目标是 C++ 实现的数学正确性，不是重写算法。

依赖：Mathlib4 (stable)，Lean 4.28.0+

## 构建指令

```bash
export PATH="$HOME/.elan/bin:$PATH"
lake build              # 编译全部
lake env lean FILE.lean # 单文件检查
```

## 验证架构：三层自顶向下

```
L3  数学基石    CLPoly/Math/         Mathlib 定理组合，纯数学
L2  算法模型    CLPoly/Algorithm/    1:1 对应 C++，证明满足 Spec
L1  表示层      CLPoly/Repr/         uint64 溢出、数组越界（可选）
```

**工作顺序**：Spec（接口规约）→ Pipeline（顶层骨架）→ Math（L3）→ Algorithm（L2）→ Repr（L1）

详见 `proof/docs/implementation-roadmap.md`

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
- `CLPoly/Algorithm/` — L2 算法模型（1:1 对应 C++）
- `CLPoly/Experiment/` — Phase 0 实验（保留参考）

### 命名
- 规约 Prop：`SquarefreeDecomp`、`DDFCorrect`、`EDFCorrect`
- 数学定理：`irreducible_dvd_X_pow_sub_X`、`gcd_X_pow_sub_X_factors`
- 算法函数：`ddfLoop`、`squarefree_Zp`（对应 C++ 函数名）

### 1:1 对应原则
Lean 算法函数的控制流必须与 C++ 一一对应，不允许重写为"更优雅的 Lean 风格"。

| C++ 构造 | Lean 4 对应 |
|---------|------------|
| `for` 循环 | 顶层递归 + `termination_by` |
| `while(true)` | `partial def` 或 fuel 参数 |
| `gcd(a, b)` | `GCDMonoid.gcd a b`（非 `Polynomial.gcd`） |
| `f /ₘ g` | `divByMonic`（需 `Monic g` 或 `Ring` 即可） |

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
