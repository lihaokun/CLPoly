# Phase 4 Recombination：nl-proof + 框架形式化

> 日期：2026-03-23
> 分支：`feature/formal-proofs`

---

## 做了什么

1. **Recombine nl-proof**（v3.1，5 轮审核通过）
   - Hensel 唯一性定理（B 首一条件精确化）
   - 因子恢复定理（lc-baking 两种情况，匹配 C++）
   - RecombineCorrect 组合

2. **Mignotte bound nl-proof**（v2，精确匹配 C++）
   - Mathlib 已有 Mahler measure 完整 API
   - 精确界 C(n,n/2)·‖f‖₂ = C++ 的 `__mignotte_bound`

3. **Lean 形式化框架**（Recombine.lean）
   - `recombine_correct` 0 sorry（Z[x] UFD 构造）
   - `mignotte_bound` 1 sorry（接口声明，Mathlib API 确认）
   - `hensel_unique` 1 sorry（接口声明）

4. **Hensel.lean 完成**（366 行 0 sorry）
   - `hensel_step` + `hensel_two_factor` 全部编译通过

## 为什么做

ZZ 因式分解管线的最后组件。Recombination 将 Hensel 因子恢复为真正的 Z[x] 不可约因子。

## 教训

### 教训 1：不要过早简化，先对齐 C++

**问题**：多次尝试用"更简单的数学"替代 C++ 的实际算法：
- 方案 A（UFD 存在性）：纯数学正确但不验证算法 → 被用户否决
- Cauchy 弱界替代 Mignotte：不能验证 C++ 的精度选择 → 被用户否决
- M(f) ≤ ‖f‖₁ 替代 M(f) ≤ ‖f‖₂：多了 √(n+1) 因子 → 审核发现

**正确做法**：
1. 先读 C++ 代码，确定 C++ 具体用了什么（`__mignotte_bound` 用 `C(n,n/2)·‖f‖₂`）
2. nl-proof 必须精确对齐 C++（同一个界，同一个精度条件）
3. 不要因为"更容易证"就换一个界——如果换了就不能验证 C++ 的正确性

**规则**：L2 证明的目标是验证 C++ 算法，不是找最简数学证明。CLAUDE.md 已更新此原则。

### 教训 2：Hensel 唯一性需要 "B 首一 in Z_m[x]"，不仅 "mod p 首一"

**问题**：第一版 nl-proof 只要求 Ā, B̄ 首一（mod p），但这不够——唯一性需要 B 在整个 Z_m[x] 中首一。

**根因**：在归纳步骤中，从 Ē = c·Ā、F̄ = -c·B̄ 推出 c = 0 需要 B₁ 的 leading coefficient 精确等于 1（不只是 ≡ 1 mod p）。仅 "mod p 首一" 允许 c ≠ 0。

**发现过程**：审核时尝试验证 deg(E) < deg(Ā) 的推导，发现 deg(E) = deg(Ā) 的情况没被排除。追查 GCL Theorem 6.1 确认需要 B_m 首一。

**教训**：经典定理的精确条件很重要。"首一"在不同环中（F_p vs Z_m）含义不同。

### 教训 3：nl-proof 中不要留探索/试错文本

**问题**：nl-proof 多次包含"等等——这不对"、"让我重新推导"、"hmm"等探索性文本，导致审核困难。

**正确做法**：写完探索后，清理为干净的证明链。每次审核修正后删除废弃路径。

### 教训 4：先查 Mathlib 再声明 sorry

**问题**：最初认为 Mignotte bound "需要复分析，很难"，声明为 sorry。后来发现 Mathlib 已有完整的 Mahler measure + Jensen 公式。

**正确做法**：在声明 sorry 之前，先搜索 Mathlib 是否已有所需引理。CLPoly 依赖的 Mathlib 版本已包含 `Analysis.Polynomial.MahlerMeasure` 和 `Analysis.Complex.JensenFormula`。

### 教训 5：审核标准要精确（√(n+1) 问题）

**问题**：nl-proof 用 M(f) ≤ ‖f‖₁ ≤ √(n+1)·‖f‖₂，比 C++ 宽了 √(n+1) 因子。第一次审核没发现（因为界"看起来差不多"）。

**正确做法**：审核时必须验证精确数值匹配：
- C++ 用 `B = C(n,n/2) · ‖f‖₂`
- 我们的界必须 ≤ B（不能有额外因子）
- √(n+1) > 1 导致 C++ 的 m 可能不满足我们的精度条件

## 涉及的文件

| 文件 | 状态 |
|------|------|
| `proof/nl-proof/phase4-recombine.md` | nl-proof v3.1（5 轮审核） |
| `proof/nl-proof/phase4-mignotte.md` | nl-proof v2（Mathlib API） |
| `proof/lean/CLPoly/Algorithm/Recombine.lean` | 框架，2 sorry |
| `proof/lean/CLPoly/Algorithm/Hensel.lean` | 完成，0 sorry |
| `proof/CLAUDE.md` | 添加 L2 原则 + 审核标准 |

## 度量

- 耗时：~10 小时（nl-proof 5 轮审核 ~6h + Hensel 形式化 ~4h）
- nl-proof 迭代：Recombine 5 版本 + Mignotte 2 版本
- Lean 新增行数：Hensel 366 行 + Recombine 94 行 = 460 行
- 放弃的方案：
  - UFD 纯存在性（不验证算法）
  - Cauchy 弱界（不匹配 C++ 精度）
  - M(f) ≤ ‖f‖₁ 路径（多 √(n+1) 因子）
  - Ā 首一条件（过度约束，且不匹配 Case 1 的非首一 Ā）

## 当前全项目状态

| 模块 | 行数 | sorry |
|------|------|-------|
| SQF | 1501 | 0 |
| DDF | 435 | 0 |
| EDF | 257 | 0 |
| Hensel | 366 | 0 |
| Recombine | 94 | 2（mignotte + hensel_unique） |
| Pipeline | 379 | 0 |
| 端到端 Zp | 75 | 0 |
| **合计** | **3107** | **2** |

`lake build` 1914 jobs 成功。
