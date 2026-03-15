# Phase 0: Lean 4 形式化验证环境搭建与 API 实验

**日期**: 2026-03-14

## 做了什么

完成 Phase 0 的 E0-E4 全部 5 个实验：

1. **E0 项目初始化**: 安装 elan/Lean 4 工具链，创建 `lean/` 项目结构，配置 Mathlib 依赖，获取 8010 个预编译 .olean 缓存文件
2. **E1 ZpPoly API 验证**: 确认 `Polynomial (ZMod p)` 的 Field、EuclideanDomain、GCDMonoid、DecidableEq 实例全部可用
3. **E2 Thm 2.1 桥接路径**: 验证从 Mathlib 有限域理论到 DDF 正确性定理的完整桥接路径
4. **E3 ZMod (p^k) 除法**: 验证 Hensel 提升所需的非整环多项式除法和层间投影
5. **E4 终止性**: 验证 DDF 循环的 `termination_by` 方案和 `.induct` 归纳原理

## 为什么做

在投入 Phase 1（DDF 端到端验证，预估 800-1200 行 Lean 代码）之前，需要确认关键技术假设成立，避免大规模返工。

## 关键决策及理由

### 1. `Fact (Nat.Prime p)` 是 ZMod p 域实例的前置条件
- Lean 不能自动从字面量 `7` 推断素性
- 需要显式声明 `instance : Fact (Nat.Prime 7) := ⟨by decide⟩`
- 泛型代码中用 `variable [Fact (Nat.Prime p)]`

### 2. EuclideanDomain / GCDMonoid 是 noncomputable 的
- 依赖域的除法算法，无法编译为可执行代码
- **不影响证明** — 只影响 `#eval`，证明只需类型和命题
- DDF 循环的 `#eval` 测试需要用纯 `Nat` 版本，不能用多项式版本

### 3. gcd 的正确路径是 `GCDMonoid.gcd`
- `Polynomial.gcd` 不存在
- `GCDMonoid.gcd` 和 `EuclideanDomain.gcd` 都可用
- `normalize_gcd` 确认 gcd 返回归一化（首一）结果

### 4. `divByMonic` 在非整环上可用
- 只需要 `Ring`，不需要 `IsDomain`
- 已在 `ZMod 343` 上实测通过
- 这为 Hensel 提升的 `ZMod (p^k)` 方案提供了关键支持

### 5. DDF 终止性：`termination_by n + 1 - 2 * d` 自动通过
- Lean 4 的 `omega` 策略自动证明递减性
- 不需要手动 `decreasing_by`
- `.induct` 自动生成 4 分支归纳原理，完美匹配 DDF 循环的 if-then-else 结构

### 6. Thm 2.1 桥接路径完整
关键 Mathlib API 已确认：
- `roots_X_pow_card_sub_X`: X^{|K|} - X 的根 = K 中所有元素
- `splits_X_pow_card_sub_X`: X^{|K|} - X 在 K 上完全分裂
- `nonempty_algHom_iff_finrank_dvd`: 域嵌入 ⟺ 维度整除
- `GaloisField.finrank`: GaloisField p n 的 finrank = n
- `expand_card`: Frobenius 幂映射

## 遇到的问题与解决方式

| 问题 | 解决方式 |
|------|---------|
| `elan` 未安装 | `curl` 安装脚本 + `export PATH` |
| `source ~/.elan/env` 失败 | 改用 `export PATH="$HOME/.elan/bin:$PATH"` |
| `Mathlib.Algebra.Squarefree.UniqueFactorizationDomain` 不存在 | 改为 `.NormalizedFactors` |
| `Mathlib.Tactic.Omega` 不存在 | `omega` 是 Lean 4 内置，无需 import |
| `ℕ` 不可用（`autoImplicit = false`） | 使用 `Nat` 或 import Mathlib |
| `let rec` 的 `termination_by` 语法不对 | 改为顶层 `def` 递归函数 |

## 涉及的文件

- `lean/lakefile.lean` — 项目配置
- `lean/CLPoly.lean` — 主入口
- `lean/CLPoly/Experiment/E1_ZpPolyAPI.lean` — API 验证
- `lean/CLPoly/Experiment/E2_TheoremBridge.lean` — Thm 2.1 桥接
- `lean/CLPoly/Experiment/E3_ZModPkDiv.lean` — ZMod p^k 除法
- `lean/CLPoly/Experiment/E4_Termination.lean` — 终止性验证
- `proof/docs/phase0-experiment-plan.md` — 实验计划
- `proof/docs/research-report.md` — 调研报告
