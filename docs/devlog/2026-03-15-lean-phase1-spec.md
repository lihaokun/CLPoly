# Lean 4 形式化验证：Phase 1 启动 — 接口规约 + 文档审计

> 日期：2026-03-15
> 分支：`feature/gcd-optimization-p1`
> 前置：Phase 0 实验完成（`b1b4627`）

---

## 1. 文档审计与修正

对验证架构的三份核心文档进行一致性审计，发现并修正 24 处不一致：

| 文档 | 修正数 | 主要问题 |
|------|--------|----------|
| `proof/docs/implementation-roadmap.md` | 10 | L2 描述混入"1:1 对应 C++"措辞 |
| `proof/docs/lean-verification-blueprint.md` | 13 | L2/L1 术语混淆 + §6 Phase 定义与 roadmap 完全脱节 |
| `proof/CLAUDE.md` | 1 | `Repr（L1）` → `Impl（L1）` |

### 关键修正

- **L2 算法模型**：统一为"算法逻辑，操作 Mathlib 类型，抽象掉 C++ 细节"
- **L1 实现模型**：统一为"1:1 对应 C++（uint64 语义、数组越界、move）"
- **Blueprint §6**：删除过时的 Phase 0-4 任务列表，替换为指向 `implementation-roadmap.md` 的参考表

## 2. Phase 0 实验重现

重新运行 4 个 Phase 0 实验，全部通过：

| 实验 | 内容 | 结果 |
|------|------|------|
| E1 | ZMod p 多项式 API 可用性 | ✓ |
| E2 | Mathlib 定理桥接（squarefree、gcd） | ✓ |
| E3 | ZMod (p^k) 除法/castHom | ✓ |
| E4 | 终止性证明模式 | ✓ |

## 3. Phase 1 T1.1：接口规约 `Spec.lean`

新建 `proof/lean/CLPoly/Spec.lean`，定义 7 个子过程的接口规约（Prop）：

### Zp[x] 子过程

| 规约 | 签名 | 核心性质 |
|------|------|----------|
| `SquarefreeDecomp` | `f → [(factor, mult)]` | Associated 还原 + 无平方 + 首一 + 互素 |
| `DDFCorrect` | `f → [(gd, d)]` | 整除 + 度数分组 + 不遗漏 + 还原 |
| `EDFCorrect` | `g, d → [factor]` | 还原 + 不可约 + 首一 + 度 = d |

### Z[x] 子过程

| 规约 | 签名 | 核心性质 |
|------|------|----------|
| `HenselCorrect` | `f, k, facs_p, facs_pk` | 模 p^k 乘积 ≡ f + 数量一致 + 度数保持 |
| `RecombineCorrect` | `f → [factor]` | Associated 还原 + 不可约 |

### 顶层

| 规约 | 签名 | 核心性质 |
|------|------|----------|
| `FactorZpCorrect` | `f → (lc, [(factor, mult)])` | 精确还原 + 不可约 + 首一 |
| `FactorZZCorrect` | `f → [factor]` | Associated 还原 + 不可约 |

### 技术要点

- `Polynomial.natDegree`（非 `natDeg`）：Mathlib 4 的完整名称
- `List.Forall₂` 替代 `Fin` 索引：避免 `omega` 无法跨假设推导长度相等的问题
- `Associated` 用于乘积还原（允许单位差异，比精确等号更合理）

## 度量
- 耗时：~4 小时（文档审计 ~1.5h、Phase 0 重现 ~0.5h、Spec.lean 编写 ~2h）
- 迭代：3 轮编译-修复循环（Spec.lean 的 typeclass / 签名调整）
- Lean 新增行数：133 行（Spec.lean）
- 对应 C++ 行数：N/A（纯规约，无对应实现）
- 放弃的方案：最初用 `Fin n → Polynomial` 索引因子列表，改为 `List.Forall₂` 避免长度推导困难

## 4. 下一步

- **T1.2**：`CLPoly/Pipeline/FactorZp.lean` — Zp[x] 顶层正确性骨架
- **T1.3**：`CLPoly/Pipeline/FactorZZ.lean` — Z[x] 顶层正确性骨架
- **T1.4**：接口一致性审查 + 锁定
