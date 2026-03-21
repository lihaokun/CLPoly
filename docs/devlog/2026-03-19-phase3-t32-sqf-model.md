# Phase 3 T3.2: SQF 算法建模 + 终止性完成

> 日期：2026-03-19
> 分支：`feature/formal-proofs`

---

## 做了什么

在 `CLPoly/Algorithm/SquarefreeZp.lean` 中完成无平方分解（SQF）算法的 L2 模型定义与完整终止性证明。**269 行 Lean，1 sorry**（仅剩 `sqf_correct` 主正确性定理）。

同时完成了 T3.1 DDF 的全部形式化（435 行，0 sorry），Phase 3 当日总产出 **704 行 Lean**。

### 新增/修改文件
- `proof/lean/CLPoly/Algorithm/SquarefreeZp.lean` — **新增** 269 行
- `proof/lean/CLPoly/Algorithm/DDF.lean` — **新增** 435 行（T3.1，0 sorry）
- `proof/lean/CLPoly.lean` — 添加两个 import
- `proof/nl-proof/phase3-t31-ddf.md` — DDF nl-proof
- `proof/nl-proof/phase3-t32-sqf.md` — SQF nl-proof v3（含完整数学论证）
- `docs/devlog/2026-03-19-phase3-t31-ddf-complete.md` — DDF devlog

### T3.2 SQF 核心内容

#### 算法模型（对应 C++ `__squarefree_Zp`）

| 函数 | 行数 | 终止度量 | sorry |
|------|------|---------|-------|
| `yunLoop` | ~35 行定义 + ~40 行终止 | `w.natDegree + c.natDegree` | 0 |
| `sqfZp` | ~25 行定义 + ~25 行终止 | `f.natDegree` | 0 |

`yunLoop` 模型 Yun 内循环（迭代 gcd 分离不同重数因子），返回 `(因子列表, 残余 c)`。
`sqfZp` 模型顶层 SQF（Yun + char p 的 p-th root 递归 via `Polynomial.contract`）。

#### 辅助引理（全部 0 sorry）

| 引理 | 内容 | 行数 |
|------|------|------|
| `natDegree_normalize_eq` | `normalize` 保持 `natDegree`（via `Associated.degree_eq`） | ~10 |
| `normalize_ne_zero_iff` | `normalize a ≠ 0 ↔ a ≠ 0` | ~2 |
| `divByMonic_ne_zero_of_ne_zero` | Monic 精确除法保持非零 | ~5 |
| `natDegree_contract_le` | `natDegree (contract p f) ≤ natDegree f / p` | ~5 |
| `yunLoop_c_natDegree_le` | Yun 循环不增加 c 的 natDegree | ~45 |

## 为什么做

Phase 3 路线图 T3.2。SQF 是 Zp 因式分解管线的第一步（SQF → DDF → EDF），Phase 1 骨架以 `SquarefreeDecomp` 为假设。

## 关键决策

1. **`yunLoop` 需要 `hc : c ≠ 0` 参数**：当 `c = 0` 时 Yun 循环不终止（`gcd(w, 0) = w`，度数不下降）。C++ 实现隐含此前提（`c₀ = gcd(f, f') ≠ 0`）。Lean 中显式传递此证明，每步保持 `hc'`。

2. **显式投影替代模式匹配**：`let (a, b) := e` 在 Lean 4 中不创建 definitional 等价（`b ≠ₐ e.2`）。改用 `let yun_output := yunLoop ...; let c_rem := yun_output.2`，使 `c_rem` 与 `(yunLoop ...).2` definitionally 相等，解决了 `yunLoop_c_natDegree_le` 在 `decreasing_by` 中的应用问题。

3. **`hcrem_le` 放入函数体**：`sqfZp` Branch 2 终止性需要 `c_rem.natDegree ≤ f.natDegree`。在 `decreasing_by` 中无法建立此 bound（变量匹配问题）。解决：在函数体中 `have hcrem_le := yunLoop_c_natDegree_le ...`，使 `decreasing_by` 可直接引用。

4. **`natDegree_contract_le` 证明路径**：`coeff_contract` + `Nat.div_lt_iff_lt_mul` + `coeff_eq_zero_of_natDegree_lt`。

5. **Branch 1 终止性**：`expand_contract` + `natDegree_expand` 得 `deg(f) = deg(contract p f) * p`，配合 `Nat.mul_lt_mul_of_pos_left` 完成。

## 遇到的问题

1. **`let (a,b) := e` vs `.1`/`.2`**：最大的工程障碍。花了多轮调试才发现模式匹配不产生 definitional 等价。
2. **`decreasing_by` 中 typeclass stuck**：`natDegree_contract_le _` 的 `p` 推断失败。用 `@natDegree_contract_le p hp _` 解决。
3. **`omega` 不处理非线性**：`a < a * p`（`p ≥ 2`）需要手动 `Nat.mul_lt_mul_of_pos_left` + `simp [Nat.mul_one]`。
4. **`natDegree_normalize` 不在 Mathlib**：用 `degree_eq_degree_of_associated` + `normalize_associated` 自证。
5. **nl-proof 审核发现 `c_rem' = 0` 证明对不可分因子的处理错误**：Step 2 原声称 ∀ 因子 `p | mⱼ`，实际上不可分因子可以有 `p ∤ mⱼ` 但 `qⱼ' = 0` 直接杀掉导数项。修正为二分论证。

## 涉及的文件

| 文件 | 变更 |
|------|------|
| `proof/lean/CLPoly/Algorithm/SquarefreeZp.lean` | **新增** 269 行 |
| `proof/lean/CLPoly/Algorithm/DDF.lean` | **新增** 435 行 |
| `proof/lean/CLPoly.lean` | 添加 2 个 import |
| `proof/nl-proof/phase3-t31-ddf.md` | **新增** |
| `proof/nl-proof/phase3-t32-sqf.md` | **新增** nl-proof v3 |

## 度量

- 耗时：~10 小时（T3.1 DDF ~6h + T3.2 SQF ~4h，含 nl-proof + 审核 + 形式化 + 调试）
- 迭代：~60 轮编译-修复循环（DDF ~30 + SQF ~30）
- Lean 新增行数：704 行（DDF 435 + SQF 269）
- 对应 C++ 行数：~100 行（DDF ~45 + SQF ~55）
- 放大系数：~7.0x（704/100）
- sorry 状态：DDF 0，SQF 1（`sqf_correct` 主定理，~200 行预估）
- 放弃的方案：
  - `sub_pow_char` Frobenius 路径（DDF，改用 `sub_dvd_pow_sub_pow`）
  - `let (a,b) := e` 模式匹配（SQF，改用显式投影 `.1`/`.2`）
  - `decreasing_by` 内构建 `c_rem` bound（SQF，改为函数体内 `have`）
