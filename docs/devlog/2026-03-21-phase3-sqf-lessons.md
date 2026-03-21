# Phase 3 T3.1-T3.2 完成：DDF 全验证 + SQF 近完成 + 形式化教训

> 日期：2026-03-19 ~ 2026-03-21
> 分支：`feature/formal-proofs`

---

## 做了什么

完成 Phase 3 L2 算法模型的两个核心定理：

| 定理 | 文件 | 行数 | sorry | 状态 |
|------|------|------|-------|------|
| T3.1 DDF | `Algorithm/DDF.lean` | 435 | **0** | ✅ 完全验证 |
| T3.2 SQF | `Algorithm/SquarefreeZp.lean` | 1034 | **4** | 🔶 ~97% |

**合计 1469 行 Lean**，对应 C++ ~100 行（`__ddf_Zp` + `__squarefree_Zp`），放大系数 ~14.7x。

### DDF（0 sorry，一次成功）

- `ddfLoop` 递归函数 + 终止性
- `ddfLoop_correct` 携带 7 个不变量的归纳证明
- `ddf_correct` 顶层定理
- 关键辅助引理：`h_cong_step`（Frobenius 替代路径）、`gd_irred_characterization`（核心桥接）

### SQF（4 sorry，结构完成）

**已验证（0 sorry）的核心引理**：
- `squarefree_div_gcd_derivative` — f/gcd(f,f') 是 squarefree（valuation 论证 + 无穷递升矛盾）
- `yunLoop_correct` — Yun 循环 10 个不变量（Y1-Y6-Y10）的完整保持证明
- `yunLoop_extracts_factor` — q | w₀ → q 被 Yun 提取到某个 acc 条目
- `yunLoop_acc_subset` — yunLoop 只追加 acc，不删除
- `yunLoop_c_natDegree_le` / `yunLoop_c_ne_zero` — yunLoop 的 c 度数单调 + 非零保持
- `sqfDecomp_pow_lift` — p-th root 下的 SquarefreeDecomp 提升
- `expand_eq_pow` — Frobenius: expand p f = f^p in F_p[X]
- `pow_dvd_derivative_of_pow_succ_dvd` — d^{k+1} | f → d^k | f'
- `sqf_correct` Case A（deg=0）、Case B（f'=0）、Case C.2（c_rem deg=0）

**4 sorry 的根因分析**：
1. `derivative(c_rem) = 0`（1 个数学 sorry）— 需要 UFD 分解 + 逐因子导数分析
2. Case C.1 的 Associated 链 + List.mem_append 匹配（3 个 plumbing sorry）— Lean 表达式匹配

---

## 关键教训（本次最大收获）

### 教训 1：normalize/Associated 是隐藏的复杂度炸弹

**现象**：`EuclideanDomain.gcd` 返回非 monic 多项式。`normalize` 使其 monic，但引入 `Associated`（∃ unit, a * u = b）而非精确等式。**每一个**涉及 gcd 的 dvd 证明都需要 `normalize_dvd_iff`、`normalize_associated` 桥接。

**代价**：SQF 中约 60% 的调试时间花在 Associated 算术上（特别是 Y1 乘积不变量的保持证明）。

**对比**：DDF 用了类似的结构但更简单（只有 1 层 normalize），SQF 有 2 层嵌套（yunLoop 的 y 和 c' 都经过 normalize），导致 Associated 链长度翻倍。

**教训**：在设计 Lean 算法模型时，**尽量减少 normalize 的使用点**。如果可能，在函数签名中直接要求 Monic 输入，避免运行时 normalize。

### 教训 2：`let (a, b) := e` 是定时炸弹

**现象**：`let (a, b) := e` 在 Lean 4 中编译为 `match`，`a` 和 `e.1` 不是 definitionally equal。这意味着 `yunLoop_c_natDegree_le` 对 `(yunLoop ...).2` 的 bound 不能直接用于 `c_rem`（从 `let` 解构得到）。

**代价**：发现这个问题花了 ~5 轮编译-调试。修复方案：改用 `let output := e; let a := output.1; let b := output.2`。

**教训**：**在递归函数定义中永远不用模式匹配 let**。始终用 `.1`/`.2` 投影。

### 教训 3：`rw` 在递归目标中是危险的

**现象**：`rw [hcrem_eq_pow]`（将 `c_rem` 替换为 `(contract p c_rem)^p`）会替换目标中**所有** `c_rem`，包括 `sqfZp(contract p c_rem)` 内部的——导致 `sqfZp` 无限展开或类型不匹配。

**代价**：Case B 和 Case C.1 的 Associated 证明反复失败，~10 轮调试。

**教训**：**对含递归函数调用的目标，永远不用 `rw`**。改用 isolated `have` 块内 `rw`，或 `Associated.trans` 手动链接。

### 教训 4：应该先攻最难的 sorry

**现象**：我按自底向上的顺序工作（函数定义 → 终止性 → 辅助引理 → 主定理 Case A/B → Case C），最后才遇到 `derivative(c_rem) = 0`（最难的数学引理）和 Case C.1 的 Associated 拼接（最复杂的 plumbing）。

**代价**：前 80% 进展很快（看起来快完成了），后 20% 花了同样多的时间。

**教训**：nl-proof 完成后，**先识别形式化的关键瓶颈**（通常是：最深的数学引理 + 最复杂的 Associated/normalize 链），从那里开始。如果瓶颈不能解决，尽早暴露。

### 教训 5：`Nat.strongRecOn` 远优于 `f.induct`

**现象**：`yunLoop.induct` 的 inaccessible name 数量不可预测（let 绑定也算，且不同 case 数量不同）。`rename_i` 需要试错。而 `Nat.strongRecOn` 完全掌控命名。

**代价**：`yunLoop.induct` 尝试花了 ~1 小时，切换到 `Nat.strongRecOn` 后 ~15 分钟完成。

**教训**：不变量多于 5 个的递归证明，**直接用 `Nat.strongRecOn`**。

### 教训 6：复杂度估计需要 3x 安全系数

| | DDF 估计 | DDF 实际 | SQF 估计 | SQF 实际 |
|---|---|---|---|---|
| 行数 | ~310 | 435 | ~340 | 1034 |
| sorry | 0 | 0 | 0-1 | 4 |
| 耗时 | ~6h | ~6h | ~4h | ~15h+ |

DDF 估计准确（1.4x 误差）。SQF 估计严重不足（3x 行数误差，sorry 目标未达成）。原因：SQF 的嵌套递归 + normalize 摩擦 + 2 层 p-th root 递归是乘法而非加法复杂度。

**教训**：有嵌套递归 + normalize 的定理，**复杂度估计乘以 3**。

---

## 涉及的文件

| 文件 | 变更 |
|------|------|
| `proof/lean/CLPoly/Algorithm/DDF.lean` | **新增** 435 行 |
| `proof/lean/CLPoly/Algorithm/SquarefreeZp.lean` | **新增** 1034 行 |
| `proof/lean/CLPoly.lean` | 添加 2 个 import |
| `proof/nl-proof/phase3-t31-ddf.md` | **新增** DDF nl-proof |
| `proof/nl-proof/phase3-t32-sqf.md` | **新增** SQF nl-proof v3（3 轮审核） |
| `proof/CLAUDE.md` | 新增"Lean 形式化经验教训"段 |

## 度量

- 耗时：~20 小时（DDF ~6h + SQF ~14h，含 nl-proof 草稿 + 3 轮审核 + 形式化 + 调试）
- 迭代：~120 轮编译-修复循环（DDF ~30 + SQF ~90）
- Lean 新增行数：1469 行（DDF 435 + SQF 1034）
- 对应 C++ 行数：~100 行（DDF ~45 + SQF ~55）
- 放大系数：14.7x（1469/100）
- 放弃的方案：
  - `sub_pow_char` Frobenius（改 `sub_dvd_pow_sub_pow`，避免 `CharP (Polynomial (ZMod p))` 实例传递）
  - `let (a,b) := e` 模式匹配（改 `.1`/`.2` 投影，避免 definitional equality 问题）
  - `yunLoop.induct` functional induction（改 `Nat.strongRecOn`，避免 inaccessible name 试错）
  - `rw [h]` 在递归目标中（改 isolated `have` + `Associated.trans`，避免无限展开）
  - `simp [sqfZp]` 展开递归（改 `unfold sqfZp; split`，避免 looping simp）
  - `squarefree_iff_emultiplicity_le_one` 路径（改直接 `d^n | f` 无穷递升 + 度数矛盾，更初等）
