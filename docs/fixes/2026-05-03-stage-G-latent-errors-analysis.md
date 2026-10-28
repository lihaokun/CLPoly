# 阶段 G — Latent 102 errors 根因分析 + cpp2lean v2 整体审视

**日期**：2026-05-03
**前置**：阶段 F（已 commit ac757f7）— 原 12 errors 全清，但 F8（FuncType IR 节点 + LetStmt ty 同步）让 Lean 不再因前 5 个 hard error 提前 abort 函数检查，暴露出 **102 个深埋的 latent errors**。

本文档：
1. 把 102 errors 按根因归类（每条都对应到 Pass）
2. 反思 cpp2lean v2 的设计/Pass 协作中导致 latent 累积的系统性问题
3. 给出 stage G 的修复路径与下一轮的优先级建议

## 0. 当前状态量化

| 指标 | 值 |
|------|-----|
| sorry | 0 |
| lake errors | 102 |
| corpus 行数 | 6816 |
| 已修原方案错误 | 9 errors + 4 sorries（原 12 全清）|

错误分布（按 Lean 报错关键字）：

| 关键字 | 数量 | 占比 |
|--------|------|------|
| Application type mismatch | 40 | 39% |
| Function expected | 22 | 22% |
| Type mismatch | 20 | 20% |
| Unknown identifier | 8 | 8% |
| Invalid field | 3 | 3% |
| failed to synthesize | 2 | 2% |
| Unknown constant | 2 | 2% |
| Invalid projection | 2 | 2% |
| 其他 | 3 | 3% |

## 1. 根因归类

我从 102 errors 中按文件位置取样了 ~10 处具体错误，做了上下文检查。共**6 个根因簇**，每个簇关联 1-3 个 Pass 的 silent miss。

### 根因 A — local-var lambda 多参调用（~15 errors）

**症状**：`Unknown identifier "row_sub" / "row_swap" / "lc_correct" / "compute_theta" / "next_p" / "upzp_coeff"` 等。

**示例** (1783)：

```lean
-- 调用点：
let __sideeff_38_0_1 := (row_sub_1 k_2 j_11 q2_1)  -- ❌ row_sub_1 在 caller scope 不可见

-- 实际定义在 outer 函数：
-- (在 __lll_reduce_ir 内)
let row_sub_1 := __hoist_lam_0_1.1
```

**根因链**：

1. C++ 用一个 IIFE 模式返回多个 lambda 的元组（典型的 LLL reduction 代码风格）：
   ```cpp
   auto [row_sub, row_swap, dot, ...] = []() {
       auto row_sub = [&](int k, int j, ZZ q) { ... };
       auto row_swap = [&](int i, int j) { ... };
       return std::tuple{row_sub, row_swap, ...};
   }();
   ```
2. Pass 3 lift **只 lift 外层 lambda**，内层 mini-lambda 留为 LambdaExpr 残留在 body 内。
3. Pass 6 SSA build 不对 `Call(callee=str)` 的 callee 做 SSA 重命名。
4. Pass 7 抽 inner loop 时，`row_sub_1` 是 outer 函数的 SSA var，在 inner loop 函数的 cap_params 列表内未被识别（callee str 不被 `_collect_var_reads_in_expr_versioned` 当作 Var read）。

**已尝试**：
- 阶段 F19/F20/F21：Pass 6 把 0-arg `Call(local_var, [])` 转 Var 引用 — **解决了 0-arg case（compute_error）但未解决多参 case（row_sub k j q）**。
- Pass 7 regex 解析 `<base>_<v>$` callee 加进 free vars — **回滚**（false positive `next_prime_64` 把全局函数误当 local）。

**正解**：

需要让 IR 的 `Call.callee` 接受 `Var` 类型（而非仅 `str | UnresolvedOp`）。Pass 1 检测 callee 是 DeclRefExpr to local var → emit `Call(callee=Var, args=...)`。然后 Pass 6 SSA 自动重命名，Pass 7 free var 自动收集。

**改动量估计**：~40 行（IR + Pass 1 + Pass 6 + Pass 8 emit_call），中等风险（IR 层改动需各 Pass 透传 isinstance 判断）。

---

### 根因 B — IIFE 元组解构 nested lambda（~12 errors）

**症状**：`__hoist_lam_X.snd has type Array (Array Int) which has only 1 field`，`Invalid projection: Index 2 is invalid`。

**示例** (1953-1959)：

```lean
let __hoist_lam_0_1 := (_lambda___lll_reduce_3_ir M U_2 n_1 M U_2 n_1)
let M_1 := __hoist_lam_0_1.2  -- 期望 4-tuple .2 但实际是 1-tuple
let U_4 := __hoist_lam_0_1.2.2  -- 同样
let row_sub_1 := __hoist_lam_0_1.1
```

**根因链**：

`_lambda___lll_reduce_3_ir` 应该返回 `(row_sub, M, U) : (Func × LLLMatrix × LLLMatrix)`（3-tuple），但 Lean 端实际推断为 `Array (Array Int)`（即 LLLMatrix）— 因为 lifted lambda 的 ret_ty 没正确 inference。

C++ source（推断）：

```cpp
auto [row_sub, M, U] = []() {
    auto row_sub = [&](int k, int j, ZZ q) {
        for (int i = 0; i < n; i++) {
            M[k][i] = M[k][i] - q * M[j][i];
            U[k][i] = U[k][i] - q * U[j][i];
        }
    };
    return std::tuple{row_sub, M, U};
}();
```

Pass 3 lift 外层 lambda 时，**ret_ty 推断**只看 ReturnStmt.value，不识别 `std::tuple{...}` 表达式（被 Pass 1 解析为某种 Call 或 ArrayLit）。结果 `_lambda_..._3_ir` 的 `ret_ty = LLLMatrix`（最后一个元素的类型），不是元组。

**与根因 A 的关系**：根因 A 是"local var 调用没 SSA bump"，根因 B 是"lifted lambda ret_ty 推断错"。两者叠加导致 `__hoist_lam_0_1.X` navigation 失效。

**正解**：

Pass 3 `_lift_lambda` 中 `_find_ret_ty` 需识别 `std::tuple` / `std::pair` / `{a, b, c}` ArrayLit 模式，构造 PairType / TupleType 作为 ret_ty。需要 ArrayLit 至 TupleType 的强制类型识别（ArrayLit 在外层 ret 位置时是 tuple，不是 array）。

**改动量估计**：~30 行（Pass 1 + Pass 3 ret_ty 推断），低风险但需测试。

---

### 根因 C — vec.assign 模板参数顺序错（~10 errors）

**症状**：`Array.replicate U ...` 中 U 是 LLLMatrix（Array），但期望 Nat。

**示例** (1947)：

```lean
let U_1 : LLLMatrix := (Array.replicate U ((n_1).toInt64.toUInt64) ...)
-- Lean Array.replicate : (n : Nat) → α → Array α
-- 但这里第一个 arg `U` 是 LLLMatrix，第二个才是 Nat —— 顺序反了
```

**根因**：

C++ `vec.assign(n, val)` 在 class_map line 109/137 映射为 `("mutate", "Array.replicate")`. Pass 5 的 mutate template 把 receiver 当作 first arg，emit `Array.replicate <receiver> <args[0]> <args[1]>`. 实际 `Array.replicate n val` 不需要 receiver，但 mutate 模板硬塞 receiver。

应该是 `Array.replicate args[0] args[1]`（去掉 receiver）。

**正解**：

class_map 的 "mutate" 语义改为：模板 substitution 时，receiver 作为 LHS（assignment target）但**不作为 args**。或者在 Pass 8 处理 mutate 模板特殊化，emit `<receiver> := <template apply args>` 时不把 receiver 注入 args。

或更直接：把 `Array.replicate` 改为一个 wrapper `Array.replicateMut` 接受 receiver + n + val（receiver 忽略）。

**改动量估计**：~10 行，低风险。

---

### 根因 D — comp_ptr / Lex / Variable 类型混乱（~8 errors）

**症状**：`comp_ptr_1 has type Lex but expected Variable × Int64` 等。

**示例** (2063)：

```lean
let var_poly_1 : Poly := (MvPolyZZ.mk (MvPolyZZ.comp gk_1))
let mono_1 : Monomial := (Monomial.empty)
let mono_2 : Monomial := (Array.push mono_1 ((var_1, (1 : Int32))))
-- 接下来：
let lc_tau_zp_1 : Array PolyZp := (Array.replicate ... (MvPolyZp.mk comp_ptr_1))
-- comp_ptr_1 : Lex, 但 MvPolyZp.mk 现在是 {α} (_ : α) → MvPolyZp，所以 α=Lex 应该 OK
-- 但 .empty / Monomial 等 stub 让类型链卡住
```

**根因**：

`comp_ptr` 在 C++ 是 `lex_<less>` 比较器对象（Pass 1 typedef alias `Lex`）。它在 polynomial 构造里以 `comp_ptr` 形式传入，但 Lean Model 把 `Lex` 定义为 `Unit`（占位），不与 Variable / Monomial 兼容。

**正解**：

无需让 `Lex` 与 Variable 兼容。但要让 `MvPolyZp.mk` / `Monomial.mk` 接受任意类型 α（已通过 `{α} (_ : α) → ...` 实现）。具体每条 site 看是 cast 错还是 stub 缺失。

**改动量估计**：~15 行 stub + 各 site 微调，低风险。

---

### 根因 E — Pass 5 cast_table 路径漏（~10 errors）

**症状**：`f has type MvPolyZZ but expected SparsePolyZZ`，`Function expected at poly_convert ...`，`HPow Zp Int64 Zp 缺失`。

**示例** (2018)：

```lean
let __mignotte_bound_upoly_ir f
-- f : MvPolyZZ 但 __mignotte_bound_upoly 期望 SparsePolyZZ
```

**根因**：

`__mignotte_bound_upoly_ir` 的 caller 是某个 multivar 函数，传入了 `MvPolyZZ` 但 callee 期望 `SparsePolyZZ`。这是 **Pass 1 类型推断错** —— 函数实例化时 Pass 1 没区分 caller 是 univar 还是 multivar instance。

类似 (2412) `Function expected at poly_convert sigma0_2[...]! result_3[...]!`：`poly_convert` 被多余地 apply 给两个 args，但 Lean 端 `poly_convert {α β} (f : α) (target : β) : β := target` 已饱和 — 第三个 arg 不能 apply。

**正解**：

需要分类逐条 trace：每个 site 看是 Pass 1 类型推断错、Pass 5 cast 漏，还是 class_map 函数签名错。

**改动量估计**：每个 5-10 行，~30-50 行总计。

---

### 根因 F — HPow / 类型类实例缺失（~5 errors）

**症状**：`failed to synthesize HPow Zp Int64 Zp`。

**正解**：补 Lean Model 缺少的类型类实例。

```lean
instance : HPow Zp Int64 Zp where
  hPow base e := ... -- base ^ e.toNat
```

**改动量**：~15 行 stub。

---

## 2. 反思：cpp2lean v2 的系统性问题

这次阶段 F 的工作让我重新审视 cpp2lean v2 的设计。三个**结构性问题**导致了 latent errors 的累积。

### 问题 1：Pass 早 abort + 单错放大

**现象**：Lean 因前 5 个 hard error abort 整个函数的 elaboration，使函数后续 stmts 中的 ~100 个 bug 被掩盖。F8 修了第一处，Lean 不再 abort，瞬间暴露 100+ 错误。

**反思**：

- v2 之前每次 commit 都看 lake errors 数，但 errors 数 ≠ 真实正确性。**前几个 hard error 修掉后，后续 stmts 的错误才会浮现**。
- 应该有更可靠的进度度量。例如：
  - 每个 top-level 函数的"独立"errors（不被其他函数影响）
  - 函数级的 elaboration 通过率
  - 模拟 Lean 的"独立 elaboration"，每个函数单独编译看错

**改进建议**：把每个 corpus 的顶层 def 写成独立 `.lean` 文件（或用 `set_option maxErrors 1000`），每次 build 报告 per-function 状态。

### 问题 2：IR 类型表达力不足导致绕路

**现象**：

- `Call.callee` 只接受 `str | UnresolvedOp`，不接受 `Var` → 局部变量调用没 SSA rename → 一类深 silent miss（根因 A）。
- 没有 `FuncType` IR 节点（直到 F8 才补）→ lifted lambda 类型用 `LambdaRef = Unit` 占位 → 大量类型链失效（根因 B、C 部分）。
- 没有 `RecordType`（结构体类型） / `TupleType` 与 `PairType` 不统一 → tuple navigation 路径混乱（n=2 用 `.fst/.snd`，n>2 自定义路径）。

**反思**：

cpp2lean v2 的 IR 是"够用就好"的增量设计——遇到新模式才补 IR 节点。这导致：

- 新 IR 节点加上后，所有 Pass 需要补 isinstance 透传（不补就炸）。
- 旧 Pass 可能在缺少 IR 节点时用 NamedType 字符串绕路，留下 silent miss。

**改进建议**：

- 提前规划 IR 全集——`FuncType`、`RecordType`、`UnionType`、`Capture`、`PartialApply`、`LocalVarCall` 等节点应该在初始设计时就定义。
- 每加一个 IR 节点，立即在所有 Pass 写"显式 unhandled"分支（assert / sorry），避免静默 fallback。

### 问题 3：Pass 间约定不严格（语义模糊地带）

**现象**：

- `LambdaRef` 占位被多个 Pass 各自解读：Pass 3 设为 caller-side var ty，Pass 4 设为 filter pred lambda ty，Pass 8 emit 为 `Unit`（让 caller param 不报错但失去类型信息）。
- "mutate" 模板的语义在 Pass 5 和 Pass 8 之间不一致（receiver 是注入还是不注入）。
- `_lambda_<host>_<n>` 的命名约定 vs 内层 mini-lambda 的命名（根因 A 中 IIFE 元组）。

**反思**：

每个 Pass 的 invariant 文档零散，且没有自动化 check。某 Pass 的修改容易在后续 Pass 触发 silent miss。

**改进建议**：

- 每个 Pass 边界写**类型级 invariant**：例如 "Pass 3 输出的 HIRFunc 中所有 LambdaExpr 已被 lift；所有 lifted lambda 的 ret_ty 都不是 LambdaRef"。
- 加 invariant check 步骤（dry-run 检查），在 corpus build 中每 Pass 结束后运行。

## 3. Stage G 修复路径

按根因簇优先级排序：

| 簇 | 数量 | 优先级 | 预估行数 | 风险 |
|----|------|-------|---------|------|
| C (vec.assign 模板) | ~10 | 🔴 高（最简单） | 10 | 低 |
| F (HPow 类型类) | ~5 | 🔴 高 | 15 | 低 |
| B (IIFE ret_ty 推断) | ~12 | 🟡 中 | 30 | 低 |
| D (类型混乱) | ~8 | 🟡 中 | 15+各 site 微调 | 中 |
| E (杂项 cast/Pass 1 推断) | ~10 | 🟡 中 | 30-50 | 中 |
| A (local-var lambda 多参) | ~15 | 🟢 低（深改动） | 40 + 多 Pass 透传 | 中 |
| 其他散点 | ~42 | 🟢 低 | 各自 trace | 不定 |

**推荐修复顺序**：C → F → B → D → E → A → 散点

理由：
- C/F 是**最简单的 Lean Model + class_map 改动**，能快速消除 ~15 errors，腾出 elaboration 空间
- B 解决 IIFE pattern，能根本性消除根因 A 的部分 case（不解 IIFE，根因 A 也不能完全修）
- A 是最深改动（IR + 多 Pass 透传），最后做风险最小（前面已修的 case 可作为参考点）

## 4. 整体审视

### 已交付（阶段 A-F，2 周累计）

| 阶段 | 主要内容 | 错误轨迹 |
|------|---------|---------|
| A (5月2日) | Pass 7 反向回归 + Pass 6 _with field 修 | sorry 90→42 |
| B | Pass 2b __refret tmp ty 推导 | sorry 42→7 |
| C | long-tail Pass 1/4/6/7 类型推断补全 | sorry 7→0 |
| D | 类型残留 7 sorry → 0 | lake errors 16→11 |
| E (5月3日上半) | typedef alias、Variable/Monomial 统一、HenselNode 结构化 | errors 102→12 |
| F (5月3日下半) | 9 个原方案修复 + FuncType IR 节点 + filter partial-app | sorry 0、errors 12→102（露 latent） |

**真实进展**：
- 原 90+ sorry → 0
- v1 时代的"打包 sorry"全部根因解决
- IR 加入 FuncType、Pass 6 SSA 改进 callee rename、Pass 4 filter capture 联动

**未完成**：
- 102 latent errors（今天暴露的 corpus 真实 bug，根因如上 6 簇）
- 即使阶段 G 全清，可能 Lean 又会暴露下一波（需要做 problem 1 的 per-function 检查机制）

### 翻译器质量评估

**优点**：

- IR 设计现已完整覆盖 C++ 主要构造（结构体、lambda、ref-out、SSA、loop lower）
- Pass 边界清晰（10 个 Pass 各自负责一种转换）
- 每个 Pass 都有 smoke test + audit doc

**缺点**：

- **正确性度量不可靠**（前述问题 1）
- **某些 IR 节点过晚补充**（如 FuncType 直到阶段 F 才有），积累了大量绕路代码
- **Lean Model 是手工补 stub 的累积**，缺乏"翻译器自动生成 stub"的机制（每次新 lake error → 手工 grep + 加 def）

### 短期 vs 长期路径

**短期（阶段 G）**：按根因 C → F → B → D → E → A 顺序攻 102 errors。预估 5-8 小时。

**长期（v3 重构）**：

- 重新设计 IR：FuncType / RecordType / LocalVarCall 等一次性补齐
- 每个 Pass 严格 invariant + check
- Lean Model 自动生成（从 C++ 头文件生成 stubs，不手工维护）

是否启动 v3 由用户决定。短期阶段 G 是为了把 v2 收尾到"lake build 通过"的可交付状态。
