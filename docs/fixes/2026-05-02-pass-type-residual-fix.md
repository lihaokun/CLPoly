# 修正方案：Pass 上游类型残留（RefType / UnknownType）系统性修复

日期：2026-05-02
状态：调研完成，方案待确认
关联：Pass 8 codegen B1 阶段后期 sorry 数 90+ 不降的根因；docs/devlog/2026-05-02-pass8-codegen-b1.md known issues

## 问题陈述

Pass 8 codegen 输出的 Corpus.lean 含 90 个 sorry，主要来自上游 Pass 1-7 的
类型残留（RefType / UnknownType / Lambda NamedType）。

继续修 Pass 8 emit 路径已属 whack-a-mole（每修一处暴露下层更多）；必须
回流上游修类型传播。

## 根因调研

### 数据：各 Pass 阶段残留计数（64 函数 corpus，不含 factorize）

| 阶段 | RefType | UnknownType | Lambda残留 |
|------|---------|-------------|------------|
| HIR0 parse | 224 | 14 | 98 |
| HIR1 ref_elim | 224 | 14 | 98 |
| HIR1b callsite_ref_elim | 227 (+3) | **252 (+238)** | 98 |
| HIR2 lambda_lift | 229 (+2) | 257 (+5) | 105 (+7) |
| HIR2b lambda_ref_elim | 229 | 303 (+46) | 93 (-12) |
| HIR3 iter_recognize | 216 (-13) | 331 (+28) | 99 (+6) |
| HIR4 operator_resolve | **131 (-85)** | 324 (-7) | **49 (-50)** |
| MIR0 ssa_build | **2 (-129)** | 5 (-319) | 1 (-48) |
| MIR1 loop_lower | **27 (+25)** | **64 (+59)** | 1 |

### 根因 #1: Pass 2b callsite_ref_elim 一次引入 238 个 UnknownType

HIR1 14 → HIR1b 252（涨 17 倍）。Pass 2b 改写 caller 端 ref 参数为
destructure，新增的 LetStmt 类型没正确传播 → 默认 UnknownType("")。

具体场景（推测，待 trace 确认）：
- `f(out_ref) → let __ret := f(out_ref); let out := __ret.snd`
- destructure stmt 的 var.ty 是 UnknownType（应该从 callee output 推）

### 根因 #2: Pass 7 loop_lower 反向引入 25 RefType + 59 UnknownType

MIR0→MIR1：本应只重组结构，但 Pass 7 反而**新增**类型残留：
- `_build_loop_func`：cap_param.ty 来自 free_var.ty（可能含 RefType 残留）
- `_splice_loop_call`：destructure LetStmt 用 `lo.ty or UnknownType("")`，
  若 live_out var 类型不全 → UnknownType
- `_make_exit_return`：`reaching_def_at` 返回的 Var ty 可能未恢复

具体修复点（推测）：
- Pass 7 emit 时 strip RefType（`lo.ty.inner if isinstance(lo.ty, RefType)`）
- live_out 类型从 phi/Cast 链继承而非默认 UnknownType

### 根因 #3: Pass 5 operator_resolve 大幅减少残留（85 RefType + 50 Lambda）

Pass 5 是当前类型清扫主力——但仍漏 131 RefType。这些是 Pass 5 cast_table /
constructor_map 未覆盖的形态（如复合 PairType<Ref<X>, Y>）。

## 影响分析

| 影响目标 | 当前 | 修复后预期 |
|---------|------|-----------|
| Corpus.lean sorry 数 | 90 | 30-50（剩余是真深 silent miss） |
| lake build errors | 10 | 5-7（少几条 Pass 1 残留导致的 type mismatch） |
| MIR1 RefType | 27 | 0 |
| MIR1 UnknownType | 64 | 10-20（Pass 1 真未识别） |

## 修复策略：分两阶段

### 阶段 A: Pass 7 反向回归修复（高 ROI，1-2 commits）

**目标**：MIR1 的 RefType/UnknownType 不应高于 MIR0。

**具体动作**：
1. `_build_loop_func` cap_param 创建处：strip RefType
   ```python
   cap_params = [HIRParam(name=..., ty=_strip_ref(v.ty), ...) for v in free_vars]
   ```
2. `_splice_loop_call` destructure：从 callee ret_ty 推导而非用 `lo.ty`
3. `_make_exit_return` reaching_def 处：保留原 Var ty（不让 RefType 漏过）

**测试**：跑 corpus 残留扫描，MIR1 RefType 应为 0，UnknownType ≤ MIR0

### 阶段 B: Pass 2b callsite_ref_elim UnknownType 修复（中 ROI，2-3 commits）

**目标**：HIR1b 新增 UnknownType 从 238 降到 < 50。

**具体动作**：
1. trace 哪条改写路径生成的 destructure LetStmt 没传 ty
2. 在 destructure 时用 callee 的 output_param 类型查表

**注意**：这是 Pass 2b 的二阶段重构，可能需要先扩展 TRANSLATION_SCOPE_OUTPUT_PARAMS
表加 output 类型信息。

### 阶段 C（必做）：Pass 1/3/4/5 long-tail 残留逐 case 修

剩 131 RefType（在 HIR4 中）→ Pass 5 cast 模板覆盖不全 + Pass 1 类型上下文丢失。

**子阶段**：
- C1: HIR1b 之外的 RefType 来源 trace（Pass 1 parse 直接产 224 个 base，
  Pass 5 reduces to 131）→ 修 Pass 5 RefType strip 路径
- C2: HIR2/HIR3 期间累积的 ~80 UnknownType（Pass 3/3b/4 改写未传 ty）
- C3: 49 Lambda residual @ HIR4 → Pass 3 lift 漏过的 case

每个 sub-stage 独立 1-3 commits。无单一修复点，但属于真 silent miss
必须收敛。预估 5-10 commits 总。

## 决策（用户确认）

A → B → C 全部做，每阶段后做 1-对-1 审视。

## 测试方案

每阶段后：
1. `bash tests/run_all_smoke.sh`（67/67 通过）
2. 残留扫描脚本（本文档对应的 scan）：MIR1 RefType / UnknownType 数应下降
3. `lake build CLPoly/Generated/Corpus.lean`：errors 不增（应少 1-3 条）
4. `proof/lean/CLPoly/Generated/Corpus.lean` sorry 数：应下降 ≥ 30

## 风险

- Pass 7 修改可能让 free_var / live_out 推断错（已经在第 11/12 轮审视
  发现过类似 silent miss）→ 修后必须再做 1-对-1 审视（Lesson A）
- Pass 2b 重构面较大，可能引入新 silent miss

## 与项目主线的关系

按 proof/CLAUDE.md "L2 禁止简化" 原则——上游残留是真 silent miss，必须
修；不能用下游 Pass 8 sorry 兜底掩盖。

## 实施顺序（已确认）

阶段 A → B → C 全部做。每阶段后做 1-对-1 审视（Lesson A 模式）。
