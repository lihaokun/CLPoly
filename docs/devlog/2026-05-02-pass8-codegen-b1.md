# Pass 8 codegen B1 阶段完成 — MIR₁ → Lean 4 (lake build 101→11)

日期：2026-05-02

## 做了什么

实现 cpp2lean v2 Pass 8（MIR₁ → Lean 4 源码），完成 B 阶段（happy-path
子集 → 全 corpus 聚合 → 错误收敛），通过 9 个 commit 把 lake build 错误从初
始 101 降到 11。

### 代码产出
- `proof/cpp2lean_v2/passes/pass8_codegen.py` ~900 行：Pass 8 核心
  - emit_type / emit_lit / emit_expr / emit_stmt / emit_param 基础 emit
  - emit_cfg + emit_merge_lambda + emit_bb_inline + emit_terminator + emit_jump_to CFG 流程
  - emit_mirfunc + codegen_pass + codegen_corpus 顶层 emit
  - EmitCtx 含 caller_instance / func_instances / merge_free_vars
  - `__ctor__<template>` / `__cast__<template>` 替换处理
- `proof/cpp2lean_v2/tests/test_pass8_emit_basics.py` 15 单测
- `proof/cpp2lean_v2/tests/test_pass8_cfg.py` 4 单测
- `proof/cpp2lean_v2/tests/smoke_pass8_full.py` + `build_pass8_corpus.py`
  全 67 函数烟测 + corpus 聚合
- `proof/lean/CLPoly/Model.lean` 拷贝 v1 + 多次扩展（abbrev、instance、stub）
- `proof/lean/CLPoly/Generated/Corpus.lean` 6814 行聚合输出（mutual block）

### 上游 Pass 修复（B1 暴露后回流）
- Pass 3 lambda_lift: lifted lambda 命名加 instance_suffix；ret_ty 从
  ReturnStmt.value/Cast.target_ty 推断 fallback
- Pass 4 iter_recognize: 同上 instance_suffix；filterMap pred 从 NamedType
  容器查表推断 elem_ty
- Pass 6 ssa_build: RangeFor cont_var ty 从 s.container.ty 推断；
  RangeFor idx 改用 Nat 类型（与 Lean Array.size 一致）；
  `__write__` 改 `Array.set!` / `_with` 语义化（auto& 写回 + record-update）
- Pass 7 loop_lower: instance_suffix；ret tuple 字段名改 Lean 投影路径
  (`1` / `2.1` / `2.2.1`)
- cast_table.py: int32→nat 用 .toNatClampNeg；nat→int32 / uint64→int32
  用 .toUInt32.toInt32 链；int32→uint64 用 .toInt64.toUInt64
- constructor_map.py: Zp:2 模板 Zp.mk → Zp.ofInt；vector:2/3 size 加 .toNat
- class_map.py: Array.mkArray → Array.replicate

### 提交序列（10 commits）
| 时间 | commit | lake errors | sorry | 关键修复 |
|------|--------|-------------|-------|---------|
| 04-30 | eb3d2c8 | (未跑) | 495 | S1+S2+S3 主体（emit_type/expr/stmt + cfg + corpus） |
| 04-30 | bd255b1 | 8 | 181 | __ctor__/Model 基础 + 第一个 lake build 通过 |
| 05-01 | 3107e2a | 8 | 177 | mutual + 命名冲突 |
| 05-01 | 7b3bd2d | 55 | 64 | 模板实例 + cast 链 + RangeFor 类型（去 sorry 暴露真错） |
| 05-01 | 3c27db3 | 42 | 64 | RangeFor idx Nat |
| 05-01 | 7b6c640 | 28 | 112 | merge BB free vars + Array.set! 语义化 |
| 05-01 | c6e8ef3 | 21 | 112 | Array.set! idx Nat + auto& 写回 set! |
| 05-02 | (一并) | 14 | 91 | Coe + HMul/HPow + Pass 3 ret_ty fallback |
| 05-02 | cdf2bcd | 14 | 91 | (commit 错位) |
| 05-02 | a1b3ffb | **11** | 91 | UMonomial.deg→Nat + HAdd 实参类型化 |

## 为什么做

按用户既定 (a+b) 计划：a = Pass 7 + Pass 8 + Lean Impl/* 桩 + lake build 端到端。
Pass 7 已完成（第 12 轮 1-对-1 审视），Pass 8 是端到端关键最后一环。

策略选择 B → C：先 happy-path 验证可行性 → 全量推 + sorry 降级 → 错误轨迹观察。

后期决定**直接全量推**（用户挑战："小子集不会触发 TupleType 等代码路径，给假信心"），事实印证：单文件 Corpus.lean 聚合 67 函数后才暴露真实跨函数协议问题（lifted lambda 命名、模板实例化、record-update 语义、merge BB free var capture）——这些 happy-path 不会触发。

## 关键决策及其理由

### 1. 单文件 mutual block 聚合，不分文件
- 调研 Agent B 报告：仅 2 个 SCC（factorize 主链 size=8 + mtshl_wmds size=2），分文件价值低
- 单文件全 mutual 让所有 forward-ref 一次性解决，省略了精细调用图分析
- **理由**：B1 阶段验证可行性优先，分块策略留给 B2/C 阶段

### 2. Sorry 降级 vs 上游修复 二选一
- 暴露真错时选择：①用 sorry/placeholder 让 codegen 通过，②回流上游修
- 大多数情况选 ②（避免 silent semantic loss）
- 例外：UnknownType 残留、RefType 残留、Lambda residual 用 sorry/comment 容错
- **理由**：sorry 仅作"未实现"信号；上游真 bug 修上游

### 3. Merge BB free vars 作为 lambda 隐式参数
- 关键 silent miss：merge BB lambda 在 entry stmts 之前 emit，body 引用
  caller-scope vars (如 f_star_2, result_2) Lean 词法作用域看不到
- 选择：emit_merge_lambda 自己分析 free vars + emit_jump_to 调用时传入
- **替代方案放弃**：emit 顺序重排（stmts 前置 + terminator 后置）— 试过引入
  6800→9512 行膨胀 + 错误反增到 48
- **理由**：free var 分析是 Pass 8 内部能力，不需改上游 Pass

### 4. Lean Coe instances + Array.set! idx 自动 .toNat
- Pass 5 cast 链不全 → Lean 端用 Coe 兜底（Int32→UInt64/Int64 等）
- Array.set! 期 Nat → Pass 6 emit 时自动 wrap idx 用 Cast(_, _, Nat)
- **理由**：上游修每条 cast 散点工作量大；Lean Coe + Pass 8 自动 cast 是
  低成本通用方案

### 5. UMonomial.deg : UInt64 → Nat
- C++ side `size_t` 在多种调用上下文下 cast 到 Int32 / Int64 / UInt64 不一致
- Lean Nat 同时支持 .toUInt64 / .toInt32 / .toInt64 等链
- 兼容多个调用点，消除 .toInt32 invalid field 错误
- **替代方案放弃**：UMonomial.deg : Int64 — 引入 SparsePolyZp.derivative 等
  Mul Nat × Int64 不匹配的链式问题

## 遇到的问题与解决方式

### 问题 1：第一次去掉 emit_lit 类型标注后错误反增
- 现象：emit_lit 不带 `(value : type)` 标注，让 Lean 推断 → HAdd Int64
  Nat 类等错误激增
- 根因：BinOp 中 lhs 类型已确定，rhs Lit 不带标注让 Lean 推为 Nat（默认）
- 解决：恢复带标注；BinOp/LetStmt 加 Lit 上下文对齐（rhs.ty 与 lhs.ty 不
  一致时重标 Lit.ty）

### 问题 2：mutual block forward-ref 失败
- 现象：聚合后 caller 调用 callee 仍 unknown
- 根因：mutual block 内某早期 def 的 parse error（空参数名 `( : ZZ)`）
  破坏了整个块解析
- 解决：emit_param 给空名兜底 `_anon_<n>`

### 问题 3：HMul/HAdd SparsePolyZZ instance 不识别
- 现象：lake build 报 HAppend SparsePolyZZ 不存在
- 根因：abbrev 不 reducible 在 instance resolver 中
- 解决：instance 用具体 `Array (UMonomial × Int)` 类型，避开 abbrev 透明度

### 问题 4：Monomial 类型 RangeFor 错误
- 试过加 `abbrev Monomial := MvMonomial` → 反而引入 5 个新错（`StdMap.get!`
  / `Array.insert` 等被解析触发新错）
- 根因：abbrev 加入会改变 Lean 类型推断的解析路径
- 解决：留 known issue，未修

### 问题 5：cross-instance 命名冲突 silent
- 现象：factorize 三实例（lex/grlex/grevlex）的 lifted lambda 重名
  `_lambda_factorize_1_ir`，单文件聚合时定义两次
- 根因：Pass 3 / 4 / 7 用 `func.base_name` 作 host name，三实例 base_name
  都是 "factorize"
- 解决：所有上游 host = `base_name + "_" + instance_suffix`

## 度量

- **耗时**：~14 小时（含调研 / S1-S4 实现 / 7 轮迭代修错 / commit 写 message）
- **迭代**：49 轮编译-修复循环（每轮 = 一次 lake env lean → 修 → 重 build）
- **Lean 新增/修改行数**：
  - Pass 8 codegen 新增：~900 行 (pass8_codegen.py)
  - Pass 8 测试新增：~400 行 (test_pass8_*.py + smoke + build)
  - Pass 3/4/6/7 修改：~80 行（各路 instance_suffix + ret_ty 推断 + Array.set! 等）
  - cast_table / constructor_map / class_map 修改：~20 行
  - Lean Model.lean 扩展：~70 行（abbrev + instance + stub）
  - 总计 ~1470 行新增/修改
- **对应 C++ 行数**：CLPoly 67 函数 ~6800 行（Lean 输出 6814 行 ≈ 1:1）
- **lake build 错误轨迹**：101 → 8 → 8 → 55 → 42 → 28 → 21 → 14 → 11
  （sorry 数 495 → 91 同步）
- **放弃的方案**：
  1. emit 顺序重排（前置非 merge BB stmts + 后置 terminator）— 引入代码膨胀
     + 错误反增；改用 merge BB free vars 隐式参数
  2. UMonomial.deg : Int64 — 与 SparsePolyZp.derivative 等冲突；改 Nat
  3. 加 abbrev Monomial = MvMonomial — 引入 5 新错；留为 known issue
  4. 各 cast 模板用一个统一 chain `({x}.toInt).toInt32` — UInt64.toInt 不存在；
     改用 cast_table 按 src 类型分支

## 涉及的文件列表

### 新增
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/cpp2lean_v2/tests/test_pass8_emit_basics.py`
- `proof/cpp2lean_v2/tests/test_pass8_cfg.py`
- `proof/cpp2lean_v2/tests/smoke_pass8_full.py`
- `proof/cpp2lean_v2/tests/build_pass8_corpus.py`
- `proof/lean/CLPoly/Model.lean`（拷贝自 v1 clpoly_model.lean，扩展）
- `proof/lean/CLPoly/Generated/Corpus.lean`
- `docs/design/l1-translation-validation/survey/pass8-smoke.md`

### 修改
- `proof/cpp2lean_v2/passes/pass3_lambda_lift.py`（host_name + ret_ty fallback）
- `proof/cpp2lean_v2/passes/pass4_iter_recognize.py`（host_name + elem_ty 推断）
- `proof/cpp2lean_v2/passes/pass6_ssa_build.py`（RangeFor + record-update + auto& 写回）
- `proof/cpp2lean_v2/passes/pass7_loop_lower.py`（host_base + tuple 投影路径）
- `proof/cpp2lean_v2/cast_table.py`（Lean 4 API 对齐）
- `proof/cpp2lean_v2/class_map.py`（Array.replicate / instance_suffix consistency）
- `proof/cpp2lean_v2/constructor_map.py`（Zp / vector 模板）
- `proof/cpp2lean_v2/tests/test_pass6_ssa.py`（__write__ → Array.set! 改测试）
- `proof/cpp2lean_v2/ir_types.py`（assert_mir1_invariant 加 TailCall arity）

## Known issues（B1 baseline）

剩 11 errors（每个独立深修，~30-60 min 成本，无主线 ROI）：
1. h_i_1 SparsePolyZZ vs SparsePolyZp（Pass 1 类型分类）
2. __upoly_random tuple vs 单值（Pass 5 output_param destructure）
3. Array.replicate () expected ZZ/Array Int32（Pass 1 vector(n,T()) 默认值）
4. Array.find? present（Pass 1 STL find 误识别）
5. Monomial RangeFor + GetElem? × 2（Pass 1 typedef → Lean 类型）
6. term_1 / g_8 unknown × 2（Pass 4/6 SSA silent miss）
7. __x.fst/snd : sorry × 2（Pass 4 filt2 elem_ty 推断）

## 下一步

1. **B1 续修**（继续修剩 11，预估 +5-7 commits 至 0 errors）
2. **B2 1-对-1 审视**累积 commits（Lesson A 模式，可能再发现 silent miss）
3. **C 转 B2B 验证**（11 errors 作 baseline，找 1 个最干净函数 lake 单独通过 → 实证 pipeline 可用）
