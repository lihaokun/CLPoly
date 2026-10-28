# cpp2lean v2 技术债清单

集中跟踪翻译器实施过程中暴露但**故意推迟**的设计 issue，按"触发条件 / 修复成本 / 影响范围"分类。Stage 3 bring-up 时可对账。

每条记录格式：
```
## TD-N: 标题
- **状态**：open / in-progress / resolved (commit)
- **触发**：什么场景下会让 lake build 失败 / 证明卡住
- **影响范围**：哪些函数、多少 sites
- **根因**：技术问题描述
- **修复方案**：具体如何修
- **修复成本**：行数 + 时间估计
- **暴露时间**：发现的 commit / 阶段
```

---

## TD-1: begin/end → toList 简化非 sentinel 模型

- **状态**：**resolved**（commit `cpp2lean-v2-iterator-model-redesign` CF-1/CF-2/CF-3 完整阶段）
- **触发**：Stage 3 lake build 暴露 STL 其他 iter 惯用法（`std::distance` / `std::remove_if` / `std::lower_bound` / 嵌套 begin+offset 复杂运算）时，Pass 5 现有 Pattern A/B 不够用
- **影响范围**：~10+ 函数（compact-erase 2、parallel-iter 1、sort/iota 7+、find lookup 4）
- **根因**：CLASS_MAP `vec.end()` → `Array.toList` 把 STL 末尾迭代器强制映射为 List，丢失"sentinel"语义；Pass 4 的 filter-loop 模式仅识别 pure filter
- **修复**：迭代器整体 desugar 为索引模型——
  - **CF-1** Pass 4 `_match_filter_loop_A` 扩展接受 mutator stmts，emit `Array.filterMap`（修 compact-erase mutate-then-filter）
  - **CF-2** Pass 6 `_collect_seq_iter_origins` + `_rewrite_seq_iter_in_stmts` 把 `<X>.toList(c)` 来源的 iterator 改写成 `__iter_idx_X : Int64`；`__deref__(it).field` → `c[idx].field`；`it != toList(c)` → `idx < c.size`；`Iterator.advance(it)` → `idx + 1`
  - **CF-3** Pass 5 `sort/iota(toList(c), toList(c), ...)` → `c := Array.sort/range_init(c, ...)`
- **修复时间**：2026-04-29

## TD-2: Lean Impl/* 工程层缺失（70+ 函数）

- **状态**：open（Stage 2 验收前必须先建）
- **触发**：Pass 8 codegen 生成 .lean 文件后，lake build 立刻 unbound identifier
- **影响范围**：CLASS_MAP / CONSTRUCTOR_MAP / FUNC_MAP 引用的 ~80 个目标 Lean 函数，其中 `proof/lean/CLPoly/Impl/Types.lean` 当前仅有 5 个方法 + 4 个 Zp 算术原语
- **根因**：v2 设计阶段集中在翻译器代码 + 数据表配置；Lean 端只做 L2/L3 的 Mathlib 数学模型，没建 L1 工程层
- **缺失清单**（代表性，详见 Pass 5 第二轮审视报告）：
  - 多项式核：`MvPolyZp/MvPolyZZ/SparsePolyZp/SparsePolyZZ.{normalization,front!,back!,size_u64,comp,empty,toList,push}` 全套
  - 单项式：`MvMonomial/UMonomial.{empty,mk,deg,size_u64,...}`、`Variable.{mk,default,find,insert,...}`
  - 数论：`Zp.{ofInt,inv,div,mk,toBool}`、`ZZ.{fdiv_q,fdiv_r,invert,sizeinbase,toBool}`
  - 容器：`StdMap.{empty,get!,insert,find,erase,size,isEmpty}`、`Iterator.{advance,deref!,fromList}`
  - RNG：`Rng.{mk,default,next}`、`UniformIntDist.mk`
  - 工具：`polynomial_GCD`、`pair_vec_div`、`leadcoeff`、`degree`、`is_number`、`get_variables` 等 15+ 函数
  - 算法主体：`__factor_squarefree_primitive_ZZ`、`__upoly_to_poly`、`__symmetric_mod` 等（67 个 TRANSLATION_SCOPE 内函数自身递归调用 → 同步 stub）
- **修复方案**：建立 `proof/lean/CLPoly/Impl/Vec.lean`、`Impl/StdMap.lean`、`Impl/MvPoly.lean`、`Impl/Wang.lean` 等模块，每个目标先放 `def X : … := sorry` stub，Stage 3 bring-up 同步填充实际定义
- **修复成本**：~500-1000 行 Lean stub + ~1500-3000 行 Lean 实现（按 Stage 3 批次推进）
- **暴露时间**：Pass 5 第二轮审视（commit `ee8e5cb`）

## TD-3: Pass 1 typectx 二次推断不完整

- **状态**：open（轻微影响 Pass 5 数据表查表精度，已通过 P0-9 让 Lean 推承担一部分）
- **触发**：嵌套 Call（如 `f.data().empty()`）的返回类型未填到 Call.ty，Pass 5 receiver 类型推断 None，依赖 Lean typeclass 兜底
- **影响范围**：之前 Pass 5 第一轮审视有 ~373 个 None receiver gap，P0-9 让 Lean 推消化大部分；剩余仍是嵌套 Call 返回类型链
- **根因**：Pass 1 `_propagate_var_types` 只处理 LetStmt/RangeFor/参数的 Var，不会推断 Call 的返回类型
- **修复方案**：扩展 `_propagate_var_types`：对 `Call(callee=str, args=...)`，若 callee 在 CLASS_MAP / FUNC_MAP / CONSTRUCTOR_MAP 中能查到返回类型，填到 Call.ty。需要数据表新增 `ret_type` 字段或硬编码常见 method 返回类型表
- **修复成本**：~80-120 行 + 数据表 ret_type 字段扩展
- **暴露时间**：Pass 5 第一轮审视（commit `c0b0b6a` 前后），P0-9 短期化解

## TD-4: HenselNode aggregate-init 误编为 ArrayLit

- **状态**：open（Stage 3 第一次涉及该 struct 时触发）
- **触发**：lake build 时 `#[Sparse, Sparse, Sparse, Sparse, 0, 0, 0, 0]` ArrayLit 不能赋给 HenselNode struct
- **影响范围**：`__hensel_tree_build`、`__hensel_tree_build_recursive` 至少 2 处；其他 Aux struct 类似（PrimeSelectionResult、WangLcResult）
- **根因**：Pass 1 解析 `InitListExpr` 一律输出 `ArrayLit`，没区分 RecordType（aggregate struct）vs container init
- **修复方案**：Pass 1 对 `InitListExpr` 的 type 是 `Record` 时输出 `Call(UnresolvedOp("construct_<TypeName>"), elems)`；Pass 5 通过 CONSTRUCTOR_MAP 解析为 struct 构造
- **修复成本**：~30 行 Pass 1 + CONSTRUCTOR_MAP 加 HenselNode/PrimeSelectionResult/WangLcResult 字段构造
- **暴露时间**：Pass 5 第二轮审视（commit `ee8e5cb`）

## TD-5: UB-6 signed 溢出 require 113 个未生成

- **状态**：open（设计文档 `type-system.md §2.4` 已说"L1 必须保守生成"但 Pass 5 实际没做）
- **触发**：Lean L1 精化定理时无法 discharge "C++ signed 溢出 == Lean Int 算术"——证明义务缺一个前提
- **影响范围**：113 个 signed 算术 site，~10 个独立 require 参数（循环不变量自动消化大部分）
- **根因**：Pass 5 `_collect_ub_requires` 只处理除零/越界/空容器三类，UB-6 留 TODO
- **修复方案**：扩展 `_collect_ub_requires` 对 `BinOp("+/-/*", _, _)` 当 ty ∈ {Int8/16/32/64} 时生成 `RequireStmt(Int.noOverflow_op a b, "h_no_ovf")`；type-system.md §2.4 已规约模板
- **修复成本**：~30 行 + 重跑 smoke 看 require 数量爆增情况
- **暴露时间**：Pass 5 设计阶段（hir-design.md §7.6 已写 TODO）

## TD-6: UB-7 unsigned→signed 截断 require 8 个未生成

- **状态**：open（同 TD-5）
- **影响范围**：8 个 cast 站点
- **根因**：CAST_TABLE `ub_kind="fits_int32/fits_int64"` 在 Pass 5 P0-5 修复后已能生成 RequireStmt；本项是检查 P0-5 是否覆盖全部 UB-7 site
- **修复方案**：实测 Pass 5 输出，若漏插则补 `_emit_cast_ub_require` 的覆盖
- **修复成本**：5 行（基本由 P0-5 覆盖，需验证）
- **暴露时间**：Pass 5 设计阶段

## TD-7: 短路 `&&`/`||` 条件 require 错误提取

- **状态**：open
- **触发**：`if (cond1 && arr[i] != 0)` 中 `arr[i]` 的 oob require 被提到 if 之前，破坏 C++ `&&` 短路语义
- **影响范围**：`__hensel_lift` line 50-52 已确认 1 处；其他 `&&`/`||` 含 ArrayAccess 的 if cond 同样问题
- **根因**：`_collect_ub_requires` 在 IfStmt cond 上一律深扫提取 require，不区分逻辑分支
- **修复方案**：identify `&&`/`||` short-circuit；右操作数的 require 只在左操作数为 true（`&&`）/ false（`||`）时才注入，需要嵌入分支 cond
- **修复成本**：~40 行
- **暴露时间**：Pass 5 第二轮审视（agent 报告）

## TD-8: `__upoly_divmod` 的 TRANSLATION_SCOPE_OUTPUT_PARAMS 漏注册

- **状态**：**resolved**（commit `cpp2lean-v2-rangefor-and-callsite-ref-elim` 修正方案 + Pass 2b 落实）
- **触发**：`(expr) __upoly_divmod(q, r, f, g)` 调用方代码不会被改写为 `(q, r) := __upoly_divmod(f, g)`
- **影响范围**：调用 `__upoly_divmod` 的所有 callsite + 同源未注册函数（`poly_convert`、`fdiv_r/q`、`pair_vec_div#4#5`、`std::swap`）
- **根因**：`class_map.py TRANSLATION_SCOPE_OUTPUT_PARAMS` dict 漏注册 + 调用点 ref-elim 整体未实现（mutation-model-design.md §3.2 设计未落实）
- **修复**：新建 Pass 2b `callsite_ref_elim_pass`（`pass2b_callsite_ref_elim.py` 200 行）+ 补 5 个未注册函数 + 添加 arity overload（`pair_vec_div#4` vs `#5`）。修后 ref-out 47 处 sideeff 全清零。
- **修复时间**：2026-04-28

## TD-9: dump_all_hir.py 不展开 LambdaExpr body

- **状态**：open（仅影响人工审视便利度，不影响正确性）
- **触发**：Pass 5+ 的人工审视看 dump 时 LambdaExpr 显示为 `LAMBDA(captures=[], params=N)` 不展开 body
- **修复方案**：`fmt_expr` 的 LambdaExpr 分支递归 `fmt_stmt` 输出 body
- **修复成本**：~10 行
- **暴露时间**：Pass 3 1-对-1 审视开始

## TD-10: STL `std::set/map<T>` 与 `std::vector<T>` 解析的 Pass 1 类型混乱

- **状态**：open
- **触发**：`__extract_monomial_content` 中 `present : std::set<Variable>` 被 Pass 1 解析为 `Array<Variable>`，但 ctor 用 `__ctor__StdMap.empty()` —— 类型与 ctor 不一致
- **影响范围**：所有用 `std::set` 的函数（CLPoly 实测约 2-3 处）
- **根因**：Pass 1 `parse_type` 对 `std::set<T>` 选择 `ArrayType` 还是 `StdMapType` 不一致；CONSTRUCTOR_MAP 的 STL pattern 把 set 走 `_SET_PATTERNS` 映射到 `StdMap.empty`
- **修复方案**：统一 set 表示——要么全 ArrayType + 改 ctor pattern，要么全 StdMapType + 改 CLASS_MAP 引用
- **修复成本**：~20 行 + 验证 set 使用场景
- **暴露时间**：Pass 5 第二轮审视

## TD-11: 模式 3 iterator/StdMap deref assign（`*it = x` / `it->second = x` / `m.get!(k) = v`）

- **状态**：**resolved**（commit `cpp2lean-v2-iterator-deref-assign` 阶段 A+B + CF-1/CF-2 阶段）
- **触发**：StdMap 迭代器 `it->second = ...` 或 `*it = ...` 在 MIR0 表现为 `let __sideeff_X := __write__(__deref__(it).second, ...)`，root 变量（it / map）未 SSA bump
- **影响范围**：9 处（`__extract_monomial_content`、`__hensel_step_linear`、`__select_eval_point`、`__si_theta_array_eval`、`__upoly_mod_coeff`）
- **根因**：B7 record-update 仅识别 ArrayAccess/FieldAccess root 是 Var；不识别 `__deref__(it)` / `StdMap.get!(m, k)` 形态；Pass 4 拒识 mutate-then-filter 体
- **修复**：分四阶段联合：
  - **阶段 A**: Pass 6 `_root_var` 识别 `Call("StdMap.get!", [m, k])` → root = m（修 P-A 模式 2 处）
  - **阶段 B**: Pass 6 `_collect_iter_origins` 扫 `it = m.find(k)` LetStmt → `__deref__(it).second` 改写为 `StdMap.get!(m, k)`（修 P-B find-source 2 处）
  - **CF-1**: Pass 4 mutate-then-filter 改写为 `Array.filterMap`（修 compact-erase 5 处）
  - **CF-2**: Pass 6 sequence iterator 索引化（修 begin-source 形态）
- **修复时间**：2026-04-29

## TD-12: lambda by-ref capture 修复后剩余 `phi_undef_ver0` (Pass 1 默认初始化局部变量)

- **状态**：open（已 workaround：smoke 阈值放宽到 ≤2）
- **触发**：`__mtshl_step_j` 的 `MvPolyZp Gi;` 与 `__wang_core` 的 `MvPolyZZ h;` 是 C++ 默认构造的局部变量，Pass 1 未生成显式 `let X := __default__(T)` → SSA phi 出现 `phi(bb_entry: X, ...)` 中 X 为 ver=0 裸名
- **影响范围**：2 处（Gi、h）；其它默认初始化局部变量未触发是因为后续都有显式赋值
- **根因**：Pass 1 解析 `VarDecl` 不带 initializer 时仅生成 `LetStmt(var=X, value=UnknownExpr("decomp_uninit"))` 或省略，没显式 emit `__default__(T)` ctor
- **修复方案**：Pass 1 检测无 initializer 的 VarDecl，按 ty 在 CONSTRUCTOR_MAP 查 `__default__` ctor，emit `LetStmt(var=X, value=Call("<ctor>"))`
- **修复成本**：~40 行 + CONSTRUCTOR_MAP 默认 ctor 字段补全
- **暴露时间**：Pass 6 lambda by-ref 修复第三/五轮审视

## TD-13: Pass 1 lambda const-ref 参数误识别为 capture

- **状态**：open（已 workaround：Pass 3 `n_caps` qual_type 编码区分 cap 与 lambda by-ref param）
- **触发**：`auto row_sub = [&](int i, int j, const ZZ& c) { … }` 中 `c` 在 Pass 3 被错误识别为 `[CONST-REF]` capture（应是 lambda 自身参数）
- **影响范围**：~3 个 lifted lambda（`__lll_reduce_3` 的 c、其他 lambda 含 const-ref param 类似）
- **根因**：Pass 1 把 `const T&` 类型的 lambda 形参识别成自由变量？或 lambda body 通过外层变量名引用了 c？需进一步调研
- **修复方案**：调研 Pass 1 lambda param 解析逻辑；确认 const-ref param 的 ParmVarDecl qualType；若解析正确，问题在 Pass 3 自由变量分析的过滤条件
- **修复成本**：~30 行 + 调研
- **暴露时间**：Pass 6 lambda by-ref 修复实施

## TD-14: 模式 9 std::sort/iota 范围迭代器形态错

- **状态**：**resolved**（commit `cpp2lean-v2-iterator-model-redesign` CF-3 阶段 3）
- **触发**：MIR0 中 `let __sideeff_X := sort(Array.toList(result), Array.toList(result), comp)` —— `Array.toList` 是 view 但 `sort` 期望返回新 Array
- **影响范围**：13 处（含 `__factor_Zp`、`__factor_multivar`、`__lll_reduce`、`__zassenhaus_recombine`、`__vanhoeij_recombine`、`__wang_leading_coeff`、`factorize_lex`、`factorize_upoly` 等）
- **根因**：Pass 5 把 STL `std::sort(begin, end, comp)` 映射到 `sort(toList(c), toList(c), comp)`，丢失"in-place 排序"语义；toList 是只读视图
- **修复**：Pass 5 ExprStmt(Call) 处理增加 `_stl_extract_toList_container` + `sort/iota` 改写：
  - `sort(toList(c), toList(c), comp)` → `c := Array.sort(c, comp)`
  - `iota(toList(c), toList(c), v)` → `c := Array.range_init(c, v)`
- **Lean 端依赖**：`Array.sort`（Mathlib 已有）、`Array.range_init`（需 Stage 3 Lean stub ~10 行；写入 TD-2 清单）
- **修复时间**：2026-04-29

## TD-15: lambda 调用作 sort comparator 的间接调用 + by-ref capture

- **状态**：**resolved**（联动 TD-14 修复）
- **触发**：`std::sort(begin, end, mv_next_combination)` 把 lambda 作为 comparator 传递；若 comparator 含 by-ref capture（如 `_lambda___wang_core_3` 的 idx），sort 内部循环修改不传播
- **影响范围**：实测仅作 comparator 的 lifted lambda 都是纯比较器（无 modified caps）。`_lambda___wang_core_3` / `_lambda___zassenhaus_recombine_1` 实际是 `next_combination` 类型——已被 Pass 3b hoist 到 do-while cond 之外（不是 sort comparator）
- **修复**：TD-14 把 `std::sort(toList, toList, comp)` 转 `Array.sort(c, comp)`，comp 作为 LambdaRef 直接传——若 comp 含 ref cap 会编译期报错（lifted lambda 是 pair-return 签名，sort 期望 Bool）；当前 corpus 实测 0 用例
- **修复时间**：2026-04-29（间接修复）

## TD-16: lambda 在 LetStmt/AssignStmt RHS 嵌入表达式深处的调用

- **状态**：partial（Pass 3b 已处理顶层 `LetStmt.value=Call`、`AssignStmt.value=Call`；表达式深处嵌套 Call 未处理）
- **触发**：`let x := f(lambda(idx, n))` —— lambda Call 不是顶层而是嵌入参数
- **影响范围**：当前 67-函数 corpus **未观测**（grep 0 命中），dormant
- **根因**：Pass 3b `_rewrite_stmt` 的 LetStmt/AssignStmt 分支只检查 RHS 顶层 Call，不递归子表达式
- **修复方案**：扩展 `_rewrite_stmt` 递归 walk 表达式找 lambda Call → hoist 到 stmt 之前
- **修复成本**：~50 行
- **暴露时间**：Pass 6 lambda by-ref 修复第五轮 agent 审视提示

---

## 跨项目教训（统一记入 memory）

### 教训 A：每个 Pass 完成必须走完"dump 审视循环"

**根因**：Pass 5 第一轮 release 后跑出"24 → 0 gap、67/67 OK、强化 invariant 全过"——但第二轮 agent 深度审视暴露 11 个 P0 真 bug。"invariant 通过 + smoke 全绿"是必要不充分。

**预防**：每完成一个 Pass，必须：
1. 跑全函数 dump 到 `/tmp/hirN_dump/`
2. **抽样 5+ 个有代表性的函数**（含最复杂、最简单、最异常的）手工读 dump
3. **派 agent 做语义模拟**：对照 C++ 源，逐函数比对"翻译后的 Lean 形态是否真等价"
4. grep 残留模式（`<method>.` / `<operator` / `__cast__` / `<callable>` / tuple-callee 等）
5. 任何"看起来对但没验证过"的数据表映射要 cross-check 目标 Lean 函数存在性

**反例**：Pass 5 第一轮我跳过了 (3) 语义模拟和 (5) cross-check，导致 NoOp delegate / Iterator++ 硬编码 / lambda apply alias 等 bug 全部潜伏。

### 教训 B：技术债必须"两步落实"——commit + 设计文档

**根因**：Pass 5 P2 标"留 Stage 3 处理"——但只在 commit message 写了，设计文档没改。这本身是另一种"承诺了不落实"。

**预防**：任何"留作 P2 / TODO / 已知 limitation" 的决定，commit 之外必须：
1. 加入 `docs/design/l1-translation-validation/tech-debt.md`（本文件）
2. 若涉及多 Pass，更新对应 §（hir-design.md / type-system.md / mir-design.md）末尾的"已知 limitations"
3. 若影响 Stage 3 bring-up 顺序，更新 `translator-v2-todo.md` Stage 3 段
