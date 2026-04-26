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

- **状态**：open（短期 Pattern A 兜底，长期需 B+ 重组 Pass 4）
- **触发**：Stage 3 lake build 暴露 STL 其他 iter 惯用法（`std::distance` / `std::remove_if` / `std::lower_bound` / 嵌套 begin+offset 复杂运算）时，Pass 5 现有 Pattern A/B 不够用
- **影响范围**：当前已识别 7+ sites（compact-erase 3、classic iter-loop 2、vector(begin±i, end±j) 6、erase(begin+j) 3、set::find==set::end 1）；CLPoly 实测覆盖完，但任何新 STL 算法引入会再发
- **根因**：`CLASS_MAP[Array]["end"] = ("method", "Array.toList")` 把 STL 末尾迭代器强制映射为 List，丢失"sentinel"语义。`set/map.find == set.end` 这种"不存在"惯用法只能靠 Pass 5 单独 pattern matching 兜底
- **修复方案 B+**：在 Pass 4 阶段统一识别**所有** begin/end 出现位置，整体重写为对应 Lean 形态（`set::find == end()` → `(s.find? k).isNone` / `for(it=begin; it!=end; ++it)` → `for x in toList c do …` / `vec.erase(begin+j)` → `Array.eraseIdx j` / `vector(begin+i, end)` → `Array.drop i` / 等）。Pass 5 之后**禁止**任何 `<method>.begin/.end` 残留（assert 强制）
- **修复成本**：~150-250 行（新 Pass 4 patterns + Pass 4/5 现有 Pattern A/B 迁移整合 + 回归）
- **暴露时间**：Pass 5 第二轮审视（commit `ee8e5cb`）

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

- **状态**：open（Stage 3 bring-up 时 callsite 改写不会发生）
- **触发**：`(expr) __upoly_divmod(q, r, f, g)` 调用方代码不会被改写为 `(q, r) := __upoly_divmod(f, g)`
- **影响范围**：调用 `__upoly_divmod` 的所有 callsite
- **根因**：`class_map.py:686 TRANSLATION_SCOPE_OUTPUT_PARAMS` dict 漏 `__upoly_divmod: [0, 1]`
- **修复方案**：补一行
- **修复成本**：1 行
- **暴露时间**：Pass 5 第二轮审视

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
