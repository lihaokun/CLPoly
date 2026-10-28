# 修正方案：cpp2lean v2 迭代器模型整体重写（TD-1 + TD-11 阶段 C）

> 状态：草稿（待用户确认）
> 对应 workflow.md §5.1 修正方案文档
> 触发：iterator-deref-assign 修正方案 §3 阶段 C 重新审视，2026-04-29
> 前置：[iterator-deref-assign](./cpp2lean-v2-iterator-deref-assign.md) 阶段 A+B 已完成

## 摘要

经 Phase 6 修复链（B7 record-update / P0-1/2 ranged-for / lambda by-ref / TD-11 阶段 A+B）后，剩 5 处 `sideeff_write` 全为 STL **compact-erase 迭代器循环**（`__upoly_mod_coeff`、`__hensel_step_linear` 各 ~3 个 sites）。这暴露了 **TD-1 STL 迭代器整体模型缺失** 的根本问题——当前 Pass 5 把 `vec.begin()/.end()` 映射成 `Array.toList`，丢失"位置 + sentinel"语义，Pass 4 的 filter-loop 识别只能覆盖纯过滤场景。

本方案重新设计**迭代器模型**：从 HIR4 阶段开始，把所有 STL iterator 用法 **desugar 为索引模型 + 函数式 Array 操作**。这一改造同步解决：
- TD-1（begin/end → toList sentinel 丢失）
- TD-11 阶段 C（compact-erase mutate-then-filter）
- TD-14（sort/iota 范围迭代器形态错）
- TD-15（lambda 作 sort comparator）

修复后 Pass 6 sideeff_write 应降至 0，sideeff_other 也大幅减少（sort/iota 转 functional 形态），Phase 6 真正 clean。

## 第一部分：实测分类

### 67-corpus 中 iterator 使用分布（grep + dump 实测）

| 模式 | 函数数 | 子模式 | 当前 sideeff |
|------|--------|--------|---|
| **C1 compact-erase**（mutate-then-filter） | 2 | `for it: it->X=mutate; if cond *out=move(*it),++out` + `f.erase(out,end)` | sideeff_write 5 |
| **C2 parallel iter**（双指针 zip-merge） | 1 (`__upoly_divmod_mod`) | `while (r_it != end \|\| g_it != end)` 两路条件推进 | sideeff_other ~3 |
| **C3 find/lookup**（point query） | 4 | `auto it = m.find(k); if (it != m.end()) ...` | 已被 TD-11 阶段 B 修 |
| **C4 sort/iota**（STL 算法） | 7 | `std::sort(c.begin(), c.end(), comp)` / `std::iota(c.begin(), c.end(), v)` | sideeff_other 13 |
| **C5 range-for** | 多 | `for (auto& x : c)` | 已被 P0-1 修 |

实测涉及函数（HIR4 grep `Iterator|toList|StdMap.find`）共 **14 个**：
- C1: `__upoly_mod_coeff`, `__hensel_step_linear`
- C2: `__upoly_divmod_mod`
- C4: `__factor_Zp`, `__factor_multivar`, `__lll_reduce`, `__zassenhaus_recombine`, `__vanhoeij_recombine`, `factorize_lex`, `factorize_upoly`
- C3 (已修): `__extract_monomial_content`, `__si_theta_array_eval`
- 其它 (含 range-for): `__wang_core`, `__wang_leading_coeff`, etc.

### 当前 Pass 5 + Pass 4 模型

```
Pass 5 CLASS_MAP:
  Array.begin / Array.end / SparsePolyZZ.begin → "Array.toList"  ← 丢位置+sentinel
  Array.erase(c, out_iter, end_iter) → "Array.erase(c, ...)"     ← out_iter 是 toList，不能直接索引
  StdMap.find(m, k) → "StdMap.find(m, k)"                        ← 返回 Iterator（位置不明）

Pass 4 _match_filter_loop_A/B:
  - 仅匹配 pure filter（_is_pure_filter_body 拒绝 mutate body）
  - 输出：Array.filter(c, pred_lambda)
```

## 第二部分：新模型设计

### 核心思路：iterator 全部 desugar 为索引

| C++ 源 | 旧模型（Pass 5 toList） | 新模型（索引 + Array.X） |
|--------|---------|---------|
| `vec.begin()` | `Array.toList(vec)` | `0`（起始索引）|
| `vec.end()` | `Array.toList(vec)`（重复！） | `Array.size(vec)` |
| `*it` | `Iterator.deref!(it)` | `vec[i]` |
| `*it = e` | `Iterator.deref!(it) = e` | `vec := vec.set i e` |
| `it->second = e` | `__deref__(it).second = e` | `vec := vec.set i { vec[i] with snd := e }` (record-update) |
| `++it` | `it = Iterator.advance(it)` | `i := i + 1` |
| `it == vec.end()` | `it == Array.toList(vec)` | `i == vec.size` |
| `vec.erase(out, vec.end())` | `vec = Array.erase(vec, out, ...)` | `vec := vec.take out_idx`（截断到 out 位置）|
| `m.find(k) != m.end()` | 已 TD-11-B 修：`StdMap.contains m k`-equiv | （保留） |
| `m.find(k)->second` | 已 TD-11-B 修：`StdMap.get! m k` | （保留） |
| `for (it=c.begin();it!=c.end();++it)` | classic iter-for-loop | `for i in [0..c.size) do let elem := c[i]; ...` |
| `std::sort(c.begin(),c.end(),comp)` | `sort(toList,toList,comp)` | `c := Array.sort comp c` |
| `std::iota(c.begin(),c.end(),v)` | `iota(toList,toList,v)` | `c := Array.range_init c v` |

### 新模型的不变量

1. **iterator 不出现在 IR**：HIR4 起，所有 iterator 已被消除（Pass 4 desugar）
2. **erase/sort/iota 等 STL 算法已 functional**：`vec := op c args`
3. **位置语义保留**：所有迭代变量是 `Nat` 索引

### 与现有 RangeForStmt 的关系

`RangeForStmt(var=x, var_ty, container, body, is_mutable_ref)` 已在 Pass 6 desugar 为索引循环（`__rangefor_idx_N`/`__rangefor_cont_N`）。新 iterator 模型可直接重用 RangeForStmt 作为目标 IR，避免引入新 stmt 类型。

## 第三部分：实施分阶段

### 阶段 1：C1 compact-erase mutate-then-filter（~120 行）

**目标**：处理 `__upoly_mod_coeff` / `__hensel_step_linear` 共 5 处 sideeff_write。

**步骤**：

1. **扩展 `_is_pure_filter_body`** 接受 body 头部的 mutator stmts：
   - 检查 body 第一组连续 stmts 是 `AssignStmt(target=*it 或 it->X, value=expr)` 形态
   - 后续 stmts 仍是 if-write-advance 模式

2. **扩展 `_build_pred_lambda`** 输出 filter+map 合一 lambda：
   - 把 mutator stmts 改写为 `let elem' := expr_substituted`
   - filter cond 改写后用 elem' 而非原 elem
   - 返回 `if cond then some elem' else none`

3. **新建 `_filter_loop_to_filterMap`**：emit `c := Array.filterMap pred_lambda c`

4. **CLASS_MAP** 注册 `Array.filterMap`（作为 functional 算法）。

5. **Lean Impl** 端 `Array.filterMap` 已是 Mathlib 标准函数，无新 stub。

**输出**：5 处 sideeff_write → 0；新增 `Array.filterMap` 调用 ~5 个。

### 阶段 2：C2 parallel-iter（~150 行）

**目标**：处理 `__upoly_divmod_mod` 双指针 zip-merge。

**步骤**：

1. 识别两路 `let r_it := c1.begin()` + `let g_it := c2.begin()` + `while (r_it != c1.end() || g_it != c2.end())`
2. desugar 为：
   ```
   let r_idx := 0; let g_idx := 0
   while (r_idx < c1.size || g_idx < c2.size) {
       body[r_it ↦ c1[r_idx], g_it ↦ c2[g_idx]]
       (advance handled by body's existing assign)
   }
   ```
3. body 内 `r_it = Iterator.advance(r_it)` 改写为 `r_idx := r_idx + 1`

**风险**：parallel-iter 的 sentinel 检测（`!=end`）在 OR 条件中——desugar 必须保持短路语义。

**输出**：1 处复杂 iter-while → 索引化 WhileStmt；后续 Pass 6 走标准 SSA。

### 阶段 3：C4 sort/iota STL 算法（~80 行）

**目标**：处理 `std::sort(begin, end, comp)` 7 处 + `std::iota(begin, end, v)` 2 处共 ~13 处 sideeff_other。

**步骤**：

1. **Pass 5** 增加 `_match_stl_algorithm` 检测：
   - `Call("sort", [Array.toList(c), Array.toList(c), comp])` → `c := Array.sort comp c`
   - `Call("iota", [Array.toList(c), Array.toList(c), v])` → `c := Array.range_init c.size v`
2. **CLASS_MAP / FUNC_MAP** 注册 `Array.sort` / `Array.range_init`
3. **Lean Impl** 端：`Array.sort` Mathlib 已有；`Array.range_init` 需小 stub（10 行）

**输出**：13 处 sideeff_other → 0；同时 TD-15 lambda comparator + by-ref capture 自然解决（comparator 是纯函数则正常，含 by-ref 则报错）。

### 阶段 4：C3 find/lookup 兼容性确认（~20 行）

**目标**：验证 TD-11 阶段 B 已完整覆盖 find-source 模式；补 begin-iteration 形态。

**步骤**：跑 smoke 看是否还有 `__deref__(it).second` 残留；若有补 `_collect_iter_origins` 覆盖。

## 第四部分：风险与依赖

### 跨阶段依赖

```
阶段 1 → 阶段 2 → 阶段 3 → 阶段 4
    ↓        ↓        ↓
   独立    独立    独立（Lean stub 依赖）
```

各阶段独立可推进；阶段 3 依赖 Lean Impl 端 `Array.sort`（已有）+ `Array.range_init`（新增 stub ~10 行）。

### 风险

1. **mutator 表达式的副作用语义**：阶段 1 把 `fdiv_r(it.second, it.second, m)` 内联进 lambda——如果 fdiv_r 不是纯函数（外部状态依赖），lambda 内副作用未定义。
   - 缓解：fdiv_r 已被 Pass 2b 改写成纯函数返回（`q := fdiv_r(q, q, m)`）；GMP fdiv_r/fdiv_q 数学上是纯函数。
2. **parallel-iter 的索引边界条件**：阶段 2 的 OR 条件 `r_idx<size || g_idx<size` 可能让 body 同时访问越界容器。
   - 缓解：body 内已有 `if (r_idx < c1.size)` 保护；保持原结构。
3. **comparator lambda**：阶段 3 把 sort comparator 视为纯函数——若 comparator 含 by-ref capture（TD-15），现行 Pass 3b 已 ref-elim 改成 pair-return，comparator 期望 Bool 但收到 (Bool, idx) → 类型不匹配。
   - 缓解：阶段 3 在 _match_stl_algorithm 时检查 comparator 不应有 modified caps（has_ref=False），违反则报 gap warning。

### 可分批 commit 计划

| commit | 阶段 | 影响 |
|--------|------|------|
| 1 | 阶段 1 mutate-then-filter | sideeff_write 5→0；__upoly_mod_coeff / __hensel_step_linear OK |
| 2 | 阶段 3 sort/iota | sideeff_other -13；7+ 函数 OK |
| 3 | 阶段 2 parallel-iter | __upoly_divmod_mod OK |
| 4 | 阶段 4 验证 + tech-debt 关闭 | TD-1 / TD-11 / TD-14 / TD-15 全 resolved |

每 commit 可独立 regression。

## 第五部分：成功判据

修复完成时：
- `bash tests/run_all_smoke.sh` 全过
- `smoke_pass6_full.py`：
  - `sideeff_write = 0`
  - `sideeff_other` ≤ 10（仅余真正 discarded return）
  - `phi_undef_ver0 = 0`
- 14 个 iterator 函数 dump 中：
  - 无任何 `Iterator.X` / `toList` 出现
  - C1: `Array.filterMap` 形态
  - C2: 索引化 WhileStmt
  - C4: `Array.sort` / `Array.range_init`
- 单测：mir_types 9 + pass6 12 + pass2b 7 + pass3b 5 + **新阶段 1/2/3 各 3-5 用例** = ~37
- TD-1 / TD-11 / TD-14 / TD-15 全 resolved
- agent 第六轮抽样审视 `__upoly_mod_coeff` / `__upoly_divmod_mod` / `__factor_Zp`：semantic OK

## 第六部分：替代方案对比

| 方案 | 工程量 | 解决问题 | 不解决 |
|------|--------|---------|--------|
| **C-Defer**（当前选项） | 0 | 无 | 5 处 sideeff_write + TD-1/11/14/15 全留 |
| **仅阶段 1**（C-Mini） | ~120 | sideeff_write 5→0 | TD-1 主体 + TD-14/15 |
| **阶段 1+3**（推荐起点） | ~200 | sideeff_write + 大部分 sideeff_other | TD-1 主体（C2 parallel-iter） |
| **完整 C-Full**（4 阶段） | ~370 | 全部 | 无（TD-1/11/14/15 全清） |

## 第七部分：建议执行顺序

**推荐**：分 commit 推进，每 commit 独立可回滚。

1. **本次 commit**：仅做阶段 1（C1 mutate-then-filter）。120 行。Phase 6 sideeff_write = 0，达成"phase 6 真正 clean"。
2. **下次 commit**：阶段 3（C4 sort/iota）。和 Pass 7 并行（sort/iota 是 Pass 4/5 工作，不阻塞 Pass 7）。
3. **可选**：阶段 2（C2 parallel-iter）。仅 1 个函数受影响，可在 Stage 3 bring-up 时再修。

理由：
- 阶段 1 是 phase 6 收尾的最低必要量
- 阶段 3 补完所有 STL 算法，但与 Pass 7 无强依赖
- 阶段 2 是 1 个特殊函数，可推迟

## 附录 A：关联文档

- `docs/design/l1-translation-validation/tech-debt.md` TD-1 / TD-11 / TD-14 / TD-15
- `docs/design/l1-translation-validation/iterators.md`（如有）
- `docs/fixes/cpp2lean-v2-iterator-deref-assign.md`（前置阶段 A+B）
- `proof/cpp2lean_v2/passes/pass4_iter_recognize.py:_match_filter_loop_A/B`
- `proof/cpp2lean_v2/passes/pass5_operator_resolve.py`（CLASS_MAP / FUNC_MAP）
