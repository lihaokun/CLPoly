# 修正方案：cpp2lean v2 iterator/StdMap deref assign 写回

> 状态：草稿（待用户确认）
> 对应 workflow.md §5.1 修正方案文档；tech-debt.md TD-11
> 触发：Pass 6 P0-1/P0-2 修复后残余 `sideeff_write = 9`，2026-04-29
> 前置：[rangefor-and-callsite-ref-elim](./cpp2lean-v2-rangefor-and-callsite-ref-elim.md) + [lambda-by-ref-capture](./cpp2lean-v2-lambda-by-ref-capture.md) 已完成

## 摘要

Pass 6 修了 `vec[i] = e` (B7) + `for (auto&)` (P0-1) + `f(out)` (P0-2) + lambda by-ref 后，剩余 9 处 `__sideeff_X := __write__(...)` 集中在 5 个函数（`__extract_monomial_content`、`__hensel_step_linear`、`__select_eval_point`、`__si_theta_array_eval`、`__upoly_mod_coeff`）。它们是**iterator/StdMap 解引用赋值**——C++ 通过 `it->second = v` 或 `m[k] = v` 写入容器，B7 record-update 链不能识别 root。

修复后 Pass 6 sideeff_write 应降至 0，**phase 6 真正 clean**。

---

## 第一部分：复现与定位

按 grep 实测分类，9 处分 3 种模式：

### 模式 P-A：`m[k] = v`（StdMap subscript assign，2 处）

```cpp
// __select_eval_point
alpha[var] = ZZ(val);

// __si_theta_array_eval
acc[e1_arr[t]] = contrib;
```

MIR0：
```
let __sideeff_28_3_1 := __write__(StdMap.get!(alpha_1, vars_2[i_5]), __ctor__(...)(val_1))
let __sideeff_40_1_1 := __write__(StdMap.get!(acc_1, e1_arr_2[t_5]), contrib_1)
```

期望（修复后）：
```
let alpha_2 := StdMap.set alpha_1 vars_2[i_5] (__ctor__(...)(val_1))
let acc_2   := StdMap.set acc_1 e1_arr_2[t_5] contrib_1
```

### 模式 P-B：`it->second = v`（map iterator deref，5 处）

```cpp
// __extract_monomial_content
auto it = min_deg.find(var);
if (it != min_deg.end()) it->second = std::min(it->second, deg);

// __upoly_mod_coeff / __hensel_step_linear
for (auto it = f.begin(); it != f.end(); ++it)
    fdiv_r(it->second, it->second, m);

// __si_theta_array_eval
auto it = acc.find(k);
if (it != acc.end()) it->second = it->second + contrib;
```

MIR0：
```
let __sideeff_X := __write__(__deref__(it_N).second, <expr>)
```

期望（修复后）：
- `it->second = v` 等价于 `m[k] = v` 其中 `it` 起源于 `m.find(k)` 或 `m.begin()` + 遍历位置 → 需要 SSA bump m
- 对 `it = m.find(k)` 形态：`it->second = v` → `m := StdMap.set m k v`
- 对 `it = m.begin() + i` 形态：复杂（需识别遍历位置 → 转 fold/map）

### 模式 P-C：`*out = move(*it)`（compact-erase 写指针，2 处）

```cpp
// __upoly_mod_coeff（典型 compact-erase 模式）
auto out = it;
for (; it != f.end(); ++it) {
    it->second = fdiv_r(it->second, it->second, m);
    if (it->second != 0) {
        if (out != it) *out = move(*it);
        ++out;
    }
}
f.erase(out, f.end());
```

MIR0：
```
let __sideeff_8_0_1 := __write__(Iterator.deref!(out_2), move(Iterator.deref!(it_2)))
```

期望（修复后）：整个循环应被 Pass 4 重写为函数式形态：
```
f := f.filterMap (fun (k, v) =>
    let new_v := fdiv_r v v m
    if new_v != 0 then some (k, new_v) else none
)
```

---

## 第二部分：根因分析

### 模式 P-A：`_root_var` 不识别 `StdMap.get!(m, k)` 形态

`pass6_ssa_build.py:_root_var` 只识别 Var/ArrayAccess/FieldAccess 链。`Call("StdMap.get!", [m, k])` 是函数调用形态（Pass 5 把 C++ `m[k]` 解析成此形）。无法回溯到 `m`。

### 模式 P-B：iterator 与 source map 的关系丢失

C++ `auto it = m.find(k)` 把 iterator 绑定到 map+key。但 Pass 5 解析后只看到 `Call("__deref__", [it])` —— 不知道 `it` 来自哪个 map。

要正确建模，需要**iterator dataflow tracking**：从 iterator 创建点（`m.find(k)` / `m.begin()` 等）记下 `(source_map, position_index)`，写回时反查。

### 模式 P-C：compact-erase 整个循环需重写

C++ 的 `*out = move(*it); ++out;` 是 erase-remove 惯用法。Pass 4 已有 `compact_erase` filter-loop 识别（详见 `pass4_iter_recognize.py:_match_filter_loop`），但**不处理 body 内对 `it->X` 的修改**（即"先变换再过滤"组合模式）。

### 为什么 Pass 4 漏

Pass 4 现有 `compact_erase` 模式要求 body 是"纯 filter"（只判断 `if (cond) write`）；本案的 body 是"先 mutate `it->second`，再判断，再 write" —— 多了一个 mutate 步骤，Pass 4 不识别。

---

## 第三部分：修复方案

按修复成本从小到大分 3 阶段：

### 阶段 A：模式 P-A 修复（2 处，~25 行）

**改 Pass 6 `_root_var`** 识别 `Call("StdMap.get!", [m, k])`：

```python
def _root_var(e: ExprIR) -> Var | None:
    if isinstance(e, Var):
        return e
    if isinstance(e, ArrayAccess):
        return _root_var(e.arr)
    if isinstance(e, FieldAccess):
        return _root_var(e.obj)
    # P-A 修复：StdMap.get!(m, k) 视为 map root 访问
    if isinstance(e, Call) and e.callee == "StdMap.get!" and len(e.args) == 2:
        return _root_var(e.args[0])
    return None
```

Pass 6 现有 B7 record-update 链自动接管：
```
m := __write__(StdMap.get!(m, k), v)   ← rename 阶段 m bump
```

Pass 8 codegen（未来）把 `__write__(StdMap.get!(m, k), v)` 模式识别为 `StdMap.set m k v`。

### 阶段 B：模式 P-B 修复（5 处，~80 行）

**Pass 6 引入 iterator dataflow（局部）**：

1. 在 `_collect_def_blocks_and_types` 阶段建 `iter_source_map`：扫所有 `LetStmt(var=it, value=Call("StdMap.find"|"StdMap.begin"|...|, [m, ...]))`，记录 `it.name → m_var.name`。
2. 扩展 `_root_var`：识别 `Call("__deref__", [Var(it)])` → 查 iter_source_map → 返回 source map Var。
3. Pass 6 B7 链生成 `m := StdMap.set m <inferred_key> v`。问题：原 C++ 用 `it->second = v`，key 来自 `it->first`。我们需把 `__deref__(it).second = v` 改写成 `m := StdMap.set m (__deref__(it).first) v`，或更简单地保留 `__write__(__deref__(it).second, v)` 形式但 bump m。

考虑到 Pass 8 codegen 复杂度，简化版：仅对**find 来源的 it**做（`it = m.find(k)`，key 已知）：
```python
# 在 LetStmt(var=it, value=Call("StdMap.find", [m, k])) 处记录 it→(m, k)
# 后续 __deref__(it).second = v 改为 m := StdMap.set m k v
```

`m.begin()` 来源的 it（遍历）暂留 P2。

### 阶段 C：模式 P-C 重写（2 处，~120 行）

**Pass 4 扩展 compact-erase 识别**：现有 `_match_filter_loop` 模式假定 body 只有 if-write；扩展为允许 body 头部的 `it->X = expr` 修改：

```cpp
for (it = ...) {
    it->X = expr;          // Pass 4 新增：识别 mutator stmts
    if (cond) *out = move(*it), ++out;
}
f.erase(out, f.end());
```

→ 重写为：
```
f := f.filterMap (fun elem => 
    let elem_mut = { elem with X = expr_substituted }
    if cond_substituted then some elem_mut else none)
```

涉及：
- 模式识别（Pass 4）
- 表达式替换（it->X → elem.X）
- 类型推导（filterMap 的 fun 签名）

依赖 Pass 4 的 filter-loop 框架，工程量较大。

---

## 第四部分：实施 + 验证

### 实施顺序

1. **阶段 A**（~25 行）：Pass 6 `_root_var` 识别 `StdMap.get!`。
2. **回归 + smoke**：sideeff_write 应从 9 → 7。
3. **阶段 B**（~80 行）：iterator dataflow（仅 `m.find(k)` 来源）。
4. **回归 + smoke**：sideeff_write 应从 7 → ~3（剩 `m.begin()` + compact-erase 共 3 处）。
5. **阶段 C**（~120 行，可分独立 ticket）：Pass 4 mutate-then-filter 模式。
6. **回归**：sideeff_write 应 → 0。
7. **agent 第六轮审视**确认。

### 验证标准

- `bash tests/run_all_smoke.sh` 全过
- `smoke_pass6_full.py` 输出 `sideeff_write = 0`（或留 `m.begin()` 来源 ~2 处作 P2 follow-up）
- `__extract_monomial_content`、`__upoly_mod_coeff`、`__hensel_step_linear` 等 dump 中 `m`/`f` 有 SSA chain 推进
- 单测：3 个新用例覆盖 P-A / P-B / P-C 各一例

### 风险

- **风险 1**：iterator dataflow 在嵌套 loop 中失效（it 跨 loop 重新赋值）。
  - 缓解：仅在直系 LetStmt 紧邻处建立映射；跨复杂控制流时不应用（保守 fallback 到 sideeff）。
- **风险 2**：阶段 C 的 mutate-then-filter 模式与现有 filter-loop 识别冲突。
  - 缓解：先做识别 + dry-run（生成警告而非改写），确认覆盖正确后再 hard-replace。
- **风险 3**：`StdMap.set` 的 Lean 语义未确认（CONSTRUCTOR_MAP / CLASS_MAP 是否注册）。
  - 缓解：Pass 8 codegen 阶段才暴露；本方案先生成 `__write__(StdMap.get!(...), ...)` 形态，让 codegen 处理。

### 推荐先做范围

**最小 commit**：仅做阶段 A + B（共 ~105 行）。预计 sideeff_write 9 → 3。剩余 3 处（m.begin() iter-loop + compact-erase）作 TD-11 follow-up 留 Pass 4 重写。

**理由**：
1. 阶段 C 涉及 Pass 4 重写，工程量大且与现有 filter-loop 框架耦合，需独立审视；
2. 阶段 A+B 已能修 7/9（78%）；
3. 剩余 3 处都在 `__upoly_mod_coeff` / `__hensel_step_linear` 的 compact-erase 循环，**不阻塞** Pass 7 推进（Pass 7 处理结构化循环 → tail-call MIR₁ 不依赖这些 sideeff 是否清零）。

---

## 第五部分：关联与依赖

### 与 P0-2 / P0-1 / lambda by-ref 的关系

本方案是 ref-elim 体系的最后一类——**iterator/map 写回**。前述：
- P0-2：调用点 ref out（`f(out)` → `out := f(out)`）
- P0-1：ranged-for `auto&`（`for (auto& x : c) x.f = v` → `c[idx] := __write__(c[idx], x_modified)` + exit `c := cont_post_loop`）
- lambda by-ref：lifted lambda 的 cap mutation 通过 pair-return 传播

这一方案补完最后一类——**iterator deref 写回**。修完后 Pass 6 输出的 MIR0 应**完整反映 C++ mutating semantics**。

### 与 TD-1 begin/end 的关系

TD-1 提到 STL `vec.end()` → `Array.toList` 丢失 sentinel 语义。本方案部分相关——iterator 写回的 codegen 需要 Pass 8 把 `__write__(StdMap.get!(m, k), v)` 映射到 `StdMap.set m k v`，依赖 CLASS_MAP 设计。

### 与 Lean Impl 层（TD-2）的关系

Lean 端需要 `StdMap.set : StdMap K V → K → V → StdMap K V` 函数（functional update）。当前 TD-2 列了 `StdMap.{empty, get!, insert, find, erase, size, isEmpty}`，缺 `set` —— 本方案需要补到 Lean stub 清单。
