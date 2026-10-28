# 修正方案：cpp2lean v2 lambda by-ref capture 写回

> 状态：草稿（待用户确认）
> 对应 workflow.md §5.1 修正方案文档
> 触发：Pass 6 第三轮 1-对-1 审视（Agent B 语义模拟），2026-04-28
> 上一方案：[rangefor-and-callsite-ref-elim.md](./cpp2lean-v2-rangefor-and-callsite-ref-elim.md)（已实施）

## 摘要

cpp2lean v2 在 Pass 3 lambda_lift 之后，**lifted lambda 的 by-ref capture 不正确传播 mutation 到外层 scope**。Pass 3 已正确：
- 检测 modified captures（body 中 AssignStmt 的 capture）
- 给 lifted lambda 的 cap params 标 `is_ref=True` / `is_const_ref=True`
- 通过 `qual_type` 字符串记录 `modified_captures`

**未做**：
1. **Lifted lambda 签名未改 pair-return**（Pass 2 ref_elim 不跑 lifted lambda）
2. **Call site 未 prepend captures**（pre-existing Pass 3 bug：lifted lambda 有 6 参数，call site 只传 3）
3. **Call site 未 destructure**（同 P0-2 根因，但 Pass 2b 不跑 lifted lambda）

修复后：6 个 lifted lambdas / 8 modified captures / 8 call sites 的 silent miss 全部修。

---

## 第一部分：复现与定位

### 最小复现（`__lll_reduce`）

C++ 源（核心 LLL 算法的 row reduce）：
```cpp
auto row_sub = [&M, &U](int n, int i, int j, ZZ c) {
    for (int k = 0; k < n; ++k) {
        M[i][k] -= c * M[j][k];
        U[i][k] -= c * U[j][k];
    }
};
auto row_swap = [&M, &U](int i, int j) {
    std::swap(M[i], M[j]);
    std::swap(U[i], U[j]);
};

// LLL 主循环中调用：
row_sub(n, k, k-1, q);     // 期望：M、U 被 row reduce
row_swap(k, k-1);          // 期望：M、U 行交换
```

### 当前 HIR3 形态（`/tmp/hir3_dump/__lll_reduce.txt:150-171`）

```
# _lambda___lll_reduce_3
qualType: lambda in __lll_reduce | modified_captures=['M', 'U']

## Params
  [0] M : LLLMatrix [REF]                    ← Pass 3 已标 is_ref=True
  [1] U : LLLMatrix [REF]                    ← Pass 3 已标 is_ref=True
  [2] n : Int32 [CONST-REF]
  [3] i : Int32
  [4] j : Int32
  [5] c : ZZ [CONST-REF]

## Return type: ?[]                          ← 未改 pair-return（应为 (LLLMatrix, LLLMatrix)）

## Body
for (init;cond=(k < n);step) {
  [init] let k : Int32 := 0
  (expr) <operator-=>(<operator[]>(<operator[]>(M, ...), ...), ...)
  (expr) <operator-=>(<operator[]>(<operator[]>(U, ...), ...), ...)
  [step] (expr) (++k)
}
```

### 当前 MIR0 调用点形态（`/tmp/mir0_dump/__lll_reduce.txt`）

```
let __sideeff_24_0_1 := _lambda___lll_reduce_3(k_2, (k_2 - 1), q_1)
let __sideeff_38_0_1 := _lambda___lll_reduce_3(k_2, j_11, q2_1)
let __sideeff_47_1_1 := _lambda___lll_reduce_4(k_2, (k_2 - 1))
```

**3 个 bug 同时显现**：
1. **bug A**：调用只传 3 args（`k_2, (k_2-1), q_1`）但签名要 6 params（`M, U, n, i, j, c`）—— Pass 3 没 prepend captures。
2. **bug B**：返回值放 `__sideeff_*`（discarded），即使 Lambda 修改了 M/U，外层永远看不到。
3. **bug C**：Lambda body 内 `M[i][k] -= ...` 在 Pass 6 里转成 `__sideeff_*=__write__(M[...], ...)`（参见 sideeff_write metric 9 处中的部分），lambda body 内 SSA 不 bump M。

### 影响范围（grep 实测，2026-04-28）

| Lambda | Modified captures | Call sites（MIR0 sideeff） |
|--------|-------------------|---------------------------|
| `_lambda___lll_reduce_3` | M, U | **2** |
| `_lambda___lll_reduce_4` | M, U | **1** |
| `_lambda___mtshl_step_j_1` | Gi | 1 |
| `_lambda___wang_core_1` | h | 4 |
| `_lambda___wang_core_3` | idx | 0（passed as comparator） |
| `_lambda___zassenhaus_recombine_1` | idx | 0（passed as comparator） |

**总计**：6 lifted lambdas，8 modified captures，**8 实际 call sites silent miss**。

`_lambda___wang_core_3` / `_lambda___zassenhaus_recombine_1` 被作为 comparator 传给 `std::sort`（第三方调用，0 直接 call），但本身仍有 `idx` 修改——也需要修，但 mechanism 不同（属于 sort 范畴）。

### 受影响算法

| 函数 | 算法 | 错误后果 |
|------|------|----------|
| **__lll_reduce** | LLL 格基约简（Van Hoeij 因式分解核心） | row reduce / row swap 全失 → LLL 无效 → 因式分解错 |
| `__mtshl_step_j` | MTSHL Newton 迭代步 | Gi 累积失 → MTSHL 收敛错 |
| `__wang_core` | Wang/EEZ 多变量因式分解核心 | h 多项式构造错 |

后果严重：影响整个 L1 翻译器对**因式分解 + LLL 子系统**的正确性。

---

## 第二部分：根因分析

### 设计 vs 实现

| `mutation-model-design.md` §位置 | 设计 | 实施 |
|---|---|---|
| §5 G3 闭包变量完备性 | "后处理 pass 把 capture 加为 lambda params" | ✓ Pass 3 落实（`is_ref` 标记） |
| §5.1 `collect_free_vars` | 自由变量收集 | ✓ Pass 3 `_collect_free_vars` |
| §3.4 函数定义签名改造（pair-return） | "TRANSLATION_SCOPE 中函数应用此规则" | ✗ **lifted lambda 不在 TRANSLATION_SCOPE 中——Pass 2 不跑** |
| §3.2 调用点变换 | "查 OUTPUT_PARAMS 改写 call site" | ✗ **lifted lambda 不在 TRANSLATION_SCOPE_OUTPUT_PARAMS 中——Pass 2b 不查** |

设计上 §3.2/§3.4 应同时覆盖 outer 函数 *和* lifted lambda，但实现把"是否 mutating ref-out"绑定到 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 静态 dict。lifted lambda 是 Pass 3 动态生成的，没有静态注册。

### Pass 3 已知 pre-existing bug（bug A）

Pass 3 `_lift_lambda` (`pass3_lambda_lift.py:280`) 创建了 lifted HIRFunc：
```python
new_params = cap_params + list(lam.params)   # 6 个 params for row_sub
```

但 **call site 重写**（`_rewrite_expr` 处理 Call）只重写 args 的内嵌 lambda，不 prepend captures：
```python
if isinstance(e, Call):
    new_args = [_rewrite_expr(a, ...) for a in e.args]
    return replace(e, args=new_args)
```

C++ `row_sub(k, k-1, q)` → `Call(callee=Var(row_sub), args=[k,k-1,q])` → Pass 3 只递归 args，不变更 callee。Pass 5 后续把 `Var(row_sub)` 解析为 lifted name `_lambda_..._3`（因为 row_sub 是 LetStmt 别名），但仍 3 args。

### 为什么 67 函数烟测当前是绿的

虽然语义错（`__lll_reduce` 输出 M/U 永远是初始值），但：
- Pass 6 invariant 检查"结构性"（CFG / SSA / phi keys）
- smoke metric 检查"残余 sideeff 调用形态"——但 `_lambda_*` 不在 `TRANSLATION_SCOPE_OUTPUT_PARAMS`，sideeff_other 计入但不报错
- 没有端到端语义测试（Lean codegen 比对 C++ 输出尚未做）

这是 Lesson A 第三次成立：**syntactic invariant + smoke + dump 抽样 ≠ semantic correct**——必须用 agent 做语义模拟才能挖出。

---

## 第三部分：修复方案

### 总体策略

**新建 Pass 3b: `lambda_ref_elim_pass`**，在 Pass 3 之后、Pass 4 之前，做三件事：

1. **签名改造**：每个 lifted lambda 中 `is_ref=True` 的 cap params → 函数返回 tuple 含 modified captures（沿用 Pass 2 ref_elim 的 `_tuple_of` + `_rewrite_returns` 逻辑）。
2. **修复 bug A：调用点 prepend captures**：每个 call to `_lambda_<name>(args)` → `_lambda_<name>(<cap_args>, args)`，其中 `<cap_args>` 是从外层作用域的 Var 名查得（与 lifted lambda 签名前 N 个 params 对应）。
3. **修复 bug B/C：调用点 destructure**：仿 Pass 2b，把 `(expr) _lambda_..._N(<caps>, <args>)` 改写为 `let __refret := ...; <cap_i> := __refret.<field_i>`。

### 步骤 P1：检测 cap params 与 modified captures

`HIRFunc` 的 lifted lambda：
- `is_ref=True` 的 params = modified captures（Pass 3 已标）
- `is_const_ref=True` 的 params = read-only captures（Pass 3 已标）
- `is_ref=False, is_const_ref=False` = 原 lambda params

提取 cap names：
```python
def _extract_cap_names(lifted: HIRFunc) -> list[str]:
    """从 Pass 3 lifted lambda 提取 capture 参数名（leading is_ref/is_const_ref）。"""
    names = []
    for p in lifted.params:
        if p.is_ref or p.is_const_ref:
            names.append(p.name)
        else:
            break  # 一旦遇到非-cap，后续都是 lam.params
    return names
```

### 步骤 P2：对 lifted lambda 应用 ref_elim_pass

复用现有 `pass2_ref_elim.ref_elim_pass`：
```python
new_aux = []
for lifted in func.aux_lambdas:
    new_lifted = ref_elim_pass(lifted)   # 复用 Pass 2 逻辑
    new_aux.append(new_lifted)
```

`ref_elim_pass` 已正确处理：
- 计算 `ref_params = [p for p in func.params if p.is_ref]`
- 改 `ret_ty = _tuple_of([orig_ret] + ref_types)`（非 void 时）或 `_tuple_of(ref_types)`（void）
- 改写 ReturnStmt 包装 tuple
- 末尾追加 ReturnStmt（void）
- 清 is_ref / is_const_ref 标记

注意：ref_elim_pass 调用前需要保存 cap names（清标记后无法再识别），因为 P3 步要在 outer body 重写时用。

### 步骤 P3：在 outer body 重写 call sites

对每个 lifted lambda L 与其 call sites（在 outer body）：

**重写规则**：
```
原：(expr) Call(callee="_lambda_<name>", args=[orig_args])
后：
  let __refret_N := Call(callee="_lambda_<name>", args=[<cap_args>, ...orig_args])
  <cap_0> := __refret_N.<field_0>     # AssignStmt for each modified cap
  <cap_1> := __refret_N.<field_1>
  ...
```

其中：
- `<cap_args>` = `[Var(name=cap_name) for cap_name in saved_cap_names_for_L]`
- modified caps 决定 destructure 的字段。原 `lifted` 的 modified caps 来自 `is_ref=True` params（清标记前）。

**关键问题**：call sites 通过什么形式标识？在 HIR3 里是 `Call(callee=Var("row_sub"), ...)`（row_sub 是绑定到 LambdaRef 的 LetStmt 别名），不是直接的 `Call(callee="_lambda_..._3", ...)`。

**解决**：先建立 alias map：在 outer body 中找 `let row_sub := Var("_lambda_..._3")` 形式的 LetStmt，建立 `{"row_sub": "_lambda_..._3"}`。然后改写 Call 时同时认 alias。

### 步骤 P4：单元测试 + smoke metric

新增 `test_pass3b_lambda_ref_elim.py`（5 用例）：
- T1: 单 modified capture（`[&x]` lambda）→ `let (x', ret') := lambda(x, args)`
- T2: 双 modified captures（`[&M, &U]`）→ pair-return + destructure
- T3: 混合 modified + read-only capture → 仅 modified 在返回值
- T4: 0 modified captures（read-only `[&]`）→ passthrough
- T5: 67 函数烟测（含 6 个有 by-ref capture 的函数）

新增 smoke metric：`lambda_ref_callsite_unwritten` = call to lifted lambda with `is_ref` params 但没有 destructure 赋回 captures 的次数。期望 0。

### 步骤 P5：文档同步

- `mutation-model-design.md` §5（G3 闭包）补充："本设计在实现层由 Pass 3 + Pass 3b 联合落实——Pass 3 检测 + 标记，Pass 3b 改写签名 + call site"。
- `tech-debt.md` 关闭对应条目（如有）。

---

## 第四部分：风险与边界

### 风险 1：alias 链解析错误

`let row_sub := Var("_lambda_..._3")` 是简单情况。若 C++ 用 `auto row_sub = [&M](){...};` + 复制：`auto row_sub2 = row_sub;`，alias 链可能多层。

**缓解**：先 grep corpus 确认无多层 alias；若有，简化方案不处理（gap warning）。

### 风险 2：lambda 在 expr 上下文调用

C++ 可能在表达式中调用 lambda：`int x = compute() + lambda(args);`。当前 Pass 2b 仅处理 stmt-level；lambda 调用同样问题。

**缓解**：本方案先做 stmt-level（覆盖 8/8 已知 call sites）；expr-level 留 follow-up。

### 风险 3：lambda 作为 comparator 传给 std::sort

`_lambda___wang_core_3` / `_lambda___zassenhaus_recombine_1` 被传给 `std::sort` 作 comparator。这种"间接调用"无法在 outer body 静态识别 call site。

但 comparator lambdas 中的 modified captures（`idx`）在 sort 内部循环修改不传播——本质是 sort 实现的责任。当前 Pass 5 把 sort 当 functional `Array.sort` 处理；comparator 只读 lambda OK，但 mutating comparator (`[&idx]`) 行为未定义。

**缓解**：本方案不处理；标 P2 follow-up（与 sort/iota 模式 9 一起处理）。

### 风险 4：lifted lambda 的 body 内 mutating call

Lambda body 内：`M[i] -= ...` (compound assign on ArrayAccess root) → 走 Pass 6 B7 record-update 链。但 Pass 6 在 SSA build 时把 lambda 当独立 MIRFunc 处理（`ssa_build_pass(aux)` 递归调用），cap param `M` 是 lambda 的 input。lambda 内部的 SSA 链应正确——M_1 (param) → M_2 (after first row sub) → M_3 (after swap) → ...

**关键**：lifted lambda 的 ret 必须是 modified captures 的 latest version，且通过 ref_elim_pass 自动处理（_rewrite_returns 包装最后的 SSA version）。

### 回退策略

本方案两个修复（lifted lambda 签名改造 + call site 改写）必须**同时上**——只改签名不改调用，类型不匹配编译报错；只改调用不改签名，运行时 dimension 不匹配。

回退：revert Pass 3b 整体引入，回到当前状态（lambda mutation silent miss）。

---

## 第五部分：实施计划

| 步骤 | 内容 | 行数估计 |
|------|------|---------|
| P1 | `passes/pass3b_lambda_ref_elim.py` 新建（含 cap 提取 + ref_elim_pass 复用 + call site 重写） | +180 |
| P2 | smoke runner 注册（运行链插入 Pass 3b） | +5 各 file |
| P3 | unit test `test_pass3b_lambda_ref_elim.py`（5 用例） | +180 |
| P4 | smoke metric `lambda_ref_callsite_unwritten` | +25 |
| P5 | `mutation-model-design.md` §5 + `tech-debt.md` 同步 | +15 |
| 验证 | run_all_smoke.sh 全过 + 第四轮 1-对-1 审视（重点 `__lll_reduce`、`__wang_core`、`__mtshl_step_j`） | – |

### 验证标准

1. `bash tests/run_all_smoke.sh` 全绿
2. `smoke_pass6_full.py` 输出：
   - `lambda_ref_callsite_unwritten` = 0
   - `sideeff_other` ≤ 30（当前 37 → 减 8 = 29 左右）
3. `__lll_reduce.txt` dump 中：
   - `_lambda___lll_reduce_3` 调用点变成 `let __refret_X := ...`
   - 后续 `M_n / U_n` 使用 destructure 后的版本
4. unit test 总数：mir_types 9 + pass6 12 + pass2b 7 + **pass3b 新 5** = 33 通过
5. Agent B 第四轮语义模拟：`__lll_reduce` 评级"OK"

### 风险控制

- **先做最小 PoC**：仅处理 `_lambda___lll_reduce_3` 一个 call site，验证设计可行；再扩展到全部 6 lambda。
- **保留旧路径**：Pass 3b 默认开启，但加 `--no-lambda-ref-elim` flag（环境变量或参数）可关闭——回归排查用。

---

## 附录 A：grep 命令清单（验证用）

```bash
# Lifted lambda call sites（修复前）
grep -h "__sideeff_.*:= _lambda" /tmp/mir0_dump/*.txt
# 修复后期望：仅 comparator-passthrough 形态（_lambda_*_3 / _lambda_*_recombine_1 等）

# 修复后期望出现：
grep -h "let __refret_.*:= _lambda" /tmp/mir0_dump/*.txt

# 验证 M/U SSA chain 在 __lll_reduce
grep -nE "M_[0-9]+|U_[0-9]+" /tmp/mir0_dump/__lll_reduce.txt | head -30
```

## 附录 B：相关文档

- `docs/design/l1-translation-validation/mutation-model-design.md` §3, §5
- `docs/fixes/cpp2lean-v2-rangefor-and-callsite-ref-elim.md`（前置 P0-1/P0-2 修复）
- `proof/cpp2lean_v2/passes/pass3_lambda_lift.py`（lift 实现，含 modified_captures 检测）
- `proof/cpp2lean_v2/passes/pass2_ref_elim.py`（可复用的 ref_elim 逻辑）
- `proof/cpp2lean_v2/passes/pass2b_callsite_ref_elim.py`（可复用的 callsite destructure 逻辑）
