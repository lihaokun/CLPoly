# cpp2lean 翻译缺陷修复方案

## 现状

66 函数全量翻译，0 Unknown，0 sorry。但逐函数审计发现 **25 PASS / 41 FAIL**。

41 个 FAIL 函数的缺陷归结为 **6 个系统性根因**，非 case-by-case 问题。

## 根因与修复

### R1. `let _ :=` 丢弃变异结果（~20 函数，最高优先级）

**根因**：`ssa_transform.py:696` 的 fallback 把未识别的 mutation 留为 `ExprStmt`，`lean_codegen.py:331` 把所有 `ExprStmt` 无条件生成为 `let _ := expr`。

**触发模式**：

```cpp
// C++: 成员方法修改对象
f.normalization();           // → ExprStmt(Call("noop", [f])) — 正确: noop 不需要捕获
result.push_back(x);         // → ExprStmt(ArrayPush(result, x)) — 应转为 LetStmt
result[i].push_back(term);   // → ExprStmt(ArrayPush(ArrayAccess(...), term)) — 嵌套
```

**修复**：在 `ssa_transform.py` 的 `transform_expr_stmt` 中，扩展 mutation 识别：

1. **ArrayPush 在 ExprStmt 位置**：当前 line 624-629 已处理顶层 ArrayPush，但嵌套 ArrayPush（对 `arr[i]` 的 push）未覆盖。需要递归检测：如果 `ArrayPush.arr` 是 `ArrayAccess`，则生成 `Array.set! arr i (arr[i]!.push elem)`。

2. **`_mutate_` 前缀 Call 在 ExprStmt**：当前 line 615-623 已处理。确认覆盖。

3. **noop 方法（normalization, comp_ptr）在 ExprStmt**：当前作为 `ExprStmt(obj)` 不需要修改。确认 `_COMMON_METHODS` 中 "noop" 类直接返回 obj，codegen 生成 `let _ := obj` 无副作用——这是正确行为（normalization 在 Lean 模型中是 noop）。

4. **关键遗漏**：`result[i] = expr` 形式的数组赋值。C++ AST 中是 `CXXOperatorCallExpr(operator=, ArrayAccess, rhs)`。检查 `clang_ast.py` 的 operator= handler 是否产生正确的 `AssignStmt(ArrayAccess(...), rhs)`，以及 `ssa_transform.py` 是否把 ArrayAccess 赋值转换为 `Array.set!`。

**验证**：修复后重新生成，grep `let _ :=`，检查每处是否确实是 noop。

---

### R2. 结构化绑定空循环体（~8 函数）

**根因**：`clang_ast.py:682-712` 的 `_handle_decomposition_decl` 产生多个 LetStmt，但 `ssa_transform.py:1238` 的 range-for 解析只期望 `children[1]` 是单个 LetStmt（循环变量声明）。结构化绑定产生的多条 LetStmt 打乱了 children 的位置映射。

**触发模式**：
```cpp
for (const auto& [var, deg] : monomial) {
    // var 和 deg 在循环体内使用
}
```
Clang AST 产生的 RangeForLoop 结构：
- `children[0]` = collection
- `children[1]` = DecompositionDecl → 多个 LetStmt
- `children[2]` = CompoundStmt（循环体）

由于 children[1] 实际是一个 DecompositionDecl（被展开为多个 LetStmt），range-for 解析器的假设（children[1] 是单个循环变量声明）不成立。

**修复**：在 `ssa_transform.py` 的 `transform_range_for` 中：

1. 检测 `children[1]` 是否包含结构化绑定（多个 LetStmt 或 UnknownStmt("DecompositionDecl")）。
2. 如果是结构化绑定，解析绑定名列表（如 `[var, deg]`）。
3. 在循环函数体的开头插入解构语句：`let var := elem.fst`，`let deg := elem.snd`。
4. 循环变量仍然是单个 `elem`（或 `_decomp`），类型为 pair/tuple。

这样循环体内的 `var`、`deg` 引用就能正确解析。

**影响函数**：`__extract_monomial_content`、`__si_theta_array_eval`、`__mtshl_sparse_int`、`__select_eval_point`、`__wang_core`、`__wang_leading_coeff` 等。

---

### R3. `assign()` 变量替换缺失（~8 函数）

**根因**：`class_map.py:315` 把 `assign` 映射为 `("id", "identity")`，`clang_ast.py:985-986` 的 identity 规则只返回第一个参数，忽略 var 和 val。

**语义**：`assign(poly, var, val)` = 在多项式 `poly` 中用 `val` 替代变量 `var`。这是多变量多项式的核心操作（partial evaluation）。

**修复**：

1. `class_map.py`：将 `assign` 的映射改为 `("MvPoly.assign", "direct")`。
2. `clpoly_model.lean`：增加 `MvPolyZZ.assign` 和 `MvPolyZp.assign` 声明（opaque 或实现为变量替换）。
3. 同时检查 CLASS_MAP 中各类的 `assign` 方法映射（MvPolyZZ、MvPolyZp 等），确保 method 版本（`poly.assign(var, val)`）和 standalone 函数版本（`assign(poly, var, val)`）都正确映射。

**影响函数**：`__taylor_coeff`、`__taylor_coeff_zp`、`__mtshl_multi_bdp`、`__mtshl_wmds`、`__mtshl_step_j`、`__assign_partial_zp`、`__wang_leading_coeff` 等。

---

### R4. continue/break 状态丢失（~6 函数）

**根因**：`ssa_transform.py:886` 创建 for/while 循环的 LoopCtx 时，`all_param_names` 只包含 `modified_names`，不包含循环控制变量（如 for-loop 的计数器 `i`、while-loop 的条件变量）。对比 range-for（line 1297）正确包含了 `[idx_name, coll_name]`。

**修复**：在 `transform_for_loop` 和 `transform_while_loop` 中，构建 `all_param_names` 时包含全部循环函数参数：

```python
# for-loop
all_param_names = [loop_var_name] + modified_names + closure_names

# while-loop
all_param_names = modified_names + closure_names
```

确保 `_make_loop_recurse` 和 `_make_loop_return` 使用的参数列表与循环函数的实际参数一致。

**验证**：检查 `__edf_Zp`、`__upoly_subtract_x`、`__select_prime` 的 continue 递归调用参数数量是否与循环函数签名匹配。

---

### R5. 迭代器 compact 模式（~3 函数）

**根因**：`ssa_transform.py:50-79` 检测到迭代器模式（变量名含 `it`/`out`）后，**整个循环体替换为** `SparsePolyZZ.compactNonzero` 硬编码调用。这对 `__upoly_mod_coeff`（纯 compact）是正确的，但对 `__hensel_step`、`__hensel_step_linear` 的双指针模式（filter + transform）是错误的。

**修复方案**：

方案 A（精确）：区分不同的迭代器模式：
- **compact-nonzero**（filter zero terms）：保留当前 `compactNonzero` 映射。
- **compact-transform**（`__hensel_step` Part 2）：需要翻译为 `filterMap` + 具体变换逻辑。

方案 B（保守）：对 `__hensel_step` 和 `__hensel_step_linear` 使用 `FUNC_BODY_OVERRIDE`，手写 Lean 等价体。这两个函数的迭代器部分是固定的 Bezout 系数更新 + compact，可以用 `Array.filterMap` 精确表达。

**推荐方案 B**：只有 2-3 个函数受影响，FUNC_BODY_OVERRIDE 更可靠。

---

### R6. Lambda 捕获体为不透明存根（~3 函数）

**根因**：`clang_ast.py:1214` 注册 lambda 时保存了 `(captures, params, body)`，但 `ssa_transform.py:505` 解包时丢弃了 captures（`_, lam_params, lam_body =`），然后用 `collect_free_vars` 重新推导。对于复杂 lambda（多语句体、嵌套调用），重新推导的 capture 列表可能不完整。

**触发条件**：lambda body 较复杂（>3 条语句），或包含闭包变量的嵌套引用。

**修复**：

1. `ssa_transform.py:505`：正确解包 captures：`captures, lam_params, lam_body = _LAMBDA_REGISTRY[lambda_id]`。
2. 用 clang_ast 提供的 captures 作为 primary source，collect_free_vars 作为 fallback validation。
3. 确保 lambda 函数体的 SSA 变换正确传播捕获变量的版本号。

**影响函数**：`__mtshl_step_j`（`lc_correct`、`product_F` lambda）。

---

## 优先级排序

| 优先级 | 根因 | 影响函数数 | 修复复杂度 | 说明 |
|--------|------|-----------|-----------|------|
| P0 | R1 | ~20 | 中 | 最广泛的缺陷，需扩展 mutation 识别 |
| P0 | R2 | ~8 | 中 | 结构化绑定是 Wang 模块的核心模式 |
| P1 | R3 | ~8 | 低 | 改 FUNC_MAP 映射 + 加 model 声明 |
| P1 | R4 | ~6 | 低 | 扩展 LoopCtx 参数列表 |
| P2 | R5 | ~3 | 中 | 2-3 个 FUNC_BODY_OVERRIDE |
| P2 | R6 | ~3 | 中 | 修复 capture 解包逻辑 |

## R1-R6 修复后状态（v2 审计）

25 → 28 PASS。R1-R6 解决了表层问题，暴露了更深层的 SSA 变量管理缺陷。

---

## 第二轮残留根因（S1-S5）

### S1. if/else 分支变量不逃逸（~8 函数）

**现象**：`if (cond) { x = a; } else { x = b; }` 后，外层用的仍是 `x` 的旧版本。

**根因**：`transform_if` 的 phi-node 逻辑（line 786-854）**本身是正确的**，但有两个 gap：

1. **`pre_existing` 过滤过严**（line 790-791）：如果变量 `delta` 是在 if 之前通过 `LetStmt` 首次声明但初始值是占位符（如 `int64_t delta;` 声明无初始化），它可能不在 `env.versions` 中，因此不在 `pre_existing` 中，导致 phi 被跳过。

2. **`_extract_stmts_for_var` 只追踪 LetStmt 依赖链**（line 801-838）：如果分支中变量修改是通过 `ExprStmt`（如 `ArrayPush`）而非 `LetStmt`，提取器找不到 final_stmt，返回 `([], None)`，phi 退化为旧变量。

**修复**：
- 在 `transform_if` 的 diff 检测中，去掉 `pre_existing` 过滤——所有分歧变量都应生成 phi。
- 在 `_extract_stmts_for_var` 中，也搜索 `ExprStmt` 产生的变量修改。

**影响函数**：`__wang_leading_coeff`（delta）、`factorize`（result）等。

---

### S2. 循环内变量不传播到递归调用（~10 函数）

**现象**：循环体内 `p_2 = next_p(p_1)` 计算了新值，但递归调用传的是 `p_1`。

**根因**：**`identify_loop_vars_from_step` 只检测 `++`/`+=`/`AssignStmt`**（line 1024-1038），不检测 `ExprStmt(BinOp("=", Var("p"), Call("next_p", ...)))` 形式的 step。

具体案例：`__select_prime` 的 for-loop step 是 `p = next_p(p)`，这是 `ExprStmt(BinOp("=", ...))` 形式。`identify_loop_vars_from_step` 没有匹配这个模式，导致 `p` 不在 `modified_names` 中，不作为循环参数传递。

**修复**：在 `identify_loop_vars_from_step` 中添加 `BinOp("=")` 赋值检测：
```python
if isinstance(step_stmt, ExprStmt):
    expr = step_stmt.expr
    # p = next_p(p) 赋值形式
    if isinstance(expr, BinOp) and expr.op == "=":
        result |= _extract_all_root_var_names(expr.lhs)
```

**影响函数**：`__select_prime`（p 不推进）。

---

### S3. bi 乘积计算缺失（3 函数）

**现象**：`if(first) bi[i]=F[l] else bi[i]*=F[l]` 只翻译了 `first` 标志，没有实际赋值/乘法。

**根因**：`bi[i] = F[l]` 和 `bi[i] *= F[l]` 是 **2D 数组赋值**（`bi` 是 vector 的 vector）。clang_ast.py 的 `operator=` handler 把 `bi[i] = expr` 翻译为 `BinOp("=", ArrayAccess(bi, i), expr)`。

但在 SSA 变换中，`ArrayAccess` 赋值（line 668-674）只处理 `arr[i] = val` 的简单情况。`bi[i] = F[l]` 中 `F[l]` 本身也是 ArrayAccess，这是正确匹配的。

实际问题更可能是：`bi` 不在 `modified_names` 中（因为 `identify_loop_vars` 在嵌套 if 内的 ArrayAccess 赋值中找不到 `bi`），所以 `bi` 不是循环参数，SSA 变换时 `bi` 的修改被丢弃。

**修复**：`identify_loop_vars` 的 ExprStmt handler 加上对 `ArrayPush(ArrayAccess(...))` 和 `ArrayPush(FieldAccess(...))` 的检测。同时确保嵌套 if 内的 ArrayAccess 赋值被检测到。

---

### S4. Lambda 引用捕获写回（~5 函数）

**现象**：`lc_correct(F[i])` 修改了 `F[i]` 但结果未写回到 `F`。`row_sub(M, U, i, j, c)` 修改了 `M` 和 `U` 但结果丢弃。

**根因**：C++ lambda 通过 `[&]` 引用捕获外层变量，lambda 内部修改直接反映到外层。翻译器把 lambda 提取为独立函数，捕获变量作为参数传入。但调用方把 lambda 调用作为 `ExprStmt`（`let _ := _lambda_N_ir ...`），修改后的值没有写回到外层变量。

**C++ 语义**：
```cpp
auto row_sub = [&](int i, int j, ZZ c) {
    for (k = 0; k < n; k++) {
        M[i][k] -= c * M[j][k];  // 修改 M
        U[i][k] -= c * U[j][k];  // 修改 U
    }
};
row_sub(k, k-1, mu_kk1);  // 调用后 M, U 已被修改
```

**Lean 当前翻译**：
```lean
partial def _lambda_5_ir (M : LLLMatrix) (U : LLLMatrix) (i j : Int) (c : ZZ) : LLLMatrix :=
  -- ... 内部修改 M, U ...

let _ := (_lambda_5_ir M_1 U_1 k_1 (k_1 - 1) mu_kk1)  -- 结果丢弃！
-- 后续仍用 M_1, U_1（未修改版本）
```

**修复方案**：在 `ssa_transform.py` 的 ExprStmt 处理中，当检测到 lambda 调用（函数名含 `_lambda_`）且 lambda body 修改了捕获变量时，将 `ExprStmt` 转换为 `LetStmt`，捕获返回的修改值。

具体实现：
1. 在 lambda 提取阶段（line 507-541），记录哪些 capture 在 lambda body 中被修改（通过检查 lambda body 的 `identify_loop_vars`）。存入 `_LAMBDA_REGISTRY` 的第 4 个元素。
2. 在 ExprStmt 处理中（line 578+），如果调用的是 lambda 且有 modified captures，将结果拆包写回对应变量：
```python
# ExprStmt(Call("_lambda_5_ir", [M, U, i, j, c]))
# → let _lam_result := _lambda_5_ir M_1 U_1 i j c
#   let M_2 := _lam_result.1
#   let U_2 := _lam_result.2
```
3. Lambda 函数的返回值从 `auto` 改为修改后的 capture 变量 tuple。

**正确性论证**：

设 lambda 函数 $\lambda(c_1, ..., c_k, p_1, ..., p_m)$ 中，$c_1, ..., c_k$ 是捕获变量，$p_1, ..., p_m$ 是形式参数。设 lambda body 修改了 $c_i, c_j$（$\{i, j\} \subseteq \{1, ..., k\}$）。

C++ 语义：调用 $\lambda$ 后，外层的 $c_i, c_j$ 被原地修改为 $c_i', c_j'$。

Lean 翻译：$\lambda$ 返回 $(c_i', c_j')$，调用方将返回值拆包写回 $c_i, c_j$。

等价性：纯函数语义下，显式返回+写回 $\equiv$ 引用捕获修改。因为：
- $\lambda$ 内部对 $c_i$ 的所有读写都作用于参数 $c_i$（SSA 保证无别名）
- 返回的 $c_i'$ 是 $\lambda$ 执行完毕后的最终值
- 写回后，外层的 $c_i$ 被更新为 $c_i'$，等同于 C++ 的引用修改

**前置条件**：lambda 的捕获变量在调用期间无其他别名访问（C++ 中 `[&]` 捕获在单线程下满足）。

**影响函数**：`__lll_reduce`（`row_sub`, `row_swap`）、`__mtshl_step_j`（`lc_correct`, `product_F`）、`__mtshl_multi_bdp`（`compute_error`）。

---

### S5. 循环内 return 传播到父函数（~4 函数）

**现象**：`__select_eval_point` 的 while(true) 循环内 `return alpha` 被翻译为循环函数的 return，但调用方把循环结果按 `modified_vars` tuple 拆包，忽略了 alpha。最终函数返回 `StdMap.empty`（循环后的默认值）。

**根因**：loop-as-function 架构不区分"循环正常退出（break/条件不满足）"和"从父函数 return"。两者都被翻译为循环函数的 `ReturnStmt`，但语义完全不同。

**C++ 模式**：
```cpp
Result func(...) {
    while (true) {
        ...
        if (found) return result;  // 退出 func，不只是循环
        ...
    }
    return default_value;  // 循环耗尽后的 fallback
}
```

**当前翻译（错误）**：
```lean
partial def loop_N_ir (vars...) :=
  ...
  if found then return result  -- 返回 result 而非 modified_vars tuple！
  ...

partial def func_ir (...) :=
  let _loop := loop_N_ir(vars...)
  let var1 := _loop.1    -- 错误拆包：result 不是 tuple
  let var2 := _loop.2
  ...
  default_value           -- 总是返回默认值
```

**修复方案：return → break + 标志位**

核心思想：将循环内的 `return val` 转换为 `_ret = true; _ret_val = val; break`。循环函数照常返回 `modified_vars` tuple（含标志位和返回值），调用方检查标志位决定是继续执行还是直接返回。

**实现步骤**：

1. **检测**：在 `transform_for_loop` / `transform_while_loop` 的循环体变换之前，扫描原始 body 是否含 `ReturnStmt`（排除 break/continue 产生的 ReturnStmt）。

2. **注入标志变量**：如果含 ReturnStmt，在循环参数中添加两个合成变量：
   - `_ret_flag : Bool`（初始 `false`）
   - `_ret_val : ReturnType`（初始 `default`）
   加入 `modified_names`。

3. **变换 ReturnStmt**：循环体内的 `ReturnStmt(val)` 替换为：
   ```
   let _ret_flag := true
   let _ret_val := val
   break  // → return modified_vars tuple（含 _ret_flag = true）
   ```

4. **调用方检查**：循环调用后，拆包 `_ret_flag`：
   ```lean
   let _loop := loop_N_ir(...)
   let _ret_flag := _loop.ret_flag_accessor
   let _ret_val := _loop.ret_val_accessor
   if _ret_flag then return _ret_val
   else
     let var1 := _loop.var1_accessor
     ... 后续代码 ...
   ```

**正确性论证**：

设 C++ 函数 $f$ 含循环 $L$，$L$ 内有 $\texttt{return}\ e$。

C++ 语义：执行 $L$ 时，若到达 $\texttt{return}\ e$，则 $f$ 立即返回 $e$，不执行循环后代码。

Lean 翻译后语义：
- 循环函数 $L_{ir}$ 执行到原 return 点时，设 $\texttt{\_ret\_flag} := \texttt{true}$，$\texttt{\_ret\_val} := e'$（$e$ 的 SSA 翻译），然后 break（返回状态 tuple）。
- 调用方检查：$\texttt{\_ret\_flag} = \texttt{true}$，直接 return $\texttt{\_ret\_val}$。
- 循环后的代码不执行（被 `if _ret_flag then return` 短路）。

等价性证明：
- **return 路径**：C++ 返回 $e$ ↔ Lean 返回 $e'$（$e' = \text{rename}(e)$，SSA 等价）。循环后代码在两种语义下都不执行。✓
- **非 return 路径**：C++ 循环正常退出后执行后续代码 ↔ Lean $\texttt{\_ret\_flag} = \texttt{false}$，跳过 return，执行后续代码。变量通过 `modified_vars` tuple 正确传播。✓
- **不变量**：$\texttt{\_ret\_flag}$ 初始为 false，仅在原 return 点被设为 true。一旦设为 true 则立即 break，不会被后续迭代重置。✓

**边界情况**：
- 循环无 return（只有 break）：不注入标志变量，行为不变。✓
- 循环有多个 return：每个 return 点都设 flag + val + break，flag 语义一致。✓
- 嵌套循环内 return：外层循环提取时，内层循环的 return 已被内层处理转换为 flag+break，外层看到的是 `if _ret_flag then break`（传播 flag 到外层）。需要递归处理。

**影响函数**：`__select_eval_point`、`__wang_core`、`__edf_Zp`、`__zassenhaus_recombine`。

---

## S1-S3 已修复状态

| 根因 | 状态 | 修复方式 |
|------|------|---------|
| S1 | ✅ 已修复 | `transform_if` 两分支一致版本同步到外层 env |
| S2 | ✅ 已修复 | `identify_loop_vars_from_step` 添加 `BinOp("=")` |
| S3 | ✅ 已修复 | `identify_loop_vars` 扩展 `_extract_all_root_var_names(expr.arr)` |

## 修复优先级

| 优先级 | 根因 | 修复复杂度 | 预计影响 |
|--------|------|-----------|---------|
| **P0** | S5 return→break+flag | 高（~50 行） | ~4 函数，解决 return 传播 |
| **P1** | S4 lambda 写回 | 中（~30 行） | ~5 函数，解决引用捕获 |

修复 S4+S5 预计 PASS 率从 ~28/66 提升至 ~40-45/66。

---

## R1-R6 + S1-S5 修复后状态（v3 审计）

42 PASS / 24 FAIL (64%)。暴露第三层机制缺陷 M1-M5。

---

## 第三轮机制缺陷（M1-M5）

### M1. 输出参数只包装最后一条 ReturnStmt（4 函数）

**现象**：`__upoly_make_monic` 有两条 return 路径（early return + 正常 return）。输出参数 `f` 的修改后值只出现在函数末尾的死代码 `Prod.mk(0, f_1)` 中，实际 return 都只返回 `lc`。级联影响所有调用者。

**根因**：`transform_func` line 292 只检查 `clean_body[-1]`：
```python
if clean_body and isinstance(clean_body[-1], ReturnStmt) and clean_body[-1].value:
    orig_ret = clean_body[-1].value
    clean_body[-1] = ReturnStmt(Call("Prod.mk", [orig_ret] + out_vars))
```
当函数体最后一条语句是 `IfStmt`（含两条 return 路径）时，`isinstance(clean_body[-1], ReturnStmt)` 为 False，走 else 分支追加死代码。

**修复方案**：遍历 clean_body 中所有 ReturnStmt（含 IfStmt 嵌套），全部包装：

```python
def _wrap_all_returns(stmts, out_vars):
    """将所有 ReturnStmt(val) 替换为 ReturnStmt(Prod.mk(val, out_vars...))。"""
    result = []
    for s in stmts:
        if isinstance(s, ReturnStmt) and s.value:
            result.append(ReturnStmt(Call("Prod.mk", [s.value] + out_vars)))
        elif isinstance(s, IfStmt):
            new_then = _wrap_all_returns(s.then_body, out_vars)
            new_else = _wrap_all_returns(s.else_body, out_vars)
            result.append(IfStmt(s.cond, new_then, new_else))
        else:
            result.append(s)
    return result
```

然后在 `transform_func` 的输出参数处理中：
```python
if func.ret_type != BaseType.VOID:
    clean_body = _wrap_all_returns(clean_body, out_vars)
```

**正确性论证**：

设 C++ 函数 $f$ 有返回类型 $R$ 和输出参数 $o_1, ..., o_k$。任意 return 路径 $p$ 返回值 $r_p$，此时 $o_1, ..., o_k$ 的值分别为 $o_1^p, ..., o_k^p$。

C++ 语义：调用方收到 $r_p$，同时 $o_i$ 通过引用被修改为 $o_i^p$。

修复后 Lean：每条路径 $p$ 返回 $(r_p, o_1^p, ..., o_k^p)$。调用方拆包获得全部值。

等价性：$\forall p$，返回值包含同一组 $(r_p, o_1^p, ..., o_k^p)$ ✓。递归包装确保嵌套 if/else 内的 return 不遗漏。

**影响函数**：`__upoly_make_monic` → `__squarefree_Zp`、`__ddf_Zp`、`__factor_Zp`

---

### M2. for-loop continue 不执行 step（~3 函数）

**现象**：C++ `for (tried = 0; tried < max; p = next_p(p))` 的 continue 应先执行 step `p = next_p(p)` 再进入下一次迭代。翻译器的 continue 直接递归循环函数，跳过 step。

**根因**：`transform_body` line 442-446 的 `ContinueStmt` 处理：
```python
if stmt.kind == "ContinueStmt":
    recurse = _make_loop_recurse(loop_ctx, env)
    result.append(ReturnStmt(recurse))
    return result
```

`_make_loop_recurse` 只用当前 env 的变量版本构造递归调用，不包含 step 语句。但 C++ for-loop 的 continue 语义是：**step → condition → body**，step 不可跳过。

**修复方案**：在 `ContinueStmt` 处理中，如果 `loop_ctx.step_stmts` 非空（for-loop），先执行 step 再递归：

```python
if stmt.kind == "ContinueStmt":
    # for-loop: continue 必须先执行 step
    if loop_ctx.step_stmts:
        for step_s in loop_ctx.step_stmts:
            result.append(step_s)
    recurse = _make_loop_recurse(loop_ctx, env)
    result.append(ReturnStmt(recurse))
    return result
```

**注意**：`step_stmts` 是在循环体变换结束后才 SSA 变换的（`transform_step` 用 `loop_env`）。在 continue 点，env 可能与 step 期望的版本不同。但对于简单的 step（如 `p = next_p(p)`、`i++`），step 只依赖自身变量的当前版本，不依赖 body 中的其他变量。

对于更复杂的 step，需要在 continue 点重新变换 step 表达式。安全方案：在 continue 时，用当前 env 重新 transform 原始 step_stmt：

```python
if stmt.kind == "ContinueStmt":
    if loop_ctx.step_stmts:
        # 用当前 env 重新变换 step（确保变量版本正确）
        fresh_step = transform_step(loop_ctx.raw_step, env)
        result.extend(fresh_step)
    recurse = _make_loop_recurse(loop_ctx, env)
    result.append(ReturnStmt(recurse))
    return result
```

这需要在 LoopCtx 中保存 `raw_step`（原始 step 语句，SSA 前）。

**正确性论证**：

C++ for-loop: `for (init; cond; step) { body }`

语义：$\text{init}; \text{while}(\text{cond}) \{ \text{body}; \text{step}; \}$

continue 语义：跳过 body 剩余部分，执行 step，重新检查 cond。

修复后 Lean：continue 点插入 step 的 SSA 变换，然后递归。递归时变量已包含 step 的修改（如 `p_2 = next_p(p_1)`）。等价于 C++ 的 step 执行。✓

while-loop 的 continue 不需要 step（while 没有 step 表达式），保持原逻辑。✓

**影响函数**：`__select_prime`（`p = next_p(p)` step 被 continue 跳过）、`__upoly_subtract_x`

---

### M3. while-loop _ret_flag 初始化缺失（~3 函数）

**现象**：`__edf_Zp` 的 while 循环调用 `while_67466_ir _ret_flag _ret_val ...`，但 `_ret_flag` 和 `_ret_val` 从未声明/初始化。

**根因**：S5 实现中，_ret_flag 初始化代码只加在了 `transform_for_loop` 的 `result_stmts` 组装处。`transform_while_loop` 虽然注入了 `_ret_flag`/`_ret_val` 到 `env.versions` 和 `modified_names`，但**没有在 `result_stmts` 中追加初始化 LetStmt**。

**修复方案**：在 `transform_while_loop` 的 `result_stmts` 组装处，与 `transform_for_loop` 对称地加初始化：

```python
result_stmts = []
# S5/M3: 初始化 _ret_flag/_ret_val
if has_parent_ret:
    result_stmts.append(LetStmt(env.current("_ret_flag"), BaseType.BOOL, Lit(False)))
    result_stmts.append(LetStmt(env.current("_ret_val"), "auto", Lit(0)))
```

**正确性论证**：

_ret_flag 初始值为 false（无 return 发生）。循环函数接收此初始值作为参数。循环体内 return 时设为 true。循环退出后检查 flag。初始化确保 flag 在循环调用前有定义。✓

**影响函数**：`__edf_Zp`、`__vanhoeij_recombine`

---

### M4. Lambda 写回检测不完整（~2 函数）

**现象**：S4 修复让 `row_swap` 的 M/U 写回生效（4 处 `_lam_out`），但 `row_sub` 的 M/U 修改仍为 `let _` 丢弃。

**根因**：S4 通过 `identify_loop_vars` 检测 lambda body 中的 modified captures。`row_swap` 的 body 直接用 `BinOp("=")` 赋值，被检测到。但 `row_sub` 的 body 通过嵌套 for-loop 修改 M/U（`M[i][k] -= c * M[j][k]`），`identify_loop_vars` 对 `CompoundStmt` 递归但不深入 `ForLoop`/`WhileLoop` 的 children。

具体来说，`identify_loop_vars` line 1783-1785：
```python
elif isinstance(stmt, UnknownStmt):
    for c in stmt.children:
        result |= identify_loop_vars(c)
```
`ForLoop` 是 `UnknownStmt(kind="ForLoop")`，其 children 包含 `[init, cond, step, body]`。递归会检测 body 中的赋值，所以**应该能检测到**。

实际问题可能是：`M[i][k] -= c * M[j][k]` 在 AST 中被解析为 `CXXOperatorCallExpr(operator-=, ArrayAccess(ArrayAccess(M, i), k), ...)` → 翻译为 `BinOp("=", ArrayAccess(ArrayAccess(M,i),k), BinOp("-", ...))` 。`_extract_all_root_var_names` 对 `ArrayAccess(ArrayAccess(M,i),k)` 应该能提取出 `M`。

如果检测确实遗漏，**保守修复**：对所有 `[&]` capture 的 lambda，假设所有 capture 都可能被修改（`modified_caps = capture_names`）。这会产生一些多余的写回（noop captures 也被返回），但不会错误。

```python
# 保守策略：[&] capture = 所有 capture 都可能被修改
modified_caps = list(capture_names)
```

**正确性论证**：

保守策略：返回所有 capture 值。如果某个 capture 未被修改，返回值 == 输入值，写回后变量不变。语义等价。✓

精确策略需要完善 `identify_loop_vars` 对所有嵌套结构的递归，且正确解析所有赋值模式（`-=`、`*=`、`[]`、`.field` 等）。

**影响函数**：`__lll_reduce`（row_sub）、`__mtshl_step_j`（lc_correct、product_F）

---

### M5. 复杂循环体逻辑丢失（~7 函数，个别排查）

这些不是单一机制缺陷，需要逐个排查 AST 解析链路。按子类分：

#### M5a. range-for body 内 if/push 逻辑丢失（2 函数）

**函数**：`__upoly_subtract_one`、`__upoly_subtract_x`

**现象**：range-for body 内的 `if (deg==0) { compute; push } else { push }` 逻辑部分丢失。

**可能根因**：
1. `parse_range_for` 的 body 提取可能遗漏某些 CompoundStmt 子节点
2. range-for body 内的 continue 使后续 push 未执行
3. SSA 变换中 if/else 分支的变量作用域不正确

**排查方法**：dump `__upoly_subtract_one` 的原始 IR（parse 后、SSA 前），检查 RangeForLoop 的 body 是否包含完整的 if/else + push 逻辑。如果 IR 完整则问题在 SSA；如果 IR 缺失则问题在 AST parser。

#### M5b. do-while 内层循环空壳（1 函数）

**函数**：`__zassenhaus_recombine`

**现象**：`do { ... } while (next_combination)` 的 body 只包含 `next_combination` 调用但无 subset-product / trial-division 逻辑。

**可能根因**：do-while 翻译为 while-loop 时，body 内容可能被 `parse_while` 或 `parse_stmt` 的 do-while 分支丢弃。

**排查方法**：检查 Clang AST 中 do-while 的 kind（`DoStmt`），确认 `parse_stmt` 是否有 DoStmt handler。

#### M5c. 迭代器 compact 模式（2 函数）

**函数**：`__hensel_step`（Part 2）、`__hensel_step_linear`

**现象**：双指针 compact 循环被错误简化或留为空壳。

**根因**：S5 修改后 `_has_iterator_ops` 更保守（只检测 it/out 双指针名），这两个函数的迭代器循环不再被 compactNonzero 覆盖，但原始翻译仍不完整。

**修复方案**：对这两个函数使用 `FUNC_BODY_OVERRIDE`，手写 Lean 等价体。迭代器 Part 2 的语义是 `Array.filterMap (fun term => let c = fdiv_r(fdiv_q(term.snd, m), p); if c != 0 then some (term.fst, c) else none)`。

#### M5d. array element write-back（2 函数）

**函数**：`__hensel_lift`（factors_adj[0] 修改未写回）、`__build_cld_matrix`（row extend 未写回）

**现象**：对数组元素的遍历修改（`for (auto& row : M) row.push_back(0)`）通过 range-for 的 `__coll` 参数传递，但 write-back 机制（R1 修复的 `coll_needs_writeback`）未正确检测到修改。

**排查方法**：检查 `loop_var_ver > 0`（write-back 触发条件）是否满足。如果 loop body 修改了 `row`（通过 push），`loop_var_ver` 应该 > 0。

---

## 修复优先级

| 优先级 | 机制 | 修复量 | 预计新增 PASS |
|--------|------|--------|-------------|
| **P0** | M1 输出参数全包装 | ~15 行 | +4（含级联） |
| **P0** | M2 continue 执行 step | ~10 行 + LoopCtx 扩展 | +3 |
| **P0** | M3 while-loop _ret_flag init | ~3 行 | +3 |
| **P1** | M4 lambda 保守写回 | ~3 行（改 1 行策略） | +2 |
| **P2** | M5a-d 逐个排查 | 各 ~10 行 | +5-7 |

**预期**：M1+M2+M3 修复后 PASS 从 42/66 → ~52/66 (79%)。全部修复后 ~57-59/66 (87-89%)。

---

## M1-M4 修复后状态（v4 审计）

v3: 42 PASS → v4: 35 PASS（**回归 7 个**）。M1 引入调用方不兼容。

---

## 第四轮机制缺陷（N1-N5）

### N1. 输出参数调用方未自动解构（~8 函数，最高优先级）

**现象**：M1 修复了 `__upoly_make_monic` 的返回值（`Prod.mk(lc, f)`），但所有调用者仍写 `let g := __upoly_make_monic_ir g`，把 `(lc, f)` pair 当作多项式使用。级联破坏 `__squarefree_Zp`、`__ddf_Zp`、`__edf_Zp`、`__factor_Zp`。

**根因**：`TRANSLATION_SCOPE_OUTPUT_PARAMS` 的调用方适配代码只处理 `ExprStmt(Call(...))` 位置（`transform_stmt` line 592-628）。但 `__upoly_make_monic` 的调用全是 `LetStmt(var, typ, Call(...))` 形式——调用结果被绑定到变量，不是独立的语句。LetStmt 的 Call 不经过 ExprStmt handler，所以 output param 解构被跳过。

**修复方案**：在 `transform_stmt` 的 LetStmt handler 中，检测 Call 目标是否在 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 中。如果是，自动解构返回值：

```python
# LetStmt: let g := __upoly_make_monic_ir(g)
# → let _out := __upoly_make_monic_ir(g)
#   let lc := _out.1    (原返回值)
#   let g := _out.2     (输出参数)
if isinstance(stmt, LetStmt) and isinstance(stmt.value, Call):
    func_base = stmt.value.func.replace("_ir", "")
    if func_base in TRANSLATION_SCOPE_OUTPUT_PARAMS:
        out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[func_base]
        # stmt.var 绑定原返回值（如 lc）
        # 输出参数从 pair.2, pair.2.2, ... 提取并写回 env
        ...
```

**正确性论证**：

设 C++ 函数 $g = f(\text{out1}, \text{out2}, ...)$ 有返回值 $g$ 和输出参数 $\text{out}_i$。

M1 将函数返回类型改为 $(g, \text{out}_1, \text{out}_2, ...)$。

N1 在调用方自动解构：
- `let _out := f(out1, out2, ...)`
- `let ret_val := _out.1` （原来 `let g := f(...)` 的语义）
- `let out1_new := _out.2` （输出参数新值）

调用方原来写 `let g := f(args)`，假设 $f$ 返回单值。修复后 `g` 绑定 `_out.1`（仍是原返回值），`out1` 通过 `_out.2` 更新。语义一致。✓

**影响函数**：`__squarefree_Zp`（2 处）、`__ddf_Zp`（6 处）、`__edf_Zp`（3 处）、`__factor_Zp`（1 处）、`__upoly_make_monic` 自身的 early-return `f_1` 引用

---

### N2. step_stmts 在 continue 处理时仍为空（~2 函数）

**现象**：M2 修复在 `ContinueStmt` 处插入 `loop_ctx.step_stmts`，但 `__select_prime` 的 continue 路径仍传 `p_1`（step 未执行）。

**根因**：`transform_for_loop` 中的执行顺序：

```python
lctx = LoopCtx(func_name, ...)         # step_stmts 默认 = []
body_stmts = transform_body(..., lctx)  # continue 在这里处理
                                        # 此时 lctx.step_stmts = [] ← 空！
step_stmts = transform_step(step, env)  # 这里才计算 step
lctx.step_stmts = step_stmts           # 太晚
```

`transform_body` 在 `transform_step` 之前执行。body 中遇到 `ContinueStmt` 时，`loop_ctx.step_stmts` 还是空列表。

**修复方案**：在 body 变换之前计算 step_stmts，但用 fork 的 env 避免 step 的副作用影响 body：

```python
# 先计算 step（用 fork env，不影响 body 的变量版本）
step_env = loop_env.fork()
step_stmts = transform_step(step_stmt, step_env)
lctx.step_stmts = step_stmts

# 再变换 body（body 中的 continue 可以正确使用 step_stmts）
body_stmts = transform_body(flatten_stmts(body_stmt), loop_env, ctx, loop_ctx=lctx)

# body 后再次计算 step（用 body 结束后的 loop_env，给正常路径用）
step_stmts_final = transform_step(step_stmt, loop_env)
```

**正确性论证**：

- continue 路径：使用 `step_stmts`（fork env 计算），step 变量版本 = 循环参数的初始版本。对于 `p = next_p(p)` 类型的 step，`p` 在 body 中未修改，初始版本 = 任意 continue 点的版本。✓
- 正常路径：使用 `step_stmts_final`（body 结束后的 env 计算），step 变量版本 = body 执行后的最新版本。✓
- 已验证 66 个函数中 0 个 conflict（step 变量不在 body 中被修改），两者版本一致。✓

**影响函数**：`__select_prime`（`p = next_p(p)` step 在 continue 时跳过）

---

### N3. Lambda 内部循环的数组 mutation 不在 modified_names（~2 函数）

**现象**：`_lambda_5_ir`（`__lll_reduce` 的 `row_sub`）body 有 `loop_9048_ir` 循环修改 `M[i][k]` 和 `U[i][k]`，但循环返回值只有 `k`（循环计数器），M 和 U 不在循环的 `modified_names` 中。Lambda 返回 `(Prod.mk M U_2 n_1)` — M 和 U 是原始输入值。

**根因**：C++ 中 `M[i][k] -= c * M[j][k]` 是 `CXXOperatorCallExpr(operator-=, ...)`。Clang AST 翻译可能产生：
1. `BinOp("=", ArrayAccess(...), BinOp("-", ...))` — `identify_loop_vars` 能检测
2. `ExprStmt(Call("operator-=", [ArrayAccess(...), ...]))` — `identify_loop_vars` **不能**检测

如果走的是路径 2，`M` 不会被加入 `modified_names`。

**修复方案**：在 `identify_loop_vars` 中扩展 compound assignment operator 检测：

```python
# operator-=, operator+=, operator*= 等 compound assignment
if isinstance(expr, Call) and expr.func.startswith("operator"):
    if any(op in expr.func for op in ["-=", "+=", "*=", "/=", "%=", "|=", "&="]):
        if expr.args:
            result |= _extract_all_root_var_names(expr.args[0])
```

**正确性论证**：`operator-=(M[i][k], val)` 等价于 `M[i][k] = M[i][k] - val`，即 `M` 被修改。将 `M` 加入 `modified_names` 确保循环函数返回修改后的 `M`。✓

**影响函数**：`__lll_reduce`（`row_sub` lambda 的内部 for 循环）

---

### N4. 重复函数定义（编译错误）

**现象**：模板实例化（`grlex` + `lex` 两个特化）产生同名的 helper 函数（如 `range_90057_ir` 出现 4 次）。Lean 不允许同名 `partial def`。

**根因**：`gen_full.py` 收集所有翻译函数并按依赖排序输出，但不去重。同一个 C++ 模板函数被 `grlex` 和 `lex` 各实例化一次，loop hash 可能相同。

**修复方案**：在 `gen_full.py` 的输出阶段去重——按函数名去重，保留第一个：

```python
seen_names = set()
for sf in ordered:
    code = gen_ssa_func(sf)
    # 提取所有 partial def 名
    names = re.findall(r'partial def (\S+)', code)
    # 去重
    if all(n not in seen_names for n in names):
        print(code)
        seen_names.update(names)
```

或更简单：在 aux_defs 收集阶段，按名称去重。

**影响**：全局，修复后可编译。

---

### N5. 类型别名缺少 struct 字段（~5 函数）

**现象**：`Factorization.factors`、`Factorization.content`、`PrimeSelectionResult.factors`、`WangLcResult.success` 等字段在 Lean 中不存在。`clpoly_model.lean` 将这些定义为 `Array Int` 别名。

**根因**：`clpoly_model.lean` 的 Factorization 等类型用 `abbrev` 定义为 `Array Int`（占位），没有真正的 struct 字段。但 C++ 中这些是 `struct` 类型有 `.factors()`、`.content()` 等方法。

**修复方案**：将 `Factorization`、`PrimeSelectionResult`、`WangLcResult` 从 `abbrev` 改为 `structure`：

```lean
structure Factorization where
  content : ZZ
  factors : Array (SparsePolyZZ × UInt64)
deriving Inhabited

structure PrimeSelectionResult where
  p : UInt64
  factors : Array SparsePolyZp
  nfactors : UInt64
deriving Inhabited

structure WangLcResult where
  success : Bool
  scaled_factors : Array MvPolyZZ
  lc_targets : Array MvPolyZZ
  delta : ZZ
deriving Inhabited
```

并在 `CLASS_MAP` 中更新对应的 method 映射为 field access。

**影响函数**：`factorize`、`__factor_Zp`、`__select_prime`、`__wang_leading_coeff`、`__factor_multivar`

---

## 修复优先级

| 优先级 | 缺陷 | 修复量 | 预计新增 PASS |
|--------|------|--------|-------------|
| **P0** | N1 调用方解构 | ~20 行 | +7（M1 回归恢复 + 级联） |
| **P0** | N2 step_stmts 时序 | ~5 行 | +2 |
| **P1** | N3 operator-= 检测 | ~5 行 | +2 |
| **P1** | N4 去重 | ~10 行 | 编译通过 |
| **P1** | N5 struct 定义 | ~30 行 model + ~10 行 map | +3-5 |

**预期**：N1+N2 修复后恢复到 v3 水平（42 PASS）并超过。全部修复后预计 ~50-55/66 (76-83%)。

---

## N1-N5 + 黑名单重构后状态（v5 审计）

60 PASS / 6 FAIL (91%)。0 sorry, 0 WARNING。

---

## 第五轮残留缺陷（P1-P3，6 函数）

### P1. S1 版本同步只更新 env 不生成代码（2 函数）

**现象**：`__wang_leading_coeff` 的 `delta_2` 在 if/else 两分支内各定义，但外部 `if (delta_2 == 0)` 引用时 `delta_2` 不在作用域。`__extract_monomial_content` 的 `var_factors_2` 在 early-return 路径引用了后面才定义的变量。

**根因**：S1 修复只做了 `env.versions[name] = v_then`（更新 env），没有生成 `LetStmt`。Lean 的 `let` 在 `if` 分支内是局部作用域——即使两分支都 `let delta_2 := ...`，外部也看不到 `delta_2`。

需要生成 phi-node CondExpr：
```lean
let delta_2 := if cond then
  let delta_2 := ... -- then branch computation
  delta_2
else
  let delta_2 := ... -- else branch computation  
  delta_2
```

**修复**：在 S1 同步逻辑中，不仅更新 env，还用 `_extract_stmts_for_var` 提取两分支的值并生成 `LetStmt(new_var, typ, CondExpr(cond, then_expr, else_expr))`。这与 diff 非空时的 phi 逻辑（line 898-912）完全相同。

```python
# 当前（只同步 env）：
if v_then == v_else and v_then > v_outer:
    env.versions[name] = v_then

# 修复（同步 env + 生成 phi）：
if v_then == v_else and v_then > v_outer:
    then_stmts, then_final = _extract_stmts_for_var(then_body, name)
    else_stmts, else_final = _extract_stmts_for_var(else_body, name)
    then_expr = BlockExpr(then_stmts, then_final) if then_stmts and then_final else Var(name, v_then)
    else_expr = BlockExpr(else_stmts, else_final) if else_stmts and else_final else Var(name, v_then)
    new_var = env.bump(name)
    result.append(LetStmt(new_var, env.get_type(name), CondExpr(cond, then_expr, else_expr)))
```

**正确性**：与 diff 非空时的 phi 逻辑一致——提取分支值、CondExpr 合并。区别是两分支版本相同（而非不同），但生成 phi 的方式完全一样。✓

**影响函数**：`__wang_leading_coeff`（delta）、`__extract_monomial_content`（var_factors）

---

### P2. 循环体逻辑缺失（2 函数）

**`__upoly_subtract_one`**：range-for body 只跟踪 `found` 标志，不做 `result.push_back(term)`（非 deg-0 项）和减法后 push（deg-0 项）。可能是 AST 解析阶段循环体的 if/else 内容丢失。

**`__edf_Zp`**：base case 和 deg≤0 case 返回 `()` 而非 `(result, rng)`。这些是 output-param 函数的非主路径 return，`_wrap_all_returns` 没覆盖到——因为它们是 if/else 表达式的值（不是 ReturnStmt）。

**排查方向**：
1. `__upoly_subtract_one`：dump raw IR 确认循环体内容是否完整
2. `__edf_Zp`：`()` 值来源——可能是 `return` 被翻译为 bare `()`（void return → Unit）

---

### P3. N1 调用方解构遗漏（2 函数）

**`__upoly_make_monic`**：自身 early-return 用 `f_1`（else 分支才定义），应用 `f`。同时 loop helper 传了未绑定的 `__coll`/`term`。

**`__factor_Zp`**：`let lc_1 := __upoly_make_monic_ir f` 没有 N1 解构。检查：该调用是 `ExprStmt(BinOp("="))` 还是 `LetStmt`？如果是 BinOp 赋值，N1 BinOp handler 应匹配。如果是 LetStmt + Cast，N1 LetStmt handler 应匹配。两者都没生效——需要 debug 具体原因。

---

## 修复优先级

| 优先级 | 缺陷 | 修复量 | 影响 |
|--------|------|--------|------|
| **P0** | P1 phi 生成代码 | ~10 行 | +2 函数 |
| **P1** | P3 N1 解构 debug | ~5 行 | +2 函数 |
| **P2** | P2 循环体排查 | 逐个 | +2 函数 |

修复 P1+P3 预计 PASS 从 60 → 64/66 (97%)。

### 调研结论：6 FAIL 只有 2 个根因

**根因 A: S1 phi 同步只改 env 不生成代码（5 函数）**

`transform_if` 的 S1 逻辑在两分支版本一致时做 `env.versions[name] = v_then` 但不生成 `LetStmt`。Lean 的 `let` 在 if 分支内是局部作用域——即使两分支都 `let delta_2 := ...`，外部也看不到。

影响函数：
- `__wang_leading_coeff`（delta_2 不可见）
- `__extract_monomial_content`（var_factors_2 不可见）
- `__upoly_subtract_one`（result push 在两分支都做了但外部看不到）
- `__upoly_make_monic`（f_1 在 else 分支定义但 early return 引用）
- `__edf_Zp`（result_1 在 if 分支定义但外部 Prod.mk 引用）

**修复**：S1 同步逻辑中，对版本一致的变量也生成 CondExpr phi——与 diff 非空时的逻辑（line 898-912）完全相同：

```python
# 当前（只同步 env）：
for name in pre_existing:
    v_then = env_then.versions.get(name, 0)
    v_else = env_else.versions.get(name, 0)
    v_outer = env.versions.get(name, 0)
    if v_then == v_else and v_then > v_outer:
        env.versions[name] = v_then  # ← 只更新 env

# 修复（同步 env + 生成 phi）：
agreed = {}
for name in pre_existing:
    v_then = env_then.versions.get(name, 0)
    v_else = env_else.versions.get(name, 0)
    v_outer = env.versions.get(name, 0)
    if v_then == v_else and v_then > v_outer:
        agreed[name] = (Var(name, v_then), Var(name, v_then))

# 与 diff 变量一起处理（共用 _extract_stmts_for_var + CondExpr 生成）
for name in sorted(set(diff.keys()) | set(agreed.keys())):
    if name in diff:
        then_var, else_var = diff[name]
    else:
        then_var, else_var = agreed[name]
    then_stmts, then_final = _extract_stmts_for_var(then_body, name)
    else_stmts, else_final = _extract_stmts_for_var(else_body, name)
    then_expr = BlockExpr(then_stmts, then_final) if then_stmts and then_final else then_var
    else_expr = BlockExpr(else_stmts, else_final) if else_stmts and else_final else else_var
    new_var = env.bump(name)
    result.append(LetStmt(new_var, env.get_type(name), CondExpr(cond, then_expr, else_expr)))
```

**正确性论证**：

两分支版本一致（都是 `v`）意味着两分支都定义了 `name_v`。在 Lean 的作用域模型中，分支内的 `let name_v := ...` 只在分支内可见。生成 `let name_{v+1} := if cond then <then_value> else <else_value>` 将两分支的值合并到外部作用域。

语义：$\text{name}_{v+1} = \begin{cases} e_{\text{then}} & \text{if cond} \\ e_{\text{else}} & \text{otherwise} \end{cases}$

两分支独立计算 `name_v`，phi 取正确的那个。等价于 C++ 的 if/else 赋值。✓

---

**根因 B: N1 LetStmt handler 未匹配 BinOp 包裹的 Call（1 函数）**

`__factor_Zp` 的 IR：`LetStmt(lc, Zp, BinOp("=", Var("f"), Call("__upoly_make_monic", [f])))`

C++ 是 `Zp lc = __upoly_make_monic(f)`——Clang 把引用参数调用翻译为赋值表达式 `f = __upoly_make_monic(f)`，然后赋值表达式的值作为 `lc` 的初始化器。

N1 LetStmt handler 检查 `isinstance(call_value, Call)` 但 `call_value` 是 `BinOp("=", ...)`，不匹配。

**修复**：在 N1 LetStmt handler 中增加 BinOp("=") 解包：

```python
call_value = stmt.value
if isinstance(call_value, Cast):
    call_value = call_value.expr
# 新增：BinOp("=") 包裹的 Call（引用参数调用）
if isinstance(call_value, BinOp) and call_value.op == "=" and isinstance(call_value.rhs, Call):
    call_value = call_value.rhs
if isinstance(call_value, Call):
    ...  # N1 解构逻辑
```

同时，`BinOp("=")` 的 LHS `Var("f")` 表示 `f` 被修改——这已经由 output-param 机制处理（`f` 是 index 0 output param），所以解构时 `_out.2` 就是修改后的 `f`。

**正确性**：`BinOp("=", Var("f"), Call("__upoly_make_monic", [f]))` 的语义是：调用 `__upoly_make_monic(f)` 修改 `f` 并返回 `lc`。解构后：`lc := _out.1`（返回值），`f := _out.2`（修改后的 f）。与 C++ 语义一致。✓

---

## 架构重构：ExprStmt 黑名单策略

### 问题

4 轮 20 个 fix 中 7 个是变异检测遗漏。根源：`transform_stmt(ExprStmt)` 的 fallthrough（line 803）默认 `return ExprStmt(rename_expr(expr, env))` → `let _ := expr`，静默丢弃任何未识别的变异。

### 重构方案

将 fallthrough 从**静默通过**改为**分类检查**：

```python
# === 已知 noop（无副作用，可安全丢弃） ===
KNOWN_NOOP_METHODS = {
    "normalization", "comp_ptr", "reserve", "toList", "data",
    "MvPolyZp.normalization", "MvPolyZZ.normalization",
    "MvMonomial.normalization", "SparsePolyZp.normalization",
    "SparsePolyZZ.normalization",
    "SparsePolyZp.toList", "SparsePolyZZ.toList",
}

KNOWN_NOOP_EXPRS = (Var, Lit)  # 变量/字面量作为语句 = noop

def _classify_expr_stmt(expr) -> str:
    """分类 ExprStmt 中的表达式。返回 "noop"/"mutation"/"unknown"。"""
    # 1. 字面量、变量引用 → noop
    if isinstance(expr, (Var, Lit)):
        return "noop"
    # 2. CondExpr 作为语句 → noop（assert 检查）
    if isinstance(expr, CondExpr):
        return "noop"
    # 3. 已知 noop 方法调用
    if isinstance(expr, Call):
        func = expr.func if isinstance(expr.func, str) else ""
        if func in KNOWN_NOOP_METHODS:
            return "noop"
        if func in ("poly_convert", "SparsePolyZZ.compactNonzero"):
            return "noop"  # 副作用已由其他机制处理
    # 4. FieldAccess 作为语句 → noop（读取但不赋值）
    if isinstance(expr, FieldAccess):
        return "noop"
    # 5. ArrayAccess 作为语句 → noop
    if isinstance(expr, ArrayAccess):
        return "noop"
    # 6. Cast 作为语句 → noop
    if isinstance(expr, Cast):
        return "noop"
    # 7. 其他 → unknown（应报警）
    return "unknown"
```

在 `transform_stmt` 的 ExprStmt fallthrough 处：

```python
# 旧代码：
return ExprStmt(rename_expr(expr, env))

# 新代码：
classification = _classify_expr_stmt(expr)
if classification == "noop":
    return ExprStmt(rename_expr(expr, env))  # 安全丢弃
else:
    import sys
    print(f"WARNING: unhandled mutation in ExprStmt: {expr}", file=sys.stderr)
    return ExprStmt(rename_expr(expr, env))  # 仍生成，但报警
```

### 同步修复 N1-N5

在重构的同时修复 N1-N5：

1. **N1**（调用方解构）：在 LetStmt handler 中也检测 `TRANSLATION_SCOPE_OUTPUT_PARAMS`
2. **N2**（step 时序）：`transform_for_loop` 中用 fork env 提前计算 step
3. **N3**（operator-=）：在 mutation handler 中加 compound assignment 检测
4. **N4**（去重）：`gen_full.py` 按函数名去重
5. **N5**（struct 字段）：`clpoly_model.lean` 中 Factorization 等改为 structure
