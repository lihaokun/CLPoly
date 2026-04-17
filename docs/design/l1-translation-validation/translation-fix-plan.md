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
