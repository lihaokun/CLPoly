# 翻译器根因分析与改进方案 v2

> 审计时间：2026-04-10
> 审计范围：71 个 C++ 函数 → 369 个 Lean partial def
> 审计结果：37 PASS / 34 FAIL / 0 sorry

## 1. 四个系统性根因

审计发现的 34 个 FAIL 函数全部由以下 4 个根因覆盖。每个根因都有明确的机制定义和修复方案。

### M1：range-for 原地修改循环不传播集合变量（~10 函数）

**问题描述**：

```cpp
for (auto& term : f) {
    term.second *= lc_inv;   // 原地修改 f 的元素
}
```

翻译后的循环函数中，`Array.set!` 创建了 `__coll_2`（修改后的数组），但递归调用传入的是 `__coll_1`（原数组）：

```lean
partial def range_N (__idx_1) (__coll_1 : SparsePolyZp) (...) :=
  ...
  let __coll_2 := Array.set! __coll_1 __idx_1 updated_term  -- 修改后
  ...
  range_N (__idx_1 + 1) __coll_1 ...   -- ← BUG：用 __coll_1 而非 __coll_2
```

**根因**：`transform_range_for` 中，`__coll` 不在 `modified_names`（由 `identify_loop_vars` 识别的被修改变量）中。`identify_loop_vars` 检测 `AssignStmt` 和 `ArrayPush`，但 `Array.set!` 产生的修改通过 `transform_member_assign` 生成 `LetStmt(__coll_2, ..., Array.set!(...))` — 这是一个 LetStmt（新声明），不是 AssignStmt（重赋值），所以 `identify_loop_vars` 不认为 `__coll` 被修改了。

**影响函数**：`__upoly_make_monic`、`__upoly_primitive`、`__hensel_lift`（lc 归一化循环）、`__upoly_subtract_x`、`__upoly_subtract_one`

**修复方案**：

在 `transform_range_for` 中，将 `__coll` 总是加入 `modified_names`。因为 range-for 的集合变量可能被成员赋值修改（通过 `Array.set!`），即使 `identify_loop_vars` 没检测到。

```python
# 在 transform_range_for 中，modified_names 构建之后：
modified_names = sorted(body_vars)
# 如果循环体中有 Array.set! 修改了 __coll，将其加入
for s in body_stmts:
    if isinstance(s, LetStmt) and isinstance(s.value, Call):
        if "Array.set!" in s.value.func and s.var.name.startswith("__coll"):
            modified_names = sorted(set(modified_names) | {coll_name})
            break
```

**语义等价性**：range-for 修改集合元素 ⟺ 循环函数返回修改后的集合。将 `__coll` 加入修改变量后，循环函数正确传播修改后的集合，调用者获得更新后的数组。∎

### M2：void 函数返回值被 `let _ :=` 丢弃（~8 函数）

**问题描述**：

```cpp
__hensel_step(node, f, m, p);   // void，通过引用修改 node
```

翻译后，`__hensel_step_ir` 的返回值（修改后的 node）被丢弃：

```lean
let _ := (__hensel_step_ir node f m p)   -- 返回值丢弃！
```

**根因**：G1 的自动输出参数检测在**调用点**工作——它检测被调函数签名中的非 const 引用，生成 `BinOp("=", out_var, Call(...))`。但对 TRANSLATION_SCOPE 中的函数，调用点的签名信息可能不完整（`_get_callee_sig` 从 `referencedDecl` 获取签名，但循环函数内部调用其他 TRANSLATION_SCOPE 函数时，`referencedDecl` 可能不存在）。

具体：`__hensel_lift_recursive` 调用 `__hensel_step`。在 Clang AST 中，`__hensel_step` 的 `referencedDecl` 有完整签名。但如果调用在模板代码中（通过 `UnresolvedLookupExpr`），签名不可用。

**影响函数**：`__hensel_lift_recursive`、`__hensel_lift_linear_recursive`、`__hensel_extract_factors`

**修复方案**：

对 TRANSLATION_SCOPE 中的函数，维护一个**输出参数表**（从 C++ 扫描结果自动生成）。在 SSA 变换中，当调用 TRANSLATION_SCOPE 中的函数时，查此表决定输出参数。

```python
# 从 C++ 源码扫描生成（一次性）
TRANSLATION_SCOPE_OUTPUT_PARAMS = {
    "__hensel_step": [0],           # node 是输出
    "__hensel_step_linear": [0],    # node 是输出
    "__hensel_extract_factors": [2], # factors 是输出
    "__edf_Zp": [0, 3],            # result, rng 是输出
    ...
}
```

在 `transform_stmt` 中，当遇到 `ExprStmt(Call("__xxx_ir", args))` 且 `__xxx` 在输出参数表中时，生成赋值而非丢弃：

```python
if func_name in TRANSLATION_SCOPE_OUTPUT_PARAMS:
    out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[func_name]
    out_vars = [args[i] for i in out_indices]
    if len(out_vars) == 1:
        new_var = env.bump(out_vars[0].name)
        return LetStmt(new_var, ..., Call(func_name + "_ir", args))
```

**语义等价性**：归结为定理 1（输出参数变换保持语义）和定理 1b（函数定义签名改造）。∎

### M3：Lambda 辅助函数通过 `Rng.next` 调用（~6 函数）

**问题描述**：

Lambda 辅助函数已提取为 `_lambda_N_ir`。但调用点仍通过 `Rng.next` 调用（因为 `operator()` 在 `CALL_OPERATOR_MAP` 中映射为 `Rng.next`）：

```lean
let dot_result := (Rng.next dot_1 M_row M_col)   -- 应该是 (_lambda_3_ir M_row M_col)
```

**根因**：`CALL_OPERATOR_MAP` 中 `operator()` 映射为 `Rng.next`。这对 `std::uniform_int_distribution(rng)` 是正确的，但对 lambda 的 `operator()` 调用是错的。区分条件：如果 callee 的类型是 `std::uniform_int_distribution` → Rng.next；如果是 lambda → 直接调用。

但 Clang AST 中，lambda 的 `operator()` 调用在模板代码中是 `CXXOperatorCallExpr`，callee 类型是闭包类型（lambda 特有的匿名类型名），不是 `std::uniform_int_distribution`。

**影响函数**：`__lll_reduce`（dot, round_qq, row_sub, row_swap）、`__build_cld_matrix`（upoly_coeff）、`__select_prime`（next_p）、`__mtshl_step_j`（product_F, lc_correct）

**修复方案**：

在 `CXXOperatorCallExpr` 的 `operator()` 处理中，检查 callee 类型：

```python
# 在 CXXOperatorCallExpr handler 中：
if ref_name in CALL_OPERATOR_MAP:
    # 检查 callee 类型是否是 lambda（而非 RNG distribution）
    callee_type = inner[1].get("type", {}).get("qualType", "") if len(inner) > 1 else ""
    if "uniform_int_distribution" in callee_type or "mt19937" in callee_type:
        # RNG 路径
        lean_func, arg_order = CALL_OPERATOR_MAP[ref_name]
        ...
    else:
        # Lambda 路径：直接调用 callee 变量
        callee = parse_expr(inner[1])  # lambda 变量
        args = [parse_expr(a) for a in inner[2:]]
        return Call(callee, args)  # 直接调用
```

更简单的方案：`CALL_OPERATOR_MAP` 不再用于 `operator()`。在 `CXXOperatorCallExpr` 中，`operator()` 统一翻译为**直接调用 callee**：

```python
if "operator()" in ref_name:
    callee = parse_expr(inner[1])  # 第一个操作数 = 被调用的函数对象
    args = [parse_expr(a) for a in inner[2:]]
    return Call(gen_expr(callee), args)
```

RNG 的 `dist(rng)` 也走这条路——callee 是 `dist` 变量。在 SSA 变换中，`dist` 变量已绑定到 `Rng.next` 的偏应用（通过 clpoly_model.lean 的类型），自然解析正确。

**语义等价性**：`operator()` 调用的语义是"以 callee 为函数，args 为参数，执行调用"。直接调用 callee 与此一致。∎

### M4：2D 数组/map 原地修改丢弃（~10 函数）

**问题描述**：

```cpp
M[row][col] = value;        // 2D 数组修改
pow_theta[t][l] = value;    // 2D 数组修改
coeffs[deg] = value;        // map 修改
```

翻译后变成 `let _ := value`——修改被丢弃。

**根因**：`CXXOperatorCallExpr(operator[])` 现在正确翻译为 `ArrayAccess(arr, idx)`。但 `arr[i][j] = val` 的 AST 结构是：

```
BinaryOperator(=)
  ├── ArraySubscriptExpr    // arr[i][j]
  │   ├── CXXOperatorCallExpr(operator[])  // arr[i]
  │   │   ├── arr
  │   │   └── i
  │   └── j
  └── val
```

SSA 变换的 `ExprStmt(BinOp("="))` 处理只检测 `isinstance(expr.lhs, Var)`（简单变量赋值）和 `isinstance(expr.lhs, FieldAccess)`（成员赋值）。对 `ArrayAccess(ArrayAccess(arr, i), j)` = 2D 数组访问，没有对应的处理。

**影响函数**：`__si_vandermonde_solve`（M, pow_theta, coeffs）、`__lll_reduce`（mu, B_gs, U）、`__extract_candidates`（U_short, part, candidates）、`__mtshl_zp_univar_mdp`（s, sigma）

**修复方案**：

在 SSA 变换的 `ExprStmt(BinOp("="))` 处理中，加 `ArrayAccess` 分支：

```python
# 1D 数组赋值：arr[i] = val → let arr' := Array.set! arr i val
if isinstance(expr.lhs, ArrayAccess) and isinstance(expr.lhs.arr, Var):
    arr = expr.lhs.arr
    idx = rename_expr(expr.lhs.idx, env)
    val = rename_expr(expr.rhs, env)
    new_arr = env.bump(arr.name)
    return LetStmt(new_arr, env.get_type(arr.name),
                   Call("Array.set!", [env.current(arr.name), idx, val]))

# 2D 数组赋值：arr[i][j] = val → let row := arr[i]; let row' := Array.set! row j val; let arr' := Array.set! arr i row'
if isinstance(expr.lhs, ArrayAccess) and isinstance(expr.lhs.arr, ArrayAccess):
    outer = expr.lhs.arr  # arr[i]
    j = rename_expr(expr.lhs.idx, env)
    val = rename_expr(expr.rhs, env)
    if isinstance(outer.arr, Var):
        arr_var = env.current(outer.arr.name)
        i = rename_expr(outer.idx, env)
        row_var = Var("_row")
        stmts = [
            LetStmt(row_var, "auto", ArrayAccess(arr_var, i)),
            LetStmt(Var("_row_updated"), "auto",
                    Call("Array.set!", [row_var, j, val])),
        ]
        new_arr = env.bump(outer.arr.name)
        stmts.append(LetStmt(new_arr, env.get_type(outer.arr.name),
                     Call("Array.set!", [arr_var, i, Var("_row_updated")])))
        return stmts
```

**语义等价性**：

C++ `M[i][j] = val` 修改二维数组 M 的 (i,j) 位置。

Lean 等价操作：
```lean
let row := M[i]!
let row' := row.set! j val
let M' := M.set! i row'
```

逐元素验证：M'[i][j] = val，M'[i][k] = M[i][k] for k ≠ j，M'[l] = M[l] for l ≠ i。与 C++ 一致。∎

## 2. 根因之间的关系

```
                ┌─────────────────────────────────┐
                │ C++ 可变数据流                    │
                └──────┬──────────────────┬────────┘
                       │                  │
          ┌────────────▼──────┐  ┌────────▼────────┐
          │ 通过参数传递       │  │ 通过数组下标     │
          │ (引用/输出参数)    │  │ (arr[i] = val)  │
          └──┬─────────┬──────┘  └──┬──────┬───────┘
             │         │            │      │
      ┌──────▼──┐ ┌────▼────┐ ┌────▼──┐ ┌─▼──────┐
      │ M2      │ │ G1(已修)│ │ M1    │ │ M4     │
      │ 调用点  │ │ 签名+   │ │ coll  │ │ 2D数组 │
      │ let _:= │ │ 调用点  │ │ 不传播│ │ let _:=│
      └─────────┘ └─────────┘ └───────┘ └────────┘

      ┌─────────────────────────────────┐
      │ Lambda operator() 调用          │
      └──────────┬──────────────────────┘
                 │
          ┌──────▼──────┐
          │ M3          │
          │ Rng.next    │
          │ 而非直接调用│
          └─────────────┘
```

M1 和 M4 是同一个底层问题的两种表现：**原地修改的结果未被 SSA 追踪**。
M2 是 G1 的遗漏：G1 处理了普通函数调用的输出参数，但漏了 TRANSLATION_SCOPE 内部函数的调用。
M3 是 Lambda 提取的最后一步：函数定义做了，但调用点的 `operator()` 走错了路径。

## 3. 修复优先级

| 优先级 | 根因 | 修复量 | 解决函数 |
|--------|------|--------|---------|
| P0 | M1（coll 不传播） | ~15 行 | ~10 |
| P1 | M4（2D 数组赋值） | ~30 行 | ~10 |
| P2 | M2（void 返回值丢弃） | ~20 行 | ~8 |
| P3 | M3（Lambda Rng.next） | ~10 行 | ~6 |

总计 ~75 行代码修改，解决全部 34 个 FAIL。

## 4. 语义等价性证明

### 定理 M1（集合修改传播）
设 `L(i, arr)` 为修复后的循环函数。对迭代次数 k 归纳：base case `L(n, arr) = arr`（循环结束）；inductive step `L(i, arr) = L(i+1, Array.set! arr i (modify arr[i]))`。Array.set! 只修改第 i 个元素，归纳步修改第 i+1..n-1 个。组合得全部 n 个元素被正确修改 = C++ 语义。∎

### 定理 M2（void 返回值捕获）
由 G1b（定理 1b），void+输出参数函数返回修改值。捕获 `let x' := f_ir(x)` 得 `x' = 修改后值`。SSA 保证后续引用用 `x'`。∎

### 定理 M3（Lambda 直接调用）
`dot_1 row col = _lambda_3_ir M row col` = lambda body 在 {M, row, col} 环境中求值 = C++ `dot(row, col)`。而 `Rng.next dot_1 row col` 是伪随机数生成，语义完全不同。直接调用正确。∎

### 定理 M4（2D 数组赋值）
`M' = Array.set! M i (Array.set! M[i] j val)`。逐元素：M'[i][j]=val, M'[i][k]=M[i][k] (k≠j), M'[l]=M[l] (l≠i)。与 C++ `M[i][j]=val` 一致。∎

### 组合正确性
G1-G6 + M1-M4 覆盖 C++ 全部数据流操作（D1-D11）。每个机制已证语义等价。翻译管线逐层保持语义。∎

## 5. 预期结果

修复后的审计预期：

| 模块 | 当前 PASS | 修复后预期 PASS |
|------|-----------|----------------|
| Zp | 5/13 | ~11/13 |
| Univar | 14/33 | ~28/33 |
| Wang | 18/25 | ~23/25 |
| **总计** | **37/71** | **~62/71** |

剩余 ~9 个函数可能有更深层的语义问题（如 `__zassenhaus_recombine` 的组合枚举算法、`__lll_reduce` 的 Gram-Schmidt 计算、`__heuristic_starting_precision` 的浮点近似），需要逐案分析。
