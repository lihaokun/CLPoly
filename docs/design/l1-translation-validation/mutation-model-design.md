# C++ → Lean 翻译器：可变数据流统一变换机制

> 状态：设计 v2（含迭代器、std::map、浮点覆盖 + 完整正确性证明）
> 覆盖：71 个 C++ 函数 → Lean IR 的全部语义鸿沟

## 1. 问题定义

### 1.1 翻译正确性目标

设 P 是 C++ 程序，T(P) 是翻译后的 Lean IR。翻译正确 ⟺

```
∀ 合法输入 x. ⟦P⟧(x) = ⟦T(P)⟧(x)
```

其中 ⟦·⟧ 是程序的指称语义（输入→输出的函数）。

### 1.2 语义鸿沟分类

审计 71 个函数发现 6 类语义鸿沟：

| 编号 | 鸿沟 | 根因 | 影响函数 |
|------|------|------|---------|
| G1 | 输出参数 | C++ 引用传递 vs Lean 值传递 | ~25 |
| G2 | 裸运算符 | 模板 ADL 延迟解析 vs Lean 即时解析 | ~20 |
| G3 | 闭包遗漏 | C++ 词法作用域 vs Lean 顶层函数 | ~10 |
| G4 | 迭代器模式 | C++ 迭代器抽象 vs Lean 无迭代器 | 4 |
| G5 | std::map | C++ 有序映射 vs Lean 无内建 map | ~10 |
| G6 | 浮点运算 | C++ double vs Lean 无浮点 | 1 |

本文为每类鸿沟提供**统一变换机制**和**语义等价性证明**。

## 2. 形式化框架

### 2.1 C++ 状态模型

C++ 程序状态 σ = (V, H)：
- V : Var → Val（栈变量）
- H : Addr → Val（堆/数组内容）

语句语义 ⟦S⟧ : State → State。

### 2.2 Lean 环境模型

Lean 程序环境 ρ : Name → Val（不可变绑定）。

表达式语义 ⟦E⟧ : Env → Val。

### 2.3 翻译等价性判据

翻译 T 正确 ⟺ 存在状态对应关系 R 使得：
- R(σ₀, ρ₀)（初始状态对应）
- ⟦S⟧(σ₀) = σ₁ ∧ ⟦T(S)⟧(ρ₀) = ρ₁ → R(σ₁, ρ₁)（每步保持对应）

R 的定义：对每个 C++ 变量 x，σ(x) = ρ(x_latest)，其中 x_latest 是 x 的最新 SSA 版本。

## 3. G1：输出参数统一变换

### 3.1 检测机制

C++ 中，非 const 引用参数 = 输出参数。Clang AST 中检测方式：

**方法 A**（从被调函数声明）：`ParmVarDecl` 的 `qualType` 含 `&` 且不含 `const`。

**方法 B**（从调用点参数）：输入参数经过 `ImplicitCastExpr(LValueToRValue)` 取值后传入；输出参数**不经过** `LValueToRValue`，直接传引用。

```python
def detect_output_params(call_node: dict) -> list[int]:
    """从 CallExpr 的参数检测输出参数位置。"""
    args = call_node.get("inner", [])[1:]  # 跳过 callee
    output_indices = []
    for i, arg in enumerate(args):
        # 输出参数：没有 LValueToRValue cast
        if not _has_lvalue_to_rvalue_cast(arg):
            # 但要排除字面量（字面量没有 LValueToRValue 是因为它不是左值）
            if _is_lvalue(arg):
                output_indices.append(i)
    return output_indices

def _has_lvalue_to_rvalue_cast(node: dict) -> bool:
    """检查节点是否经过 LValueToRValue 隐式转换。"""
    if node.get("kind") == "ImplicitCastExpr":
        if node.get("castKind") == "LValueToRValue":
            return True
    return False
```

**方法 B 比方法 A 更可靠**：不需要查找被调函数声明（在 `-ast-dump-filter` 模式下可能不可用），且对模板函数也适用。

### 3.2 变换规则

设 C++ 调用 `f(a₁, ..., aₙ)`，其中 `a_{i₁}, ..., a_{iₖ}` 是输出参数（索引集 O = {i₁,...,iₖ}）。

**规则**：

```
C++:   f(a₁, ..., aₙ);         // aᵢ (i∈O) 被修改
Lean:  let (a_{i₁}', ..., a_{iₖ}') := f_lean(a₁, ..., aₙ)
       -- 所有参数都传入（输出参数提供初始值）
       -- 返回值只含输出参数的最终值
```

当函数有非 void 返回值 r 时：
```
C++:   r = f(a₁, ..., aₙ);
Lean:  let (r', a_{i₁}', ..., a_{iₖ}') := f_lean(a₁, ..., aₙ)
```

### 3.3 语义等价性证明

**定理 1**（输出参数变换保持语义）

**形式化**：

设 C++ 函数 f 接收 n 个参数，其中 O = {i₁,...,iₖ} 是输出参数索引集。

f 的语义是一个纯函数（抽象掉引用语义后）：

```
f_pure : (v₁,...,vₙ) → (ret, w_{i₁},...,w_{iₖ})
  输入：全部 n 个参数的值（输出参数提供初始值）
  输出：返回值 ret + 输出参数的最终值 w_{iⱼ}
```

**关键**：输出参数既是输入（提供初始值，如空 vector），也是输出（接收最终值，如填充后的 vector）。f 的计算可能依赖输出参数的初始值。

C++ 调用语义：
```
σ₀ → f(σ₀(a₁),...,σ₀(aₙ)) → σ₁ = σ₀[a_{i₁} ↦ w_{i₁}, ..., a_{iₖ} ↦ w_{iₖ}]
其中 (ret, w_{i₁},...,w_{iₖ}) = f_pure(σ₀(a₁),...,σ₀(aₙ))
```

Lean 变换：
```
let (ret', a_{i₁}', ..., a_{iₖ}') := f_lean(ρ₀(a₁),...,ρ₀(aₙ))
ρ₁ = ρ₀[ret' ↦ ret, a_{i₁}' ↦ w_{i₁}, ..., a_{iₖ}' ↦ w_{iₖ}]
```

注意 f_lean **接收全部 n 个参数**（含输出参数的初始值），因为 f_pure 的计算可能依赖它们。

**前提**：调用 f 之前，SSA 对应关系 R(σ₀, ρ₀) 成立（由前序语句的翻译正确性归纳保证）。

**证明**：

需证：∀ j ∈ {1,...,k}. σ₁(a_{iⱼ}) = ρ₁(a_{iⱼ}')，且调用后 R(σ₁, ρ₁) 成立。

(1) σ₁(a_{iⱼ}) = w_{iⱼ}（C++ 赋值语义的定义）

(2) ρ₁(a_{iⱼ}') = f_lean(ρ₀(a₁),...,ρ₀(aₙ)) 的第 (j+1) 分量（Lean let 绑定的定义）

(3) 由前提 R(σ₀, ρ₀)：ρ₀(aₘ) = σ₀(aₘ) 对所有 m 成立。
    因此 f_lean(ρ₀(a₁),...,ρ₀(aₙ)) = f_lean(σ₀(a₁),...,σ₀(aₙ))

(4) 由可信基假设：f_lean 正确实现 f_pure，即
    f_lean(σ₀(a₁),...,σ₀(aₙ)) = f_pure(σ₀(a₁),...,σ₀(aₙ)) = (ret, w_{i₁},...,w_{iₖ})

(5) 由 (2)(3)(4)：ρ₁(a_{iⱼ}') = w_{iⱼ} = σ₁(a_{iⱼ})。

(6) SSA 变换保证后续代码中 a_{iⱼ} 的引用被替换为 a_{iⱼ}'（ast-dispatch-design.md §7.3 引理 2）。

(7) 调用后 R 保持：对未修改的变量 aₘ (m ∉ O)，σ₁(aₘ) = σ₀(aₘ) = ρ₀(aₘ)（C++ 不修改非输出参数）。Lean 中 aₘ 的 SSA 版本不变，ρ₁(aₘ) = ρ₀(aₘ)。因此 R(σ₁, ρ₁) 成立。

∎

### 3.4 函数定义的签名改造

§3.2 覆盖了**调用点**的变换。但当被翻译的函数**本身有输出参数**时，函数定义也需要改造。

**问题**：C++ 函数 `void __edf_Zp(vector<poly>& result, const poly& f, uint64_t d, mt19937& rng)` 翻译为 Lean 后：
- 参数 `result` 和 `rng` 是非 const 引用 → 输出参数
- C++ 返回 void，通过引用修改 result 和 rng
- 当前翻译：返回 `Unit`，修改的 result 丢失

**规则**：翻译 TRANSLATION_SCOPE 中的函数时，检测非 const 引用参数，改造返回类型。

```
C++ 签名: void f(T& out₁, const U& in₁, V& out₂)
Lean 签名: partial def f_ir (out₁ : T) (in₁ : U) (out₂ : V) : T × V :=
  -- body 中 out₁ 和 out₂ 通过 SSA 正确追踪
  -- 函数末尾返回 (out₁_latest, out₂_latest)
```

当原始返回类型非 void 时：
```
C++ 签名: R f(T& out₁, const U& in₁)
Lean 签名: partial def f_ir (out₁ : T) (in₁ : U) : R × T :=
  -- 返回 (原始返回值, out₁_latest)
```

**实现**：在 `transform_func` 中：
1. 从 `FuncIR.params` 的 `is_output` 标记识别输出参数
2. 修改 `ret_type`：`void` + 1 输出 → 输出类型；`void` + 多输出 → Prod；非 void + 输出 → Prod(ret, 输出...)
3. 在 body 末尾追加 `ReturnStmt(Prod.mk(out₁_latest, out₂_latest))`

**定理 1b**（函数定义签名改造保持语义）

设 C++ 函数 `void f(T& out, const U& in)` 的语义为：从 out 初始值 v₀ 和 in 值 u，经过函数体 S₁;S₂;...;Sₙ 执行后，out 的最终值为 v_final。

Lean 翻译：`partial def f_ir (out : T) (in_ : U) : T := let out_1 := S₁'(out, in_); ...; out_n`。

**证明**（对函数体长度 n 归纳）：

*base case (n=1)*：函数体只有一条语句 S₁。
- C++ 执行 S₁ 后 out = S₁_effect(v₀, u) = v_final。
- Lean：`out_1 = S₁'(v₀, u)`。由 S₁ 翻译正确性，`out_1 = v_final`。
- `f_ir` 返回 `out_1 = v_final`。✓

*inductive step*：假设前 k 条语句正确（Lean 中 out_k = C++ 中 out 的值），第 k+1 条 S_{k+1}：
- C++ 执行 S_{k+1} 后 out 从 v_k 变为 v_{k+1}。
- Lean：`out_{k+1} = S_{k+1}'(out_k, in_) = S_{k+1}'(v_k, u)`。
- 由 S_{k+1} 翻译正确性 + 归纳假设 out_k = v_k，得 `out_{k+1} = v_{k+1}`。

*最终*：`f_ir` 返回 `out_n = v_final`。∎

**多输出参数**：对每个输出参数独立应用上述归纳。Lean 返回 `Prod.mk(out₁_n, out₂_n)`。调用者解构。

**前提**：SSA 正确追踪每次修改（引理 2）；`is_output` 正确标识输出参数；函数体中无 C++ 生命周期副作用（透传节点引理保证）。

### 3.5 子类覆盖

| C++ 模式 | 检测 | Lean 变换 |
|---------|------|----------|
| `void f(T& out, const T& in)` | out 无 LValueToRValue | `let out' := f_lean(out, in)` |
| `R f(T& out, U in)` | out 无 LValueToRValue | `let (r', out') := f_lean(out, in)` |
| `void f(T& out1, T& out2, U in)` | out1, out2 无 LValueToRValue | `let (out1', out2') := f_lean(out1, out2, in)` |
| `void push_back(T& vec, U elem)` | CLASS_MAP 的 mutate_push | `let vec' := vec.push elem`（已有） |
| `T& operator[](size_t i)` | 返回引用，后续赋值 | ArrayAccess + Array.set!（已有） |

## 4. G2：运算符与 Lambda 统一解析

### 4.1 运算符解析规则

Clang AST 中未解析的运算符出现在两种节点中：

**情况 A**：`UnresolvedLookupExpr(name="operator-")` 作为 `CallExpr` 的 callee。

```python
# 在 CallExpr 处理中
if callee.kind == "UnresolvedLookupExpr" and "operator" in callee.name:
    op = callee.name.replace("operator", "").strip()
    args = [parse_expr(a) for a in inner[1:]]
    if len(args) == 1:
        return UnaryOp(op, args[0])
    elif len(args) == 2:
        return BinOp(op, args[0], args[1])
```

**情况 B**：`CXXDependentScopeMemberExpr(member="operator-")` 作为成员运算符。

已由现有 `handle_dependent_member` 处理（查 CLASS_MAP operators）。

### 4.2 运算符语义等价性证明

**定理 2**（运算符解析保持语义）

C++ 模板中 `UnresolvedLookupExpr(operator+)(a, b)` 在实例化时由 ADL 解析为具体类型 T 的 `T::operator+(a, b)`。

Lean 中 `BinOp("+", a, b)` 生成 `(a + b)`，由类型类 `Add T` 解析。

**证明**：分两步——(A) 解析规则正确性；(B) 运算结果等价性。

**(A) 解析规则正确性**：参数个数判断一元/二元。

Clang AST 中 `CallExpr` 的 `inner` 子节点：`inner[0]` = callee（函数/运算符引用），`inner[1:]` = 操作数。

- 一元 `operator-(a)` 的 `CallExpr`：`inner = [callee, a]`，`len(inner) - 1 = 1` → `UnaryOp`
- 二元 `operator-(a, b)` 的 `CallExpr`：`inner = [callee, a, b]`，`len(inner) - 1 = 2` → `BinOp`

Clang 保证 `inner` 的子节点个数精确等于操作数个数 + 1（callee）。因此按 `len(inner) - 1` 判断一元/二元是正确的。

**(B) 运算结果等价性**：

由可信基假设：clpoly_model.lean 中为类型 T（Zp, ZZ, SparsePolyZp, MvPolyZZ 等）定义了 `Add`/`Sub`/`Mul`/`Neg` 实例，其语义与 C++ 对应类型的运算符语义一致。

设 a : T, b : T。
- C++ 语义：`T::operator+(a, b)` = T 类型上的加法结果（由 ADL 在实例化时解析）
- Lean 语义：`Add.add a b` = `Add T` 实例的加法结果（由类型类在编译时解析）

由可信基假设，两者相等。∎

**注意**：一元运算符 `operator-` 对应 `Neg T` 实例。`UnaryOp("-", a)` 生成 `(-a)` = `Neg.neg a`。

### 4.3 Lambda 提取规则

Clang AST 中 `LambdaExpr` 的结构：

```
LambdaExpr
  ├── CXXRecordDecl
  │   ├── FieldDecl: capture x (by ref)
  │   ├── FieldDecl: capture y (by ref)
  │   └── CXXMethodDecl: operator()(ParmVarDecl a, ParmVarDecl b)
  │       └── CompoundStmt: { return a.deg < b.deg; }
  └── DeclRefExpr: x (capture init)
  └── DeclRefExpr: y (capture init)
```

**规则 A**（简单 lambda，body 是单 return）：内联。

```
[&](auto a, auto b) { return a.deg < b.deg; }
→ fun a b => a.deg < b.deg
```

**规则 B**（复杂 lambda）：提取为辅助函数。

```
[&x, &y](auto a) { S1; S2; return e; }
→ partial def _lambda_N (x y : T) (a : A) : R := S1'; S2'; e'
   -- 调用点传入 capture 变量当前值
```

### 4.4 复杂 Lambda 提取的实现规范

§4.3 规则 B 说"复杂 lambda → 提取为辅助函数"。具体实现：

**Clang AST 中 LambdaExpr 的信息提取**：

```python
def handle_lambda_expr(node):
    inner = node.get("inner", [])
    
    # 1. 提取 capture 列表
    captures = []
    for child in inner:
        if child.get("kind") == "CXXRecordDecl":
            for field in child.get("inner", []):
                if field.get("kind") == "FieldDecl":
                    captures.append(field.get("name", "_"))
    
    # 2. 提取参数列表
    params = []
    for child in inner:
        if child.get("kind") == "CXXMethodDecl":
            for p in child.get("inner", []):
                if p.get("kind") == "ParmVarDecl":
                    params.append(parse_param(p))
    
    # 3. 提取 body
    for child in inner:
        if child.get("kind") == "CXXMethodDecl":
            for sub in child.get("inner", []):
                if sub.get("kind") == "CompoundStmt":
                    body = parse_compound(sub)
    
    # 4. 简单 lambda（单 return）→ 内联
    if len(body) == 1 and isinstance(body[0], ReturnStmt):
        return body[0].value  # 内联表达式
    
    # 5. 复杂 lambda → 提取为辅助函数
    func_name = f"_lambda_{id}"
    all_params = [ParamIR(c, "auto") for c in captures] + params
    aux_func = SSAFunc(name=func_name, params=all_params, body=body, ...)
    ctx.aux_defs.append(aux_func)
    # 返回偏应用：传入 capture 变量
    return Call(func_name + "_ir", [Var(c) for c in captures])
```

**关键**：capture 列表从 `CXXRecordDecl → FieldDecl` 提取。这些是 lambda 闭包类型的字段，每个对应一个捕获变量。

**可变 capture**：如果 capture 变量在 lambda 内被修改（如 `[&idx]() { idx[s]++; }`），该 capture 同时是输出参数——归结为 §3.4 的签名改造。

### 4.5 复杂 Lambda 提取语义等价性证明

**定理 3b**（复杂 Lambda 提取保持语义）

C++ lambda `[&x₁, &x₂, y₃](A a, B b) { S₁;...;Sₙ; return e; }` 在环境 Γ = {x₁=σ(x₁), x₂=σ(x₂), y₃=y₃_captured} ∪ {a=v_a, b=v_b} 中求值 body。

Lean 提取函数 `_lambda_N(x₁, x₂, y₃, a, b)` 在环境 {x₁=ρ(x₁_cur),..., a=v_a, b=v_b} 中求值 body'。

**证明**：
(1) 由 SSA 对应关系 R：ρ(x₁_cur) = σ(x₁)（引用捕获 = 调用时的值）。
(2) 值捕获 y₃：ρ(y₃_captured) = 定义时的值。由 SSA 版本不变性一致。
(3) Lean 环境 = C++ 环境（(1)+(2)）。
(4) body 翻译正确性由归纳假设保证。
(5) 返回值在相同环境中求值 → 相等。

**可变捕获**：x₁ 在 lambda 内被修改 → x₁ 是输出参数（§3.4）→ 返回修改值 → 由定理 1b 保持语义。

**capture 完备性**：Clang AST 的 `CXXRecordDecl → FieldDecl` 穷举了全部 capture 变量（Clang 编译器保证）。∎

### 4.6 原始 Lambda 语义等价性证明（简单情况）

**定理 3**（Lambda 提取保持语义）

需要区分两种捕获模式：

**情况 A：只读捕获**（lambda 不修改捕获变量）

C++ lambda `[&x₁,...,&xₘ](a₁,...,aₖ) { body_readonly }` 的语义：
- 创建闭包 C = ⟨body, {x₁=σ(x₁),...,xₘ=σ(xₘ)}⟩
- 调用 C(v₁,...,vₖ) = eval(body, {a₁=v₁,...,aₖ=vₖ, x₁=σ(x₁),...,xₘ=σ(xₘ)})

Lean 提取函数 `_lambda_N(x₁,...,xₘ, a₁,...,aₖ) := body'` 的语义：
- 调用 `_lambda_N(ρ(x₁'),...,ρ(xₘ'), v₁,...,vₖ) = eval(body', {x₁=ρ(x₁'),..., a₁=v₁,...})`

**证明 A**：

(1) 由 SSA 对应关系 R，ρ(xᵢ') = σ(xᵢ)。
(2) 因此 Lean 参数环境 = C++ 闭包环境。
(3) body 的翻译正确性由归纳假设保证。
(4) 结果相等。∎

**情况 B：可变捕获**（lambda 修改捕获变量，如 `[&idx](){ idx[s]++; ... }`）

C++ 语义：lambda 通过引用修改 `idx`。调用后 σ'(idx) ≠ σ(idx)。

Lean 变换：可变捕获变量同时作为输入和输出（与 G1 输出参数相同机制）。

```
C++:  bool ok = lambda(args...);  // lambda 修改 idx
Lean: let (ok', idx') := _lambda_N(idx, args...)
```

**证明 B**：归结为定理 1（输出参数变换）。可变捕获变量 = lambda 函数的输出参数。由定理 1，idx' = C++ 调用后 idx 的值。∎

**情况 C：值捕获 `[x]`**

C++ 语义：捕获定义时 x 的值（拷贝）。

Lean 语义：传入定义点的 `x_current` SSA 版本。

由 SSA 构造，`x_current` = 定义时 x 的值。后续即使 x 被修改（x 的新版本），传给 lambda 的仍是旧版本。等价。∎

## 5. G3：闭包变量完备性

### 5.1 后处理 pass

在循环函数完全构造后，运行 `collect_free_vars` 覆盖整个函数体（含递归调用）。

```python
def finalize_loop_func(loop_func: SSAFunc, outer_env: VarEnv):
    """确保循环函数无自由变量。"""
    param_names = {p.name for p in loop_func.params}
    free = collect_free_vars(loop_func.body)
    
    missing = []
    for vname, vver in sorted(free):
        lean_name = Var(vname, vver).lean_name()
        if lean_name not in param_names:
            missing.append((vname, vver, lean_name))
    
    for vname, vver, lean_name in missing:
        typ = outer_env.get_type(vname)
        # 追加参数
        loop_func.params.append(ParamIR(lean_name, typ))
        # 更新所有调用点
        _add_arg_to_all_calls(loop_func, lean_name, Var(vname, vver))
```

### 5.2 `collect_free_vars` 的正确性

**引理 4a**（`collect_free_vars` 正确计算自由变量集）

`collect_free_vars(stmts)` 返回 `Referenced(stmts) \ Declared(stmts)`，其中：
- `Referenced(stmts)` = stmts 中所有 Var 引用的 (name, version) 集合
- `Declared(stmts)` = stmts 中所有 LetStmt 声明的 (name, version) 集合

**证明**：`collect_free_vars` 递归遍历所有 IR 节点类型：

| IR 节点 | 遍历方式 | 覆盖 |
|---------|---------|------|
| `Var(name, ver)` | 加入 Referenced | ✓ |
| `Lit` | 无变量引用 | ✓ |
| `BinOp(op, lhs, rhs)` | 递归 lhs, rhs | ✓ |
| `UnaryOp(op, operand)` | 递归 operand | ✓ |
| `Call(func, args)` | 递归 args | ✓ |
| `ArrayAccess(arr, idx)` | 递归 arr, idx | ✓ |
| `FieldAccess(obj, field)` | 递归 obj | ✓ |
| `ArrayPush(arr, elem)` | 递归 arr, elem | ✓ |
| `Cast(expr, _, _)` | 递归 expr | ✓ |
| `CondExpr(c, t, e)` | 递归 c, t, e | ✓ |
| `BlockExpr(stmts, value)` | 递归 stmts（via scan_stmt）+ 递归 value | ✓ |
| `LetStmt(var, typ, val)` | var 加入 Declared，递归 val | ✓ |
| `IfStmt(cond, then, else)` | 递归 cond, then 子句, else 子句 | ✓ |
| `ReturnStmt(val)` | 递归 val | ✓ |
| `ExprStmt(expr)` | 递归 expr | ✓ |

穷举覆盖了全部 IR 节点类型（ir_types.py 中的 ExprIR 和 StmtIR union）。∎

### 5.3 完备性证明

**定理 4**（闭包后处理保证无自由变量）

设 FV(B) = `collect_free_vars`(B) 为函数体 B 的自由变量集（引理 4a 保证正确），P 为参数名集。

后处理：P' = P ∪ {v ∈ FV(B) | v ∉ P}。

则 FV(B) ⊆ P'。

**证明**：设 v ∈ FV(B)。
- 若 v ∈ P，则 v ∈ P'（P ⊆ P'）。
- 若 v ∉ P，则 v ∈ FV(B) \ P。由后处理定义，FV(B) \ P ⊆ P'。故 v ∈ P'。

因此 FV(B) \ P' = ∅，循环函数体中没有自由变量。∎

### 5.3 语义等价性证明

**定理 5**（追加闭包参数保持语义）

设 `loop(x₁,...,xₙ)` 的 body 引用自由变量 z（来自外部作用域，值为 v_z）。

追加后：`loop(x₁,...,xₙ, z)` 在调用时传入 `z = v_z`。

C++ 语义：body 通过词法作用域读 z = v_z。
Lean 语义：body 通过参数 z 读 z = v_z。

两者等价（相同的值以不同方式传入）。

递归调用中 z 不变（只读），因此递归传入的 z 始终 = v_z。∎

## 6. G4：迭代器模式变换

### 6.1 问题分析

CLPoly 中迭代器模式只有一种用途：**遍历 + 原地压缩 + 尾部截断**。

```cpp
auto it = f.data().begin();
auto out = f.data().begin();
for (; it != f.data().end(); ++it) {
    fdiv_r(it->second, it->second, m);  // 修改系数
    if (it->second != 0) {
        *out = *it;
        ++out;
    }
}
f.data().erase(out, f.data().end());  // 截断尾部
```

语义：对每个元素应用变换，保留非零元素。等价于 `filter_map`。

### 6.2 变换规则

**规则**：迭代器压缩模式 → `Array.filterMap`。

```lean
let f' := f.filterMap (fun term =>
    let coeff' := ZZ.fdiv_r term.snd m
    if coeff' != 0 then some { term with snd := coeff' }
    else none)
```

### 6.3 语义等价性证明

**定理 6**（迭代器压缩 → filterMap 等价）

C++ 迭代器压缩模式的三步操作：

**步骤 1（变换）**：`fdiv_r(it->second, it->second, m)` — 原地修改 `it->second`。
设变换函数 `t(e) = fdiv_r(e.second, m)`，产生新元素 `e' = {e.first, t(e)}`。

**步骤 2（过滤）**：`if (it->second != 0)` — 保留 `e'.second ≠ 0` 的元素。
设谓词 `p(e') = (e'.second ≠ 0)`。

**步骤 3（压缩）**：`*out = *it; ++out; erase(out, end)` — 将满足条件的元素紧密排列。

三步的组合语义：
```
result = [t(e) | e ∈ arr, p(t(e))]
       = [{e.first, fdiv_r(e.second, m)} | e ∈ arr, fdiv_r(e.second, m) ≠ 0]
```

Lean `Array.filterMap`：
```
result = arr.filterMap (fun e =>
    let coeff' := ZZ.fdiv_r e.snd m
    if coeff' != 0 then some {e with snd := coeff'} else none)
```

展开 `filterMap` 定义：对每个 `e ∈ arr`，计算 `coeff' = fdiv_r(e.snd, m)`，若 `coeff' ≠ 0` 则保留 `{e with snd := coeff'}`，否则丢弃。

这与 C++ 的三步组合语义逐元素相同。

**关于原地修改 vs 新数组**：C++ 修改原数组 `f`（原地压缩），Lean 创建新数组 `f'`。由 SSA 变换，后续代码引用 `f'`（新版本），不引用 `f`（旧版本）。因此等价。

**严格证明**：对数组 `arr = [e₁, e₂, ..., eₙ]`，设 `eᵢ' = t(eᵢ)`。

C++ 结果：`[eᵢ' | i ∈ {1..n}, p(eᵢ')]`（保序——迭代器从前往后）。
Lean 结果：`arr.filterMap(t_and_p)` = `[eᵢ' | i ∈ {1..n}, p(eᵢ')]`（filterMap 保序）。

两者是相同的有序列表。∎

### 6.4 检测机制

迭代器压缩模式有固定的 AST 结构：

```
ForStmt
  ├── DeclStmt: it = f.data().begin()
  ├── BinaryOperator(!=): it != f.data().end()
  ├── UnaryOperator(++): ++it
  └── CompoundStmt
      ├── [transform step]
      ├── IfStmt(predicate)
      │   ├── [*out = *it; ++out]
      └── (no else)
```

后跟 `f.data().erase(out, f.data().end())`。

**检测规则**：连续出现 `begin/end` 迭代器初始化 + for 循环 + `erase(out, end)` → 识别为压缩模式。

```python
def detect_compact_pattern(stmts: list) -> bool:
    """检测迭代器压缩模式。"""
    # 查找 begin/end 声明 + for 循环 + erase 调用的连续序列
    ...
```

### 6.5 函数体整体替换

§6.2-6.4 覆盖了**内联压缩模式**（出现在其他函数体中的 begin/for/erase 序列）。但有些函数的**整个函数体就是迭代器模式**，逐行翻译不可行。

**问题**：`__upoly_mod_coeff` 的整个函数体是：
```cpp
auto it = f.data().begin();
auto out = it;
for (; it != f.data().end(); ++it) {
    fdiv_r(it->second, it->second, m);
    if (it->second) { if (out != it) *out = std::move(*it); ++out; }
}
f.data().erase(out, f.data().end());
```

逐行翻译产生迭代器垃圾（`operator->`、`operator++`、`operator*`）。

**规则**：对整个函数体是不可翻译模式的函数，提供**函数体覆盖表**——不翻译原始 body，用预定义的语义等价 Lean body 替换。

```python
FUNC_BODY_OVERRIDE = {
    # func_name → (param_names, lean_body_expr)
    "__upoly_mod_coeff": "SparsePolyZZ.modCoeff f m",
    "__upoly_divmod_mod": "SparsePolyZZ.divmod_mod q r f g m",
}
```

**定理 6c**（`__upoly_mod_coeff` 体替换保持语义）

设 C++ `__upoly_mod_coeff(f, m)` 的语义为：对 f 中每个 (mono, coeff)，计算 c = fdiv_r(coeff, m)，保留 c ≠ 0 的元素。

Lean `SparsePolyZZ.modCoeff f m` 定义为 `f.filterMap(fun (mono, coeff) => let c := coeff % m; if c ≠ 0 then some (mono, c) else none)`。

**证明**（逐元素）：

对第 i 个元素 (monoᵢ, coeffᵢ)：
(a) C++ 计算 cᵢ = fdiv_r(coeffᵢ, m)。当 m > 0 时，fdiv_r(a, m) = a - m·⌊a/m⌋ ∈ [0, m)。
(b) Lean 计算 cᵢ = coeffᵢ % m。当 m > 0 时，Int.emod(a, m) = a - m·⌊a/m⌋ ∈ [0, m)。
(c) 当 m > 0 时两者相等（floor division 余数 = Euclidean modulo for positive divisor）。
(d) 过滤条件相同（cᵢ ≠ 0）。
(e) 保留元素相同（(monoᵢ, cᵢ)）。
(f) 遍历顺序相同（前到后，保序）。

前提：m > 0（m 是素数的幂，由算法不变量保证）。∎

**与 L1 1:1 原则的关系**：这不是"简化"——迭代器压缩和 filterMap 是**同一算法的两种表示**（命令式 vs 函数式）。C++ 选择迭代器是因为性能（原地修改），Lean 选择 filterMap 是因为纯函数式不支持原地修改。两者计算的是同一个函数。

**影响范围**：

| 函数 | 覆盖方式 | 理由 |
|------|---------|------|
| `__upoly_mod_coeff` | `modCoeff` | 整个函数体是压缩模式 |
| `__upoly_divmod_mod` | `divmod_mod` | 整个函数体是归并模式 |

其他函数（`__hensel_step`、`__hensel_step_linear`）的内联压缩由 §6.4 的模式检测处理。

### 6.6 两种迭代器模式

CLPoly 中迭代器有两种不同的使用模式：

**模式 A：压缩（filterMap）**——遍历 + 变换 + 过滤零。定理 6 覆盖。

| 函数 | 用途 |
|------|------|
| `__upoly_mod_coeff` | 系数模约化 + 零压缩 |
| `__hensel_step` | e/ep 零压缩 |
| `__hensel_step_linear` | e 零压缩 |

3 个函数，模式完全相同。

**模式 B：归并（merge）**——两个有序列表按键归并。

`__upoly_divmod_mod` 的内层循环：
```cpp
while (r_it != r.end() || g_it != g.end()) {
    if (g_it == g.end()) { *out++ = *r_it++; }
    else if (r_it == r.end()) { ... }
    else if (deg_r > deg_g) { *out++ = *r_it++; }
    else if (deg_r < deg_g) { ... }
    else { /* deg_r == deg_g */ coeff = r + g; ... }
}
```

变换：两个有序数组的归并 → `List.merge` 或自定义递归函数。

```lean
partial def merge_polys (r g : SparsePolyZp) : SparsePolyZp :=
    if r.isEmpty && g.isEmpty then #[]
    else if g.isEmpty then r
    else if r.isEmpty then g
    else if r.front!.fst.deg > g.front!.fst.deg then
        #[r.front!] ++ merge_polys r.tail g
    else if r.front!.fst.deg < g.front!.fst.deg then
        #[g.front!] ++ merge_polys r g.tail
    else  -- equal degree: add coefficients
        let c := r.front!.snd + g.front!.snd
        if c != 0 then #[(r.front!.fst, c)] ++ merge_polys r.tail g.tail
        else merge_polys r.tail g.tail
```

**定理 6b**（归并等价）：`merge_polys` 的输出是 r 和 g 的按度数降序归并，同度数系数相加，零系数过滤。与 C++ 的迭代器归并循环逐元素等价（因为两者都按相同的度数比较规则处理，且保序）。∎

1 个函数使用此模式。

## 7. G5：std::map 变换

### 7.1 CLPoly 中的 std::map 使用模式

| 用途 | 键类型 | 值类型 | 例子 |
|------|-------|--------|------|
| 变量→值映射 | `variable` | `ZZ` | `eval_point`, `alpha`, `sub` |
| 度数→系数映射 | `int64_t` | `Zp` 或 `Poly` | `acc`, `coeffs` |
| 度数→索引分组 | `int64_t` | `vector<int>` | `groups` |

### 7.2 Lean 模型

在 clpoly_model.lean 中定义 StdMap 为 `List (K × V)`（有序对列表，按键排序）：

```lean
abbrev StdMap (K V : Type) := List (K × V)

namespace StdMap
def empty : StdMap K V := []
def find? (m : StdMap K V) (k : K) [BEq K] : Option V :=
  m.find? (fun (k', _) => k == k') |>.map Prod.snd
def insert (m : StdMap K V) (k : K) (v : V) [BEq K] : StdMap K V :=
  (k, v) :: m.filter (fun (k', _) => !(k == k'))
def erase (m : StdMap K V) (k : K) [BEq K] : StdMap K V :=
  m.filter (fun (k', _) => !(k == k'))
end StdMap
```

### 7.3 语义等价性证明

**定理 7**（StdMap 模型正确性）

std::map<K,V> 是键到值的有限映射，支持 find/insert/erase/operator[]。

**抽象语义**：将 std::map 视为数学函数 M : K → Option V。

| 操作 | 抽象语义 |
|------|---------|
| `find(k)` | 返回 M(k)（Some v 或 None） |
| `insert(k, v)` | M' = M[k ↦ Some v] |
| `erase(k)` | M' = M[k ↦ None] |
| `operator[](k)` 读 | 返回 M(k).getOrDefault |
| `operator[](k) = v` 写 | M' = M[k ↦ Some v] |

**Lean 实现**：`List (K × V)` 视为键值对列表。

```
find?(m, k) = (m.find? (fun (k', _) => k == k')).map Prod.snd
insert(m, k, v) = (k, v) :: m.filter (fun (k', _) => k ≠ k')
erase(m, k) = m.filter (fun (k', _) => k ≠ k')
```

**证明**（对每个操作验证语义一致性）：

*find?*：
- 若 (k, v) ∈ m，则 `m.find?` 返回 `some (k, v)`，`.map Prod.snd` = `some v`。与抽象语义 M(k) = Some v 一致。
- 若 k ∉ m，则 `m.find?` 返回 `none`。与 M(k) = None 一致。

*insert*：
- `filter` 移除所有旧 (k, _)。然后 `(k, v) :: ...` 添加新绑定。
- 结果：对 k 查找得到 v（新绑定在列表头），对 k' ≠ k 查找结果不变（filter 不影响）。
- 与 M[k ↦ Some v] 一致。

*erase*：
- `filter` 移除所有 (k, _)。对 k 查找得到 None，对 k' ≠ k 不变。
- 与 M[k ↦ None] 一致。

**键序问题**：C++ std::map 按键有序（支持 lower_bound、rbegin 等有序操作）。Lean List 无内在顺序。CLPoly 中 map 的有序遍历仅出现在：
- `__pseudo_remainder_x1`：按 x1 度数降序遍历系数 → 变换时显式排序
- `__si_theta_array_eval`：按度数分组 → 插入后排序

对这些场景，在需要有序遍历时调用 `List.mergeSort`。对不需要有序的场景（find/insert/erase），List 实现直接正确。∎

### 7.4 operator[] 的特殊处理

C++ `map[key]` 在键不存在时**插入默认值并返回引用**。这有副作用（修改 map）。

**读取变换**：`map[key]`（已知键存在）→ `StdMap.find! map key`

**前提条件**：键存在由算法不变量保证。若不成立，C++ 会插入默认值（改变 map 大小），Lean 会 panic。这属于 UB-freedom 假设——算法正确时不会访问不存在的键。

**写入变换**：`map[key] = value` → `let map' := StdMap.insert map key value`

写入总是安全的（insert 覆盖已有键或添加新键）。∎

## 8. G6：浮点运算近似

### 8.1 影响范围

仅 `__heuristic_starting_precision` 一个函数。C++ 代码：

```cpp
double logp  = std::log((double)p);
double a_h_d = std::ceil((2.5 * r + min_b) * std::log(2.0) / logp
                         + std::log((double)(N + 1)) / (2.0 * logp));
```

此函数计算 Hensel 提升的精度估计。结果只用于决定提升多少步，不影响因式分解的正确性（只影响性能——精度不够时有 fallback 重新提升）。

### 8.2 近似规则

| C++ | Lean 近似 | 方向 | 说明 |
|-----|----------|------|------|
| `std::log(x)` | `Nat.log 2 x` | 下界 | `Nat.log 2 x ≤ log₂(x) ≤ Nat.log 2 x + 1` |
| `std::ceil(x)` | `x + 1` | 上界 | `x + 1 ≥ ⌈x⌉` 对所有 x |
| `2.5 * r` | `(5 * r + 1) / 2` | 上界 | 向上取整：`⌈5r/2⌉ ≥ 2.5r` |
| `log(2.0)` | `1`（整数近似） | 上界 | `1 ≥ log(2) ≈ 0.693` |

### 8.3 正确性论证

**定理 8**（浮点近似保持算法正确性）

C++ 公式：
```
a_h = ⌈(2.5r + min_b) · log(2)/log(p) + log(N+1)/(2·log(p))⌉
```

Lean 近似（每一步取保守上界）：
```
a_h' = ((5r + 1)/2 + min_b) / Nat.log(2, p) + (Nat.log(2, N+1) + 1) / (2 · Nat.log(2, p)) + 1
```

**方向性分析**：

(1) `(5r + 1)/2 ≥ 2.5r`（向上取整）→ 分子更大 → 比率更大。

(2) `1/Nat.log(2,p)` vs `log(2)/log(p)`：
   - `Nat.log(2,p) = ⌊log₂(p)⌋ ≤ log₂(p) = log(p)/log(2)`
   - 因此 `1/Nat.log(2,p) ≥ log(2)/log(p)`
   - 作为乘法因子，更大 → a_h' 更大。

(3) 最外层 `+1` 代替 `ceil` → 上界。

综合：**a_h' ≥ a_h**。

**算法影响**：

a_h 决定 Hensel 提升的初始精度。a_h' ≥ a_h 意味着 Lean 版本可能多提升几步。

多提升不影响正确性，因为：
- Hensel 提升保持因子的模同余性（不管提升多少步，因子始终是 f mod p^k 的正确分解）
- 后续的 Mignotte bound 精度 a_mig 独立于 a_h 计算
- 算法取 `min(a_h, a_mig)` 作为起点——如果 a_h' 过大，a_mig 路径会被选择

**关键性质**：`__heuristic_starting_precision` 的输出仅影响**第一次 Hensel 提升的精度**。如果精度不足，van Hoeij 算法检测到后会 fallback 到 Mignotte 精度重新提升。因此即使近似不精确，算法仍然正确——只可能多做一轮提升。

因此浮点近似不影响算法的正确性保证。∎

**注意**：如果 `Nat.log(2, p) = 0`（p = 0 或 p = 1），会导致除零。这由 `p` 是素数的前提条件排除（素数 ≥ 2，`Nat.log(2, p) ≥ 1`）。

## 9. 完整性定理

### 9.1 C++ 数据流操作的穷举分类

CLPoly 因式分解代码中的全部数据流操作：

| 编号 | 操作 | 变换机制 | 状态 |
|------|------|---------|------|
| D1 | 局部变量赋值 `x = e` | SSA let 链 | 已有 |
| D2 | 函数输出参数 `f(out, in)` | **G1：自动检测 + 返回值** | 本设计 |
| D3 | 成员赋值 `obj.field = e` | SSA functional update | 已有 |
| D4 | 数组元素赋值 `arr[i] = e` | SSA Array.set! | 已有 |
| D5 | 数组元素成员赋值 `arr[i].field = e` | SSA Array.set! + functional update | 已有 |
| D6 | 自增/自减 `i++` | SSA `let i' := i + 1` | 已有 |
| D7 | 复合赋值 `a += b` | SSA `let a' := a + b` | 已有 |
| D8 | 方法突变 `vec.push_back(x)` | CLASS_MAP mutate | 已有 |
| D9 | 迭代器压缩 `begin/end/erase` | **G4：filterMap** | 本设计 |
| D10 | map 操作 `map[k] = v` | **G5：StdMap 模型** | 本设计 |
| D11 | 浮点计算 `log/ceil` | **G6：整数近似** | 本设计 |

### 9.2 表达式节点的穷举分类

| 编号 | Clang AST 节点 | 变换机制 | 状态 |
|------|---------------|---------|------|
| E1 | 字面量（Int/Bool/Float/String） | 直接翻译 | 已有 |
| E2 | 变量引用 DeclRefExpr | SSA rename | 已有 |
| E3 | 二元运算 BinaryOperator | BinOp | 已有 |
| E4 | 一元运算 UnaryOperator | UnaryOp | 已有 |
| E5 | 函数调用 CallExpr | Call + FUNC_MAP | 已有 |
| E6 | 成员调用 CXXMemberCallExpr | CLASS_MAP | 已有 |
| E7 | 运算符调用 CXXOperatorCallExpr | OPERATOR_MAP + operator[] | 已有 |
| E8 | 类型转换 ImplicitCastExpr | Cast | 已有 |
| E9 | 构造 CXXConstructExpr | CLASS_MAP constructors | 已有 |
| E10 | 模板成员 CXXDependentScopeMemberExpr | CLASS_MAP 查找 | 已有 |
| E11 | 模板函数 UnresolvedLookupExpr | **G2：参数数推断** | 本设计 |
| E12 | Lambda LambdaExpr | **G2：内联或提取** | 本设计 |
| E13 | 结构化绑定 DecompositionDecl | Prod 解构 | 已有 |
| E14 | 下标 ArraySubscriptExpr | ArrayAccess | 已有 |
| E15 | 条件 ConditionalOperator | CondExpr | 已有 |

### 9.3 控制流的穷举分类

| 编号 | 控制流 | 变换机制 | 状态 |
|------|-------|---------|------|
| C1 | if/else | IfStmt + SSA phi | 已有 |
| C2 | for 循环 | loop-as-function | 已有 |
| C3 | while 循环 | loop-as-function | 已有 |
| C4 | range-for | loop-as-function | 已有 |
| C5 | do-while | loop-as-function（同 while） | 已有 |
| C6 | break | ReturnStmt(state) | 已有 |
| C7 | continue | ReturnStmt(recurse) | 已有 |
| C8 | return | ReturnStmt | 已有 |
| C9 | throw | Throw | 已有 |

### 9.4 完整性定理

**定理 9**（翻译完整性）

CLPoly 因式分解代码中的全部数据流操作（D1-D11）、表达式节点（E1-E15）、控制流（C1-C9）均有对应的变换机制。

**证明**：

1. D1-D11 穷举了 C++ 中所有修改状态的方式。D9-D11 是本设计新增。
2. E1-E15 穷举了 Clang AST 中出现的所有表达式节点（由 §5 审计脚本验证）。E11-E12 是本设计新增。
3. C1-C9 穷举了 C++ 中所有控制流构造。

每种均有变换机制，且各机制的语义等价性已证明（定理 1-8 + ast-dispatch-design.md 中的引理 1-5）。∎

### 9.5 组合正确性定理

**定理 10**（翻译正确性）

在以下假设下，翻译 T 保持语义等价：

1. **可信基正确**：clpoly_model.lean 中类型和操作的实现与 C++ 一致
2. **AST 完整覆盖**：审计脚本确认所有节点类型已处理
3. **无 UB**：C++ 执行无未定义行为（由 require 义务保证）

**证明**：

翻译管线三层：
```
C++ AST ──AST 解析──→ IR ──SSA 变换──→ SSA IR ──codegen──→ Lean text
```

第一层：每个 AST 节点有注册的 handler（完整性由审计保证），每个 handler 保持语义（E1-E15 各自正确）。透传节点不改变语义（定理 6，ast-dispatch-design.md §7.2 引理 1）。

第二层：
- SSA 变量重命名保持语义（ast-dispatch-design.md §7.3 引理 2）
- 循环提取保持语义（loop-extraction-design.md §4）
- 输出参数变换保持语义（定理 1）
- break/continue 替换保持语义（ast-dispatch-design.md §4）
- 闭包参数追加保持语义（定理 5）

第三层：codegen 是 IR 的直接文本渲染（注入映射），不改变语义。

三层组合：每层保持语义 → 整体保持语义。∎

## 10. 实现步骤

### 已完成
1. ✓ **G1 调用点**：自动检测输出参数（从函数签名解析非 const 引用）+ 单/多输出解构
2. ✓ **G2 运算符**：CXXOperatorCallExpr 一元/二元 + CallExpr+UnresolvedLookupExpr 运算符
3. ✓ **G3 闭包**：`finalize_loop_func` 后处理 pass
4. ✓ **G4 内联压缩**：`_detect_iterator_compact` 检测 + modCoeff/compactNonzero 替换
5. ✓ **G5 StdMap**：clpoly_model.lean 定义 + CLASS_MAP 方法
6. ✓ **G6 浮点**：Nat.log + 保守 ceil

### 待完成
7. **G1 函数定义签名改造**（§3.4）：`transform_func` 中检测输出参数 → 修改返回类型 + 追加 ReturnStmt
8. **G4 函数体整体替换**（§6.5）：`FUNC_BODY_OVERRIDE` 表 → `__upoly_mod_coeff` 和 `__upoly_divmod_mod` 体替换
9. **G2 Lambda 复杂提取**（§4.4）：从 CXXRecordDecl 提取 capture + 生成辅助函数
10. **G2 残留 operator-**：检查 CXXDependentScopeMemberExpr 路径的运算符
11. 全量翻译 + 逐函数对比审计
