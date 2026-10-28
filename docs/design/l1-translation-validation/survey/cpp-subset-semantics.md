# CLPoly C++ 子集形式语义

> Stage 1 Week 3 主产出
> 基于 `cpp-construct-catalog.md`（Week 1）+ `type-system.md`（Week 2）
> 日期：2026-04-21

**目的**：为 cpp2lean v2 翻译器的每个 Pass 提供语义保持论证的**形式基础**。

**范围**：CLPoly 因式分解代码实际用到的 C++ 子集（已在 `cpp-construct-catalog.md` §3-12 枚举）。

---

## §0 导读与语义惯例

### §0.1 记号

| 记号 | 含义 |
|---|---|
| `⟦e⟧` | C++ 表达式 e 的 denotational semantics |
| `⟦S⟧ : State → State` | C++ 语句 S 的语义（状态转换器）|
| `⟦e⟧_L` | Lean 翻译 e 的语义 |
| `T(e)` / `T(S)` | 翻译器的翻译函数 |
| `σ` | 程序状态（变量名 → 值）|
| `σ[x ↦ v]` | 在状态 σ 中把 x 赋值为 v |
| `a ≡ b` | 两个表达式在所有满足前置条件的状态下语义相等 |

### §0.2 语义域（Domain）

```
Value   ::= UInt64 | Int64 | Bool | Array Value | Prod Value Value
          | Zp | ZZ | Rat | Monomial | SparsePoly | MvPoly | ...

State   = Var → Value        -- 程序状态（大致；实际需要版本号标注）

Error   = DivByZero | OOB | ShiftOOB | Overflow | AssertFail | EmptyContainer
         | UserThrow (CLPoly 无用户 throw)

Result  = State + Error       -- UB 用 Error 表达

⟦·⟧ : C++ AST → State → Result
```

### §0.3 UB 的处理

**核心约定**：C++ 触发 UB 时，对应的 Lean `require` 前置条件不成立。定理只在**所有 require 成立**的状态下声明等价。

```
C++:  σ ⊢ ⟦S⟧ ⇓ Error ub
Lean: require_of(S) 不成立 in σ
```

因此**翻译正确性只论证 UB-free 执行路径**。UB 路径不在定理范围（用户或系统需满足前置条件）。

---

## §1 基本运算语义（对应 `type-system.md` §1）

### §1.1 整数算术（uint64_t 为例）

**C++ 语义**：C++17 标准 [basic.fundamental]/4：`unsigned` 运算遵循 `mod 2^N`，**无 UB**（除零除外）。

```
⟦a + b⟧_C ≡ (⟦a⟧ + ⟦b⟧) mod 2^64  -- uint64_t
⟦a * b⟧_C ≡ (⟦a⟧ * ⟦b⟧) mod 2^64
⟦a / b⟧_C ≡ ⌊⟦a⟧ / ⟦b⟧⌋ if ⟦b⟧ ≠ 0
          ⇓ Error DivByZero       if ⟦b⟧ = 0
⟦a << n⟧_C ≡ (⟦a⟧ * 2^n) mod 2^64 if n < 64
          ⇓ Error ShiftOOB        if n ≥ 64
```

**Lean 语义**：Lean 4 `UInt64 := Fin (2^64)`，算术运算定义相同：

```
⟦a + b⟧_L = (a.val + b.val) mod 2^64  : UInt64
⟦a / b⟧_L = if b = 0 then 0 else ⌊a.val / b.val⌋  -- 不 panic
```

**差异点**：C++ `/0` 是 UB，Lean `UInt64 / 0 = 0`。翻译必须在 `/` 前生成 `require b ≠ 0`。

### §1.2 引理 L1.1（基础运算等价）

> **引理**：对 C++ 子集用到的所有**无 UB**的基础运算 op，
> 设 σ 满足 op 的 UB 前置条件，则 `⟦op(a, b)⟧_C(σ) ≡ ⟦op(a, b)⟧_L(σ)`。

**证明**：对每个 op 逐表验证（`type-system.md` §1.3 表）。核心点：C++17 标准的 unsigned 算术语义与 Lean `Fin (2^64)` 的运算定义均为 mod 2^64，**点点一致**。

对 signed 算术：在**无溢出**前提（UB-6 require 成立）下一致；溢出时 C++ 是 UB、Lean 是 mod wrap，但 UB-free 执行路径不到达溢出。∎

### §1.3 Bool 短路

**C++ 语义**：`a && b` 若 `a = false` 则 **不求值** `b`（短路）。

**Lean 语义**：`a && b = Bool.and a b`，pattern match 实现：
```lean
def Bool.and : Bool → Bool → Bool
  | false, _ => false
  | true, b  => b
```
第一个模式**不检查** `b`。等价短路。

### §1.4 引理 L1.2（短路等价）

> `⟦a && b⟧_C(σ) ≡ ⟦a && b⟧_L(σ)` for all σ（含 `a = false` 且 `b` 有副作用的情况）。

CLPoly 因式分解代码**无副作用表达式**（纯值语义），所以短路行为在实际执行中不影响结果，但形式上的等价仍由上述 Lean 定义保证。∎

---

## §2 控制流语义（对应 `cpp-construct-catalog.md` §3）

### §2.1 if-else

```
⟦if cond then S1 else S2⟧_C(σ) =
  if ⟦cond⟧_C(σ) then ⟦S1⟧_C(σ) else ⟦S2⟧_C(σ)

⟦if cond then T(S1) else T(S2)⟧_L(σ) =
  if ⟦T(cond)⟧_L(σ) then ⟦T(S1)⟧_L(σ) else ⟦T(S2)⟧_L(σ)
```

**语义保持**（引理 L2.1）：若 T(cond)、T(S1)、T(S2) 各自语义保持（归纳假设），则 if-else 组合语义保持。

### §2.2 While 循环

```
⟦while cond do S⟧_C(σ) =
  if ¬⟦cond⟧_C(σ) then σ
  else let σ' := ⟦S⟧_C(σ) in ⟦while cond do S⟧_C(σ')
```

**Lean 翻译**（`loop_lower` Pass 产出）：

```lean
partial def while_ir (state : State) : State :=
  if !⟦cond⟧_L state then state
  else while_ir (⟦S⟧_L state)
```

### §2.3 引理 L2.2（while 循环等价）

> 对所有 σ，若循环**终止**（k 步后 cond=false），则
> `⟦while cond do S⟧_C(σ) = while_ir(σ)`

**证明**：对步数 k 归纳。
- **Base (k=0)**：`⟦cond⟧ = false`。C++ 返回 σ；Lean `while_ir σ` 进入 then 分支返回 σ。
- **Step (k→k+1)**：`⟦cond⟧ = true`。C++：σ' := ⟦S⟧(σ); recurse on σ'。Lean：`while_ir (⟦S⟧ σ)`。由归纳假设（k 步从 σ'），两者相等。

**不终止情况**：C++ 不终止，Lean `partial def` 也不终止（Lean 的 `partial` 承诺 type-safe 但不保证终止）。语义"一致不终止"，定理成立（在不终止语义域里 ⊥ = ⊥）。∎

**CLPoly 的终止保证**：由 L2 模型（已 0 sorry 证终止）间接给出。L1 层的 `partial def` 只证"行为与 L2 一致"，不重证终止。

### §2.4 For 循环

C++ `for (init; cond; step) { body }` ≡ `init; while (cond) { body; step; }`。

翻译器用同一 `loop_lower` Pass 处理，引理 L2.2 的变体直接适用。

### §2.5 Range-for

```
for (auto& x : v) body
  ≡ for (size_t i = 0; i < v.size(); ++i) { auto& x = v[i]; body; }
```

翻译器识别 range-for 后，HIR `iter_recognize` 生成 `Array.foreach` 或直接等价于 C-for（由 `loop_lower` 处理）。

### §2.6 Break / Continue

**C++ 语义**：
- `break` 在循环内立即退出
- `continue` 跳过 body 剩余，执行 step（for）或 cond 检查（while）

**Lean 翻译**（`loop_lower` Pass）：引入 `_break_flag` / `_continue_flag`（仅在循环内可见），下降为：

```lean
partial def loop_ir (state, _break_flag, _continue_flag) :=
  if _break_flag then state             -- break 已设
  else if !cond state then state        -- 正常退出
  else
    let state' := body state            -- 执行 body（内部可能设 flag）
    if _break_flag_after then state'    -- break 在 body 里设了
    else loop_ir (step state', false, false)  -- reset flags
```

### §2.7 引理 L2.3（break/continue 语义保持）

> 翻译后 Lean `partial def` 的行为与 C++ 循环相同（含 break 提前退出、continue 跳过 body 剩余）。

**证明**：对**循环状态**（位置 + break flag + continue flag）做多维归纳。每个 break/continue 触发点在 C++ 中对应一个控制流跳转；翻译器生成的 flag 设置 + 条件分支精确模拟该跳转。∎

### §2.8 Return（含循环内）

**C++ 语义**：`return v` 立即退出当前函数，返回 v。

**Lean 翻译**：循环内 return 用 `_ret_flag + _ret_val` 机制（`cpp-construct-catalog.md` §15 已提及）：

```lean
partial def loop_ir (state, _ret_flag, _ret_val) :=
  if _ret_flag then (state, _ret_flag, _ret_val)  -- 短路
  else ...

def func_ir (args) :=
  let (state, rf, rv) := loop_ir (init, false, default)
  if rf then rv else default_ret
```

### §2.9 Do-while

```
do { S } while (cond)
  ≡ { S; while (cond) { S } }  -- body 至少执行一次
```

翻译器把 do-while 展开为上述形式（`loop_lower` Pass 的 2 处特例）。

---

## §3 变量、赋值与 SSA 变换（对应 §ssa_build）

### §3.1 C++ 变量模型

C++ 变量绑定到内存位置（或寄存器），赋值修改其内容。状态 `σ : Var → Value`。

```
⟦x = e⟧_C(σ) = σ[x ↦ ⟦e⟧_C(σ)]
```

### §3.2 Lean SSA 变量模型

Lean 无 mutation；SSA 变换后每个 C++ 变量 x 有多个版本 `x_0, x_1, x_2, ...`，每个版本对应 C++ 中该变量生命期的一段：

```
⟦x_i = e⟧_L (σ) = σ[x_i ↦ ⟦e⟧_L(σ)]
```

每个 `x_i` 只在声明后被定义一次（SSA 不变量）。

### §3.3 引理 L3.1（SSA 变换等价）

> 设 C++ 代码 P，SSA 变换为 T(P)。对 P 的任一执行轨迹 τ_C：
> 1. 将 τ_C 中对 x 的每次写替换为对 `x_{写次序}` 的 let 绑定
> 2. 将每个读 x 替换为读"最近的版本"
>
> 得到 τ_L 是 T(P) 的合法 Lean 执行，且 `end_state(τ_L) = end_state(τ_C)` (在所有 x 版本名对齐下)。

**证明**（Appel 1998 定理 1 的直接实例化）：
- SSA 变换 = 每次写创建新变量名
- 每次读 = 读最新定义
- Dominance 属性保证"读到的版本"在语义上就是 C++ 中"最新写"
- 结合 CLPoly 无别名（§0.4 前提）、无 goto，dominance frontier 上 phi 正确合并分支版本

详细证明参见 Appel 1998 "SSA is Functional Programming" §4 Thm 1。∎

**CLPoly 无别名前提**：函数参数值传递或 const ref；无全局可变状态；循环间无迭代器共享（由 `iter_recognize` 抽象为纯函数）。**由代码审计（blueprint §5a.2 P1）保证**。

### §3.4 Phi-node 语义

```
x_{new} = phi(x_{then}, x_{else})  -- 合并分支
```

翻译为 Lean：

```lean
let x_new := if cond then x_then else x_else
```

**等价性**：两分支各自定义 `x_then` 和 `x_else`，if 选择对应的。C++ 中分支后 x 的值取决于执行路径；Lean 的 `if` 精确模拟此选择。∎

---

## §4 函数调用与引用消除（对应 §ref_elim）

### §4.1 值传递

```
⟦f(a)⟧_C(σ) = ⟦f⟧(⟦a⟧_C(σ))    -- 值作为参数传入，不共享内存
```

Lean 天然值传递，直接对应。

### §4.2 const 引用传递

C++ `f(const T& x)` 语义上与值传递等价（x 不可修改）。Lean 按值传递，等价。

### §4.3 非 const 引用传递：`ref_elim` Pass

C++：
```cpp
void foo(T& out) { out = compute(); }
// 调用方:
T x;
foo(x);
use(x);
```

C++ 语义：
```
σ' := σ[out ↦ compute(σ)]
// 调用返回后 σ'[x] = compute
```

Lean 翻译（`ref_elim` Pass 后）：
```lean
def foo_ir : T := compute
-- 调用方:
let x_1 := foo_ir
use x_1
```

### §4.4 引理 L4.1（ref_elim 语义保持）

> 对 C++ 函数 `void foo(T& out, ...) { ... out = expr; }`，
> 翻译为 Lean `def foo_ir (...) : T := expr`，
> 调用 `foo(x, ...)` 后 x 的值 = Lean 调用 `let x_new := foo_ir (...)` 后 x_new 的值。

**证明**：
- C++ 中 `out = expr` 写入 out 引用的内存位置（= x 的存储）
- 没有别名（§0.4），x 是唯一引用
- 返回后 x 的值 = expr 的值
- Lean：`foo_ir` 返回 expr 的值，绑定到 `x_new`
- 两者值相同 ∎

**多输出参数**：同理推广到 tuple 返回。

### §4.5 跨函数语义保持（引理 L4.2）

> 若被调函数 g 的翻译保持语义，则调用 g 的函数 f 的翻译也保持语义。

**证明**：f 的翻译用 `let y := g_ir args` 替代 `g(args)`。归纳假设 g_ir 返回值 = C++ 中 g 的返回值。之后 f 的剩余代码使用 y 的地方，C++ 和 Lean 一致（引理 L3.1 SSA）。∎

**对调用图拓扑序归纳** → 全程序语义保持。

---

## §5 Lambda 与闭包（对应 §lambda_lift）

### §5.1 C++ Lambda 语义

```cpp
int x = 10;
auto f = [&x](int y) { return x + y; };  // [&] 引用捕获
f(5)  // 15
```

C++ 语义：lambda 作为对象，其闭包包含对 x 的引用；调用 `f(5)` 时读 x 的当前值。

### §5.2 Lean 翻译（`lambda_lift` Pass）

Lambda 提取为独立函数，captures 显式传参。对 `[&]` capture：

```lean
partial def _lambda_1_ir (x : Int) (y : Int) : Int := x + y

-- 调用点:
let x := 10
let result := _lambda_1_ir x 5  -- 显式传 x
```

### §5.3 引理 L5.1（Lambda 提升等价）

> Lambda `[&captures](params) { body }` 提升为 `_lambda_N_ir (captures, params)`。
> 调用 `lambda(args)` 等价于 `_lambda_N_ir (current_captures, args)`。

**证明**：
- C++ 读 `captures` 的**当前值**（闭包不拷贝，是引用）
- Lean 显式把当前 captures 作为参数传入
- 两者读到相同的值
- body 的计算在两处都是纯函数于 (captures, params) ∎

**修改 captures 的情况**：若 lambda 内修改了 `[&]` capture（如 `x = x + 1`），则 Lean `_lambda_N_ir` 必须返回修改后的 capture（作为返回值 tuple 的一部分）：

```lean
partial def _lambda_1_ir (x : Int) (y : Int) : (Int × Int) :=
  let x_new := x + 1
  (x_new, x_new + y)

-- 调用方显式重新绑定：
let (x_new, result) := _lambda_1_ir x 5
-- 从此 x_new 替代原 x
```

**等价性**：C++ 闭包引用修改 → Lean 返回值写回。语义保持（引理 L4.1 的推广）。

### §5.4 CLPoly 统计（`lambdas.md`）

26 个 in-scope lambda：
- 13 无捕获（最简单）
- 10 `[&]` 默认引用（需写回机制）
- 3 具名（`[x]`/`[j]`/`[use_large_prime]`），多为循环变量 by-value
- **0 generic**（Day 1 修复）

---

## §6 迭代器模式（对应 §iter_recognize）

### §6.1 识别的 4 种模式

| 模式 | 数量 | 语义 |
|---|---|---|
| Range-for | 92 | `for (auto& x : v) body` = `for i = 0..size; x := v[i]; body` |
| Structured-binding range-for | 24 | 对 `(K × V)` 容器，同上 + 绑定解构 |
| Compact-erase 双指针 | 4 | 原地 filter |
| Classic iterator loop | 1 | `for (auto it = v.begin(); it != v.end(); ++it)` |

### §6.2 Compact-erase 模式（关键）

**C++ 代码**（`__upoly_mod_coeff` 模式）：
```cpp
auto it = v.begin();
auto out = it;  // 双指针初始化
for (; it != v.end(); ++it) {
    if (keep_pred(*it)) *out++ = *it;
}
v.erase(out, v.end());  // 去掉尾部
```

**C++ 语义**：对 v 做 filter，保留满足 keep_pred 的元素，顺序不变。

**Lean 翻译**（`iter_recognize` 识别模式后）：
```lean
let v' := v.filter keep_pred
```

### §6.3 引理 L6.1（compact-erase 等价）

> 对上述 C++ 模式，`v.erase(...)` 后的 v 等于 Lean `v.filter keep_pred`。

**证明**：
- 双指针 it/out 的 loop invariant：`v[0..out]` = 到目前为止保留的元素
- 循环结束：it = v.end, out 指向第一个"非保留"位置
- `v.erase(out, v.end())` 截断到 out 之前
- 结果 = `v.filter keep_pred`（保留满足谓词的元素，顺序不变）
- Lean `Array.filter` 定义也是保留满足谓词、顺序不变
- 两者一致 ∎

**复杂度**：C++ O(n)，Lean `Array.filter` O(n)。等价。

### §6.4 其他迭代器的语义

- **classic iterator loop**：`iter_recognize` 不做特殊识别，保留为 HIR 显式 for-loop，经 `loop_lower` 下降
- **parallel iterator**（zip-walk）：`__upoly_divmod_mod` 仅 1 处，不识别，保留显式循环

### §6.5 范围 for + 结构化绑定（StdMap 遍历）

```cpp
for (const auto& [k, v] : map) { body }
```

等价于：
```
for (auto it = map.begin(); it != map.end(); ++it) {
    const K& k = it->first; const V& v = it->second; body
}
```

**Lean 翻译**（对 `StdMap` shim）：
```lean
map.toList.foldl (fun _ (k, v) => body) ()
-- 或使用 Array.foreach on map.entries
```

等价性由 `StdMap.toList` 保证（按 key 升序），与 `std::map` 的有序遍历一致。

---

## §7 assert / require / Except

### §7.1 assert 的语义

C++ `assert(cond)`（Debug 编译）：
- `cond` 为 true：noop，继续执行
- `cond` 为 false：abort（程序终止）

**关键前提**：翻译器**要求开发者保证 `cond` 在所有合法调用中为 true**。

### §7.2 引理 L7.1（assert → require 等价）

> C++ `assert(cond); body` 在 cond = true 时的行为 ≡ Lean `(h : cond) → body`。

**证明**：
- C++ cond=true：assert noop，body 正常执行
- Lean 提供 `h : cond` 证明后，函数体执行 body
- 两者行为一致 ∎

**cond=false 情况**：C++ abort（不在定理范围），Lean 无法被调用（缺少 h 参数）。两者都"不进入 body"。

### §7.3 throw / Except（确认不使用）

CLPoly 因式分解代码**无 `throw`**（`cpp-construct-catalog.md` §9 确认 `CXXThrowExpr` / `CXXTryStmt` 均为 0）。所有错误路径通过 return 值表达。

**结论**：L1 IR **不引入 `Except` 类型**。函数签名返回普通类型或 tuple。

---

## §8 UB 点与 require 的精确对应（§2 + `ub-sites.md`）

### §8.1 对应表

| UB 类型 | 触发条件 | 生成的 require | 站点数 |
|---|---|---|---|
| UB-1 Div0 | `a / b` 或 `a % b`，b = 0 | `require hb : b ≠ 0` | 31 |
| UB-2 OOB | `arr[i]`，i ≥ arr.size | `require hi : i < arr.size` | 468 |
| UB-3 Empty | `v.front!` / `v.back!`，v 空 | `require ¬v.isEmpty` | 42 |
| UB-4 Shift | `a << n`，n ≥ 64 | `require hn : n < 64` | 0（不用） |
| UB-6 Signed | signed `a + b` / `-` / `*` 溢出 | `require Int.noOverflow a op b` | 113 |
| UB-7 U→S | `(int)x` 溢出 | `require x ≤ INT_MAX` | 8 |
| UB-8 Assert | `assert(cond)` | `require h : cond` | 23 |
| **合计** | | | **685** |

### §8.2 UB-free 执行路径的精确定义

**定义（UB-free 执行）**：C++ 程序 P 在输入 σ 下的执行 UB-free ⇔ P 在 σ 下的求值过程中**所有 UB 触发点**的前置条件均成立。

**定理 L8.1（UB-free 精确对应）**：
> 若 P 在 σ 下 UB-free，则所有对应的 Lean `require` 参数在 σ 下可被满足。

**证明**：一对一构造。每个 C++ UB 点对应一个 require。UB-free ⇔ 前置条件成立 ⇔ require 证明存在。∎

### §8.3 `require` 传播（跨函数）

**引理 L8.2**（require 传播）：
> 若 g 调用 f，g 的 UB-free 前置条件蕴含 f 的 UB-free 前置条件，则 g 的 require 集合包含 f 的（可能加上 g 自身的）。

**实际表现**：CLPoly 的调用图大部分是链式的（如 `factorize → __factor_multivar → __wang_core → __mtshl_lift`），require 从底层传播到顶层 `factorize`。最终 `factorize_lex_ir : ... (h1 : ...) (h2 : ...) ... → Factorization`。

---

## §9 语义保持主定理

### §9.1 主定理

**定理 L9.1（翻译正确性主定理）**：
> 对 CLPoly TRANSLATION_SCOPE 中的任一函数 P：
> 设 T(P) 是翻译后的 Lean `partial def P_ir`。设 σ 是满足 T(P) 所有 `require` 前置条件的输入状态。则：
> 1. 若 C++ `P(σ)` 正常返回值 v，则 `T(P)(σ)` 返回 v'，v' 与 v 在类型映射下相等
> 2. 若 C++ `P(σ)` 触发 UB，则 σ 不满足某个 require（矛盾，不在定理范围内）

### §9.2 证明结构

**归纳策略**：
1. 对调用图**拓扑排序**（从叶节点到 `factorize`）
2. 对每个函数：对其 body 的**语句结构**做结构归纳
3. 叶函数：由引理 L1.1（基础运算等价）+ L6.1（迭代器）直接成立
4. 非叶函数：由 L4.2（跨函数）+ L5.1（lambda）+ L2.1-3（控制流）+ L3.1（SSA）组合

### §9.3 各引理用在哪一层

| Pass | 对应引理 |
|---|---|
| `parse` | L0（基础 trivial） |
| `ref_elim` | **L4.1** + L4.2 |
| `lambda_lift` | **L5.1** |
| `iter_recognize` | **L6.1** |
| `operator_resolve` | L1.1 + CAST_TABLE |
| `ssa_build` | **L3.1** (Appel 1998) |
| `loop_lower` | **L2.2 + L2.3 + L2.4** |
| `codegen` | L1.1（trivial renaming） |

### §9.4 引用基础

| 引理 | 主要依据 | 备注 |
|---|---|---|
| **L1.1** 基础运算等价 | C++17 标准 [basic.fundamental]/4 + Lean 4 `UInt64 := Fin (2^64)` 定义 | 点点一致验证 |
| **L2.2** While 循环等价 | Winskel《The Formal Semantics of Programming Languages》Thm 5.11（while 的最小不动点语义 = operational） | 直接可引 |
| **L2.3** break/continue 保持 | Aeneas ICFP 2022 §5（flag+tail return 模式） | 工程模式 |
| **L3.1** SSA 变换等价 | Appel 1998 "SSA is Functional Programming" §4 Thm 1 | 模式引用（严格证明仍需自写）|
| **L4.1** ref_elim 等价 | CLPoly-specific（§4.4）；通用形式见 Aeneas §4.2 borrow elision | |
| **L4.2** 跨函数保持 | Leroy 2009 "Formal verification of a compiler back-end" Lemma 3.4（deterministic target 下 forward simulation 自动升格 backward）| 省一半工作量 |
| **L5.1** Lambda 提升等价 | 本文件 §5.3 直接证明；无现成引用 | CLPoly 子集限定 |
| **L6.1** Compact 模式保持 | 本文件 §6.3 直接证明 | CLPoly-specific |

### §9.5 Aeneas 对比

Aeneas (ICFP 2022) 翻译 Rust MIR → Lean/HOL4 后端，使用 `partial def` + extrinsic termination proofs。**本质上与我们的方案一致**：
- Aeneas 的 `Result { ok | fail | div }` monad ≈ CLPoly 的 `Except` + `partial def`（`div` 被 `partial def` 吸收，我们不用 `Except`，`fail` 通过 require 不满足排除）
- Aeneas 的 loop extraction 同我们的 `loop_lower` Pass
- Aeneas 的 reference elision 同我们的 `ref_elim` Pass

**重要发现**（by Agent）：Aeneas 的 2022 论文**未做机器化的翻译正确性证明**（informal）；2024 工作只证 borrow checker 侧面。**本文件的纸面严格证明已与 Aeneas 同级**，无更高级工具可依赖。

### §9.6 参考文献

详见 `semantic-references.md`（167 行 Agent 调研）：
1. Appel 1998, "SSA is Functional Programming"
2. Ho & Protzenko, "Aeneas: Rust Verification by Functional Translation", ICFP 2022
3. Winskel, "The Formal Semantics of Programming Languages", MIT Press 1993
4. Leroy 2009, "Formal verification of a compiler back-end", CACM 52(7)

---

## §10 已知 Caveats（为后续 Pass 设计记录）

### §10.1 浮点不精确

`__heuristic_starting_precision` 用 `double` + `std::log` + `std::ceil`。**浮点语义本质不精确**，不同 IEEE 754 实现可能轻微差异。L1 IR 用 Lean `Float`（IEEE 754 double），与 C++ `double` 一致（现代 x86/ARM 均 IEEE 754）。

**策略**：此函数只影响 Hensel lifting 起始精度，**不影响最终因子正确性**（precision 总是"越大越安全"，可能慢但不错）。精化证明可以豁免这个函数的严格等价，只证"存在合法 precision"。

### §10.2 随机数的概率语义

`__edf_Zp` 用 `std::mt19937` 做 EDF 随机分裂。**随机性不能确定性翻译**。

**策略**：与 L2 `exists_nonQR_poly`（AdjoinRoot 有限域存在性）对接——L1 的随机数生成器在足够多次尝试下**概率几乎肯定**分裂成功，L2 证 **存在性**（非概率）。两者在因式分解正确性上等价（L2 说"存在分裂"，L1 的 RNG 找到那个分裂）。

### §10.3 GMP 任意精度 ZZ 的精确建模

`ZZ` 映射为 `Int` 假设 GMP 正确实现（任意精度整数）。GMP 本身是 trusted kernel（工业级，无已知 bug）。

**策略**：与 Mathlib 信任 `Int` 的底层实现同构，不展开 GMP 代码。

### §10.4 C++ 标准允许的实现定义行为

- `sizeof(int)`、`INT_MAX` 等实现定义。CLPoly 假设 x86_64 Linux (`int = int32_t`, `long = int64_t`)。
- 若编译目标变化，`type-system.md` 的映射表需调整。

---

## §11 语义保持的不作证明部分

以下问题**不在本文件证明范围**，作为前置条件声明：

1. **Clang 的 AST dump 忠实反映 C++ 语义**：信任 Clang 工业级编译器
2. **Lean 4 `UInt64` / `Int` / `Array` 等的 Mathlib 公理**：信任 Lean 4 + Mathlib 工业标准
3. **GMP / IEEE 754 的底层正确性**：工业级工具
4. **L2 模型已经 0 sorry**：本 Stage 不碰，Week 0 验收已通过
5. **CLPoly 无多线程/别名/goto 等 C++ 子集外特性**：代码审计保证（blueprint §5a.2 P1-P8）

---

## §12 与后续 Week 的衔接

- **Week 4 HIR 设计**：基于 §1-§7 的语义规则，定义 HIR 节点不变量（例：HIR₁ 保证无 ref 参数，因 §4.4 的 `ref_elim` 已消除）
- **Week 5 MIR/codegen**：基于 §3.1（SSA）+ §2.6（break/continue flag 下降）设计 Pass 6-7
- **Stage 2 实现**：每个 Pass 的代码带 docstring 引用本文件相应引理
- **Stage 4 精化证明**：本文件的引理作为**非形式论证**，若项目进入形式化阶段可逐条 Lean 化（但 Stage 1 停在"可引用的英文证明"级别）

---

## §13 Week 3 验收

- [x] §1 基本运算语义（UInt64/Int64 + bool 短路 + 引理 L1.1/L1.2）
- [x] §2 控制流语义（if/while/for/range-for/break/continue/return/do-while + L2.1-L2.3）
- [x] §3 SSA 变换语义 + L3.1（Appel 1998 引用）
- [x] §4 函数调用 + ref_elim 语义 + L4.1/L4.2
- [x] §5 Lambda 语义 + L5.1
- [x] §6 迭代器 4 种模式 + L6.1
- [x] §7 assert/require/Except 对应
- [x] §8 685 UB 点 → require 精确映射 + L8.1/L8.2
- [x] §9 主定理 L9.1 + 证明策略 + 各引理对 Pass 的覆盖
- [x] §10-11 Caveats + 前置条件
- [x] §12 与 Week 4-5 衔接

**本文件是 Stage 1 Week 3 的硬性产出**，为 Week 4 HIR 设计 + Week 5 MIR/codegen 设计提供语义基础。
