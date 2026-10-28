# Pass 残留修复 阶段 F — 剩 12 errors 根因 + 修复方案

**日期**：2026-05-03
**前置阶段**：A/B/C/D（2026-05-02）+ E（2026-05-03 上半场，本日）
**当前状态**：lake errors 102（baseline）→ 12，0 sorry → 4 sorry
**本阶段目标**：清零剩 12 lake errors（含 4 sorry）

## 0. 范围

阶段 E 之后剩 12 errors 分布在 7 个具体根因上。本文档逐条 trace 根因 + 给出修复路径。

## 1. 错误清单

| # | 行 | 错误概要 | 根因类别 |
|---|----|---------|---------|
| 1 | 63 | Unknown identifier `term_1` | Pass 4/6/7 find-loop 默认值 |
| 2 | 276 | Unknown identifier `g_8` | 同 1 |
| 3 | 143 | `_lambda_..._ir` has type `... → ZZ` but expected `LambdaRef` | Pass 3 lifted lambda ty propagation |
| 4 | 509 | filter lambda `Array Variable → Variable × Int64 → Bool` 多 capture 未 partial-app | Pass 3 + Pass 8 capture propagation |
| 5 | 1633×2 | `dot` LambdaRef applied as function | 同 3 |
| 6 | 835 | `MvPolyZZ × Array (Variable × Int64)` vs `Poly` | Pass 2b refret destructure caller-side |
| 7 | 972 | `(Array.replicate ((#[f]).toNat) default)` | Pass 5 vector InitList wrap 错误 |
| 8 | 1407/1408 | polynomial_GCD 4-arg form | Lean Model 重载缺失 |
| 9 | 1472 | `Int.toNat (...)` 返回 Nat 但 expected Float | cast_table 错误 path |

合计 12 errors / 4 sorry（来自 filter lambda residual）。

## 2. 逐条根因 + 修法

### 2.1 错误 #1, #2 — find-loop 默认值（2 errors）

**错误形态**：

```lean
partial def _loop__lambda___build_cld_matrix_upoly_1_0_ir
    (idx : Nat) (cont : SparsePolyZZ) (deg : Int32) : (Int64 × (UMonomial × ZZ)) :=
  if idx < cont.size then
    let term_1 := cont[idx]!
    if term_1.fst.deg.toInt32 == deg then
      ((1 : Int64), term_1)        -- found
    else
      _loop_... idx+1 cont deg     -- continue
  else
    ((0 : Int64), term_1)          -- ❌ term_1 unbound
```

**根因 trace**：

C++ 源码模式（推断）：

```cpp
ZZ result = ZZ(0);
for (auto term : poly) {
    if (term.first.deg == deg) {
        return term.second;
    }
}
return result;  // not found fallback
```

Pass 4 的 RangeFor + 早返回模式识别把 loop 转成 `(exit_kind: Int64, ...live_outs)` 元组返回。Pass 6 SSA build 把 `term_1` 作为 live_out（因为 found 分支返回它），但**没有为 outer-pred 的 phi source 提供 default 值**。

具体 silent miss 位置：`pass7_loop_lower.py:892` `init_phi_srcs = _collect_phi_sources_from(header_block, init_outer_pred)`。如果 phi target 没有从 outer pred 的 source（因为 var 是 loop body 内部新生），则 source = None。Pass 7 调用 `init_args.extend(...)` 时漏了它。

更深层：在 Pass 6 SSA build 时，进入 loop header 的 phi 应该接受 outer-pred 提供的初始值。如果 var 在 outer pred 没定义（`term_1` 在 C++ 是 `auto term : poly` — 每次 iter 重新绑定，loop 外不存在），SSA 构建器应该插入一个 `default` 初始值 phi source。

**修法**：

**Option A（推荐）**：Pass 7 `init_args` 构造时检测缺失的 outer-pred phi sources，填 `default` Lean 字面量。

```python
# pass7_loop_lower.py 大约 892 行
init_phi_srcs = _collect_phi_sources_from(header_block, init_outer_pred)
# 阶段 F 修复：phi target 在 outer pred 没 source（loop body 内部 var 被 live_out
# 化）时，用 Inhabited default 兜底，避免下游 Lean unbound identifier
for ph_idx, ph_target in enumerate(_phi_targets_at(header_block)):
    if ph_idx >= len(init_phi_srcs) or init_phi_srcs[ph_idx] is None:
        init_phi_srcs[ph_idx] = Var(name="default", version=0,
                                     ty=ph_target.ty,
                                     is_lean_keyword=True)  # 直 emit 'default'
```

需要 Pass 8 emit `Var("default", is_lean_keyword=True)` 直接输出 `default` 而不查找 var 表。

**Option B**：Pass 6 SSA build 时检测 live_out 没有 outer-pred def，插入 `LetStmt v_init := default` 在 outer-pred BB 末尾。改动小但可能影响其他 Pass 的 var 重名假设。

**优先 Option A**（局部、改动量小）。

**估计**：Pass 7 + Pass 8 emit_expr 各 ~10 行。

### 2.2 错误 #3, #5 — Pass 3 lifted lambda 类型签名传播（3 errors）

**错误形态**：

```lean
-- 错误：caller param `dot : LambdaRef` (= Unit) 不能当函数用
... (dot M[idx]!)
... (some_caller _lambda___build_cld_matrix_upoly_1_ir)
   -- _lambda_... has type `SparsePolyZZ → Int32 → ZZ` 但 caller param 期望 LambdaRef
```

**根因 trace**：

Pass 3 `pass3_lambda_lift.py:327`：

```python
replacement = Var(name=lam_name, version=0, ty=NamedType("LambdaRef"))
```

调用点的 lambda（C++ `[](...) { ... }`）被 Pass 3 lift 为顶层 `_lambda_..._ir`，调用点 Var 的类型固定标为 `NamedType("LambdaRef")`。Pass 8 的 emit_type 把 `LambdaRef` emit 为 `LambdaRef`（Lean Model 中 = `Unit`）。

但 lifted func 实际类型是 `SparsePolyZZ → Int32 → ZZ`（带具体签名）。caller 拿到这个 var 时要么作为函数应用（→ 类型错），要么传给下游期望 `α → β → γ` 的位置（→ 类型错）。

**修法**：

需要在 IR 引入 `FuncType(params: list[TypeIR], ret: TypeIR)`，Pass 3 替换 Var ty 时使用该函数类型，Pass 8 emit_type 输出 Lean 函数箭头语法 `α → β → γ`。

```python
# ir_types.py 新增
@dataclass(frozen=True)
class FuncType:
    params: tuple[TypeIR, ...]
    ret: TypeIR

# pass3_lambda_lift.py:327
replacement = Var(name=lam_name, version=0,
                  ty=FuncType(params=tuple(p.ty for p in new_params),
                              ret=orig_ret))

# pass8_codegen.py emit_type
if isinstance(ty, FuncType):
    return " → ".join([_paren(emit_type(p)) for p in ty.params] + [emit_type(ty.ret)])
```

副作用：所有 Pass（4/5/6/7）需要承认 FuncType 是合法 TypeIR。具体改动小（几乎都是透传），但 grep 一遍确保 Pass 4 isinstance 检查都不假设非 FuncType。

**估计**：~30 行新 IR 节点 + 各 Pass 透传调整。

### 2.3 错误 #4 — filter lambda 多 capture 未 partial-app（1 error / 4 sorry）

**错误形态**：

```lean
-- _lambda___extract_monomial_content_lex_filt1_ir : Array Variable → Variable × Int64 → Bool
let min_deg_8 : StdMap Variable Int64 :=
  StdMap.filter min_deg_5 _lambda___extract_monomial_content_lex_filt1_ir
-- min_deg.filter 期望 (Variable × Int64 → Bool)，但 lambda 多了 `Array Variable` 参数
```

**根因 trace**：

C++ lambda：`[present](auto x) { return ...; }`，capture `present` 一个 `Array Variable`。Pass 3 把 capture `present` 提为显式 param：

```python
# 提升后签名: (present : Array Variable, x : Variable × Int64) → Bool
```

caller 调用 `min_deg.filter(my_lambda)` 时应该 partial-apply：`my_lambda present` 后才传给 filter。但 Pass 3 / Pass 8 都没生成 partial application。

之前阶段 D 试过让 emit_var_name 自动 partial-app（lifted_caps 字段），但失败：Pass 7 没把 lambda 的 caps 加进 cap_params，所以 caller scope 没 `present` 这个名字 → 撤回。

**修法**：

需要 Pass 3 在 lift lambda 时记录原始 captures（按名字），Pass 8 emit 引用时按当前 caller scope 自动 partial-app。这需要 Pass 3 把 caps 存入 HIRFunc 元数据，Pass 8 emit_var 检索：

```python
# Pass 3
new_func.captures = [c.name for c in lam.captures]  # 记录 capture names

# Pass 8 emit_var (callee 是 _lambda_)
if name.startswith("_lambda_") and name in lifted_lambdas_with_caps:
    caps = lifted_lambdas_with_caps[name]
    # 自动 partial-app: emit `(_lambda_..._ir cap1 cap2)` 而非 `_lambda_..._ir`
    cap_args = " ".join(emit_var(Var(c, ...)) for c in caps)
    return f"({lean_name} {cap_args})"
```

但 caps 的具体 version 在 caller scope 怎么找？需要查 caller 的 SSA 名字表（Pass 6 输出）。最简：Pass 8 emit_var 时拿 caller 局部 var 表查 `present` → `present_2`（最新 SSA 版本）。

风险：如果 caps 在 caller scope 实际不存在（lambda lifted from 内层 scope），fallback 不明显。

**修法 v2**：Pass 7 在 emit cap_params 时自动加上 lambda 的 caps（如果该 lambda 在 cap_params chain 中被引用）。然后 Pass 8 emit lambda ref 时按位置 partial-app cap_params。

复杂度：高。但 4 sorry 都是这个模式，修一次清 4 个。

**估计**：Pass 3 + Pass 7 + Pass 8 联动 ~50 行。

### 2.4 错误 #6 — Pass 2b refret tuple destructure 漏（1 error）

**错误形态**：

```lean
let gk_reduced_1 : Poly := (__extract_monomial_content_lex_ir gk_1 mono_factors_1)
-- RHS 类型 (MvPolyZZ × Array (Variable × Int64))，LHS 期望 Poly = MvPolyZZ
```

**根因**：

`__extract_monomial_content_lex_ir` 的 C++ 签名：

```cpp
Poly extract_monomial_content_lex(const Poly& f, vector<pair<Variable, int64_t>>& mono_factors_out);
//                                                                         ^^^ ref out
```

Pass 2b 的 refret transform：`out` 参数被收集到返回 tuple，函数返回 `(Poly × Array(Variable × Int64))`。caller 应该：

```lean
let __refret := __extract_monomial_content_lex_ir gk_1 mono_factors_1
let gk_reduced_1 : Poly := __refret.fst
let mono_factors_2 := __refret.snd
```

但 corpus 直接 `gk_reduced_1 := __extract_...`（吞了整个 tuple）。这是 Pass 2b 在某个 caller 路径上漏 destructure。

**Trace 待补**：Pass 2b 是怎么决定要不要 destructure 的。可能是 caller 自身用 `auto&` ref-arg 形式调用，Pass 2b 走了不同路径。

**修法**：Pass 2b 把所有 ref-arg → out-param 转换的 caller 都强制 destructure，无例外。

**估计**：Pass 2b ~10 行。

### 2.5 错误 #7 — Pass 5 vector InitList wrap 错误（1 error）

**错误形态**：

```lean
(Array.replicate ((#[f]).toNat) default)
-- 应该是 #[f]
```

**根因 trace 待补**：

C++ 应该是 `return {f};` 或 `vector<Poly> r{f}; return r;`。Pass 1 把 InitListExpr 解析为 ArrayLit `#[f]`。但 Pass 5 看见 `vector` 构造（arity=2: size + value），把 `#[f]` 当 size，`default` 当 value，emit `Array.replicate (#[f].toNat) default`。

可能根因：CXXConstructExpr 节点的 arg 是 InitListExpr，arity 误判为 2。

**修法**：Pass 5 对 vector 构造检测 arity=1 + arg 是 ArrayLit 时直 emit ArrayLit（透传）。

**估计**：Pass 5 ~5 行 + 测试。

### 2.6 错误 #8 — polynomial_GCD 4-arg 形态（2 errors）

**错误形态**：

```lean
let __refret_0_1 : (SparsePolyZp × SparsePolyZp × SparsePolyZp) :=
  polynomial_GCD g_zp_2 h_zp_2 s_zp_1 t_zp_1
-- polynomial_GCD : α → α → α (2-arg) 但 caller 传 4 args 期望 3-tuple
let s_zp_2 : SparsePolyZp := __refret_0_1.snd  -- snd 是 (SparsePolyZp × SparsePolyZp)
```

**根因**：

C++ EEA 形式：

```cpp
Poly polynomial_GCD(const Poly& a, const Poly& b, Poly& s_out, Poly& t_out);
//                                              ^^^ ref out  ^^^ ref out
```

Pass 2b refret 把 (s_out, t_out) 转成 tuple → 返回 `(gcd, s, t)`。

但 Lean Model 中 `polynomial_GCD` 只有 2-arg 版返回 poly。

**修法**：

(a) Lean Model 加 3-tuple 4-arg 版本：

```lean
def polynomial_GCD_eea {α : Type} [Inhabited α] (_a _b _s _t : α) : α × α × α :=
  (default, default, default)
```

(b) Pass 2b 检测 polynomial_GCD 4-arg 调用，rename callee → `polynomial_GCD_eea`。

或者更简单：在 class_map 给 polynomial_GCD 加 output_indices = [2, 3]，Pass 2b 自动 rename。看 class_map 是否支持 callee rewrite。

**估计**：~15 行。

另：caller-side tuple destructure 的 `.snd` 需要变 `.snd.fst`（`.2.1`）— 这是 Pass 2b refret 输出 tuple navigation 错（与 Pass 7 的 n>2 tuple 投影路径同样的 bug）。可能要查 `_build_refret_tuple_ty` 的 navigation 输出。

### 2.7 错误 #9 — cast_table Float ↔ Int 路径（1 error）

**错误形态**：

```lean
let a_h_d_1 : Float := (Int.toNat ((((((2.5 : Float) * ...) ...)
-- Int.toNat 返回 Nat 但 expected Float
```

**根因 trace 待补**：

C++ 应该是 `double a_h_d = ceil(...);`。Pass 5 cast_table 在某个路径上错把 Float→Int 标记为 `Int.toNat`，应该是 `Float.toInt32` 或不 cast。

**修法**：grep cast_table 找 Float / Int 路径修。

**估计**：~5 行。

## 3. 修复顺序 + 估计

| 顺序 | 错误 | 类型 | 估计行数 | 风险 |
|------|------|------|---------|------|
| 1 | #6 (835) Pass 2b destructure | 简单 | 10 | 低 |
| 2 | #7 (972) vector InitList | 简单 | 5 | 低 |
| 3 | #8 (1407/1408) polynomial_GCD | 简单 | 15 | 低 |
| 4 | #9 (1472) cast_table | 简单 | 5 | 低 |
| 5 | #1, #2 (63/276) find-loop default | 中等 | 10 | 中（Pass 7 改） |
| 6 | #3, #5 (143/1633×2) Pass 3 LambdaRef | 中等 | 30 | 中（IR 加 FuncType 节点） |
| 7 | #4 (509) filter capture partial-app | 难 | 50 | 高（Pass 3+7+8 联动，前置失败先例） |

**预计**：1-4 完成 = 5 errors 清零；5-6 完成 = 8 errors 清零；7 完成 = 12 errors 清零（含 4 sorry）。

如果 7 攻不动，可以保留 4 sorry（filter lambda） + 1 error（509），称 stage F 部分完成。

## 4. 替代方案

如果用户希望更激进：

- **完全清零路径**：按 1→7 顺序全做，估计 2-3 小时
- **务实路径**：做 1-6（清 8 errors），剩 4 sorry 标记为 known issue
- **最小路径**：做 1-4（清 5 errors），剩 7 errors + 4 sorry 都 known

按用户指令"都修掉"——走完全清零路径。

## 5. 不确定项 / 待 trace

- 错误 #7 中 Pass 5 vector ctor 怎么把 InitList 误判为 size — 需要看 Clang AST 节点
- 错误 #6 中 Pass 2b 跳过 destructure 的具体 caller 路径条件
- 错误 #9 中 cast_table 哪个 (src, tgt) 项把 Float 加了 `.toNat`
- 错误 #4 的 partial-app 之前失败的具体调试路径（lifted_caps 撤回原因）

以上 trace 在实际修复时同步进行，发现具体路径再补充本文档。

## 6. 验证

修完每条后跑：

```bash
cd proof/cpp2lean_v2 && python3 tests/build_pass8_corpus.py
cd ../lean && lake build CLPoly.Generated.Corpus
```

目标：errors=0，sorry=0。

成功后写 devlog `proof/lean/docs/devlog/2026-05-03-pass-residual-stage-F.md`。
