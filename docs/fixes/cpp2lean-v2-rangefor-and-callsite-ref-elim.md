# 修正方案：cpp2lean v2 RangeFor `auto&` 写回 + 调用点 ref-elim

> 状态：草稿（待用户确认）
> 对应 workflow.md §5.1 修正方案文档
> 触发：Pass 6 1-对-1 审视第二轮（Agent B 语义模拟），2026-04-27

## 摘要

cpp2lean v2 Pass 6 1-对-1 审视暴露 **2 个 P0 + 1 个 P1 silent semantic bug**，且与之同根的 5 类模式可能在其他函数也存在。问题不在 Pass 6 本身，而在 Pass 1/2/4/5 协同的 mutation 模型未完整落实——`mutation-model-design.md` §3 已设计好框架（G1 输出参数变换 + ref param 检测 + 函数签名改造），但仅落实了**被调函数定义改造**（Pass 2 §3.4），**未落实调用点改写**（§3.2）；`for (auto& x : c)` 的引用语义在 Pass 6 RangeFor 真展开时也丢失。

修复后预计：
- P0-1：`__factor_multivar`、`__build_cld_matrix` 等 ~6 个函数的 ranged-for 引用变量写回正常；
- P0-2：`__edf_Zp` 递归调用、`__upoly_make_monic` 等 19 个被注册的 ref-out 函数调用点正确解构；
- 衍生：消除当前 6 个 `__sideeff_*=__write__(...)` 残余 + 14 个 mutating call sideeff。

---

## 第一部分：复现与定位

### P0-1：ranged-for `auto&` 写回丢失

**最小复现**（`__factor_multivar` bb43）：

```cpp
// C++ 源
for (auto& term : gk_pos.data()) {
    term.second = -term.second;
}
```

**当前 MIR₀ 输出**（`/tmp/mir0_dump/__factor_multivar.txt:217-218`）：
```
let term_1 := __rangefor_cont_4_1[__rangefor_idx_4_2]
let term_2 := __write__(term_1.second, (-term_1.second))
→ goto bb_latch
```

**bug**：`term_2` 是 record-update 后的 term（B7 修复正常），但 `__rangefor_cont_4_1` **没有被回写**——下一次循环迭代或 loop-exit 后读 `gk_pos` 仍是原值。C++ 通过 `auto&` reference 的写入丢失。

**等价对照**：C++ `for (auto& term : c) term.second = -term.second` 应等价于 `for (i = 0; i < c.size(); ++i) c[i].second = -c[i].second`，但当前 MIR 等价于 `for (i = 0; i < c.size(); ++i) { auto term = c[i]; term.second = -term.second; }`——`c` 不变。

**影响范围**（grep 命中）：

| 函数 | 行号 | 模式 | 影响 |
|------|------|------|------|
| `__factor_multivar` | bb43, bb64 | `term.second = -term.second` | 因子系数取负丢失，结果错误 |
| `__build_cld_matrix` | bb13 | `row.push_back(0)` | 矩阵行追加丢失 |
| `__hensel_lift_recursive` | 多处 | `__hensel_step(nodes[idx], ...)` | 节点状态更新丢失 |
| `__factor_recombine` | 待查 | `factor.coeffs[i] = ...` | 因子重组写入丢失 |
| `__lll_reduce` | 待查 | `M[i][j] = ...` (用 ranged-for) | LLL 迭代失效 |
| `__si_theta_array_eval` | 待查 | 输出累加 | 求值结果丢失 |

### P0-2：调用点未解 pair-return

**最小复现**（`__edf_Zp` bb23）：

```cpp
// C++ 源（已被 Pass 2 ref-elim 改造的签名）
__edf_Zp(result, g_8, d, rng);  // result, rng 是 ref out
```

**Pass 2 改造后的 callee 签名**：
```
__edf_Zp(result : Vec<Poly>, f : Poly, d : Nat, rng : Rng) → (Vec<Poly>, Rng)
```

**当前 MIR₀ 输出**（`/tmp/mir0_dump/__edf_Zp.txt:97-98`）：
```
let __sideeff_23_5_1 := __edf_Zp(result, g_8, d, rng)
let __sideeff_23_6_1 := __edf_Zp(result, h_part_2, d, rng)
```

**bug**：调用点把 pair-return 装到 `__sideeff_*` Var 然后丢弃。bb23 终结子返回原 `(result, rng)`（参数初值），递归调用对 result/rng 的累加全丢。

**等价对照**：应为
```
let (result_1, rng_1) := __edf_Zp(result_0, g_8, d, rng_0)
let (result_2, rng_2) := __edf_Zp(result_1, h_part_2, d, rng_1)
... return (result_2, rng_2)
```

**影响范围**：`class_map.py:684 TRANSLATION_SCOPE_OUTPUT_PARAMS` 已注册 19 个 ref-out 函数。这 19 个函数的所有调用点都是 broken。

| 注册的 ref-out 函数 | 输出 param 索引 |
|------|---|
| `__upoly_make_monic` | [0] |
| `__upoly_mod_coeff` | [0] |
| `__upoly_divmod_mod` | [0, 1] |
| `__edf_Zp` | [0, 3] |
| `__hensel_tree_build_recursive` | [0] |
| `__hensel_step` | [0] |
| `__hensel_extract_factors` | [2] |
| `__hensel_lift_recursive` | [0] |
| `__hensel_step_linear` | [0] |
| `__hensel_lift_linear_recursive` | [0] |
| `__build_cld_matrix` | [0] |
| `__lll_reduce` | [0, 1] |
| `__si_vandermonde_solve` | [2] |
| `__si_theta_array_eval` | [6] |
| `__mtshl_zp_univar_mdp` | [2] |
| `__mtshl_multi_bdp` | [3, 8] |
| `__mtshl_sparse_int` | [3, 9] |
| `__mtshl_wmds` | [3, 8] |
| `__mtshl_step_j` | [3, 5] |

### P1：mutating method/free-function call 不 bump SSA

**最小复现**（`__ddf_Zp` bb9）：
```
let __sideeff_9_0_1 := __upoly_make_monic(gd_2)   ← gd_2 应被改为 monic 形态
... 后续 line 50 仍读 gd_2 ←—— 用的是非 monic 旧值
```

**bug**：`__upoly_make_monic` 在 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 中注册为 `[0]`（mutates first arg），但调用点没生成 `gd_3 := __upoly_make_monic(gd_2)`——同 P0-2 根因。

P1 实际上是 P0-2 的子集（同根因，同修复方案）。

---

## 第二部分：根因分析

### 设计 vs 实现的对照

| 设计文档 §位置 | 设计内容 | 实现状态 |
|---|---|---|
| `mutation-model-design.md` §3.1 | 检测 ref out 参数（LValueToRValue 缺位法） | ✓ Pass 1 落实 (`is_ref`) |
| `mutation-model-design.md` §3.2 | **调用点变换**：`f(a)` → `let (a') := f(a)` | ✗ **未落实** |
| `mutation-model-design.md` §3.4 | 函数定义签名改造（pair-return） | ✓ Pass 2 落实 |
| `mutation-model-design.md` §3.5 | `push_back` 等 method 已通过 CLASS_MAP mutate 路径 | ✓ Pass 5 落实 |
| `mutation-model-design.md` §6 | RangeFor 通过 loop-as-function 改造 | 部分：Pass 6 真展开做了**索引循环**但**未处理 `auto&` 写回** |

### P0-1 根因

`mutation-model-design.md` §6（迭代器模式）讨论的是 begin/end 风格迭代器（`it != end; ++it`），未明确指出 ranged-for 中**循环变量本身是 reference 时的写回语义**。Pass 6 RangeFor 真展开（`pass6_ssa_build.py:182-244`）按"index loop"模型展开：
```
body: term := cont[idx]; <user body>
latch: idx := idx + 1
```
默认了 `term` 是 fresh 的局部副本——这等价于 C++ `for (auto term : c)`（值拷贝），**不**等价于 `for (auto& term : c)`（引用绑定）。

判定 ref vs 非 ref 的信息在 Pass 1：`RangeForStmt` AST 解析时 Clang AST 给出 loop var 的 `qualType`，含 `&` 标识。但 `RangeForStmt` IR 字段（`ir_types.py:303-309`）没有 `is_mutable_ref` 标志：
```python
@dataclass
class RangeForStmt:
    var: Var
    var_ty: TypeIR
    container: ExprIR
    body: list['StmtIR']
    decomposition: Optional[list[Var]] = None
```

Pass 1 解析时丢弃了 ref 信息，Pass 6 desugar 时不知是值还是引用，默认按值。

### P0-2 根因

Pass 2 ref_elim_pass（`pass2_ref_elim.py:38-119`）只对**当前正在翻译的函数**做了签名改造（pair-return + 末尾 ReturnStmt 包装）。但当前函数 body 里调用其他 ref-out 函数时，**调用点的 HIR 没有重写**：
- HIR 里 `Call(callee="__edf_Zp", args=[...])` 是表达式
- 它包在 `ExprStmt(Call(...))` 中
- Pass 6 把 ExprStmt 转成 `__sideeff_X := <call>`（discarded return）

`TRANSLATION_SCOPE_OUTPUT_PARAMS` 注册了哪些函数有 ref-out（19 个），但**没有任何 Pass 用这个表去改写调用点**。仅 `smoke_pass1_full.py` 用它做 ref-out 检测的 sanity check。

### P1 根因

P1 与 P0-2 同根：`__upoly_make_monic` 在 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 注册为 mutates [0]，但调用点没改写。

修 P0-2 自动修 P1。

---

## 第三部分：修复方案

### 总体策略

**两个独立修复**（无依赖关系）：

1. **P0-1**：在 Pass 1 给 `RangeForStmt` 增加 `is_mutable_ref: bool` 字段，Pass 6 RangeFor desugar 在 `is_mutable_ref=True` 时插入写回。
2. **P0-2 + P1**：在 Pass 2 ref_elim 之后添加调用点重写，或在 Pass 5 现有的 stmt-level ExprStmt 处理中按 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 改写。

两者都不需要改 `mutation-model-design.md` 设计——只是把已设计未落实的部分实施掉。

### 修复 P0-1：RangeFor `auto&` 写回

#### 步骤 A1（Pass 1 + ir_types）：给 RangeForStmt 加 `is_mutable_ref` 标志

`ir_types.py` 增加字段：
```python
@dataclass
class RangeForStmt:
    var: Var
    var_ty: TypeIR
    container: ExprIR
    body: list['StmtIR']
    decomposition: Optional[list[Var]] = None
    is_mutable_ref: bool = False  # B8 后续：auto& 写回所需
```

Pass 1 `parse_pass` 在解析 `CXXForRangeStmt`（或 `ForRangeStmt`）时检测 loop var decl 的 qualType：
```python
# 伪代码
loop_var_decl = forrange_node["loopvar"]
qual_type = loop_var_decl["type"]["qualType"]  # 例如 "int &" / "const int &" / "int"
is_mutable_ref = ("&" in qual_type) and ("const" not in qual_type)
```

#### 步骤 A2（Pass 6 RangeFor desugar）：写回回写

`pass6_ssa_build.py` `_process_one` 的 `RangeForStmt` 分支，在 latch 之前插入条件回写：

```python
# 当前实现：
latch.stmts.append(AssignStmt(
    target=idx_var,
    value=BinOp(op="+", lhs=idx_var, rhs=Lit(1, ty=BaseType.INT64),
                ty=BaseType.INT64),
))
```

改为：
```python
# 新增：mutable ref 时回写到 cont
if s.is_mutable_ref:
    latch.stmts.append(AssignStmt(
        target=cont_var,
        value=Call(callee="__write__",
                   args=[ArrayAccess(arr=cont_var, idx=idx_var, ty=s.var_ty),
                         s.var],
                   ty=UnknownType("")),
    ))
# 然后 idx ++
latch.stmts.append(AssignStmt(
    target=idx_var,
    value=BinOp(op="+", lhs=idx_var, rhs=Lit(1, ty=BaseType.INT64),
                ty=BaseType.INT64),
))
```

SSA rename 阶段会自然处理：
- `cont_var` 是 def_names 成员（`_collect_def_blocks_and_types` 已包含 root-update 检测，B7 同款），所以 cont 版本 bump
- `s.var`（`term`）的 latest version 通过 stack 取顶（rename_reads）

#### 步骤 A3（验证）

新增 unit test `test_rangefor_auto_ref_writeback`：
```python
def test_rangefor_auto_ref_writeback():
    """for (auto& x : c) x = x + 1  →  cont 必须被 bump"""
    body = [
        RangeForStmt(
            var=Var("x", ty=BaseType.INT64),
            var_ty=BaseType.INT64,
            container=Var("c", ty=ArrayType(BaseType.INT64)),
            body=[AssignStmt(target=Var("x"),
                             value=BinOp("+", Var("x"), Lit(1)))],
            is_mutable_ref=True,
        ),
        ReturnStmt(value=Var("c")),
    ]
    func = _mk_func(body, params=[HIRParam("c", ArrayType(BaseType.INT64))],
                    ret_ty=ArrayType(BaseType.INT64))
    mir = ssa_build_pass(func)
    # 验证 cont 在 latch 被 bump
    has_writeback = False
    for bb in mir.cfg.blocks.values():
        for s in bb.stmts:
            if isinstance(s, LetStmt) and s.var.name.startswith("__rangefor_cont_"):
                if isinstance(s.value, Call) and s.value.callee == "__write__":
                    has_writeback = True
    assert has_writeback
```

烟测 `smoke_pass6_full.py` 增加 metric：`rangefor_no_writeback` = ranged-for body 含 `__write__(<var>.field, ...)` 但 latch 没有对应 cont 写回的次数。期望降到 0。

### 修复 P0-2 + P1：调用点 ref-elim

#### 选项 1：新建 Pass 2.5 callsite_ref_elim

在 Pass 2 与 Pass 3 之间插入新 Pass，专做调用点改写。

```python
# pass2b_callsite_ref_elim.py（新文件）
def callsite_ref_elim_pass(func: HIRFunc) -> HIRFunc:
    """HIR₁ → HIR₁'：把 `f(out, in)` 改写为
       `let (out_new, ret_new) := f(out, in); out = out_new; return_handler(ret_new)`。
    """
    new_body = _walk_stmts(func.body)
    return replace(func, body=new_body)

def _rewrite_call_site(s: ExprStmt) -> list[StmtIR]:
    if not isinstance(s.expr, Call): return [s]
    call = s.expr
    callee = call.callee if isinstance(call.callee, str) else None
    if callee not in TRANSLATION_SCOPE_OUTPUT_PARAMS:
        return [s]
    out_indices = TRANSLATION_SCOPE_OUTPUT_PARAMS[callee]
    out_args = [call.args[i] for i in out_indices]
    # 生成 destructure
    if len(out_indices) == 1 and len(out_args) == 1:
        # `out := f(out, in)` 形式
        return [AssignStmt(target=out_args[0], value=call)]
    else:
        # `(out1, out2) := f(out1, in, out2)` 形式
        # 用 TupleExpr LHS 还是分别赋值？
        # 方案：临时变量 + field access
        tmp_var = Var(f"__refret_{counter}")
        result = [LetStmt(var=tmp_var, ty=UnknownType(""), value=call)]
        for k, idx in enumerate(out_indices):
            field = "fst" if k == 0 else f"snd" * (k)  # 简化
            result.append(AssignStmt(
                target=out_args[k],
                value=FieldAccess(obj=tmp_var, field_name=field, ty=...),
            ))
        return result
```

**优点**：职责单一，与 Pass 2 解耦。
**缺点**：增加一个 Pass，需要更新 smoke runner。

#### 选项 2：扩展 Pass 5

把调用点改写合并进 Pass 5（operator_resolve）的 ExprStmt 处理（`pass5_operator_resolve.py:976+`）：

```python
if isinstance(s, ExprStmt) and isinstance(s.expr, Call):
    callee = s.expr.callee
    if isinstance(callee, str) and callee in TRANSLATION_SCOPE_OUTPUT_PARAMS:
        # 改写为 destructure assign
        return _rewrite_ref_call_site(s, gap)
```

**优点**：不增加 Pass。
**缺点**：Pass 5 已经很复杂；ref-elim 是 Pass 2 的职责。

#### 推荐：选项 1（新 Pass 2.5）

理由：
1. 职责清晰，对应 mutation-model-design.md §3.2；
2. 修复 + 测试 ~150 行，独立成 Pass 利于回归；
3. 未来若有 expression-level call（如 `x = f(out, in)` 中 `x` 也算 return）扩展更平滑。

#### 步骤 B1：创建 `passes/pass2b_callsite_ref_elim.py`

核心实现见上方伪代码。需要处理的细节：

1. **多输出参数解构**：使用 `TupleExpr` LHS 是否在 IR 已支持？看 `ir_types.py:266 AssignStmt.target` 是 `ExprIR`，技术上 TupleExpr 可作为 lhs，但语义未明。**简化方案**：临时变量 + FieldAccess 取分量（"fst"/"snd"）+ 多个 AssignStmt：
   ```
   let __refret_N := f(...)
   out0 = __refret_N.fst
   out1 = __refret_N.snd
   ```
   注意：当 `f` 原本非 void 时，第一个分量是 return value（需丢弃或作为 expression-level result）。
2. **嵌套表达式**：若调用作为 sub-expression（如 `bool ok = f(out, in)`），更复杂。当前 67 函数 corpus 中**未观测**——需 grep 确认。
3. **transparent passthrough**：调用 `f(g(out), in)`——`g` 也可能是 ref-out 函数？grep 确认未观测。
4. **递归调用**：`__edf_Zp` 内部调用 `__edf_Zp` —— ref-out 列表是已知的。OK。

#### 步骤 B2：注册 Pass 2.5 到 smoke 运行器

更新：
- `passes/__init__.py`（如有）export
- `tests/run_all_smoke.sh` 跑顺序：Pass 1 → 2 → **2.5** → 3 → 4 → 5 → 6
- 复用现有 smoke 脚本：在 `smoke_pass2_full.py` 增加 callsite-ref-elim metric；或新建 `smoke_pass2b_full.py`

#### 步骤 B3（验证）

新增 unit test `test_pass2b_callsite_ref_elim.py`：
- T1：`f(out, in)` 调用单 out → `out := f(out, in)`
- T2：`f(out1, in, out2)` 双 out → 临时 + 分量赋值
- T3：`r = f(out, in)` 非 void + out → 临时 + return + out 赋值
- T4：递归调用同函数（`__edf_Zp` self call）
- T5：调用点不在 TRANSLATION_SCOPE_OUTPUT_PARAMS（passthrough）

烟测 `smoke_pass6_full.py` 现有 metric `sideeff_write` 应降到 0；新增 metric `unwritten_ref_call` = 调用 ref-out 函数但目标参数未 bump 的次数，期望 0。

---

## 第四部分：同类问题梳理（survey）

按"call site 不正确反映 mutation"这一根因，系统性扫描可能受影响的模式。

### Grep 实测结果（2026-04-27, 67 函数 corpus）

| 模式 | 实测命中 | 涉及函数 |
|------|----------|---------|
| **P0-1** ranged-for `auto&` 缺写回 | **40** | `__wang_leading_coeff`、`factorize_grlex`、`__hensel_step`、`__hensel_lift_recursive`、`__factor_multivar`、`__build_cld_matrix` 等 8+ 函数 |
| **P0-2** ref-out 调用点 sideeff | **47** | 跨 11 个被注册函数（`__upoly_make_monic` 9 次，`__upoly_mod_coeff` 15 次，`__edf_Zp` 4 次，等） |
| 模式 3 `*it = x` iterator deref assign | **19** | `__write__(__deref__(it).field, ...)` 形态（StdMap 迭代器 `it.second = ...`） |
| 模式 9 `std::swap(a, b)` | **2** | `__lll_reduce`（mu_9）、`__build_cld_matrix`（M_6）的行交换 |
| 模式 1 未注册 method | 待审计 | `Array.*` 多数已正确（值返回路径） |

**总计 silent semantic miss ≥ 108 处**（仅可机械识别的部分），跨 67 函数。

### 同类问题：未注册到 TRANSLATION_SCOPE_OUTPUT_PARAMS 的 mutating 调用

按 sideeff callee 出现频率 top 15（grep 实测），与 C++ 源签名核对：

| Callee | sideeff 次数 | C++ 签名（确认） | 应注册索引 | 当前状态 |
|--------|---|---|---|---|
| `poly_convert` | 26 | `void poly_convert(const X& p_in, Y& p_out)` (`upolynomial.cc:9`) | `[1]` | **未注册** |
| `__upoly_mod_coeff` | 15 | ref-out [0] | `[0]` | ✓ 已注册 |
| `fdiv_r` | 14 | `void fdiv_r(ZZ& r, const ZZ&, const ZZ&)` (`ZZ.hh:810`) | `[0]` | **未注册** |
| `pair_vec_div` | 12 | `void pair_vec_div(q&, r&, ...)` (`polynomial_factorize_wang.hh:48`) | `[0, 1]` | **未注册** |
| `__upoly_make_monic` | 9 | ref-out [0] | `[0]` | ✓ 已注册 |
| `sort` | 8 | `std::sort(begin, end)` mutates range | special（迭代器范围） | **未注册** |
| `iota` | 5 | `std::iota(begin, end, val)` mutates range | special | **未注册** |
| `fdiv_q` | 5 | `void fdiv_q(ZZ& q, ...)` (`ZZ.hh:731`) | `[0]` | **未注册** |
| `__upoly_divmod_mod` | 4 | ref-out [0,1] | `[0, 1]` | ✓ 已注册 |
| `__edf_Zp` | 4 | ref-out [0, 3] | `[0, 3]` | ✓ 已注册 |

**结论**：除已注册 19 个外，至少 **5 个高频 mutating 函数未注册**：`poly_convert`、`fdiv_r`、`fdiv_q`、`pair_vec_div`、`sort`/`iota`（迭代器范围 mutating）。

修复 P0-2 时应一并补齐这些注册。修后预计 `sideeff_*` callee 中的 mutating 调用全部消失，仅留下真正 discarded return 的副作用调用（如 `__write__` 的 root-is-param 情况、debug print 等）。

### 模式 6 lambda capture by ref（已被 Pass 3 lift 消除）

实测：`/tmp/hir3_dump/*.txt` 中 `Capture` 0 命中。Pass 3 lambda_lift 把 capture 转为参数，by-ref capture 应映射为参数 `is_ref=True`。需审计 Pass 3 lift 是否正确传递 ref 标志：lifted lambda 的参数列表中 by-ref capture 是否标 `is_ref=True`，从而被 Pass 2 ref-elim 走 pair-return 路径。

**审计命令**：
```bash
grep -h "_lambda_.*aux" /tmp/hir3_dump/*.txt   # 找 lifted lambda
# 然后检查这些 lambda 的 params 是否有 is_ref=True
```

属于"待审计"项；不在本方案修复范围。

### 模式 8 std::move 后再用

实测：`grep "move(" /tmp/mir0_dump/*.txt` 共 56 次；多数模式为 `Array.push(M_1, move(new_row_4))`、`f_star_3 := move(f_new_1)`——move 作为最后一次使用，安全。

**审计命令**：列出每个 move 的 operand，看 operand 后续是否再被读：
```python
# 在 dump 中找每个 move(X) 的 X，检查 X 后续在 stmt 中是否再被引用
```

属于"待审计"项；不在本方案修复范围。



### 模式 1：成员函数 mutating method 调用（已部分处理）

```cpp
vec.push_back(x);
m.insert({k, v});
*it = x;
it.advance();
```

**当前实现**（`pass5_operator_resolve.py`）：CLASS_MAP 标注 mutate disp 的 method 已转 `vec = vec.push x`；但仅限 CLASS_MAP 中显式注册的方法。

**风险**：未注册到 CLASS_MAP 的 mutating method（如自定义类的 `clear` / `resize` / `assign`）会被当成 query method，silent 丢 mutation。

**修复**：扩展 CLASS_MAP 覆盖率审计 + Pass 5 在遇到未识别 method 时打 gap warning（gap_log）。

### 模式 2：自由函数 ref-out（P0-2 涵盖）

已在 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 注册 19 个，本修复方案落实。

**风险**：未注册的自由函数（含 STL：`std::sort`、`std::swap`）。`std::sort(begin, end)` 通过迭代器修改 vector——当前不在注册表，但 Pass 5 可能已通过 CLASS_MAP `sort` 处理。**需 grep 确认**。

**修复**：审计所有出现在 dump 中的自由函数调用，对照 STL 文档 + CLPoly 源标注 mutating params；缺失加入 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 或 CLASS_MAP。

### 模式 3：ranged-for `auto&`（P0-1 涵盖）

本修复方案落实 `is_mutable_ref` 字段。

**风险**：iterator-style for（`for (auto it = c.begin(); it != c.end(); ++it) *it = x`）。`*it = x` 在 Pass 5 的 mutator 路径可能没覆盖。

**修复**：grep `__deref__(it).` `= ...` 模式确认；加入 Pass 5 dereference assign 的回写处理。

### 模式 4：标量 ref out（少见但存在）

```cpp
void compute(int& out) { out = 42; }
```

`TRANSLATION_SCOPE_OUTPUT_PARAMS` 已注册若干（`__upoly_divmod_mod` `[0,1]` 是 poly out）。标量 out 也走同路径。**P0-2 修复涵盖**。

### 模式 5：嵌套 record-update + ranged-for

```cpp
for (auto& row : matrix) {
    for (auto& cell : row) {
        cell *= 2;
    }
}
```

**风险**：双层 ranged-for `auto&`。修复 P0-1 后，外层 `row` 在 inner loop 结束时是否正确"回写到 matrix[i]"？需要 inner loop latch 把 `row` 写回 matrix，外层 latch 再触发……

**分析**：B7 record-update 已链式工作（root 是 row → row 自己 bump），但 row 是从 `matrix[i]` 取出的 fresh local。修 P0-1 后，外层 latch 会写 `matrix[i] := row_latest`，所以 row 的所有内层修改链正确传播。**这正是 SSA 的优势**：每层 RangeFor 独立处理，组合自然正确。

**验证**：写嵌套 ranged-for unit test 确认。

### 模式 6：lambda capture by reference

```cpp
int total = 0;
auto add = [&total](int x) { total += x; };
std::for_each(v.begin(), v.end(), add);
// total 被累加
```

**当前实现**：Pass 3 lambda_lift 把 lambda 提到顶层，capture 转为参数。**capture by ref** 应映射为"lambda 也是 ref-out"——但 lifted lambda 的签名当前没有 ref-out 标记。

**风险**：`__factor_multivar` / `__factor_Zp` 用了不少 lambda；如果 capture by ref 没正确传播，结果错。

**审计需求**：grep `Capture(by_ref=True)` 在 lifted lambda 的 corpus；确认 lifted lambda 是否在 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 自动注册。

### 模式 7：对象方法调用通过引用 self mutation

```cpp
poly.normalize();   // 修改 poly 自己
```

**当前实现**：CLASS_MAP 中 `normalize: ("mutate", "Poly.normalize")` → Pass 5 转成 `poly = Poly.normalize poly`。已正确。

**风险**：method 实现不在 TRANSLATION_SCOPE 但 mutate self（外部 STL/Boost 类）。CLASS_MAP 应覆盖。

### 模式 8：std::move 后再使用（UAF）

```cpp
auto x = std::move(y);
// y 已 moved-from
```

**当前实现**：Pass 1 显式 `move()` 调用保留为 Call。Lean 中 functional——move 仍合法（视为 alias）。但若 C++ 在 `std::move(y)` 后又用 y 当输入，Lean 会读到旧版（Lean 不区分 moved），与 C++ 可能差异。

**风险**：若任何函数依赖 moved-from 状态，翻译错误。

**审计需求**：grep `std::move` 在 corpus 中的位置 + 后续是否仍用 moved 变量。

### 模式 9：std::swap

```cpp
std::swap(a, b);
```

**当前实现**：未确认。CLASS_MAP 是否注册？

**风险**：常见函数，silent miss。

**审计需求**：grep `std::swap` 命中数 + Pass 5 处理。

### 同类问题汇总

| 模式 | 实测命中 | 优先级 | 修复路径 |
|------|----------|--------|---------|
| P0-1 ranged-for `auto&` | **40** | P0 | 本方案 §3 阶段 1 |
| P0-2 注册的 ref-out 调用点 | **47** | P0 | 本方案 §3 阶段 2 |
| P1 mutating call SSA bump | 与 P0-2 重叠 | P1 | 本方案附带修复 |
| **新增** 未注册 mutating callee | **57** (poly_convert+fdiv_*+pair_vec_div+sort+iota) | P0 | **本方案 §3 阶段 2 同时补 OUTPUT_PARAMS** |
| 模式 3 `*it = x` deref assign | **19** | P1 | 独立 ticket（与 P0-2 类似但走 iterator） |
| 模式 5 嵌套 ranged-for | 修 P0-1 后单测覆盖 | P2 | 单测 + Agent 审视 |
| 模式 6 lambda capture by ref | 0 (Pass 3 已 lift) | P1 | 审计 Pass 3 ref 标志传播 |
| 模式 7 method self mutate | partial | P2 | CLASS_MAP 持续审计 |
| 模式 8 `std::move` 后再用 | 56 出现/未确认误用 | P2 | grep + 后续读检查 |
| 模式 9 `std::swap` | **2** | P1 | 加入 OUTPUT_PARAMS（[0,1]） |

**建议**：本方案修复
1. P0-1（40 处） + P0-2（47 处） + 未注册 5 函数（57 处） + std::swap（2 处）= **146 处 silent miss**；
2. 独立 ticket 处理模式 3（19 处 iterator deref assign）+ 模式 6（lambda by-ref capture 审计）。

---

## 第五部分：实施计划

### 阶段 1：P0-1 RangeFor `auto&` 写回（~150 行）

| 步骤 | 内容 | 行数估计 |
|------|------|---------|
| A1.1 | `ir_types.py` RangeForStmt 加 `is_mutable_ref` | +1 字段 |
| A1.2 | `pass1_parse.py` 解析 loop var qualType 设置 is_mutable_ref | +5 |
| A2 | `pass6_ssa_build.py` RangeFor desugar latch 加写回 | +10 |
| A3.1 | unit test `test_rangefor_auto_ref_writeback` | +50 |
| A3.2 | smoke metric `rangefor_no_writeback` | +20 |
| A3.3 | 跑全 67 函数 smoke + 验证 dumps | – |
| 文档 | 更新 `mutation-model-design.md` §6 加 ranged-for ref 写回规则 | +30 |

### 阶段 2：P0-2 + P1 调用点 ref-elim（~250 行）

| 步骤 | 内容 | 行数估计 |
|------|------|---------|
| B1 | 新建 `passes/pass2b_callsite_ref_elim.py` | +120 |
| B2.1 | 注册到 `tests/run_all_smoke.sh` | +5 |
| B2.2 | smoke runner `smoke_pass2b_full.py` | +50 |
| B3 | unit test `test_pass2b_callsite_ref_elim.py`（5 用例） | +120 |
| 文档 | `mutation-model-design.md` §3.2 加"实现 location: pass2b_callsite_ref_elim" | +5 |
| 文档 | `tech-debt.md` 关闭 TD-8 | -10 |

### 阶段 3：回归 + 同类问题审计（~3 小时）

| 步骤 | 内容 |
|------|------|
| C1 | 跑 `bash tests/run_all_smoke.sh` 全过 |
| C2 | 重新派 3 agent 做 Pass 6 1-对-1 审视（含语义模拟） |
| C3 | 验证：`__edf_Zp` 递归调用、`__factor_multivar` 系数取负、`__build_cld_matrix` 行追加在 dump 中正确 |
| C4 | 第四部分模式 1/3/6/9 grep 审计 → 各起独立 ticket（不在本方案范围） |

### 验证标准（修复完成判据）

1. `bash tests/run_all_smoke.sh` 全绿
2. `smoke_pass6_full.py` 输出：
   - `sideeff_write` ≤ 2（仅外部 deref/iterator，不含已注册 ref-out 调用）
   - `phi_undef_ver0` = 0
   - 新增 `rangefor_no_writeback` = 0
   - 新增 `unwritten_ref_call` = 0
3. Agent B 第三轮语义模拟：`__factor_multivar`、`__edf_Zp`、`__build_cld_matrix` semantic 评级"OK"
4. unit test 总数：mir_types 9 + pass6 9 + **pass2b 新 5** = 23 通过

### 风险与回退

**风险 1**：Pass 2.5 改变 HIR₁ 结构，下游 Pass 3/4/5 可能踩坑。
- **缓解**：先跑 smoke pass1/pass2/pass2b 单独验证；通过再跑 Pass 3+。

**风险 2**：嵌套 ranged-for 在 inner break 时外层回写 corner case。
- **缓解**：T6 嵌套 unit test + Agent B 抽样审视嵌套用例。

**风险 3**：TupleExpr 作 LHS 在 IR 未支持，多输出 destructure 用临时变量+FieldAccess 字段名（fst/snd）可能与 PairType 实际字段不一致。
- **缓解**：先核查 `ir_types.py` PairType / TupleType 字段名；与 Pass 8 codegen 约定一致。

**回退**：本方案两个修复彼此独立。任一阶段失败 → revert 该阶段；另一阶段已 commit 不受影响。

---

## 附录 A：grep 命令清单（供审计）

```bash
# P0-1 残留（修复后应 = 0）
grep -h "rangefor.*term.*__write__" /tmp/mir0_dump/*.txt
grep -B2 "let term_._ := __write__" /tmp/mir0_dump/*.txt | grep -B2 -A1 "goto bb"

# P0-2 残留（修复后应 = 0）
for f in __upoly_make_monic __edf_Zp __upoly_mod_coeff __upoly_divmod_mod \
         __hensel_step __hensel_lift_recursive __build_cld_matrix \
         __lll_reduce __si_vandermonde_solve __si_theta_array_eval \
         __mtshl_zp_univar_mdp __mtshl_multi_bdp __mtshl_sparse_int \
         __mtshl_wmds __mtshl_step_j; do
    echo "=== $f ==="
    grep -h "__sideeff_.*:= $f(" /tmp/mir0_dump/*.txt | head -3
done

# 模式 1 自定义类未注册 mutating method
grep -h "__sideeff_.*\.\w*(" /tmp/mir0_dump/*.txt | sort -u | head -30

# 模式 3 iterator dereference assign
grep -h "__deref__.*=" /tmp/mir0_dump/*.txt | head -10

# 模式 6 lambda capture by ref
grep -n "Capture.*by_ref=True\|Capture(.*ref=True" /tmp/hir3_dump/*.txt | head

# 模式 8 std::move 后用 moved
# (需读 source 不是 dump)

# 模式 9 std::swap
grep -h "swap(" /tmp/mir0_dump/*.txt | head
```

## 附录 B：相关文档

- `docs/design/l1-translation-validation/mutation-model-design.md` §3 G1 输出参数（设计已存）
- `docs/design/l1-translation-validation/tech-debt.md` TD-8（已记录但未关）
- `proof/cpp2lean_v2/class_map.py:684 TRANSLATION_SCOPE_OUTPUT_PARAMS`（已注册 19 个）
- `proof/cpp2lean_v2/passes/pass2_ref_elim.py`（仅做被调函数签名）
- `proof/cpp2lean_v2/passes/pass6_ssa_build.py:182-244` RangeForStmt desugar（缺 ref 写回）
