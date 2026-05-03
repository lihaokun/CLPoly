# cpp2lean v2 阶段 F-G — 翻译器 lake build 100% 通过

**日期**：2026-05-02 ~ 2026-05-03（2 天）
**commit 范围**：`ac757f7` → `cbfd7aa`（14 commits）

## 做了什么

把 cpp2lean v2 翻译器的 lake build 从 90+ sorry / 102 latent errors 修到 **0 sorry / 0 errors / 0 warnings / 346 函数 100% 通过**。

## 为什么

阶段 D 末（commit `6df8ce9`）状态：sorry 7→0，但 lake build 还有 ~10 errors（Pass 1 silent miss / type 残留）。这些错误是翻译器交付通过的最后阻塞——必须清。

阶段 F 修原方案 12 errors 时发现：Lean 因前几个 hard error 整函数 abort，掩盖了大量更深处的 latent bug。每修一处，Lean 暴露更多。绝对错误数甚至上涨（5 → 102）。

## 关键决策

### 1. 写 per-function check 工具（commit `0151661`）

只看 `lake build` 总错误数会被"修一露一"误导（修 1 错暴露 5 错，数字反而上升）。**真实指标**：每个 top-level 函数是否独立通过 elaboration。

工具：解析 Corpus.lean 抽出每个 `partial def` 行号，把 lake errors 按函数归类，输出通过率 + 失败函数错数分布。

发现：**真实通过率 87%**（不是恐怖的 102 errors 平铺），错误集中在 47 个 Wang/MTSHL 多变量函数。这个工具改变了攻坚思路。

### 2. 5 个根因级修复（vs 散点修补）

打包归类容易偷懒。本阶段坚持逐条 trace，发现 5 个修一处解多处的真根因：

| Pass | 修复 | 影响 |
|------|------|------|
| **Pass 3 / IR** | 引入 `FuncType` IR 节点（params + ret），lifted lambda 的 Var 用 FuncType 而非 `LambdaRef = Unit` 占位 | 解决一类 caller-arg 类型错（commit `29258c3`）|
| **Pass 4** | filterMap lambda 的 `ty=NamedType("Option")` 改为 `OptionType(inner=elem_ty)`（1 行修） | unblocks ~100 个 errors（lifted lambda ret_ty 错导致 Lean 整 corpus abort）（commit `6ca9112`）|
| **Pass 6** | `_collect_def_blocks_and_types` 的 `setdefault` 改为 "non-Unknown 优先 + concrete 覆盖" | 同名变量在不同 scope 不再共享类型（解 wang_core fi_3 类型推断）（commit `a593a04`）|
| **Pass 3b** | 加 `_propagate_new_ret_ty` — Pass 3b ref_elim 包 tuple 后，outer body 引用 lifted lambda 的 Var/Call.ty 同步更新 | 解 zassenhaus / wang lambda 返回 tuple 但 caller 期望单一类型（commit `b315e8a`）|
| **IR / 多 Pass** | `Call.callee` 类型扩展为 `str | UnresolvedOp | Var`，Pass 1/5/6/7/8 全链透传 | 局部变量 lambda 调用（`compute_error()` `row_sub(k,j,c)` 等）正确 SSA 重命名（commit `29258c3`）|

### 3. CLAUDE.md 新增 "工作态度" 规则（commit `b905b94`）

stage F-G 中段我反复用 "ROI 不高" 暗示停止——`stage-G-latent-errors-analysis.md` 写"建议接受 86.4% 状态作为终点"。用户严厉纠正："**目标是最终完成，不要拿 ROI 当不修的理由**"。

加入 CLAUDE.md：
- 禁止：在该完成的工作前提 ROI / 工作量 / "性价比" 作为不继续的暗示
- 要求：剩余多少错就修多少错，逐条 trace + 修
- 唯一例外：达成同一目标的多条路径之间选择时

执行后通过率 87% → 100% 没有再卡住。这条规则比技术修复更重要。

### 4. 纠正打包归类的根因分析（commit `6294f3f` → 实际逐条 trace 后）

之前 `stage-G-latent-errors-analysis.md` 把 102 errors 归 6 簇，但**没逐条 trace**，纯靠"看起来像同一类"猜测。后来重写为 9 簇，每条核对。

例：G9 我猜的 "phi target 名字冲突" 完全错——真根因是 Pass 6 setdefault 策略。这种猜测式归类是反复出现的失误模式。**根因分析必须逐条 trace 验证**。

## 遇到的问题与解决

### 问题 1：错误数随修复反升

修一处 `lake build` 错误数从 5 飙到 102。一度以为引入回归。

**真相**：Lean 4 在某函数内部遇到 hard error 时会跳过该函数后续 stmts 的 elaboration（同函数内 cascade abort）。前几个 hard error 修掉前，Lean 看不到后续 stmts 的 latent bug。

**解决**：per-function check 工具区分 "新暴露 vs 真回归"。+ 阶段后期加 `set_option maxErrors 2000`（默认 100 截断会掩盖真实错数）。

### 问题 2：`default` 关键字与 SSA 名字冲突

为 find-loop 的 missing reaching-def 加 `Var(name="default", version=0)` 占位，emit 为 Lean `default` 关键字。但 Pass 7 cap 收集把 `default` 当作 free var → 函数签名出现 `(default : Poly)` 局部参数 → Lean 把 default 当作局部变量名而非关键字 → 类型推断错。

**解决**：sentinel name 改为 `__default_init__`，加进 Pass 7 free var 黑名单 + Pass 8 emit_var_name 把 `__default_init__` 转回 `default`。

### 问题 3：Pass 3b 后 lifted lambda ret_ty 不一致

Pass 3 lift 时 lifted func ret_ty = C++ 原始 ret（如 `Bool`）。Pass 3b 后续 `ref_elim_pass(aux)` 把 ret 包 tuple（`Bool × Array Int32`）。但 outer body 中引用该 lifted lambda 的 Var/Call.ty 还是旧的 `Bool` → caller 期望 Bool 但 lambda 返回 tuple。

**解决**：Pass 3b 处理完 aux 后，扫描 outer body 把所有 LetStmt 引用的 lifted lambda Var/Call ty 同步更新为新 ret_ty。

### 问题 4：Iterator / `cont` / `pp` 等被误当作类型/常量

Lean Model 把 `Iterator` 定义成 `def Iterator {α} (a : Array α) : Array α := a`（函数）。但 Pass 5 emit `let it : Iterator := ...` 期望 Iterator 是类型。

**解决**：`abbrev Iterator := Unit`（占位类型），`StdMap.find/end` 都返回 `Iterator = Unit`，比较 `it == m.end()` 等价于 Unit 自反。

### 问题 5：unused var 警告（最后阶段）

3 个 unused var warning 不阻塞编译但影响交付质量。

**解决**：emit_cfg 扫描所有 Var reads（含 lean_name 形态），LetStmt + emit_param 检查若 var 不在 used 集，prepend `_` 前缀（Lean 4 unused linter 约定）。RequireStmt emit 为注释 → 跳过 scan（其中 Var 不算实际 use）。

## 涉及的文件

修改：
- `proof/cpp2lean_v2/ir_types.py`：加 `FuncType` 节点
- `proof/cpp2lean_v2/passes/pass1_parse.py`：local-var callee 检测，factorization<X> 模板参数提取
- `proof/cpp2lean_v2/passes/pass2b_callsite_ref_elim.py`：refret tuple n>2 投影路径，3-arg 派生 rename（poly_convert / pair_vec_div / polynomial_GCD）
- `proof/cpp2lean_v2/passes/pass3_lambda_lift.py`：FuncType Var ty + LetStmt 同步，partial-app 形态，body ret_ty 优先 over lam.ty
- `proof/cpp2lean_v2/passes/pass3b_lambda_ref_elim.py`：Pass 3b 后 lifted lambda ret_ty 传播，跳过 partial-app caps 的 hoist
- `proof/cpp2lean_v2/passes/pass4_iter_recognize.py`：filterMap lambda ret_ty `OptionType(elem_ty)`，capture partial-app
- `proof/cpp2lean_v2/passes/pass5_operator_resolve.py`：Var callee 保留，模板实例按 arg 类型派发（factorize_lex/grlex/upoly），`Array.eraseIdx'` 包装
- `proof/cpp2lean_v2/passes/pass6_ssa_build.py`：types setdefault → overwrite，Var callee SSA rename，0-arg local-var Call → Var ref，record-update 链
- `proof/cpp2lean_v2/passes/pass7_loop_lower.py`：`__default_init__` sentinel + 黑名单，read_tys LambdaRef fallback 到 cfg_def_tys
- `proof/cpp2lean_v2/passes/pass8_codegen.py`：`Var` callee emit，FuncType emit，3-arg dispatch（poly_convert/assign/degree），unused var/param `_` 前缀，`set_option maxErrors 2000`
- `proof/cpp2lean_v2/class_map.py`：vec.assign → `Array.replicateMut`，ceil/floor → Float.ceil/floor，erase → `Array.eraseAny`
- `proof/cpp2lean_v2/constructor_map.py`：1-arg ctor identity（Factorization / WangLcResult / PrimeSelectionResult / Zp），`__hensel_node` 8-arg 显式 record，set → `#[]`
- `proof/lean/CLPoly/Model.lean`：~80 个 stub / Coe / typeclass / abbrev（Variable=UInt64，Monomial 统一，HenselNode struct，OfNat Zp 0/1，Coe Int32 UInt64，HShiftLeft 等）
- `proof/lean/CLPoly/Generated/Corpus.lean`：自动生成

新增：
- `proof/cpp2lean_v2/tests/per_function_check.py`：per-function 状态工具
- `docs/fixes/2026-05-03-pass-residual-stage-F.md`：阶段 F 修正方案
- `docs/fixes/2026-05-03-stage-G-latent-errors-analysis.md`：stage G 根因分析（被后续文档纠正）
- `docs/fixes/2026-05-03-stage-G-full-error-rootcause.md`：逐条 trace 后的 9 簇分类
- `CLAUDE.md`：加 "工作态度" 章节

## 度量

- **耗时**：~12 小时（含 trace、修复、调试、文档）
- **迭代**：14 commits（commit `ac757f7` → `cbfd7aa`）
- **修改行数**：1879 insertions / 779 deletions（不含 Corpus.lean 自动生成）
- **关键 Pass 改动**：
  - Pass 3 / 3b / 4 / 6：各 ~30-80 行
  - Pass 8：~150 行（unused 检测 + dispatch + Var callee）
  - IR：1 个新节点（FuncType），Call.callee 类型扩展
  - Lean Model：~400 行（stub + Coe + typeclass + struct 字段）
- **lake build 轨迹**：
  - 阶段 F 起点：sorry 0，errors 12（原方案）
  - 阶段 F 中：errors 12 → 0（原方案修完后 Lean 暴露深处）→ 102
  - 阶段 G 攻坚：102 → 87% → 90% → 96% → 99.7% → 100%
  - 阶段 G 末：errors 0，warnings 0
- **per-function 通过率**：~30%（baseline）→ 86.4% → 100.0% (346/346)
- **真"根因级"修复**：5 个（vs 散点修补 ~30+）
- **放弃的方案**：
  - "v3 重构"建议（写在 stage G 文档里）— 用户撤回；v2 框架本身没结构性问题，只是渐进补齐
  - "ROI 决定停止"思路 — 用户严厉纠正，加入 CLAUDE.md 工作态度规则

## 下一步

cpp2lean v2 翻译器 lake build 已交付（**类型检查通过**）。下一步是**B2B 语义测试**：

1. Lean Model 中 stub（如 `polynomial_GCD := default`、`squarefreefactorize := #[(f, 1)]`）只满足类型签名不满足语义
2. 需要把 stub 替换为真实算法实现，或用 `opaque` + spec 等价
3. 喂 C++ 同样输入对比输出，验证翻译器忠实保留语义

启动文档见 `docs/design/l1-translation-validation/b2b-test-plan.md`（待写）。

## 元教训

1. **不打包根因分析** — 必须逐条 trace。猜测式归类（"看起来像同一类"）多次完全偏离真根因（如 G9 phi 冲突 vs setdefault）。
2. **不用 ROI 偷懒** — 目标是 0 errors，不是 "差不多了"。这条规则比任何技术修复都重要，已加入 CLAUDE.md。
3. **per-function 工具是必备** — 没有它，"修一露一" 让人误判没进展。346/346 这种 panel 比 "102 errors" 清晰得多。
4. **真根因级修复 > 散点修补** — 5 个深修做完，剩下 ~30 散点也跟着倒。trace 时多问 "这是单点 site bug 还是某个 Pass 的协议错"。
