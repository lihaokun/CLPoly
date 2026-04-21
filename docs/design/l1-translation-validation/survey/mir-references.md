# MIR / Codegen / Back-to-Back 框架调研

> Stage 1 Week 5 准备资料。为 Pass 6（`ssa_build`）、Pass 7（`loop_lower`）、
> Pass 8（`codegen`）和 back-to-back 测试框架的细化设计服务。
>
> 对应 `translator-v2-plan.md` §3.2 Pass 6-8、`back2back-design.md` 整体。

---

## §1 Cytron 1991 SSA 算法要点

参考：Cytron, Ferrante, Rosen, Wegman, Zadeck. *Efficiently Computing Static
Single Assignment Form and the Control Dependence Graph*, TOPLAS 13(4), 1991.

### 1.1 四阶段流水线

```
HIR₄ ──(a)──▶ CFG ──(b)──▶ DomTree ──(c)──▶ DF(n) ──(d)──▶ SSA (phi + rename)
```

1. **(a) CFG 构造**：以基本块为节点，控制转移为边。从 HIR₄ 的结构化控制流
   （if/while/for）下降：while 展开为 header/body/exit 三块 + back-edge；
   if 展开为 cond/then/else/merge 四块。
2. **(b) 支配树 (Dominator Tree)**：节点 `d` 支配 `n` 当且仅当从 entry 到
   `n` 的所有路径都经过 `d`。`idom(n)` 是最深的严格支配者。
   - 简化算法：Cooper-Harvey-Kennedy 2001（数据流迭代法，~20 行），
     对小函数足够快（我们 67 个函数均 < 300 行）。
   - 正规算法：Lengauer-Tarjan O(E·α(E,N))，实现复杂、收益小。
3. **(c) 支配边界 (Dominance Frontier)**：`DF(n) = { m | n 支配 m 的前驱但
   不严格支配 m }`。Cytron 原文 Algorithm 从子节点向上累积：
   ```
   for each node b in post-order:
     DF(b) = {s ∈ succ(b) | idom(s) ≠ b}  -- DF_local
     for each child c of b in DomTree:
       for each w ∈ DF(c):
         if idom(w) ≠ b: DF(b) += {w}     -- DF_up
   ```
4. **(d) Phi 放置 + 变量重命名**：
   - **Phi 放置（迭代 DF⁺）**：对每个变量 `v`，从 `v` 的所有定义点 `S` 出发，
     `phi-set = DF⁺(S)`（DF 的传递闭包）。在这些块的开头插入 `v = phi(...)`。
   - **重命名（DFS on DomTree）**：为每个变量维护一个版本栈。访问块时：
     读取 = 栈顶版本；写入 = 压入新版本；phi 参数 = 前驱块出口的栈顶版本；
     回溯时弹出本块新压入的所有版本。

### 1.2 Minimal vs Pruned SSA

- **Minimal SSA（Cytron 原版）**：在 DF⁺(defs) 放 phi，不管后续是否用到。
- **Pruned SSA**：只在该变量在此块"live-in"时放 phi（要求活跃性分析）。
- **Semi-pruned SSA**：折中，不跨基本块边界的局部变量不放 phi。

**CLPoly v2 选 Minimal SSA**。理由：(1) 实现最简单（无需活跃性分析）；
(2) 67 个函数规模小，冗余 phi 在 `loop_lower` / Lean elaboration 阶段会被
`let` 折叠或 DCE 消去；(3) 正确性更容易 runtime assert。

### 1.3 可自写性评估

核心 SSA 实现量估计：
- CFG 构造（结构化控制流）：~150 行 Python
- Dominator tree（Cooper-Harvey-Kennedy）：~50 行
- DF 计算：~30 行
- Phi 放置：~40 行
- Variable renaming：~80 行
- **合计 ~350 行**，完全可控。参考：LLVM `MemorySSA` 的 Python port
  <https://github.com/llvm/llvm-project/blob/main/llvm/lib/Transforms/Utils/SSAUpdater.cpp>
  （只展开 incremental 版本时 ~400 行 C++）。

### 1.4 Lean/Coq 领域参考实现

- **CompCert SSA**（Leroy et al., *Validating dominator trees for a fast and
  efficient ssa-based computation*, 2008）：Coq 形式化的 Cytron，但它是
  **后验验证器**（外部不可信产生 SSA，Coq 验证），不是从头构造。对我们
  参考价值：如何定义 SSA invariant 用于 assertion。
- **Vellvm**（Zhao et al., *Formalizing the LLVM Intermediate Representation
  for Verified Program Transformations*, POPL 2012）：Coq 中的 LLVM IR 形式
  语义，含 SSA。
- **Lean 4**：没有现成 SSA 基础设施（见 §3）。

---

## §2 Back-to-Back 测试工程实现参考

参考：`back2back-design.md` §3-7；`proof/cpp2lean/test_back2back.cpp`（已
存在的原型，94 行，使用**空格分隔文本**而非 JSON，需升级）。

### 2.1 C++ 侧序列化

**现状**：`test_back2back.cpp` 用 `std::cout << val << " " << prime`，`gen_b2b_lean.py`
产 `eval_back2back.lean` 输出 `Zp.mk val prime`。两端格式不同，比对靠人眼。

**升级目标**：两端统一产 JSON Lines（每行一个 record）。不引入第三方
JSON 库——CLPoly 是零依赖库，测试驱动用手写 emitter 即可：

```cpp
// proof/cpp2lean_v2/b2b_emit.hh（~50 行，手写 JSON emitter）
struct JSONEmit {
  std::ostream& os;
  void record(const char* fn, const char* input_json,
              const char* output_json) {
    os << R"({"fn":")" << fn << R"(","in":)" << input_json
       << R"(,"out":)" << output_json << "}\n";
  }
};
// Zp → "{\"val\":5,\"prime\":13}"
// SparsePolyZp → "[[deg,val],[deg,val],...]"
```

Harness 架构：
1. 一个 C++ `harness.cc` include 全部 `polynomial_factorize_*.hh`。
2. 读 `vectors.jsonl`（每行 `{fn, in}`）。
3. 分派到被测函数，捕获返回，emit `{fn, in, out}` 到 stdout。
4. 输出重定向到 `cpp_out.jsonl`。

**参考 `crosscheck_flint.hh`**（test/crosscheck_flint.hh）：CLPoly 现有跨库对比都走
内存比较 + `assert`，无文件 IO。b2b 是新增数据通路，必须引入最小 JSON。

### 2.2 Lean 侧序列化

Lean 4 的 `#eval` 需要目标类型 `Repr` 或自定义 `ToString`。方案：

```lean
-- 给 b2b 目标类型加 ToJson（手写，~30 行）
instance : ToJson Zp where
  toJson z := json%{"val": $(z.val), "prime": $(z.prime)}
-- Batteries / Std 有 Lean.Json；Mathlib 有 Lean.Data.Json.
#eval IO.println (toJson (__make_zp_ir 7 13)).pretty
```

执行：
```bash
lake env lean --run harness.lean > lean_out.jsonl
```

`--run` 比 `#eval` 更适合批量：文件带 `def main : IO Unit := ...`，stderr
分开，easier CI。参考：Mathlib 自己的 `test/` 目录用 `#guard_msgs` 和
`lake env lean` 双模式。

**require 参数**：测试向量里记 `{fn, in, require_proofs: ["by decide", ...]}`，
Lean harness 生成器把证明项 splice 进调用。

### 2.3 Diff 工具

Python 脚本 `b2b_diff.py`（~80 行）：
1. 两端按 `(fn, input_hash)` join。
2. 递归 dict/list 比较，浮点用 `math.isclose(rel_tol=1e-9)`（仅
   `__heuristic_starting_precision` 涉及 double）。
3. 不一致时输出 unified diff：C++ 值 vs Lean 值，标红 key path。
4. 报告格式同 `back2back-design.md` §7（`total: N/M PASS`）。

**不用 jsondiff 库**：输入格式固定、比较规则简单（几乎全整数+数组），
手写 50 行 Python 比引入依赖清晰。

---

## §3 Lean CFG/Graph 基础设施调研

**结论：Lean 4 / Mathlib / Batteries 都没有可直接复用的编译器级 CFG 或
SSA 数据结构。** 以下是搜索结果：

| 候选 | 路径 | 是否适用 |
|------|------|--------|
| `Mathlib.Combinatorics.SimpleGraph.*` | 40+ 文件 | 数学图论（匹配、着色、度数和），无基本块概念 |
| `Mathlib.Combinatorics.Digraph.Basic` | 有向图 | 仅顶点集+邻接函数，无 entry/exit/dom 概念 |
| `Mathlib.Combinatorics.Quiver.*` | Quiver 模型 | 范畴论 quiver（带 type indexing），不适合 CFG |
| `Mathlib`中 "Dominator" | 64 文件命中 | 全部是"dominated convergence"（分析测度），**0 个编译器意义** |
| `Mathlib` 中 "SSA" / "phi node" | 0 文件 | 无 |
| Batteries `Std.Data.*` | 容器 + HashMap | 无图结构 |
| CompCert-Lean / Vellvm-Lean | — | **不存在**（均为 Coq 项目） |

**直接后果**：
1. Pass 6 的 CFG / DomTree / DF 数据结构**必须在 Python 端实现**（翻译器内部），
   不是 Lean 端。我们也没打算把 SSA 变换搬到 Lean——MIR 是翻译器的中间表示，
   不导出到证明层。
2. 生成的 Lean 代码里不会出现 phi 节点（phi 在 MIR→Lean codegen 阶段下降为
   `let` / `match` / 尾递归参数），所以 Lean 层不需要 SSA 概念。

---

## §4 对 CLPoly v2 的 Pass 6-8 + b2b 的具体建议

### 4.1 Pass 6 `ssa_build`（MIR₀ 构造）

- **数据结构**（纯 Python，~200 行）：
  ```python
  @dataclass class BasicBlock: id: int; stmts: list[HirStmt]; succs: list[int]; preds: list[int]
  @dataclass class CFG: entry: int; blocks: dict[int, BasicBlock]
  @dataclass class DomTree: idom: dict[int, int]; children: dict[int, list[int]]
  @dataclass class PhiStmt: var: str; version: int; srcs: dict[int, int]  # pred_blk -> src_version
  ```
- **实现顺序**（强制单 Pass 分步，便于测试）：
  1. HIR₄ → CFG builder：专门处理结构化控制流。`while` → header+body+exit；
     `for` 先脱糖为 `{init; while(cond){body; step}}`；`break`/`continue`/`return`
     保留为特殊 terminator，Pass 7 才下降。
  2. CFG → DomTree：Cooper-Harvey-Kennedy 迭代算法。
  3. DomTree → DF：Cytron §1.1 §(c) 的 local+up 算法。
  4. Minimal phi 放置（迭代 DF⁺）。
  5. DFS 重命名（版本栈 + 回溯）。
- **Runtime invariant 检查**（每 Pass 出口 ~30 行 assert）：
  - 每变量每静态位置至多一次赋值（遍历所有 blocks 的 stmts 计 def）。
  - 每 phi 的 `len(srcs) == len(block.preds)`。
  - 每处 use 的版本号在 entry→use 路径上有 def 覆盖。

### 4.2 Pass 7 `loop_lower`（MIR₀ → MIR₁）

- **策略**：不做通用 loop → recursion，而是**利用结构化信息**。因为 Pass 6
  保留了"这个 block 来自 while header"的标签，直接：
  - 每个 loop 提一个 `partial def _loop_N (state_tuple) : result_type`。
  - state = 循环携带的 phi 目标变量。
  - body 继续/break/return 用 `Result : | Continue state | Break value | Return value` sum type
    承载（背下降为显式 flag）。
- 这样比"先 goto-style 再 loop-recognize"可靠——我们从 HIR 进 MIR 时还
  保有高层控制流信息，丢了 structure 再捡就浪费（v1 教训）。

### 4.3 Pass 8 `codegen`

- **MIR → Lean 映射**：
  - phi 在非循环上下文下降为 `let v := match prev_blk with ...`；在循环
    上下文，phi = 递归调用的参数。
  - `TailCallStmt loop_N args` → `loop_N args`。
  - `ValueStmt e` → `e`（`partial def` 末尾表达式）。
- **类型**：MIR 节点携带 Clang 类型（`translator-v2-plan.md` §3.3），codegen
  查 type-map（`proof/cpp2lean/class_map.py`）即可；禁止任何推断。
- **输出格式**：每函数一个 `def/partial def` 块；必要时 `namespace CLPoly.L1.Zp`。

### 4.4 Back-to-back 框架搭建顺序

1. **Week 5 产出（仅设计）**：`b2b-design.md` 定稿：JSON schema、harness 骨架、
   diff 脚本接口、vector 生成策略（继承 `back2back-design.md` §4.2 的 8
   个函数表）。
2. **Stage 2 落地顺序**（不在本调研范围，但建议）：
   - Step A：JSON emitter（C++ 50 行 + Lean 30 行）。
   - Step B：harness 驱动 + 向量文件格式。
   - Step C：`b2b_diff.py`。
   - Step D：1 个 trivial 函数（`__make_zp`）跑通管线，作为 Stage 2 验收点。
3. **vectors 来源**：优先从 CLPoly 现有 test/ 套件提取实例（`test_factorize_zp.cc`
   等），随机生成作为补充；初版不追求覆盖率，只要 "每函数 ≥ 3 个向量"。
4. **关键约束**：C++ 端 harness 必须 link 真实 CLPoly（可信基），Lean 端
   原语必须是**独立正确实现**（如 `polynomial_GCD_ir` 要么用 Mathlib
   `EuclideanDomain.gcd` 桥接，要么手写欧几里得）。现版 `gen_b2b_lean.py`
   的 `polynomial_GCD_ir = fun f g => f`（`return f`）是**伪原语**，必须
   替换——否则结果一致只能证明 "两端一样错"。
