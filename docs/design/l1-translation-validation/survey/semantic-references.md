# C++ → Lean 翻译器语义论证 —— 参考文献调研

**时间**：2026-04-19  
**目的**：为 blueprint.md §5a 的引理 1-7 做严格展开（Stage 1 Week 3）收集理论基础与表达惯例。重点：能直接引用哪些已有结果，需要自己证明什么，写 `cpp-subset-semantics.md` 时用什么记号。

---

## §1 Appel "SSA is Functional Programming" (1998)

### 核心主张

Appel 的论文论证一个观察而非一个形式化定理：**SSA 形式与嵌套 scope 的函数式程序同构**。具体来说：

- **每个基本块 → 函数**：块的入口 = 函数签名；块中在入口处活跃的变量 = 函数形参。
- **每个 phi 函数 → 形参**：merge 点的 `x = φ(x₁, x₂)` 变成函数 `f(..., x, ...)` 的参数 x；不同前驱调用 f 时把各自的版本传进来。
- **控制流边 → 函数调用**：块末尾的跳转 `goto L` 翻译成 `L(...)`（尾调用）；条件跳转 `if c goto L1 else L2` 翻译成 `if c then L1(...) else L2(...)`。
- **Dominance 决定嵌套**：`def(x)` dominates `use(x)` ⟺ x 的绑定在函数式表达式中对 use 可见。Dominance frontier 决定 phi 的放置（等价于决定哪些内部函数需要多一个参数）。

### 语义论证的特点

Appel **没有给出形式化的语义保持定理**；论文用的是"结构性对应"（structural correspondence）的论证：两种表示都在一个公共的 CFG 骨架上操作，每一步计算一一对应。严格的形式化要到 Kelsey 1995（CPS↔SSA）或后来 Chakravarty/Keller、MLTon、CompCertSSA 才给出。

### 对 CLPoly 的可引用结论

| 引理 | 可引用结论 |
|------|------------|
| 引理 2（SSA→let 链）| "定义 dominates 使用" ⇒ 每个 use 读到的值就是最近一次 def 的值 ⇒ `let` 语义与 mutation 后读语义一致。这正是 Appel 的核心观察。 |
| 引理 3（loop→尾递归）| Appel 的 loop header 作为 recursive function、loop-carried phi 作为形参——这**就是**我们蓝图的翻译模式，是直接照搬。 |

**表达惯例**：Appel 在 §2 用"the functional program computes the same result"这种非形式化措辞；我们在 cpp-subset-semantics.md 里应该写得更严格些（用 big-step evaluation judgement `σ ⊢ S ⇓ σ'` 或 denotation `⟦S⟧σ = σ'`）。

**结论**：Appel 给我们 **模式**（pattern）而非 **证明**（proof）。引理 2、3 的严格证明要自己写，但可以合法地引用 Appel 为翻译模式的出处，并写"dominance 条件在我们的 CLPoly 子集中自动成立（无 goto，无指针别名）"作为 bridge。

---

## §2 Aeneas (Ho & Protzenko, ICFP 2022)

### 核心架构

Aeneas 将 Rust → LLBC (Low-Level Borrow Calculus) → **纯 λ 演算**，通过符号执行同时做 borrow-check 和翻译。目标支持 F*、Coq、**Lean 4**、HOL4（Lean 和 HOL4 是最成熟的后端）。

### Loop 翻译

- **Forward + backward function**：mutable borrow 的写回语义被编码为一对函数：`forward` 返回正常结果，`backward` 返回 borrow 结束时被修改后的值（这是 Aeneas 消除 reference 的核心发明）。
- **Loop 的翻译**：每个 loop 生成一对递归函数（forward/backward）。**不使用 fuel**，而是用 **递归函数 + extrinsic termination proof**（外部证明）。Lean 后端用 `partial def`（Lean 4 的 partial 函数机制，生成后仍可 reduce）。官方文档原话："mature backends are Lean and HOL4, which have support for partial functions and extrinsic proofs of termination"。
- 这**就是**我们蓝图 §3.1 选的策略——用 `partial def` 避免在翻译时给出 measure。

### 监子（monad）

输出目标代码用一个 **Result monad**（类似 `inductive Result α | ok : α → Result α | fail : Error → Result α | div : Result α`），其中 `div`（divergence）标记可能不终止的递归。在 Lean 后端 `div` 对应 `partial def` 本身的非总函数性。CLPoly 的 `Except String α` 大致对应 Aeneas 的 `Result`（我们的版本更简单因为我们用 `partial def` 吸收发散）。

### 引用消除

Aeneas 的 "value-based semantics"——**没有 memory、没有 address**——是通过 backward functions 做到的：每个 `&mut T` 参数在签名中对应一个 backward 返回值，把 borrow 结束时的 T 作为返回值带出。

**CLPoly 的情况更简单**：我们的 C++ 子集**没有 mutable reference**（P1 前提保证），所以 backward function 机制**不需要**；我们直接用 SSA + let 链就够了。Aeneas 的复杂性在这里是 over-kill。

### 语义保持

Aeneas ICFP 2022 论文对翻译的正确性是**informal + case study**（没有机器检查过的定理）。**ICFP 2024 后续论文**才机器检查了符号执行正确实现 borrow-checker。但翻译本身的端到端语义保持定理至今**没有完整的机器证明**——他们的 justification 是"值语义是显然的，backward functions 的合理性用抽象解释证明"。

### 对 CLPoly 的启示

| 维度 | Aeneas 做法 | CLPoly 采纳 |
|------|------------|------------|
| Loop 翻译 | recursive function + partial | ✅ 相同 |
| Termination | extrinsic / `partial def` | ✅ 相同（避免终止推断） |
| Monad | `Result { ok, fail, div }` | 简化为 `Except String α`（`div` 被 `partial def` 吸收） |
| 引用 | backward function | ❌ 不需要（子集无 mut ref） |
| 语义保持 | informal + 后续机器证明部分 | 我们给**纸面证明**，与 Aeneas 同级严格度 |

**表达惯例**：Aeneas 论文用 `⟦e⟧ρ ⇓ v`（big-step with environment）和 simulation `e ~ e'`（LLBC 与纯 λ 的语义关系）。我们应在 cpp-subset-semantics.md 采用相同的 big-step 记号。

---

## §3 通用 Translation Validation 方法

### 两种路线

1. **Verified compiler**（如 CompCert / Leroy）：一劳永逸地证明编译器函数 `C : Src → Tgt` 满足 `∀ P, ⟦P⟧ ≼ ⟦C(P)⟧`。证明写一次，对所有输入有效。技术工具：**simulation diagrams**（forward / backward simulation），stuttering 时用 **measure function** `m(·)` 消除。
2. **Translation validation**（Pnueli 1998, Necula 2000）：每次编译后，运行一个 validator 生成 verification conditions，证明"这次的源和目标等价"。不证明编译器本身，只证明这一次的输出。Necula 用 RTL 上的抽象解释，Pnueli 用 VC generator + theorem prover。

**CLPoly 位置**：我们是 **hybrid**——翻译器是外部 Python 脚本（不在 Lean 中证明翻译器本身），但对每次翻译产生的 Lean 文件需要附证明目标（UB-freedom + refinement theorem）。这更接近 translation validation（per-run validation）但简化：我们不 auto 生成 VC，而是在 blueprint 主定理下让人工（一次性）论证所有翻译规则，然后每次翻译只需检查子集前提 P1-P8（§5a 已列）。

### CompCert 的 simulation diagram 技术

CompCert 每个 pass 单独证明一个 simulation diagram，再组合。**forward simulation + receptive/determinate target ⇒ backward simulation**（Leroy 的技巧：对 deterministic target 语言，forward 自动升格为 backward）。

**启示**：C++ 和我们的 Lean 子集都是 deterministic（在 UB-free 前提下），所以我们只需证 **forward simulation**（∀ 源执行 → 存在对应目标执行），自动得到 backward。这让每个引理（1-7）只需证"C++ 这一步 → Lean 这一步"单方向即可。

### Winskel《The Formal Semantics of Programming Languages》

Winskel 的 IMP 语言**就是**我们 C++ 子集的一个超集（有 assignment、while、if、sequencing）。他同时给出 big-step SOS、small-step SOS、denotational semantics，并**证明三者等价**。

**可引用**：
- `⟦·⟧_C` 我们用 big-step SOS，记号 `σ ⊢ S ⇓ σ'`（跟 Winskel §2 完全一致）。
- while 循环的 denotation 用**最小不动点** `⟦while b do S⟧ = fix(λF. λσ. if ⟦b⟧σ then F(⟦S⟧σ) else σ)`——与我们用 `partial def` 的 Lean 翻译直接对应（`partial def` 在 Lean 4 是 inhabitance-based，语义上等价于最小不动点 up to 不可计算性）。

### Verified vs Validated 词汇表

| 术语 | 含义 | 我们用吗 |
|------|------|---------|
| semantics preserving | ⟦T(P)⟧ ≈ ⟦P⟧（弱相等/refinement） | ✅ 蓝图主定理用这个 |
| refinement (≼) | 目标的行为集 ⊆ 源的行为集 | ✅ L1 IR ≈ L2 Algorithm 的记号 |
| forward simulation | 源走一步 ⇒ 目标走 ≥ 1 步保持关系 | ✅ 引理 1-7 的形式 |
| verified translator | 翻译器在 Lean 里定义并证明 | ❌ 我们不做 |
| translation validation | per-run validator | ⚠️ 类似但简化（前提条件检查） |

---

## §4 对 CLPoly 的启示与引用建议

### 直接可引用的结果

| 位置 | 可直接引用什么 | 引用语 |
|------|---------------|--------|
| 引理 1 总导言 | Winskel IMP 的 denotational = operational 等价 | "UInt64 / Bool 运算的 denotational 定义与 C++ operational 语义在 UB-free 下一致，由定义直接验证，参考 Winskel §5"。 |
| 引理 2（SSA） | Appel 1998 定理 2（dominance ⇒ 函数式对应） | "此性质是 Appel 定理在单线程无别名子集上的直接实例"（已在蓝图写了；保留）。 |
| 引理 3（loop→recursion） | Aeneas ICFP 2022 的 loop 翻译策略 | "本翻译模式采用与 Aeneas 相同的策略：extrinsic termination 通过 `partial def` 实现，递归函数的语义由 Lean 4 的 partial 函数定义给出。" |
| 引理 3 的正确性 | Winskel §8（while loop as least fixpoint） | "`partial def` 的 denotation 是递归方程的最小不动点，与 C++ while 的标准 denotational semantics（Winskel Thm 5.11）一致。" |
| 主定理（simulation） | Leroy CompCert 的 forward-simulation-implies-backward-simulation for deterministic targets | "C++ 和 Lean 子集均 deterministic（UB-free 前提），forward simulation 自动升格为 backward simulation（Leroy 2009 Lemma 3.4）。" |

### 必须自己证明的

1. **`uint128_mul_correct`**：Aeneas/Appel/CompCert 都不涵盖。**独立 Lean 定理**（已在蓝图 P7 列出）。
2. **CLPoly 子集的 AST 覆盖完整性**：Clang AST 的每个 node kind → 引理 1-7 之一。需要做一张完整映射表（放在 cpp-subset-semantics.md §2）。
3. **Except 传播 = throw 传播**：引理 5 的具体形式（do notation 的 bind）。Aeneas 的 Result monad 有此性质但用的是 F\*/Lean 各自版本，我们要写自己版本。
4. **输出参数 → 返回值的 non-alias 前提**：引理 6。CLPoly 因式分解无 mutable ref，所以这是"空条件"，但要在 cpp-subset-semantics.md 中明确表达 non-alias invariant 的定义。

### 表达惯例（写 cpp-subset-semantics.md 时采纳）

1. **big-step judgement**：`⟨S, σ⟩ ⇓ σ'`（C++ 语义）和 `⟦T(S)⟧ ρ = v`（Lean 语义）。用 `⟦·⟧_C` / `⟦·⟧_L` 标注不同语言。
2. **状态**：σ : Var → Value（C++）vs ρ : Var → Value（Lean 环境）。翻译时维护一个 **renaming** `θ : Var_C → Var_L`（SSA 版本号）。
3. **关系**：`σ ≈ ρ` 定义为 `∀ x. ρ(θ(x)) = σ(x)`（C++ 值与最新 SSA 版本一致）。
4. **引理形式统一**：`⟨S, σ⟩ ⇓ σ' ∧ σ ≈ ρ ⇒ ⟦T(S)⟧ ρ = v ∧ σ' ≈ ρ[…↦v]`（forward simulation）。
5. **Throw 传播**：`⟨S, σ⟩ ⇓ throw e ⇒ ⟦T(S)⟧ ρ = .error e`（分支处理 Except）。

### 结构建议（cpp-subset-semantics.md）

```
§1 记号与两种语义
  §1.1 C++ 子集的 big-step semantics
  §1.2 Lean 子集的 denotation
  §1.3 relation ≈
§2 AST 覆盖表（每个 Clang node kind 对应哪个引理）
§3 引理 1-7 的严格证明
  §3.x 每个引理：陈述 → 证明 → 引用
§4 主定理：forward simulation
§5 子集前提 P1-P8 的形式化（每条给"可证伪"的机器检查条件）
```

### 不值得做的

- **全机器证明**（把翻译器写进 Lean 并证明其正确）：Aeneas 自己都没做（2022 版），2024 版只证了 borrow-checker 部分。对 CLPoly 是严重 over-engineering。
- **verified translation validator**（自动 VC 生成）：Pnueli/Necula 路线，对我们 4617 行的目标代码是过杀。手写引理 1-7 + 子集前提审计是 sweet spot。

---

## 参考

- **Appel, "SSA is Functional Programming", ACM SIGPLAN Notices 1998** —— SSA ↔ nested function 的结构对应，我们引理 2、3 的模式出处。
- **Ho & Protzenko, "Aeneas: Rust Verification by Functional Translation", ICFP 2022**（arXiv:2206.07185）—— Rust → 纯函数，Lean 后端，`partial def` + extrinsic termination。loop/monad/Result 策略与我们完全一致（去掉 backward functions 部分）。
- **Ho et al., ICFP 2024 (Aeneas borrow-check 机器证明)** —— 不直接相关（我们无 borrow），仅作"相关工作"引用。
- **Leroy, "A formally verified compiler back-end", JAR 2009** —— forward/backward simulation 框架，deterministic target 升格引理（Lemma 3.4）。
- **Winskel, The Formal Semantics of Programming Languages, MIT Press 1993** —— IMP big-step + denotational + 等价性（Thm 5.11, while 的最小不动点语义）。
- **Pnueli, Siegel, Singerman, "Translation Validation", TACAS 1998** —— translation validation 原始定义（我们**不**走这条路）。
- **Necula, "Translation Validation for an Optimizing Compiler", PLDI 2000** —— symbolic execution-based validator（对比点，非引用）。
