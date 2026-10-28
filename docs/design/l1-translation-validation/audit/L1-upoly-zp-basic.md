# L1 — `__upoly_make_monic` / `__upoly_subtract_x` / `__upoly_subtract_one` / `__upoly_mod`

**审视日期**：2026-05-03
**模块**：单变量 Zp 多项式基础操作
**状态**：4/4 ✅

---

## 1. `__upoly_make_monic` ✅

**C++ 源**：`polynomial_factorize_zp.hh:27-36`
**Lean 翻译**：`Generated/Corpus.lean:4903-4916`
**Loop**：`_loop___upoly_make_monic_0_ir` (Corpus.lean:4893-4901)

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 名 | `__upoly_make_monic(upolynomial_<Zp>& f)` | `__upoly_make_monic_ir (f : SparsePolyZp)` | ✅ |
| 返回 | `Zp` (mutates f) | `(Zp × SparsePolyZp)` | ✅（Pass 2 ref→tuple） |

### 行为对照

| C++ | Lean | 一致 |
|-----|------|------|
| `assert(!f.empty())` | `-- require (h_assert)` | ✅ |
| `Zp lc = f.front().second` | `let lc_1 := (front! f).snd` | ✅ |
| `if (lc.number() == 1) return lc` | `if lc_1.val == 1 then (lc_1, f)` | ✅ |
| `Zp lc_inv = lc.inv()` | `let lc_inv_1 := Zp.inv lc_1` | ✅ |
| `for (auto& term : f) term.second *= lc_inv` | `_loop_..._0_ir`（数组逐项 `term.snd * lc_inv` + `Array.set!`）| ✅ |
| `return lc` | `(lc_1, f_1)` | ✅ |

### 偏差

无。

---

## 2. `__upoly_subtract_x` ✅

**C++ 源**：`polynomial_factorize_zp.hh:182-214`（含尾部 `normalization()`）
**Lean 翻译**：`Generated/Corpus.lean:5136-5151`
**Loop**：`_loop___upoly_subtract_x_0_ir` (Corpus.lean:5107-5135)

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 参数 | `(const upolynomial_<Zp>& h, uint64_t p)` | `(h : SparsePolyZp) (p : UInt64)` | ✅ |
| 返回 | `upolynomial_<Zp>` | `SparsePolyZp` | ✅ |

### 行为对照

C++ 算法：从 `h` 减去 monomial `x`，结果保存到 `result`。

| C++ 行 | Lean 翻译 | 一致 |
|--------|---------|------|
| `result.reserve(h.size()+1)` | (Lean Array 自动扩容，无需 reserve) | ✅（语义保留） |
| `bool inserted = false` | `let inserted_1 : Bool := false` | ✅ |
| `for (auto& term : h)` | `_loop_..._0_ir` 循环 | ✅ |
| `if (!inserted && term.first.deg() < 1) { push (1, p-1); inserted = true; }` | line 5127-5130：相同条件 + push + inserted_3=true | ✅ |
| `if (term.first.deg() == 1) { new_c = term.second - 1; if (new_c != 0) push (1, new_c); inserted = true; continue; }` | bb_7 + bb_13：deg==1 时计算 new_c，按需 push，bb_13 强制 inserted=true，跳到下一轮（不 push 原 term） | ✅ |
| 否则 `result.push_back(term)` | else 分支 push term_1 | ✅ |
| `if (!inserted) push (1, p-1)` | 循环外 `if !inserted_2 then push (1, p-1)` | ✅ |
| `result.normalization()` | `bb_17` 内 `SparsePolyZp.normalization` | ✅（无论 inserted 与否都 normalize） |

### 偏差

- **`reserve` 缺失**：C++ 用 `result.reserve(h.size()+1)` 预分配。Lean Array 自动扩容，无 `reserve` 概念。
  - 影响：仅性能，不影响结果。
  - 处理：作为已知不影响语义的差异。

---

## 3. `__upoly_subtract_one` ✅

**C++ 源**：`polynomial_factorize_zp.hh:217-244`
**Lean 翻译**：`Generated/Corpus.lean:5090-5105`
**Loop**：`_loop___upoly_subtract_one_0_ir` (Corpus.lean:5068-5088)

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 参数 | `(const upolynomial_<Zp>& h, uint64_t p)` | `(h : SparsePolyZp) (p : UInt64)` | ✅ |
| 返回 | `upolynomial_<Zp>` | `SparsePolyZp` | ✅ |

### 行为对照

C++ 算法：从 `h` 减去常数 1。

| C++ 行 | Lean 翻译 | 一致 |
|--------|---------|------|
| `result.reserve(h.size()+1)` | (Lean Array 自动) | ✅ |
| `bool found = false` | `let found_1 : Bool := false` | ✅ |
| `for (auto& term : h)` | `_loop_..._0_ir` | ✅ |
| `if (deg == 0) { new_c = term.snd - 1; if != 0 push (0, new_c); found = true; }` | bb_10 设 found=true, conditional push | ✅ |
| `else result.push_back(term)` | else 分支 push term_1 | ✅ |
| `if (!found) { push (0, p-1); result.normalization(); }` | line 5100-5103：仅 `!found` 分支 push 后 `normalization` | ✅ |

### 关键细节确认

C++ 的 **`normalization()` 仅在 `!found` 分支内调**（line 241），**`found` 分支不调**。Lean 1:1 反映：line 5100-5103 仅在 `!found_2` 时执行 `normalization`，`else` 分支直接返回。

### 偏差

无。

---

## 4. `__upoly_mod` ✅（含小问题）

**C++ 源**：依赖 `pair_vec_div`（外部 CLPoly 原语）
**Lean 翻译**：`Generated/Corpus.lean:4918-4924`

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 名 | `__upoly_mod(...)` 等价表达：`pair_vec_div(_, _, f, g, comp)` 取 r | `__upoly_mod_ir (f : SparsePolyZp) (g : SparsePolyZp) : SparsePolyZp` | ✅ |
| 返回 | `SparsePolyZp` (即 r) | `SparsePolyZp` | ✅ |

### 行为对照

```lean
let q_1 : SparsePolyZp := (SparsePolyZp.empty)
let r_1 : SparsePolyZp := (SparsePolyZp.empty)
let __refret_0_1 := (pair_vec_div5 q_1 r_1 f g (SparsePolyZp.comp f))
let _q_2 := __refret_0_1.fst
let r_2 : SparsePolyZp := __refret_0_1.snd
r_2
```

| C++ 操作 | Lean 翻译 | 一致 |
|---------|---------|------|
| `pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp())` | `pair_vec_div5 q r f g (SparsePolyZp.comp f)` | ✅ |
| 取 `r` 作为返回 | `r_2 := __refret.snd` | ✅ |

### 偏差

- ⚠️ **stub 依赖**：`pair_vec_div5` 在 Lean Model 当前是 stub `(default, default)`（不实现真实长除法）。
  - 影响：B2B 测试时输出错（`r` 总是 default）。
  - 处理：留待 B2B 阶段补真实实现。
  - 不是翻译器 bug，是 Lean Model 原语未实现。

- ⚠️ **`SparsePolyZp.comp` 返回类型**：阶段 G2 改为 `Lex` (= Unit)。Pass 5 emit `(SparsePolyZp.comp f)` 作为 `pair_vec_div5` 第 5 参数。`pair_vec_div5` Lean 签名 `(_comp : β) : α × α` 接受任意 β。✅ 类型一致。

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| `__upoly_make_monic` | ✅ | 无 |
| `__upoly_subtract_x` | ✅ | reserve 缺（性能） |
| `__upoly_subtract_one` | ✅ | 无 |
| `__upoly_mod` | ✅ | stub 依赖（B2B 时补） |

L1 共 4 个函数，4/4 翻译忠实。无翻译器 bug。
