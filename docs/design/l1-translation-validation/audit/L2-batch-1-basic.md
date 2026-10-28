# L2 — 单变量算法基础（批 1，8 个函数）

**审视日期**：2026-05-04
**模块**：L2 单变量算法基础（mod、平方根、范数、Mignotte、squarefree、p-th 根、subset product）
**状态**：8/8 ✅

---

## 1. `__symmetric_mod` ✅

**C++**：`polynomial_factorize_univar.hh:94-103`
**Lean**：`Corpus.lean:4687-4698`

```cpp
ZZ r;            ZZ::fdiv_r(r, a, m);   // r ∈ [0, m)
ZZ half;         ZZ::fdiv_q(half, m, ZZ(2));
if (r > half) r -= m;
return r;
```

Lean：r_1 := 0; r_2 := ZZ.fdiv_r r_1 a m; half_1 := 0; half_2 := ZZ.fdiv_q half_1 m 2; if r_2 > half_2 → r_3 = r_2 - m; else r_2。
✅（fdiv_r/fdiv_q 通过 Pass 2 ref-out 重写为返回值）

---

## 2. `__binomial` ✅

**C++**：`polynomial_factorize_univar.hh:139-151`
**Lean**：`Corpus.lean:38-55`（含 loop:28-36）

```cpp
if (k < 0 || k > n) return 0;
if (k == 0 || k == n) return 1;
if (k > n - k) k = n - k;
result = 1;
for (i = 0; i < k; ++i) { result *= (n-i); result /= (i+1); }
return result;
```

Lean: 三层 if（boundaries / k=0|n / k>n-k swap）→ bb_11 → result=1, i=0, loop。
Loop: while i < k_2: result_3 := result_2 * (n-i); result_4 := result_3 / (i+1); i++; tail.
✅

---

## 3. `__isqrt_ceil` ✅

**C++**：`polynomial_factorize_univar.hh:154-173`（Newton 整数平方根）
**Lean**：`Corpus.lean:1514-1536`（含 loop:1501-1512）

```cpp
if (n <= 0) return 0;
size_t bits = n.sizeinbase(2);
ZZ x(1); x <<= (bits+1)/2;
while (true) {
    ZZ x1 = (x + n/x) / 2;
    if (x1 >= x) break;
    x = x1;
}
if (x*x < n) x += 1;
return x;
```

Lean：if n ≤ 0 → 0；else bits = sizeinbase n 2；x_1=1；x_2 = x_1 << (bits+1)/2；调 _loop_..._0_ir（while(true) + break）；后置 bb_7：if x_3*x_3 < n → x_5 = x_3+1 else x_3。
✅

⚠️ stub 依赖：`ZZ.sizeinbase` 在 Lean Model 当前是 stub（返回 0），B2B 时补真实位数计算。

---

## 4. `__mignotte_bound` ✅

**C++**：`polynomial_factorize_univar.hh:176-184`（B = C(n, ⌊n/2⌋) · ‖f‖₂）
**Lean**：`Corpus.lean:1987-1994`

```cpp
n = get_deg(f); binom = __binomial(n, n/2);
norm_sq = __upoly_norm_l2_sq(f); norm = __isqrt_ceil(norm_sq);
return binom * norm;
```

Lean：n_1 := get_deg f; binom_1 := __binomial_ir n_1 (n_1 / 2); norm_sq_1 := upoly_norm_l2_sq; norm_1 := isqrt_ceil; binom_1 * norm_1。
✅

---

## 5. `__heuristic_starting_precision` ✅

**C++**：`polynomial_factorize_univar.hh:607-631`（FLINT 启发式精度 + Mignotte 上限）
**Lean**：`Corpus.lean:1477-1499`（含 loop:1469-1475）

C++ 算法：
- logp = log(p); min_b = sizeinbase(p,2); N = f.size()-1
- a_h = max(1, ceil((2.5*r + min_b)*log(2)/logp + log(N+1)/(2*logp)))
- B_mig = mignotte_bound(f); lc_f = front!.snd; if lc_f<0: lc_f=-lc_f
- target = 2*lc_f*B_mig; a_mig=0; pa=1
- while pa <= target: pa *= p; ++a_mig
- return {min(a_mig, a_h), a_mig}

Lean：完全对应。logp/min_b/N/a_h_d/a_h 顺序计算；B_mig + lc_f 符号修正；进 bb_3 → 计算 target，loop _0_ir 推进 pa 直到 > target，返回 (min(a_mig, a_h), a_mig)。
✅

⚠️ stub 依赖：`Nat.log : Float → Float`（覆盖 Mathlib Nat.log，Lean Model 中是 stub `1.0`）；`ZZ.sizeinbase` 同上。B2B 时补真实 log。

---

## 6. `__extract_pth_root` ✅

**C++**：`polynomial_factorize_zp.hh:109-120`
**Lean**：`Corpus.lean:604-612`（含 loop:590-602）

```cpp
p = front!.snd.prime();
upolynomial_<Zp> g; g.reserve(f.size());
for (auto& term : f) {
    assert(term.first.deg() % p == 0);
    g.push_back({umonomial(term.first.deg() / p), term.second});
}
return g;
```

Lean：p_1 := front!.snd.prime; g_1 = empty; loop range-for: assert deg%p==0; push (deg/p, term.snd) into g。
✅

---

## 7. `__squarefree_Zp` ✅

**C++**：`polynomial_factorize_zp.hh:122-180`（Yun-style，p-th root 处理 derivative=0 退化）
**Lean**：`Corpus.lean:4622-4664`（含 loops 0/1/2_ir + recursive call）

C++ 控制流：
1. 若 derivative(f) 为空 → 提 p-th root → make_monic → 递归 squarefree → 每个 (s, e) push 为 (s, e·p)
2. 否则：
   - c = GCD(f, f'); w = f/c; normalize w
   - i = 1; while w 非空 && deg(w) > 0:
     - y = GCD(w, c); z = w/y; normalize z
     - 若 z 非空 && deg(z) > 0：make_monic z; push (z, i)
     - c = c/y; normalize c; w = y; ++i
   - 尾部：若 c 非空 && deg(c) > 0：提 p-th root → make_monic → 递归 squarefree → push 为 (s, e·p)

Lean：3 个 loop（_0/_1/_2 分别处理两个 sub.push 循环和主 while 循环）+ 主体 if-else + 递归 self-call。控制流逐分支匹配 C++。
✅ 包括 derivative=0 退化分支、主 while 循环、尾部高次项处理。

---

## 8. `__subset_product_mod` ✅

**C++**：`polynomial_factorize_univar.hh:725-740`
**Lean**：`Corpus.lean:4678-4685`（含 loop:4666-4676）

```cpp
prod.push_back({umonomial(0), lc_f});
for (size_t idx : subset) {
    prod = prod * factors[idx];
    prod.normalization();
    __upoly_mod_coeff(prod, m);
}
return __upoly_symmetric_mod(prod, m);
```

Lean：prod_2 = push (0, lc_f)；range-for over subset：prod_4 = prod_3 * factors[idx]；prod_5 = normalization；prod_6 = mod_coeff prod_5 m；尾部 __upoly_symmetric_mod prod_3 m。
✅

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| __symmetric_mod | ✅ | 无 |
| __binomial | ✅ | 无 |
| __isqrt_ceil | ✅ | ZZ.sizeinbase stub（B2B 补） |
| __mignotte_bound | ✅ | 无 |
| __heuristic_starting_precision | ✅ | Nat.log + sizeinbase stub（B2B 补） |
| __extract_pth_root | ✅ | 无 |
| __squarefree_Zp | ✅ | 无 |
| __subset_product_mod | ✅ | 无 |

L2 batch 1：8/8 翻译忠实。无翻译器 bug。共 3 个 stub 依赖（`ZZ.sizeinbase`、`Nat.log`(Float)、`pair_vec_div`5）需 B2B 阶段补 Lean Model 真实实现。
