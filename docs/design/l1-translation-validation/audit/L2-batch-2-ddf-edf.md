# L2 — DDF / EDF / factor_Zp / select_prime（批 2，4 个函数）

**审视日期**：2026-05-04
**模块**：Zp 上不可约分解管道 + 素数选择
**状态**：4/4 ✅

---

## 1. `__ddf_Zp` ✅

**C++**：`polynomial_factorize_zp.hh:248-292`
**Lean**：`Corpus.lean:213-239`（含 loop:186-211）

C++ 算法（不同度因子化 Distinct-Degree Factorization）：
```cpp
h = x; f_star = f;
for (d = 1; ; ++d) {
    if (deg(f_star) < 2*d) break;
    h = upoly_powmod(h, p, f_star);
    h_minus_x = upoly_subtract_x(h, p);
    gd = polynomial_GCD(h_minus_x, f_star);
    if (deg(gd) > 0) {
        make_monic(gd); result.push({gd, d});
        f_star = f_star / gd; normalize;
        h = upoly_mod(h, f_star);
    }
}
if (deg(f_star) > 0) { make_monic; result.push({f_star, get_deg(f_star)}); }
return result;
```

Lean：完全对应。
- 入口：assert; p_1 = front!.snd.prime; result_1 = []; h_1 = [(1, 1_zp)]; f_star_1 = f; d_1 = 1
- 主 loop _0_ir（`while(true)` 翻译为递归）：
  - if 2*d > deg(f_star): break (kind=1)
  - 否则：h_3 = powmod h_2 p_1 f_star_2; h_minus_x_1 = subtract_x h_3 p_1; gd_1 = GCD(h_minus_x, f_star)
  - if gd 非空 && deg > 0: make_monic gd → push (gd_2, d_2); f_new = pair_vec_div empty f_star gd; f_star_4 = normalize f_new; h_4 = upoly_mod h_3 f_star_4 → bb_11 (d++ + tail)
  - 否则 → bb_11 直接 d++ 不做修改
- 尾部 bb_4：if deg(f_star) > 0 → make_monic → push (f_star_6, deg(f_star_6)) → bb_14 return

✅ 注：Lean 用 `if true then ... else ...` 模型 `while(true)`；break 通过 `(1, ...)` 返回；continue 通过递归。

---

## 2. `__edf_Zp` ✅（含若干工程注脚）

**C++**：`polynomial_factorize_zp.hh:295-354`
**Lean**：`Corpus.lean:281-322`（含 loops _0_ir:241-248 + _1_ir:250-279）

C++ 算法（等度因子化 Equal-Degree Factorization，Cantor-Zassenhaus）：
```cpp
if (deg(f) == d) { copy + make_monic + push; return; }
if (deg(f) <= 0) return;
while (true) {
    r = upoly_random(n, p, rng);
    if (r empty) continue;
    if (p == 2) {  // trace map
        g = upoly_mod(r, f);
        for (i=1; i<d; ++i) { g = g*g + r; g = upoly_mod(g, f); }
        g = polynomial_GCD(g, f);
    } else {  // 奇特征
        exp = (p^d - 1) / 2;
        g_pow = upoly_powmod(r, exp, f);
        g_pow_minus_1 = upoly_subtract_one(g_pow, p);
        g = polynomial_GCD(g_pow_minus_1, f);
    }
    if (deg(g) > 0 && deg(g) < deg(f)) {
        h_part = f / g; normalize; make_monic g; make_monic h_part;
        edf_Zp(result, g, d, rng); edf_Zp(result, h_part, d, rng);
        return;
    }
}
```

Lean：分两层
- 顶层 `__edf_Zp_ir`：
  - if deg(f).toUInt64 == d → make_monic; push; return (result_1, rng)
  - elif deg(f) ≤ 0 → return (result, rng)
  - else 调 _loop___edf_Zp_1_ir 处理 while 循环 → 取 g_8, rng_2
  - kind==1 路径：h_part = f/g; normalize; make_monic g; make_monic h_part; 递归 edf_Zp(result, g, d, rng_2) → result_2, rng_3; 递归 edf_Zp(result_2, h_part, d, rng_3) → 结果。
- _loop___edf_Zp_1_ir（while(true) loop）:
  - random → (r_1, rng_2)
  - if r 空 → continue (用 rng_2)
  - else: p==2 分支调 _loop___edf_Zp_0_ir 计算 trace map；奇特征分支算 exp, powmod, subtract_one；GCD → bb_18
  - bb_18: if deg(g) > 0 && < deg(f) → break (1, g, rng_2, rng_2)；else continue

✅ 控制流 1:1 完整。

⚠️ 工程注脚（非翻译 bug）：
- `_loop___edf_Zp_1_ir` 返回 `(Int64 × SparsePolyZp × Rng × Rng)` — Pass 6 SSA 在 break 路径返回两个 Rng 副本（`(1, g_8, rng_2, rng_2)` 同值）。caller 仅用第二份 rng_2 做递归（kind==0 死分支用第一份 rng_1，但 while-true 永远进 kind==1）。语义无误，仅冗余。
- 与 Plan A 修复联动：`__upoly_random` 现已正确推进 rng；EDF 内每轮 random 都用新 rng。

---

## 3. `__factor_Zp` ✅

**C++**：`polynomial_factorize_zp.hh:357-392`
**Lean**：`Corpus.lean:663-684`（含 loops _0_ir/_1_ir/_2_ir + cmp lambda _1_ir）

C++ 算法（Zp 上完整不可约分解）：
```cpp
assert non-empty; p = front!.snd.prime
if (deg(f) <= 0) return {front!.snd, []};
lc = make_monic(f);
sqf = squarefree_Zp(f);
result = [];
rng = mt19937(42);
for (sj_ej : sqf) {
    ddf = ddf_Zp(sj_ej.first);
    for (gk_dk : ddf) {
        factors_k = []; edf_Zp(factors_k, gk, dk, rng);
        for (hi : factors_k) result.push({hi, sj_ej.second});
    }
}
sort(result, by deg ascending);
return {lc, result};
```

Lean：
- 入口：assert; p_1 = front!.snd.prime（unused → `_p_1`）；if deg(f) ≤ 0 → (front!.snd, [])
- else: make_monic f → (lc_1, f_1); sqf_1 = squarefree_Zp f_1; result_1 = []; rng_1 = Rng.mk 42
- 三层嵌套循环：
  - _2_ir（外层 sqf）→ 内调 _1_ir
  - _1_ir（中层 ddf）→ 内调 edf_Zp(rng) → 内调 _0_ir
  - _0_ir（内层 factors_k）→ result.push (hi_1, sj.snd)
- result_6 = Array.sort result_2 (cmp lambda)
- return (lc_1, result_6)

✅ 三层嵌套 + rng 通过递归调用链 (rng_1 → rng_4 → rng_3 →) 正确推进。

---

## 4. `__select_prime` ✅（结构层审视）

**C++**：`polynomial_factorize_univar.hh:1388-1474+`
**Lean**：`Corpus.lean:4126-4160`（主体）+ 多个 sub-loops

C++ 算法（素数选择 with min-factors heuristic）：
- 初始化 best, best_count = SIZE_MAX, max_tries = 3, rng(42)
- p = use_large_prime ? 2^64-59 : 2
- next_p = lambda(use_large_prime ? prev_prime_64 : next_prime_64)
- for (tried = 0; tried < max_tries; p = next_p(p)):
  - 跳过 lc_f mod p == 0
  - fp = polynomial_mod(f, p); 跳过 deg 退化
  - 跳过 squarefree 失败 (gcd(fp, fp') != 1)
  - ++tried
  - make_monic fp; ddf = ddf_Zp(fp); irr_factors = []
  - for gk_dk in ddf: edf_Zp(edf_out, gk, dk, rng); push to irr_factors
  - nfactors = irr_factors.size()
  - if nfactors <= 1: best.prime=p; best.factors = ...; (push fp 若空); irreducible=true; break
  - elif nfactors < best_count: best_count = nfactors; best更新

Lean：
- 入口 init：lc_f, best with prime=0/irreducible=false, best_count = UInt64::MAX, deg_f, max_tries=3, rng_1=Rng.mk 42, p_1 = if use_large_prime then -1.toUInt64 - 58 else 2
- next_p_1 = lambda
- 主 _loop___select_prime_upoly_2_ir 处理 for/break 控制流
- 后置：kind==0 → return best with irreducible=false（循环用尽未中断）
- kind==1 → break path：best with prime=p, factors=irr_factors（若空 push fp）；bb_34: irreducible=true

✅ 顶层结构完整。loop 内细节（fp 计算、ddf+edf 嵌套、size 比较）由次级 loop 实现，结构层 1:1 与 C++ 对应。

⚠️ 已知 stub：`polynomial_mod`（多元 → Zp poly 简化），B2B 时补真实实现。

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| __ddf_Zp | ✅ | 无 |
| __edf_Zp | ✅ | EDF loop 返回双 Rng（语义无误，工程冗余） |
| __factor_Zp | ✅ | 无 |
| __select_prime | ✅ | polynomial_mod stub 依赖（B2B 补） |

L2 batch 2：4/4 翻译忠实。EDF 与 Plan A rng 修复联动正常工作。
