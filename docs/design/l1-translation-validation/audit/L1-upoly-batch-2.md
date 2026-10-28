# L1 — upoly_* 批 2（11 个函数）

**审视日期**：2026-05-04
**模块**：单变量多项式辅助操作
**状态**：11/11 ✅

---

## 1. `__upoly_divmod` ✅

**C++**：`polynomial_factorize_zp.hh:49-56`（5-arg pair_vec_div 形式）
**Lean**：`Corpus.lean:4772-4776`

```cpp
void __upoly_divmod(q&, r&, f, g) {
    pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
}
```
```lean
__upoly_divmod_ir (q r f g : SparsePolyZp) : (SparsePolyZp × SparsePolyZp) :=
  let __refret_0_1 := pair_vec_div5 q r f g (SparsePolyZp.comp f)
  let q_1 := __refret_0_1.fst; let r_1 := __refret_0_1.snd
  (q_1, r_1)
```

✅ 调用 5-arg pair_vec_div5；q,r ref-out；返回 tuple。
⚠️ stub 依赖：`pair_vec_div5` 在 Lean Model 仍是 stub（与 `__upoly_mod` 同）。B2B 阶段补真实实现。

---

## 2. `__upoly_powmod` ✅

**C++**：`polynomial_factorize_zp.hh:59-85`
**Lean**：`Corpus.lean:5006-5016`（含 `_loop___upoly_powmod_0_ir:4983-5004`）

C++ 算法：标准 square-and-multiply（base^exp mod modpoly），while(e>0) 循环里：
- e 奇 → result *= b; result mod modpoly
- e /= 2
- e>0 → b *= b; b mod modpoly

Lean 翻译：
- 入口：assert; p_1 ← front.snd.prime; one_1 ← __make_zp 1 p_1; result_1 ← [(0, one_1)]; b_1 ← upoly_mod base modpoly; e_1 ← exp; 调 _loop_..._0_ir → result_2.
- Loop body（_loop_..._0_ir）：
  - if e_2 > 0：奇/偶分支 → bb_6（计算 e_3 = e/2; if e_3>0 → b 平方+mod; tail recurse）
  - else (0, result_2)

✅ 控制流 1:1（while → 递归；if 奇/偶 → if-else；嵌套 if e>0 → bb_6 + 内部 if）

---

## 3. `__upoly_const_term` ✅

**C++**：`polynomial_factorize_univar.hh:743-748`
**Lean**：`Corpus.lean:4761-4770`

```cpp
if (f.empty()) return ZZ(0);
if (f.back().first.deg() == 0) return f.back().second;
return ZZ(0);
```
```lean
if isEmpty f then 0
else if back!.fst.deg == 0 then back!.snd else 0
```
✅

---

## 4. `__upoly_norm_l1` ✅

**C++**：`polynomial_factorize_univar.hh:701-711`（绝对值之和）
**Lean**：`Corpus.lean:4958-4964`（含 loop:4942-4956）

C++：`for term: a = term.snd; if a < 0: a = -a; s += a`
Lean loop body：term.snd → a_1; if a_1 < 0 → a_2 = -a_1 else a_1; s_3 = s_2 + a。
✅

---

## 5. `__upoly_norm_l2_sq` ✅

**C++**：`polynomial_factorize_univar.hh:126-132`（系数平方和）
**Lean**：`Corpus.lean:4975-4981`（含 loop:4966-4973）

```cpp
for term: s += term.second * term.second;
```
Lean: `s_3 := s_2 + (term_1.snd * term_1.snd)`。
✅

---

## 6. `__upoly_primitive` ✅

**C++**：`polynomial_factorize_univar.hh:714-722`
**Lean**：`Corpus.lean:5029-5046`（含 loop:5018-5027）

```cpp
if (f.empty()) return {ZZ(1), f};
ZZ c = cont(f);
if (f.front().second < 0) c = -c;
for (auto& term : f) term.second /= c;
return {c, f};
```

Lean：
- `if isEmpty f then (1, id f)` else
- `c_1 := cont f`
- `if front!.snd < 0 then c_2 := -c_1; bb_7 f c_2`
- `else bb_7 f c_1`
- `bb_7`：调 _loop_..._0_ir 把每个 term.snd 除以 c_3，返回 (c_3, f_1)

✅ 一致。

---

## 7. `__upoly_to_poly` ✅

**C++**：`polynomial_factorize_univar.hh:1366-1374`
**Lean**：`Corpus.lean:5181-5184`

```cpp
polynomial_<ZZ,lex_<var_order>> result(comp_ptr);
poly_convert(up, result, var);
return result;
```
```lean
let result_1 := MvPolyZZ.mk comp_ptr
let result_2 := poly_convert3 up result_1 var
result_2
```
✅ 构造 result(comp_ptr) → MvPolyZZ.mk；poly_convert(3-arg) → poly_convert3（rename for arity disambiguation）。
⚠️ stub 依赖：`poly_convert3` 在 Lean Model 是 stub（B2B 时补真实实现）。

---

## 8. `__upoly_mod_coeff` ✅

**C++**：`polynomial_factorize_univar.hh:191-206`（in-place erase-remove pattern）
**Lean**：`Corpus.lean:4933-4935`（含 lambda filt:4928-4931）

C++：双指针 it/out，逐项 fdiv_r 后非零则保留，最后 erase 尾部。
Lean：`Array.filterMap'` + lambda：
```lean
fun __x => 
  let __m_second_0_1 := ZZ.fdiv_r __x.snd __x.snd m
  let __x_mut_1 := (__x.fst, __m_second_0_1)
  if __m_second_0_1 != 0 then Some __x_mut_1 else None
```

✅ 语义等价（filterMap 实现 erase-remove；非零保留 + 系数更新）。

---

## 9. `__upoly_divmod_mod` ✅

**C++**：`polynomial_factorize_univar.hh:209-304`（长除法 mod m，最复杂）
**Lean**：`Corpus.lean:4880-4893`（含 outer loop _1_ir:4857-4878；inner loop _0_ir:4778-4855）

外层 while（_1_ir）：r 非空 && deg(r) ≥ deg_g 时
- d = deg(r) - deg_g
- coeff = front!.snd * lc_inv（fdiv_r mod m）
- if !coeff: erase 首项 → tail recurse
- else: q.push((d, coeff)); new_r empty; iter_idx +=1（skip leading）；调 inner loop；r ← new_r；tail recurse

内层 while（_0_ir）5 分支：
1. g 用尽：c = r_it.snd, fdiv_r → push 保留 r_it
2. r 用尽：deg_term = g_it.deg + d；c = m - ((coeff*g_it.snd)%m)%m；推进 g_it
3. deg_r > deg_term：保留 r_it
4. deg_r < deg_term：从 g 减
5. 相等：c = r_it.snd - coeff*g_it.snd

Lean 5 个分支结构匹配。✅ 控制流 1:1；非线性 `(m - (x%m)%m)` 算法保留。

---

## 10. `__upoly_mul_mod` ✅

**C++**：`polynomial_factorize_univar.hh:323-331`
**Lean**：`Corpus.lean:4937-4940`

```cpp
upolynomial_<ZZ> result = a * b;
__upoly_mod_coeff(result, m);
return result;
```
```lean
let result_1 := a * b
let result_2 := __upoly_mod_coeff_upoly_ir result_1 m
result_2
```
✅

---

## 11. `__upoly_symmetric_mod` ✅

**C++**：`polynomial_factorize_univar.hh:106-119`
**Lean**：`Corpus.lean:5173-5179`（含 loop:5158-5171）

```cpp
for (auto& term : f) {
    ZZ c = __symmetric_mod(term.second, m);
    if (c) result.push_back({term.first, c});
}
```

Lean loop body：c_1 := __symmetric_mod_ir term.snd m；if c_1 ≠ 0 → push (term.fst, c_1) else 跳过。
✅

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| __upoly_divmod | ✅ | pair_vec_div5 stub（B2B 补） |
| __upoly_powmod | ✅ | 无 |
| __upoly_const_term | ✅ | 无 |
| __upoly_norm_l1 | ✅ | 无 |
| __upoly_norm_l2_sq | ✅ | 无 |
| __upoly_primitive | ✅ | 无 |
| __upoly_to_poly | ✅ | poly_convert3 stub（B2B 补） |
| __upoly_mod_coeff | ✅ | filterMap 等价于 erase-remove |
| __upoly_divmod_mod | ✅ | 无 |
| __upoly_mul_mod | ✅ | 无 |
| __upoly_symmetric_mod | ✅ | 无 |

L1 整体（含 batch-1 的 4 个 + random 1 个 + 本批 11 个）共 **16/16 ✅**，无翻译器 bug。仅 2 个 stub 依赖（pair_vec_div5 / poly_convert3），由 B2B 阶段补 Lean Model 真实实现。
