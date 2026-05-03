# L2 — Recombination + LLL（批 4，6 个函数）

**审视日期**：2026-05-04
**模块**：因子重组 (Zassenhaus / van Hoeij) + LLL 矩阵构建/约化
**状态**：6/6 ✅（结构层）

---

## 1. `__zassenhaus_recombine` ✅

**C++**：`polynomial_factorize_univar.hh:752-882`
**Lean**：`Corpus.lean:6493-6524`（含 `_loop_..._5_ir` 主循环 + 多个 sub-loop + 2 lambda）

C++ 算法（穷举子集试除）：
- 入口：if lifted size ≤ 1：deg(f) > 0 → [f] else []
- 否则：T = [0..r-1]；f_star = f；result = []
- 主 while：`while 2*s ≤ |T|`
  - 内 do-while：枚举 idx[0..s-1] 组合（next_combination lambda）
  - 剪枝 1：lc 检查 (lc_prod = lc_fstar * ∏ lifted[i].lc; lc_sq mod lc_prod == 0?)
  - 剪枝 2：常数项检查
  - 完整试除：g = subset_product；primitive(g)；divmod；r 空 → 真因子，T 删除已用，found=true, break
  - 内 while 穷举：next_combination
  - found → s=1 重启；else ++s
- 尾部：f_star 非空 → push；sort by deg ascending

Lean 结构：
- if size lifted ≤ 1 走快速分支（与 C++ 一致）
- 否则：T_2 = range_init [0..r-1]；调主 _loop_..._5_ir → (kind, f_star_2, result_2)
- bb_11：f_star 非空 + deg > 0 → push；bb_59：sort

主 _5_ir 内部嵌套子 loop 处理 next_combination 与组合枚举，每轮 lc/常数项剪枝 + subset_product + primitive + divmod 的尝试。

✅ 顶层结构 1:1，剪枝顺序与 C++ 一致。组合枚举 lambda 在 Pass 3 lambda lift 后转为外部函数。

---

## 2. `__build_cld_matrix` ✅

**C++**：`polynomial_factorize_univar.hh:937-...`
**Lean**：`Corpus.lean:143-...`

C++ 算法（构造 CLD（系数长除）矩阵的列）：
- 输入：M (LLL 矩阵)、cld 多项式数组、当前 J 宽度、目标 J、padding
- 对每个 cld[i]，从指定 deg 范围抽取系数，写入 M 的列

Lean：使用 lifted lambda `_lambda___build_cld_matrix_upoly_1_ir`（从 cld[i] 找特定 deg term），通过 _loop_..._3_ir 等遍历构建矩阵。
✅ 结构对应。

⚠️ Lean 主体里 `_upoly_coeff_1` (lambda) 标记为 unused（前缀 `_`）— 这是 LLL 构建中已 inline 到 _loop 的辅助函数，属预期。

---

## 3. `__cld_polys` ✅

**C++**：`polynomial_factorize_univar.hh:901-...`（构造 CLD 多项式数组）
**Lean**：`Corpus.lean:178-184`（含 _loop_..._0_ir）

C++：对每个 active_factor[i]：q_i = f_star quo h_i (mod m)；h_prime_i = h_i'（导数）；mod_coeff；C_i = q_i * h_prime mod m；symmetric_mod；push 到 result。

Lean：result = []；range-for over active_factors：每轮 divmod_mod q_i r_i f_star h_i m → q_i_2; h_prime_1 = derivative h_i; h_prime_2 = mod_coeff; C_i_1 = mul_mod q_i_2 h_prime_2 m; C_i_2 = symmetric_mod; push (id C_i_2)。

✅ 顺序与 C++ 完全一致。

---

## 4. `__lll_factorize` ✅

**C++**：`polynomial_factorize_univar.hh:1485-...`（vanhoeij 主入口）
**Lean**：`Corpus.lean:1538-...`

C++ 算法（用 Hensel 提升 + LLL 约化做 vanhoeij recombine）：
- 提升 factors 到 m_h（FLINT 启发式精度）
- 调 vanhoeij_recombine
- 失败回退：再提升到 m_mig（Mignotte 完整精度），再 vanhoeij_recombine

Lean：嵌套调用 hensel_lift + vanhoeij_recombine，失败检测后再 hensel_lift_to_mignotte + vanhoeij_recombine。
✅ 顶层结构对应。

---

## 5. `__lll_reduce` ✅

**C++**：`polynomial_factorize_univar.hh:993-...`（LLL 算法 with size-reduce + swap-step）
**Lean**：`Corpus.lean:1948-...`（含 lifted lambdas + 多个 sub-loops）

C++ 算法（标准 LLL）：
- 输入：基础矩阵 M、变换矩阵 U、目标界 B
- 主循环 k=1..n：size-reduce + swap or step
- 输出：short_rows 索引 + 更新的 M, U

Lean：多层 _loop（k 主循环 + 嵌套 size-reduce）+ lifted lambda（`_lambda___lll_reduce_*_ir`）处理 round-half + swap 等子步骤。lambda 内部的 fdiv_q 等通过 Pass 2/2b ref-elim 重写。
✅ 顶层主循环 + 子步骤分解结构与 C++ 对齐。

⚠️ LLL 算法核心在嵌套 lambda 与 sub-loops 中，详细数值正确性需 B2B 验证。结构层翻译忠实。

---

## 6. `__vanhoeij_recombine` ✅

**C++**：`polynomial_factorize_univar.hh:1189-...`（main vanhoeij algorithm）
**Lean**：`Corpus.lean:5424-...`（最大的单函数，含 J_target=0 pre-pass + batch extraction）

C++ 算法（经过 Maple 优化的 vanhoeij）：
- 主循环：J_cur 增加 → build_cld_matrix → lll_reduce → 检测 short rows
- short rows 解析为子集索引（batch extraction）
- 每个子集试除验证 → 真因子加入 result
- f_star 缩减；continue
- J_target=0 pre-pass：先用最小 J 拿掉 trivial 因子
- 失败回退：fall back to zassenhaus_recombine

Lean 结构：包含主 outer loop + J_target=0 special case + batch extraction loop + zassenhaus fallback。多个 lifted lambdas 处理子集解码 + 排序 + 比较。
✅ 顶层流程与 C++ 对应。

⚠️ 此函数 ~1000 行，详细逐行审视成本高，本轮做结构层验证。`lake build` 0 errors 验证类型一致；语义层留 B2B 测试覆盖。

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| __zassenhaus_recombine | ✅ | 无 |
| __build_cld_matrix | ✅ | 无 |
| __cld_polys | ✅ | 无 |
| __lll_factorize | ✅ | 无 |
| __lll_reduce | ✅ | 结构层；嵌套 lambda 数值正确性留 B2B |
| __vanhoeij_recombine | ✅ | 结构层；细节留 B2B |

L2 batch 4：6/6 顶层结构忠实。LLL 算法核心 + vanhoeij 大循环细节超出快速 1:1 审视范围，标记结构层 ✅，详细数值正确性留 B2B 测试覆盖。
