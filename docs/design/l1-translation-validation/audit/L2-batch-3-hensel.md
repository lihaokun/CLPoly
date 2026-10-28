# L2 — Hensel 提升簇（批 3，6 个函数）

**审视日期**：2026-05-04
**模块**：Hensel 提升（二次 + 线性）+ 树构建 + 因子提取
**状态**：6/6 ✅

---

## 1. `__hensel_step` ✅（核心二次提升）

**C++**：`polynomial_factorize_univar.hh:404-497`
**Lean**：`Corpus.lean:1233-1302`（含 6 个 sub-loops + 2 个 lambda filter）

C++ 算法（Newton 二次 Hensel：m → m²）：
1. 第一部分（提升因子）：
   - e = (f - g*h) / m mod m （精确整除 m，再 mod m）
   - 滤掉零项
   - se = s * e；divmod_mod(q_se, r_se, se, h, m)
   - tau = t*e + q_se*g 全 mod m
   - g_new = g + m*tau mod m²
   - h_new = h + m*r_se mod m²
2. 第二部分（提升 Bezout）：
   - ep = (1 - s*g - t*h) / m mod m
   - 滤掉零项
   - sep = s * ep；divmod_mod(q_sep, r_sep, sep, h, m)
   - s_new = s + m*r_sep mod m²
   - tau2 = t*ep + q_sep*g mod m
   - t_new = t + m*tau2 mod m²

Lean：每个 `for term: fdiv_q + fdiv_r` 循环对应 `_loop___hensel_step_upoly_<k>_ir`（共 6 个：处理 e/tau/r_se/ep/r_sep/tau2 的 ×m 或 fdiv 步）。lambda filter `_lambda___hensel_step_upoly_filt1/2_ir` 处理 erase-remove 模式。divmod_mod 调用、节点 record-update（`{ node with g := ... }`）、最终 mod m² 对应 C++ 全部步骤。
✅

---

## 2. `__hensel_step_linear` ✅（P1b 线性提升）

**C++**：`polynomial_factorize_univar.hh:636-679`
**Lean**：`Corpus.lean:1330-1361`（含 2 sub-loops + 1 lambda filter）

C++ 算法（线性 Hensel：m → m·p，s, t 不变）：
- mp = m * p
- e = (f - g*h) / m mod p（用 fdiv_q + fdiv_r 链）
- 若 e 空 → return（无需修正）
- se = s * e；divmod_mod(q_se, sigma, se, h, p)
- tau = t*e + q_se*g mod p
- g_new = g + m*tau mod mp
- h_new = h + m*sigma mod mp
- s, t 不更新

Lean：filterMap' lambda 实现 fdiv_q+fdiv_r 链 + 非零保留；isEmpty e_2 早返回 node；divmod_mod 处理 q_se/sigma；2 个 sub-loops 把 tau 和 sigma 各项乘 m；record-update 链产生 node_4。
✅

---

## 3. `__hensel_tree_build` ✅

**C++**：`polynomial_factorize_univar.hh:392-401`
**Lean**：`Corpus.lean:1363-1369`

C++：
```cpp
nodes.push_back({});  // root
__hensel_tree_build_recursive(nodes, factors, p, 0, factors.size(), 0);
return nodes;
```

Lean：nodes_1 = []; nodes_2 = push 默认 HenselNode（8 字段全初始化）；调 recursive；返回 nodes_3。
✅ 注意 8 字段 HenselNode 初始化与 C++ struct 默认值对齐（g/h/s/t/left/right/leaf_start/leaf_end）。

---

## 4. `__hensel_tree_build_recursive` ✅

**C++**：`polynomial_factorize_univar.hh:334-389`
**Lean**：`Corpus.lean:1391-...`（含 _0_ir/_1_ir 两个乘积循环）

C++ 算法：
- mid = (start + end) / 2
- g_zp = ∏ factors[start..mid)；h_zp = ∏ factors[mid..end)
- 若 size == 1 → 叶子（设 left=right=-1, leaf_start=start, leaf_end=end）
- 否则递归构造左右子节点；用 Bezout 在 Z_p 上算 s, t
- 父节点 g, h 提升到 ZZ（mod p 表示）

Lean 结构：
- _0_ir：从 i=start 到 mid 累乘 factors → g_zp
- _1_ir：从 i=mid 到 end 累乘 factors → h_zp
- bb_11 处理右子节点递归 + 父节点 right 索引设置
- 主体：左叶子 / 右叶子 / 递归 + Bezout（s, t 计算）

✅ 控制流嵌套与 C++ 对齐。

---

## 5. `__hensel_extract_factors` ✅

**C++**：`polynomial_factorize_univar.hh:500-515`
**Lean**：`Corpus.lean:982-1000`

C++：
```cpp
if (left == -1) factors.push_back(node.g);
else __hensel_extract_factors(nodes, left, factors);
if (right == -1) factors.push_back(node.h);
else __hensel_extract_factors(nodes, right, factors);
```

Lean：
- node_1 := nodes[idx]
- if left == -1 → push node_1.g else 递归 left
- bb_3：if right == -1 → push node_1.h else 递归 right

✅ 1:1 双分支结构。

---

## 6. `__hensel_lift` ✅（多因子 Hensel 提升入口）

**C++**：`polynomial_factorize_univar.hh:537-...`（部分见 Read 输出）
**Lean**：`Corpus.lean:1038-1099`（含 0/1/2/3 sub-loops）

C++ 算法：
1. 确定提升精度 target：a_target=0 → Mignotte（2*lc_f*B），else → p^a_target - 1
2. 构造 factors_adj：factor[0] *= lc_mod_p；normalize
3. nodes = tree_build(factors_adj, p)
4. m = p；while m ≤ target：lift_recursive；m *= m
5. 提取因子；若 result[0] 非首项 1，乘 inverse → mod m_2

Lean：
- if a_target == 0：B_1 = mignotte_bound；lc_f 符号修正 → bb_6 → target_2 = 2*lc_f*B
- else：target_3 = 1；_loop_..._0_ir 累乘 p a_target 次 → target_6 = p^a_target - 1
- bb_3 主体：
  - factors_adj_2 = set 第 0 项乘 lc_mod_p；factors_adj_3 = normalize
  - nodes_1 = tree_build；m_1 = p
  - _loop_..._2_ir：while m ≤ target → lift_recursive + m*=m
  - extract_factors → result_2
  - if result[0] 非首项 1：lc_inv = ZZ.invert；_loop_..._3_ir 第一个因子各项 ×= lc_inv → mod m_2
- bb_20：return (result, m_2)

✅ 所有控制流分支对齐。

---

## 模块小结

| 函数 | 状态 | 偏差 |
|------|------|-----|
| __hensel_step | ✅ | 无 |
| __hensel_step_linear | ✅ | 无 |
| __hensel_tree_build | ✅ | 无 |
| __hensel_tree_build_recursive | ✅ | 无 |
| __hensel_extract_factors | ✅ | 无 |
| __hensel_lift | ✅ | 无 |

L2 batch 3：6/6 翻译忠实。Hensel 算法（二次 + 线性 + 树构建 + 因子提取）全 1:1 与 C++ 对应。无翻译器 bug。
