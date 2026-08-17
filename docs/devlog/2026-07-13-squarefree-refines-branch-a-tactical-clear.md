# 2026-07-13 __squarefree_Zp_ir_refines：Branch A 甲类 admit 清零

## 做了什么

在 `proof/lean/CLPoly/Refinement/SquarefreeZp.lean` 中填掉了 L1→L2 精化主定理
`__squarefree_Zp_ir_refines`（导数=0 的 Branch A）及其依赖引理里的 4 处纯策略型
admit/sorry，文件内自有 admit/sorry 从 12 降到 8。每步均以
`lake build CLPoly.Refinement.SquarefreeZp` 验证（全程 1924 jobs 通过）。

具体：

1. **`Zp.inv` 的值域界（原 2013 行）** — 新增独立引理
   `modInv_val_lt (a q) (1 < q.toNat) (q.toNat ≤ UInt64.size) : (Zp.modInv a q).toNat < q.toNat`。
   关键洞察：该界**不依赖** Bezout / 逆元正确性，仅由 `modInv` 的构造得出
   （`a=0` 返回 0；否则返回 `s % q`，`Int.emod` 值域 `[0,q)`，经 toUInt64 往返保值）。
   在 `h_red_inv` 处 `calc` 调用它，配合已有 `ha_prime` 收尾。

2. **`extract_pth_root_toPoly_eq` 空数组分支（原 1755 行）** — `g=#[]` 退化情形：
   LHS 用 `loop_extract_toList` 得 toList 为空 ⇒ `toPoly = 0`；RHS `contract p 0 = 0`
   由 `coeff_contract` 逐系数得出。两侧归零后 `rw`。

3. **`h_loop_lemma` 的数组逐元素等式（原 2130 行）** — 直接用 Lean core 的
   `Array.ext' : a.toList = b.toList → a = b`，把整段 `apply Array.ext`（size + 逐元素
   admit）替换为一行 `exact Array.ext' h_toList`。

4. **`sqfZp_smul` 导数≠0 分支（原 1859 行）** — 本批最难、实为一条完整数学引理：
   `c ≠ 0 → sqfZp (C c * f) = sqfZp f`。C c 在域 `ZMod p` 上是单位，故 Yun 分支内的
   `c = normalize(gcd f f')` 与 `w = normalize(f /ₘ c)` 在单位乘法下不变 ⇒ yunLoop
   输入相同 ⇒ 整个结果相同。
   - gcd 不变：`associated_of_dvd_dvd` + `EuclideanDomain.dvd_gcd/gcd_dvd_*` +
     `associated_unit_mul_left`，再 `normalize_eq_normalize_iff_associated`。
   - w 不变：`div_modByMonic_unique` 证 `(C c*g)/ₘcc = C c*(g/ₘcc)`（cc monic 由
     `monic_normalize`），再由单位相伴 normalize 相等。
   - 收尾：`unfold sqfZp` + 4× `dif_neg` + `dsimp only` 展开 else-else 分支后，
     用 **`simp only [hcc_eq, hw_eq]`**（而非 `rw`）对齐——因为 CC1→CC2 出现在
     `yunLoop` 的依赖证明参数位置，`rw` 会报 "motive is not type correct"，simp 的
     proof-irrelevant congruence 才能处理。

## 为什么做

`__squarefree_Zp_ir_refines` 是顶点定理 `factor_ZZ_cpp_correct`（C++ 翻译代码的 Z[x]
因式分解正确性）SQF 这条腿的核心精化。L3/L2 已 0 sorry，L1 精化层是唯一前沿，其中
SQF 精化最深、最接近完成。先清导数=0 分支中不需要新前置条件的纯策略 admit，使
Branch A 逼近归零。

## 关键决策

- **`modInv_val_lt` 抽为独立轻量引理**而非复用重型的 `Zp_toZMod_inv_mul_self`：值域界
  只需构造性论证，不必牵扯 Bezout，代码更短、更可复用（Branch B 亦可能用到）。
- **`sqfZp_smul` 用 `simp only` 收尾**：`rw` 在依赖证明位破坏 motive，这是 Lean 处理
  「值出现在依赖类型/证明参数中」时的通用陷阱，记录以备后续（Branch B 的 yunLoop
  精化会反复遇到）。

## 遇到的问题与解决

- `modByMonic_add_div` 实际签名为 `(p q : R[X])` **无 Monic 前提**（q 非 monic 时
  divByMonic 返 0、modByMonic 返 p，等式仍成立），最初误传 `hcc_monic` 报类型不符，
  改为传第二个多项式实参。
- `hgcd_g_ne` 里 `h ▸ gcd_dvd_left` 产生 metavariable，改为显式 `have hdvd; rw [h] at hdvd`。

## 剩余 admit/sorry（本文件 8 处）

- **甲/乙分界**：Branch A 纯策略部分已清零。
- **乙类（需前置条件设计，待 §5.1 文档 + 用户决策）**：`h_toPolyList_specific` 的 UInt64
  溢出界（原 2242，现 2354），与 `L1.lean` 的 `h_no_overflow`/`h_deg_bound` 同源。
- **Branch B（原 2259，现 2371）**：导数≠0 → gcd/Yun 全链精化，最大剩余块。
- **GCD 链 / 终止性 / 未使用**：`get_deg_toPoly`(54)、`divmod'_wf_deg_lt`(1083)、
  `divmod_deg_decrease`(445)、`gcdAux'_wf` 终止(667/673/692)——独立模块，Branch B 会用到。

## 涉及的文件

- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`（+`modInv_val_lt` 引理；填 4 处 admit/sorry）

## 度量

- 耗时：~1.5 小时（含总目标梳理 + 调研 Mathlib API + 形式化 + 调试）
- 迭代：`modInv_val_lt`/空数组/Array.ext' 各 1 轮通过；`sqfZp_smul` 2 轮
  （首轮 3 处错误：`modByMonic_add_div` 签名、`▸` metavariable、`rw` motive）
- Lean 新增/修改行数：约 +90 行（含 `modInv_val_lt` 引理 ~22 行、`sqfZp_smul` 分支 ~55 行）
- 对应 C++ 行数：SQF `__squarefree_Zp_ir` 导数=0 分支（p 次根提取 + make_monic）
- 放弃的方案：
  - `sqfZp_smul` 结尾曾用 `rw [hcc_eq, hw_eq]`——因 yunLoop 依赖证明参数导致
    motive 不良而放弃，改 `simp only`。
  - `modInv_val_lt` 曾考虑复用 `Zp_toZMod_inv_mul_self` 的 Bezout 机制——过重，放弃。
