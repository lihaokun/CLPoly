# Phase 2 前半：T2.1-T2.3 + T2.2' 形式化完成

> 日期：2026-03-17
> 分支：`feature/formal-proofs`
> 提交：`be38bb2`

---

## 做了什么

新建 `FiniteFieldFact.lean`，完成 Phase 2（L3 数学基石）的前 4 个定理：T2.1、T2.2、T2.2'、T2.3，全部 0 sorry。同时为 Phase 1 补写了 nl-proof 草稿。

### 定理清单

| 定理 | 内容 | 行数 | 证明方法 |
|------|------|------|---------|
| T2.1 | `X^{p^d}-X` 的 Separable 性 | ~8 | `derivative` 计算 + char p 下导数 = -1 |
| T2.2 | 不可约 g, deg g \| d → g \| X^{p^d}-X | ~18 | AdjoinRoot + `pow_card_pow` 指数拼接 |
| T2.2' | g \| X^{p^d}-X → deg g \| d（逆命题） | ~55 | AlgHom ext（Frobenius = id）+ 反证法根计数 |
| T2.3 | gcd(X^{p^d}-X, f) 不可约因子刻画 | ~12 | 正向 `dvd_trans`，反向 `dvd_gcd` + T2.2 |

### 同时完成的 nl-proof 草稿

| 文件 | 内容 |
|------|------|
| `phase1-factorzp.md` | Phase 1 FactorZp 组合证明的自然语言草稿 |
| `phase1-factorzz.md` | Phase 1 FactorZZ 组合证明的自然语言草稿 |
| `phase2-t21-separable.md` | T2.1 证明草稿 |
| `phase2-t22-irred-dvd.md` | T2.2 证明草稿 |
| `phase2-t22p-converse.md` | T2.2' 证明草稿 |
| `phase2-t23-gcd-factors.md` | T2.3 证明草稿 |

## 关键技术点

1. **T2.2' 是最难的定理（~55 行）**：核心是 Step B（Frobenius 固定生成元 → 固定整个域）。最终采用 AlgHom ext 方案——构造 `iterateFrobenius` 为 AlgHom，由 `AdjoinRoot.algHom_ext` 一步完成，比备选方案（PowerBasis 展开、Polynomial.expand）简洁得多。

2. **T2.2' Step C 反证法**：假设 k ∤ d，取 r = d mod k，推出 X^{p^r}-X 在 AdjoinRoot g 中有 p^k 个根但度只有 p^r < p^k，用 `card_roots'` 矛盾。

3. **Fintype 实例构造**：AdjoinRoot g 的 Fintype 通过 `PowerBasis.basis.equivFun.toEquiv.symm` 构造，card 通过 `card_fun + card_fin + ZMod.card` 计算。这是所有涉及 AdjoinRoot 的定理的共同 boilerplate。

4. **T2.3 出乎意料地简单**（~12 行）：正向用 `dvd_trans` + `gcd_dvd`，反向用 `dvd_gcd` + T2.2。去掉了原始设计中不必要的 Monic/Squarefree 假设。

## 度量
- 耗时：~10 小时（nl-proof 草稿 ~2h、T2.1 ~0.5h、T2.2 ~1.5h、T2.2' ~4h、T2.3 ~0.5h、nl-proof Phase 1 补写 ~1.5h）
- 迭代：T2.1 1 轮；T2.2 2 轮；T2.2' ~8 轮（AlgHom commutes' 证明 + Fintype 实例 + 反证法根计数的 natDegree 计算）；T2.3 1 轮
- Lean 新增行数：167 行（FiniteFieldFact.lean）
- 对应 C++ 行数：N/A（纯数学定理；为 polynomial_factorize_zp.hh 中 DDF/EDF 的正确性提供理论基础）
- 放弃的方案：(1) T2.2' Step B 最初考虑 PowerBasis 展开法（每个元素写成 ∑aᵢαⁱ）——太繁琐；(2) 考虑 Polynomial.expand 代换法——需要额外引理；(3) 最终采用 AlgHom ext 仅 ~15 行

## 涉及文件

| 文件 | 操作 |
|------|------|
| `proof/lean/CLPoly/Math/FiniteFieldFact.lean` | 新建：T2.1-T2.3 + T2.2' |
| `proof/nl-proof/phase2-t21-separable.md` | 新建 |
| `proof/nl-proof/phase2-t22-irred-dvd.md` | 新建 |
| `proof/nl-proof/phase2-t22p-converse.md` | 新建 |
| `proof/nl-proof/phase2-t23-gcd-factors.md` | 新建 |
| `proof/nl-proof/phase1-factorzp.md` | 新建（补写） |
| `proof/nl-proof/phase1-factorzz.md` | 新建（补写） |
