# divmod_identity h_prime type mismatch fix + well-founded recursion skeleton

## 日期
2026-05-31

## 做了什么
1. 修复 `divmod_identity` 中 `h_prime` 证明的类型错误（`g[0]! ∈ g` 与 `g[0] ∈ g` 不 definitionally equal）
2. 为长期重构奠定基础：在 `SquarefreeZp.lean` 中添加了三个新结构：
   - `divmodAux'_wf` — 良基递归版 `divmodAux`（noncomputable，`decreasing_by { admit }`）
   - `divmodAux'_wf_eq` — 与 `partial def divmodAux` 等价的引理（admit）
   - `divmodAux_invariant` — divmodAux 循环不变量（admit）

## 为什么做
构建因类型错误而失败；该修复恢复了零错误状态。同时提前设计了良基递归和循环不变量的接口，这是后续消除 `gcdAux`/`divmodAux` 相关 `admit` 的前提。

## 关键决策
- `g[0]! = g[0]` 不是 definitionally equal，需要使用 `simp [hpos]` 证明；`Array.getElem_mem` 返回的是 `g[0] ∈ g`（`getElem` 形式），目标 `g[0]! ∈ g` 需要用 `simpa [h_eq]` 转换
- 新的 `divmodAux'_wf` 接口不需要 `h_2p`/`h_p2`/`h_nonempty_g` 等终止性条件（终止性用 `decreasing_by { admit }` 绕过），但 `divmodAux_invariant` 需要这些条件来确保操作的正确性

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`（新增 `divmodAux'_wf`、`divmodAux'_wf_eq`、`divmodAux_invariant`；修改 `divmod_identity` 的 `h_prime` 证明）
