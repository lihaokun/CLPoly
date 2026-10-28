# divmodAux/gcdAux 良基递归骨架 + divmodAux_invariant 归纳结构

## 日期
2026-05-31

## 做了什么
1. 解决了 `gcdAux_unfold` 和 `divmod` 相关 `admit` 的根本障碍：
   - **`partial def` 完全不透明**：`unfold`/`delta`/`simp` 均无法展开，方程引理也不可用
   - 对策：新增 `divmodAux'_wf`、`divmod'_wf`、`gcdAux'_wf` 三个良基递归（`noncomputable def` + `termination_by` + `decreasing_by { admit }`）版本，替代 `partial def` 用于证明
2. 良基递归版本的优点：生成 `.eq_1` 方程引理，可 `rw` 展开
3. 实现了 `divmodAux_invariant` 的强归纳结构：
   - 外层用 `Nat.strong_induction_on` 于 `r[0]!.fst.deg`
   - 内层用 `divmodAux'_wf.eq_1` 展开 + `simp` 化简
   - 递归步骤通过 `h_step` 展开，IH 直接适用
   - 三个子证明（度递减、代数不变量保持、q'/r' 的 AllReduced 保持）用 `admit` 暂挂
4. `divmod_identity` 已接入 `divmodAux_invariant`，编译通过

## 为什么做
`partial def` 不可推理是 `gcdAux_unfold` 和所有 `divmod` 相关证明的根本阻塞。良基递归版本的 `.eq_1` 方程引理提供了推理入口。

## 关键决策
- 不替换 Model.lean 里的 `partial def`（可能影响 C++ 翻译和运行时），而是在 SquarefreeZp.lean 里增加平行良基版本
- `divmodAux'_wf` 的终止性 measure 是 `r[0]!.fst.deg + 1`；`gcdAux'_wf` 的 measure 是 `(toPoly p g).natDegree`
- `section variable p` 在 `termination_by` 中不可见，所以 `gcdAux'_wf` 显式取 `p : ℕ` 为参数
- `let` 绑定在展开后需要用 `dsimp` 化简

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`：新增 `divmod'_wf`、`gcdAux'_wf`；重构 `divmodAux_invariant` 为强归纳结构；`divmod_identity` 保持但不依赖 `divmodAux'_wf_eq`

## 当前状态
- 构建零错误
- `divmodAux_invariant` 有 3 个 `admit` 子证明
- `divmod_deg_lt`、`divmod_snd_wellFormed`、`divmod_snd_allReduced`、`gcdAux_unfold` 保持 `admit`
- 这些 `admit` 不再有结构性阻塞——只需要填充正确的数学推理
