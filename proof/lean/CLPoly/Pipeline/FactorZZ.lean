/-
  CLPoly/Pipeline/FactorZZ.lean — Z[x] 顶层正确性：框架组合证明

  证明：若 Zp 因式分解、Hensel 提升、因子重组各自满足规约，
  则组合结果满足 FactorZZCorrect。

  证明结构（全部无 sorry）：
    1. Squarefree → f mod p ≠ 0
    2. factor_zp: f mod p = C(lc) * ∏(fᵢ^eᵢ)
    3. 构造 Hensel 输入: facs_p = [C(lc), f₁^e₁, ...]，乘积 = f mod p
    4. hensel: map f mod p^k = (hensel facs_p).prod
    5. recombine → RecombineCorrect ≡ FactorZZCorrect（定义相同）

  关键观察：RecombineCorrect 和 FactorZZCorrect 定义完全相同，
  因此组合证明只需串联三个子过程假设，无需额外辅助引理。
-/
import CLPoly.Spec
import Mathlib.RingTheory.Polynomial.Content

set_option autoImplicit false

open Polynomial

/-- Z[x] 因式分解的顶层正确性：
    假设 Zp 因式分解、Hensel 提升、因子重组各自正确，
    则组合结果是完整不可约分解。

    注：hk, hprim, hdeg 是算法前置条件（选素数/确保提升可行），
    框架组合证明不直接使用——它们通过子过程假设间接起作用。-/
theorem factor_ZZ_correct
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    -- 选定的素数和 Hensel 提升指数
    {p : ℕ} [Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    -- p 是合适的素数
    (hgood : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f))
    (hdeg : Polynomial.natDegree (Polynomial.map (Int.castRingHom (ZMod p)) f)
            = Polynomial.natDegree f)
    -- 假设子过程存在且正确
    -- 1. Zp 因式分解
    (factor_zp : Polynomial (ZMod p) → ZMod p × List (Polynomial (ZMod p) × ℕ))
    (hfzp : ∀ g, g ≠ 0 → FactorZpCorrect g (factor_zp g).1 (factor_zp g).2)
    -- 2. Hensel 提升
    (hensel : List (Polynomial (ZMod p)) → List (Polynomial (ZMod (p ^ k))))
    (hhensel : ∀ facs_p,
        Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod →
        HenselCorrect f k facs_p (hensel facs_p))
    -- 3. 因子重组
    (recombine : List (Polynomial (ZMod (p ^ k))) → List (Polynomial ℤ))
    (hrecombine : ∀ facs_pk,
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f = facs_pk.prod →
        RecombineCorrect f (recombine facs_pk))
    : ∃ result : List (Polynomial ℤ),
        FactorZZCorrect f result := by
  -- Step 1: f mod p ≠ 0（Squarefree 蕴含非零）
  set fp := Polynomial.map (Int.castRingHom (ZMod p)) f with hfp_def
  have hfp_ne : fp ≠ 0 := by
    intro h; rw [h] at hgood; exact not_squarefree_zero hgood
  -- Step 2: 对 fp 做 Zp 因式分解
  have fzp_ok := hfzp fp hfp_ne
  -- Step 3: 构造 Hensel 输入列表
  --   facs_p = [C(lc), f₁^e₁, f₂^e₂, ...]
  --   乘积 = C(lc) * ∏(fᵢ^eᵢ) = fp（由 FactorZpCorrect.1 + List.prod_cons）
  have h_facs_prod : fp =
      (Polynomial.C (factor_zp fp).1 ::
       ((factor_zp fp).2.map (fun pr => pr.1 ^ pr.2))).prod := by
    rw [List.prod_cons]; exact fzp_ok.1
  -- Step 4: Hensel 提升 mod p → mod p^k
  have hensel_ok := hhensel _ h_facs_prod
  -- Step 5: 因子重组 → FactorZZCorrect
  --   RecombineCorrect f result ≡ FactorZZCorrect f result（定义完全相同）
  exact ⟨recombine (hensel _), hrecombine _ hensel_ok.1⟩
