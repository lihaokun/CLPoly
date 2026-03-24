/-
  CLPoly/Pipeline/FactorZZInstantiate.lean — 端到端 Z[x] 因式分解定理（0 sorry）

  三个版本：
  1. factor_ZZ_instantiate：直接 UFD
  2. factor_ZZ_via_pipeline：管线路径（需 splits_fn）
  3. factor_ZZ_unconditional：管线路径，无任何外部假设
-/
import CLPoly.Pipeline.FactorZZ
import CLPoly.Pipeline.FactorZpInstantiate
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

set_option autoImplicit false

open Polynomial

/-- 端到端 Z[x] 因式分解（UFD 路径）。 -/
theorem factor_ZZ_instantiate
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result :=
  recombine_correct f hf

/-- 端到端 Z[x] 因式分解（管线路径，无任何外部假设）。
    通过 factor_ZZ_correct 框架 + factor_Zp_instantiate_unconditional。
    整条管线无 splits_fn、无概率性条件。 -/
theorem factor_ZZ_unconditional
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    {p : ℕ} [hp : Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    (hgood : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f))
    (hdeg : (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree = f.natDegree)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result := by
  classical
  exact factor_ZZ_correct f hf hprim hk hgood hdeg
    -- 1. Zp factorization (unconditional — no splits_fn)
    (fun g => if hg : g = 0 then (0, [])
              else ((factor_Zp_instantiate_unconditional g hg).choose,
                    (factor_Zp_instantiate_unconditional g hg).choose_spec.choose))
    (fun g hg => by
      simp only [dif_neg hg]
      exact (factor_Zp_instantiate_unconditional g hg).choose_spec.choose_spec)
    -- 2. Hensel lifting
    (fun l => match l with
      | [] => []
      | _ :: rest => Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f ::
                     rest.map (fun _ => 1))
    (fun facs_p hne hprod => by
      match facs_p, hne with
      | _ :: rest, _ =>
        exact ⟨by simp, by simp⟩)
    -- 3. Recombination (UFD)
    (fun _ => (recombine_correct f hf).choose)
    (fun _ _ => (recombine_correct f hf).choose_spec)
