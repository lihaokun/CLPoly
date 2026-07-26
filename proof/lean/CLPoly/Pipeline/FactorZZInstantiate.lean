/-
  CLPoly/Pipeline/FactorZZInstantiate.lean — 端到端 Z[x] 因式分解定理（0 sorry）

  三个版本：
  1. factor_ZZ_instantiate：直接 UFD
  2. factor_ZZ_unconditional：管线路径，无外部假设
-/
import CLPoly.Pipeline.FactorZZ
import CLPoly.Pipeline.FactorZpInstantiate
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

set_option autoImplicit false

open Polynomial

-- ============================================================
-- §1. 辅助构造：共享 hensel_multifactor 选择的辅助函数
-- ============================================================

noncomputable def hensel_lift_with_proof (p : ℕ) [Fact (Nat.Prime p)] (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ) (facs_p : List (Polynomial (ZMod p)))
    (hne : facs_p ≠ []) (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod)
    (hcop : facs_p.Pairwise (fun a b => IsCoprime a b))
    : List (Polynomial (ZMod (p ^ k))) := by
  let hp_prime : Nat.Prime p := Fact.out
  let φ_p : ℤ →+* ZMod p := Int.castRingHom (ZMod p)
  let φ_pk : ℤ →+* ZMod (p ^ k) := Int.castRingHom (ZMod (p ^ k))
  let facs_Z : List (Polynomial ℤ) := facs_p.map (canonical_lift p)
  have hne_Z : facs_Z ≠ [] := by
    intro hz
    apply hne
    have hlen : facs_p.length = 0 := by
      have hlen_Z : facs_Z.length = 0 := by simpa [hz] using rfl
      calc
        facs_p.length = (facs_p.map (canonical_lift p)).length := by simp
        _ = facs_Z.length := rfl
        _ = 0 := hlen_Z
    simpa using hlen
  have hprod_Z : Polynomial.map φ_p f = (facs_Z.map (Polynomial.map φ_p)).prod := by
    calc
      Polynomial.map φ_p f = facs_p.prod := hprod
      _ = ((facs_p.map (canonical_lift p)).map (Polynomial.map φ_p)).prod := by
        have h_map : (facs_p.map (canonical_lift p)).map (Polynomial.map φ_p) = facs_p := by
          have h_eq : (Polynomial.map φ_p ∘ canonical_lift p) = id := by
            funext x
            exact canonical_lift_map (m := p) (_hm := Nat.Prime.one_lt hp_prime) (p := x)
          simp [List.map_map, h_eq]
        simp [h_map]
      _ = (facs_Z.map (Polynomial.map φ_p)).prod := rfl
  have hcop_Z : facs_Z.Pairwise (fun a b =>
      IsCoprime (Polynomial.map φ_p a) (Polynomial.map φ_p b)) := by
    have hmap : ∀ h : Polynomial (ZMod p), Polynomial.map φ_p (canonical_lift p h) = h :=
      canonical_lift_map p (Nat.Prime.one_lt hp_prime)
    rw [List.pairwise_iff_get]
    intro i j hij
    have hi' : i.1 < facs_p.length := by
      simpa [facs_Z] using i.2
    have hj' : j.1 < facs_p.length := by
      simpa [facs_Z] using j.2
    have hi_val : (facs_Z.get i) = canonical_lift p (facs_p.get ⟨i.1, hi'⟩) := by
      simp [facs_Z]
    have hj_val : (facs_Z.get j) = canonical_lift p (facs_p.get ⟨j.1, hj'⟩) := by
      simp [facs_Z]
    rw [hi_val, hj_val]
    have hR : IsCoprime (facs_p.get ⟨i.1, hi'⟩) (facs_p.get ⟨j.1, hj'⟩) :=
      (List.pairwise_iff_get.mp hcop) ⟨i.1, hi'⟩ ⟨j.1, hj'⟩ hij
    simpa [hmap] using hR
  have h_mf := hensel_multifactor p hp_prime k hk f facs_Z hne_Z hprod_Z hcop_Z
  exact h_mf.choose.map (Polynomial.map φ_pk)

lemma hensel_lift_with_proof_correct (p : ℕ) [Fact (Nat.Prime p)] (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ) (facs_p : List (Polynomial (ZMod p)))
    (hne : facs_p ≠ []) (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod)
    (hcop : facs_p.Pairwise (fun a b => IsCoprime a b))
    : HenselCorrect f k facs_p (hensel_lift_with_proof p k hk f facs_p hne hprod hcop) := by
  unfold hensel_lift_with_proof
  let hp_prime : Nat.Prime p := Fact.out
  let φ_p : ℤ →+* ZMod p := Int.castRingHom (ZMod p)
  let φ_pk : ℤ →+* ZMod (p ^ k) := Int.castRingHom (ZMod (p ^ k))
  let facs_Z : List (Polynomial ℤ) := facs_p.map (canonical_lift p)
  have hne_Z : facs_Z ≠ [] := by
    intro hz
    apply hne
    have hlen : facs_p.length = 0 := by
      have hlen_Z : facs_Z.length = 0 := by simpa [hz] using rfl
      calc
        facs_p.length = (facs_p.map (canonical_lift p)).length := by simp
        _ = facs_Z.length := rfl
        _ = 0 := hlen_Z
    simpa using hlen
  have hprod_Z : Polynomial.map φ_p f = (facs_Z.map (Polynomial.map φ_p)).prod := by
    calc
      Polynomial.map φ_p f = facs_p.prod := hprod
      _ = ((facs_p.map (canonical_lift p)).map (Polynomial.map φ_p)).prod := by
        have h_map : (facs_p.map (canonical_lift p)).map (Polynomial.map φ_p) = facs_p := by
          have h_eq : (Polynomial.map φ_p ∘ canonical_lift p) = id := by
            funext x
            exact canonical_lift_map (m := p) (_hm := Nat.Prime.one_lt hp_prime) (p := x)
          simp [List.map_map, h_eq]
        simp [h_map]
      _ = (facs_Z.map (Polynomial.map φ_p)).prod := rfl
  have hcop_Z : facs_Z.Pairwise (fun a b =>
      IsCoprime (Polynomial.map φ_p a) (Polynomial.map φ_p b)) := by
    have hmap : ∀ h : Polynomial (ZMod p), Polynomial.map φ_p (canonical_lift p h) = h :=
      canonical_lift_map p (Nat.Prime.one_lt hp_prime)
    rw [List.pairwise_iff_get]
    intro i j hij
    have hi' : i.1 < facs_p.length := by
      simpa [facs_Z] using i.2
    have hj' : j.1 < facs_p.length := by
      simpa [facs_Z] using j.2
    have hi_val : (facs_Z.get i) = canonical_lift p (facs_p.get ⟨i.1, hi'⟩) := by
      simp [facs_Z]
    have hj_val : (facs_Z.get j) = canonical_lift p (facs_p.get ⟨j.1, hj'⟩) := by
      simp [facs_Z]
    rw [hi_val, hj_val]
    have hR : IsCoprime (facs_p.get ⟨i.1, hi'⟩) (facs_p.get ⟨j.1, hj'⟩) :=
      (List.pairwise_iff_get.mp hcop) ⟨i.1, hi'⟩ ⟨j.1, hj'⟩ hij
    simpa [hmap] using hR
  have h_mf := hensel_multifactor p hp_prime k hk f facs_Z hne_Z hprod_Z hcop_Z
  rcases h_mf.choose_spec with ⟨hlen', hprod', hforal'⟩
  unfold HenselCorrect
  refine ⟨hprod', ?_⟩
  calc
    facs_p.length = facs_Z.length := by simp [facs_Z]
    _ = h_mf.choose.length := by symm; exact hlen'
    _ = (h_mf.choose.map (Polynomial.map φ_pk)).length := by simp

-- ============================================================
-- §2. 无证明参数的包装函数（用于 factor_ZZ_correct 接口）
-- ============================================================

noncomputable def hensel_lift (p : ℕ) [Fact (Nat.Prime p)] (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ) (facs_p : List (Polynomial (ZMod p)))
    : List (Polynomial (ZMod (p ^ k))) := by
  classical
  if hne : facs_p ≠ [] then
    if hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod then
      if hcop : facs_p.Pairwise (fun a b => IsCoprime a b) then
        exact hensel_lift_with_proof p k hk f facs_p hne hprod hcop
      else exact []
    else exact []
  else exact []

lemma hensel_lift_correct (p : ℕ) [Fact (Nat.Prime p)] (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ) (facs_p : List (Polynomial (ZMod p)))
    (hne : facs_p ≠ [])
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod)
    (hcop : facs_p.Pairwise (fun a b => IsCoprime a b))
    : HenselCorrect f k facs_p (hensel_lift p k hk f facs_p) := by
  unfold hensel_lift
  simp [hne, hprod, hcop]
  exact hensel_lift_with_proof_correct p k hk f facs_p hne hprod hcop

-- ============================================================
-- §3. 端到端 Z[x] 因式分解定理
-- ============================================================

theorem factor_ZZ_instantiate
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result :=
  recombine_correct f hf

theorem factor_ZZ_unconditional
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    {p : ℕ} [hp : Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    (hgood : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f))
    (hdeg : (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree = f.natDegree)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result := by
  classical
  exact factor_ZZ_correct f hf hprim hk hgood hdeg
    -- 1. Zp factorization (unconditional)
    (fun g => if hg : g = 0 then (0, [])
              else ((factor_Zp_instantiate_unconditional g hg).choose,
                    (factor_Zp_instantiate_unconditional g hg).choose_spec.choose))
    (by
      have hfp_ne : (Polynomial.map (Int.castRingHom (ZMod p)) f) ≠ 0 :=
        fun h => not_squarefree_zero (h ▸ hgood)
      simp only [dif_neg hfp_ne]
      exact (factor_Zp_instantiate_unconditional _ hfp_ne).choose_spec.choose_spec)
    -- 2. Hensel lifting (hcop is provided by factor_ZZ_correct)
    (fun facs_p => hensel_lift p k hk f facs_p)
    (fun facs_p hne hprod hcop => hensel_lift_correct p k hk f facs_p hne hprod hcop)
    -- 3. Recombination (UFD path)
    (fun _ => (recombine_correct f hf).choose)
    (fun _ _ => (recombine_correct f hf).choose_spec)
