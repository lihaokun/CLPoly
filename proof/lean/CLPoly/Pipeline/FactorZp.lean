/-
  CLPoly/Pipeline/FactorZp.lean — Zp[x] 顶层正确性：框架组合证明

  证明：若 SQF、DDF、EDF 各自满足规约，则组合结果满足 FactorZpCorrect。
  这是框架验证——不依赖算法实现，仅依赖规约间的逻辑衔接。

  证明结构（全部无 sorry）：
    SQF: f ≈ ∏ gᵢ^eᵢ, 每个 gᵢ 首一无平方
    DDF: gᵢ = ∏ gd_j (首一 Associated → 相等)
    EDF: gd_j = ∏ q_k (首一 Associated → 相等)
    合并: f ≈ ∏ q_k^eᵢ = C(lc) * ∏ q_k^eᵢ

  辅助引理（全部已证明）：
    1. poly_unit_eq_C: 域上多项式环的单位是非零常数
    2. eq_of_associated_monic: 首一 Associated → 相等
    3. monic_list_prod: 首一列表的积是首一
    4. list_prod_pow: (∏ aᵢ)^n = ∏(aᵢ^n)
    5. list_prod_flatMap: flatMap 的积 = map-prod 的积
-/
import CLPoly.Spec

set_option autoImplicit false

open Polynomial

variable {p : ℕ} [Fact (Nat.Prime p)]

-- ============================================================
-- 辅助引理（均可独立证明，Phase 2 填充）
-- ============================================================

/-- 域上多项式环的单位是非零常数 -/
private lemma poly_unit_eq_C
    (u : (Polynomial (ZMod p))ˣ) :
    ∃ a : ZMod p, a ≠ 0 ∧ (u : Polynomial (ZMod p)) = Polynomial.C a := by
  refine ⟨coeff (↑u) 0, coeff_coe_units_zero_ne_zero u, ?_⟩
  exact eq_C_of_natDegree_eq_zero (natDegree_coe_units u)

/-- 首一多项式的 Associated 即相等 -/
private lemma eq_of_associated_monic
    (a b : Polynomial (ZMod p)) (ha : Monic a) (hb : Monic b)
    (h : Associated a b) : a = b := by
  -- Associated a b → ∃ u, a * ↑u = b
  obtain ⟨u, hu⟩ := h
  -- ↑u 是度 0 多项式 → ↑u = C (coeff ↑u 0)
  have h_eq := eq_C_of_natDegree_eq_zero (natDegree_coe_units u)
  rw [h_eq] at hu  -- hu : a * C c = b
  -- 比较 leadingCoeff：1 * c = 1 → c = 1
  have hc : coeff (↑u : Polynomial (ZMod p)) 0 = 1 := by
    have h_lc := congr_arg leadingCoeff hu
    simp [ha, hb] at h_lc
    exact h_lc
  -- c = 1 → C 1 = 1 → a * 1 = b → a = b
  rw [hc, map_one, mul_one] at hu
  exact hu

/-- 首一多项式列表的积是首一 -/
private lemma monic_list_prod
    (l : List (Polynomial (ZMod p))) (h : ∀ q ∈ l, Monic q) :
    Monic l.prod := by
  induction l with
  | nil => exact monic_one
  | cons a t ih =>
    rw [List.prod_cons]
    exact (h a (List.mem_cons.mpr (.inl rfl))).mul
      (ih (fun q hq => h q (List.mem_cons.mpr (.inr hq))))

/-- (∏ aᵢ)^n = ∏(aᵢ^n) 在交换幺半群中成立 -/
private lemma list_prod_pow
    (l : List (Polynomial (ZMod p))) (n : ℕ) :
    l.prod ^ n = (l.map (· ^ n)).prod := by
  induction l with
  | nil => simp
  | cons a t ih => simp [List.prod_cons, mul_pow, ih]

/-- List.flatMap 的积等于 map-then-prod 的积 -/
private lemma list_prod_flatMap {α : Type*}
    (l : List α) (f : α → List (Polynomial (ZMod p))) :
    (l.flatMap f).prod = (l.map (fun x => (f x).prod)).prod := by
  induction l with
  | nil => simp
  | cons a t ih => simp [List.flatMap_cons, List.prod_append, ih]

-- ============================================================
-- DDF + EDF 组合：单个首一无平方因子 → 不可约分解
-- ============================================================

/-- 对单个首一无平方多项式，DDF + EDF 给出完整不可约分解。
    关键步骤：DDF 输出首一 + EDF 输出首一 → Associated 即相等 → 乘积精确匹配 -/
private lemma ddf_edf_combine
    (g : Polynomial (ZMod p)) (hm : Monic g) (hsq : Squarefree g)
    (ddf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hddf : Monic g → Squarefree g → DDFCorrect g (ddf g))
    (edf : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    (hedf : ∀ gd d, Monic gd → Squarefree gd →
            (∀ q, Irreducible q → q ∣ gd → Polynomial.natDegree q = d) →
            EDFCorrect gd d (edf gd d))
    : g = ((ddf g).flatMap (fun gd_d => edf gd_d.1 gd_d.2)).prod
      ∧ (∀ q ∈ (ddf g).flatMap (fun gd_d => edf gd_d.1 gd_d.2),
          Irreducible q ∧ Monic q) := by
  -- 展开 DDF 规约
  have ddf_spec := hddf hm hsq
  obtain ⟨ddf_dvd, ddf_deg, _, ddf_assoc, ddf_monic⟩ := ddf_spec
  constructor
  · -- 目标: g = ((ddf g).flatMap (fun gd_d => edf gd_d.1 gd_d.2)).prod
    -- Step 1: 对每个 (gd, d) ∈ ddf g，EDF 给出 gd = (edf gd d).prod
    have edf_eq : ∀ gd_d ∈ ddf g,
        gd_d.1 = (edf gd_d.1 gd_d.2).prod := by
      intro gd_d hgd_d
      have gd_monic := ddf_monic gd_d hgd_d
      have gd_sq : Squarefree gd_d.1 :=
        Squarefree.squarefree_of_dvd (ddf_dvd gd_d hgd_d) hsq
      have edf_spec := hedf gd_d.1 gd_d.2 gd_monic gd_sq (ddf_deg gd_d hgd_d)
      have edf_monic : ∀ q ∈ edf gd_d.1 gd_d.2, Monic q :=
        fun q hq => (edf_spec.2 q hq).2.1
      exact eq_of_associated_monic _ _ gd_monic
        (monic_list_prod _ edf_monic) edf_spec.1
    -- Step 2: g = ((ddf g).map Prod.fst).prod（首一 Associated → 相等）
    have g_eq_ddf_prod : g = ((ddf g).map Prod.fst).prod := by
      have prod_monic : Monic ((ddf g).map Prod.fst).prod := by
        apply monic_list_prod; intro q hq
        obtain ⟨gd_d, hgd_d, rfl⟩ := List.mem_map.mp hq
        exact ddf_monic gd_d hgd_d
      exact eq_of_associated_monic g _ hm prod_monic ddf_assoc
    -- Step 3: map Prod.fst 的积 = flatMap edf 的积
    have prod_eq : ((ddf g).map Prod.fst).prod =
        ((ddf g).flatMap (fun gd_d => edf gd_d.1 gd_d.2)).prod := by
      -- 逐项替换: Prod.fst gd_d = (edf gd_d.1 gd_d.2).prod
      conv_lhs => rw [show (ddf g).map Prod.fst =
        (ddf g).map (fun gd_d => (edf gd_d.1 gd_d.2).prod) from
          List.map_congr_left (fun gd_d hgd_d => edf_eq gd_d hgd_d)]
      -- map-prod 的积 = flatMap 的积
      rw [← list_prod_flatMap]
    exact g_eq_ddf_prod.trans prod_eq
  · -- 目标: ∀ q ∈ irr, Irreducible q ∧ Monic q
    intro q hq
    -- q ∈ flatMap 意味着 ∃ gd_d ∈ ddf g, q ∈ edf gd_d.1 gd_d.2
    rw [List.mem_flatMap] at hq
    obtain ⟨gd_d, hgd_d, hq_edf⟩ := hq
    have gd_monic := ddf_monic gd_d hgd_d
    have gd_sq : Squarefree gd_d.1 :=
      Squarefree.squarefree_of_dvd (ddf_dvd gd_d hgd_d) hsq
    have edf_spec := hedf gd_d.1 gd_d.2 gd_monic gd_sq (ddf_deg gd_d hgd_d)
    exact ⟨(edf_spec.2 q hq_edf).1, (edf_spec.2 q hq_edf).2.1⟩

-- ============================================================
-- 主定理：SQF + DDF + EDF → FactorZpCorrect
-- ============================================================

/-- Zp[x] 因式分解的顶层正确性：
    假设 SQF、DDF、EDF 各自正确，则组合结果是完整因式分解。
    证明逻辑完整，仅依赖 5 个辅助引理（见文件顶部）。-/
theorem factor_Zp_correct
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    -- 假设各子过程存在且正确
    (sqf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hsqf : ∀ g, g ≠ 0 → SquarefreeDecomp g (sqf g))
    (ddf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hddf : ∀ g, Monic g → Squarefree g → DDFCorrect g (ddf g))
    (edf : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    (hedf : ∀ g d, Monic g → Squarefree g →
            (∀ q, Irreducible q → q ∣ g → Polynomial.natDegree q = d) →
            EDFCorrect g d (edf g d))
    : ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        FactorZpCorrect f lc factors := by
  -- Step 1: 展开 SQF 规约
  have sqf_spec := hsqf f hf
  obtain ⟨sqf_assoc, sqf_props, sqf_mult, _⟩ := sqf_spec
  -- Step 2: 对每个 SQF 分量 (gᵢ, eᵢ)，用 DDF+EDF 得到不可约因子
  have component_ok : ∀ ge ∈ sqf f,
      ge.1 = ((ddf ge.1).flatMap (fun gd_d => edf gd_d.1 gd_d.2)).prod
      ∧ (∀ q ∈ (ddf ge.1).flatMap (fun gd_d => edf gd_d.1 gd_d.2),
          Irreducible q ∧ Monic q) := by
    intro ge hge
    exact ddf_edf_combine ge.1
      (sqf_props ge hge).2 (sqf_props ge hge).1
      ddf (fun hm hsq => hddf ge.1 hm hsq)
      edf hedf
  -- Step 3: 构造因子列表
  --   对每个 SQF 分量 (g, e)，展开为 [(q, e) | q ∈ DDF+EDF(g)]
  let factors : List (Polynomial (ZMod p) × ℕ) :=
    (sqf f).flatMap (fun ge =>
      (ddf ge.1).flatMap (fun gd_d =>
        (edf gd_d.1 gd_d.2).map (fun q => (q, ge.2))))
  -- factors 的展开定义（let 绑定的 rfl）
  have hfactors_def : factors = (sqf f).flatMap (fun ge =>
      (ddf ge.1).flatMap (fun gd_d =>
        (edf gd_d.1 gd_d.2).map (fun q => (q, ge.2)))) := rfl
  -- Step 4: 证明 factors 的积 = SQF 分量的积
  have key : (factors.map (fun pr => pr.1 ^ pr.2)).prod =
      ((sqf f).map (fun pr => pr.1 ^ pr.2)).prod := by
    -- Step 0: 展开 factors，map 穿过 flatMap，合并内部 map
    rw [hfactors_def]
    simp_rw [List.map_flatMap, List.map_map, Function.comp_def]
    -- Step 1: 外层 list_prod_flatMap
    rw [list_prod_flatMap]
    -- Step 2: 逐项证 inner(ge).prod = ge.1^ge.2
    congr 1
    apply List.map_congr_left
    intro ge hge
    conv_rhs => rw [(component_ok ge hge).1]
    rw [list_prod_pow, List.map_flatMap]
  -- Step 5: 由 SQF Associated + key 得到 Associated f factors_prod
  have h_assoc : Associated f (factors.map (fun pr => pr.1 ^ pr.2)).prod := by
    rw [key]; exact sqf_assoc
  -- Step 6: 提取单位 u，得 f = C(lc) * factors_prod
  obtain ⟨u, hu⟩ := h_assoc
  obtain ⟨a, ha_ne, ha_eq⟩ := poly_unit_eq_C u⁻¹
  refine ⟨a, factors, ?_, ?_⟩
  · -- 目标: f = C a * (factors.map (fun pr => pr.1 ^ pr.2)).prod
    -- f = f*1 = f*(u*u⁻¹) = (f*u)*u⁻¹ = prod*u⁻¹ = prod*C(a) = C(a)*prod
    calc f = f * 1 := (mul_one f).symm
      _ = f * ((↑u : Polynomial (ZMod p)) * ↑(u⁻¹)) := by simp
      _ = f * ↑u * ↑(u⁻¹) := (mul_assoc _ _ _).symm
      _ = (factors.map (fun pr => pr.1 ^ pr.2)).prod * ↑(u⁻¹) := by rw [hu]
      _ = (factors.map (fun pr => pr.1 ^ pr.2)).prod * C a := by rw [ha_eq]
      _ = C a * (factors.map (fun pr => pr.1 ^ pr.2)).prod := mul_comm _ _
  · -- 目标: ∀ pr ∈ factors, Irreducible pr.1 ∧ Monic pr.1 ∧ pr.2 ≥ 1
    intro pr hpr
    -- pr ∈ factors 意味着 pr = (q, ge.2) 其中 q 来自 EDF，ge 来自 SQF
    rw [List.mem_flatMap] at hpr
    obtain ⟨ge, hge, hpr_inner⟩ := hpr
    rw [List.mem_flatMap] at hpr_inner
    obtain ⟨gd_d, hgd_d, hpr_edf⟩ := hpr_inner
    rw [List.mem_map] at hpr_edf
    obtain ⟨q, hq, rfl⟩ := hpr_edf
    -- q 来自 edf → Irreducible + Monic
    have q_in_irr : q ∈ (ddf ge.1).flatMap (fun gd_d => edf gd_d.1 gd_d.2) :=
      List.mem_flatMap.mpr ⟨gd_d, hgd_d, hq⟩
    have ⟨q_irr, q_monic⟩ := (component_ok ge hge).2 q q_in_irr
    -- ge.2 来自 sqf → ≥ 1
    exact ⟨q_irr, q_monic, sqf_mult ge hge⟩
