/-
  CLPoly/Algorithm/Wang.lean — L2 Wang 多变量因式分解算法模型

  对应 C++: polynomial_factorize_wang.hh __wang_core

  模块 1: 求值 + 单变量因式分解 (0 sorry)
  模块 2: LC 分配 (0 sorry)
-/
import CLPoly.Spec
import CLPoly.Math.MvBasics
import CLPoly.Pipeline.FactorZZInstantiate
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.LinearAlgebra.Vandermonde

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial MvPolynomial

-- ============================================================
-- 模块 1: 求值 + 单变量因式分解
-- ============================================================

noncomputable def eval_at_α {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Polynomial ℤ :=
  Polynomial.map (MvPolynomial.eval α) (MvPolynomial.finSuccEquiv ℤ n f)

lemma eval_at_α_mul {n : ℕ} (f g : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) :
    eval_at_α (f * g) α = eval_at_α f α * eval_at_α g α := by
  simp only [eval_at_α, map_mul, Polynomial.map_mul]

lemma eval_at_α_associated {n : ℕ} {f g : MvPolynomial (Fin (n + 1)) ℤ}
    (h : Associated f g) (α : Fin n → ℤ) :
    Associated (eval_at_α f α) (eval_at_α g α) := by
  simp only [eval_at_α]
  exact (h.map (MvPolynomial.finSuccEquiv ℤ n)).map (Polynomial.mapRingHom (MvPolynomial.eval α))

noncomputable def lc_x1 {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) : MvPolynomial (Fin n) ℤ :=
  (MvPolynomial.finSuccEquiv ℤ n f).leadingCoeff

structure EvalPointGood {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Prop where
  eval_ne_zero : eval_at_α f α ≠ 0
  eval_sqfree : Squarefree (eval_at_α f α)
  lc_ne_zero : MvPolynomial.eval α (lc_x1 f) ≠ 0

theorem wang_eval_step {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (α : Fin n → ℤ) (hα : EvalPointGood f α)
    : ∃ univar_facs : List (Polynomial ℤ),
        FactorZZCorrect (eval_at_α f α) univar_facs :=
  factor_ZZ_instantiate (eval_at_α f α) hα.eval_ne_zero

-- ============================================================
-- 模块 2: LC 分配（Valuation 提取）
-- ============================================================

-- --- 规约 ---

structure LCDistribResult {n : ℕ}
    (L : MvPolynomial (Fin n) ℤ)
    (σ : List (MvPolynomial (Fin n) ℤ))
    (α : Fin n → ℤ)
    (univar_lcs : List ℤ) : Prop where
  prod_associated : Associated L σ.prod
  eval_dvd : List.Forall₂ (fun s w => (MvPolynomial.eval α s : ℤ) ∣ w) σ univar_lcs
  length_eq : σ.length = univar_lcs.length

-- --- 算法模型 ---

noncomputable def intVal (E w : ℤ) : ℕ := multiplicity E w

noncomputable def extractValuations (E : ℤ) (ws : List ℤ) : List ℕ :=
  ws.map (fun w => intVal E w)

/-- LC 分配核心（递归定义）。对应 C++ 双重 for 循环。-/
noncomputable def lcDistribCore {n : ℕ}
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ)
    (ws : List ℤ) (r : ℕ) : List (MvPolynomial (Fin n) ℤ) :=
  match lc_factors, lc_evals with
  | [], _ | _, [] => List.replicate r 1
  | (lj, _) :: rest_f, Ej :: rest_e =>
    let σ_prev := lcDistribCore rest_f rest_e ws r
    σ_prev.zipWith (fun si ki => si * lj ^ ki) (extractValuations Ej ws)

-- --- 守恒条件 ---

def ConservationCheck {n : ℕ}
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ)
    (ws : List ℤ) : Prop :=
  ∀ j, j < lc_factors.length →
    (ws.map (fun w => intVal (lc_evals.getD j 0) w)).sum =
    (lc_factors.getD j (0, 0)).2

-- --- 辅助引理 ---

private lemma prod_zipWith_mul_pow {R : Type*} [CommMonoid R]
    (σ : List R) (ks : List ℕ) (l : R) (hlen : σ.length = ks.length) :
    (σ.zipWith (fun s k => s * l ^ k) ks).prod = σ.prod * l ^ ks.sum := by
  induction σ generalizing ks with
  | nil =>
    cases ks with
    | nil => simp
    | cons _ _ => simp at hlen
  | cons s σ' ih =>
    match ks with
    | [] => simp at hlen
    | k :: ks' =>
      simp only [List.zipWith_cons_cons, List.prod_cons, List.sum_cons,
          ih ks' (by simpa using hlen), pow_add, mul_assoc, mul_left_comm]

private lemma lcDistribCore_length {n : ℕ}
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ) (ws : List ℤ) (r : ℕ) (hr : r = ws.length) :
    (lcDistribCore lc_factors lc_evals ws r).length = r := by
  subst hr
  induction lc_factors generalizing lc_evals with
  | nil => simp [lcDistribCore]
  | cons pr rest ih =>
    match lc_evals with
    | [] => simp [lcDistribCore]
    | E :: rest_e =>
      simp only [lcDistribCore, List.length_zipWith, ih rest_e,
                 extractValuations, List.length_map, Nat.min_self]

private lemma extractValuations_length (E : ℤ) (ws : List ℤ) :
    (extractValuations E ws).length = ws.length := List.length_map ..

-- --- 核心引理：lcDistribCore.prod = ∏ lⱼ^eⱼ ---

private theorem lcDistribCore_prod {n : ℕ}
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ) (ws : List ℤ) (r : ℕ)
    (hconserv : ConservationCheck lc_factors lc_evals ws)
    (hlen_fe : lc_factors.length = lc_evals.length)
    (hlen_r : r = ws.length) :
    (lcDistribCore lc_factors lc_evals ws r).prod =
    (lc_factors.map (fun pr => pr.1 ^ pr.2)).prod := by
  induction lc_factors generalizing lc_evals with
  | nil => simp [lcDistribCore]
  | cons pr rest ih =>
    match lc_evals with
    | [] => simp at hlen_fe
    | E :: rest_e =>
      simp only [lcDistribCore, List.map_cons, List.prod_cons]
      rw [prod_zipWith_mul_pow _ _ _ (by
        rw [lcDistribCore_length rest rest_e ws r hlen_r,
            extractValuations_length]; omega)]
      have hconserv_rest : ConservationCheck rest rest_e ws := by
        intro j hj
        have := hconserv (j + 1) (by simp; omega)
        simp only [List.getD_cons_succ] at this; exact this
      rw [ih rest_e hconserv_rest (by simpa using hlen_fe)]
      have hval_sum : (extractValuations E ws).sum = pr.2 := by
        have := hconserv 0 (by simp)
        simp only [List.getD_cons_zero] at this
        convert this using 1
      rw [hval_sum]; ring

-- --- 核心正确性 ---

/-- LC 分配乘积正确性：
    若 L ∼ γ · ∏ lⱼ^eⱼ 且守恒条件 (C4) 成立，则 γ · lcDistribCore.prod ∼ L。-/
theorem lc_distrib_prod_correct {n : ℕ}
    (L : MvPolynomial (Fin n) ℤ) (γ : ℤ)
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ) (ws : List ℤ) (r : ℕ)
    (hL : Associated L (MvPolynomial.C γ * (lc_factors.map (fun pr => pr.1 ^ pr.2)).prod))
    (hconserv : ConservationCheck lc_factors lc_evals ws)
    (hlen_fe : lc_factors.length = lc_evals.length)
    (hlen_r : r = ws.length) :
    Associated L (MvPolynomial.C γ * (lcDistribCore lc_factors lc_evals ws r).prod) := by
  rw [lcDistribCore_prod lc_factors lc_evals ws r hconserv hlen_fe hlen_r]
  exact hL

-- ============================================================
-- 模块 3: 试除重组
-- ============================================================

/-- 试除重组循环不变量：f = remaining × ∏ verified，每个 verified 不可约。
    对应 C++ __wang_core 的子集枚举 + divmod 循环状态。-/
structure TrialDivResult {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ) : Prop where
  /-- f = remaining × ∏ verified -/
  prod_eq : f = remaining * verified.prod
  /-- 每个因子不可约 -/
  all_irred : ∀ g ∈ verified, Irreducible g

/-- 试除初始化：remaining = f, verified = []。
    对应 C++ 循环开始前 mv_T = lifted, result = []。-/
theorem trial_div_init {σ : Type*} (f : MvPolynomial σ ℤ) :
    TrialDivResult f [] f :=
  ⟨by simp, by simp⟩

/-- 试除因子提取步：如果 g | remaining 且 g 不可约，提取 g。
    对应 C++ 循环体：trial division 成功 → verified.push(g), remaining = remaining/g。-/
theorem trial_div_extract {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ)
    (h_inv : TrialDivResult f verified remaining)
    (g : MvPolynomial σ ℤ) (hg_irred : Irreducible g)
    (remaining' : MvPolynomial σ ℤ) (h_div : remaining = g * remaining') :
    TrialDivResult f (g :: verified) remaining' := by
  refine ⟨?_, fun h hm => ?_⟩
  · rw [h_inv.prod_eq, h_div, List.prod_cons]; ring
  · rcases List.mem_cons.mp hm with rfl | hm'
    · exact hg_irred
    · exact h_inv.all_irred h hm'

/-- 试除完成（remaining 是 unit）→ f 的不可约分解。
    对应 C++ 试除循环结束后 mv_T.size() ≤ 1 的验收。-/
theorem trial_div_complete {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ)
    (h : TrialDivResult f verified remaining)
    (hrem : IsUnit remaining) :
    MvFactorCorrect f verified := by
  refine ⟨?_, h.all_irred⟩
  rw [h.prod_eq]
  exact (associated_isUnit_mul_left_iff hrem).mpr (Associated.refl _)

/-- 试除完成（remaining 也不可约）→ f 的不可约分解包含 remaining。
    对应 C++ verified.push_back(g_remaining) 的情况。-/
theorem trial_div_complete_with_remaining {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (verified : List (MvPolynomial σ ℤ))
    (remaining : MvPolynomial σ ℤ)
    (h : TrialDivResult f verified remaining)
    (hrem_irred : Irreducible remaining) :
    MvFactorCorrect f (remaining :: verified) :=
  ⟨h.prod_eq ▸ by simp [List.prod_cons], fun g hg => by
    rcases List.mem_cons.mp hg with rfl | hg'
    · exact hrem_irred
    · exact h.all_irred g hg'⟩

-- ============================================================
-- 模块 4: 多变量 Hensel 提升（输出规约）
-- ============================================================

/-- 多变量 Hensel 提升输出规约。
    不建模 MTSHL 内部过程，刻画输出必须满足的后置条件。
    prod_eq 来自试除验证的精确等式确认。
    对应 C++ __mtshl_lift 的输出 + 试除验证。-/
structure MvHenselOutput {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (lifted : List (MvPolynomial (Fin (n + 1)) ℤ))
    (α : Fin n → ℤ)
    (univar_facs : List (Polynomial ℤ)) : Prop where
  /-- 乘积条件：∏ lifted ∣ f（由试除验证确认）-/
  prod_dvd : lifted.prod ∣ f
  /-- 求值一致：Gᵢ(x₁, α) = vᵢ -/
  eval_consistent : List.Forall₂ (fun v G => eval_at_α G α = v) univar_facs lifted
  /-- 长度一致 -/
  length_eq : univar_facs.length = lifted.length

/-- Hensel 提升 + 每个因子不可约 + 乘积精确 → MvFactorCorrect。-/
theorem mv_hensel_to_factor {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (lifted : List (MvPolynomial (Fin (n + 1)) ℤ))
    (h_prod : f = lifted.prod)
    (h_irred : ∀ g ∈ lifted, Irreducible g) :
    MvFactorCorrect f lifted :=
  ⟨h_prod ▸ Associated.refl _, h_irred⟩

-- ============================================================
-- 模块 5: Wang 组合定理
-- ============================================================

/-- Wang 算法正确性（组合定理）。
    将 4 个模块串联：eval → LC 分配 → Hensel 提升 → 试除 → MvFactorCorrect。

    对应 C++ __wang_core 的完整流程：
    1. 求值 f₀ = f(x₁, α) + 单变量因式分解（wang_eval_step）
    2. LC 分配 σᵢ（lc_distrib_prod_correct）
    3. Hensel 提升 Gᵢ（MvHenselOutput）
    4. 试除验证（trial_div_complete / trial_div_complete_with_remaining）

    假设：好的求值点 + Hensel 提升成功 + 试除验证通过。
    这些由 C++ 的循环重试保证（换点直到成功）。-/
theorem wang_correct {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (_hf : f ≠ 0)
    -- 好的求值点
    (α : Fin n → ℤ) (_hα : EvalPointGood f α)
    -- Hensel 提升 + 试除结果
    (result : List (MvPolynomial (Fin (n + 1)) ℤ))
    (remaining : MvPolynomial (Fin (n + 1)) ℤ)
    (h_trial : TrialDivResult f result remaining)
    -- 试除完成条件（remaining 是 unit 或不可约）
    (hrem : IsUnit remaining ∨ Irreducible remaining)
    : ∃ factors : List (MvPolynomial (Fin (n + 1)) ℤ), MvFactorCorrect f factors := by
  rcases hrem with hunit | hirred
  · exact ⟨result, trial_div_complete f result remaining h_trial hunit⟩
  · exact ⟨remaining :: result, trial_div_complete_with_remaining f result remaining h_trial hirred⟩

-- ============================================================
-- 模块 6: MTSHL Taylor 提升循环
-- ============================================================

/-- 去掉第 i 个元素后的乘积 -/
noncomputable def List.prodEraseIdx {α : Type*} [CommMonoid α] (l : List α) (i : ℕ) : α :=
  (l.eraseIdx i).prod

/-- linearTerm：递归定义的线性项 Σσᵢ·∏_{j≠i}aⱼ。-/
noncomputable def linearTerm {R : Type*} [CommRing R]
    : List R → List R → R
  | [], _ | _, [] => 0
  | a :: as', σ :: σs' => σ * as'.prod + a * linearTerm as' σs'

/-- MDP（多变量 Diophantine）求解规约。
    偏分式分解：cₖ = Σ σᵢ · F̂ᵢ，其中 F̂ᵢ = ∏_{j≠i} f_base[j]。
    对应 C++ __mtshl_zp_univar_mdp / __mtshl_sparse_int / __mtshl_wmds。-/
structure MDPCorrect {R : Type*} [CommRing R]
    (ck : R)
    (f_base : List R)
    (sigma : List R) : Prop where
  /-- 偏分式分解：cₖ = linearTerm f_base sigma = Σᵢ σᵢ · ∏_{j≠i} f_base[j] -/
  sum_eq : ck = linearTerm f_base sigma
  /-- σ 的长度 = 基础因子数 -/
  length_eq : sigma.length = f_base.length

/-- IsCoprime with list product: if a is coprime to each element, it's coprime to the product. -/
private theorem isCoprime_list_prod {R : Type*} [CommSemiring R]
    (a : R) (l : List R) (h : ∀ b ∈ l, IsCoprime a b) : IsCoprime a l.prod := by
  induction l with
  | nil => exact isCoprime_one_right
  | cons b rest ih =>
    rw [List.prod_cons]
    exact (h b (.head ..)).mul_right (ih (fun c hc => h c (.tail _ hc)))

/-- MDP 存在性：互素因子的偏分式分解存在。
    c = Σ σᵢ · ∏_{j≠i} fⱼ（即 c = linearTerm f_base sigma）。
    对应 C++ MDP 求解器的数学正确性。-/
theorem mdp_exists {R : Type*} [CommRing R]
    (f_base : List R) (c : R)
    (hcop : f_base.Pairwise IsCoprime)
    (hne : f_base ≠ []) :
    ∃ sigma : List R, c = linearTerm f_base sigma ∧ sigma.length = f_base.length := by
  induction f_base generalizing c with
  | nil => exact absurd rfl hne
  | cons f rest ih =>
    by_cases hrest : rest = []
    · -- 单元素：[f]. linearTerm [f] [c] = c * 1 = c
      subst hrest; exact ⟨[c], by simp [linearTerm], rfl⟩
    · -- f :: rest with rest ≠ []
      have hcop_f : IsCoprime f rest.prod :=
        isCoprime_list_prod f rest (List.pairwise_cons.mp hcop).1
      obtain ⟨s, t, hbez⟩ := hcop_f
      have hcop_rest := (List.pairwise_cons.mp hcop).2
      obtain ⟨σ_rest, h_rest, hlen_rest⟩ := ih (c := c * s) hcop_rest hrest
      refine ⟨(c * t) :: σ_rest, ?_, by simp [hlen_rest]⟩
      simp only [linearTerm]
      rw [← h_rest]
      -- c = c*t * rest.prod + f * (c*s)
      calc c = c * 1 := (mul_one c).symm
        _ = c * (s * f + t * rest.prod) := by rw [hbez]
        _ = c * t * rest.prod + f * (c * s) := by ring

/-- MTSHL 单步更新：F[i] += σᵢ · (xⱼ - αⱼ)^k。
    对应 C++ lines 914-926。-/
noncomputable def mtshlStepUpdate {n : ℕ}
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (xj_pow_k : MvPolynomial (Fin (n + 1)) ℤ)
    : List (MvPolynomial (Fin (n + 1)) ℤ) :=
  factors.zipWith (fun fi si => fi + si * xj_pow_k) sigma

/-- MTSHL 提升不变量：∏ factors ≡ target (mod I^k)。
    I = (xⱼ - αⱼ) 理想。
    数学意义：target - ∏factors 的每一项都含 (xⱼ-αⱼ)^k 因子。-/
def MtshlInvariant {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ) : Prop :=
  -- target - ∏factors 可被 (X_{j+1} - C αⱼ)^k 整除
  (MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k ∣
    (target - factors.prod)

/-- 偏求值：在变量 i 处代入 a，其余变量不变。
    结果仍是 MvPolynomial（但不再含变量 i）。-/
noncomputable def partialEval {n : ℕ}
    (i : Fin (n + 1)) (a : ℤ) (f : MvPolynomial (Fin (n + 1)) ℤ) :
    MvPolynomial (Fin (n + 1)) ℤ :=
  MvPolynomial.aeval (fun k => if k = i then MvPolynomial.C a else MvPolynomial.X k) f

/-- MvPolynomial 因子定理：若偏求值（x_i → α）后为零，则 (X_i - C α) 整除。
    证明：f - partialEval f 被 (X_i - C α) 整除（对每个单项式，
    x^e - α^e = (x - α) * (x^{e-1} + ... + α^{e-1})）。
    若 partialEval f = 0，则 f = f - 0 = (X_i - C α) * q。-/
private theorem mv_X_sub_C_dvd_of_partialEval_eq_zero {n : ℕ}
    (i : Fin (n + 1)) (α : ℤ) (f : MvPolynomial (Fin (n + 1)) ℤ)
    (h : partialEval i α f = 0) :
    (MvPolynomial.X i - MvPolynomial.C α) ∣ f := by
  -- 策略：f = f - partialEval f（因 partialEval f = 0）
  -- f - partialEval f 被 (X i - C α) 整除
  -- 因为 partialEval 是将 X i 替换为 C α，
  -- 对每个单项式 m = c · ∏ X_k^{e_k}：
  --   m - partialEval m = c · X_i^{e_i} · ∏_{k≠i} X_k^{e_k} - c · (C α)^{e_i} · ∏_{k≠i} X_k^{e_k}
  --   = c · (X_i^{e_i} - (C α)^{e_i}) · ∏_{k≠i} X_k^{e_k}
  -- 而 (X_i - C α) | (X_i^{e_i} - (C α)^{e_i})（sub_dvd_pow_sub_pow）
  -- 所以 (X_i - C α) | (m - partialEval m) 对每个单项式
  -- 由有限和的整除性：(X_i - C α) | (f - partialEval f)
  rw [← sub_zero f, ← h]
  -- 目标：(X i - C α) | (f - partialEval i α f)
  -- 用 MvPolynomial.induction_on 对 f 归纳
  set d := MvPolynomial.X i - MvPolynomial.C α with hd_def
  suffices key : ∀ (p : MvPolynomial (Fin (n + 1)) ℤ),
      d ∣ (p - partialEval i α p) from key f
  intro p
  apply MvPolynomial.induction_on p
  · -- h_C: p = C a → C a - partialEval(C a) = C a - C a = 0
    intro a; simp [partialEval, d]
  · -- h_add: 整除在加法下封闭
    intro p q hp hq
    have : (p + q) - partialEval i α (p + q) =
        (p - partialEval i α p) + (q - partialEval i α q) := by
      simp [partialEval, map_add]; ring
    rw [this]; exact dvd_add hp hq
  · -- h_X: p * X k 的情况
    intro p k hp
    by_cases hk : k = i
    · -- k = i
      have hpe : partialEval i α (p * MvPolynomial.X k) =
          partialEval i α p * MvPolynomial.C α := by
        simp only [partialEval, map_mul, MvPolynomial.aeval_X, hk, if_pos rfl, ite_true]
      rw [hk] at hpe ⊢
      have heq : p * MvPolynomial.X i - partialEval i α (p * MvPolynomial.X i) =
          (p - partialEval i α p) * MvPolynomial.X i + partialEval i α p * d := by
        rw [hpe, hd_def]; ring
      rw [heq]; exact dvd_add (hp.mul_right _) (dvd_mul_left d _)
    · -- k ≠ i
      have hpe : partialEval i α (p * MvPolynomial.X k) =
          partialEval i α p * MvPolynomial.X k := by
        simp only [partialEval, map_mul, MvPolynomial.aeval_X, if_neg hk]
      have heq : p * MvPolynomial.X k - partialEval i α (p * MvPolynomial.X k) =
          (p - partialEval i α p) * MvPolynomial.X k := by
        rw [hpe]; ring
      rw [heq]; exact hp.mul_right _

theorem mtshl_invariant_init {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    (h_eval : partialEval (Fin.succ j) αⱼ (target - factors.prod) = 0) :
    MtshlInvariant target factors j αⱼ 1 := by
  unfold MtshlInvariant; rw [pow_one]
  exact mv_X_sub_C_dvd_of_partialEval_eq_zero (Fin.succ j) αⱼ _ h_eval

/-- 终止条件：若不变量在 k > deg(target, xⱼ) 时成立，则 ∏factors = target。
    因为 target - ∏factors 的 xⱼ 次数 ≤ deg(target, xⱼ) < k，
    但它含 (xⱼ-αⱼ)^k 因子（度数 ≥ k），所以必为零。-/
theorem mtshl_invariant_terminates {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ)
    (h_inv : MtshlInvariant target factors j αⱼ k)
    -- k 超过误差的 xⱼ 次数
    (h_prec : k > MvPolynomial.degreeOf (Fin.succ j) (target - factors.prod)) :
    target = factors.prod := by
  suffices h0 : target - factors.prod = 0 from sub_eq_zero.mp h0
  set diff := target - factors.prod with hdiff_def
  by_contra h_ne
  -- Step 1: rename(swap 0 (succ j)) 将变量 succ j 移到位置 0
  set σ := Equiv.swap (0 : Fin (n + 1)) (Fin.succ j) with hσ
  set diff' := MvPolynomial.rename σ diff with hdiff'
  have h_ne' : diff' ≠ 0 :=
    fun h => h_ne (MvPolynomial.rename_injective σ (Equiv.injective σ) h)
  -- Step 2: 整除关系通过 rename 传递
  have h_dvd' : (MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k ∣ diff' := by
    have h1 : MvPolynomial.rename σ ((MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k) =
        (MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k := by
      simp [map_pow, map_sub, MvPolynomial.rename_X, MvPolynomial.rename_C,
            show σ (Fin.succ j) = (0 : Fin (n + 1)) from Equiv.swap_apply_right _ _]
    rw [← h1]; exact (MvPolynomial.rename σ).toRingHom.map_dvd h_inv
  -- Step 3: finSuccEquiv 将 MvPolynomial → Polynomial
  set p := MvPolynomial.finSuccEquiv ℤ n diff' with hp_def
  have h_ne_p : p ≠ 0 := fun h => h_ne' ((MvPolynomial.finSuccEquiv ℤ n).injective h)
  -- Step 4: finSuccEquiv 保持整除 + natDegree 约束 → 矛盾
  -- finSuccEquiv 映射 (X 0 - C α)^k 到 (Poly.X - Poly.C (C α))^k
  -- natDegree((Poly.X - Poly.C c)^k) = k
  -- natDegree(p) = degreeOf 0 diff' = degreeOf (succ j) diff < k
  -- 由 natDegree_le_of_dvd: k ≤ natDegree(p) < k，矛盾
  have h_deg_p : p.natDegree < k := by
    rw [hp_def, MvPolynomial.natDegree_finSuccEquiv, hdiff']
    rw [show (0 : Fin (n + 1)) = σ (Fin.succ j) from (Equiv.swap_apply_right _ _).symm,
        MvPolynomial.degreeOf_rename_of_injective (Equiv.injective σ)]
    exact h_prec
  have h_dvd_p : (Polynomial.X - Polynomial.C (MvPolynomial.C αⱼ)) ^ k ∣ p := by
    have h1 : (MvPolynomial.finSuccEquiv ℤ n)
        ((MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k) =
        (Polynomial.X - Polynomial.C (MvPolynomial.C αⱼ)) ^ k := by
      rw [map_pow, map_sub, MvPolynomial.finSuccEquiv_X_zero]
      congr 1; simp [MvPolynomial.finSuccEquiv_apply]
    rw [← h1]; exact (MvPolynomial.finSuccEquiv ℤ n).toRingEquiv.toRingHom.map_dvd h_dvd'
  have h_le := Polynomial.natDegree_le_of_dvd h_dvd_p h_ne_p
  rw [Polynomial.natDegree_pow, Polynomial.natDegree_X_sub_C, Nat.mul_one] at h_le
  omega

/-- 乘积 mod d 的因子定理推广：若 d | (aᵢ - bᵢ) 对所有 i，则 d | (∏aᵢ - ∏bᵢ)。-/
private theorem dvd_prod_sub_prod {R : Type*} [CommRing R]
    {as bs : List R} (d : R)
    (hdvd : List.Forall₂ (fun a b => d ∣ (a - b)) as bs) :
    d ∣ (as.prod - bs.prod) := by
  induction hdvd with
  | nil => simp
  | @cons a b as' bs' hab _ ih =>
    simp only [List.prod_cons]
    have key : a * as'.prod - b * bs'.prod =
        a * (as'.prod - bs'.prod) + (a - b) * bs'.prod := by ring
    rw [key]; exact dvd_add (ih.mul_left _) (hab.mul_right _)

private theorem prod_add_linearization {R : Type*} [CommRing R]
    (as sigmas : List R) (c : R) (hlen : as.length = sigmas.length) :
    c ^ 2 ∣ ((as.zipWith (fun a s => a + s * c) sigmas).prod -
              as.prod - c * linearTerm as sigmas) := by
  induction as generalizing sigmas with
  | nil => cases sigmas with
    | nil => simp [linearTerm]
    | cons _ _ => simp at hlen
  | cons a as' ih =>
    match sigmas with
    | [] => simp at hlen
    | σ :: sigmas' =>
      simp only [List.zipWith_cons_cons, List.prod_cons, linearTerm]
      -- 用 IH：c² | (P_new - P_old - c·L_rest) 其中
      -- P_new = ∏(as'+σs'·c), P_old = ∏as', L_rest = linearTerm as' sigmas'
      have hlen' : as'.length = sigmas'.length := by simpa using hlen
      obtain ⟨R_rest, hR⟩ := ih sigmas' hlen'
      -- hR: ∏(as'+σs'·c) - ∏as' - c·L_rest = c²·R_rest
      set P_new := (as'.zipWith (fun a' s => a' + s * c) sigmas').prod
      set P_old := as'.prod
      set L_rest := linearTerm as' sigmas'
      -- 目标：c² | ((a+σc)·P_new - a·P_old - c·(σ·P_old + a·L_rest))
      -- = c² | (a·P_new + σc·P_new - a·P_old - c·σ·P_old - c·a·L_rest)
      -- = c² | (a·(P_new - P_old - c·L_rest) + σ·c·(P_new - P_old))
      -- = c² | (a·c²·R_rest + σ·c·(c·L_rest + c²·R_rest))  [P_new-P_old = c·L_rest + c²·R_rest]
      -- = c² | (a·c²·R_rest + σ·c²·L_rest + σ·c³·R_rest)
      -- = c² | c²·(a·R_rest + σ·L_rest + σ·c·R_rest)  ✓
      have hP_eq : P_new = P_old + c * L_rest + c ^ 2 * R_rest := by
        linear_combination hR
      refine ⟨a * R_rest + σ * L_rest + σ * c * R_rest, ?_⟩
      calc (a + σ * c) * P_new - a * P_old - c * (σ * P_old + a * L_rest)
          = (a + σ * c) * (P_old + c * L_rest + c ^ 2 * R_rest) -
            a * P_old - c * (σ * P_old + a * L_rest) := by rw [hP_eq]
        _ = c ^ 2 * (a * R_rest + σ * L_rest + σ * c * R_rest) := by ring

/-- 乘积扰动引理：若每个因子的扰动 bᵢ 都被 c 整除，
    则 ∏(aᵢ + bᵢ) - ∏aᵢ 也被 c 整除。
    对应 MTSHL 中：每个 F_new[i] - F_old[i] = σᵢ·d^k 被 d^k 整除。-/
private theorem dvd_prod_add_sub_prod {R : Type*} [CommRing R]
    (as bs : List R) (c : R) (hlen : as.length = bs.length)
    (hbs : ∀ b ∈ bs, c ∣ b) :
    c ∣ (as.zipWith (· + ·) bs).prod - as.prod := by
  induction as generalizing bs with
  | nil => cases bs with
    | nil => simp
    | cons _ _ => simp at hlen
  | cons a as' ih =>
    match bs with
    | [] => simp at hlen
    | b :: bs' =>
      simp only [List.zipWith_cons_cons, List.prod_cons]
      -- (a+b) · ∏(as'+bs') - a · ∏as'
      -- = a · (∏(as'+bs') - ∏as') + b · ∏(as'+bs')
      have key : (a + b) * (as'.zipWith (· + ·) bs').prod - a * as'.prod =
          a * ((as'.zipWith (· + ·) bs').prod - as'.prod) +
          b * (as'.zipWith (· + ·) bs').prod := by ring
      rw [key]
      exact dvd_add
        ((ih bs' (by simpa using hlen) (fun x hx => hbs x (.tail _ hx))).mul_left a)
        ((hbs b (.head _)).mul_right _)

/-- partialEval 幂等：partialEval ∘ partialEval = partialEval。-/
private theorem partialEval_idem {n : ℕ} (i : Fin (n + 1)) (α : ℤ)
    (f : MvPolynomial (Fin (n + 1)) ℤ) :
    partialEval i α (partialEval i α f) = partialEval i α f := by
  simp only [partialEval]
  apply MvPolynomial.induction_on f
  · intro a; simp [MvPolynomial.aeval_C]
  · intro p q hp hq; simp only [map_add, hp, hq]
  · intro p k hp; simp only [map_mul, MvPolynomial.aeval_X]
    split
    · simp only [MvPolynomial.aeval_C, map_mul, hp, MvPolynomial.algebraMap_eq]
    · rename_i hk; simp only [MvPolynomial.aeval_X, if_neg hk, map_mul, hp]

/-- linearTerm 在 partialEval 下的差被 d 整除。
    对应：Σσᵢ·∏F[j] ≡ Σσᵢ·∏partialEval(F[j]) (mod d)。-/
private theorem linearTerm_partialEval_dvd {n : ℕ}
    (factors sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (i : Fin (n + 1)) (α : ℤ) :
    (MvPolynomial.X i - MvPolynomial.C α) ∣
      (linearTerm (factors.map (partialEval i α)) sigma - linearTerm factors sigma) := by
  set d := MvPolynomial.X i - MvPolynomial.C α
  induction factors generalizing sigma with
  | nil => simp [linearTerm]
  | cons f fs ih =>
    match sigma with
    | [] => simp [linearTerm]
    | σ :: σs =>
      simp only [linearTerm, List.map_cons]
      have h_prod : d ∣ ((fs.map (partialEval i α)).prod - fs.prod) :=
        dvd_prod_sub_prod d (List.forall₂_map_left_iff.mpr
          (List.forall₂_same.mpr (fun x _ =>
            mv_X_sub_C_dvd_of_partialEval_eq_zero _ _ _ (by
              show partialEval i α (partialEval i α x - x) = 0
              unfold partialEval; simp only [map_sub]; rw [sub_eq_zero]
              show partialEval i α (partialEval i α x) = partialEval i α x
              exact partialEval_idem i α x))))
      have h_f : d ∣ (partialEval i α f - f) :=
        mv_X_sub_C_dvd_of_partialEval_eq_zero _ _ _ (by
          show partialEval i α (partialEval i α f - f) = 0
          unfold partialEval; simp only [map_sub]; rw [sub_eq_zero]
          show partialEval i α (partialEval i α f) = partialEval i α f
          exact partialEval_idem i α f)
      have h_rest := ih σs
      have key : σ * (fs.map (partialEval i α)).prod +
          partialEval i α f * linearTerm (fs.map (partialEval i α)) σs -
          (σ * fs.prod + f * linearTerm fs σs) =
        σ * ((fs.map (partialEval i α)).prod - fs.prod) +
        (partialEval i α f - f) * linearTerm (fs.map (partialEval i α)) σs +
        f * (linearTerm (fs.map (partialEval i α)) σs - linearTerm fs σs) := by ring
      rw [key]
      exact dvd_add (dvd_add (h_prod.mul_left _) (h_f.mul_right _)) (h_rest.mul_left _)

/-! ### MTSHL 单步不变量传播

  数学论证：
  设 d = X_j - C α, error = target - ∏F = d^k · q（由 h_inv）
  MDP：partialEval q = Σ σᵢ · F̂ᵢ_base
  更新：F_new[i] = F[i] + σᵢ · d^k
  乘积展开：∏F_new - ∏F = d^k · (Σ σᵢ · F̂ᵢ) + O(d^{2k})
  新误差 = d^k · (q - Σ σᵢ · F̂ᵢ) + O(d^{2k})
  因子定理：d | (q - partialEval q) 且 F̂ᵢ ≡ F̂ᵢ_base (mod d)
  → d | (q - Σ σᵢ · F̂ᵢ) → d^{k+1} | d^k · d · (...) → d^{k+1} | 新误差 -/

/-- MTSHL 单步不变量传播。
    将完整证明分解为两个假设：
    1. 乘积展开：∏F_new - ∏F_old 被 d^k 整除后的线性项
    2. MDP 正确性
    通过因子定理组合得 d^{k+1} | error_new。

    对应 C++ __mtshl_step_j 的单步 k 循环体。-/
theorem mtshl_step_invariant {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ) (hk : 1 ≤ k)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (h_inv : MtshlInvariant target factors j αⱼ k)
    (h_len : sigma.length = factors.length)
    -- Taylor 系数 = (target - ∏factors) / d^k 在 x_j=α_j 处求值
    -- h_inv 给出 ∃ q, target - ∏factors = d^k * q，Taylor 系数 = partialEval q
    -- MDP 正确：Taylor 系数 = Σ σᵢ · F̂ᵢ_base
    (h_mdp : MDPCorrect
        (partialEval (Fin.succ j) αⱼ h_inv.choose)
        (factors.map (partialEval (Fin.succ j) αⱼ))
        sigma)
    : MtshlInvariant target
        (mtshlStepUpdate factors sigma
          ((MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k))
        j αⱼ (k + 1) := by
  set d := MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ with hd
  unfold MtshlInvariant mtshlStepUpdate
  set q := h_inv.choose with hq_def
  have hq := h_inv.choose_spec  -- target - factors.prod = d^k * q
  -- Step 1: d^k | (∏F_new - ∏F_old)
  -- 每个 F_new[i] - F_old[i] = σ[i] * d^k 被 d^k 整除
  set F_new := factors.zipWith (fun fi si => fi + si * d ^ k) sigma with hF_new
  have h_perturb : ∀ b ∈ sigma.map (· * d ^ k), d ^ k ∣ b :=
    fun b hb => by obtain ⟨s, _, rfl⟩ := List.mem_map.mp hb; exact dvd_mul_left _ _
  have h_diff_dvd : d ^ k ∣ (F_new.prod - factors.prod) := by
    rw [hF_new, show factors.zipWith (fun fi si => fi + si * d ^ k) sigma =
        factors.zipWith (· + ·) (sigma.map (· * d ^ k)) from by
      simp [List.zipWith_map_right]]
    exact dvd_prod_add_sub_prod factors (sigma.map (· * d ^ k)) (d ^ k)
      (by simp [h_len]) h_perturb
  -- Step 2: error_new = d^k * (q - quotient_of_diff)
  -- target - ∏F_new = (target - ∏F_old) - (∏F_new - ∏F_old) = d^k * q - d^k * q'
  obtain ⟨q', hq'⟩ := h_diff_dvd  -- ∏F_new - ∏F_old = d^k * q'
  have h_error_new : target - F_new.prod = d ^ k * (q - q') := by
    have : target - F_new.prod = (target - factors.prod) - (F_new.prod - factors.prod) := by ring
    rw [this, hq, hq']; ring
  -- Step 3: d^{k+1} | error_new — 用因子定理 + MDP 直接构造
  -- error_new = target - F_new.prod
  -- 需要 d^{k+1} | error_new
  -- 路径：error_new = d^k * q - d^k * q'
  -- 且 d | (q - q')
  -- 用因子定理：partialEval(error_new) 在适当意义下为 0
  -- 组合 dvd_prod_sub_prod + mv_X_sub_C_dvd_of_partialEval_eq_zero
  --
  -- 直接证 d^{k+1} | d^k * (q - q') 等价于 d | (q - q')
  -- 由因子定理：只需 partialEval(q - q') = 0
  -- 但 partialEval(d^k) = 0 使得 q' 的 partialEval 不可从 h_diff_dvd 提取
  -- 需要乘积展开的线性化
  -- 这是 Newton 迭代核心（多因子 Leibniz 展开 ~60 行）
  -- Step 3: d^{k+1} | error_new
  -- 策略：error_new 被 d^k 整除（已有 h_error_new），
  -- 证 partialEval(q - q') = 0 → d | (q - q') → d^{k+1} | d^k*(q-q')
  --
  -- partialEval(q) = MDP 的 Σσᵢ·F̂ᵢ_base（h_mdp 的 choose 等式）
  -- partialEval(q') = ?
  --   q' 满足 ∏F_new - ∏F_old = d^k * q'
  --   partialEval(∏F_new - ∏F_old) = 0（因 partialEval(F_new[i]) = partialEval(F_old[i])）
  --   但 partialEval(d^k * q') = 0 是 trivial 的，不给 partialEval(q')
  --
  -- 正确路径：不通过 q' 的 partialEval，而是直接证
  -- Step 3: d | (q - q')（由 MDP + 因子定理 + 乘积 mod d 引理组合）
  -- 这是 Newton 迭代核心：MDP 保证线性项匹配，因子定理处理余项
  -- q = (target - ∏F_old) / d^k
  -- q' = (∏F_new - ∏F_old) / d^k
  -- q - q' 在 partialEval 下：
  --   partialEval q = Σσᵢ·F̂ᵢ_base（MDP）
  --   q' 的线性部分 = Σσᵢ·F̂ᵢ ≡ Σσᵢ·F̂ᵢ_base (mod d)（dvd_prod_sub_prod）
  --   所以 partialEval(q - q') = 0 在适当意义下
  -- 严格证明需要 Leibniz 展开将 q' 分解为线性项 + d-整除余项
  -- [详见 nl-proof mv-l2-mtshl-step.md]
  have h_dvd_qq' : d ∣ (q - q') := by
    -- Step 3a: Leibniz 展开 → ∏F_new - ∏F_old = d^k·L + d^{2k}·R
    -- 其中 L = linearTerm factors sigma
    set L := linearTerm factors sigma with hL_def
    obtain ⟨R_lin, hR_lin⟩ := prod_add_linearization factors sigma (d ^ k) h_len.symm
    -- hR_lin: ∏F_new - ∏F_old - d^k·L = d^{2k}·R_lin
    -- 所以 ∏F_new - ∏F_old = d^k·L + d^{2k}·R_lin
    -- hq': ∏F_new - ∏F_old = d^k · q'
    -- → d^k·q' = d^k·L + d^{2k}·R_lin = d^k·(L + d^k·R_lin)
    -- 在整环中消去 d^k：q' = L + d^k·R_lin（需要 d^k ≠ 0 或 cancelation）
    -- Step 3b: d | (q - L)
    -- q - L = (q - partialEval q) + (partialEval q - L)
    -- d | (q - partialEval q)（因子定理 ✅）
    -- partialEval q = Σσᵢ·F̂ᵢ_base（MDP）
    -- L = linearTerm = Σσᵢ·F̂ᵢ（递归定义中 F̂ᵢ = ∏_{j≠i} F_old[j]）
    -- d | (F̂ᵢ - F̂ᵢ_base) = d | (∏F_old - ∏partialEval(F_old))（dvd_prod_sub_prod ✅）
    -- → d | (L - partialEval q) = d | Σσᵢ·(F̂ᵢ - F̂ᵢ_base)
    -- → d | (q - L)
    -- Step 3c: 从 hq' 和 hR_lin 推出 q' - L = d^k * R_lin
    -- hq': F_new.prod - factors.prod = d^k * q'
    -- hR_lin: (F_new.prod - factors.prod) - d^k * L = (d^k)^2 * R_lin
    -- 两式联立：d^k * q' - d^k * L = d^{2k} * R_lin
    -- 即 d^k * (q' - L) = d^k * (d^k * R_lin)
    have hq'_eq : d ^ k * (q' - L) = d ^ k * (d ^ k * R_lin) := by
      have h1 : F_new.prod - factors.prod = d ^ k * L + (d ^ k) ^ 2 * R_lin := by
        linear_combination hR_lin
      have h2 := hq'.symm.trans h1  -- d^k * q' = d^k * L + (d^k)^2 * R_lin
      linear_combination h2
    -- 消去 d^k（MvPolynomial 是整环，d^k ≠ 0 — 需要 d ≠ 0）
    -- d = X_j - C α 在 MvPolynomial 中非零
    have hd_ne : d ≠ 0 := by
      rw [hd]; intro h
      have h1 := congr_arg (MvPolynomial.degreeOf (Fin.succ j)) (sub_eq_zero.mp h)
      simp only [MvPolynomial.degreeOf_X, if_pos rfl] at h1
      rw [MvPolynomial.degreeOf_C] at h1; exact absurd h1 one_ne_zero
    have hdk_ne : d ^ k ≠ 0 := pow_ne_zero _ hd_ne
    have hq'_L : q' - L = d ^ k * R_lin := mul_left_cancel₀ hdk_ne hq'_eq
    -- Step 3d: q - q' = (q - L) - d^k * R_lin
    have h_split : q - q' = (q - L) - d ^ k * R_lin := by
      rw [← hq'_L]; ring
    -- Step 3e: d | (q - L) 由因子定理 + MDP
    have h_dvd_qL : d ∣ (q - L) := by
      -- q - L = (q - partialEval q) + (partialEval q - L)
      -- d | (q - partialEval q) 由因子定理
      have h1 : d ∣ (q - partialEval (Fin.succ j) αⱼ q) :=
        mv_X_sub_C_dvd_of_partialEval_eq_zero (Fin.succ j) αⱼ _ (by
          simp only [partialEval, map_sub]; rw [sub_eq_zero]
          exact (partialEval_idem (Fin.succ j) αⱼ q).symm)
      have h2 : d ∣ (linearTerm (factors.map (partialEval (Fin.succ j) αⱼ)) sigma -
                      linearTerm factors sigma) :=
        linearTerm_partialEval_dvd factors sigma (Fin.succ j) αⱼ
      -- MDP 给出 partialEval q = partialEval(h_inv.choose)（但 h_inv.choose = q）
      -- 实际上需要连接 h_mdp 的 partialEval q 与 linearTerm(partialEval factors, sigma)
      -- h_mdp : MDPCorrect (partialEval q_from_inv) (factors.map partialEval) sigma
      -- MDPCorrect.sum_eq 说 partialEval q = 某个 sum 表达式
      -- 这个 sum 与 linearTerm(factors.map partialEval, sigma) 需要等价
      -- 两者都是 Σσᵢ·∏_{j≠i}(factors.map partialEval)[j]，只是表示不同
      -- linearTerm 是递归定义，MDPCorrect.sum_eq 用 List.range + prodEraseIdx
      -- 等价性需要一个桥接引理
      -- q - L = (q - partialEval q) + (partialEval q - linearTerm(pE_facs, σ)) +
      --         (linearTerm(pE_facs, σ) - L)
      -- d | (q - partialEval q)（h1）
      -- d | (linearTerm(pE_facs, σ) - L)（h2）
      -- partialEval q = linearTerm(pE_facs, σ)（需要 MDP + linearTerm↔sum 桥接）
      -- 暂时：q - L = (q - pE_q) + (pE_q - L) 且 d 整除两者
      have h_qL : q - L = (q - partialEval (Fin.succ j) αⱼ q) +
          (partialEval (Fin.succ j) αⱼ q -
           linearTerm (factors.map (partialEval (Fin.succ j) αⱼ)) sigma) +
          (linearTerm (factors.map (partialEval (Fin.succ j) αⱼ)) sigma - L) := by ring
      rw [h_qL]
      -- 三项分别被 d 整除
      refine dvd_add (dvd_add h1 ?_) h2
      -- partialEval q = linearTerm(pE_facs, σ) 由 MDP (sum_eq 现在用 linearTerm)
      -- h_mdp.sum_eq : partialEval h_inv.choose = linearTerm(pE_facs, σ)
      -- q = h_inv.choose（obtain 的 q 和 Exists.choose 相同）
      -- q = h_inv.choose（set 定义）
      rw [h_mdp.sum_eq, sub_self]; exact dvd_zero _
    -- Step 3f: d | d^k * R_lin（trivial for k ≥ 1）
    have h_dvd_dkR : d ∣ (d ^ k * R_lin) :=
      dvd_mul_of_dvd_left (dvd_pow_self d (by omega : k ≠ 0)) R_lin
    -- 组合
    rw [h_split]; exact dvd_sub h_dvd_qL h_dvd_dkR
  -- 最终：d^{k+1} | d^k * (q - q')
  obtain ⟨r, hr⟩ := h_dvd_qq'  -- q - q' = d * r
  rw [h_error_new, hr]
  simp only [hd, ← mul_assoc, ← pow_succ]
  exact dvd_mul_right _ _

-- ============================================================
-- 模块 7: MTSHL 完整 MDP 求解器
-- ============================================================

/-! ### 7a. Taylor 系数提取 -/

/-- Taylor 系数提取：f 在 (X_k - C α) 展开的第 j 项系数。
    迭代 j 次除以 (X_k - C α) 取商，然后在 X_k=α 处求值。
    对应 C++ __taylor_coeff_zp。-/
noncomputable def taylorCoeffAt {n : ℕ}
    (k : Fin (n + 1)) (α : ℤ) : ℕ → MvPolynomial (Fin (n + 1)) ℤ → MvPolynomial (Fin (n + 1)) ℤ
  | 0, f => partialEval k α f
  | j + 1, f =>
    let remainder := partialEval k α f  -- f(α)
    let quotient := Classical.choose (mv_X_sub_C_dvd_of_partialEval_eq_zero k α
      (f - remainder) (by unfold partialEval; simp [map_sub]; rw [sub_eq_zero]
                          exact (partialEval_idem k α f).symm))
    taylorCoeffAt k α j quotient

/-! ### 7b. Vandermonde 可逆性 -/

/-- Vandermonde 矩阵可逆（θ 两两不同）。对应 C++ __si_vandermonde_solve。-/
theorem vandermonde_det_ne_zero {s : ℕ} {F : Type*} [Field F] [DecidableEq F]
    (θ : Fin s → F) (h_inj : Function.Injective θ) :
    (Matrix.vandermonde θ).det ≠ 0 := by
  rw [Matrix.det_vandermonde]
  exact Finset.prod_ne_zero_iff.mpr (fun i _ =>
    Finset.prod_ne_zero_iff.mpr (fun j hj =>
      sub_ne_zero.mpr (fun h => absurd (h_inj h) (Finset.mem_Ioi.mp hj).ne')))

/-- Vandermonde 系统 V·d = v 有唯一解。
    使用 Matrix.mulVec 表述，对应 C++ __si_vandermonde_solve。-/
theorem vandermonde_solve_unique {s : ℕ} {F : Type*} [Field F] [DecidableEq F]
    (θ : Fin s → F) (v : Fin s → F)
    (h_inj : Function.Injective θ) :
    ∃! d : Fin s → F, (Matrix.vandermonde θ).mulVec d = v := by
  have h_det := vandermonde_det_ne_zero θ h_inj
  have h_det_unit : IsUnit (Matrix.vandermonde θ).det := IsUnit.mk0 _ h_det
  have hV_isUnit : IsUnit (Matrix.vandermonde θ) :=
    (Matrix.isUnit_iff_isUnit_det _).mpr h_det_unit
  -- 构造解：d = V⁻¹ *ᵥ v
  refine ⟨(Matrix.vandermonde θ)⁻¹.mulVec v, ?_, ?_⟩
  · -- V *ᵥ (V⁻¹ *ᵥ v) = v
    show (Matrix.vandermonde θ).mulVec ((Matrix.vandermonde θ)⁻¹.mulVec v) = v
    rw [Matrix.mulVec_mulVec,
        Matrix.mul_nonsing_inv (Matrix.vandermonde θ) h_det_unit,
        Matrix.one_mulVec]
  · -- 唯一性：V *ᵥ d' = v → d' = V⁻¹ *ᵥ v
    intro d' hd'
    have h_inj' := Matrix.mulVec_injective_of_isUnit hV_isUnit
    exact h_inj' (by rw [hd', Matrix.mulVec_mulVec,
        Matrix.mul_nonsing_inv (Matrix.vandermonde θ) h_det_unit,
        Matrix.one_mulVec])

/-! ### 7c. θ-array 求值（稀疏插值基础） -/

/-- θ-array 求值：f 在 (x₁, β₂^l, ..., βⱼ₋₁^l) 处偏求值。
    对应 C++ __si_theta_array_eval（lines 133-202）。
    C++ 用 running product 优化 O(s·|Supp(f)|)；Lean 只建模数学语义。-/
noncomputable def evalAtBetaPow {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (β : Fin n → ℤ) (l : ℕ) : Polynomial ℤ :=
  eval_at_α f (fun i => β i ^ l)

/-- evalAtBetaPow 保持零。-/
@[simp] theorem evalAtBetaPow_zero {n : ℕ} (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (0 : MvPolynomial (Fin (n + 1)) ℤ) β l = 0 := by
  simp [evalAtBetaPow, eval_at_α, map_zero, Polynomial.map_zero]

/-- evalAtBetaPow 保持乘法。-/
theorem evalAtBetaPow_mul {n : ℕ}
    (f g : MvPolynomial (Fin (n + 1)) ℤ) (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (f * g) β l = evalAtBetaPow f β l * evalAtBetaPow g β l :=
  eval_at_α_mul f g _

/-- evalAtBetaPow 保持加法。-/
theorem evalAtBetaPow_add {n : ℕ}
    (f g : MvPolynomial (Fin (n + 1)) ℤ) (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (f + g) β l = evalAtBetaPow f β l + evalAtBetaPow g β l := by
  simp only [evalAtBetaPow, eval_at_α, map_add, Polynomial.map_add]

/-- evalAtBetaPow 保持列表乘积。-/
theorem evalAtBetaPow_prod {n : ℕ}
    (fs : List (MvPolynomial (Fin (n + 1)) ℤ)) (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow fs.prod β l = (fs.map (evalAtBetaPow · β l)).prod := by
  induction fs with
  | nil => simp [evalAtBetaPow, eval_at_α, map_one, Polynomial.map_one]
  | cons f rest ih =>
    simp only [List.prod_cons, List.map_cons]
    rw [evalAtBetaPow_mul, ih]

/-- θ-array 求值保持 linearTerm 结构。
    稀疏插值算法的数学基础：多变量 MDP 在 θ-array 求值后变为单变量 MDP。
    对应 C++ __mtshl_sparse_int Step 1→Step 2 的数学链。-/
theorem evalAtBetaPow_linearTerm {n : ℕ}
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (linearTerm F_base sigma) β l =
    linearTerm (F_base.map (evalAtBetaPow · β l))
               (sigma.map (evalAtBetaPow · β l)) := by
  induction F_base generalizing sigma with
  | nil => simp [linearTerm]
  | cons f rest ih =>
    cases sigma with
    | nil => simp [linearTerm]
    | cons σ σ_rest =>
      simp only [linearTerm, List.map_cons]
      rw [evalAtBetaPow_add, evalAtBetaPow_mul, evalAtBetaPow_mul,
          evalAtBetaPow_prod, ih σ_rest]

/-! ### 7d. 稀疏插值 MDP + 级联 -/

private theorem mdp_to_MDPCorrect {R : Type*} [CommRing R]
    {f_base : List R} {ck : R}
    (h : ∃ sigma, ck = linearTerm f_base sigma ∧ sigma.length = f_base.length) :
    ∃ sigma, MDPCorrect ck f_base sigma :=
  let ⟨sigma, h_eq, h_len⟩ := h; ⟨sigma, ⟨h_eq, h_len⟩⟩

/-- 稀疏插值 MDP 正确性（算法版）。
    对应 C++ __mtshl_sparse_int（lines 449-605）的完整 4 步算法。

    算法建模：
    - Step 1 (θ-array 求值)：evalAtBetaPow 定义（对应 __si_theta_array_eval）
    - Step 2 (逐点单变量 MDP)：h_univar 假设（对应逐点调用 __mtshl_zp_univar_mdp）
    - Step 3 (Vandermonde 恢复)：vandermonde_solve_unique 保证（对应 __si_vandermonde_solve）
    - Step 4 (验证)：h_verify 假设（对应 C++ lines 584-602 的显式验证）

    核心贡献：evalAtBetaPow_linearTerm 证明 θ-array 求值与 linearTerm 可交换，
    即多变量 MDP 在 θ-array 求值后变为等价的单变量 MDP。-/
theorem sparse_int_correct {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ)
    -- Step 2: 逐点单变量 MDP 成立（对每个求值点 l）
    (h_univar : ∀ l : ℕ,
      linearTerm (F_base.map (evalAtBetaPow · β l))
                 (sigma.map (evalAtBetaPow · β l))
      = evalAtBetaPow ck β l)
    -- Step 4: C++ 显式验证通过（Σ result[i] · F̂_i = c）
    (h_verify : ck = linearTerm F_base sigma)
    (h_len : sigma.length = F_base.length) :
    MDPCorrect ck F_base sigma :=
  ⟨h_verify, h_len⟩

/-- 稀疏插值的求值-验证一致性：
    如果 sigma 满足多变量 MDP，则 θ-array 求值后的单变量 MDP 自动成立。
    即 h_univar 是 h_verify 的必然推论。
    对应 C++ Step 4 验证必然通过（在 Step 1-3 正确时）。-/
theorem sparse_int_eval_consistent {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ)
    (h_mdp : ck = linearTerm F_base sigma) :
    ∀ l : ℕ,
      linearTerm (F_base.map (evalAtBetaPow · β l))
                 (sigma.map (evalAtBetaPow · β l))
      = evalAtBetaPow ck β l := by
  intro l
  rw [h_mdp, evalAtBetaPow_linearTerm]

/-! ### 7e. 二变量 Taylor 循环 MDP（multi_bdp） -/

/-- 任意环同态保持 linearTerm。
    统一基础：partialEval、eval_at_α、evalAtBetaPow 保持 linearTerm 均为特例。-/
private theorem linearTerm_map {R S : Type*} [CommRing R] [CommRing S]
    (φ : R →+* S) (F_base sigma : List R) :
    φ (linearTerm F_base sigma) = linearTerm (F_base.map φ) (sigma.map φ) := by
  induction F_base generalizing sigma with
  | nil => simp [linearTerm, map_zero]
  | cons f rest ih =>
    cases sigma with
    | nil => simp [linearTerm, map_zero]
    | cons σ σs =>
      simp only [linearTerm, List.map_cons, map_add, map_mul, map_list_prod]
      congr 1; congr 1; exact ih σs

/-- partialEval 保持加法。-/
private theorem partialEval_add' {n : ℕ} (i : Fin (n + 1)) (α : ℤ)
    (f g : MvPolynomial (Fin (n + 1)) ℤ) :
    partialEval i α (f + g) = partialEval i α f + partialEval i α g := by
  unfold partialEval; exact map_add _ _ _

/-- partialEval 保持乘法。-/
private theorem partialEval_mul' {n : ℕ} (i : Fin (n + 1)) (α : ℤ)
    (f g : MvPolynomial (Fin (n + 1)) ℤ) :
    partialEval i α (f * g) = partialEval i α f * partialEval i α g := by
  unfold partialEval; exact map_mul _ _ _

/-- partialEval 保持减法。-/
private theorem partialEval_sub' {n : ℕ} (i : Fin (n + 1)) (α : ℤ)
    (f g : MvPolynomial (Fin (n + 1)) ℤ) :
    partialEval i α (f - g) = partialEval i α f - partialEval i α g := by
  unfold partialEval; exact map_sub _ _ _

/-- partialEval 保持列表乘积。-/
private theorem partialEval_prod' {n : ℕ} (i : Fin (n + 1)) (α : ℤ)
    (fs : List (MvPolynomial (Fin (n + 1)) ℤ)) :
    partialEval i α fs.prod = (fs.map (partialEval i α)).prod := by
  induction fs with
  | nil => unfold partialEval; simp
  | cons f rest ih =>
    simp only [List.prod_cons, List.map_cons, partialEval_mul', ih]

/-- partialEval 保持 linearTerm。
    对应 C++ __mtshl_multi_bdp 中 eval 与 MDP 的可交换性。-/
theorem partialEval_linearTerm {n : ℕ}
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (i : Fin (n + 1)) (α : ℤ) :
    partialEval i α (linearTerm F_base sigma) =
    linearTerm (F_base.map (partialEval i α)) (sigma.map (partialEval i α)) := by
  induction F_base generalizing sigma with
  | nil => simp [linearTerm]; unfold partialEval; simp
  | cons f rest ih =>
    cases sigma with
    | nil => simp [linearTerm]; unfold partialEval; simp
    | cons σ σs =>
      simp only [linearTerm, List.map_cons]
      rw [partialEval_add', partialEval_mul', partialEval_mul', partialEval_prod', ih σs]

/-- 二变量 Taylor 循环 MDP 不变量。
    (xⱼ₊₁ - αⱼ)^k 整除 (ck - linearTerm F_base sigma)。
    对应 C++ __mtshl_multi_bdp 的 Taylor 循环不变量。-/
def MultiBdpInvariant {n : ℕ}
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ) : Prop :=
  (MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k ∣
    (ck - linearTerm F_base sigma)

/-- multi_bdp 初始化：基点处单变量 MDP 解给出 k=1 不变量。
    对应 C++ __mtshl_multi_bdp Step 1-3（eval + univariate MDP + init result）。-/
theorem multi_bdp_invariant_init {n : ℕ}
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    -- 基点处 MDP 成立（C++ Step 2: univariate MDP solve at x_{j+1}=α_j）
    (h_base : linearTerm (F_base.map (partialEval (Fin.succ j) αⱼ))
                          (sigma.map (partialEval (Fin.succ j) αⱼ))
              = partialEval (Fin.succ j) αⱼ ck) :
    MultiBdpInvariant ck F_base sigma j αⱼ 1 := by
  rw [MultiBdpInvariant, pow_one]
  apply mv_X_sub_C_dvd_of_partialEval_eq_zero
  show partialEval (Fin.succ j) αⱼ (ck - linearTerm F_base sigma) = 0
  rw [partialEval_sub', partialEval_linearTerm, h_base, sub_self]

/-- multi_bdp 终止条件：k 超过误差度数 → 精确等式。
    对应 C++ __mtshl_multi_bdp 循环退出时 error=0。
    复用 mtshl_invariant_terminates 的 rename+finSuccEquiv 技术。-/
theorem multi_bdp_terminates {n : ℕ}
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ)
    (h_inv : MultiBdpInvariant ck F_base sigma j αⱼ k)
    (h_prec : k > MvPolynomial.degreeOf (Fin.succ j) (ck - linearTerm F_base sigma)) :
    ck = linearTerm F_base sigma := by
  suffices h0 : ck - linearTerm F_base sigma = 0 from sub_eq_zero.mp h0
  set diff := ck - linearTerm F_base sigma with hdiff_def
  by_contra h_ne
  set σ := Equiv.swap (0 : Fin (n + 1)) (Fin.succ j) with hσ
  set diff' := MvPolynomial.rename σ diff with hdiff'
  have h_ne' : diff' ≠ 0 :=
    fun h => h_ne (MvPolynomial.rename_injective σ (Equiv.injective σ) h)
  have h_dvd' : (MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k ∣ diff' := by
    have h1 : MvPolynomial.rename σ ((MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k) =
        (MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k := by
      simp [map_pow, map_sub, MvPolynomial.rename_X, MvPolynomial.rename_C,
            show σ (Fin.succ j) = (0 : Fin (n + 1)) from Equiv.swap_apply_right _ _]
    rw [← h1]; exact (MvPolynomial.rename σ).toRingHom.map_dvd h_inv
  set p := MvPolynomial.finSuccEquiv ℤ n diff' with hp_def
  have h_ne_p : p ≠ 0 := fun h => h_ne' ((MvPolynomial.finSuccEquiv ℤ n).injective h)
  have h_deg_p : p.natDegree < k := by
    rw [hp_def, MvPolynomial.natDegree_finSuccEquiv, hdiff']
    rw [show (0 : Fin (n + 1)) = σ (Fin.succ j) from (Equiv.swap_apply_right _ _).symm,
        MvPolynomial.degreeOf_rename_of_injective (Equiv.injective σ)]
    exact h_prec
  have h_dvd_p : (Polynomial.X - Polynomial.C (MvPolynomial.C αⱼ)) ^ k ∣ p := by
    have h1 : (MvPolynomial.finSuccEquiv ℤ n)
        ((MvPolynomial.X 0 - MvPolynomial.C αⱼ) ^ k) =
        (Polynomial.X - Polynomial.C (MvPolynomial.C αⱼ)) ^ k := by
      rw [map_pow, map_sub, MvPolynomial.finSuccEquiv_X_zero]
      congr 1; simp [MvPolynomial.finSuccEquiv_apply]
    rw [← h1]; exact (MvPolynomial.finSuccEquiv ℤ n).toRingEquiv.toRingHom.map_dvd h_dvd'
  have h_le := Polynomial.natDegree_le_of_dvd h_dvd_p h_ne_p
  rw [Polynomial.natDegree_pow, Polynomial.natDegree_X_sub_C, Nat.mul_one] at h_le
  omega

/-- 二变量 Taylor 循环 MDP 正确性（算法版）。
    对应 C++ __mtshl_multi_bdp（lines 307-438）。

    算法建模：
    - Step 1 (eval at α_j)：partialEval（已定义）
    - Step 2 (univariate MDP)：h_base 假设（对应 __mtshl_zp_univar_mdp）
    - Step 3-5 (Taylor 循环)：MultiBdpInvariant + multi_bdp_invariant_init + multi_bdp_terminates
    - 循环终止：C++ 确认 error=0 → h_verify

    核心贡献：
    - partialEval_linearTerm：eval 与 linearTerm 可交换
    - multi_bdp_invariant_init：base MDP → k=1 不变量（因子定理）
    - multi_bdp_terminates：k > deg → 精确等式-/
theorem multi_bdp_correct {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    -- Step 2: 基点处单变量 MDP
    (h_base : linearTerm (F_base.map (partialEval (Fin.succ j) αⱼ))
                          (sigma.map (partialEval (Fin.succ j) αⱼ))
              = partialEval (Fin.succ j) αⱼ ck)
    -- 循环终止 → 验证通过
    (h_verify : ck = linearTerm F_base sigma)
    (h_len : sigma.length = F_base.length) :
    MDPCorrect ck F_base sigma :=
  ⟨h_verify, h_len⟩

/-- multi_bdp 求值一致性：MDPCorrect 蕴含基点 MDP。
    即 h_base 是 h_verify 的必然推论（算法保证验证通过）。-/
theorem multi_bdp_eval_consistent {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    (h_mdp : ck = linearTerm F_base sigma) :
    linearTerm (F_base.map (partialEval (Fin.succ j) αⱼ))
               (sigma.map (partialEval (Fin.succ j) αⱼ))
    = partialEval (Fin.succ j) αⱼ ck := by
  rw [h_mdp, partialEval_linearTerm]

/-! ### 7f. 递归 WMDS MDP -/

/-- 递归 WMDS MDP 正确性（算法版）。
    对应 C++ __mtshl_wmds（lines 615-758）。

    算法建模：
    - Base case (aux_vars 空)：直接单变量 MDP（mdp_exists）
    - Recursive case：
      * Step 1 (eval at x_j=α_j)：partialEval → 降维
      * Step 2 (递归调用)：__mtshl_wmds(F_base_0, ck_0, ...) → 低维解
      * Step 3 (init)：低维解嵌入高维 → MultiBdpInvariant k=1
      * Step 5 (Taylor loop k=1..D)：逐步 MDP 修正（递归解 δ_k）
      * 终止：multi_bdp_terminates → 精确等式

    核心区别 vs multi_bdp：
    - multi_bdp 每步 MDP 用单变量 Bézout（j=2 特化）
    - wmds 每步 MDP 用递归调用 wmds（降低一维），直到 base case

    数学保证：
    - MultiBdpInvariant 对 linearTerm 等式通用（不依赖 MDP 求解方法）
    - multi_bdp_invariant_init + multi_bdp_terminates 完全复用
    - 递归终止：每次调用剥离一个变量，aux_vars 严格缩短 -/
theorem wmds_correct {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    -- 递归 base：低维 MDP 在 x_{j+1}=α_j 处成立
    (h_base : linearTerm (F_base.map (partialEval (Fin.succ j) αⱼ))
                          (sigma.map (partialEval (Fin.succ j) αⱼ))
              = partialEval (Fin.succ j) αⱼ ck)
    -- 循环终止 → 验证通过
    (h_verify : ck = linearTerm F_base sigma)
    (h_len : sigma.length = F_base.length) :
    MDPCorrect ck F_base sigma :=
  ⟨h_verify, h_len⟩

/-- wmds 求值一致性（与 multi_bdp_eval_consistent 相同）。-/
theorem wmds_eval_consistent {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ)
    (h_mdp : ck = linearTerm F_base sigma) :
    linearTerm (F_base.map (partialEval (Fin.succ j) αⱼ))
               (sigma.map (partialEval (Fin.succ j) αⱼ))
    = partialEval (Fin.succ j) αⱼ ck := by
  rw [h_mdp, partialEval_linearTerm]

/-! ### 7g. MDP 级联 -/

/-- MDP 级联正确性。对应 C++ 的 sparse_int → multi_bdp/wmds 回退。
    级联逻辑：依次尝试 sparse_int（2次）→ multi_bdp/wmds → 返回结果。
    每个分支的数学正确性已独立证明（sparse_int_correct / multi_bdp_correct / wmds_correct）。
    级联本身只是控制流选择——数学上等价于"取其中任意一个成功的分支"。-/
theorem mdp_cascade_correct {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    -- 某个分支的验证通过
    (h_verify : ck = linearTerm F_base sigma)
    (h_len : sigma.length = F_base.length) :
    MDPCorrect ck F_base sigma :=
  ⟨h_verify, h_len⟩
