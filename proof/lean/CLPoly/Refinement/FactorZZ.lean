/-
  Strict end-to-end refinement assembly for the original C++
  `__factor_squarefree_primitive_ZZ` entry.

  This module is intentionally downstream of both SelectPrime and Hensel.  It
  contains cross-stage composition facts; neither component module imports the
  other or assumes a semantic result produced by a later stage.
-/
import CLPoly.Generated.StrictFactorZZ
import CLPoly.Refinement.Hensel
import CLPoly.Refinement.Recombine
import CLPoly.Refinement.SelectPrime
import Mathlib.Data.ZMod.Coprime

set_option autoImplicit false

namespace Refinement

namespace StrictFactorZZ

/-- The actual generated select-prime entry with its concrete modular
candidate callback and machine-prime iterator.  The RNG state is the local
source state supplied at entry (the C++ caller uses its fixed seed here). -/
noncomputable def concreteSelectPrime {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) :
    SparsePolyZZ → Bool → RawExec PrimeSelectionResult :=
  fun f useLargePrime =>
    Generated.StrictSelectPrime.__select_prime_raw_ir
      (StrictSelectPrime.selectPrimeRawOps
        (StrictSelectPrime.concreteTryCandidate engine provider))
      (StrictSelectPrime.selectPrimeTermination
        (StrictSelectPrime.concreteTryCandidate engine provider))
      initialRng useLargePrime f

theorem concreteSelectPrime_success {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegreeBound : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hlcSemantic : ∀ p : UInt64, Nat.Prime p.toNat →
      ((SparsePolyZZ.front! f).2 : ZMod p.toNat) =
        ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat))
    (selection : PrimeSelectionResult)
    (hrun : concreteSelectPrime engine provider initialRng f useLargePrime =
      .ok selection) :
    StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection ∧
      StrictSelectPrime.SelectionPhysical selection := by
  exact StrictSelectPrime.__select_prime_raw_ir_refines engine provider
    initialRng useLargePrime f hinitialPrimeCorrect hcanonical hnonempty hdegree hdegreeBound
      hlcSemantic selection hrun

theorem concreteSelectPrime_irreducible_size {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (selection : PrimeSelectionResult)
    (hrun : concreteSelectPrime engine provider initialRng f useLargePrime =
      .ok selection) :
    selection.irreducible = true → selection.factors.size ≤ 1 := by
  unfold concreteSelectPrime at hrun
  unfold Generated.StrictSelectPrime.__select_prime_raw_ir at hrun
  split at hrun
  · contradiction
  next hguard =>
    exact StrictSelectPrime.selectPrimeLoop_irreducible_size
      (StrictSelectPrime.concreteTryCandidate engine provider)
      f (get_deg f) (SparsePolyZZ.front! f).2 useLargePrime 3
      { tried := 0
        p := if useLargePrime then
          (18446744073709551615 : UInt64) - 58 else 2
        rng := initialRng
        bestCount := 18446744073709551615
        best := default }
      selection hrun

/-- Actual generated Hensel entry selected at the runtime prime returned by
`__select_prime`.  The dense arithmetic object and multiplication workspace
come from the same prime-indexed provider used by the modular factorization;
no Hensel output is supplied by the caller. -/
noncomputable def concreteHenselLift {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine) :
    SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
      RawExec (Array SparsePolyZZ × ZZ) :=
  fun f factors p aTarget =>
    if hp : Nat.Prime p.toNat then
      let candidate := provider.physical p hp
      letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
        (StrictHensel.strictHenselRawOps
          StrictHensel.concreteDivmodTermination)
        (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
          candidate.providers.mul)
        f factors p aTarget candidate.prime.two_le
    else .error .assertionFailure

/-- Readiness for the semantic Hensel stages.  Machine factor-count and tree
index bounds are derived separately from the successful literal execution. -/
def HenselLiftEntryReadiness
    (this : DenseUPolyZp)
    (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp)
    (aTarget : Int32) : Prop :=
  StrictHensel.HenselLiftEntryInvariant this termination mulProvider
    f factors aTarget

/-- The concrete C++ Hensel call cannot succeed unless its checked `size_t`
factor count is representable by every downstream source `int` index. -/
theorem concreteHenselLift_factorCountFits_of_success {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64)
    (aTarget : Int32) (hp : Nat.Prime p.toNat)
    (output : Array SparsePolyZZ × ZZ)
    (hrun : concreteHenselLift engine provider f factors p aTarget =
      .ok output) :
    factors.size < 2 ^ 31 := by
  let candidate := provider.physical p hp
  letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
  have hrun' :
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
        (StrictHensel.strictHenselRawOps
          StrictHensel.concreteDivmodTermination)
        (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
          candidate.providers.mul)
        f factors p aTarget candidate.prime.two_le = .ok output := by
    simpa [concreteHenselLift, hp, candidate] using hrun
  exact StrictHensel.__hensel_lift_upoly_raw_ir_factorCountFits_of_success
    (StrictHensel.strictHenselRawOps
      StrictHensel.concreteDivmodTermination)
    (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
      candidate.providers.mul)
    f factors p aTarget candidate.prime.two_le output hrun'

/-- A successful call of `concreteHenselLift` is the literal generated
Hensel execution and therefore carries `HenselLiftEntryCorrect`.  The
remaining invariant is stage readiness quantified over actual intermediate
executions; it cannot prescribe the returned factor array. -/
theorem concreteHenselLift_success {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64)
    (aTarget : Int32) (hp : Nat.Prime p.toNat)
    (hreadiness :
      let candidate := provider.physical p hp
      @HenselLiftEntryReadiness candidate.dense
        StrictHensel.concreteDivmodTermination
        candidate.providers.mul f factors aTarget)
    (output : Array SparsePolyZZ × ZZ)
    (hrun : concreteHenselLift engine provider f factors p aTarget =
      .ok output) :
    @StrictHensel.HenselLiftEntryCorrect
      StrictHensel.concreteDivmodTermination f factors p ⟨hp⟩ aTarget
      output := by
  let candidate := provider.physical p hp
  letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
  have hrun' :
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
        (StrictHensel.strictHenselRawOps
          StrictHensel.concreteDivmodTermination)
        (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
          candidate.providers.mul)
        f factors p aTarget candidate.prime.two_le = .ok output := by
    simpa [concreteHenselLift, hp, candidate] using hrun
  have hfactorFits : factors.size < 2 ^ 31 :=
    concreteHenselLift_factorCountFits_of_success engine provider f factors
      p aTarget hp output hrun
  have hfactorCount : 2 ≤ factors.size :=
    StrictHensel.__hensel_lift_upoly_raw_ir_factorCount_of_success
      (StrictHensel.strictHenselRawOps
        StrictHensel.concreteDivmodTermination)
      (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
        candidate.providers.mul)
      f factors p aTarget candidate.prime.two_le output hrun'
  rcases StrictHensel.__hensel_lift_upoly_raw_ir_refines candidate.dense
      candidate.configured
      StrictHensel.concreteDivmodTermination candidate.providers.mul
      f factors aTarget hreadiness hfactorCount hfactorFits with
        ⟨actual, hactualRun, hcorrect⟩
  have houtput : actual = output := by
    have hactualRun' :
        Generated.StrictHensel.__hensel_lift_upoly_raw_ir
          (StrictHensel.strictHenselRawOps
            StrictHensel.concreteDivmodTermination)
          (StrictHensel.strictHenselTreeBuildRawOps candidate.dense
            candidate.providers.mul)
          f factors p aTarget candidate.prime.two_le = .ok actual := by
      simpa [candidate] using hactualRun
    rw [hactualRun'] at hrun'
    exact Except.ok.inj hrun'
  subst output
  simpa [candidate] using hcorrect

/-- Concrete adapter for the generated C++ van-Hoeij callee.  The proof-only
degree guard exposes the source precondition required by the strict raw entry;
on every valid `__lll_factorize` call it erases to that exact execution. -/
def concreteVanHoeijRecombine (f : SparsePolyZZ)
    (lifted : Array SparsePolyZZ) (modulus : ZZ) :
    RawExec (Array SparsePolyZZ) :=
  if hdegree : 2 ≤ (get_deg f).toNatClampNeg then
    Generated.StrictRecombine.__vanhoeij_recombine_raw_ir
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination f lifted modulus hdegree
  else .error .assertionFailure

def concreteZassenhausRecombine (f : SparsePolyZZ)
    (lifted : Array SparsePolyZZ) (modulus : ZZ) :
    RawExec (Array SparsePolyZZ) :=
  Generated.StrictRecombine.zassenhausRecombine
    StrictRecombine.concreteZassenhausTermination f lifted modulus

/-- Concrete recombination fields plus the still-explicit generated
selection/Hensel callees.  Those two functions return only raw data; their
execution certificates are supplied separately and cannot choose results. -/
def concreteRecombineFactorZZRawOps
    (selectPrime : SparsePolyZZ → Bool → RawExec PrimeSelectionResult)
    (henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
      RawExec (Array SparsePolyZZ × ZZ)) :
    Generated.StrictFactorZZ.FactorZZRawOps where
  selectPrime := selectPrime
  henselLift := henselLift
  vanHoeijRecombine := concreteVanHoeijRecombine
  zassenhausRecombine := concreteZassenhausRecombine

theorem concreteVanHoeijRecombine_success
    (f : SparsePolyZZ) (lifted output : Array SparsePolyZZ) (modulus : ZZ)
    (hrun : concreteVanHoeijRecombine f lifted modulus = .ok output) :
    ∃ hdegree : 2 ≤ (get_deg f).toNatClampNeg,
      Generated.StrictRecombine.__vanhoeij_recombine_raw_ir
        StrictRecombine.concreteVanHoeijRawOps
        StrictRecombine.concreteVanHoeijTermination f lifted modulus hdegree =
          .ok output := by
  unfold concreteVanHoeijRecombine at hrun
  split at hrun
  next hdegree => exact ⟨hdegree, hrun⟩
  next hdegree => contradiction

/-- Below the signed 32-bit boundary, the exact C++ size conversion used by
`__lll_factorize` reflects the natural-number strict comparison. -/
theorem size_toUInt32_toInt32_lt_iff (left right : Nat)
    (hleft : left < 2 ^ 31) (hright : right < 2 ^ 31) :
    left.toUInt32.toInt32 < right.toUInt32.toInt32 ↔ left < right := by
  rw [Int32.lt_iff_toInt_lt]
  change (Int32.ofNat left).toInt < (Int32.ofNat right).toInt ↔ _
  rw [Int32.toInt_ofNat_of_lt hleft, Int32.toInt_ofNat_of_lt hright]
  omega

/-- If the concrete recombination size is already known not to exceed the
Hensel cardinality, failure of the literal signed C++ "fewer factors" test
forces exact cardinality equality. -/
theorem result_size_eq_of_not_machine_lt
    (result : Array SparsePolyZZ) (factorCount : Nat)
    (hfactorFits : factorCount < 2 ^ 31)
    (hresultLe : result.size ≤ factorCount)
    (hnot : ¬ result.size.toUInt32.toInt32 <
      factorCount.toUInt32.toInt32) :
    result.size = factorCount := by
  have hresultFits : result.size < 2 ^ 31 := lt_of_le_of_lt hresultLe hfactorFits
  rw [size_toUInt32_toInt32_lt_iff result.size factorCount hresultFits
    hfactorFits] at hnot
  omega

/-- Exact control-flow classification of the low-precision branch of the
generated C++ `__lll_factorize`.  If the first recombination is accepted, its
physical result has the full Hensel cardinality.  Otherwise the returned array
is produced by the literal second Hensel call with source argument zero and
the immediately following recombination call. -/
theorem __lll_factorize_raw_ir_low_precision_cases
    (ops : Generated.StrictFactorZZ.FactorZZRawOps)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64)
    (aH aMig : Int32) (liftedH : Array SparsePolyZZ) (mH : ZZ)
    (result output : Array SparsePolyZZ)
    (hheuristic :
      Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir f
        factors.size.toUInt32.toInt32 p = .ok (aH, aMig))
    (hlift : ops.henselLift f factors p aH = .ok (liftedH, mH))
    (hrecombine : ops.vanHoeijRecombine f liftedH mH = .ok result)
    (hliftedSize : liftedH.size = factors.size)
    (hresultLe : result.size ≤ liftedH.size)
    (hfactorFits : factors.size < 2 ^ 31)
    (hprecision : aH < aMig)
    (hrun : Generated.StrictFactorZZ.__lll_factorize_raw_ir ops f factors p =
      .ok output) :
    (output = result ∧ result.size = factors.size) ∨
      ∃ liftedMig mMig resultMig,
        ops.henselLift f factors p 0 = .ok (liftedMig, mMig) ∧
        ops.vanHoeijRecombine f liftedMig mMig = .ok resultMig ∧
        (output = resultMig ∨
          ops.zassenhausRecombine f liftedMig mMig = .ok output) := by
  unfold Generated.StrictFactorZZ.__lll_factorize_raw_ir at hrun
  dsimp only at hrun
  simp only [hheuristic, hlift, hrecombine] at hrun
  by_cases hless : Int32.ofNat result.size < Int32.ofNat factors.size
  · rw [if_pos (by simp [hless, hprecision])] at hrun
    split at hrun
    next fault hliftMig => contradiction
    next liftedMig mMig hliftMig =>
      cases hrecombineMig : ops.vanHoeijRecombine f liftedMig mMig with
      | error fault => rw [hrecombineMig] at hrun; contradiction
      | ok resultMig =>
          rw [hrecombineMig] at hrun
          change (if Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir
              resultMig.size factors.size true then
                ops.zassenhausRecombine f liftedMig mMig
              else .ok resultMig) = .ok output at hrun
          by_cases hsafety :
              Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir
                resultMig.size factors.size true = true
          · rw [if_pos hsafety] at hrun
            exact Or.inr ⟨liftedMig, mMig, resultMig, hliftMig,
              hrecombineMig, Or.inr hrun⟩
          · rw [if_neg hsafety] at hrun
            have hout := Except.ok.inj hrun
            subst output
            exact Or.inr ⟨liftedMig, mMig, resultMig, hliftMig,
              hrecombineMig, Or.inl rfl⟩
  · rw [if_neg (by simp [hless])] at hrun
    have hsafety : Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir
        result.size factors.size (aMig ≤ aH) = false := by
      simp [Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir,
        Int32.not_le.mpr hprecision]
    rw [if_neg (by simpa using hsafety)] at hrun
    have hout := Except.ok.inj hrun
    subst output
    exact Or.inl ⟨rfl, result_size_eq_of_not_machine_lt result factors.size
      hfactorFits (by omega) (by simpa using hless)⟩

/-- The first precision returned by the exact generated heuristic is always
bounded by its full Mignotte precision.  This follows from the literal final
`min`; fault branches cannot produce a successful pair. -/
theorem heuristic_starting_precision_first_le_second
    (f : SparsePolyZZ) (r : Int32) (p : UInt64) (aH aMig : Int32)
    (hrun : Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir
      f r p = .ok (aH, aMig)) :
    aH ≤ aMig := by
  unfold Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir at hrun
  split at hrun
  next hp =>
    split at hrun
    next => contradiction
    next leading hleading =>
      split at hrun
      next fault hbound => contradiction
      next bound hbound =>
        dsimp only at hrun
        split at hrun
        next hsign =>
          split at hrun
          next hfit =>
            have hout := Except.ok.inj hrun
            cases hout
            change (if _ ≤ _ then _ else _) ≤ _
            split
            next => rw [Int32.le_iff_toInt_le]
            next hnot =>
              rw [Int32.le_iff_toInt_le]
              by_contra hreverse
              apply hnot
              rw [Int32.le_iff_toInt_le]
              omega
          next hfit => contradiction
        next hsign =>
          split at hrun
          next hfit =>
            have hout := Except.ok.inj hrun
            cases hout
            change (if _ ≤ _ then _ else _) ≤ _
            split
            next => rw [Int32.le_iff_toInt_le]
            next hnot =>
              rw [Int32.le_iff_toInt_le]
              by_contra hreverse
              apply hnot
              rw [Int32.le_iff_toInt_le]
              omega
          next hfit => contradiction
  next hp => contradiction

/-- The exact well-founded precision loop stops only after the represented
prime power has crossed the concrete recovery target.  The invariant
`pa = p ^ exponent` is the literal arithmetic state of the C++ loop. -/
theorem heuristicPrecisionLoop_pow_gt
    (p target : Nat) (hp : 2 ≤ p) (pa : Nat) (hpa : 0 < pa)
    (exponent : Nat) (hpower : pa = p ^ exponent) :
    target < p ^
      Generated.StrictFactorZZ.heuristicPrecisionLoop p target hp pa hpa
        exponent := by
  rw [Generated.StrictFactorZZ.heuristicPrecisionLoop]
  split
  next hcontinue =>
    apply heuristicPrecisionLoop_pow_gt p target hp (pa * p)
      (Nat.mul_pos hpa (by omega)) (exponent + 1)
    rw [Nat.pow_succ, ← hpower, Nat.mul_comm]
  next hstop =>
    rw [← hpower]
    omega
termination_by target + 1 - pa
decreasing_by
  have hdouble : pa + pa ≤ pa * p := by
    calc
      pa + pa = 2 * pa := by omega
      _ ≤ p * pa := Nat.mul_le_mul_right pa hp
      _ = pa * p := Nat.mul_comm _ _
  have hgrow : pa < pa * p := by omega
  omega

theorem nat_toUInt32_toInt32_nonnegative_implies_fits (value : Nat)
    (hfit : value ≤ UInt32.size - 1)
    (hnonnegative : 0 ≤ value.toUInt32.toInt32) : value < 2 ^ 31 := by
  have h32 : value < 2 ^ 32 := by
    norm_num [UInt32.size] at hfit ⊢
    omega
  rw [Int32.le_iff_toInt_le] at hnonnegative
  simp only [Int32.toInt_zero] at hnonnegative
  change 0 ≤ (BitVec.ofNat 32 value).toInt at hnonnegative
  rw [BitVec.toInt_eq_toNat_cond] at hnonnegative
  have hnat : (BitVec.ofNat 32 value).toNat = value := by
    rw [BitVec.toNat_ofNat, Nat.mod_eq_of_lt h32]
  split at hnonnegative
  · rw [hnat] at *
    omega
  · rw [hnat] at *
    norm_num at hnonnegative
    omega

/-- A successful generated heuristic run exposes its exact full exponent and
the fact that the corresponding prime power crosses the concrete Mignotte
target computed by the same run. -/
theorem heuristic_starting_precision_full_pow_gt
    (f : SparsePolyZZ) (r : Int32) (p : UInt64) (aH aMig : Int32)
    (leading : UMonomial × ZZ) (bound : ZZ)
    (hleading : f[0]? = some leading)
    (hbound : Generated.StrictHensel.__mignotte_bound_upoly_raw_ir f =
      .ok bound)
    (hrun : Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir
      f r p = .ok (aH, aMig))
    (hnonnegative : 0 ≤ aMig) :
    (2 * (if leading.2 < 0 then -leading.2 else leading.2) * bound).natAbs <
      p.toNat ^ aMig.toNatClampNeg := by
  unfold Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir at hrun
  split at hrun
  next hp =>
    rw [hleading, hbound] at hrun
    dsimp only at hrun
    split at hrun
    next hsign =>
      split at hrun
      next hfit =>
        have hout := Except.ok.inj hrun
        cases hout
        have hfullNonnegative :
            0 ≤ (Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
              (2 * -leading.2 * bound).natAbs hp 1 (by omega) 0).toUInt32.toInt32 :=
          hnonnegative
        have hfullFits := nat_toUInt32_toInt32_nonnegative_implies_fits _ hfit
          hfullNonnegative
        have hroundtrip :=
          StrictRecombine.nat_toUInt32_toInt32_nonnegative_and_toNat _ hfullFits
        have hexponent :
            (Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
              (2 * -leading.2 * bound).natAbs hp 1 (by omega) 0).toUInt32.toInt32.toNatClampNeg =
              Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
                (2 * -leading.2 * bound).natAbs hp 1 (by omega) 0 := by
          rw [← Int32.toNatClampNeg_toInt64]
          exact hroundtrip.2
        rw [show (2 * (if leading.2 < 0 then -leading.2 else leading.2) *
            bound).natAbs = (2 * -leading.2 * bound).natAbs by simp [hsign]]
        rw [hexponent]
        exact heuristicPrecisionLoop_pow_gt p.toNat
          (2 * -leading.2 * bound).natAbs hp 1 (by omega) 0 (by simp)
      next => contradiction
    next hsign =>
      split at hrun
      next hfit =>
        have hout := Except.ok.inj hrun
        cases hout
        have hfullNonnegative :
            0 ≤ (Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
              (2 * leading.2 * bound).natAbs hp 1 (by omega) 0).toUInt32.toInt32 :=
          hnonnegative
        have hfullFits := nat_toUInt32_toInt32_nonnegative_implies_fits _ hfit
          hfullNonnegative
        have hroundtrip :=
          StrictRecombine.nat_toUInt32_toInt32_nonnegative_and_toNat _ hfullFits
        have hexponent :
            (Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
              (2 * leading.2 * bound).natAbs hp 1 (by omega) 0).toUInt32.toInt32.toNatClampNeg =
              Generated.StrictFactorZZ.heuristicPrecisionLoop p.toNat
                (2 * leading.2 * bound).natAbs hp 1 (by omega) 0 := by
          rw [← Int32.toNatClampNeg_toInt64]
          exact hroundtrip.2
        rw [show (2 * (if leading.2 < 0 then -leading.2 else leading.2) *
            bound).natAbs = (2 * leading.2 * bound).natAbs by simp [hsign]]
        rw [hexponent]
        exact heuristicPrecisionLoop_pow_gt p.toNat
          (2 * leading.2 * bound).natAbs hp 1 (by omega) 0 (by simp)
      next => contradiction
  next => contradiction

open CLPoly.Math

private theorem polynomial_eq_C_leadingCoeff_mul_of_associated_monic
    {K : Type*} [Field K] {f g : Polynomial K}
    (hf : f ≠ 0) (hg : g.Monic) (hassociated : Associated f g) :
    f = Polynomial.C f.leadingCoeff * g := by
  rcases hassociated with ⟨unit, hunit⟩
  rcases Polynomial.isUnit_iff.mp unit.isUnit with
    ⟨coefficient, _hcoefficientUnit, hcoefficient⟩
  have hproduct : f * Polynomial.C coefficient = g := by
    rw [hcoefficient]
    exact hunit
  have hleading := congrArg Polynomial.leadingCoeff hproduct
  rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
    hg.leadingCoeff] at hleading
  symm
  rw [← hproduct]
  calc
    Polynomial.C f.leadingCoeff *
          (f * Polynomial.C coefficient) =
        f * Polynomial.C (f.leadingCoeff * coefficient) := by
      rw [Polynomial.C_mul]
      ring
    _ = f := by rw [hleading]; simp

/-- Resolve the unit ambiguity of a modular subproduct using the leading
coefficient that the concrete Zassenhaus trial places in front of it.

If `source = divisor * quotient`, the monic modular representative associated
with `divisor` becomes exactly `C quotient.leadingCoeff * divisor` after it is
scaled by `source.leadingCoeff`.  This is the normalization needed by the
generated trial product; keeping the quotient's leading coefficient here is
essential for non-monic integer factors. -/
theorem leading_scaled_monic_associated_divisor
    (p : Nat) [Fact (Nat.Prime p)]
    (source divisor quotient : Polynomial Int)
    (candidate : Polynomial (ZMod p))
    (hfactor : source = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod p)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map (Int.castRingHom (ZMod p)) divisor).leadingCoeff =
        (divisor.leadingCoeff : ZMod p))
    (hcandidateMonic : candidate.Monic)
    (hassociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod p)) divisor) candidate) :
    Polynomial.C (source.leadingCoeff : ZMod p) * candidate =
      Polynomial.map (Int.castRingHom (ZMod p))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  have hexact := polynomial_eq_C_leadingCoeff_mul_of_associated_monic
    hdivisorModNonzero hcandidateMonic hassociated
  rw [hdivisorLeading] at hexact
  have hsourceLeading :
      source.leadingCoeff = divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [hfactor, Polynomial.leadingCoeff_mul]
  rw [hsourceLeading, Int.cast_mul, Polynomial.C_mul,
    Polynomial.map_mul, Polynomial.map_C]
  rw [hexact]
  change Polynomial.C (divisor.leadingCoeff : ZMod p) *
      Polynomial.C (quotient.leadingCoeff : ZMod p) * candidate =
    Polynomial.C (quotient.leadingCoeff : ZMod p) *
      (Polynomial.C (divisor.leadingCoeff : ZMod p) * candidate)
  ring

private theorem monic_list_product {K : Type*} [Semiring K]
    (factors : List (Polynomial K))
    (hmonic : ∀ factor ∈ factors, factor.Monic) : factors.prod.Monic := by
  induction factors with
  | nil => exact Polynomial.monic_one
  | cons factor tail ih =>
      exact (hmonic factor (by simp)).mul
        (ih fun candidate hcandidate => hmonic candidate (by simp [hcandidate]))

/-- Every divisor of a concrete finite product of irreducibles is associated
to the product of an actual sublist.  The proof recursively follows the
source-order list: if the head atom divides the divisor it is cancelled from
both sides; otherwise irreducibility makes it coprime to the divisor and it
is cancelled from the ambient product. -/
theorem divisor_associated_sublist_product
    {R : Type*} [CommRing R] [IsDomain R] [IsBezout R]
    (atoms : List R) (hirreducible : ∀ atom ∈ atoms, Irreducible atom)
    (divisor : R) (hdivisor : divisor ∣ atoms.prod) :
    ∃ chosen : List R, List.Sublist chosen atoms ∧
      Associated divisor chosen.prod := by
  induction atoms generalizing divisor with
  | nil =>
      refine ⟨[], List.Sublist.refl [], associated_one_iff_isUnit.mpr ?_⟩
      simpa using isUnit_of_dvd_one hdivisor
  | cons atom atoms ih =>
      have hatom := hirreducible atom (by simp)
      have htail : ∀ candidate ∈ atoms, Irreducible candidate := by
        intro candidate hcandidate
        exact hirreducible candidate (by simp [hcandidate])
      by_cases hheadDivisor : atom ∣ divisor
      · rcases hheadDivisor with ⟨quotient, hquotient⟩
        subst divisor
        have hquotientDivides : quotient ∣ atoms.prod := by
          exact (mul_dvd_mul_iff_left hatom.ne_zero).mp (by
            simpa [List.prod_cons] using hdivisor)
        rcases ih htail quotient hquotientDivides with
          ⟨chosen, hchosen, hassociated⟩
        refine ⟨atom :: chosen, hchosen.cons_cons atom, ?_⟩
        simpa [List.prod_cons] using
          Associated.mul_mul (Associated.refl atom) hassociated
      · have hcoprime : IsCoprime divisor atom :=
          (hatom.coprime_iff_not_dvd.mpr hheadDivisor).symm
        have hdivisorTail : divisor ∣ atoms.prod := by
          exact hcoprime.dvd_of_dvd_mul_left (by
            simpa [List.prod_cons] using hdivisor)
        rcases ih htail divisor hdivisorTail with
          ⟨chosen, hchosen, hassociated⟩
        exact ⟨chosen, hchosen.cons atom, hassociated⟩

/-- If both sides of a factorization of a concrete irreducible product are
nonunits, the occurrence-sensitive sublist representing the left factor is
proper.  Hence its cardinality is strictly smaller than the full atom list.
This is the counting fact used by the Zassenhaus minimal-subset argument. -/
theorem divisor_associated_proper_sublist_product
    {R : Type*} [CommRing R] [IsDomain R] [IsBezout R]
    (atoms : List R) (hirreducible : ∀ atom ∈ atoms, Irreducible atom)
    (left right : R) (hfactor : Associated (left * right) atoms.prod)
    (hleftNonunit : ¬IsUnit left) (hrightNonunit : ¬IsUnit right) :
    ∃ chosen : List R, List.Sublist chosen atoms ∧
      Associated left chosen.prod ∧ chosen.length < atoms.length := by
  have hatomsNe : atoms.prod ≠ 0 := by
    apply List.prod_ne_zero
    intro hzero
    exact (hirreducible 0 hzero).ne_zero rfl
  have hleftNe : left ≠ 0 := by
    intro hzero
    apply hatomsNe
    exact hfactor.eq_zero_iff.mp (by simp [hzero])
  have hleftDvd : left ∣ atoms.prod := by
    exact dvd_trans (dvd_mul_right left right) hfactor.dvd
  rcases divisor_associated_sublist_product atoms hirreducible left
      hleftDvd with ⟨chosen, hchosen, hassociated⟩
  refine ⟨chosen, hchosen, hassociated, ?_⟩
  have hlengthLe := List.Sublist.length_le hchosen
  apply lt_of_le_of_ne hlengthLe
  intro hlength
  have hchosenEq : chosen = atoms :=
    List.Sublist.eq_of_length hchosen hlength
  have hleftAtoms : Associated left atoms.prod := by
    simpa [hchosenEq] using hassociated
  have hcancel : Associated (left * 1) (left * right) := by
    simpa using hleftAtoms.trans hfactor.symm
  have honeRight : Associated (1 : R) right :=
    Associated.of_mul_left hcancel (Associated.refl left) hleftNe
  exact hrightNonunit (associated_one_iff_isUnit.mp honeRight.symm)

/-- Array-level occurrence form of the proper-sublist theorem.  A nontrivial
factorization of the modular image of an extracted candidate yields a legal
candidate over the same physical active array whose size is strictly smaller
than the extracted candidate's size. -/
theorem smaller_active_candidate_of_reducible_selected_product
    (base : Nat) [Fact (Nat.Prime base)]
    (active : Array SparsePolyZZ) (outer : Array Nat)
    (houter : StrictRecombine.LegalCombination active.size outer.size outer)
    (hirreducible : ∀ index (hindex : index < active.size),
      Irreducible (StrictHensel.toPolyMod base active[index]))
    (factor left right : Polynomial (ZMod base))
    (hfactor : factor = left * right)
    (hleftNonunit : ¬IsUnit left) (hrightNonunit : ¬IsUnit right)
    (hassociated : Associated factor
      (((StrictRecombine.selectSourceIndices active.toList outer.toList).map
        (StrictHensel.toPolyMod base)).prod)) :
    ∃ inner : Array Nat,
      StrictRecombine.LegalCombination active.size inner.size inner ∧
      inner.size < outer.size ∧
      Associated left
        (((StrictRecombine.selectSourceIndices active.toList inner.toList).map
          (StrictHensel.toPolyMod base)).prod) := by
  let selected := StrictRecombine.selectSourceIndices active.toList outer.toList
  let atoms := selected.map (StrictHensel.toPolyMod base)
  have hselectedSublist : selected.Sublist active.toList :=
    StrictRecombine.selectSourceIndices_sublist active.toList outer
      (by simpa using houter)
  have hirreducibleAtoms : ∀ atom ∈ atoms, Irreducible atom := by
    intro atom hatom
    rcases List.mem_map.mp hatom with ⟨lifted, hlifted, rfl⟩
    have hliftedActive : lifted ∈ active.toList :=
      hselectedSublist.subset hlifted
    rcases List.mem_iff_getElem.mp hliftedActive with
      ⟨index, hindex, rfl⟩
    have hindexArray : index < active.size := by simpa using hindex
    simpa [Array.getElem_toList] using hirreducible index hindexArray
  have hfactorAtoms : Associated (left * right) atoms.prod :=
    (Associated.of_eq hfactor.symm).trans (by simpa [atoms, selected] using
      hassociated)
  rcases divisor_associated_proper_sublist_product atoms hirreducibleAtoms
      left right hfactorAtoms hleftNonunit hrightNonunit with
    ⟨chosenMod, hchosenMod, hleftChosen, hlength⟩
  rcases List.sublist_map_iff.mp hchosenMod with
    ⟨chosen, hchosenSelected, rfl⟩
  have hchosenActive : chosen.Sublist active.toList :=
    hchosenSelected.trans hselectedSublist
  rcases StrictRecombine.sublist_exists_legal_combination hchosenActive with
    ⟨inner, hinner, hselect⟩
  have hinner' : StrictRecombine.LegalCombination active.size inner.size
      inner := by simpa [hinner.1] using hinner
  refine ⟨inner, hinner', ?_, ?_⟩
  · simpa [atoms, selected, hinner.1,
      StrictRecombine.selectSourceIndices] using hlength
  · rw [hselect]
    exact hleftChosen

/-- A primitive integer polynomial is irreducible when its reduction at a
prime is irreducible and the leading coefficient survives reduction.  This is
the non-monic form needed for the actual primitive factors emitted by C++
recombination.  The leading-coefficient hypothesis prevents a factor from
becoming a constant merely through degree loss modulo `p`. -/
theorem primitive_irreducible_of_irreducible_mod
    (p : Nat) [Fact (Nat.Prime p)] (g : Polynomial Int)
    (hprimitive : g.IsPrimitive)
    (hleading : (g.leadingCoeff : ZMod p) ≠ 0)
    (hmodIrreducible :
      Irreducible (g.map (Int.castRingHom (ZMod p)))) :
    Irreducible g := by
  refine ⟨fun hunit => hmodIrreducible.not_isUnit
      (hunit.map (Polynomial.mapRingHom (Int.castRingHom (ZMod p)))), ?_⟩
  intro left right hfactor
  have hleadingProduct :
      (left.leadingCoeff : ZMod p) * (right.leadingCoeff : ZMod p) =
        (g.leadingCoeff : ZMod p) := by
    rw [← Int.cast_mul]
    congr 1
    rw [← Polynomial.leadingCoeff_mul, hfactor]
  have hleftLeading : (left.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, zero_mul] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hrightLeading : (right.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hmapped :
      g.map (Int.castRingHom (ZMod p)) =
        left.map (Int.castRingHom (ZMod p)) *
          right.map (Int.castRingHom (ZMod p)) := by
    rw [← Polynomial.map_mul, hfactor]
  rcases hmodIrreducible.isUnit_or_isUnit hmapped with hleft | hright
  · have hleftPrimitive : left.IsPrimitive :=
      Polynomial.isPrimitive_of_dvd hprimitive
        (Dvd.intro right hfactor.symm)
    have hdegreeMapped :
        (left.map (Int.castRingHom (ZMod p))).degree = left.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hleftLeading
    have hdegreeZero : left.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hleft
    rw [Polynomial.eq_C_of_degree_eq_zero hdegreeZero] at hleftPrimitive ⊢
    exact Or.inl (Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hleftPrimitive
        (left.coeff 0) dvd_rfl))
  · have hrightPrimitive : right.IsPrimitive :=
      Polynomial.isPrimitive_of_dvd hprimitive
        (Dvd.intro_left left hfactor.symm)
    have hdegreeMapped :
        (right.map (Int.castRingHom (ZMod p))).degree = right.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hrightLeading
    have hdegreeZero : right.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hright
    rw [Polynomial.eq_C_of_degree_eq_zero hdegreeZero] at hrightPrimitive ⊢
    exact Or.inr (Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hrightPrimitive
        (right.coeff 0) dvd_rfl))

/-- A nontrivial factorization of a primitive integer polynomial whose leading
coefficient survives modulo `p` remains nontrivial after reduction.  This is
the converse direction needed to feed integer reducibility into the concrete
smaller-candidate construction. -/
theorem primitive_factorization_maps_nonunits
    (p : Nat) [Fact (Nat.Prime p)]
    (factor left right : Polynomial Int)
    (hprimitive : factor.IsPrimitive)
    (hleading : (factor.leadingCoeff : ZMod p) ≠ 0)
    (hfactor : factor = left * right)
    (hleftNonunit : ¬IsUnit left) (hrightNonunit : ¬IsUnit right) :
    left.IsPrimitive ∧ right.IsPrimitive ∧
      (left.leadingCoeff : ZMod p) ≠ 0 ∧
      (right.leadingCoeff : ZMod p) ≠ 0 ∧
      ¬IsUnit (left.map (Int.castRingHom (ZMod p))) ∧
      ¬IsUnit (right.map (Int.castRingHom (ZMod p))) := by
  have hleftPrimitive : left.IsPrimitive :=
    Polynomial.isPrimitive_of_dvd hprimitive
      (Dvd.intro right hfactor.symm)
  have hrightPrimitive : right.IsPrimitive :=
    Polynomial.isPrimitive_of_dvd hprimitive
      (Dvd.intro_left left hfactor.symm)
  have hleadingProduct :
      (left.leadingCoeff : ZMod p) * (right.leadingCoeff : ZMod p) =
        (factor.leadingCoeff : ZMod p) := by
    rw [← Int.cast_mul]
    congr 1
    rw [← Polynomial.leadingCoeff_mul, hfactor]
  have hleftLeading : (left.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, zero_mul] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hrightLeading : (right.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hleftMapNonunit :
      ¬IsUnit (left.map (Int.castRingHom (ZMod p))) := by
    intro hunit
    have hdegreeMapped :
        (left.map (Int.castRingHom (ZMod p))).degree = left.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hleftLeading
    have hdegreeZero : left.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hunit
    have hleftEq := Polynomial.eq_C_of_degree_eq_zero hdegreeZero
    have hconstantPrimitive := hleftPrimitive
    rw [hleftEq] at hconstantPrimitive
    apply hleftNonunit
    rw [hleftEq]
    exact Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hconstantPrimitive
        (left.coeff 0) dvd_rfl)
  have hrightMapNonunit :
      ¬IsUnit (right.map (Int.castRingHom (ZMod p))) := by
    intro hunit
    have hdegreeMapped :
        (right.map (Int.castRingHom (ZMod p))).degree = right.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hrightLeading
    have hdegreeZero : right.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hunit
    have hrightEq := Polynomial.eq_C_of_degree_eq_zero hdegreeZero
    have hconstantPrimitive := hrightPrimitive
    rw [hrightEq] at hconstantPrimitive
    apply hrightNonunit
    rw [hrightEq]
    exact Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hconstantPrimitive
        (right.coeff 0) dvd_rfl)
  exact ⟨hleftPrimitive, hrightPrimitive, hleftLeading, hrightLeading,
    hleftMapNonunit, hrightMapNonunit⟩

/-- Integer reducibility of an actually extracted primitive factor produces a
strictly smaller legal candidate over the same physical active array. -/
theorem smaller_active_candidate_of_reducible_primitive_factor
    (base : Nat) [Fact (Nat.Prime base)]
    (active : Array SparsePolyZZ) (outer : Array Nat)
    (houter : StrictRecombine.LegalCombination active.size outer.size outer)
    (hirreducible : ∀ index (hindex : index < active.size),
      Irreducible (StrictHensel.toPolyMod base active[index]))
    (factor left right : Polynomial Int)
    (hprimitive : factor.IsPrimitive)
    (hleading : (factor.leadingCoeff : ZMod base) ≠ 0)
    (hfactor : factor = left * right)
    (hleftNonunit : ¬IsUnit left) (hrightNonunit : ¬IsUnit right)
    (hassociated : Associated
      (factor.map (Int.castRingHom (ZMod base)))
      (((StrictRecombine.selectSourceIndices active.toList outer.toList).map
        (StrictHensel.toPolyMod base)).prod)) :
    ∃ inner : Array Nat,
      StrictRecombine.LegalCombination active.size inner.size inner ∧
      0 < inner.size ∧
      inner.size < outer.size ∧
      left.IsPrimitive ∧
      (left.leadingCoeff : ZMod base) ≠ 0 ∧
      Associated (left.map (Int.castRingHom (ZMod base)))
        (((StrictRecombine.selectSourceIndices active.toList inner.toList).map
          (StrictHensel.toPolyMod base)).prod) := by
  rcases primitive_factorization_maps_nonunits base factor left right
      hprimitive hleading hfactor hleftNonunit hrightNonunit with
    ⟨hleftPrimitive, _hrightPrimitive, hleftLeading, _hrightLeading,
      hleftMapNonunit, hrightMapNonunit⟩
  have hmapped : factor.map (Int.castRingHom (ZMod base)) =
      left.map (Int.castRingHom (ZMod base)) *
        right.map (Int.castRingHom (ZMod base)) := by
    rw [hfactor, Polynomial.map_mul]
  rcases smaller_active_candidate_of_reducible_selected_product base active
      outer houter hirreducible
      (factor.map (Int.castRingHom (ZMod base)))
      (left.map (Int.castRingHom (ZMod base)))
      (right.map (Int.castRingHom (ZMod base))) hmapped hleftMapNonunit
      hrightMapNonunit hassociated with
    ⟨inner, hinner, hsmall, hinnerAssociated⟩
  have hpositive : 0 < inner.size := by
    by_contra hnot
    have hempty : inner = #[] := Array.size_eq_zero_iff.mp
      (Nat.eq_zero_of_not_pos hnot)
    subst inner
    apply hleftMapNonunit
    simpa [StrictRecombine.selectSourceIndices] using
      hinnerAssociated.isUnit_iff.mpr isUnit_one
  exact ⟨inner, hinner, hpositive, hsmall, hleftPrimitive, hleftLeading,
    hinnerAssociated⟩

/-- The sparse factor array returned by a genuinely refined successful prime
selection consists pointwise of the irreducible factors recorded by that
selection's concrete candidate certificate. -/
theorem selectionFactors_irreducible
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection) :
    ∀ index (hindex : index < selection.factors.size),
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat
        selection.factors[index]) := by
  intro index hindex
  have hmember : SparsePolyZp.toPoly selection.prime.toNat
      selection.factors[index] ∈
        StrictSelectPrime.factorArrayToL2 selection.prime.toNat
          selection.factors := by
    simp only [StrictSelectPrime.factorArrayToL2, List.mem_map]
    exact ⟨selection.factors[index], Array.getElem_mem_toList hindex, rfl⟩
  exact (hselection.quality _ hmember).1

/-- The literal singleton/irreducible return branch of
`__factor_squarefree_primitive_ZZ` is sound.  Prime selection has executed the
real modular factorization; if its nonempty factor list has at most one
member, that member is associated to the complete reduction of the source
and is irreducible.  Gauss reduction then proves the primitive integer source
itself irreducible. -/
theorem selection_atMostOne_refines_singleton_FactorZZCorrect
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hcount : selection.factors.size ≤ 1) :
    FactorZZCorrect (SparsePolyZZ.toPoly f) [SparsePolyZZ.toPoly f] := by
  let decoded := StrictSelectPrime.factorArrayToL2 selection.prime.toNat
    selection.factors
  have hdecodedNonempty : decoded ≠ [] :=
    hselection.factors_nonempty hdegree
  have hdecodedLength : decoded.length = selection.factors.size := by
    simp [decoded, StrictSelectPrime.factorArrayToL2]
  have hlengthOne : decoded.length = 1 := by
    have hpositive : 0 < decoded.length := by
      apply Nat.pos_of_ne_zero
      intro hzero
      exact hdecodedNonempty (List.length_eq_zero_iff.mp hzero)
    omega
  rcases List.length_eq_one_iff.mp hlengthOne with ⟨factor, hdecoded⟩
  have hfactorIrreducible : Irreducible factor := by
    exact (hselection.quality factor (by simp [decoded, hdecoded])).1
  have hmappedAssociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f)) factor := by
    simpa [decoded, hdecoded] using hselection.productAssociated
  letI : Fact (Nat.Prime selection.prime.toNat) :=
    ⟨hselection.goodPrime.prime⟩
  have hmappedIrreducible : Irreducible
      (Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f)) :=
    hmappedAssociated.symm.irreducible hfactorIrreducible
  have hsourceIrreducible := primitive_irreducible_of_irreducible_mod
    selection.prime.toNat (SparsePolyZZ.toPoly f) hprimitive
      hselection.goodPrime.lc_nonzero hmappedIrreducible
  refine ⟨?_, ?_⟩
  · simpa using Associated.refl (SparsePolyZZ.toPoly f)
  · intro candidate hcandidate
    simp only [List.mem_singleton] at hcandidate
    subst candidate
    exact hsourceIrreducible

/-- The exact first-factor adjustment performed by Hensel preserves the
pointwise irreducibility supplied by the actual SelectPrime result. -/
theorem selectionFactors_adjusted_irreducible
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {adjusted : Array SparsePolyZp}
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hadjust : StrictHensel.HenselAdjustFirstFactorCorrect f
      selection.factors selection.prime adjusted) :
    ∀ index (hindex : index < adjusted.size),
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat
        adjusted[index]) := by
  have hcandidate : StrictSelectPrime.CandidateCorrect
      (SparsePolyZZ.toPoly f) selection.prime.toNat
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors) := hselection
  letI : Fact (Nat.Prime selection.prime.toNat) :=
    ⟨hcandidate.goodPrime.prime⟩
  have hleadingNonzero : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) ≠ 0 := by
    intro leading hleading
    rw [hleadingSemantic leading hleading]
    exact hcandidate.goodPrime.lc_nonzero
  have hrel := hadjust.unitRel hfactors hleadingNonzero
  intro index hindex
  exact hrel.irreducible (selectionFactors_irreducible hselection) index hindex

/-- The actual first-factor adjustment turns the monic finite-field factors
returned by the selected `__factor_Zp` execution into an *exact* factorization
of the integer source reduced modulo the selected prime.  This recovers the
leading unit from the concrete C++ write instead of weakening the Hensel input
to associatedness. -/
theorem selectionAdjustedFactors_product_eq_source
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {adjusted : Array SparsePolyZp}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hadjust : StrictHensel.HenselAdjustFirstFactorCorrect f
      selection.factors selection.prime adjusted) :
    (adjusted.toList.map
        (SparsePolyZp.toPoly selection.prime.toNat)).prod =
      Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f) := by
  let sourceMod := Polynomial.map
    (Int.castRingHom (ZMod selection.prime.toNat)) (SparsePolyZZ.toPoly f)
  let selected := StrictSelectPrime.factorArrayToL2
    selection.prime.toNat selection.factors
  have hselectedMonic : selected.prod.Monic := by
    apply monic_list_product
    intro factor hfactor
    exact (hselection.quality factor hfactor).2
  have hsourceLeading : sourceMod.leadingCoeff =
      ((SparsePolyZZ.toPoly f).leadingCoeff :
        ZMod selection.prime.toNat) := by
    exact Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _
      hselection.goodPrime.lc_nonzero
  have hsourceNonzero : sourceMod ≠ 0 := by
    intro hzero
    have := congrArg Polynomial.leadingCoeff hzero
    rw [hsourceLeading] at this
    exact hselection.goodPrime.lc_nonzero (by simpa using this)
  have hselectedExact :
      sourceMod = Polynomial.C sourceMod.leadingCoeff * selected.prod :=
    polynomial_eq_C_leadingCoeff_mul_of_associated_monic hsourceNonzero
      hselectedMonic hselection.productAssociated
  rcases hadjust.product_eq hfactors with
    ⟨leading, hleading, hadjustedProduct⟩
  rw [hadjustedProduct]
  change Polynomial.C (leading.2 : ZMod selection.prime.toNat) *
      selected.prod = sourceMod
  rw [hleadingSemantic leading hleading, ← hsourceLeading]
  exact hselectedExact.symm

/-- The exact selected-prime factorization now discharges the formerly
explicit source-product premise of the real Hensel loop.  Consequently the
pre-normalized factors extracted from the generated C++ tree multiply to the
integer source at the exact modulus returned by that loop. -/
theorem selectionHenselFactors_preNormalization_product
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ adjusted, ∃ extracted, ∃ outputM : Nat,
      StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors
        selection.prime adjusted ∧
      StrictHensel.HenselNormalizeCorrect extracted (outputM : Int) output.1 ∧
      output.2 = (outputM : Int) ∧
      (extracted.toList.map (StrictHensel.toPolyMod outputM)).prod =
        StrictHensel.toPolyMod outputM f := by
  apply hentry.preNormalizationProduct hcount
  intro adjusted hadjust
  have hadjustSize : adjusted.size = selection.factors.size := by
    cases hadjust
    simp [Array.set!]
  rw [← hadjustSize]
  rw [StrictHensel.henselFactorRangeList_full]
  have hexact := selectionAdjustedFactors_product_eq_source hfactors
    hleadingSemantic hselection hadjust
  simpa [StrictHensel.toPolyMod] using hexact

/-- Product invariant of the public normalized Hensel output.  The only
remaining discrepancy from the source is the concrete unit computed by the
generated normalization branch; both the modulus and the output array are the
ones returned by the actual C++-shaped entry execution. -/
theorem selectionHenselFactors_normalized_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ outputM : Nat, ∃ scale : ZMod outputM,
      output.2 = (outputM : Int) ∧ IsUnit scale ∧
      (output.1.toList.map (StrictHensel.toPolyMod outputM)).prod =
        Polynomial.C scale * StrictHensel.toPolyMod outputM f := by
  rcases selectionHenselFactors_preNormalization_product hcount hfactors
      hleadingSemantic hselection hentry with
    ⟨_adjusted, extracted, outputM, _hadjust, hnormalize, houtputM,
      hpreProduct⟩
  rcases hnormalize.product_eq_unit_mul with
    ⟨scale, hscaleUnit, hnormalizedProduct⟩
  refine ⟨outputM, scale, houtputM, hscaleUnit, ?_⟩
  rw [hnormalizedProduct, hpreProduct]

/-- Prime-power specialization of the normalized full-product invariant.  The
exponent comes from the actual well-founded Hensel execution and the scale is
the concrete normalization unit at that same returned modulus. -/
theorem selectionHenselFactors_primePower_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ exponent : Nat, 0 < exponent ∧
      ∃ scale : ZMod (selection.prime.toNat ^ exponent),
        output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
        IsUnit scale ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent))).prod =
          Polynomial.C scale *
            StrictHensel.toPolyMod
              (selection.prime.toNat ^ exponent) f := by
  rcases selectionHenselFactors_normalized_product_eq_unit_mul_source
      hcount hfactors hleadingSemantic hselection hentry with
    ⟨outputM, scale, houtputM, hscaleUnit, hproduct⟩
  rcases hentry.outputModulus_eq_prime_pow with
    ⟨exponent, hexponent, hprimePower⟩
  have hmodulus : outputM = selection.prime.toNat ^ exponent := by
    exact_mod_cast houtputM.symm.trans hprimePower
  subst outputM
  exact ⟨exponent, hexponent, scale, hprimePower, hscaleUnit, hproduct⟩

/-- Reduce the concrete prime-power normalization unit and full-product
identity back to the selected prime.  Both unit witnesses are images of the
same physical normalization scalar. -/
theorem selectionHenselFactors_prime_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ exponent : Nat, 0 < exponent ∧
      ∃ scale : ZMod (selection.prime.toNat ^ exponent),
      ∃ scaleAtPrime : ZMod selection.prime.toNat,
        output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
        IsUnit scale ∧ IsUnit scaleAtPrime ∧
        scaleAtPrime = (scale.val : ZMod selection.prime.toNat) ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent))).prod =
          Polynomial.C scale * StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent) f ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod =
          Polynomial.C scaleAtPrime *
            StrictHensel.toPolyMod selection.prime.toNat f := by
  rcases selectionHenselFactors_primePower_product_eq_unit_mul_source
      hcount hfactors hleadingSemantic hselection hentry with
    ⟨exponent, hexponent, scale, houtput, hscaleUnit, hlargeProduct⟩
  let prime := selection.prime.toNat
  have hprime : Nat.Prime prime := Fact.out
  have hpowerNe : prime ^ exponent ≠ 0 :=
    pow_ne_zero exponent hprime.ne_zero
  letI : NeZero (prime ^ exponent) := ⟨hpowerNe⟩
  let scaleAtPrime : ZMod prime := scale.val
  have hdivides : prime ∣ prime ^ exponent :=
    dvd_pow_self prime (Nat.ne_of_gt hexponent)
  have hscaleAtPrimeUnit : IsUnit scaleAtPrime := by
    have hmapped := hscaleUnit.map
      (ZMod.castHom hdivides (ZMod prime))
    rw [ZMod.castHom_apply, ZMod.cast_eq_val] at hmapped
    exact hmapped
  let scaleInt : Int := Int.ofNat scale.val
  let fullInteger : Polynomial Int :=
    (output.1.toList.map SparsePolyZZ.toPoly).prod
  let scaledSource : Polynomial Int :=
    Polynomial.C scaleInt * SparsePolyZZ.toPoly f
  have hscaleLarge : (scaleInt : ZMod (prime ^ exponent)) = scale := by
    simp [scaleInt, ZMod.natCast_zmod_val]
  have hmapFullLarge : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent))) fullInteger =
      (output.1.toList.map
        (StrictHensel.toPolyMod (prime ^ exponent))).prod := by
    simp only [fullInteger, Polynomial.map_list_prod]
    apply congrArg List.prod
    rw [List.map_map]
    apply List.map_congr_left
    intro factor hfactor
    rfl
  have hmapScaledLarge : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent))) scaledSource =
      Polynomial.C scale * StrictHensel.toPolyMod
        (prime ^ exponent) f := by
    simp only [scaledSource, Polynomial.map_mul, Polynomial.map_C]
    change Polynomial.C (scaleInt : ZMod (prime ^ exponent)) *
        StrictHensel.toPolyMod (prime ^ exponent) f = _
    rw [hscaleLarge]
  have hlargeInteger : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent)))
        fullInteger =
      Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
        scaledSource :=
    hmapFullLarge.trans (hlargeProduct.trans hmapScaledLarge.symm)
  have hprimeInteger := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime (prime ^ exponent) hdivides
    fullInteger scaledSource hlargeInteger
  have hmapFullPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) fullInteger =
      (output.1.toList.map (StrictHensel.toPolyMod prime)).prod := by
    simp only [fullInteger, Polynomial.map_list_prod]
    apply congrArg List.prod
    rw [List.map_map]
    apply List.map_congr_left
    intro factor hfactor
    rfl
  have hscalePrime : (scaleInt : ZMod prime) = scaleAtPrime := by
    simp [scaleInt, scaleAtPrime]
  have hmapScaledPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) scaledSource =
      Polynomial.C scaleAtPrime * StrictHensel.toPolyMod prime f := by
    simp only [scaledSource, Polynomial.map_mul, Polynomial.map_C]
    change Polynomial.C (scaleInt : ZMod prime) *
        StrictHensel.toPolyMod prime f = _
    rw [hscalePrime]
  have hprimeProduct :
      (output.1.toList.map (StrictHensel.toPolyMod prime)).prod =
        Polynomial.C scaleAtPrime * StrictHensel.toPolyMod prime f := by
    rw [← hmapFullPrime, hprimeInteger, hmapScaledPrime]
  exact ⟨exponent, hexponent, scale, scaleAtPrime, houtput, hscaleUnit,
    hscaleAtPrimeUnit, rfl, hlargeProduct, hprimeProduct⟩

/-- The full-precision algebraic state carried by the real Zassenhaus loop.
Both product equations refer to the current physical active array.  The two
scalars are the concrete Hensel normalization unit and its reduction, so this
certificate does not assert that recombination succeeds. -/
structure LiveHenselProduct (prime exponent : Nat) (source : SparsePolyZZ)
    (active : Array SparsePolyZZ) : Prop where
  exponentPositive : 0 < exponent
  certificate :
    ∃ scale : ZMod (prime ^ exponent), ∃ scaleAtPrime : ZMod prime,
      IsUnit scale ∧ IsUnit scaleAtPrime ∧
      (active.toList.map
        (StrictHensel.toPolyMod (prime ^ exponent))).prod =
          Polynomial.C scale *
            StrictHensel.toPolyMod (prime ^ exponent) source ∧
      (active.toList.map (StrictHensel.toPolyMod prime)).prod =
        Polynomial.C scaleAtPrime * StrictHensel.toPolyMod prime source

/-- The concrete coefficient-recovery margin carried by the real Zassenhaus
loop.  It states exactly the bounds consumed by the generated leading and
symmetric-recovery branches for every genuine factorization of the current
source. -/
structure LiveRecoveryPrecision (modulus : Nat)
    (source : SparsePolyZZ) : Prop where
  leadingBound : (SparsePolyZZ.toPoly source).leadingCoeff.natAbs * 2 < modulus
  scaledFactorBound : ∀ divisor quotient,
    SparsePolyZZ.toPoly source = divisor * quotient → divisor ≠ 0 →
    ∀ degree,
      ((Polynomial.C quotient.leadingCoeff * divisor).coeff degree).natAbs *
        2 < modulus

/-- Member-level correctness of the physical result array accumulated by the
generated Zassenhaus loop. -/
def FactorArrayIrreducible (result : Array SparsePolyZZ) : Prop :=
  ∀ factor ∈ result.toList, Irreducible (SparsePolyZZ.toPoly factor)

theorem FactorArrayIrreducible.empty : FactorArrayIrreducible #[] := by
  simp [FactorArrayIrreducible]

theorem FactorArrayIrreducible.push
    {result : Array SparsePolyZZ} {factor : SparsePolyZZ}
    (hresult : FactorArrayIrreducible result)
    (hfactor : Irreducible (SparsePolyZZ.toPoly factor)) :
    FactorArrayIrreducible (result.push factor) := by
  intro candidate hcandidate
  rw [Array.toList_push] at hcandidate
  rcases List.mem_append.mp hcandidate with hmember | hlast
  · exact hresult candidate hmember
  · have heq : candidate = factor := by simpa using hlast
    subst candidate
    exact hfactor

theorem FactorArrayIrreducible.sortFactorsByDegree
    {result : Array SparsePolyZZ}
    (hresult : FactorArrayIrreducible result) :
    FactorArrayIrreducible
      (Generated.StrictRecombine.sortFactorsByDegree result) := by
  intro factor hfactor
  unfold Generated.StrictRecombine.sortFactorsByDegree at hfactor
  have hsorted : factor ∈ result.toList.mergeSort (fun left right =>
      left[0]!.1.deg < right[0]!.1.deg) := by simpa using hfactor
  have horiginal : factor ∈ result.toList :=
    (List.Perm.mem_iff (List.mergeSort_perm _ _)).mp hsorted
  exact hresult factor horiginal

theorem FactorArrayIrreducible.finishZassenhaus
    {fStar : SparsePolyZZ} {result : Array SparsePolyZZ}
    (hresult : FactorArrayIrreducible result)
    (hfStar : 0 < fStar.size → 0 < fStar[0]!.1.deg →
      Irreducible (SparsePolyZZ.toPoly fStar)) :
    FactorArrayIrreducible
      (Generated.StrictRecombine.finishZassenhaus fStar result) := by
  unfold Generated.StrictRecombine.finishZassenhaus
  apply FactorArrayIrreducible.sortFactorsByDegree
  split
  next hnonempty =>
    split
    next hdegree =>
      exact hresult.push (hfStar hnonempty (by
        simpa [getElem!_pos fStar 0 hnonempty] using hdegree))
    next hdegree => exact hresult
  next hnonempty => exact hresult

/-- Modular irreducibility of every physical output lifts back to integer
irreducibility when the actual accumulated product is primitive and every
leading coefficient survives reduction. -/
theorem factorArrayIrreducible_of_modular
    (prime : Nat) [Fact (Nat.Prime prime)]
    (result : Array SparsePolyZZ)
    (hprimitive :
      (result.toList.map SparsePolyZZ.toPoly).prod.IsPrimitive)
    (hleading : ∀ factor ∈ result.toList,
      ((SparsePolyZZ.toPoly factor).leadingCoeff : ZMod prime) ≠ 0)
    (hmodIrreducible : ∀ factor ∈ result.toList,
      Irreducible (StrictHensel.toPolyMod prime factor)) :
    FactorArrayIrreducible result := by
  intro factor hfactor
  have hfactorMem : SparsePolyZZ.toPoly factor ∈
      result.toList.map SparsePolyZZ.toPoly :=
    List.mem_map.mpr ⟨factor, hfactor, rfl⟩
  have hfactorPrimitive : (SparsePolyZZ.toPoly factor).IsPrimitive :=
    Polynomial.isPrimitive_of_dvd hprimitive
      (List.dvd_prod hfactorMem)
  exact primitive_irreducible_of_irreducible_mod prime
    (SparsePolyZZ.toPoly factor) hfactorPrimitive (hleading factor hfactor)
    (by simpa [StrictHensel.toPolyMod] using hmodIrreducible factor hfactor)

/-- The selected-prime product relation exposed in the orientation used by
the live divisor-to-subset theorem. -/
theorem LiveHenselProduct.primeProductAssociated
    {prime exponent : Nat} {source : SparsePolyZZ}
    {active : Array SparsePolyZZ}
    (state : LiveHenselProduct prime exponent source active) :
    Associated
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source))
      (active.toList.map (StrictHensel.toPolyMod prime)).prod := by
  rcases state.certificate with
    ⟨_scale, scaleAtPrime, _hscaleUnit, hscaleAtPrimeUnit,
      _hlargeProduct, hprimeProduct⟩
  rw [hprimeProduct]
  exact ((associated_isUnit_mul_left_iff
    (Polynomial.isUnit_C.mpr hscaleAtPrimeUnit)).mpr
      (Associated.refl _)).symm

/-- Equal-cardinality output of a concrete recombination run is irreducible
when it preserves the primitive integer product and the validation facts for
every physical member.  The modular product relation is derived from the
actual Hensel certificate and the actual integer recombination product. -/
theorem factorArrayIrreducible_of_hensel_cardinality
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {source : SparsePolyZZ} {lifted output : Array SparsePolyZZ}
    (productState : LiveHenselProduct prime exponent source lifted)
    (hsourceProduct : Associated (SparsePolyZZ.toPoly source)
      (output.toList.map SparsePolyZZ.toPoly).prod)
    (hprimitive :
      (output.toList.map SparsePolyZZ.toPoly).prod.IsPrimitive)
    (hleading : ∀ factor ∈ output.toList,
      ((SparsePolyZZ.toPoly factor).leadingCoeff : ZMod prime) ≠ 0)
    (hnonunit : ∀ factor ∈ output.toList,
      ¬ IsUnit (StrictHensel.toPolyMod prime factor))
    (hirreducible : ∀ index (hindex : index < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[index]))
    (hlength : output.size = lifted.size) :
    FactorArrayIrreducible output := by
  let atoms := lifted.toList.map (StrictHensel.toPolyMod prime)
  have hatoms : ∀ atom ∈ atoms, Irreducible atom := by
    intro atom hatom
    rcases List.mem_map.mp hatom with ⟨factor, hfactor, rfl⟩
    rcases List.mem_iff_getElem.mp hfactor with ⟨index, hindex, rfl⟩
    exact hirreducible index (by simpa using hindex)
  have hsourceMapped := hsourceProduct.map
    (Polynomial.mapRingHom (Int.castRingHom (ZMod prime)))
  have houtputProduct : Associated
      ((output.toList.map (StrictHensel.toPolyMod prime)).prod)
      atoms.prod := by
    have hmapped : Associated
        (Polynomial.map (Int.castRingHom (ZMod prime))
          (SparsePolyZZ.toPoly source))
        ((output.toList.map (StrictHensel.toPolyMod prime)).prod) := by
      simpa [StrictHensel.toPolyMod, Polynomial.map_list_prod] using
        hsourceMapped
    exact hmapped.symm.trans (by
      simpa [atoms] using productState.primeProductAssociated)
  have hmodIrreducible :=
    StrictRecombine.modular_irreducible_members_of_equal_length_associated_product
      prime output atoms hatoms hleading hnonunit houtputProduct (by
        simpa [atoms] using hlength)
  exact factorArrayIrreducible_of_modular prime output hprimitive hleading
    hmodIrreducible

/-- A full-cardinality result from the literal concrete van-Hoeij entry is an
array of integer irreducibles.  Product, member quality, and cardinality all
come from the same generated execution; the Hensel certificate supplies the
modular atom product. -/
theorem concreteVanHoeij_equal_cardinality_factorArrayIrreducible
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {source : SparsePolyZZ} {lifted output : Array SparsePolyZZ}
    (productState : LiveHenselProduct prime exponent source lifted)
    (hdegree : 2 ≤ (get_deg source).toNatClampNeg)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hliftedFits : lifted.size ≤ 2 ^ 31)
    (hirreducible : ∀ index (hindex : index < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[index]))
    (hlength : output.size = lifted.size)
    (hrun : Generated.StrictRecombine.__vanhoeij_recombine_raw_ir
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination source lifted
      ((prime ^ exponent : Nat) : ZZ) hdegree = .ok output) :
    FactorArrayIrreducible output := by
  have hmodulus : 0 < prime ^ exponent :=
    pow_pos (Fact.out : Nat.Prime prime).pos exponent
  have hdivides : prime ∣ prime ^ exponent :=
    dvd_pow_self prime (Nat.ne_of_gt productState.exponentPositive)
  have hquality :=
    StrictRecombine.__vanhoeij_recombine_raw_ir_physicalFactorQuality
      source lifted output (prime ^ exponent) prime hdegree hmodulus
      (Fact.out : Nat.Prime prime).pos hdivides hcanonical hnonempty
      hprimitive hleading hliftedFits hirreducible hrun
  have hsourceProduct :=
    StrictRecombine.__vanhoeij_recombine_raw_ir_product_associated
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination source lifted
      (((prime ^ exponent : Nat) : ZZ)) hdegree hcanonical hprimitive output
      hrun
  have hsourceProduct' : Associated (SparsePolyZZ.toPoly source)
      (output.toList.map SparsePolyZZ.toPoly).prod := by
    simpa [StrictRecombine.factorArrayProduct] using hsourceProduct
  have houtputPrimitive :
      (output.toList.map SparsePolyZZ.toPoly).prod.IsPrimitive :=
    Polynomial.isPrimitive_of_dvd hprimitive hsourceProduct'.symm.dvd
  exact factorArrayIrreducible_of_hensel_cardinality productState
    hsourceProduct' houtputPrimitive
    (fun factor hfactor => (hquality factor hfactor).2.2)
    (fun factor hfactor => (hquality factor hfactor).2.1)
    hirreducible hlength

/-- A concrete integer coefficient surviving reduction modulo the selected
prime is a unit at every positive prime-power precision. -/
theorem intCast_isUnit_primePower_of_ne_zero_prime
    (prime exponent : Nat) [Fact (Nat.Prime prime)] (coefficient : Int)
    (hnonzero : (coefficient : ZMod prime) ≠ 0) :
    IsUnit (coefficient : ZMod (prime ^ exponent)) := by
  have hprimeUnit : IsUnit (coefficient : ZMod prime) :=
    isUnit_iff_ne_zero.mpr hnonzero
  have hcoprime : IsCoprime (prime : Int) coefficient :=
    (ZMod.coe_int_isUnit_iff_isCoprime coefficient prime).mp hprimeUnit
  have hpowerCoprime : IsCoprime ((prime : Int) ^ exponent) coefficient :=
    hcoprime.pow_left
  apply (ZMod.coe_int_isUnit_iff_isCoprime coefficient
    (prime ^ exponent)).mpr
  simpa using hpowerCoprime

/-- Cancellation in the prime-power polynomial ring is justified by a unit
leading coefficient, not by a false field/domain instance for `ZMod (p^k)`. -/
theorem polynomial_mul_right_cancel_of_isUnit_leadingCoeff
    {R : Type*} [CommRing R] (factor left right : Polynomial R)
    (hleading : IsUnit factor.leadingCoeff)
    (heq : factor * left = factor * right) : left = right := by
  have hzero : factor * (left - right) = 0 := by
    rw [mul_sub, heq, sub_self]
  have hdiff : left - right = 0 :=
    (Polynomial.isUnit_leadingCoeff_mul_right_eq_zero_iff hleading).mp hzero
  exact sub_eq_zero.mp hdiff

/-- Algebraic update used after a physical candidate removal.  It constructs
the next normalization unit from the four concrete units in the execution
trace and cancels the extracted factor only through its unit leading
coefficient. -/
theorem remaining_product_eq_unit_mul_quotient
    {R : Type*} [CommRing R]
    (selected remaining source factor quotient : Polynomial R)
    (scale leading content scalar : R)
    (hscale : IsUnit scale) (hleading : IsUnit leading)
    (hcontent : IsUnit content) (hscalar : IsUnit scalar)
    (hfactorLeading : IsUnit factor.leadingCoeff)
    (hactive : selected * remaining = Polynomial.C scale * source)
    (hextraction : source =
      Polynomial.C scalar * (factor * quotient))
    (hfactor : Polynomial.C content * factor =
      Polynomial.C leading * selected) :
    ∃ nextScale : R, IsUnit nextScale ∧
      remaining = Polynomial.C nextScale * quotient := by
  have hcancelEquation :
      factor * (Polynomial.C content * remaining) =
        factor * (Polynomial.C (leading * scale * scalar) * quotient) := by
    calc
      factor * (Polynomial.C content * remaining) =
          (Polynomial.C content * factor) * remaining := by ring
      _ = (Polynomial.C leading * selected) * remaining := by rw [hfactor]
      _ = Polynomial.C leading * (selected * remaining) := by ring
      _ = Polynomial.C leading * (Polynomial.C scale * source) := by
        rw [hactive]
      _ = factor *
          (Polynomial.C (leading * scale * scalar) * quotient) := by
        rw [hextraction]
        simp only [Polynomial.C_mul]
        ring
  have hcancelled : Polynomial.C content * remaining =
      Polynomial.C (leading * scale * scalar) * quotient :=
    polynomial_mul_right_cancel_of_isUnit_leadingCoeff factor _ _
      hfactorLeading hcancelEquation
  let contentInv : R := ↑(hcontent.unit⁻¹)
  have hcontentInv : contentInv * content = 1 := by
    have hspec : (↑hcontent.unit : R) = content := hcontent.unit_spec
    calc
      contentInv * content = contentInv * (↑hcontent.unit : R) := by
        exact congrArg (contentInv * ·) hspec.symm
      _ = 1 := by simp [contentInv]
  let nextUnit : Rˣ := hcontent.unit⁻¹ * hleading.unit * hscale.unit *
    hscalar.unit
  let nextScale : R := ↑nextUnit
  have hnextUnit : IsUnit nextScale := nextUnit.isUnit
  have hnextScale : nextScale = contentInv * leading * scale * scalar := by
    simp [nextScale, nextUnit, contentInv, hleading.unit_spec,
      hscale.unit_spec, hscalar.unit_spec]
  refine ⟨nextScale, hnextUnit, ?_⟩
  calc
    remaining = Polynomial.C (contentInv * content) * remaining := by
      rw [hcontentInv]
      simp
    _ = Polynomial.C contentInv * (Polynomial.C content * remaining) := by
      rw [Polynomial.C_mul]
      ring
    _ = Polynomial.C contentInv *
        (Polynomial.C (leading * scale * scalar) * quotient) := by
      rw [hcancelled]
    _ = Polynomial.C nextScale * quotient := by
      rw [hnextScale]
      simp only [Polynomial.C_mul]
      ring

/-- Variant used by the physical Zassenhaus state: cancel the selected lifted
subproduct, whose leading coefficient is a unit because every active Hensel
cell is monic.  No leading-coefficient premise is needed for the recovered
integer factor. -/
theorem remaining_product_eq_unit_mul_quotient_cancel_selected
    {R : Type*} [CommRing R]
    (selected remaining source factor quotient : Polynomial R)
    (scale leading content scalar : R)
    (hscale : IsUnit scale) (hleading : IsUnit leading)
    (hcontent : IsUnit content) (hscalar : IsUnit scalar)
    (hselectedLeading : IsUnit selected.leadingCoeff)
    (hactive : selected * remaining = Polynomial.C scale * source)
    (hextraction : source =
      Polynomial.C scalar * (factor * quotient))
    (hfactor : Polynomial.C content * factor =
      Polynomial.C leading * selected) :
    ∃ nextScale : R, IsUnit nextScale ∧
      remaining = Polynomial.C nextScale * quotient := by
  have hcancelEquation :
      selected * (Polynomial.C content * remaining) =
        selected * (Polynomial.C (scale * scalar * leading) * quotient) := by
    calc
      selected * (Polynomial.C content * remaining) =
          Polynomial.C content * (selected * remaining) := by ring
      _ = Polynomial.C content * (Polynomial.C scale * source) := by
        rw [hactive]
      _ = Polynomial.C content *
          (Polynomial.C scale *
            (Polynomial.C scalar * (factor * quotient))) := by
        rw [hextraction]
      _ = (Polynomial.C content * factor) *
          (Polynomial.C (scale * scalar) * quotient) := by
        simp only [Polynomial.C_mul]
        ring
      _ = (Polynomial.C leading * selected) *
          (Polynomial.C (scale * scalar) * quotient) := by rw [hfactor]
      _ = selected *
          (Polynomial.C (scale * scalar * leading) * quotient) := by
        simp only [Polynomial.C_mul]
        ring
  have hcancelled : Polynomial.C content * remaining =
      Polynomial.C (scale * scalar * leading) * quotient :=
    polynomial_mul_right_cancel_of_isUnit_leadingCoeff selected _ _
      hselectedLeading hcancelEquation
  let contentInv : R := ↑(hcontent.unit⁻¹)
  have hcontentInv : contentInv * content = 1 := by
    have hspec : (↑hcontent.unit : R) = content := hcontent.unit_spec
    calc
      contentInv * content = contentInv * (↑hcontent.unit : R) := by
        exact congrArg (contentInv * ·) hspec.symm
      _ = 1 := by simp [contentInv]
  let nextUnit : Rˣ := hcontent.unit⁻¹ * hscale.unit * hscalar.unit *
    hleading.unit
  let nextScale : R := ↑nextUnit
  have hnextUnit : IsUnit nextScale := nextUnit.isUnit
  have hnextScale : nextScale = contentInv * scale * scalar * leading := by
    simp [nextScale, nextUnit, contentInv, hscale.unit_spec,
      hscalar.unit_spec, hleading.unit_spec]
  refine ⟨nextScale, hnextUnit, ?_⟩
  calc
    remaining = Polynomial.C (contentInv * content) * remaining := by
      rw [hcontentInv]
      simp
    _ = Polynomial.C contentInv * (Polynomial.C content * remaining) := by
      rw [Polynomial.C_mul]
      ring
    _ = Polynomial.C contentInv *
        (Polynomial.C (scale * scalar * leading) * quotient) := by
      rw [hcancelled]
    _ = Polynomial.C nextScale * quotient := by
      rw [hnextScale]
      simp only [Polynomial.C_mul]
      ring

/-- The mapped integer polynomial itself has a unit leading coefficient at
full Hensel precision whenever its leading coefficient survives mod `p`. -/
theorem mapped_leadingCoeff_isUnit_primePower
    (prime exponent : Nat) [Fact (Nat.Prime prime)]
    (factor : Polynomial Int) (hexponent : 0 < exponent)
    (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0) :
    IsUnit ((factor.map
      (Int.castRingHom (ZMod (prime ^ exponent)))).leadingCoeff) := by
  have hcoefficientUnit := intCast_isUnit_primePower_of_ne_zero_prime prime
    exponent factor.leadingCoeff hleading
  letI : Fact (1 < prime ^ exponent) :=
    ⟨Nat.one_lt_pow hexponent.ne' (Fact.out : Nat.Prime prime).one_lt⟩
  have hcoefficientNe :
      (factor.leadingCoeff : ZMod (prime ^ exponent)) ≠ 0 :=
    hcoefficientUnit.ne_zero
  rw [Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ hcoefficientNe]
  exact hcoefficientUnit

/-- The primitive content computed by one successful physical attempt is a
unit both at the selected prime and at full Hensel precision.  The same
integer content and both exact candidate equations come from one execution
trace. -/
theorem zassenhausAttempt_extracted_content_units
    (prime exponent : Nat) [Fact (Nat.Prime prime)]
    (source factor quotient : SparsePolyZZ)
    (active : Array SparsePolyZZ) (candidate : Array Nat)
    (hactive : StrictRecombine.LiveActiveFactors prime active)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hexponent : 0 < exponent)
    (hrun : Generated.StrictRecombine.zassenhausAttempt source active
      ((prime ^ exponent : Nat) : ZZ) candidate =
        .ok (.extracted factor quotient)) :
    ∃ content : Int,
      (Polynomial.C (content : ZMod (prime ^ exponent)) *
          StrictHensel.toPolyMod (prime ^ exponent) factor =
        Polynomial.C (source[0]!.2 : ZMod (prime ^ exponent)) *
          ((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map
              (StrictHensel.toPolyMod (prime ^ exponent))).prod) ∧
      (Polynomial.C (content : ZMod prime) *
          StrictHensel.toPolyMod prime factor =
        Polynomial.C (source[0]!.2 : ZMod prime) *
          ((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map
              (StrictHensel.toPolyMod prime)).prod) ∧
      IsUnit (content : ZMod (prime ^ exponent)) ∧
      IsUnit (content : ZMod prime) := by
  let modulus := prime ^ exponent
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime prime).pos exponent
  have hprime : 0 < prime := (Fact.out : Nat.Prime prime).pos
  have hdivides : prime ∣ modulus :=
    dvd_pow_self prime (Nat.ne_of_gt hexponent)
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  rcases StrictRecombine.zassenhausAttempt_extracted_factor_mod_eq_selected_pair
      source factor quotient active modulus modulus prime candidate hmodulus
      hmodulus hprime (dvd_refl modulus) hdivides hbound hactive.fitsInt32
      (by simpa [modulus] using hrun) with
    ⟨content, hlargeEquation, hprimeEquation⟩
  have hfront : (source[0]!.2 : ZMod prime) ≠ 0 := by
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source
      hcanonical hnonempty
    simpa [getElem!_pos source 0 hnonempty] using (hhead ▸ hleading)
  have hselectedMonic := hactive.selectedToPolyModMonic candidate hbound
    (modulus := prime)
  have hselectedNe := hselectedMonic.ne_zero
  have hrhsNe : Polynomial.C (source[0]!.2 : ZMod prime) *
      ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
        (StrictHensel.toPolyMod prime)).prod ≠ 0 :=
    mul_ne_zero (Polynomial.C_ne_zero.mpr hfront) hselectedNe
  have hcontentPrimeNe : (content : ZMod prime) ≠ 0 := by
    intro hzero
    apply hrhsNe
    rw [← hprimeEquation]
    simp [hzero]
  have hcontentPrimeUnit : IsUnit (content : ZMod prime) :=
    isUnit_iff_ne_zero.mpr hcontentPrimeNe
  have hcontentLargeUnit : IsUnit (content : ZMod modulus) :=
    intCast_isUnit_primePower_of_ne_zero_prime prime exponent content
      hcontentPrimeNe
  exact ⟨content, by simpa [modulus] using hlargeEquation, hprimeEquation,
    by simpa [modulus] using hcontentLargeUnit, hcontentPrimeUnit⟩

/-- A successful literal attempt followed by the literal physical removal
preserves the full live Hensel product state for the primitive quotient that
the generated outer loop installs. -/
theorem LiveHenselProduct.extract
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {source factor quotient : SparsePolyZZ}
    {active remaining : Array SparsePolyZZ} {candidate : Array Nat}
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hattempt : Generated.StrictRecombine.zassenhausAttempt source active
      ((prime ^ exponent : Nat) : ZZ) candidate =
        .ok (.extracted factor quotient))
    (hremove : Generated.StrictRecombine.removeCombination candidate active =
      .ok remaining) :
    LiveHenselProduct prime exponent quotient remaining := by
  let modulus := prime ^ exponent
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime prime).pos exponent
  have hprime : 0 < prime := (Fact.out : Nat.Prime prime).pos
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hlegalSize : StrictRecombine.LegalCombination active.size
      candidate.size candidate := hlegal
  rcases state.certificate with
    ⟨scaleLarge, scalePrime, hscaleLarge, hscalePrime,
      hproductLarge, hproductPrime⟩
  rcases zassenhausAttempt_extracted_content_units prime exponent source factor
      quotient active candidate activeState hcanonical hnonempty hleading
      hlegal state.exponentPositive hattempt with
    ⟨content, hfactorLarge, hfactorPrime, hcontentLarge, hcontentPrime⟩
  rcases StrictRecombine.zassenhausAttempt_extracted_unit_scalar source factor
      quotient active ((modulus : Nat) : ZZ) candidate hprimitive
      (by simpa [modulus] using hattempt) with
    ⟨scalar, hscalar, hextraction⟩
  have hscalarLarge : IsUnit (scalar : ZMod modulus) :=
    hscalar.map (Int.castRingHom (ZMod modulus))
  have hscalarPrime : IsUnit (scalar : ZMod prime) :=
    hscalar.map (Int.castRingHom (ZMod prime))
  have hextractionLarge : StrictHensel.toPolyMod modulus source =
      Polynomial.C (scalar : ZMod modulus) *
        (StrictHensel.toPolyMod modulus factor *
          StrictHensel.toPolyMod modulus quotient) := by
    simpa [StrictHensel.toPolyMod, Polynomial.map_mul, Polynomial.map_C] using
      congrArg (Polynomial.map (Int.castRingHom (ZMod modulus))) hextraction
  have hextractionPrime : StrictHensel.toPolyMod prime source =
      Polynomial.C (scalar : ZMod prime) *
        (StrictHensel.toPolyMod prime factor *
          StrictHensel.toPolyMod prime quotient) := by
    simpa [StrictHensel.toPolyMod, Polynomial.map_mul, Polynomial.map_C] using
      congrArg (Polynomial.map (Int.castRingHom (ZMod prime))) hextraction
  have hpartitionLarge :=
    StrictRecombine.removeCombination_toPolyMod_product_partition modulus
      candidate active remaining hlegalSize hremove
  have hpartitionPrime :=
    StrictRecombine.removeCombination_toPolyMod_product_partition prime
      candidate active remaining hlegalSize hremove
  let selectedLarge :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      (StrictHensel.toPolyMod modulus)).prod
  let selectedPrime :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      (StrictHensel.toPolyMod prime)).prod
  let remainingLarge :=
    (remaining.toList.map (StrictHensel.toPolyMod modulus)).prod
  let remainingPrime :=
    (remaining.toList.map (StrictHensel.toPolyMod prime)).prod
  have hactiveLarge : selectedLarge * remainingLarge =
      Polynomial.C scaleLarge * StrictHensel.toPolyMod modulus source := by
    rw [hpartitionLarge]
    exact hproductLarge
  have hactivePrime : selectedPrime * remainingPrime =
      Polynomial.C scalePrime * StrictHensel.toPolyMod prime source := by
    rw [hpartitionPrime]
    exact hproductPrime
  have hselectedLargeMonic := activeState.selectedToPolyModMonic candidate
    hbound (modulus := modulus)
  have hselectedPrimeMonic := activeState.selectedToPolyModMonic candidate
    hbound (modulus := prime)
  have hselectedLargeLeading : IsUnit selectedLarge.leadingCoeff := by
    rw [hselectedLargeMonic.leadingCoeff]
    exact isUnit_one
  have hselectedPrimeLeading : IsUnit selectedPrime.leadingCoeff := by
    rw [hselectedPrimeMonic.leadingCoeff]
    exact isUnit_one
  have hfrontPrime : (source[0]!.2 : ZMod prime) ≠ 0 := by
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source
      hcanonical hnonempty
    simpa [getElem!_pos source 0 hnonempty] using (hhead ▸ hleading)
  have hfrontPrimeUnit : IsUnit (source[0]!.2 : ZMod prime) :=
    isUnit_iff_ne_zero.mpr hfrontPrime
  have hfrontLargeUnit : IsUnit (source[0]!.2 : ZMod modulus) := by
    simpa [modulus] using
      intCast_isUnit_primePower_of_ne_zero_prime prime exponent source[0]!.2
        hfrontPrime
  rcases remaining_product_eq_unit_mul_quotient_cancel_selected
      selectedLarge remainingLarge (StrictHensel.toPolyMod modulus source)
      (StrictHensel.toPolyMod modulus factor)
      (StrictHensel.toPolyMod modulus quotient) scaleLarge
      (source[0]!.2 : ZMod modulus) (content : ZMod modulus)
      (scalar : ZMod modulus) hscaleLarge hfrontLargeUnit hcontentLarge
      hscalarLarge hselectedLargeLeading hactiveLarge hextractionLarge
      hfactorLarge with
    ⟨nextScaleLarge, hnextScaleLarge, hnextProductLarge⟩
  rcases remaining_product_eq_unit_mul_quotient_cancel_selected
      selectedPrime remainingPrime (StrictHensel.toPolyMod prime source)
      (StrictHensel.toPolyMod prime factor)
      (StrictHensel.toPolyMod prime quotient) scalePrime
      (source[0]!.2 : ZMod prime) (content : ZMod prime)
      (scalar : ZMod prime) hscalePrime hfrontPrimeUnit hcontentPrime
      hscalarPrime hselectedPrimeLeading hactivePrime hextractionPrime
      hfactorPrime with
    ⟨nextScalePrime, hnextScalePrime, hnextProductPrime⟩
  exact {
    exponentPositive := state.exponentPositive
    certificate := ⟨nextScaleLarge, nextScalePrime, hnextScaleLarge,
      hnextScalePrime, by simpa [remainingLarge] using hnextProductLarge,
      by simpa [remainingPrime] using hnextProductPrime⟩ }

/-- The full live Hensel product is preserved by the candidate physically
returned by a concrete fixed-size scan and by the subsequent literal removal.
The successful attempt equation and legality certificate are recovered by
replaying that same scan; neither is supplied independently. -/
theorem LiveHenselProduct.extractScan
    {prime exponent count : Nat} [Fact (Nat.Prime prime)]
    {source factor quotient : SparsePolyZZ}
    {active remaining : Array SparsePolyZZ} {candidate : Array Nat}
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hfits : count ≤ active.size)
    (hcandidateSize : candidate.size = count)
    (hscan : Generated.StrictRecombine.scanZassenhausCombinations
      (StrictRecombine.concreteZassenhausTermination.combinations
        active.size count)
      source active (((prime ^ exponent : Nat) : ZZ))
      (Generated.StrictRecombine.initialCombination count)
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) =
        .ok (.extracted factor quotient candidate hcandidateSize))
    (hremove : Generated.StrictRecombine.removeCombination candidate active =
      .ok remaining) :
    LiveHenselProduct prime exponent quotient remaining := by
  have hlegalCount : StrictRecombine.LegalCombination active.size count
      candidate :=
    StrictRecombine.concreteScan_extracted_legal source factor quotient active
      (((prime ^ exponent : Nat) : ZZ)) candidate hcandidateSize hfits hscan
  have hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate := by
    simpa [hcandidateSize] using hlegalCount
  have hattempt : Generated.StrictRecombine.zassenhausAttempt source active
      (((prime ^ exponent : Nat) : ZZ)) candidate =
        .ok (.extracted factor quotient) :=
    StrictRecombine.scanZassenhausCombinations_extracted_attempt
      (StrictRecombine.concreteZassenhausTermination.combinations
        active.size count)
      source factor quotient active (((prime ^ exponent : Nat) : ZZ))
      (Generated.StrictRecombine.initialCombination count) candidate
      hcandidateSize
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) hscan
  exact state.extract activeState hcanonical hnonempty hprimitive hleading
    hlegal hattempt hremove

/-- Hensel uniqueness at an arbitrary recursive Zassenhaus state.  Both
factorizations are built from the current physical active array and the two
product equations carried by `LiveHenselProduct`; this is the recursive-state
counterpart of the initial-entry uniqueness theorem. -/
theorem LiveHenselProduct.selectedCongruent
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {source : SparsePolyZZ} {active : Array SparsePolyZZ}
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly source = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod prime)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map (Int.castRingHom (ZMod prime)) divisor).leadingCoeff =
        (divisor.leadingCoeff : ZMod prime))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod prime)) divisor)
      (((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
        (StrictHensel.toPolyMod prime)).prod)) :
    Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
        (Polynomial.C (SparsePolyZZ.toPoly source).leadingCoeff *
          ((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map SparsePolyZZ.toPoly).prod) =
      Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  let modulus := prime ^ exponent
  let sourcePoly := SparsePolyZZ.toPoly source
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  rcases StrictRecombine.removeCombination_succeeds candidate active hlegal with
    ⟨complement, hremove⟩
  rcases state.certificate with
    ⟨scale, scalePrime, hscaleUnit, hscalePrimeUnit,
      hfullLarge, hfullPrime⟩
  let complementInteger :=
    (complement.toList.map SparsePolyZZ.toPoly).prod
  let liftedScaled := Polynomial.C sourcePoly.leadingCoeff * selectedInteger
  let divisorScaled := Polynomial.C quotient.leadingCoeff * divisor
  let scaleInt : Int := Int.ofNat scale.val
  let complementScaled :=
    Polynomial.C (scaleInt * divisor.leadingCoeff) * quotient
  let commonSource :=
    Polynomial.C scaleInt * Polynomial.C sourcePoly.leadingCoeff * sourcePoly
  have hprime : Nat.Prime prime := Fact.out
  have hmodulusNe : modulus ≠ 0 := pow_ne_zero exponent hprime.ne_zero
  letI : NeZero modulus := ⟨hmodulusNe⟩
  have hdivides : prime ∣ modulus :=
    dvd_pow_self prime (Nat.ne_of_gt state.exponentPositive)
  have hscaleLarge : (scaleInt : ZMod modulus) = scale := by
    simp [scaleInt, modulus]
  have hsourceLeading : sourcePoly.leadingCoeff =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    dsimp [sourcePoly]
    rw [hfactor, Polynomial.leadingCoeff_mul]
  have hselectedMonic := activeState.selectedToPolyMonic candidate hbound
  have hliftedLeading : liftedScaled.leadingCoeff = sourcePoly.leadingCoeff := by
    dsimp [liftedScaled, selectedInteger]
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hselectedMonic.leadingCoeff, mul_one]
  have hdivisorLeadingInt : divisorScaled.leadingCoeff =
      sourcePoly.leadingCoeff := by
    dsimp [divisorScaled]
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hsourceLeading]
    ring
  have hleadingEqual : liftedScaled.leadingCoeff = divisorScaled.leadingCoeff :=
    hliftedLeading.trans hdivisorLeadingInt.symm
  have hleadingCoprime : ¬((prime : Int) ∣ liftedScaled.leadingCoeff) := by
    rw [hliftedLeading]
    intro hdivisor
    apply hleading
    rw [ZMod.intCast_zmod_eq_zero_iff_dvd]
    exact hdivisor
  have hpartitionLarge :=
    StrictRecombine.removeCombination_toPolyMod_product_partition modulus
      candidate active complement hlegal hremove
  have hselectedMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) selectedInteger =
      ((StrictRecombine.selectSourceIndices active.toList
        candidate.toList).map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hproductLifted : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementInteger * liftedScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [Polynomial.map_mul, liftedScaled, commonSource,
      Polynomial.map_C]
    rw [hselectedMapLarge, hcomplementMapLarge]
    change _ * (Polynomial.C (sourcePoly.leadingCoeff : ZMod modulus) * _) = _
    calc
      _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod modulus) *
          (((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map (StrictHensel.toPolyMod modulus)).prod *
            (complement.toList.map
              (StrictHensel.toPolyMod modulus)).prod) := by ring
      _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod modulus) *
          (active.toList.map (StrictHensel.toPolyMod modulus)).prod := by
            rw [hpartitionLarge]
      _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod modulus) *
          (Polynomial.C scale * StrictHensel.toPolyMod modulus source) := by
            rw [hfullLarge]
      _ = Polynomial.C (scaleInt : ZMod modulus) *
          Polynomial.C (sourcePoly.leadingCoeff : ZMod modulus) *
          Polynomial.map (Int.castRingHom (ZMod modulus)) sourcePoly := by
            rw [hscaleLarge]
            change _ = _ * _ * StrictHensel.toPolyMod modulus source
            ring
  have hfactorMapLarge := congrArg
    (Polynomial.map (Int.castRingHom (ZMod modulus))) hfactor
  have hproductDivisor : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementScaled * divisorScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [complementScaled, divisorScaled, commonSource,
      Polynomial.map_mul, Polynomial.map_C]
    rw [Polynomial.map_mul] at hfactorMapLarge
    change (Polynomial.C ((scaleInt * divisor.leadingCoeff : Int) : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) quotient) *
      (Polynomial.C (quotient.leadingCoeff : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) divisor) = _
    rw [Int.cast_mul, Polynomial.C_mul]
    rw [hsourceLeading]
    simp only [map_mul]
    change _ = _ * Polynomial.map (Int.castRingHom (ZMod modulus)) sourcePoly
    have hfactorMapLarge' :
        Polynomial.map (Int.castRingHom (ZMod modulus)) sourcePoly =
          Polynomial.map (Int.castRingHom (ZMod modulus)) divisor *
            Polynomial.map (Int.castRingHom (ZMod modulus)) quotient := by
      simpa [sourcePoly] using hfactorMapLarge
    rw [hfactorMapLarge']
    have hcast (value : Int) :
        (Int.castRingHom (ZMod modulus)) value =
          (value : ZMod modulus) := rfl
    rw [hcast scaleInt, hcast divisor.leadingCoeff,
      hcast quotient.leadingCoeff]
    ring
  have hselectedPrimeMonic :=
    activeState.selectedToPolyModMonic candidate hbound (modulus := prime)
  have hbaseB := leading_scaled_monic_associated_divisor prime sourcePoly
    divisor quotient
    (((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      (StrictHensel.toPolyMod prime)).prod)
    (by simpa [sourcePoly] using hfactor) hdivisorModNonzero hdivisorLeading
    hselectedPrimeMonic hassociated
  have hbaseB' : Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
    simpa [liftedScaled, divisorScaled, selectedInteger,
      Polynomial.map_mul, Polynomial.map_C, Polynomial.map_list_prod,
      List.map_map, StrictHensel.toPolyMod] using hbaseB
  have hpartitionPrime :=
    StrictRecombine.removeCombination_toPolyMod_product_partition prime
      candidate active complement hlegal hremove
  have hcop : IsCoprime
      (complement.toList.map (StrictHensel.toPolyMod prime)).prod
      (Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) *
        ((StrictRecombine.selectSourceIndices active.toList
          candidate.toList).map (StrictHensel.toPolyMod prime)).prod) := by
    let selected :=
      ((StrictRecombine.selectSourceIndices active.toList
        candidate.toList).map (StrictHensel.toPolyMod prime)).prod
    let remainder :=
      (complement.toList.map (StrictHensel.toPolyMod prime)).prod
    let sourceMod :=
      Polynomial.map (Int.castRingHom (ZMod prime)) sourcePoly
    have hproduct : remainder *
          (Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) * selected) =
        Polynomial.C ((sourcePoly.leadingCoeff : ZMod prime) * scalePrime) *
          sourceMod := by
      dsimp [selected, remainder, sourceMod]
      calc
        _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) *
            (((StrictRecombine.selectSourceIndices active.toList
              candidate.toList).map (StrictHensel.toPolyMod prime)).prod *
              (complement.toList.map
                (StrictHensel.toPolyMod prime)).prod) := by ring
        _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) *
            (active.toList.map (StrictHensel.toPolyMod prime)).prod := by
              rw [hpartitionPrime]
        _ = Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) *
            (Polynomial.C scalePrime * StrictHensel.toPolyMod prime source) := by
              rw [hfullPrime]
        _ = _ := by
          have hsourcePrime : StrictHensel.toPolyMod prime source =
              Polynomial.map (Int.castRingHom (ZMod prime)) sourcePoly := rfl
          rw [hsourcePrime]
          rw [Polynomial.C_mul]
          ring
    have hscalarUnit : IsUnit
        (Polynomial.C ((sourcePoly.leadingCoeff : ZMod prime) * scalePrime)) :=
      Polynomial.isUnit_C.mpr
        ((isUnit_iff_ne_zero.mpr hleading).mul hscalePrimeUnit)
    have hassociatedProduct : Associated
        (remainder *
          (Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) * selected))
        sourceMod := by
      rw [hproduct]
      exact (associated_isUnit_mul_left_iff hscalarUnit).mpr
        (Associated.refl sourceMod)
    have hproductSquarefree : Squarefree
        (remainder *
          (Polynomial.C (sourcePoly.leadingCoeff : ZMod prime) * selected)) :=
      hassociatedProduct.squarefree_iff.mpr (by simpa [sourcePoly] using hsquarefree)
    exact (IsRelPrime.of_squarefree_mul hproductSquarefree).isCoprime
  have hselectedMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) selectedInteger =
      ((StrictRecombine.selectSourceIndices active.toList
        candidate.toList).map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcop' : IsCoprime
      (Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger)
      (Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled) := by
    simpa [liftedScaled, hselectedMapPrime, hcomplementMapPrime] using hcop
  have hproductLiftedPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementInteger * liftedScaled) commonSource
    hproductLifted
  have hproductDivisorPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementScaled * divisorScaled) commonSource
    hproductDivisor
  have hbaseBNonzero :
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled ≠ 0 := by
    rw [← hbaseB']
    simp only [liftedScaled, Polynomial.map_mul, Polynomial.map_C]
    rw [hselectedMapPrime]
    exact mul_ne_zero (Polynomial.C_ne_zero.mpr hleading)
      hselectedPrimeMonic.ne_zero
  have hbaseA :
      Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger =
        Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled := by
    have heq :
        Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger *
            Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
          Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled *
            Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
      calc
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementInteger * liftedScaled) := by
                exact (Polynomial.map_mul
                  (Int.castRingHom (ZMod prime))).symm
        _ = Polynomial.map (Int.castRingHom (ZMod prime)) commonSource :=
              hproductLiftedPrime
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementScaled * divisorScaled) := hproductDivisorPrime.symm
        _ = _ := by rw [Polynomial.map_mul]
    rw [hbaseB'] at heq
    exact mul_right_cancel₀ hbaseBNonzero heq
  have hunique := hensel_unique prime hprime exponent state.exponentPositive
    commonSource complementInteger liftedScaled complementScaled divisorScaled
    hproductLifted hproductDivisor hbaseA hbaseB' hcop' hleadingEqual
    hleadingCoprime
  simpa [modulus, liftedScaled, divisorScaled] using hunique.2

/-- The quotient installed by a successful literal attempt remains squarefree
at the selected prime whenever the current source is squarefree there.  The
divisibility witness is obtained from the attempt's actual integer extraction
equation. -/
theorem zassenhausAttempt_extracted_quotient_squarefree
    (prime : Nat) [Fact (Nat.Prime prime)]
    (source factor quotient : SparsePolyZZ)
    (active : Array SparsePolyZZ) (modulus : ZZ) (candidate : Array Nat)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (hrun : Generated.StrictRecombine.zassenhausAttempt source active modulus
      candidate = .ok (.extracted factor quotient)) :
    Squarefree (Polynomial.map (Int.castRingHom (ZMod prime))
      (SparsePolyZZ.toPoly quotient)) := by
  rcases StrictRecombine.zassenhausAttempt_extracted_unit_scalar source factor
      quotient active modulus candidate hprimitive hrun with
    ⟨scalar, _hscalar, hextraction⟩
  have hmapped := congrArg
    (Polynomial.map (Int.castRingHom (ZMod prime))) hextraction
  simp only [Polynomial.map_mul, Polynomial.map_C] at hmapped
  have hscalarCast : (Int.castRingHom (ZMod prime)) scalar =
      (scalar : ZMod prime) := rfl
  rw [hscalarCast] at hmapped
  have hdivides : Polynomial.map (Int.castRingHom (ZMod prime))
      (SparsePolyZZ.toPoly quotient) ∣
      Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source) := by
    refine ⟨Polynomial.C (scalar : ZMod prime) *
      Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly factor), ?_⟩
    rw [hmapped]
    ring
  exact hsquarefree.squarefree_of_dvd hdivides

/-- The coefficient-recovery margin survives installation of the literal
primitive quotient returned by a successful attempt. -/
theorem LiveRecoveryPrecision.extract
    {modulus : Nat} {source factor quotient : SparsePolyZZ}
    {active : Array SparsePolyZZ} {candidate : Array Nat}
    (state : LiveRecoveryPrecision modulus source)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hfactorNe : SparsePolyZZ.toPoly factor ≠ 0)
    (hquotientNe : SparsePolyZZ.toPoly quotient ≠ 0)
    (hrun : Generated.StrictRecombine.zassenhausAttempt source active
      (modulus : ZZ) candidate = .ok (.extracted factor quotient)) :
    LiveRecoveryPrecision modulus quotient := by
  rcases StrictRecombine.zassenhausAttempt_extracted_unit_scalar source factor
      quotient active (modulus : ZZ) candidate hprimitive hrun with
    ⟨scalar, hscalar, hextraction⟩
  have hscalarAbs : scalar.natAbs = 1 := Int.natAbs_of_isUnit hscalar
  let factorPoly := SparsePolyZZ.toPoly factor
  let quotientPoly := SparsePolyZZ.toPoly quotient
  have hfactorLeadingAbs : 0 < factorPoly.leadingCoeff.natAbs :=
    Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hfactorNe)
  have hquotientLeadingAbs : 0 < quotientPoly.leadingCoeff.natAbs :=
    Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hquotientNe)
  have hsourceLeading : (SparsePolyZZ.toPoly source).leadingCoeff =
      scalar * (factorPoly.leadingCoeff * quotientPoly.leadingCoeff) := by
    have hlc := congrArg Polynomial.leadingCoeff hextraction
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      Polynomial.leadingCoeff_mul] at hlc
    simpa [factorPoly, quotientPoly] using hlc
  have hquotientLeadingLe : quotientPoly.leadingCoeff.natAbs ≤
      (SparsePolyZZ.toPoly source).leadingCoeff.natAbs := by
    rw [hsourceLeading, Int.natAbs_mul, Int.natAbs_mul, hscalarAbs,
      one_mul]
    exact Nat.le_mul_of_pos_left _ hfactorLeadingAbs
  have hleadingScaled := Nat.mul_le_mul_right 2 hquotientLeadingLe
  have hleadingBound : quotientPoly.leadingCoeff.natAbs * 2 < modulus :=
    lt_of_le_of_lt hleadingScaled state.leadingBound
  refine ⟨by simpa [quotientPoly] using hleadingBound, ?_⟩
  intro divisor nextQuotient hnext hdivisorNe degree
  have hnextQuotientNe : nextQuotient ≠ 0 := by
    intro hzero
    apply hquotientNe
    change SparsePolyZZ.toPoly quotient = divisor * nextQuotient at hnext
    rw [hnext, hzero, mul_zero]
  let sourceQuotient := Polynomial.C scalar * (factorPoly * nextQuotient)
  have hsourceFactorization : SparsePolyZZ.toPoly source =
      divisor * sourceQuotient := by
    dsimp [sourceQuotient]
    rw [hextraction, hnext]
    change Polynomial.C scalar *
      (factorPoly * (divisor * nextQuotient)) = _
    ring
  have holdBound := state.scaledFactorBound divisor sourceQuotient
    hsourceFactorization hdivisorNe degree
  have hsourceQuotientLeading : sourceQuotient.leadingCoeff =
      scalar * (factorPoly.leadingCoeff * nextQuotient.leadingCoeff) := by
    dsimp [sourceQuotient]
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      Polynomial.leadingCoeff_mul]
  have htargetCoeff :
      (Polynomial.C sourceQuotient.leadingCoeff * divisor).coeff degree =
        scalar * factorPoly.leadingCoeff *
          ((Polynomial.C nextQuotient.leadingCoeff * divisor).coeff degree) := by
    rw [hsourceQuotientLeading]
    change (Polynomial.C
      (scalar * (factorPoly.leadingCoeff * nextQuotient.leadingCoeff)) *
        divisor).coeff degree = _
    rw [Polynomial.coeff_C_mul]
    rw [Polynomial.coeff_C_mul]
    ring
  rw [htargetCoeff, Int.natAbs_mul, Int.natAbs_mul,
    hscalarAbs, one_mul] at holdBound
  have hdesiredLe :
      ((Polynomial.C nextQuotient.leadingCoeff * divisor).coeff degree).natAbs *
          2 ≤
        (factorPoly.leadingCoeff.natAbs *
          ((Polynomial.C nextQuotient.leadingCoeff * divisor).coeff degree).natAbs) *
            2 := by
    have hbase := Nat.le_mul_of_pos_left
      ((Polynomial.C nextQuotient.leadingCoeff * divisor).coeff degree).natAbs
      hfactorLeadingAbs
    exact Nat.mul_le_mul_right 2 hbase
  exact lt_of_le_of_lt hdesiredLe (by
    simpa [Nat.mul_assoc] using holdBound)

/-- A hypothetical nontrivial factorization of the literal extracted factor
is re-encoded as a smaller physical candidate together with every premise
needed by the later generated-attempt execution theorem. -/
theorem zassenhausAttempt_reducible_has_smaller_candidate
    (prime exponent : Nat) [Fact (Nat.Prime prime)]
    (source factor quotientSparse : SparsePolyZZ)
    (active : Array SparsePolyZZ) (outer : Array Nat)
    (subsetSize : Nat)
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (houter : StrictRecombine.LegalCombination active.size outer.size outer)
    (houterSize : outer.size = subsetSize)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hrecovery : ∀ divisor quotient,
      SparsePolyZZ.toPoly source = divisor * quotient → divisor ≠ 0 →
      ∀ degree,
        ((Polynomial.C quotient.leadingCoeff * divisor).coeff degree).natAbs *
          2 < prime ^ exponent)
    (left right : Polynomial Int)
    (hfactorization : SparsePolyZZ.toPoly factor = left * right)
    (hleftUnit : ¬IsUnit left) (hrightUnit : ¬IsUnit right)
    (hattempt : Generated.StrictRecombine.zassenhausAttempt source active
      (((prime ^ exponent : Nat) : ZZ)) outer =
        .ok (.extracted factor quotientSparse)) :
    ∃ inner : Array Nat, ∃ innerQuotient : Polynomial Int,
      StrictRecombine.LegalCombination active.size inner.size inner ∧
      0 < inner.size ∧ inner.size < subsetSize ∧
      left.IsPrimitive ∧ left ≠ 0 ∧ innerQuotient ≠ 0 ∧
      SparsePolyZZ.toPoly source = left * innerQuotient ∧
      Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
          (Polynomial.C (SparsePolyZZ.toPoly source).leadingCoeff *
            ((StrictRecombine.selectSourceIndices active.toList
              inner.toList).map SparsePolyZZ.toPoly).prod) =
        Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
          (Polynomial.C innerQuotient.leadingCoeff * left) ∧
      ∀ degree,
        ((Polynomial.C innerQuotient.leadingCoeff * left).coeff degree).natAbs *
          2 < prime ^ exponent := by
  let modulus := prime ^ exponent
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime prime).pos exponent
  have hdivides : prime ∣ modulus :=
    dvd_pow_self prime (Nat.ne_of_gt state.exponentPositive)
  have hfactorProperties :=
    StrictRecombine.zassenhausAttempt_extracted_factor_canonical_primitive
      source factor quotientSparse active (((modulus : Nat) : ZZ)) outer
      hcanonical (by simpa [modulus] using hattempt)
  have hfactorAssociated :=
    StrictRecombine.zassenhausAttempt_extracted_factor_mod_associated_selected
      source factor quotientSparse active modulus prime outer hmodulus
      (Fact.out : Nat.Prime prime).pos hdivides
      (fun position hposition => by
        simpa [getElem!_pos outer position hposition] using
          houter.2.2 position hposition)
      activeState.fitsInt32
      (by
        have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source
          hcanonical hnonempty
        simpa [getElem!_pos source 0 hnonempty] using (hhead ▸ hleading))
      activeState.irreducible (by simpa [modulus] using hattempt)
  rcases StrictRecombine.zassenhausAttempt_extracted_unit_scalar source factor
      quotientSparse active (((modulus : Nat) : ZZ)) outer hprimitive
      (by simpa [modulus] using hattempt) with
    ⟨scalar, hscalar, hextraction⟩
  have hfactorLeading :
      ((SparsePolyZZ.toPoly factor).leadingCoeff : ZMod prime) ≠ 0 := by
    have hlc := congrArg Polynomial.leadingCoeff hextraction
    simp only [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C] at hlc
    have hlcMod := congrArg (fun value : Int => (value : ZMod prime)) hlc
    simp only [Int.cast_mul] at hlcMod
    intro hzero
    apply hleading
    rw [hlcMod, hzero]
    simp
  rcases smaller_active_candidate_of_reducible_primitive_factor prime active
      outer houter activeState.irreducible (SparsePolyZZ.toPoly factor) left
      right hfactorProperties.2 hfactorLeading hfactorization hleftUnit
      hrightUnit hfactorAssociated with
    ⟨inner, hinner, hinnerPositive, hinnerSmall, hleftPrimitive,
      hleftLeading, hinnerAssociated⟩
  let quotientPoly := SparsePolyZZ.toPoly quotientSparse
  let innerQuotient := Polynomial.C scalar * (right * quotientPoly)
  have hsourceFactorization : SparsePolyZZ.toPoly source =
      left * innerQuotient := by
    dsimp [innerQuotient, quotientPoly]
    rw [hextraction, hfactorization]
    ring
  have hleftNe : left ≠ 0 := by
    intro hzero
    apply hleftLeading
    simp [hzero]
  have hinnerQuotientNe : innerQuotient ≠ 0 := by
    intro hzero
    have hsourceZero : SparsePolyZZ.toPoly source = 0 := by
      rw [hsourceFactorization, hzero, mul_zero]
    apply hleading
    rw [hsourceZero]
    simp
  have hleftMapNe : Polynomial.map (Int.castRingHom (ZMod prime)) left ≠ 0 := by
    intro hzero
    apply hleftLeading
    have hdegree := Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero
      (Int.castRingHom (ZMod prime)) hleftLeading
    rw [hzero] at hdegree
    simpa using hdegree.symm
  have hleftLeadingMap :
      (Polynomial.map (Int.castRingHom (ZMod prime)) left).leadingCoeff =
        (left.leadingCoeff : ZMod prime) :=
    Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ hleftLeading
  have hcongruent := state.selectedCongruent activeState hsquarefree hleading
    left innerQuotient hsourceFactorization hleftMapNe hleftLeadingMap inner
    hinner hinnerAssociated
  have htargetBound := hrecovery left innerQuotient hsourceFactorization
    hleftNe
  exact ⟨inner, innerQuotient, hinner, hinnerPositive,
    by simpa [houterSize] using hinnerSmall, hleftPrimitive, hleftNe,
    hinnerQuotientNe, hsourceFactorization, by simpa [modulus] using hcongruent,
    by simpa [modulus] using htargetBound⟩

/-- The normalized factors returned by the literal Hensel entry initialize a
`LiveHenselProduct` at the exact prime power returned by that execution. -/
theorem selectionHenselFactors_liveProduct
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ exponent : Nat,
      output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
        LiveHenselProduct selection.prime.toNat exponent f output.1 := by
  rcases selectionHenselFactors_prime_product_eq_unit_mul_source hcount
      hfactors hleadingSemantic hselection hentry with
    ⟨exponent, hexponent, scale, scaleAtPrime, houtput, hscaleUnit,
      hscaleAtPrimeUnit, hscaleAtPrime, hlargeProduct, hprimeProduct⟩
  exact ⟨exponent, houtput, {
    exponentPositive := hexponent
    certificate := ⟨scale, scaleAtPrime, hscaleUnit, hscaleAtPrimeUnit,
      hlargeProduct, hprimeProduct⟩ }⟩

/-- The default generated Mignotte/Hensel execution initializes the exact
recovery margin required by every later Zassenhaus candidate. -/
theorem selectionHenselFactors_liveRecoveryPrecision
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime p.toNat)]
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p 0
      output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (exponent : Nat)
    (houtput : output.2 = ((p.toNat ^ exponent : Nat) : Int)) :
    LiveRecoveryPrecision (p.toNat ^ exponent) f := by
  let modulus := p.toNat ^ exponent
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
      hcanonical hnonempty
    rw [hzero] at hhead
    exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
      (by simpa using hhead.symm)
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have honePrecision :=
    StrictRecombine.hensel_output_modulus_bounds_scaled_divisor hentry rfl
      hcanonical hnonempty hdegree leading hleading (1 : Polynomial Int)
      hsourceNe (one_dvd _)
  have hleadingLargeInt := honePrecision.2 0
  rw [houtput] at hleadingLargeInt
  have hleadingLargeInt' : (leading.2.natAbs : Int) * 2 <
      (modulus : Int) := by
    simpa [modulus] using hleadingLargeInt
  have hleadingBound : (SparsePolyZZ.toPoly f).leadingCoeff.natAbs * 2 <
      modulus := by
    rw [hsourceLeading]
    exact_mod_cast hleadingLargeInt'
  refine ⟨hleadingBound, ?_⟩
  intro divisor quotient hfactor hdivisorNe degree
  have hprecision :=
    StrictRecombine.hensel_output_modulus_bounds_scaled_divisor hentry rfl
      hcanonical hnonempty hdegree leading hleading divisor hsourceNe
      (hfactor ▸ dvd_mul_right divisor quotient)
  have hlargeInt := hprecision.2 degree
  rw [houtput] at hlargeInt
  have hlarge : (leading.2 * divisor.coeff degree).natAbs * 2 < modulus := by
    exact_mod_cast hlargeInt
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have hdivisorLeadingAbs : 0 < divisor.leadingCoeff.natAbs :=
    Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe)
  have hle : (quotient.leadingCoeff * divisor.coeff degree).natAbs * 2 ≤
      (leading.2 * divisor.coeff degree).natAbs * 2 := by
    rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_mul]
    have hbase := Nat.le_mul_of_pos_left
      (quotient.leadingCoeff.natAbs * (divisor.coeff degree).natAbs)
      hdivisorLeadingAbs
    have hscaled := Nat.mul_le_mul_right 2 hbase
    simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
  have htargetCoeff :
      (Polynomial.C quotient.leadingCoeff * divisor).coeff degree =
        quotient.leadingCoeff * divisor.coeff degree := by simp
  rw [htargetCoeff]
  exact lt_of_le_of_lt hle hlarge

/-- Any concrete modulus that exceeds the Mignotte target computed by the
generated helper supplies the same live coefficient-recovery margin.  This
form is used for the positive explicit-exponent branch of `__lll_factorize`.-/
theorem liveRecoveryPrecision_of_generated_mignotte_lt
    (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (modulus : Nat)
    (hlarge : ∀ bound,
      Generated.StrictHensel.__mignotte_bound_upoly_raw_ir f = .ok bound →
      (2 * (if leading.2 < 0 then -leading.2 else leading.2) * bound).natAbs <
        modulus) :
    LiveRecoveryPrecision modulus f := by
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
      hcanonical hnonempty
    rw [hzero] at hhead
    exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
      (by simpa using hhead.symm)
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  rcases StrictRecombine.mignotteBoundRaw_bounds_divisor f hcanonical
      hnonempty hdegree (1 : Polynomial Int) hsourceNe (one_dvd _) with
    ⟨bound, hboundRun, hboundNonnegative, hboundOne⟩
  have hboundPositive : 1 ≤ bound := by
    have hone := hboundOne 0
    simpa using hone
  have htarget := hlarge bound hboundRun
  have htargetForm : leading.2.natAbs * 2 * bound.natAbs < modulus := by
    have habs : (if leading.2 < 0 then -leading.2 else leading.2).natAbs =
        leading.2.natAbs := by
      split <;> simp_all [Int.natAbs_neg]
    rw [Int.natAbs_mul, Int.natAbs_mul, habs] at htarget
    norm_num at htarget
    simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using htarget
  have hleadingBound :
      (SparsePolyZZ.toPoly f).leadingCoeff.natAbs * 2 < modulus := by
    rw [hsourceLeading]
    have hscale : leading.2.natAbs * 2 ≤
        leading.2.natAbs * 2 * bound.natAbs := by
      have hboundNat : 1 ≤ bound.natAbs := by
        have habsInt : (bound.natAbs : Int) = bound := by
          exact Int.natAbs_of_nonneg hboundNonnegative
        rw [← habsInt] at hboundPositive
        exact_mod_cast hboundPositive
      simpa using Nat.mul_le_mul_left (leading.2.natAbs * 2) hboundNat
    exact lt_of_le_of_lt hscale htargetForm
  refine ⟨hleadingBound, ?_⟩
  intro divisor quotient hfactor hdivisorNe degree
  rcases StrictRecombine.mignotteBoundRaw_bounds_divisor f hcanonical
      hnonempty hdegree divisor hsourceNe
      (hfactor ▸ dvd_mul_right divisor quotient) with
    ⟨divisorBound, hdivisorBoundRun, hdivisorBoundNonnegative,
      hdivisorBound⟩
  have hsameBound : divisorBound = bound := by
    rw [hboundRun] at hdivisorBoundRun
    exact (Except.ok.inj hdivisorBoundRun).symm
  subst divisorBound
  have hcoefficient := hdivisorBound degree
  have hleadingScaled :
      (leading.2 * divisor.coeff degree).natAbs * 2 < modulus := by
    rw [Int.natAbs_mul]
    have hcoefficientNat : (divisor.coeff degree).natAbs ≤ bound.natAbs := by
      have habsInt : (bound.natAbs : Int) = bound := by
        exact Int.natAbs_of_nonneg hboundNonnegative
      rw [← habsInt] at hcoefficient
      exact_mod_cast hcoefficient
    have hle : leading.2.natAbs * (divisor.coeff degree).natAbs * 2 ≤
        leading.2.natAbs * 2 * bound.natAbs := by
      simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using
        Nat.mul_le_mul_left (leading.2.natAbs * 2) hcoefficientNat
    exact lt_of_le_of_lt hle htargetForm
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have hdivisorLeadingAbs : 0 < divisor.leadingCoeff.natAbs :=
    Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe)
  have hle : (quotient.leadingCoeff * divisor.coeff degree).natAbs * 2 ≤
      (leading.2 * divisor.coeff degree).natAbs * 2 := by
    rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul, Int.natAbs_mul]
    have hbase := Nat.le_mul_of_pos_left
      (quotient.leadingCoeff.natAbs * (divisor.coeff degree).natAbs)
      hdivisorLeadingAbs
    have hscaled := Nat.mul_le_mul_right 2 hbase
    simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
  have htargetCoeff :
      (Polynomial.C quotient.leadingCoeff * divisor).coeff degree =
        quotient.leadingCoeff * divisor.coeff degree := by simp
  rw [htargetCoeff]
  exact lt_of_le_of_lt hle hleadingScaled

/-- The exact selected product and the exact physical complement computed by
`removeCombination` are coprime at the selected prime after the source-leading
scaling used by Zassenhaus.  Their product differs from the squarefree source
only by the two concrete units supplied by prime selection and Hensel
normalization. -/
theorem henselCandidate_physicalComplement_coprime
    (p : Nat) [Fact (Nat.Prime p)]
    (source : Polynomial Int) (active complement : Array SparsePolyZZ)
    (candidate : Array Nat) (scale : ZMod p)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hremove : Generated.StrictRecombine.removeCombination candidate active =
      .ok complement)
    (hfull : (active.toList.map (StrictHensel.toPolyMod p)).prod =
      Polynomial.C scale *
        Polynomial.map (Int.castRingHom (ZMod p)) source)
    (hscaleUnit : IsUnit scale)
    (hleadingNonzero : (source.leadingCoeff : ZMod p) ≠ 0)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod p)) source)) :
    IsCoprime
      (complement.toList.map (StrictHensel.toPolyMod p)).prod
      (Polynomial.C (source.leadingCoeff : ZMod p) *
        ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
          (StrictHensel.toPolyMod p)).prod) := by
  let selected :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      (StrictHensel.toPolyMod p)).prod
  let remainder :=
    (complement.toList.map (StrictHensel.toPolyMod p)).prod
  let sourceMod := Polynomial.map (Int.castRingHom (ZMod p)) source
  have hpartition := StrictRecombine.removeCombination_toPolyMod_product_partition
    p candidate active complement hlegal hremove
  have hproduct : remainder *
        (Polynomial.C (source.leadingCoeff : ZMod p) * selected) =
      Polynomial.C ((source.leadingCoeff : ZMod p) * scale) * sourceMod := by
    dsimp [selected, remainder, sourceMod]
    calc
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map (StrictHensel.toPolyMod p)).prod *
            (complement.toList.map (StrictHensel.toPolyMod p)).prod) := by
              ring
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (active.toList.map (StrictHensel.toPolyMod p)).prod := by
            rw [hpartition]
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (Polynomial.C scale *
            Polynomial.map (Int.castRingHom (ZMod p)) source) := by
              rw [hfull]
      _ = _ := by rw [Polynomial.C_mul]; ring
  have hleadingUnit : IsUnit (source.leadingCoeff : ZMod p) :=
    isUnit_iff_ne_zero.mpr hleadingNonzero
  have hscalarUnit : IsUnit
      (Polynomial.C ((source.leadingCoeff : ZMod p) * scale)) :=
    Polynomial.isUnit_C.mpr (hleadingUnit.mul hscaleUnit)
  have hassociated : Associated
      (remainder * (Polynomial.C (source.leadingCoeff : ZMod p) * selected))
      sourceMod := by
    rw [hproduct]
    exact (associated_isUnit_mul_left_iff hscalarUnit).mpr
      (Associated.refl sourceMod)
  have hproductSquarefree : Squarefree
      (remainder *
        (Polynomial.C (source.leadingCoeff : ZMod p) * selected)) :=
    hassociated.squarefree_iff.mpr hsquarefree
  exact (IsRelPrime.of_squarefree_mul hproductSquarefree).isCoprime

private theorem origins_preserve_irreducible
    {p : Nat} {inputs : List SparsePolyZp} {outputs : List SparsePolyZZ}
    (horigins : List.Forall₂
      (fun input output => StrictHensel.toPolyMod p output =
        SparsePolyZp.toPoly p input) inputs outputs)
    (hirreducible : ∀ input ∈ inputs,
      Irreducible (SparsePolyZp.toPoly p input)) :
    ∀ output ∈ outputs, Irreducible (StrictHensel.toPolyMod p output) := by
  induction horigins with
  | nil => simp
  | @cons input output inputs outputs horigin horigins ih =>
      intro candidate hcandidate
      rcases List.mem_cons.mp hcandidate with rfl | htail
      · rw [horigin]
        exact hirreducible input (List.mem_cons_self)
      · exact ih
          (fun factor hfactor =>
            hirreducible factor (List.mem_cons_of_mem input hfactor))
          candidate htail

/-- The factors returned by the actual generated Hensel entry remain
irreducible after reduction modulo the selected prime.  The proof follows the
concrete adjustment, tree leaves, lift extraction, and final normalization in
source order. -/
theorem selectionHenselFactors_mod_irreducible
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∀ index (hindex : index < output.1.size),
      Irreducible (StrictHensel.toPolyMod selection.prime.toNat
        output.1[index]) := by
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjusted := selectionFactors_adjusted_irreducible
    hfactors hleadingSemantic hselection hadjust
  have hadjustSize : adjusted.size = selection.factors.size := by
    cases hadjust with
    | adjusted leading first adjusted hsource hfirst hadjustedEq =>
        have hzero : 0 < selection.factors.size := by omega
        simp [Array.set!, hzero]
  have hfull : StrictHensel.henselFactorRangeList adjusted
      selection.factors.size 0 = adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have hadjustedList : ∀ factor ∈ adjusted.toList,
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat factor) := by
    intro factor hfactor
    obtain ⟨index, hindex, hfactorEq⟩ := List.mem_iff_getElem.mp hfactor
    subst factor
    simpa using hadjusted index (by simpa using hindex)
  have hextractedList := origins_preserve_irreducible horigins hadjustedList
  have hextracted : ∀ index (hindex : index < extracted.size),
      Irreducible (StrictHensel.toPolyMod selection.prime.toNat
        extracted[index]) := by
    intro index hindex
    exact hextractedList extracted[index]
      (Array.getElem_mem_toList hindex)
  exact hnormalizeRel.irreducible hextracted

/-- The normalized array returned by the real Hensel entry establishes every
physical active-factor invariant needed by later live Zassenhaus attempts.
The size bound is the concrete bound already proved at the entry call site. -/
theorem selectionHenselFactors_liveActive
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (hfits : output.1.size ≤ 2 ^ 31) :
    StrictRecombine.LiveActiveFactors selection.prime.toNat output.1 := by
  refine {
    fitsInt32 := hfits
    canonical := hentry.outputCanonical
    nonempty := ?_
    monic := hentry.outputToPolyMonic
    irreducible := selectionHenselFactors_mod_irreducible hcount hfactors
      hleadingSemantic hselection hentry }
  intro index hindex
  rcases hentry.outputOneHead index hindex with ⟨head, tail, hlist⟩
  have hlength := congrArg List.length hlist
  simp at hlength
  omega

/-- The actual Hensel output array is pointwise associated, in the same
source order, with the factor array returned by the concrete selected-prime
execution.  The association is assembled from the generated first-factor
adjustment, exact lifted-leaf origins, and generated final normalization. -/
theorem selectionHenselFactors_pointwise_associated
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    selection.factors.size = output.1.size ∧
      ∀ index (hinput : index < selection.factors.size)
        (houtput : index < output.1.size),
        Associated
          (SparsePolyZp.toPoly selection.prime.toNat
            selection.factors[index])
          (StrictHensel.toPolyMod selection.prime.toNat output.1[index]) := by
  have hleadingNonzero : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) ≠ 0 := by
    intro leading hleading
    rw [hleadingSemantic leading hleading]
    exact hselection.goodPrime.lc_nonzero
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjustRel := hadjust.unitRel hfactors hleadingNonzero
  have hadjustSize : adjusted.size = selection.factors.size := hadjustRel.1.symm
  have hfull : StrictHensel.henselFactorRangeList adjusted
      selection.factors.size 0 = adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have horiginSize : adjusted.size = extracted.size := by
    simpa using horigins.length_eq
  have hsize : selection.factors.size = output.1.size := by
    rw [hadjustRel.1, horiginSize, hnormalizeRel.1]
  refine ⟨hsize, ?_⟩
  intro index hinput houtput
  have hadjusted : index < adjusted.size := by
    rw [← hadjustRel.1]
    exact hinput
  have hextracted : index < extracted.size := by omega
  rcases hadjustRel.2 index hinput hadjusted with
    ⟨adjustScale, hadjustScale, hadjustEq⟩
  have horigin := horigins.get (by simpa using hadjusted)
    (by simpa using hextracted)
  have horiginEq : StrictHensel.toPolyMod selection.prime.toNat
      extracted[index] = SparsePolyZp.toPoly selection.prime.toNat
        adjusted[index] := by
    simpa [Array.getElem_toList] using horigin
  rcases hnormalizeRel.2 index hextracted houtput with
    ⟨normalizeScale, hnormalizeScale, hnormalizeEq⟩
  have hadjustAssociated : Associated
      (SparsePolyZp.toPoly selection.prime.toNat selection.factors[index])
      (SparsePolyZp.toPoly selection.prime.toNat adjusted[index]) := by
    rw [hadjustEq]
    exact associated_unit_mul_right _ _
      (Polynomial.isUnit_C.mpr hadjustScale)
  have horiginAssociated : Associated
      (SparsePolyZp.toPoly selection.prime.toNat adjusted[index])
      (StrictHensel.toPolyMod selection.prime.toNat extracted[index]) :=
    Associated.of_eq horiginEq.symm
  have hnormalizeAssociated : Associated
      (StrictHensel.toPolyMod selection.prime.toNat extracted[index])
      (StrictHensel.toPolyMod selection.prime.toNat output.1[index]) := by
    rw [hnormalizeEq]
    exact associated_unit_mul_right _ _
      (Polynomial.isUnit_C.mpr hnormalizeScale)
  exact hadjustAssociated.trans
    (horiginAssociated.trans hnormalizeAssociated)

private theorem associated_prod_of_forall₂ {R : Type*} [CommMonoid R]
    {left right : List R} (hrel : List.Forall₂ Associated left right) :
    Associated left.prod right.prod := by
  induction hrel with
  | nil => exact Associated.refl 1
  | cons hhead htail ih =>
      simpa only [List.prod_cons] using hhead.mul_mul ih

/-- Product form of `selectionHenselFactors_pointwise_associated`: the
modular product of the concrete selected-prime array is associated with the
modular product of the concrete normalized Hensel output array. -/
theorem selectionHenselFactors_product_associated
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    Associated
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors).prod
      ((output.1.toList.map
        (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  have hpointwise := selectionHenselFactors_pointwise_associated hcount
    hfactors hleadingSemantic hselection hentry
  apply associated_prod_of_forall₂
  rw [List.forall₂_iff_get]
  constructor
  · simpa [StrictSelectPrime.factorArrayToL2] using hpointwise.1
  · intro index hleft hright
    have hinput : index < selection.factors.size := by
      simpa [StrictSelectPrime.factorArrayToL2] using hleft
    have houtput : index < output.1.size := by simpa using hright
    simpa [StrictSelectPrime.factorArrayToL2, Array.getElem_toList] using
      hpointwise.2 index hinput houtput

/-- Every genuine integer divisor reduces, at the actually selected prime,
to the product of an occurrence-sensitive sublist of the actual normalized
Hensel output array, up to a modular unit.  The sublist comes from recursive
irreducible divisibility over that concrete finite list; no semantic factor
oracle chooses it. -/
theorem integer_divisor_mod_associated_hensel_sublist
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f) :
    ∃ chosen : List SparsePolyZZ,
      chosen.Sublist output.1.toList ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        ((chosen.map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  let atoms := output.1.toList.map
    (StrictHensel.toPolyMod selection.prime.toNat)
  have hirreducible : ∀ atom ∈ atoms, Irreducible atom := by
    intro atom hatom
    rcases List.mem_map.mp hatom with ⟨lifted, hlifted, rfl⟩
    rcases List.mem_iff_getElem.mp hlifted with ⟨index, hindex, rfl⟩
    have hindexArray : index < output.1.size := by simpa using hindex
    simpa [Array.getElem_toList] using
      selectionHenselFactors_mod_irreducible hcount hfactors
        hleadingSemantic hselection hentry index hindexArray
  have hmapDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣
      Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f) := by
    rcases hg with ⟨quotient, hfactor⟩
    refine ⟨Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) quotient, ?_⟩
    rw [hfactor, Polynomial.map_mul]
  have hselectedDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors).prod :=
    dvd_trans hmapDvd hselection.productAssociated.dvd
  have hproductAssociated := selectionHenselFactors_product_associated hcount
    hfactors hleadingSemantic hselection hentry
  have hatomsDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣ atoms.prod :=
    dvd_trans hselectedDvd hproductAssociated.dvd
  rcases divisor_associated_sublist_product atoms hirreducible _ hatomsDvd with
    ⟨chosenMod, hchosenMod, hassociated⟩
  rcases List.sublist_map_iff.mp hchosenMod with
    ⟨chosen, hchosen, rfl⟩
  exact ⟨chosen, hchosen, hassociated⟩

/-- Live-state form of the divisor-to-combination theorem.  It applies after
successful Zassenhaus removals: if the current integer remainder reduces to
the product of the current physical active array, every genuine divisor of
that remainder is represented by a legal occurrence-sensitive combination
of that same array. -/
theorem live_divisor_mod_has_legal_candidate
    (base : Nat) [Fact (Nat.Prime base)]
    (fStar : SparsePolyZZ) (active : Array SparsePolyZZ)
    (hproduct : Associated
      (Polynomial.map (Int.castRingHom (ZMod base))
        (SparsePolyZZ.toPoly fStar))
      ((active.toList.map (StrictHensel.toPolyMod base)).prod))
    (hirreducible : ∀ index (hindex : index < active.size),
      Irreducible (StrictHensel.toPolyMod base active[index]))
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly fStar) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination active.size indices.size indices ∧
      Associated
        (Polynomial.map (Int.castRingHom (ZMod base)) g)
        (((StrictRecombine.selectSourceIndices active.toList indices.toList).map
          (StrictHensel.toPolyMod base)).prod) := by
  let atoms := active.toList.map (StrictHensel.toPolyMod base)
  have hirreducibleAtoms : ∀ atom ∈ atoms, Irreducible atom := by
    intro atom hatom
    rcases List.mem_map.mp hatom with ⟨lifted, hlifted, rfl⟩
    rcases List.mem_iff_getElem.mp hlifted with ⟨index, hindex, rfl⟩
    have hindexArray : index < active.size := by simpa using hindex
    simpa [Array.getElem_toList] using hirreducible index hindexArray
  have hmapDvd : Polynomial.map (Int.castRingHom (ZMod base)) g ∣
      Polynomial.map (Int.castRingHom (ZMod base))
        (SparsePolyZZ.toPoly fStar) := by
    rcases hg with ⟨quotient, hfactor⟩
    refine ⟨Polynomial.map (Int.castRingHom (ZMod base)) quotient, ?_⟩
    rw [hfactor, Polynomial.map_mul]
  have hatomsDvd : Polynomial.map (Int.castRingHom (ZMod base)) g ∣
      atoms.prod := dvd_trans hmapDvd hproduct.dvd
  rcases divisor_associated_sublist_product atoms hirreducibleAtoms _
      hatomsDvd with ⟨chosenMod, hchosenMod, hassociated⟩
  rcases List.sublist_map_iff.mp hchosenMod with ⟨chosen, hchosen, rfl⟩
  rcases StrictRecombine.sublist_exists_legal_combination hchosen with
    ⟨indices, hlegal, hselected⟩
  have hlegal' : StrictRecombine.LegalCombination active.size
      indices.size indices := by
    simpa [hlegal.1] using hlegal
  refine ⟨indices, hlegal', ?_⟩
  rw [hselected]
  exact hassociated

/-- Array-index form consumed by the generated Zassenhaus scanner.  Every
genuine divisor supplies a legal strictly increasing candidate over the
actual Hensel output array, and the selected occurrence-sensitive product is
associated with that divisor modulo the selected prime. -/
theorem integer_divisor_mod_has_legal_hensel_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        (((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  rcases integer_divisor_mod_associated_hensel_sublist hcount hfactors
      hleadingSemantic hselection hentry g hg with
    ⟨chosen, hchosen, hassociated⟩
  rcases StrictRecombine.sublist_exists_legal_combination hchosen with
    ⟨indices, hlegal, hselected⟩
  have hlegal' : StrictRecombine.LegalCombination output.1.size
      indices.size indices := by
    simpa [hlegal.1] using hlegal
  refine ⟨indices, hlegal', ?_⟩
  rw [hselected]
  exact hassociated

/-- Every occurrence-sensitive subproduct selected from the actual normalized
Hensel output is monic.  The proof reads the literal coefficient-one head and
canonical sparse representation carried by `HenselLiftEntryCorrect` for each
selected array cell. -/
theorem henselSelectedProduct_monic
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (modulus : Nat) (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
      (StrictHensel.toPolyMod modulus)).prod.Monic := by
  apply monic_list_product
  intro candidate hcandidate
  rcases List.mem_map.mp hcandidate with ⟨lifted, hlifted, rfl⟩
  unfold StrictRecombine.selectSourceIndices at hlifted
  rcases List.mem_map.mp hlifted with ⟨index, hindex, rfl⟩
  rcases List.mem_iff_getElem.mp hindex with
    ⟨position, hposition, hindexEq⟩
  have hpositionArray : position < indices.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  have hselected : indices[position] = index := by
    rw [← Array.getElem_toList hpositionArray]
    exact hindexEq
  rw [← hselected, getElem!_pos output.1.toList indices[position]
    (by simpa using hactive), Array.getElem_toList hactive]
  exact hentry.outputToPolyModMonic modulus indices[position] hactive

/-- Integer-polynomial version used by the concrete scalar-pruning loops. -/
theorem henselSelectedIntegerProduct_monic
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
      SparsePolyZZ.toPoly).prod.Monic := by
  apply monic_list_product
  intro candidate hcandidate
  rcases List.mem_map.mp hcandidate with ⟨lifted, hlifted, rfl⟩
  unfold StrictRecombine.selectSourceIndices at hlifted
  rcases List.mem_map.mp hlifted with ⟨index, hindex, rfl⟩
  rcases List.mem_iff_getElem.mp hindex with
    ⟨position, hposition, hindexEq⟩
  have hpositionArray : position < indices.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  have hselected : indices[position] = index := by
    rw [← Array.getElem_toList hpositionArray]
    exact hindexEq
  rw [← hselected, getElem!_pos output.1.toList indices[position]
    (by simpa using hactive), Array.getElem_toList hactive]
  exact hentry.outputToPolyMonic indices[position] hactive

/-- At the selected prime, the exact integer trial polynomial built from a
legal physical Hensel candidate is the quotient-leading-scaled genuine
divisor.  This is the base congruence supplied to `hensel_unique`. -/
theorem henselCandidate_scaled_eq_divisor_mod_prime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod p.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod p.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod p.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod p.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod p.toNat)).prod)) :
    Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
          ((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map SparsePolyZZ.toPoly).prod) =
      Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hmonic := henselSelectedProduct_monic hentry p.toNat candidate hbound
  have hscaled := leading_scaled_monic_associated_divisor p.toNat
    (SparsePolyZZ.toPoly f) divisor quotient
    (((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      (StrictHensel.toPolyMod p.toNat)).prod)
    hfactor hdivisorModNonzero hdivisorLeading hmonic hassociated
  simpa [Polynomial.map_mul, Polynomial.map_C,
    Polynomial.map_list_prod, List.map_map, StrictHensel.toPolyMod] using
      hscaled

/-- The two scaled factors compared by Hensel uniqueness have literally equal
integer leading coefficients, and that coefficient is prime to the selected
prime.  The first equality uses physical candidate monicity; the second is
exactly the `GoodPrime` leading-coefficient guard. -/
theorem henselCandidate_scaled_leadingCoeff
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate) :
    let selectedInteger :=
      ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod
    let liftedScaled :=
      Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff * selectedInteger
    let divisorScaled := Polynomial.C quotient.leadingCoeff * divisor
    liftedScaled.leadingCoeff = divisorScaled.leadingCoeff ∧
      ¬((selection.prime.toNat : Int) ∣ liftedScaled.leadingCoeff) := by
  dsimp
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hmonic := henselSelectedIntegerProduct_monic hentry candidate hbound
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [hfactor, Polynomial.leadingCoeff_mul]
  have hlifted :
      (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
        ((StrictRecombine.selectSourceIndices output.1.toList
          candidate.toList).map SparsePolyZZ.toPoly).prod).leadingCoeff =
        (SparsePolyZZ.toPoly f).leadingCoeff := by
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hmonic.leadingCoeff, mul_one]
  have hdivisor :
      (Polynomial.C quotient.leadingCoeff * divisor).leadingCoeff =
        (SparsePolyZZ.toPoly f).leadingCoeff := by
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hsourceLeading]
    ring
  refine ⟨hlifted.trans hdivisor.symm, ?_⟩
  rw [hlifted]
  intro hdivides
  apply hselection.goodPrime.lc_nonzero
  rw [ZMod.intCast_zmod_eq_zero_iff_dvd]
  exact hdivides

/-- Hensel uniqueness for the concrete occurrence-sensitive candidate.  The
first factorization is built from the literal selected product and the exact
array returned by generated reverse erasure; the second is built from the
genuine integer divisor and quotient.  All normalization scalars are explicit
representatives of the unit produced by the physical Hensel entry. -/
theorem henselCandidate_scaled_eq_divisor_mod_primePower
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ exponent : Nat, 0 < exponent ∧
      output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
      Polynomial.map
        (Int.castRingHom
          (ZMod (selection.prime.toNat ^ exponent)))
        (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
          ((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map SparsePolyZZ.toPoly).prod) =
      Polynomial.map
        (Int.castRingHom
          (ZMod (selection.prime.toNat ^ exponent)))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  rcases selectionHenselFactors_prime_product_eq_unit_mul_source
      hcount hfactors hleadingSemantic hselection hentry with
    ⟨exponent, hexponent, scale, scaleAtPrime, houtput, hscaleUnit,
      hscaleAtPrimeUnit, hscaleAtPrime, hfullLarge, hfullPrime⟩
  rcases StrictRecombine.removeCombination_succeeds candidate output.1 hlegal with
    ⟨complement, hremove⟩
  let prime := selection.prime.toNat
  let modulus := prime ^ exponent
  let source := SparsePolyZZ.toPoly f
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  let complementInteger :=
    (complement.toList.map SparsePolyZZ.toPoly).prod
  let liftedScaled := Polynomial.C source.leadingCoeff * selectedInteger
  let divisorScaled := Polynomial.C quotient.leadingCoeff * divisor
  let scaleInt : Int := Int.ofNat scale.val
  let complementScaled :=
    Polynomial.C (scaleInt * divisor.leadingCoeff) * quotient
  let commonSource :=
    Polynomial.C scaleInt * Polynomial.C source.leadingCoeff * source
  have hprime : Nat.Prime prime := Fact.out
  have hmodulusNe : modulus ≠ 0 :=
    pow_ne_zero exponent hprime.ne_zero
  letI : NeZero modulus := ⟨hmodulusNe⟩
  have hdivides : prime ∣ modulus :=
    dvd_pow_self prime (Nat.ne_of_gt hexponent)
  have hscaleLarge : (scaleInt : ZMod modulus) = scale := by
    simp [scaleInt, modulus, prime]
  have hsourceLeading : source.leadingCoeff =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    dsimp [source]
    rw [hfactor, Polynomial.leadingCoeff_mul]
  have hpartitionLarge :=
    StrictRecombine.removeCombination_toPolyMod_product_partition
      modulus candidate output.1 complement hlegal hremove
  have hselectedMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) selectedInteger =
      ((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hsourceMapLarge : Polynomial.map
      (Int.castRingHom (ZMod modulus)) source =
        StrictHensel.toPolyMod modulus f := rfl
  have hproductLifted : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementInteger * liftedScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [Polynomial.map_mul, liftedScaled, commonSource,
      Polynomial.map_C]
    rw [hselectedMapLarge, hcomplementMapLarge]
    change _ * (Polynomial.C (source.leadingCoeff : ZMod modulus) * _) = _
    calc
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map (StrictHensel.toPolyMod modulus)).prod *
            (complement.toList.map
              (StrictHensel.toPolyMod modulus)).prod) := by ring
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (output.1.toList.map
            (StrictHensel.toPolyMod modulus)).prod := by
              rw [hpartitionLarge]
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (Polynomial.C scale * StrictHensel.toPolyMod modulus f) := by
            simpa [modulus, prime] using congrArg
              (fun value => Polynomial.C (source.leadingCoeff : ZMod modulus) * value)
              hfullLarge
      _ = Polynomial.C (scaleInt : ZMod modulus) *
          Polynomial.C (source.leadingCoeff : ZMod modulus) *
          Polynomial.map (Int.castRingHom (ZMod modulus)) source := by
            rw [hsourceMapLarge, hscaleLarge]
            ring
  have hfactorMapLarge := congrArg
    (Polynomial.map (Int.castRingHom (ZMod modulus))) hfactor
  have hproductDivisor : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementScaled * divisorScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [complementScaled, divisorScaled, commonSource,
      Polynomial.map_mul, Polynomial.map_C]
    rw [Polynomial.map_mul] at hfactorMapLarge
    change (Polynomial.C ((scaleInt * divisor.leadingCoeff : Int) : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) quotient) *
      (Polynomial.C (quotient.leadingCoeff : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) divisor) = _
    rw [Int.cast_mul, Polynomial.C_mul]
    rw [hsourceLeading]
    simp only [map_mul]
    change _ = _ * Polynomial.map (Int.castRingHom (ZMod modulus)) source
    have hfactorMapLarge' :
        Polynomial.map (Int.castRingHom (ZMod modulus)) source =
          Polynomial.map (Int.castRingHom (ZMod modulus)) divisor *
            Polynomial.map (Int.castRingHom (ZMod modulus)) quotient := by
      simpa [source] using hfactorMapLarge
    rw [hfactorMapLarge']
    have hcast (value : Int) :
        (Int.castRingHom (ZMod modulus)) value =
          (value : ZMod modulus) := rfl
    rw [hcast scaleInt, hcast divisor.leadingCoeff,
      hcast quotient.leadingCoeff]
    ring
  have hbaseB := henselCandidate_scaled_eq_divisor_mod_prime hentry
    divisor quotient hfactor hdivisorModNonzero hdivisorLeading candidate
    hlegal hassociated
  have hcop := henselCandidate_physicalComplement_coprime prime source
    output.1 complement candidate scaleAtPrime hlegal hremove
    (by simpa [prime, source] using hfullPrime) hscaleAtPrimeUnit
    (by simpa [prime, source] using hselection.goodPrime.lc_nonzero)
    (by simpa [prime, source] using hselection.goodPrime.sqfree)
  have hleading := henselCandidate_scaled_leadingCoeff hselection hentry
    divisor quotient hfactor candidate hlegal
  dsimp only at hleading
  have hproductLiftedPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementInteger * liftedScaled) commonSource
    hproductLifted
  have hproductDivisorPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementScaled * divisorScaled) commonSource
    hproductDivisor
  have hselectedMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) selectedInteger =
      ((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hbaseB' : Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
    simpa [prime, source, liftedScaled, divisorScaled, selectedInteger] using hbaseB
  have hbaseBNonzero :
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled ≠ 0 := by
    rw [← hbaseB']
    simp only [liftedScaled, Polynomial.map_mul, Polynomial.map_C]
    rw [hselectedMapPrime]
    exact mul_ne_zero
      (Polynomial.C_ne_zero.mpr (by
        simpa [prime, source] using hselection.goodPrime.lc_nonzero))
      (henselSelectedProduct_monic hentry prime candidate (by
        intro position hposition
        simpa [getElem!_pos candidate position hposition] using
          hlegal.2.2 position hposition)).ne_zero
  have hbaseA :
      Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger =
        Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled := by
    have heq :
        Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger *
            Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
          Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled *
            Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
      calc
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementInteger * liftedScaled) := by
                exact (Polynomial.map_mul (Int.castRingHom (ZMod prime))).symm
        _ = Polynomial.map (Int.castRingHom (ZMod prime)) commonSource :=
              hproductLiftedPrime
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementScaled * divisorScaled) := hproductDivisorPrime.symm
        _ = _ := by rw [Polynomial.map_mul]
    rw [hbaseB'] at heq
    exact mul_right_cancel₀ hbaseBNonzero heq
  have hcop' : IsCoprime
      (Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger)
      (Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled) := by
    simpa [prime, source, liftedScaled, hselectedMapPrime,
      hcomplementMapPrime] using hcop
  have hunique := hensel_unique prime hprime exponent hexponent commonSource
    complementInteger liftedScaled complementScaled divisorScaled
    hproductLifted hproductDivisor hbaseA hbaseB' hcop'
    hleading.1 hleading.2
  exact ⟨exponent, hexponent, houtput, hunique.2⟩

/-- The first scalar-pruning branch of the literal `zassenhausAttempt` cannot
reject a bounded candidate drawn from the actual normalized Hensel output.
Every selected physical head is one, so the generated accumulator is exactly
the source leading coefficient; the generated Mignotte precision then makes
`ZZ.symmetricMod` recover that coefficient literally. -/
theorem zassenhausLeadingPrune_accepts_live_candidate
    {active : Array SparsePolyZZ} {outputM : ZZ}
    (leading : UMonomial × ZZ)
    (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < active.size)
    (hactiveCanonical : ∀ index (hindex : index < active.size),
      StrictPolynomialMod.SparsePolyZZCanonical active[index])
    (hactiveNonempty : ∀ index (hindex : index < active.size),
      0 < active[index].size)
    (hactiveMonic : ∀ index (hindex : index < active.size),
      (SparsePolyZZ.toPoly active[index]).Monic)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (houtput : outputM = (modulus : Int))
    (hleadingBound : leading.2.natAbs * 2 < modulus) :
    ∃ leadingProduct,
      Generated.StrictRecombine.selectedLeadingProductLoop indices active
        0 leading.2 = .ok leadingProduct ∧
      ¬(ZZ.symmetricMod leadingProduct outputM ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 * leading.2)
          (ZZ.symmetricMod leadingProduct outputM) ≠ 0) := by
  have hselectedCanonical : ∀ position
      (hposition : position < indices.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        active[indices[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos active indices[position] hactive] using
      hactiveCanonical indices[position] hactive
  have hselectedNonempty : ∀ position
      (hposition : position < indices.size),
      active[indices[position]]!.isEmpty = false := by
    intro position hposition
    have hactive := hbound position hposition
    have hsize := hactiveNonempty indices[position] hactive
    simpa [getElem!_pos active indices[position] hactive,
      Array.isEmpty, Nat.ne_of_gt hsize]
  have hselectedMonic :
      ((StrictRecombine.selectSourceIndices active.toList indices.toList).map
        SparsePolyZZ.toPoly).prod.Monic := by
    apply monic_list_product
    intro candidate hcandidate
    rcases List.mem_map.mp hcandidate with ⟨lifted, hlifted, rfl⟩
    unfold StrictRecombine.selectSourceIndices at hlifted
    rcases List.mem_map.mp hlifted with ⟨index, hindex, rfl⟩
    rcases List.mem_iff_getElem.mp hindex with
      ⟨position, hposition, hindexEq⟩
    have hpositionArray : position < indices.size := by simpa using hposition
    have hactive := hbound position hpositionArray
    have hselected : indices[position] = index := by
      rw [← Array.getElem_toList hpositionArray]
      exact hindexEq
    rw [← hselected, getElem!_pos active.toList indices[position]
      (by simpa using hactive), Array.getElem_toList hactive]
    exact hactiveMonic indices[position] hactive
  have hleadingValues :=
    StrictRecombine.selectedLeadingValues_prod_eq_leadingCoeff_of_canonical
      indices active hbound hselectedCanonical hselectedNonempty
  have hvaluesOne :
      (StrictRecombine.selectedLeadingValues indices active 0).prod = 1 := by
    rw [hleadingValues, hselectedMonic.leadingCoeff]
  have hloop := StrictRecombine.selectedLeadingProductLoop_succeeds indices
    active 0 leading.2 hbound hselectedNonempty
  rw [hvaluesOne, mul_one] at hloop
  refine ⟨leading.2, hloop, ?_⟩
  have hrecovered : ZZ.symmetricMod leading.2 outputM = leading.2 := by
    rw [houtput]
    exact StrictRecombine.symmetricMod_eq_of_strict_bound leading.2 modulus
      hmodulus hleadingBound
  rw [hrecovered]
  exact StrictRecombine.zassenhaus_prune_condition_false_of_dvd
    (leading.2 * leading.2) leading.2 (dvd_mul_right _ _)

/-- Initial-Hensel specialization of the live leading-prune theorem. -/
theorem zassenhausLeadingPrune_accepts_hensel_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p 0
      output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ∃ leadingProduct,
      Generated.StrictRecombine.selectedLeadingProductLoop indices output.1
        0 leading.2 = .ok leadingProduct ∧
      ¬(ZZ.symmetricMod leadingProduct output.2 ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 * leading.2)
          (ZZ.symmetricMod leadingProduct output.2) ≠ 0) := by
  have hselectedMonic := henselSelectedIntegerProduct_monic hentry indices
    hbound
  have hselectedCanonical : ∀ position
      (hposition : position < indices.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        output.1[indices[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos output.1 indices[position] hactive] using
      hentry.outputCanonical indices[position] hactive
  have hselectedNonempty : ∀ position
      (hposition : position < indices.size),
      output.1[indices[position]]!.isEmpty = false := by
    intro position hposition
    have hactive := hbound position hposition
    rcases hentry.outputOneHead indices[position] hactive with
      ⟨head, tail, hlist⟩
    have hsize : 0 < output.1[indices[position]].size := by
      have hlength := congrArg List.length hlist
      simp at hlength
      omega
    simpa [getElem!_pos output.1 indices[position] hactive,
      Array.isEmpty, Nat.ne_of_gt hsize]
  have hleadingValues :=
    StrictRecombine.selectedLeadingValues_prod_eq_leadingCoeff_of_canonical
      indices output.1 hbound hselectedCanonical hselectedNonempty
  have hvaluesOne :
      (StrictRecombine.selectedLeadingValues indices output.1 0).prod = 1 := by
    rw [hleadingValues, hselectedMonic.leadingCoeff]
  have hloop := StrictRecombine.selectedLeadingProductLoop_succeeds indices
    output.1 0 leading.2 hbound hselectedNonempty
  rw [hvaluesOne, mul_one] at hloop
  refine ⟨leading.2, hloop, ?_⟩
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
      hcanonical hnonempty
    rw [hzero] at hhead
    exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
      (by simpa using hhead.symm)
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading
      (1 : Polynomial Int) hsourceNe (one_dvd _)
  let modulus := output.2.toNat
  have hmodulusCast : (modulus : Int) = output.2 := by
    exact Int.toNat_of_nonneg hprecision.1.le
  have hmodulus : 0 < modulus := by
    exact Int.pos_iff_toNat_pos.mp hprecision.1
  have hleadingBound : leading.2.natAbs * 2 < modulus := by
    have hzeroBound := hprecision.2 0
    simpa [modulus, hmodulusCast] using hzeroBound
  have hrecovered : ZZ.symmetricMod leading.2 output.2 = leading.2 := by
    rw [← hmodulusCast]
    exact StrictRecombine.symmetricMod_eq_of_strict_bound leading.2 modulus
      hmodulus hleadingBound
  rw [hrecovered]
  exact StrictRecombine.zassenhaus_prune_condition_false_of_dvd
    (leading.2 * leading.2) leading.2 (dvd_mul_right _ _)

/-- The literal constant-coefficient pruning branch accepts the genuine
divisor candidate.  Hensel uniqueness identifies the generated accumulator
with the quotient-leading-scaled divisor at the full returned modulus; the
generated Mignotte precision then recovers its constant coefficient exactly. -/
theorem zassenhausConstantPrune_accepts_live_divisor_candidate
    {f : SparsePolyZZ} {active : Array SparsePolyZZ} {outputM : ZZ}
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hactiveCanonical : ∀ index (hindex : index < active.size),
      StrictPolynomialMod.SparsePolyZZCanonical active[index])
    (modulus : Nat) (hmodulus : 0 < modulus)
    (houtput : outputM = (modulus : Int))
    (hselectedCongruent :
      Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
            ((StrictRecombine.selectSourceIndices active.toList
              candidate.toList).map SparsePolyZZ.toPoly).prod) =
        Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C quotient.leadingCoeff * divisor))
    (htargetBound :
      ((Polynomial.C quotient.leadingCoeff * divisor).coeff 0).natAbs * 2 <
        modulus) :
    ∃ constantProduct,
      Generated.StrictRecombine.selectedConstantProductLoop candidate active
        0 leading.2 = .ok constantProduct ∧
      ZZ.symmetricMod constantProduct outputM =
        (Polynomial.C quotient.leadingCoeff * divisor).coeff 0 ∧
      ¬(ZZ.symmetricMod constantProduct outputM ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 *
          Generated.StrictRecombine.constantTerm f)
          (ZZ.symmetricMod constantProduct outputM) ≠ 0) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hselectedCanonical : ∀ position
      (hposition : position < candidate.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        active[candidate[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos active candidate[position] hactive] using
      hactiveCanonical candidate[position] hactive
  have hconstantValues :=
    StrictRecombine.selectedConstantValues_prod_eq_coeff_zero_of_canonical
      candidate active hbound hselectedCanonical
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  let target := Polynomial.C quotient.leadingCoeff * divisor
  have hloop := StrictRecombine.selectedConstantProductLoop_succeeds candidate
    active 0 leading.2 hbound
  rw [hconstantValues] at hloop
  change Generated.StrictRecombine.selectedConstantProductLoop candidate
      active 0 leading.2 =
    .ok (leading.2 * selectedInteger.coeff 0) at hloop
  have hfront : f[0] = leading := by
    rw [Array.getElem?_eq_getElem hnonempty] at hleading
    exact Option.some.inj hleading
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcoefficientCongruence :
      ((leading.2 * selectedInteger.coeff 0 : Int) : ZMod modulus) =
        (target.coeff 0 : ZMod modulus) := by
    have hcoeff := congrArg
      (fun polynomial : Polynomial (ZMod modulus) => polynomial.coeff 0)
      hselectedCongruent
    simpa [selectedInteger, target, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.coeff_map, hsourceLeading] using hcoeff
  have hrecovered := StrictRecombine.symmetricMod_eq_of_congruent_strict_bound
    (leading.2 * selectedInteger.coeff 0) (target.coeff 0) modulus hmodulus
    hcoefficientCongruence htargetBound
  rw [← houtput] at hrecovered
  refine ⟨leading.2 * selectedInteger.coeff 0, hloop, hrecovered, ?_⟩
  rw [hrecovered]
  apply StrictRecombine.zassenhaus_prune_condition_false_of_dvd
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have hsourceConstant : (SparsePolyZZ.toPoly f).coeff 0 =
      divisor.coeff 0 * quotient.coeff 0 := by
    rw [hfactor]
    simp
  have hfConstant : Generated.StrictRecombine.constantTerm f =
      (SparsePolyZZ.toPoly f).coeff 0 :=
    StrictRecombine.sparsePolyZZ_constantTerm_eq_coeff_zero f hcanonical
  rw [hfConstant, hsourceConstant, hsourceLeadingFactor]
  refine ⟨divisor.leadingCoeff * quotient.coeff 0, ?_⟩
  simp [target]
  ring

/-- Initial-Hensel specialization of the live constant-prune theorem. -/
theorem zassenhausConstantPrune_accepts_hensel_divisor_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ constantProduct,
      Generated.StrictRecombine.selectedConstantProductLoop candidate output.1
        0 leading.2 = .ok constantProduct ∧
      ZZ.symmetricMod constantProduct output.2 =
        (Polynomial.C quotient.leadingCoeff * divisor).coeff 0 ∧
      ¬(ZZ.symmetricMod constantProduct output.2 ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 *
          Generated.StrictRecombine.constantTerm f)
          (ZZ.symmetricMod constantProduct output.2) ≠ 0) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hselectedCanonical : ∀ position
      (hposition : position < candidate.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        output.1[candidate[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos output.1 candidate[position] hactive] using
      hentry.outputCanonical candidate[position] hactive
  have hconstantValues :=
    StrictRecombine.selectedConstantValues_prod_eq_coeff_zero_of_canonical
      candidate output.1 hbound hselectedCanonical
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  have hloop := StrictRecombine.selectedConstantProductLoop_succeeds candidate
    output.1 0 leading.2 hbound
  rw [hconstantValues] at hloop
  change Generated.StrictRecombine.selectedConstantProductLoop candidate
      output.1 0 leading.2 =
    .ok (leading.2 * selectedInteger.coeff 0) at hloop
  rcases henselCandidate_scaled_eq_divisor_mod_primePower hcount hfactors
      hleadingSemantic hselection hentry divisor quotient hfactor
      hdivisorModNonzero hdivisorLeading candidate hlegal hassociated with
    ⟨exponent, hexponent, houtput, hunique⟩
  let modulus := selection.prime.toNat ^ exponent
  let target := Polynomial.C quotient.leadingCoeff * divisor
  have hmodulus : 0 < modulus := by
    exact pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos exponent
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcoefficientCongruence :
      ((leading.2 * selectedInteger.coeff 0 : Int) : ZMod modulus) =
        (target.coeff 0 : ZMod modulus) := by
    have hcoeff := congrArg
      (fun polynomial : Polynomial (ZMod modulus) => polynomial.coeff 0)
      hunique
    simpa [selectedInteger, target, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.coeff_map, hsourceLeading] using hcoeff
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    apply hselection.goodPrime.lc_nonzero
    rw [hzero]
    simp
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading divisor
      hsourceNe (hfactor ▸ dvd_mul_right divisor quotient)
  have hdivisorNe : divisor ≠ 0 := by
    intro hzero
    apply hdivisorModNonzero
    simp [hzero]
  have hdivisorLeadingNe : divisor.leadingCoeff ≠ 0 := by
    exact Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have htargetBound : (target.coeff 0).natAbs * 2 < modulus := by
    have hlargeInt := hprecision.2 0
    rw [houtput] at hlargeInt
    have hlarge : (leading.2 * divisor.coeff 0).natAbs * 2 < modulus := by
      exact_mod_cast hlargeInt
    have hdivisorAbs : 0 < divisor.leadingCoeff.natAbs :=
      Int.natAbs_pos.mpr hdivisorLeadingNe
    have hle : (quotient.leadingCoeff * divisor.coeff 0).natAbs * 2 ≤
        (leading.2 * divisor.coeff 0).natAbs * 2 := by
      rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul,
        Int.natAbs_mul]
      have hbase := Nat.le_mul_of_pos_left
        (quotient.leadingCoeff.natAbs * (divisor.coeff 0).natAbs)
        hdivisorAbs
      have hscaled := Nat.mul_le_mul_right 2 hbase
      simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
    have htargetCoeff : target.coeff 0 =
        quotient.leadingCoeff * divisor.coeff 0 := by
      simp [target]
    rw [htargetCoeff]
    exact lt_of_le_of_lt hle hlarge
  have hrecovered := StrictRecombine.symmetricMod_eq_of_congruent_strict_bound
    (leading.2 * selectedInteger.coeff 0) (target.coeff 0) modulus hmodulus
    hcoefficientCongruence htargetBound
  have hmodulusInt : (modulus : Int) = output.2 := by
    simpa [modulus] using houtput.symm
  rw [hmodulusInt] at hrecovered
  refine ⟨leading.2 * selectedInteger.coeff 0, hloop, hrecovered, ?_⟩
  rw [hrecovered]
  apply StrictRecombine.zassenhaus_prune_condition_false_of_dvd
  have hsourceConstant : (SparsePolyZZ.toPoly f).coeff 0 =
      divisor.coeff 0 * quotient.coeff 0 := by
    rw [hfactor]
    simp
  have hfConstant : Generated.StrictRecombine.constantTerm f =
      (SparsePolyZZ.toPoly f).coeff 0 :=
    StrictRecombine.sparsePolyZZ_constantTerm_eq_coeff_zero f hcanonical
  rw [hfConstant, hsourceConstant, hsourceLeadingFactor]
  refine ⟨divisor.leadingCoeff * quotient.coeff 0, ?_⟩
  simp [target]
  ring

/-- After both literal scalar prunes accept a genuine Hensel candidate, the
generated checked-index conversion, modular trial product, symmetric recovery,
and primitive-part calls all execute.  The physical symmetric array denotes
exactly the quotient-leading-scaled genuine divisor. -/
theorem zassenhausCandidate_executes_through_primitive_of_live
    {f : SparsePolyZZ} {active : Array SparsePolyZZ} {outputM : ZZ}
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hactiveFits : active.size ≤ 2 ^ 31)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (houtput : outputM = (modulus : Int))
    (hselectedCongruent :
      Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
            ((StrictRecombine.selectSourceIndices active.toList
              candidate.toList).map SparsePolyZZ.toPoly).prod) =
        Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C quotient.leadingCoeff * divisor))
    (htargetBound : ∀ degree,
      ((Polynomial.C quotient.leadingCoeff * divisor).coeff degree).natAbs * 2 <
        modulus) :
    ∃ candidate32 product symmetric content recoveredFactor,
      Generated.StrictRecombine.combinationToInt32 candidate = .ok candidate32 ∧
      Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate32 active
        outputM 0 #[(⟨0⟩, leading.2)] = .ok product ∧
      Generated.StrictRecombine.symmetricModRaw product outputM = .ok symmetric ∧
      SparsePolyZZ.toPoly symmetric =
        Polynomial.C quotient.leadingCoeff * divisor ∧
      StrictPolynomialMod.SparsePolyZZCanonical symmetric ∧
      Generated.StrictRecombine.primitiveRaw symmetric =
        .ok (content, recoveredFactor) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hfits : ∀ position (hposition : position < candidate.size),
      candidate[position] < 2 ^ 31 := by
    intro position hposition
    exact lt_of_lt_of_le (hbound position hposition) hactiveFits
  rcases StrictRecombine.combinationToInt32_toList candidate hfits with
    ⟨candidate32, hconvert, _⟩
  have hvalid := StrictRecombine.combinationToInt32_candidate_valid candidate
    active.size candidate32 hbound hactiveFits hconvert
  have houtputPositive : 0 < outputM := by
    rw [houtput]
    exact_mod_cast hmodulus
  rcases StrictRecombine.trialProductLoop_complete ⟨()⟩ candidate32 active
      outputM 0 #[(⟨0⟩, leading.2)] houtputPositive.ne' hvalid with
    ⟨product, hproduct⟩
  have htrial := StrictRecombine.trialProductLoop_source_indices_refines
    modulus hmodulus candidate active candidate32
    #[(⟨0⟩, leading.2)] product hbound hactiveFits hconvert
    (by simpa [houtput] using hproduct)
  have hfront : f[0] = leading := by
    rw [Array.getElem?_eq_getElem hnonempty] at hleading
    exact Option.some.inj hleading
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcongruent : Polynomial.map (Int.castRingHom (ZMod modulus))
      (SparsePolyZZ.toPoly product) =
      Polynomial.map (Int.castRingHom (ZMod modulus))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
    change StrictHensel.toPolyMod modulus product = _
    rw [htrial]
    have hinitialMod : StrictHensel.toPolyMod modulus
        #[(⟨0⟩, leading.2)] = Polynomial.C (leading.2 : ZMod modulus) := by
      simp [StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
    rw [hinitialMod]
    simpa [StrictHensel.toPolyMod, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.map_list_prod, hsourceLeading] using hselectedCongruent
  have hinitialCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical #[(⟨0⟩, leading.2)] := by
    constructor
    · simp
    · intro term hterm
      simp at hterm
      subst term
      rw [← hfront]
      exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
  have hproductCanonical := StrictRecombine.trialProductLoop_canonical ⟨()⟩
    candidate32 active outputM 0 #[(⟨0⟩, leading.2)] product
    hinitialCanonical hproduct
  rcases StrictRecombine.symmetricModRaw_complete product outputM
      houtputPositive with ⟨symmetric, hsymmetric⟩
  have hsymmetricPoly :=
    StrictRecombine.symmetricModRaw_recovers_strictly_bounded_target product
      symmetric (Polynomial.C quotient.leadingCoeff * divisor) modulus hmodulus
      hproductCanonical (by simpa [houtput] using hsymmetric) hcongruent
      htargetBound
  have hsymmetricCanonical := StrictRecombine.symmetricModRaw_canonical product
    symmetric modulus hmodulus hproductCanonical
    (by simpa [houtput] using hsymmetric)
  rcases StrictRecombine.primitiveRaw_complete symmetric hsymmetricCanonical with
    ⟨content, recoveredFactor, hprimitive⟩
  exact ⟨candidate32, product, symmetric, content, recoveredFactor, hconvert,
    hproduct, hsymmetric, hsymmetricPoly, hsymmetricCanonical, hprimitive⟩

/-- Initial-Hensel specialization of the live-state execution lemma. -/
theorem zassenhausCandidate_executes_through_primitive
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hactiveFits : output.1.size ≤ 2 ^ 31)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ candidate32 product symmetric content recoveredFactor,
      Generated.StrictRecombine.combinationToInt32 candidate = .ok candidate32 ∧
      Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate32 output.1
        output.2 0 #[(⟨0⟩, leading.2)] = .ok product ∧
      Generated.StrictRecombine.symmetricModRaw product output.2 = .ok symmetric ∧
      SparsePolyZZ.toPoly symmetric =
        Polynomial.C quotient.leadingCoeff * divisor ∧
      StrictPolynomialMod.SparsePolyZZCanonical symmetric ∧
      Generated.StrictRecombine.primitiveRaw symmetric =
        .ok (content, recoveredFactor) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hfits : ∀ position (hposition : position < candidate.size),
      candidate[position] < 2 ^ 31 := by
    intro position hposition
    exact lt_of_lt_of_le (hbound position hposition) hactiveFits
  rcases StrictRecombine.combinationToInt32_toList candidate hfits with
    ⟨candidate32, hconvert, _⟩
  have hvalid := StrictRecombine.combinationToInt32_candidate_valid candidate
    output.1.size candidate32 hbound hactiveFits hconvert
  have houtputPositive :=
    (StrictRecombine.hensel_output_modulus_bounds_scaled_divisor hentry rfl
      hcanonical hnonempty hdegree leading hleading divisor
      (by
        intro hzero
        apply hselection.goodPrime.lc_nonzero
        rw [hzero]
        simp)
      (hfactor ▸ dvd_mul_right divisor quotient)).1
  rcases StrictRecombine.trialProductLoop_complete ⟨()⟩ candidate32 output.1
      output.2 0 #[(⟨0⟩, leading.2)] houtputPositive.ne' hvalid with
    ⟨product, hproduct⟩
  rcases henselCandidate_scaled_eq_divisor_mod_primePower hcount hfactors
      hleadingSemantic hselection hentry divisor quotient hfactor
      hdivisorModNonzero hdivisorLeading candidate hlegal hassociated with
    ⟨exponent, hexponent, houtput, hunique⟩
  let modulus := selection.prime.toNat ^ exponent
  let target := Polynomial.C quotient.leadingCoeff * divisor
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos exponent
  have houtputCast : (modulus : Int) = output.2 := by
    simpa [modulus] using houtput.symm
  have htrial := StrictRecombine.trialProductLoop_source_indices_refines
    modulus hmodulus candidate output.1 candidate32
    #[(⟨0⟩, leading.2)] product hbound hactiveFits hconvert
    (by simpa [houtputCast] using hproduct)
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcongruent : Polynomial.map (Int.castRingHom (ZMod modulus))
      (SparsePolyZZ.toPoly product) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) target := by
    change StrictHensel.toPolyMod modulus product = _
    rw [htrial]
    have hinitialMod : StrictHensel.toPolyMod modulus
        #[(⟨0⟩, leading.2)] = Polynomial.C (leading.2 : ZMod modulus) := by
      simp [StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
    rw [hinitialMod]
    simpa [modulus, StrictHensel.toPolyMod, target, Polynomial.map_mul,
      Polynomial.map_C, Polynomial.map_list_prod, hsourceLeading] using hunique
  have hinitialCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical #[(⟨0⟩, leading.2)] := by
    constructor
    · simp
    · intro term hterm
      simp at hterm
      subst term
      have hfront : f[0] = leading := by
        rw [Array.getElem?_eq_getElem hnonempty] at hleading
        exact Option.some.inj hleading
      rw [← hfront]
      exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
  have hproductCanonical := StrictRecombine.trialProductLoop_canonical ⟨()⟩
    candidate32 output.1 output.2 0 #[(⟨0⟩, leading.2)] product
    hinitialCanonical hproduct
  rcases StrictRecombine.symmetricModRaw_complete product output.2
      houtputPositive with ⟨symmetric, hsymmetric⟩
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    apply hselection.goodPrime.lc_nonzero
    rw [hzero]
    simp
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading divisor hsourceNe
      (hfactor ▸ dvd_mul_right divisor quotient)
  have hdivisorNe : divisor ≠ 0 := by
    intro hzero
    apply hdivisorModNonzero
    simp [hzero]
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have htargetBound : ∀ degree, (target.coeff degree).natAbs * 2 < modulus := by
    intro degree
    have hlargeInt := hprecision.2 degree
    rw [← houtputCast] at hlargeInt
    have hlarge : (leading.2 * divisor.coeff degree).natAbs * 2 < modulus := by
      exact_mod_cast hlargeInt
    have hdivisorLeadingAbs : 0 < divisor.leadingCoeff.natAbs :=
      Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe)
    have hle : (quotient.leadingCoeff * divisor.coeff degree).natAbs * 2 ≤
        (leading.2 * divisor.coeff degree).natAbs * 2 := by
      rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul,
        Int.natAbs_mul]
      have hbase := Nat.le_mul_of_pos_left
        (quotient.leadingCoeff.natAbs * (divisor.coeff degree).natAbs)
        hdivisorLeadingAbs
      have hscaled := Nat.mul_le_mul_right 2 hbase
      simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
    have htargetCoeff : target.coeff degree =
        quotient.leadingCoeff * divisor.coeff degree := by simp [target]
    rw [htargetCoeff]
    exact lt_of_le_of_lt hle hlarge
  have hsymmetricPoly :=
    StrictRecombine.symmetricModRaw_recovers_strictly_bounded_target product
      symmetric target modulus hmodulus hproductCanonical
      (by simpa [houtputCast] using hsymmetric) hcongruent htargetBound
  have hsymmetricCanonical := StrictRecombine.symmetricModRaw_canonical product
    symmetric modulus hmodulus hproductCanonical
    (by simpa [houtputCast] using hsymmetric)
  rcases StrictRecombine.primitiveRaw_complete symmetric hsymmetricCanonical with
    ⟨content, recoveredFactor, hprimitive⟩
  exact ⟨candidate32, product, symmetric, content, recoveredFactor, hconvert,
    hproduct, hsymmetric, hsymmetricPoly, hsymmetricCanonical, hprimitive⟩

/-- A primitive factor physically returned from a nonzero scalar multiple of
a primitive integer divisor divides that divisor over `ℤ[x]`.  The proof maps
the concrete `primitiveRaw` equation to `ℚ[x]`, cancels only proved nonzero
constant units there, and applies Gauss's lemma in the reverse direction. -/
theorem primitiveRaw_factor_dvd_scaled_primitive_divisor
    (symmetric recoveredFactor : SparsePolyZZ) (content : Int)
    (divisor quotient : Polynomial Int)
    (hsymmetric : SparsePolyZZ.toPoly symmetric =
      Polynomial.C quotient.leadingCoeff * divisor)
    (hsymmetricCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical symmetric)
    (hprimitive : Generated.StrictRecombine.primitiveRaw symmetric =
      .ok (content, recoveredFactor))
    (hdivisorPrimitive : divisor.IsPrimitive)
    (hdivisorNe : divisor ≠ 0) (hquotientNe : quotient ≠ 0) :
    SparsePolyZZ.toPoly recoveredFactor ∣ divisor := by
  have hquotientLeadingNe : quotient.leadingCoeff ≠ 0 :=
    Polynomial.leadingCoeff_ne_zero.mpr hquotientNe
  have hscaledNe : Polynomial.C quotient.leadingCoeff * divisor ≠ 0 :=
    mul_ne_zero (Polynomial.C_ne_zero.mpr hquotientLeadingNe) hdivisorNe
  have hsymmetricNe : SparsePolyZZ.toPoly symmetric ≠ 0 := by
    rw [hsymmetric]
    exact hscaledNe
  have hsymmetricNonempty : 0 < symmetric.size := by
    by_contra hnot
    have hzero : symmetric.size = 0 := Nat.eq_zero_of_not_pos hnot
    have hempty : symmetric = #[] := Array.size_eq_zero_iff.mp hzero
    apply hsymmetricNe
    simp [hempty, SparsePolyZZ.toPoly]
  have hrecoveredPrimitive := StrictRecombine.primitiveRaw_isPrimitive
    symmetric recoveredFactor content hsymmetricNonempty hsymmetricCanonical
      hprimitive
  have hsemantic := StrictRecombine.primitiveRaw_toPoly symmetric
    recoveredFactor content hprimitive
  have hcontentNe : content ≠ 0 := by
    intro hzero
    apply hsymmetricNe
    rw [hsemantic, hzero]
    simp
  have hmap := congrArg (Polynomial.map (Int.castRingHom ℚ)) hsemantic
  rw [hsymmetric] at hmap
  simp only [Polynomial.map_mul, Polynomial.map_C] at hmap
  have hcontentUnit : IsUnit (Polynomial.C (content : ℚ)) :=
    Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr
      (Int.cast_ne_zero.mpr hcontentNe))
  have hquotientUnit : IsUnit
      (Polynomial.C (quotient.leadingCoeff : ℚ)) :=
    Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr
      (Int.cast_ne_zero.mpr hquotientLeadingNe))
  let recoveredQ :=
    (SparsePolyZZ.toPoly recoveredFactor).map (Int.castRingHom ℚ)
  let divisorQ := divisor.map (Int.castRingHom ℚ)
  have hleft : Associated
      (Polynomial.C (content : ℚ) * recoveredQ) recoveredQ :=
    (associated_isUnit_mul_left_iff hcontentUnit).mpr (Associated.refl _)
  have hright : Associated
      (Polynomial.C (quotient.leadingCoeff : ℚ) * divisorQ) divisorQ :=
    (associated_isUnit_mul_left_iff hquotientUnit).mpr (Associated.refl _)
  have hassociated : Associated recoveredQ divisorQ := by
    exact hleft.symm.trans ((Associated.of_eq hmap.symm).trans hright)
  apply (Polynomial.IsPrimitive.Int.dvd_iff_map_cast_dvd_map_cast
    (SparsePolyZZ.toPoly recoveredFactor) divisor hrecoveredPrimitive
      hdivisorPrimitive).mpr
  exact hassociated.dvd

/-- A legal divisor candidate in an arbitrary live Zassenhaus state makes the
literal generated attempt return `extracted`.  The hypotheses are precisely
the representation, full-modulus congruence, and strict precision facts that
the outer loop must preserve; every computational stage below is the generated
C++-shaped operation. -/
theorem zassenhausAttempt_extracts_live_divisor_candidate
    {f : SparsePolyZZ} {active : Array SparsePolyZZ} {outputM : ZZ}
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorPrimitive : divisor.IsPrimitive)
    (hdivisorNe : divisor ≠ 0) (hquotientNe : quotient ≠ 0)
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hactiveFits : active.size ≤ 2 ^ 31)
    (hactiveCanonical : ∀ index (hindex : index < active.size),
      StrictPolynomialMod.SparsePolyZZCanonical active[index])
    (hactiveNonempty : ∀ index (hindex : index < active.size),
      0 < active[index].size)
    (hactiveMonic : ∀ index (hindex : index < active.size),
      (SparsePolyZZ.toPoly active[index]).Monic)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (houtput : outputM = (modulus : Int))
    (hleadingBound : leading.2.natAbs * 2 < modulus)
    (hselectedCongruent :
      Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
            ((StrictRecombine.selectSourceIndices active.toList
              candidate.toList).map SparsePolyZZ.toPoly).prod) =
        Polynomial.map (Int.castRingHom (ZMod modulus))
          (Polynomial.C quotient.leadingCoeff * divisor))
    (htargetBound : ∀ degree,
      ((Polynomial.C quotient.leadingCoeff * divisor).coeff degree).natAbs * 2 <
        modulus) :
    ∃ recoveredFactor recoveredQuotient,
      Generated.StrictRecombine.zassenhausAttempt f active outputM candidate =
        .ok (.extracted recoveredFactor recoveredQuotient) := by
  have hfront : f[0] = leading := by
    rw [Array.getElem?_eq_getElem hnonempty] at hleading
    exact Option.some.inj hleading
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  rcases zassenhausLeadingPrune_accepts_live_candidate leading candidate hbound
      hactiveCanonical hactiveNonempty hactiveMonic modulus hmodulus houtput
      hleadingBound with
    ⟨leadingProduct, hleadingRun, hleadingAccept⟩
  rcases zassenhausConstantPrune_accepts_live_divisor_candidate hcanonical
      hnonempty leading hleading divisor quotient hfactor candidate hlegal
      hactiveCanonical modulus hmodulus houtput hselectedCongruent
      (htargetBound 0) with
    ⟨constantProduct, hconstantRun, _hconstantRecovered, hconstantAccept⟩
  rcases zassenhausCandidate_executes_through_primitive_of_live hcanonical
      hnonempty leading hleading divisor quotient candidate hlegal hactiveFits
      modulus hmodulus houtput hselectedCongruent htargetBound with
    ⟨candidate32, product, symmetric, content, recoveredFactor, hconvert,
      hproduct, hsymmetric, hsymmetricPoly, hsymmetricCanonical, hprimitive⟩
  have hrecoveredDvdDivisor :=
    primitiveRaw_factor_dvd_scaled_primitive_divisor symmetric recoveredFactor
      content divisor quotient hsymmetricPoly hsymmetricCanonical hprimitive
      hdivisorPrimitive hdivisorNe hquotientNe
  have hrecoveredDvdSource : SparsePolyZZ.toPoly recoveredFactor ∣
      SparsePolyZZ.toPoly f := by
    exact dvd_trans hrecoveredDvdDivisor
      (hfactor ▸ dvd_mul_right divisor quotient)
  have hrecoveredCanonical := StrictRecombine.primitiveRaw_canonical symmetric
    recoveredFactor content hsymmetricCanonical hprimitive
  have hrecoveredNonempty : 0 < recoveredFactor.size := by
    by_contra hnot
    have hzero : recoveredFactor.size = 0 := Nat.eq_zero_of_not_pos hnot
    have hempty : recoveredFactor = #[] := Array.size_eq_zero_iff.mp hzero
    rw [hempty] at hrecoveredDvdDivisor
    simp [SparsePolyZZ.toPoly] at hrecoveredDvdDivisor
    exact hdivisorNe hrecoveredDvdDivisor
  rcases StrictRecombine.exactDivmodRaw_complete_of_dvd f recoveredFactor
      hcanonical hrecoveredCanonical hrecoveredNonempty hrecoveredDvdSource with
    ⟨rawQuotient, hdivmod⟩
  have hrawQuotientCanonical :=
    StrictRecombine.exactDivmodRaw_quotient_canonical f recoveredFactor
      rawQuotient #[] hcanonical.2 hdivmod
  rcases StrictRecombine.primitiveRaw_complete rawQuotient
      hrawQuotientCanonical with
    ⟨quotientContent, recoveredQuotient, hquotientPrimitive⟩
  refine ⟨recoveredFactor, recoveredQuotient, ?_⟩
  unfold Generated.StrictRecombine.zassenhausAttempt
  rw [dif_pos hnonempty]
  dsimp only
  rw [hfront]
  simp only [hleadingRun]
  rw [if_neg hleadingAccept]
  simp only [hconstantRun]
  rw [if_neg hconstantAccept]
  simp only [hconvert, hproduct, hsymmetric, hprimitive, hdivmod]
  simp only [Array.isEmpty_empty, if_true, hquotientPrimitive]

/-- Every factor returned by the literal successful attempt is irreducible
once the generated scans of all smaller positive subset sizes have physically
returned `exhausted`. -/
theorem zassenhausAttempt_extracted_irreducible
    (prime exponent : Nat) [Fact (Nat.Prime prime)]
    (source factor quotientSparse : SparsePolyZZ)
    (active : Array SparsePolyZZ) (outer : Array Nat)
    (subsetSize : Nat)
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (hhistory : StrictRecombine.SmallerZassenhausScansExhausted source active
      (((prime ^ exponent : Nat) : ZZ)) subsetSize)
    (hsubsetPositive : 0 < subsetSize)
    (houter : StrictRecombine.LegalCombination active.size outer.size outer)
    (houterSize : outer.size = subsetSize)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (leading : UMonomial × ZZ) (hleadingCell : source[0]? = some leading)
    (hleadingBound : leading.2.natAbs * 2 < prime ^ exponent)
    (hrecovery : ∀ divisor quotient,
      SparsePolyZZ.toPoly source = divisor * quotient → divisor ≠ 0 →
      ∀ degree,
        ((Polynomial.C quotient.leadingCoeff * divisor).coeff degree).natAbs *
          2 < prime ^ exponent)
    (hattempt : Generated.StrictRecombine.zassenhausAttempt source active
      (((prime ^ exponent : Nat) : ZZ)) outer =
        .ok (.extracted factor quotientSparse)) :
    Irreducible (SparsePolyZZ.toPoly factor) := by
  have hfactorProperties :=
    StrictRecombine.zassenhausAttempt_extracted_factor_canonical_primitive
      source factor quotientSparse active
      (((prime ^ exponent : Nat) : ZZ)) outer hcanonical hattempt
  have hfactorLeading :
      ((SparsePolyZZ.toPoly factor).leadingCoeff : ZMod prime) ≠ 0 := by
    rcases StrictRecombine.zassenhausAttempt_extracted_unit_scalar source factor
        quotientSparse active (((prime ^ exponent : Nat) : ZZ)) outer
        hprimitive hattempt with ⟨scalar, _hscalar, hextraction⟩
    have hlc := congrArg Polynomial.leadingCoeff hextraction
    simp only [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C] at hlc
    have hlcMod := congrArg (fun value : Int => (value : ZMod prime)) hlc
    simp only [Int.cast_mul] at hlcMod
    intro hzero
    apply hleading
    rw [hlcMod, hzero]
    simp
  have hfactorAssociated :=
    StrictRecombine.zassenhausAttempt_extracted_factor_mod_associated_selected
      source factor quotientSparse active (prime ^ exponent) prime outer
      (pow_pos (Fact.out : Nat.Prime prime).pos exponent)
      (Fact.out : Nat.Prime prime).pos
      (dvd_pow_self prime (Nat.ne_of_gt state.exponentPositive))
      (fun position hposition => by
        simpa [getElem!_pos outer position hposition] using
          houter.2.2 position hposition)
      activeState.fitsInt32
      (by
        have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source
          hcanonical hnonempty
        simpa [getElem!_pos source 0 hnonempty] using (hhead ▸ hleading))
      activeState.irreducible hattempt
  have houterPositive : 0 < outer.size := by
    simpa [houterSize] using hsubsetPositive
  have hfactorNonunit : ¬IsUnit (SparsePolyZZ.toPoly factor) := by
    intro hunit
    have hmappedUnit := hunit.map
      (Polynomial.mapRingHom (Int.castRingHom (ZMod prime)))
    have hselectedUnit := hfactorAssociated.isUnit_iff.mp hmappedUnit
    let selected :=
      (StrictRecombine.selectSourceIndices active.toList outer.toList).map
        (StrictHensel.toPolyMod prime)
    have houterZero : outer[0] < active.size := by
      simpa [getElem!_pos outer 0 houterPositive] using
        houter.2.2 0 houterPositive
    have hsourceMem : active[outer[0]] ∈
        StrictRecombine.selectSourceIndices active.toList outer.toList := by
      unfold StrictRecombine.selectSourceIndices
      apply List.mem_map.mpr
      refine ⟨outer[0], Array.getElem_mem_toList houterPositive, ?_⟩
      simpa [getElem!_pos active.toList outer[0] (by simpa using houterZero),
        Array.getElem_toList houterZero]
    have hmappedMem : StrictHensel.toPolyMod prime active[outer[0]] ∈
        selected := List.mem_map.mpr ⟨active[outer[0]], hsourceMem, rfl⟩
    have hcellUnit := List.prod_isUnit_iff.mp
      (by simpa [selected] using hselectedUnit) _ hmappedMem
    exact (activeState.irreducible outer[0] houterZero).not_isUnit hcellUnit
  refine ⟨hfactorNonunit, ?_⟩
  intro left right hfactorization
  by_cases hleftUnit : IsUnit left
  · exact Or.inl hleftUnit
  by_cases hrightUnit : IsUnit right
  · exact Or.inr hrightUnit
  rcases zassenhausAttempt_reducible_has_smaller_candidate prime exponent
      source factor quotientSparse active outer subsetSize state activeState
      houter houterSize hcanonical hnonempty hprimitive hsquarefree hleading
      hrecovery left right hfactorization hleftUnit hrightUnit hattempt with
    ⟨inner, innerQuotient, hinner, hinnerPositive, hinnerSmall,
      hleftPrimitive, hleftNe, hinnerQuotientNe, hsourceFactorization,
      hcongruent, htargetBound⟩
  rcases zassenhausAttempt_extracts_live_divisor_candidate hcanonical hnonempty
      leading hleadingCell left innerQuotient hsourceFactorization
      hleftPrimitive hleftNe hinnerQuotientNe inner hinner
      activeState.fitsInt32 activeState.canonical activeState.nonempty
      activeState.monic (prime ^ exponent)
      (pow_pos (Fact.out : Nat.Prime prime).pos exponent) rfl hleadingBound
      hcongruent htargetBound with
    ⟨recoveredFactor, recoveredQuotient, hinnerAttempt⟩
  have hrejected := hhistory.rejects inner hinnerPositive hinnerSmall hinner
  rw [hinnerAttempt] at hrejected
  simp at hrejected

/-- One literal successful fixed-size scan followed by the literal physical
removal advances every live invariant consumed by the generated outer loop. -/
theorem scanExtraction_liveStep
    (prime exponent count : Nat) [Fact (Nat.Prime prime)]
    (source factor quotient : SparsePolyZZ)
    (active : Array SparsePolyZZ) (candidate : Array Nat)
    (candidateSize : candidate.size = count)
    (state : LiveHenselProduct prime exponent source active)
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (precision : LiveRecoveryPrecision (prime ^ exponent) source)
    (history : StrictRecombine.SmallerZassenhausScansExhausted source active
      (((prime ^ exponent : Nat) : ZZ)) count)
    (hcount : 0 < count) (hfits : count ≤ active.size)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hscan : Generated.StrictRecombine.scanZassenhausCombinations
      (StrictRecombine.concreteZassenhausTermination.combinations
        active.size count)
      source active (((prime ^ exponent : Nat) : ZZ))
      (Generated.StrictRecombine.initialCombination count)
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) =
        .ok (.extracted factor quotient candidate candidateSize)) :
    ∃ remaining : Array SparsePolyZZ,
      Generated.StrictRecombine.removeCombination candidate active =
        .ok remaining ∧
      Irreducible (SparsePolyZZ.toPoly factor) ∧
      StrictPolynomialMod.SparsePolyZZCanonical quotient ∧
      (SparsePolyZZ.toPoly quotient).IsPrimitive ∧
      0 < quotient.size ∧
      ((SparsePolyZZ.toPoly quotient).leadingCoeff : ZMod prime) ≠ 0 ∧
      StrictRecombine.LiveActiveFactors prime remaining ∧
      LiveHenselProduct prime exponent quotient remaining ∧
      LiveRecoveryPrecision (prime ^ exponent) quotient ∧
      Squarefree (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly quotient)) := by
  have hlegalCount : StrictRecombine.LegalCombination active.size count
      candidate :=
    StrictRecombine.concreteScan_extracted_legal source factor quotient active
      (((prime ^ exponent : Nat) : ZZ)) candidate candidateSize hfits hscan
  have hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate := by simpa [candidateSize] using hlegalCount
  have hattempt : Generated.StrictRecombine.zassenhausAttempt source active
      (((prime ^ exponent : Nat) : ZZ)) candidate =
        .ok (.extracted factor quotient) :=
    StrictRecombine.scanZassenhausCombinations_extracted_attempt
      (StrictRecombine.concreteZassenhausTermination.combinations
        active.size count)
      source factor quotient active (((prime ^ exponent : Nat) : ZZ))
      (Generated.StrictRecombine.initialCombination count) candidate
      candidateSize
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) hscan
  rcases StrictRecombine.removeCombination_succeeds candidate active hlegal with
    ⟨remaining, hremove⟩
  have hfactorQuotient :=
    StrictRecombine.scanZassenhausCombinations_extracted_canonical_primitive
      (StrictRecombine.concreteZassenhausTermination.combinations
        active.size count)
      source factor quotient active (((prime ^ exponent : Nat) : ZZ))
      (Generated.StrictRecombine.initialCombination count) candidate
      candidateSize hcanonical hnonempty
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) hscan
  have hmodCertificate :=
    StrictRecombine.scanZassenhausCombinations_extracted_mod_certificate
      source factor quotient active (prime ^ exponent) prime
      (Generated.StrictRecombine.initialCombination count) candidate
      candidateSize
      (pow_pos (Fact.out : Nat.Prime prime).pos exponent)
      (Fact.out : Nat.Prime prime).pos
      (dvd_pow_self prime (Nat.ne_of_gt state.exponentPositive))
      activeState.fitsInt32 hcanonical hnonempty hprimitive hleading
      activeState.irreducible
      (StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size count hfits) hscan
  have hquotientLeading := hmodCertificate.2
  have hquotientNe : SparsePolyZZ.toPoly quotient ≠ 0 := by
    intro hzero
    apply hquotientLeading
    rw [hzero]
    simp
  have hquotientNonempty : 0 < quotient.size := by
    by_contra hnot
    have hempty : quotient = #[] := Array.size_eq_zero_iff.mp
      (Nat.eq_zero_of_not_pos hnot)
    apply hquotientNe
    simp [hempty, SparsePolyZZ.toPoly]
  let leading := source[0]
  have hleadingCell : source[0]? = some leading := by
    exact Array.getElem?_eq_getElem hnonempty
  have hheadLeading : leading.2 = (SparsePolyZZ.toPoly source).leadingCoeff := by
    dsimp [leading]
    exact (StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source hcanonical
      hnonempty).symm
  have hleadingBound : leading.2.natAbs * 2 < prime ^ exponent := by
    rw [hheadLeading]
    exact precision.leadingBound
  have hfactorIrreducible := zassenhausAttempt_extracted_irreducible prime
    exponent source factor quotient active candidate count state activeState
    history hcount hlegal candidateSize hcanonical hnonempty hprimitive
    hsquarefree hleading leading hleadingCell hleadingBound
    precision.scaledFactorBound hattempt
  have hfactorNe : SparsePolyZZ.toPoly factor ≠ 0 :=
    hfactorIrreducible.ne_zero
  have hactiveNext := activeState.removeCombination candidate hremove
  have hproductNext := state.extractScan activeState hcanonical hnonempty
    hprimitive hleading hfits candidateSize hscan hremove
  have hprecisionNext := precision.extract hprimitive hfactorNe hquotientNe
    hattempt
  have hsquarefreeNext := zassenhausAttempt_extracted_quotient_squarefree
    prime source factor quotient active (((prime ^ exponent : Nat) : ZZ))
    candidate hprimitive hsquarefree hattempt
  exact ⟨remaining, hremove, hfactorIrreducible, hfactorQuotient.2.1,
    hfactorQuotient.2.2, hquotientNonempty, hquotientLeading, hactiveNext,
    hproductNext, hprecisionNext, hsquarefreeNext⟩

/-- Concrete terminal state reached by the generated Zassenhaus outer loop.
It exposes the exact physical state on which `finishZassenhaus` was called. -/
structure ZassenhausTerminalCertificate (prime exponent : Nat)
    (output : Array SparsePolyZZ) : Type where
  active : Array SparsePolyZZ
  source : SparsePolyZZ
  result : Array SparsePolyZZ
  subsetSize : Nat
  output_eq : output =
    Generated.StrictRecombine.finishZassenhaus source result
  stopped : ¬(2 * subsetSize ≤ active.size)
  subsetPositive : 0 < subsetSize
  canonical : StrictPolynomialMod.SparsePolyZZCanonical source
  nonempty : 0 < source.size
  primitive : (SparsePolyZZ.toPoly source).IsPrimitive
  leading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0
  squarefree : Squarefree
    (Polynomial.map (Int.castRingHom (ZMod prime))
      (SparsePolyZZ.toPoly source))
  activeState : StrictRecombine.LiveActiveFactors prime active
  productState : LiveHenselProduct prime exponent source active
  precision : LiveRecoveryPrecision (prime ^ exponent) source
  resultIrreducible : FactorArrayIrreducible result
  history : StrictRecombine.SmallerZassenhausScansExhausted source active
    (((prime ^ exponent : Nat) : ZZ)) subsetSize

/-- Execute the complete generated outer loop to its literal terminal state
while threading every live invariant and the physical exhaustion history. -/
theorem zassenhausLoop_live_terminal
    (prime exponent : Nat) [Fact (Nat.Prime prime)]
    (active : Array SparsePolyZZ) (source : SparsePolyZZ)
    (result : Array SparsePolyZZ) (subsetSize : Nat)
    (hsubsetPositive : 0 < subsetSize)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical source)
    (hnonempty : 0 < source.size)
    (hprimitive : (SparsePolyZZ.toPoly source).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly source).leadingCoeff : ZMod prime) ≠ 0)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod prime))
        (SparsePolyZZ.toPoly source)))
    (activeState : StrictRecombine.LiveActiveFactors prime active)
    (productState : LiveHenselProduct prime exponent source active)
    (precision : LiveRecoveryPrecision (prime ^ exponent) source)
    (resultIrreducible : FactorArrayIrreducible result)
    (history : StrictRecombine.SmallerZassenhausScansExhausted source active
      (((prime ^ exponent : Nat) : ZZ)) subsetSize) :
    ∃ output,
      Generated.StrictRecombine.zassenhausLoop
        StrictRecombine.concreteZassenhausTermination
        (((prime ^ exponent : Nat) : ZZ)) active source result subsetSize
        hsubsetPositive = .ok output ∧
      Nonempty (ZassenhausTerminalCertificate prime exponent output) := by
  rw [Generated.StrictRecombine.zassenhausLoop]
  split
  next hcontinue =>
    let initial := Generated.StrictRecombine.initialCombination subsetSize
    have hfits : subsetSize ≤ active.size := by omega
    let hinitial :=
      StrictRecombine.concreteZassenhausTermination.initial_valid
        active.size subsetSize hfits
    rcases StrictRecombine.scanZassenhausCombinations_complete source active
        (prime ^ exponent) prime initial
        (pow_pos (Fact.out : Nat.Prime prime).pos exponent)
        (Fact.out : Nat.Prime prime).pos
        (dvd_pow_self prime (Nat.ne_of_gt productState.exponentPositive))
        hcanonical hnonempty hsubsetPositive activeState.fitsInt32
        (by
          have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head
            source hcanonical hnonempty
          simpa [getElem!_pos source 0 hnonempty] using (hhead ▸ hleading))
        activeState.irreducible hinitial with ⟨scanResult, hscan⟩
    cases scanResult with
    | exhausted =>
        simp only [StrictRecombine.concreteZassenhausTermination, initial,
          hinitial, hscan]
        have hcurrent : StrictRecombine.FixedSizeScanExhausted source active
            (((prime ^ exponent : Nat) : ZZ)) subsetSize :=
          ⟨hfits, hscan⟩
        exact zassenhausLoop_live_terminal prime exponent active source result
          (subsetSize + 1) (by omega) hcanonical hnonempty hprimitive hleading
          hsquarefree activeState productState precision resultIrreducible
          (history.succ hcurrent)
    | extracted factor quotient candidate candidateSize =>
        simp only [StrictRecombine.concreteZassenhausTermination, initial,
          hinitial, hscan]
        rcases scanExtraction_liveStep prime exponent subsetSize source factor
            quotient active candidate candidateSize productState activeState
            precision history hsubsetPositive hfits hcanonical hnonempty
            hprimitive hsquarefree hleading hscan with
          ⟨remaining, hremove, hfactorIrreducible, hquotientCanonical,
            hquotientPrimitive, hquotientNonempty, hquotientLeading,
            hactiveNext, hproductNext, hprecisionNext, hsquarefreeNext⟩
        rw [hremove]
        exact zassenhausLoop_live_terminal prime exponent remaining quotient
          (result.push factor) 1 (by omega) hquotientCanonical
          hquotientNonempty hquotientPrimitive hquotientLeading
          hsquarefreeNext hactiveNext hproductNext hprecisionNext
          (resultIrreducible.push hfactorIrreducible)
          (StrictRecombine.SmallerZassenhausScansExhausted.one quotient
            remaining (((prime ^ exponent : Nat) : ZZ)))
  next hcontinue =>
    refine ⟨Generated.StrictRecombine.finishZassenhaus source result, rfl, ?_⟩
    exact ⟨{
      active := active
      source := source
      result := result
      subsetSize := subsetSize
      output_eq := rfl
      stopped := hcontinue
      subsetPositive := hsubsetPositive
      canonical := hcanonical
      nonempty := hnonempty
      primitive := hprimitive
      leading := hleading
      squarefree := hsquarefree
      activeState := activeState
      productState := productState
      precision := precision
      resultIrreducible := resultIrreducible
      history := history }⟩
termination_by (active.size, active.size + 1 - subsetSize)
decreasing_by
  · exact Prod.Lex.right _ (by omega)
  · exact Prod.Lex.left _ _
      (StrictRecombine.concreteZassenhausTermination.removal_decreases active
        candidate remaining (by rw [candidateSize]; exact hsubsetPositive)
        hremove)

set_option maxHeartbeats 1000000 in
/-- The polynomial left in the literal terminal Zassenhaus state is
irreducible whenever it has positive degree.  A hypothetical factorization is
first represented by a legal candidate over the physical live Hensel array.
The generated reverse-erasure loop supplies the occurrence-sensitive
complement candidate.  One of these two candidates has size at most half of
the active array, hence is covered by the recorded literal scan history; full
recovery precision says that same generated attempt extracts it, a
contradiction. -/
theorem ZassenhausTerminalCertificate.source_irreducible
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {output : Array SparsePolyZZ}
    (certificate : ZassenhausTerminalCertificate prime exponent output)
    (hdegree : 0 < certificate.source[0]!.1.deg) :
    Irreducible (SparsePolyZZ.toPoly certificate.source) := by
  let source := certificate.source
  let active := certificate.active
  let sourcePoly := SparsePolyZZ.toPoly source
  have hsourceNonempty : 0 < source.size := by
    simpa [source] using certificate.nonempty
  have hsourceNonunit : ¬IsUnit sourcePoly := by
    have hpolyDegree : (SparsePolyZZ.toPoly source).degree =
        some source[0]!.1.deg := by
      exact StrictRecombine.sparsePolyZZ_toPoly_degree_eq_head source
        certificate.canonical hsourceNonempty
    apply Polynomial.not_isUnit_of_degree_pos sourcePoly
    change 0 < (SparsePolyZZ.toPoly source).degree
    rw [hpolyDegree]
    exact WithBot.coe_lt_coe.mpr (by simpa [source] using hdegree)
  refine ⟨hsourceNonunit, ?_⟩
  intro left right hfactorization
  by_cases hleftUnit : IsUnit left
  · exact Or.inl hleftUnit
  by_cases hrightUnit : IsUnit right
  · exact Or.inr hrightUnit
  exfalso
  rcases primitive_factorization_maps_nonunits prime sourcePoly left right
      certificate.primitive certificate.leading hfactorization hleftUnit
      hrightUnit with
    ⟨hleftPrimitive, hrightPrimitive, hleftLeading, hrightLeading,
      hleftMapNonunit, hrightMapNonunit⟩
  have hleftNe : left ≠ 0 := Polynomial.leadingCoeff_ne_zero.mp
    (fun hzero => hleftLeading (by simpa [hzero]))
  have hrightNe : right ≠ 0 := Polynomial.leadingCoeff_ne_zero.mp
    (fun hzero => hrightLeading (by simpa [hzero]))
  have hleftMapNe : left.map (Int.castRingHom (ZMod prime)) ≠ 0 := by
    intro hzero
    apply hleftLeading
    have hlead := Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero
      (Int.castRingHom (ZMod prime)) hleftLeading
    rw [hzero] at hlead
    simpa using hlead.symm
  have hrightMapNe : right.map (Int.castRingHom (ZMod prime)) ≠ 0 := by
    intro hzero
    apply hrightLeading
    have hlead := Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero
      (Int.castRingHom (ZMod prime)) hrightLeading
    rw [hzero] at hlead
    simpa using hlead.symm
  rcases live_divisor_mod_has_legal_candidate prime source active
      certificate.productState.primeProductAssociated
      certificate.activeState.irreducible left
      ⟨right, hfactorization⟩ with
    ⟨leftCandidate, hleftLegal, hleftAssociated⟩
  rcases StrictRecombine.removeCombination_complement_candidate leftCandidate
      active hleftLegal with
    ⟨remaining, rightCandidate, hremove, hrightLegal, hrightSelected,
      hsize⟩
  have hpartition :=
    StrictRecombine.removeCombination_toPolyMod_product_partition prime
      leftCandidate active remaining hleftLegal hremove
  have hrightLegal' : StrictRecombine.LegalCombination active.size
      rightCandidate.size rightCandidate := by
    simpa [hrightLegal.1] using hrightLegal
  have hrightAssociated : Associated
      (right.map (Int.castRingHom (ZMod prime)))
      (((StrictRecombine.selectSourceIndices active.toList
        rightCandidate.toList).map (StrictHensel.toPolyMod prime)).prod) := by
    have hsourceMapped :
        sourcePoly.map (Int.castRingHom (ZMod prime)) =
          left.map (Int.castRingHom (ZMod prime)) *
            right.map (Int.castRingHom (ZMod prime)) := by
      simpa [sourcePoly, source, Polynomial.map_mul] using
        congrArg (Polynomial.map (Int.castRingHom (ZMod prime)))
          hfactorization
    have htotal : Associated
        (left.map (Int.castRingHom (ZMod prime)) *
          right.map (Int.castRingHom (ZMod prime)))
        (((StrictRecombine.selectSourceIndices active.toList
          leftCandidate.toList).map (StrictHensel.toPolyMod prime)).prod *
          (remaining.toList.map (StrictHensel.toPolyMod prime)).prod) :=
      (Associated.of_eq hsourceMapped.symm).trans
        (certificate.productState.primeProductAssociated.trans
          (Associated.of_eq hpartition.symm))
    have hremainingAssociated : Associated
        (right.map (Int.castRingHom (ZMod prime)))
        ((remaining.toList.map (StrictHensel.toPolyMod prime)).prod) :=
      Associated.of_mul_left htotal hleftAssociated hleftMapNe
    simpa [hrightSelected] using hremainingAssociated
  have hleftPositive : 0 < leftCandidate.size := by
    by_contra hnot
    have hempty := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hnot)
    subst leftCandidate
    apply hleftMapNonunit
    simpa [StrictRecombine.selectSourceIndices] using
      hleftAssociated.isUnit_iff.mpr isUnit_one
  have hrightPositive : 0 < rightCandidate.size := by
    by_contra hnot
    have hempty := Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hnot)
    subst rightCandidate
    apply hrightMapNonunit
    simpa [StrictRecombine.selectSourceIndices] using
      hrightAssociated.isUnit_iff.mpr isUnit_one
  let leading := source[0]!
  have hleadingCell : source[0]? = some leading := by
    rw [Array.getElem?_eq_getElem hsourceNonempty]
    simp [leading, getElem!_pos source 0 hsourceNonempty]
  have hleadingEq : leading.2 = sourcePoly.leadingCoeff := by
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head source
      certificate.canonical hsourceNonempty
    simpa [leading, getElem!_pos source 0 hsourceNonempty] using hhead.symm
  have hleadingBound : leading.2.natAbs * 2 < prime ^ exponent := by
    simpa [hleadingEq] using certificate.precision.leadingBound
  by_cases hleftSmall : 2 * leftCandidate.size ≤ active.size
  · have hscanSmall : leftCandidate.size < certificate.subsetSize := by
      have hstopped : active.size < 2 * certificate.subsetSize := by
        simpa [active] using Nat.lt_of_not_ge certificate.stopped
      omega
    have hleftLeadingMap :
        (left.map (Int.castRingHom (ZMod prime))).leadingCoeff =
          (left.leadingCoeff : ZMod prime) :=
      Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ hleftLeading
    have hcongruent := certificate.productState.selectedCongruent
      certificate.activeState certificate.squarefree certificate.leading
      left right hfactorization hleftMapNe hleftLeadingMap leftCandidate
      hleftLegal hleftAssociated
    rcases zassenhausAttempt_extracts_live_divisor_candidate
        certificate.canonical certificate.nonempty leading hleadingCell left
        right hfactorization hleftPrimitive hleftNe hrightNe leftCandidate
        hleftLegal certificate.activeState.fitsInt32
        certificate.activeState.canonical certificate.activeState.nonempty
        certificate.activeState.monic (prime ^ exponent)
        (pow_pos (Fact.out : Nat.Prime prime).pos exponent) rfl hleadingBound
        hcongruent
        (certificate.precision.scaledFactorBound left right hfactorization
          hleftNe) with
      ⟨factor, quotient, hextracted⟩
    have hrejected := certificate.history.rejects leftCandidate hleftPositive
      hscanSmall hleftLegal
    rw [hextracted] at hrejected
    simp at hrejected
  · have hrightSmallBound : 2 * rightCandidate.size ≤ active.size := by
      have hremainingSize : remaining.size = rightCandidate.size :=
        hrightLegal.1.symm
      rw [← hremainingSize]
      omega
    have hscanSmall : rightCandidate.size < certificate.subsetSize := by
      have hstopped : active.size < 2 * certificate.subsetSize := by
        simpa [active] using Nat.lt_of_not_ge certificate.stopped
      omega
    have hrightLeadingMap :
        (right.map (Int.castRingHom (ZMod prime))).leadingCoeff =
          (right.leadingCoeff : ZMod prime) :=
      Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ hrightLeading
    have hfactorization' : sourcePoly = right * left := by
      simpa [sourcePoly, source, mul_comm] using hfactorization
    have hcongruent := certificate.productState.selectedCongruent
      certificate.activeState certificate.squarefree certificate.leading
      right left hfactorization' hrightMapNe hrightLeadingMap rightCandidate
      hrightLegal' hrightAssociated
    rcases zassenhausAttempt_extracts_live_divisor_candidate
        certificate.canonical certificate.nonempty leading hleadingCell right
        left hfactorization' hrightPrimitive hrightNe hleftNe rightCandidate
        hrightLegal' certificate.activeState.fitsInt32
        certificate.activeState.canonical certificate.activeState.nonempty
        certificate.activeState.monic (prime ^ exponent)
        (pow_pos (Fact.out : Nat.Prime prime).pos exponent) rfl hleadingBound
        hcongruent
        (certificate.precision.scaledFactorBound right left hfactorization'
          hrightNe) with
      ⟨factor, quotient, hextracted⟩
    have hrejected := certificate.history.rejects rightCandidate
      hrightPositive hscanSmall hrightLegal'
    rw [hextracted] at hrejected
    simp at hrejected

/-- Every physical factor returned by the literal terminal finish block is
irreducible: previously extracted factors use the accumulated execution
invariant, and the positive-degree remainder uses
`source_irreducible`. -/
theorem ZassenhausTerminalCertificate.output_irreducible
    {prime exponent : Nat} [Fact (Nat.Prime prime)]
    {output : Array SparsePolyZZ}
    (certificate : ZassenhausTerminalCertificate prime exponent output) :
    FactorArrayIrreducible output := by
  rw [certificate.output_eq]
  exact certificate.resultIrreducible.finishZassenhaus
    (fun _hnonempty hdegree => certificate.source_irreducible hdegree)

/-- End-to-end correctness of one concrete Hensel execution followed by an
accepted full-cardinality concrete van-Hoeij execution.  This is the exact
non-safety-net branch used by `__lll_factorize_raw_ir`. -/
theorem selectionHensel_vanHoeij_equal_cardinality_refines_FactorZZCorrect
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {henselOutput : Array SparsePolyZZ × ZZ}
    {output : Array SparsePolyZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget henselOutput)
    (hfits : henselOutput.1.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (get_deg f).toNatClampNeg)
    (hlength : output.size = henselOutput.1.size)
    (hrun : Generated.StrictRecombine.__vanhoeij_recombine_raw_ir
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination f henselOutput.1
      henselOutput.2 hdegree = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  rcases selectionHenselFactors_liveProduct hcount hfactors
      hleadingSemantic hselection hentry with
    ⟨exponent, hmodulus, productState⟩
  have hirreducible := selectionHenselFactors_mod_irreducible hcount
    hfactors hleadingSemantic hselection hentry
  have hleading := hselection.goodPrime.lc_nonzero
  have hrun' : Generated.StrictRecombine.__vanhoeij_recombine_raw_ir
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination f henselOutput.1
      (((selection.prime.toNat ^ exponent : Nat) : ZZ)) hdegree =
        .ok output := by
    simpa [hmodulus] using hrun
  have harrayIrreducible :=
    concreteVanHoeij_equal_cardinality_factorArrayIrreducible productState
      hdegree hcanonical hnonempty hprimitive hleading hfits hirreducible
      hlength hrun'
  have hproduct :=
    StrictRecombine.__vanhoeij_recombine_raw_ir_product_associated
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination f henselOutput.1
      henselOutput.2 hdegree hcanonical hprimitive output hrun
  constructor
  · simpa [StrictRecombine.factorArrayProduct] using hproduct
  · intro factor hfactor
    rcases List.mem_map.mp hfactor with ⟨physical, hphysical, rfl⟩
    exact harrayIrreducible physical hphysical

/-- End-to-end correctness of the concrete full-precision Hensel followed by
the literal generated Zassenhaus recombination entry.  Both product recovery
and irreducibility refer to the same returned physical array. -/
theorem selectionHensel_zassenhausRecombine_refines_FactorZZCorrect_of_recovery
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32}
    {henselOutput : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget henselOutput)
    (hrecovery : ∀ exponent,
      henselOutput.2 =
          ((selection.prime.toNat ^ exponent : Nat) : Int) →
        LiveRecoveryPrecision (selection.prime.toNat ^ exponent) f)
    (hfits : henselOutput.1.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.zassenhausRecombine
      StrictRecombine.concreteZassenhausTermination f henselOutput.1
        henselOutput.2 = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  have hsize := selectionHenselFactors_pointwise_associated hcount
    hfactors hleadingSemantic hselection hentry
  have hlifted : ¬henselOutput.1.size ≤ 1 := by omega
  rcases selectionHenselFactors_liveProduct hcount hfactors
      hleadingSemantic hselection hentry with
    ⟨exponent, hmodulus, productState⟩
  have activeState := selectionHenselFactors_liveActive hcount hfactors
    hleadingSemantic hselection hentry hfits
  have precision := hrecovery exponent hmodulus
  have hloopRun : Generated.StrictRecombine.zassenhausLoop
      StrictRecombine.concreteZassenhausTermination
      (((selection.prime.toNat ^ exponent : Nat) : ZZ)) henselOutput.1 f #[] 1
      (by omega) = .ok output := by
    unfold Generated.StrictRecombine.zassenhausRecombine at hrun
    rw [if_neg hlifted] at hrun
    simpa [hmodulus] using hrun
  rcases zassenhausLoop_live_terminal selection.prime.toNat exponent
      henselOutput.1 f #[] 1 (by omega) hcanonical hnonempty
      hprimitive hselection.goodPrime.lc_nonzero
      hselection.goodPrime.sqfree activeState productState precision
      FactorArrayIrreducible.empty
      (StrictRecombine.SmallerZassenhausScansExhausted.one f henselOutput.1
        (((selection.prime.toNat ^ exponent : Nat) : ZZ))) with
    ⟨terminalOutput, hterminalRun, hterminalCertificate⟩
  have houtputEq : output = terminalOutput :=
    Except.ok.inj (hloopRun.symm.trans hterminalRun)
  subst terminalOutput
  rcases hterminalCertificate with ⟨certificate⟩
  refine ⟨?_, ?_⟩
  · exact StrictRecombine.zassenhausRecombine_toPoly_product_associated
      StrictRecombine.concreteZassenhausTermination f henselOutput.1 output
      henselOutput.2 hcanonical hprimitive hrun
  · intro factor hfactor
    rcases List.mem_map.mp hfactor with ⟨physical, hphysical, rfl⟩
    exact certificate.output_irreducible physical hphysical

/-- Default-Mignotte specialization used by the literal second Hensel call
whose source argument is zero. -/
theorem selectionHensel_zassenhausRecombine_refines_FactorZZCorrect
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {henselOutput : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 henselOutput)
    (hfits : henselOutput.1.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.zassenhausRecombine
      StrictRecombine.concreteZassenhausTermination f henselOutput.1
        henselOutput.2 = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  apply selectionHensel_zassenhausRecombine_refines_FactorZZCorrect_of_recovery
    hcount hfactors hleadingSemantic hselection hentry
    (fun exponent hmodulus =>
      selectionHenselFactors_liveRecoveryPrecision hentry hcanonical
    hnonempty hdegree leading hleading exponent hmodulus)
    hfits hcanonical hprimitive hnonempty hdegree leading hleading output hrun

/-- The full-precision first-pass branch needs no semantic recovery oracle:
the generated heuristic loop crosses its own Mignotte target, and the actual
generated Hensel loop returns a prime power beyond the corresponding explicit
target. -/
theorem heuristic_full_hensel_liveRecoveryPrecision
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aH aMig : Int32} {output : Array SparsePolyZZ × ZZ}
    (r : Int32)
    (hheuristic :
      Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir f r p =
        .ok (aH, aMig))
    (hfull : aMig ≤ aH)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p aH
      output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (exponent : Nat)
    (houtput : output.2 = ((p.toNat ^ exponent : Nat) : Int)) :
    LiveRecoveryPrecision (p.toNat ^ exponent) f := by
  have hequal : aH = aMig := by
    have hle := heuristic_starting_precision_first_le_second f r p aH aMig
      hheuristic
    cases aH with
    | ofUInt32 aHUnsigned =>
      cases aMig with
      | ofUInt32 aMigUnsigned =>
        cases aHUnsigned with
        | ofBitVec aHBits =>
          cases aMigUnsigned with
          | ofBitVec aMigBits =>
            congr 3
            apply BitVec.eq_of_toInt_eq
            rw [Int32.le_iff_toInt_le] at hle hfull
            exact le_antisymm hle hfull
  by_cases hzero : aH = 0
  · have hentryZero : StrictHensel.HenselLiftEntryCorrect termination f
        factors p 0 output := by simpa [hzero] using hentry
    exact selectionHenselFactors_liveRecoveryPrecision hentryZero hcanonical
      hnonempty hdegree leading hleading exponent houtput
  · have hpositive : 0 < aH := by
      rcases hentry with
        ⟨target, adjusted, nodes, liftedNodes, outputM, extracted,
          htarget, hadjust, hsemantic, hlift, hliftedOneHead, hextract,
          hnormalize, houtputM, houtputCanonical, houtputOneHead⟩
      cases htarget with
      | mignotte sourceLeading htargetZero hsourceLeading =>
          exact False.elim (hzero htargetZero)
      | explicit htargetPositive => exact htargetPositive
    have hnonnegativeMig : 0 ≤ aMig := by
      rw [← hequal]
      exact Int32.le_of_lt hpositive
    have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
      intro hsourceZero
      have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
        hcanonical hnonempty
      rw [hsourceZero] at hhead
      exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
        (by simpa using hhead.symm)
    rcases StrictRecombine.mignotteBoundRaw_bounds_divisor f hcanonical
        hnonempty hdegree (1 : Polynomial Int) hsourceNe (one_dvd _) with
      ⟨bound, hboundRun, _hboundNonnegative, _hboundOne⟩
    have hheuristicLarge := heuristic_starting_precision_full_pow_gt
      f r p aH aMig leading bound hleading hboundRun hheuristic
        hnonnegativeMig
    have htargetNonnegative : ∀ target,
        StrictHensel.HenselLiftTargetCorrect f p aH target → 0 ≤ target := by
      intro target htarget
      cases htarget with
      | mignotte sourceLeading htargetZero hsourceLeading =>
          exact False.elim (hzero htargetZero)
      | explicit htargetPositive =>
          have honeNat : 1 ≤ p.toNat ^ aH.toNatClampNeg :=
            Nat.one_le_pow _ p.toNat (Fact.out : Nat.Prime p.toNat).pos
          have honeInt : (1 : Int) ≤
              ((p.toNat ^ aH.toNatClampNeg : Nat) : Int) := by
            exact_mod_cast honeNat
          change 0 ≤ ((p.toNat ^ aH.toNatClampNeg : Nat) : Int) - 1
          omega
    rcases hentry.outputModulus_gt_target htargetNonnegative with
      ⟨target, htarget, htargetLt⟩
    have hpowerLeOutput :
        ((p.toNat ^ aH.toNatClampNeg : Nat) : Int) ≤ output.2 := by
      cases htarget with
      | mignotte sourceLeading htargetZero hsourceLeading =>
          exact False.elim (hzero htargetZero)
      | explicit htargetPositive =>
          change ((p.toNat ^ aH.toNatClampNeg : Nat) : Int) - 1 <
            output.2 at htargetLt
          omega
    have hpowerLe : p.toNat ^ aMig.toNatClampNeg ≤ p.toNat ^ exponent := by
      rw [← hequal]
      rw [houtput] at hpowerLeOutput
      exact_mod_cast hpowerLeOutput
    apply liveRecoveryPrecision_of_generated_mignotte_lt f hcanonical hnonempty
      hdegree leading hleading (p.toNat ^ exponent)
    intro candidateBound hcandidateBoundRun
    have hsame : candidateBound = bound := by
      rw [hboundRun] at hcandidateBoundRun
      exact (Except.ok.inj hcandidateBoundRun).symm
    subst candidateBound
    exact hheuristicLarge.trans_le hpowerLe

/-- Correctness of the literal generated `__lll_factorize_raw_ir` controller
with concrete recombination callees.  The Hensel premises certify only arrays
actually returned by the supplied raw Hensel function.  Full-precision
coefficient recovery is derived from the generated heuristic and Hensel
executions rather than supplied as a semantic premise. -/
theorem __lll_factorize_raw_ir_refines_FactorZZCorrect
    (henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
      RawExec (Array SparsePolyZZ × ZZ))
    (f : SparsePolyZZ) (selection : PrimeSelectionResult)
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : ∀ aTarget henselOutput,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
        StrictHensel.HenselLiftEntryCorrect
          StrictHensel.concreteDivmodTermination f selection.factors
          selection.prime aTarget henselOutput)
    (hfactorFitsOfRun : ∀ aTarget henselOutput,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
        selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size) (hcellDegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictFactorZZ.__lll_factorize_raw_ir
      (concreteRecombineFactorZZRawOps (fun _ _ => .error .assertionFailure)
        henselLift) f selection.factors selection.prime = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  have henselSize : ∀ aTarget henselOutput,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
        henselOutput.1.size = selection.factors.size := by
    intro aTarget henselOutput hh
    exact (selectionHenselFactors_pointwise_associated hcount hfactors
      hleadingSemantic hselection (hentry aTarget henselOutput hh)).1.symm
  have vanSize : ∀ lifted modulus result,
      concreteVanHoeijRecombine f lifted modulus = .ok result →
      0 < lifted.size → result.size ≤ lifted.size := by
    intro lifted modulus result hv hlifted
    rcases concreteVanHoeijRecombine_success f lifted result modulus hv with
      ⟨hdegree, hraw⟩
    exact StrictRecombine.__vanhoeij_recombine_raw_ir_size_le
      StrictRecombine.concreteVanHoeijRawOps
      StrictRecombine.concreteVanHoeijTermination f lifted modulus hdegree
      hlifted result hraw
  have vanCorrect : ∀ aTarget henselOutput result,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
      concreteVanHoeijRecombine f henselOutput.1 henselOutput.2 = .ok result →
      result.size = selection.factors.size →
      FactorZZCorrect (SparsePolyZZ.toPoly f)
        (result.toList.map SparsePolyZZ.toPoly) := by
    intro aTarget henselOutput result hh hv hsize
    have hfactorFits := hfactorFitsOfRun aTarget henselOutput hh
    rcases concreteVanHoeijRecombine_success f henselOutput.1 result
        henselOutput.2 hv with ⟨hdegree, hraw⟩
    apply selectionHensel_vanHoeij_equal_cardinality_refines_FactorZZCorrect
      hcount hfactors hleadingSemantic hselection
      (hentry aTarget henselOutput hh) (by rw [henselSize aTarget henselOutput hh]; omega)
      hcanonical hprimitive hnonempty hdegree
      (by rw [henselSize aTarget henselOutput hh]; exact hsize) hraw
  have zassenhausCorrect : ∀ aTarget henselOutput result,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
      (∀ exponent, henselOutput.2 =
          ((selection.prime.toNat ^ exponent : Nat) : Int) →
        LiveRecoveryPrecision (selection.prime.toNat ^ exponent) f) →
      concreteZassenhausRecombine f henselOutput.1 henselOutput.2 = .ok result →
      FactorZZCorrect (SparsePolyZZ.toPoly f)
        (result.toList.map SparsePolyZZ.toPoly) := by
    intro aTarget henselOutput result hh hrecovery hz
    have hfactorFits := hfactorFitsOfRun aTarget henselOutput hh
    apply selectionHensel_zassenhausRecombine_refines_FactorZZCorrect_of_recovery
      hcount hfactors hleadingSemantic hselection
      (hentry aTarget henselOutput hh) hrecovery
      (by rw [henselSize aTarget henselOutput hh]; omega) hcanonical hprimitive
      hnonempty hcellDegree leading hleading result
    exact hz
  unfold Generated.StrictFactorZZ.__lll_factorize_raw_ir at hrun
  simp only [concreteRecombineFactorZZRawOps] at hrun
  cases hheuristic :
      Generated.StrictFactorZZ.__heuristic_starting_precision_raw_ir f
        selection.factors.size.toUInt32.toInt32 selection.prime with
  | error fault => rw [hheuristic] at hrun; contradiction
  | ok precision =>
    rcases precision with ⟨aH, aMig⟩
    rw [hheuristic] at hrun
    cases hlift : henselLift f selection.factors selection.prime aH with
    | error fault => simp [hlift] at hrun
    | ok henselH =>
      rcases henselH with ⟨liftedH, mH⟩
      simp only [hlift] at hrun
      have hfactorFits : selection.factors.size < 2 ^ 31 :=
        hfactorFitsOfRun aH (liftedH, mH) hlift
      cases hvan : concreteVanHoeijRecombine f liftedH mH with
      | error fault => simp [hvan] at hrun
      | ok result =>
        simp only [hvan] at hrun
        by_cases hlow :
            (decide (result.size.toUInt32.toInt32 <
                selection.factors.size.toUInt32.toInt32) &&
              decide (aH < aMig)) = true
        ·
          rw [if_pos hlow] at hrun
          cases hliftMig : henselLift f selection.factors selection.prime 0 with
          | error fault => simp [hliftMig] at hrun
          | ok henselMig =>
            rcases henselMig with ⟨liftedMig, mMig⟩
            simp only [hliftMig] at hrun
            cases hvanMig : concreteVanHoeijRecombine f liftedMig mMig with
            | error fault => simp [hvanMig] at hrun
            | ok resultMig =>
              simp only [hvanMig] at hrun
              by_cases hsafety :
                  Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir
                    resultMig.size selection.factors.size true = true
              · rw [if_pos hsafety] at hrun
                exact zassenhausCorrect 0 (liftedMig, mMig) output hliftMig
                  (fun exponent hmodulus =>
                    selectionHenselFactors_liveRecoveryPrecision
                      (hentry 0 (liftedMig, mMig) hliftMig) hcanonical
                      hnonempty hcellDegree leading hleading exponent hmodulus)
                  hrun
              · rw [if_neg hsafety] at hrun
                have hout := Except.ok.inj hrun
                subst output
                have hliftedSize := henselSize 0 (liftedMig, mMig) hliftMig
                have hliftedSize' : liftedMig.size = selection.factors.size := by
                  simpa using hliftedSize
                have hle := vanSize liftedMig mMig resultMig hvanMig (by
                  rw [hliftedSize']
                  omega)
                have hleFactors : resultMig.size ≤ selection.factors.size :=
                  hle.trans_eq hliftedSize'
                have hnotLess : ¬ resultMig.size < selection.factors.size := by
                  simpa [Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir]
                    using hsafety
                exact vanCorrect 0 (liftedMig, mMig) resultMig hliftMig hvanMig
                  (by omega)
        · rw [if_neg hlow] at hrun
          by_cases hsafety :
              Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir
                result.size selection.factors.size (aMig ≤ aH) = true
          · rw [if_pos hsafety] at hrun
            exact zassenhausCorrect aH (liftedH, mH) output hlift
              (fun exponent hmodulus =>
                heuristic_full_hensel_liveRecoveryPrecision
                  selection.factors.size.toUInt32.toInt32 hheuristic
                  (by
                    simp [Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir]
                      at hsafety
                    exact hsafety.1)
                  (hentry aH (liftedH, mH) hlift) hcanonical hnonempty
                  hcellDegree leading hleading exponent hmodulus) hrun
          · rw [if_neg hsafety] at hrun
            have hout := Except.ok.inj hrun
            subst output
            have hliftedSize := henselSize aH (liftedH, mH) hlift
            have hliftedSize' : liftedH.size = selection.factors.size := by
              simpa using hliftedSize
            have hle := vanSize liftedH mH result hvan (by
              rw [hliftedSize']
              omega)
            have hleFactors : result.size ≤ selection.factors.size :=
              hle.trans_eq hliftedSize'
            have hnotLess : ¬ result.size < selection.factors.size := by
              by_contra hless
              have hmachine : result.size.toUInt32.toInt32 <
                  selection.factors.size.toUInt32.toInt32 := by
                rw [size_toUInt32_toInt32_lt_iff result.size
                  selection.factors.size (by omega) hfactorFits]
                exact hless
              have hprecision := heuristic_starting_precision_first_le_second
                f selection.factors.size.toUInt32.toInt32 selection.prime aH
                aMig hheuristic
              by_cases hstrict : aH < aMig
              · have hmachine' : Int32.ofNat result.size <
                    Int32.ofNat selection.factors.size := by
                  simpa using hmachine
                apply hlow
                simpa [hmachine, hmachine', hstrict]
              · have hfull : aMig ≤ aH := Int32.not_lt.mp hstrict
                exact hsafety (by
                  simp [Generated.StrictFactorZZ.__needs_zassenhaus_safety_net_ir,
                    hfull, hless])
            exact vanCorrect aH (liftedH, mH) result hlift hvan (by omega)

/-- Exact outer control-flow composition for the original C++
`__factor_squarefree_primitive_ZZ` entry.  The selected value and every
Hensel value are constrained only after they have been returned by their raw
executions; no premise can choose the final factor array. -/
theorem __factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect
    (selectPrime : SparsePolyZZ → Bool → RawExec PrimeSelectionResult)
    (henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
      RawExec (Array SparsePolyZZ × ZZ))
    (useLargePrime : Bool) (f : SparsePolyZZ)
    (hselect : ∀ selection,
      selectPrime f useLargePrime = .ok selection →
      StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection ∧
        StrictSelectPrime.SelectionPhysical selection)
    (hirreducibleCount : ∀ selection,
      selectPrime f useLargePrime = .ok selection →
      selection.irreducible = true → selection.factors.size ≤ 1)
    (hhensel : ∀ (selection : PrimeSelectionResult)
      (hprime : Nat.Prime selection.prime.toNat),
      ∀ aTarget henselOutput,
        henselLift f selection.factors selection.prime aTarget =
            .ok henselOutput →
          @StrictHensel.HenselLiftEntryCorrect
            StrictHensel.concreteDivmodTermination f selection.factors
            selection.prime ⟨hprime⟩
            aTarget henselOutput)
    (hhenselFits : ∀ (selection : PrimeSelectionResult)
      (hprime : Nat.Prime selection.prime.toNat) aTarget henselOutput,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
        selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegreeBound : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir
      (concreteRecombineFactorZZRawOps selectPrime henselLift)
      useLargePrime f = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  unfold Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir at hrun
  have hguard : ¬(f.isEmpty || get_deg f < 2) := by
    intro hbad
    rw [if_pos hbad] at hrun
    contradiction
  rw [if_neg hguard] at hrun
  simp only [concreteRecombineFactorZZRawOps] at hrun
  cases hselectionRun : selectPrime f useLargePrime with
  | error fault => rw [hselectionRun] at hrun; contradiction
  | ok selection =>
    rw [hselectionRun] at hrun
    simp only at hrun
    have hselectionResult := hselect selection hselectionRun
    have hselection := hselectionResult.1
    have hselectedCanonical := hselectionResult.2
    letI : Fact (Nat.Prime selection.prime.toNat) :=
      ⟨hselection.goodPrime.prime⟩
    by_cases hsingle : selection.irreducible || selection.factors.size ≤ 1
    · rw [if_pos hsingle] at hrun
      have houtput : output = #[f] := Except.ok.inj hrun |>.symm
      subst output
      have hcount : selection.factors.size ≤ 1 := by
        have hsingleProp : selection.irreducible = true ∨
            selection.factors.size ≤ 1 := by
          simpa [Bool.or_eq_true] using hsingle
        rcases hsingleProp with hirreducible | hcount
        · exact hirreducibleCount selection hselectionRun hirreducible
        · exact hcount
      simpa using selection_atMostOne_refines_singleton_FactorZZCorrect
        hselection hprimitive hdegree hcount
    · rw [if_neg hsingle] at hrun
      have hcount : 2 ≤ selection.factors.size := by
        have hnot : ¬ selection.factors.size ≤ 1 := by
          intro hle
          apply hsingle
          simp [hle]
        omega
      have hleadingSemantic : ∀ sourceLeading, f[0]? = some sourceLeading →
          (sourceLeading.2 : ZMod selection.prime.toNat) =
            (SparsePolyZZ.toPoly f).leadingCoeff := by
        intro sourceLeading hsourceLeading
        have hfront : f[0] = sourceLeading := by
          rw [Array.getElem?_eq_getElem hnonempty] at hsourceLeading
          exact Option.some.inj hsourceLeading
        rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
          hnonempty, hfront]
      apply __lll_factorize_raw_ir_refines_FactorZZCorrect henselLift f
        selection hcount hselectedCanonical hleadingSemantic hselection
        (fun aTarget henselOutput hh =>
          hhensel selection hselection.goodPrime.prime
            aTarget henselOutput hh)
        (hhenselFits selection hselection.goodPrime.prime)
        hcanonical hprimitive hnonempty
        hdegreeBound leading hleading output
      exact hrun

/-- Outer refinement with the select-prime field instantiated by the actual
generated prime iterator, polynomial reduction, GCD, DDF and EDF execution.
Only downstream physical Hensel readiness remains explicit. -/
theorem concreteSelect___factor_squarefree_primitive_ZZ_raw_ir_refines
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State)
    (henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
      RawExec (Array SparsePolyZZ × ZZ))
    (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hhensel : ∀ (selection : PrimeSelectionResult)
      (hprime : Nat.Prime selection.prime.toNat),
      ∀ aTarget henselOutput,
        henselLift f selection.factors selection.prime aTarget =
            .ok henselOutput →
          @StrictHensel.HenselLiftEntryCorrect
            StrictHensel.concreteDivmodTermination f selection.factors
            selection.prime ⟨hprime⟩ aTarget henselOutput)
    (hhenselFits : ∀ (selection : PrimeSelectionResult)
      (hprime : Nat.Prime selection.prime.toNat) aTarget henselOutput,
      henselLift f selection.factors selection.prime aTarget =
          .ok henselOutput →
        selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegree62 : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hdegree63 : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir
      (concreteRecombineFactorZZRawOps
        (concreteSelectPrime engine provider initialRng) henselLift)
      useLargePrime f = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  have hlcSemantic : ∀ p : UInt64, Nat.Prime p.toNat →
      ((SparsePolyZZ.front! f).2 : ZMod p.toNat) =
        ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat) := by
    intro p hp
    have hvalue : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    have hfront : SparsePolyZZ.front! f = leading := by
      rw [SparsePolyZZ.front!, getElem!_pos f 0 hnonempty]
      exact hvalue
    rw [hfront, StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
      hcanonical hnonempty, hvalue]
  apply __factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect
    (concreteSelectPrime engine provider initialRng) henselLift useLargePrime f
    (fun selection hselectionRun =>
      concreteSelectPrime_success engine provider initialRng useLargePrime f
        hinitialPrimeCorrect hcanonical hnonempty hdegree hdegree62 hlcSemantic selection
          hselectionRun)
    (concreteSelectPrime_irreducible_size engine provider initialRng
      useLargePrime f)
    hhensel hhenselFits hcanonical hprimitive hnonempty hdegree hdegree63
      leading hleading output hrun

/-- Outer refinement with both source callees instantiated: prime selection
executes the generated modular pipeline and Hensel lifting executes the
generated tree-build/quadratic-lift/extract/normalize pipeline.  The only
Hensel premise left here is its execution-readiness invariant over actual
intermediate results; no abstract Hensel function or result oracle remains. -/
theorem concreteSelectHensel___factor_squarefree_primitive_ZZ_raw_ir_refines
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hhenselInvariant : ∀ (selection : PrimeSelectionResult)
      (hp : Nat.Prime selection.prime.toNat) (aTarget : Int32),
      let candidate := provider.physical selection.prime hp
      @HenselLiftEntryReadiness candidate.dense
        StrictHensel.concreteDivmodTermination candidate.providers.mul
        f selection.factors aTarget)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegree62 : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hdegree63 : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir
      (concreteRecombineFactorZZRawOps
        (concreteSelectPrime engine provider initialRng)
        (concreteHenselLift engine provider))
      useLargePrime f = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly) := by
  exact concreteSelect___factor_squarefree_primitive_ZZ_raw_ir_refines
    engine provider initialRng (concreteHenselLift engine provider)
    useLargePrime f hinitialPrimeCorrect
    (by
      intro selection hp aTarget henselOutput hhenselRun
      exact concreteHenselLift_success engine provider f selection.factors
        selection.prime aTarget hp
        (hhenselInvariant selection hp aTarget) henselOutput hhenselRun)
    (by
      intro selection hp aTarget henselOutput hhenselRun
      exact concreteHenselLift_factorCountFits_of_success engine provider f
        selection.factors selection.prime aTarget hp henselOutput
        hhenselRun)
    hcanonical hprimitive hnonempty hdegree hdegree62 hdegree63 leading
    hleading output hrun

/-- The literal generated Zassenhaus attempt extracts the genuine legal
Hensel candidate.  Every intermediate result is obtained from the source
execution theorems above, including exact long division and quotient
primitive normalization. -/
theorem zassenhausAttempt_extracts_hensel_divisor_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorPrimitive : divisor.IsPrimitive)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hactiveFits : output.1.size ≤ 2 ^ 31)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ recoveredFactor recoveredQuotient,
      Generated.StrictRecombine.zassenhausAttempt f output.1 output.2 candidate =
        .ok (.extracted recoveredFactor recoveredQuotient) := by
  have hfront : f[0] = leading := by
    rw [Array.getElem?_eq_getElem hnonempty] at hleading
    exact Option.some.inj hleading
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    apply hselection.goodPrime.lc_nonzero
    rw [hzero]
    simp
  have hdivisorNe : divisor ≠ 0 := by
    intro hzero
    apply hdivisorModNonzero
    simp [hzero]
  have hquotientNe : quotient ≠ 0 := by
    intro hzero
    apply hsourceNe
    rw [hfactor, hzero]
    simp
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  rcases zassenhausLeadingPrune_accepts_hensel_candidate hentry hcanonical
      hnonempty hdegree leading hleading candidate hbound with
    ⟨leadingProduct, hleadingRun, hleadingAccept⟩
  rcases zassenhausConstantPrune_accepts_hensel_divisor_candidate hcount
      hfactors hleadingSemantic hselection hentry hcanonical hnonempty hdegree
      leading hleading divisor quotient hfactor hdivisorModNonzero
      hdivisorLeading candidate hlegal hassociated with
    ⟨constantProduct, hconstantRun, _hconstantRecovered, hconstantAccept⟩
  rcases zassenhausCandidate_executes_through_primitive hcount hfactors
      hleadingSemantic hselection hentry hcanonical hnonempty hdegree leading
      hleading divisor quotient hfactor hdivisorModNonzero hdivisorLeading
      candidate hlegal hactiveFits hassociated with
    ⟨candidate32, product, symmetric, content, recoveredFactor, hconvert,
      hproduct, hsymmetric, hsymmetricPoly, hsymmetricCanonical, hprimitive⟩
  have hrecoveredDvdDivisor :=
    primitiveRaw_factor_dvd_scaled_primitive_divisor symmetric recoveredFactor
      content divisor quotient hsymmetricPoly hsymmetricCanonical hprimitive
      hdivisorPrimitive hdivisorNe hquotientNe
  have hrecoveredDvdSource : SparsePolyZZ.toPoly recoveredFactor ∣
      SparsePolyZZ.toPoly f := by
    exact dvd_trans hrecoveredDvdDivisor
      (hfactor ▸ dvd_mul_right divisor quotient)
  have hrecoveredCanonical := StrictRecombine.primitiveRaw_canonical symmetric
    recoveredFactor content hsymmetricCanonical hprimitive
  have hrecoveredNonempty : 0 < recoveredFactor.size := by
    by_contra hnot
    have hzero : recoveredFactor.size = 0 := Nat.eq_zero_of_not_pos hnot
    have hempty : recoveredFactor = #[] := Array.size_eq_zero_iff.mp hzero
    rw [hempty] at hrecoveredDvdDivisor
    simp [SparsePolyZZ.toPoly] at hrecoveredDvdDivisor
    exact hdivisorNe hrecoveredDvdDivisor
  rcases StrictRecombine.exactDivmodRaw_complete_of_dvd f recoveredFactor
      hcanonical hrecoveredCanonical hrecoveredNonempty hrecoveredDvdSource with
    ⟨rawQuotient, hdivmod⟩
  have hrawQuotientCanonical :=
    StrictRecombine.exactDivmodRaw_quotient_canonical f recoveredFactor
      rawQuotient #[] hcanonical.2 hdivmod
  rcases StrictRecombine.primitiveRaw_complete rawQuotient
      hrawQuotientCanonical with
    ⟨quotientContent, recoveredQuotient, hquotientPrimitive⟩
  refine ⟨recoveredFactor, recoveredQuotient, ?_⟩
  unfold Generated.StrictRecombine.zassenhausAttempt
  rw [dif_pos hnonempty]
  dsimp only
  rw [hfront]
  simp only [hleadingRun]
  rw [if_neg hleadingAccept]
  simp only [hconstantRun]
  rw [if_neg hconstantAccept]
  simp only [hconvert, hproduct, hsymmetric, hprimitive, hdivmod]
  simp only [Array.isEmpty_empty, if_true, hquotientPrimitive]

/-- The fixed-size source scan, started from its literal iota candidate,
reaches an actual extraction once a genuine primitive divisor supplies the
legal Hensel subset above.  All live scan invariants are projected from the
selected-prime and generated Hensel certificates. -/
theorem zassenhausScan_extracts_hensel_divisor_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorPrimitive : divisor.IsPrimitive)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hcandidate : 0 < candidate.size)
    (hactiveFits : output.1.size ≤ 2 ^ 31)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ extractedFactor extractedQuotient extractedCandidate candidateSize,
      Generated.StrictRecombine.scanZassenhausCombinations
        (StrictRecombine.concreteZassenhausTermination.combinations
          output.1.size candidate.size)
        f output.1 output.2
        (Generated.StrictRecombine.initialCombination candidate.size)
        (StrictRecombine.concreteZassenhausTermination.initial_valid
          output.1.size candidate.size
          hlegal.count_le) =
          .ok (.extracted extractedFactor extractedQuotient extractedCandidate
            candidateSize) := by
  rcases hentry.outputModulus_eq_prime_pow with
    ⟨exponent, hexponent, houtputModulus⟩
  let modulus := selection.prime.toNat ^ exponent
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos exponent
  have hbase : 0 < selection.prime.toNat :=
    (Fact.out : Nat.Prime selection.prime.toNat).pos
  have hdivides : selection.prime.toNat ∣ modulus := by
    exact dvd_pow_self _ hexponent.ne'
  have hfront : f[0] = leading := by
    rw [Array.getElem?_eq_getElem hnonempty] at hleading
    exact Option.some.inj hleading
  have hleadingMod : (f[0].2 : ZMod selection.prime.toNat) ≠ 0 := by
    rw [hfront, hleadingSemantic leading hleading]
    exact hselection.goodPrime.lc_nonzero
  have hirreducible := selectionHenselFactors_mod_irreducible hcount
    hfactors hleadingSemantic hselection hentry
  rcases zassenhausAttempt_extracts_hensel_divisor_candidate hcount
      hfactors hleadingSemantic hselection hentry hcanonical hnonempty hdegree
      leading hleading divisor quotient hfactor hdivisorPrimitive
      hdivisorModNonzero hdivisorLeading candidate hlegal hactiveFits
      hassociated with ⟨factor, recoveredQuotient, hattempt⟩
  have hfits := hlegal.count_le
  rcases StrictRecombine.zassenhausFixedSizeScan_extracts_of_candidate f
      output.1 modulus selection.prime.toNat candidate hmodulus hbase hdivides
      hcanonical hnonempty hcandidate hfits hactiveFits hleadingMod
      hirreducible hlegal factor recoveredQuotient
      (by simpa [modulus, houtputModulus] using hattempt) with
    ⟨extractedFactor, extractedQuotient, extractedCandidate, candidateSize,
      hscan⟩
  refine ⟨extractedFactor, extractedQuotient, extractedCandidate,
    candidateSize, ?_⟩
  simpa [modulus, houtputModulus] using hscan

/-- If the actual generated fixed-size Zassenhaus scan exhausts, then the
occurrence-sensitive candidate supplied by any genuine integer divisor was
not omitted: that exact index array was executed and rejected.  This is the
bridge from divisor existence to the concrete combination enumerator; the
next completeness step rules out the rejection using symmetric recovery and
exact division. -/
theorem integer_divisor_candidate_rejected_of_scan_exhausted
    {termination : Generated.StrictHensel.DivmodTermination}
    {f fStar : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f)
    (hrun : StrictRecombine.FixedSizeScanExhausted fStar output.1 output.2
      (integer_divisor_mod_has_legal_hensel_candidate hcount hfactors
        hleadingSemantic hselection hentry g hg).choose.size) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        (((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) ∧
      Generated.StrictRecombine.zassenhausAttempt fStar output.1 output.2
        indices = .ok .rejected := by
  let witness := integer_divisor_mod_has_legal_hensel_candidate hcount
    hfactors hleadingSemantic hselection hentry g hg
  let indices := witness.choose
  have hspec := witness.choose_spec
  have hrejected := hrun.rejects indices hspec.1
  exact ⟨indices, hspec.1, hspec.2, hrejected⟩

/-- Pairwise coprimality supplied for the concrete adjusted input array is
transported through the actual Hensel leaf origins, lift extraction, and
final source normalization to the returned lifted-factor array. -/
theorem henselFactors_mod_pairwise_coprime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime p.toNat)]
    (hcount : 2 ≤ factors.size)
    (hpairwise : ∀ adjusted,
      StrictHensel.HenselAdjustFirstFactorCorrect f factors p adjusted →
      ∀ i j (hi : i < adjusted.size) (hj : j < adjusted.size), i < j →
        IsCoprime
          (SparsePolyZp.toPoly p.toNat adjusted[i])
          (SparsePolyZp.toPoly p.toNat adjusted[j]))
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output) :
    ∀ i j (hi : i < output.1.size) (hj : j < output.1.size), i < j →
      IsCoprime (StrictHensel.toPolyMod p.toNat output.1[i])
        (StrictHensel.toPolyMod p.toNat output.1[j]) := by
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjustSize : adjusted.size = factors.size := by
    cases hadjust with
    | adjusted leading first adjusted hsource hfirst hadjustedEq =>
        have hzero : 0 < factors.size := by omega
        simp [Array.set!, hzero]
  have hfull : StrictHensel.henselFactorRangeList adjusted factors.size 0 =
      adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have hlength : adjusted.size = extracted.size := by
    simpa using horigins.length_eq
  intro i j hi hj hij
  have hiExtracted : i < extracted.size := by rw [hnormalizeRel.1]; exact hi
  have hjExtracted : j < extracted.size := by rw [hnormalizeRel.1]; exact hj
  have hiAdjusted : i < adjusted.size := by omega
  have hjAdjusted : j < adjusted.size := by omega
  have hiOrigin := horigins.get (by simpa using hiAdjusted)
    (by simpa using hiExtracted)
  have hjOrigin := horigins.get (by simpa using hjAdjusted)
    (by simpa using hjExtracted)
  have hiOrigin' : StrictHensel.toPolyMod p.toNat extracted[i] =
      SparsePolyZp.toPoly p.toNat adjusted[i] := by
    simpa using hiOrigin
  have hjOrigin' : StrictHensel.toPolyMod p.toNat extracted[j] =
      SparsePolyZp.toPoly p.toNat adjusted[j] := by
    simpa using hjOrigin
  have hadjustedCoprime := hpairwise adjusted hadjust i j hiAdjusted
    hjAdjusted hij
  have hextractedCoprime : IsCoprime
      (StrictHensel.toPolyMod p.toNat extracted[i]!)
      (StrictHensel.toPolyMod p.toNat extracted[j]!) := by
    rw [getElem!_pos extracted i hiExtracted,
      getElem!_pos extracted j hjExtracted]
    rw [hiOrigin', hjOrigin']
    simpa [getElem!_pos adjusted i hiAdjusted,
      getElem!_pos adjusted j hjAdjusted] using hadjustedCoprime
  simpa [getElem!_pos output.1 i hi, getElem!_pos output.1 j hj] using
    hnormalizeRel.isCoprime hi hj hextractedCoprime

end StrictFactorZZ

end Refinement
