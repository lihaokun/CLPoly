/- Strict, fuel-free model of the GMP-backed uint64 prime enumeration API. -/
import CLPoly.Model
import Mathlib.Data.Nat.Prime.Infinite

set_option autoImplicit false

namespace Generated.StrictPrimeEnumeration

/-- Scan upward inside the uint64 range.  Reaching the end is exactly the
overflow exception of C++ `next_prime_64`. -/
def nextPrimeSearch (candidate : Nat) : RawExec UInt64 :=
  if hbound : candidate < UInt64.size then
    if hprime : Nat.Prime candidate then .ok candidate.toUInt64
    else nextPrimeSearch (candidate + 1)
  else .error .arithmeticOverflow
termination_by UInt64.size - candidate
decreasing_by omega

/-- Strict model of C++ `next_prime_64`, starting strictly after `p`. -/
def next_prime_64_raw (p : UInt64) : RawExec UInt64 :=
  nextPrimeSearch (p.toNat + 1)

/-- Scan downward.  The measure is the candidate itself, so no fuel or
partial recursion is involved. -/
def prevPrimeSearch (candidate : Nat) : RawExec UInt64 :=
  if hprime : Nat.Prime candidate then .ok candidate.toUInt64
  else if candidate = 0 then .error .arithmeticDomain
  else prevPrimeSearch (candidate - 1)
termination_by candidate
decreasing_by omega

/-- Strict model of C++ `prev_prime_64`, including its `n ≤ 2` exception. -/
def prev_prime_64_raw (p : UInt64) : RawExec UInt64 :=
  if p.toNat ≤ 2 then .error .arithmeticDomain
  else prevPrimeSearch (p.toNat - 1)

theorem nextPrimeSearch_success (candidate : Nat) (result : UInt64)
    (hrun : nextPrimeSearch candidate = .ok result) :
    Nat.Prime result.toNat ∧ candidate ≤ result.toNat := by
  induction hmeasure : UInt64.size - candidate using Nat.strong_induction_on
      generalizing candidate with
  | h measure ih =>
      rw [nextPrimeSearch] at hrun
      split at hrun
      next hbound =>
        split at hrun
        next hprime =>
          have hc : candidate.toUInt64.toNat = candidate := by
            simp [UInt64.toNat_ofNat, Nat.mod_eq_of_lt hbound]
          have := Except.ok.inj hrun
          subst result
          simpa [hc] using And.intro hprime (Nat.le_refl candidate)
        next hnotPrime =>
          have hrec := ih (UInt64.size - (candidate + 1)) (by omega)
            (candidate + 1) hrun rfl
          exact ⟨hrec.1, by omega⟩
      next hbound => contradiction

theorem next_prime_64_raw_success (p result : UInt64)
    (hrun : next_prime_64_raw p = .ok result) :
    Nat.Prime result.toNat ∧ p.toNat < result.toNat := by
  have h := nextPrimeSearch_success (p.toNat + 1) result hrun
  exact ⟨h.1, by omega⟩

theorem prevPrimeSearch_success (candidate : Nat) (result : UInt64)
    (hcandidate : candidate < UInt64.size)
    (hrun : prevPrimeSearch candidate = .ok result) :
    Nat.Prime result.toNat ∧ result.toNat ≤ candidate := by
  induction candidate using Nat.strong_induction_on generalizing result with
  | h candidate ih =>
      rw [prevPrimeSearch] at hrun
      split at hrun
      next hprime =>
        have hc : candidate.toUInt64.toNat = candidate := by
          simp [UInt64.toNat_ofNat, Nat.mod_eq_of_lt hcandidate]
        have := Except.ok.inj hrun
        subst result
        simpa [hc] using And.intro hprime (Nat.le_refl candidate)
      next hnotPrime =>
        split at hrun
        next hzero => contradiction
        next hnonzero =>
          have hrec := ih (candidate - 1) (by omega) result (by omega) hrun
          exact hrec.imp_right (Nat.le_trans · (Nat.sub_le candidate 1))

theorem prev_prime_64_raw_success (p result : UInt64)
    (hrun : prev_prime_64_raw p = .ok result) :
    Nat.Prime result.toNat ∧ result.toNat < p.toNat := by
  rw [prev_prime_64_raw] at hrun
  split at hrun
  next hsmall => contradiction
  next hlarge =>
    have hbound : p.toNat - 1 < UInt64.size :=
      lt_of_le_of_lt (Nat.sub_le _ _) p.toNat_lt_size
    have h := prevPrimeSearch_success (p.toNat - 1) result hbound hrun
    exact ⟨h.1, lt_of_le_of_lt h.2 (by omega)⟩

/-- Concrete C++ branch selecting ascending or descending GMP enumeration. -/
def nextPrimeRaw (useLargePrime : Bool) (p : UInt64) : RawExec UInt64 :=
  if useLargePrime then prev_prime_64_raw p else next_prime_64_raw p

def rank (useLargePrime : Bool) (p : UInt64) : Nat :=
  if useLargePrime = true then p.toNat else UInt64.size - p.toNat

theorem nextPrimeRaw_decreases (useLargePrime : Bool) (p p' : UInt64)
    (hrun : nextPrimeRaw useLargePrime p = .ok p') :
    rank useLargePrime p' < rank useLargePrime p := by
  cases useLargePrime <;> simp [nextPrimeRaw, rank] at hrun ⊢
  · have hlt := (next_prime_64_raw_success p p' hrun).2
    have hp' := p'.toNat_lt_size
    omega
  · exact (prev_prime_64_raw_success p p' hrun).2

theorem nextPrimeRaw_prime (useLargePrime : Bool) (p p' : UInt64)
    (hrun : nextPrimeRaw useLargePrime p = .ok p') : Nat.Prime p'.toNat := by
  cases useLargePrime <;> simp only [nextPrimeRaw, if_false, if_true] at hrun
  · exact (next_prime_64_raw_success p p' hrun).1
  · exact (prev_prime_64_raw_success p p' hrun).1

end Generated.StrictPrimeEnumeration
