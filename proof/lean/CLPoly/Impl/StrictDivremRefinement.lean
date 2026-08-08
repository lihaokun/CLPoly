import CLPoly.Impl.RawPolynomialRep
import CLPoly.Impl.StrictWordArithmetic

namespace CLPoly.Impl.StrictDivremRefinement

open Generated.StrictDivrem
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictWordArithmetic

/-- Complete pointwise observation of an allocated C++ `word3` slice.  The
array is proof data only: every entry is required to be the result of the
corresponding failure-aware raw read. -/
def Word3SliceRep (heap : RawHeap) (ptr : RawPtr Word3) (length : Nat)
    (values : Array Word3) : Prop :=
  values.size = length ∧
    ∀ i (hi : i < values.size),
      heap.readWord3 ptr i = .ok values[i]

theorem word3SliceRep_exists_unique (heap : RawHeap) (ptr : RawPtr Word3)
    (length : Nat) (hvalid : heap.ValidWord3Slice ptr length) :
    ∃ values : Array Word3,
      Word3SliceRep heap ptr length values ∧
      ∀ other : Array Word3,
        Word3SliceRep heap ptr length other → other = values := by
  classical
  let observed : Fin length → Word3 := fun i =>
    Classical.choose (heap.readWord3_of_valid ptr length i hvalid i.isLt)
  let values : Array Word3 := Array.ofFn observed
  have hread (i : Nat) (hi : i < length) (hiValues : i < values.size) :
      heap.readWord3 ptr i = .ok values[i] := by
    have hchosen := Classical.choose_spec
      (heap.readWord3_of_valid ptr length i hvalid hi)
    simpa [values, observed] using hchosen
  refine ⟨values, ⟨by simp [values], ?_⟩, ?_⟩
  · intro i hi
    exact hread i (by simpa [values] using hi) hi
  intro other hother
  apply Array.ext (hother.1.trans (by simp [values]))
  intro i hiOther hiValues
  exact Except.ok.inj ((hother.2 i hiOther).symm.trans
    (hread i (by simpa [values] using hiValues) hiValues))

theorem word3SliceRep_eq (heap : RawHeap) (ptr : RawPtr Word3)
    (length : Nat) (left right : Array Word3)
    (hleft : Word3SliceRep heap ptr length left)
    (hright : Word3SliceRep heap ptr length right) :
    left = right := by
  apply Array.ext (hleft.1.trans hright.1.symm)
  intro i hiLeft hiRight
  exact Except.ok.inj ((hleft.2 i hiLeft).symm.trans (hright.2 i hiRight))

/-- Polynomial observation of a W3 accumulator array modulo the source
prime. -/
noncomputable def word3ArrayPoly (p : Nat) (values : Array Word3) :
    Polynomial (ZMod p) :=
  ∑ i : Fin values.size,
    Polynomial.monomial i.val (word3Value values[i] : ZMod p)

theorem coeff_word3ArrayPoly (p : Nat) (values : Array Word3) (degree : Nat) :
    (word3ArrayPoly p values).coeff degree =
      if h : degree < values.size then
        (word3Value values[degree] : ZMod p)
      else 0 := by
  classical
  unfold word3ArrayPoly
  rw [Polynomial.finset_sum_coeff]
  by_cases hdegree : degree < values.size
  · rw [dif_pos hdegree, Finset.sum_eq_single ⟨degree, hdegree⟩]
    · simp
    · intro index _ hne
      have hval : index.val ≠ degree := by
        intro heq
        apply hne
        exact Fin.ext heq
      simp [Polynomial.coeff_monomial, hval]
    · simp
  · rw [dif_neg hdegree]
    apply Finset.sum_eq_zero
    intro index _
    have hval : index.val ≠ degree := by
      intro heq
      apply hdegree
      simpa [heq] using index.isLt
    simp [Polynomial.coeff_monomial, hval]

def SameU64Prefix (before after : RawHeap) (ptr : RawPtr UInt64)
    (length : Nat) : Prop :=
  ∀ k value, k < length → before.readU64 ptr k = .ok value →
    after.readU64 ptr k = .ok value

/-- A valid raw coefficient slice has a safe observation of exactly the C++
declared length. -/
theorem readU64s_ok (heap : RawHeap) (ptr : RawPtr UInt64) (count : Nat)
    (hvalid : heap.ValidU64Slice ptr count) :
    ∃ coeffs, heap.readU64s ptr count = .ok coeffs ∧ coeffs.size = count := by
  cases count with
  | zero => exact ⟨#[], rfl, rfl⟩
  | succ n =>
    simp only [RawHeap.readU64s]
    rcases heap.readU64_of_valid ptr (n + 1) 0 hvalid (by omega) with
      ⟨head, hhead⟩
    simp only [hhead]
    have htail := heap.validU64Slice_add ptr (n + 1) 1 n hvalid (by omega)
    rcases readU64s_ok heap (RawPtr.add ptr 1) n htail with
      ⟨tail, hreadTail, hsizeTail⟩
    simp only [hreadTail]
    refine ⟨#[head] ++ tail, rfl, ?_⟩
    simp [hsizeTail, Nat.add_comm]

theorem readU64s_get (heap : RawHeap) (ptr : RawPtr UInt64) (count : Nat)
    (coeffs : Array UInt64) (hread : heap.readU64s ptr count = .ok coeffs)
    (hsize : coeffs.size = count) (index : Nat) (hindex : index < count) :
    heap.readU64 ptr index = .ok coeffs[index] := by
  cases count with
  | zero => omega
  | succ n =>
    simp only [RawHeap.readU64s] at hread
    cases hhead : heap.readU64 ptr 0 with
    | error fault => simp [hhead] at hread
    | ok head =>
      cases htail : heap.readU64s (RawPtr.add ptr 1) n with
      | error fault => simp [hhead, htail] at hread
      | ok tail =>
        simp only [hhead, htail] at hread
        have hcoeffs : coeffs = #[head] ++ tail := (Except.ok.inj hread).symm
        subst coeffs
        have htailSizeEq : 1 + tail.size = n + 1 := by simpa using hsize
        have htailSize : tail.size = n := by omega
        cases index with
        | zero => simpa using hhead
        | succ j =>
          have hj : j < n := by omega
          have ih := readU64s_get heap (RawPtr.add ptr 1) n tail htail
            htailSize j hj
          rw [RawHeap.readU64_add] at ih
          simpa [Nat.add_comm] using ih

theorem sliceRep_exists_unique (heap : RawHeap) (ptr : RawPtr UInt64)
    (length : Nat) (hvalid : heap.ValidU64Slice ptr length) :
    ∃ coeffs : Array UInt64,
      (heap.SliceRep ptr length coeffs ∧ coeffs.size = length) ∧
      ∀ other : Array UInt64,
        heap.SliceRep ptr length other ∧ other.size = length → other = coeffs := by
  rcases readU64s_ok heap ptr length hvalid with ⟨coeffs, hread, hsize⟩
  refine ⟨coeffs, ⟨hread, hsize⟩, ?_⟩
  intro other hother
  exact (Except.ok.inj (hread.symm.trans hother.1)).symm

theorem slicePolyRep_exists_unique (heap : RawHeap) (ptr : RawPtr UInt64)
    (length p : Nat) (hvalid : heap.ValidU64Slice ptr length) :
    ∃ poly : Polynomial (ZMod p),
      SlicePolyRep heap ptr length p poly ∧
      ∀ other : Polynomial (ZMod p),
        SlicePolyRep heap ptr length p other → other = poly := by
  rcases sliceRep_exists_unique heap ptr length hvalid with
    ⟨coeffs, hcoeffs, hunique⟩
  let poly := coeffArrayPoly p coeffs
  refine ⟨poly, ⟨coeffs, hcoeffs.1, hcoeffs.2, rfl⟩, ?_⟩
  intro other hother
  rcases hother with ⟨otherCoeffs, hreadOther, hsizeOther, hpolyOther⟩
  have heq : otherCoeffs = coeffs :=
    hunique otherCoeffs ⟨hreadOther, hsizeOther⟩
  subst otherCoeffs
  exact hpolyOther

theorem slicePolyRep_coeff (heap : RawHeap) (ptr : RawPtr UInt64)
    (length p : Nat) (poly : Polynomial (ZMod p))
    (hrep : SlicePolyRep heap ptr length p poly)
    (degree : Nat) (hdegree : degree < length) :
    ∃ value : UInt64,
      heap.readU64 ptr degree = .ok value ∧
      poly.coeff degree = (value.toNat : ZMod p) := by
  rcases hrep with ⟨coeffs, hread, hsize, hpoly⟩
  have hraw := readU64s_get heap ptr length coeffs hread hsize degree hdegree
  refine ⟨coeffs[degree], hraw, ?_⟩
  rw [hpoly, coeff_coeffArrayPoly, dif_pos]

theorem slicePolyRep_coeff_zero_of_length_le (heap : RawHeap)
    (ptr : RawPtr UInt64) (length p : Nat)
    (poly : Polynomial (ZMod p))
    (hrep : SlicePolyRep heap ptr length p poly)
    (degree : Nat) (hdegree : length ≤ degree) :
    poly.coeff degree = 0 := by
  rcases hrep with ⟨coeffs, _, hsize, rfl⟩
  rw [coeff_coeffArrayPoly, dif_neg]
  simpa [hsize] using hdegree

theorem slicePolyRep_zero_length (heap : RawHeap) (ptr : RawPtr UInt64)
    (p : Nat) (poly : Polynomial (ZMod p))
    (hrep : SlicePolyRep heap ptr 0 p poly) :
    poly = 0 := by
  ext degree
  rw [slicePolyRep_coeff_zero_of_length_le heap ptr 0 p poly hrep degree
    (Nat.zero_le _), Polynomial.coeff_zero]

/-- Extending a raw coefficient prefix by its next actual C++ cell extends
the represented polynomial by the corresponding monomial. -/
theorem slicePolyRep_succ_eq_add_monomial (heap : RawHeap)
    (ptr : RawPtr UInt64) (length p : Nat) (value : UInt64)
    (prefixPoly full : Polynomial (ZMod p))
    (hprefix : SlicePolyRep heap ptr length p prefixPoly)
    (hfull : SlicePolyRep heap ptr (length + 1) p full)
    (hread : heap.readU64 ptr length = .ok value) :
    full = prefixPoly + Polynomial.monomial length (value.toNat : ZMod p) := by
  ext degree
  rw [Polynomial.coeff_add, Polynomial.coeff_monomial]
  by_cases hdegree : degree < length
  · rcases slicePolyRep_coeff heap ptr length p prefixPoly hprefix degree
        hdegree with ⟨prefixValue, hreadPrefix, hcoeffPrefix⟩
    rcases slicePolyRep_coeff heap ptr (length + 1) p full hfull degree
        (by omega) with ⟨fullValue, hreadFull, hcoeffFull⟩
    have hvalues : fullValue = prefixValue :=
      Except.ok.inj (hreadFull.symm.trans hreadPrefix)
    rw [hcoeffFull, hcoeffPrefix, hvalues]
    rw [if_neg (show length ≠ degree by omega), add_zero]
  · by_cases heq : degree = length
    · subst degree
      rcases slicePolyRep_coeff heap ptr (length + 1) p full hfull length
          (by omega) with ⟨fullValue, hreadFull, hcoeffFull⟩
      have hvalue : fullValue = value :=
        Except.ok.inj (hreadFull.symm.trans hread)
      rw [hcoeffFull, hvalue]
      have hzero := slicePolyRep_coeff_zero_of_length_le heap ptr length p
        prefixPoly hprefix length (Nat.le_refl _)
      simp [hzero]
    · have hdegreeHigh : length + 1 ≤ degree := by omega
      rw [slicePolyRep_coeff_zero_of_length_le heap ptr (length + 1) p
          full hfull degree hdegreeHigh,
        slicePolyRep_coeff_zero_of_length_le heap ptr length p prefixPoly
          hprefix degree (by omega)]
      rw [if_neg (show length ≠ degree by omega), add_zero]

theorem slicePolyRep_extend_exists (heap : RawHeap) (ptr : RawPtr UInt64)
    (length p : Nat) (value : UInt64)
    (prefixPoly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr (length + 1))
    (hprefix : SlicePolyRep heap ptr length p prefixPoly)
    (hread : heap.readU64 ptr length = .ok value) :
    ∃ full : Polynomial (ZMod p),
      SlicePolyRep heap ptr (length + 1) p full ∧
      full = prefixPoly +
        Polynomial.monomial length (value.toNat : ZMod p) := by
  rcases slicePolyRep_exists_unique heap ptr (length + 1) p hvalid with
    ⟨full, hfull, _⟩
  exact ⟨full, hfull, slicePolyRep_succ_eq_add_monomial heap ptr length p
    value prefixPoly full hprefix hfull hread⟩

/-- Pointwise preservation of a raw UInt64 slice transports its L2
polynomial representation without re-running any L2 algorithm. -/
theorem slicePolyRep_of_same_prefix (before after : RawHeap)
    (ptr : RawPtr UInt64) (length p : Nat)
    (poly : Polynomial (ZMod p))
    (hvalidBefore : before.ValidU64Slice ptr length)
    (hvalidAfter : after.ValidU64Slice ptr length)
    (hsame : SameU64Prefix before after ptr length)
    (hrep : SlicePolyRep before ptr length p poly) :
    SlicePolyRep after ptr length p poly := by
  rcases slicePolyRep_exists_unique after ptr length p hvalidAfter with
    ⟨afterPoly, hafterRep, _⟩
  have heq : afterPoly = poly := by
    ext degree
    by_cases hdegree : degree < length
    · rcases slicePolyRep_coeff before ptr length p poly hrep degree
          hdegree with ⟨beforeValue, hreadBefore, hcoeffBefore⟩
      rcases slicePolyRep_coeff after ptr length p afterPoly hafterRep degree
          hdegree with ⟨afterValue, hreadAfter, hcoeffAfter⟩
      have hreadPreserved := hsame degree beforeValue hdegree hreadBefore
      have hvalue : afterValue = beforeValue :=
        Except.ok.inj (hreadAfter.symm.trans hreadPreserved)
      rw [hcoeffAfter, hcoeffBefore, hvalue]
    · rw [slicePolyRep_coeff_zero_of_length_le after ptr length p
          afterPoly hafterRep degree (by omega),
        slicePolyRep_coeff_zero_of_length_le before ptr length p poly hrep
          degree (by omega)]
  rw [heq] at hafterRep
  exact hafterRep

/-- `_poly_normalise` only reads within its declared prefix.  Structural
recursion on `len` proves both successful execution and the returned prefix
bound; structural recursion produces no alternate execution result. -/
theorem normaliseU64_ok (heap : RawHeap) (ptr : RawPtr UInt64) (len : Nat)
    (hvalid : heap.ValidU64Slice ptr len) :
    ∃ result, heap.normaliseU64 ptr len = .ok result ∧ result ≤ len := by
  cases len with
  | zero => exact ⟨0, rfl, Nat.le_refl 0⟩
  | succ n =>
    simp only [RawHeap.normaliseU64]
    rcases heap.readU64_of_valid ptr (n + 1) n hvalid (by omega) with
      ⟨value, hread⟩
    simp only [hread]
    split
    next hzero =>
      have hprefix := heap.validU64Slice_mono ptr (n + 1) n hvalid (by omega)
      rcases normaliseU64_ok heap ptr n hprefix with ⟨result, hresult, hle⟩
      exact ⟨result, hresult, by omega⟩
    next hnonzero => exact ⟨n + 1, rfl, Nat.le_refl _⟩

/-- Exact content specification of generated `_poly_normalise`: the returned
prefix discards only zero trailing limbs, and a nonempty returned prefix ends
in a nonzero limb. -/
theorem normaliseU64_spec (heap : RawHeap) (ptr : RawPtr UInt64) (len : Nat)
    (hvalid : heap.ValidU64Slice ptr len) :
    ∃ result, heap.normaliseU64 ptr len = .ok result ∧ result ≤ len ∧
      (∀ j, result ≤ j → j < len → heap.readU64 ptr j = .ok 0) ∧
      (result = 0 ∨ ∃ value : UInt64,
        heap.readU64 ptr (result - 1) = .ok value ∧ value ≠ 0) := by
  cases len with
  | zero =>
      refine ⟨0, rfl, Nat.le_refl 0, ?_, Or.inl rfl⟩
      intro j _ hj
      omega
  | succ n =>
      simp only [RawHeap.normaliseU64]
      rcases heap.readU64_of_valid ptr (n + 1) n hvalid (by omega) with
        ⟨value, hread⟩
      simp only [hread]
      split
      next hzero =>
        have hvalue : value = 0 := by simpa using hzero
        have hprefix := heap.validU64Slice_mono ptr (n + 1) n hvalid (by omega)
        rcases normaliseU64_spec heap ptr n hprefix with
          ⟨result, hresult, hle, hzeros, hlast⟩
        refine ⟨result, hresult, by omega, ?_, hlast⟩
        intro j hjlow hjhigh
        by_cases hjn : j < n
        · exact hzeros j hjlow hjn
        · have hjeq : j = n := by omega
          subst j
          simpa [hvalue] using hread
      next hnonzero =>
        have hvalue : value ≠ 0 := by simpa using hnonzero
        refine ⟨n + 1, rfl, Nat.le_refl _, ?_, Or.inr ?_⟩
        · intro j hjlow hjhigh
          omega
        · exact ⟨value, by simpa using hread, hvalue⟩

/-- L2 consequence of raw normalization: every coefficient at or above the
returned prefix length is zero in the represented `Polynomial (ZMod p)`. -/
theorem normaliseU64_poly_coeff_zero (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hnorm : heap.normaliseU64 ptr len = .ok result) :
    ∀ degree, result ≤ degree → poly.coeff degree = 0 := by
  rcases normaliseU64_spec heap ptr len hvalid with
    ⟨result', hnorm', hle, hzeros, hlast⟩
  have hresult : result' = result := Except.ok.inj (hnorm'.symm.trans hnorm)
  subst result'
  intro degree hdegree
  by_cases hin : degree < len
  · rcases slicePolyRep_coeff heap ptr len p poly hrep degree hin with
      ⟨value, hread, hcoeff⟩
    have hzeroRead := hzeros degree hdegree hin
    have hvalue : value = 0 := Except.ok.inj (hread.symm.trans hzeroRead)
    subst value
    simpa using hcoeff
  · rcases hrep with ⟨coeffs, hread, hsize, rfl⟩
    rw [coeff_coeffArrayPoly, dif_neg]
    simpa [hsize] using hin

/-- Shrinking a raw coefficient slice to the exact length returned by the
generated normalization routine preserves its represented polynomial. -/
theorem slicePolyRep_of_normaliseU64 (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hnorm : heap.normaliseU64 ptr len = .ok result) :
    SlicePolyRep heap ptr result p poly := by
  rcases normaliseU64_ok heap ptr len hvalid with
    ⟨result', hnorm', hresultLe⟩
  have hresult : result' = result := Except.ok.inj (hnorm'.symm.trans hnorm)
  subst result'
  have hvalidResult := heap.validU64Slice_mono ptr len result hvalid hresultLe
  rcases slicePolyRep_exists_unique heap ptr result p hvalidResult with
    ⟨prefixPoly, hprefix, _⟩
  have heq : prefixPoly = poly := by
    ext degree
    by_cases hdegree : degree < result
    · rcases slicePolyRep_coeff heap ptr result p prefixPoly hprefix degree
          hdegree with ⟨prefixValue, hreadPrefix, hcoeffPrefix⟩
      rcases slicePolyRep_coeff heap ptr len p poly hrep degree
          (by omega) with ⟨fullValue, hreadFull, hcoeffFull⟩
      have hvalue : prefixValue = fullValue :=
        Except.ok.inj (hreadPrefix.symm.trans hreadFull)
      rw [hcoeffPrefix, hcoeffFull, hvalue]
    · have hprefixZero := slicePolyRep_coeff_zero_of_length_le heap ptr
        result p prefixPoly hprefix degree (by omega)
      have hfullZero := normaliseU64_poly_coeff_zero heap ptr len p result
        poly hvalid hrep hnorm degree (by omega)
      rw [hprefixZero, hfullZero]
  rw [← heq]
  exact hprefix

theorem normaliseU64_poly_natDegree_le (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hnorm : heap.normaliseU64 ptr len = .ok result) :
    poly.natDegree ≤ result - 1 := by
  apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
  intro degree hdegree
  apply normaliseU64_poly_coeff_zero heap ptr len p result poly
    hvalid hrep hnorm degree
  omega

theorem normaliseU64_poly_eq_zero (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hnorm : heap.normaliseU64 ptr len = .ok 0) :
    poly = 0 := by
  ext degree
  have hcoeff := normaliseU64_poly_coeff_zero heap ptr len p 0 poly
    hvalid hrep hnorm degree (Nat.zero_le degree)
  simpa using hcoeff

/-- If the raw coefficients are canonical residues, the nonempty prefix
returned by C++ normalization ends in a genuinely nonzero L2 coefficient. -/
theorem normaliseU64_poly_last_coeff_ne_zero (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hreduced : ∀ i value, i < len → heap.readU64 ptr i = .ok value →
      value.toNat < p)
    (hnorm : heap.normaliseU64 ptr len = .ok result)
    (hresult : result ≠ 0) :
    poly.coeff (result - 1) ≠ 0 := by
  rcases normaliseU64_spec heap ptr len hvalid with
    ⟨result', hnorm', hle, hzeros, hlast⟩
  have heq : result' = result := Except.ok.inj (hnorm'.symm.trans hnorm)
  subst result'
  rcases hlast with hzero | ⟨value, hread, hvalue⟩
  · exact (hresult hzero).elim
  · have hindex : result - 1 < len := by omega
    rcases slicePolyRep_coeff heap ptr len p poly hrep (result - 1) hindex with
      ⟨value', hread', hcoeff⟩
    have hvalueEq : value' = value := Except.ok.inj (hread'.symm.trans hread)
    subst value'
    rw [hcoeff]
    intro hcast
    have hdvd : p ∣ value.toNat :=
      (ZMod.natCast_eq_zero_iff value.toNat p).mp hcast
    have hvalueNat : value.toNat ≠ 0 := by
      intro hzero
      apply hvalue
      apply UInt64.toNat_inj.mp
      simpa using hzero
    exact (Nat.not_dvd_of_pos_of_lt
      (Nat.pos_of_ne_zero hvalueNat)
      (hreduced (result - 1) value hindex hread)) hdvd

theorem normaliseU64_poly_natDegree_eq (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hreduced : ∀ i value, i < len → heap.readU64 ptr i = .ok value →
      value.toNat < p)
    (hnorm : heap.normaliseU64 ptr len = .ok result)
    (hresult : result ≠ 0) :
    poly.natDegree = result - 1 := by
  apply Polynomial.natDegree_eq_of_le_of_coeff_ne_zero
    (normaliseU64_poly_natDegree_le heap ptr len p result poly
      hvalid hrep hnorm)
  exact normaliseU64_poly_last_coeff_ne_zero heap ptr len p result poly
    hvalid hrep hreduced hnorm hresult

/-- Any polynomial represented by the C++ normalization buffer is either zero
or has degree strictly below the physical buffer length. -/
theorem normaliseU64_poly_degree_lt_length (heap : RawHeap)
    (ptr : RawPtr UInt64) (len p result : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr len)
    (hrep : SlicePolyRep heap ptr len p poly)
    (hnorm : heap.normaliseU64 ptr len = .ok result) :
    poly = 0 ∨ poly.natDegree < len := by
  by_cases hresult : result = 0
  · left
    subst result
    exact normaliseU64_poly_eq_zero heap ptr len p poly hvalid hrep hnorm
  · right
    rcases normaliseU64_ok heap ptr len hvalid with
      ⟨result', hnorm', hresultLe⟩
    have heq : result' = result := Except.ok.inj (hnorm'.symm.trans hnorm)
    subst result'
    have hdegree := normaliseU64_poly_natDegree_le heap ptr len p result
      poly hvalid hrep hnorm
    omega

/-- Raw-to-safe output bridge for a generated remainder buffer after its
source `_poly_normalise` call. -/
theorem generated_remainder_output_degree (heap : RawHeap)
    (R : RawPtr UInt64) (d p lenR : Nat)
    (hR : heap.ValidU64Slice R d)
    (hnorm : heap.normaliseU64 R d = .ok lenR) :
    ∃ remainder : Polynomial (ZMod p),
      SlicePolyRep heap R d p remainder ∧
      (remainder = 0 ∨ remainder.natDegree < d) := by
  rcases slicePolyRep_exists_unique heap R d p hR with
    ⟨remainder, hrep, hunique⟩
  exact ⟨remainder, hrep,
    normaliseU64_poly_degree_lt_length heap R d p lenR remainder
      hR hrep hnorm⟩

/-- Degree side of strict long division, connected to the normalized raw
divisor rather than merely its buffer length. -/
theorem generated_remainder_degree_lt_divisor (heap : RawHeap)
    (R B : RawPtr UInt64) (d p lenR : Nat)
    (remainder divisor : Polynomial (ZMod p))
    (hR : heap.ValidU64Slice R d)
    (hB : heap.ValidU64Slice B (d + 1))
    (hrepR : SlicePolyRep heap R d p remainder)
    (hrepB : SlicePolyRep heap B (d + 1) p divisor)
    (hreducedB : ∀ i value, i < d + 1 →
      heap.readU64 B i = .ok value → value.toNat < p)
    (hnormR : heap.normaliseU64 R d = .ok lenR)
    (hnormB : heap.normaliseU64 B (d + 1) = .ok (d + 1)) :
    remainder = 0 ∨ remainder.natDegree < divisor.natDegree := by
  have hdegB : divisor.natDegree = d := by
    simpa using normaliseU64_poly_natDegree_eq heap B (d + 1) p (d + 1)
      divisor hB hrepB hreducedB hnormB (by omega)
  rcases normaliseU64_poly_degree_lt_length heap R d p lenR remainder
    hR hrepR hnormR with hzero | hdegree
  · exact Or.inl hzero
  · exact Or.inr (by simpa [hdegB] using hdegree)

/-- Limb `memcpy` succeeds for valid source and destination slices.  The
returned heap has exactly the same allocation layout as the input heap, so
all caller slice invariants can be transported through the copy. -/
theorem copyU64_ok (heap : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  cases count with
  | zero =>
    refine ⟨heap, rfl, ?_⟩
    intro ptr length
    exact Iff.rfl
  | succ n =>
    simp only [RawHeap.copyU64]
    rcases heap.readU64_of_valid src (n + 1) 0 hSrc (by omega) with
      ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeU64_of_valid dst (n + 1) 0 value hDst (by omega) with
      ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst 0 value hwrite
    have hDstTail0 := heap.validU64Slice_add dst (n + 1) 1 n hDst (by omega)
    have hSrcTail0 := heap.validU64Slice_add src (n + 1) 1 n hSrc (by omega)
    have hDstTail1 : heap1.ValidU64Slice (RawPtr.add dst 1) n :=
      (hlayout1 (RawPtr.add dst 1) n).mp hDstTail0
    have hSrcTail1 : heap1.ValidU64Slice (RawPtr.add src 1) n :=
      (hlayout1 (RawPtr.add src 1) n).mp hSrcTail0
    rcases copyU64_ok heap1 (RawPtr.add dst 1) (RawPtr.add src 1) n
      hDstTail1 hSrcTail1 with ⟨heap2, hcopy, hlayout2⟩
    simp only [hcopy]
    refine ⟨heap2, rfl, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- A successful generated `memcpy` preserves any raw limb whose absolute
address is outside its destination range. -/
theorem copyU64_preserves_read (heap heap' : RawHeap)
    (dst src guard : RawPtr UInt64) (count readIndex : Nat)
    (old : UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ k, k < count →
      dst.region ≠ guard.region ∨
        dst.limbOffset + k ≠ guard.limbOffset + readIndex)
    (hcopy : heap.copyU64 dst src count = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  cases count with
  | zero =>
    simp only [RawHeap.copyU64] at hcopy
    have heq : heap' = heap := Except.ok.inj hcopy.symm
    simpa [heq] using hread
  | succ n =>
    simp only [RawHeap.copyU64] at hcopy
    rcases heap.readU64_of_valid src (n + 1) 0 hSrc (by omega) with
      ⟨value, hvalue⟩
    simp only [hvalue] at hcopy
    rcases heap.writeU64_of_valid dst (n + 1) 0 value hDst (by omega) with
      ⟨heap1, hwrite⟩
    simp only [hwrite] at hcopy
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 dst guard
      0 readIndex value old hwrite hread (houtside 0 (by omega))
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst 0 value hwrite
    have hDstTail0 := heap.validU64Slice_add dst (n + 1) 1 n hDst (by omega)
    have hSrcTail0 := heap.validU64Slice_add src (n + 1) 1 n hSrc (by omega)
    have hDstTail1 : heap1.ValidU64Slice (RawPtr.add dst 1) n :=
      (hlayout1 (RawPtr.add dst 1) n).mp hDstTail0
    have hSrcTail1 : heap1.ValidU64Slice (RawPtr.add src 1) n :=
      (hlayout1 (RawPtr.add src 1) n).mp hSrcTail0
    apply copyU64_preserves_read heap1 heap'
      (RawPtr.add dst 1) (RawPtr.add src 1) guard n readIndex old
      hDstTail1 hSrcTail1 hread1
    · intro k hk
      have hout := houtside (k + 1) (by omega)
      rcases hout with hregion | hoffset
      · exact Or.inl (by simpa [RawPtr.add] using hregion)
      · right
        dsimp [RawPtr.add]
        change dst.limbOffset + 1 + k ≠ guard.limbOffset + readIndex
        omega
    · exact hcopy
termination_by count

def CopyU64Contents (before after : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat) : Prop :=
  ∀ k value, k < count → before.readU64 src k = .ok value →
    after.readU64 dst k = .ok value

/-- Content-level semantics of the actual recursive `copyU64` used by the
short-division branch.  The C++ non-overlap precondition is explicit. -/
theorem copyU64_refines (heap : RawHeap) (dst src : RawPtr UInt64)
    (count : Nat)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hregions : dst.region ≠ src.region) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      CopyU64Contents heap heap' dst src count := by
  cases count with
  | zero =>
    exact ⟨heap, rfl, fun _ _ => Iff.rfl, by intro k _ hk; omega⟩
  | succ n =>
    simp only [RawHeap.copyU64]
    rcases heap.readU64_of_valid src (n + 1) 0 hSrc (by omega) with
      ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeU64_of_valid dst (n + 1) 0 value hDst (by omega) with
      ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst 0 value hwrite
    have hDstTail0 := heap.validU64Slice_add dst (n + 1) 1 n hDst (by omega)
    have hSrcTail0 := heap.validU64Slice_add src (n + 1) 1 n hSrc (by omega)
    have hDstTail1 : heap1.ValidU64Slice (RawPtr.add dst 1) n :=
      (hlayout1 (RawPtr.add dst 1) n).mp hDstTail0
    have hSrcTail1 : heap1.ValidU64Slice (RawPtr.add src 1) n :=
      (hlayout1 (RawPtr.add src 1) n).mp hSrcTail0
    rcases copyU64_refines heap1 (RawPtr.add dst 1) (RawPtr.add src 1) n
      hDstTail1 hSrcTail1 (by simpa [RawPtr.add] using hregions) with
      ⟨heap2, hcopy, hlayout2, hcontents⟩
    simp only [hcopy]
    refine ⟨heap2, rfl, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k old hk hold
      cases k with
      | zero =>
        have hnow := RawHeap.readU64_writeU64_same heap heap1 dst 0 value hwrite
        have hpres := copyU64_preserves_read heap1 heap2
          (RawPtr.add dst 1) (RawPtr.add src 1) dst n 0 value
          hDstTail1 hSrcTail1 hnow (by
            intro j hj
            right
            dsimp [RawPtr.add]
            change dst.limbOffset + 1 + j ≠ dst.limbOffset + 0
            omega) hcopy
        have holdEq : old = value := Except.ok.inj (hold.symm.trans hread)
        simpa [holdEq] using hpres
      | succ k =>
        have hk' : k < n := by omega
        have hold1 := RawHeap.readU64_writeU64_ne heap heap1 dst src
          0 (k + 1) value old hwrite hold (Or.inl hregions)
        have htail := hcontents k old hk'
        rw [RawHeap.readU64_add, RawHeap.readU64_add] at htail
        have hout := htail (by simpa [Nat.add_comm] using hold1)
        simpa [Nat.add_comm] using hout
termination_by count

theorem copyU64_sliceRep (heap : RawHeap) (dst src : RawPtr UInt64)
    (count : Nat) (coeffs : Array UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hregions : dst.region ≠ src.region)
    (hrep : heap.SliceRep src count coeffs)
    (hsize : coeffs.size = count) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      heap'.SliceRep dst count coeffs := by
  rcases copyU64_refines heap dst src count hDst hSrc hregions with
    ⟨heap', hcopy, hlayout, hcontents⟩
  have hDst' : heap'.ValidU64Slice dst count :=
    (hlayout dst count).mp hDst
  rcases readU64s_ok heap' dst count hDst' with
    ⟨other, hother, hotherSize⟩
  have heq : other = coeffs := by
    apply Array.ext (hotherSize.trans hsize.symm)
    intro i hiOther hiCoeffs
    have hi : i < count := by simpa [hsize] using hiCoeffs
    have hsrc := readU64s_get heap src count coeffs hrep hsize i hi
    have hdst := hcontents i coeffs[i] hi hsrc
    have hdstOther := readU64s_get heap' dst count other hother
      hotherSize i hi
    exact Except.ok.inj (hdstOther.symm.trans hdst)
  subst other
  exact ⟨heap', hcopy, hlayout, hother⟩

theorem copyU64_slicePolyRep (heap : RawHeap) (dst src : RawPtr UInt64)
    (count p : Nat) (poly : Polynomial (ZMod p))
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hregions : dst.region ≠ src.region)
    (hrep : SlicePolyRep heap src count p poly) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' dst count p poly := by
  rcases hrep with ⟨coeffs, hslice, hsize, hpoly⟩
  rcases copyU64_sliceRep heap dst src count coeffs hDst hSrc hregions
    hslice hsize with ⟨heap', hcopy, hlayout, hslice'⟩
  exact ⟨heap', hcopy, hlayout, coeffs, hslice', hsize, hpoly⟩

/-- Natural-language proof outline:

At an iteration with `i < lenA`, validity of `A[0..lenA)` gives a successful
read of `A[i]`; validity of the `lenA`-element W3 slice gives a successful
three-limb write at `W3[i]`.  `writeWord3_preserves_valid` shows that this
write changes no allocation sizes, so both slice invariants hold for the
recursive heap.  The recursive measure `lenA - (i+1)` is smaller than
`lenA-i`.  When `i ≥ lenA`, the source loop exits without a heap access. -/
theorem initW3Loop_ok (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenA i : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA) (hi : i ≤ lenA) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA ∧
      RawHeap.SameLayout heap heap' := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeWord3_of_valid W3 lenA i
      { lo := value, mid := 0, hi := 0 } hW3 hlt with ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite
        (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 i
      { lo := value, mid := 0, hi := 0 } hwrite
    rcases initW3Loop_ok heap1 A W3 lenA (i + 1) hA1 hW31 (by omega) with
      ⟨heap2, hloop, hA2, hW32, hlayout2⟩
    refine ⟨heap2, hloop, hA2, hW32, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    exact ⟨heap, rfl, hA, hW3, fun _ _ => Iff.rfl⟩
termination_by lenA - i
decreasing_by omega

/-- The already initialized W3 prefix is the exact zero-extended image of
the corresponding source coefficient prefix. -/
def InitW3Prefix (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (upto : Nat) : Prop :=
  ∀ j, j < upto → ∃ value : UInt64,
    heap.readU64 A j = .ok value ∧
    heap.readWord3 W3 j = .ok { lo := value, mid := 0, hi := 0 }

/-- Every raw coefficient in the declared C++ slice is a canonical residue.
This is a representation invariant of `dense_upoly_zp`, not an L2 polynomial
operation. -/
def CanonicalU64Prefix (heap : RawHeap) (ptr : RawPtr UInt64)
    (length : Nat) (p : UInt64) : Prop :=
  ∀ k value, k < length → heap.readU64 ptr k = .ok value →
    value.toNat < p.toNat

/-- Uniform capacity invariant for lazy C++ `word3` accumulation.  `count` is
an upper bound on the number of quotient products already added to any cell.
-/
def Word3AccumulationBudget (heap : RawHeap) (W3 : RawPtr Word3)
    (length : Nat) (p : UInt64) (count : Nat) : Prop :=
  ∀ k accum, k < length → heap.readWord3 W3 k = .ok accum →
    word3Value accum ≤ (p.toNat - 1) + count * (p.toNat - 1) ^ 2

/-- Inner-loop form of the capacity invariant.  Cells in
`[offset, offset + processed)` have received the current quotient product;
all other cells retain the previous outer-loop allowance. -/
def Word3StagedBudget (heap : RawHeap) (W3 : RawPtr Word3)
    (length : Nat) (p : UInt64) (count offset processed : Nat) : Prop :=
  ∀ k accum, k < length → heap.readWord3 W3 k = .ok accum →
    word3Value accum ≤ (p.toNat - 1) +
      (count + if offset ≤ k ∧ k < offset + processed then 1 else 0) *
        (p.toNat - 1) ^ 2

theorem accumulationBudget_to_staged_zero (heap : RawHeap)
    (W3 : RawPtr Word3) (length : Nat) (p : UInt64)
    (count offset : Nat)
    (hbudget : Word3AccumulationBudget heap W3 length p count) :
    Word3StagedBudget heap W3 length p count offset 0 := by
  intro k accum hk hread
  have hnot : ¬(offset ≤ k ∧ k < offset) := by omega
  simpa [hnot] using hbudget k accum hk hread

theorem stagedBudget_to_next_budget (heap : RawHeap)
    (W3 : RawPtr Word3) (length : Nat) (p : UInt64)
    (count offset processed : Nat)
    (hstage : Word3StagedBudget heap W3 length p count offset processed) :
    Word3AccumulationBudget heap W3 length p (count + 1) := by
  intro k accum hk hread
  have h := hstage k accum hk hread
  by_cases hin : offset ≤ k ∧ k < offset + processed
  · simpa [Word3StagedBudget, hin] using h
  · simp only [Word3StagedBudget, hin, ↓reduceIte] at h
    calc
      word3Value accum ≤
          (p.toNat - 1) + count * (p.toNat - 1) ^ 2 := by simpa using h
      _ ≤ (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by
        apply Nat.add_le_add_left
        apply Nat.mul_le_mul_right
        omega

theorem accumulationBudget_mono (heap : RawHeap) (W3 : RawPtr Word3)
    (length : Nat) (p : UInt64) (small large : Nat)
    (hbudget : Word3AccumulationBudget heap W3 length p small)
    (hle : small ≤ large) :
    Word3AccumulationBudget heap W3 length p large := by
  intro k accum hk hread
  calc
    word3Value accum ≤
        (p.toNat - 1) + small * (p.toNat - 1) ^ 2 :=
      hbudget k accum hk hread
    _ ≤ (p.toNat - 1) + large * (p.toNat - 1) ^ 2 := by
      exact Nat.add_le_add_left (Nat.mul_le_mul_right _ hle) _

/-- Zero-extension performed by the generated initialization loop establishes
the zero-product instance of the accumulation budget. -/
theorem initW3Prefix_budget (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (length : Nat) (p : UInt64)
    (hcanonical : CanonicalU64Prefix heap A length p)
    (hprefix : InitW3Prefix heap A W3 length) :
    Word3AccumulationBudget heap W3 length p 0 := by
  intro k accum hk hreadW
  rcases hprefix k hk with ⟨value, hreadA, hreadInit⟩
  have haccum : accum = { lo := value, mid := 0, hi := 0 } :=
    Except.ok.inj (hreadW.symm.trans hreadInit)
  subst accum
  have hvalue := hcanonical k value hk hreadA
  simp only [word3Value, UInt64.toNat_zero, Nat.mul_zero, Nat.add_zero,
    zero_mul]
  omega

/-- The generated W3 initialization is already the input polynomial, when
observed cell-by-cell through the raw heap.  This is the algebraic base case
for the later quotient-times-divisor accumulation proof. -/
theorem initW3Prefix_word3ArrayPoly (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (length p : Nat) (values : Array Word3)
    (dividend : Polynomial (ZMod p))
    (hprefix : InitW3Prefix heap A W3 length)
    (hvalues : Word3SliceRep heap W3 length values)
    (hdividend : SlicePolyRep heap A length p dividend) :
    word3ArrayPoly p values = dividend := by
  rcases hdividend with ⟨coeffs, hreadCoeffs, hcoeffsSize, rfl⟩
  ext degree
  rw [coeff_word3ArrayPoly, coeff_coeffArrayPoly]
  by_cases hdegree : degree < length
  · have hdegreeValues : degree < values.size := by
      simpa [hvalues.1] using hdegree
    have hdegreeCoeffs : degree < coeffs.size := by
      simpa [hcoeffsSize] using hdegree
    rw [dif_pos hdegreeValues, dif_pos hdegreeCoeffs]
    rcases hprefix degree hdegree with ⟨value, hreadA, hreadW⟩
    have hcoeffRead := readU64s_get heap A length coeffs hreadCoeffs
      hcoeffsSize degree hdegree
    have hcoeff : coeffs[degree] = value :=
      Except.ok.inj (hcoeffRead.symm.trans hreadA)
    have hword : values[degree] = { lo := value, mid := 0, hi := 0 } :=
      Except.ok.inj ((hvalues.2 degree hdegreeValues).symm.trans hreadW)
    simp [hcoeff, hword, word3Value]
  · have hdegreeValues : ¬ degree < values.size := by
      simpa [hvalues.1] using hdegree
    have hdegreeCoeffs : ¬ degree < coeffs.size := by
      simpa [hcoeffsSize] using hdegree
    rw [dif_neg hdegreeValues, dif_neg hdegreeCoeffs]

/-- Content-level refinement of the generated initialization loop.  With
the C++ non-aliasing allocation precondition, every output W3 cell is exactly
the corresponding A limb zero-extended to three limbs. -/
theorem initW3Loop_refines (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenA i : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA) (hi : i ≤ lenA)
    (hregions : W3.region ≠ A.region)
    (hprefix : InitW3Prefix heap A W3 i) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA ∧
      RawHeap.SameLayout heap heap' ∧ InitW3Prefix heap' A W3 lenA := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    let word : Word3 := { lo := value, mid := 0, hi := 0 }
    rcases heap.writeWord3_of_valid W3 lenA i word hW3 hlt with
      ⟨heap1, hwrite⟩
    dsimp [word] at hwrite ⊢
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 i word hwrite
    have hprefix1 : InitW3Prefix heap1 A W3 (i + 1) := by
      intro j hj
      by_cases hji : j = i
      · subst j
        have hreadA := RawHeap.readU64_writeWord3_region_ne heap heap1
          W3 A i i word value hwrite hread hregions
        have hreadW := RawHeap.readWord3_writeWord3_same heap heap1
          W3 i word hwrite
        exact ⟨value, hreadA, hreadW⟩
      · have hjlt : j < i := by omega
        rcases hprefix j hjlt with ⟨old, hreadA, hreadW⟩
        have hreadA1 := RawHeap.readU64_writeWord3_region_ne heap heap1
          W3 A i j word old hwrite hreadA hregions
        have hreadW1 := RawHeap.readWord3_writeWord3_ne heap heap1
          W3 i j word { lo := old, mid := 0, hi := 0 }
          hwrite hreadW (Ne.symm hji)
        exact ⟨old, hreadA1, hreadW1⟩
    rcases initW3Loop_refines heap1 A W3 lenA (i + 1) hA1 hW31
      (by omega) hregions hprefix1 with
      ⟨heap2, hloop, hA2, hW32, hlayout2, hfull⟩
    refine ⟨heap2, hloop, hA2, hW32, ?_, hfull⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    have hieq : i = lenA := by omega
    subst i
    exact ⟨heap, rfl, hA, hW3, fun _ _ => Iff.rfl, hprefix⟩
termination_by lenA - i
decreasing_by omega

/-- The initialization loop writes only W3.  Any disjoint UInt64 slice keeps
its exact contents throughout the generated recursion. -/
theorem initW3Loop_preserves_u64_region_ne (heap : RawHeap)
    (A other : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenA i otherLen : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hOther : heap.ValidU64Slice other otherLen)
    (hi : i ≤ lenA) (hregions : W3.region ≠ other.region) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA ∧
      heap'.ValidU64Slice other otherLen ∧ RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' other otherLen := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    let word : Word3 := { lo := value, mid := 0, hi := 0 }
    rcases heap.writeWord3_of_valid W3 lenA i word hW3 hlt with
      ⟨heap1, hwrite⟩
    dsimp [word] at hwrite ⊢
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite other otherLen).mp hOther
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 i word hwrite
    have hsame1 : SameU64Prefix heap heap1 other otherLen := by
      intro k old hk hreadOld
      exact RawHeap.readU64_writeWord3_region_ne heap heap1 W3 other i k
        word old hwrite hreadOld hregions
    rcases initW3Loop_preserves_u64_region_ne heap1 A other W3 lenA
      (i + 1) otherLen hA1 hW31 hOther1 (by omega) hregions with
      ⟨heap2, hloop, hA2, hW32, hOther2, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hA2, hW32, hOther2, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k old hk hreadOld
      exact hsame2 k old hk (hsame1 k old hk hreadOld)
  next hnot =>
    exact ⟨heap, rfl, hA, hW3, hOther, fun _ _ => Iff.rfl,
      fun _ _ _ hread => hread⟩
termination_by lenA - i
decreasing_by omega

/-- Natural-language proof outline:

For `j ≤ d`, the divisor invariant makes `B[j]` readable.  The bound
`i+d < lenW3` implies `i+j < lenW3`, so the accumulator is readable and
writable.  The multiplication and carry routines are total word operations.
The W3 write preserves both B and W3 allocation invariants; recurse at `j+1`
with the strictly smaller measure `d+1-j`.  At `j>d`, the source loop exits. -/
def SameWord3Above (before after : RawHeap) (ptr : RawPtr Word3)
    (top upper : Nat) : Prop :=
  ∀ k, top < k → k < upper →
    before.readWord3 ptr k = after.readWord3 ptr k

def SameWord3Below (before after : RawHeap) (ptr : RawPtr Word3)
    (lower : Nat) : Prop :=
  ∀ k, k < lower → before.readWord3 ptr k = after.readWord3 ptr k

def SameU64Above (before after : RawHeap) (ptr : RawPtr UInt64)
    (lower upper : Nat) : Prop :=
  ∀ k value, lower ≤ k → k < upper →
    before.readU64 ptr k = .ok value → after.readU64 ptr k = .ok value

theorem addMulLoop_ok (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (other : RawPtr UInt64) (otherLen : Nat)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hOther : heap.ValidU64Slice other otherLen)
    (htop : i + d < lenW3) (hj : j ≤ d + 1) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      heap'.ValidU64Slice other otherLen ∧
      RawHeap.SameLayout heap heap' := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.fst product.snd
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite other otherLen).mp hOther
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    rcases addMulLoop_ok heap1 B W3 lenW3 i d (j + 1) c other otherLen
      hB1 hW31 hOther1 htop (by omega) with
      ⟨heap2, hloop, hB2, hW32, hOther2, hlayout2⟩
    refine ⟨heap2, hloop, hB2, hW32, hOther2, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, hOther, fun _ _ => Iff.rfl⟩
termination_by d + 1 - j
decreasing_by omega

/-- The generated inner loop writes only cells `i+j`; every W3 cell strictly
below its offset `i` is byte-for-byte unchanged. -/
theorem addMulLoop_preserves_below (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (htop : i + d < lenW3) (hj : j ≤ d + 1) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameWord3Below heap heap' W3 i := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.1 product.2
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    have hsame1 : SameWord3Below heap heap1 W3 i := by
      intro k hk
      rcases heap.readWord3_of_valid W3 lenW3 k hW3 (by omega) with
        ⟨old, hreadOld⟩
      have hpreserved := RawHeap.readWord3_writeWord3_ne heap heap1 W3
        (i + j) k accum' old hwrite hreadOld (by omega)
      rw [hreadOld, hpreserved]
    rcases addMulLoop_preserves_below heap1 B W3 lenW3 i d (j + 1) c
      hB1 hW31 htop (by omega) with
      ⟨heap2, hloop, hB2, hW32, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hB2, hW32, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k hk
      exact (hsame1 k hk).trans (hsame2 k hk)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, fun _ _ => Iff.rfl, fun _ _ => rfl⟩
termination_by d + 1 - j
decreasing_by omega

/-- Every write performed by the generated inner loop targets W3.  Hence an
actual non-aliasing UInt64 slice is preserved pointwise, not merely kept
allocated. -/
theorem addMulLoop_preserves_u64_region_ne (heap : RawHeap)
    (B other : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenW3 i d j otherLen : Nat) (c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hOther : heap.ValidU64Slice other otherLen)
    (htop : i + d < lenW3) (hj : j ≤ d + 1)
    (hregions : W3.region ≠ other.region) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      heap'.ValidU64Slice other otherLen ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' other otherLen := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.1 product.2
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite other otherLen).mp hOther
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    have hsame1 : SameU64Prefix heap heap1 other otherLen := by
      intro k value hk hread
      exact RawHeap.readU64_writeWord3_region_ne heap heap1 W3 other
        (i + j) k accum' value hwrite hread hregions
    rcases addMulLoop_preserves_u64_region_ne heap1 B other W3 lenW3 i d
      (j + 1) otherLen c hB1 hW31 hOther1 htop (by omega) hregions with
      ⟨heap2, hloop, hB2, hW32, hOther2, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hB2, hW32, hOther2, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k value hk hread
      exact hsame2 k value hk (hsame1 k value hk hread)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, hOther, fun _ _ => Iff.rfl,
      fun _ _ _ hread => hread⟩
termination_by d + 1 - j
decreasing_by omega

/-- Later quotient iterations update only `W3[i..i+d]`; every allocated cell
strictly above that interval is byte-for-byte unchanged. -/
theorem addMulLoop_preserves_above (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (htop : i + d < lenW3) (hj : j ≤ d + 1) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameWord3Above heap heap' W3 (i + d) lenW3 := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.1 product.2
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    have hsame1 : SameWord3Above heap heap1 W3 (i + d) lenW3 := by
      intro k hkTop hkLen
      rcases heap.readWord3_of_valid W3 lenW3 k hW3 hkLen with
        ⟨old, hreadOld⟩
      have hpreserved := RawHeap.readWord3_writeWord3_ne heap heap1 W3
        (i + j) k accum' old hwrite hreadOld (by omega)
      rw [hreadOld, hpreserved]
    rcases addMulLoop_preserves_above heap1 B W3 lenW3 i d (j + 1) c
      hB1 hW31 htop (by omega) with
      ⟨heap2, hloop, hB2, hW32, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hB2, hW32, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k hkTop hkLen
      exact (hsame1 k hkTop hkLen).trans (hsame2 k hkTop hkLen)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, fun _ _ => Iff.rfl,
      fun _ _ _ => rfl⟩
termination_by d + 1 - j
decreasing_by omega

/-- Content semantics of the exact body executed by one generated
`addMulLoop` iteration.  It reads the real B/W3 cells, applies the generated
`_umul128` and `_add_carry3`, writes that machine result, and relates the
observed output cell to addition of `c*bj` modulo the 192-bit accumulator
width. -/
theorem addMulCell_refines (heap heap' : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (bIndex wIndex : Nat) (c bj : UInt64)
    (accum : Word3)
    (hreadB : heap.readU64 B bIndex = .ok bj)
    (hreadW : heap.readWord3 W3 wIndex = .ok accum)
    (hregions : W3.region ≠ B.region)
    (hwrite : heap.writeWord3 W3 wIndex
      (let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
       Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
         accum product.1 product.2) = .ok heap') :
    ∃ accum', heap.readWord3 W3 wIndex = .ok accum ∧
      heap'.readWord3 W3 wIndex = .ok accum' ∧
      heap'.readU64 B bIndex = .ok bj ∧
      Nat.ModEq (limbBase ^ 3) (word3Value accum')
        (word3Value accum + c.toNat * bj.toNat) := by
  let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
  let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
    accum product.1 product.2
  have hwrite' : heap.writeWord3 W3 wIndex accum' = .ok heap' := by
    simpa [accum', product] using hwrite
  have hreadW' := RawHeap.readWord3_writeWord3_same heap heap' W3 wIndex
    accum' hwrite'
  have hreadB' := RawHeap.readU64_writeWord3_region_ne heap heap'
    W3 B wIndex bIndex accum' bj hwrite' hreadB hregions
  refine ⟨accum', hreadW, hreadW', hreadB', ?_⟩
  simpa [accum', product] using addMulWord3_modEq accum c bj

/-- Exact, non-wrapping form of `addMulCell_refines`.  The additional bound is
the concrete machine-capacity obligation that will be maintained by the
quotient loop; it does not replace the generated update with mathematical
addition. -/
theorem addMulCell_refines_exact (heap heap' : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (bIndex wIndex : Nat) (c bj : UInt64) (accum : Word3)
    (hreadB : heap.readU64 B bIndex = .ok bj)
    (hreadW : heap.readWord3 W3 wIndex = .ok accum)
    (hregions : W3.region ≠ B.region)
    (hbound : word3Value accum + c.toNat * bj.toNat < limbBase ^ 3)
    (hwrite : heap.writeWord3 W3 wIndex
      (let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
       Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
         accum product.1 product.2) = .ok heap') :
    ∃ accum', heap.readWord3 W3 wIndex = .ok accum ∧
      heap'.readWord3 W3 wIndex = .ok accum' ∧
      heap'.readU64 B bIndex = .ok bj ∧
      word3Value accum' = word3Value accum + c.toNat * bj.toNat := by
  let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
  let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
    accum product.1 product.2
  have hwrite' : heap.writeWord3 W3 wIndex accum' = .ok heap' := by
    simpa [accum', product] using hwrite
  have hreadW' := RawHeap.readWord3_writeWord3_same heap heap' W3 wIndex
    accum' hwrite'
  have hreadB' := RawHeap.readU64_writeWord3_region_ne heap heap'
    W3 B wIndex bIndex accum' bj hwrite' hreadB hregions
  refine ⟨accum', hreadW, hreadW', hreadB', ?_⟩
  simpa [accum', product] using addMulWord3_exact accum c bj hbound

/-- One actual generated multiply/add write raises the uniform accumulation
budget by one.  The proof obtains the no-wrap premise from the C++ `size_t`
iteration bound and canonical residues, then uses the exact machine-update
theorem above. -/
theorem writeAddMul_preserves_budget (heap heap' : RawHeap)
    (W3 : RawPtr Word3) (length index count : Nat) (p c bj : UInt64)
    (accum : Word3)
    (hvalid : heap.ValidWord3Slice W3 length)
    (hindex : index < length)
    (hbudget : Word3AccumulationBudget heap W3 length p count)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) (hb : bj.toNat < p.toNat)
    (hread : heap.readWord3 W3 index = .ok accum)
    (hwrite : heap.writeWord3 W3 index
      (let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
       Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
         accum product.1 product.2) = .ok heap') :
    Word3AccumulationBudget heap' W3 length p (count + 1) := by
  let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
  let written := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
    accum product.1 product.2
  have hwrite' : heap.writeWord3 W3 index written = .ok heap' := by
    simpa [written, product] using hwrite
  have hbase : word3Value accum ≤
      (p.toNat - 1) + count * (p.toNat - 1) ^ 2 :=
    hbudget index accum hindex hread
  have hc' : c.toNat ≤ p.toNat - 1 := by omega
  have hb' : bj.toNat ≤ p.toNat - 1 := by omega
  have hproduct : c.toNat * bj.toNat ≤ (p.toNat - 1) ^ 2 := by
    rw [pow_two]
    exact Nat.mul_le_mul hc' hb'
  have hsum : word3Value accum + c.toNat * bj.toNat ≤
      (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          ((p.toNat - 1) + count * (p.toNat - 1) ^ 2) +
            (p.toNat - 1) ^ 2 := Nat.add_le_add hbase hproduct
      _ = (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by ring
  have hcapacity := lazyAccumulation_word3_budget p (count + 1)
    (p.toNat - 1) hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hnoWrap : word3Value accum + c.toNat * bj.toNat < limbBase ^ 3 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := hsum
      _ < p.toNat * limbBase ^ 2 := hcapacity
      _ < limbBase ^ 3 := by
        calc
          p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
            Nat.mul_lt_mul_of_pos_right hpB (pow_pos (by omega) 2)
          _ = limbBase ^ 2 * limbBase := Nat.mul_comm _ _
          _ = limbBase ^ 3 := by ring
  have hwritten : word3Value written =
      word3Value accum + c.toNat * bj.toNat := by
    simpa [written, product] using addMulWord3_exact accum c bj hnoWrap
  intro k out hk hreadOut
  by_cases hki : k = index
  · subst k
    have hsame := RawHeap.readWord3_writeWord3_same heap heap' W3 index
      written hwrite'
    have hout : out = written := Except.ok.inj (hreadOut.symm.trans hsame)
    subst out
    rw [hwritten]
    exact hsum
  · rcases heap.readWord3_of_valid W3 length k hvalid hk with
      ⟨old, hreadOld⟩
    have hpreserved := RawHeap.readWord3_writeWord3_ne heap heap' W3
      index k written old hwrite' hreadOld (Ne.symm hki)
    have hout : out = old := Except.ok.inj (hreadOut.symm.trans hpreserved)
    subst out
    have hold := hbudget k old hk hreadOld
    calc
      word3Value old ≤
          (p.toNat - 1) + count * (p.toNat - 1) ^ 2 := hold
      _ ≤ (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by
        apply Nat.add_le_add_left
        apply Nat.mul_le_mul_right
        omega

/-- The accumulation budget supplies exactly the high-limb precondition of
the generated `_lll_mod_preinv` call.  Thus quotient/remainder reduction can
consume a property derived from actual writes rather than assume `hi < p` at
each call site. -/
theorem word3_hi_lt_of_accumulation_budget (heap : RawHeap)
    (W3 : RawPtr Word3) (length count index : Nat) (p : UInt64)
    (accum : Word3)
    (hbudget : Word3AccumulationBudget heap W3 length p count)
    (hp : 1 < p.toNat) (hcount : count < limbBase)
    (hindex : index < length)
    (hread : heap.readWord3 W3 index = .ok accum) :
    accum.hi.toNat < p.toNat := by
  have hvalue := hbudget index accum hindex hread
  have hcapacity := lazyAccumulation_word3_budget p count
    (p.toNat - 1) hp hcount (by omega)
  apply word3_hi_lt_of_value_lt accum p
  exact lt_of_le_of_lt hvalue hcapacity

theorem accumulationBudget_writeU64_region_ne (heap heap' : RawHeap)
    (dst : RawPtr UInt64) (W3 : RawPtr Word3)
    (writeIndex length count : Nat) (value : UInt64) (p : UInt64)
    (hvalid : heap.ValidWord3Slice W3 length)
    (hbudget : Word3AccumulationBudget heap W3 length p count)
    (hregions : dst.region ≠ W3.region)
    (hwrite : heap.writeU64 dst writeIndex value = .ok heap') :
    Word3AccumulationBudget heap' W3 length p count := by
  intro k accum hk hread'
  rcases heap.readWord3_of_valid W3 length k hvalid hk with
    ⟨old, hread⟩
  have hpreserved := RawHeap.readWord3_writeU64_region_ne heap heap'
    dst W3 writeIndex k value old hwrite hread hregions
  have heq : accum = old := Except.ok.inj (hread'.symm.trans hpreserved)
  subst accum
  exact hbudget k old hk hread

theorem canonicalPrefix_writeU64_region_ne (heap heap' : RawHeap)
    (dst src : RawPtr UInt64) (writeIndex length : Nat)
    (value p : UInt64)
    (hvalid : heap.ValidU64Slice src length)
    (hcanonical : CanonicalU64Prefix heap src length p)
    (hregions : dst.region ≠ src.region)
    (hwrite : heap.writeU64 dst writeIndex value = .ok heap') :
    CanonicalU64Prefix heap' src length p := by
  intro k observed hk hread'
  rcases heap.readU64_of_valid src length k hvalid hk with
    ⟨old, hread⟩
  have hpreserved := RawHeap.readU64_writeU64_ne heap heap' dst src
    writeIndex k value old hwrite hread (Or.inl hregions)
  have heq : observed = old := Except.ok.inj (hread'.symm.trans hpreserved)
  subst observed
  exact hcanonical k old hk hread

/-- Exact staged transition for one inner-loop cell.  Only the cell at
`offset + processed` receives the extra allowance; every other cell keeps the
same bound. -/
theorem writeAddMul_preserves_staged_budget (heap heap' : RawHeap)
    (W3 : RawPtr Word3) (length count offset processed : Nat)
    (p c bj : UInt64) (accum : Word3)
    (hvalid : heap.ValidWord3Slice W3 length)
    (hindex : offset + processed < length)
    (hstage : Word3StagedBudget heap W3 length p count offset processed)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) (hb : bj.toNat < p.toNat)
    (hread : heap.readWord3 W3 (offset + processed) = .ok accum)
    (hwrite : heap.writeWord3 W3 (offset + processed)
      (let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
       Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
         accum product.1 product.2) = .ok heap') :
    Word3StagedBudget heap' W3 length p count offset (processed + 1) := by
  let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
  let written := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
    accum product.1 product.2
  have hwrite' : heap.writeWord3 W3 (offset + processed) written =
      .ok heap' := by simpa [written, product] using hwrite
  have hnotOld : ¬(offset ≤ offset + processed ∧
      offset + processed < offset + processed) := by omega
  have hbaseRaw := hstage (offset + processed) accum hindex hread
  have hbase : word3Value accum ≤
      (p.toNat - 1) + count * (p.toNat - 1) ^ 2 := by
    simpa [hnotOld] using hbaseRaw
  have hc' : c.toNat ≤ p.toNat - 1 := by omega
  have hb' : bj.toNat ≤ p.toNat - 1 := by omega
  have hproduct : c.toNat * bj.toNat ≤ (p.toNat - 1) ^ 2 := by
    rw [pow_two]
    exact Nat.mul_le_mul hc' hb'
  have hsum : word3Value accum + c.toNat * bj.toNat ≤
      (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          ((p.toNat - 1) + count * (p.toNat - 1) ^ 2) +
            (p.toNat - 1) ^ 2 := Nat.add_le_add hbase hproduct
      _ = (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by ring
  have hcapacity := lazyAccumulation_word3_budget p (count + 1)
    (p.toNat - 1) hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hnoWrap : word3Value accum + c.toNat * bj.toNat < limbBase ^ 3 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := hsum
      _ < p.toNat * limbBase ^ 2 := hcapacity
      _ < limbBase ^ 3 := by
        calc
          p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
            Nat.mul_lt_mul_of_pos_right hpB (pow_pos (by omega) 2)
          _ = limbBase ^ 2 * limbBase := Nat.mul_comm _ _
          _ = limbBase ^ 3 := by ring
  have hwritten : word3Value written =
      word3Value accum + c.toNat * bj.toNat := by
    simpa [written, product] using addMulWord3_exact accum c bj hnoWrap
  intro k out hk hreadOut
  by_cases hki : k = offset + processed
  · subst k
    have hsame := RawHeap.readWord3_writeWord3_same heap heap' W3
      (offset + processed) written hwrite'
    have hout : out = written := Except.ok.inj (hreadOut.symm.trans hsame)
    subst out
    have hnew : offset ≤ offset + processed ∧
        offset + processed < offset + (processed + 1) := by omega
    simpa [hnew, hwritten] using hsum
  · rcases heap.readWord3_of_valid W3 length k hvalid hk with
      ⟨old, hreadOld⟩
    have hpreserved := RawHeap.readWord3_writeWord3_ne heap heap' W3
      (offset + processed) k written old hwrite' hreadOld (Ne.symm hki)
    have hout : out = old := Except.ok.inj (hreadOut.symm.trans hpreserved)
    subst out
    have hold := hstage k old hk hreadOld
    by_cases hin : offset ≤ k ∧ k < offset + processed
    · have hin' : offset ≤ k ∧ k < offset + (processed + 1) := by omega
      simpa [hin, hin'] using hold
    · have hin' : ¬(offset ≤ k ∧ k < offset + (processed + 1)) := by
        intro hnew
        apply hin
        constructor
        · exact hnew.1
        · omega
      simpa [hin, hin'] using hold

def SameWord3PrefixAt (before after : RawHeap) (ptr : RawPtr Word3)
    (offset length : Nat) : Prop :=
  ∀ k value, k < length → before.readWord3 ptr (offset + k) = .ok value →
    after.readWord3 ptr (offset + k) = .ok value

def Word3ZeroModRange (heap : RawHeap) (ptr : RawPtr Word3)
    (lower upper p : Nat) : Prop :=
  ∀ k value, lower ≤ k → k < upper →
    heap.readWord3 ptr k = .ok value → word3Value value % p = 0

theorem zeroModRange_of_same_above (before after : RawHeap)
    (W3 : RawPtr Word3) (top upper p : Nat)
    (hsame : SameWord3Above before after W3 top upper)
    (hzero : Word3ZeroModRange before W3 (top + 1) upper p) :
    Word3ZeroModRange after W3 (top + 1) upper p := by
  intro k value hlow hhigh hreadAfter
  have hsameRead := hsame k (by omega) hhigh
  apply hzero k value hlow hhigh
  rwa [hsameRead]

theorem zeroModRange_writeU64_region_ne (before after : RawHeap)
    (dst : RawPtr UInt64) (W3 : RawPtr Word3)
    (writeIndex length lower upper p : Nat) (written : UInt64)
    (hvalid : before.ValidWord3Slice W3 length)
    (hupper : upper ≤ length)
    (hzero : Word3ZeroModRange before W3 lower upper p)
    (hregions : dst.region ≠ W3.region)
    (hwrite : before.writeU64 dst writeIndex written = .ok after) :
    Word3ZeroModRange after W3 lower upper p := by
  intro k value hlow hhigh hreadAfter
  rcases before.readWord3_of_valid W3 length k hvalid (by omega) with
    ⟨old, hreadOld⟩
  have hpreserved := RawHeap.readWord3_writeU64_region_ne before after
    dst W3 writeIndex k written old hwrite hreadOld hregions
  have hvalue : value = old :=
    Except.ok.inj (hreadAfter.symm.trans hpreserved)
  subst value
  exact hzero k old hlow hhigh hreadOld

theorem zeroModRange_extend_down (heap : RawHeap) (W3 : RawPtr Word3)
    (lower upper p : Nat)
    (hcell : ∀ value, heap.readWord3 W3 lower = .ok value →
      word3Value value % p = 0)
    (hrest : Word3ZeroModRange heap W3 (lower + 1) upper p) :
    Word3ZeroModRange heap W3 lower upper p := by
  intro k value hlow hhigh hread
  by_cases hk : k = lower
  · subst k
    exact hcell value hread
  · exact hrest k value (by omega) hhigh hread

def AddMulRangeRep (before after : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (offset start stop : Nat) (c : UInt64) : Prop :=
  ∀ k, start ≤ k → k ≤ stop → ∀ bj accum,
    before.readU64 B k = .ok bj →
    before.readWord3 W3 (offset + k) = .ok accum →
    ∃ accum', after.readWord3 W3 (offset + k) = .ok accum' ∧
      Nat.ModEq (limbBase ^ 3) (word3Value accum')
        (word3Value accum + c.toNat * bj.toNat)

def ExactAddMulRangeRep (before after : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (offset start stop : Nat) (c : UInt64) : Prop :=
  ∀ k, start ≤ k → k ≤ stop → ∀ bj accum,
    before.readU64 B k = .ok bj →
    before.readWord3 W3 (offset + k) = .ok accum →
    ∃ accum', after.readWord3 W3 (offset + k) = .ok accum' ∧
      word3Value accum' = word3Value accum + c.toNat * bj.toNat

/-- Array-observation form of the exact raw range theorem.  Both accumulator
values are obtained from the corresponding pre/post C++ heaps; the divisor
coefficient is likewise an actual raw read. -/
theorem exactAddMulRange_observed (before after : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (length offset start stop : Nat) (c : UInt64)
    (beforeValues afterValues : Array Word3)
    (hrange : ExactAddMulRangeRep before after B W3 offset start stop c)
    (hB : before.ValidU64Slice B (stop + 1))
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hafter : Word3SliceRep after W3 length afterValues)
    (k : Nat) (hstart : start ≤ k) (hstop : k ≤ stop)
    (hindex : offset + k < length) :
    ∃ bj : UInt64, ∃ accum out : Word3,
      before.readU64 B k = .ok bj ∧
      beforeValues[offset + k]? = some accum ∧
      afterValues[offset + k]? = some out ∧
      word3Value out = word3Value accum + c.toNat * bj.toNat := by
  have hbeforeIndex : offset + k < beforeValues.size := by
    simpa [hbefore.1] using hindex
  have hafterIndex : offset + k < afterValues.size := by
    simpa [hafter.1] using hindex
  rcases before.readU64_of_valid B (stop + 1) k hB (by omega) with
    ⟨bj, hreadB⟩
  have hreadBefore := hbefore.2 (offset + k) hbeforeIndex
  rcases hrange k hstart hstop bj beforeValues[offset + k]
      hreadB hreadBefore with ⟨out, hreadOut, hout⟩
  have houtEq : afterValues[offset + k] = out :=
    Except.ok.inj ((hafter.2 (offset + k) hafterIndex).symm.trans hreadOut)
  refine ⟨bj, beforeValues[offset + k], out, hreadB, ?_, ?_, hout⟩
  · simp [hbeforeIndex]
  · simp [hafterIndex, houtEq]

/-- Coefficient semantics of a complete generated multiply/add range.  This
is the pointwise convolution law needed to lift the raw loop to polynomial
addition; cells outside `offset..offset+d` are proved unchanged explicitly. -/
theorem word3ArrayPoly_coeff_addMul (before after : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (length offset d p : Nat) (c : UInt64)
    (beforeValues afterValues : Array Word3)
    (divisor : Polynomial (ZMod p))
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hafter : Word3SliceRep after W3 length afterValues)
    (hdivisor : SlicePolyRep before B (d + 1) p divisor)
    (hbelow : SameWord3Below before after W3 offset)
    (habove : SameWord3Above before after W3 (offset + d) length)
    (hrange : ExactAddMulRangeRep before after B W3 offset 0 d c)
    (htop : offset + d < length) (degree : Nat) :
    (word3ArrayPoly p afterValues).coeff degree =
      if degree < offset ∨ offset + d < degree then
        (word3ArrayPoly p beforeValues).coeff degree
      else
        (word3ArrayPoly p beforeValues).coeff degree +
          (c.toNat : ZMod p) * divisor.coeff (degree - offset) := by
  rw [coeff_word3ArrayPoly, coeff_word3ArrayPoly]
  by_cases hdegree : degree < length
  · have hbeforeDegree : degree < beforeValues.size := by
      simpa [hbefore.1] using hdegree
    have hafterDegree : degree < afterValues.size := by
      simpa [hafter.1] using hdegree
    rw [dif_pos hafterDegree, dif_pos hbeforeDegree]
    by_cases hlow : degree < offset
    · rw [if_pos (Or.inl hlow)]
      have hsame := hbelow degree hlow
      have hreadBefore := hbefore.2 degree hbeforeDegree
      have hreadAfter := hafter.2 degree hafterDegree
      rw [hreadBefore, hreadAfter] at hsame
      exact congrArg (fun value => (word3Value value : ZMod p))
        (Except.ok.inj hsame).symm
    · by_cases hhigh : offset + d < degree
      · rw [if_pos (Or.inr hhigh)]
        have hsame := habove degree hhigh hdegree
        have hreadBefore := hbefore.2 degree hbeforeDegree
        have hreadAfter := hafter.2 degree hafterDegree
        rw [hreadBefore, hreadAfter] at hsame
        exact congrArg (fun value => (word3Value value : ZMod p))
          (Except.ok.inj hsame).symm
      · rw [if_neg (by simp [hlow, hhigh])]
        let k := degree - offset
        have hk : degree = offset + k := by
          dsimp [k]
          omega
        have hkd : k ≤ d := by
          dsimp [k]
          omega
        rcases slicePolyRep_coeff before B (d + 1) p divisor hdivisor
            k (by omega) with ⟨bj, hreadB, hcoeff⟩
        have hreadBefore := hbefore.2 degree hbeforeDegree
        rcases hrange k (by omega) hkd bj beforeValues[degree]
            hreadB (by simpa [hk] using hreadBefore) with
          ⟨out, hreadOut, hout⟩
        have houtEq : afterValues[degree] = out :=
          Except.ok.inj ((hafter.2 degree hafterDegree).symm.trans
            (by simpa [hk] using hreadOut))
        subst out
        rw [hcoeff]
        simpa [k, hout, Nat.cast_add, Nat.cast_mul]
  · have hbeforeDegree : ¬ degree < beforeValues.size := by
      simpa [hbefore.1] using hdegree
    have hafterDegree : ¬ degree < afterValues.size := by
      simpa [hafter.1] using hdegree
    rw [dif_neg hafterDegree, dif_neg hbeforeDegree]
    have hhigh : offset + d < degree := by omega
    simp [hhigh]

/-- Polynomial form of the generated multiply/add loop: the observed W3
accumulator gains exactly the shifted scalar multiple written by C++.
No polynomial operation is used to execute or replace the raw loop. -/
theorem word3ArrayPoly_addMul (before after : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (length offset d p : Nat) (c : UInt64)
    (beforeValues afterValues : Array Word3)
    (divisor : Polynomial (ZMod p))
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hafter : Word3SliceRep after W3 length afterValues)
    (hdivisor : SlicePolyRep before B (d + 1) p divisor)
    (hbelow : SameWord3Below before after W3 offset)
    (habove : SameWord3Above before after W3 (offset + d) length)
    (hrange : ExactAddMulRangeRep before after B W3 offset 0 d c)
    (htop : offset + d < length) :
    word3ArrayPoly p afterValues = word3ArrayPoly p beforeValues +
      Polynomial.monomial offset (c.toNat : ZMod p) * divisor := by
  ext degree
  rw [Polynomial.coeff_add]
  have hstep := word3ArrayPoly_coeff_addMul before after B W3 length
    offset d p c beforeValues afterValues divisor hbefore hafter hdivisor
    hbelow habove hrange htop degree
  have hproduct :
      (Polynomial.monomial offset (c.toNat : ZMod p) * divisor).coeff degree =
        if offset ≤ degree then
          (c.toNat : ZMod p) * divisor.coeff (degree - offset)
        else 0 := by
    rw [← Polynomial.C_mul_X_pow_eq_monomial, mul_assoc,
      Polynomial.coeff_C_mul, Polynomial.coeff_X_pow_mul']
    by_cases hoffset : offset ≤ degree <;> simp [hoffset]
  rw [hproduct]
  by_cases hlow : degree < offset
  · rw [if_neg (by omega)]
    simpa [hlow] using hstep
  · rw [if_pos (by omega)]
    by_cases hhigh : offset + d < degree
    · have hzero := slicePolyRep_coeff_zero_of_length_le before B
        (d + 1) p divisor hdivisor (degree - offset) (by omega)
      rw [hzero, mul_zero, add_zero]
      simpa [hlow, hhigh] using hstep
    · simpa [hlow, hhigh] using hstep

theorem zmod_cast_uint64_complement (p qi : UInt64)
    (hqi : qi.toNat ≤ p.toNat) :
    ((p - qi).toNat : ZMod p.toNat) = -(qi.toNat : ZMod p.toNat) := by
  rw [UInt64.toNat_sub_of_le _ _ hqi, Nat.cast_sub hqi]
  have hpzero : (p.toNat : ZMod p.toNat) = 0 := by
    exact ZMod.natCast_self p.toNat
  rw [hpzero, zero_sub]

theorem add_complement_monomial_mul_eq_sub (p qi : UInt64)
    (i : Nat) (base divisor : Polynomial (ZMod p.toNat))
    (hqi : qi.toNat ≤ p.toNat) :
    base + Polynomial.monomial i ((p - qi).toNat : ZMod p.toNat) * divisor =
      base - Polynomial.monomial i (qi.toNat : ZMod p.toNat) * divisor := by
  rw [zmod_cast_uint64_complement p qi hqi,
    Polynomial.monomial_neg, neg_mul, sub_eq_add_neg]

theorem quotient_step_compose (p : Nat)
    (before afterStep afterRec lowerQ term divisor : Polynomial (ZMod p))
    (hstep : afterStep = before - term * divisor)
    (hrec : afterRec = afterStep - lowerQ * divisor) :
    afterRec = before - (lowerQ + term) * divisor := by
  rw [hrec, hstep]
  ring

/-- Semantic assembly used by the successor case of the generated descending
quotient loop.  The newly written raw Q cell and the recursive lower prefix
are joined before composing their two W3 subtraction steps. -/
theorem quotient_step_finalize (heap : RawHeap) (Q : RawPtr UInt64)
    (p i : Nat) (qi : UInt64)
    (lowerQ fullQ before afterStep afterRec divisor : Polynomial (ZMod p))
    (hlower : SlicePolyRep heap Q i p lowerQ)
    (hfull : SlicePolyRep heap Q (i + 1) p fullQ)
    (hread : heap.readU64 Q i = .ok qi)
    (hstep : afterStep = before -
      Polynomial.monomial i (qi.toNat : ZMod p) * divisor)
    (hrec : afterRec = afterStep - lowerQ * divisor) :
    afterRec = before - fullQ * divisor := by
  have hq := slicePolyRep_succ_eq_add_monomial heap Q i p qi lowerQ
    fullQ hlower hfull hread
  rw [hq]
  exact quotient_step_compose p before afterStep afterRec lowerQ
    (Polynomial.monomial i (qi.toNat : ZMod p)) divisor hstep hrec

/-- A non-aliasing raw Q write leaves the complete W3 polynomial observation
unchanged.  Both sides are reconstructed from their respective heaps. -/
theorem word3ArrayPoly_writeU64_region_ne (before after : RawHeap)
    (dst : RawPtr UInt64) (W3 : RawPtr Word3)
    (writeIndex length p : Nat) (written : UInt64)
    (beforeValues afterValues : Array Word3)
    (hvalid : before.ValidWord3Slice W3 length)
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hafter : Word3SliceRep after W3 length afterValues)
    (hregions : dst.region ≠ W3.region)
    (hwrite : before.writeU64 dst writeIndex written = .ok after) :
    word3ArrayPoly p afterValues = word3ArrayPoly p beforeValues := by
  ext degree
  rw [coeff_word3ArrayPoly, coeff_word3ArrayPoly]
  by_cases hdegree : degree < length
  · have hbeforeDegree : degree < beforeValues.size := by
      simpa [hbefore.1] using hdegree
    have hafterDegree : degree < afterValues.size := by
      simpa [hafter.1] using hdegree
    rw [dif_pos hafterDegree, dif_pos hbeforeDegree]
    have hreadBefore := hbefore.2 degree hbeforeDegree
    have hpreserved := RawHeap.readWord3_writeU64_region_ne before after
      dst W3 writeIndex degree written beforeValues[degree] hwrite
      hreadBefore hregions
    have hvalue : afterValues[degree] = beforeValues[degree] :=
      Except.ok.inj ((hafter.2 degree hafterDegree).symm.trans hpreserved)
    rw [hvalue]
  · have hbeforeDegree : ¬ degree < beforeValues.size := by
      simpa [hbefore.1] using hdegree
    have hafterDegree : ¬ degree < afterValues.size := by
      simpa [hafter.1] using hdegree
    rw [dif_neg hafterDegree, dif_neg hbeforeDegree]

/-- Algebraic/raw assembly of the actual `qi = 0` successor branch.  The
generated Q write is retained even though no W3 addMul is executed. -/
theorem quotient_zero_successor_finalize (before afterWrite final : RawHeap)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (length i p : Nat) (qi : UInt64)
    (lowerQ divisor : Polynomial (ZMod p))
    (beforeValues writeValues recValues finalValues : Array Word3)
    (hvalidW : before.ValidWord3Slice W3 length)
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hwriteValues : Word3SliceRep afterWrite W3 length writeValues)
    (hrecValues : Word3SliceRep afterWrite W3 length recValues)
    (hfinalValues : Word3SliceRep final W3 length finalValues)
    (hwrite : before.writeU64 Q i qi = .ok afterWrite)
    (hQW : Q.region ≠ W3.region)
    (hqi : qi = 0)
    (hlower : SlicePolyRep final Q i p lowerQ)
    (hvalidQ : final.ValidU64Slice Q (i + 1))
    (hreadFinal : final.readU64 Q i = .ok qi)
    (hrec : word3ArrayPoly p finalValues =
      word3ArrayPoly p recValues - lowerQ * divisor) :
    ∃ fullQ : Polynomial (ZMod p),
      SlicePolyRep final Q (i + 1) p fullQ ∧
      word3ArrayPoly p finalValues =
        word3ArrayPoly p beforeValues - fullQ * divisor := by
  have hwriteEq := word3ArrayPoly_writeU64_region_ne before afterWrite Q W3
    i length p qi beforeValues writeValues hvalidW hbefore hwriteValues
    hQW hwrite
  have hrecValuesEq : recValues = writeValues :=
    word3SliceRep_eq afterWrite W3 length recValues writeValues hrecValues
      hwriteValues
  have hstep : word3ArrayPoly p recValues =
      word3ArrayPoly p beforeValues -
        Polynomial.monomial i (qi.toNat : ZMod p) * divisor := by
    subst qi
    rw [hrecValuesEq, hwriteEq]
    simp
  rcases slicePolyRep_extend_exists final Q i p qi lowerQ hvalidQ hlower
      hreadFinal with ⟨fullQ, hfull, _⟩
  refine ⟨fullQ, hfull, ?_⟩
  exact quotient_step_finalize final Q p i qi lowerQ fullQ
    (word3ArrayPoly p beforeValues) (word3ArrayPoly p recValues)
    (word3ArrayPoly p finalValues) divisor hlower hfull hreadFinal hstep hrec

/-- Algebraic/raw assembly of the actual nonzero successor branch, including
the generated complement addMul between the Q write and recursive call. -/
theorem quotient_nonzero_successor_finalize
    (before afterWrite afterAdd final : RawHeap)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (length i p : Nat) (modulus qi : UInt64)
    (lowerQ divisor : Polynomial (ZMod p))
    (beforeValues writeValues addBeforeValues addValues recValues
      finalValues : Array Word3)
    (hvalidW : before.ValidWord3Slice W3 length)
    (hbefore : Word3SliceRep before W3 length beforeValues)
    (hwriteValues : Word3SliceRep afterWrite W3 length writeValues)
    (haddBefore : Word3SliceRep afterWrite W3 length addBeforeValues)
    (haddValues : Word3SliceRep afterAdd W3 length addValues)
    (hrecValues : Word3SliceRep afterAdd W3 length recValues)
    (hfinalValues : Word3SliceRep final W3 length finalValues)
    (hwrite : before.writeU64 Q i qi = .ok afterWrite)
    (hQW : Q.region ≠ W3.region)
    (hp : p = modulus.toNat)
    (hqi : qi.toNat ≤ modulus.toNat)
    (hadd : word3ArrayPoly p addValues =
      word3ArrayPoly p addBeforeValues +
        Polynomial.monomial i ((modulus - qi).toNat : ZMod p) * divisor)
    (hlower : SlicePolyRep final Q i p lowerQ)
    (hvalidQ : final.ValidU64Slice Q (i + 1))
    (hreadFinal : final.readU64 Q i = .ok qi)
    (hrec : word3ArrayPoly p finalValues =
      word3ArrayPoly p recValues - lowerQ * divisor) :
    ∃ fullQ : Polynomial (ZMod p),
      SlicePolyRep final Q (i + 1) p fullQ ∧
      word3ArrayPoly p finalValues =
        word3ArrayPoly p beforeValues - fullQ * divisor := by
  subst p
  have hwriteEq := word3ArrayPoly_writeU64_region_ne before afterWrite Q W3
    i length modulus.toNat qi beforeValues writeValues hvalidW hbefore
    hwriteValues hQW hwrite
  have haddBeforeEq : addBeforeValues = writeValues :=
    word3SliceRep_eq afterWrite W3 length addBeforeValues writeValues
      haddBefore hwriteValues
  have hrecValuesEq : recValues = addValues :=
    word3SliceRep_eq afterAdd W3 length recValues addValues hrecValues
      haddValues
  have hstepAdd : word3ArrayPoly modulus.toNat addValues =
      word3ArrayPoly modulus.toNat beforeValues -
        Polynomial.monomial i (qi.toNat : ZMod modulus.toNat) * divisor := by
    calc
      word3ArrayPoly modulus.toNat addValues =
          word3ArrayPoly modulus.toNat beforeValues +
            Polynomial.monomial i
              ((modulus - qi).toNat : ZMod modulus.toNat) * divisor := by
        rw [hadd, haddBeforeEq, hwriteEq]
      _ = word3ArrayPoly modulus.toNat beforeValues -
          Polynomial.monomial i (qi.toNat : ZMod modulus.toNat) * divisor :=
        add_complement_monomial_mul_eq_sub modulus qi i
          (word3ArrayPoly modulus.toNat beforeValues) divisor hqi
  have hstep : word3ArrayPoly modulus.toNat recValues =
      word3ArrayPoly modulus.toNat beforeValues -
        Polynomial.monomial i (qi.toNat : ZMod modulus.toNat) * divisor := by
    rw [hrecValuesEq]
    exact hstepAdd
  rcases slicePolyRep_extend_exists final Q i modulus.toNat qi lowerQ
      hvalidQ hlower hreadFinal with ⟨fullQ, hfull, _⟩
  refine ⟨fullQ, hfull, ?_⟩
  exact quotient_step_finalize final Q modulus.toNat i qi lowerQ fullQ
    (word3ArrayPoly modulus.toNat beforeValues)
    (word3ArrayPoly modulus.toNat recValues)
    (word3ArrayPoly modulus.toNat finalValues) divisor hlower hfull
    hreadFinal hstep hrec

/-- Algebraic core of the source's leading-coefficient cancellation.  It is
stated over naturals so the subsequent raw-memory theorem only has to supply
the exact generated write and canonical UInt64 observations. -/
theorem add_complement_eliminates_mod (p old qi lead r : Nat)
    (hp : 0 < p) (hqi : qi ≤ p)
    (hold : old % p = r)
    (hproduct : (qi * lead) % p = r) :
    (old + (p - qi) * lead) % p = 0 := by
  have hcongr : Nat.ModEq p old (qi * lead) := by
    exact hold.trans hproduct.symm
  have hadd := hcongr.add_right ((p - qi) * lead)
  have hrhs : qi * lead + (p - qi) * lead = p * lead := by
    calc
      qi * lead + (p - qi) * lead = (qi + (p - qi)) * lead := by ring
      _ = p * lead := by rw [Nat.add_sub_of_le hqi]
  rw [hrhs] at hadd
  have hzero : Nat.ModEq p (p * lead) 0 :=
    (show p ∣ p * lead from ⟨lead, rfl⟩).modEq_zero_nat
  exact Nat.mod_eq_of_modEq (hadd.trans hzero) hp

/-- UInt64 form used at the generated quotient-loop call site. -/
theorem word3_complement_add_eliminates_mod (p qi lead : UInt64)
    (old out : Word3) (r : Nat)
    (hp : 0 < p.toNat) (hqi : qi.toNat ≤ p.toNat)
    (hold : word3Value old % p.toNat = r)
    (hproduct : (qi.toNat * lead.toNat) % p.toNat = r)
    (hout : word3Value out =
      word3Value old + (p - qi).toNat * lead.toNat) :
    word3Value out % p.toNat = 0 := by
  rw [hout, UInt64.toNat_sub_of_le _ _ hqi]
  exact add_complement_eliminates_mod p.toNat (word3Value old)
    qi.toNat lead.toNat r hp hqi hold hproduct

/-- A modular range observation of the generated loop is an ordinary exact
addition range once the independently proved machine-capacity invariant is
available. -/
theorem addMulRangeRep_exact_of_budget (before after : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (length offset start stop count : Nat) (p c : UInt64)
    (hrange : AddMulRangeRep before after B W3 offset start stop c)
    (hbudget : Word3AccumulationBudget before W3 length p count)
    (hcanonical : CanonicalU64Prefix before B (stop + 1) p)
    (hspan : offset + stop < length)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) :
    ExactAddMulRangeRep before after B W3 offset start stop c := by
  intro k hstart hstop bj accum hreadB hreadW
  rcases hrange k hstart hstop bj accum hreadB hreadW with
    ⟨out, hreadOut, hmod⟩
  have hkLength : offset + k < length := by omega
  have hbase := hbudget (offset + k) accum hkLength hreadW
  have hbj := hcanonical k bj (by omega) hreadB
  have hc' : c.toNat ≤ p.toNat - 1 := by omega
  have hbj' : bj.toNat ≤ p.toNat - 1 := by omega
  have hproduct : c.toNat * bj.toNat ≤ (p.toNat - 1) ^ 2 := by
    rw [pow_two]
    exact Nat.mul_le_mul hc' hbj'
  have hsum : word3Value accum + c.toNat * bj.toNat ≤
      (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          ((p.toNat - 1) + count * (p.toNat - 1) ^ 2) +
            (p.toNat - 1) ^ 2 := Nat.add_le_add hbase hproduct
      _ = (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := by ring
  have hcapacity := lazyAccumulation_word3_budget p (count + 1)
    (p.toNat - 1) hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hrhs : word3Value accum + c.toNat * bj.toNat < limbBase ^ 3 := by
    calc
      word3Value accum + c.toNat * bj.toNat ≤
          (p.toNat - 1) + (count + 1) * (p.toNat - 1) ^ 2 := hsum
      _ < p.toNat * limbBase ^ 2 := hcapacity
      _ < limbBase ^ 3 := by
        calc
          p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
            Nat.mul_lt_mul_of_pos_right hpB (pow_pos (by omega) 2)
          _ = limbBase ^ 2 * limbBase := Nat.mul_comm _ _
          _ = limbBase ^ 3 := by ring
  exact ⟨out, hreadOut, hmod.eq_of_lt_of_lt (word3Value_lt out) hrhs⟩

/-- The top cell written by the actual generated quotient coefficient and
`addMulLoop` is zero modulo `p`.  This combines the generated inverse/multiply
proof with the exact raw-memory range update; no polynomial-level division is
used. -/
theorem generated_quotient_top_cell_eliminates (this : DenseUPolyZp)
    (before after : RawHeap) (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (offset d : Nat) (r lead : UInt64) (old : Word3)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat)
    (hr : r.toNat < this._p.toNat)
    (hleadPos : 0 < lead.toNat) (hlead : lead.toNat < this._p.toNat)
    (hreadLead : before.readU64 B d = .ok lead)
    (hreadOld : before.readWord3 W3 (offset + d) = .ok old)
    (hold : word3Value old % this._p.toNat = r.toNat)
    (hrange :
      let invLc := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
      let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
      ExactAddMulRangeRep before after B W3 offset 0 d (this._p - qi)) :
    ∃ out, after.readWord3 W3 (offset + d) = .ok out ∧
      word3Value out % this._p.toNat = 0 := by
  let invLc := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
  let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
  have hcoeff := quotientCoeff_eliminates_lead this r lead hcfg hprime hr
    hleadPos hlead
  change qi.toNat < this._p.toNat ∧
      (qi.toNat * lead.toNat) % this._p.toNat = r.toNat at hcoeff
  have hrange' : ExactAddMulRangeRep before after B W3 offset 0 d
      (this._p - qi) := by simpa [invLc, qi] using hrange
  rcases hrange' d (by omega) (by omega) lead old hreadLead hreadOld with
    ⟨out, hreadOut, hout⟩
  refine ⟨out, hreadOut, ?_⟩
  exact word3_complement_add_eliminates_mod this._p qi lead old out r.toNat
    hprime.pos (Nat.le_of_lt hcoeff.1) hold hcoeff.2 hout

/-- Full content invariant for the generated multiply/add loop.  B is
unchanged, cells before `j` are unchanged, and every processed cell in
`j..d` is its input accumulator plus the corresponding `c*B[k]` product. -/
theorem addMulLoop_refines (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (htop : i + d < lenW3) (hj : j ≤ d + 1)
    (hregions : W3.region ≠ B.region) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' B (d + 1) ∧
      SameWord3PrefixAt heap heap' W3 i j ∧
      AddMulRangeRep heap heap' B W3 i j d c := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.1 product.2
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    rcases addMulLoop_refines heap1 B W3 lenW3 i d (j + 1) c hB1 hW31
      htop (by omega) hregions with
      ⟨heap2, hloop, hB2, hW32, hlayout2, hsameB2, hprefix2, hrange2⟩
    refine ⟨heap2, hloop, hB2, hW32, ?_, ?_, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k value hk hread
      have hread1 := RawHeap.readU64_writeWord3_region_ne heap heap1
        W3 B (i + j) k accum' value hwrite hread hregions
      exact hsameB2 k value hk hread1
    · intro k value hk hread
      have hkj : k ≠ j := by omega
      have hread1 := RawHeap.readWord3_writeWord3_ne heap heap1 W3
        (i + j) (i + k) accum' value hwrite hread (by omega)
      exact hprefix2 k value (by omega) hread1
    · intro k hlow hhigh bk old hreadBk hreadOld
      by_cases hkj : k = j
      · subst k
        have hcell := addMulCell_refines heap heap1 B W3 j (i + j) c bj
          accum hreadB hreadW hregions (by simpa [product, accum'] using hwrite)
        rcases hcell with ⟨written, _, hreadWritten, _, hmod⟩
        have hbk : bk = bj := Except.ok.inj (hreadBk.symm.trans hreadB)
        have hold : old = accum := Except.ok.inj (hreadOld.symm.trans hreadW)
        subst bk
        subst old
        exact ⟨written, hprefix2 j written (by omega) hreadWritten, hmod⟩
      · have hjk : j + 1 ≤ k := by omega
        have hreadBk1 := RawHeap.readU64_writeWord3_region_ne heap heap1
          W3 B (i + j) k accum' bk hwrite hreadBk hregions
        have hreadOld1 := RawHeap.readWord3_writeWord3_ne heap heap1 W3
          (i + j) (i + k) accum' old hwrite hreadOld (by omega)
        exact hrange2 k hjk hhigh bk old hreadBk1 hreadOld1
  next hnot =>
    refine ⟨heap, rfl, hB, hW3, fun _ _ => Iff.rfl, ?_, ?_, ?_⟩
    · intro k value hk hread
      exact hread
    · intro k value hk hread
      exact hread
    · intro k hlow hhigh
      omega
termination_by d + 1 - j
decreasing_by omega

/-- Exact content theorem for the actual generated inner loop.  Unlike
`addMulLoop_refines`, its conclusion is equality in `Nat`, justified by the
raw capacity proof rather than merely equality modulo the word width. -/
theorem addMulLoop_refines_exact (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j count : Nat) (p c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 p count)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) p)
    (htop : i + d < lenW3) (hj : j ≤ d + 1)
    (hregions : W3.region ≠ B.region)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' B (d + 1) ∧
      SameWord3PrefixAt heap heap' W3 i j ∧
      ExactAddMulRangeRep heap heap' B W3 i j d c := by
  rcases addMulLoop_refines heap B W3 lenW3 i d j c hB hW3 htop hj
    hregions with
    ⟨heap', hrun, hB', hW3', hlayout, hsameB, hprefix, hrange⟩
  have hexact := addMulRangeRep_exact_of_budget heap heap' B W3 lenW3
    i j d count p c hrange hbudget hcanonical htop hp hcount hc
  exact ⟨heap', hrun, hB', hW3', hlayout, hsameB, hprefix, hexact⟩

/-- End-to-end semantic package for one execution of the generated inner
loop.  The returned arrays are uniquely observed from the actual pre/post
raw heaps, and their polynomial relation follows from the machine execution. -/
theorem addMulLoop_refines_polynomial (heap : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenW3 i d count : Nat) (p c : UInt64)
    (divisor : Polynomial (ZMod p.toNat))
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 p count)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) p)
    (hdivisor : SlicePolyRep heap B (d + 1) p.toNat divisor)
    (htop : i + d < lenW3)
    (hregions : W3.region ≠ B.region)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) :
    ∃ heap' beforeValues afterValues,
      addMulLoop heap B W3 i d 0 c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' B (d + 1) ∧
      Word3SliceRep heap W3 lenW3 beforeValues ∧
      Word3SliceRep heap' W3 lenW3 afterValues ∧
      word3ArrayPoly p.toNat afterValues =
        word3ArrayPoly p.toNat beforeValues +
          Polynomial.monomial i (c.toNat : ZMod p.toNat) * divisor := by
  rcases word3SliceRep_exists_unique heap W3 lenW3 hW3 with
    ⟨beforeValues, hbefore, _⟩
  rcases addMulLoop_refines_exact heap B W3 lenW3 i d 0 count p c
      hB hW3 hbudget hcanonical htop (by omega) hregions hp hcount hc with
    ⟨heap', hrun, hB', hW3', hlayout, hsameB, _, hrange⟩
  rcases addMulLoop_preserves_below heap B W3 lenW3 i d 0 c hB hW3
      htop (by omega) with
    ⟨heapBelow, hrunBelow, _, _, _, hbelow⟩
  have hheapBelow : heapBelow = heap' :=
    Except.ok.inj (hrunBelow.symm.trans hrun)
  subst heapBelow
  rcases addMulLoop_preserves_above heap B W3 lenW3 i d 0 c hB hW3
      htop (by omega) with
    ⟨heapAbove, hrunAbove, _, _, _, habove⟩
  have hheapAbove : heapAbove = heap' :=
    Except.ok.inj (hrunAbove.symm.trans hrun)
  subst heapAbove
  rcases word3SliceRep_exists_unique heap' W3 lenW3 hW3' with
    ⟨afterValues, hafter, _⟩
  refine ⟨heap', beforeValues, afterValues, hrun, hB', hW3', hlayout, hsameB,
    hbefore, hafter, ?_⟩
  exact word3ArrayPoly_addMul heap heap' B W3 lenW3 i d p.toNat c
    beforeValues afterValues divisor hbefore hafter hdivisor hbelow habove
    hrange htop

/-- Capacity refinement of the complete generated inner multiply/add loop.
The staged invariant ensures that this whole loop consumes exactly one outer
quotient allowance per affected cell, independently of the divisor degree. -/
theorem addMulLoop_preserves_staged_budget (heap : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenW3 i d j count : Nat) (p c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) p)
    (hstage : Word3StagedBudget heap W3 lenW3 p count i j)
    (htop : i + d < lenW3) (hj : j ≤ d + 1)
    (hregions : W3.region ≠ B.region)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      CanonicalU64Prefix heap' B (d + 1) p ∧
      Word3StagedBudget heap' W3 lenW3 p count i (d + 1) := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with
      ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.1 product.2
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    have hbj : bj.toNat < p.toNat := hcanonical j bj hjB hreadB
    have hstage1 : Word3StagedBudget heap1 W3 lenW3 p count i (j + 1) :=
      writeAddMul_preserves_staged_budget heap heap1 W3 lenW3 count i j
        p c bj accum hW3 hijW hstage hp hcount hc hbj hreadW
        (by simpa [product, accum'] using hwrite)
    have hcanonical1 : CanonicalU64Prefix heap1 B (d + 1) p := by
      intro k value hk hread1
      rcases heap.readU64_of_valid B (d + 1) k hB hk with
        ⟨old, hreadOld⟩
      have hpreserved := RawHeap.readU64_writeWord3_region_ne heap heap1
        W3 B (i + j) k accum' old hwrite hreadOld hregions
      have hvalue : value = old :=
        Except.ok.inj (hread1.symm.trans hpreserved)
      subst value
      exact hcanonical k old hk hreadOld
    rcases addMulLoop_preserves_staged_budget heap1 B W3 lenW3 i d
      (j + 1) count p c hB1 hW31 hcanonical1 hstage1 htop (by omega)
      hregions hp hcount hc with
      ⟨heap2, hloop, hB2, hW32, hlayout2, hcanonical2, hstage2⟩
    refine ⟨heap2, hloop, hB2, hW32, ?_, hcanonical2, hstage2⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    have hjeq : j = d + 1 := by omega
    subst j
    exact ⟨heap, rfl, hB, hW3, fun _ _ => Iff.rfl, hcanonical, hstage⟩
termination_by d + 1 - j
decreasing_by omega

/-- Public outer-step form: a complete generated `addMulLoop` raises the
per-cell accumulation count once, not once per divisor coefficient. -/
theorem addMulLoop_preserves_budget (heap : RawHeap)
    (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenW3 i d count : Nat) (p c : UInt64)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) p)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 p count)
    (htop : i + d < lenW3) (hregions : W3.region ≠ B.region)
    (hp : 1 < p.toNat) (hcount : count + 1 < limbBase)
    (hc : c.toNat < p.toNat) :
    ∃ heap', addMulLoop heap B W3 i d 0 c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      CanonicalU64Prefix heap' B (d + 1) p ∧
      Word3AccumulationBudget heap' W3 lenW3 p (count + 1) := by
  have hstage0 := accumulationBudget_to_staged_zero heap W3 lenW3 p
    count i hbudget
  rcases addMulLoop_preserves_staged_budget heap B W3 lenW3 i d 0 count
    p c hB hW3 hcanonical hstage0 htop (by omega) hregions hp hcount hc with
    ⟨heap', hrun, hB', hW3', hlayout, hcanonical', hstage'⟩
  exact ⟨heap', hrun, hB', hW3', hlayout, hcanonical',
    stagedBudget_to_next_budget heap' W3 lenW3 p count i (d + 1) hstage'⟩

/-- The descending quotient loop cannot fault when Q, B and W3 have the
capacities documented by the C++ raw API.  Its induction variable is exactly
the source `ii`: the successor case processes coefficient `i = ii-1`, then
recurses on the predecessor. -/
theorem quotientLoop_ok (this : DenseUPolyZp) (Q B : RawPtr UInt64)
    (W3 : RawPtr Word3) (qLen d lenW3 : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hii : ii ≤ qLen) (hspan : qLen + d ≤ lenW3) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' := by
  cases ii with
  | zero => exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    split
    next hnonzero =>
      rcases addMulLoop_ok heap1 B W3 lenW3 i d 0 (this._p -
          Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
            (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
              accum.hi accum.mid accum.lo this._p this._ninv this._norm) invLc)
          Q qLen hB1 hW31 hQ1 hiW (by omega) with
        ⟨heap2, hadd, hB2, hW32, hQ2, hlayout2⟩
      simp only [hadd]
      rcases quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap2 i
        hQ2 hB2 hW32 (by omega) hspan with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_⟩
      intro ptr length
      exact (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))
    next hzero =>
      rcases quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap1 i
        hQ1 hB1 hW31 (by omega) hspan with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_⟩
      intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- Full semantic base case of the generated quotient recursion. -/
theorem quotientLoop_zero_refines_polynomial (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 : Nat) (invLc : UInt64) (heap : RawHeap)
    (divisor : Polynomial (ZMod this._p.toNat))
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3) :
    ∃ heap' quotient beforeValues afterValues,
      quotientLoop this Q B W3 d invLc heap 0 = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      SlicePolyRep heap' Q 0 this._p.toNat quotient ∧
      Word3SliceRep heap W3 lenW3 beforeValues ∧
      Word3SliceRep heap' W3 lenW3 afterValues ∧
      word3ArrayPoly this._p.toNat afterValues =
        word3ArrayPoly this._p.toNat beforeValues - quotient * divisor := by
  have hQ0 : heap.ValidU64Slice Q 0 :=
    heap.validU64Slice_mono Q qLen 0 hQ (Nat.zero_le _)
  rcases slicePolyRep_exists_unique heap Q 0 this._p.toNat hQ0 with
    ⟨quotient, hquotient, _⟩
  have hquotientZero := slicePolyRep_zero_length heap Q this._p.toNat
    quotient hquotient
  rcases word3SliceRep_exists_unique heap W3 lenW3 hW3 with
    ⟨values, hvalues, _⟩
  refine ⟨heap, quotient, values, values, rfl, hQ, hB, hW3, hquotient,
    hvalues, hvalues, ?_⟩
  rw [hquotientZero, zero_mul, sub_zero]

/-- The descending generated quotient loop writes only Q indices below its
current `ii`.  All already-produced higher coefficients are preserved across
both the Q write and the optional W3 multiply/add. -/
theorem quotientLoop_preserves_Q_above (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hii : ii ≤ qLen) (hspan : qLen + d ≤ lenW3)
    (hQW : Q.region ≠ W3.region) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Above heap heap' Q ii qLen := by
  cases ii with
  | zero =>
      exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl,
        fun _ _ _ _ hread => hread⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    have hsameWrite : SameU64Above heap heap1 Q (i + 1) qLen := by
      intro k value hlow hhigh hreadOld
      exact RawHeap.readU64_writeU64_ne heap heap1 Q Q i k _ value hwrite
        hreadOld (Or.inr (by omega))
    split
    next hnonzero =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 invLc
      rcases addMulLoop_preserves_u64_region_ne heap1 B Q W3 lenW3 i d 0
          qLen (this._p - qi0) hB1 hW31 hQ1 hiW (by omega)
          (Ne.symm hQW) with
        ⟨heap2, hadd, hB2, hW32, hQ2, hlayout2, hsameAdd⟩
      simp only [qi0, r0] at hadd
      simp only [hadd]
      rcases quotientLoop_preserves_Q_above this Q B W3 qLen d lenW3
          invLc heap2 i hQ2 hB2 hW32 (by omega) hspan hQW with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3, hsameRec⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_, ?_⟩
      · intro ptr length
        exact (hlayout1 ptr length).trans
          ((hlayout2 ptr length).trans (hlayout3 ptr length))
      · intro k value hlow hhigh hreadOld
        exact hsameRec k value (by omega) hhigh
          (hsameAdd k value hhigh
            (hsameWrite k value hlow hhigh hreadOld))
    next hzero =>
      rcases quotientLoop_preserves_Q_above this Q B W3 qLen d lenW3
          invLc heap1 i hQ1 hB1 hW31 (by omega) hspan hQW with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2, hsameRec⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_, ?_⟩
      · intro ptr length
        exact (hlayout1 ptr length).trans (hlayout2 ptr length)
      · intro k value hlow hhigh hreadOld
        exact hsameRec k value (by omega) hhigh
          (hsameWrite k value hlow hhigh hreadOld)

/-- In particular, after the successor case writes `Q[i]`, the recursive
call at state `i` returns with that exact raw coefficient still present. -/
theorem quotientLoop_preserves_Q_at_current (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 i : Nat) (invLc value : UInt64)
    (heap : RawHeap)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hi : i < qLen) (hspan : qLen + d ≤ lenW3)
    (hQW : Q.region ≠ W3.region)
    (hread : heap.readU64 Q i = .ok value) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap i = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      heap'.readU64 Q i = .ok value := by
  rcases quotientLoop_preserves_Q_above this Q B W3 qLen d lenW3 invLc
      heap i hQ hB hW3 (by omega) hspan hQW with
    ⟨heap', hrun, hQ', hB', hW3', hlayout, hsame⟩
  exact ⟨heap', hrun, hQ', hB', hW3', hlayout,
    hsame i value (Nat.le_refl _) hi hread⟩

/-- The generated quotient loop only reads B.  Q writes are separated by the
raw non-aliasing precondition and W3 writes preserve B through the certified
inner-loop theorem. -/
theorem quotientLoop_preserves_B (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hii : ii ≤ qLen) (hspan : qLen + d ≤ lenW3)
    (hQB : Q.region ≠ B.region)
    (hWB : W3.region ≠ B.region) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' B (d + 1) := by
  cases ii with
  | zero =>
      exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl,
        fun _ _ _ hread => hread⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    have hsameWrite : SameU64Prefix heap heap1 B (d + 1) := by
      intro k value hk hreadOld
      exact RawHeap.readU64_writeU64_ne heap heap1 Q B i k _ value hwrite
        hreadOld (Or.inl hQB)
    split
    next hnonzero =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 invLc
      rcases addMulLoop_refines heap1 B W3 lenW3 i d 0
          (this._p - qi0) hB1 hW31 hiW (by omega) hWB with
        ⟨heap2, hadd, hB2, hW32, hlayout2, hsameAdd, _, _⟩
      simp only [qi0, r0] at hadd
      simp only [hadd]
      have hQ2 : heap2.ValidU64Slice Q qLen := (hlayout2 Q qLen).mp hQ1
      rcases quotientLoop_preserves_B this Q B W3 qLen d lenW3 invLc
          heap2 i hQ2 hB2 hW32 (by omega) hspan hQB hWB with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3, hsameRec⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_, ?_⟩
      · intro ptr length
        exact (hlayout1 ptr length).trans
          ((hlayout2 ptr length).trans (hlayout3 ptr length))
      · intro k value hk hreadOld
        exact hsameRec k value hk
          (hsameAdd k value hk (hsameWrite k value hk hreadOld))
    next hzero =>
      rcases quotientLoop_preserves_B this Q B W3 qLen d lenW3 invLc
          heap1 i hQ1 hB1 hW31 (by omega) hspan hQB hWB with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2, hsameRec⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_, ?_⟩
      · intro ptr length
        exact (hlayout1 ptr length).trans (hlayout2 ptr length)
      · intro k value hk hreadOld
        exact hsameRec k value hk (hsameWrite k value hk hreadOld)

/-- The generated descending quotient loop preserves the real lazy
accumulation capacity.  `count + ii = qLen` ties the proof counter to the
source loop state; zero quotient coefficients consume a conservative outer
allowance without performing an add, while nonzero coefficients execute the
certified generated inner loop. -/
theorem quotientLoop_preserves_budget (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 count : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) this._p)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 this._p count)
    (hstate : count + ii = qLen)
    (hspan : qLen + d ≤ lenW3) (hqLen : qLen < limbBase)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hWB : W3.region ≠ B.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      CanonicalU64Prefix heap' B (d + 1) this._p ∧
      Word3AccumulationBudget heap' W3 lenW3 this._p qLen := by
  cases ii with
  | zero =>
      have hcount : count = qLen := by omega
      subst count
      exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl, hcanonical,
        hbudget⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    have hcountLt : count < limbBase := by omega
    have hcountNext : count + 1 < limbBase := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    have hhi : accum.hi.toNat < this._p.toNat :=
      word3_hi_lt_of_accumulation_budget heap W3 lenW3 count (i + d)
        this._p accum hbudget hprime.two_le hcountLt hiW hread
    have hrEq : r.toNat = word3Value accum % this._p.toNat := by
      simpa [r] using lll_mod_preinv_ir_correct_of_configured this
        accum.hi accum.mid accum.lo hcfg hhi
    have hr : r.toNat < this._p.toNat := by
      rw [hrEq]
      exact Nat.mod_lt _ hprime.pos
    have hqiEq : qi.toNat = (r.toNat * invLc.toNat) % this._p.toNat := by
      simpa [qi] using nmod_mul_ir_correct_of_configured this r invLc hcfg hr
    have hqi : qi.toNat < this._p.toNat := by
      rw [hqiEq]
      exact Nat.mod_lt _ hprime.pos
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with
      ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    have hcanonical1 := canonicalPrefix_writeU64_region_ne heap heap1 Q B
      i (d + 1) _ this._p hB hcanonical hQB hwrite
    have hbudget1 := accumulationBudget_writeU64_region_ne heap heap1 Q W3
      i lenW3 count _ this._p hW3 hbudget hQW hwrite
    split
    next hnonzero =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 invLc
      have hqi0 : qi0.toNat < this._p.toNat := by simpa [qi0, r0] using hqi
      have hqi0ne : qi0 ≠ 0 := by simpa [qi0, r0] using hnonzero
      have hqi0pos : 0 < qi0.toNat := by
        have : qi0.toNat ≠ 0 := by
          intro hz
          apply hqi0ne
          exact UInt64.toNat_inj.mp (by simpa using hz)
        omega
      have hc : (this._p - qi0).toNat < this._p.toNat := by
        rw [UInt64.toNat_sub_of_le _ _ (Nat.le_of_lt hqi0)]
        omega
      rcases addMulLoop_preserves_budget heap1 B W3 lenW3 i d count
        this._p (this._p - qi0) hB1 hW31 hcanonical1 hbudget1 hiW hWB
        hprime.two_le hcountNext hc with
        ⟨heap2, hadd, hB2, hW32, hlayout2, hcanonical2, hbudget2⟩
      simp only [qi0, r0] at hadd
      simp only [hadd]
      have hQ2 : heap2.ValidU64Slice Q qLen := (hlayout2 Q qLen).mp hQ1
      rcases quotientLoop_preserves_budget this Q B W3 qLen d lenW3
        (count + 1) invLc heap2 i hQ2 hB2 hW32 hcanonical2 hbudget2
        (by omega) hspan hqLen hQB hQW hWB hcfg hprime with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3, hcanonical3, hbudget3⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_, hcanonical3, hbudget3⟩
      intro ptr length
      exact (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))
    next hzero =>
      have hbudgetNext : Word3AccumulationBudget heap1 W3 lenW3 this._p
          (count + 1) := accumulationBudget_mono heap1 W3 lenW3 this._p
            count (count + 1) hbudget1 (by omega)
      rcases quotientLoop_preserves_budget this Q B W3 qLen d lenW3
        (count + 1) invLc heap1 i hQ1 hB1 hW31 hcanonical1 hbudgetNext
        (by omega) hspan hqLen hQB hQW hWB hcfg hprime with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2, hcanonical2, hbudget2⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_, hcanonical2, hbudget2⟩
      intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- Full polynomial semantics of the generated descending quotient loop.
The proof follows the raw recursion itself: every successor reads the current
W3 cell, writes the computed C++ quotient limb, executes the generated
multiply/add when that limb is nonzero, and then recurses on `i`. -/
theorem quotientLoop_refines_polynomial (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 count : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (divisor : Polynomial (ZMod this._p.toNat))
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) this._p)
    (hdivisor : SlicePolyRep heap B (d + 1) this._p.toNat divisor)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 this._p count)
    (hstate : count + ii = qLen)
    (hspan : qLen + d ≤ lenW3) (hqLen : qLen < limbBase)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hWB : W3.region ≠ B.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    ∃ heap' quotient beforeValues afterValues,
      quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      CanonicalU64Prefix heap' B (d + 1) this._p ∧
      Word3AccumulationBudget heap' W3 lenW3 this._p qLen ∧
      SlicePolyRep heap' Q ii this._p.toNat quotient ∧
      Word3SliceRep heap W3 lenW3 beforeValues ∧
      Word3SliceRep heap' W3 lenW3 afterValues ∧
      word3ArrayPoly this._p.toNat afterValues =
        word3ArrayPoly this._p.toNat beforeValues - quotient * divisor := by
  cases ii with
  | zero =>
      have hcount : count = qLen := by omega
      subst count
      have hQ0 : heap.ValidU64Slice Q 0 :=
        heap.validU64Slice_mono Q qLen 0 hQ (Nat.zero_le _)
      rcases slicePolyRep_exists_unique heap Q 0 this._p.toNat hQ0 with
        ⟨quotient, hquotient, _⟩
      have hquotientZero := slicePolyRep_zero_length heap Q this._p.toNat
        quotient hquotient
      rcases word3SliceRep_exists_unique heap W3 lenW3 hW3 with
        ⟨values, hvalues, _⟩
      refine ⟨heap, quotient, values, values, rfl, hQ, hB, hW3,
        fun _ _ => Iff.rfl, hcanonical, hbudget, hquotient, hvalues,
        hvalues, ?_⟩
      rw [hquotientZero, zero_mul, sub_zero]
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    have hcountLt : count < limbBase := by omega
    have hcountNext : count + 1 < limbBase := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    have hhi : accum.hi.toNat < this._p.toNat :=
      word3_hi_lt_of_accumulation_budget heap W3 lenW3 count (i + d)
        this._p accum hbudget hprime.two_le hcountLt hiW hread
    have hrEq : r.toNat = word3Value accum % this._p.toNat := by
      simpa [r] using lll_mod_preinv_ir_correct_of_configured this
        accum.hi accum.mid accum.lo hcfg hhi
    have hr : r.toNat < this._p.toNat := by
      rw [hrEq]
      exact Nat.mod_lt _ hprime.pos
    have hqiEq : qi.toNat = (r.toNat * invLc.toNat) % this._p.toNat := by
      simpa [qi] using nmod_mul_ir_correct_of_configured this r invLc hcfg hr
    have hqi : qi.toNat < this._p.toNat := by
      rw [hqiEq]
      exact Nat.mod_lt _ hprime.pos
    rcases word3SliceRep_exists_unique heap W3 lenW3 hW3 with
      ⟨beforeValues, hbeforeValues, _⟩
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with
      ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    have hcanonical1 := canonicalPrefix_writeU64_region_ne heap heap1 Q B
      i (d + 1) _ this._p hB hcanonical hQB hwrite
    have hbudget1 := accumulationBudget_writeU64_region_ne heap heap1 Q W3
      i lenW3 count _ this._p hW3 hbudget hQW hwrite
    have hsameB1 : SameU64Prefix heap heap1 B (d + 1) := by
      intro k value hk hreadOld
      exact RawHeap.readU64_writeU64_ne heap heap1 Q B i k _ value hwrite
        hreadOld (Or.inl hQB)
    have hdivisor1 := slicePolyRep_of_same_prefix heap heap1 B (d + 1)
      this._p.toNat divisor hB hB1 hsameB1 hdivisor
    have hreadQ1 := RawHeap.readU64_writeU64_same heap heap1 Q i qi hwrite
    split
    next hnonzero =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 invLc
      have hqi0 : qi0.toNat < this._p.toNat := by simpa [qi0, r0] using hqi
      have hc : (this._p - qi0).toNat < this._p.toNat := by
        rw [UInt64.toNat_sub_of_le _ _ (Nat.le_of_lt hqi0)]
        have hqi0ne : qi0 ≠ 0 := by simpa [qi0, r0] using hnonzero
        have : qi0.toNat ≠ 0 := by
          intro hz
          apply hqi0ne
          exact UInt64.toNat_inj.mp (by simpa using hz)
        omega
      rcases addMulLoop_refines_polynomial heap1 B W3 lenW3 i d count
        this._p (this._p - qi0) divisor hB1 hW31 hbudget1 hcanonical1
        hdivisor1 hiW hWB hprime.two_le hcountNext hc with
        ⟨heap2, writeValues, addValues, hadd, hB2, hW32, hlayout2,
          hsameB2, hwriteValues, haddValues, haddPoly⟩
      simp only [qi0, r0] at hadd
      simp only [hadd]
      have hQ2 : heap2.ValidU64Slice Q qLen := (hlayout2 Q qLen).mp hQ1
      have hcanonical2 : CanonicalU64Prefix heap2 B (d + 1) this._p := by
        intro k value hk hread2
        rcases heap1.readU64_of_valid B (d + 1) k hB1 hk with
          ⟨old, hread1⟩
        have : value = old := Except.ok.inj (hread2.symm.trans
          (hsameB2 k old hk hread1))
        subst value
        exact hcanonical1 k old hk hread1
      have hdivisor2 := slicePolyRep_of_same_prefix heap1 heap2 B (d + 1)
        this._p.toNat divisor hB1 hB2 hsameB2 hdivisor1
      rcases addMulLoop_preserves_budget heap1 B W3 lenW3 i d count
        this._p (this._p - qi0) hB1 hW31 hcanonical1 hbudget1 hiW hWB
        hprime.two_le hcountNext hc with
        ⟨heapBudget, haddBudget, _, _, _, _, hbudget2⟩
      have hheapBudget : heapBudget = heap2 :=
        Except.ok.inj (haddBudget.symm.trans hadd)
      subst heapBudget
      rcases addMulLoop_preserves_u64_region_ne heap1 B Q W3 lenW3 i d 0
        qLen (this._p - qi0) hB1 hW31 hQ1 hiW (by omega)
        (Ne.symm hQW) with
        ⟨heapQ, haddQ, _, _, _, _, hsameQ⟩
      have hheapQ : heapQ = heap2 := Except.ok.inj (haddQ.symm.trans hadd)
      subst heapQ
      have hreadQ2 := hsameQ i qi (by omega) hreadQ1
      rcases quotientLoop_refines_polynomial this Q B W3 qLen d lenW3
        (count + 1) invLc heap2 i divisor hQ2 hB2 hW32 hcanonical2
        hdivisor2 hbudget2 (by omega) hspan hqLen hQB hQW hWB hcfg hprime with
        ⟨heap3, lowerQ, recValues, finalValues, hloop, hQ3, hB3, hW33,
          hlayout3, hcanonical3, hbudget3, hlowerQ, hrecValues,
          hfinalValues, hrecPoly⟩
      rcases quotientLoop_preserves_Q_at_current this Q B W3 qLen d lenW3
        i invLc qi heap2 hQ2 hB2 hW32 hiQ hspan hQW hreadQ2 with
        ⟨heapRead, hloopRead, _, _, _, _, hreadFinal⟩
      have hheapRead : heapRead = heap3 :=
        Except.ok.inj (hloopRead.symm.trans hloop)
      subst heapRead
      rcases quotient_nonzero_successor_finalize heap heap1 heap2 heap3 Q W3
        lenW3 i this._p.toNat this._p qi lowerQ divisor beforeValues
        writeValues writeValues addValues recValues finalValues hW3
        hbeforeValues hwriteValues hwriteValues haddValues hrecValues
        hfinalValues hwrite hQW rfl (Nat.le_of_lt hqi) haddPoly hlowerQ
        (heap3.validU64Slice_mono Q qLen (i + 1) hQ3 (by omega)) hreadFinal
        hrecPoly with ⟨fullQ, hfullQ, halgebra⟩
      refine ⟨heap3, fullQ, beforeValues, finalValues, hloop, hQ3, hB3,
        hW33, ?_, hcanonical3, hbudget3, hfullQ, hbeforeValues,
        hfinalValues, halgebra⟩
      intro ptr length
      exact (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))
    next hzero =>
      have hbudgetNext : Word3AccumulationBudget heap1 W3 lenW3 this._p
          (count + 1) := accumulationBudget_mono heap1 W3 lenW3 this._p
            count (count + 1) hbudget1 (by omega)
      rcases word3SliceRep_exists_unique heap1 W3 lenW3 hW31 with
        ⟨writeValues, hwriteValues, _⟩
      rcases quotientLoop_refines_polynomial this Q B W3 qLen d lenW3
        (count + 1) invLc heap1 i divisor hQ1 hB1 hW31 hcanonical1
        hdivisor1 hbudgetNext (by omega) hspan hqLen hQB hQW hWB hcfg
        hprime with
        ⟨heap2, lowerQ, recValues, finalValues, hloop, hQ2, hB2, hW32,
          hlayout2, hcanonical2, hbudget2, hlowerQ, hrecValues,
          hfinalValues, hrecPoly⟩
      rcases quotientLoop_preserves_Q_at_current this Q B W3 qLen d lenW3
        i invLc qi heap1 hQ1 hB1 hW31 hiQ hspan hQW hreadQ1 with
        ⟨heapRead, hloopRead, _, _, _, _, hreadFinal⟩
      have hheapRead : heapRead = heap2 :=
        Except.ok.inj (hloopRead.symm.trans hloop)
      subst heapRead
      rcases quotient_zero_successor_finalize heap heap1 heap2 Q W3 lenW3
        i this._p.toNat qi lowerQ divisor beforeValues writeValues recValues
        finalValues hW3 hbeforeValues hwriteValues hrecValues hfinalValues
        hwrite hQW (by simpa [qi, r] using hzero) hlowerQ
        (heap2.validU64Slice_mono Q qLen (i + 1) hQ2 (by omega)) hreadFinal
        hrecPoly with ⟨fullQ, hfullQ, halgebra⟩
      refine ⟨heap2, fullQ, beforeValues, finalValues, hloop, hQ2, hB2,
        hW32, ?_, hcanonical2, hbudget2, hfullQ, hbeforeValues,
        hfinalValues, halgebra⟩
      intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- Strengthening of `quotientLoop_preserves_budget` with the descending
leading-cell invariant.  Every processed high W3 cell is proved zero modulo
the source prime, including the generated `qi = 0` branch. -/
theorem quotientLoop_zeroes_high_cells (this : DenseUPolyZp)
    (Q B : RawPtr UInt64) (W3 : RawPtr Word3)
    (qLen d lenW3 count : Nat) (lead invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hcanonical : CanonicalU64Prefix heap B (d + 1) this._p)
    (hreadLead : heap.readU64 B d = .ok lead)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 this._p count)
    (hzero : Word3ZeroModRange heap W3 (ii + d) (qLen + d)
      this._p.toNat)
    (hstate : count + ii = qLen)
    (hspan : qLen + d ≤ lenW3) (hqLen : qLen < limbBase)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hWB : W3.region ≠ B.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat)
    (hleadPos : 0 < lead.toNat) (hlead : lead.toNat < this._p.toNat)
    (hinv : invLc =
      Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      CanonicalU64Prefix heap' B (d + 1) this._p ∧
      Word3AccumulationBudget heap' W3 lenW3 this._p qLen ∧
      Word3ZeroModRange heap' W3 d (qLen + d) this._p.toNat := by
  subst invLc
  cases ii with
  | zero =>
      have hcount : count = qLen := by omega
      subst count
      exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl, hcanonical,
        hbudget, by simpa using hzero⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    have hcountLt : count < limbBase := by omega
    have hcountNext : count + 1 < limbBase := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let inv := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r inv
    have hhi : accum.hi.toNat < this._p.toNat :=
      word3_hi_lt_of_accumulation_budget heap W3 lenW3 count (i + d)
        this._p accum hbudget hprime.two_le hcountLt hiW hread
    have hrEq : r.toNat = word3Value accum % this._p.toNat := by
      simpa [r] using lll_mod_preinv_ir_correct_of_configured this
        accum.hi accum.mid accum.lo hcfg hhi
    have hr : r.toNat < this._p.toNat := by
      rw [hrEq]
      exact Nat.mod_lt _ hprime.pos
    have hcoeff := quotientCoeff_eliminates_lead this r lead hcfg hprime hr
      hleadPos hlead
    change qi.toNat < this._p.toNat ∧
      (qi.toNat * lead.toNat) % this._p.toNat = r.toNat at hcoeff
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with
      ⟨heap1, hwrite⟩
    dsimp [r, inv, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    have hcanonical1 := canonicalPrefix_writeU64_region_ne heap heap1 Q B
      i (d + 1) _ this._p hB hcanonical hQB hwrite
    have hbudget1 := accumulationBudget_writeU64_region_ne heap heap1 Q W3
      i lenW3 count _ this._p hW3 hbudget hQW hwrite
    have hzero1 := zeroModRange_writeU64_region_ne heap heap1 Q W3 i lenW3
      (i + 1 + d) (qLen + d) this._p.toNat _ hW3 hspan hzero hQW hwrite
    have hreadLead1 := RawHeap.readU64_writeU64_ne heap heap1 Q B i d _ lead
      hwrite hreadLead (Or.inl hQB)
    have hreadW1 := RawHeap.readWord3_writeU64_region_ne heap heap1 Q W3
      i (i + d) _ accum hwrite hread hQW
    split
    next hnonzero =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let inv0 := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 inv0
      have hqi0 : qi0.toNat < this._p.toNat := by
        simpa [qi0, inv0, r0] using hcoeff.1
      have hqi0ne : qi0 ≠ 0 := by simpa [qi0, inv0, r0] using hnonzero
      have hqi0pos : 0 < qi0.toNat := by
        have : qi0.toNat ≠ 0 := by
          intro hz
          apply hqi0ne
          exact UInt64.toNat_inj.mp (by simpa using hz)
        omega
      have hc : (this._p - qi0).toNat < this._p.toNat := by
        rw [UInt64.toNat_sub_of_le _ _ (Nat.le_of_lt hqi0)]
        omega
      rcases addMulLoop_preserves_budget heap1 B W3 lenW3 i d count
        this._p (this._p - qi0) hB1 hW31 hcanonical1 hbudget1 hiW hWB
        hprime.two_le hcountNext hc with
        ⟨heap2, hadd, hB2, hW32, hlayout2, hcanonical2, hbudget2⟩
      rcases addMulLoop_refines_exact heap1 B W3 lenW3 i d 0 count
        this._p (this._p - qi0) hB1 hW31 hbudget1 hcanonical1 hiW
        (by omega) hWB hprime.two_le hcountNext hc with
        ⟨heapExact, haddExact, _, _, _, hsameB, _, hrange⟩
      have hheapExact : heapExact = heap2 :=
        Except.ok.inj (haddExact.symm.trans hadd)
      subst heapExact
      rcases generated_quotient_top_cell_eliminates this heap1 heap2 B W3
        i d r0 lead accum hcfg hprime (by simpa [r0] using hr) hleadPos
        hlead hreadLead1 hreadW1 (by simpa [r0] using hrEq.symm) (by
          simpa [qi0, inv0, r0] using hrange) with
        ⟨topOut, hreadTop, htopZero⟩
      rcases addMulLoop_preserves_above heap1 B W3 lenW3 i d 0
        (this._p - qi0) hB1 hW31 hiW (by omega) with
        ⟨heapAbove, haddAbove, _, _, _, hsameAbove⟩
      have hheapAbove : heapAbove = heap2 :=
        Except.ok.inj (haddAbove.symm.trans hadd)
      subst heapAbove
      have hrest : Word3ZeroModRange heap1 W3 (i + d + 1)
          (qLen + d) this._p.toNat := by
        convert hzero1 using 1 <;> omega
      have hsameRestricted : SameWord3Above heap1 heap2 W3 (i + d)
          (qLen + d) := by
        intro k hkTop hkUpper
        exact hsameAbove k hkTop (by omega)
      have hrest2 := zeroModRange_of_same_above heap1 heap2 W3 (i + d)
        (qLen + d) this._p.toNat hsameRestricted hrest
      have hzero2 : Word3ZeroModRange heap2 W3 (i + d) (qLen + d)
          this._p.toNat := zeroModRange_extend_down heap2 W3 (i + d)
            (qLen + d) this._p.toNat (by
              intro value hreadValue
              have : value = topOut :=
                Except.ok.inj (hreadValue.symm.trans hreadTop)
              subst value
              exact htopZero) hrest2
      have hQ2 : heap2.ValidU64Slice Q qLen := (hlayout2 Q qLen).mp hQ1
      have hreadLead2 := hsameB d lead (by omega) hreadLead1
      simp only [qi0, inv0, r0] at hadd
      simp only [hadd]
      rcases quotientLoop_zeroes_high_cells this Q B W3 qLen d lenW3
        (count + 1) lead inv0 heap2 i hQ2 hB2 hW32 hcanonical2
        hreadLead2 hbudget2 hzero2 (by omega) hspan hqLen hQB hQW hWB
        hcfg hprime hleadPos hlead rfl with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3, hcanonical3,
          hbudget3, hzero3⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_, hcanonical3,
        hbudget3, hzero3⟩
      intro ptr length
      exact (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))
    next hzeroBranch =>
      let r0 := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let inv0 := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
      let qi0 := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r0 inv0
      have hqi0zero : qi0 = 0 := by simpa [qi0, inv0, r0] using hzeroBranch
      have hr0zero : r0.toNat = 0 := by
        have hprod0 : (qi0.toNat * lead.toNat) % this._p.toNat =
            r0.toNat := by simpa [qi0, inv0, r0] using hcoeff.2
        rw [hqi0zero] at hprod0
        simpa using hprod0.symm
      have hcell : ∀ value, heap1.readWord3 W3 (i + d) = .ok value →
          word3Value value % this._p.toNat = 0 := by
        intro value hreadValue
        have hvalue : value = accum :=
          Except.ok.inj (hreadValue.symm.trans hreadW1)
        subst value
        rw [← hrEq]
        simpa [r0] using hr0zero
      have hrest : Word3ZeroModRange heap1 W3 (i + d + 1)
          (qLen + d) this._p.toNat := by
        convert hzero1 using 1 <;> omega
      have hzeroNext := zeroModRange_extend_down heap1 W3 (i + d)
        (qLen + d) this._p.toNat hcell hrest
      have hbudgetNext : Word3AccumulationBudget heap1 W3 lenW3 this._p
          (count + 1) := accumulationBudget_mono heap1 W3 lenW3 this._p
            count (count + 1) hbudget1 (by omega)
      rcases quotientLoop_zeroes_high_cells this Q B W3 qLen d lenW3
        (count + 1) lead inv0 heap1 i hQ1 hB1 hW31 hcanonical1
        hreadLead1 hbudgetNext hzeroNext (by omega) hspan hqLen hQB hQW
        hWB hcfg hprime hleadPos hlead rfl with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2, hcanonical2,
          hbudget2, hzero2⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_, hcanonical2,
        hbudget2, hzero2⟩
      intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- The final C++ remainder loop reads exactly W3[0..d) and writes R[0..d).
Both layouts are preserved at each iteration, and `d-i` decreases. -/
theorem remainderLoop_ok (this : DenseUPolyZp) (R : RawPtr UInt64)
    (W3 : RawPtr Word3) (d lenW3 i : Nat) (heap : RawHeap)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hdW : d ≤ lenW3) (hi : i ≤ d) :
    ∃ heap', remainderLoop this R W3 d i heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' := by
  rw [remainderLoop]
  split
  next hlt =>
    have hiW : i < lenW3 := by omega
    rcases heap.readWord3_of_valid W3 lenW3 i hW3 hiW with ⟨accum, hread⟩
    simp only [hread]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid R d i value hR hlt with ⟨heap1, hwrite⟩
    dsimp [value] at hwrite ⊢
    simp only [hwrite]
    have hR1 : heap1.ValidU64Slice R d :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i _ hwrite R d).mp hR
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 R i _ hwrite
    rcases remainderLoop_ok this R W3 d lenW3 (i + 1) heap1
      hR1 hW31 hdW (by omega) with
      ⟨heap2, hloop, hR2, hW32, hlayout2⟩
    refine ⟨heap2, hloop, hR2, hW32, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot => exact ⟨heap, rfl, hR, hW3, fun _ _ => Iff.rfl⟩
termination_by d - i
decreasing_by omega

/-- The generated remainder loop writes only R.  Under the raw non-aliasing
precondition it preserves every W3 cell, not merely the allocation layout. -/
theorem remainderLoop_preserves_W3 (this : DenseUPolyZp)
    (R : RawPtr UInt64) (W3 : RawPtr Word3)
    (d lenW3 i : Nat) (heap : RawHeap)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hdW : d ≤ lenW3) (hi : i ≤ d)
    (hregions : R.region ≠ W3.region) :
    ∃ heap', remainderLoop this R W3 d i heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      SameWord3PrefixAt heap heap' W3 0 lenW3 := by
  rw [remainderLoop]
  split
  next hlt =>
    have hiW : i < lenW3 := by omega
    rcases heap.readWord3_of_valid W3 lenW3 i hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid R d i value hR hlt with
      ⟨heap1, hwrite⟩
    dsimp [value] at hwrite ⊢
    simp only [hwrite]
    have hR1 : heap1.ValidU64Slice R d :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite R d).mp hR
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 R i value hwrite
    have hsame1 : SameWord3PrefixAt heap heap1 W3 0 lenW3 := by
      intro k old hk hreadOld
      simpa using RawHeap.readWord3_writeU64_region_ne heap heap1 R W3 i k
        value old hwrite (by simpa using hreadOld) hregions
    rcases remainderLoop_preserves_W3 this R W3 d lenW3 (i + 1) heap1
      hR1 hW31 hdW (by omega) hregions with
      ⟨heap2, hloop, hR2, hW32, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hR2, hW32, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k old hk hreadOld
      exact hsame2 k old hk (hsame1 k old hk hreadOld)
  next hnot =>
    exact ⟨heap, rfl, hR, hW3, fun _ _ => Iff.rfl,
      fun _ _ _ hread => hread⟩
termination_by d - i
decreasing_by omega

/-- The generated remainder loop also preserves every cell of an arbitrary
UInt64 slice disjoint from R; this transports the quotient representation to
the heap returned by `_poly_divrem`. -/
theorem remainderLoop_preserves_u64_region_ne (this : DenseUPolyZp)
    (R other : RawPtr UInt64) (W3 : RawPtr Word3)
    (d lenW3 i otherLen : Nat) (heap : RawHeap)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hOther : heap.ValidU64Slice other otherLen)
    (hdW : d ≤ lenW3) (hi : i ≤ d)
    (hregions : R.region ≠ other.region) :
    ∃ heap', remainderLoop this R W3 d i heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      heap'.ValidU64Slice other otherLen ∧ RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' other otherLen := by
  rw [remainderLoop]
  split
  next hlt =>
    have hiW : i < lenW3 := by omega
    rcases heap.readWord3_of_valid W3 lenW3 i hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid R d i value hR hlt with
      ⟨heap1, hwrite⟩
    dsimp [value] at hwrite ⊢
    simp only [hwrite]
    have hR1 : heap1.ValidU64Slice R d :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite R d).mp hR
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite
        other otherLen).mp hOther
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 R i value hwrite
    have hsame1 : SameU64Prefix heap heap1 other otherLen := by
      intro k old hk hreadOld
      exact RawHeap.readU64_writeU64_ne heap heap1 R other i k value old
        hwrite hreadOld (Or.inl hregions)
    rcases remainderLoop_preserves_u64_region_ne this R other W3 d lenW3
      (i + 1) otherLen heap1 hR1 hW31 hOther1 hdW (by omega) hregions with
      ⟨heap2, hloop, hR2, hW32, hOther2, hlayout2, hsame2⟩
    refine ⟨heap2, hloop, hR2, hW32, hOther2, ?_, ?_⟩
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
    · intro k old hk hreadOld
      exact hsame2 k old hk (hsame1 k old hk hreadOld)
  next hnot =>
    exact ⟨heap, rfl, hR, hW3, hOther, fun _ _ => Iff.rfl,
      fun _ _ _ hread => hread⟩
termination_by d - i
decreasing_by omega

/-- The completed remainder prefix contains exactly the outputs of the
generated three-limb reduction on the corresponding W3 cells. -/
def RemainderPrefix (this : DenseUPolyZp) (heap : RawHeap)
    (R : RawPtr UInt64) (W3 : RawPtr Word3) (upto : Nat) : Prop :=
  ∀ j, j < upto → ∃ accum : Word3,
    heap.readWord3 W3 j = .ok accum ∧
    heap.readU64 R j = .ok
      (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm)

def RemainderModPrefix (this : DenseUPolyZp) (heap : RawHeap)
    (R : RawPtr UInt64) (W3 : RawPtr Word3) (upto : Nat) : Prop :=
  ∀ j, j < upto → ∃ accum value,
    heap.readWord3 W3 j = .ok accum ∧
    heap.readU64 R j = .ok value ∧
    value.toNat = word3Value accum % this._p.toNat

theorem remainderPrefix_to_mod (this : DenseUPolyZp) (heap : RawHeap)
    (R : RawPtr UInt64) (W3 : RawPtr Word3) (upto : Nat)
    (hprefix : RemainderPrefix this heap R W3 upto)
    (hn : this._norm.toNat < 64)
    (hpn : (this._p <<< this._norm).toNat =
      this._p.toNat * 2 ^ this._norm.toNat)
    (hpnB : this._p.toNat * 2 ^ this._norm.toNat < limbBase)
    (hnorm : limbBase ≤ 2 * (this._p <<< this._norm).toNat)
    (hmul : (limbBase + this._ninv.toNat) *
      (this._p <<< this._norm).toNat < limbBase ^ 2)
    (hlower : limbBase ^ 2 ≤
      (limbBase + this._ninv.toNat + 1) *
        (this._p <<< this._norm).toNat)
    (hhi : ∀ j accum, j < upto → heap.readWord3 W3 j = .ok accum →
      accum.hi.toNat < this._p.toNat) :
    RemainderModPrefix this heap R W3 upto := by
  intro j hj
  rcases hprefix j hj with ⟨accum, hreadW, hreadR⟩
  let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
    accum.hi accum.mid accum.lo this._p this._ninv this._norm
  refine ⟨accum, value, hreadW, ?_, ?_⟩
  · simpa [value] using hreadR
  · exact lll_mod_preinv_ir_correct accum.hi accum.mid accum.lo
      this._p this._ninv this._norm hn (hhi j accum hj hreadW)
      hpn hpnB hnorm hmul hlower

theorem remainderPrefix_to_mod_of_configured (this : DenseUPolyZp)
    (heap : RawHeap) (R : RawPtr UInt64) (W3 : RawPtr Word3) (upto : Nat)
    (hprefix : RemainderPrefix this heap R W3 upto)
    (hcfg : DensePreinvConfigured this)
    (hhi : ∀ j accum, j < upto → heap.readWord3 W3 j = .ok accum →
      accum.hi.toNat < this._p.toNat) :
    RemainderModPrefix this heap R W3 upto := by
  rcases densePreinvConfigured_conditions this hcfg with
    ⟨hn, hpn, hpnB, hnorm, hmul, hlower⟩
  exact remainderPrefix_to_mod this heap R W3 upto hprefix hn hpn hpnB
    hnorm hmul hlower hhi

/-- Consume the quotient-loop capacity invariant directly when interpreting
the generated remainder writes. -/
theorem remainderPrefix_to_mod_of_budget (this : DenseUPolyZp)
    (heap : RawHeap) (R : RawPtr UInt64) (W3 : RawPtr Word3)
    (upto length count : Nat)
    (hprefix : RemainderPrefix this heap R W3 upto)
    (hbudget : Word3AccumulationBudget heap W3 length this._p count)
    (hupto : upto ≤ length) (hcount : count < limbBase)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    RemainderModPrefix this heap R W3 upto := by
  apply remainderPrefix_to_mod_of_configured this heap R W3 upto hprefix hcfg
  intro j accum hj hread
  exact word3_hi_lt_of_accumulation_budget heap W3 length count j this._p
    accum hbudget hprime.two_le hcount (by omega) hread

/-- A completed raw remainder prefix is exactly the polynomial represented by
the W3 accumulator when all cells above the physical remainder buffer are
zero modulo the source prime. -/
theorem remainderModPrefix_poly_eq (this : DenseUPolyZp) (heap : RawHeap)
    (R : RawPtr UInt64) (W3 : RawPtr Word3) (d length : Nat)
    (values : Array Word3) (remainder : Polynomial (ZMod this._p.toNat))
    (hprefix : RemainderModPrefix this heap R W3 d)
    (hvalues : Word3SliceRep heap W3 length values)
    (hrep : SlicePolyRep heap R d this._p.toNat remainder)
    (hd : d ≤ length)
    (hzero : Word3ZeroModRange heap W3 d length this._p.toNat) :
    remainder = word3ArrayPoly this._p.toNat values := by
  ext degree
  rw [coeff_word3ArrayPoly]
  by_cases hdegree : degree < d
  · have hdegreeValues : degree < values.size := by
      simpa [hvalues.1] using (show degree < length by omega)
    rw [dif_pos hdegreeValues]
    rcases slicePolyRep_coeff heap R d this._p.toNat remainder hrep degree
        hdegree with ⟨value, hreadR, hcoeff⟩
    rcases hprefix degree hdegree with
      ⟨accum, reduced, hreadW, hreadReduced, hreduced⟩
    have hvalue : value = reduced :=
      Except.ok.inj (hreadR.symm.trans hreadReduced)
    have haccum : values[degree] = accum :=
      Except.ok.inj ((hvalues.2 degree hdegreeValues).symm.trans hreadW)
    rw [hcoeff, hvalue, haccum, hreduced]
    exact ZMod.natCast_mod _ _
  · have hcoeffZero := slicePolyRep_coeff_zero_of_length_le heap R d
        this._p.toNat remainder hrep degree (by omega)
    rw [hcoeffZero]
    by_cases hlength : degree < length
    · have hdegreeValues : degree < values.size := by
        simpa [hvalues.1] using hlength
      rw [dif_pos hdegreeValues]
      have hread := hvalues.2 degree hdegreeValues
      have hz := hzero degree values[degree] (by omega) hlength hread
      rw [← ZMod.natCast_mod (word3Value values[degree]) this._p.toNat, hz]
      simp
    · have hdegreeValues : ¬degree < values.size := by
        simpa [hvalues.1] using hlength
      rw [dif_neg hdegreeValues]

/-- Content-level refinement of the generated final remainder loop. -/
theorem remainderLoop_refines (this : DenseUPolyZp) (R : RawPtr UInt64)
    (W3 : RawPtr Word3) (d lenW3 i count : Nat) (heap : RawHeap)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hdW : d ≤ lenW3) (hi : i ≤ d)
    (hregions : R.region ≠ W3.region)
    (hprefix : RemainderPrefix this heap R W3 i)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 this._p count) :
    ∃ heap', remainderLoop this R W3 d i heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      RemainderPrefix this heap' R W3 d ∧
      Word3AccumulationBudget heap' W3 lenW3 this._p count := by
  rw [remainderLoop]
  split
  next hlt =>
    have hiW : i < lenW3 := by omega
    rcases heap.readWord3_of_valid W3 lenW3 i hW3 hiW with ⟨accum, hread⟩
    simp only [hread]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid R d i value hR hlt with ⟨heap1, hwrite⟩
    dsimp [value] at hwrite ⊢
    simp only [hwrite]
    have hR1 : heap1.ValidU64Slice R d :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite R d).mp hR
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i value hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 R i value hwrite
    have hbudget1 := accumulationBudget_writeU64_region_ne heap heap1 R W3
      i lenW3 count value this._p hW3 hbudget hregions hwrite
    have hprefix1 : RemainderPrefix this heap1 R W3 (i + 1) := by
      intro j hj
      by_cases hji : j = i
      · subst j
        have hreadW := RawHeap.readWord3_writeU64_region_ne heap heap1
          R W3 i i value accum hwrite hread hregions
        have hreadR := RawHeap.readU64_writeU64_same heap heap1 R i value hwrite
        exact ⟨accum, hreadW, hreadR⟩
      · have hjlt : j < i := by omega
        rcases hprefix j hjlt with ⟨old, hreadW, hreadR⟩
        have hreadW1 := RawHeap.readWord3_writeU64_region_ne heap heap1
          R W3 i j value old hwrite hreadW hregions
        have hreadR1 := RawHeap.readU64_writeU64_ne heap heap1 R R i j
          value _ hwrite hreadR (Or.inr (by omega))
        exact ⟨old, hreadW1, hreadR1⟩
    rcases remainderLoop_refines this R W3 d lenW3 (i + 1) count heap1
      hR1 hW31 hdW (by omega) hregions hprefix1 hbudget1 with
      ⟨heap2, hloop, hR2, hW32, hlayout2, hfull, hbudget2⟩
    refine ⟨heap2, hloop, hR2, hW32, ?_, hfull, hbudget2⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    have hieq : i = d := by omega
    subst i
    exact ⟨heap, rfl, hR, hW3, fun _ _ => Iff.rfl, hprefix, hbudget⟩
termination_by d - i
decreasing_by omega

/-- End-to-end polynomial semantics of the actual generated remainder loop.
The returned polynomial is reconstructed from the raw R buffer; its equality
with W3 follows from the generated three-limb reductions and the quotient
loop's zero-high-cell invariant. -/
theorem remainderLoop_refines_polynomial (this : DenseUPolyZp)
    (R : RawPtr UInt64) (W3 : RawPtr Word3)
    (d lenW3 count : Nat) (heap : RawHeap) (values : Array Word3)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hvalues : Word3SliceRep heap W3 lenW3 values)
    (hdW : d ≤ lenW3) (hregions : R.region ≠ W3.region)
    (hbudget : Word3AccumulationBudget heap W3 lenW3 this._p count)
    (hcount : count < limbBase)
    (hzero : Word3ZeroModRange heap W3 d lenW3 this._p.toNat)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    ∃ heap' remainder,
      remainderLoop this R W3 d 0 heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' ∧
      Word3AccumulationBudget heap' W3 lenW3 this._p count ∧
      SlicePolyRep heap' R d this._p.toNat remainder ∧
      Word3SliceRep heap' W3 lenW3 values ∧
      remainder = word3ArrayPoly this._p.toNat values := by
  have hempty : RemainderPrefix this heap R W3 0 := by
    intro _ h
    omega
  rcases remainderLoop_refines this R W3 d lenW3 0 count heap hR hW3
    hdW (Nat.zero_le _) hregions hempty hbudget with
    ⟨heap', hrun, hR', hW3', hlayout, hprefix, hbudget'⟩
  rcases remainderLoop_preserves_W3 this R W3 d lenW3 0 heap hR hW3
    hdW (Nat.zero_le _) hregions with
    ⟨heapSame, hrunSame, _, _, _, hsame⟩
  have hheapSame : heapSame = heap' :=
    Except.ok.inj (hrunSame.symm.trans hrun)
  subst heapSame
  have hvalues' : Word3SliceRep heap' W3 lenW3 values := by
    refine ⟨hvalues.1, ?_⟩
    intro k hk
    have hkLen : k < lenW3 := by simpa [hvalues.1] using hk
    simpa using hsame k values[k] hkLen (by simpa using hvalues.2 k hk)
  have hzero' : Word3ZeroModRange heap' W3 d lenW3 this._p.toNat := by
    intro k value hdk hkLen hreadFinal
    have hkValues : k < values.size := by simpa [hvalues.1] using hkLen
    have hreadOld := hvalues.2 k hkValues
    have hreadPreserved := hsame k values[k] hkLen (by simpa using hreadOld)
    have hreadPreserved' : heap'.readWord3 W3 k = .ok values[k] := by
      simpa using hreadPreserved
    have hvalue : value = values[k] :=
      Except.ok.inj (hreadFinal.symm.trans hreadPreserved')
    subst value
    exact hzero k values[k] hdk hkLen hreadOld
  have hmod := remainderPrefix_to_mod_of_budget this heap' R W3 d lenW3
    count hprefix hbudget' hdW hcount hcfg hprime
  rcases slicePolyRep_exists_unique heap' R d this._p.toNat hR' with
    ⟨remainder, hrep, _⟩
  refine ⟨heap', remainder, hrun, hR', hW3', hlayout, hbudget', hrep,
    hvalues', ?_⟩
  exact remainderModPrefix_poly_eq this heap' R W3 d lenW3 values
    remainder hmod hvalues' hrep hdW hzero'

/-- Genuine semantic refinement of the long branch of generated C++
`_poly_divrem`.  The theorem follows the emitted initialization, quotient,
remainder, and normalization calls and reconstructs both result polynomials
from the final raw heap. -/
theorem polyDivrem_long_refines (this : DenseUPolyZp)
    (Q R A B : RawPtr UInt64) (lenA d : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlong : d < lenA)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B (d + 1))
    (hQ : heap.ValidU64Slice Q (lenA - d))
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hcanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hcanonicalB : CanonicalU64Prefix heap B (d + 1) this._p)
    (hdividend : SlicePolyRep heap A lenA this._p.toNat dividend)
    (hdivisor : SlicePolyRep heap B (d + 1) this._p.toNat divisor)
    (hnormB : heap.normaliseU64 B (d + 1) = .ok (d + 1))
    (hqLen : lenA - d < limbBase)
    (hWA : W3.region ≠ A.region) (hWB : W3.region ≠ B.region)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hRW : R.region ≠ W3.region) (hRQ : R.region ≠ Q.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    ∃ heap' lenQ lenR quotient remainder,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B (d + 1) W3 heap =
        .ok (heap', lenQ, lenR) ∧
      SlicePolyRep heap' Q lenQ this._p.toNat quotient ∧
      SlicePolyRep heap' R lenR this._p.toNat remainder ∧
      dividend = quotient * divisor + remainder ∧
      (remainder = 0 ∨ remainder.natDegree < d) ∧
      lenQ ≤ lenA - d ∧ lenR ≤ d := by
  rcases heap.readU64_of_valid B (d + 1) d hB (by omega) with
    ⟨lead, hreadLead⟩
  have hleadLt : lead.toNat < this._p.toNat :=
    hcanonicalB d lead (by omega) hreadLead
  have hleadNe : lead ≠ 0 := by
    rcases normaliseU64_spec heap B (d + 1) hB with
      ⟨result, hnorm, _, _, hlast⟩
    have hresult : result = d + 1 := Except.ok.inj (hnorm.symm.trans hnormB)
    subst result
    rcases hlast with hzero | ⟨last, hreadLast, hlastNe⟩
    · omega
    · have hlastEq : last = lead :=
        have hreadLast' : heap.readU64 B d = .ok last := by
          simpa using hreadLast
        Except.ok.inj (hreadLast'.symm.trans hreadLead)
      simpa [hlastEq] using hlastNe
  have hleadPos : 0 < lead.toNat := by
    have : lead.toNat ≠ 0 := by
      intro hz
      apply hleadNe
      exact UInt64.toNat_inj.mp (by simpa using hz)
    omega
  have hempty : InitW3Prefix heap A W3 0 := by
    intro _ h
    omega
  rcases initW3Loop_refines heap A W3 lenA 0 hA hW3 (Nat.zero_le _)
    hWA hempty with
    ⟨heap1, hinit, hA1, hW31, hlayout1, hprefix1⟩
  rcases initW3Loop_preserves_u64_region_ne heap A B W3 lenA 0 (d + 1)
    hA hW3 hB (Nat.zero_le _) hWB with
    ⟨heapB, hinitB, _, _, hB1, _, hsameB⟩
  have hheapB : heapB = heap1 := Except.ok.inj (hinitB.symm.trans hinit)
  subst heapB
  rcases initW3Loop_preserves_u64_region_ne heap A A W3 lenA 0 lenA
    hA hW3 hA (Nat.zero_le _) hWA with
    ⟨heapA, hinitA, _, _, _, _, hsameA⟩
  have hheapA : heapA = heap1 := Except.ok.inj (hinitA.symm.trans hinit)
  subst heapA
  have hcanonicalA1 : CanonicalU64Prefix heap1 A lenA this._p := by
    intro k value hk hread1
    rcases heap.readU64_of_valid A lenA k hA hk with ⟨old, hread0⟩
    have hvalue : value = old :=
      Except.ok.inj (hread1.symm.trans (hsameA k old hk hread0))
    subst value
    exact hcanonicalA k old hk hread0
  have hcanonicalB1 : CanonicalU64Prefix heap1 B (d + 1) this._p := by
    intro k value hk hread1
    rcases heap.readU64_of_valid B (d + 1) k hB hk with ⟨old, hread0⟩
    have hvalue : value = old :=
      Except.ok.inj (hread1.symm.trans (hsameB k old hk hread0))
    subst value
    exact hcanonicalB k old hk hread0
  have hdividend1 := slicePolyRep_of_same_prefix heap heap1 A lenA
    this._p.toNat dividend hA hA1 hsameA hdividend
  have hdivisor1 := slicePolyRep_of_same_prefix heap heap1 B (d + 1)
    this._p.toNat divisor hB hB1 hsameB hdivisor
  have hQ1 : heap1.ValidU64Slice Q (lenA - d) :=
    (hlayout1 Q (lenA - d)).mp hQ
  have hR1 : heap1.ValidU64Slice R d := (hlayout1 R d).mp hR
  have hbudget1 := initW3Prefix_budget heap1 A W3 lenA this._p
    hcanonicalA1 hprefix1
  rcases word3SliceRep_exists_unique heap1 W3 lenA hW31 with
    ⟨initialValues, hinitialValues, _⟩
  have hinitialPoly := initW3Prefix_word3ArrayPoly heap1 A W3 lenA
    this._p.toNat initialValues dividend hprefix1 hinitialValues hdividend1
  have hreadLead1 := hsameB d lead (by omega) hreadLead
  let invLc := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
  rcases quotientLoop_refines_polynomial this Q B W3 (lenA - d) d lenA 0
    invLc heap1 (lenA - d) divisor hQ1 hB1 hW31 hcanonicalB1 hdivisor1
    hbudget1 (by omega) (by omega) hqLen hQB hQW hWB hcfg hprime with
    ⟨heap2, quotient, beforeValues, quotientValues, hquot, hQ2, hB2,
      hW32, hlayout2, hcanonicalB2, hbudget2, hquotient,
      hbeforeValues, hquotientValues, hquotientPoly⟩
  have hbeforeEq : beforeValues = initialValues :=
    word3SliceRep_eq heap1 W3 lenA beforeValues initialValues
      hbeforeValues hinitialValues
  have hzeroStart : Word3ZeroModRange heap1 W3
      ((lenA - d) + d) ((lenA - d) + d) this._p.toNat := by
    intro _ _ hlow hhigh _
    omega
  rcases quotientLoop_zeroes_high_cells this Q B W3 (lenA - d) d lenA 0
    lead invLc heap1 (lenA - d) hQ1 hB1 hW31 hcanonicalB1 hreadLead1
    hbudget1 hzeroStart (by omega) (by omega) hqLen hQB hQW hWB hcfg
    hprime hleadPos hleadLt rfl with
    ⟨heapZero, hquotZero, _, _, _, _, _, _, hzero2⟩
  have hheapZero : heapZero = heap2 :=
    Except.ok.inj (hquotZero.symm.trans hquot)
  subst heapZero
  have hzero2' : Word3ZeroModRange heap2 W3 d lenA this._p.toNat := by
    convert hzero2 using 1 <;> omega
  have hR2 : heap2.ValidU64Slice R d := (hlayout2 R d).mp hR1
  rcases remainderLoop_refines_polynomial this R W3 d lenA (lenA - d)
    heap2 quotientValues hR2 hW32 hquotientValues (by omega) hRW hbudget2
    hqLen hzero2' hcfg hprime with
    ⟨heap3, remainder, hrem, hR3, hW33, hlayout3, _, hremainder,
      _, hremainderPoly⟩
  rcases remainderLoop_preserves_u64_region_ne this R Q W3 d lenA 0
    (lenA - d) heap2 hR2 hW32 hQ2 (by omega) (Nat.zero_le _) hRQ with
    ⟨heapQ, hremQ, _, _, hQ3, _, hsameQ⟩
  have hheapQ : heapQ = heap3 := Except.ok.inj (hremQ.symm.trans hrem)
  subst heapQ
  have hquotient3 := slicePolyRep_of_same_prefix heap2 heap3 Q (lenA - d)
    this._p.toNat quotient hQ2 hQ3 hsameQ hquotient
  rcases normaliseU64_ok heap3 Q (lenA - d) hQ3 with
    ⟨lenQ, hnormQ, hlenQ⟩
  rcases normaliseU64_ok heap3 R d hR3 with ⟨lenR, hnormR, hlenR⟩
  have halgebra : dividend = quotient * divisor + remainder := by
    rw [hremainderPoly, hquotientPoly, hbeforeEq, hinitialPoly]
    ring
  have hdegree := normaliseU64_poly_degree_lt_length heap3 R d
    this._p.toNat lenR remainder hR3 hremainder hnormR
  have hquotientNorm := slicePolyRep_of_normaliseU64 heap3 Q (lenA - d)
    this._p.toNat lenQ quotient hQ3 hquotient3 hnormQ
  have hremainderNorm := slicePolyRep_of_normaliseU64 heap3 R d
    this._p.toNat lenR remainder hR3 hremainder hnormR
  refine ⟨heap3, lenQ, lenR, quotient, remainder, ?_, hquotientNorm,
    hremainderNorm, halgebra, hdegree, hlenQ, hlenR⟩
  simp [dense_upoly_zp__poly_divrem_ir, hlong, hreadLead, hinit, hquot,
    hrem, hnormQ, hnormR, invLc]

/-- Complete content semantics of the C++ short-division branch: the
generated function returns quotient length zero and copies the represented
dividend, unchanged, into the remainder buffer. -/
theorem polyDivrem_short_refines (this : DenseUPolyZp)
    (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (p : Nat) (dividend : Polynomial (ZMod p))
    (hlenB : 0 < lenB) (hshort : lenA < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hR : heap.ValidU64Slice R lenA)
    (hregions : R.region ≠ A.region)
    (hrep : SlicePolyRep heap A lenA p dividend) :
    ∃ heap',
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', 0, lenA) ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' R lenA p dividend := by
  cases lenB with
  | zero => omega
  | succ d =>
    have hbranch : lenA < d + 1 := by omega
    rcases copyU64_slicePolyRep heap R A lenA p dividend hR hA hregions
      hrep with ⟨heap', hcopy, hlayout, hrep'⟩
    refine ⟨heap', ?_, hlayout, hrep'⟩
    simp [dense_upoly_zp__poly_divrem_ir, hbranch, hcopy]

/-- The short branch also satisfies the mathematical remainder-degree side
condition when the caller-provided C++ lengths are normalized canonical
residue lengths. -/
theorem polyDivrem_short_remainder_degree (heap : RawHeap)
    (A B : RawPtr UInt64) (lenA lenB p : Nat)
    (dividend divisor : Polynomial (ZMod p))
    (hshort : lenA < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hrepA : SlicePolyRep heap A lenA p dividend)
    (hrepB : SlicePolyRep heap B lenB p divisor)
    (hreducedA : ∀ i value, i < lenA → heap.readU64 A i = .ok value →
      value.toNat < p)
    (hreducedB : ∀ i value, i < lenB → heap.readU64 B i = .ok value →
      value.toNat < p)
    (hnormA : heap.normaliseU64 A lenA = .ok lenA)
    (hnormB : heap.normaliseU64 B lenB = .ok lenB) :
    dividend = 0 ∨ dividend.natDegree < divisor.natDegree := by
  by_cases hlenA : lenA = 0
  · left
    subst lenA
    exact normaliseU64_poly_eq_zero heap A 0 p dividend hA hrepA hnormA
  · right
    have hlenB : lenB ≠ 0 := by omega
    rw [normaliseU64_poly_natDegree_eq heap A lenA p lenA dividend
        hA hrepA hreducedA hnormA hlenA,
      normaliseU64_poly_natDegree_eq heap B lenB p lenB divisor
        hB hrepB hreducedB hnormB hlenB]
    omega

/-- Unified strict refinement of generated C++ `_poly_divrem`.  Both source
branches return raw slices at their actual normalized lengths, satisfy the
division identity, and produce a remainder strictly smaller than the
normalized divisor. -/
theorem polyDivrem_refines (this : DenseUPolyZp)
    (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlenB : 0 < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hcanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hcanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hdividend : SlicePolyRep heap A lenA this._p.toNat dividend)
    (hdivisor : SlicePolyRep heap B lenB this._p.toNat divisor)
    (hnormA : heap.normaliseU64 A lenA = .ok lenA)
    (hnormB : heap.normaliseU64 B lenB = .ok lenB)
    (hqCapacity : lenA - (lenB - 1) < limbBase)
    (hRA : R.region ≠ A.region)
    (hWA : W3.region ≠ A.region) (hWB : W3.region ≠ B.region)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hRW : R.region ≠ W3.region) (hRQ : R.region ≠ Q.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat) :
    ∃ heap' lenQ lenR quotient remainder,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      SlicePolyRep heap' Q lenQ this._p.toNat quotient ∧
      SlicePolyRep heap' R lenR this._p.toNat remainder ∧
      dividend = quotient * divisor + remainder ∧
      (remainder = 0 ∨ remainder.natDegree < divisor.natDegree) ∧
      lenQ ≤ lenA - (lenB - 1) ∧ lenR < lenB := by
  cases lenB with
  | zero => omega
  | succ d =>
    by_cases hshort : lenA < d + 1
    · have hlenAd : lenA ≤ d := by omega
      have hRfull : heap.ValidU64Slice R lenA := by
        simpa [Nat.min_eq_left hlenAd] using hR
      rcases polyDivrem_short_refines this Q R A B lenA (d + 1) W3 heap
        this._p.toNat dividend (by omega) hshort hA hRfull hRA hdividend with
        ⟨heap', hrun, hlayout, hremainder⟩
      have hQ' : heap'.ValidU64Slice Q (lenA - d) :=
        (hlayout Q (lenA - d)).mp (by simpa using hQ)
      have hQ0 := heap'.validU64Slice_mono Q (lenA - d) 0 hQ'
        (Nat.zero_le _)
      rcases slicePolyRep_exists_unique heap' Q 0 this._p.toNat hQ0 with
        ⟨quotient, hquotient, _⟩
      have hquotientZero := slicePolyRep_zero_length heap' Q this._p.toNat
        quotient hquotient
      subst quotient
      have hdegree := polyDivrem_short_remainder_degree heap A B lenA
        (d + 1) this._p.toNat dividend divisor hshort hA hB hdividend
        hdivisor hcanonicalA hcanonicalB hnormA hnormB
      refine ⟨heap', 0, lenA, 0, dividend, hrun, hquotient, hremainder,
        ?_, hdegree, by omega, by omega⟩
      simp
    · have hlong : d < lenA := by omega
      have hdleA : d ≤ lenA := Nat.le_of_lt hlong
      have hRfull : heap.ValidU64Slice R d := by
        simpa [Nat.min_eq_right hdleA] using hR
      rcases polyDivrem_long_refines this Q R A B lenA d W3 heap dividend
        divisor hlong hA hB (by simpa using hQ) hRfull hW3 hcanonicalA
        hcanonicalB hdividend hdivisor hnormB (by simpa using hqCapacity)
        hWA hWB hQB hQW hRW hRQ hcfg hprime with
        ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient,
          hremainder, halgebra, hdegree, hlenQ, hlenR⟩
      have hdivisorDegree : divisor.natDegree = d := by
        simpa using normaliseU64_poly_natDegree_eq heap B (d + 1)
          this._p.toNat (d + 1) divisor hB hdivisor hcanonicalB hnormB
          (by omega)
      have hdegree' : remainder = 0 ∨
          remainder.natDegree < divisor.natDegree := by
        rcases hdegree with hzero | hlt
        · exact Or.inl hzero
        · exact Or.inr (by simpa [hdivisorDegree] using hlt)
      exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient,
        hremainder, halgebra, hdegree', by simpa using hlenQ, by omega⟩

/-- Under exactly the capacities documented on the C++ raw API,
`_poly_divrem` cannot take `RawFault`.  This theorem is only the termination
and memory-safety bridge; the quotient/remainder algebraic invariant is the
next refinement obligation. -/
theorem polyDivrem_ok (this : DenseUPolyZp) (Q R A B : RawPtr UInt64)
    (lenA lenB : Nat) (W3 : RawPtr Word3) (heap : RawHeap)
    (hlenB : 0 < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA) :
    ∃ heap' lenQ lenR,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      lenQ ≤ lenA - (lenB - 1) ∧
      lenR ≤ Nat.min lenA (lenB - 1) := by
  cases lenB with
  | zero => omega
  | succ d =>
    simp only [dense_upoly_zp__poly_divrem_ir]
    split
    next hshort =>
      have hlenAd : lenA ≤ d := by omega
      have hRfull : heap.ValidU64Slice R lenA := by
        simpa [Nat.min_eq_left hlenAd] using hR
      rcases copyU64_ok heap R A lenA hRfull hA with
        ⟨heap1, hcopy, hlayout⟩
      simp only [hcopy]
      exact ⟨heap1, 0, lenA, rfl, by omega, by simpa [Nat.min_eq_left hlenAd]⟩
    next hlong =>
      have hdA : d < lenA := by omega
      have hdleA : d ≤ lenA := Nat.le_of_lt hdA
      have hRfull : heap.ValidU64Slice R d := by
        simpa [Nat.min_eq_right hdleA] using hR
      rcases heap.readU64_of_valid B (d + 1) d hB (by omega) with
        ⟨lead, hlead⟩
      simp only [hlead]
      let invLc := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
      rcases initW3Loop_ok heap A W3 lenA 0 hA hW3 (by omega) with
        ⟨heap1, hinit, hA1, hW31, hlayout1⟩
      simp only [hinit]
      have hQ1 : heap1.ValidU64Slice Q (lenA - d) :=
        (hlayout1 Q (lenA - d)).mp (by simpa using hQ)
      have hB1 : heap1.ValidU64Slice B (d + 1) :=
        (hlayout1 B (d + 1)).mp hB
      have hR1 : heap1.ValidU64Slice R d := (hlayout1 R d).mp hRfull
      rcases quotientLoop_ok this Q B W3 (lenA - d) d lenA invLc heap1
        (lenA - d) hQ1 hB1 hW31 (Nat.le_refl _) (by omega) with
        ⟨heap2, hquot, hQ2, hB2, hW32, hlayout2⟩
      dsimp [invLc] at hquot ⊢
      simp only [hquot]
      have hR2 : heap2.ValidU64Slice R d := (hlayout2 R d).mp hR1
      rcases remainderLoop_ok this R W3 d lenA 0 heap2 hR2 hW32
        (by omega) (by omega) with
        ⟨heap3, hrem, hR3, hW33, hlayout3⟩
      simp only [hrem]
      have hQ3 : heap3.ValidU64Slice Q (lenA - d) :=
        (hlayout3 Q (lenA - d)).mp ((hlayout2 Q (lenA - d)).mp hQ1)
      rcases normaliseU64_ok heap3 Q (lenA - d) hQ3 with
        ⟨lenQ, hnormQ, hlenQ⟩
      simp only [hnormQ]
      rcases normaliseU64_ok heap3 R d hR3 with ⟨lenR, hnormR, hlenR⟩
      simp only [hnormR]
      refine ⟨heap3, lenQ, lenR, rfl, ?_, ?_⟩
      · simpa using hlenQ
      · simpa [Nat.min_eq_right hdleA] using hlenR

end CLPoly.Impl.StrictDivremRefinement
