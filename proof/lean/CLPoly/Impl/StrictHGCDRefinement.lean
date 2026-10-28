import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRefinement

/-- Two polynomial pairs have exactly the same common divisors. -/
def SameCommonDivisors {R : Type*} [CommRing R]
    (a b c d : R) : Prop :=
  ∀ x : R, (x ∣ a ∧ x ∣ b) ↔ (x ∣ c ∧ x ∣ d)

/-- Algebraic orientation of the matrix produced by C++ HGCD: the original
pair is the matrix applied to the current Euclidean pair. -/
def HgcdTransform {R : Type*} [CommRing R]
    (left right currentA currentB m00 m01 m10 m11 : R) : Prop :=
  left = m00 * currentA + m01 * currentB ∧
    right = m10 * currentA + m11 * currentB

/-- One real Euclidean division followed by the two C++ row updates preserves
the HGCD transform orientation.  Each source row `[u,v]` becomes
`[v + quotient*u, u]`. -/
theorem hgcdTransform_euclid_step {R : Type*} [CommRing R]
    (left right currentA currentB remainder quotient
      m00 m01 m10 m11 : R)
    (htransform : HgcdTransform left right currentA currentB
      m00 m01 m10 m11)
    (hdivision : currentA = quotient * currentB + remainder) :
    HgcdTransform left right currentB remainder
      (m01 + quotient * m00) m00
      (m11 + quotient * m10) m10 := by
  rcases htransform with ⟨hleft, hright⟩
  constructor
  · rw [hleft, hdivision]
    ring
  · rw [hright, hdivision]
    ring

/-- The exact two C++ row updates negate the tracked 2×2 determinant. -/
theorem hgcdRowUpdate_determinant {R : Type*} [CommRing R]
    (quotient m00 m01 m10 m11 : R) :
    (m01 + quotient * m00) * m10 -
        m00 * (m11 + quotient * m10) =
      -(m00 * m11 - m01 * m10) := by
  ring

/-- Meaning of the integer sign returned by C++ HGCD: it records whether the
tracked polynomial matrix has determinant `1` or `-1`. -/
def HgcdSignedDet {R : Type*} [CommRing R]
    (sgn : Int) (m00 m01 m10 m11 : R) : Prop :=
  (sgn = 1 ∧ m00 * m11 - m01 * m10 = 1) ∨
    (sgn = -1 ∧ m00 * m11 - m01 * m10 = -1)

/-- Mathematical four-entry matrix state after the two C++ row updates in
one HGCD Euclidean iteration. -/
def hgcdStepEntries {R : Type*} [CommRing R]
    (quotient : R) (entries : Fin 4 → R) : Fin 4 → R
  | ⟨0, _⟩ => entries ⟨1, by omega⟩ + quotient * entries ⟨0, by omega⟩
  | ⟨1, _⟩ => entries ⟨0, by omega⟩
  | ⟨2, _⟩ => entries ⟨3, by omega⟩ + quotient * entries ⟨2, by omega⟩
  | ⟨3, _⟩ => entries ⟨2, by omega⟩

/-- Entry-array form of `hgcdTransform_euclid_step`, matching the exact
flat matrix layout used by C++. -/
theorem hgcdStepEntries_preserves_transform {R : Type*} [CommRing R]
    (left right currentA currentB remainder quotient : R)
    (entries : Fin 4 → R)
    (htransform : HgcdTransform left right currentA currentB
      (entries ⟨0, by omega⟩) (entries ⟨1, by omega⟩)
      (entries ⟨2, by omega⟩) (entries ⟨3, by omega⟩))
    (hdivision : currentA = quotient * currentB + remainder) :
    HgcdTransform left right currentB remainder
      (hgcdStepEntries quotient entries ⟨0, by omega⟩)
      (hgcdStepEntries quotient entries ⟨1, by omega⟩)
      (hgcdStepEntries quotient entries ⟨2, by omega⟩)
      (hgcdStepEntries quotient entries ⟨3, by omega⟩) := by
  simpa [hgcdStepEntries] using hgcdTransform_euclid_step left right
    currentA currentB remainder quotient (entries ⟨0, by omega⟩)
    (entries ⟨1, by omega⟩) (entries ⟨2, by omega⟩)
    (entries ⟨3, by omega⟩) htransform hdivision

/-- The two generated row updates and the source assignment `sgn = -sgn`
preserve the signed-determinant invariant. -/
theorem hgcdSignedDet_euclid_step {R : Type*} [CommRing R]
    (sgn : Int) (quotient m00 m01 m10 m11 : R)
    (hdet : HgcdSignedDet sgn m00 m01 m10 m11) :
    HgcdSignedDet (-sgn)
      (m01 + quotient * m00) m00
      (m11 + quotient * m10) m10 := by
  rcases hdet with ⟨hsgn, hdet⟩ | ⟨hsgn, hdet⟩
  · right
    constructor
    · simp [hsgn]
    · rw [hgcdRowUpdate_determinant, hdet]
  · left
    constructor
    · simp [hsgn]
    · rw [hgcdRowUpdate_determinant, hdet]
      simp

/-- Entry-array form of the synchronized matrix/sign update. -/
theorem hgcdStepEntries_preserves_signedDet {R : Type*} [CommRing R]
    (sgn : Int) (quotient : R) (entries : Fin 4 → R)
    (hdet : HgcdSignedDet sgn
      (entries ⟨0, by omega⟩) (entries ⟨1, by omega⟩)
      (entries ⟨2, by omega⟩) (entries ⟨3, by omega⟩)) :
    HgcdSignedDet (-sgn)
      (hgcdStepEntries quotient entries ⟨0, by omega⟩)
      (hgcdStepEntries quotient entries ⟨1, by omega⟩)
      (hgcdStepEntries quotient entries ⟨2, by omega⟩)
      (hgcdStepEntries quotient entries ⟨3, by omega⟩) := by
  simpa [hgcdStepEntries] using hgcdSignedDet_euclid_step sgn quotient
    (entries ⟨0, by omega⟩) (entries ⟨1, by omega⟩)
    (entries ⟨2, by omega⟩) (entries ⟨3, by omega⟩) hdet

/-- Equality of normalized Euclidean gcds follows solely from equality of
common divisors.  This does not execute either gcd algorithm. -/
theorem normalize_gcd_eq_of_sameCommonDivisors
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (a b c d : R) (hcommon : SameCommonDivisors a b c d) :
    normalize (EuclideanDomain.gcd a b) =
      normalize (EuclideanDomain.gcd c d) := by
  rw [normalize_eq_normalize_iff_associated]
  apply associated_of_dvd_dvd
  · have h := (hcommon (EuclideanDomain.gcd a b)).mp
      ⟨EuclideanDomain.gcd_dvd_left a b,
        EuclideanDomain.gcd_dvd_right a b⟩
    exact EuclideanDomain.dvd_gcd h.1 h.2
  · have h := (hcommon (EuclideanDomain.gcd c d)).mpr
      ⟨EuclideanDomain.gcd_dvd_left c d,
        EuclideanDomain.gcd_dvd_right c d⟩
    exact EuclideanDomain.dvd_gcd h.1 h.2

/-- A pair of forward and inverse 2×2 linear transformations preserves all
common divisors.  These are precisely the algebraic obligations supplied by
the generated HGCD matrix updates.  No determinant shortcut or external
specification is assumed. -/
theorem sameCommonDivisors_of_bidirectional_linear_combinations
    {R : Type*} [CommRing R]
    (a b c d m00 m01 m10 m11 n00 n01 n10 n11 : R)
    (hc : c = m00 * a + m01 * b)
    (hd : d = m10 * a + m11 * b)
    (ha : a = n00 * c + n01 * d)
    (hb : b = n10 * c + n11 * d) :
    SameCommonDivisors a b c d := by
  intro x
  constructor
  · rintro ⟨hxa, hxb⟩
    constructor
    · rw [hc]
      exact dvd_add (dvd_mul_of_dvd_right hxa m00)
        (dvd_mul_of_dvd_right hxb m01)
    · rw [hd]
      exact dvd_add (dvd_mul_of_dvd_right hxa m10)
        (dvd_mul_of_dvd_right hxb m11)
  · rintro ⟨hxc, hxd⟩
    constructor
    · rw [ha]
      exact dvd_add (dvd_mul_of_dvd_right hxc n00)
        (dvd_mul_of_dvd_right hxd n01)
    · rw [hb]
      exact dvd_add (dvd_mul_of_dvd_right hxc n10)
        (dvd_mul_of_dvd_right hxd n11)

/-- Core HGCD semantic invariant: an actual matrix transformation together
with its actual inverse preserves the normalized L2 gcd. -/
theorem normalize_gcd_eq_of_hgcd_transform
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (a b c d m00 m01 m10 m11 n00 n01 n10 n11 : R)
    (hc : c = m00 * a + m01 * b)
    (hd : d = m10 * a + m11 * b)
    (ha : a = n00 * c + n01 * d)
    (hb : b = n10 * c + n11 * d) :
    normalize (EuclideanDomain.gcd a b) =
      normalize (EuclideanDomain.gcd c d) := by
  exact normalize_gcd_eq_of_sameCommonDivisors a b c d
    (sameCommonDivisors_of_bidirectional_linear_combinations
      a b c d m00 m01 m10 m11 n00 n01 n10 n11 hc hd ha hb)

/-- Determinant-one form used by one parity of the generated HGCD matrix. -/
theorem normalize_gcd_eq_of_det_one_transform
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (a b c d m00 m01 m10 m11 : R)
    (hc : c = m00 * a + m01 * b)
    (hd : d = m10 * a + m11 * b)
    (hdet : m00 * m11 - m01 * m10 = 1) :
    normalize (EuclideanDomain.gcd a b) =
      normalize (EuclideanDomain.gcd c d) := by
  apply normalize_gcd_eq_of_hgcd_transform a b c d m00 m01 m10 m11
      m11 (-m01) (-m10) m00 hc hd
  · rw [hc, hd]
    symm
    calc
      m11 * (m00 * a + m01 * b) + -m01 * (m10 * a + m11 * b) =
          (m00 * m11 - m01 * m10) * a := by ring
      _ = a := by rw [hdet, one_mul]
  · rw [hc, hd]
    symm
    calc
      -m10 * (m00 * a + m01 * b) + m00 * (m10 * a + m11 * b) =
          (m00 * m11 - m01 * m10) * b := by ring
      _ = b := by rw [hdet, one_mul]

/-- Determinant-minus-one form used by the other parity of the generated
HGCD matrix. -/
theorem normalize_gcd_eq_of_det_neg_one_transform
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (a b c d m00 m01 m10 m11 : R)
    (hc : c = m00 * a + m01 * b)
    (hd : d = m10 * a + m11 * b)
    (hdet : m00 * m11 - m01 * m10 = -1) :
    normalize (EuclideanDomain.gcd a b) =
      normalize (EuclideanDomain.gcd c d) := by
  apply normalize_gcd_eq_of_hgcd_transform a b c d m00 m01 m10 m11
      (-m11) m01 m10 (-m00) hc hd
  · rw [hc, hd]
    symm
    calc
      -m11 * (m00 * a + m01 * b) + m01 * (m10 * a + m11 * b) =
          -(m00 * m11 - m01 * m10) * a := by ring
      _ = a := by rw [hdet]; ring
  · rw [hc, hd]
    symm
    calc
      m10 * (m00 * a + m01 * b) + -m00 * (m10 * a + m11 * b) =
          -(m00 * m11 - m01 * m10) * b := by ring
      _ = b := by rw [hdet]; ring

/-- A represented HGCD transform with its actual signed determinant preserves
the normalized gcd of the original and current polynomial pairs. -/
theorem normalize_gcd_eq_of_hgcd_signed_transform
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (sgn : Int) (left right currentA currentB m00 m01 m10 m11 : R)
    (htransform : HgcdTransform left right currentA currentB
      m00 m01 m10 m11)
    (hdet : HgcdSignedDet sgn m00 m01 m10 m11) :
    normalize (EuclideanDomain.gcd left right) =
      normalize (EuclideanDomain.gcd currentA currentB) := by
  rcases htransform with ⟨hleft, hright⟩
  rcases hdet with ⟨_, hdet⟩ | ⟨_, hdet⟩
  · exact (normalize_gcd_eq_of_det_one_transform currentA currentB left
      right m00 m01 m10 m11 hleft hright hdet).symm
  · exact (normalize_gcd_eq_of_det_neg_one_transform currentA currentB left
      right m00 m01 m10 m11 hleft hright hdet).symm

end CLPoly.Impl.StrictHGCDRefinement
