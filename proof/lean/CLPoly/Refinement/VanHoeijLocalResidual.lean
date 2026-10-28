import CLPoly.Refinement.VanHoeijResidualPolynomial
import CLPoly.Algorithm.VanHoeijLocalCld

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The actual modular CLD certificates identify the mapped physical row
combination with its weighted logarithmic derivatives. -/
theorem rowCldCombination_modular_cld
    (f : SparsePolyZZ) (active original : Array SparsePolyZZ) (size modulus : Nat)
    (coefficients : Nat → Int)
    (hcld : ∀ i, i < original.size → ModularCldValue f active[i]! original[i]! modulus) :
    Polynomial.map (Int.castRingHom (ZMod modulus)) (rowCldCombination original size coefficients) =
      ∑ i : Fin size, if i.val < original.size then
        Polynomial.C ((coefficients i.val : Int) : ZMod modulus) *
          ((StrictHensel.toPolyMod modulus f /ₘ StrictHensel.toPolyMod modulus active[i.val]!) *
            Polynomial.derivative (StrictHensel.toPolyMod modulus active[i.val]!)) else 0 := by
  rw [rowCldCombination, Polynomial.map_sum]
  apply Finset.sum_congr rfl
  intro i _
  by_cases hi : i.val < original.size
  · simp only [hi, ↓reduceIte, Polynomial.map_mul, Polynomial.map_C]
    change Polynomial.C ((coefficients i.val : Int) : ZMod modulus) *
      StrictHensel.toPolyMod modulus original[i.val]! = _
    rw [(hcld i.val hi).semantic]
  · simp [hi]

/-- The constructed residual has the needed local divisibility relation.
The monic/divisor/coprimality premises concern the lifted input factors;
they must still be discharged from the actual Hensel/good-prime execution. -/
theorem rowResidualPolynomial_local_dvd
    (f : SparsePolyZZ) (active original : Array SparsePolyZZ) (size width modulus : Nat)
    (coefficients : Nat → Int)
    (hcld : ∀ i, i < original.size → ModularCldValue f active[i]! original[i]! modulus)
    (selected : Fin size) (hselected : selected.val < original.size)
    (hmonic : ∀ i, i < original.size → (StrictHensel.toPolyMod modulus active[i]!).Monic)
    (hdiv : ∀ i, i < original.size →
      StrictHensel.toPolyMod modulus active[i]! ∣ StrictHensel.toPolyMod modulus f)
    (hcop : ∀ i, i < original.size → i ≠ selected.val →
      IsCoprime (StrictHensel.toPolyMod modulus active[selected.val]!)
        (StrictHensel.toPolyMod modulus active[i]!)) :
    StrictHensel.toPolyMod modulus active[selected.val]! ∣
      Polynomial.map (Int.castRingHom (ZMod modulus))
        (rowResidualPolynomial original size width coefficients (modulus : Int)) -
        Polynomial.C ((coefficients selected.val : Int) : ZMod modulus) *
          Polynomial.derivative (StrictHensel.toPolyMod modulus f) := by
  rw [rowResidualPolynomial_modular,
    rowCldCombination_modular_cld f active original size modulus coefficients hcld]
  let support := Finset.univ.filter (fun i : Fin size => i.val < original.size)
  have hs : selected ∈ support := by simp [support, hselected]
  have hm : ∀ i ∈ support, (StrictHensel.toPolyMod modulus active[i.val]!).Monic := by
    intro i hi
    exact hmonic i.val (Finset.mem_filter.mp hi).2
  have hd : ∀ i ∈ support,
      StrictHensel.toPolyMod modulus active[i.val]! ∣ StrictHensel.toPolyMod modulus f := by
    intro i hi
    exact hdiv i.val (Finset.mem_filter.mp hi).2
  have hc : ∀ i ∈ support, i ≠ selected →
      IsCoprime (StrictHensel.toPolyMod modulus active[selected.val]!)
        (StrictHensel.toPolyMod modulus active[i.val]!) := by
    intro i hi hne
    exact hcop i.val (Finset.mem_filter.mp hi).2 (fun he => hne (Fin.ext he))
  have hl := CLPoly.VanHoeij.local_factor_dvd_weighted_cld_residual support
    (StrictHensel.toPolyMod modulus f) (fun i : Fin size => StrictHensel.toPolyMod modulus active[i.val]!)
    (fun i => ((coefficients i.val : Int) : ZMod modulus)) selected hs hm hd hc
  simpa only [support, Finset.sum_filter] using hl

end Refinement.StrictRecombine
