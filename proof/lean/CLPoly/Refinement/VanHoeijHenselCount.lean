import CLPoly.Refinement.Hensel

set_option autoImplicit false

namespace Refinement.StrictHensel
open _root_.Generated.StrictHensel

theorem henselTreeInternalNodeCount_add_one
    (start stop : Nat) (hpositive : start < stop) :
    henselTreeInternalNodeCount start stop + 1 = stop - start := by
  rw [henselTreeInternalNodeCount]
  split
  next hlength =>
    have hmidLt := henselTreeMidpoint_lt_stop start stop hlength
    have hstartLt := henselTreeMidpoint_gt_start start stop hlength
    have hleft := henselTreeInternalNodeCount_add_one
      start ((start + stop) / 2) hstartLt
    have hright := henselTreeInternalNodeCount_add_one
      ((start + stop) / 2) stop hmidLt
    dsimp only
    omega
  next hshort => omega
termination_by stop - start

theorem extract_factors_success_size (tree : HenselLiftTree)
    (nodes : Array HenselNode) (input output : Array SparsePolyZZ)
    (hrun : __hensel_extract_factors_raw_ir tree nodes input = .ok output) :
    output.size = input.size + tree.nodeCount + 1 := by
  induction hmeasure : tree.nodeCount using Nat.strong_induction_on
      generalizing tree input output with
  | h measure ih =>
    cases tree with
    | node index left right =>
      cases hn : nodes[index]? with
      | none => simp [__hensel_extract_factors_raw_ir, hn] at hrun
      | some node =>
        cases left with
        | none =>
          cases right with
          | none =>
            simp only [__hensel_extract_factors_raw_ir, hn, bind, Except.bind] at hrun
            have he := Except.ok.inj hrun
            rw [← he]
            change 1 = measure at hmeasure
            simp [← hmeasure]
          | some child =>
            have hc : __hensel_extract_factors_raw_ir child nodes (input.push node.g) = .ok output := by
              simpa only [__hensel_extract_factors_raw_ir, hn, bind, Except.bind] using hrun
            have hlt := HenselLiftTree.right_nodeCount_lt index none child
            rw [hmeasure] at hlt
            have hs := ih child.nodeCount hlt child (input.push node.g) output hc rfl
            change 1 + 0 + child.nodeCount = measure at hmeasure
            simp only [Array.size_push] at hs
            omega
        | some child =>
          cases hl : __hensel_extract_factors_raw_ir child nodes input with
          | error fault => simp [__hensel_extract_factors_raw_ir, hn, hl, bind, Except.bind] at hrun
          | ok middle =>
            have hlt := HenselLiftTree.left_nodeCount_lt index child right
            rw [hmeasure] at hlt
            have hs := ih child.nodeCount hlt child input middle hl rfl
            cases right with
            | none =>
              simp only [__hensel_extract_factors_raw_ir, hn, hl, bind, Except.bind] at hrun
              have he := Except.ok.inj hrun
              rw [← he]
              change 1 + child.nodeCount + 0 = measure at hmeasure
              simp only [Array.size_push]
              omega
            | some rightChild =>
              have hc : __hensel_extract_factors_raw_ir rightChild nodes middle = .ok output := by
                simpa only [__hensel_extract_factors_raw_ir, hn, hl, bind, Except.bind] using hrun
              have hrt := HenselLiftTree.right_nodeCount_lt index (some child) rightChild
              rw [hmeasure] at hrt
              have hr := ih rightChild.nodeCount hrt rightChild middle output hc rfl
              change 1 + child.nodeCount + rightChild.nodeCount = measure at hmeasure
              omega

theorem normalize_result_success_size (input output : Array SparsePolyZZ) (modulus : ZZ)
    (hrun : __hensel_normalize_result_raw_ir input modulus = .ok output) :
    output.size = input.size := by
  unfold __hensel_normalize_result_raw_ir at hrun
  split at hrun
  next hn => cases Except.ok.inj hrun; rfl
  next first hf =>
    split at hrun
    next hn => contradiction
    next leading hl =>
      split at hrun
      next hnotOne =>
        dsimp only at hrun
        split at hrun
        next hinvert =>
          cases hm : __upoly_mod_coeff_raw_ir (scaleCoeffs first (ZZ.invert 0 leading.2 modulus).2) modulus with
          | error fault => simp only [hm, bind, Except.bind] at hrun; contradiction
          | ok normalized =>
            simp only [hm, bind, Except.bind] at hrun
            have he := Except.ok.inj hrun
            rw [← he]
            simp
        next hinvert => contradiction
      next hone => cases Except.ok.inj hrun; rfl

theorem hensel_entry_success_count (step : HenselStepRawOps) (treeOps : HenselTreeBuildRawOps)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64) (target : Int32)
    (hp : 2 ≤ p.toNat) (output : Array SparsePolyZZ × ZZ)
    (hrun : __hensel_lift_upoly_raw_ir step treeOps f factors p target hp = .ok output) :
    output.1.size = factors.size ∧ output.1.size < 2 ^ 31 := by
  have hcount := __hensel_lift_upoly_raw_ir_factorCount_of_success step treeOps f factors p target hp output hrun
  have hfits := __hensel_lift_upoly_raw_ir_factorCountFits_of_success step treeOps f factors p target hp output hrun
  have hguard : __hensel_factor_count_fits_ir factors.size = true := by
    simpa only [__hensel_factor_count_fits_ir, decide_eq_true_eq] using hfits
  have hnodes : (henselTreeBuildTopologyRawIR 0 factors.size 0).nodeCount + 1 = factors.size := by
    rw [henselTreeBuildTopologyRawIR_eq,
      henselTreeBuildTopology_nodeCount 0 factors.size 0 (by omega)]
    simpa only [Nat.sub_zero] using henselTreeInternalNodeCount_add_one 0 factors.size (by omega)
  unfold __hensel_lift_upoly_raw_ir at hrun
  simp only [hcount, hguard, ↓reduceIte, pure, Except.pure, bind, Except.bind] at hrun
  split at hrun
  next fault he => contradiction
  next targetZZ ht =>
    split at hrun
    next htarget =>
      split at hrun
      next fault he => contradiction
      next adjusted ha =>
        split at hrun
        next fault he => contradiction
        next nodes hn =>
          split at hrun
          next fault he => contradiction
          next lifted hl =>
            split at hrun
            next fault he => contradiction
            next extracted he =>
              split at hrun
              next fault he => contradiction
              next normalized hnorm =>
                have heq := Except.ok.inj hrun
                rw [← heq]
                have hex := extract_factors_success_size
                  (henselTreeBuildTopologyRawIR 0 factors.size 0) lifted.1 #[] extracted he
                have hns := normalize_result_success_size extracted normalized lifted.2 hnorm
                simp only [Array.size_empty] at hex
                dsimp only
                constructor <;> omega
    next htarget => contradiction

end Refinement.StrictHensel
