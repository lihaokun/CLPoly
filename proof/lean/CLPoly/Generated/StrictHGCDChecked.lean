import CLPoly.Generated.StrictHGCD

set_option autoImplicit false

namespace Generated.StrictHGCD


/-- Well-founded raw execution of the generated `_hgcd_recursive` body.
The two checks are the executable form of the source recursion contract:
both operands must remain ordered and the first length must strictly decrease.
They produce a raw assertion fault on an invalid call and never substitute a
specification or an L2 implementation. -/
def hgcdRecursiveCallChecked (this : DenseUPolyZp) : HgcdRecursiveCall
  | M, hM, computeM, A, B, a, b, lenA, lenB, W, scratch, heap =>
      hgcdRecursiveBody this
        (fun childM childHM childCompute childA childB childa childb
            childLenA childLenB childW childScratch childHeap =>
          if horder : childLenB < childLenA then
            if hdecrease : childLenA < lenA then
              hgcdRecursiveCallChecked this childM childHM childCompute childA
                childB childa childb childLenA childLenB childW childScratch
                childHeap
            else
              .error .assertionFailure
          else
            .error .assertionFailure)
        M hM computeM A B a b lenA lenB W scratch heap
termination_by M hM computeM A B a b lenA lenB W scratch heap => lenA
decreasing_by exact hdecrease

end Generated.StrictHGCD
