-- CLPoly's public Lean verification surface.

-- Specifications and mathematical factorization pipelines.
import CLPoly.Spec
import CLPoly.Pipeline.FactorZp
import CLPoly.Pipeline.FactorZZ

-- Mathematical foundations and executable L2 algorithms.
import CLPoly.Math.FiniteFieldFact
import CLPoly.Math.MvBasics
import CLPoly.Algorithm.DDF
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.EDF
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

-- L2 existence/correctness instances retained for downstream algorithms.
import CLPoly.Pipeline.FactorZpInstantiate
import CLPoly.Pipeline.FactorZZInstantiate

-- Strict generated C++ L1 → Lean L2 refinement and public final contracts.
import CLPoly.Refinement.Basic
import CLPoly.Refinement.DDF
import CLPoly.Refinement.SquarefreeZp
import CLPoly.Refinement.Generated
import CLPoly.Refinement.ZZArith
import CLPoly.Refinement.ZpArith
import CLPoly.Pipeline.L1
