-- CLPoly: Lean 4 formal verification of polynomial factorization
-- See proof/docs/ for blueprint and research report

-- Phase 1: 接口规约 + 顶层骨架
import CLPoly.Spec
import CLPoly.Pipeline.FactorZp
import CLPoly.Pipeline.FactorZZ

-- Phase 2: L3 数学基础
import CLPoly.Math.FiniteFieldFact
import CLPoly.Math.MvBasics

-- Phase 3: L2 算法模型
import CLPoly.Algorithm.DDF
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.EDF
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

-- 端到端实例化
import CLPoly.Pipeline.FactorZpInstantiate
import CLPoly.Pipeline.FactorZZInstantiate

-- L1 实例化（翻译代码管线）
import CLPoly.Pipeline.L1

-- Phase 0: 实验（保留参考，E4 被 Algorithm.DDF 取代）
import CLPoly.Experiment.E1_ZpPolyAPI
import CLPoly.Experiment.E2_TheoremBridge
import CLPoly.Experiment.E3_ZModPkDiv
-- import CLPoly.Experiment.E4_Termination  -- 与 Algorithm.DDF 冲突（ddfLoop 重定义）

-- Phase 9: L1 → L2 精化定理
import CLPoly.Refinement.Basic
import CLPoly.Refinement.DDF
import CLPoly.Refinement.SquarefreeZp
import CLPoly.Refinement.Generated.SquarefreeZp
import CLPoly.Refinement.ZZArith
import CLPoly.Refinement.ZpArith
