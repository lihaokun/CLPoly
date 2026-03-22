-- CLPoly: Lean 4 formal verification of polynomial factorization
-- See proof/docs/ for blueprint and research report

-- Phase 1: 接口规约 + 顶层骨架
import CLPoly.Spec
import CLPoly.Pipeline.FactorZp
import CLPoly.Pipeline.FactorZZ

-- Phase 3: L2 算法模型
import CLPoly.Algorithm.DDF
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.EDF
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

-- 端到端实例化
import CLPoly.Pipeline.FactorZpInstantiate

-- Phase 0: 实验（保留参考，E4 被 Algorithm.DDF 取代）
import CLPoly.Experiment.E1_ZpPolyAPI
import CLPoly.Experiment.E2_TheoremBridge
import CLPoly.Experiment.E3_ZModPkDiv
-- import CLPoly.Experiment.E4_Termination  -- 与 Algorithm.DDF 冲突（ddfLoop 重定义）
