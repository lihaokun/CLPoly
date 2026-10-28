/-
  E4: 递归函数终止性与函数归纳
  Phase 0 实验 — 确认 DDF 循环的 termination_by / decreasing_by 方案

  DDF 循环特征：
  - n = deg(f_star) 单调不增
  - d 单调递增
  - 终止条件：n < 2*d
  - 递减度量：n + 1 - 2*d
-/
open Nat List

-- ============================================================
-- 1. Toy DDF：顶层递归版本
-- ============================================================

/-- 模拟 DDF 循环（顶层递归函数）：
    - n 模拟 deg(f_star)
    - d 从 1 开始递增
    - 每步：若 n % 3 = 0 则模拟提取度为 d 的因子（n 减小）
    - 当 n < 2*d 时终止 -/
def ddfLoop (n d : Nat) (acc : List (Nat × Nat)) : List (Nat × Nat) :=
  if h : n < 2 * d then
    if 0 < n then acc ++ [(n, n)] else acc
  else
    if n % 3 = 0 then
      ddfLoop (n - d) (d + 1) (acc ++ [(d, d)])
    else
      ddfLoop n (d + 1) acc
termination_by n + 1 - 2 * d

def toyDDF (n : Nat) : List (Nat × Nat) := ddfLoop n 1 []

-- 测试
#eval toyDDF 10
#eval toyDDF 6
#eval toyDDF 0

-- ============================================================
-- 2. 函数归纳原理
-- ============================================================

-- Lean 4 是否自动生成 .induct？
#check @ddfLoop.induct

-- ============================================================
-- 3. 用归纳原理证明简单性质
-- ============================================================

-- 先检查 .induct 签名
-- 用 True 占位，主要目的是验证 .induct 可用
theorem ddfLoop_true (n d : Nat) (acc : List (Nat × Nat)) :
    True := by
  trivial

-- ============================================================
-- 4. 手动 decreasing_by 版本
-- ============================================================

def ddfLoop2 (n d : Nat) (acc : List (Nat × Nat)) : List (Nat × Nat) :=
  if h : n < 2 * d then
    acc
  else
    if hmod : n % 3 = 0 then
      ddfLoop2 (n - d) (d + 1) (acc ++ [(d, d)])
    else
      ddfLoop2 n (d + 1) acc
termination_by n + 1 - 2 * d
decreasing_by all_goals omega

-- ============================================================
-- 5. partial def 备选方案（EDF 概率循环）
-- ============================================================

/-- EDF 的 Cantor-Zassenhaus：概率性，不保证终止 -/
partial def toyEDF (n : Nat) : Nat :=
  if n ≤ 1 then n
  else toyEDF (n / 2)

#eval toyEDF 100

-- ============================================================
-- 6. fuel 参数备选方案
-- ============================================================

/-- fuel 参数保证终止 -/
def toyEDF_fuel (fuel : Nat) (n : Nat) : Option Nat :=
  match fuel with
  | 0 => none
  | fuel' + 1 =>
    if n ≤ 1 then some n
    else toyEDF_fuel fuel' (n / 2)

#eval toyEDF_fuel 100 1000
