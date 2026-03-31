-- 测试桩：替代 opaque 声明，使翻译函数可 #eval
-- 仅用于背靠背测试，不用于证明

-- ============================================================
-- 多项式求导（朴素实现）
-- ============================================================

def derivative_stub (f : SparsePolyZp) : SparsePolyZp :=
  f.filterMap (fun (m, c) =>
    if m.deg == 0 then none
    else some (⟨m.deg - 1⟩, ⟨c.val * m.deg.toNat.toUInt64 % c.prime, c.prime⟩))

-- ============================================================
-- 多项式 GCD（欧几里得算法，朴素版）
-- ============================================================

-- 简化：返回第一个非零参数（不是真 GCD，仅测试用）
partial def polynomial_GCD_stub (f g : SparsePolyZp) : SparsePolyZp :=
  if g.isEmpty then f else g

-- ============================================================
-- 多项式除法（不实现，返回空）
-- ============================================================

def pair_vec_div_stub (_a _b _c _d : SparsePolyZp) : Unit := ()

-- ============================================================
-- 其他桩
-- ============================================================

def get_deg_stub (f : SparsePolyZp) : UInt64 :=
  if f.isEmpty then 0 else (f[0]!).fst.deg

def normalize_stub (f : SparsePolyZp) : SparsePolyZp := f

def comp_stub (_f : SparsePolyZp) : UInt64 := 0

def inv_stub (a : Zp) : Zp := a

def number_stub (a : UInt64) (p : UInt64) : Zp := ⟨a % p, p⟩
