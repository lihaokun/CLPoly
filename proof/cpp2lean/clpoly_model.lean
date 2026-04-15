/-
  CLPoly 可信基类 Lean 模型

  定义因式分解模块依赖的所有类型和操作。
  翻译器生成的代码 import 此文件。
  背靠背测试中此文件提供可执行的原语实现。
-/

set_option autoImplicit false

-- ============================================================
-- §1. Zp：Z/pZ 系数
-- ============================================================

structure Zp where
  val : UInt64
  prime : UInt64
deriving Repr, Inhabited, BEq

namespace Zp

def ofInt (v : Int) (p : UInt64) : Zp :=
  let pn : Int := p.toNat
  let r := v % pn
  let r := if r < 0 then r + pn else r
  ⟨r.toNat.toUInt64, p⟩

def ofUInt64 (v p : UInt64) : Zp := ⟨v % p, p⟩

instance : Add Zp where add a b := ⟨(a.val + b.val) % a.prime, a.prime⟩
instance : Sub Zp where sub a b := ⟨(a.val + a.prime - b.val) % a.prime, a.prime⟩
instance : Mul Zp where mul a b := ⟨(a.val * b.val) % a.prime, a.prime⟩
instance : Neg Zp where neg a := ⟨(a.prime - a.val) % a.prime, a.prime⟩

-- 扩展欧几里得：gcd(a, b) = a*x + b*y，返回 (gcd, x)
partial def extGcdAux (old_r r old_s s : Int) : Int × Int :=
  if r == 0 then (old_r, old_s)
  else
    let q := old_r / r
    extGcdAux r (old_r - q * r) s (old_s - q * s)

def modInv (a p : UInt64) : UInt64 :=
  if a == 0 then 0
  else
    let (_, s) := extGcdAux (p.toNat : Int) (a.toNat : Int) 0 1
    let r := s % (p.toNat : Int)
    let r := if r < 0 then r + p.toNat else r
    r.toNat.toUInt64

def inv (a : Zp) : Zp := ⟨modInv a.val a.prime, a.prime⟩
def div (a b : Zp) : Zp := a * b.inv

end Zp

-- Zp 隐式转换（对应 C++ 的 implicit conversion operators）
instance : Coe Zp UInt64 where coe z := z.val
instance : Coe Zp Int where coe z := z.val.toNat

-- ============================================================
-- §2. ZZ：大整数
-- ============================================================

abbrev ZZ := Int

-- ============================================================
-- §3. UMonomial：单变量单项式
-- ============================================================

structure UMonomial where
  deg : UInt64
deriving Repr, Inhabited, BEq

-- ============================================================
-- §4. SparsePolyZp：Z/pZ 上稀疏多项式
-- ============================================================

abbrev SparsePolyZp := Array (UMonomial × Zp)

namespace SparsePolyZp

def empty : SparsePolyZp := #[]

def front! (f : SparsePolyZp) : UMonomial × Zp := f[0]!
def back! (f : SparsePolyZp) : UMonomial × Zp := f[f.size - 1]!
def getDeg (f : SparsePolyZp) : UInt64 := if f.isEmpty then 0 else f[0]!.fst.deg
def size_u64 (f : SparsePolyZp) : UInt64 := f.size.toUInt64

def normalization (f : SparsePolyZp) : SparsePolyZp :=
  f.filter (fun (_, c) => c.val != 0)

def reserve (f : SparsePolyZp) (_n : UInt64) : SparsePolyZp := f
def data (f : SparsePolyZp) : SparsePolyZp := f

-- 求导：d/dx (Σ c_i x^{d_i}) = Σ d_i * c_i * x^{d_i - 1}
def derivative (f : SparsePolyZp) : SparsePolyZp :=
  if f.isEmpty then #[]
  else
    let p := f[0]!.snd.prime
    f.filterMap (fun (m, c) =>
      if m.deg == 0 then none
      else some (⟨m.deg - 1⟩, ⟨c.val * m.deg % p, p⟩))

-- GCD：欧几里得算法（需要多项式除法，暂用简化版）
-- 完整实现需要 divmod，此处仅供背靠背测试框架编译
partial def gcd (f g : SparsePolyZp) : SparsePolyZp :=
  if g.isEmpty then f else gcd g (f)  -- TODO: 需要 poly mod

-- 多项式除法：f = q * g + r
-- 完整实现需要逐项消去，暂返回 (f, #[])
def divmod (f _g : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  (f, #[])  -- TODO: 实现多项式长除法

-- 比较器：单变量多项式的单项式序（降幂排列）
-- 在 Lean 模型中不需要显式比较器（约定高次在前），返回 0 占位
def comp (_f : SparsePolyZp) : UInt64 := 0

end SparsePolyZp

-- ============================================================
-- §5. Array 辅助
-- ============================================================

def Array.front! {α : Type} [Inhabited α] (a : Array α) : α := a[0]!
def Array.size_u64 {α : Type} (a : Array α) : UInt64 := a.size.toUInt64

-- ============================================================
-- §5a. Rng：伪随机数生成器模型
-- ============================================================

-- 对应 C++ 的 std::mt19937 + std::uniform_int_distribution。
-- 状态 = UInt64 种子。next 返回 [0, upper) 的伪随机值。
-- 简化实现（xorshift64）：足够背靠背测试的确定性，不求密码学安全。

namespace Rng

def next (seed upper : UInt64) : UInt64 :=
  if upper == 0 then 0
  else
    -- xorshift64 步进
    let s := seed ^^^ (seed <<< 13)
    let s := s ^^^ (s >>> 7)
    let s := s ^^^ (s <<< 17)
    s % upper

def step (seed : UInt64) : UInt64 :=
  let s := seed ^^^ (seed <<< 13)
  let s := s ^^^ (s >>> 7)
  s ^^^ (s <<< 17)

end Rng

-- ============================================================
-- §5a2. 迭代器压缩模式的函数式等价（设计 §6）
-- ============================================================

-- __upoly_mod_coeff(f, m) 的语义：对每个系数做 fdiv_r，过滤掉零
def SparsePolyZZ.modCoeff (f : SparsePolyZZ) (m : Int) : SparsePolyZZ :=
  f.filterMap (fun (mono, coeff) =>
    let c := coeff % m
    if c != 0 then some (mono, c) else none)

-- __hensel_step 中的零压缩：只过滤零系数项
def SparsePolyZZ.compactNonzero (f : SparsePolyZZ) : SparsePolyZZ :=
  f.filter (fun (_, coeff) => coeff != 0)

-- ============================================================
-- §5b. StdMap：有序映射模型
-- ============================================================

-- 对应 C++ 的 std::map<K, V>。
-- 实现为 List (K × V)（有序对列表）。
-- 语义正确（find/insert/erase），性能 O(n)（翻译目标是正确性不是性能）。

abbrev StdMap (K V : Type) := List (K × V)

namespace StdMap

def empty : StdMap K V := []

def find? [BEq K] (m : StdMap K V) (k : K) : Option V :=
  match m.find? (fun (k', _) => k == k') with
  | some (_, v) => some v
  | none => none

def find! [BEq K] [Inhabited V] (m : StdMap K V) (k : K) : V :=
  match find? m k with
  | some v => v
  | none => default

def insert [BEq K] (m : StdMap K V) (k : K) (v : V) : StdMap K V :=
  (k, v) :: m.filter (fun (k', _) => !(k == k'))

def erase [BEq K] (m : StdMap K V) (k : K) : StdMap K V :=
  m.filter (fun (k', _) => !(k == k'))

def size (m : StdMap K V) : Nat := m.length

def isEmpty (m : StdMap K V) : Bool := m.isEmpty

-- end/begin 占位（迭代器语义，翻译中用 List 遍历替代）
def end_ (m : StdMap K V) : Nat := m.length
def begin_ (m : StdMap K V) : Nat := 0

end StdMap

-- ============================================================
-- §5c. 其他辅助类型
-- ============================================================

-- 多变量多项式（占位）
abbrev MvPolyZZ := Array (Array (UInt64 × UInt64) × Int)
abbrev MvPolyZp := Array (Array (UInt64 × UInt64) × Zp)
abbrev MvMonomial := Array (UInt64 × UInt64)
abbrev Variable := Array (String × Int)
abbrev SparsePolyZZ := Array (UMonomial × Int)
abbrev LLLMatrix := Array (Array Int)
abbrev HenselNode := Array Int  -- 占位
abbrev Factorization := Array Int  -- 占位
abbrev PrimeSelectionResult := Array Int  -- 占位
abbrev WangLcResult := Array Int  -- 占位

-- ============================================================
-- §6. 验证测试
-- ============================================================

-- Zp.ofInt 对负数
#eval Zp.ofInt (-1) 5        -- 应为 ⟨4, 5⟩
#eval Zp.ofInt 0 3           -- 应为 ⟨0, 3⟩
#eval Zp.ofInt 7 13          -- 应为 ⟨7, 13⟩
#eval Zp.ofInt 100 17        -- 应为 ⟨15, 17⟩

-- Zp 算术
#eval (Zp.ofUInt64 3 7) + (Zp.ofUInt64 5 7)   -- 应为 ⟨1, 7⟩ (3+5=8≡1)
#eval (Zp.ofUInt64 3 7) * (Zp.ofUInt64 5 7)   -- 应为 ⟨1, 7⟩ (3*5=15≡1)

-- 模逆
#eval Zp.modInv 3 7          -- 应为 5 (3*5=15≡1 mod 7)

-- SparsePolyZp
#eval SparsePolyZp.getDeg #[(⟨3⟩, ⟨2, 5⟩), (⟨1⟩, ⟨1, 5⟩)]  -- 应为 3

-- derivative: d/dx(2x^3 + x) = 6x^2 + 1 over F_5 = x^2 + 1
#eval SparsePolyZp.derivative #[(⟨3⟩, ⟨2, 5⟩), (⟨1⟩, ⟨1, 5⟩)]
