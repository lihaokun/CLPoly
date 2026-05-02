/-
  CLPoly 可信基类 Lean 模型

  定义因式分解模块依赖的所有类型和操作。
  翻译器生成的代码 import 此文件。
  背靠背测试中此文件提供可执行的原语实现。
-/

-- Model.lean 来自 v1 cpp2lean/clpoly_model.lean，借用其隐式变量风格。
-- proof/lean lakefile 全局禁用 autoImplicit，本文件局部启用。
set_option autoImplicit true

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

-- UMonomial.deg 用 Nat（C++ size_t；Lean 4 Nat 同时支持 .toInt32 等转换，
-- 解决 (term_1.fst.deg).toInt32 等场景的 UInt64.toInt32 invalid field 问题）
structure UMonomial where
  deg : Nat
deriving Repr, Inhabited, BEq

-- ============================================================
-- §4. SparsePolyZp：Z/pZ 上稀疏多项式
-- ============================================================

abbrev SparsePolyZp := Array (UMonomial × Zp)

namespace SparsePolyZp

def empty : SparsePolyZp := #[]

def front! (f : SparsePolyZp) : UMonomial × Zp := f[0]!
def back! (f : SparsePolyZp) : UMonomial × Zp := f[f.size - 1]!
def getDeg (f : SparsePolyZp) : UInt64 := if f.isEmpty then 0 else f[0]!.fst.deg.toUInt64
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
      else some (⟨m.deg - 1⟩, ⟨c.val * m.deg.toUInt64 % p, p⟩))

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
-- §5a2. 迭代器压缩模式的函数式等价（设计 §6）—— SparsePolyZZ 操作
-- 移到 §5c（abbrev SparsePolyZZ）之后；保留 forward declaration 以 import 兼容
-- ============================================================

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
  match List.find? (fun (k', _) => k == k') m with
  | some (_, v) => some v
  | none => none

def find! [BEq K] [Inhabited V] (m : StdMap K V) (k : K) : V :=
  match StdMap.find? m k with
  | some v => v
  | none => default

-- get! 别名（Pass 5 把 C++ map[k] 解析为 StdMap.get!，与 find! 等价语义）
def get! [BEq K] [Inhabited V] (m : StdMap K V) (k : K) : V := find! m k

def insert [BEq K] (m : StdMap K V) (k : K) (v : V) : StdMap K V :=
  (k, v) :: m.filter (fun (k', _) => !(k == k'))

def erase [BEq K] (m : StdMap K V) (k : K) : StdMap K V :=
  List.filter (fun (k', _) => !(k == k')) m

-- StdMap.filter / filterMap：predicate / map 风格的过滤（与 Array.filter 对齐）
def filter (m : StdMap K V) (p : K × V → Bool) : StdMap K V := List.filter p m
def filterMap (m : StdMap K V) (f : K × V → Option (K × V)) : StdMap K V :=
  List.filterMap f m

def size (m : StdMap K V) : Nat := m.length

def isEmpty (m : StdMap K V) : Bool := List.isEmpty m

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
-- Monomial 是 C++ Poly::monomial_type 的 typedef
-- 实际 corpus: 迭代 Monomial 得到 (Variable, Int64) pairs（指数表达）
-- 与 MvMonomial 不同（MvMonomial 是 (UInt64, UInt64) 即 (var_id, deg)）
abbrev Monomial := Array (Variable × Int64)
-- C++ side polynomial<Zp> 的 Lean 别名（语义相同 MvPolyZp）
abbrev PolyZp := MvPolyZp
abbrev PolyZZ := MvPolyZZ
abbrev PolyQQ := Array (Array (UInt64 × UInt64) × Rat)

-- MvPolyZp 操作（stub；实际语义留 Pass 上游 / B2B 测试细化）
def MvPolyZp.normalization (f : MvPolyZp) : MvPolyZp := f
def MvPolyZp.mk (f : MvPolyZp) : MvPolyZp := f
def MvPolyZZ.normalization (f : MvPolyZZ) : MvPolyZZ := f
def MvPolyZZ.mk (f : MvPolyZZ) : MvPolyZZ := f
-- 通用 stub（与 SparsePolyZZ 解耦，无前向引用）
def __write__ (_old : α) (new : α) : α := new
def polynomial_GCD [Inhabited α] (_a _b : α) : α := default
-- pair_vec_div: 4 参数版本（C++ side: pair_vec_div(f, g, q, comp) → 返回 quotient）
-- 占位实现，B2B 测试细化（comp 通常是比较器/函数对象，签名宽松）
def pair_vec_div [Inhabited α] (_f _g _q : α) (_comp : β) : α := default

-- Array.insert: C++ STL set::insert / vec.push_back 的占位（push 到末尾）
def Array.insert (a : Array α) (v : α) : Array α := a.push v

-- Array.findVal: C++ std::find(vec, x) 的语义（按值找）
-- Lean 4 Array.find? 期望 predicate；这里 wrap 成"按 BEq 找值"
def Array.findVal [BEq α] (a : Array α) (x : α) : Option α := a.find? (· == x)

-- comp 方法占位（已存在于 namespace SparsePolyZp 之内为 UInt64）；
-- 补 SparsePolyZZ / MvPolyZp / MvPolyZZ
def SparsePolyZZ.comp (_f : SparsePolyZZ) : UInt64 := 0
def MvPolyZp.comp (_f : MvPolyZp) : UInt64 := 0
def MvPolyZZ.comp (_f : MvPolyZZ) : UInt64 := 0

-- HMul / HAdd / HSub / HPow 等 Lean 类型类 stub（B2B 测试时细化）
instance : HMul SparsePolyZp SparsePolyZp SparsePolyZp where
  hMul a b := Array.append a b
instance : HAdd SparsePolyZp SparsePolyZp SparsePolyZp where
  hAdd a b := Array.append a b
instance : HSub SparsePolyZp SparsePolyZp SparsePolyZp where
  hSub a b := Array.append a b
-- 用具体 array element 类型避免 abbrev 透明度问题
instance instHMulSparsePolyZZ :
    HMul (Array (UMonomial × Int)) (Array (UMonomial × Int)) (Array (UMonomial × Int)) where
  hMul a b := a ++ b
instance instHAddSparsePolyZZ :
    HAdd (Array (UMonomial × Int)) (Array (UMonomial × Int)) (Array (UMonomial × Int)) where
  hAdd a b := a ++ b
instance instHSubSparsePolyZZ :
    HSub (Array (UMonomial × Int)) (Array (UMonomial × Int)) (Array (UMonomial × Int)) where
  hSub a b := a ++ b
instance : HPow Int UInt64 Int where
  hPow base e := base ^ e.toNat
instance : HPow ZZ UInt64 ZZ where
  hPow base e := base ^ e.toNat

-- Coe Int32 → UInt64 / Int64（Pass 1 把 C++ 字面量识别为 Int32，Lean 端
-- 函数参数常需 UInt64/Int64；自动 Coe 解决 ~5 处 Application mismatch）
instance : Coe Int32 UInt64 where coe n := n.toInt64.toUInt64
instance : Coe Int32 Int64 where coe n := n.toInt64
instance : Coe UInt64 Nat where coe n := n.toNat
instance : Coe Int64 Nat where coe n := n.toNatClampNeg
abbrev SparsePolyZZ := Array (UMonomial × Int)

-- §5a2 迁移：SparsePolyZZ 操作（filterMap 等需要 SparsePolyZZ 已定义）
def SparsePolyZZ.modCoeff (f : SparsePolyZZ) (m : Int) : SparsePolyZZ :=
  f.filterMap (fun (mono, coeff) =>
    let c := coeff % m
    if c != 0 then some (mono, c) else none)

def SparsePolyZZ.compactNonzero (f : SparsePolyZZ) : SparsePolyZZ :=
  f.filter (fun (_, coeff) => coeff != 0)

def SparsePolyZZ.empty : SparsePolyZZ := #[]
def SparsePolyZZ.front! (f : SparsePolyZZ) : UMonomial × Int := f[0]!
def SparsePolyZZ.back! (f : SparsePolyZZ) : UMonomial × Int := f[f.size - 1]!
def SparsePolyZZ.getDeg (f : SparsePolyZZ) : UInt64 := if f.isEmpty then 0 else f[0]!.fst.deg.toUInt64

-- get_deg: 泛型化（C++ side 多模板实例化共用同一 Lean 实现）
-- 适用 SparsePolyZZ / SparsePolyZp 两种容器（结构相同：Array (UMonomial × _)）
-- 返回 Int64（多数 Pass 1 调用点把 get_deg 视为 int64_t / signed comparison 上下文）
def get_deg {α : Type} [Inhabited α] (f : Array (UMonomial × α)) : Int64 :=
  if f.isEmpty then 0 else (f[0]!).fst.deg.toUInt64.toInt64

abbrev LLLMatrix := Array (Array Int)
abbrev HenselNode := Array Int  -- 占位

structure Factorization where
  content : ZZ
  factors : Array (SparsePolyZZ × UInt64)
deriving Inhabited

structure PrimeSelectionResult where
  p : UInt64
  factors : Array SparsePolyZp
  nfactors : UInt64
deriving Inhabited

structure WangLcResult where
  success : Bool
  scaled_factors : Array MvPolyZZ
  lc_targets : Array MvPolyZZ
  delta : ZZ
deriving Inhabited

-- assign：多项式变量代入 poly[var := val]
-- C++ assign(poly, var, val) = 用 val 替代 poly 中的变量 var
opaque assign (poly : α) (var : Variable) (val : β) : α := poly

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
