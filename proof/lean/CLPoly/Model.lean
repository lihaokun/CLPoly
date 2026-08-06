/-
  CLPoly 可信基类 Lean 模型

  定义因式分解模块依赖的所有类型和操作。
  翻译器生成的代码 import 此文件。
  背靠背测试中此文件提供可执行的原语实现。
-/

-- Model.lean 来自 v1 cpp2lean/clpoly_model.lean，借用其隐式变量风格。
-- proof/lean lakefile 全局禁用 autoImplicit，本文件局部启用。
-- 不 import Mathlib：Mathlib 的 `Nat.log` 与 cpp2lean Pass 5 emit 的 `Nat.log`
-- 冲突（已弃用，改 emit `Float.log`）。需要的数学引理见 CLPoly.Math.Bigint。
set_option autoImplicit true

-- ============================================================
-- §1. Zp：Z/pZ 系数
-- ============================================================

structure Zp where
  val : UInt64
  prime : UInt64
deriving Repr, Inhabited, BEq

-- OfNat Zp 0：字面量 0 作为 Zp 时 prime=1 占位（== 比较时只看 val）
-- C++ side `if (zp != 0)` 直 emit 时，0 elaborate 为 Zp（OfNat） + BEq Zp 比较
instance : OfNat Zp 0 where ofNat := { val := 0, prime := 1 }
instance : OfNat Zp 1 where ofNat := { val := 1, prime := 1 }

namespace Zp

def ofInt (v : Int) (p : UInt64) : Zp :=
  let pn : Int := p.toNat
  let r := v % pn
  let r := if r < 0 then r + pn else r
  ⟨r.toNat.toUInt64, p⟩

def ofUInt64 (v p : UInt64) : Zp := ⟨v % p, p⟩

-- Zp 算术：UInt64 溢出安全。a.val + b.val 或 a.val * b.val 可能 > UInt64.size
-- （尤其对大素数 p > 2^32 的情况），统一用 Nat 中间值再 toUInt64
instance : Add Zp where add a b :=
  ⟨((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64, a.prime⟩
instance : Sub Zp where sub a b :=
  ⟨((a.val.toNat + a.prime.toNat - b.val.toNat) % a.prime.toNat).toUInt64, a.prime⟩
instance : Mul Zp where mul a b :=
  ⟨((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64, a.prime⟩
-- Neg: a.val < a.prime 已保证，p - val 不溢出；外层 mod p 处理 val=0 → 0 的边界
instance : Neg Zp where neg a := ⟨(a.prime - a.val) % a.prime, a.prime⟩

-- 扩展欧几里得：gcd(a, b) = a*x + b*y，返回 (gcd, x)
def extGcdAux (old_r r old_s s : Int) : Int × Int :=
  if h : r == 0 then (old_r, old_s)
  else
    let q := old_r / r
    have hr_ne_zero : r ≠ 0 := by
      intro hzero
      apply h
      simp [hzero]
    have h_decr : (old_r - q * r).natAbs < r.natAbs := by
      have hq : old_r - q * r = old_r % r := by
        dsimp [q]
        have h := Int.ediv_add_emod old_r r
        -- h: r * (old_r / r) + old_r % r = old_r
        calc
          old_r - (old_r / r) * r = (r * (old_r / r) + old_r % r) - (old_r / r) * r := by rw [h]
          _ = (r * (old_r / r) + old_r % r) - r * (old_r / r) := by
            rw [Int.mul_comm (old_r / r) r, Int.mul_comm r (old_r / r)]
          _ = old_r % r := by
            rw [Int.sub_eq_add_neg, Int.add_comm (r * (old_r / r)) (old_r % r), Int.add_assoc,
              Int.add_right_neg (r * (old_r / r)), Int.add_zero]
      rw [hq]
      have h_nonneg : 0 ≤ old_r % r := Int.emod_nonneg old_r hr_ne_zero
      have h_mod_lt : old_r % r < (r.natAbs : Int) := Int.emod_lt old_r hr_ne_zero
      have h_natAbs_lt : (old_r % r).natAbs < r.natAbs := by
        have h_eq : ((old_r % r).natAbs : Int) = old_r % r :=
          Int.natAbs_of_nonneg h_nonneg
        have h' : ((old_r % r).natAbs : Int) < (r.natAbs : Int) := by
          rw [h_eq]; exact h_mod_lt
        exact (Int.ofNat_lt.mp h')
      exact h_natAbs_lt
    extGcdAux r (old_r - q * r) s (old_s - q * s)
termination_by r.natAbs

def modInv (a p : UInt64) : UInt64 :=
  if a == 0 then 0
  else
    let (_, s) := extGcdAux (p.toNat : Int) (a.toNat : Int) 0 1
    let r := s % (p.toNat : Int)
    let r := if r < 0 then r + p.toNat else r
    r.toNat.toUInt64

def inv (a : Zp) : Zp := ⟨modInv a.val a.prime, a.prime⟩
def div (a b : Zp) : Zp := a * b.inv

-- 注：Zp.pow 在 §HPow 区块定义（Nat 输入，简单递归）；HPow Zp Int64 Zp
-- 走 toNatClampNeg → 负 exponent 退化为 0 → 返 1，与 C++ pow(Zp,int64_t)
-- (while i>0 模式) 一致。B2B __pow_zp 通过 `a ^ i` 调用。

end Zp

-- Zp 隐式转换（对应 C++ 的 implicit conversion operators）
instance : Coe Zp UInt64 where coe z := z.val
instance : Coe Zp Int where coe z := z.val.toNat
-- 阶段 G+ 修复：Pass 5 emit `Zp != (0 : Int32)` / `assign x (- alpha_j)` 等
-- 字面量保持 Int32 而非 Zp。提供 Coe 让 Lean 自动桥接（占位 prime=1）。
instance : Coe Int32 Zp where coe i := { val := i.toInt64.toUInt64, prime := 1 }
instance : Coe Int Zp where coe i := { val := i.toNat.toUInt64, prime := 1 }

-- ============================================================
-- §2. ZZ：大整数
-- ============================================================

abbrev ZZ := Int

instance : Coe Int64 ZZ where
  coe x := x.toInt

-- ============================================================
-- §3. UMonomial：单变量单项式
-- ============================================================

-- UMonomial.deg 用 Nat（C++ size_t；Lean 4 Nat 同时支持 .toInt32 等转换，
-- 解决 (term_1.fst.deg).toInt32 等场景的 UInt64.toInt32 invalid field 问题）
structure UMonomial where
  deg : Nat
deriving Repr, Inhabited, BEq

/-- Exact unsigned 128-bit machine word used by the dense C++ arithmetic. -/
abbrev UInt128 := BitVec 128

/-- Low 64 bits of an unsigned 128-bit word. -/
def uint128_lo (x : UInt128) : UInt64 := UInt64.ofNat x.toNat

/-- Exact 64-bit count-leading-zeros operation used by `__builtin_clzll`.
The C++ call is guarded against zero; this total definition also assigns the
conventional bit width `64` to zero. -/
def uint64_clz (x : UInt64) : Int32 :=
  Int32.ofNat x.toBitVec.clz.toNat

/-- Clang keeps the dense arithmetic shift count as `uint32_t`; Lean's
machine-word shift instances use `Nat`, so expose the exact conversion. -/
instance : HShiftLeft UInt64 UInt32 UInt64 where
  hShiftLeft x n := x <<< n.toUInt64

instance : HShiftRight UInt64 UInt32 UInt64 where
  hShiftRight x n := x >>> n.toUInt64

/-- Termination fact for the C++ descending `int64_t` loop.

Proof sketch: the loop guard gives `0 ≤ i`.  Hence subtracting one stays in
the signed 64-bit interval (including the endpoint `-1`), so machine
subtraction agrees with integer subtraction.  The nonnegative measure
`max (i + 1) 0` therefore drops by exactly one. -/
theorem int64_sub_one_measure_lt (i : Int64) (h : i >= (0 : Int64)) :
    (((i - (1 : Int64)).toInt + 1).toNat < (i.toInt + 1).toNat) := by
  rw [Int64.toInt_sub]
  have hone : (1 : Int64).toInt = 1 := by decide
  rw [hone]
  rw [Int.bmod_eq_of_le]
  · have hlo : 0 ≤ i.toInt := by
      simpa [Int64.le_iff_toInt_le] using h
    omega
  · have hlo : 0 ≤ i.toInt := by
      simpa [Int64.le_iff_toInt_le] using h
    omega
  · have hi := Int64.toInt_lt i
    omega

/-- Termination fact for the C++ binary-exponentiation loop.

Proof sketch: a positive integer has positive absolute value; division by the
positive constant two commutes with `natAbs`, and natural-number division by
two is strictly smaller than every positive dividend. -/
theorem int_natAbs_ediv_two_lt (e : Int) (h : 0 < e) :
    (e / 2).natAbs < e.natAbs := by
  rw [Int.natAbs_ediv_of_nonneg (Int.le_of_lt h)]
  have he : 0 < e.natAbs := Int.natAbs_pos.mpr (by omega)
  simpa using Nat.div_lt_self he (by decide : 1 < 2)

/-- State carried by the C++ `dense_upoly_zp` implementation.  This is only
the L1 representation type; its constructors and algorithms are translated
from the corresponding C++ member bodies into generated namespaces. -/
structure DenseUPolyZp where
  _coeffs : Array UInt64 := #[]
  _p : UInt64 := 0
  _ninv : UInt64 := 0
  _norm : UInt32 := 0
deriving Repr, Inhabited, BEq

/-- Exact L1 representation of `dense_upoly_zp::word3`, the three-limb lazy
accumulator used by the C++ multiplication and division routines. -/
structure Word3 where
  lo : UInt64 := 0
  mid : UInt64 := 0
  hi : UInt64 := 0
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
-- Bug fix (B2B A3)：modular reduction 后 new_val 可能为 0（如 d/dx(x^p)=p*x^{p-1}=0
-- 在 F_p 上），需 filter 掉 0-coef 项，否则 SparsePolyZp 失常规化
def derivative (f : SparsePolyZp) : SparsePolyZp :=
  if f.isEmpty then #[]
  else
    let p := f[0]!.snd.prime
    f.filterMap (fun (m, c) =>
      if m.deg == 0 then none
      else
        let new_val := c.val * m.deg.toUInt64 % p
        if new_val = 0 then none
        else some (⟨m.deg - 1⟩, ⟨new_val, p⟩))

-- divmod / gcd 实现移到 §5d.1 之后（需要 SparsePolyZp 算术 instance 已注册）

-- 比较器：单变量多项式的单项式序（降幂排列）
-- C++ `.comp()` 返回比较器（lex_<less> 类对象）；Lean 用 Unit 占位
-- (MonomialOrder = Unit abbrev 见 §5d，此处先用 Unit 避免 forward ref)
def comp (_f : SparsePolyZp) : Unit := ()

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

-- 同时返回随机值 + 推进后的 seed（C++ side `dist(rng)` 的 ref-out 语义）
-- Pass 5 把 `dist(rng)` translate 为 Rng.next_advance；Pass 2b（重跑）
-- 把 idx=0 的 rng 当作 ref-out，destructure 为 (val, new_rng) tuple
def next_advance (seed upper : UInt64) : UInt64 × UInt64 :=
  (next seed upper, step seed)

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
def begin_ (_m : StdMap K V) : Nat := 0

end StdMap

-- ============================================================
-- §5c. 其他辅助类型
-- ============================================================

-- 多变量多项式（占位）
-- 阶段 E：Monomial / MvMonomial / Variable 类型统一
-- C++ side: `Monomial = basic_monomial<lex_<less>>` = vector<pair<variable, exponent>>
-- variable 是 UInt64 ID（CLPoly clpoly::variable handle），exponent 是 Int64
abbrev Variable := UInt64
abbrev Monomial := Array (Variable × Int64)
-- MvMonomial 与 Monomial 同义（Pass 1 不同 context 偶尔产生 MvMonomial 别名）
abbrev MvMonomial := Monomial

-- §3a Monomial 比较 / normalize（lex order，与 C++ lex_<less> 一致）
namespace Monomial

-- 规范化：合并同 var 项 + 按 var id 升序 + 剔除 exp = 0
def normalize (m : Monomial) : Monomial :=
  let grouped : Monomial := m.foldl (fun acc (v, e) =>
    match acc.findIdx? (fun t => t.fst == v) with
    | some idx => acc.modify idx (fun (v', e') => (v', e' + e))
    | none => acc.push (v, e)) #[]
  let nonZero : Monomial := grouped.filter (fun (_, e) => e ≠ 0)
  nonZero.qsort (fun a b => a.fst < b.fst)

-- Monomial 乘法：exp 相加（合并同 var）
def mul (m1 m2 : Monomial) : Monomial :=
  Monomial.normalize (m1 ++ m2)

-- Lex 比较 (small var id = high priority)
-- 假设输入已 normalize（按 var id 升序、无 0 exp）
partial def ltAux : List (Variable × Int64) → List (Variable × Int64) → Bool
  | [], [] => false
  | [], _ :: _ => true   -- m1 = 0 < m2 在该 var 上 > 0
  | _ :: _, [] => false  -- m1 > m2
  | (v1, e1) :: rest1, (v2, e2) :: rest2 =>
    if v1 == v2 then
      if e1 == e2 then ltAux rest1 rest2
      else e1 < e2
    else if v1 < v2 then
      false  -- m1 在更小 var 上有 > 0 exp，m2 那位 = 0 → m1 > m2
    else
      true   -- 反之

def lt (m1 m2 : Monomial) : Bool :=
  Monomial.ltAux m1.toList m2.toList

def eq (m1 m2 : Monomial) : Bool :=
  m1 == m2

-- 总度数（所有 exp 之和的 toNat）
def totalDeg (m : Monomial) : Nat :=
  m.foldl (fun acc (_, e) => acc + e.toNatClampNeg) 0

end Monomial
def Monomial.empty : Monomial := #[]
-- HasPolyMk：1-arg ctor 派发
-- Pass 5 把 C++ 多种 1-arg ctor（copy / Poly(comp_ptr) / Poly(elem) / Variable("x")）
-- 都映射到 .mk。用 typeclass 区分：
--   - 输入类型 = 输出类型（copy ctor） → 返回输入（避免 silent identity bug）
--   - 否则（comp_ptr/MonomialOrder/String 等元素） → 返回 default（空 Array / 0）
class HasPolyMk (β α : Type) where
  mkImpl : α → β
-- 同型 copy ctor（高优先级）
instance {β : Type} : HasPolyMk β β where mkImpl x := x
-- 异型兜底：返回 default
instance (priority := 0) {β α : Type} [Inhabited β] : HasPolyMk β α where mkImpl _ := default
def Monomial.mk {α : Type} (x : α) [HasPolyMk Monomial α] : Monomial := HasPolyMk.mkImpl x
def MvMonomial.empty : MvMonomial := #[]
-- 多变量多项式：内部元素的第一槽是 Monomial，与 Pass 1 推断的 (Monomial × ZZ)
-- 一致；ZZ = Int / Zp 同样为系数类型
abbrev MvPolyZZ := Array (Monomial × Int)
abbrev MvPolyZp := Array (Monomial × Zp)
abbrev PolyZp := MvPolyZp
abbrev PolyZZ := MvPolyZZ
abbrev PolyQQ := Array (Monomial × Rat)

-- §3b MvPoly normalize：normalize 每个 monomial、合并同 monomial、排序
-- MvPolyZp normalization
def MvPolyZp.normalization (f : MvPolyZp) : MvPolyZp :=
  let normMonos : MvPolyZp := f.map (fun (m, c) => (Monomial.normalize m, c))
  let grouped : MvPolyZp := normMonos.foldl (fun acc term =>
    match acc.findIdx? (fun t => Monomial.eq t.fst term.fst) with
    | some idx => acc.modify idx (fun (m, c) => (m, c + term.snd))
    | none => acc.push term) #[]
  let nonZero := grouped.filter (fun t => t.snd.val ≠ 0)
  nonZero.qsort (fun a b => Monomial.lt b.fst a.fst)
def MvPolyZZ.empty : MvPolyZZ := #[]
def MvPolyZp.empty : MvPolyZp := #[]
-- C++ 的 Poly(comp_t) / Poly(const Poly&) ctor 都映射到 .mk；
-- 用泛型 input → 空 Poly 占位（语义在 B2B 层细化）
def MvPolyZp.mk {α : Type} (x : α) [HasPolyMk MvPolyZp α] : MvPolyZp := HasPolyMk.mkImpl x
def MvPolyZZ.normalization (f : MvPolyZZ) : MvPolyZZ :=
  let normMonos : MvPolyZZ := f.map (fun (m, c) => (Monomial.normalize m, c))
  let grouped : MvPolyZZ := normMonos.foldl (fun acc term =>
    match acc.findIdx? (fun t => Monomial.eq t.fst term.fst) with
    | some idx => acc.modify idx (fun (m, c) => (m, c + term.snd))
    | none => acc.push term) #[]
  let nonZero := grouped.filter (fun t => t.snd ≠ 0)
  nonZero.qsort (fun a b => Monomial.lt b.fst a.fst)
def MvPolyZZ.mk {α : Type} (x : α) [HasPolyMk MvPolyZZ α] : MvPolyZZ := HasPolyMk.mkImpl x
-- 通用 stub（与 SparsePolyZZ 解耦，无前向引用）
def __write__ (_old : α) (new : α) : α := new

-- ZZ.toBool: ZZ → Bool（C++ side `if (zz)` 的语义：非零为 true）
def ZZ.toBool (z : ZZ) : Bool := z != 0

-- LambdaRef: Pass 3 lifted lambda 在 caller-arg 位置的 placeholder 类型
-- (Pass 3 不保留具体函数签名)。语义占位为 Unit。
abbrev LambdaRef := Unit

-- 模板/typedef 残留 alias（Pass 1 未替换为具体实例化）
-- Rng: std::mt19937 — 状态用 UInt64 种子（与 §5a 命名一致）
abbrev Rng := UInt64
-- UniformIntDist: std::uniform_int_distribution<> — 占位（仅承载 upper bound）
abbrev UniformIntDist := UInt64
-- Poly: typedef polynomial_<ZZ, lex_<less>> 的别名（PolyZp/PolyZZ 已在 §5c 声明）
abbrev Poly := MvPolyZZ
-- QQ: typedef rational — 用 Rat 占位（C++ 用任意精度有理数）
abbrev QQ := Rat
-- MonomialOrder: 单项式序 tag —— 仅类型层占位（C++ side `lex_<var_order>`）
-- 命名避开 Mathlib `Lex` (Order.Synonym) 名冲突。
abbrev MonomialOrder := Unit
-- §5d.0 多项式 GCD / divmod typeclass dispatch
-- 通用兜底（Inhabited α）+ SparsePolyZp 特化（在 SparsePolyZp.gcd / divmod 之后注册）
class HasPolyGCD (α : Type) where
  polyGCD : α → α → α

class HasPolyDivmod (α : Type) where
  polyDivmod : α → α → α × α

-- Generic dispatch（与 cpp2lean Pass 5 emit 一致）
def polynomial_GCD {α : Type} [HasPolyGCD α] (a b : α) : α := HasPolyGCD.polyGCD a b

-- 4-arg pair_vec_div: (q_out, v1, v2, comp) → q（C++ basic.hh:568）
def pair_vec_div {α : Type} [HasPolyDivmod α] {β : Type}
    (_q_out v1 v2 : α) (_comp : β) : α :=
  (HasPolyDivmod.polyDivmod v1 v2).fst

-- 5-arg pair_vec_div5: (q_out, r_out, v1, v2, comp) → (q, r)（C++ basic.hh:698）
def pair_vec_div5 {α : Type} [HasPolyDivmod α] {β : Type}
    (_q_out _r_out v1 v2 : α) (_comp : β) : α × α :=
  HasPolyDivmod.polyDivmod v1 v2

-- 兜底 instance（任意 Inhabited 类型 → default）
instance (priority := 0) {α : Type} [Inhabited α] : HasPolyGCD α where
  polyGCD _ _ := default

instance (priority := 0) {α : Type} [Inhabited α] : HasPolyDivmod α where
  polyDivmod _ _ := (default, default)

-- 4-arg Bezout EEA 形式：(a, b, s_out, t_out) → (gcd, s, t) 满足 a*s + b*t = g
-- typeclass dispatch；SparsePolyZp 特化在 SparsePolyZp.gcd / divmod 之后注册
class HasPolyGCDEEA (α : Type) where
  polyGCDEEA : α → α → α × α × α

def polynomial_GCD_eea {α : Type} [HasPolyGCDEEA α]
    (a b _s_old _t_old : α) : α × α × α :=
  HasPolyGCDEEA.polyGCDEEA a b

instance (priority := 0) {α : Type} [Inhabited α] : HasPolyGCDEEA α where
  polyGCDEEA _ _ := (default, default, default)

-- Array.insert: C++ STL set::insert / vec.push_back 的占位（push 到末尾）
def Array.insert (a : Array α) (v : α) : Array α := a.push v

-- Pass 5 vec.erase(idx) 偶尔 emit `Array.erase a idx`（idx 是 Int64/Int32），
-- Lean Array.erase 是 by-value (a : Array α) (v : α) → Array α，类型不匹配。
-- 提供 fallback Array.erase 处理（Pass 5 emit 时若 v 不是 α 类型，Lean 会找
-- 这个 instance）
def Array.eraseAny {α : Type} [Inhabited α] (a : Array α) {β : Type} (_v : β) : Array α :=
  if a.size > 0 then a.eraseIdxIfInBounds 0 else a

-- Array.findVal: C++ std::find(vec, x) 的语义（按值找）
-- Lean 4 Array.find? 期望 predicate；这里 wrap 成"按 BEq 找值"
def Array.findVal [BEq α] (a : Array α) (x : α) : Option α := a.find? (· == x)

-- comp 方法占位（已存在于 namespace SparsePolyZp 之内为 UInt64）；
-- 补 SparsePolyZZ / MvPolyZp / MvPolyZZ
-- C++ `Poly.comp()` 返回比较器对象（lex_<less>）而非 UInt64
def SparsePolyZZ.comp (_f : SparsePolyZZ) : MonomialOrder := ()
def MvPolyZp.comp (_f : MvPolyZp) : MonomialOrder := ()
def MvPolyZZ.comp (_f : MvPolyZZ) : MonomialOrder := ()

-- §5e. C++ 容器 mutate 占位 + degree typeclass
-- ============================================================

-- vec.clear() 的占位：忽略 receiver，返回空 Array
def Array.clearVec {α : Type} (_ : Array α) : Array α := #[]

-- Pass 5 vec.erase(begin + j) 转为 `Array.eraseIdx arr idx` —— 包装版本
-- 接受任意 int-like 索引（Int32/Int64/UInt64）并 toNat 转为 Lean Nat
class ToIdxNat (α : Type) where toIdxNat : α → Nat
instance : ToIdxNat Nat where toIdxNat n := n
instance : ToIdxNat Int32 where toIdxNat i := i.toNatClampNeg
instance : ToIdxNat Int64 where toIdxNat i := i.toNatClampNeg
instance : ToIdxNat UInt32 where toIdxNat u := u.toNat
instance : ToIdxNat UInt64 where toIdxNat u := u.toNat
instance : ToIdxNat Int where toIdxNat i := i.toNat

def Array.eraseIdx' {α : Type} (a : Array α) {β : Type} [ToIdxNat β]
    (i : β) : Array α :=
  a.eraseIdxIfInBounds (ToIdxNat.toIdxNat i)

-- vec.assign(n, val) 的占位：忽略 receiver，按 (n, val) 复制
def Array.replicateMut {α : Type} (_ : Array α) (n : Nat) (v : α) : Array α :=
  Array.replicate n v

-- vec.empty() 谓词占位（Mv* 别名转发到 Array.isEmpty）
def MvPolyZZ.isEmpty (f : MvPolyZZ) : Bool := Array.isEmpty f
def MvPolyZp.isEmpty (f : MvPolyZp) : Bool := Array.isEmpty f
def MvMonomial.isEmpty (m : MvMonomial) : Bool := Array.isEmpty m

-- Rng.mk: 用整数种子构造 RNG（C++: std::mt19937(seed)）
def Rng.mk (seed : Int32) : Rng := seed.toInt64.toUInt64

-- Array.sort: C++ std::sort with comparator
-- cmp a b 返回 true 表示 a 应排在 b 之前（与 std::less 一致）
def Array.sort {α : Type} (a : Array α) (cmp : α → α → Bool) : Array α := a.qsort cmp

-- C++ 自由函数 degree(poly) — 多态（lambda 比较器里用）。
-- 用 typeclass 解决"未限定 degree"调用的多重 receiver 类型问题。
class HasDegree (α : Type) where
  degree : α → UInt64

-- HasDegree MvPoly：取所有 monomial 的 max total degree
instance : HasDegree MvPolyZZ where
  degree (f : MvPolyZZ) : UInt64 :=
    f.foldl (fun acc (mono, _) =>
      max acc (Monomial.totalDeg mono).toUInt64) 0

instance : HasDegree MvPolyZp where
  degree (f : MvPolyZp) : UInt64 :=
    f.foldl (fun acc (mono, _) =>
      max acc (Monomial.totalDeg mono).toUInt64) 0
-- SparsePolyZZ HasDegree instance 在 abbrev 后再加（见 §5c）

def degree {α : Type} [HasDegree α] (a : α) : UInt64 := HasDegree.degree a
-- 2-arg degree: poly + main var → 关于该 var 的次数（C++ degree(poly, var)）
def degree2 {α : Type} [HasDegree α] (a : α) (_var : Variable) : Int64 :=
  (HasDegree.degree a).toInt64

-- C++ free functions: is_number(poly) / get_variables(poly)
class IsNumber (α : Type) where
  isNumber : α → Bool

-- IsNumber：常数多项式（空 / 单项 0 度 / 所有 monomial 都为空）
instance : IsNumber MvPolyZZ where
  isNumber (f : MvPolyZZ) : Bool :=
    f.all (fun (mono, _) => mono.all (fun (_, e) => e = 0))

instance : IsNumber MvPolyZp where
  isNumber (f : MvPolyZp) : Bool :=
    f.all (fun (mono, _) => mono.all (fun (_, e) => e = 0))
-- SparsePolyZZ instance 在其 abbrev 之后再定义（见 §5c）

def is_number {α : Type} [IsNumber α] (a : α) : Bool := IsNumber.isNumber a

class GetVariables (α : Type) where
  vars : α → Array (Variable × Int64)

-- GetVariables：收集 f 中出现的所有 (var, max_exp) 对（去重）
instance : GetVariables MvPolyZZ where
  vars (f : MvPolyZZ) : Array (Variable × Int64) :=
    f.foldl (fun acc (mono, _) =>
      mono.foldl (fun a (v, e) =>
        match a.findIdx? (fun t => t.fst == v) with
        | some idx => a.modify idx (fun (v', e') => (v', max e' e))
        | none => a.push (v, e)) acc) #[]

instance : GetVariables MvPolyZp where
  vars (f : MvPolyZp) : Array (Variable × Int64) :=
    f.foldl (fun acc (mono, _) =>
      mono.foldl (fun a (v, e) =>
        match a.findIdx? (fun t => t.fst == v) with
        | some idx => a.modify idx (fun (v', e') => (v', max e' e))
        | none => a.push (v, e)) acc) #[]

def get_variables {α : Type} [GetVariables α] (a : α) : Array (Variable × Int64) :=
  GetVariables.vars a

-- MvPolyZZ.front! / MvPolyZp.front!：取首项（mono × coeff）
def MvPolyZZ.front! (f : MvPolyZZ) : (Monomial × Int) := f[0]!
def MvPolyZp.front! (f : MvPolyZp) : (Monomial × Zp) := f[0]!

-- ============================================================
-- §5d.1 SparsePolyZp 算术：真实多项式加减乘
--
-- 不变量（与 C++ upolynomial_<Zp> 一致）：
--   1. 按 deg 降序排序
--   2. 无重复 deg
--   3. 无 0 系数
-- 算法基于归并；spec 在 CLPoly.Math.Univariate（Mathlib Polynomial bridge）。
-- ============================================================

namespace SparsePolyZp

-- 归并加：合并两个降序列表，同 deg 项相加，删 0
-- termination_by: xs.length + ys.length 严格递减
def mergeAdd : List (UMonomial × Zp) → List (UMonomial × Zp) → List (UMonomial × Zp)
  | [], gs => gs
  | f :: fs, [] => f :: fs
  | f :: fs, g :: gs =>
    if f.fst.deg > g.fst.deg then
      f :: mergeAdd fs (g :: gs)
    else if f.fst.deg < g.fst.deg then
      g :: mergeAdd (f :: fs) gs
    else
      let s := f.snd + g.snd
      if s.val = 0 then mergeAdd fs gs
      else (f.fst, s) :: mergeAdd fs gs
termination_by xs ys => xs.length + ys.length

-- 加法：实现 +
def addImpl (f g : SparsePolyZp) : SparsePolyZp :=
  (mergeAdd f.toList g.toList).toArray

-- 取负：每项系数取负
def negImpl (f : SparsePolyZp) : SparsePolyZp :=
  f.map (fun (m, c) => (m, -c))

-- 减法：f + (-g)
def subImpl (f g : SparsePolyZp) : SparsePolyZp :=
  addImpl f (negImpl g)

-- 标量乘 monomial：c * x^d * f
def scaleByMonomial (m : UMonomial) (c : Zp) (f : SparsePolyZp) : SparsePolyZp :=
  if c.val = 0 then #[]
  else f.filterMap (fun (mf, cf) =>
    let prod := c * cf
    if prod.val = 0 then none
    else some (UMonomial.mk (m.deg + mf.deg), prod))

-- 乘法：朴素 O(n*m) 累加
-- f * g = Σ_i (f[i].snd * x^f[i].fst) * g
def mulImpl (f g : SparsePolyZp) : SparsePolyZp :=
  f.foldl (fun acc (mf, cf) => addImpl acc (scaleByMonomial mf cf g)) #[]

end SparsePolyZp

instance : HAdd SparsePolyZp SparsePolyZp SparsePolyZp where
  hAdd a b := SparsePolyZp.addImpl a b
instance : HSub SparsePolyZp SparsePolyZp SparsePolyZp where
  hSub a b := SparsePolyZp.subImpl a b
instance : Neg SparsePolyZp where
  neg a := SparsePolyZp.negImpl a
instance : HMul SparsePolyZp SparsePolyZp SparsePolyZp where
  hMul a b := SparsePolyZp.mulImpl a b

namespace SparsePolyZp

-- §5d.2 多项式长除法 + GCD（须在 HMul / HSub instance 之后定义）
-- 不变量：g 非空（除数 ≠ 0），所有 Zp 共享 prime（WellFormed）

-- 稀疏多项式「按 deg 严格降序」不变量：首项在 index 0。
-- mergeAdd/subImpl 等运算隐含假定输入降序；此谓词把该不变量显式化，
-- 用于长除法主循环的良基终止（首项抵消 ⇒ deg 严格下降）。
-- Model.lean 不 import Mathlib，故用 Bool 递归实现（天然可判定，供 divmod 依赖 if）。
def sortedListB : List (UMonomial × Zp) → Bool
  | [] => true
  | [_] => true
  | a :: b :: rest => (b.fst.deg < a.fst.deg) && sortedListB (b :: rest)

def Sorted (f : SparsePolyZp) : Prop := sortedListB f.toList = true

instance (f : SparsePolyZp) : Decidable (Sorted f) :=
  decEq (sortedListB f.toList) true

-- ── 有序性基础设施：sortedListB 的结构性质 + 各数组运算保持有序 ──

/-- 严格降序 ⇔ 头严格大于尾中所有元素 ∧ 尾亦有序。 -/
theorem sortedListB_iff (a : UMonomial × Zp) (rest : List (UMonomial × Zp)) :
    sortedListB (a :: rest) = true ↔
      (∀ x ∈ rest, x.fst.deg < a.fst.deg) ∧ sortedListB rest = true := by
  induction rest generalizing a with
  | nil => simp [sortedListB]
  | cons b rest' ih =>
    simp only [sortedListB, Bool.and_eq_true, decide_eq_true_eq]
    constructor
    · rintro ⟨hba, hsb⟩
      refine ⟨?_, hsb⟩
      intro x hx
      rcases List.mem_cons.mp hx with rfl | hxr
      · exact hba
      · exact Nat.lt_trans (((ih b).mp hsb).1 x hxr) hba
    · rintro ⟨hall, hsb⟩
      exact ⟨hall b (by simp), hsb⟩

/-- mergeAdd 的元素度数被两输入的公共上界所界。 -/
theorem mergeAdd_lt_all (d : Nat) : ∀ (xs ys : List (UMonomial × Zp)),
    (∀ x ∈ xs, x.fst.deg < d) → (∀ y ∈ ys, y.fst.deg < d) →
    (∀ z ∈ mergeAdd xs ys, z.fst.deg < d) := by
  intro xs ys
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys => intro _ hy z hz; rw [mergeAdd] at hz; exact hy z hz
  | case2 f fs => intro hx _ z hz; rw [mergeAdd] at hz; exact hx z hz
  | case3 f fs g gs hfg ih =>
    intro hx hy z hz
    rw [mergeAdd, if_pos hfg] at hz
    rcases List.mem_cons.mp hz with rfl | hz'
    · exact hx z (by simp)
    · exact ih (fun x hx' => hx x (List.mem_cons_of_mem _ hx')) hy z hz'
  | case4 f fs g gs hfg hfg2 ih =>
    intro hx hy z hz
    rw [mergeAdd, if_neg hfg, if_pos hfg2] at hz
    rcases List.mem_cons.mp hz with rfl | hz'
    · exact hy z (by simp)
    · exact ih hx (fun y hy' => hy y (List.mem_cons_of_mem _ hy')) z hz'
  | case5 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy z hz
    have h_eq : mergeAdd (f :: fs) (g :: gs) = mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]
      exact if_pos hs
    rw [h_eq] at hz
    exact ih (fun x hx' => hx x (List.mem_cons_of_mem _ hx'))
      (fun y hy' => hy y (List.mem_cons_of_mem _ hy')) z hz
  | case6 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy z hz
    have h_eq : mergeAdd (f :: fs) (g :: gs) = (f.fst, f.snd + g.snd) :: mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]
      exact if_neg hs
    rw [h_eq] at hz
    rcases List.mem_cons.mp hz with rfl | hz'
    · exact hx f (by simp)
    · exact ih (fun x hx' => hx x (List.mem_cons_of_mem _ hx'))
        (fun y hy' => hy y (List.mem_cons_of_mem _ hy')) z hz'

/-- mergeAdd 保持严格降序有序。 -/
theorem mergeAdd_sorted : ∀ (xs ys : List (UMonomial × Zp)),
    sortedListB xs = true → sortedListB ys = true → sortedListB (mergeAdd xs ys) = true := by
  intro xs ys
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys => intro _ hy; rw [mergeAdd]; exact hy
  | case2 f fs => intro hx _; rw [mergeAdd]; exact hx
  | case3 f fs g gs hfg ih =>
    intro hx hy
    rw [mergeAdd, if_pos hfg, sortedListB_iff]
    have hx' := (sortedListB_iff f fs).mp hx
    refine ⟨?_, ih hx'.2 hy⟩
    apply mergeAdd_lt_all f.fst.deg
    · exact hx'.1
    · intro y hy'; rcases List.mem_cons.mp hy' with rfl | hy''
      · exact hfg
      · exact Nat.lt_trans (((sortedListB_iff g gs).mp hy).1 y hy'') hfg
  | case4 f fs g gs hfg hfg2 ih =>
    intro hx hy
    rw [mergeAdd, if_neg hfg, if_pos hfg2, sortedListB_iff]
    have hy' := (sortedListB_iff g gs).mp hy
    refine ⟨?_, ih hx hy'.2⟩
    apply mergeAdd_lt_all g.fst.deg
    · intro x hx'; rcases List.mem_cons.mp hx' with rfl | hx''
      · exact hfg2
      · exact Nat.lt_trans (((sortedListB_iff f fs).mp hx).1 x hx'') hfg2
    · exact hy'.1
  | case5 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]
      exact if_pos hs
    rw [h_eq]
    exact ih ((sortedListB_iff f fs).mp hx).2 ((sortedListB_iff g gs).mp hy).2
  | case6 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = (f.fst, f.snd + g.snd) :: mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]
      exact if_neg hs
    rw [h_eq, sortedListB_iff]
    have hx' := (sortedListB_iff f fs).mp hx
    have hy' := (sortedListB_iff g gs).mp hy
    have hfe : f.fst.deg = g.fst.deg := by omega
    refine ⟨?_, ih hx'.2 hy'.2⟩
    apply mergeAdd_lt_all f.fst.deg
    · exact hx'.1
    · intro y hy''; rw [hfe]; exact hy'.1 y hy''

/-- 保持 .fst 的 map 不改变有序性（negImpl 用：只改系数、不改单项式）。 -/
theorem sortedListB_map_fst (F : (UMonomial × Zp) → (UMonomial × Zp))
    (hF : ∀ x, (F x).fst = x.fst) (l : List (UMonomial × Zp)) :
    sortedListB (l.map F) = sortedListB l := by
  induction l with
  | nil => rfl
  | cons a t iha =>
    cases t with
    | nil => rfl
    | cons b t' =>
      simp only [List.map_cons, sortedListB, hF]
      rw [← List.map_cons, iha]

/-- subImpl 保持严格降序有序（addImpl=mergeAdd, negImpl 只改系数）。 -/
theorem subImpl_sorted (f g : SparsePolyZp) (hf : Sorted f) (hg : Sorted g) :
    Sorted (subImpl f g) := by
  have hng : sortedListB (negImpl g).toList = true := by
    unfold negImpl
    rw [Array.toList_map]
    exact (sortedListB_map_fst _ (by rintro ⟨m, c⟩; rfl) g.toList).trans hg
  -- (mergeAdd ...).toArray.toList = mergeAdd ... 由 toList_toArray (= rfl) 定义式相等
  show sortedListB (subImpl f g).toList = true
  unfold subImpl addImpl
  exact mergeAdd_sorted _ _ hf hng

/-- filterMap 且输出度数 = sh + 输入度数 时：保持「度数 < 界」与有序性。
    （scaleByMonomial 用：单项式乘法把每个度数平移 sh = m.deg，保序。） -/
theorem sortedListB_filterMapShift (sh : Nat) (G : (UMonomial × Zp) → Option (UMonomial × Zp))
    (hG : ∀ x y, G x = some y → y.fst.deg = sh + x.fst.deg) (l : List (UMonomial × Zp)) :
    (∀ d, (∀ x ∈ l, x.fst.deg < d) → (∀ z ∈ l.filterMap G, z.fst.deg < sh + d)) ∧
    (sortedListB l = true → sortedListB (l.filterMap G) = true) := by
  induction l with
  | nil =>
    refine ⟨?_, ?_⟩
    · intro d _ z hz; simp at hz
    · intro _; rfl
  | cons a t ih =>
    obtain ⟨ih_lt, ih_sorted⟩ := ih
    refine ⟨?_, ?_⟩
    · intro d hlt z hz
      cases hGa : G a with
      | none =>
        rw [List.filterMap_cons_none hGa] at hz
        exact ih_lt d (fun x hx => hlt x (List.mem_cons_of_mem _ hx)) z hz
      | some b =>
        rw [List.filterMap_cons_some hGa] at hz
        rcases List.mem_cons.mp hz with hzb | hz'
        · rw [hzb, hG a b hGa]; have := hlt a (by simp); omega
        · exact ih_lt d (fun x hx => hlt x (List.mem_cons_of_mem _ hx)) z hz'
    · intro hsort
      cases hGa : G a with
      | none =>
        rw [List.filterMap_cons_none hGa]
        exact ih_sorted ((sortedListB_iff a t).mp hsort).2
      | some b =>
        rw [List.filterMap_cons_some hGa, sortedListB_iff]
        have hst := (sortedListB_iff a t).mp hsort
        refine ⟨?_, ih_sorted hst.2⟩
        intro z hz
        rw [hG a b hGa]
        exact ih_lt a.fst.deg hst.1 z hz

/-- scaleByMonomial（单项式乘法）保持有序：度数整体平移 m.deg。 -/
theorem scaleByMonomial_sorted (m : UMonomial) (c : Zp) (f : SparsePolyZp) (hf : Sorted f) :
    Sorted (scaleByMonomial m c f) := by
  unfold scaleByMonomial
  split
  · rfl
  · show sortedListB (Array.filterMap _ f).toList = true
    rw [Array.toList_filterMap]
    refine (sortedListB_filterMapShift m.deg _ ?_ f.toList).2 hf
    intro x y hxy
    obtain ⟨mf, cf⟩ := x
    simp only at hxy
    split at hxy
    · exact absurd hxy (by simp)
    · injection hxy with hxy'; rw [← hxy']

/-- addImpl 保持有序（= mergeAdd）。 -/
theorem addImpl_sorted (f g : SparsePolyZp) (hf : Sorted f) (hg : Sorted g) :
    Sorted (addImpl f g) := by
  show sortedListB (addImpl f g).toList = true
  unfold addImpl
  exact mergeAdd_sorted _ _ hf hg

/-- 单项式（单元素数组）乘多项式保持有序：mulImpl 单元素折叠 = scaleByMonomial。 -/
theorem single_mul_sorted (m : UMonomial) (c : Zp) (g : SparsePolyZp) (hg : Sorted g) :
    Sorted ((#[(m, c)] : SparsePolyZp) * g) := by
  have h_eq : (#[(m, c)] : SparsePolyZp) * g = addImpl #[] (scaleByMonomial m c g) := rfl
  rw [h_eq]
  exact addImpl_sorted _ _ rfl (scaleByMonomial_sorted m c g hg)

-- ── 约化+同素数不变量：所有系数 val < q 且 prime = q（长除法首项抵消所需） ──

def reducedListB (q : UInt64) : List (UMonomial × Zp) → Bool
  | [] => true
  | a :: rest => (a.snd.val < q) && (a.snd.prime == q) && reducedListB q rest

def ReducedB (f : SparsePolyZp) (q : UInt64) : Prop := reducedListB q f.toList = true

instance (f : SparsePolyZp) (q : UInt64) : Decidable (ReducedB f q) :=
  decEq (reducedListB q f.toList) true

theorem reducedListB_cons (q : UInt64) (a : UMonomial × Zp) (rest : List (UMonomial × Zp)) :
    reducedListB q (a :: rest) = true ↔
      (a.snd.val < q ∧ a.snd.prime = q) ∧ reducedListB q rest = true := by
  simp only [reducedListB, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, and_assoc]

/-- Zp 加法结果的 val 界（同素数 q>0 时 < q）。 -/
theorem Zp_add_reduced (a b : Zp) (hq : 0 < a.prime.toNat) :
    (a + b).val < a.prime ∧ (a + b).prime = a.prime := by
  refine ⟨?_, rfl⟩
  show ((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64 < a.prime
  have hlt : (a.val.toNat + b.val.toNat) % a.prime.toNat < a.prime.toNat := Nat.mod_lt _ hq
  have h2 : (a.val.toNat + b.val.toNat) % a.prime.toNat < UInt64.size :=
    Nat.lt_trans hlt a.prime.toNat_lt_size
  rw [UInt64.lt_iff_toNat_lt,
    show ((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64.toNat
      = (a.val.toNat + b.val.toNat) % a.prime.toNat from by simp [h2]]
  exact hlt

/-- mergeAdd 保持约化+同素数不变量。 -/
theorem reducedListB_mergeAdd (q : UInt64) (hq : 0 < q.toNat) : ∀ (xs ys : List (UMonomial × Zp)),
    reducedListB q xs = true → reducedListB q ys = true → reducedListB q (mergeAdd xs ys) = true := by
  intro xs ys
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys => intro _ hy; rw [mergeAdd]; exact hy
  | case2 f fs => intro hx _; rw [mergeAdd]; exact hx
  | case3 f fs g gs hfg ih =>
    intro hx hy
    rw [mergeAdd, if_pos hfg, reducedListB_cons]
    have hx' := (reducedListB_cons q f fs).mp hx
    exact ⟨hx'.1, ih hx'.2 hy⟩
  | case4 f fs g gs hfg hfg2 ih =>
    intro hx hy
    rw [mergeAdd, if_neg hfg, if_pos hfg2, reducedListB_cons]
    have hy' := (reducedListB_cons q g gs).mp hy
    exact ⟨hy'.1, ih hx hy'.2⟩
  | case5 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]; exact if_pos hs
    rw [h_eq]
    exact ih ((reducedListB_cons q f fs).mp hx).2 ((reducedListB_cons q g gs).mp hy).2
  | case6 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = (f.fst, f.snd + g.snd) :: mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]; exact if_neg hs
    rw [h_eq, reducedListB_cons]
    have hx' := (reducedListB_cons q f fs).mp hx
    have hprime : f.snd.prime = q := hx'.1.2
    have hqp : 0 < f.snd.prime.toNat := by rw [hprime]; exact hq
    have hadd := Zp_add_reduced f.snd g.snd hqp
    refine ⟨⟨?_, ?_⟩, ih hx'.2 ((reducedListB_cons q g gs).mp hy).2⟩
    · rw [← hprime]; exact hadd.1
    · rw [hadd.2]; exact hprime

/-- Zp 乘法结果的 val 界（同素数 q>0 时 < q）。 -/
theorem Zp_mul_reduced (a b : Zp) (hq : 0 < a.prime.toNat) :
    (a * b).val < a.prime ∧ (a * b).prime = a.prime := by
  refine ⟨?_, rfl⟩
  show ((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64 < a.prime
  have hlt : (a.val.toNat * b.val.toNat) % a.prime.toNat < a.prime.toNat := Nat.mod_lt _ hq
  have h2 : (a.val.toNat * b.val.toNat) % a.prime.toNat < UInt64.size :=
    Nat.lt_trans hlt a.prime.toNat_lt_size
  rw [UInt64.lt_iff_toNat_lt,
    show ((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64.toNat
      = (a.val.toNat * b.val.toNat) % a.prime.toNat from by simp [h2]]
  exact hlt

/-- Zp 取负结果的 val 界。 -/
theorem Zp_neg_reduced (c : Zp) (hq : 0 < c.prime.toNat) :
    (-c).val < c.prime ∧ (-c).prime = c.prime := by
  refine ⟨?_, rfl⟩
  show ((c.prime - c.val) % c.prime) < c.prime
  rw [UInt64.lt_iff_toNat_lt, UInt64.toNat_mod]
  exact Nat.mod_lt _ hq

/-- filterMap 若每个输出都约化，则结果约化（无需输入约化）。 -/
theorem reducedListB_filterMap (q : UInt64) (G : (UMonomial × Zp) → Option (UMonomial × Zp))
    (hG : ∀ x y, G x = some y → y.snd.val < q ∧ y.snd.prime = q) :
    ∀ (l : List (UMonomial × Zp)), reducedListB q (l.filterMap G) = true := by
  intro l
  induction l with
  | nil => rfl
  | cons a rest ih =>
    cases hGa : G a with
    | none => rw [List.filterMap_cons_none hGa]; exact ih
    | some b => rw [List.filterMap_cons_some hGa, reducedListB_cons]; exact ⟨hG a b hGa, ih⟩

/-- negImpl（映射取负）保持约化。 -/
theorem reducedListB_neg (q : UInt64) (hq : 0 < q.toNat) :
    ∀ (l : List (UMonomial × Zp)), reducedListB q l = true →
    reducedListB q (l.map (fun x => (x.1, -x.2))) = true := by
  intro l
  induction l with
  | nil => intro _; rfl
  | cons a rest ih =>
    intro h
    have h' := (reducedListB_cons q a rest).mp h
    rw [List.map_cons, reducedListB_cons]
    have hqp : 0 < a.snd.prime.toNat := by rw [h'.1.2]; exact hq
    have hneg := Zp_neg_reduced a.snd hqp
    exact ⟨⟨h'.1.2 ▸ hneg.1, hneg.2.trans h'.1.2⟩, ih h'.2⟩

/-- addImpl 保持约化。 -/
theorem ReducedB_addImpl (f g : SparsePolyZp) (q : UInt64) (hq : 0 < q.toNat)
    (hf : ReducedB f q) (hg : ReducedB g q) : ReducedB (addImpl f g) q := by
  show reducedListB q (addImpl f g).toList = true
  unfold addImpl
  exact reducedListB_mergeAdd q hq _ _ hf hg

/-- subImpl 保持约化。 -/
theorem ReducedB_subImpl (f g : SparsePolyZp) (q : UInt64) (hq : 0 < q.toNat)
    (hf : ReducedB f q) (hg : ReducedB g q) : ReducedB (subImpl f g) q := by
  have hng : reducedListB q (negImpl g).toList = true := by
    unfold negImpl
    rw [Array.toList_map]
    exact reducedListB_neg q hq g.toList hg
  show reducedListB q (subImpl f g).toList = true
  unfold subImpl addImpl
  exact reducedListB_mergeAdd q hq _ _ hf hng

/-- scaleByMonomial 结果约化（系数素数 = q）。 -/
theorem ReducedB_scaleByMonomial (m : UMonomial) (c : Zp) (f : SparsePolyZp) (q : UInt64)
    (hc : c.prime = q) (hq : 0 < q.toNat) : ReducedB (scaleByMonomial m c f) q := by
  unfold scaleByMonomial
  split
  · exact rfl
  · show reducedListB q (Array.filterMap _ f).toList = true
    rw [Array.toList_filterMap]
    apply reducedListB_filterMap q
    intro x y hxy
    obtain ⟨mf, cf⟩ := x
    simp only at hxy
    split at hxy
    · exact absurd hxy (by simp)
    · injection hxy with hxy'
      rw [← hxy']
      have hqp : 0 < c.prime.toNat := by rw [hc]; exact hq
      have hmul := Zp_mul_reduced c cf hqp
      exact ⟨hc ▸ hmul.1, hmul.2.trans hc⟩

/-- 单项式乘多项式结果约化。 -/
theorem ReducedB_single_mul (m : UMonomial) (c : Zp) (g : SparsePolyZp) (q : UInt64)
    (hc : c.prime = q) (hq : 0 < q.toNat) : ReducedB ((#[(m, c)] : SparsePolyZp) * g) q := by
  have h_eq : (#[(m, c)] : SparsePolyZp) * g = addImpl #[] (scaleByMonomial m c g) := rfl
  rw [h_eq]
  exact ReducedB_addImpl _ _ q hq rfl (ReducedB_scaleByMonomial m c g q hc hq)

-- ── 非零系数不变量：所有系数 val ≠ 0（首项抵消所需——保证 term*g 首项存活） ──

def nonzeroListB : List (UMonomial × Zp) → Bool
  | [] => true
  | a :: rest => (a.snd.val != 0) && nonzeroListB rest

def NonZeroB (f : SparsePolyZp) : Prop := nonzeroListB f.toList = true

instance (f : SparsePolyZp) : Decidable (NonZeroB f) := decEq (nonzeroListB f.toList) true

theorem nonzeroListB_cons (a : UMonomial × Zp) (rest : List (UMonomial × Zp)) :
    nonzeroListB (a :: rest) = true ↔ a.snd.val ≠ 0 ∧ nonzeroListB rest = true := by
  simp only [nonzeroListB, Bool.and_eq_true, bne_iff_ne, ne_eq]

/-- mergeAdd 保持非零系数（case6 的 s.val≠0 条件 + 保留原元素）。 -/
theorem nonzeroListB_mergeAdd : ∀ (xs ys : List (UMonomial × Zp)),
    nonzeroListB xs = true → nonzeroListB ys = true → nonzeroListB (mergeAdd xs ys) = true := by
  intro xs ys
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys => intro _ hy; rw [mergeAdd]; exact hy
  | case2 f fs => intro hx _; rw [mergeAdd]; exact hx
  | case3 f fs g gs hfg ih =>
    intro hx hy
    rw [mergeAdd, if_pos hfg, nonzeroListB_cons]
    exact ⟨((nonzeroListB_cons f fs).mp hx).1, ih ((nonzeroListB_cons f fs).mp hx).2 hy⟩
  | case4 f fs g gs hfg hfg2 ih =>
    intro hx hy
    rw [mergeAdd, if_neg hfg, if_pos hfg2, nonzeroListB_cons]
    exact ⟨((nonzeroListB_cons g gs).mp hy).1, ih hx ((nonzeroListB_cons g gs).mp hy).2⟩
  | case5 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]; exact if_pos hs
    rw [h_eq]
    exact ih ((nonzeroListB_cons f fs).mp hx).2 ((nonzeroListB_cons g gs).mp hy).2
  | case6 f fs g gs hfg hfg2 s hs ih =>
    intro hx hy
    have h_eq : mergeAdd (f :: fs) (g :: gs) = (f.fst, f.snd + g.snd) :: mergeAdd fs gs := by
      rw [mergeAdd, if_neg hfg, if_neg hfg2]; exact if_neg hs
    rw [h_eq, nonzeroListB_cons]
    exact ⟨hs, ih ((nonzeroListB_cons f fs).mp hx).2 ((nonzeroListB_cons g gs).mp hy).2⟩

/-- Zp 取负保持非零（需约化 val<prime）。 -/
theorem Zp_neg_nonzero (c : Zp) (hred : c.val < c.prime) (hnz : c.val ≠ 0)
    (hq : 0 < c.prime.toNat) : (-c).val ≠ 0 := by
  show ((c.prime - c.val) % c.prime) ≠ 0
  have h1 : c.val.toNat < c.prime.toNat := (UInt64.lt_iff_toNat_lt).mp hred
  have h2 : c.val.toNat ≠ 0 := fun h => hnz (UInt64.toNat_inj.mp (by simp [h]))
  have hsub : (c.prime - c.val).toNat = c.prime.toNat - c.val.toNat :=
    UInt64.toNat_sub_of_le _ _ (Nat.le_of_lt h1)
  intro hzero
  have hz : ((c.prime - c.val) % c.prime).toNat = 0 := by simp [hzero]
  rw [UInt64.toNat_mod, hsub, Nat.mod_eq_of_lt (by omega)] at hz
  omega

/-- filterMap 若每个输出都非零，则结果非零。 -/
theorem nonzeroListB_filterMap (G : (UMonomial × Zp) → Option (UMonomial × Zp))
    (hG : ∀ x y, G x = some y → y.snd.val ≠ 0) :
    ∀ (l : List (UMonomial × Zp)), nonzeroListB (l.filterMap G) = true := by
  intro l
  induction l with
  | nil => rfl
  | cons a rest ih =>
    cases hGa : G a with
    | none => rw [List.filterMap_cons_none hGa]; exact ih
    | some b => rw [List.filterMap_cons_some hGa, nonzeroListB_cons]; exact ⟨hG a b hGa, ih⟩

/-- negImpl 保持非零（需约化）。 -/
theorem nonzeroListB_neg (q : UInt64) (hq : 0 < q.toNat) :
    ∀ (l : List (UMonomial × Zp)), reducedListB q l = true → nonzeroListB l = true →
    nonzeroListB (l.map (fun x => (x.1, -x.2))) = true := by
  intro l
  induction l with
  | nil => intro _ _; rfl
  | cons a rest ih =>
    intro hred hnz
    have hred' := (reducedListB_cons q a rest).mp hred
    have hnz' := (nonzeroListB_cons a rest).mp hnz
    rw [List.map_cons, nonzeroListB_cons]
    have hqp : 0 < a.snd.prime.toNat := by rw [hred'.1.2]; exact hq
    exact ⟨Zp_neg_nonzero a.snd (hred'.1.2 ▸ hred'.1.1) hnz'.1 hqp, ih hred'.2 hnz'.2⟩

/-- addImpl 保持非零。 -/
theorem NonZeroB_addImpl (f g : SparsePolyZp) (hf : NonZeroB f) (hg : NonZeroB g) :
    NonZeroB (addImpl f g) := by
  show nonzeroListB (addImpl f g).toList = true
  unfold addImpl
  exact nonzeroListB_mergeAdd _ _ hf hg

/-- subImpl 保持非零（neg 需约化）。 -/
theorem NonZeroB_subImpl (f g : SparsePolyZp) (q : UInt64) (hq : 0 < q.toNat)
    (hredg : ReducedB g q) (hf : NonZeroB f) (hg : NonZeroB g) : NonZeroB (subImpl f g) := by
  have hng : nonzeroListB (negImpl g).toList = true := by
    unfold negImpl
    rw [Array.toList_map]
    exact nonzeroListB_neg q hq g.toList hredg hg
  show nonzeroListB (subImpl f g).toList = true
  unfold subImpl addImpl
  exact nonzeroListB_mergeAdd _ _ hf hng

/-- scaleByMonomial 结果非零（filterMap 丢弃零系数）。 -/
theorem NonZeroB_scaleByMonomial (m : UMonomial) (c : Zp) (f : SparsePolyZp) :
    NonZeroB (scaleByMonomial m c f) := by
  unfold scaleByMonomial
  split
  · exact rfl
  · show nonzeroListB (Array.filterMap _ f).toList = true
    rw [Array.toList_filterMap]
    apply nonzeroListB_filterMap
    intro x y hxy
    obtain ⟨mf, cf⟩ := x
    simp only at hxy
    split at hxy
    · exact absurd hxy (by simp)
    · rename_i hcond
      injection hxy with hxy'
      rw [← hxy']
      exact hcond

/-- 单项式乘多项式结果非零。 -/
theorem NonZeroB_single_mul (m : UMonomial) (c : Zp) (g : SparsePolyZp) :
    NonZeroB ((#[(m, c)] : SparsePolyZp) * g) := by
  have h_eq : (#[(m, c)] : SparsePolyZp) * g = addImpl #[] (scaleByMonomial m c g) := rfl
  rw [h_eq]
  exact NonZeroB_addImpl _ _ rfl (NonZeroB_scaleByMonomial m c g)

-- ── divmodAux 终止性基础设施：首项真抵消 ⇒ 度数严格下降 ──
-- nl-proof: docs/design/divmodAux-termination-nlproof.md

/-- Zp 乘法的 val 用 Nat 表示（供模运算链）。 -/
theorem Zp_mul_val_toNat (x y : Zp) (hq : 0 < x.prime.toNat) :
    (x * y).val.toNat = (x.val.toNat * y.val.toNat) % x.prime.toNat := by
  show ((x.val.toNat * y.val.toNat) % x.prime.toNat).toUInt64.toNat = _
  have hlt : (x.val.toNat * y.val.toNat) % x.prime.toNat < x.prime.toNat := Nat.mod_lt _ hq
  have h2 : (x.val.toNat * y.val.toNat) % x.prime.toNat < UInt64.size :=
    Nat.lt_trans hlt x.prime.toNat_lt_size
  simp [h2]

/-- 引理 B：三重积逆化。(lc*gh).val = 1 ⇒ ((a*lc)*gh).val = a.val（同素数 pm，a 约化）。 -/
theorem Zp_mul3_lc (a lc gh : Zp) (pm : UInt64)
    (hap : a.prime = pm) (hlp : lc.prime = pm) (hgp : gh.prime = pm)
    (hav : a.val < pm) (hlcgh : (lc * gh).val = 1) (hq : 0 < pm.toNat) :
    ((a * lc) * gh).val = a.val := by
  have hq_a : 0 < a.prime.toNat := hap ▸ hq
  have hq_lc : 0 < lc.prime.toNat := hlp ▸ hq
  have hav' : a.val.toNat < pm.toNat := (UInt64.lt_iff_toNat_lt).mp hav
  have hlcgh_nat : (lc.val.toNat * gh.val.toNat) % pm.toNat = 1 := by
    have hm := Zp_mul_val_toNat lc gh hq_lc
    rw [hlp] at hm
    rw [← hm, hlcgh]; rfl
  apply UInt64.toNat_inj.mp
  have halc_prime : (a * lc).prime = pm := hap
  have hq_alc : 0 < (a * lc).prime.toNat := halc_prime ▸ hq
  rw [Zp_mul_val_toNat (a * lc) gh hq_alc, halc_prime, Zp_mul_val_toNat a lc hq_a, hap,
      Nat.mod_mul_mod, Nat.mul_assoc, Nat.mul_mod a.val.toNat (lc.val.toNat * gh.val.toNat),
      hlcgh_nat, Nat.mul_one, Nat.mod_mod, Nat.mod_eq_of_lt hav']

/-- 引理 C：a + (-c) = 0（当 c 与 a 同 val 同 prime，a 约化）。 -/
theorem Zp_add_neg_cancel (a c : Zp) (hval : c.val = a.val) (hprime : c.prime = a.prime)
    (hred : a.val < a.prime) (hq : 0 < a.prime.toNat) : (a + (-c)).val = 0 := by
  show ((a.val.toNat + (-c).val.toNat) % a.prime.toNat).toUInt64 = 0
  have h1 : a.val.toNat < a.prime.toNat := (UInt64.lt_iff_toNat_lt).mp hred
  have hnegval : (-c).val.toNat = (a.prime.toNat - a.val.toNat) % a.prime.toNat := by
    show ((c.prime - c.val) % c.prime).toNat = _
    rw [UInt64.toNat_mod, hprime, hval, UInt64.toNat_sub_of_le _ _ (Nat.le_of_lt h1)]
  rw [hnegval]
  have hkey : (a.val.toNat + (a.prime.toNat - a.val.toNat) % a.prime.toNat) % a.prime.toNat = 0 := by
    rcases Nat.eq_zero_or_pos a.val.toNat with h0 | hpos
    · simp [h0, Nat.mod_self]
    · rw [Nat.mod_eq_of_lt (show a.prime.toNat - a.val.toNat < a.prime.toNat by omega)]
      rw [show a.val.toNat + (a.prime.toNat - a.val.toNat) = a.prime.toNat by omega, Nat.mod_self]
  rw [hkey]; rfl

/-- 引理 D：首项等度且系数和为 0 ⇒ mergeAdd 结果全元素度数 < d。 -/
theorem mergeAdd_cancel_lead (d : Nat) (a b : UMonomial × Zp)
    (ra rb : List (UMonomial × Zp))
    (ha : a.fst.deg = d) (hb : b.fst.deg = d)
    (hsum : (a.snd + b.snd).val = 0)
    (hra : ∀ x ∈ ra, x.fst.deg < d) (hrb : ∀ y ∈ rb, y.fst.deg < d) :
    ∀ z ∈ mergeAdd (a :: ra) (b :: rb), z.fst.deg < d := by
  have h_eq : mergeAdd (a :: ra) (b :: rb) = mergeAdd ra rb := by
    rw [mergeAdd, if_neg (show ¬ a.fst.deg > b.fst.deg by omega),
        if_neg (show ¬ a.fst.deg < b.fst.deg by omega)]
    exact if_pos hsum
  rw [h_eq]
  exact mergeAdd_lt_all d ra rb hra hrb

/-- 非空数组的 toList = 首项 :: 尾。 -/
theorem toList_cons_of_ne_empty (r : SparsePolyZp) (h : ¬ r.isEmpty) :
    r.toList = (r[0]!) :: r.toList.tail := by
  match hrl : r.toList with
  | [] =>
    exfalso; apply h
    simp [Array.isEmpty, Array.size_eq_length_toList, hrl]
  | hd :: tl =>
    have hsz : 0 < r.size := by rw [Array.size_eq_length_toList, hrl]; simp
    have hhead : r[0]! = hd := by
      rw [getElem!_pos r 0 hsz]
      have h2 : r[0] = r.toList[0]'(by rw [Array.size_eq_length_toList] at hsz; exact hsz) :=
        (Array.getElem_toList (xs := r) (i := 0) hsz).symm
      rw [h2]; simp [hrl]
    simp only [hhead, hrl, List.tail_cons]

/-- addImpl #[] x 的 toList = x 的 toList（mergeAdd [] = id）。 -/
theorem addImpl_nil_toList (x : SparsePolyZp) : (addImpl #[] x).toList = x.toList := by
  unfold addImpl; simp [mergeAdd]

/-- scaleByMonomial 内部 filterMap 的 lambda（提取为 def 以便 unfold + 复用）。 -/
def scaleLambda (m : UMonomial) (coeff : Zp) : (UMonomial × Zp) → Option (UMonomial × Zp) :=
  fun x => let prod := coeff * x.2
    if prod.val = 0 then none else some (UMonomial.mk (m.deg + x.1.deg), prod)

/-- 引理 A：scaleByMonomial 的首项 = (m.deg + dg, coeff·gh)，尾全部度数 < m.deg + dg。 -/
theorem scaleByMonomial_head_drop (m : UMonomial) (coeff : Zp) (g : SparsePolyZp)
    (gh : UMonomial × Zp) (gt : List (UMonomial × Zp)) (dgv : Nat)
    (hgl : g.toList = gh :: gt) (hcoeff_nz : coeff.val ≠ 0)
    (hhead_nz : (coeff * gh.snd).val ≠ 0)
    (hghdeg : gh.fst.deg = dgv) (hgt_lt : ∀ x ∈ gt, x.fst.deg < dgv) :
    ∃ T, (scaleByMonomial m coeff g).toList
          = (UMonomial.mk (m.deg + dgv), coeff * gh.snd) :: T
      ∧ (∀ z ∈ T, z.fst.deg < m.deg + dgv) := by
  have hF_shift : ∀ x y, scaleLambda m coeff x = some y → y.fst.deg = m.deg + x.fst.deg := by
    intro x y hxy
    unfold scaleLambda at hxy; simp only at hxy; split at hxy
    · exact absurd hxy (by simp)
    · injection hxy with h'; rw [← h']
  have hF_gh : scaleLambda m coeff gh = some (UMonomial.mk (m.deg + dgv), coeff * gh.snd) := by
    unfold scaleLambda; simp only; rw [if_neg hhead_nz, hghdeg]
  refine ⟨gt.filterMap (scaleLambda m coeff), ?_, ?_⟩
  · show (scaleByMonomial m coeff g).toList = _
    unfold scaleByMonomial
    rw [if_neg hcoeff_nz, Array.toList_filterMap]
    show g.toList.filterMap (scaleLambda m coeff) = _
    rw [hgl, List.filterMap_cons_some hF_gh]
  · intro z hz
    exact (sortedListB_filterMapShift m.deg (scaleLambda m coeff) hF_shift gt).1 dgv hgt_lt z hz

/-- subImpl 的 toList = mergeAdd（addImpl 定义展开）。 -/
theorem subImpl_toList (f g : SparsePolyZp) :
    (subImpl f g).toList = mergeAdd f.toList (negImpl g).toList := by
  unfold subImpl addImpl; simp

/-- negImpl 的 toList = 逐项取负 map。 -/
theorem negImpl_toList (f : SparsePolyZp) :
    (negImpl f).toList = f.toList.map (fun x => (x.1, -x.2)) := by
  unfold negImpl; rw [Array.toList_map]

/-- 核心：长除法一步后余式 r' = r - term*g 的所有项度数严格 < 原首项度数 dr。
    首项真抵消（coeff·lc(g) = r 首项系数，相减为 0）⇒ 度量严格下降。 -/
theorem divmod_step_drop (g r : SparsePolyZp) (dg : Nat) (lc_g_inv : Zp) (pm : UInt64)
    (hq : 0 < pm.toNat)
    (hg_ne : ¬ g.isEmpty) (hg_sorted : Sorted g) (hg_red : ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg)
    (hlp : lc_g_inv.prime = pm) (h_lc : (lc_g_inv * (g[0]!).snd).val = 1)
    (hr_ne : ¬ r.isEmpty) (hr_sorted : Sorted r) (hr_red : ReducedB r pm) (hr_nz : NonZeroB r)
    (hdg_le : dg ≤ (r[0]!).fst.deg) :
    ∀ z ∈ (r - (#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g).toList,
      z.fst.deg < (r[0]!).fst.deg := by
  -- r 分解
  have hr_red' : reducedListB pm r.toList = true := hr_red
  have hr_nz' : nonzeroListB r.toList = true := hr_nz
  have hr_sorted' : sortedListB r.toList = true := hr_sorted
  have hg_red' : reducedListB pm g.toList = true := hg_red
  have hg_sorted' : sortedListB g.toList = true := hg_sorted
  have hrl : r.toList = (r[0]!) :: r.toList.tail := toList_cons_of_ne_empty r hr_ne
  have ha_val_lt : (r[0]!).snd.val < pm := ((reducedListB_cons pm _ _).mp (hrl ▸ hr_red')).1.1
  have ha_prime : (r[0]!).snd.prime = pm := ((reducedListB_cons pm _ _).mp (hrl ▸ hr_red')).1.2
  have ha_nz : (r[0]!).snd.val ≠ 0 := ((nonzeroListB_cons _ _).mp (hrl ▸ hr_nz')).1
  have hra_lt : ∀ x ∈ r.toList.tail, x.fst.deg < (r[0]!).fst.deg :=
    ((sortedListB_iff _ _).mp (hrl ▸ hr_sorted')).1
  -- g 分解
  have hgl : g.toList = (g[0]!) :: g.toList.tail := toList_cons_of_ne_empty g hg_ne
  have hgh_prime : (g[0]!).snd.prime = pm := ((reducedListB_cons pm _ _).mp (hgl ▸ hg_red')).1.2
  have hgt_lt : ∀ x ∈ g.toList.tail, x.fst.deg < dg := by
    intro x hx
    have hh := ((sortedListB_iff _ _).mp (hgl ▸ hg_sorted')).1 x hx
    rwa [h_dg] at hh
  -- coeff/c 性质
  have hcoeff_prime : ((r[0]!).snd * lc_g_inv).prime = pm := ha_prime
  have hc_val : (((r[0]!).snd * lc_g_inv) * (g[0]!).snd).val = (r[0]!).snd.val :=
    Zp_mul3_lc (r[0]!).snd lc_g_inv (g[0]!).snd pm ha_prime hlp hgh_prime ha_val_lt h_lc hq
  have hc_prime : (((r[0]!).snd * lc_g_inv) * (g[0]!).snd).prime = pm := hcoeff_prime
  have hc_nz : (((r[0]!).snd * lc_g_inv) * (g[0]!).snd).val ≠ 0 := by rw [hc_val]; exact ha_nz
  have hcoeff_nz : ((r[0]!).snd * lc_g_inv).val ≠ 0 := by
    intro h0
    have hz0 : (((r[0]!).snd * lc_g_inv) * (g[0]!).snd).val.toNat = 0 := by
      rw [Zp_mul_val_toNat _ _ (by rw [hcoeff_prime]; exact hq), h0]; simp
    exact hc_nz (UInt64.toNat_inj.mp (by rw [hz0]; rfl))
  -- dr - dg + dg = dr
  have hdrdg : (r[0]!).fst.deg - dg + dg = (r[0]!).fst.deg := by omega
  -- term*g 分解
  obtain ⟨T, hT_eq, hT_lt⟩ := scaleByMonomial_head_drop ⟨(r[0]!).fst.deg - dg⟩
    ((r[0]!).snd * lc_g_inv) g (g[0]!) (g.toList.tail) dg hgl hcoeff_nz hc_nz h_dg hgt_lt
  -- 转 (term*g).toList
  have hterm_eq : ((#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g).toList
      = (scaleByMonomial ⟨(r[0]!).fst.deg - dg⟩ ((r[0]!).snd * lc_g_inv) g).toList := by
    rw [show ((#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g)
          = addImpl #[] (scaleByMonomial ⟨(r[0]!).fst.deg - dg⟩ ((r[0]!).snd * lc_g_inv) g) from rfl,
        addImpl_nil_toList]
  -- 头 monomial deg 化简
  have heqdeg : (⟨(r[0]!).fst.deg - dg⟩ : UMonomial).deg + dg = (r[0]!).fst.deg := by
    show (r[0]!).fst.deg - dg + dg = (r[0]!).fst.deg; omega
  have hheadmono : (UMonomial.mk ((⟨(r[0]!).fst.deg - dg⟩ : UMonomial).deg + dg))
      = (⟨(r[0]!).fst.deg⟩ : UMonomial) := by rw [heqdeg]
  -- (term*g).toList = (⟨dr⟩, c) :: T
  have hTG : ((#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g).toList
      = (⟨(r[0]!).fst.deg⟩, ((r[0]!).snd * lc_g_inv) * (g[0]!).snd) :: T := by
    rw [hterm_eq, hT_eq, hheadmono]
  have hT_lt' : ∀ z ∈ T, z.fst.deg < (r[0]!).fst.deg := by
    intro z hz; have hh := hT_lt z hz; rwa [heqdeg] at hh
  -- 目标转 mergeAdd 形式
  rw [show r - (#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g
        = subImpl r ((#[(⟨(r[0]!).fst.deg - dg⟩, (r[0]!).snd * lc_g_inv)] : SparsePolyZp) * g) from rfl,
      subImpl_toList, negImpl_toList, hTG, hrl]
  rw [List.map_cons]
  apply mergeAdd_cancel_lead (r[0]!).fst.deg (r[0]!)
    (⟨(r[0]!).fst.deg⟩, -((r[0]!).snd * lc_g_inv * (g[0]!).snd))
  · rfl
  · rfl
  · exact Zp_add_neg_cancel (r[0]!).snd (((r[0]!).snd * lc_g_inv) * (g[0]!).snd)
      hc_val (hc_prime.trans ha_prime.symm) (ha_prime ▸ ha_val_lt) (by rw [ha_prime]; exact hq)
  · exact hra_lt
  · intro y hy
    rw [List.mem_map] at hy
    obtain ⟨x, hxT, hxy⟩ := hy
    rw [← hxy]
    exact hT_lt' x hxT

-- 长除法主循环：从 r 中持续减去 g 的倍数，累加商到 q。
-- 良基递归（真 WF，0 admit）：measure = if r.isEmpty then 0 else r[0].deg+1。
-- 终止依赖「首项真抵消 ⇒ deg 严格下降」（divmod_step_drop），需以下不变量进签名：
--   pm/hq/hlp/h_lc : 系数同素数 pm、lc_g_inv 为 g 首项之逆
--   hg_ne/hg_red/h_dg/h_sorted_g : 除数 g 非空、约化、首项度 dg、降序
--   h_sorted_r/hr_red/hr_nz : 余式 r 降序、约化、非零系数（r' 递归保持）
-- 全部在 `divmod` 包装处经可判定的依赖 if 提供。
def divmodAux (g : SparsePolyZp) (dg : Nat) (lc_g_inv : Zp) (pm : UInt64)
    (hq : 0 < pm.toNat) (hg_ne : ¬ g.isEmpty) (hg_red : ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg) (hlp : lc_g_inv.prime = pm)
    (h_lc : (lc_g_inv * (g[0]!).snd).val = 1) (h_sorted_g : Sorted g)
    (q r : SparsePolyZp) (h_sorted_r : Sorted r) (hr_red : ReducedB r pm)
    (hr_nz : NonZeroB r) : SparsePolyZp × SparsePolyZp :=
  if hr : r.isEmpty then (q, r)
  else
    let dr := r[0]!.fst.deg
    if hd : dr < dg then (q, r)
    else
      let coeff := r[0]!.snd * lc_g_inv
      let d := dr - dg
      let term : SparsePolyZp := #[(⟨d⟩, coeff)]
      let r' := r - (term * g)
      let q' := q.push (⟨d⟩, coeff)
      have h_sorted_r' : Sorted r' :=
        subImpl_sorted r (term * g) h_sorted_r (single_mul_sorted ⟨d⟩ coeff g h_sorted_g)
      have hcoeff_prime : coeff.prime = pm := by
        have hr_red' : reducedListB pm r.toList = true := hr_red
        have hrl : r.toList = r[0]! :: r.toList.tail := toList_cons_of_ne_empty r hr
        exact ((reducedListB_cons pm _ _).mp (hrl ▸ hr_red')).1.2
      have h_tg_red : ReducedB (term * g) pm := ReducedB_single_mul ⟨d⟩ coeff g pm hcoeff_prime hq
      have h_tg_nz : NonZeroB (term * g) := NonZeroB_single_mul ⟨d⟩ coeff g
      have hr'_red : ReducedB r' pm := ReducedB_subImpl r (term * g) pm hq hr_red h_tg_red
      have hr'_nz : NonZeroB r' := NonZeroB_subImpl r (term * g) pm hq h_tg_red hr_nz h_tg_nz
      divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc h_sorted_g q' r'
        h_sorted_r' hr'_red hr'_nz
-- measure 用 `if r.isEmpty then 0 else r[0].deg+1`：空余式度量为 0，非空为首项度数+1，
-- 保证「r' 空 或 r'首项度<r首项度」两种情形都严格递减（含 constant÷constant 边界）。
termination_by (if r.isEmpty then 0 else r[0]!.fst.deg + 1)
decreasing_by
  have hr_ne : ¬ r.isEmpty := hr
  have hdg_le : dg ≤ r[0]!.fst.deg := Nat.le_of_not_lt hd
  have hdrop := divmod_step_drop g r dg lc_g_inv pm hq hg_ne h_sorted_g hg_red h_dg hlp h_lc
    hr_ne h_sorted_r hr_red hr_nz hdg_le
  rw [if_neg hr_ne]
  generalize hRdef : (r - (#[(⟨r[0]!.fst.deg - dg⟩, r[0]!.snd * lc_g_inv)] : SparsePolyZp) * g) = R'
    at hdrop ⊢
  split
  · omega
  · rename_i hR'ne
    have hr'l := toList_cons_of_ne_empty R' hR'ne
    have hmem : R'[0]! ∈ R'.toList := by rw [hr'l]; simp
    have hlt := hdrop R'[0]! hmem
    omega

-- divmodAux 输出余式：空 或 首项度 < dg（欧几里得余式界，供 gcdAux WF 化）。
theorem divmodAux_snd_deg_lt (g : SparsePolyZp) (dg : Nat) (lc_g_inv : Zp) (pm : UInt64)
    (hq : 0 < pm.toNat) (hg_ne : ¬ g.isEmpty) (hg_red : ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg) (hlp : lc_g_inv.prime = pm)
    (h_lc : (lc_g_inv * (g[0]!).snd).val = 1) (h_sorted_g : Sorted g) :
    ∀ (q r : SparsePolyZp) (h_sorted_r : Sorted r) (hr_red : ReducedB r pm) (hr_nz : NonZeroB r),
      (divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc h_sorted_g
        q r h_sorted_r hr_red hr_nz).snd.isEmpty = true
      ∨ (divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc h_sorted_g
        q r h_sorted_r hr_red hr_nz).snd[0]!.fst.deg < dg := by
  suffices H : ∀ n, ∀ (q r : SparsePolyZp) (h_sorted_r : Sorted r) (hr_red : ReducedB r pm)
      (hr_nz : NonZeroB r), (if r.isEmpty then 0 else r[0]!.fst.deg + 1) = n →
      (divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc h_sorted_g
        q r h_sorted_r hr_red hr_nz).snd.isEmpty = true
      ∨ (divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc h_sorted_g
        q r h_sorted_r hr_red hr_nz).snd[0]!.fst.deg < dg by
    intro q r h_sorted_r hr_red hr_nz; exact H _ q r h_sorted_r hr_red hr_nz rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro q r h_sorted_r hr_red hr_nz hn
    rw [divmodAux]
    split
    · rename_i hr; left; exact hr
    · rename_i hr
      dsimp only
      split
      · rename_i hd; right; exact hd
      · rename_i hd
        have hr_ne : ¬ r.isEmpty := hr
        have hdg_le : dg ≤ r[0]!.fst.deg := Nat.le_of_not_lt hd
        have hdrop := divmod_step_drop g r dg lc_g_inv pm hq hg_ne h_sorted_g hg_red h_dg hlp h_lc
          hr_ne h_sorted_r hr_red hr_nz hdg_le
        have hmeas : (if (r - (#[(⟨r[0]!.fst.deg - dg⟩, r[0]!.snd * lc_g_inv)] : SparsePolyZp) * g).isEmpty
            then 0
            else (r - (#[(⟨r[0]!.fst.deg - dg⟩, r[0]!.snd * lc_g_inv)] : SparsePolyZp) * g)[0]!.fst.deg + 1)
            < n := by
          rw [← hn, if_neg hr_ne]
          generalize hRdef : (r - (#[(⟨r[0]!.fst.deg - dg⟩, r[0]!.snd * lc_g_inv)] : SparsePolyZp) * g) = R'
            at hdrop ⊢
          split
          · omega
          · rename_i hR'ne
            have hr'l := toList_cons_of_ne_empty R' hR'ne
            have hmem : R'[0]! ∈ R'.toList := by rw [hr'l]; simp
            have hlt := hdrop R'[0]! hmem
            omega
        exact ih _ hmeas _ _ _ _ _ rfl

-- 多项式长除法：f = q * g + r, deg(r) < deg(g)
-- 退化 / 前提不满足（g 空、lc 非逆、未降序、未约化、有零系数）→ 返回 (#[], f) 占位。
-- 合法输入（域上、降序、g≠0、同素数约化）下依赖 if 走真分支 = WF 主循环。
def divmod (f g : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  if h : ¬ g.isEmpty ∧ ((g[0]!.snd.inv * (g[0]!).snd).val = 1) ∧ Sorted g ∧ Sorted f
      ∧ 0 < (g[0]!.snd.prime).toNat ∧ ReducedB g (g[0]!.snd.prime)
      ∧ ReducedB f (g[0]!.snd.prime) ∧ NonZeroB f then
    divmodAux g (g[0]!.fst.deg) (g[0]!.snd.inv) (g[0]!.snd.prime)
      h.2.2.2.2.1 h.1 h.2.2.2.2.2.1 rfl rfl h.2.1 h.2.2.1 #[] f
      h.2.2.2.1 h.2.2.2.2.2.2.1 h.2.2.2.2.2.2.2
  else (#[], f)

-- 标量乘（用于首一化）：multiply all coefs by const Zp
def scalarMul (c : Zp) (f : SparsePolyZp) : SparsePolyZp :=
  f.filterMap (fun (m, x) =>
    let new_val := x * c
    if new_val.val = 0 then none
    else some (m, new_val))

-- 首一化：multiply by inv(lc) — 与 C++ polynomial_GCD 输出约定一致
def makeMonic (f : SparsePolyZp) : SparsePolyZp :=
  if f.isEmpty then f
  else
    let lc_inv := f[0]!.snd.inv
    scalarMul lc_inv f

-- 欧几里得 GCD：gcd(f, g) = gcd(g, f mod g)；最后首一化（与 C++ 一致）。
--
-- 合法的规范形输入上，`divmodAux_snd_deg_lt` 保证余式为空或首项次数严格下降。
-- `divmod` 对不满足其表示不变量的输入会走保守回退分支；这里也显式检查下降条件，
-- 使实现对所有 Array 输入总化，同时在合法输入上保持原 Euclid 控制流不变。
def gcdAux (f g : SparsePolyZp) : SparsePolyZp :=
  if hg : g.isEmpty then f
  else
    let r := (divmod f g).snd
    if hr : r.isEmpty then g
    else if hd : r[0]!.fst.deg < g[0]!.fst.deg then gcdAux g r
    else g
termination_by if g.isEmpty then 0 else g[0]!.fst.deg + 1
decreasing_by
  simp_wf
  have hg' : g ≠ #[] := by
    intro h
    subst g
    simp at hg
  have hr' : r ≠ #[] := by
    intro h
    apply hr
    simp [h]
  change (if r = #[] then 0 else r[0]!.fst.deg + 1) <
    (if g = #[] then 0 else g[0]!.fst.deg + 1)
  rw [if_neg hg', if_neg hr']
  omega

def gcd (f g : SparsePolyZp) : SparsePolyZp :=
  makeMonic (gcdAux f g)

-- 扩展欧几里得：返回 (g, s, t) 满足 a*s + b*t = g
-- 类比 Nat.extGcd：
--   if b = 0 then (a, 1, 0)
--   else (q, r) := divmod a b
--        (g, s', t') := extGcd b r
--        a*t' + b*(s' - q*t') = g
-- C++ 约定：最后 g, s, t 三者同乘 inv(lc(g))，使 g 首一
partial def extGcdAux (a b : SparsePolyZp) : SparsePolyZp × SparsePolyZp × SparsePolyZp :=
  if b.isEmpty then
    if a.isEmpty then (#[], #[], #[])
    else
      let p := a[0]!.snd.prime
      let one_poly : SparsePolyZp := #[(⟨0⟩, Zp.ofUInt64 1 p)]
      (a, one_poly, #[])
  else
    let (q, r) := divmod a b
    let (g, s', t') := extGcdAux b r
    -- a*t' + b*(s' - q*t') = g
    let new_t := s' - q * t'
    (g, t', new_t)

def extGcd (a b : SparsePolyZp) : SparsePolyZp × SparsePolyZp × SparsePolyZp :=
  let (g, s, t) := extGcdAux a b
  if g.isEmpty then (g, s, t)
  else
    let lc_inv := g[0]!.snd.inv
    (scalarMul lc_inv g, scalarMul lc_inv s, scalarMul lc_inv t)

end SparsePolyZp

-- SparsePolyZp 特化 instance（高优先级，覆盖兜底 default）
instance : HasPolyGCD SparsePolyZp where
  polyGCD := SparsePolyZp.gcd

instance : HasPolyDivmod SparsePolyZp where
  polyDivmod := SparsePolyZp.divmod

instance : HasPolyGCDEEA SparsePolyZp where
  polyGCDEEA := SparsePolyZp.extGcd

-- #eval 数值验证（小例 over F_5）
-- 注：divmodAux 现为真 WF 递归（decreasing_by 完整证明，0 admit，无 sorryAx），
-- divmod 不再 sorry-tainted，#eval 已恢复。
-- (x^2 - 1) / (x - 1) = x + 1, remainder 0
#eval SparsePolyZp.divmod
  (#[(⟨2⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)  -- x^2 - 1
  (#[(⟨1⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)  -- x - 1
-- 期望: (#[(1, 1), (0, 1)], #[])  — q = x+1, r = 0

-- gcd(x^2 - 1, x - 1) = x - 1（因为 x-1 整除；F_5 下 -1 = 4）
#eval SparsePolyZp.gcd
  (#[(⟨2⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)
  (#[(⟨1⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)

-- #eval 数值验证（手算）
-- f = 2x^2 + 3x，g = x^2 + 4 (over F_7)
-- f + g = 3x^2 + 3x + 4
#eval (#[(⟨2⟩, Zp.ofInt 2 7), (⟨1⟩, Zp.ofInt 3 7)] : SparsePolyZp) +
      (#[(⟨2⟩, Zp.ofInt 1 7), (⟨0⟩, Zp.ofInt 4 7)] : SparsePolyZp)
-- 期望: #[(⟨2⟩, ⟨3, 7⟩), (⟨1⟩, ⟨3, 7⟩), (⟨0⟩, ⟨4, 7⟩)]

-- f - f = 0（删空）
#eval (#[(⟨2⟩, Zp.ofInt 2 7), (⟨1⟩, Zp.ofInt 3 7)] : SparsePolyZp) -
      (#[(⟨2⟩, Zp.ofInt 2 7), (⟨1⟩, Zp.ofInt 3 7)] : SparsePolyZp)
-- 期望: #[]

-- (x + 1)(x - 1) = x^2 - 1 (over F_7)
#eval (#[(⟨1⟩, Zp.ofInt 1 7), (⟨0⟩, Zp.ofInt 1 7)] : SparsePolyZp) *
      (#[(⟨1⟩, Zp.ofInt 1 7), (⟨0⟩, Zp.ofInt (-1) 7)] : SparsePolyZp)
-- 期望: #[(⟨2⟩, ⟨1, 7⟩), (⟨0⟩, ⟨6, 7⟩)] (因为 -1 ≡ 6 mod 7)

-- f * 0 = 0
#eval (#[(⟨2⟩, Zp.ofInt 2 7), (⟨1⟩, Zp.ofInt 3 7)] : SparsePolyZp) *
      (#[] : SparsePolyZp)
-- 期望: #[]
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
-- 阶段 G-F：补各种 exponent 类型（C++ 用 int64_t / uint64_t / int 作 exponent）
def Zp.pow (base : Zp) : Nat → Zp
  | 0 => ⟨1 % base.prime, base.prime⟩
  | k+1 => base * (Zp.pow base k)
instance : HPow Zp Nat Zp where hPow := Zp.pow
instance : HPow Zp Int64 Zp where hPow base e := Zp.pow base e.toNatClampNeg
instance : HPow Zp UInt64 Zp where hPow base e := Zp.pow base e.toNat
instance : HPow Int Int64 Int where hPow base e := base ^ e.toNatClampNeg
instance : HPow ZZ Int64 ZZ where hPow base e := base ^ e.toNatClampNeg

-- MvPolyZZ / MvPolyZp 的算术 stub
-- HShiftLeft Int UInt64：C++ side `1 << bits` 等
instance : HShiftLeft Int UInt64 Int where hShiftLeft a b := a <<< b.toNat
instance : HShiftLeft Int Int Int where hShiftLeft a b := a <<< b.toNat
-- MvPoly 算术：append 后 normalize（merge 同 monomial、排序、剔零）
-- MvPolyZZ
def MvPolyZZ.addImpl (f g : MvPolyZZ) : MvPolyZZ :=
  MvPolyZZ.normalization (f ++ g)

def MvPolyZZ.negImpl (f : MvPolyZZ) : MvPolyZZ :=
  f.map (fun (m, c) => (m, -c))

def MvPolyZZ.subImpl (f g : MvPolyZZ) : MvPolyZZ :=
  MvPolyZZ.addImpl f (MvPolyZZ.negImpl g)

def MvPolyZZ.mulImpl (f g : MvPolyZZ) : MvPolyZZ :=
  let prods : MvPolyZZ := f.foldl (fun acc (mf, cf) =>
    g.foldl (fun a (mg, cg) =>
      a.push (Monomial.mul mf mg, cf * cg)) acc) #[]
  MvPolyZZ.normalization prods

-- MvPolyZp
def MvPolyZp.addImpl (f g : MvPolyZp) : MvPolyZp :=
  MvPolyZp.normalization (f ++ g)

def MvPolyZp.negImpl (f : MvPolyZp) : MvPolyZp :=
  f.map (fun (m, c) => (m, -c))

def MvPolyZp.subImpl (f g : MvPolyZp) : MvPolyZp :=
  MvPolyZp.addImpl f (MvPolyZp.negImpl g)

def MvPolyZp.mulImpl (f g : MvPolyZp) : MvPolyZp :=
  let prods : MvPolyZp := f.foldl (fun acc (mf, cf) =>
    g.foldl (fun a (mg, cg) =>
      a.push (Monomial.mul mf mg, cf * cg)) acc) #[]
  MvPolyZp.normalization prods

instance : HMul MvPolyZZ MvPolyZZ MvPolyZZ where hMul := MvPolyZZ.mulImpl
instance : HAdd MvPolyZZ MvPolyZZ MvPolyZZ where hAdd := MvPolyZZ.addImpl
instance : HSub MvPolyZZ MvPolyZZ MvPolyZZ where hSub := MvPolyZZ.subImpl
instance : Neg MvPolyZZ where neg := MvPolyZZ.negImpl
instance : HMul MvPolyZp MvPolyZp MvPolyZp where hMul := MvPolyZp.mulImpl
instance : HAdd MvPolyZp MvPolyZp MvPolyZp where hAdd := MvPolyZp.addImpl
instance : HSub MvPolyZp MvPolyZp MvPolyZp where hSub := MvPolyZp.subImpl
instance : Neg MvPolyZp where neg := MvPolyZp.negImpl

-- derivative typeclass：C++ free function `derivative(poly)` 对所有 poly 类型多态
class HasDerivative (α : Type) where
  derivative : α → α

-- HasDerivative SparsePolyZp 真实实现已在 §4 namespace SparsePolyZp（line ~108）
-- 但 typeclass instance 还要绑定到那个 def
instance : HasDerivative SparsePolyZp where derivative := SparsePolyZp.derivative

-- HasDerivative SparsePolyZZ 真实实现见 §5c（abbrev SparsePolyZZ 之后）
-- MvPoly 的 1-arg derivative：对应 C++ polynomial.hh:880 derivative<lex_<...>>
-- 关于第一个变量（主变量，lex order 的最高 var）求导
def MvPolyZZ.derivativeMv (f : MvPolyZZ) : MvPolyZZ :=
  if h : 0 < f.size then
    let firstMono := f[0].fst
    if hm : 0 < firstMono.size then
      let mainVar : Variable := firstMono[0].fst
      f.filterMap (fun term =>
        let mono := term.fst
        let c := term.snd
        if hm' : 0 < mono.size then
          let firstVar := mono[0].fst
          let firstExp := mono[0].snd
          if firstVar = mainVar then
            -- 该 term 含主变量：新系数 = 原系数 * exp；新 mono 去掉 mainVar 或降一次
            let newCoef : Int := c * firstExp.toInt
            let rest : Monomial := mono.extract 1 mono.size
            let newMono : Monomial :=
              if firstExp = 1 then rest
              else #[(mainVar, firstExp - 1)] ++ rest
            some (newMono, newCoef)
          else
            -- 该 term 首变量不是主变量（lex 在主变量后），无主变量项 → 求导为 0
            none
        else
          -- 空 mono（常数项）→ 求导为 0
          none)
    else
      #[]  -- 首项是常数（特殊情形）
  else
    #[]  -- 空多项式

instance : HasDerivative MvPolyZZ where derivative := MvPolyZZ.derivativeMv

-- 兜底（MvPolyZp 等仍走 priority := 0 fallback）
instance (priority := 0) {α : Type} : HasDerivative α where derivative f := f

def derivative {α : Type} [HasDerivative α] (a : α) : α := HasDerivative.derivative a

-- squarefreefactorize 现在由 cpp2lean v2 翻译进 Generated/Corpus.lean
-- （`squarefreefactorize_lex_ir`），不在 Model.lean 提供 stub。
-- caller 经 Pass 8 codegen 派发到 _lex_ir 后缀，不会调到此处。

-- poly_convert: 跨域多项式转换 typeclass dispatch
class HasPolyConvert (α β : Type) where
  convert : α → β → β

class HasPolyConvert3 (α β γ : Type) where
  convert3 : α → β → γ → β

-- 通用兜底：返回 target 不变
instance (priority := 0) {α β : Type} : HasPolyConvert α β where
  convert _ target := target

instance (priority := 0) {α β γ : Type} : HasPolyConvert3 α β γ where
  convert3 _ target _ := target

def poly_convert {α β : Type} [HasPolyConvert α β] (f : α) (target : β) : β :=
  HasPolyConvert.convert f target

def poly_convert3 {α β γ : Type} [HasPolyConvert3 α β γ] (f : α) (target : β) (ctx : γ) : β :=
  HasPolyConvert3.convert3 f target ctx

-- 具体实例（特化高优先级覆盖兜底）见 §5c 末尾（abbrev SparsePolyZZ 之后）

-- SparsePolyZZ 的 OfNat 0 实例：见 §5c（abbrev 定义之后）

-- C++ 全局常量 / 宏：占位（B2B 时填实际值）
def ZASSENHAUS_THRESHOLD : Int32 := 10  -- C++: clpoly/polynomial_factorize_univar.hh:889
def __g_use_large_prime : Bool := false

-- ZZ.invert: GMP `mpz_invert(out, num, mod)` → 模逆元，true 表存在。
-- 数学定义：返回 inv 使得 num * inv ≡ 1 (mod mod)，仅当 gcd(num, mod) = 1。
-- 退化：mod ≤ 1 → false（mod 0/1 域无意义）。

-- 扩展欧几里得（Nat 版，well-founded on b）
-- 返回 (g, x, y) 满足 (a : Int) * x + (b : Int) * y = g, g = gcd(a, b)
-- x, y 用 Int 因可能为负
-- spec proof: 见 CLPoly.Math.Bigint (Bezout)
def Nat.extGcd (a b : Nat) : Nat × Int × Int :=
  if h : b = 0 then (a, 1, 0)
  else
    let q := a / b
    let r := a % b
    let (g, x, y) := Nat.extGcd b r
    (g, y, x - (q : Int) * y)
termination_by b
decreasing_by
  exact Nat.mod_lt a (Nat.pos_of_ne_zero h)

-- 内部计算：返回 (success, inv)
-- spec proof: 见 CLPoly.Math.Bigint (invert_correct)
def ZZ.invertImpl (a m : Int) : Bool × Int :=
  if m ≤ 1 then (false, 0)
  else
    -- 把 a 折回 [0, m)
    let am : Int := ((a % m) + m) % m
    if am = 0 then (false, 0)
    else
      -- extGcd am.natAbs m.natAbs：返回 (g, x, y) 满足 am * x + m * y = g
      -- 当 g = 1 时 am * x ≡ 1 (mod m)，即 x 是 am 的逆
      let (g, x, _) := Nat.extGcd am.natAbs m.natAbs
      if g = 1 then
        let inv := ((x % m) + m) % m
        (true, inv)
      else
        (false, 0)

-- 公开签名（与 Pass 5 emit 一致 + Pass 2b ref-elim）：返回 (Bool × ZZ) tuple。
-- class_map 的 OUTPUT_PARAMS 把 _out 当 ref-out → Pass 2b 把调用点 destructure
-- 为 `let __refret := ZZ.invert old_out num mod; lc_inv := __refret.snd`。
def ZZ.invert (_out num mod : ZZ) : Bool × ZZ := ZZ.invertImpl num mod

-- #eval 数值验证
#eval ZZ.invertImpl 3 7      -- (true, 5)：3 * 5 = 15 ≡ 1 mod 7
#eval ZZ.invertImpl 2 4      -- (false, 0)：gcd(2, 4) = 2 ≠ 1
#eval ZZ.invertImpl 5 13     -- (true, 8)：5 * 8 = 40 ≡ 1 mod 13
#eval ZZ.invertImpl 7 13     -- (true, 2)：7 * 2 = 14 ≡ 1 mod 13
#eval ZZ.invertImpl 0 5      -- (false, 0)
#eval ZZ.invertImpl 1 5      -- (true, 1)

-- ZZ.fdiv_q / ZZ.fdiv_r: 向下取整除法（floor division，对应 GMP mpz_fdiv_*）
-- 数学定义：q = ⌊a/b⌋，r = a - b·q（即 b·q + r = a，0 ≤ r < |b| 当 b > 0）
-- Lean stdlib 的 Int.fdiv / Int.fmod 即此语义；C++ 自带 `a / b` (即 Int.tdiv) 是
-- 截断除法，对负被除数与 GMP 不一致，故必须用 fdiv/fmod。
-- _out 参数：Pass 2 ref-elim 的"前一 SSA out 值"占位，未使用。
def ZZ.fdiv_q (_out a b : ZZ) : ZZ := Int.fdiv a b
def ZZ.fdiv_r (_out a b : ZZ) : ZZ := Int.fmod a b

-- spec：除法恒等式 b·q + r = a（Lean stdlib 已证）
theorem ZZ.fdiv_add_fmod (out a b : ZZ) :
    b * ZZ.fdiv_q out a b + ZZ.fdiv_r out a b = a :=
  Int.mul_fdiv_add_fmod a b

-- spec：余数非负界（仅 b > 0）
theorem ZZ.fdiv_r_nonneg (out a b : ZZ) (hb : 0 < b) :
    0 ≤ ZZ.fdiv_r out a b :=
  Int.fmod_nonneg_of_pos a hb

theorem ZZ.fdiv_r_lt (out a b : ZZ) (hb : 0 < b) :
    ZZ.fdiv_r out a b < b :=
  Int.fmod_lt_of_pos a hb

-- #eval 数值验证
#eval ZZ.fdiv_q 0 7 2          -- 3
#eval ZZ.fdiv_q 0 (-7) 2       -- -4 (floor, GMP 一致；非 -3)
#eval ZZ.fdiv_r 0 7 2          -- 1
#eval ZZ.fdiv_r 0 (-7) 2       -- 1 (在 [0, 2) 内)

-- ZZ = Int alias 时 `(x : ZZ).toInt` 不合法（Int 没 .toInt）。
-- Pass 5 cast_table 在某些 ZZ → Int 路径仍 emit `.toInt`；提供 identity 兜底。
def Int.toInt (x : Int) : Int := x

-- C++ log(x : double) → Lean Float.log（自然对数；Float stdlib 内置，无需 def）
-- cpp2lean Pass 5 已改 emit `Float.log` 直接调用；旧 `Nat.log` 占位被移除以避免
-- 与 Mathlib `Nat.log : Nat → Nat → Nat` 名冲突。

#eval Float.log (1.0 : Float)        -- 应为 0.0
#eval Float.log (2.71828182845904523536 : Float)  -- 应接近 1.0
#eval Float.log (10.0 : Float)       -- 应接近 2.302585...
#eval Float.log (100.0 : Float)      -- 应接近 4.605170...
def Int.toFloat (n : Int) : Float := Float.ofInt n

-- ZZ.sizeinbase: GMP `mpz_sizeinbase(z, base)` — |z| 在 base 进制下的位数。
-- 数学定义（base ≥ 2, z ≠ 0）：sizeinbase = ⌊log_base |z|⌋ + 1
-- GMP 约定：z = 0 时返回 1。
-- 退化 case（base ≤ 1）数学未定义；此处返回 1 兜底。
--
-- 不复用 Mathlib `Nat.log` 因为名字与 cpp2lean Pass 5 emit 冲突。本地实现 + 证明。

-- 整数对数：最大 k 使得 b^k ≤ n（要求 b ≥ 2, n ≥ 1）
-- 退化：b ≤ 1 或 n = 0 → 0
def Nat.intLog (b n : Nat) : Nat :=
  if h : 2 ≤ b ∧ b ≤ n then Nat.intLog b (n / b) + 1
  else 0
termination_by n
decreasing_by
  simp_wf
  exact Nat.div_lt_self (by omega) (by omega)

-- spec: b^(intLog b n) ≤ n （n ≥ 1, b ≥ 2 时）
theorem Nat.pow_intLog_le {b n : Nat} (hb : 2 ≤ b) (hn : 0 < n) :
    b ^ (Nat.intLog b n) ≤ n := by
  induction n using Nat.strongRecOn with
  | _ n ih =>
    rw [Nat.intLog]
    by_cases h : 2 ≤ b ∧ b ≤ n
    · simp [h]
      have hbpos : 0 < b := by omega
      have hpos_n : 0 < n := hn
      have hdiv_lt : n / b < n := Nat.div_lt_self hpos_n h.1
      have hdiv_pos : 0 < n / b := Nat.div_pos h.2 hbpos
      have ih' : b ^ Nat.intLog b (n / b) ≤ n / b := ih (n / b) hdiv_lt hdiv_pos
      calc b ^ (Nat.intLog b (n / b) + 1)
          = b ^ Nat.intLog b (n / b) * b := by rw [Nat.pow_succ]
        _ ≤ (n / b) * b := Nat.mul_le_mul_right b ih'
        _ ≤ n := Nat.div_mul_le_self n b
    · simp only [h, ↓reduceDIte, Nat.pow_zero]
      -- b^0 = 1 ≤ n （由 0 < n 得）
      exact hn

-- spec: n < b^(intLog b n + 1) （b ≥ 2 时；n=0 也成立 因 0 < b）
theorem Nat.lt_pow_succ_intLog {b n : Nat} (hb : 2 ≤ b) :
    n < b ^ (Nat.intLog b n + 1) := by
  induction n using Nat.strongRecOn with
  | _ n ih =>
    rw [Nat.intLog]
    by_cases h : 2 ≤ b ∧ b ≤ n
    · simp [h]
      have hbpos : 0 < b := by omega
      have hpos : 0 < n := by omega
      have hdiv_lt : n / b < n := Nat.div_lt_self hpos h.1
      have ih' : n / b < b ^ (Nat.intLog b (n / b) + 1) := ih (n / b) hdiv_lt
      -- n < (n / b + 1) * b （由 div_add_mod + mod_lt 得）
      have hmod : n % b < b := Nat.mod_lt n hbpos
      -- 注意 div_add_mod 是 b * (n/b) + n%b = n；交换乘法顺序
      have hdivmod : (n / b) * b + n % b = n := by
        rw [Nat.mul_comm]; exact Nat.div_add_mod n b
      have hub : n < (n / b + 1) * b := by
        rw [Nat.add_mul, Nat.one_mul]
        omega
      calc n < (n / b + 1) * b := hub
        _ ≤ b ^ (Nat.intLog b (n / b) + 1) * b := Nat.mul_le_mul_right b ih'
        _ = b ^ (Nat.intLog b (n / b) + 1 + 1) := (Nat.pow_succ b _).symm
    · simp only [h, ↓reduceDIte, Nat.pow_succ, Nat.pow_zero, Nat.one_mul]
      -- 目标 n < b。h ¬(2 ≤ b ∧ b ≤ n)，hb 2 ≤ b → ¬(b ≤ n) → n < b
      exact Nat.lt_of_not_ge (fun hbn => h ⟨hb, hbn⟩)

-- _nat 中间层：避开 UInt64 转换噪音（spec 都在 Nat 层）
def ZZ.sizeinbase_nat (z : ZZ) (base : Nat) : Nat :=
  if z.natAbs = 0 then 1
  else Nat.intLog base z.natAbs + 1

def ZZ.sizeinbase (z : ZZ) (base : Int32) : UInt64 :=
  (ZZ.sizeinbase_nat z base.toInt64.toNatClampNeg).toUInt64

/-- 对称模运算 a mod m → [-m/2, m/2] -/
def ZZ.symmetricMod (a m : ZZ) : ZZ :=
  let r := a.fmod m
  if r * 2 ≤ m then r else r - m

/-- 二项式系数 C(n, k) — 由于 Lean 标准库可能无 Nat.choose，手动实现 -/
def Nat.myChoose : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ => 0
  | n+1, k+1 => Nat.myChoose n (k+1) + Nat.myChoose n k

/-- 整数二项式系数 -/
def ZZ.binomial (n k : ZZ) : ZZ :=
  if 0 ≤ k ∧ k ≤ n then
    (Nat.myChoose n.natAbs k.natAbs : ZZ)
  else 0

/-- 整数向上取整平方根 -/
noncomputable def ZZ.isqrtCeil (n : ZZ) : ZZ :=
  if n ≤ 0 then 0
  else
    -- TODO: implement Nat.sqrt for ZZ
    n

-- spec 1: 总是 ≥ 1
theorem ZZ.sizeinbase_nat_pos (z : ZZ) (b : Nat) :
    1 ≤ ZZ.sizeinbase_nat z b := by
  unfold ZZ.sizeinbase_nat
  split <;> omega

-- spec 2: 当 z ≠ 0 且 base ≥ 2，base^(sizeinbase - 1) ≤ |z|
theorem ZZ.pow_pred_sizeinbase_nat_le {z : ZZ} {b : Nat}
    (hz : z ≠ 0) (hb : 2 ≤ b) :
    b ^ (ZZ.sizeinbase_nat z b - 1) ≤ z.natAbs := by
  have hn : z.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr hz
  have hpos : 0 < z.natAbs := Nat.pos_of_ne_zero hn
  unfold ZZ.sizeinbase_nat
  simp [hn]
  exact Nat.pow_intLog_le hb hpos

-- spec 3: 当 z ≠ 0 且 base ≥ 2，|z| < base^sizeinbase
theorem ZZ.lt_pow_sizeinbase_nat {z : ZZ} {b : Nat}
    (hz : z ≠ 0) (hb : 2 ≤ b) :
    z.natAbs < b ^ ZZ.sizeinbase_nat z b := by
  have hn : z.natAbs ≠ 0 := Int.natAbs_ne_zero.mpr hz
  unfold ZZ.sizeinbase_nat
  simp [hn]
  exact Nat.lt_pow_succ_intLog hb

-- #eval 数值验证
#eval ZZ.sizeinbase 0 2          -- 1
#eval ZZ.sizeinbase 1 2          -- 1
#eval ZZ.sizeinbase 2 2          -- 2
#eval ZZ.sizeinbase 7 2          -- 3
#eval ZZ.sizeinbase 8 2          -- 4
#eval ZZ.sizeinbase (-100) 10    -- 3
-- SparsePolyZZ.size_u64 见 §5c (abbrev 之后)
-- QQ = Rat 的 .num / .den 别名（Pass 5 emit `QQ.num q`，需要显式 const）
def QQ.num (q : QQ) : Int := Rat.num q
def QQ.den (q : QQ) : Int := (Rat.den q : Int)
def QQ.mk (n : Int) (d : Int) : QQ :=
  if d = 0 then 0 else (n : Rat) / (d : Rat)
def QQ.ofInt (n : Int) : QQ := (n : Rat)

-- 阶段 F #3 后续 — Lean 端 cast / API 占位（此前被 LambdaRef 错误屏蔽）
def Int64.toNat (i : Int64) : Nat := i.toNatClampNeg
def UInt64.toInt (u : UInt64) : Int := u.toNat
def Nat.toNat (n : Nat) : Nat := n  -- identity（cast_table 偶尔多余加的）

-- 兼容 vec.resize(n) 和 vec.resize(n, val) 两种 overload：默认 v 用 Inhabited
def Array.resize {α : Type} [Inhabited α] (a : Array α) (n : Nat) (v : α := default) : Array α :=
  if n ≤ a.size then a.extract 0 n
  else a ++ Array.replicate (n - a.size) v
def Array.getLast! {α : Type} [Inhabited α] (a : Array α) : α := a.back!
def Array.head! {α : Type} [Inhabited α] (a : Array α) : α := a[0]!

def Variable.mk {α : Type} (x : α) [HasPolyMk Variable α] : Variable := HasPolyMk.mkImpl x
def UniformIntDist.mk (_lo _hi : Int32) : UniformIntDist := 0
def Rng.default : Rng := 42
-- 阶段 G+：Rng = UInt64 abbrev，Pass 5/8 偶尔在某些上下文 emit `.toUInt64`
-- 让 abbrev Rng → UInt64 自动通过；定义 identity 避免 invalid field 错
def Rng.toUInt64 (r : Rng) : UInt64 := r
def UInt64.toUInt64 (u : UInt64) : UInt64 := u

-- Iterator: C++ STL iterator 的 Lean 抽象。语义层只区分"指向有效元素"与"end"。
-- true = valid iterator（指向有效元素），false = end sentinel（"未找到"）。
-- 这样 `it == m.end()` 在 Lean 端 ↔ found == false ↔ 未找到，与 C++ 语义一致。
abbrev Iterator := Bool
def Iterator.fromList {α : Type} (_a : Array α) : Iterator := true
-- MvMonomial 与 Monomial 同构（都是 Array (Variable × Int64)）；规范化复用 Monomial.normalize
def MvMonomial.normalization (m : MvMonomial) : MvMonomial := Monomial.normalize m
def gcd (a b : Int) : Int := Int.gcd a b
-- polynomial_mod(f : SparsePolyZZ, p : UInt64) : SparsePolyZp
-- 系数 mod p 把 ZZ 多项式变成 Zp 多项式
-- polynomial_mod: SparsePolyZZ + p → SparsePolyZp（实现移到 abbrev SparsePolyZZ 之后）
-- 见下方 §5c 末尾

-- next_prime_64: 返回 > n 的最小素数（C++ clpoly/number/ZZ.hh:1036, GMP mpz_nextprime 包装）
-- Lean 端用 trial division（O(√n) 单次），足够 B2B 测试小素数场景
def Nat.isPrime64 (n : Nat) : Bool :=
  if n < 2 then false
  else if n = 2 then true
  else if n % 2 = 0 then false
  else
    let rec loop (d : Nat) (fuel : Nat) : Bool :=
      match fuel with
      | 0 => true
      | fuel'+1 => if d * d > n then true
                   else if n % d = 0 then false
                   else loop (d + 2) fuel'
    loop 3 n
def next_prime_64 (p : UInt64) : UInt64 :=
  let rec scan (cand : Nat) (fuel : Nat) : Nat :=
    match fuel with
    | 0 => cand
    | fuel'+1 => if Nat.isPrime64 cand then cand else scan (cand + 1) fuel'
  (scan (p.toNat + 1) 100000).toUInt64
def prev_prime_64 (p : UInt64) : UInt64 := if p > 0 then p - 1 else 0
-- leadcoeff: 1-arg / 2-arg overload (Pass 5 emit 都用同一名)
-- 1-arg `leadcoeff p` 返回 ZZ；2-arg `leadcoeff p var` 返回 Poly
-- Lean 端：2-arg 版本（多变量主用），1-arg 用 leadcoeff1 区分
-- leadcoeff (2-arg): 关于 var 的首项系数（剥离 var 后的剩余多项式）
class HasLeadcoeff (α : Type) where
  leadcoeffImpl : α → Variable → α

-- 通用兜底：返回 default
instance (priority := 0) {α : Type} [Inhabited α] : HasLeadcoeff α where
  leadcoeffImpl _ _ := default

def leadcoeff {α : Type} [HasLeadcoeff α] (p : α) (var : Variable) : α :=
  HasLeadcoeff.leadcoeffImpl p var

-- MvPoly 实例见 §5c 末尾（abbrev 之后）
def leadcoeff1 {α : Type} [Inhabited α] (_p : α) : ZZ := 0
-- ZZ.fdiv_ui: GMP mpz_fdiv_ui(a, b) — 返回 a mod b 的非负残余（≤ b - 1）
-- 数学定义：与 fdiv_r 相同语义，被除数同 fdiv_r；只是除数是 UInt64 而非 ZZ。
-- C++ 调用点：clpoly::Zp(ZZ, prime) 用此把 ZZ 折回 Zp 域。
def ZZ.fdiv_ui (a : ZZ) (b : UInt64) : ZZ := Int.fmod a (b.toNat : Int)

-- spec：当 b > 0 时 fdiv_ui 落在 [0, b)
theorem ZZ.fdiv_ui_nonneg (a : ZZ) (b : UInt64) (hb : (b.toNat : Int) > 0) :
    0 ≤ ZZ.fdiv_ui a b :=
  Int.fmod_nonneg_of_pos a hb

theorem ZZ.fdiv_ui_lt (a : ZZ) (b : UInt64) (hb : (b.toNat : Int) > 0) :
    ZZ.fdiv_ui a b < (b.toNat : Int) :=
  Int.fmod_lt_of_pos a hb

-- #eval 数值验证
#eval ZZ.fdiv_ui 100 (7 : UInt64)     -- 2
#eval ZZ.fdiv_ui (-1) (7 : UInt64)    -- 6
#eval ZZ.fdiv_ui 0 (13 : UInt64)      -- 0
-- StdMap.find / .end 返回 Iterator（Bool 语义）
-- find m k = isSome (find? m k)（true 表已存在，false 表未找到）
-- end _   = false（end sentinel）
-- 故 `it == m.end()` ↔ 未找到，与 C++ STL `std::map::find` 语义一致
def StdMap.find {κ ν : Type} [BEq κ] [Inhabited ν] (m : StdMap κ ν) (k : κ) : Iterator :=
  (StdMap.find? m k).isSome
def StdMap.end {κ ν : Type} (_ : StdMap κ ν) : Iterator := false

-- 阶段 G-E：补 corpus 还需要的 stub 占位
def MvPolyZp.size_u64 (f : MvPolyZp) : UInt64 := (Array.size f).toUInt64
-- 真实 SparsePolyZZ.normalization 见 §5c.cont（abbrev 之后），不在此处放占位
-- （以前此处放 `def SparsePolyZZ.normalization (f) := f` 是 autoImplicit silent bug）
-- Array.range_init: 多 arity overload，C++ 写法 `iota(arr.begin(), arr.end(), start)`
-- Pass 5 emit 通常是 (arr, start) 2-arg；arr 决定大小，start 是初值
def Array.range_init {α : Type} (a : Array α) (_start : Int32) : Array Int32 :=
  (Array.range a.size).map (·.toUInt32.toInt32)

-- 注：旧 corpus 中曾出现 `compute_theta` / `upzp_coeff` / `next_p` 裸名引用（Pass 3
-- lift 不完整的产物）。当前 corpus 全部通过 `_1` lambda-arg 注入；不再需要全局占位。
-- cont(poly) → output：多项式整数系数的 content
-- C++ 端 `cont` 是 template，输入决定输出：
--   - cont(upolynomial_<ZZ>) → ZZ  （单变量，对应 SparsePolyZZ）
--   - cont(polynomial_<ZZ, lex_<var_order>>) → polynomial_<ZZ, lex_<var_order>>
--                                              （多变量，对应 MvPolyZZ；cont 在剩余变量中）
-- Lean 端用 associated output type 模拟：HasCont.Out 决定输出类型，默认 ZZ。
-- 注意：SparsePolyZZ.{cont,pp}Impl 实现移到 abbrev SparsePolyZZ 之后（见下方）
class HasCont (α : Type) where
  Out : Type := ZZ
  cont : α → Out

class HasPP (α : Type) where
  pp : α → α

-- 全局调度入口（与 cpp2lean Pass 5 emit 一致）
def cont {α : Type} [c : HasCont α] (p : α) : c.Out := c.cont p
def pp {α : Type} [HasPP α] (p : α) : α := HasPP.pp p

-- MvPolyZZ: cont = non-negative gcd of all int coeffs (issue #17 / C++ PR #18); pp = f / cont
def MvPolyZZ.contImpl (f : MvPolyZZ) : ZZ :=
  if f.isEmpty then 0
  else
    let c_nat : Nat := f.foldl (fun (acc : Nat) (term : Monomial × Int) =>
      Nat.gcd acc term.snd.natAbs) 0
    (c_nat : Int)

def MvPolyZZ.ppImpl (f : MvPolyZZ) : MvPolyZZ :=
  let c : ZZ := MvPolyZZ.contImpl f
  if c = 0 then f
  else f.map (fun term => (term.fst, term.snd / c))

-- C++ 多变量 cont 返回 polynomial（关于剩余变量的 GCD）。
-- 对 univariate-as-mvpoly 输入：剩余变量集合为空，"系数多项式 GCD" 退化为
-- 整数 GCD，故"整数 cont 包成常数 polynomial"是数学正确的。
-- 对真正多变量输入：是不完整占位（应返回剩余变量的多项式）。Phase F-impl-v2 续。
def MvPolyZZ.contMv (f : MvPolyZZ) : MvPolyZZ :=
  let c : ZZ := MvPolyZZ.contImpl f
  if c = 0 then #[]
  else #[(#[], c)]

instance : HasCont MvPolyZZ where
  Out := MvPolyZZ
  cont := MvPolyZZ.contMv
instance : HasPP MvPolyZZ where pp := MvPolyZZ.ppImpl

-- HDiv MvPolyZZ MvPolyZZ MvPolyZZ：真实现见 §5e 末尾（abbrev 之后的 MvPolyZZ.divExact）

-- MvPolyZp / SparsePolyZp（Zp 是域）: cont 类型为 ZZ，仅起 unit 标记作用
-- 约定：cont = (isEmpty ? 0 : 1)，pp = f（无需提主因子，仅在 ZZ 上有意义）
instance : HasCont MvPolyZp where
  cont f := if f.isEmpty then 0 else 1
instance : HasPP MvPolyZp where pp f := f
instance : HasCont SparsePolyZp where
  cont f := if f.isEmpty then 0 else 1
instance : HasPP SparsePolyZp where pp f := f
-- 依赖 SparsePolyZZ / LLLMatrix abbrev：见 §5c (abbrev 之后)
-- C++ std::swap(a, b)：值语义返回 (b, a) 元组（ref-elim 已转 SSA）
def swap {α β : Type} (a : α) (b : β) : β × α := (b, a)

-- Pass 4 filter-loop 转 `Array.filter(arr, pred)` —— Lean 4 Array.filter 期望
-- (pred, arr) 顺序。提供 (arr, pred) 包装。
namespace Array
def filter' {α : Type} (a : Array α) (p : α → Bool) : Array α := Array.filter p a
def filterMap' {α β : Type} (a : Array α) (f : α → Option β) : Array β := Array.filterMap f a
end Array

-- Coe Int32 → UInt64 / Int64（Pass 1 把 C++ 字面量识别为 Int32，Lean 端
-- 函数参数常需 UInt64/Int64；自动 Coe 解决 ~5 处 Application mismatch）
instance : Coe Int32 UInt64 where coe n := n.toInt64.toUInt64
-- 阶段 G+ 修复：C++ `{g, 1}` 嵌入 tuple，1 是 int 但目标 uint64
-- Lean tuple 元素 Coe 不自动传播，定义具体 pair Coe
instance : Coe (MvPolyZZ × Int32) (MvPolyZZ × UInt64) where
  coe p := (p.fst, p.snd.toInt64.toUInt64)
-- 同样为 Array (MvPolyZZ × Int32) → Array (MvPolyZZ × UInt64)
instance : Coe (Array (MvPolyZZ × Int32)) (Array (MvPolyZZ × UInt64)) where
  coe a := a.map (fun (m, i) => (m, i.toInt64.toUInt64))
-- Pass 1 在 Zp context 误把 `Poly` 映射到 MvPolyZZ；某些 site 实际期望 MvPolyZp。
-- 历史上有 lossy `coe _ := #[]` escape hatch — 但那是 silent dropping bug。
-- 改为 panic：让任何隐式触发点立刻爆出来，逼 Pass 1 修类型推断 / 调用方加 explicit cast。
-- 若 corpus 编译失败，说明真有触发点；需修 Pass 1 而不是回退到丢值的 Coe。
instance : Coe MvPolyZZ MvPolyZp where
  coe _ := panic! "Coe MvPolyZZ→MvPolyZp triggered: Pass 1 type inference broken; add explicit conversion at call site"
instance : Coe MvPolyZp MvPolyZZ where
  coe _ := panic! "Coe MvPolyZp→MvPolyZZ triggered: Pass 1 type inference broken; add explicit conversion at call site"
instance : Coe Int32 Int64 where coe n := n.toInt64
instance : Coe UInt64 Nat where coe n := n.toNat
instance : Coe Int64 Nat where coe n := n.toNatClampNeg
-- 阶段 G+：Pass 5 cast 漏的 site 用 Lean Coe 自动桥（uni-directional safe casts）
instance : Coe Nat Int64 where coe n := n.toUInt64.toInt64
instance : Coe Nat Int32 where coe n := n.toUInt32.toInt32
instance : Coe UInt64 Int64 where coe u := u.toInt64
instance : Coe UInt32 UInt64 where coe u := u.toUInt64
instance : Coe ZZ Nat where coe z := z.toNat

/-- For UInt64 values, the roundtrip through Int64 to Nat (via toNatClampNeg) is bounded above
    by the direct UInt64.toNat. This is universally true because toNatClampNeg clamps negative
    Int64 values to 0, and returns the original Nat for non-negative values.
    Used in SquarefreeZp to bound degrees after UInt64 division. -/
theorem UInt64_toInt64_toNatClampNeg_le_toNat (u : UInt64) : u.toInt64.toNatClampNeg ≤ u.toNat := by
  -- toNatClampNeg = toInt.toNat；toInt64 保位，故 toInt = (toNat).bmod 2^64。
  show (u.toInt64.toInt).toNat ≤ u.toNat
  have h1 : u.toInt64.toInt = ((u.toNat : Int)).bmod (2 ^ 64) := by
    show u.toBitVec.toInt = _
    rw [BitVec.toInt_eq_toNat_bmod, UInt64.toNat_toBitVec]
  rw [h1, Int.bmod_def]; omega

/-- For UInt64 values less than 2^63, the roundtrip through Int64 to Nat preserves the value.
    This is because toInt64 maps the value identically and toNatClampNeg returns the Nat value
    for non-negative Int64. -/
theorem UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt {u : UInt64} (h : u.toNat < 2 ^ 63) : u.toInt64.toNatClampNeg = u.toNat := by
  -- u.toNat < 2^63 ⇒ bmod 落在正区间，等于原值。
  show (u.toInt64.toInt).toNat = u.toNat
  have h1 : u.toInt64.toInt = ((u.toNat : Int)).bmod (2 ^ 64) := by
    show u.toBitVec.toInt = _
    rw [BitVec.toInt_eq_toNat_bmod, UInt64.toNat_toBitVec]
  rw [h1, Int.bmod_def]; omega
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

-- IsNumber / HasDegree instance：SparsePolyZZ abbrev 之后才能定义
instance : IsNumber SparsePolyZZ where isNumber f := (f : Array _).isEmpty
instance : IsNumber SparsePolyZp where isNumber f := (f : Array _).isEmpty
-- 单变量稀疏多项式：约定首项已是 leading（按 deg 降序），degree = front.deg
instance : HasDegree SparsePolyZZ where
  degree f := if f.isEmpty then 0 else f[0]!.fst.deg.toUInt64
instance : HasDegree SparsePolyZp where
  degree f := if f.isEmpty then 0 else f[0]!.fst.deg.toUInt64

-- OfNat 0 实例：C++ `SparsePolyZZ x = 0` → 空多项式
instance : OfNat SparsePolyZZ 0 where ofNat := #[]
instance : OfNat SparsePolyZp 0 where ofNat := #[]
instance : OfNat MvPolyZZ 0 where ofNat := #[]
instance : OfNat MvPolyZp 0 where ofNat := #[]

def SparsePolyZZ.size_u64 (f : SparsePolyZZ) : UInt64 := f.size.toUInt64

-- §5c 续：SparsePolyZZ.cont / pp 实现（在 abbrev 之后才能用 .map / .foldl 等）
-- cont = gcd 所有系数（始终非负，issue #17 / C++ PR #18 对齐）
def SparsePolyZZ.contImpl (f : SparsePolyZZ) : ZZ :=
  if f.isEmpty then 0
  else
    let c_nat : Nat := f.foldl (fun (acc : Nat) (term : UMonomial × Int) =>
      Nat.gcd acc term.snd.natAbs) 0
    (c_nat : Int)

-- pp = f / cont(f)
def SparsePolyZZ.ppImpl (f : SparsePolyZZ) : SparsePolyZZ :=
  let c : ZZ := SparsePolyZZ.contImpl f
  if c = 0 then f
  else f.map (fun term => (term.fst, term.snd / c))

instance : HasCont SparsePolyZZ where cont := SparsePolyZZ.contImpl
instance : HasPP SparsePolyZZ where pp := SparsePolyZZ.ppImpl

-- §5c.cont SparsePolyZZ 真实 derivative + normalization（abbrev 之后才能用 Array API）
def SparsePolyZZ.derivativeImpl (f : SparsePolyZZ) : SparsePolyZZ :=
  f.filterMap (fun term =>
    if term.fst.deg = 0 then none
    else
      let new_c : Int := term.snd * term.fst.deg
      if new_c = 0 then none
      else some (⟨term.fst.deg - 1⟩, new_c))

instance : HasDerivative SparsePolyZZ where
  derivative := SparsePolyZZ.derivativeImpl

-- normalization：按 deg 降序、合并同 deg、剔零
def SparsePolyZZ.normalization (f : SparsePolyZZ) : SparsePolyZZ :=
  -- step 1: group by deg, summing coefs (O(n²) but simple)
  let grouped : SparsePolyZZ := f.foldl (fun acc term =>
    match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
    | some idx => acc.modify idx (fun (m, c) => (m, c + term.snd))
    | none => acc.push term) #[]
  -- step 2: drop zero coefficients
  let nonZero : SparsePolyZZ := grouped.filter (fun t => t.snd ≠ 0)
  -- step 3: sort descending by deg
  nonZero.qsort (fun a b => a.fst.deg > b.fst.deg)

-- polynomial_mod: SparsePolyZZ + p → SparsePolyZp
-- 数学定义：每个系数 mod p（用 Zp.ofInt 折回 [0, p)），剔除 0 系数
def polynomial_mod (f : SparsePolyZZ) (p : UInt64) : SparsePolyZp :=
  f.filterMap (fun term =>
    let zp := Zp.ofInt term.snd p
    if zp.val = 0 then none else some (term.fst, zp))

-- ============================================================
-- SparsePolyZZ.gcd（单变量 ZZ primitive PRS GCD）—— Phase F-impl-X.4
-- ============================================================

-- 标量乘
def SparsePolyZZ.scalarMul (c : Int) (f : SparsePolyZZ) : SparsePolyZZ :=
  if c = 0 then #[]
  else f.map (fun (m, x) => (m, x * c))

-- 移位乘：f * x^k
def SparsePolyZZ.shiftMul (k : Nat) (f : SparsePolyZZ) : SparsePolyZZ :=
  f.map (fun (m, c) => (⟨m.deg + k⟩, c))

-- 真加法（合并同 deg + 排序）
def SparsePolyZZ.addReal (f g : SparsePolyZZ) : SparsePolyZZ :=
  SparsePolyZZ.normalization (f ++ g)

-- 真减法
def SparsePolyZZ.subReal (f g : SparsePolyZZ) : SparsePolyZZ :=
  let neg_g := g.map (fun (m, c) => (m, -c))
  SparsePolyZZ.normalization (f ++ neg_g)

-- 伪余数：lc(G)^k * F mod G, k = deg(F) - deg(G) + 1
-- 通过迭代单步消首项实现（每步乘 lc(G) 然后消首项）
partial def SparsePolyZZ.pseudoRem (F G : SparsePolyZZ) : SparsePolyZZ :=
  if F.isEmpty then #[]
  else if G.isEmpty then F  -- div by 0 placeholder
  else
    let dF := F[0]!.fst.deg
    let dG := G[0]!.fst.deg
    if dF < dG then F
    else
      let lcG := G[0]!.snd
      let lcF := F[0]!.snd
      let k := dF - dG
      let F_scaled := SparsePolyZZ.scalarMul lcG F
      let shifted := SparsePolyZZ.scalarMul lcF (SparsePolyZZ.shiftMul k G)
      let newF := SparsePolyZZ.subReal F_scaled shifted
      SparsePolyZZ.pseudoRem newF G

-- 基本 PRS：交替伪余数 + 取 primitive part
partial def SparsePolyZZ.primitivePRS (F G : SparsePolyZZ) : SparsePolyZZ :=
  if G.isEmpty then F
  else
    let R := SparsePolyZZ.pseudoRem F G
    if R.isEmpty then G
    else
      let Rp := SparsePolyZZ.ppImpl R
      SparsePolyZZ.primitivePRS G Rp

-- 单变量 ZZ GCD
def SparsePolyZZ.gcd (F G : SparsePolyZZ) : SparsePolyZZ :=
  if F.isEmpty then G
  else if G.isEmpty then F
  else
    let cF := SparsePolyZZ.contImpl F
    let cG := SparsePolyZZ.contImpl G
    -- 整数 gcd（取绝对值 + 符号约定为首项正）
    let d : Int := (Int.gcd cF cG : Int)
    let F0 := SparsePolyZZ.ppImpl F
    let G0 := SparsePolyZZ.ppImpl G
    let primGCD := SparsePolyZZ.primitivePRS F0 G0
    -- 首项规范化（issue #14 一致：factor 永远首项正）
    if h : 0 < primGCD.size then
      let lead := primGCD[0].snd
      if lead < 0 then
        SparsePolyZZ.scalarMul (-d) primGCD
      else
        SparsePolyZZ.scalarMul d primGCD
    else
      #[]

-- #eval 验证 derivative + normalization
-- d/dx (3x² + 5x + 7) = 6x + 5
#eval (SparsePolyZZ.derivativeImpl
  (#[(⟨2⟩, (3 : Int)), (⟨1⟩, (5 : Int)), (⟨0⟩, (7 : Int))] : SparsePolyZZ))
-- 期望: #[(1, 6), (0, 5)]

-- normalization: [(0,1), (2,3), (1,2), (2,4)] → [(2,7), (1,2), (0,1)]
#eval (SparsePolyZZ.normalization
  (#[(⟨0⟩, (1 : Int)), (⟨2⟩, (3 : Int)), (⟨1⟩, (2 : Int)), (⟨2⟩, (4 : Int))] : SparsePolyZZ))
-- 期望: #[(2, 7), (1, 2), (0, 1)]

-- #eval 验证：(2x² + 3x + 5) mod 5 = 2x² + 3x （5 mod 5 = 0 剔除）
#eval polynomial_mod
  (#[(⟨2⟩, (2 : Int)), (⟨1⟩, (3 : Int)), (⟨0⟩, (5 : Int))] : SparsePolyZZ) 5
-- 期望: #[(2, [2, 5]), (1, [3, 5])]

-- #eval 验证：cont(2x² + 4x + 6) = 2，pp = x² + 2x + 3
#eval SparsePolyZZ.contImpl
  (#[(⟨2⟩, (2 : Int)), (⟨1⟩, (4 : Int)), (⟨0⟩, (6 : Int))] : SparsePolyZZ)
-- 期望: 2

#eval SparsePolyZZ.ppImpl
  (#[(⟨2⟩, (2 : Int)), (⟨1⟩, (4 : Int)), (⟨0⟩, (6 : Int))] : SparsePolyZZ)
-- 期望: #[(2, 1), (1, 2), (0, 3)]

-- cont(-2x² - 4) = 2 (issue #17 / PR #18：始终非负)
#eval SparsePolyZZ.contImpl
  (#[(⟨2⟩, (-2 : Int)), (⟨0⟩, (-4 : Int))] : SparsePolyZZ)
-- 期望: 2

-- issue #17 回归用例：cont 在单项 / 常数 / 负 leading 情形始终非负
-- cont(-x) = 1
#eval SparsePolyZZ.contImpl
  (#[(⟨1⟩, (-1 : Int))] : SparsePolyZZ)
-- 期望: 1

-- cont(-6 const) = 6
#eval SparsePolyZZ.contImpl
  (#[(⟨0⟩, (-6 : Int))] : SparsePolyZZ)
-- 期望: 6

-- pp(-2x² - 4) = -x² - 2 （Maple/SymPy/FLINT 约定：pp 首项可能为负）
#eval SparsePolyZZ.ppImpl
  (#[(⟨2⟩, (-2 : Int)), (⟨0⟩, (-4 : Int))] : SparsePolyZZ)
-- 期望: #[(2, -1), (0, -2)]

-- 阶段 F 后续：依赖 SparsePolyZZ 的 stub（LLLMatrix.size 见 abbrev 之后）
-- get_first_deg: 多变量 / 单变量两态。Lean 端泛型占位（语义层 B2B 细化）
-- get_first_deg：典型用法是取多变量首项关于第一个变量的 deg
-- 对 MvPoly：f.front!.fst.front!.snd（首单项的首变量 exp）
-- 对其他类型：兜底 0
class HasFirstDeg (α : Type) where
  firstDeg : α → Int64

instance : HasFirstDeg MvPolyZZ where
  firstDeg (f : MvPolyZZ) : Int64 :=
    if f.isEmpty then 0
    else
      let mono := f[0]!.fst
      if mono.isEmpty then 0
      else mono[0]!.snd

instance : HasFirstDeg MvPolyZp where
  firstDeg (f : MvPolyZp) : Int64 :=
    if f.isEmpty then 0
    else
      let mono := f[0]!.fst
      if mono.isEmpty then 0
      else mono[0]!.snd

instance (priority := 0) {α : Type} : HasFirstDeg α where
  firstDeg _ := 0

def get_first_deg {α : Type} [HasFirstDeg α] (f : α) : Int64 := HasFirstDeg.firstDeg f

-- get_deg: 泛型化（C++ side 多模板实例化共用同一 Lean 实现）
-- 适用 SparsePolyZZ / SparsePolyZp 两种容器（结构相同：Array (UMonomial × _)）
-- 返回 Int64（多数 Pass 1 调用点把 get_deg 视为 int64_t / signed comparison 上下文）
def get_deg {α : Type} [Inhabited α] (f : Array (UMonomial × α)) : Int64 :=
  if f.isEmpty then 0 else (f[0]!).fst.deg.toUInt64.toInt64

abbrev LLLMatrix := Array (Array Int)
def LLLMatrix.size (m : LLLMatrix) : UInt64 := (Array.size m).toUInt64
def LLLMatrix.empty : LLLMatrix := #[]
def LLLMatrix.replicate (n : UInt64) (row : Array Int) (_ : Unit := ()) : LLLMatrix :=
  Array.replicate n.toNat row

-- HenselNode: Hensel 提升二叉节点（C++ __hensel_node）
-- left / right: 子节点 Int32 索引（-1 表叶节点）
-- g, h, s, t: 多项式因子 / Bezout 系数（C++ side 用 SparsePolyZZ —— 模 m 的整数表示）
structure HenselNode where
  left : Int32 := -1
  right : Int32 := -1
  g : SparsePolyZZ := #[]
  h : SparsePolyZZ := #[]
  s : SparsePolyZZ := #[]
  t : SparsePolyZZ := #[]
  leaf_start : Int32 := 0
  leaf_end : Int32 := 0
deriving Inhabited

-- Pass 1 把 C++ aggregate init `HenselNode{g,h,s,t,left,right,ls,le}` emit 为
-- Array 字面量。提供 lossy coercion：取默认 HenselNode（语义层 B2B 测试细化）
instance : CoeHTCT (Array SparsePolyZZ) HenselNode where
  coe _ := default

-- ValueType: Pass 1 把 `typename Container::value_type` 在某些 corpus 路径上
-- 简化为 NamedType("ValueType")。在 Hensel 上下文 = HenselNode（Array.value_type）。
abbrev ValueType := HenselNode

-- A 方案：Factorization 参数化为 (PolyT : Type)，C++ `factorization<X>` 直接
-- emit 为 `Factorization X`。`Factorization.empty` 用泛型 inhabit 默认值
-- （PolyT 必须 [Inhabited]）。
structure Factorization (PolyT : Type) where
  content : ZZ := 0
  factors : Array (PolyT × UInt64) := #[]

instance {PolyT : Type} [Inhabited PolyT] : Inhabited (Factorization PolyT) where
  default := { }

def Factorization.empty {PolyT : Type} : Factorization PolyT :=
  { content := 0, factors := #[] }

structure PrimeSelectionResult where
  p : UInt64 := 0
  prime : UInt64 := 0
  factors : Array SparsePolyZp := #[]
  nfactors : UInt64 := 0
  irreducible : Bool := false
deriving Inhabited

-- 阶段 G9 续修：对照 C++ __wang_lc_result（polynomial_factorize_wang.hh:1314）
-- 5 字段 + delta（CLPoly 给 result.delta 的工程附加）
structure WangLcResult where
  success : Bool := false
  f_scaled : MvPolyZZ := #[]
  lc_assignments : Array MvPolyZZ := #[]
  lc_targets : Array MvPolyZZ := #[]
  scaled_factors : Array SparsePolyZZ := #[]
  delta : ZZ := 0
deriving Inhabited

-- assign：多项式变量代入 poly[var := val]
-- typeclass dispatch：MvPolyZp + Zp 真实实现，其他类型兜底
class HasAssign (α β : Type) where
  assignImpl : α → Variable → β → α

-- MvPolyZp + Zp: 真实 substitution + normalize
instance : HasAssign MvPolyZp Zp where
  assignImpl (f : MvPolyZp) (var : Variable) (val : Zp) : MvPolyZp :=
    let p := val.prime
    let one_zp : Zp := Zp.ofUInt64 1 p
    let substituted : MvPolyZp := f.map (fun (mono, c) =>
      let acc0 : Monomial × Zp := (#[], one_zp)
      let result := mono.foldl (fun (acc : Monomial × Zp) (entry : Variable × Int64) =>
        if entry.fst == var then
          (acc.fst, acc.snd * (val ^ entry.snd.toNatClampNeg))
        else
          (acc.fst.push entry, acc.snd)) acc0
      (result.fst, c * result.snd))
    MvPolyZp.normalization substituted

-- MvPolyZZ + ZZ: 真实 substitution + normalize
instance : HasAssign MvPolyZZ ZZ where
  assignImpl (f : MvPolyZZ) (var : Variable) (val : ZZ) : MvPolyZZ :=
    let substituted : MvPolyZZ := f.map (fun (mono, c) =>
      let acc0 : Monomial × ZZ := (#[], (1 : Int))
      let result := mono.foldl (fun (acc : Monomial × ZZ) (entry : Variable × Int64) =>
        if entry.fst == var then
          (acc.fst, acc.snd * (val ^ entry.snd.toNatClampNeg))
        else
          (acc.fst.push entry, acc.snd)) acc0
      (result.fst, c * result.snd))
    MvPolyZZ.normalization substituted

-- 通用兜底（identity）
instance (priority := 0) {α β : Type} : HasAssign α β where
  assignImpl p _ _ := p

def assign {α β : Type} [HasAssign α β] (poly : α) (var : Variable) (val : β) : α :=
  HasAssign.assignImpl poly var val

-- 2-arg overload：assign(poly, eval_point) 用 map 一次性代入多个变量
class HasAssign2 (α β : Type) where
  assign2Impl : α → StdMap Variable β → α

instance : HasAssign2 MvPolyZZ ZZ where
  assign2Impl (poly : MvPolyZZ) (eval : StdMap Variable ZZ) : MvPolyZZ :=
    eval.foldl (fun acc (v, val) => HasAssign.assignImpl acc v val) poly

instance : HasAssign2 MvPolyZp Zp where
  assign2Impl (poly : MvPolyZp) (eval : StdMap Variable Zp) : MvPolyZp :=
    eval.foldl (fun acc (v, val) => HasAssign.assignImpl acc v val) poly

instance (priority := 0) {α β : Type} : HasAssign2 α β where
  assign2Impl p _ := p

def assign2 {α : Type} (poly : α) (eval_point : StdMap Variable ZZ)
    [HasAssign2 α ZZ] : α :=
  HasAssign2.assign2Impl poly eval_point

-- §5e. leadcoeff (2-arg) / poly_convert / poly_convert3 具体实例
-- 这些实例必须在 abbrev MvPolyZZ/MvPolyZp/SparsePolyZZ 之后定义

-- 工具：取 mono 中关于 var 的 exp（缺失视作 0）
def Monomial.expOf (m : Monomial) (var : Variable) : Int64 :=
  match m.find? (fun e => e.fst == var) with
  | some e => e.snd
  | none => 0

-- 工具：strip mono 中的 var 条目
def Monomial.dropVar (m : Monomial) (var : Variable) : Monomial :=
  m.filter (fun e => e.fst != var)

-- m1 divides m2 ⇔ ∀ var v ∈ m1, m1[v].exp ≤ m2[v].exp（m2 缺失项视为 0）
def Monomial.divides (m1 m2 : Monomial) : Bool :=
  m1.all (fun (v, e) => e ≤ Monomial.expOf m2 v)

-- m2 / m1 = monomial 商：var-by-var exp 相减（假设 m1 divides m2）
def Monomial.div (m2 m1 : Monomial) : Monomial :=
  Monomial.normalize (m2.map (fun (v, e) => (v, e - Monomial.expOf m1 v)))

-- HDiv MvPolyZZ MvPolyZZ MvPolyZZ：exact division via leading-mono reduction
-- 假设 g | f（除尽）；否则返回不完整商（leading 不能整除时停）
-- 对应 C++ basic_polynomial::operator/(p1, p) = pair_vec_div(...)
-- C++ 是 heap-based 优化版；Lean 用经典 leading-mono 归约（语义等价）
partial def MvPolyZZ.divExact (f g : MvPolyZZ) : MvPolyZZ :=
  if g.isEmpty then f  -- div by 0 placeholder（数学未定义）
  else
    let rec loop (rem : MvPolyZZ) (acc : MvPolyZZ) : MvPolyZZ :=
      if rem.isEmpty then acc
      else
        let lc_r : Monomial × Int := rem[0]!
        let lc_g : Monomial × Int := g[0]!
        if Monomial.divides lc_g.fst lc_r.fst ∧ lc_r.snd % lc_g.snd = 0 then
          let q_mono : Monomial := Monomial.div lc_r.fst lc_g.fst
          let q_coef : Int := lc_r.snd / lc_g.snd
          let q_term : MvPolyZZ := #[(q_mono, q_coef)]
          let prod : MvPolyZZ := MvPolyZZ.mulImpl q_term g
          let new_rem : MvPolyZZ := MvPolyZZ.normalization (rem - prod)
          loop new_rem (acc.push (q_mono, q_coef))
        else
          acc  -- leading 无法整除（非 exact），返回当前商
    loop f #[]

-- 覆盖之前的 placeholder HDiv（identity）
instance : HDiv MvPolyZZ MvPolyZZ MvPolyZZ where
  hDiv := MvPolyZZ.divExact

-- ============================================================
-- MvPolyZZ.polynomialGCD（多变量 ZZ GCD）—— Phase F-impl-X.6
-- ============================================================

-- 整数 cont（所有系数 gcd 绝对值，始终非负，issue #17 / C++ PR #18 对齐）
def MvPolyZZ.contInt (f : MvPolyZZ) : Int :=
  if 0 < f.size then
    let nat_gcd : Nat := f.foldl (fun (acc : Nat) (t : Monomial × Int) =>
      Nat.gcd acc t.snd.natAbs) 0
    (nat_gcd : Int)
  else 0

-- 是否常数（所有 mono 为空）
def MvPolyZZ.isConstantMv (f : MvPolyZZ) : Bool :=
  f.all (fun t => t.fst.isEmpty)

-- 取首项的主变量（lex 最高）；首项是常数返回 none
def MvPolyZZ.getMainVar (f : MvPolyZZ) : Option Variable :=
  if h : 0 < f.size then
    let firstMono := f[0].fst
    if hm : 0 < firstMono.size then some firstMono[0].fst
    else none
  else none

-- 是否仅使用变量 v
def MvPolyZZ.isUnivariateInVar (f : MvPolyZZ) (v : Variable) : Bool :=
  f.all (fun t =>
    t.fst.all (fun e => e.fst == v))

-- univariate-in-v → SparsePolyZZ
def MvPolyZZ.mvToSparseZZ (v : Variable) (f : MvPolyZZ) : SparsePolyZZ :=
  let s : SparsePolyZZ := f.map (fun term =>
    let deg : Nat := (Monomial.expOf term.fst v).toNatClampNeg
    (⟨deg⟩, term.snd))
  SparsePolyZZ.normalization s

-- SparsePolyZZ → MvPolyZZ (单变量 v)
def MvPolyZZ.sparseZZToMv (v : Variable) (f : SparsePolyZZ) : MvPolyZZ :=
  f.map (fun term =>
    let mono : Monomial :=
      if term.fst.deg = 0 then #[]
      else #[(v, (term.fst.deg.toUInt64.toInt64 : Int64))]
    (mono, term.snd))

-- 多变量 polynomial_GCD：
--   - 任一空 → 另一个
--   - 任一常数 → 常数 GCD（整数 gcd 包成常数 poly）
--   - 同主变量且都 univariate → 转 SparsePolyZZ 调单变量 gcd 再转回
--   - 否则（真多变量）→ 退化到整数 cont gcd（已知不完整，TODO Phase F-impl-v2）
def MvPolyZZ.polynomialGCD (F G : MvPolyZZ) : MvPolyZZ :=
  if F.isEmpty then G
  else if G.isEmpty then F
  else if MvPolyZZ.isConstantMv F then
    let c : Int := F[0]!.snd
    let cInt := MvPolyZZ.contInt G
    let g : Int := (Nat.gcd c.natAbs cInt.natAbs : Int)
    if g = 0 then #[] else #[(#[], g)]
  else if MvPolyZZ.isConstantMv G then
    let c : Int := G[0]!.snd
    let cInt := MvPolyZZ.contInt F
    let g : Int := (Nat.gcd c.natAbs cInt.natAbs : Int)
    if g = 0 then #[] else #[(#[], g)]
  else
    match MvPolyZZ.getMainVar F, MvPolyZZ.getMainVar G with
    | some vF, some vG =>
      if vF = vG ∧ MvPolyZZ.isUnivariateInVar F vF ∧ MvPolyZZ.isUnivariateInVar G vG then
        -- univariate-as-mvpoly：转 SparsePolyZZ 调 univariate GCD
        let sF := MvPolyZZ.mvToSparseZZ vF F
        let sG := MvPolyZZ.mvToSparseZZ vG G
        let sGcd := SparsePolyZZ.gcd sF sG
        MvPolyZZ.sparseZZToMv vF sGcd
      else
        -- 真多变量退化：整数 cont gcd（不完整占位，TODO Phase F-impl-v2）
        let cF := MvPolyZZ.contInt F
        let cG := MvPolyZZ.contInt G
        let g : Int := (Nat.gcd cF.natAbs cG.natAbs : Int)
        if g = 0 then #[] else #[(#[], g)]
    | _, _ => #[]

instance : HasPolyGCD MvPolyZZ where polyGCD := MvPolyZZ.polynomialGCD

-- leadcoeff(p, var)：找 var 的最大 exp，收集这些项，剥离 var
instance : HasLeadcoeff MvPolyZZ where
  leadcoeffImpl (f : MvPolyZZ) (var : Variable) :=
    let exps : Array Int64 := f.map (fun t => Monomial.expOf t.fst var)
    let maxExp : Int64 := exps.foldl (fun acc e => if e > acc then e else acc) 0
    let lc : MvPolyZZ := f.filterMap (fun (m, c) =>
      if Monomial.expOf m var = maxExp then some (Monomial.dropVar m var, c) else none)
    MvPolyZZ.normalization lc

instance : HasLeadcoeff MvPolyZp where
  leadcoeffImpl (f : MvPolyZp) (var : Variable) :=
    let exps : Array Int64 := f.map (fun t => Monomial.expOf t.fst var)
    let maxExp : Int64 := exps.foldl (fun acc e => if e > acc then e else acc) 0
    let lc : MvPolyZp := f.filterMap (fun (m, c) =>
      if Monomial.expOf m var = maxExp then some (Monomial.dropVar m var, c) else none)
    MvPolyZp.normalization lc

-- SparsePolyZp → SparsePolyZZ: Zp.val.toInt → Int（Hensel 升提用）
instance : HasPolyConvert SparsePolyZp SparsePolyZZ where
  convert (f : SparsePolyZp) (_target : SparsePolyZZ) : SparsePolyZZ :=
    f.map (fun (m, c) => (m, (c.val.toNat : Int)))

-- MvPolyZp → SparsePolyZp: 投影到单变量（假设 mono 仅含一个 var 条目）
instance : HasPolyConvert MvPolyZp SparsePolyZp where
  convert (f : MvPolyZp) (_target : SparsePolyZp) : SparsePolyZp :=
    SparsePolyZp.normalization (f.map (fun (m, c) =>
      let deg : Nat := if h : 0 < m.size then m[0].snd.toNatClampNeg else 0
      (⟨deg⟩, c)))

-- SparsePolyZp + var → MvPolyZp: 单变量提升到多变量
instance : HasPolyConvert3 SparsePolyZp MvPolyZp Variable where
  convert3 (f : SparsePolyZp) (_target : MvPolyZp) (var : Variable) : MvPolyZp :=
    f.map (fun (umono, c) =>
      let deg64 : Int64 := (umono.deg : Int64)
      let mono : Monomial := if deg64 = 0 then #[] else #[(var, deg64)]
      (mono, c))

-- SparsePolyZZ + var → MvPolyZZ: 同上对 ZZ
instance : HasPolyConvert3 SparsePolyZZ MvPolyZZ Variable where
  convert3 (f : SparsePolyZZ) (_target : MvPolyZZ) (var : Variable) : MvPolyZZ :=
    f.map (fun (umono, c) =>
      let deg64 : Int64 := (umono.deg : Int64)
      let mono : Monomial := if deg64 = 0 then #[] else #[(var, deg64)]
      (mono, c))

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

-- Zp.pow 验证（HPow Zp Int64 → 走 toNatClampNeg → Zp.pow Nat）
#eval ((⟨3, 7⟩ : Zp) ^ (4 : Int64))    -- 期望: 3^4 mod 7 = 4
#eval ((⟨2, 13⟩ : Zp) ^ (10 : Int64))  -- 期望: 2^10 mod 13 = 10
#eval ((⟨5, 11⟩ : Zp) ^ (0 : Int64))   -- 期望: 1
#eval ((⟨5, 11⟩ : Zp) ^ (1 : Int64))   -- 期望: 5
