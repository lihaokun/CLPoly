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
-- 阶段 G+ 修复：Pass 5 emit `Zp != (0 : Int32)` / `assign x (- alpha_j)` 等
-- 字面量保持 Int32 而非 Zp。提供 Coe 让 Lean 自动桥接（占位 prime=1）。
instance : Coe Int32 Zp where coe i := { val := i.toInt64.toUInt64, prime := 1 }
instance : Coe Int Zp where coe i := { val := i.toNat.toUInt64, prime := 1 }

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
-- C++ `.comp()` 返回比较器（lex_<less> 类对象）；Lean 用 Unit 占位
-- (Lex = Unit abbrev 见 §5d，此处先用 Unit 避免 forward ref)
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
-- 阶段 E：Monomial / MvMonomial / Variable 类型统一
-- C++ side: `Monomial = basic_monomial<lex_<less>>` = vector<pair<variable, exponent>>
-- variable 是 UInt64 ID（CLPoly clpoly::variable handle），exponent 是 Int64
abbrev Variable := UInt64
abbrev Monomial := Array (Variable × Int64)
-- MvMonomial 与 Monomial 同义（Pass 1 不同 context 偶尔产生 MvMonomial 别名）
abbrev MvMonomial := Monomial
def Monomial.empty : Monomial := #[]
-- Pass 5 把 C++ 多种 1-arg Monomial 构造（comp_ptr / element 等）都映射到 .mk；
-- 用泛型 input 统一占位（Variable×Int64 输入仍可隐式 Coe）
def Monomial.mk {α : Type} (_ : α) : Monomial := #[]
def MvMonomial.empty : MvMonomial := #[]
-- 多变量多项式：内部元素的第一槽是 Monomial，与 Pass 1 推断的 (Monomial × ZZ)
-- 一致；ZZ = Int / Zp 同样为系数类型
abbrev MvPolyZZ := Array (Monomial × Int)
abbrev MvPolyZp := Array (Monomial × Zp)
abbrev PolyZp := MvPolyZp
abbrev PolyZZ := MvPolyZZ
abbrev PolyQQ := Array (Monomial × Rat)

-- MvPolyZp 操作（stub；实际语义留 Pass 上游 / B2B 测试细化）
def MvPolyZp.normalization (f : MvPolyZp) : MvPolyZp := f
def MvPolyZZ.empty : MvPolyZZ := #[]
def MvPolyZp.empty : MvPolyZp := #[]
-- C++ 的 Poly(comp_t) / Poly(const Poly&) ctor 都映射到 .mk；
-- 用泛型 input → 空 Poly 占位（语义在 B2B 层细化）
def MvPolyZp.mk {α : Type} (_ : α) : MvPolyZp := #[]
def MvPolyZZ.normalization (f : MvPolyZZ) : MvPolyZZ := f
def MvPolyZZ.mk {α : Type} (_ : α) : MvPolyZZ := #[]
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
-- Lex: 单项式序 tag —— 仅类型层占位
abbrev Lex := Unit
def polynomial_GCD [Inhabited α] (_a _b : α) : α := default
-- 4-arg Bezout EEA 形式：(a, b, s_out, t_out) → (gcd, s, t)
-- Pass 2b refret transform 把 ref-out 收成 tuple 返回，Pass 2b 同步 rename
-- callee → polynomial_GCD_eea（避免与 2-arg 版本签名冲突）
def polynomial_GCD_eea [Inhabited α] (_a _b _s _t : α) : α × α × α :=
  (default, default, default)
-- pair_vec_div: 4 参数版本（C++ side: pair_vec_div(f, g, q, comp) → 返回 quotient）
-- 占位实现，B2B 测试细化（comp 通常是比较器/函数对象，签名宽松）
def pair_vec_div [Inhabited α] (_f _g _q : α) (_comp : β) : α := default
-- 5-arg overload: (new_v, R, v1, v2, comp) — basic.hh:698 形态
-- Pass 2b refret 把 R 收成 tuple → return (q, R)
def pair_vec_div5 [Inhabited α] (_f _g _q _r : α) (_comp : β) : α × α := (default, default)

-- Array.insert: C++ STL set::insert / vec.push_back 的占位（push 到末尾）
def Array.insert (a : Array α) (v : α) : Array α := a.push v

-- Array.findVal: C++ std::find(vec, x) 的语义（按值找）
-- Lean 4 Array.find? 期望 predicate；这里 wrap 成"按 BEq 找值"
def Array.findVal [BEq α] (a : Array α) (x : α) : Option α := a.find? (· == x)

-- comp 方法占位（已存在于 namespace SparsePolyZp 之内为 UInt64）；
-- 补 SparsePolyZZ / MvPolyZp / MvPolyZZ
-- C++ `Poly.comp()` 返回比较器对象（lex_<less>）而非 UInt64
def SparsePolyZZ.comp (_f : SparsePolyZZ) : Lex := ()
def MvPolyZp.comp (_f : MvPolyZp) : Lex := ()
def MvPolyZZ.comp (_f : MvPolyZZ) : Lex := ()

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

-- Array.sort: 占位（C++ std::sort with comparator）
def Array.sort {α : Type} (a : Array α) (_cmp : α → α → Bool) : Array α := a

-- C++ 自由函数 degree(poly) — 多态（lambda 比较器里用）。
-- 用 typeclass 解决"未限定 degree"调用的多重 receiver 类型问题。
class HasDegree (α : Type) where
  degree : α → UInt64

instance : HasDegree MvPolyZZ where degree _ := 0
instance : HasDegree MvPolyZp where degree _ := 0
-- SparsePolyZZ HasDegree instance 在 abbrev 后再加（见 §5c）

def degree {α : Type} [HasDegree α] (a : α) : UInt64 := HasDegree.degree a
-- 2-arg degree: poly + main var → 关于该 var 的次数（C++ degree(poly, var)）
def degree2 {α : Type} [HasDegree α] (a : α) (_var : Variable) : Int64 :=
  (HasDegree.degree a).toInt64

-- C++ free functions: is_number(poly) / get_variables(poly)
class IsNumber (α : Type) where
  isNumber : α → Bool

instance : IsNumber MvPolyZZ where isNumber f := (f : Array _).isEmpty
instance : IsNumber MvPolyZp where isNumber f := (f : Array _).isEmpty
-- SparsePolyZZ instance 在其 abbrev 之后再定义（见 §5c）

def is_number {α : Type} [IsNumber α] (a : α) : Bool := IsNumber.isNumber a

class GetVariables (α : Type) where
  vars : α → Array (Variable × Int64)

instance : GetVariables MvPolyZZ where vars _ := #[]
instance : GetVariables MvPolyZp where vars _ := #[]

def get_variables {α : Type} [GetVariables α] (a : α) : Array (Variable × Int64) :=
  GetVariables.vars a

-- MvPolyZZ.front! / MvPolyZp.front!：取首项（mono × coeff）
def MvPolyZZ.front! (f : MvPolyZZ) : (Monomial × Int) := f[0]!
def MvPolyZp.front! (f : MvPolyZp) : (Monomial × Zp) := f[0]!

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
instance : HMul MvPolyZZ MvPolyZZ MvPolyZZ where hMul a b := a ++ b
instance : HAdd MvPolyZZ MvPolyZZ MvPolyZZ where hAdd a b := a ++ b
instance : HSub MvPolyZZ MvPolyZZ MvPolyZZ where hSub a b := a ++ b
instance : HMul MvPolyZp MvPolyZp MvPolyZp where hMul a b := a ++ b
instance : HAdd MvPolyZp MvPolyZp MvPolyZp where hAdd a b := a ++ b
instance : HSub MvPolyZp MvPolyZp MvPolyZp where hSub a b := a ++ b

-- derivative typeclass：C++ free function `derivative(poly)` 对所有 poly 类型多态
class HasDerivative (α : Type) where
  derivative : α → α

instance : HasDerivative SparsePolyZp where derivative f := f
instance : HasDerivative SparsePolyZZ where derivative f := f
instance : HasDerivative MvPolyZZ where derivative f := f
instance : HasDerivative MvPolyZp where derivative f := f

def derivative {α : Type} [HasDerivative α] (a : α) : α := HasDerivative.derivative a

-- squarefreefactorize 占位（多变量 ZZ 默认；其他实例需要时再加）
def squarefreefactorize (f : MvPolyZZ) : Array (MvPolyZZ × UInt64) := #[(f, 1)]

-- poly_convert: 跨域多项式系数转换占位（C++ 模板函数）
-- 2-arg 版本（C++ side `poly_convert(p, target)`）
def poly_convert {α β : Type} (_f : α) (target : β) : β := target
-- 3-arg 版本（C++ side `poly_convert(p, target, ctx)`，ctx 是某 lex/var 标识）
def poly_convert3 {α β γ : Type} (_f : α) (target : β) (_ctx : γ) : β := target

-- SparsePolyZZ 的 OfNat 0 实例：见 §5c（abbrev 定义之后）

-- C++ 全局常量 / 宏：占位（B2B 时填实际值）
def ZASSENHAUS_THRESHOLD : Int32 := 8
def __g_use_large_prime : Bool := false

-- ZZ.invert: 模逆元（C++ mpz_invert(result, op, mod) → 0/1 success）
-- 3 参数版：(out_dummy, num, mod) → Bool
def ZZ.invert (_out _num _mod : ZZ) : Bool := true

-- ZZ.fdiv_q / ZZ.fdiv_r: 向下取整除法（3 参数版：result, dividend, divisor → unit）
def ZZ.fdiv_q (_out a b : ZZ) : ZZ := a / b
def ZZ.fdiv_r (_out a b : ZZ) : ZZ := a % b

-- ZZ = Int alias 时 `(x : ZZ).toInt` 不合法（Int 没 .toInt）。
-- Pass 5 cast_table 在某些 ZZ → Int 路径仍 emit `.toInt`；提供 identity 兜底。
def Int.toInt (x : Int) : Int := x

-- C++ 数学/IO 内置占位（B2B 测试时细化语义）
-- C++ log(x : double) : double — 提供同名 Lean 占位（屏蔽 Mathlib Nat.log 名冲突）
namespace Nat
def log (_x : Float) : Float := 1.0
end Nat
def Int.toFloat (n : Int) : Float := Float.ofInt n
def ZZ.sizeinbase (_z : ZZ) (_base : Int32) : UInt64 := 0
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

def Variable.mk {α : Type} (_ : α) : Variable := 0
def UniformIntDist.mk (_lo _hi : Int32) : UniformIntDist := 0
def Rng.default : Rng := 42
-- 阶段 G+：Rng = UInt64 abbrev，Pass 5/8 偶尔在某些上下文 emit `.toUInt64`
-- 让 abbrev Rng → UInt64 自动通过；定义 identity 避免 invalid field 错
def Rng.toUInt64 (r : Rng) : UInt64 := r
def UInt64.toUInt64 (u : UInt64) : UInt64 := u

def Iterator {α : Type} (a : Array α) : Array α := a
def MvMonomial.normalization (m : MvMonomial) : MvMonomial := m
def gcd (a b : Int) : Int := Int.gcd a b
def polynomial_mod {α : Type} [Inhabited α] (_a _b : α) : α := default
def next_prime_64 (p : UInt64) : UInt64 := p + 1
def prev_prime_64 (p : UInt64) : UInt64 := if p > 0 then p - 1 else 0
-- leadcoeff: 1-arg / 2-arg overload (Pass 5 emit 都用同一名)
-- 1-arg `leadcoeff p` 返回 ZZ；2-arg `leadcoeff p var` 返回 Poly
-- Lean 端：2-arg 版本（多变量主用），1-arg 用 leadcoeff1 区分
def leadcoeff {α : Type} [Inhabited α] (_p : α) (_var : Variable) : α := default
def leadcoeff1 {α : Type} [Inhabited α] (_p : α) : ZZ := 0
def ZZ.fdiv_ui (_a : ZZ) (_b : UInt64) : ZZ := 0
def StdMap.find {κ ν : Type} [BEq κ] [Inhabited ν] (m : StdMap κ ν) (k : κ) : ν :=
  StdMap.get! m k
def StdMap.end {κ ν : Type} (_ : StdMap κ ν) : Unit := ()
def rd {α : Type} [Inhabited α] (_ : α) : α := default

-- 阶段 G-E：补 corpus 还需要的 stub 占位
def MvPolyZp.size_u64 (f : MvPolyZp) : UInt64 := (Array.size f).toUInt64
def SparsePolyZZ.normalization (f : SparsePolyZZ) : SparsePolyZZ := f
def Array.range_init (n : UInt64) : Array Int32 :=
  (Array.range n.toNat).map (·.toUInt32.toInt32)

-- 这些是 C++ 局部 lambda（Pass 3 lift 后理论上以 _lambda_<host>_<n>_ir 形式
-- 出现，但 corpus 中存在裸名引用，疑似 Pass 3 漏 lift 或别名重置失败）。
-- 占位让 Lean 通过类型检查；正确语义留 stage G-A 整改。
def compute_theta {α : Type} [Inhabited α] : α := default
def upzp_coeff {α : Type} [Inhabited α] : α := default
def next_p : UInt64 := 2
-- cont(poly) → ZZ: 多项式整数系数的 content (gcd)
def cont {α : Type} (_p : α) : ZZ := 0
-- pp(poly) → poly: primitive part (poly / cont)
def pp {α : Type} [Inhabited α] (_p : α) : α := default
def all_div : Bool := false
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
instance : Coe Int32 Int64 where coe n := n.toInt64
instance : Coe UInt64 Nat where coe n := n.toNat
instance : Coe Int64 Nat where coe n := n.toNatClampNeg
-- 阶段 G+：Pass 5 cast 漏的 site 用 Lean Coe 自动桥（uni-directional safe casts）
instance : Coe Nat Int64 where coe n := n.toUInt64.toInt64
instance : Coe Nat Int32 where coe n := n.toUInt32.toInt32
instance : Coe UInt64 Int64 where coe u := u.toInt64
instance : Coe UInt32 UInt64 where coe u := u.toUInt64
instance : Coe ZZ Nat where coe z := z.toNat
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
instance : HasDegree SparsePolyZZ where degree _ := 0
instance : HasDegree SparsePolyZp where degree _ := 0

-- OfNat 0 实例：C++ `SparsePolyZZ x = 0` → 空多项式
instance : OfNat SparsePolyZZ 0 where ofNat := #[]
instance : OfNat SparsePolyZp 0 where ofNat := #[]
instance : OfNat MvPolyZZ 0 where ofNat := #[]
instance : OfNat MvPolyZp 0 where ofNat := #[]

def SparsePolyZZ.size_u64 (f : SparsePolyZZ) : UInt64 := f.size.toUInt64

-- 阶段 F 后续：依赖 SparsePolyZZ 的 stub（LLLMatrix.size 见 abbrev 之后）
-- get_first_deg: 多变量 / 单变量两态。Lean 端泛型占位（语义层 B2B 细化）
def get_first_deg {α : Type} (_f : α) : Int64 := 0

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
-- C++ assign(poly, var, val) = 用 val 替代 poly 中的变量 var
opaque assign (poly : α) (var : Variable) (val : β) : α := poly
-- 2-arg overload：assign(poly, eval_point) 用 map 一次性代入多个变量
def assign2 {α : Type} (poly : α) (_eval_point : StdMap Variable ZZ) : α := poly

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
