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

-- 4-arg Bezout EEA 形式：(a, b, s_out, t_out) → (gcd, s, t)
-- 留 stub（真实 EEA for SparsePolyZp 是 phase 2A.10）
-- Pass 2b refret transform 把 ref-out 收成 tuple 返回，Pass 2b 同步 rename
-- callee → polynomial_GCD_eea（避免与 2-arg 版本签名冲突）
def polynomial_GCD_eea [Inhabited α] (_a _b _s _t : α) : α × α × α :=
  (default, default, default)

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

-- 长除法主循环：从 r 中持续减去 g 的倍数，累加商到 q
-- partial def 因 Lean 终止性证明依赖 deg(r') < deg(r) 的严格递减（数学正确）
partial def divmodAux (g : SparsePolyZp) (dg : Nat) (lc_g_inv : Zp)
    (q r : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  if r.isEmpty then (q, r)
  else
    let dr := r[0]!.fst.deg
    if dr < dg then (q, r)
    else
      let coeff := r[0]!.snd * lc_g_inv
      let d := dr - dg
      let term : SparsePolyZp := #[(⟨d⟩, coeff)]
      let r' := r - (term * g)
      let q' := q.push (⟨d⟩, coeff)
      divmodAux g dg lc_g_inv q' r'

-- 多项式长除法：f = q * g + r, deg(r) < deg(g)
-- 退化：g 为空（除以 0）→ 返回 (#[], f) 占位
def divmod (f g : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  if g.isEmpty then (#[], f)
  else
    let dg := g[0]!.fst.deg
    let lc_g_inv := g[0]!.snd.inv
    divmodAux g dg lc_g_inv #[] f

-- 欧几里得 GCD：gcd(f, g) = gcd(g, f mod g)
-- partial def 因终止性依赖 deg(f mod g) < deg(g) 严格递减
partial def gcd (f g : SparsePolyZp) : SparsePolyZp :=
  if g.isEmpty then f
  else gcd g (divmod f g).snd

end SparsePolyZp

-- SparsePolyZp 特化 instance（高优先级，覆盖兜底 default）
instance : HasPolyGCD SparsePolyZp where
  polyGCD := SparsePolyZp.gcd

instance : HasPolyDivmod SparsePolyZp where
  polyDivmod := SparsePolyZp.divmod

-- #eval 数值验证（小例 over F_5）
-- (x^2 - 1) / (x - 1) = x + 1, remainder 0
-- (1 - 1 = 0; x^2 ≡ x^2 mod 5)
#eval SparsePolyZp.divmod
  (#[(⟨2⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)  -- x^2 - 1
  (#[(⟨1⟩, Zp.ofInt 1 5), (⟨0⟩, Zp.ofInt (-1) 5)] : SparsePolyZp)  -- x - 1
-- 期望: (#[(1, 1), (0, 1)], #[])  — q = x+1, r = 0

-- gcd(x^2 - 1, x - 1) = x - 1（因为 x-1 整除）
-- 实际返回的可能是 normalized 形式，待 normalization
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

def Variable.mk {α : Type} (_ : α) : Variable := 0
def UniformIntDist.mk (_lo _hi : Int32) : UniformIntDist := 0
def Rng.default : Rng := 42
-- 阶段 G+：Rng = UInt64 abbrev，Pass 5/8 偶尔在某些上下文 emit `.toUInt64`
-- 让 abbrev Rng → UInt64 自动通过；定义 identity 避免 invalid field 错
def Rng.toUInt64 (r : Rng) : UInt64 := r
def UInt64.toUInt64 (u : UInt64) : UInt64 := u

-- Iterator: 占位类型（C++ STL iterator 的 Lean 抽象）
-- Pass 1 把 std::map iterator 等都映射到 Iterator（无具体 elem 类型）
-- Lean 端用 Unit 占位，配合 BEq 实例支持 it == m.end() 比较
abbrev Iterator := Unit
def Iterator.fromList {α : Type} (_a : Array α) : Iterator := ()
def MvMonomial.normalization (m : MvMonomial) : MvMonomial := m
def gcd (a b : Int) : Int := Int.gcd a b
-- polynomial_mod(f : SparsePolyZZ, p : UInt64) : SparsePolyZp
-- 系数 mod p 把 ZZ 多项式变成 Zp 多项式
-- polynomial_mod: SparsePolyZZ + p → SparsePolyZp（实现移到 abbrev SparsePolyZZ 之后）
-- 见下方 §5c 末尾
def next_prime_64 (p : UInt64) : UInt64 := p + 1
def prev_prime_64 (p : UInt64) : UInt64 := if p > 0 then p - 1 else 0
-- leadcoeff: 1-arg / 2-arg overload (Pass 5 emit 都用同一名)
-- 1-arg `leadcoeff p` 返回 ZZ；2-arg `leadcoeff p var` 返回 Poly
-- Lean 端：2-arg 版本（多变量主用），1-arg 用 leadcoeff1 区分
def leadcoeff {α : Type} [Inhabited α] (_p : α) (_var : Variable) : α := default
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
-- StdMap.find / .end 返回 Iterator（C++ side iterator 语义）
-- Lean 端 Iterator = Unit，find 与 end 比较等价于"成员存在"判断
def StdMap.find {κ ν : Type} [BEq κ] [Inhabited ν] (_m : StdMap κ ν) (_k : κ) : Iterator := ()
def StdMap.end {κ ν : Type} (_ : StdMap κ ν) : Iterator := ()
def rd {α : Type} [Inhabited α] (_ : α) : α := default

-- 阶段 G-E：补 corpus 还需要的 stub 占位
def MvPolyZp.size_u64 (f : MvPolyZp) : UInt64 := (Array.size f).toUInt64
def SparsePolyZZ.normalization (f : SparsePolyZZ) : SparsePolyZZ := f
-- Array.range_init: 多 arity overload，C++ 写法 `iota(arr.begin(), arr.end(), start)`
-- Pass 5 emit 通常是 (arr, start) 2-arg；arr 决定大小，start 是初值
def Array.range_init {α : Type} (a : Array α) (_start : Int32) : Array Int32 :=
  (Array.range a.size).map (·.toUInt32.toInt32)

-- 这些是 C++ 局部 lambda（Pass 3 lift 后理论上以 _lambda_<host>_<n>_ir 形式
-- 出现，但 corpus 中存在裸名引用，疑似 Pass 3 漏 lift 或别名重置失败）。
-- 占位让 Lean 通过类型检查；正确语义留 stage G-A 整改。
def compute_theta {α : Type} [Inhabited α] : α := default
def upzp_coeff {α : Type} [Inhabited α] : α := default
def next_p : UInt64 := 2
-- cont(poly) → ZZ: 多项式整数系数的 content (gcd)
-- ContPP typeclass：cont(poly) = gcd 系数（带符号），pp(poly) = poly / cont(poly)
-- 注意：SparsePolyZZ.{cont,pp}Impl 实现移到 abbrev SparsePolyZZ 之后（见下方）
class HasCont (α : Type) where
  cont : α → ZZ

class HasPP (α : Type) where
  pp : α → α

-- 全局调度入口（与 cpp2lean Pass 5 emit 一致）
def cont {α : Type} [HasCont α] (p : α) : ZZ := HasCont.cont p
def pp {α : Type} [HasPP α] (p : α) : α := HasPP.pp p

-- MvPolyZZ / MvPolyZp 暂保 stub（多变量留 phase 2B）
instance : HasCont MvPolyZZ where cont _ := 0
instance : HasPP MvPolyZZ where pp f := f
instance : HasCont MvPolyZp where cont _ := 0
instance : HasPP MvPolyZp where pp f := f
instance : HasCont SparsePolyZp where cont _ := 0  -- Zp 是域，cont 总是 1 / unit；占位
instance : HasPP SparsePolyZp where pp f := f
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
-- 阶段 G+ 修复：C++ `{g, 1}` 嵌入 tuple，1 是 int 但目标 uint64
-- Lean tuple 元素 Coe 不自动传播，定义具体 pair Coe
instance : Coe (MvPolyZZ × Int32) (MvPolyZZ × UInt64) where
  coe p := (p.fst, p.snd.toInt64.toUInt64)
-- 同样为 Array (MvPolyZZ × Int32) → Array (MvPolyZZ × UInt64)
instance : Coe (Array (MvPolyZZ × Int32)) (Array (MvPolyZZ × UInt64)) where
  coe a := a.map (fun (m, i) => (m, i.toInt64.toUInt64))
-- Pass 1 在 Zp context 误把 `Poly` 映射到 MvPolyZZ；某些 site 实际期望 MvPolyZp
-- 提供 lossy Coe（语义层 B2B 测试时再细化）
instance : Coe MvPolyZZ MvPolyZp where coe _ := #[]
instance : Coe MvPolyZp MvPolyZZ where coe _ := #[]
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

-- §5c 续：SparsePolyZZ.cont / pp 实现（在 abbrev 之后才能用 .map / .foldl 等）
-- cont = gcd 所有系数（符号匹配 leading coeff）
def SparsePolyZZ.contImpl (f : SparsePolyZZ) : ZZ :=
  if f.isEmpty then 0
  else
    let c_nat := f.foldl (fun (acc : Nat) (term : UMonomial × Int) =>
      Nat.gcd acc term.snd.natAbs) 0
    let c_int : Int := c_nat
    if f[0]!.snd < 0 then -c_int else c_int

-- pp = f / cont(f)
def SparsePolyZZ.ppImpl (f : SparsePolyZZ) : SparsePolyZZ :=
  let c : ZZ := SparsePolyZZ.contImpl f
  if c = 0 then f
  else f.map (fun term => (term.fst, term.snd / c))

instance : HasCont SparsePolyZZ where cont := SparsePolyZZ.contImpl
instance : HasPP SparsePolyZZ where pp := SparsePolyZZ.ppImpl

-- polynomial_mod: SparsePolyZZ + p → SparsePolyZp
-- 数学定义：每个系数 mod p（用 Zp.ofInt 折回 [0, p)），剔除 0 系数
def polynomial_mod (f : SparsePolyZZ) (p : UInt64) : SparsePolyZp :=
  f.filterMap (fun term =>
    let zp := Zp.ofInt term.snd p
    if zp.val = 0 then none else some (term.fst, zp))

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

-- cont(-2x² - 4) = -2 (sign matches leading coeff)
#eval SparsePolyZZ.contImpl
  (#[(⟨2⟩, (-2 : Int)), (⟨0⟩, (-4 : Int))] : SparsePolyZZ)
-- 期望: -2

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
