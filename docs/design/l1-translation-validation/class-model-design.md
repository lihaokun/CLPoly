# 可信基类建模设计

## 1. 问题

翻译器当前对 CLPoly 类的处理是 ad-hoc 的——每个方法/构造函数/运算符都是硬编码特判。没有统一模型，导致：
- 遗漏（某些方法没处理 → 错误翻译而非 sorry）
- 不一致（同一个类的方法在不同地方用不同机制处理）
- 不可扩展（加新类需要改翻译器多处代码）

## 2. 方案

### 2.1 三步走

**步骤 1**：在 Lean 中定义可信基类的**完整模型**（类型 + 所有公共操作）
**步骤 2**：构建 C++ → Lean **映射表**（类操作 + 独立函数，每个 C++ 操作对应哪个 Lean 函数）
**步骤 3**：翻译器**只查表**。在表中 → 翻译；不在表中 → sorry

### 2.2 可信基清单

CLPoly 因式分解模块依赖的：

**类**：

| C++ 类 | 用途 |
|-------|------|
| `Zp` | Z/pZ 系数 |
| `ZZ` | 大整数（仅 `__upoly_powmod` 和 `__edf_Zp` 使用） |
| `umonomial` | 单变量单项式 |
| `upolynomial_<Zp>` | Z/pZ 上稀疏多项式 |
| `std::vector<T>` | 动态数组 |
| `std::pair<A,B>` | 二元组 |

**独立函数**（非类方法）：

| C++ 函数 | 用途 |
|---------|------|
| `derivative(f)` | 多项式求导 |
| `polynomial_GCD(f, g)` | 多项式 GCD |
| `pair_vec_div(q, r, f, g, comp)` | 多项式除法（5 参数，前 2 个是输出） |
| `get_deg(f)` | 取多项式首项度数 |
| `std::move(x)` | 移动语义（纯函数式中 = identity） |
| `std::make_pair(a, b)` | 构造 pair |

## 3. Lean 类模型

### 3.1 Zp

```lean
structure Zp where
  val : UInt64
  prime : UInt64
deriving Repr, Inhabited, BEq

namespace Zp
-- 构造函数
-- C++ Zp(int64_t val, uint64_t p)：val 可能为负
def ofInt (v : Int) (p : UInt64) : Zp :=
  let r := v % (p.toNat : Int)
  let r := if r < 0 then r + p.toNat else r
  ⟨r.toNat.toUInt64, p⟩

-- C++ Zp(uint64_t val, uint64_t p)
def ofUInt64 (v p : UInt64) : Zp := ⟨v % p, p⟩

-- getter（.number() → .val, .prime() → .prime）
-- val 和 prime 是结构体字段，自动有 accessor

-- 算术（mod p）
instance : Add Zp where add a b := ⟨(a.val + b.val) % a.prime, a.prime⟩
instance : Sub Zp where sub a b := ⟨(a.val + a.prime - b.val) % a.prime, a.prime⟩
instance : Mul Zp where mul a b := ⟨(a.val * b.val) % a.prime, a.prime⟩
instance : Neg Zp where neg a := ⟨(a.prime - a.val) % a.prime, a.prime⟩

-- 模逆（扩展欧几里得）
partial def modInvAux (a b x0 x1 : Int) : Int :=
  if a <= 1 then x1
  else modInvAux b (a % b) (x1 - (a / b) * x0) x0

def modInv (a p : UInt64) : UInt64 :=
  let r := modInvAux (a.toNat : Int) (p.toNat : Int) 0 1
  let r := if r < 0 then r + p.toNat else r
  r.toNat.toUInt64

def inv (a : Zp) : Zp := ⟨modInv a.val a.prime, a.prime⟩
def div (a b : Zp) : Zp := a * b.inv
end Zp
```

### 3.2 ZZ

```lean
-- CLPoly 的 ZZ 是 GMP 大整数。Lean 的 Int 是任意精度。
-- 因式分解中 ZZ 仅用于 __upoly_powmod 的指数和 __edf_Zp 的幂计算。
abbrev ZZ := Int

namespace ZZ
def ofInt (v : Int) : ZZ := v
def div (a b : ZZ) : ZZ := a / b
end ZZ
```

### 3.3 UMonomial

```lean
structure UMonomial where
  deg : UInt64
deriving Repr, Inhabited, BEq

namespace UMonomial
def mk (d : UInt64) : UMonomial := ⟨d⟩
end UMonomial
```

### 3.4 SparsePolyZp

```lean
abbrev SparsePolyZp := Array (UMonomial × Zp)

namespace SparsePolyZp
-- 构造
def empty : SparsePolyZp := #[]

-- 查询
def isEmpty (f : SparsePolyZp) : Bool := f.isEmpty
def front! (f : SparsePolyZp) : UMonomial × Zp := f[0]!
def back! (f : SparsePolyZp) : UMonomial × Zp := f[f.size - 1]!
def getDeg (f : SparsePolyZp) : UInt64 := if f.isEmpty then 0 else f[0]!.fst.deg
def size_u64 (f : SparsePolyZp) : UInt64 := f.size.toUInt64

-- 修改（返回新数组——纯函数式）
def pushBack (f : SparsePolyZp) (t : UMonomial × Zp) : SparsePolyZp := f.push t

-- normalization：移除零系数项（C++ 原地修改，Lean 返回新数组）
def normalization (f : SparsePolyZp) : SparsePolyZp :=
  f.filter (fun (_, c) => c.val != 0)

-- no-op 操作
def reserve (f : SparsePolyZp) (_n : UInt64) : SparsePolyZp := f
def data (f : SparsePolyZp) : SparsePolyZp := f  -- .data() → identity

-- comp：比较函数（CLPoly 内部使用，翻译中不传递）
-- C++ 中 f.comp() 返回比较函数对象，在 pair_vec_div 中使用
-- Lean 中不需要——pair_vec_div 的 Lean 实现不使用比较函数
end SparsePolyZp
```

### 3.5 std::pair

```lean
-- Lean Prod (α × β) 已有 .fst / .snd
-- C++ .first → .fst, .second → .snd
-- std::make_pair(a, b) → Prod.mk a b
```

### 3.6 std::vector

```lean
-- Lean Array α 已有所有操作
-- C++ 操作通过 CLASS_MAP 映射到 Lean Array 方法
```

## 4. 独立函数模型（FUNC_MAP）

CLASS_MAP 覆盖类方法。独立函数用 **FUNC_MAP** 覆盖。

```python
FUNC_MAP = {
    # C++ 函数名 → (Lean 函数名, 参数处理规则)
    # 参数处理：
    #   "direct" = 直接传参
    #   "out_param(n)" = 第 n 个参数是输出参数，变为返回值
    #   "drop(n)" = 丢弃第 n 个参数（如 comp）

    "derivative": ("SparsePolyZp.derivative", "direct"),
    "polynomial_GCD": ("SparsePolyZp.gcd", "direct"),
    "pair_vec_div": ("SparsePolyZp.divmod", "out_param(0,1),drop(4)"),
    # C++: pair_vec_div(q, r, f, g, comp) → Lean: SparsePolyZp.divmod f g
    # 返回 (q, r)，丢弃 comp 参数

    "get_deg": ("SparsePolyZp.getDeg", "direct"),
    "std::move": ("identity", "pass_through"),  # move(x) → x
    "std::make_pair": ("Prod.mk", "direct"),
}
```

### 4.1 pair_vec_div 的特殊处理

C++ 签名：`pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp())`
- 参数 0, 1 是输出（q, r）
- 参数 2, 3 是输入（f, g）
- 参数 4 是比较函数（Lean 中不需要）

Lean 翻译：
```lean
let (q, r) := SparsePolyZp.divmod f g
```

翻译器规则：遇到 `pair_vec_div` → 提取参数 2,3 作为输入 → 生成 `let (out0, out1) := divmod arg2 arg3` → SSA 更新 out0, out1 的版本。

## 5. C++ → Lean 映射表（CLASS_MAP）

翻译器使用此表。**不在表中的操作 → sorry。**

映射表的 key 使用 `TypeIR`（不是字符串），与翻译器内部类型一致。

Clang 的 `type.qualType` 可能包含命名空间和 const/引用（如 `"const clpoly::Zp &"`）。`parse_type` 已经负责清理这些 → 得到 `TypeIR`。映射表基于清理后的 `TypeIR`。

```python
from ir_types import BaseType, StructType

CLASS_MAP = {
    # key = StructType.name（parse_type 的输出）

    "Zp": {
        "lean_type": StructType("Zp", []),
        "constructors": {
            # key = 参数类型 tuple
            (BaseType.INT64, BaseType.UINT64): "Zp.ofInt",
            (BaseType.UINT64, BaseType.UINT64): "Zp.ofUInt64",
            (): "default",
        },
        "methods": {
            "number": ("field", "val"),       # .number() → .val
            "prime": ("field", "prime"),       # .prime() → .prime
            "val": ("field", "val"),           # .val 直接访问
        },
        "operators": {
            "+": None,   # Lean Add instance 处理
            "-": None,   # Lean Sub instance 处理
            "*": None,   # Lean Mul instance 处理
            "/": "Zp.div",
            "==": None,  # Lean BEq instance 处理
            "!=": None,
        },
    },

    "UMonomial": {
        "lean_type": StructType("UMonomial", []),
        "constructors": {
            (BaseType.UINT64,): "UMonomial.mk",
            (): "default",
        },
        "methods": {
            "deg": ("field", "deg"),
        },
    },

    "SparsePolyZp": {
        "lean_type": StructType("SparsePolyZp", []),
        "constructors": {
            (): "SparsePolyZp.empty",
        },
        "methods": {
            # (类别, Lean 名)
            # field: 字段访问 → obj.lean_name
            # method: 函数调用 → lean_name obj
            # mutate: 返回新值 → let obj' := lean_name obj（SSA 新版本）
            # mutate_push: push → let obj' := obj.push elem
            # noop: 忽略 → obj
            # identity: 透传 → obj
            "empty": ("method", "SparsePolyZp.isEmpty"),
            "front": ("method", "SparsePolyZp.front!"),
            "back": ("method", "SparsePolyZp.back!"),
            "size": ("method", "SparsePolyZp.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "SparsePolyZp.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
        },
    },
}
```

### 5.1 类别说明

| 类别 | 含义 | Lean 生成 | SSA 处理 |
|------|------|----------|---------|
| `field` | 字段访问 | `obj.field_name` | 无（纯读取） |
| `method` | 无副作用方法调用 | `lean_name obj` | 无 |
| `mutate` | 修改对象的方法 | `lean_name obj` | 生成 SSA 新版本 `let obj' := lean_name obj` |
| `mutate_push` | push_back | `obj.push elem` | 生成 SSA 新版本 |
| `noop` | 无操作（如 reserve） | `obj`（透传） | 无 |
| `identity` | 返回自身（如 .data()） | `obj` | 无 |

`mutate` 类别的处理在 `ssa_transform.py` 中：当翻译器通过 CLASS_MAP 识别出 `mutate` 方法时，自动生成 SSA 新版本。这替代了当前 ssa_transform 中的硬编码特判（`_mutate_normalize`、`__upoly_make_monic` 等）。

## 6. 未知操作处理

翻译器遇到的每个操作，按以下顺序查找：

1. **CLASS_MAP**：`obj_type.name → methods[method_name]`
2. **FUNC_MAP**：`func_name → (lean_name, param_rule)`
3. **LEAN_STDLIB**：`func_name → lean_name`（Prod.mk, Array.empty 等）
4. **CAST_TABLE**：`(source, target) → lean_expr`
5. **都不匹配 → `sorry`** + 注释标注 C++ 原始操作

## 7. 实施步骤

1. **`clpoly_model.lean`**：Lean 类模型文件（§3 的全部定义）→ 先编译验证
2. **`class_map.py`**：CLASS_MAP + FUNC_MAP 字典
3. **改造 `clang_ast.py`**：`CXXMemberCallExpr`、`CXXConstructExpr`、`CXXOperatorCallExpr`、`MemberExpr` → 查 CLASS_MAP
4. **改造 `ssa_transform.py`**：`mutate` 类别统一处理 → 删除 `_mutate_normalize` 等特判
5. **改造 `gen_full.py`**：prelude 用 `clpoly_model.lean` 替代手写 header
6. **删除硬编码**：`CLPOLY_TYPE_MAP`、`METHOD_MAP`、`LEAN_STDLIB`、`NOOP_FUNCS` → 全部合并入 CLASS_MAP/FUNC_MAP
7. **验证**：重新翻译 13 函数 → 编译通过 → 背靠背测试

## 8. 前提条件

| 前提 | 验证方式 |
|------|---------|
| `clpoly_model.lean` 编译通过 | `lake env lean clpoly_model.lean` |
| CLASS_MAP 覆盖因式分解模块用到的所有类操作 | 对 13 函数翻译，sorry 数 = 0 |
| FUNC_MAP 覆盖所有独立函数 | 同上 |
| Zp.ofInt 对负数正确 | 单元测试 `Zp.ofInt (-1) 5 = ⟨4, 5⟩` |
| SparsePolyZp.normalization 行为与 C++ 一致 | 背靠背测试 |
