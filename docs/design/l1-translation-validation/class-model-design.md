# 可信基类建模设计

## 1. 问题

翻译器当前对 CLPoly 类的处理是 ad-hoc 的——每个方法/构造函数/运算符都是硬编码特判。没有统一模型，导致：
- 遗漏（某些方法没处理 → 错误翻译而非 sorry）
- 不一致（同一个类的方法在不同地方用不同机制处理）
- 不可扩展（加新类需要改翻译器多处代码）

## 2. 方案

### 2.1 三步走

**步骤 1**：在 Lean 中定义可信基类的**完整模型**（类型 + 所有公共操作）
**步骤 2**：构建 C++ → Lean **映射表**（每个 C++ 操作对应哪个 Lean 函数）
**步骤 3**：翻译器**只查表**。在表中 → 翻译；不在表中 → sorry

### 2.2 可信基类清单

CLPoly 因式分解模块依赖的类：

| C++ 类 | 用途 | 头文件 |
|-------|------|--------|
| `Zp` | Z/pZ 系数 | `number.hh` |
| `ZZ` | 大整数 | `number.hh` |
| `umonomial` | 单变量单项式 | `monomial.hh` |
| `upolynomial_<Zp>` | Z/pZ 上稀疏多项式 | `upolynomial.hh` |
| `std::vector<T>` | 动态数组 | 标准库 |
| `std::pair<A,B>` | 二元组 | 标准库 |

## 3. Lean 类模型

### 3.1 Zp

```lean
structure Zp where
  val : UInt64
  prime : UInt64
deriving Repr, Inhabited, BEq

namespace Zp
-- 构造
def ofInt (v : Int) (p : UInt64) : Zp := ⟨(v.toNat % p.toNat).toUInt64, p⟩
def ofUInt64 (v p : UInt64) : Zp := ⟨v % p, p⟩

-- getter（对应 C++ .number(), .prime()）
-- val 和 prime 是结构体字段，自动有 accessor

-- 算术（mod p）
instance : Add Zp where add a b := ⟨(a.val + b.val) % a.prime, a.prime⟩
instance : Sub Zp where sub a b := ⟨(a.val + a.prime - b.val) % a.prime, a.prime⟩
instance : Mul Zp where mul a b := ⟨(a.val * b.val) % a.prime, a.prime⟩
instance : Neg Zp where neg a := ⟨(a.prime - a.val) % a.prime, a.prime⟩

-- 模逆
partial def modInv (a p : UInt64) : UInt64 := sorry -- 扩展欧几里得
def inv (a : Zp) : Zp := ⟨modInv a.val a.prime, a.prime⟩

-- 比较
instance : BEq Zp where beq a b := a.val == b.val && a.prime == b.prime
end Zp
```

### 3.2 UMonomial

```lean
structure UMonomial where
  deg : UInt64
deriving Repr, Inhabited, BEq

namespace UMonomial
def mk (d : UInt64) : UMonomial := ⟨d⟩
-- deg 是结构体字段
end UMonomial
```

### 3.3 SparsePolyZp

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

-- 修改（返回新数组）
def pushBack (f : SparsePolyZp) (t : UMonomial × Zp) : SparsePolyZp := f.push t
def normalization (f : SparsePolyZp) : SparsePolyZp := f -- 移除零项+排序
def reserve (f : SparsePolyZp) (_n : UInt64) : SparsePolyZp := f -- no-op

-- .data() → identity（C++ 返回裸指针，Lean 中无意义）
def data (f : SparsePolyZp) : SparsePolyZp := f

-- .comp() → 比较函数（CLPoly 内部）
def comp (_f : SparsePolyZp) : UInt64 := 0
end SparsePolyZp
```

### 3.4 std::pair

```lean
-- Lean 的 Prod (α × β) 已有 .fst / .snd
-- C++ .first → .fst, .second → .snd
-- std::make_pair(a, b) → (a, b) 或 Prod.mk a b
```

### 3.5 std::vector

```lean
-- Lean 的 Array α 已有所有操作
-- .empty() → .isEmpty
-- .front() → 自定义 Array.front!
-- .push_back(x) → .push x
-- .size() → .size_u64（返回 UInt64 而非 Nat）
-- operator[] → arr[i.toNat]!
```

## 4. C++ → Lean 映射表

翻译器使用此表。**不在表中的操作 → sorry。**

```python
CLASS_MAP = {
    # ============================================================
    # Zp
    # ============================================================
    "Zp": {
        "lean_type": StructType("Zp"),
        "constructors": {
            # (参数类型列表) → Lean 构造函数名
            ("int64_t", "uint64_t"): "Zp.ofInt",
            ("uint64_t", "uint64_t"): "Zp.ofUInt64",
            (): "default",  # Zp() → ⟨0, 0⟩
        },
        "methods": {
            # 方法名 → (类别, Lean 名)
            # 类别: field=字段访问, method=函数调用, mutate=返回新值, noop=忽略
            "number": ("field", "val"),
            "prime": ("field", "prime"),
        },
        "operators": {
            "+": "HAdd.hAdd",
            "-": "HSub.hSub",
            "*": "HMul.hMul",
            "/": "Zp.div",
            "==": "BEq.beq",
            "!=": "bne",
        },
    },

    # ============================================================
    # umonomial
    # ============================================================
    "umonomial": {
        "lean_type": StructType("UMonomial"),
        "constructors": {
            ("uint64_t",): "UMonomial.mk",
            (): "default",
        },
        "methods": {
            "deg": ("field", "deg"),
        },
    },

    # ============================================================
    # upolynomial_<Zp> (= SparsePolyZp)
    # ============================================================
    "upolynomial_<Zp>": {
        "lean_type": StructType("SparsePolyZp"),
        "constructors": {
            (): "SparsePolyZp.empty",  # 默认构造 → #[]
        },
        "methods": {
            "empty": ("method", "SparsePolyZp.isEmpty"),
            "front": ("method", "SparsePolyZp.front!"),
            "back": ("method", "SparsePolyZp.back!"),
            "size": ("method", "SparsePolyZp.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "normalization": ("mutate", "SparsePolyZp.normalization"),
            "reserve": ("noop", None),
            "data": ("identity", None),
            "comp": ("method", "SparsePolyZp.comp"),
        },
    },

    # ============================================================
    # std::vector<T>
    # ============================================================
    "std::vector": {
        "methods": {
            "empty": ("method", "Array.isEmpty"),
            "front": ("method", "Array.front!"),
            "back": ("method", "Array.back!"),
            "size": ("method", "Array.size_u64"),
            "push_back": ("mutate_push", "Array.push"),
            "reserve": ("noop", None),
        },
    },

    # ============================================================
    # std::pair<A,B>
    # ============================================================
    "std::pair": {
        "constructors": {
            2: "Prod.mk",  # make_pair(a, b)
        },
        "fields": {
            "first": "fst",
            "second": "snd",
        },
    },
}
```

## 5. 翻译器改造

### 5.1 parse_expr 改造

当前：每种方法都是 if/elif 特判。

改后：查 `CLASS_MAP`。

```python
def parse_member_call(node):
    obj_type = get_type(obj_node)  # 从 Clang AST 获取 this 对象类型
    method_name = member.get("name")

    class_info = CLASS_MAP.get(obj_type)
    if class_info is None:
        return sorry("unknown class")

    method_info = class_info["methods"].get(method_name)
    if method_info is None:
        return sorry(f"unknown method {obj_type}.{method_name}")

    category, lean_name = method_info
    if category == "field":
        return FieldAccess(obj, lean_name)
    elif category == "method":
        return Call(lean_name, [obj])
    elif category == "mutate":
        ...  # 返回新 SSA 版本
    elif category == "noop":
        return obj  # 透传
    elif category == "identity":
        return obj
```

### 5.2 CXXConstructExpr 改造

```python
def parse_construct(node):
    target_type = get_type(node)
    class_info = CLASS_MAP.get(target_type)
    if class_info is None:
        return sorry("unknown constructor")

    arg_types = tuple(get_type(arg) for arg in args)
    ctor = class_info["constructors"].get(arg_types)
    if ctor is None:
        return sorry(f"unknown constructor {target_type}({arg_types})")

    return Call(ctor, args)
```

### 5.3 未知操作 → sorry

任何不在 `CLASS_MAP` 中的操作：

```python
return UnknownExpr(f"unmapped: {class_name}.{method_name}")
# gen_expr 输出: /- unmapped: Zp.some_method -/ sorry
```

这比当前的行为（静默错误翻译）安全得多。

## 6. 实施步骤

1. **Lean 类模型文件** `tools/cpp2lean/clpoly_model.lean`：所有类型定义 + 方法实现
2. **映射表** `tools/cpp2lean/class_map.py`：CLASS_MAP 字典
3. **翻译器改造**：`clang_ast.py` 的 `CXXMemberCallExpr`、`CXXConstructExpr`、`CXXOperatorCallExpr`、`MemberExpr` 全部改为查表
4. **删除硬编码**：`CLPOLY_TYPE_MAP`、`METHOD_MAP`、所有 if/elif 特判
5. **验证**：重新翻译 13 函数，对比结果

## 7. 预期效果

- 翻译器代码大幅简化（删除 ~100 行特判代码）
- 新增类只需扩展 CLASS_MAP（不改翻译器）
- 未知操作显式 sorry（而非静默错误）
- Lean 类模型可独立编译验证
- 背靠背测试中，类模型就是原语的正确实现
