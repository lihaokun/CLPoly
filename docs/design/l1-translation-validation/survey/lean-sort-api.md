# Lean 4 排序 API 综合目录

## 目标映射

CLPoly C++ 的 `std::sort(v.begin(), v.end(), comparator)` 应映射为 Lean 的 **`Array.qsort arr cmp`** 或 **`List.mergeSort l cmp`**。

---

## 1. Lean 核心库 (Init)

### 1.1 Array.qsort

**位置**: `/home/haokun/.elan/toolchains/leanprover--lean4---v4.28.0/src/lean/Init/Data/Array/QSort/Basic.lean`

**完整签名**:
```lean
@[inline] def Array.qsort (as : Array α) (lt : α → α → Bool := by exact (· < ·))
    (lo := 0) (hi := as.size - 1) : Array α
```

**详细分析**:
- **参数**:
  - `as : Array α` — 待排序数组
  - `lt : α → α → Bool` — 自定义比较函数（less-than 风格，返回 `Bool`）
    - 默认：`(· < ·)` （需要 `LT` 实例）
    - **支持自定义比较器**：`fun x y => x.field ≤ y.field`，或自定义谓词
  - `lo, hi : Nat` — 排序范围 `[lo, hi]`（默认整个数组）
  
- **返回值**: `Array α` — 新排序后的数组（返回新副本，非原地）
  
- **比较器格式**: **Bool 式** (less-than)
  - 接受 `α → α → Bool`
  - 返回 `true` 当 `a < b`，`false` 否则
  - **直接对应 C++ `comparator(a, b)` 返回 bool 的风格**
  
- **Computable**: ✅ 是 (可 `#eval`)
  
- **Partial**: ❌ 否 (fully terminating)
  
- **稳定性**: ❌ 否 (不稳定——quicksort 的特性)
  
- **Typeclass 需求**:
  - 如使用默认 `(· < ·)`：需要 `[LT α]` 且 `(· < ·)` 可决定
  - 如自定义 `lt`：**无 typeclass 需求**（完全自由）
  
- **In-place**: ❌ 返回新 Array （Lean 值语义，无原地修改）

**使用示例**:
```lean
-- 用自定义比较器按降序排序
def descending (a b : Nat) : Bool := b < a
#eval Array.qsort #[3, 1, 4, 1, 5] descending
-- [5, 4, 3, 1, 1]

-- 按字符串长度排序
#eval Array.qsort #["a", "bcd", "ab"] (fun s t => s.length < t.length)
-- ["a", "ab", "bcd"]
```

**关键优势**:
- 最直接对应 C++ `std::sort` 的语义
- 比较器是自由的，不依赖 `Ord` 实例
- Computable，可直接执行和验证

---

### 1.2 Array.qsortOrd

**位置**: 同上 (`QSort/Basic.lean`)

**完整签名**:
```lean
def Array.qsortOrd [ord : Ord α] (xs : Array α) : Array α
```

**详细分析**:
- **参数**:
  - `xs : Array α` — 待排序数组
  - `[ord : Ord α]` — **隐式**，需要定义 `Ord` 实例
  
- **返回值**: `Array α` — 排序后的数组

- **比较器格式**: 使用 `Ord.compare : α → α → Ordering`
  - 内部实现：`xs.qsort fun x y => compare x y |>.isLT`
  - `Ordering.lt` → 输出 `true`，其他 → `false`
  - 不直接接受用户自定义比较器，比较器由 `Ord` 实例固定
  
- **Computable**: ✅ 是 (如果 `Ord` 实例 computable)
  
- **Partial**: ❌ 否
  
- **稳定性**: ❌ 否
  
- **Typeclass 需求**: `[Ord α]` (必需)
  
- **in-place**: ❌ 返回新 Array

**使用示例**:
```lean
-- 对 Nat 排序 (使用默认 Ord 实例)
#eval Array.qsortOrd #[3, 1, 4, 1, 5]
-- #[1, 1, 3, 4, 5]

-- 需要为自定义类型定义 Ord 实例
structure Person where
  name : String
  age : Nat
deriving Ord  -- 自动生成按字典序比较

#eval Array.qsortOrd #[⟨"Alice", 30⟩, ⟨"Bob", 25⟩]
```

**与 qsort 的关系**: `qsortOrd` 是 `qsort` 的高级版本，用于有 `Ord` 实例的类型。

---

### 1.3 Array.insertionSort

**位置**: `/home/haokun/.elan/toolchains/leanprover--lean4---v4.28.0/src/lean/Init/Data/Array/InsertionSort.lean`

**完整签名**:
```lean
@[inline] def Array.insertionSort (xs : Array α) (lt : α → α → Bool := by exact (· < ·)) : Array α
```

**详细分析**:
- **参数**:
  - `xs : Array α` — 待排序数组
  - `lt : α → α → Bool` — less-than 比较函数（同 `qsort`）
  - 默认：`(· < ·)`
  
- **返回值**: `Array α` — 排序后的数组
  
- **比较器格式**: Bool 式 (less-than)，与 `qsort` 完全相同
  
- **Computable**: ✅ 是
  
- **Partial**: ❌ 否
  
- **稳定性**: ✅ **是** (insertion sort 天然稳定)
  
- **Typeclass 需求**: 同 `qsort`（如果使用自定义 `lt`，则无需求）
  
- **in-place**: ❌ 返回新 Array
  
- **时间复杂度**: O(n²) 平均和最坏情况（vs qsort O(n log n) 平均）

**使用示例**:
```lean
#eval Array.insertionSort #[3, 1, 4, 1, 5]
-- #[1, 1, 3, 4, 5]

-- 自定义比较器（按降序）
#eval Array.insertionSort #[3, 1, 4] (fun a b => b < a)
-- #[4, 3, 1]
```

**使用场景**: 数组较小或已近似排序时更高效；需要**稳定性**时必须用此。

---

## 2. List 排序 (Mathlib)

### 2.1 List.insertionSort

**位置**: `/home/haokun/projects/CLPoly/proof/lean/.lake/packages/mathlib/Mathlib/Data/List/Sort.lean`

**完整签名**:
```lean
def List.insertionSort (r : α → α → Prop) [DecidableRel r] : List α → List α
```

**详细分析**:
- **参数**:
  - `r : α → α → Prop` — 二元关系（**命题**，不是 `Bool`！）
  - `[DecidableRel r]` — **必须** 有可决定实例
  - 输入 List
  
- **返回值**: 排序后的 List
  
- **比较器格式**: **命题式** (`α → α → Prop`)
  - 例如 `(· ≤ ·)`、`(· < ·)`（可从 `LT`/`LE` 获取）
  - **与 `Array.insertionSort` 的 Bool 比较器不兼容**
  
- **Computable**: ✅ 是（因为有 `DecidableRel`）
  
- **Partial**: ❌ 否
  
- **稳定性**: ✅ **是**
  
- **Typeclass 需求**: `[DecidableRel r]`
  
- **in-place**: N/A (Lists 都是链表)
  
- **复杂度**: O(n²)

**使用示例**:
```lean
#eval List.insertionSort (· ≤ ·) [3, 1, 4, 1, 5]
-- [1, 1, 3, 4, 5]

-- 自定义关系：按降序
#eval List.insertionSort (· ≥ ·) [3, 1, 4, 1, 5]
-- [5, 4, 3, 1, 1]
```

**关键限制**: 
- **不能直接用 Bool 函数** — Mathlib 版本要求 `DecidableRel`
- 需要显式命题 + 可决定证明
- 对 Lean 初学者较陡峭

---

### 2.2 List.mergeSort

**位置**: `/home/haokun/.elan/toolchains/leanprover--lean4---v4.28.0/src/lean/Init/Data/List/Sort/Basic.lean` (核心库)
以及 Mathlib 补充

**完整签名** (核心库):
```lean
def List.mergeSort : ∀ (xs : List α) (le : α → α → Bool := by exact fun a b => a ≤ b), List α
```

**详细分析**:
- **参数**:
  - `xs : List α` — 待排序 List
  - `le : α → α → Bool` — **Bool 式** less-or-equal 比较（与 `Array.qsort` 的 `lt` 略有不同）
  - 默认：`(· ≤ ·)`
  
- **返回值**: 排序后的 List
  
- **比较器格式**: Bool 式，但约定为 `≤` 而非 `<`（细微差异）
  
- **Computable**: ✅ 是
  
- **Partial**: ❌ 否
  
- **稳定性**: ✅ **是** (merge sort 天然稳定)
  
- **Typeclass 需求**: 如使用默认 `(· ≤ ·)`，需要 `[LE α]`；自定义 `le` 则无需求
  
- **in-place**: N/A (链表)
  
- **复杂度**: O(n log n)

**使用示例**:
```lean
#eval List.mergeSort [3, 1, 4, 1, 5]
-- [1, 1, 3, 4, 5]

-- 自定义比较：按降序
#eval List.mergeSort (fun a b => b ≤ a) [3, 1, 4, 1, 5]
-- [5, 4, 3, 1, 1]
```

**vs `insertionSort`**: 时间复杂度更好 O(n log n)，且同样稳定。

---

## 3. Batteries 扩展

### 3.1 Array.merge

**位置**: `/home/haokun/projects/CLPoly/proof/lean/.lake/packages/batteries/Batteries/Data/Array/Merge.lean`

**完整签名**:
```lean
def Array.merge (lt : α → α → Bool) (xs ys : Array α) : Array α
```

**详细分析**:
- **用途**: 归并两个已排序的数组
- **参数**:
  - `lt : α → α → Bool` — less-than 比较器
  - `xs, ys : Array α` — 两个输入数组（必须已按 `lt` 排序）
  
- **返回值**: 归并后的数组（已排序）
  
- **复杂度**: O(|xs| + |ys|)
  
- **稳定性**: ✅ **是** (保留原相对顺序)
  
- **注意**: 仅用于 *已排序* 数组，不能用作通用排序
  
**使用场景**: 分治排序的 merge 阶段，或实现自定义高级排序（如 Tim sort）。

---

### 3.2 Array.sortDedup

**位置**: 同上 (`Merge.lean`)

**完整签名**:
```lean
def Array.sortDedup [ord : Ord α] (xs : Array α) : Array α
```

**详细分析**:
- **用途**: 排序 + 去重（一步到位）
- **参数**: `xs : Array α`，需要 `[Ord α]` 实例
- **返回值**: 排序且去重的数组
- **复杂度**: O(n log n)
- **实现**: `dedupSorted (xs.qsort (compare · · |>.isLT))`
  - 先调用 `Array.qsort`，再调用 `dedupSorted`
  
**使用场景**: 需要排序集合（无重复元素）时。

---

## 对比速查表

| 特性 | `Array.qsort` | `Array.insertionSort` | `Array.qsortOrd` | `List.insertionSort` | `List.mergeSort` |
|-----|-------|---------|-------|-------|-------|
| **位置** | Init | Init | Init | Mathlib | Init |
| **比较器类型** | `Bool` 式 (`lt`) | `Bool` 式 (`lt`) | 需 `Ord` 实例 | `Prop` 式 + `DecidableRel` | `Bool` 式 (`le`) |
| **接受自定义比较** | ✅ 直接 | ✅ 直接 | ❌ 仅 `Ord` 实例 | ⚠ 需 `DecidableRel` | ✅ 直接 |
| **Computable** | ✅ | ✅ | ✅ | ✅ | ✅ |
| **In-place** | ❌ | ❌ | ❌ | N/A | N/A |
| **稳定** | ❌ | ✅ | ❌ | ✅ | ✅ |
| **时间复杂度** | O(n log n) avg | O(n²) | O(n log n) avg | O(n²) | O(n log n) |
| **最坏复杂度** | O(n²) | O(n²) | O(n²) | O(n²) | O(n log n) |
| **适用场景** | 通用快速排序，需自定 cmp | 小数组/近序，需稳定 | 有 `Ord` 类型 | 理论验证（Prop 式） | 稳定排序/大数据 |

---

## 4. C++ → Lean 映射建议

### 4.1 直接映射

**C++ 代码**:
```cpp
std::vector<T> v = {...};
std::sort(v.begin(), v.end(), [](const T& a, const T& b) { 
    return cmp(a, b);  // 返回 bool，true 表示 a < b
});
```

**Lean 映射（推荐）**:
```lean
def cmp (a b : T) : Bool := ...  -- 自定义比较，返回 true/false
let sorted := Array.qsort v cmp   -- 直接对应 std::sort
```

**优点**:
- 最直接的对应
- 完全保留 C++ 的比较器语义
- Computable，可直接 `#eval` 验证

### 4.2 有 Ord 实例的情况

**C++ 代码** (使用 std::less):
```cpp
std::sort(v.begin(), v.end());  // 默认 < 排序
```

**Lean 映射**:
```lean
-- 如果已有 [Ord T] 实例
let sorted := Array.qsortOrd v   -- 使用 Ord 实例
-- 或
let sorted := Array.qsort v (fun a b => compare a b |>.isLT)
```

### 4.3 需要稳定性

**如 C++ 代码依赖稳定排序**:
```lean
-- 用 insertionSort (小数组) 或 mergeSort (List)
let sorted := Array.insertionSort v cmp    -- Array, 稳定
let sorted := List.mergeSort cmp v.toList  -- List, 稳定，需转换
```

---

## 5. 常见问题

### Q1: CLPoly 应该用哪个 API？

**答**: 
- **首选**: `Array.qsort cmp arr`
  - 直接对应 C++ `std::sort` 的语义
  - 自定义比较器最灵活
  - Computable，可验证
  
- **次选**: `Array.insertionSort` (如果需要稳定性或数组很小)

### Q2: CLPoly 的比较器返回 bool，对应哪个签名？

**答**: **Bool 式** (less-than)，对应：
- `Array.qsort : Array α → (α → α → Bool) → Array α`
- `Array.insertionSort : Array α → (α → α → Bool) → Array α`
- `List.mergeSort : List α → (α → α → Bool) → List α`

**不是** Mathlib 的 `List.insertionSort` (需要 `Prop` + `DecidableRel`)。

### Q3: 排序是否必须稳定？

**答**: CLPoly 的算法不依赖稳定性。如果偶发用到（例如分层排序），可用 `insertionSort` 或 `mergeSort`。

### Q4: 是否有 `Array.sortBy : Array α → (α → β) → Array α`（按字段排序）？

**答**: 无直接 `sortBy`。替代方案：
```lean
-- 按字段 f 排序
let sorted := Array.qsort arr (fun a b => f a < f b)
-- 或用 Ord derive 在结构体上
deriving Ord  -- 自动按字段字典序
```

### Q5: Computable 意味着什么？

**答**: 能直接 `#eval` 执行和验证结果。对 L1 验证很重要（可在 Lean 中执行 C++ 算法的形式化模型）。

---

## 附录：Lean 类型类速查

| 类型类 | 签名 | 何时使用 |
|-------|------|--------|
| `LT α` | `(· < ·) : α → α → Prop` | 定义小于关系，通常结合 `DecidableLE` |
| `LE α` | `(· ≤ ·) : α → α → Prop` | 定义小于等于关系 |
| `Ord α` | `compare : α → α → Ordering` | 全序关系（返回 `.lt`/`.eq`/`.gt`） |
| `DecidableEq α` | 相等性可决定 | 需要判断 `a = b` 时 |
| `DecidableRel (· ≤ ·)` | 关系可决定 | Mathlib `insertionSort` 必需 |

---

## 总结建议

**CLPoly L1 翻译**应采用：
1. **主要排序**: `Array.qsort arr cmp`（直接对应 C++ `std::sort`）
2. **自定义类型 + 需稳定**: `Array.insertionSort arr cmp`
3. **List 上下文**: `List.mergeSort cmp l`（稳定 + O(n log n)）
4. **有 Ord 实例**: 可选用 `Array.qsortOrd`（但不如显式 `cmp` 灵活）

所有 API 都是 **computable**，可直接 `#eval` 验证，符合 L1 的验证需求。

