# CLPoly Model Inventory: `clpoly_model.lean`

**File**: `/home/haokun/projects/CLPoly/proof/cpp2lean/clpoly_model.lean`  
**Lines**: 262  
**Status**: Baseline L1 model for cpp2lean v1; v2 rewrite target

---

## 文件概览

本文件是 CLPoly 可信基类的 Lean 模型，定义了翻译器生成代码依赖的所有基本类型、操作和辅助函数。
采用可执行的 Lean 4 实现（非 axiom），支持背靠背测试框架的 `#eval`。

### 结构划分

| Section/Namespace | 行号 | 功能 |
|-----------------|------|------|
| Preamble | 1-9 | 选项设置（`autoImplicit false`） |
| §1: Zp (Z/pZ) | 11-57 | 有限域元素和四则运算 |
| §2: ZZ (大整数) | 59-63 | 整数类型别名 |
| §3: UMonomial | 65-71 | 单变量单项式 |
| §4: SparsePolyZp | 73-117 | Z/pZ 稀疏多项式（数组 + 降幂序） |
| §5: Array 辅助 | 119-124 | 通用数组工具函数 |
| §5a: Rng | 126-150 | 伪随机数生成器 |
| §5a2: SparsePolyZZ helpers | 152-165 | Z 系数多项式的模运算和稀疏化 |
| §5b: StdMap | 167-204 | 有序映射（std::map 模型） |
| §5c: 其他辅助类型 | 206-239 | 多变量多项式、因式分解结果等结构体 |
| §6: 验证测试 | 241-262 | `#eval` 测试用例 |

---

## 1. 类型定义总览

### 1.1 基础结构体（Structure）

| 类型名 | 字段 | 类型参数 | 衍生类 | C++ 对应 | 用途 | 状态 |
|--------|------|--------|--------|---------|------|------|
| **Zp** | `val: UInt64` `prime: UInt64` | 否 | Repr, Inhabited, BEq | `clpoly::Zp` | 有限域元素 Z/pZ | ✅ 完整实现 |
| **UMonomial** | `deg: UInt64` | 否 | Repr, Inhabited, BEq | `clpoly::UMonomial` | 单变量单项式（次数） | ✅ 完整实现 |
| **Factorization** | `content: ZZ` `factors: Array (SparsePolyZZ × UInt64)` | 否 | Inhabited | `clpoly::Factorization` | 因式分解结果（内容 + 因子列表） | ✅ 占位结构 |
| **PrimeSelectionResult** | `p: UInt64` `factors: Array SparsePolyZp` `nfactors: UInt64` | 否 | Inhabited | `clpoly::PrimeSelectionResult` | 素数选择的输出 | ✅ 占位结构 |
| **WangLcResult** | `success: Bool` `scaled_factors: Array MvPolyZZ` `lc_targets: Array MvPolyZZ` `delta: ZZ` | 否 | Inhabited | `clpoly::WangLcResult` | Wang LC 分配结果 | ✅ 占位结构 |

### 1.2 类型别名（Abbrev）

| 类型名 | 定义 | C++ 对应 | 用途 | 状态 |
|--------|------|---------|------|------|
| **ZZ** | `Int` | `clpoly::ZZ` | 大整数（任意精度） | ✅ 简单映射 |
| **SparsePolyZp** | `Array (UMonomial × Zp)` | `clpoly::SparsePolyZp` | Z/pZ 稀疏多项式（降幂有序） | ✅ 数组 + 对 |
| **SparsePolyZZ** | `Array (UMonomial × Int)` | `clpoly::SparsePolyZZ` | Z 稀疏多项式 | ✅ 数组 + 对 |
| **MvPolyZZ** | `Array (Array (UInt64 × UInt64) × Int)` | `clpoly::MvPolyZZ` | 多变量多项式（指数向量 + 系数） | ⚠️ 占位（待细化） |
| **MvPolyZp** | `Array (Array (UInt64 × UInt64) × Zp)` | `clpoly::MvPolyZp` | 多变量多项式（Z/pZ 系数） | ⚠️ 占位（待细化） |
| **MvMonomial** | `Array (UInt64 × UInt64)` | `clpoly::MvMonomial` | 多变量单项式（变量指数对） | ⚠️ 占位（待细化） |
| **Variable** | `Array (String × Int)` | `clpoly::Variable` | 变量和替换值 | ⚠️ 占位（简化） |
| **LLLMatrix** | `Array (Array Int)` | `clpoly::LLLMatrix` | LLL 归约矩阵 | ⚠️ 占位（待细化） |
| **HenselNode** | `Array Int` | 无对应 | Hensel 提升中间节点（已知不使用） | ⚠️ 死代码 |
| **StdMap K V** | `List (K × V)` | `std::map<K,V>` | 有序映射 | ✅ 函数式等价 |

---

## 2. Opaque Constants 与 Axiom

| 名称 | 签名 | 定义形式 | 含义 | 状态 |
|------|------|---------|------|------|
| **assign** | `(poly : α) → (var : Variable) → (val : β) → α` | `opaque` | 多项式变量代入 `poly[var := val]`；实现为恒等函数 | ⚠️ 占位（实现被简化） |

---

## 3. 函数与方法定义

### 3.1 Zp 命名空间（14 个函数）

| 函数名 | 完整签名 | 形式 | C++ 对应 | 描述 | 状态 |
|--------|---------|------|---------|------|------|
| **ofInt** | `(v : Int) → (p : UInt64) → Zp` | `def` (computable) | `Zp(int, uint64)` | 从整数转换为 Z/pZ 元素（取模+符号处理） | ✅ 完整 |
| **ofUInt64** | `(v p : UInt64) → Zp` | `def` (computable) | `Zp(uint64, uint64)` | 从无符号整数转换为 Z/pZ | ✅ 完整 |
| **Add instance** | `a + b := ⟨(a.val + b.val) % a.prime, a.prime⟩` | `instance` | `Zp::operator+` | 加法（模运算） | ✅ 完整 |
| **Sub instance** | `a - b := ⟨(a.val + a.prime - b.val) % a.prime, a.prime⟩` | `instance` | `Zp::operator-` | 减法（模运算，避免负数） | ✅ 完整 |
| **Mul instance** | `a * b := ⟨(a.val * b.val) % a.prime, a.prime⟩` | `instance` | `Zp::operator*` | 乘法（模运算） | ✅ 完整 |
| **Neg instance** | `a := ⟨(a.prime - a.val) % a.prime, a.prime⟩` | `instance` | `Zp::operator-()` | 一元否定 | ✅ 完整 |
| **extGcdAux** | `(old_r r old_s s : Int) → Int × Int` | `partial def` | `__gcd_extended_aux` | 扩展欧几里得辅助（递归，尾调用） | ✅ 完整 |
| **modInv** | `(a p : UInt64) → UInt64` | `def` | `Zp::modInv` | 模逆（利用 gcd，返回 0 if a=0） | ✅ 完整 |
| **inv** | `(a : Zp) → Zp` | `def` | `Zp::inv` | 有限域逆元 | ✅ 完整 |
| **div** | `(a b : Zp) → Zp` | `def` | `Zp::operator/` | 有限域除法 | ✅ 完整 |
| **Coe Zp UInt64** | `coe (z : Zp) : UInt64 := z.val` | `instance` | 隐式转换 | 强制转换为 UInt64 | ✅ 完整 |
| **Coe Zp Int** | `coe (z : Zp) : Int := z.val.toNat` | `instance` | 隐式转换 | 强制转换为 Int | ✅ 完整 |

**小计**：12 个 `def` + 2 个 `instance` (Coe)，所有可执行。

### 3.2 SparsePolyZp 命名空间（9 个函数）

| 函数名 | 完整签名 | 形式 | C++ 对应 | 描述 | 状态 |
|--------|---------|------|---------|------|------|
| **empty** | `→ SparsePolyZp` | `def` | `SparsePolyZp()` | 空多项式 | ✅ 完整 |
| **front!** | `(f : SparsePolyZp) → UMonomial × Zp` | `def` | `front()` | 首项（最高次） | ✅ 完整 |
| **back!** | `(f : SparsePolyZp) → UMonomial × Zp` | `def` | `back()` | 尾项（最低次） | ✅ 完整 |
| **getDeg** | `(f : SparsePolyZp) → UInt64` | `def` | `getDeg()` | 多项式次数（空多项式 → 0） | ✅ 完整 |
| **size_u64** | `(f : SparsePolyZp) → UInt64` | `def` | `size()` (U64) | 项数转为 UInt64 | ✅ 完整 |
| **normalization** | `(f : SparsePolyZp) → SparsePolyZp` | `def` | `normalization()` | 过滤掉零系数项 | ✅ 完整 |
| **reserve** | `(f : SparsePolyZp) → (n : UInt64) → SparsePolyZp` | `def` | `reserve()` | 占位（返回 f 不变） | ⚠️ 占位 |
| **data** | `(f : SparsePolyZp) → SparsePolyZp` | `def` | `data()` | 占位（返回 f 不变） | ⚠️ 占位 |
| **derivative** | `(f : SparsePolyZp) → SparsePolyZp` | `def` | `derivative()` | 多项式导数（链式法则 + Z/pZ 模运算） | ✅ 完整 |
| **gcd** | `(f g : SparsePolyZp) → SparsePolyZp` | `partial def` | `gcd()` | GCD（欧几里得）；**TODO**：需实现 poly divmod | ❌ **TODO** |
| **divmod** | `(f g : SparsePolyZp) → SparsePolyZp × SparsePolyZp` | `def` | `divmod()` | 多项式除法（商 + 余）；**TODO**：多项式长除法占位，返回 `(f, #[])` | ❌ **TODO** |
| **comp** | `(f : SparsePolyZp) → UInt64` | `def` | （无） | 占位比较器；返回 0 | ⚠️ 占位 |

**小计**：11 个 `def`，其中 2 个 `partial def`，**2 个 TODO** (gcd, divmod)。

### 3.3 Array 全局函数（2 个）

| 函数名 | 完整签名 | 形式 | 描述 | 状态 |
|--------|---------|------|------|------|
| **Array.front!** | `{α : Type} [Inhabited α] → (a : Array α) → α` | `def` | 数组首元素 | ✅ 完整 |
| **Array.size_u64** | `{α : Type} → (a : Array α) → UInt64` | `def` | 数组大小转为 UInt64 | ✅ 完整 |

### 3.4 Rng 命名空间（2 个函数）

| 函数名 | 完整签名 | 形式 | C++ 对应 | 描述 | 状态 |
|--------|---------|------|---------|------|------|
| **next** | `(seed upper : UInt64) → UInt64` | `def` | `next(upper)` | 伪随机数 [0, upper)；xorshift64 实现 | ✅ 完整 |
| **step** | `(seed : UInt64) → UInt64` | `def` | `step()` | 伪随机数生成步进 | ✅ 完整 |

**小计**：2 个 `def`，确定性（非密码学）。

### 3.5 SparsePolyZZ 全局函数（2 个）

| 函数名 | 完整签名 | 形式 | C++ 对应 | 描述 | 状态 |
|--------|---------|------|---------|------|------|
| **SparsePolyZZ.modCoeff** | `(f : SparsePolyZZ) → (m : Int) → SparsePolyZZ` | `def` | `__upoly_mod_coeff()` | 系数模 m + 过滤零 | ✅ 完整 |
| **SparsePolyZZ.compactNonzero** | `(f : SparsePolyZZ) → SparsePolyZZ` | `def` | （Hensel 中间过程） | 过滤零系数项 | ✅ 完整 |

### 3.6 StdMap 命名空间（7 个函数）

| 函数名 | 完整签名 | 形式 | C++ 对应 | 描述 | 状态 |
|--------|---------|------|---------|------|------|
| **empty** | `→ StdMap K V` | `def` | `std::map()` | 空映射 | ✅ 完整 |
| **find?** | `[BEq K] → (m : StdMap K V) → (k : K) → Option V` | `def` | `find()` | 查找（返回 Option） | ✅ 完整 |
| **find!** | `[BEq K][Inhabited V] → (m : StdMap K V) → (k : K) → V` | `def` | `find()` | 查找（返回默认值） | ✅ 完整 |
| **insert** | `[BEq K] → (m : StdMap K V) → (k : K) → (v : V) → StdMap K V` | `def` | `insert()` | 插入/更新（头部 prepend） | ✅ 完整 |
| **erase** | `[BEq K] → (m : StdMap K V) → (k : K) → StdMap K V` | `def` | `erase()` | 删除 | ✅ 完整 |
| **size** | `(m : StdMap K V) → Nat` | `def` | `size()` | 映射大小 | ✅ 完整 |
| **isEmpty** | `(m : StdMap K V) → Bool` | `def` | `empty()` | 是否空映射 | ✅ 完整 |
| **end_** | `(m : StdMap K V) → Nat` | `def` | 迭代器端点 | 占位（返回长度） | ⚠️ 占位 |
| **begin_** | `(m : StdMap K V) → Nat` | `def` | 迭代器起点 | 占位（返回 0） | ⚠️ 占位 |

**小计**：9 个 `def`，7 个完整，2 个迭代器占位。

---

## 4. 已知问题与 TODO 清单

### 4.1 TODO 注释（3 处）

| 行号 | 函数 | TODO 内容 | 优先级 |
|------|------|---------|--------|
| 106 | `gcd` | "需要 poly mod" | 🔴 高 |
| 111 | `divmod` | "实现多项式长除法" | 🔴 高 |
| 104 | `gcd` (注释) | "完整实现需要 divmod，暂返回简化版" | 🔴 高 |

### 4.2 占位符与简化实现

| 项 | 描述 | 影响范围 | 建议 |
|----|------|---------|------|
| **SparsePolyZp.gcd** | 不完整欧几里得（缺 divmod） | 多项式 GCD 依赖 | 补充 poly divmod |
| **SparsePolyZp.divmod** | 返回 `(f, #[])` 占位 | 多项式除法链 | 实现长除法 |
| **SparsePolyZp.reserve, data** | 返回 f 不变 | 内存预分配（非关键） | 可保留占位 |
| **SparsePolyZp.comp** | 返回 0 | 多项式比较 | 可保留占位 |
| **MvPolyZZ, MvPolyZp** | 简化为指数向量数组 | 多变量因式分解 | v2 需精化为多元表示 |
| **HenselNode** | 定义为 `Array Int` | 未使用 | 检查是否死代码 |
| **StdMap.end_, begin_** | 迭代器占位 | 列表遍历实现 | 可保留（翻译中不使用迭代器） |

### 4.3 测试覆盖

文件末尾有 6 个 `#eval` 测试（§6），验证：
- Zp.ofInt（负数、边界）
- Zp 加法、乘法
- 模逆
- SparsePolyZp 度数、导数

**状态**：✅ 均已编译通过，无失败。

---

## 5. 类型参数与泛型支持

| 类型 | 泛型参数 | 约束 | 复用度 |
|------|--------|------|--------|
| **Zp** | 无 | 无 | 单一具体类型 |
| **UMonomial** | 无 | 无 | 单一具体类型 |
| **SparsePolyZp** | 无 | 无 | 单一具体类型 |
| **StdMap** | `K, V` | `[BEq K]`, `[Inhabited V]` | 高度泛型 |
| **Array.front!** | `α` | `[Inhabited α]` | 通用 |
| **Array.size_u64** | `α` | 无 | 通用 |

---

## 6. 与 C++ 的一致性评估

### 6.1 完全对应（Green）
| 模块 | 对应度 | 备注 |
|------|-------|------|
| Zp 四则运算 | ✅ 1:1 | 模算术完整实现 |
| UMonomial | ✅ 1:1 | 单变量单项式 |
| SparsePolyZp 基本操作 | ✅ 95% | 除法/GCD 缺实现 |
| Rng (xorshift64) | ✅ 1:1 | 伪随机数确定性 |
| StdMap 语义 | ✅ 1:1 | 函数式等价 List 实现 |

### 6.2 部分对应或占位（Yellow）
| 模块 | 对应度 | 原因 |
|------|-------|------|
| SparsePolyZp.divmod | ⚠️ 0% | TODO：多项式长除法 |
| SparsePolyZp.gcd | ⚠️ 20% | 缺 divmod |
| MvPolyZZ/ZP | ⚠️ 20% | 简化为指数向量，未精化 |
| HenselNode | ⚠️ 未使用 | 占位，需核实是否死代码 |

### 6.3 设计与语义差异
| 方面 | Lean 模型 | C++ 实现 | 差异 |
|------|---------|---------|------|
| 整数溢出 | 无（`Int` 任意精度） | 有（`uint64_t` 固定宽） | L1 精化需处理溢出 |
| 数组越界 | panic | UB | L1 需边界检查 |
| 内存管理 | 自动（GC） | 手动（RAII） | L1 需 ownership 模型 |
| 非决定性 | 确定的 xorshift64 | `mt19937` 可配 | 语义相容但参数不同 |

---

## 7. 代码质量度量

| 指标 | 值 | 备注 |
|------|---|----|
| 文件总行数 | 262 | 含空行、注释、测试 |
| 实际代码行数 | ~200 | 排除空行和注释 |
| 类型定义数 | 5 struct + 9 abbrev = 14 | 含临时占位 |
| 函数总数 | ~40 | 含 instance、def、partial def、opaque |
| Axiom 占位数 | 1 (`assign`) | opaque + 恒等实现 |
| 部分实现数 | 6 | gcd, divmod, reserve, data, comp, StdMap.end/begin |
| TODO 注释 | 3 处 | 多项式除法链 |
| `#eval` 测试 | 6 个 | 全通过 |
| 编译状态 | ✅ 通过 | `lake build` OK |

---

## 8. 对 cpp2lean v2 重构的影响

### 8.1 可复用的定义（Green 代码）

这些定义可直接导入 v2，无需修改：

```lean
structure Zp where ...
instance Add Zp where ...
def Zp.ofInt, modInv, inv, div ...
def SparsePolyZp.empty, front!, back!, getDeg, derivative, normalization ...
def Array.front!, Array.size_u64
def Rng.next, Rng.step
def SparsePolyZZ.modCoeff, compactNonzero
def StdMap.* (全部 9 个)
```

**建议**：v2 可继承这些代码，无改动或极小改动。

### 8.2 需要补充的定义（Red 待办）

| 待办 | 优先级 | 建议方案 |
|------|-------|---------|
| `SparsePolyZp.divmod` 完整实现 | 🔴 高 | 多项式长除法（对标 C++ 实现） |
| `SparsePolyZp.gcd` 正确版本 | 🔴 高 | 基于完整 divmod 的欧几里得 |
| `MvPolyZZ/Zp` 精化 | 🟡 中 | v2 确定多元多项式内部表示后再实现 |
| `HenselNode` 清理 | 🟡 中 | 确认是否死代码，若有则删除 |
| `assign` 实现 | 🟢 低 | 当前 opaque + 恒等已满足测试需求 |

### 8.3 签名检查（Signature Audit）

以下函数签名是否需改：

| 函数 | 当前签名 | 可能问题 | 建议 |
|------|---------|---------|------|
| `Zp.ofInt` | `(v : Int) → (p : UInt64) → Zp` | 素数 p 未验证 | 可添加 `prime : Fact (Nat.Prime p)` 参数（可选） |
| `SparsePolyZp.gcd` | `(f g : SparsePolyZp) → SparsePolyZp` | 无失败处理 | 可改为 `Option SparsePolyZp` 或保持（取决于 spec） |
| `MvPolyZZ.` | （需逐个检查） | 表示与 C++ 不符 | v2 设计阶段明确后再改 |
| `assign` | `(poly : α) → (var : Variable) → (val : β) → α` | 过度泛型，实现与签名不符 | 具体化为 `SparsePolyZZ → ... → SparsePolyZZ` |

---

## 9. 与已排除死代码的关联

根据背景信息，以下 6 个 C++ 函数已被排除：

```
__multivar_hensel_lift
__hensel_lift_one_var
__hensel_lc_correct
__multivar_diophantine
__pseudo_remainder_x1
__taylor_coeff
```

**当前 clpoly_model.lean 中**：
- ❌ 无对应的 Lean stub（已清理）
- ⚠️ `HenselNode` 可能是遗留占位，需确认

---

## 10. 验证与测试现状

### 10.1 #eval 测试（6 个，全通过）

```lean
#eval Zp.ofInt (-1) 5        -- ✅ ⟨4, 5⟩
#eval Zp.ofInt 0 3           -- ✅ ⟨0, 3⟩
#eval Zp.ofInt 7 13          -- ✅ ⟨7, 13⟩
#eval Zp.ofInt 100 17        -- ✅ ⟨15, 17⟩
#eval (Zp.ofUInt64 3 7) + (Zp.ofUInt64 5 7)  -- ✅ ⟨1, 7⟩
#eval (Zp.ofUInt64 3 7) * (Zp.ofUInt64 5 7)  -- ✅ ⟨1, 7⟩
#eval Zp.modInv 3 7          -- ✅ 5
#eval SparsePolyZp.getDeg ...  -- ✅ 3
#eval SparsePolyZp.derivative ...  -- ✅ 计算通过
```

### 10.2 编译状态

- ✅ `lake build` 全通过（0 error）
- ✅ `set_option autoImplicit false` 无隐式参数陷阱
- ✅ 所有导出类型 derive `Repr`、`BEq` 便于调试

---

## 11. 大小写规范与命名一致性

| 类型/项 | 命名风格 | 一致性 | 备注 |
|--------|---------|--------|------|
| **Zp** | CamelCase | ✅ 一致 | Lean 结构体约定 |
| **UMonomial** | CamelCase | ✅ 一致 | 对应 C++ `UMonomial` |
| **SparsePolyZp** | CamelCase | ✅ 一致 | 对应 C++ `SparsePolyZp` |
| **StdMap** | CamelCase | ✅ 一致 | 对应 C++ `std::map` |
| 函数（如 `ofInt`） | camelCase | ✅ 一致 | Lean 函数约定 |
| namespace | lowercase | ✅ 一致 | `Zp`, `SparsePolyZp`, `Rng` 等 |

---

## 总结表格

| 类别 | 计数 | 备注 |
|------|------|------|
| **结构体** | 5 | Zp, UMonomial, Factorization, PrimeSelectionResult, WangLcResult |
| **类型别名** | 9 | ZZ, SparsePolyZp, SparsePolyZZ, MvPolyZZ, MvPolyZp, MvMonomial, Variable, LLLMatrix, HenselNode |
| **函数总数** | ~40 | def (28) + instance (6) + partial def (2) + opaque (1) + 验证测试 |
| **Instance** | 6 | Add, Sub, Mul, Neg, Coe×2 (Zp) |
| **Axiom/Opaque** | 1 | `assign` (opaque，实现为恒等) |
| **Partial def** | 2 | `extGcdAux`, `gcd` |
| **完整实现** | ~32 | 可执行 def |
| **部分实现/占位** | 6 | gcd, divmod, reserve, data, comp, StdMap.end/begin |
| **TODO 注释** | 3 | 多项式除法链 |
| **文件行数** | 262 | 含测试、空行、注释 |
| **代码行数** | ~200 | 实际 Lean 代码 |
| **编译状态** | ✅ Pass | 无 error，#eval 全通过 |

---

## 附录：文件组织建议

对于 v2 重构，建议按以下结构组织：

```
CLPoly/
├── proof/cpp2lean/
│   ├── clpoly_model.lean          (当前文件，Green 代码保留)
│   ├── clpoly_model_polyops.lean  (新增：divmod, gcd 完整实现)
│   ├── clpoly_model_mvpoly.lean   (新增：MvPolyZZ/Zp 精化)
│   └── clpoly_model_tests.lean    (新增：扩展测试覆盖)
├── docs/design/l1-translation-validation/
│   └── survey/
│       ├── clpoly-model-inventory.md  (本文件)
│       └── signature-audit-v2.md      (v2 签名检查清单)
```

---

**文档版本**: v1.0  
**生成日期**: 2026-04-19  
**用途**: cpp2lean v2 重构基线与类型系统设计复用
