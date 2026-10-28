# Mathlib Gap Report — CLPoly 因式分解形式化验证

> Phase 0 交付物。基于 E1-E4 实验结果，评估 Mathlib4 对 CLPoly DDF/EDF/Hensel 形式化验证的支持程度。

**Mathlib 版本**: stable (2026-03, Lean 4.28.0)
**评估日期**: 2026-03-14

---

## 1. 可直接使用的 Mathlib API

### 1.1 基础代数结构

| API | 签名 | 状态 | 备注 |
|-----|------|------|------|
| `Field (ZMod p)` | 需要 `Fact (Nat.Prime p)` | ✅ | |
| `EuclideanDomain (Polynomial (ZMod p))` | noncomputable | ✅ | 不影响证明 |
| `GCDMonoid (Polynomial (ZMod p))` | noncomputable | ✅ | |
| `DecidableEq (Polynomial (ZMod p))` | | ✅ | |
| `UniqueFactorizationMonoid (Polynomial (ZMod p))` | | ✅ | |
| `CommRing (ZMod (p^k))` | 非 IsDomain | ✅ | Hensel 提升需要 |

### 1.2 多项式运算

| API | 签名 | 用途 |
|-----|------|------|
| `divByMonic` | `Ring R → R[X] → R[X] → R[X]` | DDF 主循环、Hensel 提升 |
| `modByMonic` | `Ring R → R[X] → R[X] → R[X]` | powmod 实现 |
| `modByMonic_add_div` | `p %ₘ q + q * (p /ₘ q) = p` (需 Monic q) | 除法恒等式 |
| `Monic.map` | `p.Monic → (map f p).Monic` | 层间投影保持首一 |
| `GCDMonoid.gcd` | 归一化 gcd | DDF gcd 步骤 |
| `normalize_gcd` | `normalize (gcd a b) = gcd a b` | gcd 首一性 |
| `derivative` | 标准导数 | Separable 证明 |

### 1.3 有限域理论

| API | 签名 | 用途 |
|-----|------|------|
| `roots_X_pow_card_sub_X` | `(X^{\|K\|} - X).roots = univ.val` | Thm 2.1 核心 |
| `splits_X_pow_card_sub_X` | `(map ... (X^{\|K\|} - X)).Splits` | 分裂域 |
| `isSplittingField_sub` | `IsSplittingField F K (X^{\|K\|} - X)` | 域构造 |
| `nonempty_algHom_iff_finrank_dvd` | `Nonempty (K →ₐ[F] L) ↔ finrank F K ∣ finrank F L` | 域嵌入 |
| `GaloisField` | `(p : ℕ) → [Fact (Nat.Prime p)] → ℕ → Type` | 有限域构造 |
| `GaloisField.card` | `Nat.card (GaloisField p n) = p ^ n` | 基数 |
| `GaloisField.finrank` | `finrank (ZMod p) (GaloisField p n) = n` | 维度 |
| `expand_card` | `(expand K \|K\|) f = f ^ \|K\|` | Frobenius |

### 1.4 因式分解相关

| API | 签名 | 用途 |
|-----|------|------|
| `Irreducible` | `¬IsUnit p ∧ ∀ a b, p = a * b → IsUnit a ∨ IsUnit b` | DDF/EDF 正确性 |
| `Squarefree` | `∀ d, d * d ∣ a → IsUnit d` | 预处理条件 |
| `Separable` | `IsCoprime f (derivative f)` | 等价于 squarefree（域上） |
| `Separable.squarefree` | `f.Separable → Squarefree f` | |
| `squarefree_iff_nodup_normalizedFactors` | UFD 中 squarefree ⟺ normalizedFactors 无重复 | |

### 1.5 ZMod (p^k) 基础设施

| API | 签名 | 用途 |
|-----|------|------|
| `ZMod.castHom` | `m ∣ n → ZMod n →+* ZMod m` | Hensel 层间投影 |
| `ZMod.val` | `ZMod n → ℕ` | 值提取 |
| `ZMod.valMinAbs` | `ZMod n → ℤ` | 对称模表示 |
| `ZMod.unitOfCoprime` | `x.Coprime n → (ZMod n)ˣ` | 单位构造 |
| `Polynomial.map` | 系数环同态推广到多项式环 | |

---

## 2. 需要自建的部分

### 2.1 DDF 算法模型 (~200-300 行)

Mathlib 没有 DDF 算法的形式化。需要自建：

- `ddfLoop` 递归函数定义
- 循环不变量（gcd 条件、因子收集完整性）
- 正确性定理：`∀ g, Irreducible g → natDeg g = d → g ∣ gd`

**技术方案确认**: `termination_by n + 1 - 2 * d` + 自动 `.induct` 归纳原理

### 2.2 EDF 算法模型 (~100-200 行)

Cantor-Zassenhaus 随机分裂。需要自建：

- `partial def` 或 `fuel` 参数方案
- 正确性条件化证明（假设找到非平凡因子时的正确性）

### 2.3 Hensel 提升 (~1000-1500 行)

**这是最大的缺口。** Mathlib 只有 p-adic 标量根提升，没有多项式因子提升。

需要自建：
- Hensel 提升引理：模 p^k 的因式分解 → 模 p^{k+1} 的因式分解
- 唯一性定理
- 迭代构造

**可参考**: Isabelle AFP 的 Berlekamp-Zassenhaus 形式化

### 2.4 X^{p^d} - X 的 Separable 性 (~20-50 行)

实验中用 `sorry` 占位。需要证明：

```
derivative(X^{p^d} - X) = p^d · X^{p^d-1} - 1 = -1  (char p)
gcd(X^{p^d} - X, -1) = 1 → Separable
```

需要的引理：`CharP.cast_eq_zero`、`derivative_X_pow`、`derivative_sub`，全部已确认可用。

### 2.5 "不可约 g, deg g = d → g ∣ (X^{p^d} - X)" 桥接引理 (~50-100 行)

这是 Thm 2.1 的核心桥接步骤。思路：
1. `deg g = d` → g 的分裂域是 `GaloisField p d`
2. `GaloisField p d` 中所有元素是 `X^{p^d} - X` 的根
3. g 的所有根都是 `X^{p^d} - X` 的根
4. g ∣ (X^{p^d} - X)（因为域上多项式由根决定）

所有中间步骤的 API 已确认可用。

### 2.6 IsUnit (a : ZMod (p^k)) 的构造方法 (~10-20 行)

`ZMod.unitOfCoprime` 已确认可用，但需要方便的 wrapper。实验中 `sorry` 占位的 `IsUnit (3 : ZMod 343)` 需要通过 `Nat.Coprime` 证明。

---

## 3. 风险评估

| 风险 | 级别 | 缓解措施 |
|------|------|---------|
| Hensel 提升工作量超预期 | 🟡 中 | 可先用 `axiom` 占位，Phase 3 再填充 |
| Thm 2.1 桥接引理细节困难 | 🟢 低 | 所有中间步骤 API 已确认，主要是组合工作 |
| `noncomputable` 阻碍测试 | 🟢 低 | 用 `Nat` 版本做 `#eval` 测试，多项式版本只做证明 |
| 编译时间过长 | 🟢 低 | 增量编译 ~1 秒，预编译缓存已下载 |
| `termination_by` 不够表达力 | ✅ 已消除 | 实验确认自动接受 |
| `.induct` 不可用 | ✅ 已消除 | 实验确认自动生成，4 分支 |

---

## 4. 结论

**Mathlib 覆盖率: ~70%**。基础代数结构、有限域理论、多项式运算、因子谓词全部就绪。主要缺口是算法模型层（DDF/EDF/Hensel），这是预期的——Mathlib 是数学库而非算法库。

**建议**：按路线图继续执行 Phase 1（DDF 端到端验证），无需调整技术方案。所有 Phase 0 假设均已验证通过。
