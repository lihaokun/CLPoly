# L0 — `__make_zp`

**审视日期**：2026-05-03
**状态**：✅ 一致

## 1. `__make_zp`

**C++ 源**：`clpoly/polynomial_factorize_zp.hh:20`
**Lean 翻译**：`Generated/Corpus.lean:1984-1985`
**关联**：依赖 `Zp.ofInt`（Lean Model）+ C++ `Zp(int64_t, uint64_t)` 构造

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 函数名 | `__make_zp` | `__make_zp_ir` | ✅（`_ir` 后缀是 cpp2lean v2 约定）|
| 参数 1 | `int64_t val` | `(val : Int64)` | ✅ |
| 参数 2 | `uint64_t p` | `(p : UInt64)` | ✅ |
| 返回 | `Zp` | `Zp` | ✅ |

### 行为对照

**C++**：
```cpp
inline Zp __make_zp(int64_t val, uint64_t p) { return Zp(val, p); }
```

→ 调用 `Zp(int64_t, uint64_t)` 构造（`Zp.hh:142-148`）：

```cpp
Zp(int64_t i, uint64_t p) : _p(p) {
    __precompute(p, _ninv, _norm);
    uint64_t abs_i = (i >= 0) ? (uint64_t)i : -(uint64_t)i;
    uint64_t r = abs_i % p;
    _i = (r == 0 || i > 0) ? r : p - r;
}
```

行为：标准 mathematical mod，结果 ∈ [0, p)。

**Lean**：
```lean
partial def __make_zp_ir (val : Int64) (p : UInt64) : Zp :=
  (Zp.ofInt (val).toInt p)
```

→ 调用 `Zp.ofInt`（Model.lean:29-33）：

```lean
def ofInt (v : Int) (p : UInt64) : Zp :=
  let pn : Int := p.toNat
  let r := v % pn
  let r := if r < 0 then r + pn else r
  ⟨r.toNat.toUInt64, p⟩
```

行为：标准 mathematical mod（Int 取模可能负 → +pn），结果 ∈ [0, p)。

### 偏差分析

| 维度 | 评估 |
|------|------|
| 控制流 | ✅ 1 对 1（无控制流） |
| 数据流 | ✅ Int64 → Int (Lean 内置) → mod p → Zp |
| 调用 | ✅ 同一个 ctor 语义（C++ 内联 mod，Lean 通过 `Zp.ofInt`） |
| 边界 | ⚠️ `p = 0` 时 C++ `__precompute` 是 UB，Lean `v % 0 = 0` 返回 `Zp(0, 0)` —— 微小差异，但 caller 应保证 `p > 0`（数学约定 prime > 0）|
| `__precompute` 缺失 | ⚠️ Lean 端没保留 `_ninv` / `_norm`（C++ 用于 Barrett 归约加速）— 不影响数学结果，性能层面差异 |

### 结论

✅ **语义 1-对-1**。C++ 与 Lean 对所有有效输入（`p > 0`）返回相同 Zp 值。

仅有 2 个非语义层面差异：
- 性能（Barrett 归约缺失）
- UB 边界（`p=0` 时行为不同，但 caller 保证 `p > 0` 排除）

无需修复。

---

## 模块小结

| 函数 | 状态 | 备注 |
|------|------|------|
| `__make_zp` | ✅ | 翻译忠实 |

L0 共 1 个函数，全部审视通过。
