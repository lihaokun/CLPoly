# P3-HGCD 架构文档：dense_upoly_zp Half-GCD

> 状态：待确认
> 调研依据：`docs/research/hgcd-research.md`

---

## 1. 核心问题

P3a（稠密表示）和 P3-mul（Karatsuba + lazy divrem）完成后，CLPoly 的 Euclid GCD 在 deg <= 500 已接近 FLINT。但 FLINT 在 deg >= 340 自动切换 HGCD，高次差距随 degree 增长：

| 用例 | CLPoly (Euclid) | FLINT (HGCD) | ratio |
|------|-----------------|--------------|-------|
| gcd deg500+common250 | 0.57ms | 0.40ms | 1.44x |
| gcd deg1000+common500 | 2.27ms | 0.85ms | 2.67x |
| gcd deg3000+common1500 | 20.0ms | 4.70ms | 4.28x |
| gcd deg5000+common2500 | 55.5ms | 10.5ms | 5.28x |
| gcd deg10000+common5000 | 221ms | 32.2ms | 6.87x |

根本原因：Euclid GCD 是 O(n^2)，HGCD 是 O(M(n) log n)。CLPoly 已有 Karatsuba（M(n) = O(n^1.585)），具备 HGCD 收益的前提。

## 2. 设计概览

### 2.1 三层架构（对标 FLINT）

```
Layer 1: gcd 分派
  dense_upoly_zp::gcd(G, A, B)
    if len(B) < GCD_CUTOFF(340):
      _gcd_euclid(...)                  // 现有 Euclid
    else:
      _gcd_hgcd(...)                    // 新增：HGCD-based GCD

Layer 2: HGCD-based GCD 循环（对照 FLINT gcd_hgcd.c:54-98）
  _gcd_hgcd(G, A, B):             // 前置：len(A) >= len(B)
    divrem(_, R, A, B)
    if R == 0: G = B; return       // B 整除 A
    _hgcd(G, J, B, R)             // 首次 HGCD 缩减 (B, R) → (G, J)
    while J != 0:                  // J==0 表示 HGCD 内部已找到 GCD
      divrem(_, R, G, J)
      if R == 0: G = J; break      // J 整除 G
      if len(J) < GCD_CUTOFF:
        _gcd_euclid(G, J, R)      // fallback Euclid
        break
      _hgcd(G, J, J, R)           // 注意：输出 J 别名输入 J（见 §3 Decision 7）

Layer 3: HGCD 递归
  _hgcd_recursive(M, A, B, a, b, W):
    m = len(a) / 2
    if len(b) < m+1: base case
    if len < HGCD_CUTOFF(100):
      _hgcd_iter(...)                   // 迭代 base case
    else:
      递归 → 矩阵合并
```

### 2.2 新增操作

| 操作 | 说明 | 层级 |
|------|------|------|
| `_poly_add` | 原始数组逐系数模加 | 内部 raw API |
| `_poly_sub` | 原始数组逐系数模减 | 内部 raw API |
| `_poly_divrem` | 原始数组长除法（商+余式）| 内部 raw API |
| `_mul` | 原始数组不等长乘法 | 内部 raw API |
| `add` / `sub` | `dense_upoly_zp` 对象级加减 | 公开 API |

### 2.3 不修改的部分

- `dense_upoly_zp` 的现有 API（mul, divrem, gcd 签名不变）
- `polynomial_gcd.hh` 的集成代码（`__polynomial_GCD` 调用 `gcd` 不变）
- 所有上层调用方

## 3. 关键设计决策

### Decision 1: 内部使用 raw pointer API

HGCD 递归过程中需要频繁操作多项式子数组（右移取高位、截断取低位）和工作区切片。使用 `dense_upoly_zp` 对象会引入不必要的拷贝和分配。

**方案**：内部 HGCD 函数全部使用 `uint64_t*` + `size_t len`，与现有 `_classical_mul`、`_kar_mul` 风格一致。仅在 `gcd` 公开入口做对象 <-> raw 转换。

```cpp
// 内部 raw API 示例
size_t _poly_add(uint64_t* C, const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const;
```

返回值为实际输出长度（已 normalize）。

### Decision 2: 预分配工作区 + 指针切片

FLINT 的做法：预分配一大块连续内存，通过指针偏移切分给各层递归使用。

**方案**：在 `_gcd_hgcd` 入口一次性分配 `std::vector<uint64_t>` 工作区，大小 `22 * n + 16 * (ceil_log2(n) + 1)`，传入 raw pointer。内部递归函数自行切片、推进指针。

```cpp
void _gcd_hgcd(uint64_t* G, size_t& len_g,
               const uint64_t* A, size_t len_a,
               const uint64_t* B, size_t len_b) const
{
    size_t alloc = 22 * len_a + 16 * (ceil_log2(len_a) + 1);
    // + gcd 外层循环的临时空间
    std::vector<uint64_t> W(alloc + 4 * len_a);
    // ...
}
```

避免递归过程中反复分配/释放，零额外 heap 操作。

### Decision 3: 2x2 矩阵用 4 个 (pointer, length) 对表示

FLINT 用 `gr_ptr M[4]` + `slong lenM[4]`。我们用相同方式：

```cpp
struct hgcd_mat {
    uint64_t* poly[4];  // [0]=M00, [1]=M01, [2]=M10, [3]=M11
    size_t len[4];
};
```

矩阵的 4 个多项式均指向工作区内的预分配空间（不独立分配）。最大长度 `(n+1)/2`（HGCD 过程中矩阵元素的度数上界）。

### Decision 4: 初始仅实现 classical 矩阵乘法

FLINT 在矩阵最小元素度数 >= 20 时用 Strassen（7 mul vs 8 mul）。对 CLPoly 当前规模（deg <= 10000），矩阵元素度数通常不大，Strassen 收益有限。

**方案**：第一版仅实现 classical（8 次 `_mul` + 4 次 `_poly_add`）。后续若 benchmark 显示矩阵乘法是瓶颈，再加 Strassen。

### Decision 5: cutoff 参数

| 参数 | 初始值 | 说明 |
|------|--------|------|
| `GCD_CUTOFF` | 340 | `gcd` 分派：Euclid vs HGCD |
| `HGCD_CUTOFF` | 100 | HGCD 递归 base case：迭代 vs 递归 |

沿用 FLINT 3.0.1 对 >8-bit 素数的统一值。后续可通过 benchmark 微调。

### Decision 6: `_mul` 在 HGCD 内部的调用方式

HGCD 矩阵乘法需要调用多项式乘法。两个操作数可能不等长，且常为非 2 的幂。

**方案**：新增一个 `_mul` 内部接口，支持不等长输入：

```cpp
// 不等长乘法：len_a >= len_b > 0，C 预分配 len_a + len_b - 1 元素
void _mul(uint64_t* C, const uint64_t* A, size_t len_a,
          const uint64_t* B, size_t len_b, uint64_t* scratch) const;
```

内部根据 `min(len_a, len_b)` 选择 schoolbook 或 Karatsuba。若短操作数 length < KARATSUBA_THRESHOLD，使用 schoolbook（O(len_a * len_b)）；否则零填充到等长后 `_kar_mul`。

**为什么用 min 而非 max**：HGCD 迭代 base case 中 `Q * M[i]`，Q 可能很短（常数商，length=1），M[i] 可能很长。若按 max 判断会零填充 Q 到等长后 Karatsuba，复杂度 O(max^1.585)，远差于 schoolbook O(len_a * len_b)。scratch 由调用方提供。

### Decision 7: `_hgcd_recursive` 必须支持 {B, a} aliasing

GCD 外层循环中：`_hgcd(G, J, J, R)` — 输出 B（= J 缓冲区）与输入 a（= J 缓冲区）指向同一块内存。

**安全条件**：`_hgcd_recursive` 在 base case 中先执行 `A ← a`（拷贝），再执行 `B ← b`。当 B 别名 a 时，第二步覆盖 a，但 a 已在第一步拷贝到 A。递归 case 中，a 仅通过指针偏移（`a >> m`、`a mod x^m`）读取，所有读取发生在写出 A, B 之前。

**实现约束**：`_hgcd_recursive` 内部写入 A 的操作必须在任何可能覆盖 a 的操作之前完成。FLINT 的 `__set(A, a)` 在 `__set(B, b)` 之前（hgcd.c:396-397），必须保持此顺序。

### Decision 8: 重构步必须零填充

重构 `b2 = b2_lo + x^m * b3` 时，`b2_lo` 的 length 可能 < m。加高位 b3 之前，`b2[lenb2 .. m-1]` 必须清零：

```cpp
// FLINT hgcd.c:475 对应操作
memset(b2 + lenb2, 0, (m + lenb3 - lenb2) * sizeof(uint64_t));
// 然后 b2[m:] += b3
```

遗漏零填充会导致高位加到未初始化的内存上，产出错误的重构多项式。a2 的重构同理。

## 4. 模块分解

### M0: 基础操作（add / sub / divrem / _mul）

新增 `dense_upoly_zp` 的基础操作，为 HGCD 提供子操作。

#### M0.1: `_poly_add` / `_poly_sub`（raw API）

```cpp
// 逐系数模加，C 预分配 max(len_a, len_b) 元素，返回 normalize 后长度
size_t _poly_add(uint64_t* C, const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const;

// 逐系数模减，同上
size_t _poly_sub(uint64_t* C, const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const;
```

实现：短多项式用 0 补齐。add 用 `nmod_add`，sub 用 `nmod_sub`。返回前 normalize。

支持 aliasing（C 可与 A 或 B 重叠，因为是逐元素原地操作）。

#### M0.2: `_poly_divrem`（raw API）

```cpp
// 完整除法：Q 预分配 len_a - len_b + 1 元素，R 预分配 len_b - 1 元素
// 返回 {len_q, len_r}（均已 normalize）
// 使用 lazy 3-word 归约（同现有 divrem）
std::pair<size_t, size_t> _poly_divrem(
    uint64_t* Q, uint64_t* R,
    const uint64_t* A, size_t len_a,
    const uint64_t* B, size_t len_b) const;
```

HGCD 迭代 base case（M1）和递归步（M2）都需要商 Q（用于矩阵更新），因此必须是 `_poly_divrem` 而非仅 `_poly_rem`。

GCD 外层循环不需要 Q，但可直接调用 `_poly_divrem` 并忽略 Q（FLINT 也是如此，标注了 `/* todo: only rem */`）。后续可选加 `_poly_rem` 优化。

#### M0.3: `_mul`（不等长 raw API）

```cpp
// 不等长乘法，C 预分配 len_a + len_b - 1，scratch 预分配 6*max(len_a,len_b)
void _mul(uint64_t* C, const uint64_t* A, size_t len_a,
          const uint64_t* B, size_t len_b, uint64_t* scratch) const;
```

内部分派：`min(len_a, len_b) < KARATSUBA_THRESHOLD` → `_classical_mul`，否则零填充到等长后 `_kar_mul`。

#### M0.4: 对象级 `add` / `sub`（公开 API）

```cpp
static void add(dense_upoly_zp& C, const dense_upoly_zp& A, const dense_upoly_zp& B);
static void sub(dense_upoly_zp& C, const dense_upoly_zp& A, const dense_upoly_zp& B);
```

注：对象级 `rem` 不在第一版范围内。GCD 外层循环直接用 `_poly_divrem` 忽略商即可（FLINT 也是如此）。

### M1: HGCD 迭代 base case

当 `len(a) < HGCD_CUTOFF(100)` 时使用的迭代 Euclid + 矩阵追踪。

```cpp
// 迭代 base case，返回符号 sgn
// A, B, T：3 个多项式缓冲区（各 len_a 大小），通过指针旋转复用
// M：矩阵的 4 个多项式槽 + t 临时槽（各 (len_a+1)/2 大小）
// 所有 swap 都是指针交换（O(1)），不拷贝数据
int _hgcd_iter(hgcd_mat& M,
               uint64_t** A, size_t* len_a_out,
               uint64_t** B, size_t* len_b_out,
               const uint64_t* a, size_t len_a,
               const uint64_t* b, size_t len_b,
               uint64_t* Q, uint64_t** T, uint64_t** t,
               uint64_t* scratch) const;
```

**指针交换机制**（对标 FLINT `gr_ptr *`）：

FLINT 用 `gr_ptr *A`（指针的指针），swap 交换指针值而非拷贝数据。三个缓冲区 {*A, *B, *T} 通过指针旋转复用：
- divrem 后：`*T = *A mod *B`
- swap(*B, *T)：*B 指向余式，*T 指向旧 *B
- swap(*A, *T)：*A 指向旧 *B，*T 可被下次 divrem 覆盖

矩阵项 {M.poly[i]} 和临时 {*t} 同理——每步 Euclid 通过 swap 旋转，避免 O(n) 拷贝。

算法：
```
m = len_a / 2
M = identity
*A = a, *B = b
while len(*B) >= m + 1:
  divrem(Q, *T, *A, *B)
  swap(*A, *B); swap(*B, *T)            // 指针交换，O(1)
  // 更新矩阵（正确性推导见 correctness-proof.md §2.2）：
  //   M_new = M_old * [[Q, 1], [1, 0]]
  *T = Q * M[2]; *t = M[3] + *T; swap(M[3], M[2]); swap(M[2], *t)
  *T = Q * M[0]; *t = M[1] + *T; swap(M[1], M[0]); swap(M[0], *t)
  sgn = -sgn
```

### M2: HGCD 递归

核心递归算法。A, B 为输出缓冲区（各 len_a 大小），W 为工作区。

```cpp
// 递归 HGCD，返回符号 sgn
// W 为工作区（每层消耗 6*len_a + 10*(len_a+1)/2 元素，递归传递剩余）
// A, B 不需要指针的指针（递归步中通过 __set 拷贝到输出，不做指针旋转）
int _hgcd_recursive(hgcd_mat& M,
                    uint64_t* A, size_t& len_a_out,
                    uint64_t* B, size_t& len_b_out,
                    const uint64_t* a, size_t len_a,
                    const uint64_t* b, size_t len_b,
                    uint64_t* W, uint64_t* scratch) const;
```

注：与 M1 不同，递归步中 A, B 是单指针。FLINT 的 `_hgcd_recursive` 也是 `gr_ptr A`（非 `gr_ptr *A`）。递归步内部用工作区 a2, b2 等做中间计算，最终拷贝到 A, B。只有 M1 迭代 base case 需要指针的指针（因为多次 swap 避免拷贝）。

**aliasing**：支持 {B, a} 别名（见 Decision 7）。base case 中 `A ← a` 必须在 `B ← b` 之前。

算法（见调研文档 §3.2 Layer 2，正确性推导见 correctness-proof.md §3）：

1. `m = len_a / 2`
2. 若 `len_b < m + 1`：`A ← a`，然后 `B ← b`（顺序不可交换），M = identity，返回
3. 取高位 `a0 = a >> m`，`b0 = b >> m`（指针偏移，零拷贝）
4. 若 `len(a0) < HGCD_CUTOFF`：`_hgcd_iter(R, a3, b3, a0, b0, ...)`
   否则：`_hgcd_recursive(R, a3, b3, a0, b0, W', ...)`
5. 重构完整 `a2`, `b2`（公式见 correctness-proof.md §3.3）：
   - sgnR > 0: `b2_lo = R[0]*b_lo - R[2]*a_lo`
   - sgnR < 0: `b2_lo = R[2]*a_lo - R[0]*b_lo`
   - **零填充 `b2[len(b2_lo) .. m-1] = 0`**（Decision 8）
   - `b2[m:] += b3`
   - 类似重构 `a2`（同样需零填充）
6. 若 `len(b2) < m + 1`：M = R，A = a2，B = b2，返回
7. `_poly_divrem(q, d, a2, b2)`（需要商 q 用于步骤 11 的矩阵合并）
8. 取新高位 `c0 = b2 >> k`，`d0 = d >> k`（k = 2m - len(b2) + 1）
9. 递归或迭代得矩阵 S
10. 重构 A, B（同步骤 5，用 S 和 sgnS）
11. 合并矩阵：`S_mod = [[q,1],[1,0]] * S`，然后 `M = R * S_mod`（用 `_mat_mul`）
12. 符号：`sgn = -(sgnR * sgnS)`

### M3: 矩阵乘法

```cpp
// classical 2x2 矩阵乘法：C = A * B
// T 为临时空间，至少 max_entry_len(A) + max_entry_len(B) - 1 元素
void _mat_mul(hgcd_mat& C, const hgcd_mat& A, const hgcd_mat& B,
              uint64_t* T, uint64_t* scratch) const;
```

Classical 实现（8 次 `_mul` + 4 次 `_poly_add`）：
```
C[0] = A[0]*B[0] + A[1]*B[2]
C[1] = A[0]*B[1] + A[1]*B[3]
C[2] = A[2]*B[0] + A[3]*B[2]
C[3] = A[2]*B[1] + A[3]*B[3]
```

### M4: GCD 分派与 HGCD-based GCD

```cpp
// HGCD-based GCD 循环
// 前置：len_a >= len_b >= GCD_CUTOFF
void _gcd_hgcd(uint64_t* G, size_t& len_g,
               const uint64_t* A, size_t len_a,
               const uint64_t* B, size_t len_b) const;
```

**完整流程**（对照 FLINT gcd_hgcd.c:54-98，三个零检查缺一不可）：

```
_gcd_hgcd(G, A, B):
    // 分配工作区：J[n], R[n], Q[n] + HGCD 工作区 + scratch
    _poly_divrem(Q, R, A, B)
    if len(R) == 0:                    // 检查 ①：B 整除 A
        G ← B; return
    _hgcd_recursive(_, G, J, B, R, W)  // 首次 HGCD 缩减 (B, R) → (G, J)
    while len(J) != 0:                 // 检查 ②：J==0 表示 GCD 已在 HGCD 内找到
        _poly_divrem(Q, R, G, J)
        if len(R) == 0:                // 检查 ③：J 整除 G
            G ← J; break
        if len(J) < GCD_CUTOFF:
            _gcd_euclid(G, J, R)       // fallback Euclid
            break
        _hgcd_recursive(_, G, J, J, R, W)  // aliasing: 输出 J = 输入 a (Decision 7)
```

修改 `dense_upoly_zp::gcd`，加入分派：

```cpp
static void gcd(dense_upoly_zp& G,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    // ... 规范化：确保 deg(a) >= deg(b)
    if (len_b < GCD_CUTOFF) {
        // 现有 Euclid
    } else {
        // _gcd_hgcd
    }
}
```

## 5. 工作区设计

### 5.1 HGCD 递归工作区

每层递归消耗的空间（FLINT 同款布局）：

```
a2[n]  b2[n]  a3[n]  b3[n]  q[(n+1)/2]  d[n]  T0[n]  T1[(n+1)/2]
R[0..3] 各 (n+1)/2    S[0..3] 各 (n+1)/2
```

每层：`6n + 10 * (n+1)/2` 元素。

总工作区（含递归深度 O(log n)）：
```
W_hgcd = 22 * n + 16 * (ceil_log2(n) + 1)
```

### 5.2 GCD 外层工作区

`_gcd_hgcd` 循环额外需要：

```
A_buf[n]  B_buf[n]  R_buf[n]  Q_buf[n]  scratch[6n]
```

约 `10n` 元素。

### 5.3 总分配

```
total = 32 * n + 16 * ceil_log2(n)    (约 32n，对 n=10000 约 2.5MB)
```

一次 `std::vector<uint64_t>` 分配，整个 GCD 过程零额外 heap 操作。

## 6. 文件改动

| 文件 | 改动 | 估计行数 |
|------|------|---------|
| `clpoly/dense_upoly_zp.hh` | M0-M4 全部新增 | +400~500 行 |

所有改动集中在 `dense_upoly_zp.hh` 一个文件中，保持 header-only 风格（与现有 Karatsuba/divrem 一致）。

无需修改其他文件：`gcd` 签名不变，`polynomial_gcd.hh` 集成代码不变。

## 7. 实施顺序

```
M0: 基础操作 (add/sub/divrem/_mul)
 |  独立可测试
 ↓
M1: _hgcd_iter (迭代 base case)
 |  可用小 case 测试矩阵追踪正确性
 ↓
M2: _hgcd_recursive
 |  组合 M1，可测试 HGCD 缩减正确性
 ↓
M3: _mat_mul (classical)
 |  M2 依赖此操作
 ↓
M4: _gcd_hgcd + gcd 分派
 |  组合 M0-M3，全量回归测试
 ↓
性能验证: bench_clpoly + bench_comparative
```

注：M2 和 M3 有循环依赖（递归中用矩阵乘法，矩阵乘法用多项式乘法），实际实现时 M3 先于 M2 完成，但测试可一起进行。

## 8. 测试方案

1. **M0 单元测试**：
   - add/sub：随机多项式对比 `nmod_add`/`nmod_sub` 逐系数结果
   - `_poly_divrem`：验证 A = Q*B + R 且 deg(R) < deg(B)
   - `_mul`：不等长乘法对比现有 `mul` 结果

2. **M1-M3 集成测试**：
   - HGCD 缩减正确性：验证 `M * (A, B)^T = (a, b)^T`（M 映射缩减后 → 原始）
   - 矩阵行列式：`det(M) = (-1)^k`

3. **M4 全量回归**：
   - `bash test/run_all_tests.sh`（全部 279 单元测试）
   - `make crosscheck`（258 FLINT 交叉验证）
   - 高次专项：deg500, 1000, 3000, 5000, 10000 的 GCD 正确性

4. **性能验证**：
   - `make bench-clpoly` + `make bench-comparative`
   - 关注 deg1000+ 的 ratio 变化

## 9. 预期收益

| 度数 | 当前 ratio vs FLINT | HGCD 后预估 | 说明 |
|------|--------------------|-----------|----|
| deg200 | 1.25x | ~1.25x | 低于 cutoff，不变 |
| deg500 | 1.44x | ~1.2-1.3x | 刚过 cutoff |
| deg1000 | 2.67x | ~1.3-1.5x | HGCD 主力区间 |
| deg3000 | 4.28x | ~1.5-2x | 同上 |
| deg5000 | 5.28x | ~1.5-2x | 同上 |
| deg10000 | 6.87x | ~2-2.5x | 剩余差距来自 mul 层级（Karatsuba vs KS4/NTT） |

HGCD 后的剩余差距主要来自乘法算法层级：CLPoly Karatsuba O(n^1.585) vs FLINT KS4→GMP O(n^~1.4)。需 P3-KS（Kronecker Substitution → GMP mpn_mul）进一步消除。
