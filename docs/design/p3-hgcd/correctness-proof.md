# P3-HGCD 正确性推导

> 本文档对架构文档中的 HGCD 算法进行数学推导，验证每一层的正确性。
> 对照 FLINT 3.0.1 源码 `src/gr_poly/hgcd.c`。

---

## 1. 记号与约定

### 1.1 基本记号

- 多项式 f 的 **length** = deg(f) + 1（FLINT 惯例，CLPoly `_coeffs.size()`）
- `len(f) = 0` 表示零多项式
- 升幂排列：`f = f[0] + f[1]x + ... + f[n-1]x^{n-1}`，length = n

### 1.2 多项式分割

对多项式 f（length = n）和位置 m：
- **高位**（shift）：`f_hi = f >> m`，即 `f_hi[i] = f[m+i]`，length = max(n - m, 0)
- **低位**（truncate）：`f_lo = f mod x^m`，即 `f_lo[i] = f[i]`（i < m），length = min(n, m)
- 恒等式：**f = f_lo + x^m * f_hi**

### 1.3 矩阵约定

2x2 矩阵 M = [[M[0], M[1]], [M[2], M[3]]]，索引按 FLINT 布局。

**FLINT 的矩阵约定**（关键！）：

> **M 映射缩减后 → 原始：M * (A, B)^T = (a, b)^T**

即 M 是*逆变换*，将 HGCD 输出的缩减多项式映射回原始输入。

这与直觉的 "正向变换" 相反，但简化了递归中的矩阵合并。

### 1.4 Euclid 步的矩阵表示

一步 Euclid：`A_new = B, B_new = A - Q*B`（其中 Q = A div B）。

正向变换：`(A_new, B_new) = [[0, 1], [1, -Q]] * (A, B)`

**逆变换**：`(A, B) = [[Q, 1], [1, 0]] * (A_new, B_new)`

验证：`[[Q, 1], [1, 0]] * (B, A-QB) = (QB + A-QB, B) = (A, B)` ✓

det([[Q, 1], [1, 0]]) = -1，故 k 步 Euclid 后 det(M) = (-1)^k。

---

## 2. 迭代 base case 正确性（Layer 3）

**函数**：`_hgcd_iter(M, A, B, a, b)`

**输入**：多项式 a（length = n），b（length < n），m = n/2

**输出**：矩阵 M，多项式 A, B，符号 sgn

**不变量**：每次循环迭代后，`M * (A, B)^T = (a, b)^T` 且 `det(M) = sgn`

### 2.1 初始化

```
M = I（单位矩阵），A = a，B = b，sgn = 1
```

- `I * (a, b)^T = (a, b)^T` ✓
- `det(I) = 1 = sgn` ✓

### 2.2 循环体

**前置条件**：`M_old * (A_old, B_old)^T = (a, b)^T`，`det(M_old) = sgn_old`

**操作**（FLINT hgcd.c:309-354）：
```
divrem(Q, T, A, B)              // T = A mod B
swap(B, T); swap(A, old_value)  // A_new = B_old, B_new = T = A_old mod B_old
```

Euclid 步矩阵 E = [[Q, 1], [1, 0]]：
```
E * (A_new, B_new) = E * (B_old, A_old - Q*B_old)
                   = (Q*B_old + A_old - Q*B_old, B_old)
                   = (A_old, B_old)
```

**矩阵更新**（FLINT hgcd.c:344-352）：

FLINT 代码（忽略 res 逻辑）：
```
T = Q * M[2];  t = M[3] + T;  swap(M[3], M[2]);  swap(M[2], t)
T = Q * M[0];  t = M[1] + T;  swap(M[1], M[0]);  swap(M[0], t)
```

逐步追踪（行 2 为例）：
1. `T = Q * M[2]_old`
2. `t = M[3]_old + Q * M[2]_old`
3. swap(M[3], M[2]): `M[3] = M[2]_old`, `M[2] = M[3]_old`
4. swap(M[2], t): `M[2] = t = M[3]_old + Q * M[2]_old`

最终：
```
M_new = [[M[1]_old + Q*M[0]_old,  M[0]_old],
         [M[3]_old + Q*M[2]_old,  M[2]_old]]
```

**验证 M_new = M_old * E**：
```
M_old * E = [[M[0], M[1]], [M[2], M[3]]] * [[Q, 1], [1, 0]]
          = [[M[0]*Q + M[1], M[0]],
             [M[2]*Q + M[3], M[2]]]
```

与 M_new 完全一致 ✓

**后置条件**：
```
M_new * (A_new, B_new) = (M_old * E) * (A_new, B_new)
                       = M_old * (E * (A_new, B_new))
                       = M_old * (A_old, B_old)
                       = (a, b)  ✓
```

**行列式**：`det(M_new) = det(M_old) * det(E) = sgn_old * (-1) = -sgn_old`

FLINT 用 `sgn = -sgn`（hgcd.c:354），一致 ✓

### 2.3 终止条件

循环在 `len(B) < m + 1` 时终止。

**输出保证**：`len(A) >= m + 1 > len(B)`（因为 A 来自上一轮的 B，而上一轮 B 满足 `len >= m+1`）

### 2.4 终止性

每步 divrem 严格减少 len(B)（因为 len(T) = len(A mod B) < len(B)），所以循环必终止。

---

## 3. 递归步正确性（Layer 2）

**函数**：`_hgcd_recursive(M, A, B, a, b, W)`

**输入**：a（length = n），b（length <= n），m = n/2

**输出**：M, A, B, sgn 满足 `M * (A, B)^T = (a, b)^T`

### 3.1 Base case（len(b) < m + 1）

```
M = I, A = a, B = b, sgn = 1
```

`I * (a, b) = (a, b)` ✓。但注意此时 HGCD 不保证 `len(A) >= m+1 > len(B)`——因为 b 本身就不够长，没有 Euclid 步可执行。

### 3.2 第一次递归

**分割**（FLINT hgcd.c:436-437）：
```
a0 = a >> m    (a 的高位部分，length = n - m)
b0 = b >> m    (b 的高位部分，length = len(b) - m)
```

**递归**：`HGCD(a0, b0) → R, a3, b3, sgnR`

**递归后不变量**：`R * (a3, b3)^T = (a0, b0)^T`

### 3.3 重构公式推导

已知 `R * (a3, b3)^T = (a0, b0)^T`，即 `(a3, b3)^T = R^{-1} * (a0, b0)^T`。

我们要求 `(a2, b2)` 使得 `R * (a2, b2)^T = (a, b)^T`。

由 `a = a_lo + x^m * a_hi = a_lo + x^m * a0`（a_lo = a mod x^m）：

```
(a, b)^T = (a_lo, b_lo)^T + x^m * (a0, b0)^T
         = (a_lo, b_lo)^T + x^m * R * (a3, b3)^T
```

所以：
```
R * (a2, b2)^T = (a_lo, b_lo)^T + x^m * R * (a3, b3)^T
(a2, b2)^T = R^{-1} * (a_lo, b_lo)^T + x^m * (a3, b3)^T
```

设 R = [[R0, R1], [R2, R3]]，det(R) = sgnR。

R^{-1} = (1/sgnR) * [[R3, -R1], [-R2, R0]]

所以：
```
a2 = (1/sgnR) * (R3 * a_lo - R1 * b_lo) + x^m * a3
b2 = (1/sgnR) * (-R2 * a_lo + R0 * b_lo) + x^m * b3
```

分两种情况：

**sgnR = +1**（偶数步）：
```
a2_lo = R3 * a_lo - R1 * b_lo      (低 m 项)
b2_lo = R0 * b_lo - R2 * a_lo      (低 m 项)
a2 = a2_lo + x^m * a3
b2 = b2_lo + x^m * b3
```

**sgnR = -1**（奇数步）：
```
a2_lo = R1 * b_lo - R3 * a_lo
b2_lo = R2 * a_lo - R0 * b_lo
a2 = a2_lo + x^m * a3
b2 = b2_lo + x^m * b3
```

### 3.4 FLINT 代码验证

FLINT 计算 b2（hgcd.c:467-481）：
```c
__mul(b2, R[2], s);              // b2 = R2 * a_lo
__mul(T0, R[0], t);              // T0 = R0 * b_lo
if (sgnR < 0)
    __sub(b2, b2, T0);           // b2 = R2*a_lo - R0*b_lo
else
    __sub(b2, T0, b2);           // b2 = R0*b_lo - R2*a_lo
// 高位：b2[m:] += b3
```

- sgnR > 0: `b2_lo = R0*b_lo - R2*a_lo` ✓（与公式一致）
- sgnR < 0: `b2_lo = R2*a_lo - R0*b_lo` ✓（与公式一致）

FLINT 计算 a2（hgcd.c:483-496）：
```c
__mul(a2, R[3], s);              // a2 = R3 * a_lo
__mul(T0, R[1], t);              // T0 = R1 * b_lo
if (sgnR < 0)
    __sub(a2, T0, a2);           // a2 = R1*b_lo - R3*a_lo
else
    __sub(a2, a2, T0);           // a2 = R3*a_lo - R1*b_lo
// 高位：a2[m:] += a3
```

- sgnR > 0: `a2_lo = R3*a_lo - R1*b_lo` ✓
- sgnR < 0: `a2_lo = R1*b_lo - R3*a_lo` ✓

**重构后的验证**：`R * (a2, b2)^T = (a, b)^T` ✓

### 3.5 提前终止（len(b2) < m + 1）

FLINT hgcd.c:498-512：若 `len(b2) < m + 1`，则 M = R，A = a2，B = b2。

- `M * (A, B) = R * (a2, b2) = (a, b)` ✓
- `sgn = sgnR` ✓

### 3.6 中间 divrem 步

若 `len(b2) >= m + 1`，需继续缩减。

FLINT hgcd.c:552：
```c
__divrem(q, d, a2, b2)      // d = a2 mod b2, q = a2 div b2
```

Euclid 步矩阵 E = [[q, 1], [1, 0]]：
```
E * (b2, d)^T = (q*b2 + d, b2) = (a2, b2)
```

合并到此：
```
R * E * (b2, d)^T = R * (a2, b2)^T = (a, b)^T
```

### 3.7 第二次递归

**新的分割位置**：`k = 2m - len(b2) + 1`

```
c0 = b2 >> k     (b2 的高位)
d0 = d >> k      (d 的高位)
```

**递归**：`HGCD(c0, d0) → S, a3', b3', sgnS`

**递归后不变量**：`S * (a3', b3')^T = (c0, d0)^T`

### 3.8 第二次重构

与第一次重构完全类似，将 S^{-1} 应用到 (b2, d) 的低 k 项，加上 x^k * (a3', b3')，得到最终的 (A, B)。

```
S * (A, B)^T = (b2, d)^T
```

### 3.9 矩阵合并

FLINT hgcd.c:635-644：

```c
// 将 q 合并进 S：S_modified = [[q,1],[1,0]] * S
__swap(S[0], S[2]); __swap(S[1], S[3]);
S[0] += q * S[2];   S[1] += q * S[3];

// 最终矩阵：M = R * S_modified
__mat_mul(M, R, S_modified);
```

**验证**：
```
M * (A, B)^T = R * [[q,1],[1,0]] * S * (A, B)^T
             = R * [[q,1],[1,0]] * (b2, d)^T
             = R * (q*b2 + d, b2)^T
             = R * (a2, b2)^T
             = (a, b)^T  ✓
```

### 3.10 符号合并

FLINT hgcd.c:647：`sgn = -(sgnR * sgnS)`

- R 贡献 sgnR（其 det）
- 中间 divrem 步贡献 -1（一步 Euclid 的 det）
- S 贡献 sgnS

总 det = sgnR * (-1) * sgnS = -(sgnR * sgnS) ✓

---

## 4. GCD 外层循环正确性（Layer 1）

**函数**：`_gcd_hgcd(G, A, B)`

**FLINT 实现**（gcd_hgcd.c:54-105，简化）：

```
divrem(Q, R, A, B)            // 初始 A mod B
if R == 0: G = B; return

hgcd(_, _, G, J, B, R)        // 缩减 (B, R) → (G, J)
while J != 0:
    divrem(Q, R, G, J)
    if R == 0: G = J; break
    if len(J) < cutoff:
        gcd_euclidean(G, J, R)
        break
    hgcd(_, _, G, J, J, R)    // 继续缩减
```

### 4.1 GCD 保持

**引理**：若 `M * (A, B)^T = (a, b)^T` 且 det(M) = +/-1，则 `gcd(A, B) = gcd(a, b)`。

**证明**：M 可逆（det = +/-1，在 Zp[x] 中系数可逆），所以 (a, b) 和 (A, B) 互相是线性组合。任何整除 (a, b) 的多项式必整除 (A, B)，反之亦然。□

每步操作：
- `rem(R, A, B)`: gcd(A, B) = gcd(B, R) ✓
- `hgcd(G, J, A, B)`: gcd(G, J) = gcd(A, B)（由引理）✓
- `gcd_euclidean`: 标准 Euclid 正确性 ✓

### 4.2 终止性

HGCD 保证缩减：输出 len(B) 严格小于输入 len(b)（至少缩减到 m+1 以下）。外层循环每轮 HGCD + divrem 至少减少 len，最终 len < cutoff 时 fallback Euclid。

### 4.3 不需要矩阵 M

GCD 外层循环只需要缩减后的 (A, B)，不需要矩阵 M。FLINT 传 NULL 给 M 参数。我们的实现同理——HGCD 递归内部需要 M 做重构，但 `_gcd_hgcd` 入口不需要。

---

## 5. 边界条件与特殊情况

### 5.1 m 的定义

`m = len(a) / 2`（整数除法，向下取整）。

- len(a) = 1: m = 0, 条件 `len(b) >= m+1 = 1` 几乎总成立（b 非零即成立），但 divrem(a, b) 对常数多项式就是 R = 0，循环立即终止。
- len(a) = 2: m = 1, 需要 len(b) >= 2，即 b 至少线性。

### 5.2 k 的范围

第二次递归的移位量 `k = 2m - len(b2) + 1`。

由条件 `len(b2) >= m + 1`（否则已提前终止）：
```
k = 2m - len(b2) + 1 <= 2m - (m+1) + 1 = m
```

由 `len(b2) <= len(a2) < len(a) = n`（divrem 后长度递减）且 `m = n/2`：
```
k = 2m - len(b2) + 1 >= 2m - n + 1 = 2*(n/2) - n + 1
```

对偶数 n: k >= 1；对奇数 n: k >= 0。但 FLINT 的 `lena / 2` 对奇数 n 给出 m = (n-1)/2，所以 k >= 2*((n-1)/2) - n + 1 = n - 1 - n + 1 = 0。

实际上 k > 0 总成立，因为 `len(b2) <= n - 1`（b2 是 HGCD 缩减后的 B，比 a 短）。

### 5.3 零多项式

HGCD 内部**可以**遇到零多项式，但**不需要额外的零检查**——现有控制流自然覆盖：

1. **迭代 base case**：divrem 余式 R=0 时，swap 后 B=0（length=0）。`while (len(B) >= m+1)` 条件排除 B=0，循环自然退出。
2. **递归步内部**：`divrem(q, d, a2, b2)` 若 b2 整除 a2，则 d=0，第二次递归收到 d0 = empty（length=0）。base case `len(b) < m+1` 包含 b=0，返回 identity 矩阵。
3. **GCD 外层循环**：HGCD 返回 B=0 时，`while (J != 0)` 退出循环，G 即为 GCD。`R == 0` 检查处理 divrem 余式为零的情况。

三处均由现有条件判断自然处理，无需单独的零多项式分支。

### 5.4 工作区大小

每层递归消耗 `6*n + 10*(n+1)/2` 元素：
```
a2[n] + b2[n] + a3[n] + b3[n] + q[(n+1)/2] + d[n] + T0[n] + T1[(n+1)/2]
+ R[0..3] 各 (n+1)/2 + S[0..3] 各 (n+1)/2
= 6n + 10*(n+1)/2
```

递归深度 = O(log n)（每层 n 减半），总计：
```
sum_{i=0}^{log n} (6*(n/2^i) + 10*(n/2^i + 1)/2)
= 12n + 10*log(n) + ...
≈ 22n + 16*(ceil_log2(n) + 1)     (FLINT 的精确公式)
```

---

## 6. 正确性小结

| 性质 | 保证 | 依据 |
|------|------|------|
| `M * (A, B) = (a, b)` | 每层递归/迭代的不变量 | §2.2, §3.3-3.9 |
| `det(M) = sgn = (-1)^k` | Euclid 步矩阵 det = -1 的累积 | §2.2, §3.10 |
| `gcd(A, B) = gcd(a, b)` | M 可逆 → 线性组合保持 GCD | §4.1 |
| `len(A) >= m+1 > len(B)` | 迭代终止条件 | §2.3 |
| 终止性 | 每步严格缩减 len(B) | §2.4, §4.2 |
| 工作区充足 | 逐层计算 + 递归深度 O(log n) | §5.4 |

### 对我们实现的要求

1. **矩阵约定**：必须使用 FLINT 同款约定 `M * (A, B) = (a, b)`，不可反过来
2. **符号追踪**：每步 Euclid 翻转 sgn，两次递归的符号用 `-(sgnR * sgnS)` 合并
3. **重构公式**：低位部分的符号取决于 sgnR/sgnS，必须正确分支
4. **零填充**：重构时低位部分可能短于 m/k，必须零填充到合并位置
5. **k 计算**：`k = 2*m - len(b2) + 1`，不可用其他公式

这些都是"差一个符号就全错"的硬约束，实现时必须逐条对照验证。
