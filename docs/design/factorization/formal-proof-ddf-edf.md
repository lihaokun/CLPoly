# DDF/EDF 形式化数学论证

## 1. 符号约定

- $\mathbb{F}_p = \mathbb{Z}/p\mathbb{Z}$，$p$ 为素数
- $\mathbb{F}_{p^d}$ = 含 $p^d$ 个元素的有限域 = $\mathbb{F}_p[x]/(g)$ 对任意 $d$ 次不可约 $g$
- $f \in \mathbb{F}_p[x]$ 为首一无平方多项式
- $f = f_1 f_2 \cdots f_k$ 为不可约分解

## 2. 核心定理

**定理 2.1 (有限域的多项式分解)**
$$x^{p^d} - x = \prod_{\substack{g \in \mathbb{F}_p[x] \\ g \text{ 首一不可约} \\ \deg(g) \mid d}} g$$

**证明**：$x^{p^d} - x$ 的根恰为 $\mathbb{F}_{p^d}$ 的全部 $p^d$ 个元素。每个 $\alpha \in \mathbb{F}_{p^d}$ 的极小多项式 $m_\alpha(x)$ 不可约且 $\deg(m_\alpha) \mid d$。反之，若 $g$ 不可约且 $\deg(g) = k \mid d$，则 $\mathbb{F}_{p^k} \hookrightarrow \mathbb{F}_{p^d}$，$g$ 的所有根在 $\mathbb{F}_{p^d}$ 中，因此 $g \mid (x^{p^d} - x)$。

**推论 2.2**
$$\gcd(x^{p^d} - x,\; f) = \prod_{\substack{f_i \mid f \\ f_i \text{ 不可约} \\ \deg(f_i) \mid d}} f_i$$

**证明**：由定理 2.1 和 $f$ 无平方直接得到。

## 3. DDF 算法正确性

### 3.1 算法描述

对应 `polynomial_factorize_zp.hh` 第 247-292 行 `__ddf_Zp`：

```
输入: 首一无平方 f ∈ F_p[x]
输出: {(g_d, d)} 其中 g_d = ∏{f 的 d 次不可约因子}

h ← x,  f* ← f
for d = 1, 2, ... :
    if deg(f*) < 2d: break
    h ← h^p mod f*
    g_d ← gcd(h - x, f*)
    if deg(g_d) > 0:
        output (g_d, d)
        f* ← f* / g_d
        h ← h mod f*
if deg(f*) > 0: output (f*, deg(f*))
```

### 3.2 循环不变量

**引理 3.1 (h 不变量)**：在第 $d$ 次迭代开始时（计算 $h \leftarrow h^p$ 之前），$h \equiv x^{p^{d-1}} \pmod{f^*}$。计算后，$h \equiv x^{p^d} \pmod{f^*}$。

**证明**（归纳法）：
- **基础**：$d = 1$ 前，$h = x = x^{p^0}$。
- **归纳**：假设第 $d$ 次开始前 $h \equiv x^{p^{d-1}} \pmod{f^*}$，则 $h^p \equiv x^{p^d} \pmod{f^*}$。

**关键**：提取 $g_d$ 后更新 $f^*_{\text{new}} = f^* / g_d$，$h \leftarrow h \bmod f^*_{\text{new}}$。

由于 $f^*_{\text{new}} \mid f^*_{\text{old}}$，$h \equiv x^{p^d} \pmod{f^*_{\text{old}}}$ 蕴含 $h \equiv x^{p^d} \pmod{f^*_{\text{new}}}$。显式取模仅减小代表元次数。不变量成立。$\square$

### 3.3 分裂正确性

**定理 3.2**：在第 $d$ 次迭代中，$g_d = \gcd(h - x, f^*)$ 恰为 $f^*$ 中所有 $d$ 次不可约因子之积。

**证明**：
1. 由引理 3.1，$h - x \equiv x^{p^d} - x \pmod{f^*}$。
2. 因此 $\gcd(h - x, f^*) = \gcd(x^{p^d} - x, f^*)$。
3. 由推论 2.2，$\gcd(x^{p^d} - x, f^*) = \prod\{f_i \mid f^* : \deg(f_i) \mid d\}$。
4. 当前 $f^*$ 仅包含度 $\geq d$ 的不可约因子（度 $< d$ 的因子在先前迭代已提取）。
5. 故满足 $\deg(f_i) \mid d$ 且 $\deg(f_i) \geq d$ 的因子，恰为 $\deg(f_i) = d$ 的因子。
6. 因此 $g_d$ 恰为 $f^*$ 中所有 $d$ 次不可约因子之积。$\square$

### 3.4 提前终止正确性

**命题 3.3**：当 $\deg(f^*) < 2d$ 时，$f^*$ 至多有一个不可约因子。

**证明**：$f^*$ 无平方，其所有不可约因子度 $\geq d$。若有两个以上不可约因子，则 $\deg(f^*) \geq 2d$，矛盾。

若 $\deg(f^*) > 0$，则 $f^*$ 本身就是唯一的不可约因子，度为 $\deg(f^*)$。输出 $(f^*, \deg(f^*))$ 正确。$\square$

### 3.5 终止性

**命题 3.4**：算法必然终止。

**证明**：$d$ 单调递增。无论 $g_d$ 是否非平凡，$2d$ 终将超过 $\deg(f^*)$（后者有界），循环退出。$\square$

### 3.6 代码实现对照

| 代码行 | 数学操作 | 正确性依据 |
|--------|---------|-----------|
| L264: `h = __upoly_powmod(h, ZZ(p), f_star)` | $h \leftarrow h^p \bmod f^*$ | 引理 3.1 |
| L267: `__upoly_subtract_x(h, p)` | $h - x$ | 需要 $-1 \equiv p-1 \pmod{p}$ 正确表示 |
| L269: `polynomial_GCD(h_minus_x, f_star)` | $\gcd(h-x, f^*)$ | 定理 3.2 |
| L277-278: `pair_vec_div(f_new, f_star, gd, ...)` | $f^* \leftarrow f^*/g_d$ | 精确整除（$g_d \mid f^*$） |
| L281: `h = __upoly_mod(h, f_star)` | 减小 $h$ 的次数 | 引理 3.1 中 div 不变量 |
| L260: `deg(f*) < 2d` → break | 提前终止 | 命题 3.3 |
| L288: output `(f_star, deg(f_star))` | 残余单不可约因子 | 命题 3.3 |

## 4. EDF (Cantor-Zassenhaus) 正确性

### 4.1 算法描述

对应 `polynomial_factorize_zp.hh` 第 293-353 行 `__edf_Zp`：

```
输入: 首一无平方 f ∈ F_p[x]，所有不可约因子度为 d，p 为奇素数
输出: f 的全部不可约因子

if deg(f) = d: return {f}
repeat:
    r ← random polynomial, deg(r) < deg(f)
    g ← gcd(r^{(p^d-1)/2} - 1, f)
    if 0 < deg(g) < deg(f):
        递归 edf(g, d) 和 edf(f/g, d)
        return
```

### 4.2 CRT 分解

**引理 4.1 (中国剩余定理)**：设 $f = f_1 \cdots f_k$，$f_i$ 两两互素不可约。则
$$\mathbb{F}_p[x]/(f) \cong \mathbb{F}_{p^d} \times \cdots \times \mathbb{F}_{p^d} \quad (k \text{ 份})$$
其中同构将 $r \bmod f$ 映为 $(r \bmod f_1, \ldots, r \bmod f_k)$。

### 4.3 Legendre 符号论证

**定理 4.2**：对 $r \in \mathbb{F}_p[x]/(f)$，记 $r_i = r \bmod f_i \in \mathbb{F}_{p^d}$。当 $r_i \neq 0$ 时，
$$r_i^{(p^d-1)/2} \in \{1, -1\}$$
（Fermat 小定理在 $\mathbb{F}_{p^d}^*$ 中的推论）。

**证明**：$\mathbb{F}_{p^d}^* = \langle \alpha \rangle$ 是 $p^d - 1$ 阶循环群。$r_i^{p^d-1} = 1$，故 $r_i^{(p^d-1)/2}$ 是 $1$ 的平方根，即 $\pm 1$。$\square$

**推论 4.3**：$r^{(p^d-1)/2} \bmod f$ 的 CRT 分量各为 $\pm 1$（或 $0$ 若 $r_i = 0$，概率为 $1/p^d$，可忽略）。

### 4.4 分裂概率

**定理 4.4**：设 $f$ 有 $k \geq 2$ 个不可约因子。对随机 $r$，$\gcd(r^{(p^d-1)/2} - 1, f)$ 为非平凡因子的概率 $\geq 1 - 2^{1-k} \geq 1/2$。

**证明**：
1. $g = \gcd(r^{(p^d-1)/2} - 1, f) = \prod_{i : r_i^{(p^d-1)/2} = 1} f_i$
2. $g$ 平凡 $\Leftrightarrow$ 全部 $r_i^{(p^d-1)/2}$ 同号：全 $+1$ 或全 $-1$
3. 每个 $r_i^{(p^d-1)/2}$ 独立取 $\pm 1$，各概率 $1/2$（对随机非零 $r_i$）
4. 全同号概率 $= 2 \cdot (1/2)^k = 2^{1-k}$
5. 非平凡分裂概率 $= 1 - 2^{1-k} \geq 1/2$（$k \geq 2$）

期望尝试次数：$k = 2$ 时 $\leq 2$ 次，$k \geq 3$ 时 $\leq 4/3$ 次。$\square$

### 4.5 代码实现对照

| 代码行 | 数学操作 | 正确性依据 |
|--------|---------|-----------|
| L334: `ZZ exp = (pow(ZZ(p), d) - 1) / 2` | $(p^d-1)/2$ | $p$ 奇数 ⇒ $p^d-1$ 偶数 |
| L335: `__upoly_powmod(r, exp, f)` | $r^{(p^d-1)/2} \bmod f$ | 二进制快速幂 + 多项式取模 |
| L336: `__upoly_subtract_one(g_pow, p)` | $g - 1$ | **需 $-1 \equiv p-1$ 正确表示** |
| L337: `polynomial_GCD(g_pow_minus_1, f)` | $\gcd(g-1, f)$ | 定理 4.4 |
| L340: `deg(g) > 0 && deg(g) < deg(f)` | 非平凡分裂检测 | — |
| L344: `pair_vec_div(h_part, f, g, ...)` | $f / g$ | 精确整除 |
| L348-349: 递归 `__edf_Zp(g, d)` 和 `__edf_Zp(h_part, d)` | 分治 | 每次分裂减小次数 |

## 5. 已修复的 Bug 及其数学根因

### 5.1 Bug 描述

`__upoly_subtract_x` (L193, L209) 和 `__upoly_subtract_one` (L239) 中：
```cpp
Zp neg_one = __make_zp((int64_t)(p - 1), p);
```

### 5.2 数学根因

需要计算 $-1 \equiv p - 1 \pmod{p}$。当 $p > 2^{63}$ 时（例如 $p = 2^{64} - 59$）：

$$p - 1 = 2^{64} - 60 > 2^{63} - 1 = \texttt{INT64\_MAX}$$

强制转换 `(int64_t)(p - 1)` 产生补码回绕：$(2^{64} - 60) \bmod 2^{64}$ 作为有符号整数 $= -60$。

构造函数 `Zp(-60, p)` 计算 $|-{60}| \bmod p = 60$，再取负 $p - 60$。结果为 $p - 60 \neq p - 1$。

### 5.3 影响分析

- **DDF**：$\gcd(h - x, f^*)$ 中 $h - x$ 计算错误 → 因子分组错误
- **EDF**：$\gcd(g - 1, f)$ 中 $g - 1$ 计算错误 → 分裂始终失败 → 死循环

### 5.4 修复

将 `__make_zp((int64_t)(p - 1), p)` 替换为 `Zp(p - 1, p)`。

`Zp(uint64_t i, uint64_t p)` 构造函数执行 `_i = i % p`：
- $i = p - 1$，$i \% p = p - 1$
- 对任何 $p \leq 2^{64} - 1$ 均正确

### 5.5 Bug 触发条件

仅当同时满足以下条件时触发：
1. $p > 2^{63}$（大素数模式）
2. `__upoly_subtract_x/one` 进入"无对应项"分支（即 $h$ 无 $x^1$ 项或无 $x^0$ 项）

条件 2 在数学上等价于：Legendre 符号分组结果的某些系数恰为零。这对随机多项式以约 $1/p$ 的概率发生，但在 EDF 的 Cantor-Zassenhaus 中（恰好需要分裂成功的随机元素），零常数项出现概率约 $1/2$（当两个根的 Legendre 符号相反时）。因此该 bug 在 EDF 中几乎必然触发。
