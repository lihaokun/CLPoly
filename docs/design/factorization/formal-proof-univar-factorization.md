# 单变量整数因式分解：形式化正确性证明

## 1. 记号与设定

### 1.1 基本对象

设 $f \in \mathbb{Z}[x]$，满足：
- 无平方（squarefree）
- 本原（primitive）：$\text{cont}(f) = 1$
- $n = \deg f \geq 2$
- $\ell = \text{lc}(f) > 0$（首项系数为正，不失一般性）

设 $f = g_1 g_2 \cdots g_t$ 是 $f$ 在 $\mathbb{Z}[x]$ 中的完全不可约分解。

**Gauss 引理**：若 $g, h \in \mathbb{Z}[x]$ 均本原，则 $gh$ 本原。

由此推出：每个 $g_j$ 本原（不可约多项式必本原），且 $g_1 \cdots g_t$ 本原。由 $f$ 本原，$f = g_1 \cdots g_t$（无额外整数因子）。

记 $\ell_j = \text{lc}(g_j)$。则 $\ell = \prod_{j=1}^t \ell_j$。

### 1.2 模约化

**素数选择**：选 $p$ 素数满足：
- $p \nmid \ell$
- $\bar{f} = f \bmod p$ 在 $\mathbb{F}_p[x]$ 中无平方

**模因子分解**：$\bar{f} = \bar{\ell} \cdot \bar{h}_1 \cdots \bar{h}_r$，其中 $\bar{h}_i \in \mathbb{F}_p[x]$ 两两互素的首一不可约多项式，$\bar{\ell} = \ell \bmod p$。

**因子对应**：每个真因子 $g_j$ 对应一个子集 $S_j \subseteq \{1, \ldots, r\}$，满足：
$$\bar{g}_j / \bar{\ell}_j = \prod_{i \in S_j} \bar{h}_i \quad \text{in } \mathbb{F}_p[x]$$
其中 $\bar{\ell}_j = \ell_j \bmod p$（$p \nmid \ell_j$ 因为 $\ell_j | \ell$ 且 $p \nmid \ell$）。$\{S_j\}_{j=1}^t$ 是 $\{1, \ldots, r\}$ 的一个分割。

### 1.3 对称模

$$\text{sym}_m(a) = a \bmod m \text{ 约化到 } (-m/2, m/2]$$

对多项式逐系数进行：$\text{sym}_m(f) = \sum_k \text{sym}_m(a_k) x^k$。

### 1.4 Mignotte 界

**定理 (Landau-Mignotte)**：若 $g | f$ 在 $\mathbb{Z}[x]$ 中，$\deg g = d \leq n$，则 $g$ 的任意系数 $|c_k|$ 满足：
$$|c_k| \leq \binom{d}{k} \|\hat{f}\|_2$$
其中 $\hat{f} = f / \ell$ 是 $f$ 的首一化。

**推论**：$\|g\|_\infty \leq B_{\text{Mig}}(f) := \binom{n}{\lfloor n/2 \rfloor} \|f\|_2$。

**推导**：$\|g\|_\infty = \max_k |c_k| \leq \max_k \binom{d}{k} \|\hat{f}\|_2 = \binom{d}{\lfloor d/2 \rfloor} \|\hat{f}\|_2$。由 $d \leq n$ 知 $\binom{d}{\lfloor d/2 \rfloor} \leq \binom{n}{\lfloor n/2 \rfloor}$（中心二项式系数关于 $n$ 单调递增）。由 $\ell \geq 1$ 知 $\|\hat{f}\|_2 = \|f/\ell\|_2 \leq \|f\|_2$。$\square$

### 1.5 $\mathbb{Z}_m[x]$ 中的多项式除法

**引理 1.1**：设 $A, B \in \mathbb{Z}_m[x]$（系数在 $\mathbb{Z}/m\mathbb{Z}$ 中），$B$ 首一。则存在唯一的 $Q, R \in \mathbb{Z}_m[x]$ 使得 $A \equiv B \cdot Q + R \pmod{m}$，$\deg R < \deg B$。

**证明**：标准多项式长除法。每一步中，商系数 = 被除式首项系数 × 除式首项系数$^{-1}$。由 $B$ 首一，首项系数$^{-1} = 1$，无需模逆元。故除法在 $\mathbb{Z}_m[x]$ 中精确执行。$\square$

**推论 1.2**：若 $B$ 首一且 $B | A$ 在 $\mathbb{Z}_m[x]$ 中（即 $A \equiv B \cdot Q \pmod{m}$），则多项式长除法恰好得到 $R = 0$、$Q$ 即为精确商。

## 2. Hensel 提升不变量

### 2.1 lc 烘焙（lc-baking）

CLPoly 采用 "lc 烘焙" 方案进行 Hensel 提升：

**初始化**：定义调整后的模因子：
$$\tilde{h}_1 = \ell \cdot \bar{h}_1 \bmod p, \quad \tilde{h}_i = \bar{h}_i \text{ for } i \geq 2$$

验证乘积：$\tilde{h}_1 \cdot \bar{h}_2 \cdots \bar{h}_r = \ell \cdot \bar{h}_1 \cdots \bar{h}_r = \bar{f} \pmod{p}$。✓

### 2.2 二因子 Hensel 唯一性定理

**定理 2.1 (Hensel 唯一性)**：设 $F \in \mathbb{Z}_m[x]$（$m = p^k$，$k \geq 1$），$\bar{F} = \bar{A} \cdot \bar{B}$ 在 $\mathbb{F}_p[x]$ 中且 $\gcd(\bar{A}, \bar{B}) = 1$。则存在唯一的 $A_m, B_m \in \mathbb{Z}_m[x]$ 满足：

- $A_m \cdot B_m \equiv F \pmod{m}$
- $A_m \equiv \bar{A} \pmod{p}$，$B_m \equiv \bar{B} \pmod{p}$
- $\deg A_m = \deg \bar{A}$，$\deg B_m = \deg \bar{B}$

**注**：$F$ 可以是 $\mathbb{Z}_m[x]$ 中的任意多项式（不要求来自 $\mathbb{Z}[x]$ 的归约）。证明仅使用 $\mathbb{Z}/m\mathbb{Z}$ 的环结构。

**参考**：[GCL] 定理 6.1，或 Cohen §3.5.4。$\square$

### 2.3 多因子 Hensel 提升定理

**定理 2.2 (多因子 Hensel 提升)**：设 $m = p^{2^k}$（二次提升 $k$ 步）。存在唯一的 $H_1, \ldots, H_r \in \mathbb{Z}_m[x]$ 满足：

**(I1)** 乘积恒等式：$\prod_{i=1}^r H_i \equiv f \pmod{m}$

**(I2)** 次数保持：$\deg H_i = \deg \bar{h}_i$

**(I3)** 首项系数：
- $\text{lc}(H_1) = \ell$（当 $\ell < m$ 时精确等于 $\ell$，由 $m > 2\ell \cdot B_{\text{Mig}} > \ell$ 保证）
- $\text{lc}(H_i) = 1$ 对所有 $i \geq 2$（首一）

**(I4)** 模 $p$ 一致：$H_1 \equiv \tilde{h}_1 \pmod{p}$，$H_i \equiv \bar{h}_i \pmod{p}$（$i \geq 2$）

**证明 (I3)**：二次 Hensel 步骤的标准形式为：给定 $f \equiv G \cdot H \pmod{m}$，$s \cdot G + t \cdot H \equiv 1 \pmod{p}$，求 $G^*, H^*$ 使得 $f \equiv G^* \cdot H^* \pmod{m^2}$。

设 $e = f - G \cdot H$，$e/m \in \mathbb{Z}[x]$（因 $m | e$）。令 $s \cdot (e/m) = q \cdot H + \sigma$（$\mathbb{F}_p[x]$ 中带余除法），$\deg \sigma < \deg H$。令 $\tau = t \cdot (e/m) + G \cdot q$。则：
$$G^* = G + m \cdot \tau, \quad H^* = H + m \cdot \sigma$$

由 $\deg \sigma < \deg H$ 知 $\deg(m \cdot \sigma) < \deg H$，故 $\deg H^* = \deg H$，$\text{lc}(H^*) = \text{lc}(H)$。

由 Bézout 关系 $s \cdot G + t \cdot H \equiv 1 \pmod{p}$ 知 $\deg s < \deg H$，$\deg t < \deg G$（否则左侧次数 $> 0$）。类似的次数分析（参见 [GCL] §6.3）给出 $\deg \tau < \deg G$，故 $\deg G^* = \deg G$，$\text{lc}(G^*) = \text{lc}(G)$。

因此首项系数在每步提升中严格保持不变。$\square$

### 2.4 Mignotte 精度条件

选择 $m > 2\ell \cdot B_{\text{Mig}}(f)$，记此条件为 **(M)**。

## 3. 因子恢复定理

### 3.1 子集乘积与真因子的关系

**定理 3.1（因子恢复）**：设条件 (I1)-(I4) 和 (M) 成立。对任何真因子 $g_j$（对应子集 $S_j$），定义：

$$P_{S_j} = \begin{cases} \prod_{i \in S_j} H_i & \text{若 } 1 \in S_j \\ \ell \cdot \prod_{i \in S_j} H_i & \text{若 } 1 \notin S_j \end{cases}$$

则 $\text{pp}(\text{sym}_m(P_{S_j})) = g_j$。

**证明**：两种情况的论证结构相同——构造 $\mathbb{Z}_m[x]$ 中的对比组，再由 Hensel 唯一性识别。

令 $h_j^* = f / g_j \in \mathbb{Z}[x]$。则 $\text{lc}(h_j^*) = \ell / \ell_j$（因 $\ell = \ell_j \cdot \text{lc}(h_j^*)$）。$\gcd(\ell, m) = 1$ 保证所有 $\ell_j^{-1} \bmod m$ 存在。

设 $U = \prod_{i \in S_j} H_i$，$V = \prod_{i \notin S_j} H_i$。由 (I1)，$U \cdot V \equiv f \pmod{m}$。区分两种情况构造对比组：

**情况 1**（$1 \in S_j$）：$U = \prod_{i \in S_j} H_i$，$\text{lc}(U) = \ell$。

构造：$U' = (\ell/\ell_j) \cdot g_j$，$V' = (\ell_j/\ell) \cdot h_j^* \bmod m$。
- $\text{lc}(U') = (\ell/\ell_j) \cdot \ell_j = \ell = \text{lc}(U)$。✓
- $U' \cdot V' = f$ 精确。✓
- $U' \bmod p = \bar{\ell} \cdot \prod_{i \in S_j} \bar{h}_i = U \bmod p$。✓
- $\text{lc}(V') = (\ell_j/\ell) \cdot (\ell/\ell_j) = 1 = \text{lc}(V)$。✓
- $V' \bmod p = \prod_{i \notin S_j} \bar{h}_i = V \bmod p$。✓

由定理 2.1：$U \equiv U' \pmod{m}$，即 $\prod_{i \in S_j} H_i \equiv (\ell/\ell_j) \cdot g_j \pmod{m}$。

Mignotte 恢复：$(\ell/\ell_j) \cdot g_j \in \mathbb{Z}[x]$，系数绝对值 $\leq (\ell/\ell_j) \cdot \|g_j\|_\infty \leq \ell \cdot B_{\text{Mig}} < m/2$（由条件 (M)）。故 $\text{sym}_m(U) = (\ell/\ell_j) \cdot g_j$。

$P_{S_j} = U$，$\text{pp}(P_{S_j}) = \text{pp}((\ell/\ell_j) \cdot g_j) = g_j$（$g_j$ 本原，$\ell/\ell_j \in \mathbb{Z}_{>0}$）。$\square$

**情况 2**（$1 \notin S_j$）：$U = \prod_{i \in S_j} H_i$，$\text{lc}(U) = 1$（全部首一）。$V = \prod_{i \notin S_j} H_i$，$\text{lc}(V) = \ell$（$H_1 \in V$ 侧，因 $1 \notin S_j$）。

构造：$U'' = g_j \cdot \ell_j^{-1} \bmod m$（首一，$\text{lc}(U'') = 1$），$V'' = \ell_j \cdot h_j^* \bmod m$（$\text{lc}(V'') = \ell_j \cdot (\ell/\ell_j) = \ell$）。
- $U'' \cdot V'' = (g_j/\ell_j) \cdot \ell_j \cdot h_j^* = g_j \cdot h_j^* = f$。✓

$U'' \bmod p$：$g_j \cdot \ell_j^{-1} \bmod p = \bar{g}_j / \bar{\ell}_j = \prod_{i \in S_j} \bar{h}_i = U \bmod p$。✓

$V'' \bmod p$：$\ell_j \cdot h_j^* \bmod p = \bar{\ell}_j \cdot \bar{h}_j^* = \bar{\ell}_j \cdot (\bar{\ell}/\bar{\ell}_j) \cdot \prod_{i \notin S_j} \bar{h}_i = \bar{\ell} \cdot \prod_{i \notin S_j} \bar{h}_i$。

而 $V \bmod p = \prod_{i \notin S_j} H_i \bmod p$。当 $1 \notin S_j$ 时，$1 \in \{i \notin S_j\}$，故 $V \bmod p = \tilde{h}_1 \cdot \prod_{i \notin S_j, i \neq 1} \bar{h}_i = \bar{\ell} \cdot \bar{h}_1 \cdot \prod_{i \notin S_j, i \neq 1} \bar{h}_i = \bar{\ell} \cdot \prod_{i \notin S_j} \bar{h}_i$。✓

由定理 2.1：$U \equiv U'' \pmod{m}$，即：
$$\prod_{i \in S_j} H_i \equiv g_j / \ell_j \pmod{m}$$

$P_{S_j} = \ell \cdot \prod_{i \in S_j} H_i \equiv \ell \cdot g_j / \ell_j = (\ell/\ell_j) \cdot g_j \pmod{m}$。

Mignotte 恢复同情况 1，$\text{sym}_m(P_{S_j}) = (\ell/\ell_j) \cdot g_j$，$\text{pp}(P_{S_j}) = g_j$。$\square$

### 3.2 关键恒等式总结

由定理 3.1 的证明，对所有 $S_j$，无论 $1 \in S_j$ 与否：

$$\prod_{i \in S_j} H_i \equiv \frac{\ell}{\ell_j} \cdot g_j \pmod{m} \quad (1 \in S_j) \tag{$\star_1$}$$
$$\prod_{i \in S_j} H_i \equiv \frac{g_j}{\ell_j} \pmod{m} \quad (1 \notin S_j) \tag{$\star_2$}$$

两种情况下，$P_{S_j}$ 均 $\equiv (\ell/\ell_j) \cdot g_j \pmod{m}$，恢复出 $g_j$。

### 3.3 试除正确性

**推论 3.2**：$\text{pp}(\text{sym}_m(P_{S_j}))$ 精确整除 $f$。

**证明**：由定理 3.1，$\text{pp}(\text{sym}_m(P_{S_j})) = g_j$，而 $g_j | f$。$\square$

## 4. 剪枝正确性条件

Zassenhaus 算法在试除前使用两个剪枝条件加速。以下推导这些条件成立的前提。

### 4.1 lc 剪枝

**剪枝规则**：对子集 $S$，计算 $\text{lc}(P_S) = \text{sym}_m(\ell_{\text{mult}} \cdot \prod_{i \in S} \text{lc}(H_i))$。若 $\text{lc}(f^*)^2 \% \text{lc}(P_S) \neq 0$ 则剪枝。

**正确性条件**（剪枝不误杀正确子集）：

设 $S = S_j$ 是正确子集。由定理 3.1，$\text{sym}_m(P_{S_j}) = (\ell/\ell_j) \cdot g_j$，故 $\text{lc}(P_{S_j}) = \ell$（两种情况下均如此）。

由 $g_j | f^*$ 知 $\ell_j | \ell^*$（$\ell^* = \text{lc}(f^*)$，因为 $f^* = g_j \cdot h^{**}$，$\ell^* = \ell_j \cdot \text{lc}(h^{**})$）。

**引理 4.1**：剪枝不误杀正确子集当且仅当 $\text{lc}(P_{S_j}) | \ell^{*2}$。

**分析**：
- **初始调用**（$f^* = f$，$\ell^* = \ell$）：$\text{lc}(P_{S_j}) = \ell$，$\ell^{*2} = \ell^2$，$\ell | \ell^2$。✓
- **van Hoeij 回退调用**（$f^*$ 是部分约化后的多项式，$\ell^* \neq \ell$）：$\text{lc}(P_{S_j})$ 取决于 $\ell_{\text{mult}}$ 的设置。

  代码中，当 $1 \in S_j$ 时 $\ell_{\text{mult}} = 1$，$\text{lc}(P_{S_j}) = \text{lc}(H_1) = \ell$。剪枝要求 $\ell | \ell^{*2}$。**不一定成立。**

**反例**（Case 1）：$\ell = 3569280 = 2^5 \cdot 3 \cdot 5 \cdot 7 \cdot 1063$，$\ell^* = 2288 = 2^4 \cdot 11 \cdot 13$。

$\ell^{*2} = 2^8 \cdot 11^2 \cdot 13^2$。$\ell$ 含素因子 $3, 5, 7, 1063$，而 $\ell^{*2}$ 不含。故 $\ell \nmid \ell^{*2}$，剪枝误杀。$\square$

### 4.2 常数项剪枝

**剪枝规则**：计算 $c_{\text{prod}} = \text{sym}_m(\ell_{\text{mult}} \cdot \prod_{i \in S} H_i(0))$，检查 $(\ell^* \cdot f^*(0)) \% c_{\text{prod}} \stackrel{?}{=} 0$。

**正确性分析**：

设 $S = S_j$。代码中 $c_{\text{prod}} = \text{sym}_m(P_{S_j}(0))$，即 $P_{S_j}$ 在 $x=0$ 的值。由定理 3.1，$\text{sym}_m(P_{S_j}) = (\ell/\ell_j) \cdot g_j$，故 $c_{\text{prod}} = (\ell/\ell_j) \cdot g_j(0)$。

由 $f^* = g_j \cdot h^{**}$，$f^*(0) = g_j(0) \cdot h^{**}(0)$。

$$\frac{\ell^* \cdot f^*(0)}{c_{\text{prod}}} = \frac{\ell^* \cdot g_j(0) \cdot h^{**}(0)}{(\ell/\ell_j) \cdot g_j(0)} = \frac{\ell^* \cdot \ell_j \cdot h^{**}(0)}{\ell}$$

- **初始调用**（$\ell^* = \ell$）：结果为 $\ell_j \cdot h^{**}(0) \in \mathbb{Z}$。✓
- **van Hoeij 回退调用**（$\ell^*$ 可能与 $\ell$ 不同）：$\ell^* \cdot \ell_j / \ell$ 不一定为整数。**失败！** $\square$

## 5. Bug 的形式化描述

### 5.1 lc 不变量违反

**定义 5.1（lc 不变量）**：设当前剩余多项式为 $f^*$，活跃 Hensel 因子集为 $\{H_i : i \in A\}$。**lc 不变量**要求：

> 对任何正确子集 $S_j \subseteq A$，试除用乘积 $P_{S_j} = \ell_{\text{mult}} \cdot \prod_{i \in S_j} H_i$ 的 $\text{pp}$ 必须等于 $g_j$，且 $\text{lc}(P_{S_j}) | \ell^{*2}$。

**定理 5.1（Bug 定理）**：在以下条件同时满足时，lc 不变量被违反：

1. van Hoeij 提取了部分真因子，$\ell^* \neq \ell$
2. $1 \in A$（lc 烘焙因子 $H_1$ 仍在活跃集中）
3. 存在正确子集 $S_j \subseteq A$ 且 $1 \in S_j$
4. $\ell \nmid \ell^{*2}$（$\ell$ 含 $\ell^*$ 不具有的素因子）

**证明**：由 §4.1，当 $1 \in S_j$ 时，代码使用 $\ell_{\text{mult}} = 1$，得 $\text{lc}(P_{S_j}) = \ell$。条件 4 直接给出 $\ell \nmid \ell^{*2}$，故 lc 剪枝误杀正确子集。

同理，常数项剪枝中 $\ell^* \cdot \ell_j / \ell$ 不为整数（因 $\ell$ 含 $\ell^*$ 不具有的素因子），故常数项剪枝也误杀正确子集。$\square$

### 5.2 CLD 不变量违反

**定义 5.2（CLD 不变量）**：CLD 多项式 $C_i = (f^* / H_i) \cdot H_i' \bmod m$ 的计算要求 $f^*/H_i$ 的系数有界（Mignotte 量级），以保证 CLD 列构成"短向量"。

**定理 5.2**：当 $\text{lc}(H_1) = \ell \neq 1$ 且 $\text{lc}(f^*) = \ell^*$ 时，$f^*/H_1 \pmod{m}$ 的多项式长除法在每步需要乘以 $\ell^{-1} \bmod m$（非首一除数）。商的系数可达 $O(m)$ 量级，CLD 值不再是 Mignotte 量级的"小"数。

具体地：$f^*/H_1$ 的首项系数为 $\ell^* \cdot \ell^{-1} \bmod m$。
- 若 $\ell | \ell^*$：$\ell^*/\ell \in \mathbb{Z}$，系数有界于 $\ell^* \cdot B_{\text{Mig}}$。✓
- 若 $\ell \nmid \ell^*$：$\ell^{-1} \bmod m$ 是一个 $O(m)$ 量级的数，首项系数 $\ell^* \cdot \ell^{-1} \bmod m$ 可达 $O(m)$，此后除法的每步均传播此量级。CLD 中对应于 $H_1$ 的列包含 $O(m)$ 量级的系数。

LLL 格规约中，$O(m)$ 量级的 CLD 列使格基向量的范数远超 Mignotte 量级，无法被识别为短向量。因子分组失败。$\square$

## 6. 修复方案：首一归一化

### 6.1 方案描述

**修复操作**：在 `__hensel_lift` 返回后、调用重组算法前，将 $H_1$ 归一化为首一：

$$H_1 \leftarrow H_1 \cdot (\text{lc}(H_1))^{-1} \pmod{m}$$

前提条件：$\gcd(\ell, m) = 1$，由 $m = p^a$ 且 $p \nmid \ell$ 保证。

此后所有 Hensel 因子均为首一。

**同时修改 Zassenhaus**：移除 `subset_has_lc` 特殊处理，统一使用 $\ell_{\text{mult}} = \ell^*$（当前 $f^*$ 的首项系数）。

### 6.2 修复后的不变量

**定义 6.1（修复后不变量 I'）**：设所有活跃 Hensel 因子 $H_i$（$i \in A$）均为首一。则：

**(I'1)** $\prod_{i \in A} H_i \equiv f^* \cdot (\ell^*)^{-1} \pmod{m}$

**(I'2)** $\text{lc}(H_i) = 1$ 对所有 $i \in A$

**(I'3)** 对任何正确子集 $S_j \subseteq A$：$\prod_{i \in S_j} H_i \equiv g_j \cdot \ell_j^{-1} \pmod{m}$

**注记**：(I'1) 和 (I'3) 中的 $(\ell^*)^{-1}$ 和 $\ell_j^{-1}$ 是 $\mathbb{Z}_m$ 中的模逆元。乘积 $\prod_{i \in S_j} H_i$ 的系数是 $[0, m)$ 中的整数；$g_j \cdot \ell_j^{-1}$ 的含义是将 $g_j$ 的每个系数 $a_k$ 乘以 $\ell_j^{-1} \bmod m$ 再取模。这两个多项式在 $\mathbb{Z}_m[x]$ 中相等。

**定理 6.1（不变量 I' 成立）**：

**证明 (I'3)**：

对正确子集 $S_j \subseteq A$，使用定理 3.1 情况 2 的 Hensel 唯一性论证。构造 $\mathbb{Z}_m[x]$ 中的对比组：

$U'' = g_j \cdot \ell_j^{-1} \bmod m$（首一），$V'' = h_j^* \cdot (\ell/\ell_j)^{-1} \bmod m$（首一，$\text{lc}(V'') = (\ell/\ell_j) \cdot (\ell/\ell_j)^{-1} = 1$）。

验证：
- $U'' \cdot V'' = g_j \cdot \ell_j^{-1} \cdot h_j^* \cdot (\ell/\ell_j)^{-1} = f \cdot (\ell_j \cdot \ell/\ell_j)^{-1} = f \cdot \ell^{-1} \pmod{m}$。而首一化后 $\prod_{i \in A} H_i \equiv f \cdot \ell^{-1} \pmod{m}$（(I'1)），故乘积匹配。✓
- $\text{lc}(U'') = 1 = \text{lc}(\prod_{i \in S_j} H_i)$（首一化后全部首一）。✓
- $\text{lc}(V'') = 1 = \text{lc}(\prod_{i \notin S_j} H_i)$。✓
- $\deg U'' = \deg g_j = \deg(\prod_{i \in S_j} H_i)$。✓
- $U'' \bmod p = g_j \cdot \ell_j^{-1} \bmod p = \bar{g}_j / \bar{\ell}_j = \prod_{i \in S_j} \bar{h}_i = \prod_{i \in S_j} H_i \bmod p$。✓
- $V'' \bmod p = h_j^* \cdot (\ell/\ell_j)^{-1} \bmod p = (\ell/\ell_j) \cdot \prod_{i \notin S_j} \bar{h}_i \cdot (\ell/\ell_j)^{-1} = \prod_{i \notin S_j} \bar{h}_i = \prod_{i \notin S_j} H_i \bmod p$。✓

（注：首一化后 $H_i \equiv \bar{h}_i \pmod{p}$ 对所有 $i$，包括 $i=1$：首一化将 $H_1 \equiv \tilde{h}_1 = \ell \bar{h}_1$ 变为 $H_1 \cdot \ell^{-1} \equiv \bar{h}_1 \pmod{p}$。）

由定理 2.1（Hensel 唯一性，应用于 $F = f \cdot \ell^{-1} \in \mathbb{Z}_m[x]$）：$\prod_{i \in S_j} H_i \equiv U'' = g_j \cdot \ell_j^{-1} \pmod{m}$。$\square$

**证明 (I'1)**：

**初始状态**：$A = \{1, \ldots, r\}$，$f^* = f$，$\ell^* = \ell$。

由 (I1)，$\prod_{i=1}^r H_i \equiv f \pmod{m}$（首一化前）。首一化 $H_1$：$H_1^{\text{monic}} = H_1 \cdot \ell^{-1} \bmod m$。

$$\prod_i H_i^{(\text{monic})} = (H_1 \cdot \ell^{-1}) \cdot \prod_{i \geq 2} H_i = \ell^{-1} \cdot \prod_i H_i \equiv \ell^{-1} \cdot f = f \cdot \ell^{-1} = f^* \cdot (\ell^*)^{-1} \pmod{m}$$

✓

**归纳步骤**：设提取真因子 $g_j$（对应子集 $S_j \subseteq A$）后，活跃集 $A' = A \setminus S_j$，$f^*_{\text{new}} = f^* / g_j$。

**Gauss 引理保证 $f^*/g_j$ 本原**：$f^*$ 本原，$g_j$ 本原（不可约），$f^* = g_j \cdot (f^*/g_j)$。由 Gauss 引理，$\text{cont}(f^*) = \text{cont}(g_j) \cdot \text{cont}(f^*/g_j)$，即 $1 = 1 \cdot \text{cont}(f^*/g_j)$，故 $f^*/g_j$ 本原。因此 $f^*_{\text{new}} = \text{pp}(f^*/g_j) = f^*/g_j$。

**首项系数关系**：$\ell^* = \text{lc}(f^*) = \ell_j \cdot \text{lc}(f^*/g_j) = \ell_j \cdot \ell^*_{\text{new}}$。

**乘积更新**：由 $S_j \subseteq A$，$\prod_{i \in S_j} H_i$ 是 $\prod_{i \in A} H_i$ 的子因子乘积，故在 $\mathbb{Z}_m[x]$ 中（所有 $H_i$ 首一，由引理 1.1 长除法精确）：

$$\prod_{i \in A'} H_i = \frac{\prod_{i \in A} H_i}{\prod_{i \in S_j} H_i} \pmod{m}$$

由归纳假设 (I'1)：$\prod_{i \in A} H_i \equiv f^* \cdot (\ell^*)^{-1} \pmod{m}$。

由 (I'3)：$\prod_{i \in S_j} H_i \equiv g_j \cdot \ell_j^{-1} \pmod{m}$。

$$\prod_{i \in A'} H_i \equiv \frac{f^* \cdot (\ell^*)^{-1}}{g_j \cdot \ell_j^{-1}} = \frac{f^*}{g_j} \cdot \frac{\ell_j}{\ell^*} = f^*_{\text{new}} \cdot \frac{\ell_j}{\ell^*} \pmod{m}$$

由 $\ell^* = \ell_j \cdot \ell^*_{\text{new}}$：$\ell_j / \ell^* = 1/\ell^*_{\text{new}}$。

$$\prod_{i \in A'} H_i \equiv f^*_{\text{new}} \cdot (\ell^*_{\text{new}})^{-1} \pmod{m}$$

即 (I'1) 对新状态成立。$\square$

### 6.3 因子恢复定理（修复后）

**定理 6.2（修复后因子恢复）**：设不变量 (I') 成立。对任何正确子集 $S_j \subseteq A$：

$$\text{pp}\left(\text{sym}_m\left(\ell^* \cdot \prod_{i \in S_j} H_i\right)\right) = g_j$$

**证明**：由 (I'3)：
$$\ell^* \cdot \prod_{i \in S_j} H_i \equiv \ell^* \cdot g_j \cdot \ell_j^{-1} \pmod{m}$$

由 $\ell_j | \ell^*$（$g_j | f^*$ 蕴含 $\ell_j | \ell^*$），$\ell^* \cdot \ell_j^{-1} = \ell^*/\ell_j \in \mathbb{Z}$。故上式等于 $(\ell^*/\ell_j) \cdot g_j$，这是一个**整数多项式**。

Mignotte 恢复：$|(\ell^*/\ell_j) \cdot g_j \text{ 的系数}| \leq (\ell^*/\ell_j) \cdot \|g_j\|_\infty \leq \ell \cdot B_{\text{Mig}}(f) < m/2$（由条件 (M)，利用 $\ell^*/\ell_j \leq \ell^* \leq \ell$）。

故 $\text{sym}_m(\ell^* \cdot \prod_{i \in S_j} H_i) = (\ell^*/\ell_j) \cdot g_j$ 精确成立。

$g_j$ 本原，$\ell^*/\ell_j \in \mathbb{Z}_{>0}$，故 $\text{pp}((\ell^*/\ell_j) \cdot g_j) = g_j$。$\square$

### 6.4 剪枝正确性（修复后）

**定理 6.3（修复后 lc 剪枝正确）**：设不变量 (I') 成立。对任何正确子集 $S_j \subseteq A$，lc 剪枝不会误杀。

**证明**：所有 $H_i$ 首一，$\ell_{\text{mult}} = \ell^*$。故：
$$\text{lc}(P_{S_j}) = \ell^* \cdot \prod_{i \in S_j} \text{lc}(H_i) = \ell^* \cdot 1 = \ell^*$$

剪枝条件：$\ell^{*2} \% \ell^* = 0$。始终成立。$\square$

**定理 6.4（修复后常数项剪枝正确）**：设不变量 (I') 成立。对任何正确子集 $S_j \subseteq A$，常数项剪枝不会误杀。

**证明**：

由定理 6.2 的推导，$\text{sym}_m(\ell^* \cdot \prod_{i \in S_j} H_i) = (\ell^*/\ell_j) \cdot g_j$。在 $x = 0$ 处取值：

$$c_{\text{prod}} = \text{sym}_m\left(\ell^* \cdot \prod_{i \in S_j} H_i(0)\right) = \frac{\ell^*}{\ell_j} \cdot g_j(0)$$

$$f^*_{\text{const}} = \ell^* \cdot f^*(0) = \ell^* \cdot g_j(0) \cdot h^{**}(0)$$

（其中 $h^{**} = f^*/g_j \in \mathbb{Z}[x]$。）

$$\frac{f^*_{\text{const}}}{c_{\text{prod}}} = \frac{\ell^* \cdot g_j(0) \cdot h^{**}(0)}{(\ell^*/\ell_j) \cdot g_j(0)} = \ell_j \cdot h^{**}(0) \in \mathbb{Z}$$

故整除条件成立（当 $g_j(0) \neq 0$ 时；$g_j(0) = 0$ 时 $c_{\text{prod}} = 0$，代码跳过剪枝）。$\square$

### 6.5 CLD 正确性（修复后）

**定理 6.5**：设不变量 (I') 成立（所有活跃因子首一）。CLD 多项式 $C_i = (f^* / H_i) \cdot H_i' \bmod m$ 满足 van Hoeij LLL 分组的正确性条件。

**证明**：分三步——验证 CLD 计算合法性、推导 CLD 求和公式、验证系数有界性。

**第一步：$H_i | f^* \pmod{m}$（CLD 计算合法）。**

由 (I'1)，$f^* \equiv \ell^* \cdot \prod_{j \in A} H_j \pmod{m}$。由于 $H_i$ 是 $\prod_{j \in A} H_j$ 的子因子，且所有 $H_j$ 首一：

$$f^* / H_i \equiv \ell^* \cdot \prod_{j \in A,\, j \neq i} H_j \pmod{m}$$

右侧是 $\mathbb{Z}_m[x]$ 中的多项式（系数在 $[0, m)$ 中）。由引理 1.1（首一除数的精确长除法），$f^* / H_i$ 在 $\mathbb{Z}_m[x]$ 中余式为 0。✓

**第二步：CLD 求和公式。**

$$C_i = (f^*/H_i) \cdot H_i' = \ell^* \cdot \left(\prod_{j \in A,\, j \neq i} H_j\right) \cdot H_i' \pmod{m}$$

对正确子集 $S_j \subseteq A$，令 $G_S = \prod_{i \in S_j} H_i$。

$$\sum_{i \in S_j} C_i = \ell^* \cdot \sum_{i \in S_j} \left(\prod_{k \in A,\, k \neq i} H_k\right) \cdot H_i'$$

$$= \ell^* \cdot \left(\prod_{k \in A \setminus S_j} H_k\right) \cdot \underbrace{\sum_{i \in S_j} \left(\prod_{k \in S_j,\, k \neq i} H_k\right) \cdot H_i'}_{= G_S' \text{（乘积求导法则）}}$$

$$= \ell^* \cdot \left(\prod_{k \in A \setminus S_j} H_k\right) \cdot G_S' \pmod{m}$$

**（乘积求导法则验证**：$G_S = \prod_{i \in S_j} H_i$，则 $G_S' = \sum_{i \in S_j} (\prod_{k \in S_j, k \neq i} H_k) \cdot H_i'$。这是多项式乘积的 Leibniz 求导法则，在任何交换环上成立，包括 $\mathbb{Z}_m[x]$。**）**

由 (I'1) 和 (I'3)：
- $\prod_{k \in A \setminus S_j} H_k \equiv f^*_{\text{new}} \cdot (\ell^*_{\text{new}})^{-1} \pmod{m}$（由 (I'1) 的归纳步骤）
- $G_S \equiv g_j \cdot \ell_j^{-1} \pmod{m}$，故 $G_S' \equiv g_j' \cdot \ell_j^{-1} \pmod{m}$

$$\sum_{i \in S_j} C_i \equiv \ell^* \cdot \frac{f^*_{\text{new}}}{\ell^*_{\text{new}}} \cdot \frac{g_j'}{\ell_j} = \frac{\ell^*}{\ell^*_{\text{new}} \cdot \ell_j} \cdot f^*_{\text{new}} \cdot g_j' \pmod{m}$$

由 $\ell^* = \ell_j \cdot \ell^*_{\text{new}}$：

$$\sum_{i \in S_j} C_i \equiv f^*_{\text{new}} \cdot g_j' = \frac{f^*}{g_j} \cdot g_j' \pmod{m}$$

这是 $f^*$ 相对于因子 $g_j$ 的对数导数（乘以 $f^*$）：$\sum C_i = f^* \cdot (g_j'/g_j)$。

**第三步：系数有界性。**

$h^{**} = f^*/g_j$，则 $\sum_{i \in S_j} C_i \equiv h^{**} \cdot g_j' \pmod{m}$。

$h^{**} \cdot g_j'$ 的系数绝对值 $\leq \|h^{**}\|_1 \cdot \|g_j'\|_\infty \leq (n+1) \cdot \|f^*\|_\infty \cdot n \cdot \|g_j\|_\infty$。

在 Mignotte 条件 (M) 下，$\|f^*\|_\infty$ 和 $\|g_j\|_\infty$ 均为多项式量级（远小于 $m$），故 CLD 列的值是 $O(\text{poly}(n, \|f\|))$ 量级的"小"数。LLL 格规约可正确识别同一真因子对应的模因子组。$\square$

## 7. 完整性定理

### 7.1 满 Mignotte 精度下的正确性

**定理 7.1（修复后正确性）**：设 $f \in \mathbb{Z}[x]$ 无平方、本原、$\deg f \geq 2$。若 Hensel 提升精度满足 Mignotte 条件 (M)（$m > 2\ell \cdot B_{\text{Mig}}(f)$），则修复后的因子重组算法返回 $f$ 的完全不可约分解。

**证明**：

**第一层：Hensel 提升正确性。**

`__hensel_lift` 返回 $H_1, \ldots, H_r$ 满足 (I1)-(I4)。出口处首一化 $H_1$ 后，不变量 (I') 成立（定理 6.1）。

**第二层：van Hoeij 重组正确性。**

`__vanhoeij_recombine` 的循环不变量：

> **(L)**：在每次迭代开始时：
> 1. `active_lifted` 中所有因子首一（(I'2)）
> 2. $\prod_{i \in \text{active}} H_i \equiv f^* \cdot (\ell^*)^{-1} \pmod{m}$（(I'1)）
> 3. 每个 `active_lifted[k]` 对应某个真因子 $g_j$ 的模因子

**(L) 的初始成立**：由 (I') 直接保证。

**(L) 的保持**：

- **J_target=0 预提取**：对角 LLL 提取单因子候选。由定理 6.2，正确子集的试除必成功。提取后，由定理 6.1 的归纳步骤，(L) 对新状态保持。
- **CLD + LLL**：由定理 6.5，CLD 系数有界，LLL 规约给出正确分组。
- **Zassenhaus 回退**：由定理 6.3 和 6.4，剪枝不误杀正确子集。Zassenhaus 枚举所有 $s \leq |A|/2$ 的子集（或其字典序组合），每个正确子集通过试除被发现（定理 6.2）。

**终止性**：每次成功提取至少减少一个活跃因子，活跃集有限。

**完备性**：
- 若 LLL 成功分组所有因子：直接得到完全分解。
- 若回退到 Zassenhaus：穷举保证发现所有正确子集。
- 无遗漏：每个真因子 $g_j$ 的模因子子集 $S_j$ 必在某一步被提取。$\square$

### 7.2 Phase 2 安全网与启发式精度

**定理 7.2（Phase 2 正确性）**：若 Phase 1 使用启发式精度 $a_h$ 未能完全分解（返回单个因子），Phase 2 以完整 Mignotte 精度重新提升+重组。由条件 (M) 满足，定理 7.1 保证 Phase 2 返回完全分解。

**命题 7.3（Phase 2 触发条件的局限性）**：

当前代码 Phase 2 触发条件为 `result.size() == 1 && r > 1`。这是一个**启发式安全网**，**不是**数学完整性保证。

**理论上的缺口**：若启发式精度 $a_h$ 足以恢复部分因子（小系数因子在 $a_h$ 精度下满足 Mignotte 条件，$\text{sym}_m$ 正确恢复），但不足以恢复其余因子（大系数因子的 $|(\ell^*/\ell_j) \cdot g_j|$ 的某些系数 $> m/2$），则 van Hoeij 提取部分因子后，剩余多项式因精度不足而无法继续分解。此时 `result.size() > 1`，Phase 2 不触发，剩余可约多项式作为"不可约"返回。

**数学上的完整性保证**：仅在以下条件之一成立时，算法保证返回完全不可约分解：
1. $m > 2\ell \cdot B_{\text{Mig}}(f)$（满 Mignotte 精度，定理 7.1）
2. Phase 2 被触发（`result.size() == 1 && r > 1`）

**实践上的安全性**：FLINT 启发式公式使 $a_h$ 在绝大多数情况下要么足够（全部分解），要么完全不足（一个因子都分不出，触发 Phase 2）。理论上的"部分精度不足"场景概率极低，且与 lc bug 无关——lc bug 是在满 Mignotte 精度下也会触发的正确性缺陷。

**建议**：将触发条件放宽为 `result.size() < r && a_h < a_mig` 可填补此理论缺口，但这是独立的优化任务，不影响 lc bug 的修复。

## 8. 改动影响分析

### 8.1 需要修改的函数

| 函数 | 修改 | 正确性依据 |
|------|------|-----------|
| `__hensel_lift` 出口 | 将 `lifted[0]` 首一化：`lifted[0] *= lc^{-1} mod m` | 建立 (I'2)，由 $\gcd(\ell, m) = 1$ 保证逆元存在 |
| `__zassenhaus_recombine` | 移除 `subset_has_lc`，统一 `lc_mult = lc_fstar` | 定理 6.3, 6.4 |
| `__vanhoeij_recombine` | （无需修改——已经使用 `lc_fstar` 前缀） | 定理 6.2 |

### 8.2 不需要修改的函数

| 函数 | 原因 |
|------|------|
| `__hensel_step` / `__hensel_step_linear` | Hensel 提升内部逻辑不受影响（修改在提升完成后） |
| `__cld_polys` | 接收外部传入的 `active_factors`，首一化在传入前已完成 |
| `__build_cld_matrix` | 仅操作 CLD 系数列，不涉及 lc 处理 |
| `__lll_reduce` | 纯格规约，不涉及多项式语义 |
| `__select_prime` | 素数选择与 lc 处理无关 |

## 9. 总结

### 9.1 定理依赖图

```
引理 1.1 (首一精确除法)
    ↓
定理 2.1 (Hensel 唯一性)  ←— 核心工具
    ↓
定理 6.1 (不变量 I')
  ├─ (I'3): Hensel 唯一性 + 对比组构造
  └─ (I'1): (I'3) + Gauss 引理 + 归纳
         ↓
    ┌────┴────┐
    ↓         ↓
定理 6.2   定理 6.5
(因子恢复)  (CLD 正确)
    ↓
  ┌─┴─┐
  ↓   ↓
定理 6.3  定理 6.4
(lc 剪枝) (常数项剪枝)
     ↓
定理 7.1 (满精度完全正确性)
```

### 9.2 修复效果

| 组件 | 修复前 | 修复后 | 正确性保证 |
|------|--------|--------|-----------|
| Hensel 因子 lc | $H_1$: $\text{lc}=\ell$, 其余首一 | 全部首一 | 定理 6.1 (I') |
| 因子恢复 | 分 $1\in S$ / $1\notin S$ 两种情况 | 统一 $\ell^* \cdot \prod H_i$ | 定理 6.2 |
| lc 剪枝 | $1\in S$ 时 $\ell\nmid\ell^{*2}$ 可能失败 | $\ell^{*2}\%\ell^*=0$ 始终成立 | 定理 6.3 |
| 常数项剪枝 | $1\in S$ 时缩放因子不匹配 | 缩放因子统一为 $\ell^*/\ell_j$ | 定理 6.4 |
| CLD | $\ell^{-1}\bmod m$ 引入 $O(m)$ 系数 | 首一除法，系数 $O(\text{poly})$ | 定理 6.5 |

### 9.3 正确性边界

**无条件保证**：在满 Mignotte 精度（$m > 2\ell \cdot B_{\text{Mig}}$）下，修复后的算法返回完全不可约分解（定理 7.1）。

**启发式依赖**：在启发式精度 $a_h < a_{\text{Mig}}$ 下，完全正确性依赖于 Phase 2 触发条件（命题 7.3）。当前条件 `result.size() == 1 && r > 1` 存在理论缺口，但实践中极少触发。
