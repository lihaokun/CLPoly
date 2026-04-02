# ZZ[x] 单变量因式分解：正确性分析

## §1 算法流程概述

输入：$f \in \mathbb{Z}[x]$，无平方、本原、$\deg(f) \ge 2$、$\text{lc}(f) > 0$。

```
factorize(f):
  1. __select_prime(f)           → (p, {h₁,...,hᵣ} ⊂ Fₚ[x], irreducible?)
  2. __lll_factorize(f, {hᵢ}, p)
     2a. Phase 1: __hensel_lift(f, {hᵢ}, p, a_h)    → {H₁,...,Hᵣ} mod m_h
         __vanhoeij_recombine(f, {Hᵢ}, m_h)          → result₁
     2b. Phase 2 (条件触发):
         __hensel_lift(f, {hᵢ}, p)                   → {H₁,...,Hᵣ} mod m_mig
         __vanhoeij_recombine(f, {Hᵢ}, m_mig)        → result₂
  3. 返回 result
```

记号约定：
- $s$：$f$ 在 $\mathbb{Z}[x]$ 中的不可约因子个数
- $r$：$f \bmod p$ 在 $\mathbb{F}_p[x]$ 中的不可约因子个数
- $B_{\text{Mig}}(f) = \binom{n}{\lfloor n/2 \rfloor} \cdot \|f\|_2$：Mignotte 界
- $a_h$：启发式精度（FLINT 公式）
- $a_{\text{mig}}$：Mignotte 精度，满足 $p^{a_{\text{mig}}} > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$
- $m_h = p^{2^{\lceil \log_2 a_h \rceil}}$：Phase 1 实际模数（二次提升取整）
- $m_{\text{mig}} = p^{2^{\lceil \log_2 a_{\text{mig}} \rceil}}$：Phase 2 实际模数

---

## §2 各阶段正确性

### §2.1 素数选择 (`__select_prime`)

**前置条件**：$f$ 无平方、本原、$\deg(f) \ge 2$。

**后置条件**：返回素数 $p$ 和 $f \bmod p$ 的不可约分解 $\{h_1, \ldots, h_r\}$，满足：
1. $p \nmid \text{lc}(f)$
2. $\deg(f \bmod p) = \deg(f)$（首项系数不消失）
3. $f \bmod p$ 无平方（$\gcd(f \bmod p,\, f' \bmod p) = 1$）
4. $h_i$ 两两互素（由无平方性保证）
5. $f \equiv c \cdot h_1 \cdots h_r \pmod{p}$，其中 $c = \text{lc}(f) \bmod p$

**命题 2.1.1**（$r \ge s$）：$f$ 的每个 $\mathbb{Z}[x]$-不可约因子 $g_j$ 对应 $f \bmod p$ 的一个或多个 $\mathbb{F}_p[x]$-不可约因子。因此 $r \ge s$。

*证明*：$f = g_1 \cdots g_s$ 在 $\mathbb{Z}[x]$ 中。模 $p$ 后 $f \equiv \prod_j (g_j \bmod p)$。每个 $g_j \bmod p$ 在 $\mathbb{F}_p[x]$ 中至少有一个不可约因子。由 $f \bmod p$ 无平方，各 $g_j \bmod p$ 之间无公共因子（否则 $f \bmod p$ 有重因子）。因此 $r = \sum_j r_j \ge s$，其中 $r_j$ 是 $g_j \bmod p$ 的不可约因子个数。$\square$

**命题 2.1.2**（不可约检测）：若某素数 $p$ 使得 $r = 1$，则 $f$ 在 $\mathbb{Z}[x]$ 中不可约。

*证明*：$r \ge s \ge 1$，$r = 1$ 蕴含 $s = 1$。$\square$

**实现说明**：`__select_prime` 至多尝试 3 个合法素数，取 $r$ 最小者。若任一素数给出 $r = 1$，立即返回不可约。

### §2.2 Hensel 提升 (`__hensel_lift`)

**前置条件**：
- $f \equiv \text{lc}(f) \cdot h_1 \cdots h_r \pmod{p}$
- $h_i$ 两两互素于 $\mathbb{F}_p[x]$

**后置条件**：返回 $\{H_1, \ldots, H_r\}$ 和模数 $m$，满足：
1. $f \equiv H_1 \cdots H_r \pmod{m}$
2. $H_i \equiv h_i \pmod{p}$（对 $i \ge 2$）；$H_1 \equiv \text{lc}(f) \cdot h_1 \pmod{p}$
3. $\deg(H_i) = \deg(h_i)$
4. $m \ge p^{a_{\text{target}}}$（若指定）或 $m > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$（默认）

**定理 2.2（Hensel 引理）**：由于 $h_i$ 两两互素，提升唯一确定。二次提升每步 $m \to m^2$，$O(\log a)$ 步达到精度 $p^a$。

*正确性*：这是 Hensel 引理的标准推论，此处不再证明。CLPoly 实现使用二叉树结构的多因子 Hensel 提升（GCL §6.4），每层递归处理 $(g, h)$ 对。

### §2.3 Mignotte 界与对称模恢复

**定理 2.3（Mignotte 界）**：设 $f \in \mathbb{Z}[x]$，$\deg(f) = n$，$g | f$ 在 $\mathbb{Z}[x]$ 中。则 $g$ 的每个系数满足：

$$|g_k| \le \binom{n}{\lfloor n/2 \rfloor} \cdot \|f\|_2 =: B_{\text{Mig}}(f)$$

**推论 2.3.1（对称模恢复）**：设 $m > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$。对 $f$ 的任意真因子 $g$，设 $S \subseteq \{1,\ldots,r\}$ 使得 $\text{lc}(f) \cdot g \equiv \text{lc}(f) \cdot \prod_{i \in S} H_i \pmod{m}$。则：

$$\text{sym}_m\!\bigl(\text{lc}(f) \cdot \textstyle\prod_{i \in S} H_i \bmod m\bigr) = \text{lc}(f) \cdot g$$

即对称约化精确恢复 $\text{lc}(f) \cdot g$。

*证明*：$\text{lc}(f) \cdot g$ 的每个系数绝对值 $\le |\text{lc}(f)| \cdot B_{\text{Mig}} < m/2$。因此 $\text{lc}(f) \cdot g$ 的每个系数在 $(-m/2, m/2]$ 中，对称约化无损。$\square$

**定义**：称模数 $m$ 为 $f$ 的**充分精度**，当且仅当 $m > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$。

### §2.4 重组 (`__zassenhaus_recombine` / `__vanhoeij_recombine`)

重组算法接收 $f$、提升因子 $\{H_i\}$、模数 $m$，返回 $f$ 的因子列表。

两种重组策略共享相同的核心步骤：
1. 从 $\{H_1, \ldots, H_r\}$ 中选取子集 $S$
2. 计算候选因子 $\tilde{g} = \text{pp}\bigl(\text{sym}_m(\text{lc}(f^*) \cdot \prod_{i \in S} H_i \bmod m)\bigr)$
3. 试除：$\tilde{g} \mid f^*$ 在 $\mathbb{Z}[x]$ 中？

**命题 2.4.1（试除正确性）**：试除永不引入假阳性。即：若 $\tilde{g} \mid f^*$ 在 $\mathbb{Z}[x]$ 中，则 $\tilde{g}$ 是 $f^*$ 的一个真因子（$\tilde{g}$ 是一个或多个 $\mathbb{Z}[x]$-不可约因子的乘积）。

*证明*：$\tilde{g} \mid f^*$ 是在 $\mathbb{Z}[x]$ 中精确除法（余式为零）。$f^*$ 无平方（$f$ 无平方且 $f^* \mid f$），因此 $\tilde{g}$ 是 $f^*$ 的某些不可约因子的乘积。$\square$

**命题 2.4.2（试除可能引入假阴性）**：当 $m$ 不是充分精度时，对称模恢复可能不精确（推论 2.3.1 不适用），导致正确子集 $S$ 的候选 $\tilde{g}$ 与真因子 $g$ 不同，从而 $\tilde{g} \nmid f^*$（试除失败）。

*证明*：反证法。设 $m < 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}$。$\text{lc}(f) \cdot g$ 的某些系数可能超过 $m/2$，对称约化截断，$\tilde{g} \ne \text{pp}(\text{lc}(f) \cdot g)$。一个随机多项式整除 $f^*$ 的概率为零，因此 $\tilde{g}$ 几乎必然不整除 $f^*$。$\square$

**定理 2.4.3（充分精度下重组的完全性）**：当 $m$ 为充分精度时，Zassenhaus 或 van Hoeij 重组返回恰好 $s$ 个因子，且每个因子均为 $f$ 的 $\mathbb{Z}[x]$-不可约因子。

*证明概要*：
- **Zassenhaus**：穷举子集，子集大小从 1 到 $r/2$。由推论 2.3.1，每个真因子对应的子集 $S$ 的候选 $\tilde{g}$ 精确等于 $\text{pp}(\text{lc} \cdot g)$，试除必定成功。
- **van Hoeij**：CLD 矩阵在充分精度下正确计算，LLL 约化保证在多项式时间内找到所有真因子子集（van Hoeij 2002 定理）。安全网 Zassenhaus 在 LLL 列耗尽时接管，此时同样在充分精度下工作。$\square$

---

## §3 `__lll_factorize` 的正确性分析

### §3.1 当前实现

```
__lll_factorize(f, {hᵢ}, p):
  a_h = __heuristic_starting_precision(f, r, p)   // a_h = min(a_mig, a_h_FLINT)
  Phase 1:
    {Hᵢ}, m_h = __hensel_lift(f, {hᵢ}, p, a_h)
    result = __vanhoeij_recombine(f, {Hᵢ}, m_h)
  Phase 2 (当前触发条件):
    if result.size() == 1 AND r > 1:
      {Hᵢ}, m_mig = __hensel_lift(f, {hᵢ}, p)     // 完整 Mignotte 精度
      result = __vanhoeij_recombine(f, {Hᵢ}, m_mig)
  return result
```

### §3.2 启发式精度

`__heuristic_starting_precision` 返回 $a_h = \min(a_{\text{mig}},\, a_{h,\text{FLINT}})$，其中：

$$a_{h,\text{FLINT}} = \biggl\lceil \frac{(2.5r + \text{min\_b}) \ln 2}{\ln p} + \frac{\ln(N+1)}{2 \ln p} \biggr\rceil$$

两种情况：
- **$a_{h,\text{FLINT}} \ge a_{\text{mig}}$**：$a_h = a_{\text{mig}}$，Phase 1 即充分精度，重组保证正确（定理 2.4.3）。
- **$a_{h,\text{FLINT}} < a_{\text{mig}}$**：$a_h = a_{h,\text{FLINT}} < a_{\text{mig}}$，Phase 1 精度不足，重组可能不完全。

### §3.3 精度不足时的返回值分析

**引理 3.3.1**：重组（Zassenhaus 或 van Hoeij）的返回值 $\text{result}$ 始终满足：

$$\text{content}(\text{result}) \cdot \prod_{g \in \text{result}} g = f$$

即返回因子的乘积等于 $f$（不论精度是否充分）。

*证明*：
- Zassenhaus：每找到一个因子 $\tilde{g}$，$f^* \leftarrow \text{pp}(f^* / \tilde{g})$。循环结束时，剩余 $f^*$ 加入结果。因此 $\prod g_i = f$（到常数倍）。
- van Hoeij：同理，每找到一个因子后更新 $f^*$，最终剩余 $f^*$ 加入结果。$\square$

**引理 3.3.2**：返回的每个因子都是 $f$ 的 $\mathbb{Z}[x]$-因子（一个或多个不可约因子的乘积）。

*证明*：直接由命题 2.4.1（试除正确性）。$\square$

**引理 3.3.3**（$|\text{result}| \le s$）：返回因子个数不超过 $f$ 的不可约因子个数 $s$。

*证明*：每个返回因子 $g_i$ 是一个或多个不可约因子的乘积（引理 3.3.2），且由引理 3.3.1，$\prod g_i = f$。因此 $\{g_i\}$ 对应 $f$ 的不可约因子集 $\{g_1^*, \ldots, g_s^*\}$ 的一个**划分的粗化**（coarsening）。划分的类数 $\le s$。$\square$

**定义**：称返回结果 $\text{result}$ **完全**，当且仅当 $|\text{result}| = s$（每个因子均不可约）。

**引理 3.3.4**（不完全当且仅当存在合成因子）：$|\text{result}| < s$ 当且仅当至少一个返回因子是两个或以上不可约因子的乘积（合成因子）。

*证明*：引理 3.3.3 的直接推论。$\square$

### §3.4 Phase 2 触发条件分析

**目标**：Phase 2 触发条件应满足：
1. **完备性（Soundness）**：每当 Phase 1 返回不完全结果时，Phase 2 触发。
2. **无害性**：Phase 2 触发但结果已正确时，Phase 2 仍返回正确结果（仅浪费时间）。

#### §3.4.1 当前条件 `result.size() == 1 && r > 1` 的缺陷

**命题 3.4.1**（当前条件不完备）：存在输入使得 Phase 1 返回不完全结果（$1 < |\text{result}| < s$），但当前条件不触发。

*构造*：取 $f = (x+1)(x-800)(x+800)(96000x + 38368001)$，$s = 4$。
- `__select_prime` 选择 $p = 7$，$r = 4$。
- Phase 1：$a_h = 6$，$m_h = 7^{2^3} = 7^8 = 5764801$。
- 对称模恢复需要 $m > 2 \cdot 96000 \cdot B_{\text{Mig}}$。由 Mignotte 界计算，$a_{\text{mig}} \approx 90/\log_2 7 \approx 32$。
- Phase 1 精度远不足（$a_h = 6 \ll 32 = a_{\text{mig}}$）。
- van Hoeij 返回 2 因子：$(x+1)$ 和 $(96000x^3 + \cdots)$。
- 当前触发条件：$|\text{result}| = 2 \ne 1$，不触发 Phase 2。
- 返回 2 因子（不完全：三次因子可进一步分解为 3 个不可约因子）。$\square$

#### §3.4.2 修正条件 `result.size() < r`

**命题 3.4.2**（修正条件的完备性）：若 Phase 1 返回不完全结果（$|\text{result}| < s$），则 $|\text{result}| < r$。

*证明*：$|\text{result}| < s \le r$（由命题 2.1.1，$r \ge s$）。因此 $|\text{result}| < r$。$\square$

**命题 3.4.3**（修正条件可能过度触发）：存在输入使得 Phase 1 返回完全结果（$|\text{result}| = s$），但 $s < r$，导致 Phase 2 不必要地触发。

*构造*：取 $f = x^4 + 1$（在 $\mathbb{Z}[x]$ 中不可约，$s = 1$）。
- $f \bmod 3 = (x^2 + x + 2)(x^2 + 2x + 2)$，$r = 2$。
- `__select_prime` 尝试 3 个素数均给出 $r \ge 2$（$x^4+1$ 对所有素数 $p$ 分裂）。
- Phase 1（充分精度）：正确返回 1 因子（$f$ 不可约）。
- $|\text{result}| = 1 < 2 = r$：Phase 2 触发。
- Phase 2 同样返回 1 因子。结果正确，仅浪费 Phase 2 的时间。$\square$

**命题 3.4.4**（修正条件的无害性）：当 $|\text{result}| < r$ 但 Phase 1 已正确时，Phase 2 返回相同的正确结果。

*证明*：Phase 2 使用充分精度 $m_{\text{mig}} > 2 \cdot |\text{lc}| \cdot B_{\text{Mig}}$。由定理 2.4.3，重组返回恰好 $s$ 个不可约因子。若 Phase 1 已正确（$|\text{result}| = s$），Phase 2 返回相同的 $s$ 个因子。$\square$

**命题 3.4.5**（$|\text{result}| = r$ 蕴含完全）：若 $|\text{result}| = r$，则结果完全且正确。

*证明*：由引理 3.3.3，$|\text{result}| \le s \le r$。$|\text{result}| = r$ 蕴含 $s = r$，且 $|\text{result}| = s$。每个返回因子恰对应一个不可约因子。$\square$

### §3.5 正确性定理

**定理 3.5（修正后 `__lll_factorize` 的正确性）**：将 Phase 2 触发条件从 `result.size() == 1 && r > 1` 改为 `result.size() < r` 后，`__lll_factorize` 对所有合法输入返回完全正确的不可约分解。

*证明*：分两种情况。

**情况 1**：Phase 1 返回完全结果（$|\text{result}| = s$）。
- 若 $s = r$：不触发 Phase 2，返回正确结果（命题 3.4.5）。
- 若 $s < r$：触发 Phase 2，Phase 2 返回正确结果（命题 3.4.4）。

**情况 2**：Phase 1 返回不完全结果（$|\text{result}| < s$）。
- 由命题 3.4.2，$|\text{result}| < r$，Phase 2 触发。
- Phase 2 使用充分精度 $m_{\text{mig}}$，由定理 2.4.3 返回 $s$ 个不可约因子。$\square$

### §3.6 过度触发的影响分析

修正条件 `result.size() < r` 的唯一代价是：当 $s < r$ 时 Phase 2 不必要地触发。

**常见情况**：
- $s = r$（每个 $\mathbb{Z}[x]$-不可约因子在 $\mathbb{F}_p[x]$ 中也不可约）：不触发。**这是最常见的情况**（随机素数下，大多数不可约因子 mod $p$ 仍不可约）。
- $s = 1$（$f$ 不可约）：`__select_prime` 大概率检测到某素数 $r = 1$ 并提前返回。仅当 $f$ 对所有尝试的素数都分裂时（如 $x^4 + 1$、分圆多项式），Phase 2 触发。这类多项式很少见。

**进一步优化（可选）**：若 `__heuristic_starting_precision` 同时返回 $a_h$ 和 $a_{\text{mig}}$，可加入条件 `a_h < a_mig`，避免 Phase 1 已用充分精度时触发 Phase 2：

```cpp
if ((int)result.size() < r && a_h < a_mig) { /* Phase 2 */ }
```

**命题 3.6.1**（$a_h = a_{\text{mig}}$ 时 Phase 2 不必要）：当 $a_h = a_{\text{mig}}$ 时，Phase 1 返回完全正确的分解，Phase 2 无需触发。

*证明*：
1. $a_h = \min(a_{\text{mig}},\, a_{h,\text{FLINT}})$。$a_h = a_{\text{mig}}$ 蕴含 $a_{h,\text{FLINT}} \ge a_{\text{mig}}$。
2. Phase 1 调用 `__hensel_lift(f, {hᵢ}, p, a_{\text{mig}})`。内部设 $\text{target} = p^{a_{\text{mig}}} - 1$，二次提升产生 $m = p^{2^k}$（$2^k$ 是 $\ge a_{\text{mig}}$ 的最小 2 的幂），因此 $m \ge p^{a_{\text{mig}}}$。
3. 由 `__heuristic_starting_precision` 中 $a_{\text{mig}}$ 的定义：循环 `while (pa <= target) { pa *= p; ++a_mig; }` 结束时 $p^{a_{\text{mig}}} > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$。
4. 综合 (2)(3)：$m \ge p^{a_{\text{mig}}} > 2 \cdot |\text{lc}(f)| \cdot B_{\text{Mig}}(f)$，即 $m$ 是充分精度（§2.3 定义）。
5. 由定理 2.4.3，充分精度下重组返回 $s$ 个不可约因子。Phase 1 结果已完全正确。$\square$

因此增加 `a_h < a_mig` 守卫可安全消除所有 $a_h = a_{\text{mig}}$ 场景下的过度触发。

---

## §4 `__select_prime` 的性质

### §4.1 因子数的单调性与最优素数

`__select_prime` 尝试至多 3 个合法素数，取 $r$ 最小者。

**注意**：最小 $r$ 不一定是 $s$。例如 $f = (x^2+1)(x^2+x+1)$ 在 $\mathbb{Z}[x]$ 中有 $s=2$ 个不可约因子，但 $f \bmod 2 = (x^2+1)(x^2+x+1)$ 给出 $r=2$，而 $f \bmod 13 = (x-5)(x+5)(x-4)(x+4+\cdots)$ 可能给出 $r=4$。选最小 $r$ 减少了 Hensel 提升和重组的代价。

### §4.2 `max_tries = 3` 的局限

仅尝试 3 个素数可能选不到 $r = s$ 的最优素数。但这不影响正确性：无论选哪个合法素数，Hensel 提升 + 充分精度重组（Phase 2）都保证正确。`max_tries` 仅影响性能。

---

## §5 端到端正确性总结

**定理 5.1**（端到端正确性）：修正 Phase 2 触发条件后，`factorize(f)` 对所有 $f \in \mathbb{Z}[x]$ 返回正确的完全因式分解。

*证明链*：
1. Squarefree 分解正确（标准算法）。
2. 对每个无平方因子 $f_i$：
   a. `__select_prime` 返回合法素数 + 正确的 $\mathbb{F}_p[x]$ 因子（§2.1）。
   b. `__hensel_lift` 正确提升（Hensel 引理，§2.2）。
   c. Phase 1 重组返回 $|\text{result}| \le s$ 个因子（引理 3.3.1-3.3.3）。
   d. 若 $|\text{result}| < r$，Phase 2 用充分精度重新提升+重组，返回 $s$ 个不可约因子（定理 3.5）。
   e. 若 $|\text{result}| = r$，Phase 1 结果已正确（命题 3.4.5）。
3. 乘积验证断言（debug 模式 line 1627）确认 $\text{content} \cdot \prod f_i^{e_i} = F$。$\square$

---

## §6 B2 Bug 总结

| 项目 | 内容 |
|------|------|
| **位置** | `polynomial_factorize_univar.hh:1488` |
| **当前条件** | `result.size() == 1 && r > 1` |
| **缺陷** | 不完备：当 $1 < |\text{result}| < s$ 时不触发 Phase 2 |
| **修正条件** | `result.size() < r` |
| **正确性** | 定理 3.5 |
| **副作用** | 当 $s < r$ 时过度触发 Phase 2（无害，仅浪费时间） |
| **可选优化** | 增加 `a_h < a_mig` 条件，消除 Phase 1 已充分精度时的过度触发 |
| **复现** | $f = 96000x^4 + 38464001x^3 - 61401631999x^2 - 24616960640000x - 24555520640000$ |
