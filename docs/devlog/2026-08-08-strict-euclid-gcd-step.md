# 严格 Euclid 单步 GCD 不变量

日期：2026-08-08

## 做了什么

- 新增 `StrictEuclidRefinement.lean`。
- 从除法恒等式 `A = Q * B + R` 证明 Euclid 轮换 `(A, B) → (B, R)`
  保持规范化的数学 GCD。
- 将该代数引理直接接到生成的 C++ `_poly_divrem` RawHeap 精化定理，得到
  带有真实执行、内存表示、余式降阶和 GCD 保持性质的单步结果。
- 证明 normalization 返回长度是原始扫描的固定点。
- 分别证明短分支 `memcpy` 与长分支三字 limb 约减都保持规范系数、标准
  长度、物理容量及 heap layout。
- 定义 `RawDensePolyRep`，并证明真实除法输出的 `R` 满足下一轮 Euclid
  可直接消费的完整 raw 状态不变量。
- 逐段证明除法不改写旧 `B`：短分支检查 `memcpy` 目标区域，长分支组合
  W3 初始化、quotient loop 与 remainder loop 的逐 cell 保持定理。
- `polyDivrem_next_state` 同时返回轮换后 `A = old B` 与 `B = R` 的完整
  `RawDensePolyRep`，不再把源除数内容保持留作递归假设。

## 为什么做

SQF 调用 C++ GCD。要证明真正的 L1→L2 精化，必须从生成的除法执行逐步建立
GCD 不变量，不能调用 L2 `%`、规格 oracle 或回退实现替代 C++ 执行。

## 关键决策

- 证明只使用已由生成代码建立的除法恒等式，不重新计算 L2 余式。
- 先证明两组参数具有相同共同因子，再由互相整除得到关联性并规范化；这与
  C++ 的 `a = move(b); b = move(r)` 轮换完全对应。
- 保留 `Fact (Nat.Prime p)` 作为多项式 Euclidean domain 实例的显式前提。
- 不从 `Polynomial (ZMod p)` 的值反推 raw limb；规范性来自实际 `memcpy`
  内容关系或生成三字约减的模等式，避免丢失机器表示信息。

## 遇到的问题与解决方式

- 直接改写整除见证会误改写 GCD 参数内部的同名多项式。改用 `dvd_sub`、
  `dvd_add` 和 `dvd_mul_of_dvd_right` 构造共同因子证明，避免脆弱的全局改写。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictEuclidRefinement.lean`
- `docs/devlog/2026-08-08-strict-euclid-gcd-step.md`
