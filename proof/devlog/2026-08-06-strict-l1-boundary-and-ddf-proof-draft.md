# 严格 L1 边界清理与 DDF 证明草案

## 结论

本阶段撤销所有不能证明 cpp2lean 生成执行树精化 L2 的公开入口。删除的内容包括：

- EDF/Hensel 的“候选结果通过后验规约则采用，否则调用 L2 算法”包装；
- 独立手写的 DDF/EDF/Hensel safe 实现；
- 用 `Classical.choose` 或 UFD 存在性见证冒充 L1 算法的 Pipeline 实例；
- Recombine 中把生成结果等同于选择见证的占位定理；
- SQF 在常数输入上绕开生成入口的 wrapper 及其 Pipeline 出口；
- ZZ 算术中尚未完成的占位精化定理。

因此当前 Pipeline 不宣称任何尚未直接闭合的算法级 L1 正确性。已有数学 L2
定理仍可作为目标规约，但不能作为 L1 的实现或 fallback。

## StrictDDF 自然语言证明草案

目标是直接证明 `Generated.StrictDDF.__ddf_Zp_ir` 的输出满足 L2 `DDFCorrect`，
中间不定义另一份可执行 DDF。

1. 对生成的 `make_zp`、`make_monic`、`mod`、`powmod`、`subtract_x`、GCD
   和精确除法逐一建立表示保持与多项式语义定理。每个定理都以生成函数本身为
   主语，不经由独立实现。
2. `powmod` 内循环按生成函数的良基归纳原理证明。度量是指数绝对值；在 C++
   循环条件 `e > 0` 下，除以二得到严格更小的非负指数。循环不变量是
   `result * base^e` 在模多项式下与初始目标同余，同时所有稀疏多项式保持规范形。
3. DDF 主循环按生成的
   `ddfWellFoundedMeasure fStar d = deg(fStar)+1-2*d` 归纳。循环不变量包括：
   已输出块的乘积乘剩余 `fStar` 等于输入；输出块两两互素、首一且只含指定次数
   的不可约因子；`h` 等于 `X^(p^d)` 模 `fStar`；所有表示均规范且整数转换无溢出。
4. 若 GCD 非平凡，生成代码首一化 GCD、push `(gd,d)`、精确除去 `gd`，再把
   `h` 对新余因子取模。利用 GCD 的有限域因子次数刻画保持上述不变量，并由
   `fStar` 次数下降证明生成递归调用确实走 `hdec` 分支。
5. 若 GCD 平凡，只递增 `d`；`2*d` 增加二，因此良基度量严格下降。终止分支
   `deg(fStar) < 2*d` 表明剩余非恒定因子的所有不可约因子次数相同，生成入口的
   尾部 push 正好补齐最后一块。
6. 最后展开 `__ddf_Zp_ir`，证明其初始化 `h=X, fStar=f, d=1, result=[]`
   满足循环不变量，并将数组观察映射转换为 L2 的列表规约。

这个证明只允许引用生成定义、表示层引理和 L2 数学定理；不得引用 `*Safe`、
候选验证器、L2 算法函数作为实现，或任何 fuel 版本。

## 度量

- 删除不合格代码：2057 行（本次工作树统计，提交前值）。
- 删除不合格生成模块：2 个。
- 当前严格扫描：`sorry/admit` 0 个；fallback/后验选择定义 0 个。
- 受影响 Lean 目标：9 个全部构建成功。
- StrictDDF 递归：良基递归 2 处（powmod、DDF），fuel 0 处。
