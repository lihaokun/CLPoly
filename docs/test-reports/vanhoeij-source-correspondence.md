# van Hoeij 修复的源码与严格模型核对

日期：2026-10-16。研究分支 `research/vanhoeij-lattice-certificate`。
本记录核对最终执行路径，不能替代尚在运行的最终测试。

## 执行路径

当前公共入口依次执行选素数、单模因子分支或纯 LLL 重试。纯重试
执行实际启发式、实际 Hensel、实际完整重组证书；失败则提高精度。
没有固定次数预算或终端 Zassenhaus 分支，也不再比较真因子数与
模因子数决定是否丢弃证书。

| C++ 源操作 | 严格模型 | 已检查的对应要点 |
|---|---|---|
| `__heuristic_starting_precision` 及两个精度 helper | `StrictVanHoeijHeuristic` | 浮点保护先于转换；实际 Mignotte 调用；先比较幂与目标，再检查 INT_MAX，再乘 p、自增；返回 min |
| `__vanhoeij_next_precision`、`__lll_factorize` | `StrictVanHoeijRetry` | 超过半上限时饱和；仅最大指数失败时报错；每轮重新 Hensel；仅完整证书返回 |
| `__vanhoeij_try_recombine` | `StrictVanHoeijEntry` | singleton 验证、参数选择、CLD、截断、宽度、列循环顺序相同 |
| `__vanhoeij_choose_parameters` | `StrictVanHoeijParameters` | D₀ 倍增、D₀² 上限、同步 floor 模数、残差界及所有拒绝分支 |
| 规范矩阵及 `__build_cld_matrix` | `StrictVanHoeij` | 每轮重新建立 32I；螺旋系数次序；增加零列、CLD 条目和模数行 |
| `__lll_reduce` 与缓存 GS 范数 | `StrictRecombine`、`StrictVanHoeij.lllReduceWithNorms` | 相同 QQ 规约、舍入、Lovász 比较、行交换、变换矩阵更新及最终范数 |
| `__lll_short_span_from_norms` | `StrictVanHoeij.shortSpanFromNorms` | 只删除连续超界后缀；不是普通行范数过滤 |
| `__extract_candidates` | `StrictRecombine.extractCandidates` | 比较变换矩阵前 r 列在所有保留行的坐标，按首次出现建立等价类 |
| `__vanhoeij_validate_partition` | `StrictVanHoeij.validatePartition` | 当前首项系数乘模因子、对称恢复、本原化、次数检查、精确除法；全部组验证后余项必须为 1 |
| 列增长控制器 | `StrictVanHoeijController` | 获证立即返回；全宽失败返回未完成；否则安全倍增至宽度 |
| `__factor_squarefree_primitive_ZZ` | `StrictVanHoeijFactorZZ` | 使用同一次选素数结果；单因子直接返回，其他调用纯重试 |

源码集中在 [单变量实现](../../clpoly/polynomial_factorize_univar.hh)，
对应模板在 [严格生成器目录](../../proof/cpp2lean_v2/tests)，
模型在 [Generated](../../proof/lean/CLPoly/Generated)。模板的 `--check`
检查文件可复现性；模板是人工维护的控制流翻译，不能将该检查或
B2B 通过称为通用 Clang AST 翻译正确性或任意 C++ 程序等价证明。

## 数值与运行边界

- ZZ/QQ 对应任意精度整数/有理数；GMP floor、对称约化及精确
  除法使用已有严格算术模型和已证 callee 契约。
- 因子数、请求指数和 CLD 宽度保留 Int32。正指数自增及倍增的
  signed 边界在实际分支中检查，未以 Nat 无限整数替代。
- 格维数是 r+j，不要求它也落在 signed int 内。LLL 行号、内积
  下标、保留行数组及 iota 计数均改为 size_t；r+j 的加法在转为
  size_t 后执行。由各自的 Int32 范围与实际新增列数界，实际格
  维数小于 2^32；`VanHoeijIndexWidth` 连接实际矩阵大小、USize
  加法及安全循环递增/递减。额外规约使用正剩余数倒数，避免
  unsigned 的 −1 哨兵。
- 浮点估计不用于判定因子不可约。完整正确性和成功存在性只需
  clamp 结果为正，后续可提高精度；B2B 另核对具体浮点返回值。
- 计时和 profile 字段不参与决策；move 对应同一值的转移。资源
  分配失败不在 Lean 算法成功定理的范围内。
- 最终因子按次数排序；相同次数的相对顺序不是分解正确性契约。
  数学结论验证因子乘积与各因子不可约性，不将排序内部实现当作
  新的代数假设。

## 证明依据与明确前提

由实际规范格和 CLD 构造真因子短向量，经实际 LLL 变换和 GS
后缀界保证保留，再由实际列等价类和整个分区验证推出输出因子
不可约。充分精度结式分离界依赖源输入，D 上限也是输入固定值。
实际入口和实际纯重试的成功执行均已构造，不给出重组成功或输出
不可约性 oracle。

公共定理仍明确要求规范非空源、正首项、本原条件、合法选素数/
模因子、物理次数及因子数范围、恢复目标与分离阈值可用 signed
请求指数表示，以及既有 Hensel 运行就绪契约。它不是对任意非法
物理数组、任意超出机器接口范围输入或资源不足状态的无条件定理。

主要证据是 [入口完整执行](../../proof/lean/CLPoly/Refinement/VanHoeijEntryTotal.lean)、
[实际前缀](../../proof/lean/CLPoly/Refinement/VanHoeijHeuristicComplete.lean)、
[纯重试总正确性](../../proof/lean/CLPoly/Refinement/VanHoeijRetryTotal.lean)、
[格下标范围](../../proof/lean/CLPoly/Refinement/VanHoeijIndexWidth.lean)。

固定精度兼容入口 `__vanhoeij_recombine` 在证书不足时仍有显式
Zassenhaus 回退，因为它没有提高 Hensel 精度所需的原模因子。
这个兼容入口不在上述公共纯重试调用路径内，不能用其回退证明
纯 van Hoeij 的终止性。旧 StrictFactorZZ 中的旧循环只为历史
证明/共享契约保留；当前注册、公共流水线和 B2B 使用新入口。

## 验收状态

完整 `lake build`（3612 jobs）、前缀及纯重试审计 34/34、公共
审计 340/340 已通过；格下标范围
的 6 项新增证明及扩大分离审计 212/212 已通过。最新格下标源码
的重组回归 401/401，加上兼容入口计数检查后 profiling 回归 406/406，
全套原生 27 文件通过；18 文件、3107 用例
的最新两端对照尚待完成。固定 FLINT 及原 PR 基线性能结果随后记录到
[参考报告](vanhoeij-flint-reference.md)，当前尚不可据此合并。
