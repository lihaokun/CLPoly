# van Hoeij 格构造与完成证书修正方案

日期：2026-10-15；进度更新：2026-10-16。状态：用户已确认修复。实际入口、完整分组证书、检查式启发式及纯 van Hoeij 重试的成功存在性和不可约性已连接；明确保留机器精度可表示性和既有 Hensel 契约。重试审计 34/34、公共审计 340/340 已通过。源码对应审查进一步修复格维数的 signed 窄化，6 项范围定理已编译；最新全套回归、3107 项两端对照及三版本分阶段性能验收仍在进行，**当前不可合并**。下文按时间保留研究过程，以本段和末尾进展为当前状态。

研究分支：`research/vanhoeij-lattice-certificate`，起点 `428a27fe2c22dcbeee100d1100e310d8a08f1e77`（`feature/l1-l2-refinement-codex`）。起点工作区干净。本地 Git 配置关联 lihaokun/CLPoly PR #21，在线 PR 工具未能核实。

## 1. 复现与定位

调用链：[__lll_factorize → __vanhoeij_recombine → __build_cld_matrix / __lll_reduce / __extract_candidates](../../clpoly/polynomial_factorize_univar.hh)。研究时对应行号分别为 1506、1195、943、1008、1144。

- 外层 safety net（1497、1546、1551 行）在满精度且输出数小于模因子数时，对原始整份输入重新执行 Zassenhaus。模因子需要组合时，这个条件是正常现象，不能表示 LLL 失败。
- 内层控制器（1340 行）可能已对剩余输入执行 Zassenhaus；外层仍会重算。现有 profiling 的 recombine 计时结束于外层 safety net 之前，测量时须单独计入这部分。
- 格构造忽略参数 m，实际为 `[s I | 0; C^T | I]`，其中 s=2^U_exp。旧测试明确断言这个错误结构。

最小数学输入：f=x²+1，p=13，m=169，提升因子 h₁=x−70、h₂=x+70。70²≡−1 mod 169，因此 CLD 是 x+70、x−70。首个常数系数约束应表达 `70e₁−70e₂ mod 169`。

当前构造与未截断的标准行格构造分别为：

```
当前                       标准模约束行格
[32   0  0]                [32   0   70]
[ 0  32  0]                [ 0  32  -70]
[70 -70  1]                [ 0   0  169]
```

标准格的向量为 `(s e, e C + m z)`；当前格没有表达这个同余关系。仅把末尾 1 改成 m 仍然不正确，CLD 的方向也须修正。

可直接编译的定位程序（使用项目已有 GMP 配置；不是生产修复）：

```cpp
#include <clpoly/polynomial_factorize_univar.hh>
#include <iostream>
using namespace clpoly;
int main() {
    upolynomial_<ZZ> f, h1, h2;
    f.push_back({umonomial(2), ZZ(1)});
    f.push_back({umonomial(0), ZZ(1)});
    h1.push_back({umonomial(1), ZZ(1)});
    h1.push_back({umonomial(0), ZZ(-70)});
    h2.push_back({umonomial(1), ZZ(1)});
    h2.push_back({umonomial(0), ZZ(70)});
    auto cld = __cld_polys(f, {h1,h2}, ZZ(169));
    LLLMatrix M(2, std::vector<ZZ>(2, ZZ(0)));
    M[0][0] = M[1][1] = ZZ(32);
    __build_cld_matrix(M, cld, 0, 1, ZZ(169));
    for (const auto& row : M) {
        for (const auto& v : row) std::cout << v << ' ';
        std::cout << '\n';
    }
}
```

将代码保存为 /tmp/vanhoeij_probe.cc 后，从仓库根目录执行：

```
c++ -std=c++17 -I. /tmp/vanhoeij_probe.cc clpoly/upolynomial.cc clpoly/variable.cc -lgmpxx -lgmp -o /tmp/vanhoeij_probe
/tmp/vanhoeij_probe
```

本机系统 clang++ 因 Command Line Tools 缺失无法启动；未运行该程序、完整 C++/Lean 测试或 r=24/32 基准。上面的矩阵是逐步源码推演。用户报告的约 50% 退化尚未独立复测；“历史上从未通过 CLD 成功分组”也未通过运行记录核实。

## 2. 根因与证明边界

1. **同余格错误。** 当前格始终包含 `(s e_i,0)`，它们均小于现有阈值 `(r+1)s²`。因此若正确保留所有短向量，其因子坐标空间仍包含全部单位向量，CLD 无法提供所需的分组约束。
2. **短向量完备性缺失。** 当前仅保留规约后行范数不超过 B 的行。一般不能由此推出所有长度不超过界的格向量均落在保留行张成空间内。需要 Gram–Schmidt 尾部移除等可证明的保留规则；不能直接照搬 FLINT 的 B 而省略其截断误差与缩放分析。
3. **完成证书缺失。** 精确试除证明“是因子”，不证明“不可约”。需要从实际格状态推出真因子指示向量全部被保留，再结合列等价类和完整试除证明不可约性。
4. **精度概念混淆。** Mignotte 精度保证因子系数恢复，不自动保证 LLL 已能分离组合空间。列耗尽后仍可能需要提升精度。
5. **旧算法层证明没有排除该实现。** [Recombine.lean](../../proof/lean/CLPoly/Algorithm/Recombine.lean) 中 `vanHoeij_correct`（700 行）与 `vanHoeij_algo_correct'`（1037 行）直接接收 `h_irred`。同文件 `vanHoeijMatrix` 的 CLD 方向甚至与当前 C++ 不同，且仍用单位矩阵代替模数块。

当前 safety net 的证明动机见[已有修正记录](../devlog/2026-08-17-vanhoeij-coordinate-and-safety-fix-plan.md)和[等基数不可约性记录](../../proof/docs/devlog/2026-08-17-vanhoeij-equal-cardinality-irreducible.md)。这些记录支持“用回退覆盖未认证路径”的解释；本次未重新构建全部证明，不能将静态审查表述为重新验证了整个 PR。

应区分控制程序最终回退并停止，与 van Hoeij 在足够精度下独立完成不可约分解。前者不能代替后者。

## 3. FLINT 一手实现对照

2026-10-15 读取以下 master 源码（尚未固定上游 commit，正式实现前应固定版本）：

- [factor_van_hoeij.c](https://github.com/flintlib/flint/blob/master/src/fmpz_poly_factor/factor_van_hoeij.c)
- [next_col_van_hoeij.c](https://github.com/flintlib/flint/blob/master/src/fmpz_mat/next_col_van_hoeij.c)
- [van_hoeij_check_if_solved.c](https://github.com/flintlib/flint/blob/master/src/fmpz_poly_factor/van_hoeij_check_if_solved.c)

具体差异：

| 阶段 | 当前实现 | FLINT |
|---|---|---|
| 新约束 | 新行前 r 项放 CLD，末项 1，m 不使用 | 从矩阵前 r 列取组合坐标 U，计算变换后的 CLD 列，加入 P_trunc 模数行 |
| 信息不足 | 不做界过滤 | 比较模数位数、CLD 界和实际列信息，必要时跳过 |
| LLL 后保留 | 根据各行普通范数筛选 | with-removal 并缩小基矩阵 |
| 分组 | 读取局部幺模变换的前 r 列 | 读取当前保留矩阵的前 r 列 |
| 完成 | 接受可试除分组，外层按数量回退 | 检查列分组；单一组判不可约，否则验证完整候选分解 |
| 未完成 | 列数超界回退 | 增加 CLD 覆盖，耗尽后继续 Hensel 提升、指数翻倍 |

对上述同一输入，FLINT 的初始 sI 给出两个单元素候选；x±70 均不能整除 x²+1，故检查不成功。m=169 的位数不足以通过其首个信息量检查（bit_r 至少 20），因此不会直接加入 CLD 列，而是继续提升。当前实现则直接加入错误的 `[70,-70,1]` 行，列预算耗尽便回退。这里是源码路径对照，不是双实现运行 trace；动态对照仍是实现前待办。

[旧调研](../research/vanhoeij-lll-research.md) 的矩阵示意、removal 解释及“按局部变换矩阵分组”的推断不可继续作为依据。

## 4. 推荐修正方案与验收

### 推荐路径

保留当前 PR 的保障机制，在研究分支实现有证书的 van Hoeij，再合回原分支。不要先删 safety net，也不要把“只有返回一个因子才回退”作为修复：多个可约分组同样可能漏分解。

1. 修订设计：先用未截断的 `[sI | C; 0 | mI]` 建立准确语义，给出真因子 CLD 界及短向量界。实际高效实现可采用有误差界的截断版本；保留每轮规范重建的坐标约定，增量优化另作证明。
2. 改造保留空间：证明真因子指示向量在投影空间 W 中；LLL 行操作保持格，安全移除保持全部相关短向量。不能只修改 `__build_cld_matrix`。
3. 引入可检查的完成状态。设列等价类为 B₁,…,Bₖ：W 中任意向量在同一类内坐标相等。因真因子 indicator 属于 W，同一类不能跨越两个真因子的支持；若每个类的模乘积都恢复为非平凡整数因子、支持匹配且完整乘积通过验证，则每类恰为一个真因子。该论证不要求 k=r，也不必先证明 W 精确等于真因子空间。
4. 未获证书时增加约束/精度。证明足够精度可分离所需空间，并给出精度增长的终止论证；若保留工程回退，应返回明确的未认证状态。小规模残余可走 Zassenhaus，但最终结果必须覆盖已提取分组的不可约性，不能只认证剩余项。
5. `__lll_factorize` 依据证书完成状态决定返回、继续提升或回退，替换 `result_count < r`。将所有回退计入 profiling。
6. 同步生成器、严格 L1、refinement、B2B 与 axiom audit，从具体运行状态推出证书，不以输出不可约性为新增前提。

### 验收条件

- 非平凡分组：至少两个整数不可约因子各分裂为多个模因子，CLD 路径确实成功，Zassenhaus 调用数为零。
- 整数不可约但模可约：通过保留空间证书判定不可约，不枚举子集。
- 非首一、多轮 CLD、提高精度、真短向量是多个基行组合等案例。
- 每例验证确切因子/不可约性，不能仅检查乘积；记录初始模因子数与路径计数。
- r=24、32 固定输入、素数、编译参数，对比原基线、当前 safety net、修复版；拆分 Hensel、LLL、内外回退耗时，重复测量。
- C++ 回归、参考实现动态对照、生成一致性、B2B、Lean 构建及 axiom audit 全通过后合回。

当前待办：恢复可运行工具链、固定 FLINT 版本并完成同输入动态对照、确认方案后细化证书接口与数学界。此文不是已完成修复或性能报告。

## 5. 确认后的实现进展

已找到并使用本机 Nix clang 19.1.7、GMP、Boost 与 Lean 4.29，工具链问题不再阻挡工作。上文的“待确认”和“尚未运行”描述保留为研究时点记录。

- C++ 已改为正确的模约束行格，加入显式 CLD 界与共同截断误差界，只移除超界的 GS 后缀；主控制器依据完整分区验证结果返回，不再调用数量判据。
- 为保持控制器有限终止，当前采用至多三次精度尝试，最终未完成时执行一次完整精度 Zassenhaus。这是明确的工程回退预算，不能宣称纯 van Hoeij 的多项式终止性已完成证明。
- 新增 [CLD 界证明](../../proof/lean/CLPoly/Algorithm/VanHoeijBound.lean)和[完成证书基础引理](../../proof/lean/CLPoly/Refinement/VanHoeijCertificate.lean)，均已构建通过；审计仅含标准 Lean/Mathlib 公理。
- 原生重组回归 382/382；既有因式分解测试在临时缩窄头文件包含后 256/256。
- 旧数量 helper 暂为旧 B2B 接口保留，当前生产控制器不调用它；旧 B2B 通过不能作为新控制器的验证证据。
- 继续迁移时新增独立 [StrictVanHoeij](../../proof/lean/CLPoly/Generated/StrictVanHoeij.lean)，覆盖矩阵追加、实际保留行数组和 LLL 范数导出。安全保留证明已连接实际输出，新增相关 B2B 190/190；完整分组与顶层尚未切换。详见[桥接日志](../../proof/docs/devlog/2026-10-15-vanhoeij-strict-helper-bridge.md)。
- 已进一步证明实际候选组在任意初始短整数格向量上坐标相等，并给出实际 CLD 的模环语义，去掉旧桥接中不适用于一般 Hensel 因子的整数整除前提。尚需构造真因子对应的实际短向量并连接候选整数验证，详见[分组坐标日志](../../proof/docs/devlog/2026-10-15-vanhoeij-physical-group-coordinates.md)。
- 规范格模块已从实际初始化/追加列证明完整矩阵公式和满秩，并由整数除余推导共同截断残差，组合出支持 indicator 与负进位对应的短格向量。仍须由实际 Hensel/CLD/截断执行产生其支持、同余和系数对应前提，详见[规范格日志](../../proof/docs/devlog/2026-10-15-vanhoeij-canonical-lattice.md)。
- 实际截断 helper、严格循环和规范化后的取系数桥接已实现；新增 62 项 B2B 通过，原 190 项 B2B 与 382 项原生回归通过。该系数定理待与顶层流程组合，详见[截断执行日志](../../proof/docs/devlog/2026-10-15-vanhoeij-truncation-execution.md)。

2026-10-16 更新：上述格、CLD、短向量保留、实际完整分组验证已连接为
新单次尝试和新重试入口的正确性定理。新重试定理覆盖完成与 Zassenhaus
回退两支，同一实际返回数组满足乘积还原及不可约性；输出计数、输出
次数、回退精度不再独立假设。它仍要求源输入条件、正确的素数选择、
算术配置和底层递归 Hensel 执行就绪条件。完整重试 B2B 27/27 通过。
详见[重试正确性日志](../../proof/docs/devlog/2026-10-16-vanhoeij-retry-correctness.md)。

外层 squarefree FactorZZ 及公开 Pipeline 已迁移到新重试路径，删除了
不再调用的旧数量 helper、其四项原生断言和独立 B2B 接口。新的原生
回归为 378/378（比此前少的四项仅是已删除 helper 的断言）。公开
Pipeline、最终公理审计及 B2B 驱动构建通过，共 3615 jobs；新增链条
21 项公理审计仅有标准公理。详见[外层迁移日志](../../proof/docs/devlog/2026-10-16-vanhoeij-factorzz-migration.md)。

外层及重试 B2B 50/50 全部通过。随后修正无关实根头文件的定义顺序，
完整库构建及仓库全部 27 个原生测试文件通过，原样因式分解测试为
256/256；不再依赖缩窄头文件的临时测试副本。

充分精度方向新增十项结式基础定理和五项实际 LLL 保留 GS 范数界，
十五项构建与公理审计通过。尚须连接物理基行坐标、截断残差和实际
分组验证的执行完备性；参数策略仍未改动。详见[分离基础日志](../../proof/docs/devlog/2026-10-16-vanhoeij-separation-foundations.md)。

随后已通过实际 Gram 等式连接保留行平方范数与物理整数坐标上界，
新增四项定理；合并十九项公理审计通过。待办进一步缩小到跨精度的
统一参数界、实际截断残差、局部 CLD 支持及成功验证的完备性。

进一步补齐任意有符号组合的截断残差界，并接到同一次实际 LLL 的
transform、物理坐标和保留扫描；新增八项定理。实际参数选择又补四项
界，给定加入列数 J≤n 时，平方 cutoff 有输入统一上界。合并分离
审计 31/31 通过。残差界仍含 D，且完整次数覆盖、实际 J≤n 和局部
支持条件仍待连接，不能据此宣布充分精度终止。详见[残差日志](../../proof/docs/devlog/2026-10-16-vanhoeij-truncated-row-residual.md)。

完整列遍历与宽度扫描、全系数残差多项式及局部 CLD 整除关系又新增
23 项证明，合并分离审计 54/54 通过。残差高度界已连接同一次实际
LLL 保留行及已有 CLD 证书，覆盖截断后消失的高次项。实际 CLD 的
输入次数界、Hensel 局部条件、充分精度策略和分组验证完备性仍待完成。
详见[完整残差日志](../../proof/docs/devlog/2026-10-16-vanhoeij-residual-polynomial.md)。

实际 Hensel/CLD 执行现已导出上述 CLD 次数界、局部因子的整除、
两两互素、导数互素和正次数，统一接入分离输入证书，新增十项定理；
分离审计 64/64 通过。剩余需组合统一输入界和真实因子支持、证明
保留行在支持上恒定，并完成实际分组完备性及精度策略。
详见[执行前提日志](../../proof/docs/devlog/2026-10-16-vanhoeij-separation-input.md)。

进一步新增 33 项定理：统一实际维数、cutoff 和坐标界，给出含 D 的
显式结式阈值，并证明超过阈值时实际保留行在真实 Hensel 支持上恒定。
实际分组的闭包与无故障返回已证明，组合目标 3580 jobs 通过，合并
分离公理审计 97/97 通过且仅含标准公理；不再假设 extractCandidates
成功运行。仍须消去保留行非空前提、连接完整
支持与组的对应及恢复验证成功，并解决 D/精度策略的阈值可达性。
详见[充分精度分组日志](../../proof/docs/devlog/2026-10-16-vanhoeij-sufficient-precision-grouping.md)。

随后新增 17 项定理：由实际 Hensel 短向量消去保留行非空假设；用同一
支持的短向量与结式分离，得到实际输出组恰为完整真实支持；再从实际
收集顺序证明索引无重复，将组的 SelectedProductMod 精确连接到该
Hensel 支持乘积。最终组合目标 3600 jobs 通过，合并公理审计 114/114
通过且仅含标准公理。剩余重点是全部组的
对应、逐组恢复验证必成功，以及能达到显式阈值的实际精度策略。
详见[精确支持乘积日志](../../proof/docs/devlog/2026-10-16-vanhoeij-exact-support-product.md)。

恢复验证方向再新增 19 项定理，覆盖实际 Int64 次数、非首一恢复、
本原化、精确除法和完整验证循环；同一结式阈值已证明覆盖恢复精度。
给定与所有实际组对齐的精确数学因子目标时，实际 validatePartition
入口必成功，目标构建 3597 jobs 通过。目标记录没有不可约性或成功
执行字段，但仍须从真实分组构造全部目标，不能据此宣布端到端终止。
详见[恢复验证日志](../../proof/docs/devlog/2026-10-16-vanhoeij-validation-completeness.md)。

恢复验证合并审计已完成：133/133，仅标准公理，3608 jobs。再新增
6 项实际类覆盖与支持唯一性定理，证明每组有物理来源、同列类唯一，
以及逐因子精确组覆盖全部实际组；最终模块目标 3601 jobs 通过。
尚须实例化数学因子族并构造完整有序目标，不能将这些中间定理视为
终止性完成。详见[组覆盖日志](../../proof/docs/devlog/2026-10-16-vanhoeij-family-coverage.md)。

再新增 8 项组合证明，已从输入本原、正首项和好素数条件构造无重复
真实因子族，证明因子到实际组的双射并构造完整有序目标。充分精度
下，同一次实际 Hensel/LLL 的组提取和 validatePartition 必成功，
不再要求完整目标数组、输出不可约性或验证成功为前提。最终目标
3614 jobs 通过；合并公理审计 147/147，仅标准公理，3616 jobs。
现在剩余的主要证明义务是把全宽单次成功接入实际
列循环，并修复随 m 增大的 D 与重试控制器，使显式阈值确实可达。
详见[充分精度验证日志](../../proof/docs/devlog/2026-10-16-vanhoeij-sufficient-validation.md)。

精度策略已采用输入固定的截断上界：令原高度倍增除数为 D₀，
把位数截断产生的除数限制为 D≤D₀²。原先的快速位数规则在未触及
上界时保持；上界以后有效模数可以继续增长。参数选择最终成功、
含 D 的结式阈值单调性、统一输入阈值及充分素数幂存在性均已证明；
合并审计 154/154、3654 jobs。原生重组 383/383，完整原生回归
27 个测试文件通过。B2B 与本轮性能记录见
[截断上界日志](../../proof/docs/devlog/2026-10-16-vanhoeij-bounded-divisor.md)。
此修改不使现有三次重试自动成为纯 van Hoeij 终止证明。

实际全宽轮次与正目标列循环已连接到统一充分精度门槛：正常返回必
完成，且最终组合定理不保留“全宽成功”oracle。新增 4 项定理，合并
审计 158/158、3625 jobs。五组新参数 B2B 1738/1738 已通过。
零宽度入口、内部执行完备性与真实外层重试可达性仍待完成，详见
[列循环记录](../../proof/docs/devlog/2026-10-16-vanhoeij-column-loop-completion.md)。

零列宽的恢复验证义务也已接通：实际截断宽度为零时，充分精度迫使
每个真实因子的实际 Hensel 支持为单元素；从源输入构造因子族并与
实际 singletonGroups 建立双射，进而证明实际完整验证成功。新增
7 项定理、4 模块 374 行，组合目标 3617 jobs，合并审计 165/165、
3631 jobs。该分支未假设任何未执行的 LLL。下一步将正目标/零宽度
两条链统一接到真实入口，消去 trace 条件，并处理内部错误不可达
与外层重试精度可达性。详见
[零列宽记录](../../proof/docs/devlog/2026-10-16-vanhoeij-zero-width-validation.md)。

充分精度结论现已统一接入真实 tryRecombine：无需调用者给出 CLD、
截断、LLL、中间完整目标或 complete 条件，即可从入口正常返回推出
complete、乘积精确及所有实际输出不可约。新增 3 项定理、185 行，
合并审计 168/168、3643 jobs。保留的正常返回前提仍须通过内部执行
完备性消去，外层指数增长与机器范围义务未完成。详见
[入口记录](../../proof/docs/devlog/2026-10-16-vanhoeij-entry-sufficient.md)。

LLL 内部执行义务已完成：从非空、方阵且满秩的实际输入构造初始化、
行消元、Lovász 推进/交换、主循环和短行收集的正常返回，最终实际
lllReduceWithNorms 无需成功执行假设。新增 14 项定理、357 行，合并
审计 182/182、3648 jobs。主循环沿用真实行列式/索引下降，不使用
燃料或回退。仍须消去 CLD、构格和整个分区验证的错误分支，再接
外层精度增长。详见[LLL 执行记录](../../proof/docs/devlog/2026-10-16-vanhoeij-lll-total-execution.md)。

固定参考版本已落实为本机 FLINT 3.1.2，源码包与实际库均记录 SHA-256。
同 p=97、同模因子的直接 van Hoeij 结果对照 18/18 通过，覆盖首一、
非首一二次族及高次 Eisenstein 族的 r=24/32；CLPoly 均完整获证且零
回退。内部格状态仅作源码核对，详见[参考报告](../test-reports/vanhoeij-flint-reference.md)。

实际重组入口的正常返回前提现已消去：从同一次 Hensel 输出与源输入
条件，构造 CLD（包含空余式断言）、截断、规范格、LLL、任意合法
候选分区验证及列数循环的实际返回。达到统一输入阈值时存在实际
入口返回，complete=true、乘积精确且各输出不可约。新增 24 项定理、
8 模块 755 行，合并审计 206/206、3656 jobs，仅标准公理；严格边界、
登记和补丁检查通过。没有改动 C++，沿用已通过的运行回归。
外层固定三次尝试与机器范围仍须处理，不能将该入口结论当成整个
重试控制器的成功终止。详见
[入口完整执行记录](../../proof/docs/devlog/2026-10-16-vanhoeij-entry-total-execution.md)。

**仍须完成**：完整生成对应关系与最终审计；纯 van Hoeij 精度可达性与终止；
性能验收。当前控制流模型由手写模板生成，
不能以 B2B 通过声称 AST 翻译或普遍实现等价；也不能用旧模型的构建成功
代替这些义务。

## 外层重试实施决定（2026-10-16）

入口充分精度成功已完成后，下一修改取消三次预算及其 Zassenhaus
终端回退。精度正常增长取 max(full,max(1,2*current))；超过半整数
上限时饱和至 INT_MAX，仅在最大指数仍未获证时报告溢出。不会执行
溢出的乘法，也不靠回退取得分解结果。严格模型使用实际指数距机器
上限的下降度量；成功存在性另由输入阈值及可表示精度条件证明。
同步迁移成功来源、公共 FactorZZCorrect、B2B 与原生测试，完成前
不以旧三次控制器的通过结果作为新实现验收。

上述生产控制器和严格模型已实现，成功来源仅有实际完整证书。
从正请求指数开始，在 T<p^INT_MAX 的明确可表示条件及既有 Hensel
就绪契约下，实际纯重试必成功；公共包装已连接精确分解与不可约性，
仍保留实际启发式前缀和正起始指数条件。新审计 25/25、3667 jobs；
原生 390/390，完整 27 测试文件通过，相关 245/245 条 B2B 已通过。
前缀核对另发现 fullNat 的模板 UInt32 范围与 C++ signed int 不一致，
必须继续修正计数溢出和浮点转换对应，再完成整体验收和性能测试。
详见[纯重试记录](../../proof/docs/devlog/2026-10-16-vanhoeij-certificate-only-retry.md)。

## 检查式启发式前缀（2026-10-16）

实际前缀已改为先检查浮点转换、在循环内检查 signed 指数自增，
项数直接转 double。新模板保留真实分支顺序，并显式使用有符号
浮点转换。已构造实际 Mignotte/完整指数/启发式的正常正返回，
`concrete_lllFactorize_succeeds` 已编译通过，不再以启发式成功或
正起始指数作假设。保留两项输入阈值可表示性和既有 Hensel 契约。
最新原生重组 401/401；扩大审计、公共流水线、最终全套回归及
最新两端 B2B 尚在收尾，源对应审查和安静环境性能验收仍未完成。
详见[前缀记录](../../proof/docs/devlog/2026-10-16-vanhoeij-checked-heuristic.md)。

前缀/纯重试审计已完成 34/34，公共流水线及公理审计 340/340 通过。
随后源码审查发现 r+j 格维数的 signed 窄化问题；已改为容器下标，
没有收紧源输入范围。新增 6 项范围定理编译通过，最新 401 项原生
重组通过；212 项分离审计、全套原生和最新 3107 项 B2B 继续运行。
源码对应关系记录见[审查](../test-reports/vanhoeij-source-correspondence.md)，
具体版本和度量见[下标记录](../../proof/docs/devlog/2026-10-16-vanhoeij-index-width.md)。
