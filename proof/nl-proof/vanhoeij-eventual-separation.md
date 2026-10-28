# van Hoeij 充分精度分离：结式路径

日期：2026-10-16。状态：数学基础正在形式化；尚未连接生产控制器。

## 与文献、当前代码的关系

[正式发表版本](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.655.pdf)
§5.1 的完整格含模数块，§5.2 的截断量由输入系数界确定，§5.3 用
有界非零结式不能被任意大模数整除来证明终止。不能把早期预印本的
矩阵文字抽取结果当作最终公式。

当前代码还使用随模数增长的 `r+16` 位截断上限，并每轮重建格。
因此文献终止定理不能直接实例化。下面先建立独立的结式基础，再
决定如何加入不丢失充分精度信息的尝试；不以数学定理的假设冒充
已验证的执行性质。

## 1. 模共同因子迫使结式被模数整除

设 F,G∈Z[X]，次数界 u≥deg F、v≥deg G，且 u+v>0。
设 M>1，q∈(Z/MZ)[X] 首一、次数为正，同时整除 F mod M 和 G mod M。

Mathlib 的 `exists_mul_add_mul_eq_C_resultant` 在任意交换环成立，
所以在 Z/MZ 上给出 A、B 及 F·A+G·B=C(Res(F,G;u,v))。
q 整除等式左端，故整除右端常数。首一正次多项式除常数的余式
就是该常数；整除又使余式为零，故常数为零。结式与环同态交换，
于是 M 整除整数结式。不要求 Z/MZ 为域；这是模素数幂情形所必需的。

## 2. 初等、显式的结式系数界

若 F、G 所有系数绝对值≤B（B≥0），Sylvester 矩阵每个元素也≤B。
其维数 d=u+v。行列式的 Leibniz 展开有 d! 项，每项绝对值≤B^d，
由三角不等式得 |Res(F,G;u,v)|≤d!·B^d。
这比 Hadamard 界松，但足以证明充分精度存在，且其位长可显式控制。

如果 M>d!·B^d，则第一节的整除关系强迫整数结式为零。
取实际次数 u=deg F、v=deg G，把结式映到 Q。
次数在 Z→Q 下不变；当 F 在 Q[X] 不可约时，零结式说明 F、G
不互素，再由不可约性得 F|G。这里的不可约前提只描述一个待分析的
真实整数因子，绝不假设算法的输出不可约。

转回整数整除时不能要求 H−e_i f' 本原：它可能为零或有非单位内容。
零多项式直接处理；非零时分成内容乘本原部分，内容映到 Q 后可逆，
先消去内容。对两个本原多项式用 Gauss 引理，再乘回原内容。

同一支持上两个局部残差的比较也在 Z/MZ 完成：q_j 同时整除
H−e_i f' 和 H−e_j f'，相减得到 q_j | C(e_j−e_i)f'。
`IsCoprime q_j f'` 消去 f'，首一正次 q_j 整除常数迫使该常数为零。
由 |e_i|、|e_j|≤K 和 M>2K 得整数差为零。局部互素条件目前仅是
这个代数引理的条件，尚须由实际 Hensel/好素数的证书推导。

边界：F 非零且正次由 Q 上不可约性给出；G=0 允许；M>1
排除平凡模环。若一般行列式界的 d=0，空乘积与 0!=1 仍满足公式。

## 3. 接到 CLD 分离的后续义务

对保留基行的整数因子坐标 e 和对应的小残差多项式 H，局部因子 q_i
整除 H−e_i f' mod M。若真实因子 g 的支持含 i，q_i 也整除 g mod M。
前两节在足够大的 M 下推出 g | H−e_i f'。对同一支持中的 j 再约化，
由局部互素性得到 M | (e_i−e_j)；若 M>2 max|e_i|，即有 e_i=e_j。
这将证明所有保留行在每个真因子支持上恒定。结合已证的真 indicator
保留及完整恢复验证，才能推出单次尝试成功。

仍需逐项证明：

1. 从实际 LLL 的 size-reduced 与 Lovász 性质和实际保留末行 GS 界，
   推出每个保留行的整数坐标界，界不能依赖不断增长的模数。
2. 从实际共同截断格恢复上述 H，证明模 CLD 恒等式及系数界。
   要求该证明所用截断除数最终独立于模数，或另证保留信息的增量策略。
3. 从实际 Hensel 支持和 CLD 证明两个局部整除条件，不能新增此类 oracle。
4. 证明实际分组、恢复验证在空间完全分离且恢复精度足够时返回完成；
   现有成功结果的正确性定理不能反向用作执行完备性。
5. 接到实际精度策略；有限次数后 Zassenhaus 回退仍不是这个终止结论。

2026-10-16 进度：§4–10 已补实际保留坐标界、完整残差多项式及
执行导出的 CLD 次数/局部整除/互素条件。跨精度的统一界仍需组合，
真实因子支持尚未代入最终结式比较，§3 的分组执行完备性和精度策略
仍未完成。后面的基础定理不得替代这几个剩余结论。

## 4. 实际 LLL 的 GS 范数上界

从 `lllReduceWithNorms` 的实际初始化及主循环运行取出同一最终状态。
初始化的 k=1，所以 size-reduced 与 Lovász 的已处理前缀都真空成立。
现有主循环不变量给出最终所有 |mu[i,i−1]|≤1/2，以及
`(3/4−mu²)·norm[i−1]≤norm[i]`。真实 GS 范数为正，故
`norm[i−1]≤2·norm[i]`。这里引用的是实际执行定理，不是假设一个
抽象 LLL 黑盒已返回规约基。

对实际后缀扫描递归归纳：若保留数量 t>0，则最后保留的
`norm[t−1]≤B`；递归删除分支沿归纳传递，首次不删除分支正好给出
该不等式。把上述相邻范数界迭代 t−1−i 次，即得所有 i<t 满足
`norm[i]≤2^(t−1−i)·B`。空保留前缀不索引 t−1。

后续还要用实际 Gram 等式将 GS 范数界转换为每个物理基行坐标的界。
若 N 是完整格维数，粗略可用行范数平方≤N·2^N·B，整数坐标
绝对值再用 `N·2^N·B+1` 控制。随后才可把这些界代入结式阈值。

物理坐标界的展开：实际 Gram 恒等式 `G=L D Lᵀ` 的对角项给出
`sum_c row[c]^2=sum_j L[i,j]^2 norm[j]`。当 j>i 时 L[i,j]=0；
当 j=i 时 L[i,j]=1；当 j<i 时 size reduction 给出 |L[i,j]|≤1/2。
因此对保留的 i<t，每个右侧非零项≤`2^N B`，共 N 项，得到上述
行平方范数界。任一物理坐标的平方≤整行平方和。若 x²≤C 且 C≥0，
由非负平方 `(abs x−1)²` 可得 `2 abs x≤C+1`，从而 `abs x≤C+1`。
最后把有理数不等式映回整数；不能直接把 GS 范数当作普通行范数。

形式化路径已检查：Mathlib Resultant.Basic 的通用结式/Bézout/映射/零判据，
Polynomial.Div 的首一余式性质，PrincipalIdealDomain 的不可约互素判据，
以及项目现有实际 LLL size-reduced/Lovász 证明。先完成 §1–2，不宣称
尚未完成的 §3 已被这些基础引理替代。

## 5. 从截断坐标恢复任意保留行的有界残差

这里的行组合系数可以为负，不能使用仅适用于真因子 0/1 indicator
的已有截断误差引理。设 D>0、m'=floor(m/D)≥1，原 CLD 系数
满足 |C_i|≤m；e_i 是整数行组合系数，|e_i|≤E，项数至多 r。
设一个截断数据坐标 y=Σ e_i floor(C_i/D)+z m' 满足 |y|≤Y。

首先，欧几里得恒等式 C_i=D floor(C_i/D)+(C_i%D)、m=Dm'+(m%D)，
及 0≤余数<D 给出 |floor(C_i/D)|≤m'+1。负 C_i 必须保留这个 +1，
不能把负数的 floor 当作向零截断。三角不等式和乘积绝对值给出
|Σ e_i floor(C_i/D)|≤rE(m'+1)，故
|z|m'≤Y+rE(m'+1)≤(Y+2rE)m'。因 m'≥1，得 |z|≤Y+2rE。

取 H=Σ e_i C_i+z m。逐项展开欧几里得恒等式可得
H=D y+Σ e_i(C_i%D)+z(m%D)。余数非负且<D，所以
|H|≤D Y+rED+|z|D≤D(2Y+3rE)。同时 m 整除 H−Σ e_i C_i。
这给出同一个 H 的模 CLD 关系和整数高度界，而不是分别假设它们。
当 e_i 的物理因子坐标为 s e_i 且整数 s≥1 时，|e_i|≤|s e_i|；
因此已有保留坐标界 K 可同时取作 E、Y。

形式化拆成 floor 绝对值界、带权有限和界、模数行系数界、残差
恒等式和残差高度界。使用 Int.mul_ediv_add_emod、Int.emod_nonneg、
Int.emod_lt_of_pos、Finset.abs_sum_le_sum_abs，最后由正数乘法消去
m'。空支持允许；E、Y 非负明确列为前提；m'≥1 排除退化模数行。
这个界仍含 D，只有进一步证明 D 最终由输入固定，才可将它作为
充分精度论证中的统一高度界。当前按 m 增长的截断策略尚不满足该条件。

执行桥接：取实际 LLL 返回的 transform 第 k 行为整数组合系数。
已有 `LLLTransformRel` 说明这同一行乘原始规范格正好等于返回的物理
matrix 第 k 行。前 r 列规范公式为 s e_i，用 s≥1 和物理坐标界
控制 |e_i|。数据列规范公式由实际 truncateCld 的系数定理改写为
上面的 y 等式，因而直接应用残差界。有限和中非因子坐标置零，
项数用完整维数 N 粗界，得到 D(2+3N)K；不额外假定 transform
本身的小系数，也不假定待构造残差已满足高度界。此阶段只覆盖
实际加入的每个 CLD 列，未加入的次数和完整残差多项式留待下一步。

## 6. 实际参数选择的统一短向量 cutoff

当前 divisorLoop 从 1 翻倍至 32D≥height，后续位数限制只会增大 D。
所以实际 selectParameters 返回时始终有 height≤32D。对 D>0，
由整数除余恒等式，若 height≤D C，则
floor((height+D−1)/D)≤C：反设商至少 C+1，乘 D 后与分子上界
D C+D−1 矛盾。取 C=32，得到实际 residual≤32+2r。
在 height≥0、r≥0 的合法调用下 residual 也非负。

若实际加入列数 J≤n，则 cutoff B=r·32²+J·residual²
≤r·32²+n·(32+2r)²。因此 cutoff 不依赖 m，虽然残差高度界中的 D
尚会依赖 m。本阶段证明参数部分以及给定 J≤n 下的不等式；实际
CLD 宽度、列循环和输入次数之间的 J≤n 仍要从执行中导出，不能
把这一条件悄悄当成已经完成的顶层结果。

## 7. 完整次数覆盖与残差多项式

实际列遍历的次数是偶数列 k/2、奇数列 W−1−(k−1)/2。
对任意 d<W，取逆索引
`k=2d`（2d<W 时），否则取 `k=2(W−1−d)+1`。
分别展开奇偶分支可得 k<W 且遍历次数恰为 d；W=0 时没有此类 d。
这些都是自然数线性算术与除以 2 的取余关系，可以直接由 omega 检查。

从 current=0 开始且 target≤W 时，实际追加循环的成功执行恰好加入
target 列。对 target−added 归纳：added<target 时也有 added<W，
故必须执行追加并递归；若追加报错则不可能满足成功执行前提。
停止支给出 added=target。取 target=W 得完整覆盖，不能仅用现有
“加入列数≤target”反推已经覆盖所有次数。

宽度扫描保持 maximum 的上界，同时最终 maximum≥每个尚待扫描的
非空多项式的首项次数+1。对数组尾长归纳，逐项用 max 的两侧界。
规范多项式首项次数等于 natDegree，因此所有次数≥最终宽度的系数
均为零；空多项式由其 toPoly=0 单独处理。给定每个原 CLD 的次数
小于 n，截断不会引入新次数，故实际宽度 W≤n。这个原 CLD 次数界
仍须从真实模除法和导数执行得到。

完整残差取次数 d<n 的有限单项式和。d<W 时使用上述逆索引对应的
实际模数行系数 z_d；W≤d<n 时取 z_d=0。后一段的截断系数为零，
由实际 floor 桥接知原系数在 [0,D) 内，因此对应残差≤N K D，仍受
D(2+3N)K 控制。d≥n 时原 CLD 和构造的残差系数都为零。逐系数可
得同一个残差多项式的模同余与全系数高度界；之后才接局部支持和结式。

形式化中先采用等价且不需要预先传入 n 的表达式
`H=P+C(m) Z`，其中 P 是原 CLD 的整数组合，
`Z=Σ(d<W) monomial(d,z_d)`。这样模同余直接由 m 映到 Z/mZ 为零
得到，原 CLD 的高次项也完整保留。全系数高度界不需要原次数上界；
只有后续给 H 的次数作输入统一界时才使用原 CLD 次数<n。对 d<W
用实际列残差界，对 d≥W 用 floor 为零的余数界。W=0 时 Z=0，
所有次数均走后一支。该构造与前述 n 截止的有限和在原次数界成立时相同。

## 8. 局部 CLD 整除关系

在任意交换系数环上，首一 q_i 整除 f 时写 f=q_i g_i。首一除法
返回 g_i，由乘积求导得 `q_i | f'−g_i q_i'`。若 j≠i 且 q_i、q_j
互素，q_j 首一并整除 f，则由 f=q_j g_j 和互素消去得到 q_i|g_j，
因此 q_i 整除第 j 个 CLD。有限加权和中，i 项减去 e_i f' 后可整除，
其余每项也可整除，得 `q_i | Σ e_j CLD_j−e_i f'`。这在 Z/p^aZ
上有效，不使用“系数环是域”的错误前提。

与实际残差多项式连接时使用 §7 已证的模同余和实际 CLD semantic 字段。
此处首一、整除、互素条件先作为局部代数引理的前提。现有
Algorithm/Hensel.isCoprime_lift_pk 可把模 p 的互素性提升到模 p^a；
仍须组合实际 Hensel 输出与好素数性质，不能把这些前提算作已连接。

## 9. 从实际 CLD 语义导出次数界

设 n=deg f>0，模环中 q 首一。若 deg q=0，q'=0，CLD 为零。
若 deg q>deg(f mod m)，首一除法的商为零。其余情形下，商次数
是 deg(f mod m)−deg q，导数次数≤deg q−1，因此 CLD 次数<n。
这里只用一般交换环上的首一除法次数定理，不依赖模环为域。

实际 CLD 的模语义给出次数≥n 的整数系数 c 映到 Z/mZ 后为零，
故 m|c。对称约化已有 `2 natAbs(c)≤m`，且 m>1，所以 |c|<m。
若 c≠0，整除关系反而给出 m≤|c|，矛盾。故原整数 CLD 在次数≥n
也为零；由 Polynomial.natDegree_le_iff_coeff_eq_zero 得次数<n。
再接已证的实际截断宽度与残差次数定理，消去此前独立传入的 CLD
次数前提。CLD 可能为零，必须使用 n>0 处理 natDegree(0)=0。

## 10. 实际 Hensel 输出的局部条件

实际 Hensel 入口给出同一个返回模数 m=p^a、a>0，以及规范化后
所有物理因子的乘积等于 C(u)f mod m，其中 u 可逆。每个实际因子
整除此乘积，消去单位 C(u) 后整除 f。首一性由实际物理首项为 1
的证书取得，不需要另行指定首一代表。

模 p 的各物理因子互素性已有执行证明。将 i<j 和 j<i 两支合并后，
对同一对整数多项式调用已有 isCoprime_lift_pk，得到它们在实际
返回模数 m 上的互素性。

好素数证书说明 f mod p 无平方；有限域是完美域，所以它可分，
即与自己的导数互素。物理因子 q_i 整除 f mod p，因此 q_i 与
f' mod p 互素。再次对 q_i 和整数 f' 调用相同提升定理，再用
导数与 map 交换，得到实际模数下 q_i 与 f' 的互素性。所有条件
来自同一选素数和 Hensel 执行，没有假设算法最终输出不可约。

组合执行时把这些已导出性质放入只含 Prop 的分离输入证书。证书
描述同一 Hensel 因子数组和同一 CLD 返回数组：尺寸相同、模数>1、
每个 CLD 的真实语义/系数界/次数界，以及局部因子的首一、整除、
两两互素和与 f' 互素。构造该证书的定理以实际 Hensel 执行和实际
cldPolys 成功返回为前提，不能让调用者凭空给出另一套因子。
随后由证书与真实 truncateCld 运行统一导出残差次数和局部整除。
局部因子的正次数也由同一执行导出：其模 p 约化由既有选素数/Hensel
细化证明为不可约，故次数为正；整数物理因子首一，因此映到模 p
及任意 m>1 的次数均保持。这提供结式共同因子所需的正次数条件。

## 11. 将维数与 cutoff 统一为输入界

设 r 是实际提升因子数，n 是源多项式次数。实际截断数组的大小仍为 r，
宽度≤n。目标列数≤宽度时，已有循环列数界给出 J≤n；实际规范格
维数为 N=r+J≤r+n。选定参数满足 residual≤32+2r，因而
`B≤B0=r·32²+n·(32+2r)²`。B0 非负。

对非负 cutoff，函数 `N·2^N·B+1` 在 N、B 上单调：分别用自然数
嵌入、2 的幂单调性、非负乘法和加法，得实际物理坐标统一界
`K0=(r+n)·2^(r+n)·B0+1`。原 CLD 残差的统一高度界随之为
`S0=D·(2+3(r+n))·K0`。这里 D 仍明确出现，不将它误称为独立于
精度的输入常数；只有采用最终固定 D 的策略后，S0 才是统一常数。

## 12. 统一结式阈值与支持上的相等坐标

令 L=Σ|f_k|，n=deg f。任意源系数≤L；求导后的第 k 个系数是
(k+1)f_(k+1)，非零时 k+1≤n，故所有导数系数≤nL。真实因子
F|f 的系数由现有 Mignotte 定理控制在 A=choose(n,n/2)L。
对保留行的残差 H 和 |e_i|≤K0，G=H−e_i f' 的系数≤S0+K0 nL。
取统一 B=1+A+S0+K0 nL，则 B≥1，F 与 G 的所有系数均≤B。

真实因子次数≤n；H 次数<n，f' 次数≤n−1，因此 G 次数≤n。
结式尺寸 d=deg F+deg G≤2n，阶乘和 B 的幂均单调，故
`d! B^d≤(2n)! B^(2n)`。若 m 大于右端，且局部因子 q_i 同时整除
F mod m 和 G mod m，现有结式定理推出整数 F|G。对同一支持内的
q_j，再约化得到 q_j|H−e_i f'；实际 CLD 局部关系还给 q_j|H−e_j f'。
两式相减、由导数互素消去 f'，得 m|(e_i−e_j)。若同时 m>2K0，
则 e_i=e_j。这只是充分精度下的分离，仍要连接实际分组成功和
能达到此阈值的精度策略。B 中的 D 必须最终固定才能证明阈值可达到。

实际支持连接不另作存在性假设：使用已有
`hensel_divisor_physical_support`，对同一次 Hensel 运行和真实分解
f=F·complement 取出合法物理索引数组及对应模乘积等式。成员 q_i
整除该乘积，因而整除 F mod m。将同一数组的任意两个成员和同一次
LLL 的任意保留行代入上面的结式比较，得到实际 transform 的两列
逐行相等。得到的是单次完整格在显式阈值下的分离性质，不是仅证明
某个抽象格或任选支持存在。

## 13. 从逐行相等接到实际比较与分组循环

实际 shortSpanFromNorms 是从 0 开始的连续索引数组，长度就是
retainedPrefixCount；iotaLoop_spec 给出长度和每项值。因此保留行上的
相等直接满足 candidateColumnsEqual_true_complete 的逐项条件，实际
比较器返回 true，索引安全由同一次 LLL 的 transform 尺寸给出。

记 E(i,j) 为实际比较器对列 i,j 返回 true。已有证明给出 E 的自反、
对称、传递性。内部 assignCandidateClass 从 start 扫描到 r，其对
任意有效索引 j 的最终值恰为：若 start≤j<r、初始 classes[j]=none
且 E(rep,j)，则 some classId；否则保留原值。按剩余扫描长度归纳，
访问 j 的分支直接写入，其他分支不改变 j；已有非空标号不会被覆盖。

外层保持两个不变量：(a) 小于当前 column 的项均已赋值；(b) E(i,j)
蕴含 classes[i]=classes[j]（含两者都为 none 的情况）。初始全 none
满足二者。已赋值的当前列直接跳过。否则先播种当前列，再运行内部
扫描。对所有 j<r，更新结果可写为：若 classes[j]=none 且 E(column,j)
则 some classCount，否则原值。j=column 由自反性；j<column 因 (a)
不会是 none；j>column 用内部精确公式。E(i,j) 使两者原值相同，且
E(column,i)↔E(column,j)，故 (b) 保持；(a) 由旧值保留和播种保持。
终止时同等列标号相同，再由已证明的全部赋值、实际 collect 的成员
覆盖，得到两者出现在同一输出组。此层先证明实际成功运行的分组
完备性；无故障返回及后续验证成功仍须另接，不能暗含为前提。

分组自身的无故障返回可独立证明。CandidateColumnsValid 保证比较器
每次读到的行与两列均有效；按剩余行数归纳，它返回某个 Bool。内部
扫描的 class 数组长度≥r，且代表列<r，因此比较器和每次写入都有效；
按剩余列数归纳，它返回数组并保持长度。外层同样按剩余列数归纳，
并保持所有已赋标号<classCount；播种的新标号仅等于旧 classCount，
故下一步的新界是 classCount+1。已有 ids_bound 引理保证内部扫描
保持这一界。终止后全部列已赋值，collect 的每个标号都小于输出
数组长度，所以按剩余来源列数归纳得到实际输出。空保留行分支直接
返回空数组；此分支无故障，但不宣称验证完整分解成功。

组合时先从真实 LLL 执行取得 CandidateColumnsValid，再用无故障定理
产生真实 groups 数组及运行等式。给定保留行非空，同一 Hensel 支持
上的两个有效物理索引由 §12 逐行相等，再由比较器完备性返回 true，
最后由分组闭包得到同一输出组。这个组合不假设 extractCandidates
成功；保留行非空目前仍是显式前提，下一层需由非零真实因子短向量
保留性消去。当前结论对同一支持内任意两成员成立，尚不替代后续
分组等于完整不可约支持及恢复验证必成功的证明。

## 14. 消去保留行非空假设

先对同一次实际 LLL 证明：若初始格有非零整数坐标 c，且其物理
向量平方范数≤cutoff，则保留前缀非空。用返回的幺模 transform 的
整数逆把 c 转为约化基坐标 z。已有重构定理保持物理向量和范数，
已有非零定理给出 z≠0。取 z 最后一个非零位置 t，已有 GS 下界
给出 norm[t]≤cutoff。若保留前缀为空，实际 suffix 扫描的性质却给
每个 norm[t]>cutoff，矛盾。这使用真实 GS 状态和真实返回索引数组。

随后证明真实正次数因子的 Hensel 支持非空。把支持等式投影到原素数
p；好素数条件及 F|f 保证 F 的首项系数模 p 非零，所以映射保持次数。
若支持为空，其模像与 1 相伴，成为单位，次数为零，矛盾。
这里投影、首项非零和模数为 p 的正幂都从已有执行证书导出。

为证明保留行非空，不必预先找到一个不可约因子：对整个 f 与补因子 1
调用已有 hensel_factor_truncated_short。f 正次数，故其实际支持有成员 i；
短向量的第 i 个因子坐标为 1，因此非零。原 CLD 高度取实际 cldBound
的数学值，且 m>height，实际 cutoff 至少为该短向量定理的界。
将此向量代入第一段即可消去非空前提。cutoff 的下界与此前结式证明
使用的上界必须同时保留：只有上界不能排除 cutoff=0 的空前缀。

## 15. 输出组恰为一个完整真实支持

首先加强实际分组的结论。对任意固定物理列 i，取外层扫描返回的
classes[i]=some c。全部赋值定理保证这样的 c 存在。实际 collect 的
成员来源定理和覆盖定理给出：输出组 c 包含物理列 j，当且仅当
classes[j]=some c。这里源索引经 UInt32/Int32 转换，必须使用 r≤2^31
及已证明的往返性质，排除不同 Nat 索引转换成同一机器整数。
外层已有同组比较正确性和上一节闭包性质，分别给出
classes[i]=classes[j] ⇒ E(i,j) 与 E(i,j) ⇒ classes[i]=classes[j]。
因此组 c 的成员恰是与 i 列比较相等的所有有效物理列。

对真实因子 F，用同一次 hensel_factor_truncated_short 同时取得支持
indices 与初始短向量 coefficients；不分别选取两个支持然后假设它们
相同。支持非空，取成员 i。若比较器 E(i,j) 为真，实际保留行逐项相等，
已有短向量保留定理推出 coefficients[i]=coefficients[j]。前者为 1，
后者为支持的 0/1 指示值，故 j 在支持内。反向若 j 在支持内，直接把
该同一支持等式代入结式分离定理，得到全部实际保留行相等，再由
比较器完备性得到 E(i,j)。结合上段，实际输出组 c 恰为该完整支持。
这个结论仍不暗含其他组都对应因子或恢复验证已成功；这些是后续
partition 完备性所需单独组合的性质。

## 16. 实际收集顺序、重复次数与支持乘积

仅有成员集合相等不能推出多项式乘积相等。定义某 class 的剩余来源
列表为 `[start, …, classes.size−1]` 中满足 `classes[i]=some classId`
的索引，等价于 range/drop/filter。对实际 collect 运行按剩余长度归纳，
证明输出该行恰为初始该行附加上剩余来源列表的 Int32 转换。当前来源
属于该 class 时实际 push 一项，否则该行不变；与 filter 的头部展开
逐项一致。完成时剩余列表为空。

从空数组开始，该行因而是有界、无重复的来源列表的映射。对每个成员
用 Nat→UInt32→Int32→Int64→Nat 的有界往返，解码后的组恰为上述来源
列表，所以解码组 Nodup。结合上一节完整支持的双向成员关系和合法
Hensel 支持的 Nodup，用 List.perm_ext_iff_of_nodup 得到解码组与实际
支持索引的排列等价。乘积按该排列不变，故真实输出组的
SelectedProductMod 恰好等于同一个 Hensel 支持乘积。该证明不依赖
验证已经成功，不借助“已验证因子模 p 无平方”来反推无重复性。

## 17. 实际组乘积执行与机器次数累加

partitionProductLoop 的每个来源索引由 CandidateIndicesValid 保证非负
且在 lifted 内。非零模数使 multiplyNormalizeModRaw 必定返回；按剩余
组长度归纳，实际组乘积循环正常返回 product 和 Int64 degree。
乘积语义复用已有该循环到 trialProductLoop 的执行投影，不引入新乘法
oracle。次数值精确等于从输入 degree 开始，对剩余组条目执行
`acc + get_deg(lifted[index])` 的 foldl；按循环递归逐项对应可得。

规范且非空的物理多项式满足 get_deg=Int64.ofNat(natDegree)，这条
等式本身允许机器截断，不能直接据此做有序比较。用 Int64.ofNat_add
按列表归纳，实际 foldl 等于对自然数次数总和作一次 ofNat。随后再
由该组真实支持乘积的次数等式，证明总和等于真实因子次数；只有在
真实因子次数<2^63 后才用 toInt_ofNat_of_lt 排除符号溢出并通过正次数
检查。完整证明不能把机器累加直接当成无限精度整数加法。

组的次数等式可在原素数 p 上证明：把已证明的模 m 支持等式投影到
p（p|m）。真实因子的首项模 p 非零，保持次数；物理 Hensel 因子首一，
映射到 p 后同样保持次数。首一多项式的列表乘积次数等于次数之和，
前面的非零常数首项系数不改变次数。这样实际组次数和等于真实因子
次数，无需假设整数提升因子本身整除 f，也不在有零因子的环上误用
域的次数乘法定理。

## 18. 非首一因子的实际恢复与本原化

设当前剩余 R=F·Q，F 本原且首项为正，Q 首项也为正，实际组的模支持
满足 F mod m = lc(F)·selected。验证代码乘积初值为 lc(R)，因此实际
模乘积并非直接恢复 F，而是恢复 T=lc(Q)·F。由首项乘法等式和支持
等式，实际 product 与 T 同余。若每个 T 系数绝对值严格小于 m/2，
已有对称恢复定理给出实际 symmetric 的整数多项式恰为 T。

对于规范物理输入 T，contentLoop 等于数学 content。F 本原，且
lc(Q)>0，因此 content(T)=lc(Q)；T 首项为 lc(Q)lc(F)>0，primitiveRaw
不会翻转 content 的符号。它返回的实际 content 正是 lc(Q)。将此
代入已有 primitiveRaw_toPoly 并消去非零常数多项式，得到实际恢复
因子恰为 F。这同时处理非首一 F，不把本原化假设成恒等操作。

随后 F|R 使实际精确除法返回空余数。返回商的语义由除法恒等式和
F≠0 消去得到 Q。规范表示保持性提供下一次验证迭代的物理条件；
§17 的次数等式与机器界使当前因子为空、次数不匹配或非正次数三个
拒绝条件均为假。这里先完成单组的真实执行，整个列表的因子对应和
最终剩余 1 仍需逐组归纳，不能在单组定理中预设。

恢复所需的系数精度不另加一个不可达 oracle：若 F 和 Q 都整除原输入
f，既有 Mignotte 界及统一系数界 B 给出 |lc(Q)|≤B、|F_k|≤B，所以
`2|(lc(Q)F)_k|≤2B²`。f 正次数时 2≤2n，且 B≥1，既有阶乘/幂单调
定理给出 `2B²=2!B²≤(2n)!B^(2n)`。因此同一结式阈值已经保证所有
这些缩放恢复目标都严格落在半模数区间内，不需要把输出恢复成功
放入前提。

## 19. 给定精确有序支持分解时，完整验证循环成功

为组合单组步骤，记录一份数学目标数组 targets：它与实际 groups 等长，
乘积是原输入，每项本原、正首项、正次数，且每项的模支持正是对应
实际组。这里不假设输出不可约，也不假设任何验证运行已成功；从
实际分离/分组证明构造全部 targets 仍是后续义务，不能把此条件记录
当作最终端到端证明。

循环不变量为实际 remaining 的整数多项式等于 targets 从 index 开始
的后缀乘积，并保持规范表示。正首项的列表乘积仍有正首项。当前
目标与下一后缀均整除原输入，因此同一阈值控制缩放恢复系数；当前
剩余次数也不超过原输入次数。§18 给出真实恢复步骤并返回该下一后缀。
按 groups.size−index 归纳，得到正常完成整个循环。

终止时后缀为空，remaining 的整数多项式为 1。不能只用数学等式跳过
实际机器检查：规范表示和头项次数定理给出头次数为 0；严格降序的
尾部不可能有次数小于 0 的项，所以数组恰有一项，且头系数为 1。
于是末尾的 size/get_deg/系数检查全部通过，返回实际排序后的 factors。
原输入正次数又排除空 groups，使公开 validatePartition 的空组分支
不可达，故整个入口返回 complete=true。

## 20. 实际输出类非空及类的唯一性

为排除“真因子都有对应组，但还存在额外空组”的缺口，在实际外层
partitionCandidateColumns 循环维护：每个小于 classCount 的类编号，
都有某个物理列的 classes 槽保存该编号。初始 classCount=0，性质为空。
遇到已分配列时编号不变；遇到未分配列时，把该列写成 classCount，
再运行内层 assign。旧类的见证列不可能是这个原先未分配的列，故
set 保持其槽；assignCandidateClass_preserves_some 保持所有已有槽。
新类以刚写入的当前列为见证，同样被内层保持。由 count-column 的
强归纳，终态每个编号仍有见证。此证明只使用真实成功执行及数组界，
不要求列比较器的数学性质。

collectCandidateClasses_size 给出输出数组长度等于最终 classCount；
collectCandidateClasses_includes_source 将每个槽的见证列放入对应实际
组，所以每个实际输出组都有物理来源列。空 rows 分支直接返回空数组，
不存在合法组编号；因而非空性结论也覆盖此退化分支。

对唯一性，若同一合法物理列属于两个组，机器索引界 count≤2^31
保证 Int32 编码不混淆不同来源。已有 collector 的成员充要条件把两次
成员关系分别化为同一 classes[source] 槽等于 some c 和 some d；
Option.some 的单射性给出 c=d。这里只证明实际类的结构，后续仍需
用所有模因子的真因子覆盖来构造完整目标数组。

## 21. 同一精确模支持确定唯一正首项本原因子

设 F₁、F₂ 均整除非零正次数输入 f，本原且首项为正，并且
F₁ mod m=lc(F₁)S、F₂ mod m=lc(F₂)S。交叉乘首项系数后，
lc(F₂)F₁ 与 lc(F₁)F₂ 同余。§18 的统一阈值分别以 F₁、F₂
作为因子和缩放来源，保证两边所有系数严格小于半模数；对称区间
唯一性给出整数多项式等式。两边取 content，本原性和正首项使
content 分别等于 lc(F₂) 和 lc(F₁)，所以两个首项相等；再消去
这个非零常数多项式，得到 F₁=F₂。此结论不要求不可约性，避免
在构造目标数组时以待证的输出不可约性或恢复成功来证明唯一性。

## 22. 因子族的支持覆盖实际输出组

设数学因子列表的乘积为输入 f，每个因子已有 §16 的精确实际组及
Hensel 支持。任取实际输出组，§20 提供组内来源列 i。它在原素数 p
下的模因子 qᵢ 是素元且整除 f mod p，因此整除因子列表映射后的
乘积。对列表归纳，素元不能整除空乘积 1；非空列表用素性把整除
分解为头项或尾部，从而找出 qᵢ 整除的某个真实因子 F mod p。

将 F 的模 m 支持投影到 p，得到 F mod p 与支持内模因子乘积相伴。
各个物理模因子是两两不相伴的素元，prime_atom_dvd_support_iff
因此给出 i 属于 F 的支持。F 的精确组包含 i；§20 的实际类唯一性
将这个组与原先任取的输出组识别。所以每个实际组都对应因子列表
中的某一项。这个桥接步骤以数学因子族和逐因子的精确支持结论为
输入，尚不替代从原输入构造该族及证明完整排列/乘积的后续工作。

## 23. 按实际组顺序构造完整验证目标

把无重复的数学因子族表示为数组 factors。逐因子的精确组结论给出
函数 assigned：每个因子下标映到一个实际组下标，并有该因子的合法
Hensel 支持及精确成员充要关系。这是逐因子定理的选择，不是预设
每个输出组已有目标。

assigned 满射：对任意输出组，用 §20 取来源列；重复 §22 的素元
整除列表乘积论证，找到某个 factors[j]，其支持包含该列，故其
assigned(j) 与任取组相等。这里保留找到的真实数组下标，避免只得到
因子值存在而丢失排列所需的下标信息。

assigned 单射：若两个因子下标映到同组，两个因子的完整模支持相同。
每个因子整除原输入，正首项且本原；§21 在同一显式阈值下给出因子
相等。原 factors 的 Nodup 将元素相等提升为下标相等。由双射构造
有限类型等价 e，并定义 targets[i]=factors[e⁻¹(i)]。

targets 与实际 groups 等长。有限乘积沿等价重排保持乘积，所以
targets 乘积仍是 source；本原、正首项、正次数逐项从 factors 继承。
e(e⁻¹(i))=i 把逐因子支持等式改写成第 i 个实际组的支持，构造完整
ExactPartitionTargets。空数组同样适用有限等价与空乘积；实际正次数
输入会由前面的 source 乘积约束排除空因子族。

这些有限重排仅用于证明实际验证循环的数学目标，不替换实际循环、
不生成另一套算法，也不要求或假设输出不可约。尚须由原输入得到
无重复正首项本原因子族，并实例化实际 Hensel/LLL 的逐因子分组数据。

## 24. 从输入构造用于分析实际算法的数学因子族

输入 f 本原、首项为正，且好素数 p 保持其首项非零、f mod p 无平方
因子。先证明 f 在整数多项式环中无平方因子：若 g²|f，则 g|f，
本原性传递给 g，好素数条件保持 g 的次数。映射到 p 后，平方整除
与 f mod p 无平方因子使 g mod p 为单位，因此 g 的次数为零。本原
常数多项式必为单位，所以得到 Squarefree f。

取 normalizedFactors f 的列表为数学因子族。非零 f 的规范因子乘积
是 normalize f；整数多项式规范化只乘首项的 normUnit，而正首项的
normUnit=1，故乘积恰是 f。无平方因子等价于规范因子多重集 Nodup。
每项是不可约因子并整除 f，所以本原；Gauss 引理给出映射到 Q 后
不可约，其次数为正，并由到 Q 的单射映射保持次数。每项已规范化，
首项也已规范化，是非负整数；不可约项非零使首项非零，故首项为正。

这个存在性因子族只提供 §23 重排证明和真实短向量/结式证明的分析
对象；实际输出仍须由实际 Hensel、LLL、组提取、恢复和除法链得到。
不直接从 UFD 存在性断言算法返回，也不把其输出不可约性写入前提。

## 25. 充分精度下实际分组与验证的组合成功

在同一实际 Hensel 输出、CLD 截断数组、全宽规范格和实际 LLL 执行
之上组合前述定理。由规范格维数和 CLD 数组大小得到模因子列数不
超过矩阵大小，实际 LLL 的列形状定理与 extractCandidates_succeeds
给出一份实际 groups。§24 从原输入本原、正首项、好素数条件构造
数学因子族；每项整除 f，所以取其补因子并应用 hensel_exact_support_group。
该定理返回的组数组与刚取得的 groups 都是同一次确定性函数的返回值，
由 Except.ok 单射识别。这里没有假设某个独立选取的支持恰相同。

实际 Hensel 定理给出原素数下模因子的素性与不相伴性；实际乘积相伴
定理给出每个模因子整除原输入的映射。§23 因而构造完整目标数组。
输入正次数与目标乘积排除空 groups，也排除产生空 groups 的空保留行
分支。实际索引有效性定理、Hensel 物理头/规范表示/首一性提供验证
循环的表示条件；同一阈值与输入机器次数界给出 validatePartition
入口成功。最终结论同时返回实际 extract 的成功结果和实际验证成功。

此定理仍以实际 Hensel、CLD 和 LLL 的前序执行事实、规范全宽格及
精度不等式为输入；下一阶段必须把这些事实接入实际列循环和精度
控制器，证明可达，不能将单次充分精度成功冒充完整终止性。

## 26. 输入固定的截断上界与参数成功阈值

当前 D=max(D₀,2^(bits(m)−r−16)) 无上界，不能从 m 增长推出充分精度。
修正为 D=min(max(D₀,2^(bits(m)−r−16)),D₀²)，未进入位数截断分支
时仍取 D₀。D₀ 是原来从 1 倍增到 32D₀≥H 的正整数，只依赖输入。
D₀≥1 给出 D₀≤D₀²；max 和 min 的序关系保证 D₀≤D≤D₀²，所以
旧有 32D≥H、正除数、残差界 T≤32+2r 均保持。

取固定 Dmax=D₀²，R=32+2r，r≥0。若 m>H 且 m>Dmax(8rR+2)，
则 m>D(8rR+2)。整数除法给出 floor(m/D)≥8rR+2，因此它既大于 1，
也严格大于 8rT；参数选择的所有拒绝分支均不可达，返回实际参数。

separationResidualHeight 对非负 D 单调，其余系数界项与 D 无关；
所以 separationPolynomialHeight 及其非负整数幂、结式阈值也对 D
单调。用 Dmax 计算固定结式阈值后，m 超过它就超过实际 D 对应阈值。
再同时要求 m>H、m>Dmax(8rR+2)、m>2K，即得到只依赖输入的统一
充分精度门槛。该阈值不作为生产运算的大整数求值，而用于证明；
后续仍须连接实际提升指数增长、机器溢出条件、全宽列循环与重试。

位数规则仍由实际 C++ sizeinbase 和严格模型的相同运算决定；本证明
只用 min/max 夹界，不假设或重写机器位数运算为数学对数。

对 p>1，整数幂无界定理给出 N，使 p^N>max(统一门槛,1)；后者
同时排除 N=0。整数幂对指数单调，所以每个 a≥N 都超过同一门槛，
实际参数选择成功且实际 D 的结式与坐标不等式成立。这只证明所需
素数幂精度存在，尚未将机器指数范围或当前有限重试当作无界增长。

## 27. 全宽 columnRound 的正常返回不能拒绝

对实际 columnRound 的正常返回使用现有 trace 定理，取得它真实构造
的初始格、追加列结果、LLL 输出、组数组及验证返回。实际目标 width
对应截断 CLD 的完整扫描宽度；buildCldMatrix_full_width 给出 added
恰等于此宽度。实际 target 已是 Int32，added≤target.toInt.toNat 与
Int32 的有符号上界保证 added 的 Int32 转换精确，无需另加次数 2^31 界。
实际 count 等于 Hensel 因子数，截断保持数组长度。

由真实初始/追加格执行得到 CanonicalCldBasis，由 CLD 数据得到宽度
不超过输入次数、矩阵维数不超过 r+n。参数选择的实际 residual 等式
把 columnRound 使用的 cutoff 识别为真因子短向量界，原 cutoff 上界
定理提供结式分离所需的统一界。统一输入门槛给出实际 D 的两个
分离不等式；所有前提因而实例化 §25 的充分精度验证定理。

该定理产生的组和验证返回与 trace 的实际调用参数相同。函数确定性
先识别组数组，再比较验证返回的 Bool，得到 accepted=true。这里仍
只分析正常返回；不能由此忽略前序 RawExec 的错误分支或声称整轮
无条件返回，后续另接执行完备性和入口条件。

## 28. 正目标列循环的正常退出必须完成

先由 columnRound_trace 与 validatePartition_success_frame 证明每个
accepted=true 的实际轮次都令返回结果 complete=true。随后对实际
columnLoop 的剩余宽度 width−target 作强归纳。目标为正且不超过
width；调用出错与已有正常返回条件矛盾。若接受，使用上述 frame；
若拒绝且 target=width，与 §27 的全宽不可拒绝结论矛盾；若拒绝但
未到全宽，实际 nextCldTarget_progress 保证下一目标仍正且严格增长，
ControllerExecution 的已证明下降条件使归纳成立。

通用循环引理以“全宽正常返回不可拒绝”为条件，随后必须用实际
充分精度定理实例化，不能把该条件留作最终算法的成功 oracle。
结论分析实际正常返回，尚不证明所有内部错误不可达；target=0 的
入口也没有被正目标前提覆盖，须另外处理。

## 29. 零扫描宽度时，一个真实不可约因子最多含一个模因子

先证明零宽度意味着每个原始 CLD 的每个系数绝对值小于 D。实际
truncateCld_spec 给出正除数、长度保持和规范化，cldSpiralWidth_sparseCoeff_zero
把所有截断系数置零，truncateCld_sparseCoeff 识别它们为原系数的
向下取整商；floor_zero_abs_lt 给出所需严格界，包括负系数情况。

对任意两个模因子索引 i、j，选择系数向量 e=δᵢ。宽度零意味着
rowCarryPolynomial=0，而有限和的唯一非零项使 rowResidualPolynomial
恰等于第 i 个原始 CLD。K₀ 的定义有末项 +1，故每个 e 坐标均满足
|e|≤K₀。D>0、K₀≥1、2+3(r+n)≥1 给出 D≤S₀，因此上述 CLD 系数界
也满足已有结式分离定理的残差高度条件。

若同一真实本原不可约因子的完整模支持含 i、j，由
CldSeparationData.scalars_equal_of_sufficient_precision 得 eᵢ=eⱼ。
若 i≠j，则左端为 1、右端为 0，矛盾，所以 i=j。该结论直接使用
实际 CLD 截断及源因子，没有引入零宽度分支实际未执行的 LLL。
这是零宽度入口的支持结构结论，仍须接完整单元素分区恢复成功。

进一步将此结论接到同一次 Hensel 提升：hensel_divisor_physical_support
提供合法递增索引数组与完整模支持；素数幂及输入好素数保证真实因子
首项模 p 非零，divisor_modular_support_nonempty 排除空支持。若支持
至少含两个索引，取数组第 0、1 项，前述零宽度定理迫使它们相等，
与 LegalCombination 的严格递增距离界冲突。所以数组长度恰为 1，
其模乘积等式化为 F mod m = lc(F)·qᵢ。此步骤从真实 Hensel 执行
构造索引，不把单元素支持作为调用者假设。

## 30. 用真实单元素分组构造完整恢复目标

给定由输入导出的无重复真实因子族及 §29 的逐因子单元素支持，为
每个真实因子选择唯一模因子索引。singletonGroups_size/entries
说明实际入口生成的第 i 组恰为 [i]，Int32 无溢出事实说明它的解码
索引仍为 i。因此逐因子支持与实际第 i 组的 SelectedProductMod 一致。

映射的单射性直接使用正首项本原因子的模支持唯一性及因子族 Nodup。
满射性逐实际组证明：该组对应的模 p 素因子整除输入乘积，所以整除
某个真实因子的模 p 像。将该真实因子的单元素完整模支持投影到 p，
其首项为非零标量，得到它与所选模 p 因子 Associated。于是前一素
因子整除所选素因子；不可约性给出 Associated，实际 Hensel 原子的
互异性迫使两个索引相同。这样每个实际单元素组都被命中。

最后调用已证的双射重排定理构造 ExactPartitionTargets。这里完整
目标数组是证明产物，不是算法成功或输出不可约性假设；所用分组
是入口实际 singletonGroups 的返回，不引入虚拟 extractCandidates。

## 31. 零列宽时实际单元素分组验证成功

由输入本原、正首项和好素数条件，按 §25 构造无重复真实不可约因子
族。对每个因子应用 §29 的真实 Hensel 单元素支持定理，再把实际
Hensel 原子的素性、互异性、输入整除关系送入 §30，产生与入口实际
singletonGroups 返回对齐的全部恢复目标。

同一次 Hensel 输出提供规范多项式、非空数组和首一性；实际单元素
组的已证索引合法性及 count 长度等式给出恢复循环物理前提。于是
validatePartition_complete_of_targets 产生实际验证返回 (true,after)。
验证函数确定性意味着任何已知该次返回都接受且 complete=true，
可用于排除充分精度入口“单元素验证失败后又遇零列宽”的路径。
CLD 数据目前仍由显式执行给定；后续入口组合必须从入口 trace 构造。

## 32. 真实重组入口在充分精度下正常返回必完成

先连接实际 Int32 宽度扫描：实际 CLD 数据的次数界与截断不增次数
给出自然扫描宽度≤输入次数。输入次数<2^31 时，每个非空截断多项式
的首项次数加一≤自然宽度<2^31，所以 cldWidth_eq_natural 可用于实际
扫描结果，排除有符号溢出。

直接展开 tryRecombine 的正常返回执行。断言失败和任何 RawExec
错误分支都与正常返回条件矛盾。cldBound_polynomial 将实际高度
识别为输入固定的数学高度；选参充分精度定理排除 none 返回。
若最初单元素分组验证接受，使用实际成功 frame 得到 complete=true。
否则保留实际拒绝返回，并从后续真实 cldPolys、truncateCld 执行
构造 CldSeparationData 和宽度等式。

若自然宽度为零，§31 给出同一单元素分组、同一初态的验证返回 true，
与入口记录的 false 矛盾。若自然宽度为正，实际 min(2,width)>0，
调用 §28 的已实例化列循环完成定理。所有中间执行事实来自入口，
调用者只提供实际 Hensel 执行、源输入合法性、机器大小界、统一
充分精度不等式及入口正常返回，不提供 CLD、截断、LLL 或成功证书。
该结论仍是正常返回完成性，内部无错执行与外层充分指数可达性另证。

将入口得到的 complete=true 送入已有
hensel_tryRecombine_complete_irreducible，可进一步得到实际输出因子
乘积等于输入且每个输出因子在 Z[x] 中不可约；complete 不再由调用者
假设。这个组合仍明确保留入口正常返回前提，不能声称无条件终止。

## 33. 实际初始矩阵与 LLL 初始化必正常执行

先证明 setInitialDiagonalLoop 在 size×size 方阵上不会越界。按
size−index 归纳；index<size 时两次读取由方阵形状保证，内外 set
保持矩阵大小和每行大小，下一 index 使度量下降。结束时返回当前
矩阵。zeroMatrix_size/get 与 zeroMatrixRow_size 为 makeInitialMatrix
提供初始形状，因此该函数对任意 scale 均正常返回，不需要先假设
执行成功或 scale≠0（后者只涉及格的非退化性）。

对非空方阵 initializeLLL，应用上述实际单位矩阵构造，再用
dotRows_eq_fin_sum 取得第 0 行内积的正常返回。zeroQQMatrix 的大小
和行大小、Array.replicate/set 的长度构成 GramStorageShape。
initializeGramSchmidtLoop_succeeds 从下标 1 完成后续实际初始化。
只证明无越界和返回存在；矩阵非退化时的正范数与执行不变量继续
使用原有 initializeLLL_concrete_valid，不添加数值成功假设。

## 34. 实际行消元正常返回

subtractRowsLoop 在四个输入行长度都至少为 size 时，每个 index<size
的读取都合法；按 size−index 归纳直接构造返回，累积数组只 push，
不需要额外输出数组长度假设。subtractMatrixRows 在矩阵与变换矩阵
均为同阶方阵且 target/source 均有效时，四次行访问及上述行循环
全部正常返回。

对 sizeReduceAt，实际执行不变量给出 mu 行/列大小，j<k<n 排除
两次入口越界。roundQQ_eq_round 已给出无条件正常返回；若 q=0，
源程序直接返回原状态。否则应用上述矩阵行消元执行定理，以及
reduceMuPrefixLoop_eq_array。后者的总函数结果保持 mu 大小和每行
长度，因此最终写入 mu[k][j] 仍合法。这里只构造源操作返回，数学
不变量保持继续使用原有已证的 preserves_execution_valid 定理。

## 35. LLL 单步的两条实际分支都能正常执行

对额外降序消元循环按 remaining 归纳，前提 remaining≤k<n；每次
调用 §34 的 sizeReduceAt，已有不变量、k/norms 不变及变换关系保持
给出递归前提。变换关系包含两矩阵同阶与变换矩阵每行大小，不需要
假设返回状态的物理尺寸。

lllStep 首先用 k>0、k<n 正常完成对 k−1 的消元。消元后的实际
不变量保证 norms 和 mu 的四次访问合法。Lovász 成立分支应用上述
额外消元总执行；不成立分支的 matrix/transform/mu 两行交换均由
合法 k、k−1 保证正常。QQ 行交换和后续一次 set 保持行列尺寸，故
updateMuAfterSwapLoop_eq_array 的全部索引前提成立，完成实际更新。
该证明同时覆盖零或非零 newNorm 分支；范数函数本身不返回 RawExec
错误，其数学正确性由既有不变量保持证明承担。

两分支还直接给出下次 k>0：推进为 k+1，交换为 max(k−1,1)。
因此以 concreteLLLRank 作强归纳可以构造 lllMainLoop 的正常返回。
k≥n 时直接返回；k<n 时上述单步正常返回，已有 step_valid、变换
关系保持和严格 rank 下降给出递归调用全部前提。初始 k=1，故没有
遗漏 k=0 的断言错误。这里的循环终止来自既有行列式/索引度量，
不是人为燃料或有限重试。

## 36. 实际 LLL 含最终范数接口的完整执行

insertShortRow 的输入行和已收集索引均合法时，所有内积均为同一行
与自身，dotRows_eq_fin_sum 排除错误；按 sorted.size−position 归纳
完成插入。collectShortRows 按 matrix.size−row 归纳，已证插入保持
索引合法性为下一轮提供前提，故从空数组收集必正常返回。

对非空、方阵且行列式非零的实际输入，§33 构造 initializeLLL 的
返回，既有定理得到 ConcreteLLLExecutionValid 及真实初始变换关系。
§35 构造主循环返回，再应用上述收集执行定理，直接展开实际
lllReduceWithNorms 得到其返回。最终范数取同一实际 reduced.norms，
没有把结果、LLL 成功或缓存范数换成外部假设。

## 37. 实际规范格构造正常返回

fillCldColumnLoop 在 cld.size≤matrix.size 且每行含有待写列时不会越界。
按 cld.size−index 归纳，每次内外 Array.set 保持行数与各行长度，
故递归前提成立。appendZeroColumn 已有无条件执行等式；给输入方阵
每行追加零以后，新列 cld.size+existingColumns 恰为旧行长，严格
小于扩展行长。上述填写循环正常返回；零行长度为新列+1，末尾的
模数写入也合法，构成 appendCldColumn 正常返回。

buildCldMatrixLoop 的每轮以 CanonicalCldBasis 为不变量。实际追加
函数正常返回后，已有 appendCldColumn_canonical 得到同一实际输出
的下一规范格，按 target−added 归纳。达到 target 或 width 均直接
返回，所有目标/宽度边界原样保留。最终 buildCldMatrix 因而从规范
输入格正常返回，不以目标是否足够分离因子作为构格的执行前提。

## 38. 模除法实际余式的系数区间与空余式

实际 divmodMergeLoop 写入的每个系数都经过正模数的 floor remainder，
故位于 [0,m)。先证明 pushNonzero 保持该区间，再按两个剩余列表
长度之和归纳覆盖 merge 的全部分支。DivmodTrace 的 done 分支
保留当前余式，vanished 分支仅删除首项，subtract 分支应用上述
merge 结论，因此实际除法返回的余式系数仍全部约化。

对规范且非空的约化数组，首项系数非零且在 [0,m)，其模 m 像非零。
由实际首项与数学 leadingCoeff/natDegree 的对应，模映射不降低
首项次数。实际模除法若除数物理首项为一，源输入的模像被除数整除，
则返回等式 R+QG=F 推出 G|R。若 R 非空，上述次数保持与实际停止
时每项次数严格低于 G 首项给出 deg R<deg G。这里 Z/mZ 可以有零
因子，不能使用一般整环的整除次数界：物理首项为一还给出模像 G
monic。写 R=G·S，R≠0 推出 S≠0，Monic.natDegree_mul' 在任意
交换环上给出 deg R=deg G+deg S，才得矛盾。故实际余式数组为空；
不把源断言检查当作额外成功前提。

## 39. 实际 CLD 生成与截断的正常执行

先取 concreteDivmodTermination，要求源与每个活动因子数组规范、次数
可表示，每个活动因子非空且首项系数为一，并且它的模像整除源模像。
按 active.size−index 归纳 CLD 外循环。§38 给出实际模除法返回以及
空余式，所以源 isEmpty 断言成立。导数函数直接返回其纯循环结果；
非零模数保证 modCoeffLoop、multiplyNormalizeModRaw 正常返回，正
模数保证 symmetricModRaw 正常返回。将这些实际返回代入源控制流，
递归剩余活动因子即可。没有以预先生成的 CLD 数组作为执行前提。

对于实际 Hensel 输出，从 outputOneHead 推导非空与首项为一；
outputCanonical 给出规范数组。既有 Hensel 首项次数界，加上规范
数组所有次数不超过首项，推出每项的机器次数界。源次数界由同一
规范首项等于输入 natDegree 推出。选素数和实际 Hensel 乘积恒等式
给出每个提升因子的模整除性，正的素数幂给出 modulus>1。因此上段
全部条件从输入及 Hensel 执行事实导出，适用于同一输出数组。

truncateTermsLoop 在 divisor≠0 时只执行保持数组大小的 set；以
size−index 归纳正常返回。truncateCldLoop 使用这一结论，然后
normalization 是总函数，外层 set 仍保持数组大小，同样按剩余项数
归纳。入口 divisor>0 排除唯一算术错误分支，得到实际 truncateCld
正常返回。空数组、越过末尾的起始 index 均由源直接返回分支覆盖。

## 40. 任意合法候选分组的验证正常返回

入口还会验证未能分离真因子的候选分组，因此不能只证明正确目标
分组的接受执行。对正 modulus、规范源数组及全部合法分组索引，
按 groups.size−index 归纳实际 validatePartitionLoop 的正常返回。
空 remaining、次数非正、候选次数不符及除法余式非空均是正常拒绝
分支，不是错误。余下分支由 partitionProductLoop_succeeds 得到
真实模乘结果；其规范性来自实际乘积循环的规范性保持，初始单项
系数为当前 remaining 的非零首项。实际对称恢复与 primitiveRaw
各自正常返回并保持规范性。检查通过保证 factor 非空，故已有
exactDivmodRaw_complete 适用，不需要因子事先整除 remaining。
若余式为空，实际 quotient 仍规范，递归前提由原来的索引合法性
以及次数无关的商规范性给出。循环结束仅返回 none 或排序后的
some 数组，两者均正常。入口空 groups 直接拒绝，其余将这一实际
循环结果映射为 Bool/VanHoeijResult。此定理仅保证执行，不宣称
任意分组都会被接受，也不把拒绝结果当成不可约性证据。

## 41. 任意列数回合及列数增长循环的正常执行

取规范源、非空活动因子数组、count 与该数组/CLD 数量相等且索引
可表示，原模数为正、格的 scale 和模数非零。§33 的单位格执行
与 §37 的规范格构造执行给出实际 built，已有构格规范性与尺寸
保持保证 built 非空，ControllerExecution.buildValid 给出真实
ConcreteLLLInputValid。§36 给出同一 ready 证明下的实际 LLL 返回。
真实变换矩阵与 retained 行索引满足 CandidateColumnsValid，故
实际 extractCandidates 正常返回。retained 为空时返回空 groups，
索引合法性为空真；否则已有成员来源定理与 count 的 int32 界给出
每个实际分组索引合法。§40 于是保证验证正常返回，构成整个
columnRound 正常返回，不需要完整列数或足够精度。

columnLoop 以实际 (width−target).toNat 强归纳。target≤0 直接返回；
正常回合的接受分支直接返回，拒绝且 target=width 也直接返回；
其余分支由 ControllerExecution.targetDecreases 严格减小度量，
targetBound 提供下一轮前提。归纳构造实际返回，而不假设循环已
成功执行。该结论仅保证正常返回；是否 complete 仍由充分精度
定理判定，不能混为一谈。

## 42. 从实际 Hensel 输入构造入口返回与充分精度分解

对具体模除法终止实例，选素数条件、至少两个模因子、实际 Hensel
输出正确性及源规范性/正次数/机器次数界蕴涵 tryRecombine 正常
返回。选素数与 Hensel 因子逐项相伴定理给出输出数量等于模因子
数量，素数正幂给出 modulus>1，故最外层断言通过。规范首项等于
natDegree，Int64→Int32 的转换在源次数界内精确，cldBound 正常
返回。singletonGroups 正常返回且索引合法，§40 构造第一次验证。
接受时直接返回；参数选择 none 时正常返回未完成结果。参数选择
some 时得到 divisor>0，§39 构造实际 CLD 和截断数组，既有大小
保持证明供给 columnLoop 的 countSize，§41 给出实际列数循环返回。
所有 dependent match 的证明参数只用于这些相同执行值，不引入
可替代计算结果的回调。

把构造出的实际入口执行应用于 §32 的充分精度定理，得到存在
after，实际入口返回 after，after.complete=true、因子乘积为源且
各因子不可约。这里不再假设入口已正常返回，也不假设 CLD/LLL/
分组验证的任何中间执行。Hensel 调用本身的正确执行仍为输入事实；
外层如何增长精度并再次获得这一事实，保留为独立控制器义务。

## 43. 外层精度增长：取消次数预算，显式机器上限

修复外层只在完整证书出现时返回分解；否则继续提高 Hensel 请求
指数，不再在第三次失败后转 Zassenhaus。保留 full_precision 作为
下一指数下界。新增源 helper：precision=INT_MAX 时抛出溢出异常；
precision>INT_MAX/2 时下一指数取 INT_MAX；其余取
max(full_precision,max(1,2*precision))。这样可到达全部可表示的正
请求指数，避免在半上限处提前报错，且不执行溢出的整数乘法。

实际 Hensel 成功蕴涵请求指数为零或正数，复用其已有源入口检查
定理作为 erased 事实。对 precision≥0 的任何正常 nextPrecision
返回，指数严格增加：饱和分支由 precision≠INT_MAX 与 signed 上界
得严格增加；倍增分支在 precision≤INT_MAX/2 下精确，precision>0
时 2*precision>precision，precision=0 时 max(1,...)>0。full 的
max 不会降低结果。因此实际 retryLoop 的
(INT_MAX−precision.toInt).toNat 严格下降，构造无燃料的完整生成模型。
错误分支原样返回，成功分解只能来自实际 complete=true 的入口。

此机器度量排除控制器无限运行，但单独不证明成功分解。还必须用
§42 的充分精度阈值排除正常输入在最大指数处仍未完成的情况。若
阈值超出机器可表示指数，则准确保留溢出结果，不冒充成功终止。
profile 仅用 first 布尔量区分初次/后续调用，避免无界 attempt 计数。

## 44. 重试成功结果的纯证书来源

按 §43 的实际精度剩余量归纳 retryLoop 的正常分解返回。Hensel 或
tryRecombine 错误不能产生此返回。入口 complete=true 时直接提供
同一次 Hensel、同一次入口调用及其实际因子数组；complete=false
时 nextPrecision 的错误也不能返回分解，正常值严格提高指数并可
递归应用归纳假设。不存在 Zassenhaus 来源分支。将此唯一来源接入
既有 Hensel 语义证书与入口不可约性定理，保持公共 FactorZZCorrect。
该正确性证明仍区别于“存在成功返回”的终止定理，不混用两者。

## 45. 实际 Hensel 模数达到请求指数

当请求 a>0，实际 HenselLiftTargetCorrect 不可能走 a=0 分支；
显式分支的目标精确为 p^a−1，其非负性由 p≥2 得到。已有实际
Hensel 循环终止时 outputModulus>target，故 outputModulus≥p^a。
这将输入固定的分离阈值与机器请求指数直接关联，不假设实际提升
模数恰等于 p^a（实际二次提升可能超过它）。

## 46. 为任意合法请求构造同一次实际 Hensel 返回

从选素数正确性、模因子物理规范性、源规范性/非空/次数界、模因子
数量至少二且可表示，以及该请求的 HenselLiftRuntimeReadiness，
构造已有 HenselLiftEntryInvariant。显式正请求给出 exponent 条件；
源首项与 leadingCoeff 对应给出首项模 p 非零，已证明的选素数因子
调整规范性、非空、首一和互素性逐项提供其余代数条件。readiness
字段只描述真实 Hensel 递归/归一化的执行前提，沿用公共原有契约。
调用已有 __hensel_lift_upoly_raw_ir_refines，直接得到实际运行与语义
证书，不再先假设 Hensel 已成功返回。所有 source/target/provider
参数保持同一实例，不添加能指定输出因子的替代函数。

## 47. 在可表示充分精度范围内，实际纯重试必成功

令 H 为输入的实际 CLD 高度公式，r 为选择的模因子数，T 为已有
vanHoeijInputThreshold(f,r,H)。明确机器可行条件 T<p^INT_MAX。
对正初始请求，以实际 (INT_MAX−precision).toNat 强归纳。§46
构造当前 Hensel 调用，实际大小保持给出合法的入口 countFits；
§42 构造当前 tryRecombine 返回。若 complete=true，实际控制器
直接返回其因子。若未完成且 precision=INT_MAX，§45 给出实际
modulus≥p^INT_MAX>T；Hensel 输出大小等于 r，将同一阈值代入
充分精度完成定理，得到 complete=true 的矛盾。因此未完成时指数
不是最大值，实际 nextPrecision 正常返回且严格增加，新的指数
仍为正并且 rank 严格下降，归纳构造后续实际返回。

最后与 §44 的实际成功来源/不可约性定理组合，得到实际纯重试
返回的分解正确性。这里没有 Zassenhaus 分支、固定尝试次数、
输出不可约性或重组成功 oracle。T<p^INT_MAX 是有界机器接口的
明确可行范围，不声称任意大小的数学输入都可用 int 请求指数处理。

## 48. 实际启发式的有符号范围与浮点转换

生产启发式改为两个可单独核对的实际 helper。浮点 ceil 值先检查
NaN 或 ≥INT_MAX，返回 INT_MAX；≤1 返回 1；剩下范围才转换为 int，
并保留 max(1,converted)。这样 NaN、两种无穷和所有越界浮点值均
不进入 C++ 未定义的转换，真实结果始终为正。Lean 对 positivity
只使用最后的整数 max 和两个常量分支，不为 Float 算术添加公理。
Float intrinsic 与 C++ 对应另以边界及有限值 B2B 核对。

项数直接以 size_t 转 double 后取 log，删除无必要的 signed int
中间值 N=size−1，避免 size=2^31 时的减一溢出。严格模型的 r 和
min_b 显式使用 Int32.toInt 转 Float，不经公共模型的 unsigned
隐式转换链；项数使用 Nat.toFloat。此变更不改变可表示项数下的公式。

完整指数 helper 保留 C++ while(pa≤target) 的先后顺序：条件成立
后先检查 exponent==INT_MAX，溢出则立即抛异常；否则 pa*=p，再
++exponent。严格模型使用同一 ZZ pa、Int32 exponent 和分支顺序，
不再用纯 Nat 循环结束后的 UInt32 检查替代它。正 pa、p≥2 保证
(target+1−pa).toNat 严格下降；计数自增的精确性来自 exponent<INT_MAX，
不是假设它永不溢出。零或负 target 在初始 pa=1 时直接返回 0。

## 49. 真实有界指数循环的成功与正性

循环保持 pa=p^exponent、exponent≥0。若 target<p^INT_MAX，则在
exponent==INT_MAX 且 pa≤target 的错误分支矛盾；其余继续分支
先证自增精确且非负，再保持 pa=p^next，按实际 pa 距离归纳构造
正常返回。另独立证明任意正常返回的 exponent≥初始指数；若
target≥1，从初始 pa=1、exponent=0 必执行一次安全自增，故最终
指数至少 1。源 guard 返回零指数的退化输入仍原样建模，不作无效
普遍正性声明。

## 50. 构造真实启发式前缀并连接纯重试

规范非空且物理次数可表示的非零源，由实际 Mignotte helper 的
已证执行得到同一 bound；取数学除子 1 的系数界可知 bound≥1。
源首项非零，故实际恢复 target=2*abs(lc)*bound≥2。若此实际 target
小于 p^INT_MAX，§49 构造 full helper 的正返回，§48 的浮点 helper
恒为正，最终 min 为正。所有函数为实际源函数，不给出额外的
启发式输出或成功运行 oracle。再连接 §47 的实际纯重试，消去
公共 lllFactorize 包装中的前缀正常返回与初始指数正性假设；保留
明确的恢复精度和分离精度可表示条件，以及已有 Hensel 就绪契约。

## 51. 格维数与容器下标的机器范围

新增模数行后，格维数为 r+j；即使 r 和 CLD 宽度分别可用 int
表示，r+j 也可能超过 INT_MAX。原 C++ 的 LLL 行数、短行编号及
内积循环将容器大小窄化成 int，与严格模型的 Nat 下标不一致。
修复使用 size_t 保存格行数、行号及列下标；因子编号和 CLD 螺旋
位置仍是 int，但在计算 r+j 和分配 r+j+1 前先转 size_t。

额外大小规约原 j=k−2,...,0 改为 remaining=k−1,...,1，实际访问
j=remaining−1，顺序逐项相同。k≥1、remaining>0 分别保证减法不
下溢。短行数组及其 iota 的起始类型均为 size_t，避免输出类型虽
扩大而 iota 内部计数仍以 int 溢出。候选分组只将格行下标扩大，
列/因子编号保持原来有符号接口。

由 r<2^31、CLD 宽度<2^31 及 j≤width，得到 r+j<2^32<2^64。
故本任务的 64 位目标平台上，容器下标的所有实际加法均可表示；
Nat 模型不再静默绕过 signed 行维数溢出。对 k<n 的递增，k+1≤n；
对剩余计数的递减，remaining>0 排除下溢。此修复保留全部实际
LLL 运算、顺序、候选分组和终止判据，不增加输出成功假设，也
不缩小已经声明的源次数/因子数输入范围。
