# van Hoeij 格下标的容器宽度

日期：2026-10-16。先写自然语言证明
[§51](../../nl-proof/vanhoeij-eventual-separation.md)，再修改实现和证明。

## 问题与修复

规范格增加模数行后实际维数为 r+j，两个 signed int 的可表示性
不能推出和也可表示。原生 LLL 将矩阵/行长度转成 int，而严格
模型使用 Nat，存在真实的范围差异。将格维数、行下标、短行数组、
内积索引、iota 起始值改为 size_t；计算 r+j 前先扩大操作数。
候选因子编号仍为 int，遵循原来的 Hensel 因子数约束。

降序大小规约改为正剩余计数，对应严格模型已经采用的
`extraSizeReduceLoop`：remaining=k−1,...,1 时访问 remaining−1。
不使用 unsigned −1 哨兵，不改变实际访问顺序，也不增加新的
输入范围限制。纯重试和充分精度成功定理的数学内容无需改变。

## 证明及对应审查

新增 [VanHoeijIndexWidth](../../lean/CLPoly/Refinement/VanHoeijIndexWidth.lean)
的 6 项定理：r+j<2^32<2^64；UInt64 和 USize 的实际加法精确；
由实际初始化、实际建格输出及新增列数界，得到实际矩阵和每行
大小的界；循环自增和正剩余数递减对应自然数运算。

源调用路径、参数块、截断/GS 后缀、完整验证、机器边界、现存
callee 契约及旧兼容入口范围记录于
[源码对应审查](../../../docs/test-reports/vanhoeij-source-correspondence.md)。
此记录明确区分人工控制流翻译、内核证明与 B2B，不声称通用
AST 翻译器或任意编译器产物的普遍等价性。

## 验证

- 新增 6 项定理与扩大分离审计 212/212 已通过，仅标准公理，3657 jobs：
  `/tmp/clpoly-vh-index-width-audit2.log`。
- 最新 C++ 重组回归 401/401：`/tmp/clpoly-vh-index-width-native.log`。
- 最新全套原生回归 27 文件通过、0 失败：`/tmp/clpoly-vh-index-width-all-native.log`。
- 最新独立 C++ B2B 驱动 `_build/cpp_driver_index_width` 已重新编译。
  18 个向量文件共 3107 项正在与最新 Lean 驱动对照，日志
  `/tmp/clpoly-vh-final-b2b.log`。不把前一版两端的运行结果混记为
  最新版本通过。
- 严格边界、已证登记、B2B 运行时隔离及补丁检查通过。完整
  `lake build` 也已通过，3612 jobs，日志
  `/tmp/clpoly-vh-final-full-lake-build.log`。
- 三版本固定基线已准备：eb45ca48（safety net 之前）、428a27fe
  （当前 PR）、修复工作区。
  `test/prepare_vanhoeij_benchmark.py` 仅对隔离副本添加统一 RAII
  计时，区分 Hensel、LLL、内部和外部 Zassenhaus；
  `test/run_vanhoeij_benchmark.py` 记录相同编译参数、二进制哈希、
  每次确切因子/FLINT 校验及完整样本，明确记录旧版超时或失败。
  三个基准均已编译，命令与二进制 SHA-256 记录在隔离构建目录的
  `build-manifest.json`；必须等后台测试完成后再计时。当前没有合并。

## 度量

- 耗时：未记录精确墙钟起止，不估填。
- 迭代：2 次显式 Lake 构建；第一次需修正 UInt64.ofNat 与
  numeral 形式的引理匹配，第二次新增定理及审计均已通过。
- Lean 新增：68 行、6 项定理；现有严格执行模型的 Nat 行索引
  与循环顺序保持不变。
- C++：涉及建格索引、LLL、候选行读取、GS 后缀及独立 GS 重算；
  两端注册的行类型同步修改。算法决策与算术公式不变。
- 原生编译：第一次因命令中的 GMP include 路径拼写错误失败；
  修正路径后的编译和 401 项回归通过。不是实现失败。
- 放弃的方案：收紧输入度数范围以隐藏 r+j 溢出；仅扩大数组
  元素类型而继续用 int iota 计数；使用 unsigned 负哨兵。
