# HGCD 递归基础分支统一不变量

日期：2026-08-09

## 做了什么

- 保留 `_hgcd_recursive` 矩阵基础分支真实 `_mat_one` 初始化已经建立的恒等变换、带符号行列式和矩阵长度界。
- 将 `computeM = true` 的真实矩阵初始化及 A/B 拷贝封装为 `HgcdRecursiveRawInvariant`。
- 将 `computeM = false` 的真实 A/B 拷贝封装为同一不变量；矩阵语义与长度义务仅在矩阵实际被请求时成立。
- 增加按真实 `computeM` 分支统一入口的基础情形定理，供完整递归主体直接调用。
- 将第一递归调用后的真实 early-return 拷贝出口封装为统一递归不变量；矩阵描述符长度通过实际四项矩阵拷贝保持。
- 证明带符号行列式控制的两套 C++ 低半重构公式互为实际 HGCD 矩阵的逆，并在加回移位高半后恢复完整多项式 transform。
- 将第一子调用的统一 raw invariant、真实输入低/高切分与四调用重构执行组合，导出完整输入的 transform、GCD 保持及重构长度界。
- 为最终尾部加入纯物理 workspace，显式 frame 重构后的 A/B 以及 R、S、商缓冲区，并闭合真实重构到可选商更新/完整矩阵乘法的 raw 执行桥。
- 证明真实商矩阵的行列式取反、真实矩阵乘法的行列式乘法，以及两次递归矩阵围绕中间 divrem 的 transform/带符号行列式复合。
- 从真实 normalization 字段证明 raw 描述符零长度当且仅当零多项式，并证明非零描述符长度精确等于 `natDegree + 1`，建立最终矩阵次数界回到 C++ 长度界的桥。
- 从四个真实、规范化 raw 输入描述符推出单个矩阵乘法输出 `P*Q + R*S` 的精确规范化长度上界，为最终四项矩阵长度不变量提供可直接实例化的算术桥。
- 将该界实例化到真实 `_mat_mul` 返回堆；一个对 `Fin 4` 的统一定理同时约束四个实际 C++ 输出描述符，并保留左右输入矩阵的真实 raw 表示来源。
- 为源码商矩阵更新 `top + quotient*bottom` 建立规范化 raw 长度界；零多项式分支和乘积次数分支都从描述符表示推出，不引入规格计算或备用实现。
- 证明规范化 raw 描述符长度由其 L2 多项式唯一决定，并将商矩阵块的四项结果闭合为两条更新上界和两条精确保留等式，覆盖真实 row swap 与两次 guarded update。
- 将商矩阵四项界与真实 `hgcdRecursiveCombineMatrix` 执行组合，导出实际中间矩阵及最终 `_mat_mul` 四项描述符界；中间对象来自生成执行返回值而非规格侧构造。
- 保留真实 guarded multiplication 的精确 `lenLeft + lenRight - 1` 容量界，并让成对重构同时导出 A/B 两侧 normalization 返回长度上界；由源码恒等式 `k + lenC0 = reconstructed.lenB` 闭合最终两项均不超过外层输入长度。
- 从规范化 low/high/output 三个 raw 描述符证明：当 low 严格位于 shift 以下且 high 非空时，真实 `liftHigh` normalization 返回长度精确为 `shift + highLength`；成对重构现已保留可实例化的 A/B 精确长度条件。
- 直接按真实 `hgcdIterLoop` 良基执行证明停止时 `inputLength/2 < lenA`：继续分支的下一 A 就是通过 source guard 的旧 B；该事实加入递归长度归纳结果，并在 iterator/base/early-copy 路径中真实传递。
- 把两个统一基础分支定理加入严格 HGCD 源码门禁。

## 为什么做

良基递归总精化需要所有返回分支提供同一种语义归纳结果。此前基础分支只导出了部分 raw 表示，丢弃了矩阵初始化已经证明的语义和长度证据，无法直接作为递归主体的基础情形。

## 关键决策

- 矩阵启用分支的所有语义均来自真实 `hgcdRecursiveBase` 执行以及已有的 `hgcdIterInit_refines`，没有另行计算规格结果。
- 矩阵禁用分支不声称输入矩阵具有任何代数语义，只保留源程序确实执行的两次拷贝及 GCD 不变性。
- 长度归纳信息来自真实 `_mat_one` 描述符长度和拷贝后的返回长度。

## 验证

- `lake build CLPoly.Impl.StrictHGCDRawRefinement`
- 严格 HGCD 源码门禁（在提交前运行）
- 新增定理公理审计（在提交前运行）

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-hgcd-recursive-base-invariant.md`
