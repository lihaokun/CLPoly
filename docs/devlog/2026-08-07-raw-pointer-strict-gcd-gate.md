# Raw pointer 语义与 StrictGCD 准入收紧

## 日期

2026-08-07

## 做了什么

- 将 C++ 原始指针从引用类型中分离，保留分配区身份与以 `uint64_t` limb 为单位的偏移。
- 增加 `RawHeap`、`RawPtr`、`RawFault` 与 `RawExec`；越界、无效区和未初始化指针均显式失败。
- 为 `word3 *` 的重解释与三 limb 读写建立模型。
- 增加 `memcpy` 的逐 limb `copyU64` 语义和 `_poly_normalise` 的
  `normaliseU64` 语义；两者按剩余长度结构递归并传播 `RawFault`。
- 将尚使用 `Array.get!`、`Array.set!` 和默认填充值的 `divrem`、构造器及数组例程撤出 StrictGCD 准入集合。
- StrictGCD 生成器现在拒绝 `sorry`、`partial def`、未解析调用、默认数组访问、fuel、Safe/oracle/fallback 和 axiom。
- 删除 `size_t` 未初始化局部变量自动取 0 的翻译规则；在证明第一次读取前存在支配写之前，相关 HGCD 函数会生成失败而非产生伪语义。

## 为什么做

Lean 的 `Array.get!`/`set!` 在越界时具有默认行为，而 C++ 原始指针越界是未定义行为。若把前者直接作为 L1 执行语义，错误路径可能返回普通值并被误认为 L2 结果，不能构成真正的 C++ L1 到 Lean L2 精化。

## 自然语言证明草稿

严格生成链中的函数必须来自唯一的 Clang AST 函数体。每次原始指针访问在线性 heap 上执行：有效 region 且 limb 范围足够时，其读取或写入与 C++ 对应地址一致；否则返回 `RawFault`，不产生多项式结果。精化定理随后以布局、边界、非别名和初始化不变量排除 `RawFault`。源循环保持原控制流，并以剩余迭代区间为自然数度量证明每条递归边严格下降；不引入 fuel 或替代算法。

## 关键决策

- 指针偏移统一用 limb 计量，以保持 `uint64_t *` 与 `word3 *` 的地址重解释。
- 未初始化指针用 `region = none` 表示；任何读取立即产生 `invalidPointer`。
- 未初始化整数不赋予任意可执行值。必须先完成 def-use/支配写证明，函数才能进入严格链。
- 暂时缩小 StrictGCD 的已准入范围，明确撤销不满足内存语义标准的旧 `divrem`，而不是保留名义上的完成状态。

## 遇到的问题与解决方式

- `hgcd_mat` 的 C 数组固定为四项：模型保留数组布局，并用 `HgcdMat.Valid` 记录长度不变量。
- 项目根目录没有 Lake 配置：Lean 验证从 `proof/lean/` 执行。
- Pass 1 全量脚本有三个既存 fixture AST 选择失败；本次新增的指针与 `hgcd_mat` 类型测试均通过，Lean 模型与 StrictGCD 均编译通过。

## 涉及文件

- `proof/cpp2lean_v2/ir_types.py`
- `proof/cpp2lean_v2/passes/pass1_parse.py`
- `proof/cpp2lean_v2/passes/pass5_operator_resolve.py`
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/cpp2lean_v2/cast_table.py`
- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/constructor_map.py`
- `proof/cpp2lean_v2/tests/build_strict_gcd.py`
- `proof/cpp2lean_v2/tests/test_pass1_parse.py`
- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Generated/StrictGCD.lean`

## 度量

- 耗时：约 1 小时（语义设计、AST/IR 探测、生成与编译验证）
- 迭代：4 轮生成/编译/审计
- Lean 新增/修改行数：约 150 行；StrictGCD 撤回约 260 行未准入定义
- 对应 C++ 行数：本阶段为底层内存语义基础，未新增已完成精化的 C++ 算法行
- 放弃的方案：把未初始化 `size_t` 解释为 0；原因是它会把缺失的支配写证明变成可执行后备值
