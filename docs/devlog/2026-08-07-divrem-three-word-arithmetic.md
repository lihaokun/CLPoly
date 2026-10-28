# `divrem` 三字算术严格生成

## 日期

2026-08-07

## 做了什么

- 将 `_umul128`、`_add_carry3`、`_lll_mod_preinv` 加入严格 GCD 的 C++ AST 生成闭包。
- 为当前 x86_64 分支的 `GCCAsmStmt` 增加严格形状识别和三条 add/adc 指令的可执行 Lean 语义。
- 修正 `UInt64 → UInt128` 为显式零扩展，避免把 Lean 类型标注误当作机器整数转换。
- 生成后的 `StrictGCD.lean` 通过 Lean 构建、可复现检查与零占位审计。

## 为什么做

稠密长除法不是抽象域除法：C++ 使用 128 位乘法、三 limb 带进位累加和预逆元归约。若用 L2 多项式除法替换任一辅助函数，就不再是 C++ L1 执行的语义精化。

## 自然语言证明草稿

`_add_carry3` 的 x86_64 分支依次执行 `addq b0, lo`、`adcq b1, mid`、`adcq 0, hi`。第一步的 128 位和低位写回 `lo`、高位是第一进位；第二步把 `mid+b1+第一进位` 的低位写回 `mid`、高位作为第二进位；第三步将第二进位加入 `hi`，所有字段均按 64 位回绕。Lean 定义逐条实现这一状态变换。

## 关键决策

- 汇编识别只接受 `lo, mid, hi, b0, b1` 五个操作数的精确 AST 形状；其他汇编仍产生未知节点并被严格门拒绝。
- `UInt128` 转换使用 `BitVec.ofNat 128 x.toNat`，明确表达 C++ 无符号零扩展。
- 辅助函数仍保持 C++ 中间 limb 状态，不提前替换成模素数的数学结果。

## 问题与解决

Clang JSON 不携带汇编模板文本，但保留完整输出/输入操作数，并且源码函数与目标架构分支已确定。生成器对该唯一操作数形状映射到命名为 x86 的指令语义；不匹配时拒绝生成。原有 `(x : UInt128)` 不是有效的零扩展，因此新增显式转换函数并统一更新 cast 表。

## 涉及文件

- `proof/cpp2lean_v2/passes/pass1_parse.py`
- `proof/cpp2lean_v2/cast_table.py`
- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass8_codegen.py`
- `proof/cpp2lean_v2/tests/build_strict_gcd.py`
- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Generated/StrictGCD.lean`

## 度量

- 耗时：约 0.6 小时
- 迭代：4 轮探测、生成与 Lean 编译
- Lean 新增/修改行数：模型约 18 行，生成产物约 70 行
- 对应 C++ 行数：约 70 行（三个辅助函数）
- 放弃的方案：未选择预处理器强制走 portable `#else`，因为那不是当前 x86_64 实际执行分支；未选择数学模运算 oracle，因为会跳过 limb 级实现。
