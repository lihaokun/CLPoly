# Pass 1 parse 完成 — Stage 2 Week 1 Day 1

日期：2026-04-21

## 做了什么

实现 cpp2lean v2 的 Pass 1 parse，并通过完整三路 1:1 等价性审核。

### 代码产出（commit bb4cf66 + 7696a4f）
- 目录骨架 `proof/cpp2lean_v2/`：复用 v1 的 `clang_hybrid.py` (192) + `class_map.py` (531)
- `ir_types.py` 474 行：HIR + MIR 数据类（40 节点）
- `passes/pass1_parse.py` 713 行：AST → HIR₀
- `tests/fixtures/*.cc` 5 个 fixture + `test_pass1_parse.py` 9 测试
- `tests/smoke_pass1_full.py` + `tests/dump_all_hir.py`
- `docs/.../survey/pass1-smoke.md` 烟测报告

### 审核方式
三路 Agent 并发逐函数 1:1 审核：
- Zp agent: 13 函数
- Univar agent: 34 函数（含 factorize 3 实例）
- Wang agent: 18 函数

## 为什么做

Stage 1 完成了完整设计（13000 行文档），Stage 2 进入实施。Pass 1 是管道入口，必须牢靠——否则后续 7 个 Pass 都建在沙子上。

## 关键决策及其理由

### 1. Pass 1 允许 UnresolvedOp / LambdaExpr 占位
- 所有 `CXXOperatorCallExpr` 保留为 `Call(UnresolvedOp(op_name))`，由 Pass 5 operator_resolve 解析
- Lambda 的 captures 留空，由 Pass 3 lambda_lift 推导
- **理由**：分层 Pass 设计要求每 Pass 只做一件事；在 Pass 1 解析运算符会重复 class_map 查表，违反单一职责

### 2. 类型映射保守：未识别的 STL 内部类型映射到占位 NamedType
- 如 `allocator_type` → `NamedType("Allocator")`、`_Bit_reference` → `NamedType("StlInternal")`
- **理由**：这些类型在 L1 精化证明中不重要，占位保留信息即可，Pass 5 再查表

### 3. 指针类型用 `RefType(is_pointer=True)` 而非新节点
- CLPoly 仅 1 处指针参数 (`__upoly_to_poly` 的 `comp_ptr`)，不值得独立 `PointerType` 节点
- **理由**：YAGNI；若未来指针使用扩展可再拆

### 4. assert 识别两种模式
- v1 的 `if (!cond) __assert_fail(...)` 形式
- 以及 Clang 宏展开的三元 `(cond ? (void)0 : __assert_fail(...))` 形式
- **理由**：实测 CLPoly 的 assert 大多是三元形式（23 个 site）——如果只识别 if 形式会漏掉全部

## 遇到的问题与解决方式

### 问题 1：初次烟测 55 函数有 Unknown 节点
- **根因**：`parse_type` 对 STL typedef、CLPoly 内部 typedef、lambda closure 类型等覆盖不全；parser 对 CXXForRangeStmt 的 inner 索引错位
- **解决**：~165 行增补修复，烟测 55 WARN → 0 WARN

### 问题 2：三路审核 11 函数 FAIL
- 修复 4 项：
  * assert 三元识别（7 函数）
  * 指针类型（1 函数）
  * factorize 多实例 dump（1 函数 → 3 实例）
  * Lambda params 填充（25+ 函数可读性）
- **结果**：重审 65/65 PASS

### 问题 3：RecursionError in assert_hir0_invariant
- **根因**：初版用 `__dict__` 递归走所有节点，走到 Enum 时循环
- **解决**：用 dataclass_fields 显式枚举 + 只递归 IR 节点（不追 TypeIR）

## 量化结果

| 指标 | 数值 |
|---|---|
| 单元测试 | 9/9 通过 |
| 烟测覆盖 | 65 函数，0 WARN / 0 FAIL |
| 人工审核 | 65/65 PASS（Zp 13 + Univar 34 + Wang 18）|
| UB 识别 | 21 函数 25 条 `require` |
| Lambda 参数填充 | ~9 处 lambda 正确 params 数 |
| Unknown AST 节点 | 0 |

## 涉及的文件

### 新建
- `proof/cpp2lean_v2/README.md`
- `proof/cpp2lean_v2/ir_types.py`（474 行）
- `proof/cpp2lean_v2/passes/pass1_parse.py`（713 行）
- `proof/cpp2lean_v2/passes/__init__.py`
- `proof/cpp2lean_v2/tests/__init__.py`
- `proof/cpp2lean_v2/tests/fixtures/{5个.cc,README.md}`
- `proof/cpp2lean_v2/tests/test_pass1_parse.py`（244 行）
- `proof/cpp2lean_v2/tests/smoke_pass1_full.py`
- `proof/cpp2lean_v2/tests/dump_all_hir.py`
- `docs/.../survey/pass1-smoke.md`

### 复用（从 v1 100%）
- `proof/cpp2lean_v2/clang_hybrid.py`
- `proof/cpp2lean_v2/class_map.py`

### 修改
- `.gitignore`（`__pycache__/`）

## 下一步

**Pass 2 ref_elim**（HIR₀ → HIR₁）：
- 消除 `T&` / `T&&` 参数 → tuple 返回值
- 依据 `hir-design.md` §4、`cpp-subset-semantics.md` §4.3-4.4 引理 L4.1
- 预估 ~200 行 Python + ~100 行单元测试
- Pass 2 完成后：HIR₁ 不变量 = 无 `is_ref` 参数

Pass 2 后续：Pass 3 lambda_lift / Pass 4 iter_recognize / Pass 5 operator_resolve / Pass 6 ssa_build（Cytron）/ Pass 7 loop_lower / Pass 8 codegen。

## 度量

- 耗时：~4-5 小时（写 + 修 + 审 + re-audit）
- 迭代：3 轮（初版 → 修 Unknown → 修 11 FAIL）
- Python 代码：~1500 行（含复用）
- 对应 C++ 行数：4852（65 函数）
- Agent 调用：4 次（1 次 fixture 选择 + 第一轮 3 agent 审 + 第二轮 3 agent 重审）
- 放弃的方案：
  - 在 Pass 1 解析运算符（留给 Pass 5）
  - 独立 PointerType 节点（用 RefType 的 is_pointer 标志足够）
  - `__dict__` 递归遍历 IR（改用 dataclass_fields + _IR_NODE_TYPES 白名单）
