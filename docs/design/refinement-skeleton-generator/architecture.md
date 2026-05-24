# 精化定理骨架生成器 — 架构设计 v2

> 对应 workflow.md §2.2 架构阶段
> 日期：2026-05-24

## 1. 设计变更

**v1 方案**：独立的 Python 生成器 + 手写映射表
**v2 方案（优化后）**：利用 cpp2lean v2 现有管线，加 **Pass 9 (`refine_gen_pass`)**，复用全部基础设施

## 2. 核心流程

```
TRANSLATION_SCOPE (67 函数)
  ↓ Pass 1-8
MIRFunc 列表 ──→ Pass 9: refine_gen_pass ──→ Refinement/*.lean
                      ↑                                (全 sorry 骨架)
              REFINEMENT_MAP (新表, 在 class_map.py)
```

**改动点**：
- `class_map.py`：新增 `REFINEMENT_MAP` 表
- `pass9_refine_gen.py`：新 Pass 9，接收 `MIRFunc` 列表 → 输出 `.lean` 骨架
- `build_pass8_corpus.py`：末尾增加一行调用 Pass 9

**零改动**：Corpus.lean、Pass 1-8、ir_types.py

## 3. REFINEMENT_MAP 设计

`class_map.py` 中新增映射表，每条记录：

```python
REFINEMENT_MAP = {
    # L1 base_name → {L2 映射信息}
    "__squarefree_Zp": {
        "l2_name": "sqfZp",
        "l2_import": "CLPoly.Algorithm.SquarefreeZp",
        "refinement_file": "SquarefreeZp",
        "bridge": "toPoly",        # SparsePolyZp.toPoly
        "result_kind": "map_eq",    # L1 Array → L2 List
    },
    "__ddf_Zp": {
        "l2_name": "ddf",
        "l2_import": "CLPoly.Algorithm.DDF",
        "refinement_file": "DDF",
        "bridge": "toPoly",
        "result_kind": "map_eq",
    },
    "__edf_Zp": {
        "l2_name": "edf",
        "l2_import": "CLPoly.Algorithm.EDF",
        "refinement_file": "EDF",
        "bridge": "toPoly",
        "result_kind": "map_eq",
    },
    ...
}
```

**命名规则**：key 是 `base_name`（无 `_ir` 后缀），因为 Pass 9 从 `MIRFunc.base_name` 查表，再用 `MIRFunc.lean_name` 拼出定理名。

## 4. Pass 9 设计

### 4.1 功能规约

```
模块名称：pass9_refine_gen

功能描述：读取已翻译的 MIRFunc 列表 + REFINEMENT_MAP → 生成精化定理 Lean 骨架

前置条件：
  - MIRFunc 列表已完整（Pass 1-8 已跑完）
  - REFINEMENT_MAP 覆盖所有需要精化的入口函数

后置条件：
  - 输出 Refinement/<module>.lean 文件（每个文件含多个定理）
  - 每个定理为 `theorem <l1_name>_refines ... := by sorry`
  - 所有文件在 lake build 下可编译

副作用：写入文件系统（proof/lean/CLPoly/Refinement/）
```

### 4.2 输出文件结构

```
proof/lean/CLPoly/Refinement/
  Basic.lean                      # 手写桥引用理（不被 Pass 9 覆盖）
  SquarefreeZp.lean               # 自动生成
  DDF.lean
  EDF.lean
  Hensel.lean
  Recombine.lean
  Wang.lean
  FactorZp.lean
  FactorZZ.lean
  FactorMv.lean
```

每个自动生成文件的结构：

```lean
import CLPoly.Generated.Corpus
import CLPoly.Algorithm.<Module>
import CLPoly.Refinement.Basic
import CLPoly.Math.Univariate

set_option autoImplicit false

open Polynomial

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- L1 <func> → L2 <func> -/
theorem <l1_name>_refines (f : SparsePolyZp)
    (hwf : SparsePolyZp.WellFormed p f)
    (hred : SparsePolyZp.AllReduced p f)
    (hp_size : 2 * p ≤ UInt64.size)
    : (Generated.<l1_lean_name> f).map (·.toPoly p) = <l2_name> (f.toPoly p) :=
by
  sorry

end Refinement
```

### 4.3 四类 result_kind

| result_kind | 定理结论模板 | 适用场景 |
|---|---|---|
| `direct_eq` | `x = y` | 纯数值函数（`__make_zp`, `__binomial`） |
| `map_eq` | `(l1_ir ...).map (·.toPoly p) = l2 (...)` | 多项式运算（`__squarefree_Zp`, `__ddf_Zp`） |
| `conditional` | `l1_ir ... = .ok r → ...` | 有 Except/throw 的函数 |
| `pair_eq` | `(l1_ir ...).1.toPoly = (l2 ...).1 ∧ ...` | 多返回值函数 |

### 4.4 实现细节

Pass 9 需要从 `MIRFunc` 提取的信息：
- `lean_name`：L1 函数名（如 `__squarefree_Zp_ir`）
- `base_name`：查 REFINEMENT_MAP 的 key
- `params`, `ret_ty`：函数签名（用于骨架签名精度）
- `requires`：require 条件（用于 conditional 形态）

## 5. 与现有代码的关系

### 5.1 修改文件

| 文件 | 改动 |
|---|---|
| `proof/cpp2lean_v2/class_map.py` | 新增 `REFINEMENT_MAP` 表 |
| `proof/cpp2lean_v2/tests/build_pass8_corpus.py` | 末尾加 Pass 9 调用 |

### 5.2 新建文件

| 文件 | 用途 |
|---|---|
| `proof/cpp2lean_v2/passes/pass9_refine_gen.py` | Pass 9 实现 |
| `proof/lean/CLPoly/Refinement/Basic.lean` | 手写 bridge 引理 |
| `proof/lean/CLPoly/Refinement/*.lean` | 自动生成骨架 |

## 6. 与 v1 方案对比

| 维度 | v1（独立生成器） | v2（Pass 9 集成） |
|---|---|---|
| 基础设施复用 | 需从零写签名提取 | 直接复用 MIRFunc |
| 类型一致性 | 需手动同步 Corpus 类型 | 自动一致（同一源） |
| 维护成本 | 两份映射表 | 一份 REFINEMENT_MAP |
| 编译验证 | 需独立验证 | 随 lake build 自动验证 |
| 复杂度 | 中（~200 行） | 低（~100 行） |
