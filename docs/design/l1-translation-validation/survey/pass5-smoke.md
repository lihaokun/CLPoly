# Pass 1 + 2 + 3 + 4 + 5 全量烟测

- 目标函数数：**65**（factorize 3 实例 → 67 HIRs）
- OK: **67** / FAIL: **0**
- Gap 总数: **2**（B 策略残余，Pass 8 codegen 时输出 sorry）

## Gap 分类

| 类别 | 数量 |
|---|---|
| `method_miss` | 2 |

## Top 20 Gap 详情

| 次数 | 类别 | 详情 |
|---|---|---|
| 2 | method_miss | `('BaseType.UNIT', 'operator bool')` |

## Per-function Gap 数

| 函数 | gap |
|---|---|
| `__hensel_step_linear` | 1 |
| `__upoly_mod_coeff` | 1 |
