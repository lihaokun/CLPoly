# Pass 1 + 2 + 3 + 4 全量烟测

- 目标函数数（TRANSLATION_SCOPE）：**65**（factorize 展开 3 实例 → 67 HIRs）
- OK: **67** / FAIL: **0**

- 识别 filter-loop 总数: **4**（预期 4）
- 残留 `RangeForStmt.decomposition`: **0**（预期 0）
- 有 filter-loop 的宿主: **2**

## Filter-loop 识别详情

| 宿主函数 | 识别数 |
|---|---|
| `__extract_monomial_content` | 2 |
| `__hensel_step` | 2 |

## 核对

✅ filter-loop = 预期 & decomposition 全清空