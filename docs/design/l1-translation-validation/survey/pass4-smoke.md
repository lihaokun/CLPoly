# Pass 1 + 2 + 3 + 4 全量烟测

- 目标函数数（TRANSLATION_SCOPE）：**66**（factorize 展开 3 实例 → 68 HIRs）
- OK: **68** / FAIL: **0**

- 识别 filter-loop 总数: **0**（预期 4）
- 残留 `RangeForStmt.decomposition`: **0**（预期 0）
- 有 filter-loop 的宿主: **0**

## 核对

⚠️ 差异：filter 0 vs 预期 4；残留 decomp 0