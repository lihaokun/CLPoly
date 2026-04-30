# Pass 1-8 codegen 全量烟测

- 目标 mir1：**67**（factorize 展开 3 实例）
- OK: **67** / FAIL: **0**
- 总 sorry：**181**（avg 2.7/func）
- 总输出行数：7217（avg 107.7/func）
- 输出目录：`/tmp/v2_lean_dump/`

## sorry 数量 Top 15（残留 hot-spot）
- `__factor_multivar`: sorry=18, lines=285
- `__hensel_step`: sorry=14, lines=145
- `__mtshl_lift`: sorry=10, lines=340
- `__vanhoeij_recombine`: sorry=10, lines=292
- `__wang_core`: sorry=9, lines=451
- `__extract_monomial_content`: sorry=8, lines=140
- `__select_eval_point`: sorry=8, lines=282
- `factorize_upoly`: sorry=8, lines=140
- `factorize_lex`: sorry=8, lines=144
- `__factor_Zp`: sorry=7, lines=78
- `__hensel_step_linear`: sorry=6, lines=65
- `__lll_reduce`: sorry=6, lines=432
- `__mtshl_sparse_int`: sorry=6, lines=410
- `__mtshl_step_j`: sorry=6, lines=312
- `__build_cld_matrix`: sorry=5, lines=108