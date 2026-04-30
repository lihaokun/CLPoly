# Pass 1-8 codegen 全量烟测

- 目标 mir1：**67**（factorize 展开 3 实例）
- OK: **67** / FAIL: **0**
- 总 sorry：**495**（avg 7.4/func）
- 总输出行数：7217（avg 107.7/func）
- 输出目录：`/tmp/v2_lean_dump/`

## sorry 数量 Top 15（残留 hot-spot）
- `__vanhoeij_recombine`: sorry=32, lines=292
- `__wang_leading_coeff`: sorry=30, lines=392
- `__factor_multivar`: sorry=29, lines=285
- `__lll_reduce`: sorry=27, lines=432
- `__mtshl_lift`: sorry=27, lines=340
- `__wang_core`: sorry=27, lines=451
- `__select_eval_point`: sorry=19, lines=282
- `__mtshl_sparse_int`: sorry=17, lines=410
- `__mtshl_step_j`: sorry=17, lines=312
- `__zassenhaus_recombine`: sorry=17, lines=237
- `factorize_upoly`: sorry=17, lines=140
- `__hensel_lift`: sorry=16, lines=105
- `__hensel_step`: sorry=15, lines=145
- `factorize_lex`: sorry=15, lines=144
- `__mtshl_wmds`: sorry=13, lines=310