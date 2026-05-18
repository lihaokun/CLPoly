# Pass 1 + 2 + 3 + 4 + 5 全量烟测

- 目标函数数：**66**（factorize 3 实例 → 60 HIRs）
- OK: **60** / FAIL: **8**
- Gap 总数: **0**（B 策略残余，Pass 8 codegen 时输出 sorry）

## FAIL 列表

- `__lll_reduce`: TranslationError: aux(_lambda___lll_reduce_5)[8]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='dot', version=0, ty=NamedType(name='Lambda'))
- `__mtshl_multi_bdp`: TranslationError: body[20]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='compute_error', version=0, ty=NamedType(name='Lambda'))
- `__mtshl_sparse_int`: TranslationError: body[8]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='rd', version=0, ty=NamedType(name='Rng'))
- `__mtshl_step_j`: TranslationError: body[15]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='product_F', version=0, ty=NamedType(name='Lambda'))
- `__mtshl_wmds`: TranslationError: body[31]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='compute_error', version=0, ty=NamedType(name='Lambda'))
- `__select_prime`: TranslationError: body[13]/step[0]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='next_p', version=0, ty=NamedType(name='Lambda'))
- `__wang_core`: TranslationError: body[11]/body[0]/body[5]/then[4]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='all_div', version=0, ty=NamedType(name='Lambda'))
- `__wang_leading_coeff`: TranslationError: body[9]: P0-2 violation — Call.callee is Var (must be str or UnresolvedOp), got Var(name='make_const', version=0, ty=NamedType(name='Lambda'))

## Gap 分类

| 类别 | 数量 |
|---|---|

## Top 20 Gap 详情

| 次数 | 类别 | 详情 |
|---|---|---|

## Per-function Gap 数

✅ 0 gap，所有节点都已 resolve
