# L1 — `__upoly_random`（含 rng 状态前进 bug 修复）

**审视日期**：2026-05-04
**模块**：单变量 Zp 多项式随机生成
**状态**：✅（修复后）

---

## C++ 源（`polynomial_factorize_zp.hh:38-50`）

```cpp
inline void __upoly_random(upolynomial_<Zp>& result, int max_deg, uint64_t p, std::mt19937& rng) {
    std::uniform_int_distribution<uint64_t> dist(0, p-1);
    for (int d = max_deg - 1; d >= 0; --d) {
        uint64_t c = dist(rng);  // dist(rng) 改 rng 内部状态（mt19937）
        if (c != 0) result.push_back({d, Zp(c, p)});
    }
}
```

## 发现的偏差

**初始翻译**（修前）：

```lean
partial def _loop___upoly_random_0_ir (..) (rng : Rng) :=
  if (d_2 >= 0) then
    let c_1 : UInt64 := (Rng.next rng dist_1)   -- ❌ 纯函数，rng 不变
    if c_1 != 0 then
      ... (Array.push ... c_1 ...)
      _loop___upoly_random_0_ir d_2 rng ...   -- ❌ rng 没变，下轮 c 完全相同
```

**问题根因**：Lean 的 `Rng.next : Rng → UInt64 → UInt64` 是纯函数，不返回更新后的 seed；C++ 的 `dist(rng)` 改写 rng 状态。所有循环迭代里 `c` 取相同值，结果不是真正的随机多项式。

## 修复（Plan A，方案 1）

1. **Lean Model**（`Model.lean`）：新增 `Rng.next_advance`：
   ```lean
   def next_advance (seed upper : UInt64) : UInt64 × UInt64 :=
     (next seed upper, step seed)  -- 返回 (随机值, 推进后的 seed)
   ```

2. **Pass 5**（`pass5_operator_resolve.py`）：`dist(rng)` → `Rng.next_advance rng dist`，emit ty=UInt64（orig_ret），让 Pass 2b 自动 wrap 成 PairType(UInt64, Rng)

3. **class_map**（`class_map.py`）：`TRANSLATION_SCOPE_OUTPUT_PARAMS["Rng.next_advance"] = [0]`（rng 是 ref-out at idx 0）

4. **Pass 2b**（`pass2b_callsite_ref_elim.py`）：
   - `_hoist_ref_call` 升级为递归 walker，处理嵌套表达式中的 ref-out Call（如 `Array.set! arr k (Zp.ofInt (Rng.next_advance gen dist).toInt p)`）
   - 新增 `callee_filter` 参数，rerun 模式只 hoist 指定 callee

5. **build_pass8_corpus.py**：在 Pass 5 之后加 Pass 2b rerun，传 `callee_filter={"Rng.next_advance"}`，避免破坏首轮已 destructure 的 `fdiv_q`/`__upoly_divmod` 等 callsite

## 修复后翻译

```lean
partial def _loop___upoly_random_0_ir (..) (rng_1 : Rng) :=
  if (d_2 >= 0) then
    let __refret_0_1 : (UInt64 × Rng) := (Rng.next_advance rng_1 dist_1)
    let rng_2 : Rng := __refret_0_1.snd                   -- ✅ 推进后的 rng
    let c_1 : UInt64 := __refret_0_1.fst                  -- 随机值
    if (c_1 != 0) then
      ...
      _loop___upoly_random_0_ir d_2 rng_2 ...             -- ✅ 传新 rng
```

## 行为对照

| C++ | Lean | 一致 |
|-----|------|------|
| `std::uniform_int_distribution<uint64_t> dist(0, p-1)` | `let dist_1 := UniformIntDist.mk 0 (p-1)` | ✅ |
| `for (int d = max_deg-1; d >= 0; --d)` | for-loop lower 为 `_loop_..._0_ir` 递归 + `d_3 := d_2 - 1` | ✅ |
| `uint64_t c = dist(rng)`（mut rng）| `Rng.next_advance rng_1 dist_1 → (c_1, rng_2)` | ✅ |
| `if (c != 0) result.push_back({d, Zp(c, p)})` | `if c_1 != 0 then Array.push ... (UMonomial.mk d_2, Zp.ofInt c_1.toInt p)` | ✅ |

## 偏差

无（修复后）。

## 关联修复影响

`Rng.next_advance` 同样应用于 `_loop___mtshl_sparse_int_lex_1_ir`（line 2700 nested 形式），同样修复成功。

## 测试结果

`lake build CLPoly.Generated.Corpus`：✅ **0 errors / 0 warnings**（all 346 functions pass）

修复前：18 errors（Type mismatch: PairType vs UInt64 等）
