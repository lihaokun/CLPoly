# P1 度数界穿线 — 轻量接口设计

日期：2026-07-26
目标：清 `Pipeline/L1.lean` 的 2 个 sqf 溢出 admit（h_no_overflow / h_deg_bound），使 sqf L1 路径无 sorry。

## 核心洞察

`sqfZp_l1_correct` 需前提 `hdeg : f.natDegree * p < 2^64`（C++ 模型对 coeff·deg 做 UInt64 乘）。
但 L2 管线胶水 `factor_Zp_correct.hsqf` / `factor_ZZ_correct.hfzp` 是 `∀ g` 全称——大度数下 C++ 真溢出，全称不可证。

**关键事实**：这两个 ∀-g 假设各自**只在一个具体多项式上被调用**：
- `factor_Zp_correct` 证明体只用 `hsqf f hf`（f 自身）。
- `factor_ZZ_correct` 证明体只用 `hfzp fp hfp_ne`，`fp = f.map (castRingHom (ZMod p))`，且 `hdeg : fp.natDegree = f.natDegree` 在手。

## 设计：把子过程假设弱化为「仅对实际使用的多项式」

数值度数界不进胶水接口，只在 L1 叶子调用点消解。L2 无条件定理不受影响。

| 定理 | 改动 |
|------|------|
| `sqfZp_l1_correct` (L1) | 加 `(hdeg : f.natDegree * p < 2^64)`，证 2 admit（helper 已备） |
| `factor_Zp_correct` (FactorZp) | `hsqf : ∀ g, g≠0 → SquarefreeDecomp g (sqf g)` → `hsqf : SquarefreeDecomp f (sqf f)`；证明体 `hsqf f hf` → `hsqf` |
| `factor_Zp_instantiate` (×2) | `(fun g hg => sqf_correct g hg)` → `(sqf_correct f hf)` |
| `factor_Zp_l1` (L1) | 加 `(hdeg : f.natDegree * p < 2^64)`；hsqf 参数 `(sqfZp_l1_correct hp_size hp2)` → `(sqfZp_l1_correct hp_size hp2 f hf hdeg)` |
| `factor_zp_l1_func` (L1) | 定义加度数界 if：`else if hb : g.natDegree*p<2^64 then (factor_Zp_l1 … g hg hb).choose… else (0,[])` |
| `factor_zp_l1_func_correct` (L1) | 加前提 `(hb : g.natDegree*p<2^64)`；证明 `simp [hg, hb]` 后取 choose_spec |
| `factor_ZZ_correct` (FactorZZ) | `hfzp : ∀ g, g≠0 → FactorZpCorrect g …` → `hfzp : FactorZpCorrect fp (factor_zp fp).1 (factor_zp fp).2`（fp 具体）；证明体 `hfzp fp hfp_ne` → `hfzp` |
| `factor_ZZ_instantiate` (FactorZZInstantiate) | hfzp 参数改为具体 `factor_zp_correct fp (map_ne_zero…)` |
| `factor_ZZ_cpp_correct` (L1) | 加 `(hfbound : f.natDegree * p < 2^64)`；hfzp 参数 `factor_zp_l1_func_correct … fp hfp_ne (由 hfbound + hdeg 得 fp 界)` |

ddf/edf 的 hddf/hedf 仍 `∀ g`（用于多个 sqf 分量），本次不动——ddf_l1/edf_l1 仍 sorry（P2），
其 ∀-g 类型满足未改的 hddf/hedf。P2 的 ddf/edf 精化将复用同一 `toSparsePolyZp_deg_le` + hdeg 模式。

## 验证
- `#print axioms sqfZp_l1_correct` 应无 sorryAx（条件版）。
- `factor_Zp_instantiate_unconditional` / `factor_ZZ_instantiate` 保持无 sorryAx（弱化不破坏）。
- `factor_ZZ_cpp_correct` 仍 sorry（ddf/edf/hensel/recombine 未完成），但 sqf 分量不再贡献 sorry。
- `lake build` 全绿。
