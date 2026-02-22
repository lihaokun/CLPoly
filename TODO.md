# CLPoly TODO

## 性能优化

### P1: 单变量因式分解优化

#### P1a: van Hoeij LLL 重组（**已完成**）

**实现结果**（2026-02-22 benchmark）：

| 用例 | 实现前 | 实现后 | 提升 |
|------|--------|--------|------|
| 单变量 70 因子 | 80 ms | 77 ms | 无回退 ✓ |
| Wilkinson W(15) | 13.9 ms | 4.4 ms | 3.2x |
| x^24-1 (cyclotomic) | 24.6 ms | 9.2 ms | 2.7x |

架构文档：`docs/design/vanhoeij/architecture.md`
调研报告：`docs/research/vanhoeij-lll-research.md`

实现任务（全部在 `clpoly/polynomial_factorize_univar.hh`）：

- [x] **M1** `__cld_polys`：CLD 多项式 C_i = (f/h_i)·h_i' mod m
- [x] **M2** `__build_cld_matrix`：内螺旋列喂入格矩阵（初版无列过滤）
- [x] **M3** `__lll_reduce`：整数 LLL（Cohen §2.6，QQ Gram-Schmidt）
- [x] **M4** `__extract_candidates`：U_short 列等价类分组
- [x] **M5** `__vanhoeij_recombine`：主控循环（批量提取 + 对角 LLL 预通道）
- [x] `__lll_factorize`：入口函数（二次 Hensel 提升 + van Hoeij LLL，全 r 路由）

**关键设计**：
- 初始 J_target=0（对角 LLL）免费执行 s=1 单因子检验，避免全分裂多项式触发昂贵 CLD 计算
- 批量提取：单次 LLL 调用中处理所有候选，将 r 次重复 LLL 降为 O(1) 次

#### P1b: 线性 Hensel 提升（调研完成，实现保留为基础设施）

调研与设计（Monagan 2019）已完成：`docs/research/linear-hensel-research.md`、`docs/design/vanhoeij/detailed-design-p1b.md`

- [x] 线性 Hensel 基础设施：`__hensel_step_linear`、`__hensel_lift_linear_recursive`（保留供未来真正线性交织优化）
- [x] `__heuristic_starting_precision`：FLINT 启发式精度估算

**结论**：van Hoeij LLL 需在满精度（Mignotte 界）下运行才能高效收敛；低精度 LLL 导致 O(n^4) 规约退化。P1a 的二次 Hensel + van Hoeij LLL 已达到实用性能，P1b 完整实现（精度自适应 + 早期因子检测）可作为后续优化项。

参考：van Hoeij 2002, Belabas 2009, Hart 2011（P1a）；Monagan 2019（P1b）

### P2: 多变量因式分解算法选择与稀疏路径

多变量因式分解当前使用 Wang EEZ + Zassenhaus 重组，多因子时性能极差。

| 用例 | CLPoly | FLINT | 比率 |
|------|--------|-------|------|
| 双变量 70 因子 | 71 s | 10 ms | 7200x |
| 三变量 60 因子 | 6.7 s | 1.3 ms | 5100x |

当前 CLPoly 无论密度一律走 Wang + Zassenhaus。需要：1) 实现稀疏算法路径；2) 根据密度/变量数自动分流。

**稀疏算法候选：**

- **Zippel 稀疏插值**（FLINT 默认）：在多个随机点求值 → 单变量因式分解 → Vandermonde 插值重建多变量因子。完全绕开重组问题，复杂度 O(t³)，t = 项数。参考：Zippel 1979/1993, FLINT `irred_zippel.c`
- **MTSHL 稀疏 Hensel lift**（Maple 2019+）：计算大量双变量像 → 双变量 Hensel lift → 稀疏插值重建多变量因子。避免多变量 Diophantine 方程，高变量数时可能更优。参考：Monagan-Tuncer 2019 (ISSAC)

**算法分流策略**（参考 FLINT/Maple）：稠密 → Wang EEZ，稀疏 → Zippel/MTSHL，回退 → Zassenhaus。

- [ ] 评估 Zippel vs MTSHL 的适用场景
- [ ] 实现选定的稀疏算法作为 Wang EEZ 的备选路径
- [ ] 实现密度计算 + 根据密度/变量数自动选择算法路径

详见：`docs/research/polynomial-systems-comparison.md` §4, `docs/design/multivariate-factorization-design.md`

### P3: 多项式乘法优化

- [ ] FFT / Karatsuba 大次数单变量乘法
- [ ] 多变量乘法堆排序优化

### P4: GCD 性能

当前多变量 GCD 比 FLINT 慢 3-9x。

- [ ] 调研 FLINT 的模 GCD 实现
- [ ] 考虑稀疏模 GCD（SPMOD）或 Zippel GCD

## 文档完善

- [ ] API 文档：为核心头文件添加 Doxygen 注释
- [ ] 算法文档：补充各模块采用的算法描述和参考文献
- [ ] 示例代码：在 `examples/` 目录添加典型用例
- [ ] 性能指南：记录各操作的复杂度和适用场景
