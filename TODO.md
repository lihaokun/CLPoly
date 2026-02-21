# CLPoly TODO

## 性能优化

### P1: van Hoeij LLL 重组（单变量因式分解）

当前单变量 `__factor_recombine` 使用 Zassenhaus 子集枚举，O(2^r) 复杂度。
所有主流 CAS（FLINT、NTL、Maple、Magma）均已采用 van Hoeij LLL 方案，复杂度 O(poly(r))。

当前差距（多因子场景）：

| 用例 | CLPoly | FLINT | 比率 |
|------|--------|-------|------|
| 单变量 70 因子 | 80 ms | 3.6 ms | 22x |
| Wilkinson W(15) | 4.5 ms | 0.4 ms | 10x |

- [ ] 实现或引入 LLL 格规约（可用 fplll 库或基于 GMP 自写）
- [ ] Newton 幂和迹计算（从模因子）
- [ ] 格矩阵构造 + 短向量提取 → 因子分组
- [ ] 替换 `__factor_recombine`

参考：van Hoeij 2002, Hart-van Hoeij-Novocin 2011

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
