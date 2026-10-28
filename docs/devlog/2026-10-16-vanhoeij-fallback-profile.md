# van Hoeij 兼容入口回退计数

日期：2026-10-16。

最终验收核对发现，纯 van Hoeij 主路径已取消 Zassenhaus 回退，但固定精度兼容入口仍保留回退，其调用未增加 profiling 计数。这会使路径统计漏报。

在兼容入口实际调用 Zassenhaus 前增加一次计数；增加固定精度不足、恢复界已满足的二次因子用例，验证拒绝证书、准确恢复两个因子且回退恰好一次。纯主路径继续要求证书成功且零回退。

修改仅影响 `CLPOLY_PROFILE` 编译配置。最终 profiling 原生测试 406/406 通过（日志 `/tmp/clpoly-vh-final-profile-native.log`）；不带该宏的 27 文件原生回归与正在运行的最终 B2B 对应同一算法。三版本性能快照需重新生成，以包含这次统计修正。

涉及 [生产实现](../../clpoly/polynomial_factorize_univar.hh)、[回归测试](../../test/test_recombine.cc)及[性能验收脚本](../../test/run_vanhoeij_benchmark.py)。性能脚本同时检查独立计时器记录的 LLL 调用大于零、内外 Zassenhaus 调用均为零，避免仅依赖生产计数。

## 度量

- 新增原生断言：5；总计 406，失败 0。
- Lean 变更：0；计数不参与算法控制流。
- 放弃方案：无。
