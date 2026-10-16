# van Hoeij 格下标范围修复

日期：2026-10-16。

源码核对发现格维数 r+j 可能超出 int，而因子数与列数分别仍合法。
将格行下标和计数改为 size_t，并使用正剩余数执行降序规约，保留
因子编号接口和算法顺序。6 项机器范围定理及扩大审计 212/212 已通过；最新重组回归
401/401，全套原生 27 文件通过。最新 3107 项 B2B 仍在运行。

详见[证明日志](../../proof/docs/devlog/2026-10-16-vanhoeij-index-width.md)与
[源码对应审查](../test-reports/vanhoeij-source-correspondence.md)。
