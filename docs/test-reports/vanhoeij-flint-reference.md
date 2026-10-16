# 固定 FLINT 3.1.2 的同输入 van Hoeij 动态对照

日期：2026-10-16。研究分支：`research/vanhoeij-lattice-certificate`。

## 固定基线

- 运行库：本机 Nix `flint3-3.1.2`，不可变 store 路径
  `/nix/store/h4xmcb94mg0mjdiic8p9adfy2hcfdf7v-flint3-3.1.2`。
- 实际 `libflint.19.0.dylib` 的 SHA-256：
  `792ec384cc6497b83f0cad6c9ed638eb456b0d32fd8c0abfec53e144c09a5870`。
- [官方源码包](https://flintlib.org/download/flint-3.1.2.tar.gz)已下载，并与
  本机该版本 Nix 源码 derivation 的 SHA-256 核对一致：
  `fdb3a431a37464834acff3bdc145f4fe8d0f951dd5327c4c6f93f4cbac5c2700`。
- 库采用本机 Nix 构建；没有声称本轮重新从 tarball 构建了该二进制。

已核对该包的 `factor_van_hoeij.c`、`next_col_van_hoeij.c` 和
`van_hoeij_check_if_solved.c`。固定版本保留已有格状态、加入变换后的
约束列、使用 with-removal，并在未完成时继续提升精度。当前 CLPoly
每轮重建及固定虚拟位上限仍须自己的充分精度证明。

## 动态测试方法

[基准](../../test/bench_vanhoeij_certificate.cc) 增加可选
`CLPOLY_FLINT_REFERENCE` 模式。向两端直接提供同一个整数多项式、
同一素数 p=97 和同一组已分解模因子。FLINT 直接调用
`fmpz_poly_factor_van_hoeij`，不另行选择素数。

两端分别核对每个确切已知不可约因子和因子总数；FLINT 还检查重数
均为一。CLPoly 每次要求证书完成计数为一、Zassenhaus 回退为零。
覆盖六个不同输入，每个重复三次，共 **18/18 对照通过**：

- 首一正二次因子族，r=24/32，分别有 12/16 个真因子。
- 非首一正二次因子族，r=24/32；构造后取本原部分，负判别式保证不可约。
- 两个 Eisenstein(2) 高次因子，r=24/32，真因子次数分别为 12/16。

前两族均为两次提升、一次 LLL、两列、零回退。高次族 r=24 为两次
LLL、累计六列；r=32 为四次 LLL、累计三十列，均两次提升、零回退。
动态对照覆盖输入和最终因子；内部格状态只作源码核对，未采集 FLINT
的逐步内部执行 trace。

## 时间样本

同一进程依次运行 CLPoly 与 FLINT，单次调用前后计时；输入转换及
结果核验在计时外。单位毫秒，各项为三次中位数：

| 输入族 | r | CLPoly | FLINT 3.1.2 |
|---|---:|---:|---:|
| 首一二次 | 24 | 82.5105 | 1.10992 |
| 首一二次 | 32 | 201.599 | 2.32350 |
| 非首一二次 | 24 | 83.7097 | 1.16483 |
| 非首一二次 | 32 | 199.757 | 2.92313 |
| 两个高次因子 | 24 | 635.446 | 1.90833 |
| 两个高次因子 | 32 | 21072.4 | 8.52313 |

这些样本显示仍有很大的参考性能差距。高次 r=32 本轮 CLPoly 样本
为 21072.4、28348.9、19620.4，波动明显；参考模式还增加了动态库
及头文件，因此不将此表替代此前同编译配置的旧/新基线对比。
原始用户输入尚未取得，不能宣称已经复测其约 50% 退化。

## 复跑

本轮使用 Nix clang 19.1.7、GMP 6.3.0、Boost 1.87.0、MPFR 4.2.1，
`-std=c++17 -O2`，启用 `CLPOLY_PROFILE`、`CLPOLY_CERTIFICATE_BENCH`、
`CLPOLY_FLINT_REFERENCE`。编译基准并链接 `upolynomial.cc`、`variable.cc`、
FLINT 和 GMP；依赖 include/lib 路径使用上述已安装 Nix 输出。

```sh
_build/bench_vanhoeij_flint
_build/bench_vanhoeij_flint nonmonic
_build/bench_vanhoeij_flint hard
```

三份原始日志分别是 `/tmp/clpoly-vh-flint-default.log`、
`/tmp/clpoly-vh-flint-nonmonic.log`、`/tmp/clpoly-vh-flint-hard.log`。

固定版本的结果对照已通过。充分精度终止证明、实际源码与模型的最终
对应审查以及性能验收继续进行，当前仍不可合并。

## 2026-10-16：输入固定截断上界的诊断对照

将 D 限制为 D₀² 后，用相同编译选项比较修改前研究分支二进制与新版。
三族各 r=24/32、每项三次，新版 18/18 与 FLINT 3.1.2 及已知因子一致，
每次 certificate=1、fallback=0。此处旧版是增加截断上界之前的研究
实现，并非最初 PR 的 safety-net 基线。

| 输入族 | r | 加上界前中位数 ms | 加上界后中位数 ms |
|---|---:|---:|---:|
| 首一二次 | 24 | 84.6982 | 83.9992 |
| 首一二次 | 32 | 214.278 | 219.935 |
| 非首一二次 | 24 | 85.5645 | 87.1101 |
| 非首一二次 | 32 | 198.967 | 202.506 |
| 两个高次因子 | 24 | 599.879 | 599.865 |
| 两个高次因子 | 32 | 13252.1 | 12974.7 |

运行期间 B2B 的 Lean 解释器仍在后台计算，不能将上述时间用作最终
性能验收。原始汇总为 `/tmp/clpoly-vh-bounded-divisor-benchmark.json`；
六份日志为 `/tmp/clpoly-vh-bounded-divisor-{quadratic,nonmonic,hard}-{before,after}.log`。
