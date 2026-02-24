# CLPoly

CLPoly 是一个 C++ 多项式计算库，提供多变量多项式的算术运算、GCD、因式分解、结式、特征列、实根隔离等功能。

CLPoly is a C++ library for polynomial computation, supporting multivariate polynomial arithmetic, GCD, factorization, resultants, characteristic sets, and real root isolation.

## 特性概览

| 能力 | 说明 |
|------|------|
| 多项式算术 | Z[x₁,…,xₙ] 和 Q[x₁,…,xₙ] 上的加减乘、伪除、求导、代入 |
| GCD | 多变量多项式 GCD，无平方分解与无平方基 |
| 因式分解 | 单变量 + 多变量不可约分解（Hensel 提升 + Zassenhaus / Wang EEZ） |
| 结式 | 结式、判别式、子结式链 |
| 特征列 | 吴方法特征列、正则三角分解（开源 C++ 库中少见） |
| 实根隔离 | 基于 Uspensky 算法的单/多多项式实根隔离 |
| 数值类型 | `ZZ`（任意精度整数，小整数内联优化）、`QQ`（有理数）、`Zp`（素域） |
| 单项式序 | grlex（默认）、lex、自定义变量优先级 |
| 设计 | 泛型模板，同一代码工作于不同单项式序和系数域 |

## 快速开始

```cpp
#include <clpoly/clpoly.hh>
using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");
    polynomial_ZZ f = -6*pow(x,4) - 7*pow(x,2)*y - x*pow(y,3)*z - 3*x*y - 9;
    polynomial_ZZ g = 10*pow(x,4)*y + 2*pow(x,4) - 7*pow(x,3)*pow(y,2)
                      + 9*pow(x,2) + 9*pow(z,4) - 9*y + 2*z;

    std::cout << "gcd(f^2, g*f) = " << polynomial_GCD(f*f, g*f) << std::endl;
    std::cout << "resultant(f,g,x) = " << resultant(f, g, x) << std::endl;

    // 因式分解
    polynomial_ZZ h = pow(x,2) - pow(y,2);
    auto fac = factorize(h);
    for (auto& [fi, ei] : fac.factors)
        std::cout << "  (" << fi << ")^" << ei << std::endl;
    // 输出: (x-y)^1  (x+y)^1
}
```

## 依赖

- [GMP](https://gmplib.org/)（含 gmpxx）
- [Boost](https://www.boost.org/)（math/prime）

## 构建与测试

```bash
# 构建全部库（debug/release × .a/.so → lib/）
make

# 编译并运行单个测试
make test/test_factorize && _build/debug/bin/test_factorize

# 运行全部测试（23 个测试文件，~1950 assertions）
bash test/run_all_tests.sh

# 性能基准（release 构建，fork 子进程隔离，同时测量时间和内存）
make bench              # 运行全部基准（bench-clpoly + bench-comparative）
make bench-clpoly       # 纯 CLPoly 基准
make bench-comparative  # CLPoly vs FLINT/NTL 对比基准

# 交叉验证测试（需要 FLINT 和 NTL）
make crosscheck
```

构建产物布局：

```
lib/                        # 公开库输出（-Llib -lclpoly）
├── libclpoly.a/.so         # release
└── debug/libclpoly.a/.so   # debug
_build/                     # 内部构建产物
├── debug/obj/              # debug 对象文件
├── debug/bin/              # 测试二进制
└── release/obj|bin/        # release 对象文件与基准二进制
```

### 在项目中使用

```bash
g++ -O3 -I/path/to/CLPoly your_code.cc -Llib -lclpoly -lgmpxx -lgmp
```

```cpp
#include <clpoly/clpoly.hh>
```

## 与其他库的对比

| | CLPoly | FLINT | NTL | Singular | SymPy |
|--|--------|-------|-----|----------|-------|
| 语言 | C++ | C | C++ | C/C++ | Python |
| 多变量多项式 | 支持 | 支持 | 不支持 | 支持 | 支持 |
| 单项式序 | grlex/lex/自定义 | deglex/degrevlex | — | dp/lp/Dp 等 | lex/grlex 等 |
| 系数域 | Z, Q, Z/pZ | Z, Q, Z/pZ, Z/nZ, Fq | Z, Z/pZ, GF(p^k) | Z, Q, Z/pZ, GF(p^k) | Z, Q, Z/pZ, 代数数域 |
| 因式分解 | 单变量+多变量 Z/Q | 单变量+多变量 | 单变量 Z | 单变量+多变量 | 单变量+多变量 |
| GCD | 多变量 | 多变量 | 单变量 | 多变量 | 多变量 |
| 结式/子结式 | 多变量 | 单变量 | 不支持 | 不支持 | 多变量 |
| 特征列 | 支持 | 不支持 | 不支持 | 不支持 | 不支持 |
| 实根隔离 | 支持 | 不支持 | 不支持 | 不支持 | 支持 |
| Gröbner 基 | 支持 | 支持 | 不支持 | 支持 | 支持 |
| 内存布局 | AoS（交织） | SoA（分离） | 稠密数组 | AoS | Python 对象树 |
| 设计取向 | 泛型模板 | 极致性能 | 数论专用 | 交换代数 | 通用 CAS |

## 核心功能与 API 参考

### 数值类型

```cpp
// ZZ — 任意精度整数（小整数 ≤ 63 bit 内联存储，自动升级到 GMP）
ZZ a(42), b("-99999999999999999999");
ZZ c = a * b;                       // 算术：+  -  *  /  %
bool s = a.is_zero();               // 查询：is_zero()  is_one()  is_odd()  sgn()
ZZ g = gcd(a, b);                   // gcd / lcm / abs / pow

// QQ — 有理数（自动约分）
QQ q(3, 7);                         // 3/7
QQ r = q + QQ(1, 2);               // 算术：+  -  *  /
const ZZ& num = q.get_num();        // 分子/分母访问

// Zp — 素域 Z/pZ
Zp x(5, 13);                        // 5 mod 13
Zp y = x.inv();                     // 模逆元
```

### 变量与单项式序

```cpp
variable x("x"), y("y"), z("z");

// 支持的单项式序
polynomial_<ZZ, grlex> f;            // grlex（默认）：先比总次数，再反字典序
polynomial_<ZZ, lex>   g;            // lex：纯字典序
// 类型别名
polynomial_ZZ  h;                    // = polynomial_<ZZ, grlex>
polynomial_QQ  q;                    // = polynomial_<QQ, grlex>
upolynomial_ZZ u;                    // 单变量多项式

// 序之间可转换
poly_convert(f, g);                  // grlex → lex
```

### 多项式算术

```cpp
variable x("x"), y("y");
polynomial_ZZ f = pow(x,3) - 2*x*y + 1;
polynomial_ZZ g = pow(x,2) + y;

// 基本运算
auto h = f * g;                      // 加、减、乘
auto p = pow(f, 5);                  // 幂

// 伪除法
auto r = prem(f, g, x);             // 伪余式
polynomial_ZZ q_out, r_out;
pquo(q_out, r_out, f, g, x);        // 伪商 + 伪余式

// 查询
auto d = degree(f);                  // 总次数
auto dx = degree(f, x);             // 关于 x 的次数
auto lc = leadcoeff(f, x);          // 关于 x 的首项系数
auto cf = coeff(f, x);              // 提取 x 的各次系数
bool is_c = is_number(f);           // 是否为常数
auto vars = get_variables(f);       // 获取变量列表 [(var, max_deg), ...]

// 求导与求值
auto df = derivative(f, x);         // 偏导数 ∂f/∂x
auto v1 = assign(f, x, ZZ(3));      // 代入 x=3
auto v2 = assign(f, {{x, ZZ(1)}, {y, ZZ(2)}}); // 多变量代入
```

### GCD 与无平方分解

```cpp
auto g = gcd(f1, f2);               // 多项式 GCD（多变量 / 单变量均可）
auto sqf = squarefreefactorize(f);   // 无平方分解: [(因子, 重数), ...]
auto basis = squarefreebasis({f1, f2}); // 无平方基
bool sf = is_squarefree(f);          // 判断是否无平方
```

### 因式分解

支持 `polynomial_<ZZ, lex>`、`polynomial_<ZZ, grlex>`、`polynomial_QQ`、`upolynomial_ZZ`。

返回 `factorization<Poly>`，包含 `content`（内容）和 `factors`（`[(因子, 重数), ...]`）。

```cpp
// 单变量
upolynomial_ZZ f({{6, ZZ(1)}, {0, ZZ(-1)}});  // x^6 - 1
auto fac = factorize(f);
// fac.content = 1
// fac.factors = [(x-1,1), (x+1,1), (x²-x+1,1), (x²+x+1,1)]

// 多变量（lex 或 grlex 均可）
polynomial_ZZ g = pow(x,2) - pow(y,2);
auto fac2 = factorize(g);
// fac2.factors = [(x-y,1), (x+y,1)]

// 有理数域
polynomial_QQ q = ...;
auto fac3 = factorize(q);           // 内容为 QQ，因子为 polynomial_QQ
```

**算法**：单变量 Berlekamp + Hensel + Zassenhaus；多变量 Wang EEZ（LC 分配 + 多变量 Hensel 提升 + 子集枚举试除）。

### 结式与判别式

```cpp
auto res  = resultant(f, g, x);      // 结式 res_x(f, g)
auto disc = discriminant(f, x);      // 判别式 disc_x(f)
auto src  = subresultant(f, g, x);   // 子结式链 [S_0, S_1, ...]
```

### Gröbner 基

支持 Z[x₁,…,xₙ] 和 Q[x₁,…,xₙ]。算法：Buchberger + Sugar 策略 + Gebauer-Möller 准则。

```cpp
variable x("x"), y("y");
polynomial_ZZ f1 = pow(x,2) + y - 1;
polynomial_ZZ f2 = x + pow(y,2) - 1;
std::vector<polynomial_ZZ> gens = {f1, f2};

// 直接计算
auto gb = groebner_basis(gens);      // 返回约化 Gröbner 基

// 面向对象接口
Ideal<ZZ> I({f1, f2});
auto& gb2 = I.groebner_basis();     // 惰性计算，缓存结果
bool mem  = I.contains(f1 * f2);    // 理想成员判定
auto nf   = I.reduce(f1 + f2);      // 对 Gröbner 基求标准形

// S-多项式与标准形
auto sp = s_polynomial(f1, f2);      // S(f1, f2)
auto nf2 = normal_form(sp, gb);     // 对基约化
```

### 特征列

```cpp
auto cs  = charset(polys);           // 基本特征列
auto wcs = Wucharset(polys);         // 吴特征列
auto rcs = GRDforZD(polys);          // 正则三角分解
```

### 实根隔离

基于 Uspensky 算法（Descartes 符号规则 + 二分）。

```cpp
// 单多项式
upolynomial_ZZ f = ...;
auto intervals = uspensky(f);        // 返回隔离区间 [(left, right), ...]

// 多多项式联合实根隔离
auto roots = realroot({f1, f2, f3});

// 实根辅助
ZZ bound = RealRootBound(f);         // Cauchy 根界
int sc = coeffsignchanges(f);        // 系数符号变化数（Descartes 规则）
```

## 许可证

MIT License
