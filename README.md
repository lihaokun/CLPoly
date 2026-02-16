# CLPoly

CLPoly 是一个 C++ 多项式计算库，提供多变量多项式的算术运算、GCD、因式分解、结式、特征列、实根隔离等功能。

CLPoly is a C++ library for polynomial computation, supporting multivariate polynomial arithmetic, GCD, factorization, resultants, characteristic sets, and real root isolation.

## 功能

### 多项式算术

支持 Z[x₁,...,xₙ] 和 Q[x₁,...,xₙ] 上的多项式运算，提供自然的运算符重载：

```cpp
#include <clpoly/clpoly.hh>
using namespace clpoly;

variable x("x"), y("y");
polynomial_ZZ f = pow(x,3) - 2*x*y + 1;
polynomial_ZZ g = pow(x,2) + y;

auto h = f * g;              // 乘法
auto r = prem(f, g, x);      // 伪余式
auto lc = leadcoeff(f, x);   // 首项系数
auto d = derivative(f, x);   // 导数
auto v = assign(f, x, ZZ(3)); // 代入求值
```

### GCD 与无平方分解

```cpp
auto g = gcd(f1, f2);                  // 多项式 GCD
auto sqf = squarefreefactorize(f);     // 无平方分解: [(因子, 重数), ...]
auto basis = squarefreebasis({f1, f2}); // 无平方基
```

### 因式分解

Z[x] 和 Q[x] 上的单变量不可约分解（Hensel 提升 + Zassenhaus 重组）：

```cpp
polynomial_ZZ f = pow(x,6) - 1;
auto fac = factorize(f);
// fac.content = 1
// fac.factors = [(x-1,1), (x+1,1), (x²-x+1,1), (x²+x+1,1)]

// 也支持 upolynomial_ZZ 和 polynomial_QQ
auto fac_u = factorize(uf);
auto fac_q = factorize(fq);
```

### 结式与判别式

```cpp
auto res  = resultant(f, g, x);       // 结式
auto disc = discriminant(f, x);        // 判别式
auto src  = subresultant(f, g, x);     // 子结式链
```

### 特征列

基于吴方法的特征列与正则三角分解：

```cpp
auto cs = charset(polys);              // 基本特征列
auto wcs = Wucharset(polys);           // 吴特征列
auto rcs = GRDforZD(polys);            // 正则三角分解
```

### 实根隔离

基于 Uspensky 算法的实数根隔离：

```cpp
upolynomial_ZZ f = ...;
auto intervals = uspensky(f);          // 返回隔离区间 [(left,right), ...]
auto roots = realroot({f1, f2, f3});   // 多多项式联合实根隔离
```

### 数值类型

| 类型 | 说明 |
|------|------|
| `ZZ` | 任意精度整数（基于 GMP，小整数内联优化） |
| `QQ` | 有理数 |
| `Zp` | 素域 Z/pZ 元素 |

### 单项式序

| 序 | 类型 | 说明 |
|----|------|------|
| grlex | `grlex` | 分次字典序（默认） |
| lex | `lex` | 字典序 |
| 自定义 | `lex_<custom_var_order>` | 自定义变量优先级 |

## 依赖

- [GMP](https://gmplib.org/)（含 gmpxx）
- [Boost](https://www.boost.org/)（math/prime）

## 构建与测试

```bash
# 构建库
make

# 编译并运行单个测试
make test/test_factorize && test/test_factorize

# 运行全部测试（18 个测试文件，~1000 assertions）
bash test/run_all_tests.sh

# 交叉验证测试（需要 FLINT 和 NTL）
make crosscheck
```

### 在项目中使用

```bash
# 编译
g++ -O3 -I/path/to/CLPoly your_code.cc -lgmpxx -lgmp -Llib/clpoly -lclpoly

# 单头文件引入
#include <clpoly/clpoly.hh>
```

## 示例

```cpp
#include <clpoly/clpoly.hh>
using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");
    polynomial_ZZ f = -6*pow(x,4) - 7*pow(x,2)*y - x*pow(y,3)*z - 3*x*y - 9;
    polynomial_ZZ g = 10*pow(x,4)*y + 2*pow(x,4) - 7*pow(x,3)*pow(y,2)
                      + 9*pow(x,2) + 9*pow(z,4) - 9*y + 2*z;

    std::cout << "f = " << f << std::endl;
    std::cout << "g = " << g << std::endl;
    std::cout << "prem(f,g,x) = " << prem(f, g, x) << std::endl;
    std::cout << "resultant(f,g,x) = " << resultant(f, g, x) << std::endl;
    std::cout << "gcd(f^2, g*f) = " << polynomial_GCD(f*f, g*f) << std::endl;

    // 因式分解
    polynomial_ZZ h = pow(x,6) - 1;
    auto fac = factorize(h);
    std::cout << "factorize(x^6-1): content=" << fac.content << std::endl;
    for (auto& [fi, ei] : fac.factors)
        std::cout << "  (" << fi << ")^" << ei << std::endl;
}
```

```bash
make example
```

## 与其他库的对比

| | CLPoly | FLINT | NTL | Singular |
|--|--------|-------|-----|----------|
| 语言 | C++ | C | C++ | C/C++ |
| 多变量多项式 | 支持 | 支持 | 不支持 | 支持 |
| 单项式序 | grlex/lex/自定义 | deglex/degrevlex | — | dp/lp/Dp 等 |
| 系数域 | Z, Q, Z/pZ | Z, Q, Z/pZ, Z/nZ, Fq | Z, Z/pZ, GF(p^k) | Z, Q, Z/pZ, GF(p^k) |
| 因式分解 | 单变量 Z/Q | 单变量+多变量 | 单变量 Z | 单变量+多变量 |
| GCD | 多变量 | 多变量 | 单变量 | 多变量 |
| 结式/子结式 | 多变量 | 单变量 | 不支持 | 不支持 |
| 特征列 | 支持 | 不支持 | 不支持 | 不支持 |
| 实根隔离 | 支持 | 不支持 | 不支持 | 不支持 |
| Gröbner 基 | 不支持 | 支持 | 不支持 | 支持 |
| 内存布局 | AoS（交织） | SoA（分离） | 稠密数组 | AoS |
| 设计取向 | 泛型模板 | 极致性能 | 数论专用 | 交换代数 |

**CLPoly 的特点**：特征列方法与实根隔离的组合在开源 C++ 库中较为少见。泛型模板设计使得同一代码可工作于不同单项式序和系数域。因式分解目前限于单变量，多变量因式分解是后续工作。

## 许可证

MIT License
