# 类型转换（Cast）枚举

总计 **5670** 次 cast（4 种 CastExpr kind）。按 (source, target, castKind) 三元组去重统计。

## castKind 直方图

| castKind | 次数 | 含义 |
|---|---|---|
| `FunctionToPointerDecay` | 1786 | 函数名到函数指针 |
| `NoOp` | 1630 | 无操作（类型同质，可消除） |
| `LValueToRValue` | 1263 | 左值转右值（读取） |
| `IntegralCast` | 704 | 整数类型转换（需检查位宽/符号） |
| `ConstructorConversion` | 142 | 调构造函数 |
| `BuiltinFnToFnPtr` | 53 |  |
| `ArrayToPointerDecay` | 46 | 数组到指针 |
| `ToVoid` | 25 |  |
| `UserDefinedConversion` | 16 | 用户自定义转换 |
| `IntegralToFloating` | 4 | 整数转浮点 |
| `FloatingToIntegral` | 1 | 浮点转整数（截断） |

## (source, target, castKind) 三元组 — 前 80

| Source | Target | CastKind | 次数 | 宿主数 |
|---|---|---|---|---|
| `int` | `int` | `LValueToRValue` | 814 | 25 |
| `int` | `size_type` | `IntegralCast` | 465 | 25 |
| `reference (size_type) noexcept` | `reference (*)(size_type) noexcept` | `FunctionToPointerDecay` | 389 | 20 |
| `ZZ` | `ZZ` | `NoOp` | 308 | 29 |
| `basic_polynomial<...>` | `basic_polynomial<...>` | `FunctionToPointerDecay` | 169 | 25 |
| `std::vector<...>` | `std::vector<...>` | `NoOp` | 133 | 25 |
| `iterator` | `__normal_iterator<...>` | `NoOp` | 132 | 20 |
| `uint64_t` | `uint64_t` | `LValueToRValue` | 125 | 35 |
| `__normal_iterator<...>` | `__normal_iterator<...>` | `FunctionToPointerDecay` | 114 | 36 |
| `size_t` | `size_t` | `LValueToRValue` | 104 | 9 |
| `reference () noexcept` | `reference (*)() noexcept` | `FunctionToPointerDecay` | 100 | 35 |
| `bool (__normal_iterator<...>` | `bool (*)(__normal_iterator<...>` | `FunctionToPointerDecay` | 99 | 36 |
| `upolynomial_<...>` | `basic_polynomial<...>` | `NoOp` | 90 | 18 |
| `iterator` | `__gnu_cxx::__normal_iterator<...>` | `NoOp` | 87 | 20 |
| `upolynomial_<...>` | `upolynomial_<...>` | `NoOp` | 84 | 27 |
| `ZZ` | `ZZ` | `ConstructorConversion` | 82 | 22 |
| `int` | `int` | `NoOp` | 71 | 19 |
| `const_iterator` | `__normal_iterator<...>` | `NoOp` | 66 | 20 |
| `const_reference (size_type) noexcept` | `const_reference (*)(size_type) noexcept` | `FunctionToPointerDecay` | 65 | 18 |
| `ZZ (ZZ &, ZZ &)` | `ZZ (*)(ZZ &, ZZ &)` | `FunctionToPointerDecay` | 63 | 18 |
| `int` | `int64_t` | `IntegralCast` | 61 | 28 |
| `typename tuple_element<...>` | `typename tuple_element<...>` | `FunctionToPointerDecay` | 60 | 13 |
| `int64_t` | `int64_t` | `LValueToRValue` | 58 | 14 |
| `<...>` | `typename std::remove_reference<...>` | `BuiltinFnToFnPtr` | 53 | 22 |
| `bool (ZZ &, ZZ &)` | `bool (*)(ZZ &, ZZ &)` | `FunctionToPointerDecay` | 49 | 17 |
| `lex_<...>` | `lex_<...>` | `LValueToRValue` | 48 | 8 |
| `value_type` | `basic_polynomial<...>` | `NoOp` | 45 | 9 |
| `size_type` | `int` | `IntegralCast` | 44 | 19 |
| `Zp` | `Zp` | `NoOp` | 44 | 12 |
| `ZZ &(ZZ &)` | `ZZ &(*)(ZZ &)` | `FunctionToPointerDecay` | 43 | 19 |
| `int64_t (upolynomial_<...>` | `int64_t (*)(upolynomial_<...>` | `FunctionToPointerDecay` | 41 | 12 |
| `QQ` | `QQ` | `NoOp` | 40 | 1 |
| `const_iterator` | `__gnu_cxx::__normal_iterator<...>` | `NoOp` | 39 | 20 |
| `int` | `uint64_t` | `IntegralCast` | 37 | 18 |
| `bool` | `bool` | `LValueToRValue` | 35 | 18 |
| `PolyZp` | `basic_polynomial<...>` | `NoOp` | 35 | 5 |
| `upolynomial_<...>` | `upolynomial_<...>` | `FunctionToPointerDecay` | 32 | 12 |
| `polynomial_<...>` | `polynomial_<...>` | `FunctionToPointerDecay` | 31 | 10 |
| `value_type` | `value_type` | `LValueToRValue` | 31 | 5 |
| `pointer () noexcept` | `pointer (*)() noexcept` | `FunctionToPointerDecay` | 29 | 6 |
| `ZZ &(ZZ &&) noexcept` | `ZZ &(*)(ZZ &&) noexcept` | `FunctionToPointerDecay` | 28 | 14 |
| `void (upolynomial_<...>` | `void (*)(upolynomial_<...>` | `FunctionToPointerDecay` | 28 | 13 |
| `void (std::vector<...>` | `void (*)(std::vector<...>` | `FunctionToPointerDecay` | 27 | 18 |
| `int` | `size_t` | `IntegralCast` | 26 | 7 |
| `bool` | `bool` | `NoOp` | 25 | 21 |
| `value_type` | `QQ` | `NoOp` | 25 | 1 |
| `PolyZp` | `PolyZp` | `NoOp` | 23 | 7 |
| `int` | `void` | `ToVoid` | 23 | 19 |
| `void (char *, char *, unsigned int, char *) __attribute__((noreturn)) noexcept(true)` | `void (*)(char *, char *, unsigned int, char *) __attribute__((noreturn)) noexcept(true)` | `FunctionToPointerDecay` | 23 | 19 |
| `std::pair<...>` | `std::pair<...>` | `NoOp` | 23 | 8 |
| `pair<...>` | `pair<...>` | `FunctionToPointerDecay` | 22 | 12 |
| `basic_polynomial<...>` | `basic_polynomial<...>` | `NoOp` | 22 | 9 |
| `void (ZZ &, ZZ &, ZZ &)` | `void (*)(ZZ &, ZZ &, ZZ &)` | `FunctionToPointerDecay` | 21 | 9 |
| `QQ (QQ &, QQ &)` | `QQ (*)(QQ &, QQ &)` | `FunctionToPointerDecay` | 21 | 1 |
| `Poly` | `basic_polynomial<...>` | `NoOp` | 20 | 5 |
| `PolyZp` | `polynomial_<...>` | `NoOp` | 19 | 5 |
| `ZZ (ZZ &)` | `ZZ (*)(ZZ &)` | `FunctionToPointerDecay` | 19 | 10 |
| `iterator` | `iterator` | `NoOp` | 18 | 6 |
| `std::vector<...>` | `std::vector<...>` | `FunctionToPointerDecay` | 18 | 10 |
| `std::tuple_element<...>` | `std::tuple_element<...>` | `LValueToRValue` | 18 | 4 |
| `PolyZp` | `PolyZp` | `ConstructorConversion` | 18 | 5 |
| `umonomial` | `umonomial` | `ConstructorConversion` | 17 | 12 |
| `Poly` | `polynomial_<...>` | `NoOp` | 17 | 5 |
| `bool (polynomial_<...>` | `bool (*)(polynomial_<...>` | `FunctionToPointerDecay` | 16 | 5 |
| `bool` | `bool` | `UserDefinedConversion` | 16 | 8 |
| `value_type` | `std::vector<...>` | `NoOp` | 16 | 1 |
| `value_type` | `Zp` | `NoOp` | 15 | 4 |
| `bool (std::vector<...>` | `bool (*)(std::vector<...>` | `FunctionToPointerDecay` | 14 | 5 |
| `iterator` | `_Self` | `NoOp` | 13 | 4 |
| `void (polynomial_<...>` | `void (*)(polynomial_<...>` | `FunctionToPointerDecay` | 13 | 7 |
| `polynomial_<...>` | `polynomial_<...>` | `NoOp` | 13 | 6 |
| `Zp (Zp, Zp &)` | `Zp (*)(Zp, Zp &)` | `FunctionToPointerDecay` | 13 | 5 |
| `value_type` | `ZZ` | `NoOp` | 12 | 4 |
| `const_iterator` | `const_iterator` | `ConstructorConversion` | 12 | 7 |
| `Zp &(Zp &&) noexcept` | `Zp &(*)(Zp &&) noexcept` | `FunctionToPointerDecay` | 12 | 3 |
| `bool (_Self &, _Self &) noexcept` | `bool (*)(_Self &, _Self &) noexcept` | `FunctionToPointerDecay` | 11 | 5 |
| `Poly` | `Poly` | `NoOp` | 11 | 6 |
| `unsigned long` | `unsigned long` | `LValueToRValue` | 11 | 5 |
| `void (__gnu_cxx::__normal_iterator<...>` | `void (*)(__gnu_cxx::__normal_iterator<...>` | `FunctionToPointerDecay` | 11 | 7 |
| `std::pair<...>` | `std::pair<...>` | `FunctionToPointerDecay` | 11 | 5 |

## 基础数值类型间的转换（UB 分析重点）

| Source | Target | Kind | 次数 |
|---|---|---|---|
| `int` | `int` | `LValueToRValue` | 814 |
| `int` | `size_type` | `IntegralCast` | 465 |
| `reference (size_type) noexcept` | `reference (*)(size_type) noexcept` | `FunctionToPointerDecay` | 389 |
| `uint64_t` | `uint64_t` | `LValueToRValue` | 125 |
| `size_t` | `size_t` | `LValueToRValue` | 104 |
| `bool (__normal_iterator<...>` | `bool (*)(__normal_iterator<...>` | `FunctionToPointerDecay` | 99 |
| `int` | `int` | `NoOp` | 71 |
| `const_reference (size_type) noexcept` | `const_reference (*)(size_type) noexcept` | `FunctionToPointerDecay` | 65 |
| `int` | `int64_t` | `IntegralCast` | 61 |
| `int64_t` | `int64_t` | `LValueToRValue` | 58 |
| `bool (ZZ &, ZZ &)` | `bool (*)(ZZ &, ZZ &)` | `FunctionToPointerDecay` | 49 |
| `size_type` | `int` | `IntegralCast` | 44 |
| `int64_t (upolynomial_<...>` | `int64_t (*)(upolynomial_<...>` | `FunctionToPointerDecay` | 41 |
| `int` | `uint64_t` | `IntegralCast` | 37 |
| `bool` | `bool` | `LValueToRValue` | 35 |
| `pointer () noexcept` | `pointer (*)() noexcept` | `FunctionToPointerDecay` | 29 |
| `int` | `size_t` | `IntegralCast` | 26 |
| `bool` | `bool` | `NoOp` | 25 |
| `void (char *, char *, unsigned int, char *) __attribute__((noreturn)) noexcept(true)` | `void (*)(char *, char *, unsigned int, char *) __attribute__((noreturn)) noexcept(true)` | `FunctionToPointerDecay` | 23 |
| `bool (polynomial_<...>` | `bool (*)(polynomial_<...>` | `FunctionToPointerDecay` | 16 |
| `bool` | `bool` | `UserDefinedConversion` | 16 |
| `bool (std::vector<...>` | `bool (*)(std::vector<...>` | `FunctionToPointerDecay` | 14 |
| `bool (_Self &, _Self &) noexcept` | `bool (*)(_Self &, _Self &) noexcept` | `FunctionToPointerDecay` | 11 |
| `unsigned long` | `unsigned long` | `LValueToRValue` | 11 |
| `long` | `long` | `LValueToRValue` | 10 |
| `int64_t (polynomial_<...>` | `int64_t (*)(polynomial_<...>` | `FunctionToPointerDecay` | 9 |
| `bool (variable &) const` | `bool (*)(variable &) const` | `FunctionToPointerDecay` | 8 |
| `int64_t` | `int` | `IntegralCast` | 7 |
| `char[11]` | `char` | `ArrayToPointerDecay` | 7 |
| `uint64_t` | `uint64_t` | `NoOp` | 7 |
| `int &(int &, int &)` | `int &(*)(int &, int &)` | `FunctionToPointerDecay` | 6 |
| `int64_t` | `uint64_t` | `IntegralCast` | 6 |
| `int` | `unsigned long` | `IntegralCast` | 6 |
| `reference (size_type)` | `reference (*)(size_type)` | `FunctionToPointerDecay` | 6 |
| `size_t` | `int` | `IntegralCast` | 5 |
| `uint64_t (uint64_t)` | `uint64_t (*)(uint64_t)` | `FunctionToPointerDecay` | 5 |
| `Zp (int64_t, uint64_t)` | `Zp (*)(int64_t, uint64_t)` | `FunctionToPointerDecay` | 4 |
| `double (double) noexcept(true)` | `double (*)(double) noexcept(true)` | `FunctionToPointerDecay` | 4 |
| `bool (QQ &, QQ &)` | `bool (*)(QQ &, QQ &)` | `FunctionToPointerDecay` | 4 |
| `unsigned long long` | `unsigned long long` | `NoOp` | 4 |
| `uint64_t` | `unsigned long long` | `IntegralCast` | 4 |
| `bool (Zp &, int64_t)` | `bool (*)(Zp &, int64_t)` | `FunctionToPointerDecay` | 4 |
| `int` | `std::size_t` | `IntegralCast` | 4 |
| `ZZ (ZZ &&, unsigned long)` | `ZZ (*)(ZZ &&, unsigned long)` | `FunctionToPointerDecay` | 4 |
| `long &(long &, long &)` | `long &(*)(long &, long &)` | `FunctionToPointerDecay` | 3 |
| `char[104]` | `char` | `ArrayToPointerDecay` | 3 |
| `char[147]` | `char` | `ArrayToPointerDecay` | 3 |
| `int` | `double` | `IntegralToFloating` | 3 |
| `double` | `double` | `LValueToRValue` | 3 |
| `Zp (Zp &, int64_t)` | `Zp (*)(Zp &, int64_t)` | `FunctionToPointerDecay` | 3 |
