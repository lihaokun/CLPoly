# CLPoly 自定义数值类型设计方案

> 目标：用自定义 ZZ/QQ 类型替换 `mpz_class`/`mpq_class`，引入 small integer optimization，
> 同时解决 macOS 编译问题 (GitHub issue #3)。

---

## 1. 动机

### 1.1 当前问题

CLPoly 当前直接使用 GMP C++ 包装：

```cpp
typedef mpz_class ZZ;
typedef mpq_class QQ;
```

存在以下问题：

- **性能**：`mpz_class` 对每个值都进行堆分配，即使值为 0 或 1
- **跨平台**：GMP C++ 接口不支持 `long long` 构造，macOS 上 `int64_t = long long`，导致
  `int64_t → ZZ` 隐式转换歧义（issue #3）
- **缓存效率**：多项式系数数组中每个 `mpz_class` 都是独立的堆对象，缓存不友好

### 1.2 参考：FLINT 的 fmpz 方案

FLINT 用单个 `slong`（8 字节）实现双模式存储，通过位编码区分小整数和 mpz 指针。
小整数场景下比 `mpz_class` 快约 3 倍，但依赖平台地址空间布局假设，编码复杂。

### 1.3 CLPoly 的方案选择

采用 16 字节双字段方案：`int64_t` 存值 + `mpz_ptr` 做判别（`nullptr` = 小整数）。

| | FLINT 8B 位编码 | CLPoly 16B 双字段 |
|---|---|---|
| 大小 | 8 字节 | 16 字节 |
| 小整数范围 | 62 位 | **完整 64 位** |
| 平台假设 | 依赖指针高位为 0 | **无** |
| 判别方式 | 位移 + 比较 | 指针判空 |
| 小整数运算速度 | 相同 | 相同 |
| 代码复杂度 | 高 | **低** |

对于 CLPoly 的符号计算场景（多项式项数通常 <10000），16B 与 8B 的缓存差异可忽略。
如果未来性能成为瓶颈，可仅改 ZZ 内部实现切换到 8B，公共 API 不变。

---

## 2. 设计方案

### 2.1 新 ZZ 类

> **GMP 隔离策略**：头文件需 `#include <gmp.h>`（私有成员 `mpz_ptr` 需要），
> 但**公共 API 中不暴露任何 GMP 类型**。用户代码不需要直接使用 GMP。

```cpp
class ZZ {
private:
    int64_t _val;     // 小整数时存值；大整数时未使用
    mpz_ptr _mpz;     // nullptr = 小整数模式；非空 = 大整数模式

    // ---- 内部辅助 ----
    bool _is_small() const { return _mpz == nullptr; }
    void _promote();             // 小 → 大：分配 mpz，拷贝 _val
    void _demote_if_small();     // 大 → 小：若值在 int64_t 范围内，释放 mpz

public:
    // ---- 构造/析构 ----
    ZZ();                        // = 0，_mpz = nullptr，无堆分配
    ZZ(int v);
    ZZ(long v);
    ZZ(long long v);             // 解决 issue #3
    ZZ(unsigned long v);
    ZZ(unsigned long long v);
    ZZ(const char* str, int base = 10);  // 大数字符串构造
    ZZ(const ZZ& other);
    ZZ(ZZ&& other) noexcept;
    ~ZZ();

    ZZ& operator=(const ZZ& other);
    ZZ& operator=(ZZ&& other) noexcept;
    ZZ& operator=(long v);
    ZZ& operator=(long long v);

    // ---- 算术 ----
    friend ZZ operator+(const ZZ&, const ZZ&);
    friend ZZ operator-(const ZZ&, const ZZ&);
    friend ZZ operator*(const ZZ&, const ZZ&);
    friend ZZ operator/(const ZZ&, const ZZ&);  // 整除
    friend ZZ operator%(const ZZ&, const ZZ&);
    friend ZZ operator-(const ZZ&);              // 取负

    ZZ& operator+=(const ZZ&);
    ZZ& operator-=(const ZZ&);
    ZZ& operator*=(const ZZ&);
    ZZ& operator/=(const ZZ&);
    ZZ& operator%=(const ZZ&);

    // ---- 比较 ----
    friend bool operator==(const ZZ&, const ZZ&);
    friend bool operator!=(const ZZ&, const ZZ&);
    friend bool operator<(const ZZ&, const ZZ&);
    friend bool operator>(const ZZ&, const ZZ&);
    friend bool operator<=(const ZZ&, const ZZ&);
    friend bool operator>=(const ZZ&, const ZZ&);

    // ---- 状态查询 ----
    explicit operator bool() const;   // 非零判定
    bool is_zero() const;
    bool is_one() const;
    int sign() const;                 // -1, 0, 1

    // ---- 数值特化操作（替代 GMP C API） ----
    void addmul(const ZZ& a, const ZZ& b);   // *this += a * b
    void submul(const ZZ& a, const ZZ& b);   // *this -= a * b
    uint64_t fdiv_ui(uint32_t p) const;       // 取模（给 Zp 用）
    size_t sizeinbase(int base) const;

    // ---- 友元特殊函数 ----
    friend ZZ pow(const ZZ& base, uint64_t exp);
    friend ZZ gcd(const ZZ& a, const ZZ& b);
    friend ZZ lcm(const ZZ& a, const ZZ& b);
    friend ZZ abs(const ZZ& a);

    // ---- IO ----
    friend std::ostream& operator<<(std::ostream&, const ZZ&);
    friend std::istream& operator>>(std::istream&, ZZ&);

    // ---- Hash ----
    std::size_t hash() const;
};

// 全局 hash 支持
inline std::size_t hash_value(const ZZ& v) { return v.hash(); }

namespace std {
    template<> struct hash<ZZ> {
        std::size_t operator()(const ZZ& v) const { return v.hash(); }
    };
}
```

**存储布局**（16 字节 / 个）：

```
小整数模式：[ _val: int64_t (8B) | _mpz: nullptr (8B) ]  → 零堆分配
大整数模式：[ _val: unused  (8B) | _mpz: mpz_ptr (8B) ]  → _mpz 指向堆上的 mpz_t
```

**关键设计原则**：

- 小/大判别：`_mpz == nullptr`，无位编码，无平台假设
- 小整数范围：完整 64 位（`int64_t`）
- 公共接口中**不暴露任何 GMP 类型**
- GMP 仅作为内部实现细节，处理溢出到大整数的情况
- 所有 `int64_t` / `long long` 构造原生支持，彻底解决 issue #3

### 2.2 新 QQ 类

内部用两个 ZZ 表示分子/分母：

```cpp
class QQ {
private:
    ZZ _num;   // 分子
    ZZ _den;   // 分母（始终 > 0）

    void _canonicalize();  // 约分 + 分母正规化

public:
    QQ();                        // = 0/1
    QQ(int v);
    QQ(long v);
    QQ(long long v);
    QQ(const ZZ& v);             // 整数 → 有理数
    QQ(const ZZ& num, const ZZ& den);
    QQ(long num, long den);
    QQ(const QQ&);
    QQ(QQ&&) noexcept;
    ~QQ();

    QQ& operator=(const QQ&);
    QQ& operator=(QQ&&) noexcept;
    QQ& operator=(long v);

    // ---- 访问 ----
    const ZZ& num() const;
    const ZZ& den() const;

    // ---- 算术 ----
    friend QQ operator+(const QQ&, const QQ&);
    friend QQ operator-(const QQ&, const QQ&);
    friend QQ operator*(const QQ&, const QQ&);
    friend QQ operator/(const QQ&, const QQ&);
    friend QQ operator-(const QQ&);

    QQ& operator+=(const QQ&);
    QQ& operator-=(const QQ&);
    QQ& operator*=(const QQ&);
    QQ& operator/=(const QQ&);

    // ---- 比较 ----
    friend bool operator==(const QQ&, const QQ&);
    friend bool operator!=(const QQ&, const QQ&);
    friend bool operator<(const QQ&, const QQ&);
    friend bool operator>(const QQ&, const QQ&);
    friend bool operator<=(const QQ&, const QQ&);
    friend bool operator>=(const QQ&, const QQ&);

    // ---- 状态 ----
    explicit operator bool() const;
    bool is_zero() const;
    int sign() const;

    // ---- IO / Hash ----
    friend std::ostream& operator<<(std::ostream&, const QQ&);
    std::size_t hash() const;
};
```

**优势**：小分数（如 1/2, 3/4）的分子和分母都是小整数，完全零堆分配。

### 2.3 Zp 类适配

Zp 当前有**两处** GMP 依赖：

```cpp
// 构造函数（number.hh:132）
// 旧：Zp(ZZ i, uint32_t p) : _i(mpz_fdiv_ui(i.get_mpz_t(), p)) {}
// 新：Zp(ZZ i, uint32_t p) : _i(i.fdiv_ui(p)) {}

// 赋值运算符（number.hh:149）
// 旧：this->_i = mpz_fdiv_ui(i.get_mpz_t(), this->_p);
// 新：this->_i = i.fdiv_ui(this->_p);
```

改动极小。

---

## 3. 影响范围

### 3.1 当前 GMP 直接访问清单

GMP 内部结构访问**高度集中**于 `number.hh` 和 `interval.hh`：

#### `number.hh`（需完全重写）

| 行号 | 代码 | 访问类型 |
|------|------|----------|
| 9 | `#include <gmpxx.h>` | 头文件依赖 |
| 26 | `p.get_mpz_t()->_mp_d` | 深层内部结构 — 直接读取 limb 数组 |
| 26-27 | `p.get_mpz_t()->_mp_size` | 深层内部结构 — 直接读取 size 字段 |
| 38-39 | `typedef mpz_class ZZ / mpq_class QQ` | 类型别名 |
| 44 | `mpz_addmul(op.get_mpz_t(), ...)` | GMP C API |
| 49 | `mpz_submul(op.get_mpz_t(), ...)` | GMP C API |
| 54 | `mpz_fdiv_q(op.get_mpz_t(), ...)` | GMP C API |
| 59 | `mpz_fdiv_qr(op.get_mpz_t(), ...)` | GMP C API |
| 64, 69 | `mpq_class(op1, op2)` | GMP C++ 构造 |
| 75 | `mpz_set_si(op.get_mpz_t(), 0)` | GMP C API |
| 95 | `mpz_pow_ui(x.get_mpz_t(), ...)` | GMP C API |
| 100-101 | `mpz_pow_ui(x.get_num_mpz_t(), ...)` | GMP C API (QQ 内部) |
| 106 | `mpz_sizeinbase(x.get_mpz_t(), i)` | GMP C API |
| 132 | `mpz_fdiv_ui(i.get_mpz_t(), p)` | GMP C API (Zp 构造) |
| 149 | `mpz_fdiv_ui(i.get_mpz_t(), this->_p)` | GMP C API (Zp 赋值) |

#### `polynomial_convert.hh`（4 处 — QQ 访问器改名）

| 行号 | 代码 | 迁移方式 |
|------|------|----------|
| 125 | `i.second.get_den()` | → `.den()` |
| 131 | `i.second.get_num()` / `.get_den()` | → `.num()` / `.den()` |
| 143 | `i.second.get_den()` | → `.den()` |
| 149 | `i.second.get_num()` / `.get_den()` | → `.num()` / `.den()` |

#### `realroot.cc`（4 处 — QQ 访问器改名）

| 行号 | 代码 | 迁移方式 |
|------|------|----------|
| 226 | `a.get_den()`, `b.get_den()` | → `.den()` |
| 227 | `a.get_den()`, `b.get_den()` | → `.den()` |
| 229 | `a.get_num()`, `b.get_den()` | → `.num()` / `.den()` |
| 230 | `b.get_num()`, `a.get_den()` | → `.num()` / `.den()` |

#### `upolynomial.cc`（2 处 — QQ 访问器改名）

| 行号 | 代码 | 迁移方式 |
|------|------|----------|
| 29 | `i.second.get_den()` | → `.den()` |
| 33 | `i.second.get_num()` / `.get_den()` | → `.num()` / `.den()` |

#### `interval.hh`（~30 处 — 暂缓）

深度耦合 `mpq_t`/`mpq_srcptr`/`mpria_*` C API，暂缓处理。

### 3.2 需要修改的文件

| 文件 | 改动量 | 说明 |
|------|--------|------|
| `number.hh` | **重写** | typedef → class 定义；template 特化改为 ZZ 成员方法或新特化 |
| `polynomial_.hh` | **简化** | int64_t→ZZ 的 operator 重载自动工作，可能可以删减 |
| `polynomial_convert.hh` | 少量 | `.get_num()` → `.num()`，`.get_den()` → `.den()`（4 处） |
| `realroot.cc` | 少量 | `.get_num()` → `.num()`，`.get_den()` → `.den()`（4 处） |
| `upolynomial.cc` | 少量 | `.get_num()` → `.num()`，`.get_den()` → `.den()`（2 处） |
| `polynomial_gcd.hh` | 少量 | ZZ 作为模板参数透明传递 |
| `resultant.hh` | 少量 | 同上 |
| `realroot.hh` | 中等 | QQ 作为区间端点 |
| `upolynomial.hh` | 少量 | 类型别名 + eval |
| `interval.hh/cc` | 暂缓 | 深度耦合 mpria（当前不可用），待后续单独处理 |
| 测试文件 | 几乎无 | `ZZ(42)` 等构造自动兼容 |

### 3.3 不需要改的

- `basic_monomial.hh` — 不涉及系数类型
- `basic_polynomial.hh` — 纯模板，系数类型透明
- `variable.hh` — 不涉及
- `associatedgraph.hh` — 不涉及
- `charset.hh` — 不涉及（通过模板参数使用）

### 3.4 template 特化迁移

旧代码中的 template 特化：

```cpp
// 旧
template<> void addmul(mpz_class& op, const mpz_class& a, const mpz_class& b) {
    mpz_addmul(op.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}
template<> struct zore_check<mpz_class> { ... };
```

迁移为 ZZ 的成员或新特化：

```cpp
// 新
template<> void addmul(ZZ& op, const ZZ& a, const ZZ& b) {
    op.addmul(a, b);  // ZZ 内部处理 small/large 两种路径
}
template<> struct zore_check<ZZ> {
    bool operator()(const ZZ& op) { return op.is_zero(); }
};
```

---

## 4. 实施计划

```
Phase 1: 实现 ZZ 类
         ├── int64_t + mpz_ptr 双字段存储（16 字节）
         ├── promote（小→大）/ demote_if_small（大→小）
         ├── 全套运算符和比较（含 %=）
         ├── addmul / submul / pow / gcd / lcm / fdiv_ui / sizeinbase
         ├── hash（小整数: 直接哈希 _val；大整数: 哈希 mpz limb 数组）
         └── 单元测试

Phase 2: 替换 number.hh
         ├── 删除 typedef mpz_class ZZ
         ├── 更新所有 template 特化（addmul, submul, __div, set_zero, zore_check）
         ├── 更新 Zp 构造函数和赋值运算符（2 处）
         └── 编译测试

Phase 3: 实现 QQ 类
         ├── 两个 ZZ 存储 + _canonicalize() 自动约分
         ├── num() / den() 接口
         ├── 完整比较运算符（含 <=, >=）
         ├── pow(QQ, uint64_t)
         └── 单元测试

Phase 4: 替换 QQ 相关代码
         ├── polynomial_convert.hh: .get_num() → .num(), .get_den() → .den()（4 处）
         ├── realroot.cc: 同上（4 处）
         ├── upolynomial.cc: 同上（2 处）
         └── 编译测试

Phase 5: 全量回归测试
```

Phase 1-2 完成后即可解决 issue #3 并获得 ZZ 性能提升。

---

## 5. 性能预期

| 场景 | 当前 (mpz_class) | 预期 (自定义 ZZ) |
|------|-----------------|-----------------|
| 小整数构造 | 堆分配 | 零分配，直接赋值 |
| 小整数加法 | GMP 函数调用 | 一条 ADD 指令 + 溢出检查 |
| 多项式系数存储 | 每系数 ~24 字节 + 堆 | 小系数 16 字节，零堆分配 |
| 缓存命中率 | 低（指针追踪） | 高（连续内存） |
| macOS 编译 | 失败 (issue #3) | 正常 |

---

## 6. 参考

- FLINT fmpz 实现: https://flintlib.org/doc/fmpz.html
- GMP C++ 接口限制: https://gmplib.org/manual/C_002b_002b-Interface-General
- GitHub issue #3: int64_t → mpz_class 歧义
