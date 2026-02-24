# P2b 细化文档：Zp 标量算术优化

> 状态：已实现
> 依据：`docs/design/p2b-zp-arithmetic/architecture.md`

---

## 1. 修改范围

| 文件 | 改动类型 |
|------|----------|
| `clpoly/number.hh` | `Zp` 类完全重写（字段、构造函数、运算符、辅助函数） |
| `clpoly/number.hh` | `addmul<Zp>`、`submul<Zp>`、`set_zero<Zp>` 模板特化修复 |
| `clpoly/polynomial_gcd.hh` | 第 918、1031 行 `number()` 直接赋值修复 |
| `clpoly/polynomial_gcd.cc` | 无需修改（`number()` 只读，隐式拓宽；编译验证） |
| `clpoly/polynomial_factorize_univar.hh` | 无需修改（`static_cast<int64_t>(number())` 安全拓宽） |
| `clpoly/number.hh` `__div<Zp>` 特化 | 无需修改（自动受益于新 `operator/`） |

---

## 2. `clpoly/number.hh`：`Zp` 类

### 2.1 新结构定义

```cpp
class Zp {
private:
    uint32_t _i;     // 值 ∈ [0, p)；类型与语义一致
    uint32_t _p;     // 素数 p（0 表示默认构造的未初始化状态）
    uint64_t _ninv;  // Barrett 常数：UINT64_MAX / p（0 表示未初始化）

    static uint32_t _s_prime;  // 类级别当前素数，默认 0（未设置）；可选调用 set_prime()

    // 预计算 Barrett 常数
    // p==0：未初始化哨兵，返回 0；保持旧 Zp() 行为（_p=0 = "uninitialized"），
    //        addmul/submul 的 |= trick 依赖此哨兵，不需要改动
    // p>=2：UINT64_MAX/p 等价于 floor(2^64/p)（两值差 δ≤1，product<2^64，整数误差=0）
    static uint64_t __barrett_ninv(uint32_t p) {
        if (p == 0) return 0;
        assert(p >= 2 && p < (1u << 31));  // p < 2^31：保证加法 _i+_i 不溢出 uint32_t
        return UINT64_MAX / p;
    }

    // Barrett 归约：将 product（< p²，即 < 2^64）归约到 [0, p)
    uint32_t __barrett_reduce(uint64_t product) const {
        assert(_p != 0);  // 防止 _p=0 的未初始化对象进入乘法
        uint64_t q = (unsigned __int128)product * _ninv >> 64;
        uint64_t r = product - q * _p;
        return (uint32_t)(r >= _p ? r - _p : r);
    }

public:
    // 类级别素数管理（类似 NTL ZZ_p::init）
    static void     set_prime(uint32_t p) { assert(p >= 2 && p < (1u << 31)); _s_prime = p; }
    static uint32_t cur_prime()           { return _s_prime; }

    // ... 见 2.2~2.6
};

inline uint32_t Zp::_s_prime = 0;  // 默认 0（未设置）；与旧 Zp() 的 _p=0 行为一致
```

### 2.2 构造函数

```cpp
// 默认构造：_p=_s_prime（默认 0，即 uninitialized），与旧行为一致
// 若先调用 set_prime(p)，则 Zp() 直接可用；否则需由 addmul/submul |= trick 延迟初始化
Zp() : _i(0), _p(_s_prime), _ninv(__barrett_ninv(_s_prime)) {}

explicit Zp(uint32_t p) : _i(0), _p(p), _ninv(__barrett_ninv(p)) {}
// assert(p>=2) 不在此处：__barrett_ninv 已覆盖（初始化列表先于构造函数体执行）

Zp(uint64_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) {
    _i = (uint32_t)(i % p);
}
Zp(uint32_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) {
    _i = i % p;
}
Zp(int64_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) {
    // 修复原始 bug：i=-p 时 p - (-p)%p = p - 0 = p（越界）
    // 正确做法：C++ 截断式 % 可能给负余数，加 p 再取模
    int64_t r = i % (int64_t)p;
    _i = (uint32_t)(r >= 0 ? r : r + p);
}
Zp(int i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) {
    int64_t r = (int64_t)i % (int64_t)p;
    _i = (uint32_t)(r >= 0 ? r : r + p);
}
Zp(const ZZ& i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) {
    _i = (uint32_t)i.fdiv_ui(p);
}
```

### 2.3 赋值运算符

```cpp
constexpr Zp& operator=(int64_t i) {
    if (_p == 0) { assert(i == 0); _i = 0; return *this; }
    _i = i >= 0 ? (uint32_t)(i % _p) : (uint32_t)(_p - (uint32_t)((-i) % _p));
    return *this;
}
inline Zp& operator=(const ZZ& i) {
    if (_p == 0) { assert(!i); _i = 0; return *this; }
    _i = (uint32_t)i.fdiv_ui(_p);
    return *this;
}
```

### 2.4 访问接口

```cpp
constexpr uint32_t prime()  const { return _p; }
// uint32_t& prime() 删除：_p 和 _ninv 必须同步，不能暴露可变引用
constexpr uint32_t number() const { return _i; }    // 返回类型改为 uint32_t
constexpr uint32_t& number()      { return _i; }    // 返回类型改为 uint32_t&
constexpr void prime(uint32_t p)  { _p = p; _ninv = __barrett_ninv(p); }  // 唯一修改 _p 的入口

explicit constexpr operator uint64_t() const { return _i; }  // 隐式扩展，OK
explicit constexpr operator bool()     const { return _i != 0; }

constexpr void normalization() {
    assert(_p);
    _i %= _p;  // 通常 _i 已在 [0,p)，此处保险用，频率低
}
```

**注**：`prime(uint32_t p)` 在修改 `_p` 时必须同步更新 `_ninv`，这是新增的约束。

### 2.5 算术运算符

```cpp
// 一元负号
Zp operator-() const {
    return Zp(_i == 0 ? 0u : _p - _i, _p);
}

// 加法：a + b < 2p，一次条件减，无除法
friend inline Zp operator+(Zp op1, const Zp& op2) {
    assert(op1._p == op2._p);
    uint32_t r = op1._i + op2._i;
    op1._i = (r >= op1._p) ? r - op1._p : r;
    return op1;
}
inline Zp& operator+=(const Zp& op2) {
    assert(_p == op2._p);
    uint32_t r = _i + op2._i;
    _i = (r >= _p) ? r - _p : r;
    return *this;
}

// 减法：a - b = a + (p - b)，同样一次条件减
friend inline Zp operator-(Zp op1, const Zp& op2) {
    assert(op1._p == op2._p);
    op1._i = op1._i >= op2._i ? op1._i - op2._i : op1._p - op2._i + op1._i;
    return op1;
}
inline Zp& operator-=(const Zp& op2) {
    assert(_p == op2._p);
    _i = _i >= op2._i ? _i - op2._i : _p - op2._i + _i;
    return *this;
}

// 乘法：Barrett 归约
friend inline Zp operator*(Zp op1, const Zp& op2) {
    assert(op1._p == op2._p);
    op1._i = op1.__barrett_reduce((uint64_t)op1._i * op2._i);
    return op1;
}
friend inline Zp operator*(Zp op1, int64_t op2) {
    op1._i = op1.__barrett_reduce((uint64_t)op1._i * (uint64_t)Zp(op2, op1._p)._i);
    return op1;
}
friend inline Zp operator*(int64_t op1, const Zp& op2) { return op2 * op1; }
inline Zp& operator*=(const Zp& op2) {
    assert(_p == op2._p);
    _i = __barrett_reduce((uint64_t)_i * op2._i);
    return *this;
}

// 除法：乘逆元，inv_prime 不变
friend inline Zp operator/(Zp op1, const Zp& op2) {
    assert(op1._p == op2._p);
    op1._i = op1.__barrett_reduce((uint64_t)op1._i * inv_prime(op2._i, op1._p));
    return op1;
}
inline Zp& operator/=(const Zp& op2) {
    assert(_p == op2._p);
    _i = __barrett_reduce((uint64_t)_i * inv_prime(op2._i, _p));
    return *this;
}
```

**inv 方法**（不变，`inv_prime` 依赖 `uint64_t` 参数，仍可用）：
```cpp
inline Zp inv() const {
    Zp r(_p);
    r._i = (uint32_t)inv_prime(_i, _p);
    return r;
}
```

### 2.6 比较运算符（不变，仅类型调整）

`_i` 从 `uint64_t` 改为 `uint32_t` 后，`== != < >` 等比较逻辑不变，参数类型随之调整即可。

---

## 3. `clpoly/number.hh`：模板特化修复

### 3.1 `set_zero<Zp>`

```cpp
// 改前（直接操作 number() 引用，绕过接口）：
template<>
inline void set_zero(Zp& op) {
    op.number() = 0;
}

// 改后：删除特化，走泛型模板即可
// 泛型模板调用 op = 0（即 operator=(int64_t)），结果同样为 _i=0，无需特化
```

### 3.2 `addmul<Zp>`（`op = op + op1 * op2`）

```cpp
// 改前（直接操作 number()，依赖 uint64_t 防溢出）：
template<>
inline void addmul(Zp& op, const Zp& op1, const Zp& op2) {
    assert(...);
    op.prime() |= op1.prime();
    op.number() += op1.number() * op2.number();  // uint64_t 乘，unsafe with uint32_t
    op.number() %= op.prime();
}

// 改后（使用 Zp 运算符，安全且正确）：
template<>
inline void addmul(Zp& op, const Zp& op1, const Zp& op2) {
    assert((op.prime() == op1.prime() || op.prime() == 0)
           && op1.prime() == op2.prime() && op1.prime());
    if (op.prime() == 0) op.prime(op1.prime());  // prime() setter 同步更新 _p 和 _ninv
    op += op1 * op2;
}
```

### 3.3 `submul<Zp>`（`op = op - op1 * op2`）

```cpp
// 改前：
template<>
inline void submul(Zp& op, const Zp& op1, const Zp& op2) {
    assert(...);
    op.prime() |= op1.prime();
    op.number() += op.prime() - (op1.number() * op2.number()) % op.prime();
    op.number() %= op.prime();
}

// 改后：
template<>
inline void submul(Zp& op, const Zp& op1, const Zp& op2) {
    assert((op.prime() == op1.prime() || op.prime() == 0)
           && op1.prime() == op2.prime() && op1.prime());
    if (op.prime() == 0) op.prime(op1.prime());  // prime() setter 同步更新 _p 和 _ninv
    op -= op1 * op2;
}
```

---

## 4. `clpoly/polynomial_gcd.hh`：直接赋值修复

### 4.1 第 918 行

```cpp
// 改前：
v_bool[j_tmp-1] = 0; p_.number() = j_tmp - 1;

// 改后（j_tmp 是 uint64_t，j_tmp-1 < prime < 2^32，显式截断消除 narrowing 警告）：
v_bool[j_tmp-1] = 0; p_.number() = (uint32_t)(j_tmp - 1);
```

### 4.2 第 1031 行

```cpp
// 改前：
p_.number() = 0;

// 改后（消除 signed→unsigned 警告）：
p_.number() = 0u;
```

---

## 5. `clpoly/polynomial_gcd.hh/.cc`：`number()` 作为 ZZ 标量（只读，无需修改）

CRT 合并步骤中，`number()` 仅用于读取值参与 ZZ 运算：

```cpp
// polynomial_gcd.cc:90 — 构造 ZZ 元素
Pout_.push_back({i.first, i.second.number()});   // uint32_t → ZZ（隐式拓宽）

// polynomial_gcd.cc:119,127,132,141,146 和 polynomial_gcd.hh:387,395,400,409,414
// number() 作为 ZZ 乘法中的标量：
ZZ * tmp_inv.number() * Pout_prime   // uint32_t → ZZ * uint32_t
```

`uint32_t` 拓宽到 `int64_t` 再到 `ZZ` 是安全的，语义不变，**无需修改**。
编译时若出现 `ZZ * uint32_t` 重载歧义（ZZ 无 `uint32_t` 重载），添加显式转型 `(int64_t)number()` 即可。**实现时须编译验证**。

---

## 6. `clpoly/polynomial_factorize_univar.hh`：`number()` 读取（无需修改）

```cpp
// 第 261 行
ZZ(static_cast<int64_t>(term.second.number()))
```

`uint32_t` → `int64_t` 是拓宽转换，始终正确，**无需修改**。

---

## 7. 加法溢出约束（隐性前置条件）

新加法实现：
```cpp
uint32_t r = op1._i + op2._i;   // 若 p >= 2^31，两个 (p-1) 相加可能溢出 uint32_t
```

**约束**：`p < 2^31` 时 `(p-1)+(p-1) < 2^32`，加法不溢出 `uint32_t`。
`__barrett_ninv` 中已加 `assert(p < (1u << 31))`，所有构造路径均覆盖。
因式分解所用素数远小于 2^31，此 assert 在正常使用中不会触发。

---

## 8. 调用关系与复用

```
Zp(构造)          → __barrett_ninv(p)            [一次预计算，_p/_ninv 同步]
void prime(p)     → __barrett_ninv(p)            [唯一修改 _p 的入口，_ninv 同步]
operator*         → __barrett_reduce             [每次乘法]
operator/         → __barrett_reduce + inv_prime [inv_prime 不变]
addmul            → prime(p) setter（首次）+ operator+= + operator*
submul            → prime(p) setter（首次）+ operator-= + operator*
```

---

## 9. 测试要求

按实现顺序：

1. **构造函数与基本访问**：各构造路径下 `number()` 和 `prime()` 返回预期值
2. **Barrett 正确性**：对 p ∈ {2, 3, 17, 65537, 999983}，对所有 a, b ∈ [0, p) 验证 `(Zp(a,p) * Zp(b,p)).number() == a*b%p`（小 p 可穷举，大 p 随机抽样）
3. **加减法**：同上，覆盖边界（a=0, b=0, a=p-1, b=p-1）
4. **`addmul` / `submul`**：与改前实现对相同输入结果一致
5. **全量测试**：`bash test/run_all_tests.sh`，278/278 通过
6. **性能验证**：`profile_factorize` 确认 DDF per trial 时间降低 ~3x
