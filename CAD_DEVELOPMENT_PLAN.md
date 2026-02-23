# CAD (Cylindrical Algebraic Decomposition) 开发计划
## 项目概述

实现柱形代数分解（CAD）算法，用于多项式系统的符号计算和求解。本项目基于现有的 CLPoly 库，已具备多项式基础运算能力，将分阶段实现 CAD 的核心组件。

## 开发阶段

### 第一阶段：基础实现 

#### 投影函数


**目标**：实现 CAD 算法的投影阶段核心函数

##### 函数签名
```cpp
enum class projection_method {
    MCCALLUM,    // McCallum 投影算子（默认）
    LAZARD      // Lazard 投影算子
};

std::vector<clpoly::polynomial_ZZ> project(
    const std::vector<clpoly::polynomial_ZZ>& polys, 
    const clpoly::variable& x, 
    projection_method method = projection_method::MCCALLUM
);
```

##### 功能描述
- **输入参数**：
  - `polys`: 多项式集合 `{f₁, f₂, ..., fₙ}`
  - `x`: 需要消去的投影变量
  - `method`: 投影算子类型

- **输出**：投影后的多项式集合（不含变量 `x`）

####  提升

```cpp
std::vector<clpoly::QQ> sample_open_intervals(
    const std::vector<clpoly::upolynomial_ZZ>& polys
);
```
#### 功能描述

**sample_open_intervals** - 开区间采样点提升

- **输入参数**：
  - `polys`: 一元多项式集合（关于提升变量的多项式，系数为有理数）
  
- **输出**：
  - 开区间内的采样点（有理数），保证每个polys隔离出来的开区间各一个采样点


## 参考资料