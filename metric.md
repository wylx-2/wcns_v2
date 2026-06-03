# SCMM 曲线坐标度量系数计算模块

---

# 一、总体任务

实现三维曲线坐标网格下的：

1. SCMM 度量系数计算
2. Jacobian 计算
3. 单元中心 → 半节点（面）度量插值
4. 周期/并行 connectivity 处理
5. 通用高阶插值—差分算子封装

注意：

所有导数计算必须统一采用：

$$
\text{Interpolation} \rightarrow \text{Differentiation}
$$

即：

* 先插值到半节点
* 再利用半节点值求导

禁止直接在中心点上差分替代该流程。

---

# 二、输入数据

已知三维单元中心坐标：

$$
x(i,j,k),\quad y(i,j,k),\quad z(i,j,k)
$$

数组大小：

$$
(nx+5)\times(ny+5)\times(nz+5)
$$

其中包含 ghost layers。

输入包括：

* 网格坐标
* 网格尺寸

$$
\Delta \xi,\Delta \eta,\Delta \zeta
$$

* 周期/并行 connectivity 信息

---

# 三、SCMM 度量计算

采用 SCMM 方法。

度量系数定义如下：

## 1. ξ 方向度量

$$
{{{\hat{\xi }}}_{x}}=\frac{1}{2}\left[ {{\left( z{{y}_{\eta }} \right)}_{\zeta }}+{{\left( y{{z}_{\zeta }} \right)}_{\eta }}-{{\left( z{{y}_{\zeta }} \right)}_{\eta }}-{{\left( y{{z}_{\eta }} \right)}_{\zeta }} \right] 
$$

$$
{{{\hat{\xi }}}_{y}}=\frac{1}{2}\left[ {{\left( x{{z}_{\eta }} \right)}_{\zeta }}+{{\left( z{{x}_{\zeta }} \right)}_{\eta }}-{{\left( x{{z}_{\zeta }} \right)}_{\eta }}-{{\left( z{{x}_{\eta }} \right)}_{\zeta }} \right] 
$$

$$
{{{\hat{\xi }}}_{z}}=\frac{1}{2}\left[ {{\left( y{{x}_{\eta }} \right)}_{\zeta }}+{{\left( x{{y}_{\zeta }} \right)}_{\eta }}-{{\left( y{{x}_{\zeta }} \right)}_{\eta }}-{{\left( x{{y}_{\eta }} \right)}_{\zeta }} \right] 
$$

## 2. η 方向度量

$$
{{{\hat{\eta }}}_{x}}=\frac{1}{2}\left[ {{\left( z{{y}_{\zeta }} \right)}_{\xi }}+{{\left( y{{z}_{\xi }} \right)}_{\zeta }}-{{\left( z{{y}_{\xi }} \right)}_{\zeta }}-{{\left( y{{z}_{\zeta }} \right)}_{\xi }} \right] 
$$

$$
{{{\hat{\eta }}}_{y}}=\frac{1}{2}\left[ {{\left( x{{z}_{\zeta }} \right)}_{\xi }}+{{\left( z{{x}_{\xi }} \right)}_{\zeta }}-{{\left( x{{z}_{\xi }} \right)}_{\zeta }}-{{\left( z{{x}_{\zeta }} \right)}_{\xi }} \right] 
$$

$$
{{{\hat{\eta }}}_{z}}=\frac{1}{2}\left[ {{\left( y{{x}_{\zeta }} \right)}_{\xi }}+{{\left( x{{y}_{\xi }} \right)}_{\zeta }}-{{\left( y{{x}_{\xi }} \right)}_{\zeta }}-{{\left( x{{y}_{\zeta }} \right)}_{\xi }} \right] 
$$

---

## 3. ζ 方向度量

$$
 {{{\hat{\zeta }}}_{x}}=\frac{1}{2}\left[ {{\left( z{{y}_{\xi }} \right)}_{\eta }}+{{\left( y{{z}_{\eta }} \right)}_{\xi }}-{{\left( z{{y}_{\eta }} \right)}_{\xi }}-{{\left( y{{z}_{\xi }} \right)}_{\eta }} \right] 
$$

$$
 {{{\hat{\zeta }}}_{y}}=\frac{1}{2}\left[ {{\left( x{{z}_{\xi }} \right)}_{\eta }}+{{\left( z{{x}_{\eta }} \right)}_{\xi }}-{{\left( x{{z}_{\eta }} \right)}_{\xi }}-{{\left( z{{x}_{\xi }} \right)}_{\eta }} \right] 
$$

$$
{{{\hat{\zeta }}}_{z}}=\frac{1}{2}\left[ {{\left( y{{x}_{\xi }} \right)}_{\eta }}+{{\left( x{{y}_{\eta }} \right)}_{\xi }}-{{\left( y{{x}_{\eta }} \right)}_{\xi }}-{{\left( x{{y}_{\xi }} \right)}_{\eta }} \right] 
$$

---

# 四、Jacobian 计算

Jacobian 按 SCMM 形式计算：

$$
 \frac{1}{J}=\frac{1}{3}\left[ {{\left( x{{{\hat{\xi }}}_{x}}+y{{{\hat{\xi }}}_{y}}+z{{{\hat{\xi }}}_{z}} \right)}_{\xi }}+{{\left( x{{{\hat{\eta }}}_{x}}+y{{{\hat{\eta }}}_{y}}+z{{{\hat{\eta }}}_{z}} \right)}_{\eta }}+{{\left( x{{{\hat{\zeta }}}_{x}}+y{{{\hat{\zeta }}}_{y}}+z{{{\hat{\zeta }}}_{z}} \right)}_{\zeta }} \right] 
$$

即必须按照以下工作流：

例如先计算
$$
y_\eta
$$
再计算
$$
(z y_\eta)_\zeta
$$
最终组装：
$$
\hat{\xi}_x
$$
计算 Jacobian 时继续求：
$$
{{\left( x{{{\hat{\xi }}}_{x}}+y{{{\hat{\xi }}}_{y}}+z{{{\hat{\xi }}}_{z}} \right)}_{\xi }}
$$
必须严格保持此计算顺序。

---

# 五、统一高阶插值—差分算法

所有导数统一采用：

## Step 1：中心 → 半节点插值

例如 ζ 方向：

$$
{{a}_{i,j,k+1/2}}=\frac{75}{128}\left( {{a}_{i,j,k}}+{{a}_{i,j,k+1}} \right)-\frac{25}{256}\left( {{a}_{i,j,k-1}}+{{a}_{i,j,k+2}} \right)+\frac{3}{256}\left( {{a}_{i,j,k-2}}+{{a}_{i,j,k+3}} \right)
$$

---

## Step 2：半节点差分

ζ 导数：

$$
{{a}_{\zeta (i,j,k)}}=\frac{75}{64\Delta \zeta }\left( {{a}_{i,j,k+1/2}}-{{a}_{i,j,k-1/2}} \right)-\frac{25}{384\Delta \zeta }\left( {{a}_{i,j,k+3/2}}-{{a}_{i,j,k-3/2}} \right)+\frac{3}{640\Delta \zeta }\left( {{a}_{i,j,k+5/2}}-{{a}_{i,j,k-5/2}} \right)
$$

---

# 六、非周期边界处理

对于非周期边界，采用单边插值 + 单边差分。

以左边界为例。

---

## 1. 单边插值

### 第一半节点

$$
{{a}_{1/2}}=\frac{1}{128}\left( 315{{a}_{1}}-420{{a}_{2}}+378{{a}_{3}}-180{{a}_{4}}+35{{a}_{5}} \right) 
$$

### 第二半节点

$$
{{a}_{3/2}}=\frac{1}{128}\left( 35{{a}_{1}}+140{{a}_{2}}-70{{a}_{3}}+28{{a}_{4}}-5{{a}_{5}} \right)
$$

---

## 2. 单边差分

第一点：

$$
{{a}_{\xi ,1}}=\frac{1}{24\Delta \xi }\left( -22{{a}_{1/2}}+17{{a}_{3/2}}+9{{a}_{5/2}}-5{{a}_{7/2}}+{{a}_{9/2}} \right) 
$$

第二点：

$$
{{a}_{\xi ,2}}=\frac{1}{24\Delta \xi }\left( {{a}_{1/2}}-27{{a}_{3/2}}+27{{a}_{5/2}}-{{a}_{7/2}} \right)
$$

要求：

* 单边处理独立实现
* 不使用简单 ghost extrapolation 替代
* 自动根据边界位置切换 stencil

---

# 七、周期边界与并行 Connectivity

对于：

* periodic boundary
* MPI connectivity
* block connectivity

注意：

高阶 stencil 可能需要 **多次 halo 数据交换**。不假设一次通信足够。同时周期边界交换中间量时可能注意叠加坐标平移，如交换$(xy_\eta)$等数据时可能涉及$x$的再次平移。


---

# 八、度量存储结构

## 1. 中心度量

所有 SCMM 度量首先存储于单元中心点处，大小为

$$
(nx+5)(ny+5)(nz+5)
$$

包括：

* $\hat{\xi}_x,\hat{\xi}_y,\hat{\xi}_z$
* $\hat{\eta}_x,\hat{\eta}_y,\hat{\eta}_z$
* $\hat{\zeta}_x,\hat{\zeta}_y,\hat{\zeta}_z$
* $J$

储存在相应的类中，如grid

---

## 2. 面度量

随后将单元中心处的度量系数插值至半节点，插值采用前面相同方法，包括六点中心插值和边界单边插值

于是X面数组的大小为

$$
(nx+6)(ny+5)(nz+5)
$$

Y/Z 类似。储存在相应的类中，如grid

---

# 九、程序架构要求

由于插值—差分过程与单边处理会反复调用，必须写为通用功能模块。包括

* 中心插值
* 单边插值
* 中心差分
* 单边差分
* 周期交换或者并行通信

这一部分按需求新建类或函数实现
