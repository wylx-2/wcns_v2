# WCNS_v2 —— 三维可压缩 Navier-Stokes 方程并行求解器

基于 WCNS 格式的三维可压缩 NS 方程 CFD 求解器，采用多块结构网格与 MPI 并行计算。

## 构建与运行

```bash
# 配置（含 MPI）
cmake -S . -B build -DENABLE_MPI=ON
cmake --build build

# 运行（N 为 MPI 进程数）
mpirun -np N ./build/wcns_v2 input.cfg
```

- C++17，编译选项 `-Wall -Wextra -Wpedantic`

---

## 一、文件结构

```
wcns_v2/
├── CMakeLists.txt                # 顶层 CMake 构建文件
├── README.md                     # 项目说明与设计文档
├── CLAUDE.md                     # AI 辅助开发指引
├── cmake/                        # CMake 模块
│   ├── FindMPI.cmake
│   └── CompilerFlags.cmake
├── doc/                          # 详细文档
│   ├── theory.md                 # 控制方程与数值格式理论
│   └── user_guide.md             # 用户手册
├── input/                        # 算例输入文件
│   └── example.cfg               # 示例配置文件
├── src/
│   ├── main.cpp                  # 入口，解析命令行并启动求解器
│   │
│   ├── core/                     # 核心框架
│   │   ├── solver.h / .cpp       # 求解器总控，驱动整个计算流程
│   │   ├── config.h / .cpp       # 配置文件解析与参数管理
│   │   └── timer.h / .cpp        # 计时与性能统计
│   │
│   ├── grid/                     # 网格模块
│   │   ├── grid.h / .cpp         # 多块网格管理（Grid 类）
│   │   ├── block.h / .cpp        # 单块结构网格（Block 类，含度量系数）
│   │   ├── metric.h / .cpp       # 网格导数和 Jacobian 计算
│   │   └── grid_reader.h / .cpp  # 网格文件读取（Plot3D / CGNS / 自定义格式）
│   │
│   ├── field/                    # 流场数据模块
│   │   ├── field.h / .cpp        # 流场变量存储（保守/原始变量）
│   │   └── thermodynamics.h / .cpp # 热力学关系（状态方程、声速等）
│   │
│   ├── scheme/                   # 数值格式模块
│   │   ├── wcns.h / .cpp         # WCNS 非线性插值核心（WCNS5/7, WCNS-Z）
│   │   ├── wcns_coeff.h          # WCNS 线性/非线性权重预计算系数
│   │   ├── compact_diff.h / .cpp # 紧凑/显式差分求导（三对角/五对角）
│   │   ├── inviscid.h / .cpp     # 无粘通量：WCNS 插值 + Riemann 求解器
│   │   ├── riemann.h / .cpp      # Riemann 求解器（Roe / HLLC / AUSM+）
│   │   ├── viscous.h / .cpp      # 粘性通量：中心差分
│   │   └── characteristic.h / .cpp # 特征空间投影与反投影
│   │
│   ├── bc/                       # 边界条件模块
│   │   ├── bc_manager.h / .cpp   # 边界条件管理器
│   │   ├── bc.h                  # 边界条件基类
│   │   ├── farfield.h / .cpp     # 远场 / 特征边界条件
│   │   ├── wall.h / .cpp         # 物面边界（无滑移绝热 / 等温）
│   │   ├── symmetry.h / .cpp     # 对称面边界
│   │   ├── inflow.h / .cpp       # 亚/超声速入口
│   │   ├── outflow.h / .cpp      # 亚/超声速出口
│   │   └── periodic.h / .cpp     # 周期边界条件
│   │
│   ├── time/                     # 时间推进模块
│   │   ├── time_integrator.h     # 时间推进基类
│   │   ├── explicit_rk.h / .cpp  # 显式 Runge-Kutta（RK3-TVD / RK4 / RK5）
│   │   └── implicit_lusgs.h / .cpp # 隐式 LU-SGS
│   │
│   ├── parallel/                 # 并行模块
│   │   ├── mpi_env.h / .cpp      # MPI 环境初始化与销毁
│   │   ├── domain_decomp.h / .cpp # 多块区域分解（负载均衡）
│   │   ├── halo_exchange.h / .cpp # 虚拟网格数据交换
│   │   └── parallel_compact.h / .cpp # 并行三对角求解器（Pipelined Thomas）
│   │
│   ├── io/                       # 输入输出模块
│   │   ├── cfg_reader.h / .cpp   # 配置文件解析器
│   │   ├── solution_io.h / .cpp  # 解文件输出（Tecplot / VTK / Plot3D）
│   │   └── restart.h / .cpp      # 断点续算文件读写
│   │
│   └── utils/                    # 通用工具
│       ├── types.h               # 全局类型定义（Real, Int, Vec3, 数组别名）
│       ├── constants.h           # 物理/数学常量（气体常数、比热比等）
│       ├── math_utils.h / .cpp   # 数学辅助函数
│       └── logger.h / .cpp       # 日志系统（控制台 + 文件输出）
│
└── test/                         # 单元测试
    ├── test_wcns.cpp             # WCNS 插值精度测试
    └── test_riemann.cpp          # Riemann 求解器测试
```

---

## 二、程序结构（类层次与模块关系）

### 2.1 顶层驱动

```
main.cpp
  │
  └─ Solver (求解器总控)
       ├── Config           (参数配置)
       ├── Grid             (多块网格)
       │     └── Block[0..M]  (各网格块)
       │           ├── Metric (度量系数)
       │           └── Field  (该块的流场)
       ├── WCNS             (WCNS 插值方案)
       ├── CompactDiff      (差分求导)
       ├── InviscidFlux     (无粘通量)
       │     ├── WCNS (重构)
       │     ├── Characteristic (特征投影)
       │     └── Riemann
       ├── ViscousFlux      (粘性通量)
       ├── TimeIntegrator   (时间推进)
       ├── BCManager        (边界条件)
       ├── Parallel         (MPI 并行)
       │     ├── DomainDecomp
       │     ├── HaloExchange
       │     └── ParallelCompact
       └── SolutionIO       (输入输出)
```

### 2.2 核心类职责

| 类 | 职责 |
|---|---|
| `Solver` | 总控，持有所有子模块，驱动初始化→时间循环→输出→终止的完整流程 |
| `Config` | 解析 `input.cfg`，提供键值对访问，管理物理模型、数值格式、时间步、输出策略等参数 |
| `Grid` | 管理全部网格块，提供网格块查询、全局→局部索引映射 |
| `Block` | 单个三维结构网格块：维度 `(ni,nj,nk)`，含坐标 `(x,y,z)`、度量系数、虚拟网格层数 `ng` |
| `Field` | 流场数据：保守变量 `Q`、原始变量 `V`、通量 `RHS`、特征变量，支持内存分配与拷贝 |
| `WCNS` | WCNS 非线性插值核心：生成候选模板、计算光滑因子和权重，输出半节点处插值 |
| `CompactDiff` | 紧凑/显式差分求导，支持三对角/五对角系统并行求解 |
| `InviscidFlux` | 组合 WCNS 重构与 Riemann 求解器，计算无粘通量 |
| `Riemann` | Roe / HLLC / AUSM+ 等迎风型 Riemann 近似求解器 |
| `ViscousFlux` | 中心差分计算粘性通量，含 Sutherland 粘性律 |
| `BCManager` | 遍历所有边界块，按面/类型施加虚拟网格值 |
| `TimeIntegrator` | 虚基类，派生 `ExplicitRK`（RK3-TVD/RK4）和 `ImplicitLUSGS` |
| `DomainDecomp` | 按块/网格线分块分配 MPI 进程，实现负载均衡 |
| `HaloExchange` | 相邻块间交换虚拟网格层（SendRecv），支持面/棱/角区域 |

### 2.3 多态/策略模式

- **数值通量**：`FluxScheme` 基类 → `RoeFlux`, `HLLCFlux`, `AUSMPlusFlux`
- **时间推进**：`TimeIntegrator` 基类 → `ExplicitRK`, `ImplicitLUSGS`
- **边界条件**：`BC` 基类 → `FarfieldBC`, `WallBC`, `SymmetryBC`, `InflowBC`, `OutflowBC`

---

## 三、数据结构

### 3.1 全局类型（`utils/types.h`）

```cpp
using Real    = double;                    // 浮点精度（可切换 float）
using Int     = int;                       // 整数索引
using Vec3    = std::array<Real, 3>;       // 三维向量

// 三维数组：使用连续内存分配，方便 MPI 打包/解包
// MultiArray3D<Real> 封装 std::vector<Real> + (ni,nj,nk) 维度，支持 (i,j,k) 索引
template<typename T>
using MultiArray3D = /* 连续存储的 3D 数组，手动索引: idx = k + nk*(j + nj*i) */;
```

### 3.2 网格块（`grid/block.h`）

```cpp
struct Block {
    Int id;            // 全局块编号
    Int ni, nj, nk;    // 各方向网格节点数
    Int ng;            // 虚拟网格层数（WCNS-5 需 ng=3, WCNS-7 需 ng=4）

    // 网格坐标（含虚拟网格）
    MultiArray3D<Real> x, y, z;     // 形状: (ni+2*ng) × (nj+2*ng) × (nk+2*ng)

    // 网格度量系数（Jacobi 逆矩阵分量，各方向）
    MultiArray3D<Real> xi_x, xi_y, xi_z;   // ∂ξ/∂x, ∂ξ/∂y, ∂ξ/∂z
    MultiArray3D<Real> eta_x, eta_y, eta_z;
    MultiArray3D<Real> zeta_x, zeta_y, zeta_z;
    MultiArray3D<Real> jac;                 // Jacobian |∂(ξ,η,ζ)/∂(x,y,z)| = 1/|J|

    // 相邻块信息（用于 MPI 通信）
    std::array<Int, 6> neighbor_id;   // -1 表示物理边界（IMIN, IMAX, JMIN, JMAX, KMIN, KMAX）
    std::array<Int, 6> neighbor_rank; // 相邻块所在进程号
};
```

### 3.3 流场变量（`field/field.h`）

```cpp
struct Field {
    Int ni, nj, nk, ng;

    // 保守变量 Q[5]: [ρ, ρu, ρv, ρw, ρE]
    // 形状: (ni+2*ng) × (nj+2*ng) × (nk+2*ng)
    MultiArray3D<Real> rho, rhou, rhov, rhow, rhoE;

    // 原始变量（按需计算和同步）
    MultiArray3D<Real> u, v, w, p, T;

    // 右端项（残差 RHS），与 Q 同形状
    MultiArray3D<Real> rhs_rho, rhs_rhou, rhs_rhov, rhs_rhow, rhs_rhoE;

    // 当地时间步长
    MultiArray3D<Real> dt_local;

    // 守恒 → 原始变量转换
    void conservative_to_primitive();
    // 原始 → 守恒变量转换
    void primitive_to_conservative();
};
```

### 3.4 WCNS 插值数据（`scheme/wcns.h`）

```cpp
// WCNS 参数集（编译期/运行期常量）
struct WCNSParams {
    Int order;           // 精度阶数：5 或 7
    Real eps;            // 小量（防除零），通常 1e-6 ~ 1e-40
    Real epsilon;        // WCNS-Z 的 ε 参数（p=2 推荐）
    Int weno_type;       // JS / Z / CU 型权重
    bool characteristic; // 是否使用特征空间插值
};

// 模板线性权重 d_k（5阶：k=0,1,2 → d=[1/16, 10/16, 5/16] 左偏）
// 候选模板插值系数 c_{k,r}（k=模板索引, r=模板内点偏移）

// 运行时：在块内每条网格线上依次调用
//   wcns_interpolate_xi(dir, field, block) → 输出左右状态 QL, QR 在半节点处
```

### 3.5 半节点/面通量存储

```cpp
// 块内面通量（用于紧凑差分或积分）
struct FaceFlux {
    // F_{i+1/2,j,k} 等，存储在 半整数索引 对应的面上
    // 由于结构网格，面数 ≈ 体单元数，使用 (ni+1) × nj × nk 存储
    MultiArray3D<Vec3> flux_i;  // i 方向面（F̂ 或 F̂-F̂v）
    MultiArray3D<Vec3> flux_j;  // j 方向面（Ĝ 或 Ĝ-Ĝv）
    MultiArray3D<Vec3> flux_k;  // k 方向面（Ĥ 或 Ĥ-Ĥv）
};
```

### 3.6 MPI 通信缓冲区

```cpp
struct HaloBuffer {
    // 六个面各需要发送 ng 层数据，包含 5 个保守变量
    // 按面打包为连续缓冲区（IMIN: ng×nj×nk×5, IMAX: ng×nj×nk×5, ...）
    std::array<std::vector<Real>, 6> send_buf;
    std::array<std::vector<Real>, 6> recv_buf;
    std::array<MPI_Request, 12>   requests;   // 6发 + 6收 (非阻塞)
};
```

---

## 四、算法过程

### 4.1 控制方程

曲线坐标系(ξ,η,ζ)下的无量纲化三维可压缩 Navier-Stokes 方程：

$$\frac{\partial \hat{\mathbf{Q}}}{\partial t} + \frac{\partial \hat{\mathbf{F}}}{\partial \xi} + \frac{\partial \hat{\mathbf{G}}}{\partial \eta} + \frac{\partial \hat{\mathbf{H}}}{\partial \zeta}
= \frac{\partial \hat{\mathbf{F}}_v}{\partial \xi} + \frac{\partial \hat{\mathbf{G}}_v}{\partial \eta} + \frac{\partial \hat{\mathbf{H}}_v}{\partial \zeta}$$

其中：
- $\hat{\mathbf{Q}} = \mathbf{Q}/J$，$\mathbf{Q} = [\rho, \rho u, \rho v, \rho w, \rho E]^T$
- $\hat{\mathbf{F}}, \hat{\mathbf{G}}, \hat{\mathbf{H}}$ 为曲线坐标系下无粘通量
- $\hat{\mathbf{F}}_v, \hat{\mathbf{G}}_v, \hat{\mathbf{H}}_v$ 为曲线坐标系下粘性通量
- $J = |\partial(\xi,\eta,\zeta)/\partial(x,y,z)|$ 为坐标变换 Jacobian

### 4.2 求解流程总览

```
┌──────────────────────────────────────────────┐
│  1. 初始化                                    │
│     ├ 读取配置文件 (input.cfg)                 │
│     ├ MPI_Init，获取 rank / nprocs            │
│     ├ 读取网格文件（Plot3D / CGNS）            │
│     ├ 区域分解（block → rank 映射）            │
│     ├ 分配各块内存                             │
│     ├ 计算网格度量系数                         │
│     ├ 初始化流场（自由流 / 重启文件）            │
│     └ 建立块间通信拓扑（邻居关系）              │
└───────────────┬──────────────────────────────┘
                ▼
┌──────────────────────────────────────────────┐
│  2. 主时间循环 (n = 0, 1, ..., max_iter)      │
│                                                │
│  ┌──────────────────────────────────────┐     │
│  │ 2.1 施加边界条件 (BCManager)          │     │
│  └──────────┬───────────────────────────┘     │
│             ▼                                  │
│  ┌──────────────────────────────────────┐     │
│  │ 2.2 计算当地时间步长 Δt_local (CFL)   │     │
│  └──────────┬───────────────────────────┘     │
│             ▼                                  │
│  ┌──────────────────────────────────────┐     │
│  │ 2.3 显式 RK 子步循环 (m = 1..N_stage) │     │
│  │                                      │     │
│  │  ┌─────────────────────────────────┐ │     │
│  │  │ a. 交换虚拟网格数据 (HaloExchange)│ │     │
│  │  │    → Q 在 ng 层虚网格处更新       │ │     │
│  │  ├─────────────────────────────────┤ │     │
│  │  │ b. 计算无粘通量 RHS_inv:         │ │     │
│  │  │    ① WCNS 插值：在特征空间内对  │ │     │
│  │  │      每个方向分别进行非线性插值   │ │     │
│  │  │    ② Riemann 求解：在界面处计算 │ │     │
│  │  │      数值通量                     │ │     │
│  │  │    ③ 紧凑/显式差分：对界面通量   │ │     │
│  │  │      沿网格线求导                 │ │     │
│  │  ├─────────────────────────────────┤ │     │
│  │  │ c. 计算粘性通量 RHS_vis:         │ │     │
│  │  │    ① 重构界面处的原始变量        │ │     │
│  │  │    ② 计算速度梯度与应力张量      │ │     │
│  │  │    ③ 计算热通量                  │ │     │
│  │  │    ④ 组装粘性通量                │ │     │
│  │  ├─────────────────────────────────┤ │     │
│  │  │ d. 累积 RHS = RHS_inv + RHS_vis │ │     │
│  │  ├─────────────────────────────────┤ │     │
│  │  │ e. 更新: Q^(m) =                 │ │     │
│  │  │    α_m * Q^(0) + β_m * (Q^(m-1)  │ │     │
│  │  │    + Δt * RHS^(m-1))             │ │     │
│  │  └─────────────────────────────────┘ │     │
│  └──────────┬───────────────────────────┘     │
│             ▼                                  │
│  ┌──────────────────────────────────────┐     │
│  │ 2.4 收敛判断 / 残差计算               │     │
│  └──────────┬───────────────────────────┘     │
│             ▼                                  │
│  ┌──────────────────────────────────────┐     │
│  │ 2.5 按间隔输出 (solution, restart)    │     │
│  └──────────────────────────────────────┘     │
└───────────────┬──────────────────────────────┘
                ▼
┌──────────────────────────────────────────────┐
│  3. 终止                                      │
│     ├ 写入最终解文件                          │
│     ├ 写入重启文件                            │
│     └ MPI_Finalize                           │
└──────────────────────────────────────────────┘
```

### 4.3 WCNS 非线性插值详解（核心算法）

WCNS 的核心思想是将 **非线性插值**（WENO 型，在半节点处）与 **线性差分**（紧凑型，在节点处）解耦，从而同时获得高精度和无振荡的激波捕捉能力。

#### 步骤 1：选择插值方向与变量

对 ξ 方向（η、ζ 方向同理），需在 i+1/2 面处获得左右状态。插值可作用于：
- **守恒变量** Q（简单，但激波处易振荡）
- **特征变量** W = L·Q（推荐，激波捕捉鲁棒性更强）：先投影到特征空间，插值后再反投影
- **通量分量** F（经济，但跨激波时可能不守恒）

#### 步骤 2：构造候选插值模板

以 **WCNS-5**（五阶）为例，对于 i+1/2 面的左状态，取 3 个三点模板：

| 模板 k | 所用节点 | 理想权重 d_k |
|--------|---------|-------------|
| 0 | i-2, i-1, i | 1/16 |
| 1 | i-1, i, i+1 | 10/16 |
| 2 | i, i+1, i+2 | 5/16 |

每个模板给出一个插值 $q_{i+1/2}^{(k)}$：

$$q_{i+1/2}^{(0)} = \tfrac{3}{8}q_{i-2} - \tfrac{5}{4}q_{i-1} + \tfrac{15}{8}q_i$$

$$q_{i+1/2}^{(1)} = -\tfrac{1}{8}q_{i-1} + \tfrac{3}{4}q_i + \tfrac{3}{8}q_{i+1}$$

$$q_{i+1/2}^{(2)} = \tfrac{3}{8}q_i + \tfrac{3}{4}q_{i+1} - \tfrac{1}{8}q_{i+2}$$

#### 步骤 3：计算光滑因子 β_k

Jiang-Shu 型（WCNS-JS）：

$$\beta_k = \sum_{l=1}^{2} \Delta\xi^{2l-1} \int_{\xi_{i-1/2}}^{\xi_{i+1/2}} \left(\frac{d^l p_k(\xi)}{d\xi^l}\right)^2 d\xi$$

对均匀网格，有显式公式（以模板 1 为例）：

$$\beta_1 = \tfrac{13}{12}(q_{i-1} - 2q_i + q_{i+1})^2 + \tfrac{1}{4}(q_{i+1} - q_{i-1})^2$$

#### 步骤 4：计算非线性权重

**WCNS-JS**：

$$\omega_k = \frac{\alpha_k}{\sum \alpha_l}, \quad \alpha_k = \frac{d_k}{(\beta_k + \varepsilon)^p}$$

其中 p=2（推荐），ε ≈ 1e-6 ~ 1e-40。

**WCNS-Z**（更高分辨率）：

$$\tau_5 = |\beta_0 - \beta_2|$$

$$\omega_k = \frac{\alpha_k}{\sum \alpha_l}, \quad \alpha_k = d_k\left[1 + \left(\frac{\tau_5}{\beta_k + \varepsilon}\right)^p\right]$$

#### 步骤 5：加权得到半节点值

$$q_{i+1/2}^L = \sum_{k=0}^{2} \omega_k \cdot q_{i+1/2}^{(k)}$$

右状态 $q_{i+1/2}^R$ 同理（对称地取 i-1,i,i+1,i+2,i+3）。

#### 步骤 6：Riemann 求解

在界面处用左右状态求解 Riemann 问题：

$$\hat{\mathbf{F}}_{i+1/2} = \text{Riemann}(\mathbf{Q}_{i+1/2}^L, \mathbf{Q}_{i+1/2}^R)$$

支持的求解器：
- **Roe**（含熵修正）
- **HLLC**（Harten-Lax-van Leer Contact）
- **AUSM+**（Advection Upstream Splitting Method）

#### 步骤 7：差分求导（节点处）

用线性紧凑差分对界面通量求导得到节点处的 $\partial\hat{\mathbf{F}}/\partial\xi$：

**六阶三对角紧凑格式**：

$$\frac{1}{3}\left(\frac{\partial\hat{\mathbf{F}}}{\partial\xi}\right)_{i-1} + \left(\frac{\partial\hat{\mathbf{F}}}{\partial\xi}\right)_i + \frac{1}{3}\left(\frac{\partial\hat{\mathbf{F}}}{\partial\xi}\right)_{i+1} = \frac{14}{9}\frac{\hat{\mathbf{F}}_{i+1/2} - \hat{\mathbf{F}}_{i-1/2}}{\Delta\xi} + \frac{1}{9}\frac{\hat{\mathbf{F}}_{i+3/2} - \hat{\mathbf{F}}_{i-3/2}}{\Delta\xi}$$

对每条网格线求解三对角系统（Thomas 算法或并行 Pipelined Thomas）。

**WCNS-E**（显式替代，并行友好）：

直接用显式中心差分公式，省去三对角求解。精度略低于紧凑格式但完全无全局依赖，天然并行。

### 4.4 粘性通量计算

在曲线坐标系下，粘性通量包含应力张量和热通量：

1. **重构面心原始变量**：用中心/迎风插值获得 $\mathbf{V}_{i+1/2}$
2. **计算面心速度梯度**：$\partial u_i/\partial x_j$（中心差分 + 度量系数链式法则）
3. **应力张量**：$\tau_{ij} = \mu(\partial u_i/\partial x_j + \partial u_j/\partial x_i) - \tfrac{2}{3}\mu\delta_{ij}\partial u_k/\partial x_k$
4. **热通量**：$q_i = -\kappa \cdot \partial T/\partial x_i$
5. **粘性系数** μ：Sutherland 公式
6. **组装**：粘性通量 = 应力功率 + 热通量

### 4.5 时间推进格式

#### 显式 RK3-TVD（三阶 TVD Runge-Kutta）

$$\mathbf{Q}^{(1)} = \mathbf{Q}^{(0)} + \Delta t \cdot \mathbf{R}(\mathbf{Q}^{(0)})$$

$$\mathbf{Q}^{(2)} = \tfrac{3}{4}\mathbf{Q}^{(0)} + \tfrac{1}{4}\mathbf{Q}^{(1)} + \tfrac{1}{4}\Delta t \cdot \mathbf{R}(\mathbf{Q}^{(1)})$$

$$\mathbf{Q}^{(3)} = \tfrac{1}{3}\mathbf{Q}^{(0)} + \tfrac{2}{3}\mathbf{Q}^{(2)} + \tfrac{2}{3}\Delta t \cdot \mathbf{R}(\mathbf{Q}^{(2)})$$

#### 隐式 LU-SGS（定常计算加速）

求解 $\left(I + \Delta t (D_\xi A + D_\eta B + D_\zeta C)\right) \delta \mathbf{Q} = -\Delta t \cdot \mathbf{R}(\mathbf{Q}^n)$

通过前后扫掠避免大矩阵求逆，CFL 数可达 1e4 以上。

### 4.6 边界条件

对物理边界面的虚拟网格赋值：

| 类型 | 赋值方式 |
|------|---------|
| 远场 | 基于一维特征线理论（Riemann 不变量） |
| 无滑移壁面（绝热） | u=v=w=0，∂T/∂n=0，∂p/∂n=0 |
| 无滑移壁面（等温） | u=v=w=0，T=给定，∂p/∂n=0 |
| 对称面 | 法向速度反向，其余镜像对称 |
| 亚声速入口 | 指定总温总压 + 特征关系 |
| 超声速入口 | 全部变量给定 |
| 亚声速出口 | 指定背压 + 特征关系 |
| 超声速出口 | 全部变量外推 |
| 周期边界 | 配对面上虚拟网格直接拷贝对应内点值：Q_ghost(i) = Q_interior(N-i) |

### 4.7 MPI 并行策略

1. **区域分解**：按网格块分配（1块/进程 或 多块/进程），按格点量均衡负载
2. **虚拟网格交换**：每个 RK 子步前，非阻塞 SendRecv 交换各面 ng 层数据
3. **紧凑差分并行**：沿网格线方向的 Thomas 算法使用流水线并行（Pipelined Thomas）
4. **全局归约**：收敛判断时使用 `MPI_Allreduce` 归约最大/平均残差

### 4.8 推荐的 WCNS 组合方案

| 场景 | WCNS 变种 | 通量格式 | 时间推进 | 说明 |
|------|----------|---------|---------|------|
| 定常亚/跨声速 | WCNS5-Z | Roe/M | LU-SGS | 高分辨率 + 大 CFL 加速收敛 |
| 定常超声速 | WCNS5-JS | HLLC | LU-SGS | 强激波鲁棒性 |
| 非定常（DNS/LES） | WCNS7-Z | Roe | RK4(显式) | 高阶低耗散 |
| 激波-湍流干扰 | WCNS5-CU | HLLC | RK3-TVD | 兼顾分辨率和鲁棒性 |
| 粗网格快速验证 | WCNS5-JS | AUSM+ | RK3-TVD | 最鲁棒的组合 |

---

## 五、开发路线图

1. **Phase 1**：单块、单进程原型
   - 网格读取、度量计算、自由流初始化
   - WCNS5-JS + Roe + RK3-TVD
   - 基本边界条件（远场+壁面）
   - Tecplot 输出

2. **Phase 2**：多块与 MPI 并行
   - 多块网格支持
   - 区域分解 + 虚网格交换
   - 并行紧凑差分（Pipelined Thomas）
   - 残差全局归约

3. **Phase 3**：高阶与高级功能
   - WCNS7 / WCNS-Z
   - LU-SGS 隐式推进
   - 特征空间插值
   - 湍流模型（SA / SST）
   - 重启/断点续算

4. **Phase 4**：验证与优化
   - 制造解精度测试
   - 经典算例验证（平板边界层、RAE2822、ONERA M6、双椭球等）
   - 性能优化（向量化、内存对齐、通信/计算重叠）
   - 负载均衡优化
