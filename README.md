# WCNS v2 —— 三维可压缩 Navier-Stokes 方程并行求解器

基于 WCNS（Weighted Compact Nonlinear Scheme）格式的三维可压缩 Navier-Stokes 方程 CFD 求解器，采用多块结构网格与 MPI 并行计算。

## 构建与运行

```bash
# 配置（需 MPI 和 CGNS 库）
cmake -S . -B build
cmake --build build

# 运行（N 为 MPI 进程数）
mpirun -np N ./build/wcns_v2 <grid.cgns> [config.ini]
```

- C++17，编译选项 `-Wall -Wextra -Wpedantic`
- 依赖：MPI、CGNS（读取网格文件）

---

## 一、文件结构

```
wcns_v2/
├── CMakeLists.txt
├── README.md
├── CLAUDE.md                          # AI 辅助开发指引
│
├── include/wcns_v2/                   # 公共头文件（header-only 库）
│   ├── bc/
│   │   └── bc_applier.h / .hxx       # 边界条件施加器
│   ├── core/
│   │   ├── config.h / .hxx           # 物理参数、计算控制、格式选项
│   │   ├── residual.h / .hxx         # 残差计算、流场监控、收敛判据
│   │   └── history_monitor.h / .hxx  # 截面平均流速历史监控
│   ├── field/
│   │   └── field.h / .hxx            # 流场变量存储（保守/原始/通量/梯度）
│   ├── grid/
│   │   ├── grid.h / .hxx             # 结构网格块，含度量系数与 Jacobian
│   │   ├── boundary_condition.h / .hxx # BC 面片定义
│   │   └── connectivity.h / .hxx     # 块间 1-to-1 连接（周期性/界面）
│   ├── init/
│   │   └── flow_initializer.h / .hxx # 多种初始流场（均匀/泊肃叶/Riemann/涡）
│   ├── io/
│   │   ├── config_reader.h           # INI 配置文件解析器
│   │   ├── cgns_reader.h             # CGNS 网格文件读取器
│   │   ├── solution_writer.h / .hxx  # 解文件输出（Tecplot）
│   │   └── restart_writer.h / .hxx   # 断点续算文件输出
│   ├── parallel/
│   │   ├── parallel_env.h / .hxx     # MPI 初始化/销毁封装
│   │   ├── parallel_manager.h / .hxx # 并行总控（分解、通信、全局归约）
│   │   ├── local_block.h / .hxx      # 本地子域 = Grid + Field + 邻接信息
│   │   ├── halo_exchange.h / .hxx    # 虚拟网格 MPI 交换
│   │   ├── flux_halo_exchange.h / .hxx # 面通量/面值 MPI 交换
│   │   └── decomposer.h / .hxx       # 区域分解算法
│   ├── scheme/
│   │   ├── wcns_interp.h / .hxx      # WCNS 插值（WENO-JS / MDCD Linear / MDCD Hybrid）
│   │   ├── interp_diff.h / .hxx      # 6 阶插值/差分核心算子
│   │   ├── riemann_solver.h / .hxx   # Riemann 求解器（Roe / Rusanov / HLL / HLLC）
│   │   ├── inviscid_rhs.h / .hxx     # 无粘 RHS：通量差分求导
│   │   ├── viscid_rhs.h / .hxx       # 粘性 RHS：梯度→应力→通量→差分
│   │   └── body_force.h / .hxx       # 体积力源项
│   ├── time/
│   │   ├── time_integrator.h / .hxx  # 显式 RK（RK3-TVD / RK4）
│   │   ├── time_step.h / .hxx        # CFL 时间步长计算
│   │   └── lu_sgs.h / .hxx           # 隐式 LU-SGS 时间推进
│   └── utils/
│       └── types.h / .hxx            # 全局类型（Real, Int, MultiArray3D）
│
├── src/
│   ├── main.cpp                      # 入口：配置→初始化→主时间循环→输出
│   └── io/
│       ├── config_reader.cpp         # INI 解析实现
│       └── cgns_reader.cpp           # CGNS 读取实现
│
└── input/
    ├── isentropic/                   # 等熵涡算例（精度验证）
    ├── Riemann_2d/                   # 2D Riemann 问题（激波捕捉）
    ├── poiseuille/                   # Poiseuille 槽道流
    ├── poiseuille_2/                 # 静止启动槽道流（等温壁）
    └── channel_turbulence/           # 湍流槽道初始化
```

---

## 二、程序结构

### 2.1 顶层流程（`main.cpp`）

```
main()
 │
 ├─ 1. ParallelEnv::init()              MPI 初始化
 ├─ 2. ConfigReader::read()             读取 INI 配置文件
 ├─ 3. ParallelManager::initialize()    读 CGNS 网格 → 区域分解 → 构建 LocalBlock → 度量计算
 ├─ 4. FlowInitializer::initialize()    初始流场赋值
 ├─ 5. BC + Halo Exchange               边界条件 → 虚网格交换 → 棱/角修正
 │
 └─ 6. 主时间循环 (iter ≤ max_iter)
      │
      ├─ TimeStep::compute()           计算全局 Δt（CFL 条件）
      │
      ├─ [显式 RK 子步循环 (stage = 0..n_stages-1)]
      │   ├─ WCNS 插值 (ξ,η,ζ)         cell-center → face Q_L/Q_R
      │   ├─ Riemann 求解 (ξ,η,ζ)       Q_L/Q_R → inv_xi/eta/zeta
      │   ├─ 通量虚网格交换            MPI 相邻块界面通量一致
      │   ├─ InviscidRHS::compute()     6 阶中心差分 → rhs
      │   ├─ [粘性 RHS 流水线]
      │   │   ├─ 5a. 插值到面           (u,v,w,T) 面心值
      │   │   ├─ 5b. 梯度计算           链式法则 + 6 阶差分
      │   │   ├─ 5c. 粘性通量           应力张量 + 热通量
      │   │   ├─ 5d. 面通量组装         插值 + 度量点乘
      │   │   └─ 5e. RHS 差分           6 阶差分 → rhs（累加）
      │   ├─ BodyForce::add_to_rhs()    体积力源项（累加到 rhs）
      │   ├─ TimeIntegrator / LuSgs     时间推进，更新 cons
      │   └─ 后处理                     cons→prim → BC → halo → prim→cons
      │
      ├─ 残差 / 监控计算                L2 残差、min/max 物理量
      ├─ 收敛判据 / 发散检测
      ├─ 按间隔输出解文件 / 重启文件
      └─ 进度打印到 stdout
```

### 2.2 模块依赖关系

```
                   ┌─────────────┐
                   │   Config    │  全局参数（所有模块共享）
                   └──────┬──────┘
                          │
    ┌─────────────────────┼─────────────────────┐
    │                     │                     │
    ▼                     ▼                     ▼
┌────────┐  ┌──────────────────────┐  ┌────────────────┐
│  Grid  │  │      LocalBlock      │  │  FlowInit      │
│ 度量   │  │  Grid + Field + 邻接  │  │  初始流场       │
└───┬────┘  └──────────┬───────────┘  └────────────────┘
    │                  │
    ▼                  ▼
┌──────────┐  ┌─────────────────────────────────────┐
│  Field   │  │            Scheme 层                 │
│ prim/cons│  │  WcnsInterp → Riemann → InviscidRHS │
│ rhs/flux │  │  ViscidRHS (5a→5b→5c→5d→5e)        │
│ gradients│  │  BodyForce                           │
└──────────┘  └──────────────┬──────────────────────┘
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
        ┌──────────┐  ┌──────────┐  ┌──────────────┐
        │TimeInteg │  │  LuSgs   │  │    BCApplier │
        │RK3-TVD   │  │  隐式推进 │  │  边界条件     │
        │ RK4      │  └──────────┘  └──────────────┘
        └──────────┘
              │
              ▼
    ┌─────────────────────────────────────┐
    │         Parallel 层                  │
    │  ParallelManager (总控)              │
    │  HaloExchange / FluxHaloExchange    │
    │  Decomposer (区域分解)               │
    │  ParallelEnv (MPI 封装)              │
    └─────────────────────────────────────┘
              │
              ▼
    ┌─────────────────────────────────────┐
    │           I/O 层                     │
    │  SolutionWriter / RestartWriter     │
    │  Residual / HistoryMonitor          │
    └─────────────────────────────────────┘
```

### 2.3 多态/策略模式

| 基类 | 派生类 | 选择方式 |
|------|--------|---------|
| `WcnsInterpBase` | `WcnsNonlinearInterp`, `WcnsMdcdLinearInterp`, `WcnsMdcdHybridInterp` | `cfg.interp_type` |
| `RiemannSolverBase` | `RiemannSolverRoe`, `RiemannSolverRusanov`, `RiemannSolverHLL`, `RiemannSolverHLLC` | `cfg.riemann_type` |

时间推进为静态方法（`TimeIntegrator::advance_stage`, `LuSgs::advance`），通过 `cfg.time_scheme` 分发。

---

## 三、数据结构

### 3.1 全局类型（`utils/types.h`）

```cpp
using Real = double;       // 浮点精度
using Int  = int;          // 整数索引

// 三维连续存储数组：索引 idx = i + ni*(j + nj*k)
template<typename T>
class MultiArray3D {
    // operator()(i,j,k) —— 读写元素
    // allocate(ni,nj,nk) —— 分配内存
    // fill(val) —— 填充常值
    // ni(), nj(), nk() —— 查询维度
};
```

### 3.2 网格块（`grid/grid.h`）

```cpp
class Grid {
    // --- 节点空间（含 ghost）---
    Int ni, nj, nk;                 // 顶点数
    MultiArray3D<Real> node_x, node_y, node_z;

    // --- 单元空间（含 ghost）---
    Int nci, ncj, nck;              // 单元数 (= ni-1)
    MultiArray3D<Real> cell_x, cell_y, cell_z;  // 单元中心坐标
    MultiArray3D<Real> cell_vol;                 // 单元体积

    // --- 物理核心 ---
    Int ni_core, nj_core, nk_core;  // 原始（未扩展）顶点数

    // --- Ghost 层 ---
    Int ng;                          // ghost 层数（WCNS-5 需 ng=3）

    // --- SCMM 度量系数（单元中心）---
    MultiArray3D<Real> met_xi_x, met_xi_y, met_xi_z;    // ξ̂ = J·∇ξ
    MultiArray3D<Real> met_eta_x, met_eta_y, met_eta_z; // η̂ = J·∇η
    MultiArray3D<Real> met_zeta_x, met_zeta_y, met_zeta_z; // ζ̂ = J·∇ζ
    MultiArray3D<Real> jacobian;     // J = |∂(ξ,η,ζ)/∂(x,y,z)|

    // --- 面度量（半节点）---
    // ξ-face(i+1/2): (nci+1)×ncj×nck, etc.
    MultiArray3D<Real> face_xi_x, face_xi_y, face_xi_z;
    MultiArray3D<Real> face_eta_x, face_eta_y, face_eta_z;
    MultiArray3D<Real> face_zeta_x, face_zeta_y, face_zeta_z;

    // --- BC 与连接性 ---
    BoundaryCondition bc;            // BC 面片列表
    ConnectivityList connections;    // 1-to-1 连接（周期性/块间界面）
};
```

**关键索引约定**：物理核心单元索引范围为 `[ng, ng+nci_core-1]` 等。ghost 层在 `[0, ng-1]` 和 `[ng+nci_core, nci-1]`。

### 3.3 流场变量（`field/field.h`）

```cpp
// 原始变量（单元中心，nci×ncj×nck）
struct PrimitiveVars {
    MultiArray3D<Real> rho;   // 密度
    MultiArray3D<Real> u, v, w; // 速度分量
    MultiArray3D<Real> p;     // 压力
};

// 守恒变量（单元中心）
struct ConservativeVars {
    MultiArray3D<Real> rho;   // ρ
    MultiArray3D<Real> rhou;  // ρu
    MultiArray3D<Real> rhov;  // ρv
    MultiArray3D<Real> rhow;  // ρw
    MultiArray3D<Real> rhoE;  // ρE（总能）
};

// 面通量（5 分量）
struct FluxVars {
    MultiArray3D<Real> f1;    // 密度通量
    MultiArray3D<Real> f2;    // x-动量通量
    MultiArray3D<Real> f3;    // y-动量通量
    MultiArray3D<Real> f4;    // z-动量通量
    MultiArray3D<Real> f5;    // 能量通量
};

class Field {
    // --- 核心状态 ---
    PrimitiveVars    prim;
    ConservativeVars cons;
    ConservativeVars Q0;       // Q^(0) — 时间步起始状态（RK 多级用）

    // --- 面插值左右状态 ---
    ConservativeVars ql_xi, qr_xi;      // ξ-face (i+1/2)
    ConservativeVars ql_eta, qr_eta;    // η-face (j+1/2)
    ConservativeVars ql_zeta, qr_zeta;  // ζ-face (k+1/2)

    // --- 无粘/粘性面通量 ---
    FluxVars inv_xi, inv_eta, inv_zeta; // 无粘通量（ξ/η/ζ 面）
    FluxVars vis_xi, vis_eta, vis_zeta; // 粘性通量（ξ/η/ζ 面）

    // --- 单元中心粘性通量 ---
    FluxVars vis_x, vis_y, vis_z;       // Cartesian 粘性通量 F_vis, G_vis, H_vis

    // --- 面插值物理量（粘性用）---
    MultiArray3D<Real> u_face_xi, v_face_xi, w_face_xi, T_face_xi;
    // ... η, ζ 方向同理（共 12 个数组）

    // --- 速度/温度梯度（粘性用，共 12 个数组）---
    MultiArray3D<Real> du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz;
    MultiArray3D<Real> dw_dx, dw_dy, dw_dz, dT_dx, dT_dy, dT_dz;

    // --- 右端项 ---
    ConservativeVars rhs;      // RHS = -(∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ)

    // --- 扩展字段（湍流模型等）---
    std::map<std::string, MultiArray3D<Real>> ext_fields_;
};
```

**面数组维度**（以 nci×ncj×nck 单元为例）：
- ξ-face `(i+1/2)`：`(nci+1) × ncj × nck`
- η-face `(j+1/2)`：`nci × (ncj+1) × nck`
- ζ-face `(k+1/2)`：`nci × ncj × (nck+1)`

### 3.4 本地块与邻接信息（`parallel/local_block.h`）

```cpp
struct NeighborInfo {
    int  target_rank;      // 邻居 MPI rank（-1 = 同一进程）
    int  target_block;     // 邻居全局 block ID
    int  target_face;      // 邻居的连接面（0..5）
    int  transform[3];     // 方向映射（CGNS transform）
    bool active;           // false = 物理边界（无需 halo 交换）
    bool is_periodic;      // true = 周期性连接
};

class LocalBlock {
    int block_id;          // 全局唯一 ID
    int zone_id;           // 原始 zone 编号
    int sub_id;            // zone 内子块编号
    Grid  grid;            // 网格与度量
    Field field;           // 流场
    NeighborInfo neighbors[6]; // 6 个面：IMIN/IMAX/JMIN/JMAX/KMIN/KMAX
};
```

### 3.5 MPI 通信缓冲（`parallel/halo_exchange.h`）

```cpp
class HaloExchange {
    struct FaceBuffer {
        vector<Real> send_buf, recv_buf;  // 打包/解包缓冲
        int  target_rank;
        bool is_remote;                    // 跨进程通信 vs 本地拷贝
        MPI_Request send_req, recv_req;   // 非阻塞通信句柄
    };
    FaceBuffer bufs_[6];                   // 每个面一个缓冲
};
```

### 3.6 参数配置（`core/config.h`）

全部配置参数集中在 `Config` 结构体中，按 section 分类：`[physical]`, `[reference]`, `[control]`, `[scheme]`, `[initialization]`, `[output]`。关键字段见下文各节说明。

---

## 四、数学实现细节

### 4.1 控制方程

曲线坐标系 (ξ, η, ζ) 下的无量纲化三维可压缩 Navier-Stokes 方程：

$$\frac{\partial \hat{\mathbf{Q}}}{\partial t} + \frac{\partial \hat{\mathbf{F}}}{\partial \xi} + \frac{\partial \hat{\mathbf{G}}}{\partial \eta} + \frac{\partial \hat{\mathbf{H}}}{\partial \zeta}
= \frac{\partial \hat{\mathbf{F}}_v}{\partial \xi} + \frac{\partial \hat{\mathbf{G}}_v}{\partial \eta} + \frac{\partial \hat{\mathbf{H}}_v}{\partial \zeta} + \hat{\mathbf{S}}$$

其中：
- $\hat{\mathbf{Q}} = \mathbf{Q}/J$，$\mathbf{Q} = [\rho, \rho u, \rho v, \rho w, \rho E]^T$
- $\hat{\mathbf{F}} = J^{-1}(\xi_x\mathbf{F} + \xi_y\mathbf{G} + \xi_z\mathbf{H})$，无粘通量（ξηζ 方向同理）
- $\hat{\mathbf{F}}_v = J^{-1}(\xi_x\mathbf{F}_v + \xi_y\mathbf{G}_v + \xi_z\mathbf{H}_v)$，粘性通量
- $\hat{\mathbf{S}} = J^{-1}\mathbf{S}$，体积力源项

### 4.2 无量纲化

| 物理量 | 参考量 | 无量纲形式 |
|--------|--------|-----------|
| 密度 | ρ_ref | ρ* = ρ / ρ_ref |
| 速度 | U_ref | u* = u / U_ref |
| 温度 | T_ref | T* = T / T_ref |
| 压力 | ρ_ref·U_ref² | p* = p / (ρ_ref·U_ref²) |
| 长度 | L_ref | L* = L / L_ref |
| 时间 | L_ref/U_ref | t* = t·U_ref / L_ref |
| 粘性 | μ_ref | μ_ref = ρ_ref·U_ref·L_ref / Re |

无量纲状态方程：**p\* = ρ\* · T\* / (γ · Mach²)**（其中 $Mach^2 = U_{ref}^2 / (\gamma R_{gas} T_{ref})$）

### 4.3 WCNS 插值算法

WCNS（Weighted Compact Nonlinear Scheme）的核心思想是将非线性插值（WENO 型，半节点处）与线性差分（节点处）解耦。

#### 4.3.1 半节点约定

半节点 `h` 位于单元 `h-1` 和 `h` 之间（h ∈ [0, n]），Q_L 为左偏插值，Q_R 为右偏（镜像对称）插值。

#### 4.3.2 WENO-JS（5 阶）

**左偏 3 个子模板**（到 i+1/2）：

| 模板 k | 所用节点 | 理想权重 d_k |
|--------|---------|-------------|
| 0 | i-2, i-1, i | 1/16 |
| 1 | i-1, i, i+1 | 10/16 |
| 2 | i, i+1, i+2 | 5/16 |

**光滑因子**（Jiang-Shu）：

$$\beta_k = \frac{13}{12}(q_{k,0} - 2q_{k,1} + q_{k,2})^2 + \frac{1}{4}(q_{k,0} - 4q_{k,1} + 3q_{k,2})^2$$

（k=0 时为 (q_{i-2}, q_{i-1}, q_i)，k=1 时为 (q_{i-1}, q_i, q_{i+1})，k=2 时为 (q_i, q_{i+1}, q_{i+2})）

**非线性权重**（p=2, ε=1e-6）：

$$\omega_k = \frac{\alpha_k}{\sum \alpha_l}, \quad \alpha_k = \frac{d_k}{(\beta_k + \varepsilon)^2}$$

**加权插值**：$q_{i+1/2}^L = \sum \omega_k \cdot v_k$

#### 4.3.3 MDCD Linear（6 阶）

6 点线性插值，含可调耗散/色散系数：

$$q_{i+1/2}^L(\text{diss}, \text{disp}) = \sum_{k=0}^{5} c_k(\text{diss}, \text{disp}) \cdot a[i-2+k]$$

当 diss = disp = 0 时退化为标准 6 阶线性插值。diss > 0 增加耗散，disp 调节色散。

#### 4.3.4 MDCD Hybrid（6 阶）

对每个半节点：
1. 计算间断探测器 $\sigma$（基于相邻二阶差分）
2. 若 $\sigma > \text{sai\_ref}$（间断）：回退到 MDCD Linear（耗散抑制振荡）
3. 若光滑：使用 4 个子模板的 WENO-Z 型非线性重构，含 6 阶全局光滑指示器 $\tau_6$

### 4.4 Riemann 求解器

#### Roe 格式

$$\hat{\mathbf{F}}_{i+1/2} = \frac{1}{2}\left[\hat{\mathbf{F}}(\mathbf{Q}_L) + \hat{\mathbf{F}}(\mathbf{Q}_R)\right] - \frac{1}{2}|\tilde{\mathbf{A}}|(\mathbf{Q}_R - \mathbf{Q}_L)$$

- $\tilde{\mathbf{A}}$ 为面法向通量 Jacobian 在 Roe 平均态下的值
- $|\tilde{\mathbf{A}}| = \mathbf{R}|\mathbf{\Lambda}|\mathbf{R}^{-1}$，$\mathbf{\Lambda} = \text{diag}(\lambda_1, \ldots, \lambda_5)$
- **熵修正**（Harten）：$|\lambda|^* = \max(|\lambda|, \delta/2 + \lambda^2/\delta)$, $\delta = \text{eps} \cdot (|U_n| + c)$
- 特征值：$\lambda_{1,2,3} = U_n$, $\lambda_4 = U_n - c|S|$, $\lambda_5 = U_n + c|S|$

#### Rusanov（Local Lax-Friedrichs）

$$\hat{\mathbf{F}} = \frac{1}{2}(\hat{\mathbf{F}}_L + \hat{\mathbf{F}}_R) - \frac{1}{2}S_{max}(\mathbf{Q}_R - \mathbf{Q}_L), \quad S_{max} = \max(|U_{n,L}| + c_L|S|, |U_{n,R}| + c_R|S|)$$

#### HLL

两波模型，波速估计 $S_L = \min(U_{n,L} - c_L|S|, U_{n,R} - c_R|S|)$，$S_R = \max(U_{n,L} + c_L|S|, U_{n,R} + c_R|S|)$

$$\hat{\mathbf{F}} = \begin{cases} \hat{\mathbf{F}}_L & S_L > 0 \\ \frac{S_R\hat{\mathbf{F}}_L - S_L\hat{\mathbf{F}}_R + S_L S_R(\mathbf{Q}_R - \mathbf{Q}_L)}{S_R - S_L} & S_L \leq 0 \leq S_R \\ \hat{\mathbf{F}}_R & S_R < 0 \end{cases}$$

#### HLLC

在 HLL 基础上恢复接触波（三波模型），精确分辨接触间断，保证正密度性。

### 4.5 SCMM 度量系数

采用对称守恒度量方法（Symmetric Conservative Metric Method）：

$$\hat{\xi}_x = \frac{1}{2}\left[(z y_\eta)_\zeta + (y z_\zeta)_\eta - (z y_\zeta)_\eta - (y z_\eta)_\zeta\right]$$

各方向同理（循环置换），Jacobian 由 $J = \frac{1}{3}\left[(S_\xi)_\xi + (S_\eta)_\eta + (S_\zeta)_\zeta\right]$ 计算，其中 $S_d = x\hat{d}_x + y\hat{d}_y + z\hat{d}_z$。6 阶中心差分求导保证度量系数的自由流保持特性。

### 4.6 无粘 RHS 差分

面通量 → 节点导数的 6 阶中心差分格式：

$$\left(\frac{\partial \hat{F}}{\partial \xi}\right)_i = \frac{1}{\Delta\xi}\sum_{k=-3}^{3} d_k \hat{F}_{i+k-1/2}$$

其中 `d = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60]`（7 点模板，跨越 6 个面通量）。

### 4.7 粘性通量计算（5 步流水线）

| 步骤 | 操作 | 输入 | 输出 |
|------|------|------|------|
| 5a | 6 阶插值 (u,v,w,T) 到面 | cell-center prim | u_face_*/v_face_*/w_face_*/T_face_* |
| 5b | 链式法则求导 | 面值 × 面度量 | du/dx, du/dy, ..., dT/dz（12 数组）|
| 5c | 应力/热通量 → Cartesian 粘性通量 | 梯度 | vis_x/vis_y/vis_z（单元中心）|
| 5d | Cartesian 通量插值到面 → 组装曲线通量 | vis_x/y/z + 面度量 | vis_xi/vis_eta/vis_zeta |
| 5e | 6 阶差分求导 | vis_xi/eta/zeta | rhs（累加到无粘项）|

粘性系数模型：
- `constant`: μ* = mu_const / Re
- `sutherland`: μ* = T^{1.5}·(1+S*)/(T+S*) / Re（S* = 110.4/T_ref）
- `none`: 跳过全部粘性计算

### 4.8 时间推进

#### RK3-TVD（三阶 TVD Runge-Kutta）

| Stage | α | β | 公式 |
|-------|---|---|------|
| 0 | 0 | 1 | Q⁽¹⁾ = Q⁽⁰⁾ + Δt·R(Q⁽⁰⁾) |
| 1 | 3/4 | 1/4 | Q⁽²⁾ = ¾Q⁽⁰⁾ + ¼(Q⁽¹⁾ + Δt·R(Q⁽¹⁾)) |
| 2 | 1/3 | 2/3 | Q⁽³⁾ = ⅓Q⁽⁰⁾ + ⅔(Q⁽²⁾ + Δt·R(Q⁽²⁾)) |

更新式：Q_new = α·Q0 + β·(Q_old + Δt·J·rhs)，仅更新内部单元 [3..nci-4]。

#### RK4（经典四阶）

4 级 4 阶，Butcher 系数：α=(1,1,1,0)，β=(0.5,0.5,1,1/6)。

#### LU-SGS（标量对角近似，隐式）

求解 (D+L)D⁻¹(D+U)ΔQ = RHS：

1. **对角**：D = Ω/Δt + κ·(σ_ξ + σ_η + σ_ζ) + 粘性贡献
2. **前向扫描**（i↑, j↑, k↑）：ΔQ* = D⁻¹[Ω·RHS + ½max(σ_L,σ_R)·ΔQ*_{neighbor}]
3. **后向扫描**（i↓, j↓, k↓）：ΔQ = ΔQ* − D⁻¹[½max(σ_L,σ_R)·ΔQ_{neighbor}]
4. **更新**：Q^{n+1} = Q^n + ΔQ

其中 σ 为谱半径，κ 为超松弛因子（κ > 1 增强对角占优，6 阶中心差分需 κ ≥ 4 维持稳定性）。

#### 时间步长（CFL 条件）

$$\Delta t = \text{CFL} \cdot \min_i \frac{\Omega_i}{\sum_d (\Lambda_c^d + \Lambda_v^d)}$$

对流谱半径：$\Lambda_c^d = |U_{contra}^d| + c \cdot |S^d|$
粘性谱半径：$\Lambda_v^d = \beta \cdot (\mu/\rho)/Pr \cdot |S^d|^2 / \Omega$

### 4.9 边界条件

施加顺序：**面 ghost → halo 交换 → 棱 ghost → 角 ghost**

| BC 类型 | 施加方式 | 适用场景 |
|---------|---------|---------|
| adiabatic_noslip | u=v=w=0, ∂T/∂n=0, ∂p/∂n=0 | 绝热固壁 |
| isothermal_noslip | u=v=w=0, T=给定, ∂p/∂n=0 | 等温固壁 |
| slip | 法向速度反射, 其余对称 | 无粘壁面/对称面 |
| farfield | 基于 Riemann 不变量的特征 BC | 远场边界 |
| symmetry | 法向速度反向, 其余镜像 | 对称面 |
| inflow | 指定自由流值的亚/超声速入口 | 入口 |
| outflow | 全部外推的亚/超声速出口 | 出口 |

### 4.10 MPI 并行策略

1. **区域分解**：按块分配（或块内 Cartesian 分割），负载近似均衡
2. **虚拟网格交换**：每个 RK 子步后，非阻塞 SendRecv 交换 6 个面 ng 层数据
3. **通量交换**：Riemann 求解后交换 ng+1 层面通量，确保跨块界面通量一致
4. **全局归约**：残差/监控计算使用 MPI_Allreduce

---

## 五、算法过程（简洁版）

```
主循环 (iter):
│
├─ 时间步长: Δt = CFL × min(Ω / ΣΛ)
│
├─ [RK 子步 stage:]
│   ├─ [stage=0] 快照 Q0 = cons
│   │
│   ├─ WCNS 插值:  cell-center → face Q_L/Q_R
│   ├─ Riemann:    Q_L/Q_R → inv_xi/eta/zeta
│   ├─ [仅多块]    通量 halo 交换
│   ├─ 无粘 RHS:  6 阶中心差分 → rhs
│   │
│   ├─ [粘性]:
│   │   ├─ 5a: 插值 (u,v,w,T) 到面
│   │   ├─ 5b: 计算梯度 (12 个梯度分量)
│   │   ├─ 5c: 单元中心 Cartesian 粘性通量
│   │   ├─ 5d: 面通量组装
│   │   └─ 5e: 粘性 RHS → rhs (累加)
│   │
│   ├─ 体积力: rhs += J⁻¹·S
│   ├─ 时间推进: cons = α·Q0 + β·(cons + dt·J·rhs)
│   │            或 LU-SGS: forward → backward → update
│   │
│   └─ 后处理:
│       cons→prim → BC (面 ghost) → halo 交换 → BC (棱/角) → prim→cons
│
├─ 残差 & 发散检测
├─ 输出 (按间隔)
└─ 进度打印
```

---

## 六、可扩展模块

### 6.1 已预留 / 计划中

| 扩展项 | 状态 | 说明 |
|--------|------|------|
| 特征空间插值 | 计划已就绪 | `interp_vars="characteristic"` — 左/右特征向量投影 + 标量插值 + 反投影，激波处鲁棒性更优 |
| 湍流模型 | 预留接口 | `Field::add_extension("k_turb")` 等可动态添加湍流变量，`FlowMonitor` 预留 extras |
| CGNS 输出 | 预留 | `output_format="cgns"` 在 config 中已定义 |
| VTK 输出 | 预留 | `output_format="vtk"` |
| 更多初始条件 | 易扩展 | `FlowInitializer` 使用静态方法 dispatch 模式，添加新 `init_type` 仅需一个新方法 |
| 更多 BC 类型 | 易扩展 | `BCApplier` 按 `BCType` 分发，添加新类型仅需实现 `apply_*()` 方法 |
| 更多 Riemann 求解器 | 易扩展 | `RiemannSolverBase` 工厂模式，添加新求解器仅需继承基类 |

### 6.2 建议扩展方向

| 方向 | 难度 | 影响 |
|------|------|------|
| 双时间步（Dual Time-Stepping） | 中等 | LU-SGS 作为内迭代求解器 → 非定常问题大 CFL + 高阶时间精度 |
| 矩阵值 LU-SGS | 较大 | 5×5 块对角替代标量对角 → 更好隐式耦合 → 更高精度 |
| WENO-Z / WENO-CU 变体 | 较小 | 替代现有 WENO-JS 权重公式 → 更高分辨率 |
| 高阶紧凑差分（8 阶三对角/五对角） | 中等 | 更陡峭的修正波数特性 |
| GPU 加速（CUDA/HIP） | 大 | 大规模计算提速 |
| 自适应网格加密（AMR） | 大 | 激波/涡结构局部高分辨率 |
| RANS 湍流模型（SA / k-ω SST） | 中等 | 工程实用定常问题 |
| LES 亚格子模型 | 中等 | 非定常大涡模拟 |

---

## 七、可用算例

### 7.1 等熵涡 (Isentropic Vortex) — `input/isentropic/`

Euler 方程的精确解：涡被平均流无失真对流输运。用于验证空间格式的形式精度。

| 文件 | 说明 |
|------|------|
| `input.ini` | mdcd_hybrid, CFL=0.5, t_max=1.0 |
| `input_uniform.ini` | mdcd_linear, CFL=0.5, uniform metrics |
| `input_lu_sgs.ini` | LU-SGS (κ=4.0), CFL=0.5 |
| `input_rk3_debug.ini` | RK3 调试用（少迭代） |

关键参数：`isentropic_vortex_strength=1.0, radius=0.05, xc=0.5, yc=0.5, u_inf=1, v_inf=1`

### 7.2 二维 Riemann 问题 — `input/Riemann_2d/`

Lax & Liu (1998) 配置 3：四激波相互作用。验证激波捕捉和接触间断分辨能力。

| 文件 | 说明 |
|------|------|
| `input.ini` | weno_js + roe, CFL=0.3, t_max=0.3 |

### 7.3 Poiseuille 槽道流 — `input/poiseuille/`

体积力驱动的定常层流：抛物速度型与理论解对比。验证粘性项和固壁 BC。

| 参数 | 值 |
|------|-----|
| Re | 100 |
| Mach | 0.1 |
| 时间格式 | RK3-TVD |
| 粘性 | constant (μ*=1/Re) |
| BC | adiabatic_noslip |

### 7.4 Poiseuille 静止启动 — `input/poiseuille_2/`

从静止态开始，体积力加速至定常抛物型。验证等温壁 BC 和能量方程耦合。

| 参数 | 值 |
|------|-----|
| Re | 20 |
| Mach | 0.2 |
| 初始 | rest（零速） |
| BC | isothermal_noslip |

### 7.5 湍流槽道初始化 — `input/channel_turbulence/`

Spalding 壁面律平均剖面 + 正弦扰动叠加。用于湍流模拟的初始条件。

| 参数 | 值 |
|------|-----|
| Re | 180（摩擦雷诺数 Re_τ） |
| Mach | 0.1 |
| BC | isothermal_noslip |

---

## 八、配置文件快速参考

INI 格式，section `[name]`，key=value。`#` 或 `;` 开头为注释。

```ini
[physical]          # 物理参数
gamma   = 1.4       # 比热比
prandtl = 0.72      # Prandtl 数
re      = 1.0e6     # Reynolds 数
mach    = 0.5       # 马赫数

[reference]         # 参考量（有量纲，用于无量纲化）
length = 1.0        # 参考长度
rho    = 1.225      # 参考密度
temp   = 288.15     # 参考温度

[control]           # 计算控制
cfl           = 0.5     # CFL 数
max_iter      = 10000   # 最大迭代步数
max_time      = 1.0     # 最大物理时间（0=不限）
time_scheme   = rk3-tvd # rk3-tvd | rk4 | lu-sgs
lu_sgs_kappa  = 4.0     # LU-SGS 超松弛因子
ng            = 3       # Ghost 层数

[scheme]            # 数值格式
interp_type     = weno_js     # weno_js | mdcd_linear | mdcd_hybrid
interp_vars     = conservative # conservative（characteristic 计划中）
riemann_type    = roe         # roe | rusanov | hll | hllc
viscous_type    = none        # none | constant | sutherland
mdcd_diss       = 0.02        # MDCD 耗散系数
mdcd_disp       = 0.046       # MDCD 色散系数

[initialization]    # 初始条件
init_type = uniform | poiseuille | rest | riemann_2d | isentropic_vortex | channel_turbulence
# ... 各类型特定参数见 config.h

[output]            # 输出控制
output_format = tecplot
output_dir    = ./output
```

---

## 九、未来改进

1. **特征空间插值**：`interp_vars="characteristic"` 实现（计划已就绪），显著提升激波捕捉鲁棒性
2. **双时间步**：LU-SGS 内迭代 + 二阶 BDF 外迭代 → 非定常问题兼顾大 CFL 和高阶时间精度
3. **矩阵值 LU-SGS**：5×5 块对角替代标量对角 → 隐式耦合更强，非定常精度显著提高
4. **湍流模型**：SA（一方程）和 k-ω SST（两方程）RANS 模型，利用 Field 扩展字段接口
5. **WENO-Z/CU 变体**：更高分辨率的非线性权重公式
6. **8 阶紧凑差分**：三对角/五对角格式的 8 阶版本
7. **VTK/CGNS 输出**：支持更多后处理格式
8. **GPU 加速**：CUDA/HIP 后端用于大规模计算
9. **自适应网格加密（AMR）**：局部高分辨率激波/涡结构
10. **重启功能完善**：从 restart_*.bin 恢复计算状态
11. **二阶精度验证**：制造解方法（Method of Manufactured Solutions）验证格式收敛阶
12. **负载均衡优化**：加权区域分解，考虑粘性计算开销的异构性
