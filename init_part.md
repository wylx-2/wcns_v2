# WCNS v2 — 初始化部分功能与算法说明

## 目录

1. [总体概述](#1-总体概述)
2. [源文件组织架构](#2-源文件组织架构)
3. [初始化执行流程](#3-初始化执行流程)
4. [各模块详细说明](#4-各模块详细说明)
   - [4.1 MPI 环境初始化](#41-mpi-环境初始化)
   - [4.2 配置文件读取](#42-配置文件读取)
   - [4.3 CGNS 网格文件读取](#43-cgns-网格文件读取)
   - [4.4 Ghost 层扩展](#44-ghost-层扩展)
   - [4.5 网格度量系数计算](#45-网格度量系数计算)
   - [4.6 区域分解](#46-区域分解)
   - [4.7 本块构建](#47-本块构建)
   - [4.8 Halo 交换建立](#48-halo-交换建立)
   - [4.9 流场初始化](#49-流场初始化)
   - [4.10 边界条件应用](#410-边界条件应用)
   - [4.11 Halo 交换执行](#411-halo-交换执行)
5. [数据结构总览](#5-数据结构总览)
6. [网格索引布局](#6-网格索引布局)
7. [当前功能完备性检查](#7-当前功能完备性检查)
8. [后续待实现模块](#8-后续待实现模块)

---

## 1. 总体概述

WCNS v2 是一个基于 WCNS（Weighted Compact Nonlinear Scheme）格式的 3D 可压缩 Navier-Stokes 并行求解器，采用结构化多块网格和 MPI 并行。求解器使用有限体积法/有限差分法混合格式，物理量存储在单元中心（cell-centered），边界条件通过 ghost cell 施加。

本文档描述求解器的**初始化阶段**（即进入时间推进循环之前的所有步骤），包括从启动到完成所有前处理的完整流程。

### 无量纲化约定

所有内部计算均使用无量纲量：

| 量 | 无量纲化 | 符号 |
|----|----------|------|
| 密度 | ρ* = ρ / ρ_ref | rho |
| 速度 | u* = u / U_ref | u, v, w |
| 温度 | T* = T / T_ref | T |
| 压力 | p* = p / (ρ_ref · U_ref²) | p |
| 长度 | L* = L / L_ref | x, y, z |
| 粘度 | μ* = μ / μ_ref | mu |

状态方程（无量纲）：

```
p* = ρ* · T* / (γ · Mach²) = ρ* · T* · eos_factor
```

其中 `eos_factor = 1 / (γ · Mach²)`。

---

## 2. 源文件组织架构

```
wcns_v2/
├── CMakeLists.txt                          # CMake 构建（C++17 + MPI + CGNS）
├── input.ini                               # 算例配置文件
├── init_part.md                            # 本文档
│
├── include/wcns_v2/
│   ├── utils/
│   │   ├── types.h / .hxx                  # 基础类型：Real, Int, MultiArray3D
│   │
│   ├── core/
│   │   ├── config.h / .hxx                 # 物理参数/控制参数的 Config 结构体
│   │
│   ├── grid/
│   │   ├── grid.h / .hxx                   # 结构化网格 Grid 类
│   │   ├── boundary_condition.h / .hxx     # BCType 枚举 & BCPatch/BC 类
│   │   ├── connectivity.h / .hxx           # 1-to-1 连接（块接口/周期边界）
│   │
│   ├── field/
│   │   ├── field.h / .hxx                  # 流场存储 Field 类（原始/守恒/通量）
│   │
│   ├── scheme/
│   │   ├── interp_diff.h / .hxx            # 高阶插值与差分算子（WCNS 6阶中心+5阶边界）
│   │
│   ├── io/
│   │   ├── config_reader.h                 # INI 文件解析器
│   │   ├── cgns_reader.h                   # CGNS 网格文件读取器
│   │
│   ├── parallel/
│   │   ├── parallel_env.h / .hxx           # MPI 环境封装（单例）
│   │   ├── decomposer.h / .hxx             # 区域分解（贪婪+笛卡尔拆分）
│   │   ├── local_block.h / .hxx            # 本块：Grid + Field + NeighborInfo
│   │   ├── halo_exchange.h / .hxx          # Halo 数据交换（MPI 阻塞/非阻塞）
│   │   ├── parallel_manager.h / .hxx       # 并行管理器（顶层调度）
│   │
│   ├── init/
│   │   ├── flow_initializer.h / .hxx       # 流场初始化（均匀流/泊肃叶流）
│   │
│   ├── bc/
│   │   ├── bc_applier.h / .hxx             # 边界条件应用（face→edge→corner）
│   │
│   └── time/                               # （待实现）时间推进
│
├── src/
│   ├── main.cpp                            # 主程序入口
│   ├── io/
│   │   ├── config_reader.cpp               # ConfigReader 实现
│   │   ├── cgns_reader.cpp                 # CGNSReader 实现
│   │
└── build/                                  # CMake 构建目录
```

### 文件归类

| 类别 | 文件数 | 说明 |
|------|--------|------|
| 基础工具 | 1 | types.h |
| 核心配置 | 1 | config.h |
| 网格 | 3 | grid.h, boundary_condition.h, connectivity.h |
| 流场 | 1 | field.h |
| 数值格式 | 1 | interp_diff.h |
| 输入输出 | 2 | config_reader, cgns_reader |
| 并行 | 5 | env, decomposer, local_block, halo_exchange, manager |
| 初始化 | 1 | flow_initializer.h |
| 边界条件 | 1 | bc_applier.h |
| 主程序 | 1 | main.cpp |
| **总计** | **17 个模块** | |

---

## 3. 初始化执行流程

```
main()
 ├── [1] MPI_Init → ParallelEnv::init()
 ├── [2] 解析命令行参数（网格文件、配置文件）
 ├── [3] ConfigReader::read() → Config
 │        ├── 解析 [physical]    (gamma, Mach, Re, Pr, AoA, beta)
 │        ├── 解析 [reference]   (L_ref, rho_ref, T_ref, R_gas)
 │        ├── 解析 [control]     (CFL, max_iter, time_scheme, ng, ...)
 │        ├── 解析 [initialization] (init_type, body_force, wall_type, ...)
 │        └── cfg.finalize() → 计算导出量 (U_ref, p_ref, mu_ref)
 ├── [4] Field 转换测试（prim ⇄ cons 精度验证）
 ├── [5] ParallelManager::initialize()
 │        ├── [5.1] CGNSReader::read_all() — 所有进程独立读取网格文件
 │        │         ├── read_zone() × Nzones
 │        │         │    ├── cg_zone_read() → 维数
 │        │         │    ├── grid.allocate(ni, nj, nk, ng)
 │        │         │    ├── read_coords() → node_x/y/z
 │        │         │    ├── compute_cell_centers() → 8顶点平均
 │        │         │    ├── compute_cell_volumes() → 均匀直角近似
 │        │         │    ├── read_bc() → BCPatch 列表
 │        │         │    └── read_1to1() → Connectivity 列表
 │        ├── [5.2] extend_ghost_layers() — 扩展 ghost 节点
 │        │         ├── 保存原始节点坐标
 │        │         ├── 分配扩展数组 (ni+2ng)×(nj+2ng)×(nk+2ng)
 │        │         ├── 拷贝核心节点至偏移 ng 位置
 │        │         ├── fill_ghost_nodes()
 │        │         │    ├── 周期面 → fill_ghost_face_periodic()
 │        │         │    └── 非周期面 → fill_ghost_face_extrapolate()
 │        │         ├── compute_cell_centers() → 扩展后重算
 │        │         └── compute_cell_volumes()
 │        ├── [5.3] compute_metrics() → SCMM 度量系数
 │        │         ├── 计算坐标导数 (x_ξ, x_η, x_ζ 等)
 │        │         ├── 计算 9 个度量项 (met_xi_x...met_zeta_z)
 │        │         └── 计算 Jacobian
 │        ├── [5.4] compute_face_metrics() → 半节点面度量
 │        ├── [5.5] BlockDecomposer::decompose() → 区域分解
 │        │         ├── 贪婪分配 (nprocs ≤ nzones)
 │        │         └── 笛卡尔拆分 (nprocs > nzones)
 │        ├── [5.6] 构建 LocalBlock
 │        │         ├── from_full_zone() — 完整 zone
 │        │         └── from_sub_zone() — 拆分 zone（重算度量）
 │        ├── [5.7] HaloExchange::setup() — 预分配通信缓冲区
 │        └── [5.8] 打印并行设置摘要
 ├── [6] FlowInitializer::initialize() × Nblocks
 │        ├── init_type = "uniform"    → init_uniform()
 │        └── init_type = "poiseuille" → init_poiseuille()
 │        └── prim_to_cons() — 同步守恒量
 ├── [7] BoundaryConditionApplier::apply_all() × Nblocks
 │        ├── Stage 1: apply_face_ghost() — 6面 ghost 填充
 │        │    ├── 确定每面 BC 类型（从 BCPatch）
 │        │    ├── 对整个面应用 BC（全切向范围）
 │        │    │    ├── Farfield → 1D Riemann 不变量
 │        │    │    ├── Wall (no-slip) → 偶反射
 │        │    │    ├── Wall (slip) → 法向反射
 │        │    │    ├── Symmetry → 法向反射
 │        │    │    ├── Inflow → 指定来流
 │        │    │    └── Outflow → 零阶外推
 │        │    └── 自周期面：从对面内点拷贝
 │        ├── Stage 2: apply_edge_ghost() — 12 条边平均
 │        └── Stage 3: apply_corner_ghost() — 8 个角平均
 ├── [8] ParallelManager::exchange_all_halos() — MPI halo 交换
 ├── [9] verify_initialization() — 验证内点/ghost 值
 └── [10] Metric 验证 (GCL 残差 + Jacobian 范围)
```

### 主程序伪代码

```cpp
int main(int argc, char* argv[]) {
    // 1. MPI 初始化
    ParallelEnv::init(argc, argv);

    // 2. 读取配置
    Config cfg = ConfigReader::read(argv[2]);   // 或使用默认参数
    cfg.finalize();
    print_config(cfg);

    // 3. 场转换测试（正确性验证）
    test_field_conversion(cfg);

    // 4. 并行初始化（网格读取 + 分解 + LocalBlock 构建）
    ParallelManager pm;
    std::vector<LocalBlock> local_blocks;
    pm.initialize(argv[1], cfg, local_blocks);

    // 5. 流场初始化（仅内点）
    for (auto& lb : local_blocks)
        FlowInitializer::initialize(lb, cfg);

    // 6. 边界条件应用（ghost cells）
    for (auto& lb : local_blocks)
        BoundaryConditionApplier::apply_all(lb, cfg);

    // 7. Halo 交换（并行 ghost cells）
    pm.exchange_all_halos(local_blocks);

    // 8. 验证
    verify_initialization(local_blocks, cfg);
    // ... GCL 验证 ...

    ParallelEnv::finalize();
}
```

---

## 4. 各模块详细说明

### 4.1 MPI 环境初始化

**文件**: [parallel_env.h](include/wcns_v2/parallel/parallel_env.h) / [.hxx](include/wcns_v2/parallel/parallel_env.hxx)

**类**: `ParallelEnv`（全静态方法，单例模式）

**功能**:
- 调用 `MPI_Init_thread` 初始化 MPI（请求 `MPI_THREAD_MULTIPLE` 支持）
- 获取 rank / nprocs
- 提供 `barrier()`, `is_master()`, `is_parallel()` 等便捷查询
- `finalize()` 调用 `MPI_Finalize`

**关键特性**: 单进程模式（无 mpirun）下也能正常工作，rank=0, size=1。

### 4.2 配置文件读取

**文件**: [config_reader.h](include/wcns_v2/io/config_reader.h) / [config_reader.cpp](src/io/config_reader.cpp)

**类**: `ConfigReader`

**功能**: 解析 INI 格式配置文件，支持 `#` / `;` 注释和行内注释。

**Section 分派**:

| Section | 方法 | 解析参数 |
|---------|------|---------|
| `[physical]` | `set_physical()` | gamma, Prandtl, Re, Mach, AoA, beta |
| `[reference]` | `set_reference()` | L_ref, rho_ref, T_ref, R_gas |
| `[control]` | `set_control()` | CFL, max_iter, output_freq, time_scheme, ng, converge_tol |
| `[initialization]` | `set_initialization()` | init_type, body_force_x/y/z, poiseuille_umax/y_min/y_max, wall_type |

**Config 导出量** (由 `cfg.finalize()` 计算):

```
c_ref  = sqrt(γ · R_gas · T_ref)         — 参考声速
U_ref  = Mach · c_ref                     — 参考速度
p_ref  = rho_ref · U_ref²                 — 参考压力
mu_ref = rho_ref · U_ref · L_ref / Re     — 参考粘度
```

**自由来流速度** (由 `cfg.free_stream_velocity()` 计算):

```
u_inf = cos(AoA) · cos(beta)              — 无量纲，模为 1
v_inf = sin(beta)
w_inf = sin(AoA) · cos(beta)
```

### 4.3 CGNS 网格文件读取

**文件**: [cgns_reader.h](include/wcns_v2/io/cgns_reader.h) / [cgns_reader.cpp](src/io/cgns_reader.cpp)

**类**: `CGNSReader`

**功能**: 读取 CGNS 格式的结构化网格文件。所有 MPI 进程独立读取同一文件（无需 MPI I/O）。

**read_zone() 流程**:

1. **读取 zone 信息**: `cg_zone_read()` → 节点维数 (ni, nj, nk)
2. **分配数组**: `grid.allocate(ni, nj, nk, ng)` — 此时为核心（物理）维数
3. **读取坐标**: `cg_coord_read()` → `node_x`, `node_y`, `node_z` (vertex-centered)
   - 先尝试 `RealDouble`，失败则回退 `RealSingle`
4. **计算单元中心**: `grid.compute_cell_centers()` — 8 顶点平均
5. **计算单元体积**: `grid.compute_cell_volumes()` — 均匀直角近似
6. **读取边界条件**: `read_bc()` → BCType 枚举 + 面检测
   - CGNS BC 类型 → BCType 枚举映射（Wall, Farfield, Inflow, Outflow, Symmetry, Periodic）
   - `detect_face()`: 通过缩并维度确定面编号（0=IMIN, 1=IMAX, ..., 5=KMAX）
7. **读取 1-to-1 连接**: `read_1to1()` → Connectivity 列表
   - 周期连接信息（平移向量、旋转中心/角度）
   - Transform 数组（方向映射）

### 4.4 Ghost 层扩展

**文件**: [grid.hxx](include/wcns_v2/grid/grid.hxx) `extend_ghost_layers()`

**功能**: 将核心网格的所有数组扩展 ng 层 ghost。

**步骤**:

1. 保存原始节点坐标 (`old_nx/ny/nz`)，记 `ni_src = ni_core`
2. 分配扩展数组: `(ni_src+2·ng) × (nj_src+2·ng) × (nk_src+2·ng)`
3. 拷贝核心节点至偏移 ng 位置: `new(i+ng, j+ng, k+ng) = old(i, j, k)`
4. 更新维数: `ni = ni_src+2ng`, `nci = ni-1`, 等
5. **填充 ghost 节点**: `fill_ghost_nodes()`
   - 周期面 → `fill_ghost_face_periodic()`: 从对面 donor 拷贝并加平移向量
   - 非周期面 → `fill_ghost_face_extrapolate()`: **逐层线性外推**
     ```
     ghost(i) = 2 · base(i) - prev(i)
     ```
     其中 base 是比 ghost 近一层的已知值，prev 是近两层的已知值。逐层从内向外进行。
6. 重新计算扩展后的 `cell_x/y/z` 和 `cell_vol`

**维数关系**:
```
扩展前（核心）:  ni_core, nj_core, nk_core       (vertex)
                nci_core = ni_core - 1            (cell)
扩展后（总计）:  ni = ni_core + 2·ng
                nci = ni - 1 = nci_core + 2·ng
```

### 4.5 网格度量系数计算

**文件**: [grid.hxx](include/wcns_v2/grid/grid.hxx) `compute_metrics()` / `compute_face_metrics()`

**算法**: SCMM（Symmetric Conservative Metric Method — 对称守恒度量方法）

#### 4.5.1 单元中心度量 `compute_metrics()`

使用**插值→微分**两步法（与 WCNS 非线性格式的线性部分一致）：

**第一步**: 用 6 阶中心差分计算坐标导数（物理空间对计算空间）:

```
x_ξ, x_η, x_ζ,  y_ξ, y_η, y_ζ,  z_ξ, z_η, z_ζ
```
调用 `InterpDiff::derivative()`，在 ξ/η/ζ 方向分别计算。

**第二步**: 用 SCMM 公式计算度量项。以 ξ̂_x 为例:

```
ξ̂_x = ½[(z·y_η)_ζ + (y·z_ζ)_η - (z·y_ζ)_η - (y·z_η)_ζ]
```

类似计算其余 8 个分量（η̂_x/y/z, ζ̂_x/y/z）。每个分量的计算模式为：两个正项 + 两个负项的差分组合。

**第三步**: 计算 Jacobian（守恒形式的 1/J）:

```
S_ξ  = x·ξ̂_x + y·ξ̂_y + z·ξ̂_z
S_η  = x·η̂_x + y·η̂_y + z·η̂_z
S_ζ  = x·ζ̂_x + y·ζ̂_y + z·ζ̂_z

1/J = ⅓[(S_ξ)_ξ + (S_η)_η + (S_ζ)_ζ]
J   = 3 / (1/J)
```

**边界处理**: 非周期边界使用 5 阶单边差分，周期边界使用 6 阶中心差分（利用 ghost 中的周期数据）。

#### 4.5.2 面度量 `compute_face_metrics()`

将单元中心度量插值到半节点（面中心）:

- ξ 面 (i+1/2): `(nci+1) × ncj × nck`
- η 面 (j+1/2): `nci × (ncj+1) × nck`
- ζ 面 (k+1/2): `nci × ncj × (nck+1)`

调用 `InterpDiff::interp_to_faces()`，内部边界同样区分中心/单边格式。

#### 4.5.3 插值/差分格式 (`InterpDiff`)

**文件**: [interp_diff.h](include/wcns_v2/scheme/interp_diff.h)

| 位置 | 插值 | 差分 |
|------|------|------|
| 内部 | 6 阶中心 (7 点) | 6 阶中心 (7 点半节点) |
| 左边界 1,2 层 | 5 阶单边 | 5 阶单边 |
| 右边界 1,2 层 | 5 阶单边 | 5 阶单边 |

第 3 层及之后的边界层退化到由 ghost 值支撑的中心格式。

### 4.6 区域分解

**文件**: [decomposer.h](include/wcns_v2/parallel/decomposer.h) / [.hxx](include/wcns_v2/parallel/decomposer.hxx)

**类**: `BlockDecomposer`

**算法**:

1. **收集各 zone 单元数** (core cells): `(ni_core-1) × (nj_core-1) × (nk_core-1)`

2. **Phase 1 — 贪婪分配** (`assign_greedy`):
   - 按 cell 数从大到小排序
   - 贪心将每个 zone 分配给当前负载最小的进程
   - 适用于 `nprocs ≤ nzones` 情况

3. **Phase 2 — 拆分** (当 `nprocs > nzones`):
   - 计算所需额外块数: `n_extra = nprocs - nzones`
   - 按 zone 大小降序分配拆分份数
   - `find_split_factors()`: 寻找使子块形状最接近立方的分解因子
     - 约束: 每个子块 ≥ `2·ng` cells 每方向
     - 评分: 各方向尺寸与平均值偏差 + 过度拆分的惩罚
   - 轮转分配子块到各进程

**SubBlock 结构**: 记录在**原始 zone 坐标**中的 cell 范围（0-based），用于后续子块提取。

### 4.7 本块构建

**文件**: [local_block.h](include/wcns_v2/parallel/local_block.h) / [.hxx](include/wcns_v2/parallel/local_block.hxx)

**类**: `LocalBlock`

**两种构建路径**:

#### from_full_zone() — zone 未被拆分

- 直接深拷贝整个 Grid（含所有度量、BC、连接）
- 分配 Field: `(nci, ncj, nck)` 匹配 grid
- 构建 `NeighborInfo[6]`

#### from_sub_zone() — zone 被拆分

- 从完整 grid 提取子区域节点坐标
- 子块坐标映射: `sub(i,j,k) = full(i + offset_i, j + offset_j, k + offset_k)`
  其中 `offset = (ci_min, cj_min, ck_min)`（对应 full zone 的 node 索引）
- **重算** cell centers、cell volumes、metrics、face metrics
- **重映射** BC patches: 只保留与子块重叠的 patch，转换到子块局部坐标
- 分配 Field

#### NeighborInfo 构建 (`build_neighbors()`)

对每个面（0=IMIN ... 5=KMAX）:

| 条件 | NeighborInfo | 说明 |
|------|-------------|------|
| 面在原始 zone 边界上 + 周期连接 | active=true, is_periodic=true, rank=-1 | 自周期面 |
| 面在原始 zone 边界上 + 无周期连接 | active=false | BC 面，由 BCApplier 处理 |
| 内部拆分面 | active=true, is_periodic=false | 查找相邻子块，记录 rank |

### 4.8 Halo 交换建立

**文件**: [halo_exchange.h](include/wcns_v2/parallel/halo_exchange.h) / [.hxx](include/wcns_v2/parallel/halo_exchange.hxx)

**类**: `HaloExchange`

**`setup()` 流程**:
- 对 6 个面各初始化 `FaceBuffer`
- 远程邻居面：预分配 send_buf / recv_buf
  - I-face: `ng × ncj × nck`
  - J-face: `nci × ng × nck`
  - K-face: `nci × ncj × ng`

**`exchange_multi()` 流程**（阻塞式）:
1. 将所有数组（rho, u, v, w, p）顺序 pack 到 send buffer
2. 从**内点**区域提取数据（IMIN: i=ng..2ng-1, IMAX: i=nci-2ng..nci-ng-1）
3. MPI_Isend + MPI_Irecv → MPI_Wait
4. Unpack 到 ghost 区域（IMIN: i=0..ng-1, IMAX: i=nci-ng..nci-1）

**关键点**: send buffer 从内点（非 ghost）提取数据，确保发送的是已正确计算的值。

### 4.9 流场初始化

**文件**: [flow_initializer.h](include/wcns_v2/init/flow_initializer.h) / [.hxx](include/wcns_v2/init/flow_initializer.hxx)

**类**: `FlowInitializer`

**原则**: **仅填充内点**（cell 索引 `[ng, ng+nci_core-1]` 范围）

**主入口** `initialize()`:
1. 根据 `cfg.init_type` 分派到具体方法
2. 调用 `lb.field.prim_to_cons(cfg.gamma)` 同步守恒量

#### 4.9.1 均匀流 `init_uniform()`

```
对所有内点 (i,j,k):
    prim.rho = 1.0
    prim.u   = u_inf  (由 Mach + AoA + beta 确定)
    prim.v   = v_inf
    prim.w   = w_inf
    prim.p   = 1.0 / (gamma · Mach²)   (无量纲来流静压)
```

#### 4.9.2 泊肃叶流 `init_poiseuille()`

**物理背景**: 两无限大平板间的压力驱动层流。理论解为抛物线速度剖面。

**流速方向判定**: 比较 `body_force_x/y/z` 的绝对值，取最大分量方向为流向。

```
若 |bfx| 最大 → 流向为 x，剖面在 y 方向:
    yc  = cell_y(i,j,k)             — 单元中心 y 坐标
    yn  = (yc - y_min) / (y_max - y_min)  — 归一化 [0, 1]
    u   = 4 · umax · yn · (1 - yn)  — 抛物线剖面
    v   = 0, w = 0

若 |bfy| 最大 → 流向为 y，剖面在 x 方向
若 |bfz| 最大 → 流向为 z，剖面在 x 方向
```

其他量（ρ, p, 非流向速度分量）与均匀流一致。

### 4.10 边界条件应用

**文件**: [bc_applier.h](include/wcns_v2/bc/bc_applier.h) / [.hxx](include/wcns_v2/bc/bc_applier.hxx)

**类**: `BoundaryConditionApplier`

**设计**: 仅处理 `neighbors[f].active == false` 的面（BC 面），并行邻居面由 HaloExchange 填充。

**三阶段流程**:

```
Stage 1: apply_face_ghost  → 6 个面 ghost 填充（含边/角区域）
Stage 2: apply_edge_ghost  → 12 条边 ghost 修正（两面平均）
Stage 3: apply_corner_ghost → 8 个角 ghost 修正（三面平均）
```

#### Stage 1: 面 Ghost 填充

1. 从 BCPatch 列表确定每面的 BC 类型（取第一个匹配的 patch 类型）
2. 对**整个面的切向范围**（而非 patch 范围）应用 BC，因为边/角修正需要完整的面 ghost 数据
3. 对每个 BC 类型调用对应的填充函数
4. 自周期面（target_rank==-1）：从对面内点逐层拷贝
   - IMIN ghost(d) ← IMAX interior(nci-1-ng-d)
   - 类似处理 J/K 面

#### Stage 2: 边 Ghost 修正

12 条边 = I 面 ghost × J 面 ghost 的 4 组合 + I 面 ghost × K 面 ghost 的 4 组合 + J 面 ghost × K 面 ghost 的 4 组合。

对每条边:
- 计算两面 ghost 范围的交集
- 重叠区域 cell 值:
  - 两面都是 BC → 取平均
  - 仅一面是 BC → 取该面值

#### Stage 3: 角 Ghost 修正

8 个角 = I 面 × J 面 × K 面的所有组合（每面可取 MIN 或 MAX）。

对每个角:
- 计算三面 ghost 范围的交集
- 收集所有 BC 面的值取平均

#### 各 BC 类型算法

##### Farfield（远场 / 特征边界条件）

基于 1D Riemann 不变量理论:

```
对 ghost 层 d = 0..ng-1:
    取对应内点值 (rho_int, u_int, v_int, w_int, p_int)
    c_int = sqrt(γ · p_int / ρ_int)
    c_inf = sqrt(γ · p_inf / rho_inf)

    // 面法向速度（正 = 出域）
    u_n_int = sign · u_int  (sign = +1 for MIN, -1 for MAX)

    // Riemann 不变量
    R⁺ = u_n_int + 2·c_int/(γ-1)       (出射)
    R⁻ = u_n_inf - 2·c_inf/(γ-1)       (入射)

    // 边界值
    u_n_bnd = ½(R⁺ + R⁻)
    c_bnd   = ¼(γ-1)(R⁺ - R⁻)

    // 熵（依据流动方向）
    若 u_n_bnd > 0 (出流): s = p_int / ρ_int^γ
                          v_t_bnd = v_t_int
    若 u_n_bnd < 0 (入流): s = p_inf / ρ_inf^γ
                          v_t_bnd = v_t_inf

    ρ_bnd = (c_bnd² / (γ·s))^(1/(γ-1))
    p_bnd = s · ρ_bnd^γ
```

##### Wall (no-slip, adiabatic)

均匀反射（ghost cell 与 interior cell 关于壁面对称）:

```
ghost(d)     ← mirror = interior(ng + (ng-1-d))   [IMIN]
ghost(gi)    ← mirror = interior(nci-1-ng - (ng-1-d))  [IMAX]

ρ_ghost = ρ_mirror      (零梯度)
p_ghost = p_mirror      (零梯度 / 绝热)
u_ghost = -u_mirror     (无滑移：全反射)
v_ghost = -v_mirror
w_ghost = -w_mirror
```

##### Wall (slip)

与 Symmetry 相同：仅法向速度分量反射。

##### Symmetry / Slip Wall

```
ρ_ghost = ρ_mirror
p_ghost = p_mirror
法向速度: u_n_ghost = -u_n_mirror
切向速度: u_t_ghost = u_t_mirror
```

##### Inflow

所有 ghost cell 赋来流值:
```
ρ_ghost = 1.0
u/v/w_ghost = 自由来流速度
p_ghost = 1/(γ·Mach²)
```

##### Outflow

零阶外推（所有 ghost 层取第一层内点值）:
```
prim_ghost(d) = prim_interior(ng)  [IMIN]
```

### 4.11 Halo 交换执行

**方法**: `ParallelManager::exchange_all_halos()`

对每个本块的 5 个原始变量数组（rho, u, v, w, p）进行 MPI 通信:
- 向远程邻居发送内点数据
- 接收邻居内点数据填充本块 ghost

**注意**: 自周期面（同进程 target_rank==-1）的 ghost 已在 BC 阶段通过自拷贝完成；`copy_local()` 目前为 no-op（从完整 zone 提取的 ghost 已含正确数据）。

---

## 5. 数据结构总览

### 类型别名

| 别名 | 基础类型 | 说明 |
|------|----------|------|
| `Real` | `double` | 浮点精度 |
| `Int` | `int` | 整数类型 |

### 核心类关系

```
Config                  — 物理参数/控制参数（值对象）
    │
Grid                    — 结构化网格（节点坐标、单元中心、度量系数）
    ├── node_x/y/z      — 顶点坐标 (ni×nj×nk)
    ├── cell_x/y/z      — 单元中心坐标 (nci×ncj×nck)
    ├── met_xi_x..met_zeta_z — SCMM 度量项 (nci×ncj×nck)
    ├── face_xi_x..face_zeta_z — 面度量
    ├── jacobian        — Jacobian (nci×ncj×nck)
    ├── bc: BoundaryCondition — BCPatch 列表
    └── connections: ConnectivityList — 1-to-1 连接列表
    │
Field                   — 流场变量
    ├── prim: PrimitiveVars    — (rho, u, v, w, p) × (nci×ncj×nck)
    ├── cons: ConservativeVars — (rho, rhou, rhov, rhow, rhoE)
    ├── inv_xi/eta/zeta: FluxVars — 无粘通量（面中心）
    ├── vis_xi/eta/zeta: FluxVars — 粘性通量（面中心）
    └── rhs: PrimitiveVars     — 右端项
    │
LocalBlock              — 本进程块（Grid + Field + 邻居信息）
    ├── grid: Grid
    ├── field: Field
    └── neighbors[6]: NeighborInfo
    │
ParallelManager         — 顶层并行管理器
    ├── initialize()    — 驱动整个前处理流程
    ├── exchange_all_halos() — MPI 通信
    └── global_max/min/sum() — 全局归约
```

### BCType 枚举

| 值 | 对应 CGNS 类型 | 说明 |
|----|---------------|------|
| `Inflow` | BCInflow, etc. | 入流 |
| `Outflow` | BCOutflow, etc. | 出流 |
| `Wall` | BCWall, BCWallInviscid | 壁面 |
| `Farfield` | BCFarfield | 远场 |
| `Symmetry` | BCSymmetryPlane, BCSymmetryPolar | 对称面 |
| `Periodic` | — (来自 1-to-1 连接) | 周期 |
| `Unknown` | 其他 | 未知（退化外推） |

---

## 6. 网格索引布局

详见记忆文件: [grid-index-layout.md](https://)

以 `ng=3`, i 方向为例：

```
Cell 索引 (0-based):
┌──────────┬───────────────────────────┬──────────┐
│ ghost    │       interior             │ ghost    │
│ [0, 2]   │       [3, nx+2]           │ [nx+3,   │
│          │                            │  nx+5]   │
└──────────┴───────────────────────────┴──────────┘
  3 layers      nci_core cells            3 layers

Node 索引 (0-based):
┌──────────┬───────────────────────────────┬──────────┐
│ ghost    │        core nodes             │ ghost    │
│ [0, 2]   │        [3, nx+3]             │ [nx+4,   │
│          │                               │  nx+6]   │
└──────────┴───────────────────────────────┴──────────┘
```

**维数关系**:
```
ni = ni_core + 2·ng          (总节点数)
nci = ni - 1 = nci_core + 2·ng  (总单元数)
nci_core = ni_core - 1        (核心单元数)

内点范围: [ng, ng + nci_core - 1]  即 [3, nx+2] (当 ng=3)
虚拟网格: [0, ng-1] 和 [ng+nci_core, nci-1]
```

**关键结论**:
1. 物理量存储在单元中心（cell center），不在网格节点上
2. 物理边界在网格节点上，不在单元中心上
3. Ghost cells 在物理域外，用于支持格式模板和边界条件
4. 边界条件通过 ghost cell 值影响通量计算，而非直接修正边界上的物理量

---

## 7. 当前功能完备性检查

### 已完成模块

| # | 模块 | 状态 | 文件 |
|---|------|------|------|
| 1 | MPI 环境初始化 | ✅ 完成 | parallel_env.h |
| 2 | 配置文件读取 | ✅ 完成 | config_reader.h/.cpp |
| 3 | CGNS 网格读取 | ✅ 完成 | cgns_reader.h/.cpp |
| 4 | Ghost 层扩展 | ✅ 完成 | grid.hxx |
| 5 | 网格度量系数 (SCMM) | ✅ 完成 | grid.hxx |
| 6 | 面度量插值 | ✅ 完成 | grid.hxx |
| 7 | 区域分解 | ✅ 完成 | decomposer.h |
| 8 | LocalBlock 构建 | ✅ 完成 | local_block.h |
| 9 | Halo 交换建立 | ✅ 完成 | halo_exchange.h |
| 10 | 流场初始化 | ✅ 完成 | flow_initializer.h |
| 11 | 边界条件应用 | ✅ 完成 | bc_applier.h |
| 12 | BC 类型: Farfield | ✅ 完成 | Riemann 不变量 |
| 13 | BC 类型: Wall (no-slip) | ✅ 完成 | 偶反射 |
| 14 | BC 类型: Wall (slip) | ✅ 完成 | 法向反射 |
| 15 | BC 类型: Symmetry | ✅ 完成 | 法向反射 |
| 16 | BC 类型: Inflow | ✅ 完成 | 指定值 |
| 17 | BC 类型: Outflow | ✅ 完成 | 零阶外推 |
| 18 | Ghost 边/角修正 | ✅ 完成 | 平均 |
| 19 | 自周期面填充 | ✅ 完成 | 对面内点拷贝 |
| 20 | Halo MPI 交换 | ✅ 完成 | 阻塞式 |
| 21 | prim ⇄ cons 转换 | ✅ 完成 | field.hxx |
| 22 | GCL 验证 | ✅ 完成 | main.cpp |
| 23 | 初始化验证 | ✅ 完成 | main.cpp |

### 已验证的测试场景

| 场景 | 进程 | 结果 |
|------|------|------|
| 单 zone 均匀流 | 1 | GCL ~10⁻¹⁷, ghost ρ∈[1,1], corner ρ=1 |
| 单 zone 拆分均匀流 | 2 | GCL ~10⁻¹⁶, ghost ρ∈[1,1] |
| 双 zone 均匀流 | 2 | GCL ~10⁻¹⁷, ghost ρ∈[1,1] |
| 单 zone 泊肃叶流 | 1 | GCL ~10⁻¹⁷, ghost ρ∈[1,1] |

### 已知限制

| 限制 | 说明 |
|------|------|
| 周期连接 transform | 目前仅支持 transform={1,2,3}（无方向反转） |
| copy_local() | 同进程块间直接拷贝为 no-op（依赖初始提取） |
| 单元体积 | 当前为均匀直角近似，非一般曲线网格的真实体积 |
| non-blocking halo | start_exchange/wait_exchange 的 unpack 未完成 |
| 无 CGNS 写出 | 仅能读取网格，不能输出 CGNS 格式结果 |

---

## 8. 后续待实现模块

按自然的求解器开发顺序:

| 顺序 | 模块 | 说明 |
|------|------|------|
| 1 | **WCNS 非线性插值** | 特征变量到半节点的 WCNS 插值（5/6 阶） |
| 2 | **Riemann 求解器** | 面通量计算（Roe, AUSM+, HLLC 等） |
| 3 | **粘性通量** | 速度/温度梯度 + 应力张量 + 热通量 |
| 4 | **时间推进** | RK3-TVD / RK4 / LU-SGS 的实现 |
| 5 | **源项** | 体积力（body force）添加到动量方程 |
| 6 | **结果输出** | Tecplot/VTK/Plot3D 格式写出 |
| 7 | **残差监控** | 收敛历史、L2 范数全局归约 |
| 8 | **重启文件** | 读写中间状态 |
| 9 | **湍流模型** | SA / k-ω SST 等 |
