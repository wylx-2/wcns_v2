# MPI 跨进程通讯详细文档

## 目录

1. [概述与通讯总览](#一概述与通讯总览)
2. [MPI 环境管理](#二mpi-环境管理)
3. [Cell 中心量 Halo 交换 (HaloExchange)](#三cell-中心量-halo-交换-haloexchange)
4. [面通量 Halo 交换 (FluxHaloExchange)](#四面通量-halo-交换-fluxhaloexchange)
5. [全局归约操作 (ParallelManager)](#五全局归约操作-parallelmanager)
6. [History Monitor 归约](#六history-monitor-归约)
7. [Solution Writer 收集](#七solution-writer-收集)
8. [Restart Writer 通讯](#八restart-writer-通讯)
9. [主时间循环通讯时序](#九主时间循环通讯时序)
10. [初始化通讯时序](#十初始化通讯时序)
11. [通讯数据汇总表](#十一通讯数据汇总表)
12. [已知问题与注意事项](#十二已知问题与注意事项)

---

## 一、概述与通讯总览

本求解器使用 MPI 进行多进程并行计算。通讯分为三大类：

| 类别 | 模块 | 通讯模式 | 频率 |
|------|------|----------|------|
| **运行时数据交换** | `HaloExchange` | 点对点 (P2P) | 每 RK 阶段 |
| **运行时数据交换** | `FluxHaloExchange` | 点对点 (P2P) | 每 RK 阶段 |
| **全局归约** | `ParallelManager` | 集体 (Collective) | 每时间步 |
| **全局归约** | `HistoryMonitor` | 集体 (Collective) | 每 `residual_freq` 步 |
| **I/O 收集** | `SolutionWriter` | 集体 (Gather) | 每 `output_freq` 步 |
| **I/O 收集** | `RestartWriter` | 集体 (Gather/Bcast) | 每 `restart_freq` 步 / 启动时 |

### 关键索引约定

- **ng = 3**: 虚拟网格层数
- **Cell 索引**: `[0..ng-1]` 虚拟, `[ng..nci-ng-1]` 内部, `[nci-ng..nci-1]` 虚拟
- **Face 索引**: `[0..nci]`，face h 位于 cell h-1 与 cell h 之间
- **六面编号**: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX

### 通讯协议

所有点对点通讯采用**三段式（3-Phase）无死锁协议**：

1. **Phase 1**: 先 Post 所有 `MPI_Irecv`
2. **Phase 2**: 打包数据 + `MPI_Isend`
3. **Phase 3**: `MPI_Waitall` + 解包

此顺序保证即使周期边界跨 face 索引（如 rank0 IMIN ↔ rank1 IMAX）也不会死锁。

---

## 二、MPI 环境管理

**文件**: [include/parallel/parallel_env.hxx](include/parallel/parallel_env.hxx)

| 调用 | MPI 函数 | 时机 | 说明 |
|------|----------|------|------|
| 初始化 | `MPI_Init_thread(THREAD_MULTIPLE)` | 程序启动 | 请求多线程支持 |
| 获取 rank | `MPI_Comm_rank` | 初始化 | 存入 `rank_` |
| 获取 size | `MPI_Comm_size` | 初始化 | 存入 `size_` |
|  Barrier | `MPI_Barrier` | 初始化打印时 | 同步终端输出 |
| 终结 | `MPI_Finalize` | 程序退出 | |

- 通讯域统一使用 `MPI_COMM_WORLD`
- 单进程模式 (`size_ == 1`) 时，所有归约操作跳过 MPI 调用直接返回本地值

---

## 三、Cell 中心量 Halo 交换 (HaloExchange)

**文件**: [include/parallel/halo_exchange.h](include/parallel/halo_exchange.h), [include/parallel/halo_exchange.hxx](include/parallel/halo_exchange.hxx)

这是最核心的运行时通讯模块，负责在相邻 block 之间交换 cell 中心的 ghost 层数据。

### 3.1 数据结构

```cpp
struct FaceBuffer {
    vector<Real> send_buf;     // 发送缓冲区
    vector<Real> recv_buf;     // 接收缓冲区
    int  target_rank;          // 目标 MPI rank
    int  target_block;         // 目标全局 block ID
    bool is_remote;            // 是否为跨进程通讯
    bool active;               // 该面是否激活
    MPI_Request send_req;      // 非阻塞发送请求
    MPI_Request recv_req;      // 非阻塞接收请求
    Int ni_send, nj_send, nk_send; // 发送/接收子体积维度
};
```

每个 block 有 6 个 `FaceBuffer`（每面一个）。

### 3.2 `is_remote` 判定

**文件**: [include/parallel/halo_exchange.hxx:29-30](include/parallel/halo_exchange.hxx#L29-L30)

```cpp
fb.is_remote = ni.active && ni.target_rank >= 0 && ni.target_rank != -1;
```

- `target_rank == -1` → 同一 block 自周期（如 KMIN↔KMAX 周期边界），不通过 MPI
- `target_rank == my_rank` → 同一进程内的不同 block（内部剖分边界），**当前仍走 MPI 自发送**
- `target_rank != my_rank` → 真正的跨进程通讯

**⚠️ 已知 Bug 1**: 同进程多 block 使用 `target_rank == my_rank`，`is_remote` 为 `true`，导致 MPI 自发送。`exchange_multi` 中 Phase 2 的 `!fb.is_remote` 分支不被执行，而是走远程分支进行 MPI 通讯——这在本地上是可行的（MPI 会处理自消息），但低效且绕过了更直接的内存拷贝路径。

### 3.3 缓冲区大小

| 面方向 | 缓冲区大小 (每个数组) | 说明 |
|--------|----------------------|------|
| IMIN / IMAX | `ng × ncj × nck` | I-face，变维 i 方向取 ng 层 |
| JMIN / JMAX | `nci × ng × nck` | J-face，变维 j 方向取 ng 层 |
| KMIN / KMAX | `nci × ncj × ng` | K-face，变维 k 方向取 ng 层 |

多数组时总大小 = 单数组大小 × 数组数量。

### 3.4 pack_face — 发送数据打包

**位置**: [include/parallel/halo_exchange.hxx:80-134](include/parallel/halo_exchange.hxx#L80-L134)

从当前 block 的**内部区域**复制 ng 层数据到发送缓冲区：

| 面 | 发送源索引范围 | 说明 |
|----|---------------|------|
| IMIN (face=0) | `i = [ng, 2*ng-1]` | 内部区域靠近 IMIN 的 ng 层 |
| IMAX (face=1) | `i = [nci-2*ng, nci-ng-1]` | 内部区域靠近 IMAX 的 ng 层 |
| JMIN (face=2) | `j = [ng, 2*ng-1]` | 内部区域靠近 JMIN 的 ng 层 |
| JMAX (face=3) | `j = [ncj-2*ng, ncj-ng-1]` | 内部区域靠近 JMAX 的 ng 层 |
| KMIN (face=4) | `k = [ng, 2*ng-1]` | 内部区域靠近 KMIN 的 ng 层 |
| KMAX (face=5) | `k = [nck-2*ng, nck-ng-1]` | 内部区域靠近 KMAX 的 ng 层 |

**发送数据的位置语义**:
- **IMIN 面**: 发送 `i ∈ [ng, 2*ng-1]` → 邻居接收后放入 `i ∈ [nci_nbr-ng, nci_nbr-1]` (邻居的 IMAX ghost)
- **IMAX 面**: 发送 `i ∈ [nci-2*ng, nci-ng-1]` → 邻居接收后放入 `i ∈ [0, ng-1]` (邻居的 IMIN ghost)

### 3.5 unpack_face — 接收数据解包

**位置**: [include/parallel/halo_exchange.hxx:140-195](include/parallel/halo_exchange.hxx#L140-L195)

从接收缓冲区复制到当前 block 的 **ghost 区域**：

| 面 | 接收目标索引范围 | 说明 |
|----|-----------------|------|
| IMIN (face=0) | `i = [0, ng-1]` | IMIN 侧 ghost 层 |
| IMAX (face=1) | `i = [nci-ng, nci-1]` | IMAX 侧 ghost 层 |
| JMIN (face=2) | `j = [0, ng-1]` | JMIN 侧 ghost 层 |
| JMAX (face=3) | `j = [ncj-ng, ncj-1]` | JMAX 侧 ghost 层 |
| KMIN (face=4) | `k = [0, ng-1]` | KMIN 侧 ghost 层 |
| KMAX (face=5) | `k = [nck-ng, nck-1]` | KMAX 侧 ghost 层 |

### 3.6 位置一致性验证

**发送方 (IMIN 面为例)**: 发送 `[ng, 2*ng-1]` — 这是本 block 紧邻 IMIN ghost 区的内部 cell。
**接收方 (IMAX 面)**: 接收后填入 `[nci-ng, nci-1]` — IMAX ghost 区。

这两个区域在物理空间上对应**同一个位置**：block 边界面两侧的 cell。发送方发送其内部值，接收方用它填充 ghost。对于标准邻居（face A 的 IMIN ↔ face B 的 IMAX），位置严格匹配。

**同理**:
- 发送方 IMAX `[nci-2*ng, nci-ng-1]` ↔ 接收方 IMIN `[0, ng-1]`
- 发送方 JMIN `[ng, 2*ng-1]` ↔ 接收方 JMAX `[ncj-ng, ncj-1]`
- 等等

**结论**: pack/unpack 的索引位置是**严格对称**的，发送方和接收方的空间位置一致。

### 3.7 exchange_multi — 阻塞式多数组交换

**位置**: [include/parallel/halo_exchange.hxx:225-345](include/parallel/halo_exchange.hxx#L225-L345)

**签名**: `exchange_multi(arrays, block, all_blocks)` — `all_blocks` 为本进程所有 block 的列表，用于同进程多 block 间的直接内存拷贝。

流程：
1. **Phase 1 (行 246-266)**: 遍历 6 面，跳过非活跃或非远程面，为每个活跃远程面 Post `MPI_Irecv`
2. **Phase 2 (行 269-325)**:
   - 若 `!fb.is_remote`：
     - 同 block 自周期 (`target_block == -1 || target_block == block_id`)：`pack_face` 从自身 `target_face` 侧读取 + `unpack_face` 写入自身 `face` 侧 ghost
     - 同进程不同 block：从 `all_blocks` 中查找邻居，调用 `copy_local` 直接从邻居 block 的内部区域拷贝
   - 若远程：打包到 `send_buf`，然后 `MPI_Isend`
3. **Phase 3 (行 328-345)**: `MPI_Waitall` 等待所有通讯完成 → 对每个远程面调用 `unpack_face`

**Tag 计算** (行 258-261):
```cpp
tag = min(block_id, target_block) * 1000
    + min(face, target_face) * 10
    + (is_periodic ? 5 : 0)
    + 1
```

- `min(block_id, target_block)` 保证对称：无论 A→B 还是 B→A，发送方和接收方算出的 tag 相同
- `min(face, target_face)` 同理保证对称
- 周期标记参与 tag 以确保不会与非周期通讯混淆

### 3.8 copy_local — 同进程本地拷贝 **[已实现]**

**位置**: [include/parallel/halo_exchange.hxx:205-222](include/parallel/halo_exchange.hxx#L205-L222)

直接从邻居 block 的内部区域拷贝数据到当前 block 的 ghost 区域：

```cpp
void copy_local(arr, face, block, neighbor) {
    pack_face(arr, ni.target_face, neighbor, ...);  // 从邻居的 target_face 侧内部读取
    unpack_face(arr, face, block, ...);              // 写入当前 block 的 face 侧 ghost
}
```

- 邻居打包面 = `ni.target_face`（邻居与当前 block 连接的面）
- 当前解包面 = `face`（当前 block 需要 ghost 数据的面）
- 这两个面对应物理空间中的同一边界，数据位置一致

### 3.9 start_exchange / wait_exchange — 非阻塞交换

**位置**: [include/parallel/halo_exchange.hxx:336-399](include/parallel/halo_exchange.hxx#L336-L399)

- `start_exchange`: 为每个活跃面打包数据并启动 `MPI_Isend`/`MPI_Irecv`
- `wait_exchange`: 调用 `MPI_Wait` 等待完成
- **注意**: 非阻塞模式不处理同进程邻居（行 354-358 直接 `continue`），且 `wait_exchange` 不执行 `unpack`（行 395-398 注释）
- 当前代码仅使用阻塞式 `exchange_multi`，非阻塞接口未使用

### 3.10 调用的数据

通过 `exchange_multi` 交换的数据（按调用频率排序）：

| 调用方法 | 数组数量 | 具体数据 | 每 RK 阶段? |
|----------|---------|---------|-------------|
| `exchange_all_halos` | 5 | `prim.rho, prim.u, prim.v, prim.w, prim.p` | ✅ |
| `exchange_gradient_halos` | 12 | `du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dT_dx, dT_dy, dT_dz` | ✅ (若有粘性) |
| `exchange_viscous_flux_halos` | 15 | `vis_x.f1..f5, vis_y.f1..f5, vis_z.f1..f5` | ✅ (若有粘性) |

---

## 四、面通量 Halo 交换 (FluxHaloExchange)

**文件**: [include/parallel/flux_halo_exchange.h](include/parallel/flux_halo_exchange.h), [include/parallel/flux_halo_exchange.hxx](include/parallel/flux_halo_exchange.hxx)

### 4.1 目的与时机

在 Riemann 求解器计算完无粘面通量 (`inv_xi/eta/zeta`) 后，相邻 block 共享的 connectivity 边界处的 ghost face 通量必须与邻居的内部 face 通量一致。这是因为 6 点中心差分格式需要 ±3 个 face 的通量值，而 cell 中心的 ghost 交换不能保证 face 通量的正确性。

**时机**: Riemann 求解器之后，InviscidRHS 差分之前。

**交换量**: ng+1 = 4 个 face 切片（比 cell ghost 多 1，因为差分需要更多 face）

### 4.2 数据结构

```cpp
struct FluxFaceInfo {
    bool active;           // 是否激活
    bool is_remote;        // 是否跨进程
    int  target_rank;      // 目标 rank
    int  target_block;     // 目标全局 block ID
    Int send_begin, send_end;  // 发送 face 范围
    Int recv_begin, recv_end;  // 接收 face 范围
    Int n_faces;           // face 数量 (= ng+1 = 4)
    Int dim1, dim2;        // 单个 face 切片的两个横向维度
    vector<Real> send_buf; // 发送缓冲 (5分量 × n_faces × dim1 × dim2)
    vector<Real> recv_buf; // 接收缓冲
};
```

### 4.3 setup — 连通性解析

**位置**: [include/parallel/flux_halo_exchange.hxx:24-106](include/parallel/flux_halo_exchange.hxx#L24-L106)

通过 `Grid::find_face_connection(face)` 查找 1-to-1 连接，再通过 `zone_to_block` 映射解析目标 block。对于块内周期连接（同 block、同 zone 的内部周期），不激活 FluxHaloExchange。

**Face 索引范围** (以 JMAX ↔ JMIN 为例，ncj 个 cell):

| | 发送范围 | 接收范围 |
|----|---------|---------|
| Block A (JMAX) | `[ncj-2*ng, ncj-ng]` (内部侧) | `[ncj-ng, ncj]` (ghost 侧) |
| Block B (JMIN) | `[ng, 2*ng]` (内部侧) | `[0, ng]` (ghost 侧) |

对应关系：
- A 发送 `[ncj-2*ng, ncj-ng]` → B 接收 `[0, ng]`
- B 发送 `[ng, 2*ng]` → A 接收 `[ncj-ng, ncj]`

### 4.4 pack_flux — 打包面通量

**位置**: [include/parallel/flux_halo_exchange.hxx:112-213](include/parallel/flux_halo_exchange.hxx#L112-L213)

从 `FluxVars` (f1..f5) 中提取 `send_begin` 到 `send_end` 的 face 切片。每个分量独立打包，5 个分量顺序排列。

### 4.5 unpack_flux — 解包面通量

**位置**: [include/parallel/flux_halo_exchange.hxx:219-313](include/parallel/flux_halo_exchange.hxx#L219-L313)

将接收缓冲写入 `recv_begin` 到 `recv_end` 的 face 位置。

### 4.6 copy_local_flux — 同进程面通量拷贝

**位置**: [include/parallel/flux_halo_exchange.hxx:319-403](include/parallel/flux_halo_exchange.hxx#L319-L403)

直接访问邻居 block 的 `FluxVars` 进行内存拷贝：

- 接收方 (MIN face): 邻居是 MAX → 从邻居 `[ndim-2*ng-1, ndim-ng-1]` 拷贝
- 接收方 (MAX face): 邻居是 MIN → 从邻居 `[ng, 2*ng]` 拷贝

注意 MAX→MIN 时使用了 `-1` 修正 (`nbr_s0 = nbr_ndim - 2*ng - 1`)，因为 face 索引的特殊性（face 0 对应 ghost 最外侧）。

### 4.7 exchange — 无粘通量交换

**位置**: [include/parallel/flux_halo_exchange.hxx:409-518](include/parallel/flux_halo_exchange.hxx#L409-L518)

- 通过 `all_blocks` 查找邻居是否在同进程，决定走本地拷贝还是 MPI
- 自周期 block 跳过（`neighbor_ptr->block_id != block.block_id`）——因为 Riemann 求解器已利用 BC applier 的 ghost 数据计算了正确的通量
- **Tag 计算** (行 446):
  ```cpp
  tag = min(block_id, target_block) * 100
      + min(face, face^1) * 10 + dir + 200
  ```
  - `face^1` 翻转 MIN↔MAX（0↔1, 2↔3, 4↔5），保证发送方/接收方算出相同 tag
  - `+200` 偏移与 cell halo 区分

**⚠️ 已知 Bug 3**: Tag 使用 `face^1` 而非 `target_face`。对于标准几何对边连接（IMIN↔IMAX），`face^1 == target_face`，tag 正确。但对非标准连接（通过 `transform` 数组指定的），`face^1` 可能不等于 `target_face`，导致 tag 不匹配。

### 4.8 exchange_face_arrays — 单分量面数组交换

**位置**: [include/parallel/flux_halo_exchange.hxx:524-646](include/parallel/flux_halo_exchange.hxx#L524-L646)

用于粘性通量管道中的面插值数组交换。一次处理一个方向（dir=0,1,2），交换 N 个单分量面数组。

- buffer 布局: N 个分量顺序排列，每个分量 `n_faces × dim1 × dim2`
- Tag: `min(block_id, target_block)*100 + min(face, face^1)*10 + dir + 300`
- 仅处理远程 MPI 邻居，不处理同进程邻居

**调用场景**（在 `ViscidRHS` 中，对每个 dir）:
1. 面插值后的原始变量 (`u_face, v_face, w_face, T_face` 等)
2. 面插值后的 Cartesian 粘性通量

---

## 五、全局归约操作 (ParallelManager)

**文件**: [include/parallel/parallel_manager.hxx](include/parallel/parallel_manager.hxx)

| 方法 | MPI 调用 | 操作 | 用途 |
|------|----------|------|------|
| `global_max(local_val)` | `MPI_Allreduce(MPI_MAX)` | 全局最大值 | 速度/温度/Mach 数极值 |
| `global_min(local_val)` | `MPI_Allreduce(MPI_MIN)` | 全局最小值 | 时间步长 dt，密度/温度最小值 |
| `global_sum(local_val)` | `MPI_Allreduce(MPI_SUM)` | 全局求和 | 残差平方和，cell 总数 |
| `broadcast(val)` | `MPI_Bcast(from=0)` | 广播 | 将 rank 0 的值广播到所有 rank |

所有归约操作在单进程模式下跳过 MPI 直接返回本地值。

### 调用时机与数据

| 调用位置 | 方法 | 数据 | 频率 |
|----------|------|------|------|
| [main.cpp:322](src/main.cpp#L322) | `global_min` | `dt` (时间步长) | 每迭代 |
| `Residual::compute` | `global_sum` | 残差平方和 (5 分量) + cell 数 | 每 `residual_freq` 步 |
| `Residual::monitor` | `global_min/max` | rho, T, p, u, v, w, Mach 极值 | 每 `residual_freq` 步 |
| `HistoryMonitor` | `global_min` | `dx` (网格间距估计) | 首次调用 |

---

## 六、History Monitor 归约

**文件**: [include/core/history_monitor.hxx](include/core/history_monitor.hxx)

**位置**: [include/core/history_monitor.hxx:112-117](include/core/history_monitor.hxx#L112-L117)

```cpp
MPI_Allreduce(local_sum_u.data(), global_sum_u.data(), n_targets, MPI_DOUBLE, MPI_SUM, comm);
MPI_Allreduce(local_vol.data(),  global_vol.data(),  n_targets, MPI_DOUBLE, MPI_SUM, comm);
```

- **时机**: 每次 residual 输出时
- **数据**: `n_targets` 个截面位置的体积加权速度积分和总体积
- **操作**: `MPI_SUM` 归约→ rank 0 计算加权平均
- **结果**: 仅 rank 0 写入 `history.dat`

---

## 七、Solution Writer 收集

**文件**: [include/io/solution_writer.hxx](include/io/solution_writer.hxx)

### 并行模式下的 MPI 收集（写到 rank 0）

| 步骤 | MPI 调用 | 数据 | 类型 |
|------|----------|------|------|
| 1 | `MPI_Gather` | 每个 rank 的 block 数量 | `MPI_INT` |
| 2 | `MPI_Gatherv` | 每个 block 的元数据 `{block_id, ni_core, nj_core, nk_core, data_size}` | `MPI_INT` |
| 3 | `MPI_Gather` | 每个 rank 的总数据量 | `MPI_INT` |
| 4 | `MPI_Gatherv` | 结构化解数据 `{x,y,z,rho,u,v,w,p,T}` 每内部 cell | `MPI_DOUBLE` |

- rank 0 收集后按 `block_id` 排序，写入 Tecplot 文件
- 其他 rank 在 `MPI_Gatherv` 后直接返回

---

## 八、Restart Writer 通讯

**文件**: [include/io/restart_writer.hxx](include/io/restart_writer.hxx)

### 写入

每个 rank 独立写自己的 `restart_XXXXXX_rN.bin` 文件（并行 I/O，**不通过 MPI 传输解数据**）。

| MPI 调用 | 数据 | 用途 |
|----------|------|------|
| `MPI_Gather` | block 数量 | 收集到 rank 0 写索引文件 |
| `MPI_Gatherv` | block ID 列表 | 收集到 rank 0 写索引文件 |

### 读取

每个 rank 独立读自己的 `restart_XXXXXX_rN.bin` 文件（并行 I/O）。

| MPI 调用 | 数据 | 用途 |
|----------|------|------|
| `MPI_Bcast` | `saved_iter` (Int) | 从 rank 0 广播迭代数 |
| `MPI_Bcast` | `saved_time` (Real) | 从 rank 0 广播时间 |

---

## 九、主时间循环通讯时序

每个 RK 阶段的完整通讯序列（[src/main.cpp](src/main.cpp) 行 344-557）：

### 9.1 无粘管道

```
  [已从上一阶段得到: prim ghost 正确, cons ghost 正确]
  │
  ├─ interp_xi/eta/zeta     (WCNS 插值，读取 cons ghost)
  ├─ solve_xi/eta/zeta      (Riemann 求解器，写 inv_xi/eta/zeta)
  │
  ├─ ★ exchange_flux_halos  (FluxHaloExchange::exchange)
  │    ├─ 同进程: copy_local_flux (从邻居 block 内存拷贝)
  │    └─ 跨进程: MPI_Isend/Irecv/Waitall (ng+1 面通量切片)
  │    交换数据: inv_xi/eta/zeta 的 f1..f5 (5分量)
  │
  ├─ InviscidRHS::compute   (6 点中心差分)
```

### 9.2 粘性管道

```
  ├─ interp_to_faces           (cell → face 插值 u,v,w,T)
  ├─ compute_gradients + exchange_and_assemble (面插值量交换)
  │    └─ ★ exchange_face_arrays (FluxHaloExchange)
  │         交换数据: 面插值原始变量 (u_face, v_face, w_face, T_face)
  │
  ├─ ★ exchange_gradient_halos (HaloExchange::exchange_multi)
  │    交换数据: 12 个梯度数组 (du_dx..dT_dz, 每个 nci×ncj×nck)
  │
  ├─ compute_cell_viscous_flux
  │
  ├─ ★ exchange_viscous_flux_halos (HaloExchange::exchange_multi)
  │    交换数据: 15 个粘性通量数组 (vis_x/y/z 各 f1..f5, 每个 nci×ncj×nck)
  │
  ├─ interp_cart_flux_to_faces + exchange_and_assemble_face_flux
  │    └─ ★ exchange_face_arrays (FluxHaloExchange, 每方向)
  │         交换数据: 面插值 Cartesian 粘性通量 (每方向多分量)
  │
  ├─ ViscidRHS::compute_rhs   (6 点中心差分)
```

### 9.3 时间推进后

```
  ├─ cons_to_prim           (更新后的守恒量 → 原始变量)
  ├─ apply_face_ghost       (填充 BC 面 ghost + 自周期 ghost)
  │
  ├─ ★ exchange_all_halos   (HaloExchange::exchange_multi)
  │    交换数据: 5 个原始变量数组 (rho, u, v, w, p)
  │
  ├─ apply_edge_ghost       (修正棱边 ghost，依赖 MPI ghost 数据)
  ├─ apply_corner_ghost     (修正角点 ghost，依赖 MPI ghost 数据)
  ├─ prim_to_cons           (同步 cons ghost)
```

### 9.4 非 RK 通讯（每时间步）

```
  ├─ ★ global_min(dt)       (MPI_Allreduce MIN)
  │
  ├─ (每 residual_freq 步)
  │    ├─ ★ global_sum/max/min  (残差归约)
  │    └─ ★ MPI_Allreduce SUM   (history monitor)
  │
  ├─ (每 output_freq 步)
  │    └─ ★ MPI_Gather/Gatherv  (solution writer)
  │
  └─ (每 restart_freq 步)
       └─ ★ MPI_Gather/Gatherv  (restart writer)
```

---

## 十、初始化通讯时序

### 10.1 首次初始化 (main.cpp:77-210)

```
  ★ MPI_Init_thread         (所有 rank)
  ★ MPI_Comm_rank/size      (所有 rank)
  │
  ├─ (所有 rank 独立读 CGNS, 构建 Grid, 分解)
  │
  ├─ (每个 rank 构建自己的 LocalBlock)
  │    └─ build_neighbors: 判定每个面的邻居 rank/block
  │
  ├─ HaloExchange::setup    (预分配缓冲区, 无 MPI 调用)
  ├─ FluxHaloExchange::setup(预分配缓冲区, 无 MPI 调用)
  │
  ├─ FlowInitializer         (所有 rank 独立初始化流场)
  ├─ apply_face_ghost        (所有 rank, 本地操作)
  │
  ├─ ★ exchange_all_halos    (HaloExchange, 首次 MPI 数据交换)
  │    交换数据: prim.rho, u, v, w, p (5 数组)
  │
  ├─ apply_edge/corner_ghost (所有 rank, 本地操作)
  ├─ prim_to_cons            (所有 rank, 本地操作)
  │
  └─ ★ MPI_Gather/Gatherv   (SolutionWriter, 写到 rank 0)
```

### 10.2 重启动初始化 (main.cpp:122-162)

```
  ├─ (每个 rank 独立读 restart_XXX_rN.bin)
  │
  ├─ ★ MPI_Bcast saved_iter (rank 0 → all)
  ├─ ★ MPI_Bcast saved_time (rank 0 → all)
  │
  ├─ cons_to_prim            (所有 rank)
  ├─ apply_face_ghost        (重新应用 BC)
  │
  ├─ ★ exchange_all_halos    (HaloExchange, 重新填充 MPI ghost)
  │    交换数据: prim.rho, u, v, w, p
  │
  ├─ apply_edge/corner_ghost
  ├─ prim_to_cons
  │
  └─ ★ MPI_Gather/Gatherv   (SolutionWriter, 写重启动状态)
```

---

## 十一、通讯数据汇总表

### 11.1 点对点通讯 (P2P)

| 序号 | 模块 | 方法 | 数据类型 | 每面数据量 | 分量数 | 总数据量(最大) | 频率 |
|------|------|------|----------|-----------|--------|---------------|------|
| 1 | HaloExchange | `exchange_all_halos` | Cell prim | `ng×dim2×dim3` | 5 | `5 × ng × 2*(ncj*nck + nci*nck + nci*ncj)` | 每 RK 阶段 |
| 2 | HaloExchange | `exchange_gradient_halos` | Cell grad | `ng×dim2×dim3` | 12 | `12 × ng × ...` | 每 RK 阶段(粘性) |
| 3 | HaloExchange | `exchange_viscous_flux_halos` | Cell vis flux | `ng×dim2×dim3` | 15 | `15 × ng × ...` | 每 RK 阶段(粘性) |
| 4 | FluxHaloExchange | `exchange` | Face inv flux | `(ng+1)×dim2×dim3` | 5×3方向 | `15 × (ng+1) × ...` | 每 RK 阶段 |
| 5 | FluxHaloExchange | `exchange_face_arrays` | Face misc | `(ng+1)×dim2×dim3` | N | `N × (ng+1) × ...` | 每 RK 阶段(粘性) |

### 11.2 集体通讯 (Collective)

| 序号 | 模块 | MPI 调用 | 数据 | 数据量 | 频率 |
|------|------|----------|------|--------|------|
| C1 | ParallelManager | `MPI_Allreduce(MIN)` | `dt` | 1 Real | 每迭代 |
| C2 | Residual | `MPI_Allreduce(SUM)` | 残差平方和 | 5 Real | 每 `residual_freq` 步 |
| C3 | Residual | `MPI_Allreduce(SUM)` | cell 总数 | 1 Real | 每 `residual_freq` 步 |
| C4 | Residual | `MPI_Allreduce(MIN/MAX)` | 流场极值 | 9 Real | 每 `residual_freq` 步 |
| C5 | HistoryMonitor | `MPI_Allreduce(SUM)` | 加权速度 | N Real | 每 `residual_freq` 步 |
| C6 | SolutionWriter | `MPI_Gather` | block 数量 | 1 Int | 每 `output_freq` 步 |
| C7 | SolutionWriter | `MPI_Gatherv` | block 元数据 | 5×Nblocks Int | 每 `output_freq` 步 |
| C8 | SolutionWriter | `MPI_Gatherv` | 解数据 | 9×Ncells Real | 每 `output_freq` 步 |
| C9 | RestartWriter | `MPI_Gather` | block 数量 | 1 Int | 每 `restart_freq` 步 |
| C10 | RestartWriter | `MPI_Gatherv` | block ID | Nblocks Int | 每 `restart_freq` 步 |
| C11 | RestartWriter | `MPI_Bcast` | iter, time | 1 Int + 1 Real | 启动时 |

### 11.3 生命期通讯

| 序号 | MPI 调用 | 时机 |
|------|----------|------|
| L1 | `MPI_Init_thread` | 程序启动 |
| L2 | `MPI_Comm_rank/size` | 程序启动 |
| L3 | `MPI_Barrier` | 初始化打印（多次） |
| L4 | `MPI_Finalize` | 程序退出 |

---

## 十二、已知问题与注意事项

### 12.1 ✅ Bug 1: HaloExchange 同进程自发送 **[已修复]**

**位置**: [include/parallel/halo_exchange.hxx:29-31](include/parallel/halo_exchange.hxx#L29-L31)

**问题**: 同进程内多 block 的内部剖分边界 `target_rank == my_rank`，`is_remote` 判定为 `true`，走 MPI 自发送路径。

**修复**: 添加 `ni.target_rank != ParallelEnv::rank()` 条件，使同进程 block 走本地内存拷贝路径（更高效且语义更清晰）。

### 12.2 ✅ Bug 2: HaloExchange 同进程内部剖分拷贝数据源错误 **[已修复]**

**位置**: [include/parallel/halo_exchange.hxx:277-312](include/parallel/halo_exchange.hxx#L277-L312)

**问题**: Phase 2 的 `!fb.is_remote` 分支中，对同进程多 block 的拷贝从当前 block 打包数据（使用 `target_face`），而非从邻居 block 打包。

**修复**: 
- 区分两种同进程场景：同 block 自周期 (`target_block == -1`) vs 同进程不同 block
- 同 block 自周期：保持原有行为（从自身对面拷贝）
- 同进程不同 block：通过 `all_blocks` 查找邻居 block，调用 `copy_local` 从邻居的内部区域拷贝
- 实现 `copy_local`：从邻居 block 的 `target_face` 侧打包，写入当前 block 的 `face` 侧 ghost

### 12.3 ✅ Bug 3: FluxHaloExchange Tag **[确认为非问题]**

**位置**: [include/parallel/flux_halo_exchange.hxx:446-447](include/parallel/flux_halo_exchange.hxx#L446-L447)

**分析**: Tag 使用 `face ^ 1`（几何对边）。在 CGNS 结构化网格的 1-to-1 连接中，连接面始终是几何对边（IMIN↔IMAX, JMIN↔JMAX, KMIN↔KMAX），`face^1` 等价于 `target_face`。故 tag 始终对称，不会不匹配。无需修改。

### 12.4 设计注意事项

1. **自周期 block (单 block 周期边界)**: 不通过 MPI 通讯，由 BC Applier 的 `apply_face_ghost` 直接拷贝内部数据到对面 ghost。FluxHaloExchange 也正确跳过（`neighbor_ptr->block_id != block.block_id` 守卫）。

2. **FluxHaloExchange 不处理同进程 face 数组交换**: `exchange_face_arrays` 方法只处理远程 MPI 邻居。对于同进程的邻居（如内部剖分后的邻接 block），粘性通量管道中的 `ViscidRHS::exchange_and_assemble_face_flux` 需要额外的处理逻辑。

3. **非阻塞 HaloExchange 不完整**: `start_exchange`/`wait_exchange` 不支持同进程邻居，且 `wait_exchange` 不执行 `unpack`。当前代码仅使用阻塞式 `exchange_multi`。

4. **所有点对点通讯使用 `MPI_DOUBLE`**: `Real` 类型与 `double` 等价，所有浮点通讯使用 `MPI_DOUBLE`。

5. **所有集体操作使用 `MPI_COMM_WORLD`**: 不使用子通讯域。

6. **Buffer 预分配**: `HaloExchange` 和 `FluxHaloExchange` 在 `setup` 阶段预分配 buffer，避免运行时动态分配。

---

*文档生成时间: 2026-06-11*
*相关文档: [ghost.md](ghost.md)*
