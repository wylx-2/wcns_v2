# Ghost Layer 完整文档

## 概述

本项目使用 `ng = 3` 层 ghost cell 支持 WCNS 6 阶格式（需要 ±3 个模板点）。Ghost 层同时存在于**节点坐标**和**单元中心变量**上，且分为 cell 值和 face 值两大类。

### 索引布局

```
Ghost cells (0..ng-1) | Interior cells (ng..nci-ng-1) | Ghost cells (nci-ng..nci-1)
                       |←——— nci_core cells ———→|
   total: nci = nci_core + 2*ng cells

Faces (0..nci): face h 位于 cell h-1 与 cell h 之间
   Ghost faces: 0..ng 和 nci-ng..nci
   Interior faces: ng..nci-ng
```

---

## 一、节点坐标 Ghost（初始化阶段，仅执行一次）

### 1.1 Grid::extend_ghost_layers()

- **文件**: [include/grid/grid.hxx:80-129](include/grid/grid.hxx#L80-L129)
- **调用**: `ParallelManager::initialize()` 中，在 CGNS 读取后（[parallel_manager.hxx:31](include/parallel/parallel_manager.hxx#L31)）
- **作用**: 
  1. 保存原始 core 节点坐标
  2. 分配 `(ni_core+2*ng) × (nj_core+2*ng) × (nk_core+2*ng)` 扩展数组
  3. 将 core 数据复制到扩展数组的 interior（偏移 ng）
  4. 更新网格维度 (`ni, nj, nk, nci, ncj, nck`)
  5. 调用 `fill_ghost_nodes()` 填充 ghost 节点
  6. 重新分配并计算 cell 中心坐标和体积

### 1.2 Grid::fill_ghost_nodes()

- **文件**: [include/grid/grid.hxx:131-140](include/grid/grid.hxx#L131-L140)
- **调用**: `extend_ghost_layers()` 内部
- **作用**: 遍历 6 个面，对每个面如果有周期连接则调用 `fill_ghost_face_periodic()`，否则调用 `fill_ghost_face_extrapolate()`

### 1.3 Grid::fill_ghost_face_periodic(face, conn)

- **文件**: [include/grid/grid.hxx:341-496](include/grid/grid.hxx#L341-L496)
- **修饰**: **节点坐标 ghost**，周期拷贝
- **作用**: 从同一 zone 的对面 interior 拷贝节点坐标到 ghost 区域，减去周期平移向量。支持全部 6 个面，`transform={1,2,3}` 简单周期情形。

### 1.4 Grid::fill_ghost_face_extrapolate(face)

- **文件**: [include/grid/grid.hxx:498-563](include/grid/grid.hxx#L498-L563)
- **修饰**: **节点坐标 ghost**，线性外推
- **作用**: 对非周期面，逐层向外线性外推：`x(-d) = 2*x(-(d-1)) - x(-(d-2))`

### 1.5 Grid::fix_interface_ghost(face, donor, conn)

- **文件**: [include/grid/grid.hxx:207-339](include/grid/grid.hxx#L207-L339)
- **调用**: `ParallelManager::initialize()` 中，在 `extend_ghost_layers()` 之后（[parallel_manager.hxx:51](include/parallel/parallel_manager.hxx#L51)）
- **修饰**: **节点坐标 ghost**，inter-zone 接口修正
- **作用**: 对于跨 zone 的 1-to-1 连接（非同 zone），用 donor zone 的 interior 节点坐标覆盖之前线性外推的 ghost 节点。处理方向匹配的连接，方向不匹配的暂时跳过（TODO）。

### 1.6 Grid::extract_metrics_from(donor, ci0, cj0, ck0)

- **文件**: [include/grid/grid.hxx:942-1021+](include/grid/grid.hxx#L942)
- **调用**: `LocalBlock::from_sub_zone()` 中（[local_block.hxx:186](include/parallel/local_block.hxx#L186)）
- **修饰**: **cell 中心度量系数 + 面度量系数**（包含 ghost 区域）
- **作用**: 从 full zone 的已计算度量系数中提取子块对应的区域。Full zone 的度量系数是在正确的 ghost 模板支持下计算的，因此子块内部 ghost 区域的度量系数也是正确的。提取范围包括所有 cell-center metric（含 ghost）和 face metric（含 ghost face）。

---

## 二、Cell 中心 Primitive 变量 Ghost（每 RK 阶段执行）

### 2.1 BoundaryConditionApplier::apply_face_ghost(lb, cfg)

- **文件**: [include/bc/bc_applier.hxx:160-279](include/bc/bc_applier.hxx#L160-L279)
- **调用时序**: 初始化阶段 + 每个 RK 阶段的 post-update
- **修饰**: **Cell prim ghost**，面区域（整面，包括边角重叠区）
- **处理两类面**:

#### A. BC 面（物理边界，`neighbors[face].active == false`）
  1. 确定面的 BC 类型（从第一个匹配的 BC patch）
  2. 在整个切向范围 `[0, n_tangential-1]` 上调用 `apply_bc_on_face()`
  3. 如果没有 BC patch，回退到 outflow 外推

#### B. 自周期面（同进程，`target_rank == -1`）
  直接拷贝对面 interior cell 的 prim 值到 ghost cell:
  - IMIN (face 0): ghost `gi = d` ← source `si = nci - 1 - ng - d`
  - IMAX (face 1): ghost `gi = nci - 1 - d` ← source `si = ng + d`
  - JMIN/JMAX/KMIN/KMAX 同理

### 2.2 BoundaryConditionApplier::apply_bc_on_face(face, bc_type, ...)

- **文件**: [include/bc/bc_applier.hxx:285-316](include/bc/bc_applier.hxx#L285-L316)
- **修饰**: **Cell prim ghost**，分派器
- **分派到**:

| BC 类型 | 函数 | 策略 |
|---------|------|------|
| Farfield | `apply_farfield` ([L322-491](include/bc/bc_applier.hxx#L322-L491)) | 1D Riemann 不变量 |
| Wall (slip) | `apply_wall_slip` → `apply_symmetry` ([L679-685](include/bc/bc_applier.hxx#L679-L685)) | 法向速度反射，切向零梯度 |
| Wall (noslip) | `apply_wall_noslip` ([L497-574](include/bc/bc_applier.hxx#L497-L574)) | 所有速度反射 + 偶数镜像 |
| Wall (noslip iso) | `apply_wall_noslip_isothermal` ([L580-673](include/bc/bc_applier.hxx#L580-L673)) | 速度反射 + 等温壁面 T |
| Symmetry | `apply_symmetry` ([L691-779](include/bc/bc_applier.hxx#L691-L779)) | 法向速度反射 |
| Inflow | `apply_inflow` ([L785-838](include/bc/bc_applier.hxx#L785-L838)) | 全部设为自由来流 |
| Outflow | `apply_outflow` ([L844-928](include/bc/bc_applier.hxx#L844-L928)) | 零阶外推 |

所有 BC 函数都使用 ghost mirror 索引:
- IMIN: `ghost(d) ↔ mirror(ng + (ng-1-d))`
- IMAX: `ghost(nci-1-d) ↔ mirror(nci-1-ng-(ng-1-d))`

### 2.3 HaloExchange::exchange_multi(arrays, block)

- **文件**: [include/parallel/halo_exchange.hxx:225-330](include/parallel/halo_exchange.hxx#L225-L330)
- **调用时序**: 每次 `apply_face_ghost` 之后（初始化 + 每个 RK 阶段）
- **修饰**: **Cell prim ghost** + **Cell gradient ghost** + **Cell viscous flux ghost**
- **协议**: 3 阶段无死锁 MPI 交换
  1. **Phase 1**: 对所有 `is_remote == true` 的面 post `MPI_Irecv`
  2. **Phase 2**: 
     - 非 remote 面（自周期/同进程）: `pack_face(对面 interior) → unpack_face(本面 ghost)`
     - Remote 面: `pack_face(本面 interior) → MPI_Isend`
  3. **Phase 3**: `MPI_Waitall` + `unpack_face(recv_buf → ghost)`

- **pack_face 发送范围** (ng 层 interior):
  | 面 | 发送索引 |
  |----|---------|
  | IMIN (0) | `[ng, 2*ng-1]` |
  | IMAX (1) | `[nci-2*ng, nci-ng-1]` |
  | JMIN (2) | `[ng, 2*ng-1]` |
  | JMAX (3) | `[ncj-2*ng, ncj-ng-1]` |
  | KMIN (4) | `[ng, 2*ng-1]` |
  | KMAX (5) | `[nck-2*ng, nck-ng-1]` |

- **unpack_face 接收范围** (ng 层 ghost):
  | 面 | 接收索引 |
  |----|---------|
  | IMIN (0) | `[0, ng-1]` |
  | IMAX (1) | `[nci-ng, nci-1]` |
  | JMIN (2) | `[0, ng-1]` |
  | JMAX (3) | `[ncj-ng, ncj-1]` |
  | KMIN (4) | `[0, ng-1]` |
  | KMAX (5) | `[nck-ng, nck-1]` |

- **is_remote 判断** ([halo_exchange.hxx:29-31](include/parallel/halo_exchange.hxx#L29-L31)):
  ```cpp
  fb.is_remote = ni.active && ni.target_rank >= 0 && ni.target_rank != -1
              && ni.target_rank != ParallelEnv::rank();
  ```
  ✅ **已修复**: 同进程多 block 现在正确走本地内存拷贝路径。

### 2.4 HaloExchange::copy_local(arr, face, block, neighbor)

- **文件**: [include/parallel/halo_exchange.hxx:205-222](include/parallel/halo_exchange.hxx#L205-L222)
- **状态**: ✅ **已实现** — 从邻居 block 的 `target_face` 侧内部区域拷贝到当前 block 的 `face` 侧 ghost
- 同 block 自周期保持原有行为（从自身对面拷贝）

### 2.5 BoundaryConditionApplier::apply_edge_ghost(lb)

- **文件**: [include/bc/bc_applier.hxx:934-1045](include/bc/bc_applier.hxx#L934-L1045)
- **调用时序**: **必须在 HaloExchange 之后**（需要 MPI 传递的 ghost 数据）
- **修饰**: **Cell prim ghost**，12 条边的重叠区域（两个 ghost 面的交集）
- **策略**: 
  - 两个面都是 BC → 取两面值的平均
  - 只有一个面是 BC → 取该面的值
  - 两个面都是 parallel neighbor → 跳过

### 2.6 BoundaryConditionApplier::apply_corner_ghost(lb)

- **文件**: [include/bc/bc_applier.hxx:1051-1126](include/bc/bc_applier.hxx#L1051-L1126)
- **调用时序**: **必须在 HaloExchange 之后**（在 apply_edge_ghost 之后）
- **修饰**: **Cell prim ghost**，8 个角点区域（三个 ghost 面的交集）
- **策略**: 对所有 BC 面的值取算术平均

---

## 三、Cell 中心 Conservative 变量 Ghost（通过转换同步）

### 3.1 Field::prim_to_cons(gamma)

- **文件**: [include/field/field.hxx:144-164](include/field/field.hxx#L144-L164)
- **调用时序**: BC + HaloExchange 完成 prim ghost 填充后，**必须调用**
- **修饰**: **Cell cons ghost**（间接——通过转换 prim ghost 值）
- **作用**: 对**所有 cell**（包括 ghost cell）执行 prim → cons 转换。Ghost cell 的 cons 值由已正确填充的 prim ghost 值转换而来。

### 3.2 Field::cons_to_prim(gamma)

- **文件**: [include/field/field.hxx:166-190](include/field/field.hxx#L166-L190)
- **调用时序**: 每个 RK 阶段的 post-update（在 BC 应用之前）
- **修饰**: **Cell prim**（间接——从更新后的 cons 值转换）
- **作用**: 对**所有 cell**（包括 ghost cell）执行 cons → prim 转换

### 3.3 完整 prim↔cons 同步流程

```
时间推进更新 cons (interior only)
      │
      ▼
cons_to_prim (全部 cell)     ← 同步 prim
      │
      ▼
apply_face_ghost             ← 填充 prim ghost
      │
      ▼
exchange_all_halos           ← MPI 交换 prim ghost
      │
      ▼
apply_edge_ghost             ← 修正 prim edge ghost
      │
      ▼
apply_corner_ghost           ← 修正 prim corner ghost
      │
      ▼
prim_to_cons (全部 cell)     ← ★ 同步 cons ghost
      │
      ▼
WCNS 插值 (使用 cons ghost)
```

---

## 四、Cell 中心梯度 Ghost（粘性项 5b 阶段）

### 4.1 ParallelManager::exchange_gradient_halos(blocks)

- **文件**: [include/parallel/parallel_manager.hxx:207-217](include/parallel/parallel_manager.hxx#L207-L217)
- **调用时机**: main.cpp 中 `ViscidRHS::compute_gradients()` 之后（[main.cpp:411](src/main.cpp#L411)）
- **修饰**: **Cell gradient ghost**，12 个梯度数组
- **交换内容**: `du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dT_dx, dT_dy, dT_dz`
- **实现**: 调用 `HaloExchange::exchange_multi(arrays, block)` 将所有 12 个数组打包到单个 MPI 消息中

---

## 五、Cell 中心粘性通量 Ghost（粘性项 5c 阶段）

### 5.1 ParallelManager::exchange_viscous_flux_halos(blocks)

- **文件**: [include/parallel/parallel_manager.hxx:223-238](include/parallel/parallel_manager.hxx#L223-L238)
- **调用时机**: main.cpp 中 `ViscidRHS::compute_cell_viscous_flux()` 之后（[main.cpp:419](src/main.cpp#L419)）
- **修饰**: **Cell viscous flux ghost**，15 个数组
- **交换内容**: `vis_x.{f1,f2,f3,f4,f5}, vis_y.{f1,f2,f3,f4,f5}, vis_z.{f1,f2,f3,f4,f5}`
- **实现**: 调用 `HaloExchange::exchange_multi(arrays, block)`

---

## 六、面通量 Ghost（Face Flux Exchange）

### 6.1 ParallelManager::exchange_flux_halos(blocks)

- **文件**: [include/parallel/parallel_manager.hxx:189-193](include/parallel/parallel_manager.hxx#L189-L193)
- **调用时机**: main.cpp 中 Riemann 求解器之后、InviscidRHS 之前（[main.cpp:385](src/main.cpp#L385)）
- **修饰**: **无粘面通量 ghost face**，3 方向 × 5 分量 = 15 个 face 数组
- **交换面数**: `ng + 1 = 4` 个面（比 cell ghost 多 1 层，因为 6 点中心差分需要 ±3 面的通量差）

### 6.2 FluxHaloExchange::setup(block, zone_to_block)

- **文件**: [include/parallel/flux_halo_exchange.hxx:24-106](include/parallel/flux_halo_exchange.hxx#L24-L106)
- **调用**: `ParallelManager::initialize()` 中（[parallel_manager.hxx:140](include/parallel/parallel_manager.hxx#L140)）
- **面索引映射**:

| 面类型 | send 范围 | recv 范围 |
|--------|----------|----------|
| MIN (0,2,4) | `[ng, 2*ng]` = 4 faces | `[0, ng]` = 4 faces |
| MAX (1,3,5) | `[ndim-2*ng, ndim-ng]` = 4 faces | `[ndim-ng, ndim]` = 4 faces |

- 缓冲区大小: `5 分量 × (ng+1) × dim1 × dim2`

### 6.3 FluxHaloExchange::pack_flux(fv, dir, info)

- **文件**: [include/parallel/flux_halo_exchange.hxx:112-213](include/parallel/flux_halo_exchange.hxx#L112-L213)
- **修饰**: 将 `ng+1` 个面切片的 5 个通量分量从 FluxVars 打包到连续发送缓冲区

### 6.4 FluxHaloExchange::unpack_flux(fv, dir, info)

- **文件**: [include/parallel/flux_halo_exchange.hxx:219-313](include/parallel/flux_halo_exchange.hxx#L219-L313)
- **修饰**: 从接收缓冲区解包到 ghost face 位置

### 6.5 FluxHaloExchange::copy_local_flux(block, face, dir, info, neighbor)

- **文件**: [include/parallel/flux_halo_exchange.hxx:319-403](include/parallel/flux_halo_exchange.hxx#L319-L403)
- **修饰**: **同进程 direct copy** 面通量 ghost
- MIN 面: 从 neighbor 的 `[nbr_ndim-2*ng-1, nbr_ndim-ng-1]` 拷贝
- MAX 面: 从 neighbor 的 `[ng, 2*ng]` 拷贝

### 6.6 FluxHaloExchange::exchange(block, all_blocks)

- **文件**: [include/parallel/flux_halo_exchange.hxx:409-518](include/parallel/flux_halo_exchange.hxx#L409-L518)
- **3 阶段协议**:
  1. **Phase 1**: 对 remote 面 post `MPI_Irecv`
  2. **Phase 2**: 
     - 同进程邻居: `copy_local_flux`（跳过自周期 `block_id == neighbor->block_id`）
     - Remote: `pack_flux → MPI_Isend`
  3. **Phase 3**: `MPI_Waitall → unpack_flux`
- **自周期面正确处理**: `block_id == neighbor->block_id` 时跳过拷贝（Riemann 求解器已用正确的周期 ghost 数据计算了所有面的通量）

### 6.7 FluxHaloExchange::exchange_face_arrays(face_arrs, dir, block)

- **文件**: [include/parallel/flux_halo_exchange.hxx:524-646](include/parallel/flux_halo_exchange.hxx#L524-L646)
- **修饰**: **通用单分量面数组** MPI 交换（不限于 FluxVars）
- **调用方**: `ViscidRHS::compute_gradients()` (5b) 和 `ViscidRHS::exchange_and_assemble_face_flux()` (5d)
- **交换面数**: `ng + 1` 个面
- **仅处理 MPI remote 邻居**；同进程 copy 由调用方单独处理

### 6.8 ViscidRHS::copy_local_face_array(my_arr, nbr_arr, dir, ...)

- **文件**: [include/scheme/viscid_rhs.hxx:73-133](include/scheme/viscid_rhs.hxx#L73-L133)
- **修饰**: **同进程 face ghost 拷贝**（粘性项面数组）
- **适用范围**: 粘性管道中所有面数组（速度/温度插值面、面乘积、Cartesian 通量面）
- **复制数量**: `ng + 1` 个面切片
- **面索引**:
  - MIN 面接收: ghost face `[0, ng]` ← neighbor face `[nbr_ndim-2*ng, nbr_ndim-ng]`
  - MAX 面接收: ghost face `[ndim-ng, ndim]` ← neighbor face `[ng, 2*ng]`

---

## 七、粘性通量管道的 Ghost 操作序列

粘性项需要在每个 RK 阶段执行以下子步骤，每个步骤都涉及 ghost 操作：

### 7.1 5a: ViscidRHS::interp_to_faces(lb, cfg)

- **文件**: [include/scheme/viscid_rhs.hxx:162-195](include/scheme/viscid_rhs.hxx#L162-L195)
- **修饰**: **Face 值**（u, v, w, T 的 face 插值），依赖 cell cons ghost 有效
- **策略**: 周期面使用需要 ghost 数据的中心模板；非周期面使用 5 阶单侧模板

### 7.2 5b: ViscidRHS::compute_gradients(lb, all_blocks, flux_ex, cfg)

- **文件**: [include/scheme/viscid_rhs.hxx:201-390+]
- **Ghost 操作**:
  1. 构成面乘积 `phi_face × face_metric`
  2. 调用 `flux_ex.exchange_face_arrays()` 交换面乘积 ghost
  3. 调用 `copy_local_face_array()` 进行同进程面乘积 ghost 拷贝
  4. 拷贝后重新构成面乘积
  5. 调用 `InterpDiff::derivative_from_faces()` 求导（使用 ghost face 值）

### 7.3 5c 后: exchange_gradient_halos

- 交换 12 个 cell-center 梯度数组的 ghost cell → 5d 需要

### 7.4 5c 后: exchange_viscous_flux_halos

- 交换 15 个 cell-center Cartesian 粘性通量的 ghost cell → 5d 需要

### 7.5 5d-step1: ViscidRHS::interp_cart_flux_to_faces(lb, dir, cfg)

- **文件**: [include/scheme/viscid_rhs.hxx:830-858]
- **修饰**: **Face 值**（12 个 Cartesian 粘性通量分量插值到面）

### 7.6 5d-step2+3: ViscidRHS::exchange_and_assemble_face_flux(...)

- **文件**: [include/scheme/viscid_rhs.hxx:864-942]
- **Ghost 操作**:
  1. `flux_ex.exchange_face_arrays()` 交换 12 个面插值 Cartesian 通量 ghost
  2. `copy_local_face_array()` × 12 同进程拷贝
  3. 使用面插值值和面度量系数组装 `vis_xi/eta/zeta`

### 7.7 5e: ViscidRHS::compute_rhs(lb)

- 6 阶中心差分粘性面通量 → 累加到 rhs
- 需要 ghost face 通量正确

---

## 八、Ghost 初始化阶段完整时序

### 8.1 Grid 层面（所有进程同时执行）

```
CGNSReader::read_all(zones, ng)
  └─ Grid::allocate(..., ng)          ← 分配 core 尺寸，ng 暂未生效
      │
      ▼
Grid::extend_ghost_layers()           ← 扩展并填充 ghost 节点
  ├─ fill_ghost_face_periodic()       ← 周期 ghost 节点
  └─ fill_ghost_face_extrapolate()    ← 非周期 ghost 节点
      │
      ▼
Grid::compute_cell_centers()          ← 使用 ghost 节点计算 cell 中心
Grid::compute_cell_volumes()          ← 包括 ghost cell 体积
      │
      ▼
Grid::fix_interface_ghost()           ← 修正 inter-zone ghost 节点
  └─ 再次 compute_cell_centers()
  └─ 再次 compute_cell_volumes()
      │
      ▼
Grid::compute_metrics()               ← 使用 ghost 节点计算度量系数
Grid::compute_face_metrics()          ← 使用 cell 度量系数计算面度量
```

### 8.2 LocalBlock 层面（每个进程只构建自己的块）

```
LocalBlock::from_full_zone()          ← 深拷贝整个 zone（含 ghost）
  └─ build_neighbors()
      ├─ 周期面: target_rank=-1, active=true
      └─ BC 面: active=false

LocalBlock::from_sub_zone()           ← 从 full zone 提取子区域
  ├─ 拷贝节点坐标 (含 ghost offset)
  ├─ compute_cell_centers()
  ├─ compute_cell_volumes()
  ├─ 拷贝 full zone 的周期连接信息
  ├─ extract_metrics_from()           ← ★ 从 full zone 提取（含 ghost 度量系数）
  ├─ 拷贝 BC patches
  └─ build_neighbors()
      ├─ 周期面 (同 rank): target_rank=-1
      ├─ 周期面 (跨 rank): target_rank=donor_rank
      ├─ 内部切分面: target_rank=neighbor_rank
      └─ BC 面: active=false
```

### 8.3 Field 初始化（per-block）

```
FlowInitializer::initialize(lb, cfg)
  └─ 只填 interior cells [ng, ng+nci_core-1]
  └─ prim_to_cons()                  ← 转换 interior cons
      │
      ▼
apply_face_ghost(lb, cfg)            ← 填充 BC 面和自周期面的 prim ghost
      │
      ▼
exchange_all_halos(blocks)           ← MPI 交换 prim ghost
      │
      ▼
apply_edge_ghost(lb)                 ← 修正 prim edge ghost
apply_corner_ghost(lb)               ← 修正 prim corner ghost
      │
      ▼
prim_to_cons(gamma)                  ← ★ 同步 cons ghost
```

---

## 九、每个 RK 阶段的 Ghost 处理完整时序

```
┌─ Stage 0 ──────────────────────────────────────────────────────┐
│ Save Q0 (全部 cell, 包括 ghost)                                 │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ 无粘 RHS ─────────────────────────────────────────────────────┐
│ WCNS 插值 (cons cell → face L/R)      [需要 cons ghost 正确]    │
│ Riemann 求解器 (face L/R → face flux)  [所有 face 都计算]       │
│ exchange_flux_halos                    [交换 ghost face 通量]   │
│ InviscidRHS::compute                   [使用 ghost face 通量]   │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ 粘性 RHS (如果启用) ──────────────────────────────────────────┐
│ 5a: interp_to_faces (u,v,w,T)          [需要 cons ghost]       │
│ 5b: compute_gradients                  [交换 face 乘积 ghost]  │
│     └─ exchange_face_arrays + copy_local_face_array            │
│ exchange_gradient_halos                 [交换梯度 ghost cell]   │
│ 5c: compute_cell_viscous_flux           [需要梯度 ghost]        │
│ exchange_viscous_flux_halos             [交换粘性通量 ghost]    │
│ 5d: interp_cart_flux_to_faces          [需要通量 ghost cell]    │
│     exchange_and_assemble_face_flux     [交换 face 通量 ghost]  │
│     └─ exchange_face_arrays + copy_local_face_array            │
│ 5e: compute_rhs                         [使用 ghost face 通量] │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ 时间推进 ──────────────────────────────────────────────────────┐
│ advance_stage / LuSgs::advance          [只更新 interior cells]  │
└────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─ Post-stage Ghost 刷新 ────────────────────────────────────────┐
│ cons_to_prim (全部 cell)                ← 同步 prim             │
│ apply_face_ghost                        ← BC 面 + 自周期 prim   │
│ exchange_all_halos                      ← MPI 交换 prim ghost   │
│ apply_edge_ghost                        ← prim edge 修正        │
│ apply_corner_ghost                      ← prim corner 修正      │
│ prim_to_cons (全部 cell)                ← ★ 同步 cons ghost     │
└────────────────────────────────────────────────────────────────┘
```

---

## 十、所有涉及 Ghost 操作的数据类型汇总

| 数据类型 | 存储位置 | Ghost 维度 | 填充者 | 时机 |
|---------|---------|-----------|--------|------|
| **Node 坐标** (x,y,z) | `Grid::node_*` | `(ni_tot, nj_tot, nk_tot)` = `(ni_core+2*ng, ...)` | `fill_ghost_nodes` | 初始化 |
| **Cell 中心坐标** | `Grid::cell_*` | `(nci, ncj, nck)` | `compute_cell_centers` | 初始化 |
| **Cell 体积** | `Grid::cell_vol` | `(nci, ncj, nck)` | `compute_cell_volumes` | 初始化 |
| **Cell 度量系数** (met_xi/eta/zeta) | `Grid::met_*` | `(nci, ncj, nck)` | `compute_metrics` / `extract_metrics_from` | 初始化 |
| **面度量系数** (face_xi/eta/zeta) | `Grid::face_*` | `(nci+1,ncj,nck)` 等 | `compute_face_metrics` / `extract_metrics_from` | 初始化 |
| **Prim 变量** (rho,u,v,w,p) | `Field::prim.*` | `(nci, ncj, nck)` — 全部 cell | BC Applier + HaloExchange | 初始化 + 每 RK 阶段 |
| **Cons 变量** (rho,rhou,rhov,rhow,rhoE) | `Field::cons.*` | `(nci, ncj, nck)` — 全部 cell | `prim_to_cons` | 每次 prim ghost 更新后 |
| **RHS** | `Field::rhs.*` | `(nci, ncj, nck)` — interior only | `InviscidRHS::compute` + `ViscidRHS::compute_rhs` | 每 RK 阶段 |
| **Q0 快照** | `Field::Q0.*` | `(nci, ncj, nck)` — 全部 cell | main loop (stage 0 begin) | 每时间步 |
| **面插值 L/R 状态** (ql/qr) | `Field::ql_*`, `Field::qr_*` | face 尺寸 | WCNS 插值 | 每 RK 阶段 |
| **无粘面通量** (inv_xi/eta/zeta) | `Field::inv_*` | face 尺寸 | Riemann 求解器 + FluxHaloExchange | 每 RK 阶段 |
| **粘性面通量** (vis_xi/eta/zeta) | `Field::vis_*` | face 尺寸 | ViscidRHS (5d 组装) | 每 RK 阶段 |
| **速度/温度面插值** (u/v/w/T_face) | `Field::*_face_*` | face 尺寸 | ViscidRHS (5a 插值) | 每 RK 阶段 |
| **梯度** (du_dx .. dT_dz, 12 arrays) | `Field::d*` | `(nci, ncj, nck)` — 全部 cell | ViscidRHS (5b) + HaloExchange | 每 RK 阶段 |
| **Cartesian 粘性通量** (vis_x/y/z, 15 arrays) | `Field::vis_*` | `(nci, ncj, nck)` — 全部 cell | ViscidRHS (5c) + HaloExchange | 每 RK 阶段 |
| **Cartesian 通量面插值** (12 arrays) | `Field::vis_*_face` | face 尺寸 | ViscidRHS (5d-step1) + FluxHaloExchange | 每 RK 阶段 |

---

## 十一、已知问题

### ✅ Bug 1: HaloExchange 同进程多块误用 MPI 自发送 **[已修复]**

- **文件**: [halo_exchange.hxx:29-31](include/parallel/halo_exchange.hxx#L29-L31)
- **修复**: 添加 `ni.target_rank != ParallelEnv::rank()` 条件排除同进程

### ✅ Bug 2: HaloExchange 同进程跨块周期拷贝来源错误 **[已修复]**

- **文件**: [halo_exchange.hxx:277-312](include/parallel/halo_exchange.hxx#L277-L312)
- **修复**: 区分同 block 自周期 vs 同进程不同 block，后者通过 `copy_local` 从邻居 block 拷贝

### ✅ Bug 3: Non-blocking 路径不处理同进程通信 **[部分修复]**

- **文件**: [halo_exchange.hxx:380-398](include/parallel/halo_exchange.hxx#L380-L398)
- **状态**: `start_exchange` 现在处理同 block 自周期；同进程多 block 仍只在 `exchange_multi` 中支持
- **影响**: 当前未使用非阻塞模式

### InterpDiff 模块的 ghost 处理说明

- **文件**: [interp_diff.hxx:107-108](include/scheme/interp_diff.hxx#L107-L108)
- `interp_line` 和 `diff_line` 的 `left_periodic/right_periodic` 参数**当前未使用**（标记为 `(void)`）
- 一阶单侧格式在绝对数组端点始终应用，不论是否周期
- 周期 ghost 数据的正确性依仗于外部填充（BC Applier + HaloExchange）

---

## 十二、关键设计原则

1. **初始值只填 interior**: `FlowInitializer` 只写 `[ng, ng+nci_core-1]`，ghost 留待 BC/Halo 填充
2. **时间推进只改 interior**: `TimeIntegrator::advance_stage` 和 `LuSgs::advance` 只更新 interior cell
3. **BC → Halo → Edge → Corner → prim_to_cons**: 严格的执行顺序保证 ghost cons 正确
4. **Face flux 自周期无操作**: FluxHaloExchange 跳过自 block 拷贝，依赖 BC 周期 ghost 使 Riemann 求解器面通量正确
5. **子块度量系数从 full zone 提取**: 避免内部切分面上的外推 ghost 导致度量系数错误
6. **ng+1 面交换**: 面通量交换比 cell ghost 多 1 层，因为 6 点差分需要 ±3 面的通量
