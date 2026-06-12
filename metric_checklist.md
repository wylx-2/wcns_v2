# 度量系数计算 — 完整流程检查清单

> 对应规范文档: [metric.md](metric.md) | 生成日期: 2026-06-11

本文档以**每一步变量初始化、更新、计算**的粒度，追踪从 CGNS 网格文件读取到最终面度量系数的完整流程。每步标注精确的代码位置和变量状态。

---

## 目录

1. [阶段 0: 配置读取](#阶段-0-配置读取)
2. [阶段 1: CGNS 网格读取](#阶段-1-cgns-网格读取)
3. [阶段 2: Ghost 层扩展](#阶段-2-ghost-层扩展)
4. [阶段 3: Ghost 节点填充](#阶段-3-ghost-节点填充)
5. [阶段 4: 跨 Zone 界面 Ghost 修正](#阶段-4-跨-zone-界面-ghost-修正)
6. [阶段 5: Cell Center 重算](#阶段-5-cell-center-重算)
7. [阶段 6: 域分解与子块提取](#阶段-6-域分解与子块提取)
8. [阶段 7: SCMM 度量系数计算](#阶段-7-scmm-度量系数计算)
9. [阶段 8: 面度量插值](#阶段-8-面度量插值)
10. [阶段 9: Uniform 网格路径](#阶段-9-uniform-网格路径)
11. [完整流水线时序图](#完整流水线时序图)
12. [变量生命周期汇总表](#变量生命周期汇总表)
13. [人工检查清单](#人工检查清单)

---

## 阶段 0: 配置读取

**文件**: [src/io/config_reader.cpp](src/io/config_reader.cpp)

**关键变量**:
```cpp
Config cfg;                          // 默认构造
cfg.grid_metrics = "auto";           // 第 168 行 config.h — 默认 SCMM
```

- `"auto"` → SCMM 从网格几何计算
- `"uniform"` → 解析均匀 Cartesian 度量

**检查项**:
- [ ] `grid_metrics` 是否按预期设置（`auto` 或 `uniform`）

---

## 阶段 1: CGNS 网格读取

**文件**: [src/io/cgns_reader.cpp](src/io/cgns_reader.cpp)

### 1.1 文件打开

```cpp
// cgns_reader.cpp:14-31
CGNSReader reader;
reader.open(grid_file);
```

### 1.2 逐 Zone 读取 → Grid 对象

**入口**: `CGNSReader::read_all()` → `read_zone()` (行 296-303)

```cpp
for (Int z = 1; z <= num_zones(1); ++z) {
    Grid g;                          // 默认构造 — 所有指针为空, ni=nj=nk=nci=ncj=nck=0
    read_zone(1, z, g, ng);          // 传入 ng (ghost 层数)
    grids.push_back(std::move(g));
}
```

### 1.3 read_zone — 逐步变量初始化

**文件**: [src/io/cgns_reader.cpp:52-90](src/io/cgns_reader.cpp#L52-L90)

#### Step 1.3a: 读取 zone 尺寸

```cpp
// 行 57-65
cgsize_t sizes[9];
cg_zone_read(file_id_, B, Z, zone_name, sizes);
Int ni = sizes[0];    // 顶点数 (i 方向)
Int nj = sizes[1];    // 顶点数 (j 方向)
Int nk = sizes[2];    // 顶点数 (k 方向)
```

**此时 ni/nj/nk 含义**: **物理域顶点数（不含 ghost）**

#### Step 1.3b: 分配 Grid 数组

```cpp
// 行 73 — 调用 Grid::allocate(ni, nj, nk, ng)
grid.allocate(ni, nj, nk, ng);
```

**Grid::allocate 内部**: [include/grid/grid.hxx:13-27](include/grid/grid.hxx#L13-L27)

```cpp
void Grid::allocate(Int ni_, Int nj_, Int nk_, Int ng_) {
    ni = ni_; nj = nj_; nk = nk_; ng = ng_;          // 顶点维度 = 物理域核心
    nci = ni - 1; ncj = nj - 1; nck = nk - 1;       // cell 维度 = 顶点-1
    ni_core = ni; nj_core = nj; nk_core = nk;        // 保存核心尺寸

    node_x.allocate(ni, nj, nk);   // ← 仅为核心域分配！此时无 ghost
    node_y.allocate(ni, nj, nk);
    node_z.allocate(ni, nj, nk);

    cell_x.allocate(nci, ncj, nck);  // 核心域 cell 坐标
    cell_y.allocate(nci, ncj, nck);
    cell_z.allocate(nci, ncj, nck);

    cell_vol.allocate(nci, ncj, nck);  // 核心域 cell 体积
}
```

**变量状态检查**:
| 变量 | 值 | 含义 |
|------|-----|------|
| `ni` | CGNS 顶点数 | 核心域顶点 (无 ghost) |
| `nj, nk` | CGNS 顶点数 | 核心域顶点 |
| `nci` | `ni-1` | 核心域 cell 数 |
| `ncj, nck` | `nj-1, nk-1` | 核心域 cell 数 |
| `ng` | 配置值 (3) | ghost 层数 (尚未扩展) |
| `ni_core` | `ni` | 记录核心尺寸 |
| `node_x/y/z` | `(ni,nj,nk)` | 仅核心域，无 ghost |
| `cell_x/y/z` | `(nci,ncj,nck)` | 仅核心域 |
| `cell_vol` | `(nci,ncj,nck)` | 仅核心域 |

#### Step 1.3c: 读取节点坐标

```cpp
// 行 79 — read_coords(B, Z, grid)
// 行 92-119
// 对 d=0,1,2 分别读取 CoordinateX/Y/Z
cg_coord_read(file_id_, B, Z, "CoordinateX", RealDouble, rmin, rmax, grid.node_x.data());
cg_coord_read(file_id_, B, Z, "CoordinateY", RealDouble, rmin, rmax, grid.node_y.data());
cg_coord_read(file_id_, B, Z, "CoordinateZ", RealDouble, rmin, rmax, grid.node_z.data());
```

**CGNS 读取范围**: `rmin=[1,1,1], rmax=[ni,nj,nk]` — 1-based

**存储**: `node_x(0..ni-1, 0..nj-1, 0..nk-1)` — 0-based

#### Step 1.3d: 计算 cell 中心坐标

```cpp
// 行 82 — grid.compute_cell_centers()
```

**文件**: [include/grid/grid.hxx:29-47](include/grid/grid.hxx#L29-L47)

```cpp
for (Int k = 0; k < nck; ++k)     // nck = 核心 cell 数
for (Int j = 0; j < ncj; ++j)
for (Int i = 0; i < nci; ++i) {
    Real cx=0, cy=0, cz=0;
    for (Int dk=0; dk<=1; ++dk)
    for (Int dj=0; dj<=1; ++dj)
    for (Int di=0; di<=1; ++di) {
        cx += node_x(i+di, j+dj, k+dk);   // 8 个角点求和
        cy += node_y(i+di, j+dj, k+dk);
        cz += node_z(i+di, j+dj, k+dk);
    }
    cell_x(i,j,k) = cx * 0.125;           // ÷8 得中心
    cell_y(i,j,k) = cy * 0.125;
    cell_z(i,j,k) = cz * 0.125;
}
```

**约定**: cell (i,j,k) 被顶点 (i,j,k) 到 (i+1,j+1,k+1) 包围

#### Step 1.3e: 计算 cell 体积 (初始近似)

```cpp
// 行 83 — grid.compute_cell_volumes()
```

**文件**: [include/grid/grid.hxx:49-61](include/grid/grid.hxx#L49-L61)

```cpp
for (Int k = 0; k < nck; ++k)
for (Int j = 0; j < ncj; ++j)
for (Int i = 0; i < nci; ++i) {
    Real dx = node_x(i+1,j,k) - node_x(i,j,k);   // 粗略 Δx
    Real dy = node_y(i,j+1,k) - node_y(i,j,k);   // 粗略 Δy
    Real dz = node_z(i,j,k+1) - node_z(i,j,k);   // 粗略 Δz
    cell_vol(i,j,k) = std::abs(dx * dy * dz);     // 均匀 Cartesian 近似
}
```

**注意**: 这是初始近似值。SCMM 计算后会通过 `jacobian` 给出精确体积。后续 `extract_metrics_from` 也会用 donor 的 `cell_vol` 覆盖。

#### Step 1.3f: 读取边界条件

```cpp
// 行 86 — read_bc(B, Z, grid)
// 每个 BC patch 记录: name, type, face, imin..kmax (1-based, 节点索引)
```

#### Step 1.3g: 读取 1-to-1 连通性

```cpp
// 行 89 — read_1to1(B, Z, grid)
// 每个 connection 记录:
//   name, donor_name
//   imin..kmax (当前 zone, 1-based 节点索引)
//   donor_imin..donor_kmax (donor zone, 1-based 节点索引)
//   transform[3]
//   is_periodic, translation[3], rotation_center[3], rotation_angle[3]
```

---

## 阶段 2: Ghost 层扩展

**文件**: [include/grid/grid.hxx:80-129](include/grid/grid.hxx#L80-L129)

**触发**: `Grid::extend_ghost_layers()` — 在 `ParallelManager::initialize` 中对所有 zone 调用

### Step 2.1: 保存旧节点数据

```cpp
Int ni_src = ni_core;  // 保存核心尺寸 (与当前 ni 相同)
Int nj_src = nj_core;
Int nk_src = nk_core;

MultiArray3D<Real> old_nx = std::move(node_x);  // 窃取旧数据
MultiArray3D<Real> old_ny = std::move(node_y);
MultiArray3D<Real> old_nz = std::move(node_z);
```

**此时**: `old_nx` 的维度为 `(ni_src, nj_src, nk_src)` — 即核心域，无 ghost

### Step 2.2: 分配扩展节点数组

```cpp
Int nie = ni_src + 2 * ng;   // 核心 + 2×ng ghost (每侧 ng 层)
Int nje = nj_src + 2 * ng;
Int nke = nk_src + 2 * ng;

node_x.allocate(nie, nje, nke);   // 现在 node_x 有 ghost 空间
node_y.allocate(nie, nje, nke);
node_z.allocate(nie, nje, nke);
```

### Step 2.3: 复制核心节点到 interior

```cpp
// 带 ng 偏移
for (Int k = 0; k < nk_src; ++k)
for (Int j = 0; j < nj_src; ++j)
for (Int i = 0; i < ni_src; ++i) {
    node_x(i+ng, j+ng, k+ng) = old_nx(i,j,k);
    node_y(i+ng, j+ng, k+ng) = old_ny(i,j,k);
    node_z(i+ng, j+ng, k+ng) = old_nz(i,j,k);
}
```

**节点索引约定** (扩展后):

```
  i 范围:    [0 .. ng-1]           ghost (IMIN 侧)
             [ng .. ng+ni_src-1]   核心节点
             [ng+ni_src .. nie-1]  ghost (IMAX 侧)
```

### Step 2.4: 更新维度

```cpp
ni = nie;  nj = nje;  nk = nke;           // 顶点总数 (含 ghost)
nci = ni - 1;  ncj = nj - 1;  nck = nk - 1;  // cell 总数 (含 ghost)
```

**此时 `ni_core, nj_core, nk_core` 保持不变** — 仍为原始核心顶点数。

### Step 2.5: 填充 ghost 节点

```cpp
fill_ghost_nodes();  // 见阶段 3
```

### Step 2.6: 分配扩展 cell 数组并重算

```cpp
cell_x.allocate(nci, ncj, nck);   // 重新分配 — 现在包含 ghost cell
cell_y.allocate(nci, ncj, nck);
cell_z.allocate(nci, ncj, nck);
cell_vol.allocate(nci, ncj, nck);

compute_cell_centers();   // 用扩展后的节点坐标 (含 ghost 节点) 重算
compute_cell_volumes();   // 用扩展后的节点坐标重算
```

**变量状态检查** (扩展后):

| 变量 | 值 | 说明 |
|------|-----|------|
| `ni, nj, nk` | `nie, nje, nke` | 顶点总数 = 核心 + 2×ng |
| `nci, ncj, nck` | `ni-1, nj-1, nk-1` | cell 总数 = 核心 + 2×ng |
| `ng` | 3 | ghost 层数 (不变) |
| `ni_core, nj_core, nk_core` | 原始核心顶点数 | 不变 |
| `node_x/y/z` | `(ni,nj,nk)` | 含 ghost 节点 (已填充) |
| `cell_x/y/z` | `(nci,ncj,nck)` | 含 ghost cell (从扩展后节点计算) |
| `cell_vol` | `(nci,ncj,nck)` | 含 ghost cell 的粗略体积 |

---

## 阶段 3: Ghost 节点填充

**文件**: [include/grid/grid.hxx:131-140](include/grid/grid.hxx#L131-L140), [341-563](include/grid/grid.hxx#L341-L563)

### Step 3.1: 分面调度

```cpp
void Grid::fill_ghost_nodes() {
    for (int face = 0; face < 6; ++face) {
        const Connectivity* conn = find_periodic_connection(face);
        if (conn)
            fill_ghost_face_periodic(face, *conn);   // 周期面 → 拷贝对面 interior
        else
            fill_ghost_face_extrapolate(face);        // 其他面 → 线性外推
    }
}
```

### Step 3.2: 周期面 Ghost 填充 (fill_ghost_face_periodic)

**文件**: [include/grid/grid.hxx:341-496](include/grid/grid.hxx#L341-L496)

**逻辑**: 将对面 interior 的节点坐标拷贝到当前 ghost，减去平移向量。

以 IMIN (face=0) 为例:

```cpp
// 判断 donor 面是高侧还是低侧
bool donor_is_high = (conn.donor_imin == ni_core);  // 1-based

// Ghost 循环: i = 0 .. ng-1
for (Int i = 0; i < ng; ++i) {
    Int d = ng - i;    // ghost 层深度: 1(最近) .. ng(最远)

    Int donor_index;
    if (donor_is_high)
        donor_index = ng + ni_core - 1 - d;   // 从 MAX 面向内 d 层
    else
        donor_index = ng + d;                  // 从 MIN 面向内 d 层

    // 拷贝并减平移
    node_x(i, j, k) = node_x(donor_index, j, k) - tx;
    node_y(i, j, k) = node_y(donor_index, j, k) - ty;
    node_z(i, j, k) = node_z(donor_index, j, k) - tz;
}
```

**关键**: 周期平移为 `ghost = donor - translation`。如果平移向量正确（由 CGNS 提供），ghost 节点坐标应精确对应物理空间中的正确位置。

### Step 3.3: 非周期面 Ghost 填充 (fill_ghost_face_extrapolate)

**文件**: [include/grid/grid.hxx:498-563](include/grid/grid.hxx#L498-L563)

**逻辑**: 逐层线性外推: `x(-d) = 2*x(-(d-1)) - x(-(d-2))`。

以 IMIN (face=0) 为例:

```cpp
for (Int layer = 1; layer <= ng; ++layer) {
    Int ghost_i = ng - layer;        // 当前 ghost 层
    Int base_i  = ng - layer + 1;    // 内一层 (更靠近面)
    Int prev_i  = ng - layer + 2;    // 内两层

    node_x(ghost_i, j, k) = 2.0 * node_x(base_i, j, k) - node_x(prev_i, j, k);
    node_y(ghost_i, j, k) = 2.0 * node_y(base_i, j, k) - node_y(prev_i, j, k);
    node_z(ghost_i, j, k) = 2.0 * node_z(base_i, j, k) - node_z(prev_i, j, k);
}
```

**注意**: 对于非周期 BC 面（固壁、远场等），外推的 ghost 节点仅用于支撑内部 stencil，不反映真实几何。这对度量系数在边界附近的精度有影响。BC 面的物理处理由 BC Applier 负责。

**检查项**:
- [ ] 周期面的 `translation` 向量是否从 CGNS 正确读取
- [ ] 非周期面的外推 ghost 是否在后续 `fix_interface_ghost` 中被覆盖

---

## 阶段 4: 跨 Zone 界面 Ghost 修正

**文件**: [include/grid/grid.hxx:207-339](include/grid/grid.hxx#L207-L339), [include/parallel/parallel_manager.hxx:37-54](include/parallel/parallel_manager.hxx#L37-L54)

**触发**: `ParallelManager::initialize` 中，对所有 zone 检查 1-to-1 连接，若 donor 为不同 zone，则调用 `fix_interface_ghost`。

```cpp
// parallel_manager.hxx:37-54
for (auto& z : zones) {
    for (int face = 0; face < 6; ++face) {
        const Connectivity* conn = z.find_face_connection(face);
        if (!conn) continue;
        if (conn->donor_name == z.name) continue;   // 同 zone (周期或自连接)

        Grid* donor = nullptr;
        for (auto& d : zones) {
            if (d.name == conn->donor_name) { donor = &d; break; }
        }
        if (donor) {
            z.fix_interface_ghost(face, *donor, *conn);
        }
    }
}
```

### Step 4.1: 确定 donor face

```cpp
// grid.hxx:212-219
int donor_face = -1;
if (conn.donor_imin == conn.donor_imax)
    donor_face = (conn.donor_imin == 1) ? 0 : 1;
else if (conn.donor_jmin == conn.donor_jmax)
    donor_face = (conn.donor_jmin == 1) ? 2 : 3;
else if (conn.donor_kmin == conn.donor_kmax)
    donor_face = (conn.donor_kmin == 1) ? 4 : 5;
```

### Step 4.2: 方向匹配检查

```cpp
if (cur_dim != donor_dim) return;  // 方向不匹配 → 保持外推值 (TODO)
```

### Step 4.3: 拷贝 donor interior 节点到 ghost

以 IMIN ghost (`cur_dim==0, cur_is_high==false`) 且 `donor_is_high==true` 为例:

```cpp
for (Int i = 0; i < ng; ++i) {
    Int d = ng - i;                          // ghost 层深度
    Int donor_i = donor.ng + donor.ni_core - 1 - d;  // donor MAX 面向内 d 层

    node_x(i, j, k) = donor.node_x(donor_i, j, k) - tx;
    node_y(i, j, k) = donor.node_y(donor_i, j, k) - ty;
    node_z(i, j, k) = donor.node_z(donor_i, j, k) - tz;
}
```

**关键**: 这会覆盖阶段 3 中外推的 ghost 节点坐标。对于跨 zone 界面，ghost 节点现在精确对应相邻 zone 的 interior 节点，保证度量系数在界面处的连续性。

**检查项**:
- [ ] `cur_dim == donor_dim` 是否对所有连接成立？（不成立时被跳过）
- [ ] 平移向量 `translation` 是否正确？

---

## 阶段 5: Cell Center 重算

**文件**: [include/parallel/parallel_manager.hxx:61-64](include/parallel/parallel_manager.hxx#L61-L64)

在 `fix_interface_ghost` 修改了 ghost 节点坐标后，必须重算 cell center:

```cpp
// parallel_manager.hxx:61-64
for (auto& z : zones) {
    z.compute_cell_centers();   // 用修正后的 ghost 节点重算
    z.compute_cell_volumes();   // 用修正后的 ghost 节点重算
}
```

**这很重要**: 阶段 2 中 `extend_ghost_layers` 末尾已调用过一次 `compute_cell_centers` 和 `compute_cell_volumes`，但那时 ghost 节点可能还是外推值。现在 ghost 节点已被 interface correction 修正，cell center 必须重算以保证 ghost cell 的坐标正确。

**检查项**:
- [ ] 确认 `compute_cell_centers` 在 `fix_interface_ghost` 之后被调用
- [ ] 确认 `compute_cell_volumes` 在 `fix_interface_ghost` 之后被调用

---

## 阶段 6: 域分解与子块提取

**文件**: [include/parallel/parallel_manager.hxx:74-129](include/parallel/parallel_manager.hxx#L74-L129)

### 完整 Zone (未剖分)

```cpp
// local_block.hxx:39-66
LocalBlock lb = LocalBlock::from_full_zone(zones[zone_id], block_id, ...);
// Grid 直接 deep-copy: lb.grid = full_zone;
// 所有度量数组通过 Grid 拷贝构造函数复制
```

### 子 Zone (Cartesian 剖分)

```cpp
// local_block.hxx:72-251
LocalBlock lb = LocalBlock::from_sub_zone(full_zone, sub, block_id, ...);
```

#### 节点坐标提取:

```cpp
// local_block.hxx:128-134
for (Int k = 0; k < nk_tot; ++k)
for (Int j = 0; j < nj_tot; ++j)
for (Int i = 0; i < ni_tot; ++i) {
    lb.grid.node_x(i,j,k) = full_zone.node_x(i + offset_i, j + offset_j, k + offset_k);
    // offset_i = sub.ci_min = 子块在完整 zone 中的起始 cell 索引
    // 这包括了 ghost 节点 — 因为 full_zone 已有 ghost 层
}
```

#### 度量系数提取 (避免子块重算):

```cpp
// local_block.hxx:185-186
lb.grid.metrics_type = full_zone.metrics_type;
lb.grid.extract_metrics_from(full_zone, offset_i, offset_j, offset_k);
```

**为什么提取而不重算**: 子块的内部剖分边界 ghost 节点来自相邻子块的 interior 区域，但在子块构建时这些 ghost 节点并未被正确填充（没有做子块间的 ghost exchange）。如果在子块上直接调用 `compute_metrics`，内部剖分边界处的 stencil 会读到错误的 ghost cell 坐标，导致边界附近的度量系数错误。

**解决方法**: 从完整 zone 的已有度量数组中按偏移提取。完整 zone 的度量是在所有 ghost 层正确填充后计算的，内部剖分边界处的度量已使用正确的邻居数据。

**检查项**:
- [ ] 子块提取的 `offset_i/j/k` 是否与完整 zone 中的位置对应
- [ ] `extract_metrics_from` 是否正确复制了所有 9 个中心度量 + Jacobian + 3×3 面度量
- [ ] `cell_vol` 是否也从 donor 萃取

---

## 阶段 7: SCMM 度量系数计算

**文件**: [include/grid/grid.hxx:584-774](include/grid/grid.hxx#L584-L774)

**入口**: `Grid::compute_metrics()` → 非 uniform 时走 SCMM 路径

### 变量初始化

```cpp
// 行 591-592
bool fp[6];
build_face_periodic(fp);     // 标记 6 个面是否为周期面
const Real dh = 1.0;         // 计算空间网格间距
```

`fp[i]` 决定 InterpDiff 在对应面上使用何种边界处理 — 但 `interp_line`/`diff_line` 的一侧闭合始终从绝对数组端应用，`fp` 仅作为 API 参数传递。

### Step 7.1: 分配度量数组

```cpp
// 行 597-606
met_xi_x.allocate(nci, ncj, nck);
met_xi_y.allocate(nci, ncj, nck);
met_xi_z.allocate(nci, ncj, nck);
// ... met_eta_* , met_zeta_* (共 9 个)
jacobian.allocate(nci, ncj, nck);
```

**数组维度**: `(nci, ncj, nck)` — 全部 cell (含 ghost)。每个数组在分配后未显式初始化。

### Step 7.2: 计算坐标导数

```cpp
// 行 610-616
MultiArray3D<Real> x_xi(nci,ncj,nck), x_eta(...), x_zeta(...);
MultiArray3D<Real> y_xi(nci,ncj,nck), y_eta(...), y_zeta(...);
MultiArray3D<Real> z_xi(nci,ncj,nck), z_eta(...), z_zeta(...);

compute_coord_deriv(0, dh, fp, x_xi,  y_xi,  z_xi);   // ∂/∂ξ
compute_coord_deriv(1, dh, fp, x_eta, y_eta, z_eta);  // ∂/∂η
compute_coord_deriv(2, dh, fp, x_zeta,y_zeta,z_zeta); // ∂/∂ζ
```

**`compute_coord_deriv` 内部**: [grid.hxx:575-582](include/grid/grid.hxx#L575-L582)

```cpp
void compute_coord_deriv(int dir, Real dh, const bool fp[6],
                         MultiArray3D<Real>& dx_d, ...) {
    InterpDiff::derivative(cell_x, dx_d, dir, dh, ng, fp);
    InterpDiff::derivative(cell_y, dy_d, dir, dh, ng, fp);
    InterpDiff::derivative(cell_z, dz_d, dir, dh, ng, fp);
}
```

### Step 7.3: InterpDiff::derivative — 插值→差分流水线

**文件**: [include/scheme/interp_diff.hxx:195-271](include/scheme/interp_diff.hxx#L195-L271)

以 ξ 方向 (`dir=0`) 为例:

```cpp
// 1. 分配半节点数组: (ni+1, nj, nk) — 每 cell 线多 1 个 face
MultiArray3D<Real> ah;  ah.allocate(ni + 1, nj, nk);

// 2. 逐线处理
for (Int k = 0; k < nk; ++k)
for (Int j = 0; j < nj; ++j) {
    const Real* a_line = &a(0, j, k);   // cell_x 的 i-line (含 ghost)
    Real*      ah_line = &ah(0, j, k);  // 半节点数组
    Real*      da_line = &da(0, j, k);  // 导数输出

    interp_line(a_line, ah_line, ni, lp, rp);  // cell → half-node
    diff_line(ah_line, da_line, ni, dh, lp, rp); // half-node → derivative
}
```

**`interp_line` 逐半节点处理**: [interp_diff.hxx:105-148](include/scheme/interp_diff.hxx#L105-L148)

```
  半节点 h:  0        1        2      3  ... n-3    n-2     n-1     n
            a_{-1/2} a_{1/2} a_{3/2}              a_{n-5/2} a_{n-3/2} a_{n-1/2}

  h=0:     interp_left_1st(a[0..4])       = (315a0 - 420a1 + 378a2 - 180a3 + 35a4)/128
  h=1:     interp_left_2nd(a[0..4])        = (35a0 + 140a1 - 70a2 + 28a3 - 5a4)/128
  h=2:     interp_left_3rd(a[0..4])        = (-5a0 + 60a1 + 90a2 - 20a3 + 3a4)/128
  h=3..n-3: interp_center_6pt(a, h-1)     = 75/128*(a[h-1]+a[h]) - 25/256*(a[h-2]+a[h+1]) + 3/256*(a[h-3]+a[h+2])
  h=n-2:   interp_right_3rd(a[n-5..n-1])  (对称)
  h=n-1:   interp_right_2nd(a[n-5..n-1])  (对称)
  h=n:     interp_right_1st(a[n-5..n-1])  (对称)
```

**`diff_line` 逐 cell 处理**: [interp_diff.hxx:150-189](include/scheme/interp_diff.hxx#L150-L189)

```
  cell i:  0        1        2  ... n-3    n-2     n-1

  i=0:     diff_left_1st(ah, 0, dh)   = (-22ah[0] + 17ah[1] + 9ah[2] - 5ah[3] + ah[4])/(24dh)
  i=1:     diff_left_2nd(ah, 0, dh)   = (ah[0] - 27ah[1] + 27ah[2] - ah[3])/(24dh)
  i=2..n-3: diff_center_6pt(ah, i, dh) = 75/(64dh)*(ah[i+1]-ah[i]) - 25/(384dh)*(ah[i+2]-ah[i-1]) + 3/(640dh)*(ah[i+3]-ah[i-2])
  i=n-2:   diff_right_2nd(ah, n, dh)   (对称)
  i=n-1:   diff_right_1st(ah, n, dh)   = (22ah[n] - 17ah[n-1] - 9ah[n-2] + 5ah[n-3] - ah[n-4])/(24dh)
```

**关键依赖**: 中心插值 `interp_center_6pt(a, h-1)` 需要 `a[h-3..h+2]`，即 6 个 cell 值。中心差分 `diff_center_6pt(ah, i, dh)` 需要 `ah[i-2..i+3]`，即 6 个半节点值。这些 stencil 要求 ghost cell 的坐标值正确——由阶段 2-5 保证。

**检查项**:
- [ ] dh = 1.0（计算空间步长）
- [ ] 确认 `interp_left_1st` 使用的 a[0..4] 在 ghost 区域内正确
- [ ] 确认 `diff_left_1st` 使用的 ah[0..4] 正确
- [ ] 对周期面，ghost cell 坐标必须来自对面 interior（阶段 3 period copy）

### Step 7.4: SCMM 度量组装

**文件**: [grid.hxx:618-724](include/grid/grid.hxx#L618-L724)

**临时数组**: `t_prod(nci,ncj,nck)` 和 `t_deriv(nci,ncj,nck)`

**Lambda**: `mul(a, b)` — 逐点乘积 → `t_prod`
**Lambda**: `accum_deriv(met, coef, dir)` — `derivative(t_prod, t_deriv, dir)` 然后 `met += coef * t_deriv`

以 **ξ̂_x** 为例 (4 步):

```cpp
// ξ̂_x = 0.5*[(z·y_η)_ζ + (y·z_ζ)_η - (z·y_ζ)_η - (y·z_η)_ζ]

met_xi_x.fill(0.0);

// 第1项: +(z·y_η)_ζ * 0.5
mul(cell_z, y_eta);          // t_prod(i,j,k) = cell_z(i,j,k) * y_eta(i,j,k)
accum_deriv(met_xi_x, 0.5, 2);  // ∂(t_prod)/∂ζ → met_xi_x += 0.5 * ∂_ζ(z·y_η)

// 第2项: +(y·z_ζ)_η * 0.5
mul(cell_y, z_zeta);         // t_prod = cell_y * z_zeta
accum_deriv(met_xi_x, 0.5, 1);  // met_xi_x += 0.5 * ∂_η(y·z_ζ)

// 第3项: -(z·y_ζ)_η * 0.5
mul(cell_z, y_zeta);         // t_prod = cell_z * y_zeta
accum_deriv(met_xi_x, -0.5, 1); // met_xi_x -= 0.5 * ∂_η(z·y_ζ)

// 第4项: -(y·z_η)_ζ * 0.5
mul(cell_y, z_eta);          // t_prod = cell_y * z_eta
accum_deriv(met_xi_x, -0.5, 2); // met_xi_x -= 0.5 * ∂_ζ(y·z_η)
```

**每项 `accum_deriv` 内部**:
1. `InterpDiff::derivative(t_prod, t_deriv, dir, dh, ng, fp)` — t_prod 的导数 (又走一遍 interp→diff)
2. `met(i,j,k) += coef * t_deriv(i,j,k)` — 累加到输出

**其他 8 个度量分量**按相同模式计算，共 36 次 `accum_deriv` 调用。每次 `accum_deriv` 都执行一次完整的 `InterpDiff::derivative`。

### Step 7.5: Jacobian 计算

**文件**: [grid.hxx:726-771](include/grid/grid.hxx#L726-L771)

```cpp
// J = (1/3) * [ ∂_ξ(S_ξ) + ∂_η(S_η) + ∂_ζ(S_ζ) ]
// 其中 S_dir = x·met_dir_x + y·met_dir_y + z·met_dir_z

// S_ξ = x*ξ̂_x + y*ξ̂_y + z*ξ̂_z
t_prod(i,j,k) = cell_x*met_xi_x + cell_y*met_xi_y + cell_z*met_xi_z;
InterpDiff::derivative(t_prod, jacobian, 0, dh, ng, fp);  // jacobian ← ∂_ξ(S_ξ)

// S_η → + ∂_η(S_η)
t_prod(i,j,k) = cell_x*met_eta_x + cell_y*met_eta_y + cell_z*met_eta_z;
InterpDiff::derivative(t_prod, t_deriv, 1, dh, ng, fp);
jacobian += t_deriv;

// S_ζ → + ∂_ζ(S_ζ)
t_prod(i,j,k) = cell_x*met_zeta_x + cell_y*met_zeta_y + cell_z*met_zeta_z;
InterpDiff::derivative(t_prod, t_deriv, 2, dh, ng, fp);
jacobian += t_deriv;

// J = sum / 3
jacobian(i,j,k) /= 3.0;
```

**检查项**:
- [ ] `jacobian` 的正负号（右手系网格应为正）
- [ ] Jacobian 初值是否被 `fill(0)` 后正确累加
- [ ] `S_dir` 累加时使用了已计算的 `met_*` 数组
- [ ] `cell_vol` 是否在 SCMM 后需要更新为 `abs(jacobian)`（当前未更新——`cell_vol` 保持阶段 5 的 Cartesian 近似，但历史监控使用时取 `cell_vol`）

---

## 阶段 8: 面度量插值

**文件**: [include/grid/grid.hxx:848-935](include/grid/grid.hxx#L848-L935)

**入口**: `Grid::compute_face_metrics()`

**目的**: 将 cell 中心的度量系数插值到半节点（面），因为 Riemann 求解器在面上需要面法向量。

### Step 8.1: 分配面度量数组

```cpp
// X 面 (i+1/2): (nci+1) × ncj × nck
face_xi_x.allocate(nci + 1, ncj, nck);
face_xi_y.allocate(nci + 1, ncj, nck);
face_xi_z.allocate(nci + 1, ncj, nck);

// Y 面 (j+1/2): nci × (ncj+1) × nck
face_eta_x.allocate(nci, ncj + 1, nck);
face_eta_y.allocate(nci, ncj + 1, nck);
face_eta_z.allocate(nci, ncj + 1, nck);

// Z 面 (k+1/2): nci × ncj × (nck+1)
face_zeta_x.allocate(nci, ncj, nck + 1);
face_zeta_y.allocate(nci, ncj, nck + 1);
face_zeta_z.allocate(nci, ncj, nck + 1);
```

**维度**: 面数组在插值方向多 1 个元素（nci+1 个面对应 nci 个 cell）

### Step 8.2: 调用 interp_to_faces

```cpp
InterpDiff::interp_to_faces(met_xi_x, face_xi_x, 0, ng, fp);
InterpDiff::interp_to_faces(met_xi_y, face_xi_y, 0, ng, fp);
InterpDiff::interp_to_faces(met_xi_z, face_xi_z, 0, ng, fp);
// 同理 η, ζ 方向
```

**`interp_to_faces` 内部**: [interp_diff.hxx:273-321](include/scheme/interp_diff.hxx#L273-L321)

以 ξ 方向为例:

```cpp
for (Int k = 0; k < nk; ++k)
for (Int j = 0; j < nj; ++j) {
    const Real* a_line = &a(0, j, k);   // met_xi_x 的 i-line
    Real*      af_line = &af(0, j, k);  // face_xi_x 的 i-line

    interp_line(a_line, af_line, ni, lp, rp);  // 重用相同的 interp_line!
}
```

**关键**: 面度量插值使用与坐标导数相同的 `interp_line` 函数，保证离散一致性。

**检查项**:
- [ ] 面数组维度是否正确：`(nci+1, ncj, nck)` 等
- [ ] 对 `uniform` 网格是否跳过插值直接赋常值
- [ ] 子块情况下是否从 donor `extract_metrics_from` 萃取面度量（跳过插值）

---

## 阶段 9: Uniform 网格路径

**文件**: [include/grid/grid.hxx:780-903](include/grid/grid.hxx#L780-L903)

**触发**: `metrics_type == "uniform"`

### Cell 中心度量

```cpp
// 从 interior 推断 dx, dy, dz
Int i0 = ng + nci / 4;  // 取远离边界的 interior cell
Real dx = cell_x(i0+1, j0, k0) - cell_x(i0, j0, k0);
Real dy = cell_y(i0, j0+1, k0) - cell_y(i0, j0, k0);
Real dz = cell_z(i0, j0, k0+1) - cell_z(i0, j0, k0);

// 常值赋值
Real xi_area  = dy * dz;    // |S_ξ|
Real eta_area = dx * dz;    // |S_η|
Real zeta_area = dx * dy;   // |S_ζ|
Real J = dx * dy * dz;

FOR_ALL_CELLS:
    met_xi_x = xi_area;    met_xi_y = 0;    met_xi_z = 0;
    met_eta_x = 0;         met_eta_y = eta_area; met_eta_z = 0;
    met_zeta_x = 0;        met_zeta_y = 0;   met_zeta_z = zeta_area;
    jacobian = J;
```

### 面度量

直接赋常值，不插值:
```cpp
face_xi_x = xi_area;   // 所有面相同
// ...
```

**检查项**:
- [ ] dx, dy, dz 是否正确（等间距网格应处处相等）
- [ ] Jacobian 是否 = dx*dy*dz

---

## 完整流水线时序图

```
main.cpp / ParallelManager::initialize
│
├─ [Phase A] 所有 rank 独立读 CGNS
│   ├─ A1. CGNSReader::open(grid_file)
│   ├─ A2. read_all(zones, ng)
│   │   └─ for each zone:
│   │       ├─ Grid::allocate(ni, nj, nk, ng)     → node_x/y/z: (ni,nj,nk) 核心
│   │       ├─ read_coords(B,Z,grid)               → node_x/y/z 填充核心节点
│   │       ├─ grid.compute_cell_centers()          → cell_x/y/z: (nci,ncj,nck) 核心
│   │       ├─ grid.compute_cell_volumes()          → cell_vol: Cartesian 近似
│   │       ├─ read_bc(B,Z,grid)                    → grid.bc
│   │       └─ read_1to1(B,Z,grid)                  → grid.connections
│   │
│   └─ zones[] = [Grid0, Grid1, ...]
│
├─ [Phase B] Ghost 层扩展 (所有 zone)
│   └─ for each zone:
│       └─ z.extend_ghost_layers()
│           ├─ 保存旧 node_x/y/z → old_nx/y/z
│           ├─ 分配扩展 node_x/y/z: (nie, nje, nke) = 核心 + 2*ng
│           ├─ 核心节点带 ng 偏移拷贝到 interior
│           ├─ ni/nj/nk/nci/ncj/nck ← 更新为扩展后维度
│           ├─ fill_ghost_nodes()
│           │   └─ for each face:
│           │       ├─ find_periodic_connection(face)?
│           │       │   ├─ yes → fill_ghost_face_periodic  (拷贝对面 interior)
│           │       │   └─ no  → fill_ghost_face_extrapolate (线性外推)
│           ├─ compute_cell_centers() → cell_x/y/z 扩展后维度
│           └─ compute_cell_volumes() → cell_vol 扩展后维度
│
├─ [Phase C] 跨 zone 界面 ghost 修正 (所有 zone)
│   └─ for each zone:
│       └─ for each face:
│           └─ find_face_connection(face) → donor zone?
│               ├─ donor_name == self? → skip
│               └─ donor_name != self? → fix_interface_ghost(face, donor, conn)
│                   └─ 用 donor interior 节点覆盖 ghost 节点 (减平移)
│
├─ [Phase D] 界面修正后重算 cell center (所有 zone)
│   └─ for each zone:
│       ├─ z.compute_cell_centers()    ← 使用修正后的 ghost 节点！
│       └─ z.compute_cell_volumes()
│
├─ [Phase E] 计算度量系数 (所有 zone)
│   └─ for each zone:
│       ├─ z.compute_metrics()
│       │   ├─ metrics_type == "uniform"?
│       │   │   ├─ yes → compute_metrics_uniform()
│       │   │   └─ no  → SCMM:
│       │   │       ├─ build_face_periodic(fp[6])
│       │   │       ├─ 分配 met_xi/eta/zeta_x/y/z + jacobian: (nci,ncj,nck)
│       │   │       ├─ 计算 9 个坐标导数 (x/y/z_xi/eta/zeta)
│       │   │       │   └─ compute_coord_deriv → InterpDiff::derivative → interp_line + diff_line
│       │   │       ├─ 组装 9 个度量分量 (36 次 accum_deriv)
│       │   │       │   └─ mul → accum_deriv → InterpDiff::derivative → interp_line + diff_line
│       │   │       └─ 计算 Jacobian (3 次 S_dir → derivative → 组装)
│       │   └─ z.compute_face_metrics()
│       │       ├─ metrics_type == "uniform"?
│       │       │   ├─ yes → 常值面度量
│       │       │   └─ no  → interp_to_faces (9 次, 每个分量每方向)
│       │       │            └─ 重用 interp_line, 与坐标导数相同 stencil
│       │       └─ 分配 face_xi/eta/zeta_x/y/z
│       │
│       └─ (结果存入每个 zone 的 Grid 对象中)
│
├─ [Phase F] 域分解 (所有 rank 相同结果)
│   └─ BlockDecomposer::decompose(zones, nprocs, ng)
│       → std::vector<SubBlock> decomp
│
├─ [Phase G] 构建本地 LocalBlock (仅本 rank 的 block)
│   └─ for each SubBlock assigned to my_rank:
│       ├─ zone 未剖分:
│       │   └─ LocalBlock::from_full_zone(zone, block_id, ...)
│       │       └─ lb.grid = full_zone (deep copy — 包含已计算的度量)
│       │
│       └─ zone 已剖分:
│           └─ LocalBlock::from_sub_zone(zone, sub, block_id, ...)
│               ├─ 从 full_zone 按偏移提取 node_x/y/z
│               ├─ 从 full_zone 提取 periodic connectivity
│               ├─ compute_cell_centers() / compute_cell_volumes()
│               └─ ★ extract_metrics_from(full_zone, ci0, cj0, ck0)
│                   ├─ 复制 9 个 cell 中心度量
│                   ├─ 复制 Jacobian
│                   ├─ 复制 9 个面度量
│                   └─ 复制 cell_vol
│
└─ [Phase H] 设置 HaloExchange / FluxHaloExchange
    └─ setup(block) → 预分配通讯 buffer
```

---

## 变量生命周期汇总表

| 变量 | 首次创建 | 维度 | 更新节点 | 最终状态 |
|------|---------|------|---------|---------|
| `node_x/y/z` | Phase A2: alloc (核心) | (ni,nj,nk) | Phase B: 扩展 + ghost填充; Phase C: interface修正; Phase G: 子块提取 | Z: 完整zone全部正确; L: 子块依赖full_zone ghost |
| `cell_x/y/z` | Phase A2: 核心 | (nci,ncj,nck) | Phase B: 扩展后重算; Phase D: interface修正后重算; Phase G: 子块重算 | Z: 含正确ghost; L: 从full_zone节点计算 |
| `cell_vol` | Phase A2: Cartesian近似 | (nci,ncj,nck) | Phase B/D: 重算; Phase G: extract覆盖 | Z: Cartesian近似; L: 从full_zone萃取 |
| `met_xi/eta/zeta_x/y/z` | Phase E: SCMM | (nci,ncj,nck) | Phase G: extract覆盖(子块) | 全部cell(含ghost) |
| `jacobian` | Phase E: SCMM | (nci,ncj,nck) | Phase G: extract覆盖(子块) | 全部cell |
| `face_xi/eta/zeta_x/y/z` | Phase E: 面插值 | (nci+1,ncj,nck)等 | Phase G: extract覆盖(子块) | 全部面(含ghost面) |
| `x_xi/eta/zeta` | Phase E: 临时 (stack) | (nci,ncj,nck) | 仅Phase E内使用，用完析构 | — |
| `t_prod` | Phase E: 临时 | (nci,ncj,nck) | 每项mul/accum_deriv复用 | — |
| `t_deriv` | Phase E: 临时 | (nci,ncj,nck) | 每项accum_deriv复用 | — |

---

## 人工检查清单

### A. 网格读取

- [ ] **A1** CGNS 文件中的 `CoordinateX/Y/Z` 是否全部读取成功（无 fallback 到 `RealSingle` 的警告）
- [ ] **A2** `cell_x/y/z` 初始值：对均匀 Cartesian 网格，相邻 cell 中心间距应为常数
- [ ] **A3** `cell_vol` 初始 Cartesian 近似：全部为正值
- [ ] **A4** 1-to-1 connectivity 的 `transform[3]` 是否为 `{1,2,3}`（无旋转）
- [ ] **A5** 周期 connectivity 的 `translation[3]` 是否与几何一致

### B. Ghost 扩展

- [ ] **B1** `extend_ghost_layers` 后 `ni/nj/nk` = `ni_core + 2*ng`
- [ ] **B2** `nci/ncj/nck` = `ni-1` 等
- [ ] **B3** 核心节点拷贝偏移 = `ng` — 确认 interior 节点坐标与原始核心节点一致
- [ ] **B4** 周期面 ghost 节点：验证 `ghost = donor_interior - translation` 给出合理的物理坐标
- [ ] **B5** 非周期面 ghost：线性外推生成的节点不会与真实几何偏差过大（对近正交网格通常 OK）

### C. Interface 修正

- [ ] **C1** 所有跨 zone 的 1-to-1 连接都被 `fix_interface_ghost` 处理（无 `donor_name != z.name` 遗漏）
- [ ] **C2** `cur_dim == donor_dim` 对所有连接成立（否则被跳过，ghost 为外推值）
- [ ] **C3** interface 修正后 `compute_cell_centers` 和 `compute_cell_volumes` 都被调用

### D. 坐标导数

- [ ] **D1** `dh = 1.0`（计算空间步长）
- [ ] **D2** `x_xi` 等 9 个导数数组全部非零（无全零导数）
- [ ] **D3** 对均匀 Cartesian 网格，`x_xi` 应为常数 ≈ dx, `x_eta` ≈ 0, `x_zeta` ≈ 0
- [ ] **D4** 对曲面网格，ghost cell 附近的导数应与 interior 连续（无突变）

### E. SCMM 度量

- [ ] **E1** 9 个度量分量 + Jacobian 全部为有限值（无 NaN/Inf）
- [ ] **E2** Jacobian 全部为正值（右手系网格）
- [ ] **E3** 度量满足 GCL：对常数场，`∂_ξ(ξ̂_x) + ∂_η(η̂_x) + ∂_ζ(ζ̂_x)` 应 ≈ 0（机器精度）
- [ ] **E4** 对均匀 Cartesian 网格，SCMM 结果应与 uniform 路径一致
- [ ] **E5** 周期边界附近度量连续（周期面 ghost 数据正确 → SCMM stencil 正确）

### F. 面度量

- [ ] **F1** 面数组维度：`face_xi_*` 为 `(nci+1, ncj, nck)` 等
- [ ] **F2** 面度量在 cell 内点与周围 cell 中心度量的插值一致（`interp_line` 与 `derivative` 使用相同算子）
- [ ] **F3** 面度量在周期边界处连续（周期面 ghost cell 度量正确 → 插值正确）

### G. 子块度量萃取

- [ ] **G1** `extract_metrics_from` 的 `ci0/cj0/ck0` 偏移正确对应子块在完整 zone 中的位置
- [ ] **G2** 子块的度量、Jacobian、面度量、cell_vol 全部从 donor 萃取
- [ ] **G3** 萃取后子块的 interior 度量与完整 zone 对应位置一致（逐 cell 对比）
- [ ] **G4** 子块的 ghost 度量与完整 zone 对应位置一致

### H. 均匀自由流验证

- [ ] **H1** 在曲面网格上用均匀场 (`rho=const, u=v=w=const, p=const`) 初始化
- [ ] **H2** 残差初始值应为 ≈ 0（机器精度） — 度量 GCL 满足则对流项无伪源项
- [ ] **H3** 多步运行后密度/压力的最大偏差在 1e-12 量级

---

*相关文档: [metric.md](metric.md), [ghost.md](ghost.md), [mpi.md](mpi.md)*
