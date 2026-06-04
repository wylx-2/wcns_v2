# WCNS v2 — 主时间循环实现计划

## 概述

本文档描述主时间循环的实现方案，涵盖时间步长计算、显式 Runge-Kutta 子步循环（含无粘/粘性通量计算）、收敛判断、以及输出功能。

### 当前已完成模块（上下文）

| 模块 | 状态 | 说明 |
|------|------|------|
| Grid + 度量系数 + GCL | ✅ 完成 | cell-center 度量、face 度量、Jacobian，GCL 验证通过 |
| Field (prim/cons/flux/RHS) | ✅ 完成 | 原始/守恒变量转换、通量数组、RHS 数组、扩展字段 |
| InterpDiff | ✅ 完成 | 6阶中心+5阶单边插值/差分，`ah[i]=a_{i-1/2}` 约定 |
| FlowInitializer | ✅ 完成 | uniform、poiseuille 初始场 |
| BoundaryConditionApplier | ✅ 完成 | farfield/wall/symmetry/inflow/outflow + 边/角 ghost |
| HaloExchange | ✅ 完成 | MPI 并行 halo 交换 |
| Config + ConfigReader | ✅ 完成 | 参数管理 |

### 总体数据流（单时间步）

```
compute_dt()  →  for m in 1..N_stage:
                    inviscid_rhs()    (插值 → Riemann → 差分)
                  + viscid_rhs()      (插值 → 梯度 → 粘性通量 → 插值 → 差分)
                  + body_force()       (体积力源项 → 累加入 RHS)
                  → RHS = inv + vis + body_force
                  → Q^(m) = α_m*Q^(0) + β_m*(Q^(m-1) + Δt*RHS)
                  → BC apply + Halo exchange
                  → [final stage] residual(5-comp L2) + monitor(T/p/ρ min/max)
                  → convergence/divergence check → log residual.dat
              → solution output (Tecplot .plt, interior prim vars) / restart (binary .bin, full cons vars) (every N steps)
```

---

## 子任务拆解

### 子任务 1: 时间步长计算

**目标**：根据 CFL 条件和当地谱半径计算全局/当地时间步长。

**内容**：
- 对每个 cell 计算各方向的谱半径（对流项 + 粘性项贡献）
- 对流谱半径：`Λ_c = |u·n| + c·|n|`（各方向面法向）
- 粘性谱半径：`Λ_v = max(4/3,γ)/ρ * μ/Pr * |n|² / V`
- 取各方向最小 Δt_local，全局 min 得到 Δt_global
- `Δt = CFL * min(Δt_global)` 或 local time stepping

**涉及文件**：
- 新文件 `include/wcns_v2/time/time_step.h` / `.hxx`

**依赖**：Grid（度量系数、Jacobian）、Field（primitive vars）、Config（CFL, Re, Pr）

---

### 子任务 2: WCNS 非线性插值类

**目标**：实现 WCNS 原始非线性插值，将 cell-center 处的守恒变量（或特征变量）插值到半节点（face），得到左右值 `Q_L`, `Q_R`。这是 Riemann 求解器的前置步骤。

**核心算法**（WCNS 5 阶 WENO 非线性插值）：

对每个半节点 `i+1/2`（face），沿方向逐线对每个守恒/特征分量独立进行插值。

**左值 `Q_L(i+1/2)` — 左偏模板 `[i-2, i-1, i, i+1, i+2]`**：

3 个子模板及其 3 阶插值（到半节点 `i+1/2`）：
```
S0 = {i-2, i-1, i}:    v0 = ( 3·a[i-2] - 10·a[i-1] + 15·a[i]  ) / 8
S1 = {i-1, i, i+1}:    v1 = ( -1·a[i-1] +  6·a[i]   +  3·a[i+1]) / 8
S2 = {i, i+1, i+2}:    v2 = (  3·a[i]   +  6·a[i+1] -  1·a[i+2]) / 8
```

光滑度指示器（Jiang–Shu 1996）：
```
β0 = 13/12·(a[i-2] - 2·a[i-1] + a[i])² + 1/4·(a[i-2] - 4·a[i-1] + 3·a[i])²
β1 = 13/12·(a[i-1] - 2·a[i] + a[i+1])² + 1/4·(a[i-1] - a[i+1])²
β2 = 13/12·(a[i] - 2·a[i+1] + a[i+2])² + 1/4·(3·a[i] - 4·a[i+1] + a[i+2])²
```

非线性权（WENO-JS，p=2，ε=1e-6）：
```
α_k = d_k / (ε + β_k)²      线性最优权: d0=1/16, d1=10/16, d2=5/16
ω_k = α_k / (α_0 + α_1 + α_2)
```

左值：`Q_L = ω_0·v0 + ω_1·v1 + ω_2·v2`

**右值 `Q_R(i+1/2)` — 右偏模板 `[i-1, i, i+1, i+2, i+3]`**（对称）：

3 个子模板（对称翻转）：
```
S0' = {i+3, i+2, i+1}:  v0' = ( 3·a[i+3] - 10·a[i+2] + 15·a[i+1]) / 8
S1' = {i+2, i+1, i}:    v1' = (-1·a[i+2] +  6·a[i+1] +  3·a[i]  ) / 8
S2' = {i+1, i, i-1}:    v2' = ( 3·a[i+1] +  6·a[i]   -  1·a[i-1]) / 8
```

光滑度指示器和权计算与左值对称，得 `Q_R`。

**边界处理**（半节点 h=0,1,2 和 h=n-2,n-1,n）：

左右端点附近模板不足 5 点，无法构造完整 WENO 插值。此时：
- 使用 `InterpDiff` 中已实现的 5 阶单边插值公式
- **将单边插值结果同时赋给 Q_L 和 Q_R**（`Q_L = Q_R = InterpDiff::interp_*`），确保 Riemann 求解器在边界处正常工作

| 半节点 | 索引 h | 使用公式 | 模板 |
|--------|--------|---------|------|
| `a_{-1/2}` | 0 | `interp_left_1st(&a[0])` | a[0..4] |
| `a_{1/2}` | 1 | `interp_left_2nd(&a[0])` | a[0..4] |
| `a_{3/2}` | 2 | `interp_left_3rd(&a[0])` | a[0..4] |
| ... | 3..n-3 | WENO 非线性插值 | (Q_L, Q_R 分别计算) |
| `a_{n-5/2}` | n-2 | `interp_right_3rd(&a[n-5])` | a[n-5..n-1] |
| `a_{n-3/2}` | n-1 | `interp_right_2nd(&a[n-5])` | a[n-5..n-1] |
| `a_{n-1/2}` | n | `interp_right_1st(&a[n-5])` | a[n-5..n-1] |

**输出数组大小**（以 x 方向为例）：
- cell-center 数组：`(0 ~ nx+5) × (0 ~ ny+5) × (0 ~ nz+5)`，即 `(nci) × (ncj) × (nck)`，其中 `nci = nx+6`
- x 方向 face 数组（Q_L, Q_R）：**`(0 ~ nx+6) × (0 ~ ny+5) × (0 ~ nz+5)`**，即 `(nci+1) × ncj × nck`
- 类似，y 方向 face 数组：`nci × (ncj+1) × nck`；z 方向 face 数组：`nci × ncj × (nck+1)`

**设计要点**：

1. **类层次**（预留扩展）：
   ```
   WcnsInterpBase          ← 抽象基类，定义虚接口
      ├── WcnsNonlinearInterp   ← 当前实现：WCNS 原始非线性插值（WENO-JS）
      ├── (预留) WcnsLinearInterp     ← 线性 6 阶中心插值（无 WENO 权）
      ├── (预留) WcnsWenoZInterp      ← WENO-Z 改进格式
      └── (预留) WcnsWenoJSInterp     ← WENO-JS 标准格式（可参数化 p, ε）
   ```
   - 当前 `WcnsNonlinearInterp` 使用**静态方法**（无状态，所有参数来自 Config）
   - 基类 `WcnsInterpBase` 定义纯虚接口，后续派生类可覆盖

2. **Config 扩展**：新增 `interp_vars` 选项控制插值变量空间
   - `interp_vars = "conservative"`（默认）— 对守恒变量 `(rho, rhou, rhov, rhow, rhoE)` 逐分量插值
   - `interp_vars = "characteristic"`（预留）— 投影到特征空间插值后投影回来
   - 当前仅实现 `conservative` 模式，`characteristic` 选择时抛出异常提示未实现

3. **分量循环**：对 5 个守恒分量分别进行相同的插值操作（分量间独立）

4. **逐方向接口**：
   - `interp_xi(lb, cfg, Q_L, Q_R)` — 沿 ξ(i) 方向对所有 j,k 线进行插值
   - `interp_eta(lb, cfg, Q_L, Q_R)` — 沿 η(j) 方向对所有 i,k 线进行插值
   - `interp_zeta(lb, cfg, Q_L, Q_R)` — 沿 ζ(k) 方向对所有 i,j 线进行插值
   - 每个函数内：提取 1D 线 → 对每个分量调用 1D WENO 插值 → 写回 Q_L/Q_R

5. **Q_L/Q_R 存储**：采用方案 A，在 `Field` 中预分配 6 组 face `ConservativeVars`：
   ```
   ql_xi, qr_xi   — ξ-face (nci+1)×ncj×nck
   ql_eta, qr_eta — η-face nci×(ncj+1)×nck
   ql_zeta, qr_zeta — ζ-face nci×ncj×(nck+1)
   ```
   与 `inv_xi/vis_xi` 等通量数组尺寸一致，在 `Field::allocate()` 中统一分配。避免每时间步重复申请内存。

**涉及文件**：
- 新文件 `include/wcns_v2/scheme/wcns_interp.h` — 类声明（基类 + 非线性派生类）
- 新文件 `include/wcns_v2/scheme/wcns_interp.hxx` — 内联实现（WENO 权、1D 插值核、3D 逐方向接口）
- 修改 `include/wcns_v2/field/field.h` — 新增 `ql_xi/qr_xi`、`ql_eta/qr_eta`、`ql_zeta/qr_zeta`（6 组 face ConservativeVars）
- 修改 `include/wcns_v2/field/field.hxx` — 新增 Q_L/Q_R 分配逻辑
- 修改 `include/wcns_v2/core/config.h` — 新增 `interp_vars` 参数
- 修改 `include/wcns_v2/io/config_reader.h` — 新增 `set_scheme` 方法声明
- 修改 `src/io/config_reader.cpp` — 新增 `[scheme]` section 和 `interp_vars` 解析
- 修改 `input.ini` — 新增 `[scheme]` section
- 修改 `include/wcns_v2/scheme/interp_diff.h` — 新增 `friend class WcnsNonlinearInterp`（访问单边插值系数）

**1D 插值核伪代码**：
```cpp
// 对一条长度为 n 的线进行 WCNS 非线性插值
// a[0..n-1] : 输入的 cell-center 值（含 ghost cell）
// ql[0..n]  : 输出的左值（n+1 个半节点）
// qr[0..n]  : 输出的右值（n+1 个半节点）
void wcns_interp_1d(const Real* a, Real* ql, Real* qr, Int n) {
    for (Int h = 0; h <= n; ++h) {
        // 左边界：单边插值，Q_L = Q_R
        if (h == 0) { ql[h] = qr[h] = InterpDiff::interp_left_1st(&a[0]); continue; }
        if (h == 1) { ql[h] = qr[h] = InterpDiff::interp_left_2nd(&a[0]); continue; }
        if (h == 2) { ql[h] = qr[h] = InterpDiff::interp_left_3rd(&a[0]); continue; }

        // 右边界：单边插值，Q_L = Q_R
        if (h == n)     { ql[h] = qr[h] = InterpDiff::interp_right_1st(&a[n-5]); continue; }
        if (h == n - 1) { ql[h] = qr[h] = InterpDiff::interp_right_2nd(&a[n-5]); continue; }
        if (h == n - 2) { ql[h] = qr[h] = InterpDiff::interp_right_3rd(&a[n-5]); continue; }

        // 内点：WENO 非线性插值
        // Q_L: 左偏模板 (i = h-1, 使用 a[h-3..h+1])
        ql[h] = weno_left(a, h - 1);

        // Q_R: 右偏模板 (i = h-1, 使用 a[h-2..h+2])
        qr[h] = weno_right(a, h - 1);
    }
}
```

**依赖**：InterpDiff（边界单边插值系数）、Field（cons 数组）、Grid（索引信息）、Config（interp_vars）

---

### 子任务 3: Riemann 求解器类

**目标**：在 face 处根据左右状态 `Q_L`, `Q_R` 计算数值通量。

**核心算法**（Roe 格式）：
- 对每个 face 计算 Roe 平均状态
- 构造特征值和特征向量
- 通量：`F = 1/2 * (F_L + F_R) - 1/2 * |A| * (Q_R - Q_L)`
- 其中 `|A| = R * |Λ| * R^{-1}`，`Λ` 为特征值对角阵
- 需处理特征值过零时的熵修正（Harten 校正）

**设计要点**：
- 独立类 `RiemannSolver`，便于后续添加 HLLC、AUSM+ 等格式
- 基类 `RiemannSolverBase` 定义虚接口
- 各方向（ξ, η, ζ）的 face-normal 投影速度处理方式相同，通过传入度量系数区分

**涉及文件**：
- 新文件 `include/wcns_v2/scheme/riemann_solver.h` / `.hxx`

**依赖**：Field（通量数组 face 维度）、Grid（face 度量系数）

---

### 子任务 3.5: 界面通量通讯（Face Flux Halo Exchange）

**目标**：在 Riemann 求解器完成后，对周期边界和并行连接等 connectivity 界面处的半节点通量进行通讯，将内部计算得到的正确通量值传递到相邻 block 对应的“虚拟网格”（ghost face）位置，确保交界面附近差分使用的通量值一致且正确。

**背景**：Riemann 求解器在各 block 上独立计算所有 face 的通量，包括 ghost 区域的 face。但 ghost 区域 face 的左右插值依赖于 ghost cell 数据（由 cell-centered halo exchange 填充），如果 cell-centered halo exchange 未能正确填充 ghost cell（例如同进程多 block 的 local_copy 未实现），则：
- 近交界面 face 的 `Q_L`/`Q_R` 插值不准确 → 通量值错误
- 进而导致差分时产生 NaN/Inf

**解决方案**：不做 cell-centered halo exchange，而是直接将 Riemann 求解后得到的、位于 block 内部（interior）的正确 face 通量值，通过 MPI 或直接内存拷贝，覆盖到相邻 block 的 ghost face 位置上。

**核心对应关系**（以 J 方向 JMAX↔JMIN 界面为例）：

```
Block 0 (左)                          Block 1 (右)
cell: ... 31  32  33  34 | 35  36  37    cell:  0   1   2 |  3   4   5 ...
            [--- interior ---] [-- ghost]       [-- ghost --] [--- interior ---]
face: ... 32  33  34  35  36  37  38     face:  0   1   2   3   4   5   6 ...
            [Block 0 内部面] [ghost面]          [ghost面] [Block 1 内部面]
```

共享界面 face：**Block 0 face 35 = Block 1 face 3**（同一物理位置的半节点）

通量对应关系（Block 0 内部 face → Block 1 ghost face）：
```
Block 0 face 32  ←→  Block 1 face 0
Block 0 face 33  ←→  Block 1 face 1
Block 0 face 34  ←→  Block 1 face 2
Block 0 face 35  ←→  Block 1 face 3  (共享面)
```

反向（Block 1 内部 face → Block 0 ghost face）：
```
Block 1 face 3  ←→  Block 0 face 35  (共享面)
Block 1 face 4  ←→  Block 0 face 36
Block 1 face 5  ←→  Block 0 face 37
Block 1 face 6  ←→  Block 0 face 38
```

**通用公式**（以 JMAX/JMIN 连接为例，`ncj = ny_core + 2*ng`）：
- 需交换的 face 数量：`ng + 1`（共享面 + 每侧各 ng 个面）
- Block 0 (JMAX) 发送：inv_eta 的 face `[ncj - 2*ng, ncj - ng]` 共 ng+1 个 → 对应 Block 1 的 face `[0, ng]`
- Block 1 (JMIN) 发送：inv_eta 的 face `[ng, 2*ng]` 共 ng+1 个 → 对应 Block 0 的 face `[ncj - ng, ncj]`

**三个方向的总表**（以 `nc = nc_dim + 2*ng` 表示该方向的 cell 总数，`nc = core + 2*ng`）：

| 方向 | 连接类型 | 发送方 (face 范围) | 接收方 (face 范围) | 数量 |
|------|---------|-------------------|-------------------|------|
| ξ (IMIN/IMAX) | 周期性 / MPI 内部连接 | src: `[0, ng]` 或 `[nci-ng, nci]` | dst: 对应 ghost face | `ng+1` |
| η (JMIN/JMAX) | 周期性 / MPI 内部连接 | src: `[0, ng]` 或 `[ncj-2*ng, ncj-ng]` | dst: 对应 ghost face | `ng+1` |
| ζ (KMIN/KMAX) | 周期性 / MPI 内部连接 | src: `[0, ng]` 或 `[nck-2*ng, nck-ng]` | dst: 对应 ghost face | `ng+1` |

**涉及的通量数组**（每个方向 5 个分量）：
- `inv_xi.f1, f2, f3, f4, f5`（ξ 方向）
- `inv_eta.f1, f2, f3, f4, f5`（η 方向）
- `inv_zeta.f1, f2, f3, f4, f5`（ζ 方向）

**实现方式**：

1. **同进程多 block（local copy）**：直接内存拷贝，face 之间有固定的偏移关系
2. **跨进程（MPI）**：使用 `MPI_Isend`/`MPI_Irecv` 发送 face 子集（ng+1 个 face × 横向维度 × 5 分量）

**设计要点**：
- 新增 `FluxHaloExchange` 类（或在现有 `HaloExchange` 中扩展 face 维度支持）
- 复用 `LocalBlock::neighbors` 中的连通性信息（target_rank, target_block, target_face）
- 实现时需要正确处理 face 尺寸差异（ξ/η/ζ 方向 face 数组维度互不相同）

**涉及文件**：
- 新文件 `include/wcns_v2/parallel/flux_halo_exchange.h` / `.hxx`（或扩展现有 `halo_exchange.h`）
- 修改 `include/wcns_v2/parallel/parallel_manager.h` / `.hxx`（新增 `exchange_flux_halos` 方法）

**依赖**：子任务 3 (RiemannSolver)、子任务 4 (InviscidRHS)、LocalBlock（neighbors 信息）

---

### 子任务 4: 无粘通量 RHS 计算

**目标**：串联插值 → Riemann 求解 → 显式差分，完成一个方向的 RHS_inv 计算。

**流程**（以 ξ 方向为例）：
1. **插值**：`WcnsInterp::interp_xi(lb, Q_L, Q_R)` — 得到面左右值，数组大小 `(0~nx+6) × (0~ny+5) × (0~nz+5)`
2. **Riemann**：`RiemannSolver::solve_xi(lb, Q_L, Q_R, inv_xi_flux)` — 计算面通量
3. **差分**：对每个通量分量，调用 `InterpDiff::derivative(inv_xi_flux.f1, rhs_inv.rho, dir=0, dh=1, ...)` — 注意差分只在内点 `(3~nx+2)×(3~ny+2)×(3~nz+2)` 上计算，边界不需要

**要点**：
- Q_L/Q_R 临时数组大小：x方向面 `(nci+1)×ncj×nck`（x方向总cell为nx+6=ni，面数量为ni+1，即nx+7）
- inv_xi_flux 已预分配为 `(nci+1, ncj, nck)` 在 Field 中
- 差分使用 `InterpDiff::derivative`（内部调用 `interp_line` + `diff_line`），仅需内点值
- 三个方向依次处理后累积到 RHS_inv

**涉及文件**：
- 新文件 `include/wcns_v2/scheme/inviscid_flux.h` / `.hxx`

**依赖**：子任务 2 (WcnsInterp)、子任务 3 (RiemannSolver)、InterpDiff

---

### 子任务 5: 粘性通量 RHS 计算

**目标**：计算粘性应力张量和热通量贡献的 RHS_vis。

**流程**：

**5a. 半节点物理量插值**：
- 将 cell-center 的 `(u, v, w, T)` 插值到各方向的半节点
- 调用 `InterpDiff::interp_to_faces`，包括单边插值
- 温度从 `T = p/(rho * R_gas_nd)` 计算（无量纲下 `R_gas_nd = 1/(gamma*Mach^2)`，即 `T = gamma*Mach^2 * p/rho`）
- **输出存储**：4 物理量 × 3 方向 = 12 个 face 数组，在 `Field::allocate()` 中预分配，命名与尺寸：

  | 数组 | ξ-face `(nci+1)×ncj×nck` | η-face `nci×(ncj+1)×nck` | ζ-face `nci×ncj×(nck+1)` |
  |------|--------------------------|--------------------------|--------------------------|
  | `u_face_*` | `u_face_xi` | `u_face_eta` | `u_face_zeta` |
  | `v_face_*` | `v_face_xi` | `v_face_eta` | `v_face_zeta` |
  | `w_face_*` | `w_face_xi` | `w_face_eta` | `w_face_zeta` |
  | `T_face_*` | `T_face_xi` | `T_face_eta` | `T_face_zeta` |

  （每个 face 数组为单分量 `MultiArray3D<Real>`，与 `ql_xi` 等共存于 `Field` 类中）

**5b. 单元中心物理量导数计算**（5 步流程）：

**步骤 1 — 形成半节点乘积**：
- 将 5a 插值得到的 `(u,v,w,T)` 半节点值与对应方向的 face 度量系数相乘
- 以 ξ 方向为例，对速度分量 u：
  ```
  prod_u_xi_x(i+1/2) = u_face_xi * face_xi_x    （u * ξ̂_x 在 ξ-face 半节点）
  prod_u_xi_y(i+1/2) = u_face_xi * face_xi_y    （u * ξ̂_y 在 ξ-face 半节点）
  prod_u_xi_z(i+1/2) = u_face_xi * face_xi_z    （u * ξ̂_z 在 ξ-face 半节点）
  ```
- 对 v, w, T 同理。每物理量每方向产生 3 个乘积数组
- **逐方向处理**：一次只处理一个方向（先 ξ，再 η，再 ζ），同时仅需 **12 个临时 face 数组**（4 物理量 × 3 度量分量），处理完一个方向后在下一方向复用
- 临时乘积数组大小：ξ 方向 `(nci+1)×ncj×nck`，η 方向 `nci×(ncj+1)×nck`，ζ 方向 `nci×ncj×(nck+1)`

**步骤 2 — 半节点乘积数据交换**：
- 对周期边界或并行连接（connectivity）区域，将内部正确的半节点乘积值交换到 ghost face 位置
- 逻辑与子任务 3.5 `FluxHaloExchange` 相同：`ng+1` 个 face 切片在连通边界上交换
- **需要扩展 `FluxHaloExchange`**：当前只支持 5 分量 `FluxVars`（`inv_xi/eta/zeta`），需新增单分量 face 数组交换方法：
  ```cpp
  void exchange_face_array(MultiArray3D<Real>& face_arr, int dir,
                           const LocalBlock& block,
                           const std::vector<LocalBlock>& all_blocks);
  ```
- 每方向交换 12 个单分量 face 数组，建议将同一方向的 12 次交换数据合并为一次 MPI 通信以减少延迟

**步骤 3 — 半节点值差分为 cell-center 导数**：
- **需要新增 `InterpDiff::derivative_from_faces` 公开方法**

  原因：`InterpDiff::derivative` 接受 cell-center 输入、内部执行 `interp_line` → `diff_line` 两步，无法直接消费步骤 1-2 已产生的半节点乘积
- 新方法签名：
  ```cpp
  /// 从半节点值直接差分为 cell-center 导数（跳过插值步骤）。
  /// @param af   输入半节点数组（face，大小取决于 dir：ξ→(nci+1)×nj×nk 等）
  /// @param da   输出导数数组（cell-center, nci×ncj×nck），须预分配
  /// @param dir  方向: 0=ξ, 1=η, 2=ζ
  /// @param dh   计算空间网格间距（=1.0）
  /// @param face_is_periodic  各面周期性标志 [IMIN,IMAX,JMIN,JMAX,KMIN,KMAX]
  static void derivative_from_faces(const MultiArray3D<Real>& af,
                                     MultiArray3D<Real>& da,
                                     int dir, Real dh,
                                     const bool face_is_periodic[6]);
  ```
- 内部实现：仅执行 `diff_line`（6 阶中心 + 5 阶单边差分），不执行 `interp_line`
- 3D 遍历逻辑与现有 `derivative` 相同：ξ 方向按 `(j,k)` 遍历 1D 线；η/ζ 方向逐线 gather → `diff_line` → scatter
- 计算范围：全 cell `(nci)×(ncj)×(nck)`，即 `0..nci-1, 0..ncj-1, 0..nck-1`（含 ghost cell）

**步骤 4 — 组合链式法则**：
- 对每个物理量 φ ∈ {u,v,w,T}，三个方向导数累加后乘以 `ja_inv`：
  ```
  dφ/dx = ja_inv * ( d(φ*ξ̂_x)/dξ + d(φ*η̂_x)/dη + d(φ*ζ̂_x)/dζ )
  dφ/dy = ja_inv * ( d(φ*ξ̂_y)/dξ + d(φ*η̂_y)/dη + d(φ*ζ̂_y)/dζ )
  dφ/dz = ja_inv * ( d(φ*ξ̂_z)/dξ + d(φ*η̂_z)/dη + d(φ*ζ̂_z)/dζ )
  ```
- 其中 `ja_inv(i,j,k) = 1.0 / jacobian(i,j,k)`
- 得到 **12 个梯度数组**：`du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dT_dx, dT_dy, dT_dz`，大小均为 `(nci, ncj, nck)`
- 梯度数组作为 `ViscidRHS` 的成员变量（首次调用时分配，后续复用，避免每时间步重复 new/delete）

**步骤 5 — 梯度数组 cell-center 数据交换**：
- 使用现有 `HaloExchange::exchange_multi` 对 12 个梯度数组进行 cell-center halo 交换
- **目的**：填充 connectivity 边界 ghost cell 位置的梯度值。5c 应力计算需要所有 cell（含 ghost）的梯度，因为 5d 插值应力到 face 时需要 ghost cell 的应力值作为插值模板支撑
- 实现示例：
  ```cpp
  std::vector<MultiArray3D<Real>*> grads = {
      &du_dx, &du_dy, &du_dz, &dv_dx, &dv_dy, &dv_dz,
      &dw_dx, &dw_dy, &dw_dz, &dT_dx, &dT_dy, &dT_dz
  };
  halo_exch.exchange_multi(grads, lb);
  ```

**5c. 计算 cell-center 粘性通量（Cartesian 坐标系）**：

**设计决策**：τ（应力张量）和 q（热通量）是中间量，5d/5e 后续操作的对象是粘性通量本身。直接在 cell center 计算并存储 Cartesian 坐标系下的三个粘性通量向量 `F_vis`, `G_vis`, `H_vis`（各 5 分量 = `FluxVars`）。τ 和 q 不单独存储。

**粘性系数配置**（`Config` 新增字段）：

```cpp
// 粘性计算模式
// "none"       — 无粘模式，跳过所有粘性计算，RHS_vis = 0
// "constant"   — 固定粘性系数，μ = μ_const / Re
// "sutherland" — Sutherland 公式（温度相关）
std::string viscous_type = "constant";

// 固定粘性系数乘子（仅 viscous_type="constant" 时使用，通常=1.0）
Real mu_const = 1.0;

// Sutherland 常数 [K]（仅 viscous_type="sutherland" 时使用）
Real sutherland_S = 110.4;   // 维度值，finalize() 中转为无量纲
Real sutherland_S_nd;        // S* = S / T_ref（由 finalize() 计算）
```

**粘性系数 μ 计算**（无量纲，在所有 cell 上进行）：

| 模式 | 公式 |
|------|------|
| `none` | μ = 0 |
| `constant` | μ = μ_const / Re |
| `sutherland` | μ = T^(1.5) · (1 + S*) / (T + S*) / Re |

温度 `T = γ · Mach² · p / ρ`（与 5a 中 `compute_temperature` 公式一致）。

**导热系数 k**（无量纲 Fourier 热通量 `q_i = -k · ∂T/∂x_i`）：

无量纲下气体常数 `R_gas_nd = 1/(γ·Mach²)`，等压比热 `cp* = γ/(γ-1) · R_gas_nd = 1/((γ-1)·Mach²)`。因此：

```
k* = μ* · cp* / Pr = μ* / ((γ-1) · Mach² · Pr)
```

**应力张量**（中间量，不存储，直接在通量组装时使用）：

```
divV = ∂u/∂x + ∂v/∂y + ∂w/∂z

τ_xx = μ · (2·∂u/∂x - 2/3·divV)
τ_yy = μ · (2·∂v/∂y - 2/3·divV)
τ_zz = μ · (2·∂w/∂z - 2/3·divV)
τ_xy = μ · (∂u/∂y + ∂v/∂x)
τ_xz = μ · (∂u/∂z + ∂w/∂x)
τ_yz = μ · (∂v/∂z + ∂w/∂y)
```

**热通量**（中间量，不存储）：

```
q_x = -k · ∂T/∂x,   q_y = -k · ∂T/∂y,   q_z = -k · ∂T/∂z
```

**Cartesian 粘性通量**（`FluxVars` × 3，cell-center，大小 nci×ncj×nck）：

```
// x-direction: F_vis = [f1, f2, f3, f4, f5]^T
F_vis.f1 = 0                          (质量通量 = 0)
F_vis.f2 = τ_xx                       (x-动量通量)
F_vis.f3 = τ_xy                       (y-动量通量)
F_vis.f4 = τ_xz                       (z-动量通量)
F_vis.f5 = u·τ_xx + v·τ_xy + w·τ_xz - q_x   (能量通量)

// y-direction: G_vis
G_vis.f1 = 0
G_vis.f2 = τ_yx
G_vis.f3 = τ_yy
G_vis.f4 = τ_yz
G_vis.f5 = u·τ_yx + v·τ_yy + w·τ_yz - q_y

// z-direction: H_vis
H_vis.f1 = 0
H_vis.f2 = τ_zx
H_vis.f3 = τ_zy
H_vis.f4 = τ_zz
H_vis.f5 = u·τ_zx + v·τ_zy + w·τ_zz - q_z
```

**Field 新增存储**（cell-center，nci×ncj×nck）：

```cpp
FluxVars vis_x;   // F_vis — Cartesian x-direction viscous flux (f1..f5, f1≡0)
FluxVars vis_y;   // G_vis — Cartesian y-direction viscous flux
FluxVars vis_z;   // H_vis — Cartesian z-direction viscous flux
```

**ViscidRHS 新增方法**：

```cpp
/// 5c — Compute Cartesian viscous flux vectors at cell centers.
///
/// Prerequisites: compute_gradients() + exchange_gradient_halos()
///   must have been called (fills du_dx, du_dy, ..., dT_dz).
///
/// Reads gradients + prim.{u,v,w,p,rho} from Field.
/// Writes vis_x, vis_y, vis_z (FluxVars) in Field at nci×ncj×nck.
///
/// If cfg.viscous_type == "none": sets all components to 0.
static void compute_cell_viscous_flux(LocalBlock& lb, const Config& cfg);
```

**应力/热通量 halo 交换**：

5c 计算完 vis_x/vis_y/vis_z 后，需要交换 ghost cell 的粘性通量值（供 5d 插值到 face 使用）：

```cpp
// ParallelManager 新增方法
void exchange_viscous_flux_halos(std::vector<LocalBlock>& blocks) {
    for (size_t i = 0; i < blocks.size(); ++i) {
        std::vector<MultiArray3D<Real>*> arrays = {
            &blocks[i].field.vis_x.f1, &blocks[i].field.vis_x.f2,
            &blocks[i].field.vis_x.f3, &blocks[i].field.vis_x.f4,
            &blocks[i].field.vis_x.f5,
            &blocks[i].field.vis_y.f1, &blocks[i].field.vis_y.f2,
            &blocks[i].field.vis_y.f3, &blocks[i].field.vis_y.f4,
            &blocks[i].field.vis_y.f5,
            &blocks[i].field.vis_z.f1, &blocks[i].field.vis_z.f2,
            &blocks[i].field.vis_z.f3, &blocks[i].field.vis_z.f4,
            &blocks[i].field.vis_z.f5
        };
        halo_ex_[i].exchange_multi(arrays, blocks[i]);
    }
}
```

**5d. 半节点粘性通量组装（方案 C）**：

**背景分析**：经过对 `HaloExchange` 源码审查确认——`HaloExchange::copy_local` 目前为空实现（TODO），即 cell-center halo exchange 对同进程多 block 的 ghost cell 不执行 local copy。因此，5c 中 `exchange_viscous_flux_halos` 对同进程邻居的 ghost cell 无法保证 `vis_x/y/z` 的正确性。5d 在对 ghost cell 做 `interp_to_faces` 时，ghost cell 的 `vis_x/y/z` 可能错误 → 插值到 ghost face 的 Cartesian 通量错误 → 必须在 face 级别进行交换纠正。

**方案 C 核心思想**：先交换不依赖 face metric 的 Cartesian 通量（物理量，跨 block 一致），再用各 block 自身的 face metric 合成 `vis_xi/eta/zeta`。这与 5b 的 "交换乘积 → 差分 → 组合" 逻辑完全对应。

**实施步骤**（以 ξ 方向为例，η/ζ 对称）：

```
5d-step1（所有 block）: 将 vis_x.{f2,f3,f4,f5}, vis_y.{f2,f3,f4,f5}, vis_z.{f2,f3,f4,f5}
                       （共 12 个 scalar arrays）从 cell-center 插值到 ξ-face
                       → 12 个临时 face 数组，存储在 Field 中
5d-step2（所有 block）: MPI 交换 + 同进程 local copy 交换 12 个 face 数组
5d-step3（所有 block）: 逐 face 点组装 vis_xi:
   vis_xi.f1 = 0
   vis_xi.f2 = vis_x_f2_face * face_xi_x + vis_y_f2_face * face_xi_y + vis_z_f2_face * face_xi_z
   vis_xi.f3 = vis_x_f3_face * face_xi_x + vis_y_f3_face * face_xi_y + vis_z_f3_face * face_xi_z
   vis_xi.f4 = vis_x_f4_face * face_xi_x + vis_y_f4_face * face_xi_y + vis_z_f4_face * face_xi_z
   vis_xi.f5 = vis_x_f5_face * face_xi_x + vis_y_f5_face * face_xi_y + vis_z_f5_face * face_xi_z
```

**关键设计点**：

1. **两步调用模式**：step1 必须在所有 block 上完成后才能开始 step2。原因：step2 的 local copy 需要从邻居 block 的 Field 中读取临时 face 数组——邻居的临时 face 数组必须已填充。这与 5a→5b 的关系相同（5a 在所有 block 上计算 face-interpolated 物理量并存入 Field，5b 再读取）。

   ```cpp
   // main.cpp 调用顺序（每个方向 dir ∈ {0,1,2}）:
   for (auto& lb : blocks)
       ViscidRHS::interp_cart_flux_to_faces(lb, dir, cfg);     // step1
   for (auto& lb : blocks) {
       ViscidRHS::exchange_and_assemble_face_flux(lb, dir,      // step2+3
           blocks, flux_halo_ex, cfg);
   }
   ```

2. **临时 face 数组存储位置**：存储在 `Field` 中，共 12 个 `MultiArray3D<Real>`，按当前方向分配尺寸。这些是纯临时量（仅 5d 使用），步进完成后释放。命名：
   ```cpp
   // Field 新增 — 5d 临时 face 数组（单方向，12 个 scalars）
   MultiArray3D<Real> vis_x_f2_face, vis_x_f3_face, vis_x_f4_face, vis_x_f5_face;
   MultiArray3D<Real> vis_y_f2_face, vis_y_f3_face, vis_y_f4_face, vis_y_f5_face;
   MultiArray3D<Real> vis_z_f2_face, vis_z_f3_face, vis_z_f4_face, vis_z_f5_face;
   ```
   逐方向分配和释放，峰值内存 = 12 个 face 数组（约 `12 × max(ni,nj,nk)² × ng` 量级）。对 128³ 网格约 ~220 MB，可接受。

3. **同进程 local copy**：复用 5b 的 `ViscidRHS::copy_local_face_array` 方法。对 12 个临时 face 数组逐一调用，从邻居 block 的 Field 中读取对应临时 face 数组，写入本 block ghost face 区域。逻辑与 5b 中 `u_face_xi` 的 local copy 完全一致。

4. **face_metric 读取**：组装时使用的 `face_xi_x/y/z` 等 face metric 来自 `lb.grid`。每个 block 使用自己的 face metric——这是方案 C 相对于交换已组装 `vis_xi` 的关键优势（避免了 metric 不一致问题）。

**ViscidRHS 新增方法**：

```cpp
/// 5d-step1 — Interpolate Cartesian viscous flux components to faces.
///
/// Interpolates vis_x.{f2,f3,f4,f5}, vis_y.{f2,f3,f4,f5}, vis_z.{f2,f3,f4,f5}
/// (12 scalar arrays) from cell centers to face half-nodes in direction @p dir.
/// Results stored in lb.field temporary face arrays (vis_x_f2_face, ...).
///
/// @param dir  0=ξ, 1=η, 2=ζ
static void interp_cart_flux_to_faces(LocalBlock& lb, int dir, const Config& cfg);

/// 5d-step2+3 — Exchange face-interpolated Cartesian fluxes and assemble vis_xi/eta/zeta.
///
/// Prerequisites: interp_cart_flux_to_faces() must have been called on ALL blocks
/// for the same @p dir (the temporary face arrays in Field must be populated).
///
/// 1. MPI exchange via FluxHaloExchange::exchange_face_arrays (12 arrays)
/// 2. Same-process local copy via copy_local_face_array (12 arrays)
/// 3. Assemble vis_xi/eta/zeta: f1=0, f2..f5 = Σ vis_*_face_comp × face_metric_comp
///
/// The 12 temporary face arrays are deallocated after assembly.
///
/// @param dir  0=ξ, 1=η, 2=ζ
static void exchange_and_assemble_face_flux(LocalBlock& lb, int dir,
                                             const std::vector<LocalBlock>& all_blocks,
                                             FluxHaloExchange& flux_ex,
                                             const Config& cfg);
```

**5e. 显式差分 — 粘性 RHS 计算**：

**目标**：将 face 粘性通量 `vis_xi/eta/zeta` 差分回 cell-center，累加入 `lb.field.rhs`。

**方案**：使用与 `InviscidRHS::compute` 相同的 inline 6 阶中心差分系数，直接对 `vis_xi/eta/zeta` 的 face-flux 进行差分。不使用 `InterpDiff::derivative_from_faces`（该方法对全 grid cell 差分，而我们仅需内点；且需额外分配临时数组）。

**差分公式**（dh=1.0，计算空间）：

```
∂F/∂ξ(i) = c0·(F[i+1] - F[i]) + c1·(F[i+2] - F[i-1]) + c2·(F[i+3] - F[i-2])
c0 = 75/64,  c1 = -25/384,  c2 = 3/640
```

其中 `F[h] = vis_xi(h,j,k)` 是 face 数组（`h` 为 face index，`F[h] = F_{h-1/2}`）。

**计算范围**：仅内点 `(3..nci-4)×(3..ncj-4)×(3..nck-4)`，与 `InviscidRHS` 完全一致。

**累加方式**（注意符号——粘性项 RHS 与无粘项符号相同）：

```
RHS = RHS - (∂F_vis/∂ξ + ∂G_vis/∂η + ∂H_vis/∂ζ)
```

即每个方向的 flux 导数都**减去**。若 `InviscidRHS::compute` 已在 `rhs` 中写入了 `-∂F_inv/∂ξ` 等，则 5e 直接在此基础上再减去粘性项导数。

**ViscidRHS 新增方法**：

```cpp
/// 5e — Compute viscous RHS contribution from face fluxes.
///
/// Prerequisites: assemble_face_fluxes() must have been called (fills vis_xi/eta/zeta).
///
/// Differentiates vis_xi, vis_eta, vis_zeta at interior cells using 6th-order
/// centered differences.  Subtracts each flux derivative from lb.field.rhs
/// (accumulating with the inviscid RHS).
///
/// Only interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] are updated.
static void compute_rhs(LocalBlock& lb);
```

**实现模式**：完全复用 `InviscidRHS::compute` 的结构，将 `inv_xi/eta/zeta` 替换为 `vis_xi/eta/zeta`。f1 分量（密度通量）对粘性项恒为 0，其导数和减去操作仍执行以保证代码一致性。

**涉及文件**（子任务 5c/5d/5e 汇总更新）：

| 文件 | 修改内容 |
|------|---------|
| `include/wcns_v2/field/field.h` | ✅ 5c `vis_x/y/z` (FluxVars)；✅ 5d 12 临时 face 数组 |
| `include/wcns_v2/field/field.hxx` | ✅ 5c `vis_x/y/z.allocate()` |
| `include/wcns_v2/scheme/viscid_rhs.h` | ✅ 5c/5d/5e 全部方法声明 |
| `include/wcns_v2/scheme/viscid_rhs.hxx` | ✅ 5c/5d/5e 全部实现 |
| `include/wcns_v2/parallel/parallel_manager.h` | ✅ `exchange_viscous_flux_halos` |
| `include/wcns_v2/parallel/parallel_manager.hxx` | ✅ `exchange_viscous_flux_halos` 实现 |
| `src/io/config_reader.cpp` | ✅ 解析 viscous 参数 |
| `src/main.cpp` | ✅ 5a/5b/5c/5d/5e 测试段 |
| `input.ini` | ✅ `viscous_type = constant` |
| `main_cycle_plan.md` | ✅ 5a-5e 设计文档 |

**依赖**：子任务 5c（vis_x/y/z + halo 交换）、FluxHaloExchange（face 级别 MPI 交换）、Grid（face metric）

### 无粘模式行为

当 `viscous_type = "none"` 时：
- 5a/5b 仍需执行（梯度计算可能有其他用途，如湍流模型）
- 5c: `compute_cell_viscous_flux` 将 `vis_x/vis_y/vis_z` 所有分量设为 0
- 5d-5e: 后续粘性通量组装和差分结果为 0，`RHS_vis ≡ 0`
- 最终 `RHS = RHS_inv + RHS_vis = RHS_inv`

---

### 子任务 5.5: 体积力源项

**目标**：实现体积力（body force）源项计算，将物理源项以 `J⁻¹·S` 形式累加入 RHS，为后续重力、离心力、Coriolis 力等体积力提供统一接口。

**背景 — SCMM 方程中的源项**：

带体积力源项的 NS 方程（物理空间）：
```
∂Q/∂t + ∂F/∂x + ∂G/∂y + ∂H/∂z = S_physical
```
其中物理源项 `S_physical = [0, ρ·f_x, ρ·f_y, ρ·f_z, ρ·(f_x·u + f_y·v + f_z·w)]^T`，`f = (f_x, f_y, f_z)` 为体积力加速度矢量（单位质量受力）。

SCMM 变换后：
```
J⁻¹·∂Q/∂t + ∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ = J⁻¹·S_physical
```

子任务 4/5 计算的 flux RHS 为：
```
rhs_flux = -(∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ) = J⁻¹·∂Q/∂t - J⁻¹·S_physical
```

为使完整 RHS 满足 `rhs = J⁻¹·∂Q/∂t`，需将源项补入：
```
rhs = rhs_flux + J⁻¹·S_physical
```

时间推进器中 `Q_new = α·Q0 + β·(Q_old + dt·J·rhs)` 将自动得到正确的物理更新：
```
dt·J·rhs = dt·J·rhs_flux + dt·J·(J⁻¹·S_physical) = dt·J·rhs_flux + dt·S_physical
```

**物理源项公式**（cell-center，5 分量）：

```
S_rho  = 0
S_rhou = ρ · f_x
S_rhov = ρ · f_y
S_rhow = ρ · f_z
S_rhoE = ρ · (f_x·u + f_y·v + f_z·w)
```

其中 ρ 为密度，`(u,v,w)` 为速度分量，`(f_x,f_y,f_z)` 为体积力加速度。

**Config 新增参数**（`[physics]` section）：

```cpp
// 体积力类型
// "none"     — 无体积力（默认）
// "constant" — 固定体积力加速度矢量
std::string body_force_type = "none";

// 体积力加速度分量（仅 body_force_type="constant" 时使用）
Real body_force_x = 0.0;
Real body_force_y = 0.0;
Real body_force_z = 0.0;
```

**类设计**：

```cpp
/// @file body_force.h
/// @brief Volumetric body force source term for NS equations.

class BodyForce {
public:
    /// Compute body force source term and add to RHS.
    ///
    /// Physical source: S = [0, ρ·f_x, ρ·f_y, ρ·f_z, ρ·f·v]^T
    /// Adds J⁻¹·S to lb.field.rhs so that the time integrator's
    /// "multiply by J" step correctly recovers S_physical:
    ///   Q_new = α·Q0 + β·(Q_old + dt·J·rhs)
    ///         = α·Q0 + β·(Q_old + dt·J·rhs_flux + dt·S_physical)
    ///
    /// @param lb  Local block (reads prim.{rho,u,v,w}, grid.jacobian;
    ///             adds to field.rhs)
    /// @param cfg Configuration (body_force_type, body_force_{x,y,z})
    ///
    /// Only interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] are updated.
    /// If body_force_type == "none", returns immediately (no-op).
    static void add_to_rhs(LocalBlock& lb, const Config& cfg);
};
```

**实现核心**（仅 interior cells）：

```cpp
if (cfg.body_force_type == "none") return;

Real fx = cfg.body_force_x;
Real fy = cfg.body_force_y;
Real fz = cfg.body_force_z;

for (Int k = k0; k <= k1; ++k) {
for (Int j = j0; j <= j1; ++j) {
for (Int i = i0; i <= i1; ++i) {
    Real rho = f.prim.rho(i,j,k);
    Real u   = f.prim.u(i,j,k);
    Real v   = f.prim.v(i,j,k);
    Real w   = f.prim.w(i,j,k);
    Real invJ = 1.0 / J(i,j,k);

    // S_physical → J⁻¹·S_physical → rhs
    f.rhs.rho(i,j,k)  += 0.0;                               // = 0
    f.rhs.rhou(i,j,k) += invJ * rho * fx;
    f.rhs.rhov(i,j,k) += invJ * rho * fy;
    f.rhs.rhow(i,j,k) += invJ * rho * fz;
    f.rhs.rhoE(i,j,k) += invJ * rho * (fx*u + fy*v + fz*w);
}}}
```

**调用时机**（主循环中，见子任务 9 伪代码）：

```
2e. compute_rhs (粘性 RHS 差分)  ← rhs 此时含 inv+vis 两部分
2f. BodyForce::add_to_rhs         ← rhs 再加入体积力源项 = inv+vis+body
3.  TimeIntegrator::advance_stage ← 时间推进
```

即在 inviscid RHS 和 viscid RHS 都已累加入 `rhs` 之后、时间推进之前调用 `BodyForce::add_to_rhs`。

**涉及文件**：

| 文件 | 操作 | 说明 |
|------|------|------|
| `include/wcns_v2/scheme/body_force.h` | **新建** | `BodyForce` 类声明 |
| `include/wcns_v2/scheme/body_force.hxx` | **新建** | `BodyForce::add_to_rhs` 实现 |
| `include/wcns_v2/core/config.h` | 修改 | 新增 `body_force_type`, `body_force_{x,y,z}` |
| `src/io/config_reader.cpp` | 修改 | 新增 `[physics]` section 解析 |
| `input.ini` | 修改 | 新增 `[physics]` section |
| `src/main.cpp` | 修改 | 新增 5.5 测试段 |

**依赖**：子任务 4 (InviscidRHS)、子任务 5 (ViscidRHS) — 仅依赖 RHS 数组已分配、prim 变量有效、Jacobian 可访问。与子任务 4/5 无数据流耦合，可与子任务 4/5 并行开发。

**测试验证**：

1. **none 模式**：`body_force_type = "none"` → RHS 不变（与未调用一致）
2. **零力矢量**：`f = (0,0,0)` → RHS 不变
3. **均匀流 + 恒定 f_x**：仅 `rhs.rhou` 和 `rhs.rhoE` 有非零增量，验证 `J⁻¹·ρ·f_x` 数值正确
4. **无 NaN/Inf**：在单区/双区网格上运行，确保源项计算无异常值

---

### 子任务 6: 时间推进器

**目标**：实现显式 Runge-Kutta 子步循环，时间推进格式在 input 中切换，当前实现 RK3-TVD。

**核心设计 — Jacobian 修正**：

SCMM 形式的 Euler/NS 方程：
```
J⁻¹ · ∂Q/∂t + ∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ = 0
```
其中 `F̂ = ξ̂_x·F_phys + ξ̂_y·G_phys + ξ̂_z·H_phys`（即 `face_metric · physical_flux`）。

各 RHS Compute（InviscidRHS::compute / ViscidRHS::compute_rhs）差分 face flux 得到的 rhs 是：
```
rhs = -(∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ) = J⁻¹ · ∂Q/∂t
```
在时间推进时，需要转换为 `∂Q/∂t`：
```
∂Q/∂t = J · rhs
Q_new = Q_old + Δt · J · rhs   （乘以 Jacobian）
```

**RK3-TVD 格式**（当前实现）：
| m | α_m | β_m |
|---|-----|-----|
| 1 | 0   | 1   |
| 2 | 3/4 | 1/4 |
| 3 | 1/3 | 2/3 |

更新公式：`Q^(m) = α_m * Q^(0) + β_m * (Q^(m-1) + Δt · J · RHS^(m-1))`

**RK4 格式**（经典四阶，预留）：
| m | α_m | β_m |
|---|-----|-----|
| 1 | 1.0 | 1/2 |
| 2 | 1.0 | 1/2 |
| 3 | 1.0 | 1.0 |
| 4 | 0.0 | 1/6 |  (需要前三阶段加权)

**设计要点**：
- 类 `TimeIntegrator`，纯静态方法（无状态）
- `n_stages(cfg)` — 返回 RK 子步数
- `rk_coeffs(cfg, stage, alpha, beta)` — 返回给定 stage 的 α, β 系数
- `advance_stage(lb, cfg, dt, stage, Q0)` — 执行单个 RK 子步：
  1. 读取 `lb.field.rhs`（需预先计算 inv+vis RHS）
  2. 读取 `lb.field.cons`（当前 Q^(m-1)）
  3. 读取 `Q0`（初始 Q^(0)）
  4. 计算 `Q_new(i,j,k) = α*Q0(i,j,k) + β*(Q_old(i,j,k) + dt*J(i,j,k)*rhs(i,j,k))`
  5. 仅更新 interior cells `[3..nci-4]×[3..ncj-4]×[3..nck-4]`
- `Q0` 存储为 `Field::Q0` (`ConservativeVars`)，在每时间步开始时保存
- 预留 `lu_sgs_step()` 空接口（隐式方法）

**涉及文件**：
- 修改 `include/wcns_v2/field/field.h` — 新增 `ConservativeVars Q0`
- 修改 `include/wcns_v2/field/field.hxx` — 新增 `Q0.allocate()`
- 新文件 `include/wcns_v2/time/time_integrator.h` — 类声明
- 新文件 `include/wcns_v2/time/time_integrator.hxx` — 实现
- 修改 `src/main.cpp` — 新增 6a/6b 测试段

**依赖**：子任务 4 (InviscidRHS)、子任务 5 (ViscidRHS)、子任务 1 (TimeStep)、Grid (Jacobian)

---

### 子任务 7: 残差计算与物理量追踪

**目标**：计算全部 5 个守恒量残差的 L2 范数，统计全场物理量极值（温度、压力、密度），定期输出到 `residual.dat`，并判断收敛。

#### 7.1 残差计算

对每个 block 的内点 `[3..nci-4]×[3..ncj-4]×[3..nck-4]`（与 RHS 计算范围一致），累加各守恒量分量的平方和：

```
local_sq[0] = Σ rhs.rho(i,j,k)²
local_sq[1] = Σ rhs.rhou(i,j,k)²
local_sq[2] = Σ rhs.rhov(i,j,k)²
local_sq[3] = Σ rhs.rhow(i,j,k)²
local_sq[4] = Σ rhs.rhoE(i,j,k)²
local_N    = 内点总数 (nci_core * ncj_core * nck_core)
```

跨 block + MPI 全局归约（`ParallelManager::global_sum`）：

```
res[c] = sqrt(global_sum(local_sq[c]) / global_sum(local_N))
```

**收敛判据**：以密度残差 `res_rho` 为基准，`res_rho < converge_tol` 时视为收敛。

#### 7.2 物理量追踪

同一次内点遍历中统计全场极值（避免重复循环）：

| 物理量 | 符号 | 无量纲公式 | 说明 |
|--------|------|-----------|------|
| 密度 | ρ | `prim.rho` | 直接读取 |
| 温度 | T | `prim.p * γ * Mach² / prim.rho` | 由无量纲 EOS 反推 |
| 压力 | p | `prim.p` | 直接读取 |
| 速度幅值 | \|U\| | `sqrt(u²+v²+w²)` | |
| 当地 Mach 数 | M_local | `|U| / sqrt(γ*p/ρ)` | 可选，检测超声速区 |

每个 block 计算 local min/max，跨 block 用 `ParallelManager::global_min` / `global_max` 归约。

**发散检测**：若任一监测量出现 NaN 或 Inf，立即终止并报错。

#### 7.3 数据结构

```cpp
/// 5 分量守恒残差 L2 范数
struct ResidualNorms {
    Real rho;   ///< ||RHS_rho||_L2
    Real rhou;  ///< ||RHS_rhou||_L2
    Real rhov;  ///< ||RHS_rhov||_L2
    Real rhow;  ///< ||RHS_rhow||_L2
    Real rhoE;  ///< ||RHS_rhoE||_L2
};

/// 全场物理量极值（可扩展）
struct FlowMonitor {
    Real rho_min, rho_max;   ///< 密度范围
    Real T_min, T_max;       ///< 温度范围
    Real p_min, p_max;       ///< 压力范围
    Real vel_min, vel_max;   ///< 速度幅值范围
    Real Mach_max;           ///< 最大当地 Mach 数

    /// 仅 rank 0 有效，全部为 false 表示正常
    bool has_nan = false;
    bool has_inf = false;

    /// 预留扩展：后续可增加湍流 k/ω、组分浓度等
    // std::map<std::string, Real> extras;
};
```

#### 7.4 类接口

```cpp
class Residual {
public:
    /// 计算全部 5 守恒量残差 L2 范数（跨 block 全局归约，所有 rank 结果一致）
    static ResidualNorms compute(const std::vector<LocalBlock>& blocks);

    /// 统计全场物理量极值（跨 block 全局归约，所有 rank 结果一致）
    static FlowMonitor monitor(const std::vector<LocalBlock>& blocks, const Config& cfg);

    /// 写入 residual.dat 表头（仅在文件开头调用一次）
    static void write_header(std::ostream& os);

    /// 写入一行残差 + 监测量记录
    static void log(std::ostream& os, Int iter, Real dt,
                    const ResidualNorms& res, const FlowMonitor& mon);

    /// 收敛判断：密度残差 < 收敛容差
    static bool converged(const ResidualNorms& res, Real tol);

    /// 发散检测：任一监测量出现 NaN/Inf
    static bool diverged(const FlowMonitor& mon);
};
```

**设计要点**：
- `compute()` 和 `monitor()` 分开调用（但可共享同一次循环的临时结果以优化性能）
- 两者均通过 `ParallelManager::global_sum/min/max` 归约，所有 rank 结果一致
- 仅 rank 0 写文件，避免竞争

#### 7.5 残差文件格式 (`residual.dat`)

表头 + 数据行（空格分隔，科学计数法）：

```
# iter       dt             res_rho        res_rhou       res_rhov       res_rhow       res_rhoE        rho_min        rho_max        T_min          T_max          p_min          p_max          vel_max        Mach_max
     1  1.234567e-03  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00  1.000000e+00
```

#### 7.6 Config 新增字段

```cpp
Int residual_freq = 1;   // 残差/监测计算频率（默认每步，开销极小：仅一次内点遍历 + 归约）
```

若 `residual_freq = 1`，每时间步的最终 RK stage 后计算并记录。

#### 7.7 涉及文件

| 文件 | 操作 | 说明 |
|------|------|------|
| `include/wcns_v2/core/residual.h` | **新建** | `ResidualNorms`、`FlowMonitor` 结构体 + `Residual` 类声明 |
| `include/wcns_v2/core/residual.hxx` | **新建** | `compute()`、`monitor()`、`log()`、`converged()`、`diverged()` 实现 |
| `include/wcns_v2/core/config.h` | 修改 | 新增 `residual_freq` 字段（`[control]` section） |
| `src/io/config_reader.cpp` | 修改 | 新增 `residual_freq` 解析 |
| `input.ini` | 修改 | 新增 `residual_freq = 1` |

#### 7.8 依赖

- Field（prim、rhs）：内点遍历读取
- Config（gamma、Mach、converge_tol、residual_freq）：参数
- ParallelManager（global_sum、global_min、global_max）：MPI 归约
- 无 Block 间顺序依赖，可独立于子任务 8 开发

---

### 子任务 8: 输出模块（Solution + Restart）

**目标**：每间隔 N 步或达到指定物理时刻时，输出流场文件用于后处理和续算。Solution 输出仅含 primitive 变量（`rho, u, v, w, p`）和 cell-center 坐标；Restart 输出完整的 conservative 变量（含 ghost cell），读取 restart 时自行由 cons vars 计算 prim vars。

**输出格式**：通过 `output_format` 参数选择，当前实现 Tecplot（ASCII `.plt`），CGNS 和 VTK 预留接口。

---

#### 8.1 Config 新增参数

新增 `[output]` section：

```cpp
// =========================================================================
// Output control
// =========================================================================

/// Solution output format: "tecplot" (implemented) | "cgns" (reserved) | "vtk" (reserved)
std::string output_format = "tecplot";

/// Output directory ("" = current working directory)
std::string output_dir = "";

/// Physical time interval for solution output (0 = disabled, use step-based only)
/// Reserved — requires time accumulation in main loop.
Real output_time_interval = 0.0;
```

已有字段 `output_freq` / `restart_freq` 在 `[control]` section 中，无需移动。

---

#### 8.2 内点输出范围

Solution 输出**内点** `[ng .. nci-1-ng]`，即 cell index `[3 .. nci-4]`：

- ghost cell `[0..2]` 和 `[nci-3..nci-1]` 不输出（由 BC/Halo 填充，非物理）
- 输出维度 = `nci_core × ncj_core × nck_core`，其中 `nci_core = nci - 2*ng`

Restart 输出**全范围** cell `[0 .. nci-1]`（含 ghost），保证续算时状态精确还原。

---

#### 8.3 类设计

```
                    ┌─────────────────────────┐
                    │     SolutionWriter        │  ← 统一入口，根据 output_format 分发
                    │    (静态方法)              │
                    └────────────┬──────────────┘
                                 │
                ┌────────────────┼────────────────┐
                ▼                ▼                ▼
       ┌────────────┐   ┌────────────┐   ┌────────────┐
       │TecplotWriter│   │ CgnsWriter │   │ VtkWriter  │
       │  (已实现)    │   │  (预留)    │   │  (预留)    │
       └────────────┘   └────────────┘   └────────────┘

       ┌─────────────────────────┐
       │     RestartWriter        │  ← 二进制 restart I/O
       │   write() + read()      │
       └─────────────────────────┘
```

**SolutionWriter**：

```cpp
class SolutionWriter {
public:
    /// Write solution at current time step.
    /// Dispatches to format-specific writer based on cfg.output_format.
    /// Only rank 0 writes the file; all ranks participate in data gathering.
    static void write(const std::vector<LocalBlock>& blocks,
                      const Config& cfg, Int iter, Real time);
};
```

**TecplotWriter**（内部实现，用户不直接调用）：

```cpp
class TecplotWriter {
public:
    /// Write multi-zone Tecplot ASCII file.
    ///
    /// Gathers interior primitive vars and cell-center coordinates from all
    /// blocks to rank 0 via MPI_Gather, then writes a single .plt file with
    /// one ZONE per block.
    ///
    /// @param filename  Full path to output file (e.g. "output/sol_000100.plt")
    /// @param iter      Current iteration number (for title annotation)
    /// @param time      Current physical time (for title annotation)
    static void write(const std::vector<LocalBlock>& blocks,
                      const std::string& filename, Int iter, Real time);
};
```

**RestartWriter**：

```cpp
class RestartWriter {
public:
    /// Write restart files (binary).
    ///
    /// Each rank writes its own block data to restart_XXXXXX_rN.bin.
    /// Rank 0 additionally writes an index file restart_XXXXXX.bin mapping
    /// global_block_id → rank so the decomposition can be reconstructed on
    /// a different number of processes.
    ///
    /// Saved per block: nci, ncj, nck, ng, block_id, zone_id,
    ///   cons.rho[], cons.rhou[], cons.rhov[], cons.rhow[], cons.rhoE[]
    ///   (full arrays including ghost cells, flat binary dump in (k,j,i) order)
    static void write(const std::vector<LocalBlock>& blocks,
                      const Config& cfg, Int iter, Real time);

    /// Read restart files and populate fields.
    ///
    /// Each rank reads its own restart_XXXXXX_rN.bin, loads conservative vars,
    /// then computes primitive vars via Field::cons_to_prim(gamma).
    /// The index file restart_XXXXXX.bin is used to verify block distribution.
    ///
    /// @param basename  Restart base name (e.g. "restart_000500")
    /// @param blocks    [out] Local blocks with fields to populate
    /// @param cfg       Configuration (provides gamma for cons→prim conversion)
    /// @return          {iteration, time} from the restart point
    static std::pair<Int, Real> read(const std::string& basename,
                                      std::vector<LocalBlock>& blocks,
                                      const Config& cfg);
};
```

---

#### 8.4 Tecplot 文件格式 (`.plt`)

多 block 结构化网格，每个 LocalBlock 为一个 ZONE：

```
TITLE = "WCNS v2 Solution — iter=000100, time=1.234567e-01"
VARIABLES = "X", "Y", "Z", "rho", "u", "v", "w", "p"
ZONE T="block_0", I=64, J=32, K=16, DATAPACKING=POINT
 1.0000e+00  0.0000e+00  0.0000e+00  1.2250e+00  1.0000e+00  0.0000e+00  0.0000e+00  1.0000e+00
 1.0156e+00  0.0000e+00  0.0000e+00  1.2250e+00  1.0000e+00  0.0000e+00  0.0000e+00  1.0000e+00
 ...
ZONE T="block_1", I=64, J=32, K=16, DATAPACKING=POINT
 ...
```

**要点**：
- `DATAPACKING=POINT`：每行一个点的全部变量（与代码逐点写出顺序一致）
- 坐标输出 **cell-center 坐标**（`grid.cell_x/y/z`），非 node 坐标
- ZONE name 用 `"block_" + block_id` 区分
- 内点遍历顺序：`k → j → i`（最内层 i，与内存布局一致）

---

#### 8.5 Restart 二进制格式 (`.bin`)

**索引文件**（rank 0 写入，记录 block→rank 映射）：

| 字段 | 类型 | 字节 |
|------|------|------|
| magic | `uint32_t` | 4 (`0x57434E53` = "WCNS") |
| version | `uint32_t` | 4 (`1`) |
| iter | `Int` | 4 |
| time | `Real` | 8 |
| nprocs | `Int` | 4 |
| n_blocks_total | `Int` | 4 |
| for each block: block_id (Int), rank (Int) | — | `n_blocks_total × 8` |

**每 rank 数据文件** (`restart_XXXXXX_rN.bin`)：

| 字段 | 类型 | 字节 |
|------|------|------|
| n_blocks_local | `Int` | 4 |
| **for each local block:** | | |
| — block_id | `Int` | 4 |
| — zone_id | `Int` | 4 |
| — nci | `Int` | 4 |
| — ncj | `Int` | 4 |
| — nck | `Int` | 4 |
| — ng | `Int` | 4 |
| — cons.rho[] | `Real × nci×ncj×nck` | … |
| — cons.rhou[] | `Real × nci×ncj×nck` | … |
| — cons.rhov[] | `Real × nci×ncj×nck` | … |
| — cons.rhow[] | `Real × nci×ncj×nck` | … |
| — cons.rhoE[] | `Real × nci×ncj×nck` | … |

每个 `MultiArray3D<Real>` 按 `(k,j,i)` 顺序写入（与内存 `data()` 布局一致），一次 `ostream::write` 调用：
```cpp
os.write(reinterpret_cast<const char*>(arr.data()), arr.size() * sizeof(Real));
```

**读取时**：按相同顺序 `read` 回 `MultiArray3D::data()`，然后调用 `Field::cons_to_prim(gamma)` 恢复 prim 变量。

---

#### 8.6 MPI 数据收集策略

**Solution 输出 — Gather-to-rank-0 方案**：

1. 每个 rank 将其本地 blocks 的内点 primitive vars + 坐标打包为 buffer
2. `MPI_Gather` 各 rank 的 buffer size 到 rank 0
3. `MPI_Gatherv` 收集实际数据到 rank 0
4. Rank 0 按 block_id 顺序写入 Tecplot 文件（单进程 I/O，无竞争）

**优点**：输出为单文件，后处理直接加载；实现简单
**内存代价**（rank 0 峰值，128³ 网格 × 5 blocks × 8 变量）≈ 400 MB，可接受

**Restart 输出 — 每 rank 独立文件**：

1. 每个 rank 将其本地 blocks 的 cons vars 写入 `restart_XXXXXX_rN.bin`
2. Rank 0 额外写索引文件 `restart_XXXXXX.bin`
3. 无需 MPI 通信

**优点**：无 gather 开销；可处理大规模网格；支持进程数变化的续算

---

#### 8.7 文件名生成

```cpp
/// 生成带 6 位补零迭代号的文件名
/// make_filename("output", "sol", 100, ".plt") → "output/sol_000100.plt"
inline std::string make_filename(const std::string& dir,
                                  const std::string& prefix,
                                  Int iter,
                                  const std::string& suffix) {
    std::ostringstream oss;
    if (!dir.empty()) oss << dir << "/";
    oss << prefix << "_" << std::setw(6) << std::setfill('0') << iter << suffix;
    return oss.str();
}
```

调用示例：
```cpp
auto filename = make_filename(cfg.output_dir, "sol", iter, ".plt");
TecplotWriter::write(blocks, filename, iter, time);
```

---

#### 8.8 CGNS / VTK 预留接口

```cpp
/// CGNS 输出（预留，暂未实现）
class CgnsWriter {
public:
    static void write(const std::vector<LocalBlock>& blocks,
                      const std::string& filename, Int iter, Real time);
    // → throw std::runtime_error("CGNS output not yet implemented")
};

/// VTK 输出（预留，暂未实现）
class VtkWriter {
public:
    static void write(const std::vector<LocalBlock>& blocks,
                      const std::string& filename, Int iter, Real time);
    // → throw std::runtime_error("VTK output not yet implemented")
};
```

`SolutionWriter::write` 分发逻辑：

```cpp
if (cfg.output_format == "tecplot") {
    TecplotWriter::write(blocks, filename, iter, time);
} else if (cfg.output_format == "cgns") {
    CgnsWriter::write(blocks, filename, iter, time);
} else if (cfg.output_format == "vtk") {
    VtkWriter::write(blocks, filename, iter, time);
} else {
    throw std::runtime_error("Unknown output_format: " + cfg.output_format);
}
```

---

#### 8.9 涉及文件

| 文件 | 操作 | 说明 |
|------|------|------|
| `include/wcns_v2/io/solution_writer.h` | **新建** | `SolutionWriter`、`TecplotWriter` 类声明 + `CgnsWriter`/`VtkWriter` 预留声明 |
| `include/wcns_v2/io/solution_writer.hxx` | **新建** | Tecplot 写文件、MPI gather、文件名生成、分发逻辑 |
| `include/wcns_v2/io/restart_writer.h` | **新建** | `RestartWriter` 类声明 |
| `include/wcns_v2/io/restart_writer.hxx` | **新建** | 二进制 write/read 实现 |
| `include/wcns_v2/core/config.h` | 修改 | 新增 `output_format`、`output_dir`、`output_time_interval` 字段 |
| `src/io/config_reader.cpp` | 修改 | 新增 `[output]` section 解析 |
| `input.ini` | 修改 | 新增 `[output]` section |
| `src/main.cpp` | 修改 | 调用 `SolutionWriter::write` 和 `RestartWriter::write` |

---

#### 8.10 依赖

- Field（prim、cons 访问）
- Grid（cell_x/y/z 坐标、`nci/ncj/nck` 维度信息、`ng`）
- Config（`output_format`、`output_freq`、`output_dir`）
- ParallelEnv（`rank`、`is_master`、`MPI_Gather`/`MPI_Gatherv`）
- 无 Block 间数据依赖，可与子任务 1、7 并行开发

---

### 子任务 9: 主循环集成

**目标**：在 `main.cpp` 中串联所有模块，形成完整的时间推进循环。

**伪代码**（展示每个 RK 子步的完整 RHS 计算流程）：
```cpp
// 初始化（已完成）
FlowInitializer::initialize(lb, cfg);
BoundaryConditionApplier::apply_all(lb, cfg);
pm.exchange_all_halos(blocks);

// 保存 Q^(0) 副本
std::vector<ConservativeVars> Q0(blocks.size());

// 打开残差日志文件 (rank 0 only)
std::ofstream res_file;
if (ParallelEnv::is_master()) {
    res_file.open("residual.dat");
    Residual::write_header(res_file);
}

bool converged_flag = false;

// 主循环
for (Int iter = 1; iter <= cfg.max_iter; ++iter) {
    Real dt = TimeStep::compute(blocks, cfg);

    for (Int stage = 0; stage < n_stages; ++stage) {
        // 保存初始守恒量（stage 0）
        if (stage == 0) save_Q0();

        // ============================================================
        // 1. 无粘 RHS: 插值 → Riemann → 通量交换 → 差分
        // ============================================================
        // 1a. WCNS 插值 (conservative vars → face Q_L/Q_R)
        for (auto& lb : blocks)
            WcnsInterp::interp_all(lb, cfg);

        // 1b. Riemann 求解 (Q_L/Q_R → inv_xi/eta/zeta face fluxes)
        for (auto& lb : blocks)
            RiemannSolver::solve_all(lb, cfg);

        // 1c. Face 通量 halo 交换 (ng+1 faces at connectivity boundaries)
        pm.exchange_flux_halos(blocks);

        // 1d. 差分 → RHS_inv (interior cells, accumulated in rhs)
        for (auto& lb : blocks)
            InviscidRHS::compute(lb);

        // ============================================================
        // 2. 粘性 RHS: 插值 → 梯度 → 粘性通量 → 组装 → 差分
        // ============================================================
        // 2a. 物理量插值到 face (u,v,w,T cell-center → faces, 3 dirs)
        for (auto& lb : blocks)
            ViscidRHS::interp_to_faces(lb, cfg);

        // 2b. 梯度计算 (chain rule, 含 face-product exchange)
        for (auto& lb : blocks)
            ViscidRHS::compute_gradients(lb, blocks, pm.flux_halo_ex(idx), cfg);

        // 2b-exchange: 梯度 ghost cell 交换
        pm.exchange_gradient_halos(blocks);

        // 2c. Cell-center Cartesian 粘性通量 (τ, q 中间量 → vis_x/y/z)
        for (auto& lb : blocks)
            ViscidRHS::compute_cell_viscous_flux(lb, cfg);

        // 2c-exchange: vis_x/y/z ghost cell 交换
        pm.exchange_viscous_flux_halos(blocks);

        // 2d. 半节点通量组装 (方案C: 插值到face → 交换 → metric组装)
        for (int dir = 0; dir < 3; ++dir) {
            // step1: 所有 block 插值 Cartesian 通量到 faces
            for (auto& lb : blocks)
                ViscidRHS::interp_cart_flux_to_faces(lb, dir, cfg);
            // step2+3: 交换 + Metric 组装 → vis_xi/eta/zeta
            for (auto& lb : blocks)
                ViscidRHS::exchange_and_assemble_face_flux(lb, dir, blocks,
                    pm.flux_halo_ex(lb_idx), cfg);
        }

        // 2e. 粘性 RHS 差分 (interior cells, accumulated in rhs)
        for (auto& lb : blocks)
            ViscidRHS::compute_rhs(lb);

        // 2f. 体积力源项 (累加入 rhs)
        for (auto& lb : blocks)
            BodyForce::add_to_rhs(lb, cfg);

        // ============================================================
        // 3. 时间更新
        // ============================================================
        for (auto& lb : blocks) {
            // Q^(m) = α*Q^(0) + β*(Q^(m-1) + dt*RHS)
            TimeIntegrator::advance(lb, cfg, dt, stage, Q0);

            // 转换为原始变量
            lb.field.cons_to_prim(cfg.gamma);

            // 施加 BC
            BoundaryConditionApplier::apply_all(lb, cfg);
        }
        // Primitive variable halo exchange
        pm.exchange_all_halos(blocks);

        // 最终 stage: 残差计算 + 物理量追踪（遍历一次内点，开销极小）
        if (stage == n_stages - 1) {
            auto res = Residual::compute(blocks);         // 5 分量 L2 范数
            auto mon = Residual::monitor(blocks, cfg);    // 全场 T/p/ρ 极值

            // 发散检测（NaN/Inf 即时终止）
            if (Residual::diverged(mon)) {
                if (ParallelEnv::is_master()) {
                    std::cerr << "FATAL: Divergence detected at iter="
                              << iter << "\n";
                }
                break;  // 或 throw
            }

            // 按频率写入 residual.dat
            if (iter % cfg.residual_freq == 0) {
                if (ParallelEnv::is_master())
                    Residual::log(res_file, iter, dt, res, mon);
            }

            // 收敛判断
            if (Residual::converged(res, cfg.converge_tol)) {
                if (ParallelEnv::is_master())
                    std::cout << "Converged at iter=" << iter << "\n";
                converged_flag = true;
            }
        }
    }
    if (converged_flag) break;

    // 定频输出
    if (iter % cfg.output_freq == 0) {
        SolutionWriter::write(blocks, cfg, iter);
    }
    if (iter % cfg.restart_freq == 0) {
        RestartWriter::write(blocks, cfg, iter);
    }
}
```

**涉及文件**：
- 修改 `src/main.cpp`
- 可能需要 `include/wcns_v2/core/solver_driver.h` / `.hxx`（可选，将主循环逻辑从 main 抽离）

**依赖**：所有子任务

---

## 依赖关系图

```
子任务 1 (时间步长)
    │
    ├──→ 子任务 4 (无粘通量 RHS)
    │        ├── 依赖 子任务 2 (WCNS 插值)
    │        ├── 依赖 子任务 3 (Riemann 求解器)
    │        └── 依赖 子任务 3.5 (界面通量通讯)
    │
    ├──→ 子任务 3.5 (界面通量通讯)
    │        └── 依赖 子任务 3 (Riemann 求解器)
    │
    ├──→ 子任务 5 (粘性通量 RHS)
    │        ├── 5a (物理量插值到 face) — 依赖 InterpDiff
    │        ├── 5b (梯度计算) — 依赖 5a、FluxHaloExchange、InterpDiff
    │        │     └── exchange_gradient_halos (ParallelManager)
    │        ├── 5c (cell-center Cartesian 粘性通量) — 依赖 5b、Config
    │        │     └── exchange_viscous_flux_halos (ParallelManager)
    │        ├── 5d (face 通量组装, 方案C) — 依赖 5c、FluxHaloExchange、Grid metrics
    │        │     ├── interp_cart_flux_to_faces (step1: 插值 → Field 临时 face 数组)
    │        │     └── exchange_and_assemble_face_flux (step2+3: 交换 + 组装)
    │        └── 5e (粘性 RHS 差分) — 依赖 5d (vis_xi/eta/zeta)
    │              └── compute_rhs (6阶中心差分, interior cells)
    │
    ├──→ 子任务 5.5 (体积力源项)
    │        └── 依赖 Field (prim, rhs, Jacobian) — 可并行于 4/5 开发
    │
    ├──→ 子任务 6 (时间推进器)
    │        ├── 依赖 子任务 4
    │        ├── 依赖 子任务 5
    │        ├── 依赖 子任务 5.5
    │        └── 依赖 子任务 1
    │
    ├──→ 子任务 7 (残差/收敛)
    │
    ├──→ 子任务 8 (输出)
    │
    └──→ 子任务 9 (主循环集成)
             └── 依赖 所有子任务
```

子任务 2、3、1、7、8 可并行开发。子任务 3.5 依赖 3，子任务 4 依赖 2+3+3.5，子任务 5 依赖 3.5（FluxHaloExchange），子任务 5.5 可并行于 4/5 开发（仅依赖 Field 已有分配），子任务 6 依赖 4+5+5.5+1，子任务 9 收尾。

## 实现顺序建议

| 阶段 | 子任务 | 状态 | 说明 |
|------|--------|------|------|
| 1 | 子任务 1 | ✅ 已完成 | 时间步长计算 |
| 2 | 子任务 3 | ✅ 已完成 | Roe Riemann 求解器 |
| 3 | 子任务 2 | ✅ 已完成 | WCNS 非线性插值 |
| 3.5 | 子任务 3.5 | ✅ 已完成 | 界面通量通讯 (FluxHaloExchange) |
| 4 | 子任务 4 | ✅ 已完成 | 无粘 RHS (InviscidRHS::compute) |
| 5a | 子任务 5a | ✅ 已完成 | 物理量 face 插值 (u,v,w,T) |
| 5b | 子任务 5b | ✅ 已完成 | 梯度计算 + halo 交换 |
| 5c | 子任务 5c | ✅ 已完成 | Cell-center Cartesian 粘性通量 |
| **5d** | **子任务 5d** | **✅ 已完成** | **Face 粘性通量组装 (方案C)** |
| **5e** | **子任务 5e** | **✅ 已完成** | **粘性 RHS 差分** |
| **5.5** | **子任务 5.5** | **✅ 已完成** | **体积力源项 (BodyForce)** |
| 6 | 子任务 6 | ✅ 已完成 | 时间推进器 (RK3-TVD) |
| 7 | 子任务 7 | ✅ 已完成 | 残差计算与物理量追踪 |
| 8 | 子任务 8 | ✅ 已完成 | 输出模块 (Solution + Restart) |
| 9 | 子任务 9 | ✅ 已完成 | 主循环集成 |

## 边界与 ghost cell 处理汇总

关键索引约定（以 x 方向为例，`nci = nx+6`）：

| 区域 | cell index | 半节点 index | 操作 |
|------|-----------|-------------|------|
| 左 ghost | 0..2 | — | BC/Halo 填充 |
| 左边界区 | 0..2 | 0..2 | InterpDiff 单边插值/差分 / WCNS 单边 (Q_L=Q_R) |
| 内点区 | 3..nci-4 | 3..nci-3 | WCNS 非线性插值 (Q_L≠Q_R) / 6阶中心差分 |
| 右边界区 | nci-3..nci-1 | nci-2..nci | InterpDiff 单边插值/差分 / WCNS 单边 (Q_L=Q_R) |
| 右 ghost | nci-3..nci-1 | — | BC/Halo 填充 |

**RHS 计算范围**：仅内点 `(3~nx+2)`（即 cell index `3..nci-4`）。边界点的 RHS 不需要精确计算（会被 BC 覆盖），但 ghost cells 必须在外插值中提供正确的支撑值。

**半节点数组大小**：
- x 方向 face：`(nci+1) × ncj × nck`，半节点下标 `0..nci`
- y 方向 face：`nci × (ncj+1) × nck`，半节点下标 `0..ncj`
- z 方向 face：`nci × ncj × (nck+1)`，半节点下标 `0..nck`

半节点约定：`ah[i]` 对应 `a_{i-1/2}`（half-node between cell i-1 and cell i）。
WCNS 插值输出 Q_L/Q_R 使用相同尺寸和约定。
