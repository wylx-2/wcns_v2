# MPI 通讯修复 — 人工检查清单

> 生成日期: 2026-06-11 | 对应修复: Bug 1, Bug 2 (HaloExchange 同进程通讯)

---

## 修复摘要

| Bug | 描述 | 修复内容 | 影响文件 |
|-----|------|----------|----------|
| 1 | 同进程多 block 误走 MPI 自发送 | `is_remote` 排除同 rank | [halo_exchange.hxx](include/parallel/halo_exchange.hxx#L29-L31) |
| 2 | 同进程多 block 拷贝数据源错误 | 实现 `copy_local`，从邻居 block 读取 | [halo_exchange.h](include/parallel/halo_exchange.h), [halo_exchange.hxx](include/parallel/halo_exchange.hxx#L205-L222), [halo_exchange.hxx:277-312](include/parallel/halo_exchange.hxx#L277-L312) |
| 3 | FluxHaloExchange tag | 确认为非问题，无需修改 | — |
| API | `exchange_multi` 新参数 | 增加 `all_blocks` 参数 | [halo_exchange.h](include/parallel/halo_exchange.h), [parallel_manager.hxx](include/parallel/parallel_manager.hxx) |

---

## 一、编译验证

- [ ] **1.1** `cmake --build build` 零错误通过
- [ ] **1.2** 无新增编译警告（`exchange_multi` 中 `face_recv_idx`/`face_send_idx` 的 `-Wunused-but-set-variable` 是原有的）

---

## 二、单进程功能测试

### 2.1 单 block 周期边界

- [ ] **2.1.1** 单 block 周期算例（如 `Isentropic_curl.cgns` 1 进程），验证正常运行
- [ ] **2.1.2** 检查初始 ghost 与多步运行后结果与修复前一致（同 block 自周期路径未改变）

### 2.2 多 block（同进程，内部剖分）

- [ ] **2.2.1** 多 block 算例（如 `channel_turbulence` 1 进程多 block），验证正常运行
- [ ] **2.2.2** **关键验证**: 在 `exchange_multi` 的 Phase 2 `!fb.is_remote` 分支（同进程不同 block）中添加临时日志，确认：
  - 进入了 Case 2 分支（`else` 分支，即 `copy_local` 路径）
  - `neighbor` 不为 `nullptr`
  - `neighbor->block_id` 与预期的邻居 block ID 一致
- [ ] **2.2.3** 对比修复前后的残差历史，确认残差一致或改善（不应变差）

---

## 三、多进程功能测试

### 3.1 跨进程 halo 交换

- [ ] **3.1.1** 多 block 算例 2+ 进程运行，确认无 MPI 错误
- [ ] **3.1.2** 对比修复前（如果之前跑过）和修复后的残差历史，确认一致
- [ ] **3.1.3** 确认跨进程的 `is_remote == true` 路径未受影响

### 3.2 重启动

- [ ] **3.2.1** 用修复前产生的 restart 文件重启动，确认能正常读取并继续计算
- [ ] **3.2.2** 用不同进程数重启动（如 4 进程 → 2 进程），确认 halo exchange 重建正确

### 3.3 输出一致性

- [ ] **3.3.1** 同一算例，分别用 1 进程和 2 进程运行 10 步，对比第 10 步输出的 tecplot 文件
  - 两文件应逐字节一致（或差异在机器精度内）
  - 这是验证 MPI 分解不影响结果的关键测试

---

## 四、粘性管道验证

- [ ] **4.1** 粘性算例（如 `poiseuille_2` 或 `channel_turbulence`）2+ 进程运行，确认：
  - 梯度 exchange 正常
  - 粘性通量 exchange 正常
  - 面通量 exchange 正常（`FluxHaloExchange`，未修改但需确认未被间接影响）
- [ ] **4.2** 对比多进程与单进程的 U_avg 历史，确认一致

---

## 五、压力测试

- [ ] **5.1** 多 block 多进程（4+ 进程）运行 100+ 步，确认：
  - 无死锁
  - 无 segfault
  - 残差单调下降（不发散）
- [ ] **5.2** 不同进程数（1, 2, 4）运行同一算例相同步数，用 `diff` 比较最终输出的 `residual.dat`，确认差异在 1e-12 以内

---

## 六、代码审查

- [ ] **6.1** `exchange_multi` 新增的 `all_blocks` 参数在所有 3 个调用位置（`exchange_all_halos`、`exchange_gradient_halos`、`exchange_viscous_flux_halos`）都正确传入
- [ ] **6.2** `copy_local` 中 `pack_face` 使用 `ni.target_face` 和 `neighbor` block，确认：
  - 当 block A (IMAX, face=1) 连接 block B (IMIN, face=0) 时，`ni.target_face = 0`
  - `pack_face` 从 block B 的 `[ng, 2*ng-1]`（IMIN 侧内部）读取 ✓
  - `unpack_face` 写入 block A 的 `[nci-ng, nci-1]`（IMAX ghost） ✓
- [ ] **6.3** 同 block 自周期路径（`target_block == -1`）未被修改，行为与修复前一致
- [ ] **6.4** `start_exchange` / `wait_exchange` 签名已更新但未被调用，确认对现有功能无影响

---

## 七、性能验证

- [ ] **7.1** 多进程运行 50 步，对比修复前后的 wall-clock time，确认无性能退化
  - 预期：修复后同进程多 block 场景应略有加速（减少 MPI 自发送开销）
  - 跨进程场景性能应不变

---

## 八、已知未修复问题（非阻塞路径）

- [ ] **8.1** 确认: `start_exchange` / `wait_exchange` 仍不完整——`wait_exchange` 不执行 `unpack`，同进程多 block 也不支持。这些接口当前未被调用，不影响正常运行。若将来需要非阻塞模式，需进一步实现。

---

## 九、边界情况

- [ ] **9.1** 单 cell 厚度的 block（如 `ncj_core == 1`）多进程运行，确认 halo exchange 不越界
- [ ] **9.2** `exchange_multi` 传入空的 `all_blocks` 时，若邻居在同一进程，应安全降级（`neighbor == nullptr` 时跳过，不崩溃）
- [ ] **9.3** 同时存在同进程邻居和跨进程邻居的 block，确认两种路径都正确执行

---

*相关文档: [mpi.md](mpi.md), [ghost.md](ghost.md)*
