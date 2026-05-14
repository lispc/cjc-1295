# 活跃任务跟踪

## 任务 1：Rosetta FlexPepDock — WT GHRH(1-29) vs DPP-IV
**状态**：✅ ab-initio 完成（960 models），判定为方法学错误；🔄 refine-only 重跑中（10 models）
**阻塞项**：无
**下一步**：refine-only 完成后提取 score 与 RMSD
**备注**：Claude review 指出原协议误用 ab-initio 模式（`-lowres_preoptimize`）。起始结构来自 7CZ5 晶体对齐，应使用 refine-only（`-pep_refine` 无 `-lowres_preoptimize`）。当前 refine-only 正在运行。

---

## 任务 2：Rosetta FlexPepDock — D-Ala2 突变体 vs DPP-IV
**状态**：🔄 refine-only 运行中（10 models）
**阻塞项**：无
**下一步**：与 WT refine 结果对比 I_sc、BSA、interface HB
**备注**：PDB 中 `ALA` → `DAL`，Rosetta 自动应用 `D_AA` patch。已完成第 1 个 model（281s），观察到 `Using simple Rotamer generation logic for DALA`。

---

## 任务 3：GROMACS MD — WT 复合物（200 ns）
**状态**：🔄 生产 MD 运行中（预计 ~2.3 天）
**阻塞项**：无
**下一步**：完成后运行分析脚本（RMSF、SASA、几何时序、氢键）
**备注**：EM/NVT/NPT 均已完成。起始结构直接来自 7CZ5（通过 catalytic 4/4 验证）。

---

## 任务 4：GROMACS MD — D-Ala2 突变体（50 ns）
**状态**：🔄 NPT 运行中 → 随后启动生产 MD
**阻塞项**：无
**下一步**：50 ns 生产 MD → 与 WT 对比分析
**备注**：EM 已完成（1721 步，Epot = -5.219 MJ）。EM 后 φ=174°, ψ=177°（backbone 被拉直到 trans 构象），χ1 = -115.68°（手性翻转确认）。

---

## 任务 5：Step 4 数据分析与可视化
**状态**：⏳ 脚本已就绪，等 MD 完成
**阻塞项**：需 Task 3、4 完成
**下一步**：
  - RMSF、SASA 分析
  - 关键距离/角度时序（Ser630-Ala2、N-term-Glu205/206）
  - 氢键占有率统计
  - WT vs D-Ala2 对比图表
**预计完成时间**：MD 完成后 1–2 小时

---

## 任务 6：三肽几何对比可视化
**状态**：✅ 已完成
**产出**：
  - `workspace/figures/tripeptide_wt_vs_DAla_comparison.png`（3D 镜面 + Fischer 投影）
  - `workspace/figures/catalytic_geometry_WT_vs_DAla.png`（距离柱状图 + 攻击角仪表盘）

---

## 紧急修复（来自 Claude review）
- ✅ P0: `cc.sh` 已删除，`.gitignore` 已更新（API key 安全）
- ✅ P0: 根目录散文件已清理（`io.mc`、`posre.itp`、`#test_dala.top.1#`）
- ⏳ P0: refine-only FlexPepDock 重跑中（WT + D-Ala2）
- ⏳ P0: score-only 待执行（WT + D-Ala2 起始 pose）
- ⏳ P2: `run_flexpepdock_parallel.py` 的 `combine_silent` 命令待修复
- ⏳ P2: step2 无用 silent/log 待压缩清理

---

*最后更新：2026-05-13 21:30*
*维护者：Kimi Code CLI*
