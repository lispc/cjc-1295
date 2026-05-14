# 活跃任务跟踪

## 任务 1：Rosetta FlexPepDock — WT GHRH(1-29) vs DPP-IV
**状态**：✅ ab-initio 完成（960 models），判定为方法学错误；✅ refine-only 已完成（6 models）
**阻塞项**：refine 未传 `-native`，RMS 列自指无意义；需要带约束重跑
**下一步**：带 `-native` + CoordinateConstraint 的 refine-only 重跑
**备注**：Claude review 指出原协议误用 ab-initio 模式。当前 refine-only 仍缺 `-native`，需修复。

---

## 任务 2：Rosetta FlexPepDock — D-Ala2 突变体 vs DPP-IV
**状态**：✅ refine-only 已完成（7 models）
**阻塞项**：同 Task 1，未传 `-native`
**下一步**：与 WT 同步带约束重跑
**备注**：PDB 中 `ALA` → `DAL`，Rosetta 自动应用 `D_AA` patch。

---

## 任务 3：GROMACS MD — WT 复合物（200 ns）
**状态**：🔄 生产 MD 运行中（~40 ns / 200 ns，GPU 0）
**阻塞项**：无
**下一步**：完成后运行分析脚本（RMSF、SASA、几何时序、氢键）
**备注**：从 9.44 ns checkpoint 续跑成功。dt=0.002，~101 ns/day。

---

## 任务 4：GROMACS MD — D-Ala2 突变体（200 ns）
**状态**：🔄 生产 MD 运行中（~12.9 ns / 200 ns，GPU 1）
**阻塞项**：无
**下一步**：继续运行到 200 ns
**备注**：前 ~12 ns 为弛豫期（从 phi=-6° 弛豫到 phi=-106°，beta-sheet 区域）。
  - dt=0.001（从 EM 开始稳定）
  - 起始结构 backbone 压平 bug 已确认根因（`fab()` + `pair_fit()`）
  - 系统已自我纠正，无需重启

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
- ✅ P0: 根目录散文件已清理
- ✅ P0: step3/ autosave 文件已清理
- ✅ P0: step2/ worker silent 已清理（保留 combined + archive）
- ✅ P0: `prepare_structures.py` / `prepare_structures_v2.py` `fab()` bug 已修复
- ✅ P0: `prepare_docking_start.py` 已添加 phi/psi 验证
- ✅ P1: 文档已更新（Phase1 v1.2, CJC-1295 Known Issues）
- ⏳ P0: refine-only 需带 `-native` + CoordinateConstraint 重跑
- ⏳ P1: `run_flexpepdock_parallel.py` 需添加 `-native` 支持
- ⏳ P2: `DAla2_MD_status_report.md` Finding 1 需改写（同义反复）

---

*最后更新：2026-05-14*
*维护者：Kimi Code CLI*
