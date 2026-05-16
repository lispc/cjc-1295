# 活跃任务跟踪

*最后更新：2026-05-16 上午*

---

## 任务 1：Rosetta FlexPepDock — WT GHRH(1-29) vs DPP-IV
**状态**：✅ 已完成
**结论**：长肽无约束 refine 会漂移。起始结构（7CZ5 对齐）4/4 PASS，直接用于 MD。

---

## 任务 2：Rosetta FlexPepDock — D-Ala2 突变体 vs DPP-IV
**状态**：✅ 已完成
**结论**：同 WT，无约束 refine 破坏催化几何。起始结构用于 MD。

---

## 任务 3：GROMACS MD — CJC-1295 (D-Ala2, 1-29) 全长
**状态**：🔄 GPU 0 运行中（**184.5 ns** / 200 ns）
**性能**：74 ns/day | **预计今天完成**
**关键发现**：Ser630–Ala2C ≈ 3.15 Å, 93.8% <4Å。紧密结合，攻击角错误 → 抑制机制

---

## 任务 4：GROMACS MD — 天然 GHRH (L-Ala, 1-29) 全长
**状态**：🔄 GPU 1 运行中（**53.7 ns** / 200 ns）
**性能**：19 ns/day（dt=0.001）| 预计 **5/24** 完成
**注意**：dt=0.002 加速尝试失败（GPU 利用率降至 13%，性能更差）。dt=0.001 是唯一稳定配置。
**关键发现**：Ser630–Ala2C ≈ 6.31 Å，肽漂离口袋 → 底物解离

---

## 任务 5：D-Ala2 短肽 MD (GHRH 1-10)
**状态**：🔄 GPU 3 运行中（**184.1 ns** / 200 ns）
**性能**：69 ns/day | **预计今天完成**
**关键发现**：Glu205–Tyr1 盐桥 92.4% 稳定，催化几何 34% <4Å

---

## 任务 6：L-Ala 短肽 MD (GHRH 1-10)
**状态**：🔄 GPU 2 运行中（**25.9 ns** / 200 ns, v2 优化版）
**性能**：60 ns/day（优化后，原 25 ns/day）| 预计 **5/19** 完成
**优化历史**：rc=0.9 + nstlist=50 + rlist=1.1 + ntmp=6，较原配置 +140%
**目的**：与 D-Ala2 短肽形成直接对照，孤立 Ala2 手性效应
**注意**：v1（10.7 ns）已废弃重置；dt=0.004 尝试 crash（需 HMR）

---

## 任务 7：FEP (Alchemical Free Energy, D-Ala2 ↔ L-Ala2)
**状态**：✅ 已完成
**结果**：ΔG(D→L) = +0.83 ± 0.31 kJ/mol，≪ kT → 结合亲和力无显著差异
**输出**：`bar_output.txt`, `fep_summary.xvg`

---

## 任务 8：手性标签审计与文档化
**状态**：✅ 完成
**文档**：`docs/CHIRALITY_CORRECTION.md`, `workspace/step3/CHIRALITY_LABELS.txt`

---

## 任务 9：催化几何初步分析
**状态**：✅ 初步完成，待全部 MD 完成后完善
**关键产出**：
- D-Ala2 (CJC-1295): 3.15±0.52 Å, 93.8% <4Å, **0% productive** (angle wrong)
- L-Ala (native): 6.31±1.48 Å, 8.3% <4Å, **3.58% productive** (rare correct angle)
- FEP ΔG ≈ 0 → 抑制机制纯属几何效应

---

## 任务 10：可视化
**状态**：✅ 已完成
**产出**：
- `workspace/results/mechanism_composite.png` — PyMOL 结构图 + 数据图三面板
- `workspace/step3/panel_DALA.png` / `panel_LALA.png` — 干净 PyMOL 渲染

---

## 2×2 矩阵进度

| | **全长 (1-29)** | **短肽 (1-10)** |
|---|---|---|
| **D-Ala2** | 184.5 ns ✅ **今天完** | 184.1 ns ✅ **今天完** |
| **L-Ala** | 53.7 ns 🔄 ~8天 | 25.9 ns 🔄 ~3天 |

---

## 待办

- GPU0/3 完成后立即做 D-Ala2 全长 vs 短肽对比分析
- 全部完成后四系统统一分析（RMSF/SASA/氢键/盐桥/MM-PBSA）
- GPU1 是长杆，无法加速，耐心等

---

*维护者：Claude Code*
