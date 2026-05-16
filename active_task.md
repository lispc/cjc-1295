# 活跃任务跟踪

*最后更新：2026-05-16 晚间*

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
**状态**：✅ 已完成（~200 ns，有 9.4-66.7 ns 空缺）
**性能**：GPU 1, ~33 ns/day（优化后）
**关键发现（修正后）**：
- Ser630–Ala2C = 4.29±0.56 Å, 39% <4Å（66.7-200ns）
- 100-150ns 瞬态窗口：3.96Å, 66% <4Å
- 150ns 后脱结合：4.71Å, 0% <4Å
- **结论：不形成稳定催化几何，抑制是空间/变构效应**

---

## 任务 4：GROMACS MD — 天然 GHRH (L-Ala, 1-29) 全长
**状态**：🔄 已跑 ~80 ps（LAla_full_opt.xtc），**未续跑**
**注意**：此系统实际进度远低于预期。之前上下文中的 "62ns" 数据属于标签混淆时期，可能实际对应 D-Ala2。
**当前真实的 L-Ala 全长轨迹仅 80 ps**，远不足以做统计。

---

## 任务 5：D-Ala2 短肽 MD (GHRH 1-10)
**状态**：✅ 已完成（~200 ns）
**性能**：GPU 3, 69 ns/day
**关键发现（修正后）**：
- Ser630–Ala2C = 5.41±0.52 Å, **3.0% <4Å**
- **结论：短肽完全不结合，GHRH 11-29 对口袋锚定必不可少**

---

## 任务 6：L-Ala 短肽 MD (GHRH 1-10)
**状态**：🔄 GPU 2 运行中（LAla_short_v2, ~200 ns 目标）
**性能**：60 ns/day
**进度**：从 xtc 大小估算已跑较大量

---

## 任务 7：FEP (Alchemical Free Energy, D-Ala2 ↔ L-Ala2)
**状态**：✅ 已完成
**结果**：ΔG(D→L) = +0.83 ± 0.31 kJ/mol ≈ 0 → 结合亲和力无显著差异
**输出**：`bar_output.txt`, `fep_summary.xvg`

---

## 任务 8：手性标签审计与文档化
**状态**：✅ 完成
**文档**：`docs/CHIRALITY_CORRECTION.md`, `workspace/step3/CHIRALITY_LABELS.txt`

---

## 任务 9：催化几何分析（已修正）
**状态**：✅ 修正完成，GROMACS 数据已重新分析
**关键产出**：
- `workspace/step3/CATALYTIC_GEOMETRY_CORRECTED.md`
- `workspace/step3/analysis_corrected.ndx`
- **所有旧 "3.15Å" 结论已废弃**

---

## 任务 10：原子索引灾难修正
**状态**：✅ 完成
**根因**：`make_ndx r N` 使用 PDB 编号而非 resind，导致所有 GHRH 原子选错
**修正后原子**：
- Ser630 OG = 9562
- GHRH Tyr1 N = 11655
- GHRH Ala2 C = 11686
- GHRH Ala2 O = 11687
- Glu205 OE1 = 2787, OE2 = 2788

---

## 任务 11：OpenMM 复刻自然 GHRH (L-Ala2, 1-29)
**状态**：✅ 已完成（1 ns 基准测试 + 稳定性验证）
**GPU**：GPU 0
**速度**：NPT **44.4 ns/day**（vs GROMACS ~33 ns/day，快 34%）
**稳定性**：NVT 100 ps + NPT 100 ps 无 NaN
**文件**：`openmm_production.py`, `OPENMM_REPLICA.md`
**注意**：
- OpenMM thermostat = 300 K，GROMACS ref_t = 310 K（10K 差异）
- 溶剂化条件不完全相同（OpenMM 从 PDB 重建 vs GROMACS 原始溶剂化）
- 对于速度基准和物理等价性验证已足够

---

## 任务 12：可视化
**状态**：✅ 已完成
**产出**：`workspace/results/mechanism_composite.png` 等

---

## 2×2 矩阵进度（修正标签后）

| | **全长 (1-29)** | **短肽 (1-10)** |
|---|---|---|
| **D-Ala2** | ✅ ~200 ns（有空缺） | ✅ ~200 ns |
| **L-Ala** | ⚠️ 仅 ~80 ps | 🔄 ~200 ns 进行中 |

---

## 待办

1. **OpenMM 1 ns 轨迹分析**：等生产跑完后分析温度/压力/RMSD，与 GROMACS 80 ps 对比
2. **GROMACS L-Ala 全长续跑**：当前仅 80 ps，需决定是否继续跑至 200 ns
3. **四系统统一分析**：等 L-Ala 短肽完成后做 RMSF/SASA/氢键/盐桥/MM-PBSA
4. **Tyr1-N → Glu205 盐桥重新分析**：用正确原子（11655）重新计算
5. **GHRH 构象聚类**：分析 D-Ala2 全长 100-150ns 瞬态窗口的构象特征

---

*维护者：Claude Code*
