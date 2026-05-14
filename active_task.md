# 活跃任务跟踪

## 任务 1：Rosetta FlexPepDock — WT GHRH(1-29) vs DPP-IV
**状态**：✅ 已完成（refine-only，960 models 废弃；refine 6 models 保留）
**结论**：长肽无约束 refine 会漂移 1.6–7 Å。起始结构（来自 7CZ5 对齐）4/4 PASS，直接用于 MD。

---

## 任务 2：Rosetta FlexPepDock — D-Ala2 突变体 vs DPP-IV
**状态**：✅ 已完成（refine-only，7 models）
**结论**：同 WT，无约束 refine 破坏催化几何。起始结构用于 MD。

---

## 任务 3：GROMACS MD — WT 复合物（200 ns，无约束）
**状态**：🔄 GPU 0 运行中（~55 ns / 200 ns）
**关键发现**：P0a 分析显示 WT 有 3.58% 帧处于 2D productive pose（距离+角度同时满足）。
**下一步**：继续跑完 200 ns。

---

## 任务 4：GROMACS MD — D-Ala2 突变体（200 ns，无约束）
**状态**：🔄 GPU 1 运行中（~18 ns / 200 ns）
**关键发现**：P0a 分析显示 D-Ala2 **0.00%** 帧处于 2D productive pose。攻击角系统性偏低（77.7°）。
**下一步**：继续跑完 200 ns。

---

## 任务 5：WT + N 端 POSRES MD（200 ns）
**状态**：⏹️ 已停止（GPU 2 已释放）
**结论**：CA-only POSRES 不足以维持催化几何。约束/截断把距离锁定在 ~4.9 Å（productive zone 外）。

---

## 任务 6：短肽 GHRH(1-10) 验证 MD（200 ns，无约束）
**状态**：🔄 GPU 3 运行中（~6 ns / 200 ns）
**关键发现**：N 端盐桥保持，但催化距离锁定在 ~4.9 Å。角度很好（83.3°）。
**下一步**：继续跑完 200 ns。

---

## 任务 7：P0a 数据分析 ✅
**状态**：已完成
**产出**：
- WT 21.96% distance-only productive → **3.58% 2D productive**
- D-Ala2 1.10% distance-only → **0.00% 2D productive**
- 攻击角数据证实 D-Ala2 手性翻转同时破坏距离和角度

---

## 任务 8：P0b Rosetta Constrained Refine 🔄
**状态**：运行中（16 workers，各 200 models）
**预计完成**：~40-60 分钟
**配置**：CoordinateConstraint GHRH 1-5 CA, σ=0.5 Å, weight=10.0

---

## 任务 9：P1-FEP 筹备 ⏳
**状态**：准备中
**目标**：计算 D-Ala2 vs WT 的结合自由能差 ΔΔG_bind
**体系**：DPP-IV + GHRH(1-10) 短肽（~100k 原子）
**方法**：GROMACS alchemical FEP（L-Ala2 ↔ D-Ala2 手性翻转）
**难点**：手性翻转不是常规 vdW 关闭，需特殊处理 dual-topology 或 coordinate transformation

---

## 当前阻塞项

1. **P0b 等待中**：Rosetta refine 完成后才能分析
2. **P1-FEP 方案未定**：需要确认 GROMACS 2026 中手性翻转 FEP 的具体实现路径

---

*最后更新：2026-05-14 晚间*
*维护者：Claude Code*
