# AGENTS.md — CJC-1295 Phase 1 项目指引

## 项目概述

本项目是 CJC-1295（长效 GHRH 类似物）计算机辅助设计的第一阶段复现：
**在硅基环境中复现天然 GHRH(1-29) 的 L-Ala2 与 DPP-IV 催化口袋的"完美契合"现象**。

核心方案文档：`docs/Phase1_DPP4_GHRH_Perfect_Fit_Protocol.md`

---

## 重要文档位置

| 文档 | 路径 | 说明 |
|------|------|------|
| **活跃任务跟踪** | `/home/scroll/personal/cjc-1295/active_task.md` | 记录正在运行的任务状态和阻塞项 |
| **项目日志** | `/home/scroll/personal/cjc-1295/project_log.md` | 记录关键发现、决策和结果 |
| **第一阶段方案** | `/home/scroll/personal/cjc-1295/docs/Phase1_DPP4_GHRH_Perfect_Fit_Protocol.md` | 原始实验方案文档 |
| **复现总览** | `/home/scroll/personal/cjc-1295/docs/CJC-1295 Computational Design Reproduction.md` | 完整复现指南（含后续阶段） |

---

## 重要软件与环境

### Conda 环境（miniforge3）

| 环境名 | 路径 | 用途 | 关键工具 |
|--------|------|------|----------|
| `base` | `/home/scroll/miniforge3` | 默认环境 | PyMOL, OpenMM, Biopython |
| `rosetta` | `/home/scroll/miniforge3/envs/rosetta` | Rosetta 2026.15 | FlexPepDocking, docking_protocol, combine_silent |
| `gmx` | `/home/scroll/miniforge3/envs/gmx` | GROMACS 2026.0 | gmx, GPU 加速 MD |
| `cgas-md` | `/home/scroll/miniforge3/envs/cgas-md` | 已有完整 MD 工具链 | tleap, pdb4amber, AmberTools 24.8 |
| `boltz` | `/home/scroll/miniforge3/envs/boltz` | Boltz-2 结构预测 | boltz |

### 关键软件二进制

| 软件 | 路径 |
|------|------|
| **Rosetta FlexPepDock** | `/home/scroll/miniforge3/envs/rosetta/bin/FlexPepDocking` |
| **Rosetta 数据库** | `/home/scroll/miniforge3/envs/rosetta/database` |
| **GROMACS** | `/home/scroll/miniforge3/envs/gmx/bin/gmx` |
| **PyMOL** | `/home/scroll/miniforge3/bin/pymol` |
| **AmberTools tleap** | `/home/scroll/miniforge3/envs/cgas-md/bin/tleap` |

### GPU 资源

- 4× NVIDIA RTX 3090 (24 GB VRAM)
- CUDA Driver 580.95.05, CUDA Runtime 12.9
- GROMACS 已编译 CUDA 支持

---

## 工作目录结构

```
/home/scroll/personal/cjc-1295/
├── docs/                          # 原始方案文档
├── workspace/                     # 所有计算数据和结果
│   ├── structures/                # 下载的 PDB 文件（1NU8, 7CZ5, 7V9M）
│   ├── step0/                     # Step 0: YAD 三肽快速验证
│   ├── step1/                     # Step 1: 结构准备
│   ├── step2/                     # Step 2: Rosetta 对接
│   ├── step3/                     # Step 3: GROMACS MD 配置
│   ├── results/                   # 结果输出（图、报告）
│   ├── scripts/                   # 所有脚本
│   └── logs/                      # 日志
├── active_task.md                 # 活跃任务跟踪
├── project_log.md                 # 项目日志
└── AGENTS.md                      # 本文件
```

### 关键数据文件

| 文件 | 路径 | 说明 |
|------|------|------|
| DPP-IV 受体（clean） | `workspace/step1/DPP4_clean.pdb` | 1NU8 链 B，去除溶剂/NAG |
| DPP-IV + Diprotin A | `workspace/step1/DPP4_with_diprotinA.pdb` | 含配体参考 |
| GHRH(1-29) WT | `workspace/step1/GHRH_1-29.pdb` | 基于 7CZ5 构建的 29 残基 |
| GHRH(1-29) D-Ala2 | `workspace/step1/GHRH_1-29_DAla2.pdb` | 手性翻转的突变体 |
| YAD 三肽 | `workspace/step0/YAD_tripeptide.pdb` | 用于 Step 0 验证 |
| YAD 对齐到 DPP-IV | `workspace/step0/YAD_aligned_to_DPP4.pdb` | 已对齐到活性口袋 |
| 静电势图 | `workspace/results/step0_electrostatics.png` | APBS 电势 + YAD 叠加 |
| Rosetta 预处理结构 | `workspace/step2/prepacked_DPP4_GHRH_start_0001.pdb` | 对接起始构象 |

---

## 已知问题与注意事项

1. **GROMACS 2026.0 conda 包的 amber14sb aminoacids.hdb 有上游 bug**：文件内有多处损坏（ALA 条目重复、CSER 混入 THR 规则）。已用 amber19sb 的干净版本替换。
2. **Rosetta 种子参数**：FlexPepDocking 不认 `-seed`，需用 `-constant_seed` + `-jran`。
3. **1NU8 是二聚体**：Diprotin A（chain D）结合在 **chain B** 上，不是 chain A。预处理时必须保留 chain B。
4. **D-Ala 力场**：Amber14sb 无预装 D-氨基酸。已在 `aminoacids.rtp` 和 `aminoacids.hdb` 中手动添加 `DALA` 条目（复制自 ALA），验证通过。

---

## 已知陷阱（Known Traps）

> 以下陷阱已被踩过至少一次，写入此处防止重复犯错。

### Trap 1: PyMOL `cmd.fab("YAD", ss=0)` 压平 backbone
- **症状**：GHRH N-端 φ≈0°, ψ≈0°（Ramachandran 禁区），GROMACS segfault
- **根因**：`fab()` 生成线性肽，`pair_fit()` 强制 backbone 匹配
- **修复**：从晶体结构提取天然 YAD 构象（见 `prepare_structures.py`）
- **检测**：`prepare_docking_start.py` 现在自动验证 phi/psi

### Trap 2: Rosetta `-native` 参数漏传
- **症状**：silent 文件中 `startRMSall = 0.000` 对所有模型，RMS 列完全无意义
- **根因**：`-native` 不是可选参数；缺省时 Rosetta 用输入结构自身作为参考
- **修复**：所有 FlexPepDock 调用必须包含 `-native <ref.pdb>`
- **检测**：检查 silent 文件 header 中 `startRMSall` 是否全为 0

### Trap 3: Rosetta refine-only 无约束时长肽漂移
- **症状**：`-pep_refine` 后 Ser630→Ala2 距离从 2.8 Å 增加到 5–10 Å
- **根因**：29 残基长肽在无约束 refine 中 backbone RMS 漂移 1.6–7 Å
- **修复**：对 GHRH 1–5 CA 加 CoordinateConstraint（sd=0.5 Å）
- **替代**：直接用未 refine 的起始结构做 MD（催化几何更好）

### Trap 4: D-Ala2 segfault 误判为 GROMACS bug
- **症状**：D-Ala2 MD dt=0.002 时 segfault，WT 正常
- **根因**：起始结构 backbone 在禁区（Trap 1），不是 improper dihedral
- **验证**：dt=0.001 允许系统缓慢弛豫；10+ ns checkpoint 后 dt=0.002 也稳定
- **教训**：先检查起始结构几何，再怀疑软件 bug

---

*维护者：Kimi Code CLI*
*最后更新：2026-05-14*
