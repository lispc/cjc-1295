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
| `gmx` (OpenMM) | `/home/scroll/miniforge3/envs/gmx` | OpenMM 8.5.1 已安装 | `openmm`, `parmed` |

### 关键软件二进制

| 软件 | 路径 |
|------|------|
| **Rosetta FlexPepDock** | `/home/scroll/miniforge3/envs/rosetta/bin/FlexPepDocking` |
| **Rosetta 数据库** | `/home/scroll/miniforge3/envs/rosetta/database` |
| **GROMACS** | `/home/scroll/miniforge3/envs/gmx/bin/gmx` |
| **PyMOL** | `/home/scroll/miniforge3/bin/pymol` |
| **AmberTools tleap** | `/home/scroll/miniforge3/envs/cgas-md/bin/tleap` |
| **OpenMM** | Python: `import openmm` (conda env `gmx`) |

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

### Trap 5: 29 残基长肽在无约束 MD 中从口袋漂移
- **症状**：WT ~40 ns 后 0/4 PASS，D-Ala2 ~12.9 ns 后 1/4 PASS（起始均为 4/4）
- **根因**：GHRH(1-29) C 端完全暴露，热涨落拉动 N 端离开口袋
- **验证**：提取 MD 最后帧，Ser630→Ala2 C 从 2.82 Å 漂到 4.84 Å（WT）/ 6.40 Å（D-Ala2）
- **教训**：长肽-受体复合物 MD 必须加 N 端约束（POSRES）或 PLUMED，否则 catalytic geometry 无法维持

### Trap 6: CA-only POSRES 不足以维持催化几何
- **症状**：POSRES on GHRH 1-5 CA (100 kJ/mol/nm²) 跑 2.8 ns，Ser630→Ala2 C = 4.93 Å（起始 2.82 Å）
- **根因**：CA 约束只固定 Cα 位置，骨架可绕 CA 旋转，C=O 取向不受控
- **验证**：短肽（无 C 端拖动）Tyr1N→Glu205 保持在 2.76 Å（盐桥完好），但 Ser630→Ala2C 仍 4.87 Å — 说明即使 N 端锚定了，催化攻击距离也不自动满足
- **修复**：使用全骨架重原子约束（N/CA/C/O, 500-1000 kJ）或 PLUMED 距离约束 `dist(Ser630 OG, Ala2 C) = 3.0 ± 0.3 Å`

### Trap 7: 混合 POSRES + C-rescale/P-R 压力耦合需特殊处理
- **症状**：grompp 警告 "combining position restraints with pressure coupling can lead to instabilities"
- **根因**：GROMACS 2026 要求使用 `-r` 指定参考坐标、使用 C-rescale（非 P-R）、加 `refcoord_scaling = com`
- **修复**：grompp 加 `-r ref.gro`，mdp 用 `pcoupl = C-rescale`

### Trap 8: D/L-氨基酸手性标签容易混淆
- **症状**：项目中所有 WT/D-Ala2 标签与实际手性完全相反，跑了一周才发现
- **根因**：`build_DAla_mutant.py` 的镜像翻转逻辑写反；7CZ5 模板本身已是 D-Ala2，翻转后生成 L-Ala 但标签写 D-Ala2
- **检测**：检查 Ala2 的 N-CA-CB-C 二面角符号：χ > 0 → L-Ala，χ < 0 → D-Ala
- **教训**：任何手性操作后必须验证 χ 二面角符号，不能信任标签。详见 `docs/CHIRALITY_CORRECTION.md`

### Trap 9: FEP dual-topology λ=0 需真实 bonded 参数
- **症状**：FEP λ₀₀ 窗口 crash，A-state (dummy atom) 的 bonded 参数全为 0
- **根因**：GROMACS 2026 不允许 bonds/angles/dihedrals 力常数为 0，即使 dummy atom 也不接受
- **修复**：给 dummy atom 的 bonds/angles/dihedrals 赋予与真实原子相同的力常数

### Trap 10: 不能通过提取帧+新建 TPR 来改变 dt
- **症状**：提取 checkpoint 坐标，用新 dt=0.002 建 TPR 重启 → GPU 利用率 13%，性能 9 ns/day
- **根因**：提取的帧不含速度，系统需重新热化，且盒子/压力状态与 P-R barostat 不匹配
- **正确做法**：用 `gmx convert-tpr` 修改已有 TPR 的 dt；或继续用 checkpoint 续跑原 dt
- **教训**：别试图从 checkpoint 中提取帧再重启来改 dt——几乎总会变慢

### Trap 11: Amber14sb + TIP3P 不能用 4 fs timestep
- **症状**：dt=0.004 时 CUDA illegal memory access crash
- **根因**：Amber14sb 未做 Hydrogen Mass Repartitioning (HMR)，氢原子质量太轻，4 fs 步长导致积分不稳定
- **修复**：如要 4 fs，需先对拓扑做 HMR 处理；或在 mdp 中用 dt=0.002（安全上限）

### Trap 12: 97k 原子体系在 GROMACS 2026 + RTX 3090 上性能上限约 70 ns/day
- **观察**：短肽体系（97k 原子）无论如何调参数，性能上限 ~60-70 ns/day
- **根因**：PME FFT 网格（96³）相对 compute 占比过高（PME load ~0.25），GPU 大部分时间在等 FFT
- **对比**：325k 体系（144³ 网格）也能跑 74 ns/day——PME/PP 比例对中等体系不友好
- **经验**：小体系不一定更快；该瓶颈无法通过调参突破

### Trap 13: `make_ndx r N` 使用 PDB 残基编号而非系统索引
- **症状**：所有催化几何分析完全错误（短肽 "3.15Å" 实际是 5.41Å）
- **根因**：`gmx make_ndx` 的 `r N` 使用 PDB 残基编号；`-merge all` 后 GHRH 保留 chain B 编号 `r 1-29`，DPP-IV 保留 chain A 编号 `r 39-766`
- **检测**：`gmx dump -s md.tpr | grep resind` 查看系统索引；`make_ndx` 中用 `l` 列出残基确认编号
- **教训**：合并链拓扑中，PDB 编号和系统索引完全脱耦。必须用 `gmx dump` 确认原子映射，不能信任 `make_ndx` 的直觉

### Trap 14: GROMACS → OpenMM 坐标转换时的包装问题
- **症状**：OpenMM 能量 34.8 亿 kJ/mol，键长 16.6 nm
- **根因**：`gmx editconf` 输出 PDB 时坐标被包装到盒子边界；OpenMM 直接读取 wrapped 坐标计算键能
- **修复**：`gmx trjconv -pbc whole` 恢复分子完整性后再转 PDB
- **教训**：任何跨引擎格式转换都必须做 `whole` 恢复，不能直接转 raw 轨迹

---

*维护者：Claude Code*
*最后更新：2026-05-16 上午*

*维护者：Claude Code*
*最后更新：2026-05-15 下午*
*最后更新：2026-05-14 下午*
