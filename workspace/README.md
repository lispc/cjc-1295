# Workspace 目录说明

本目录存放 CJC-1295 Phase 1 复现的全部计算数据、脚本和结果。

## 目录结构

```
workspace/
├── structures/          # 原始 PDB 下载文件
├── step0/               # Step 0: YAD 三肽快速验证
├── step1/               # Step 1: 大分子结构准备
├── step2/               # Step 2: Rosetta 柔性对接
├── step3/               # Step 3: GROMACS MD 配置与输入
├── results/             # 结果输出（图、表格、报告）
├── scripts/             # 所有自动化脚本
└── logs/                # 运行日志
```

---

## 各目录详情

### `structures/`

| 文件 | 来源 | 说明 |
|------|------|------|
| `1NU8.pdb` | RCSB PDB | DPP-IV 与 Diprotin A 共晶结构（2.00 Å） |
| `7CZ5.pdb` | RCSB PDB | GHRH-GHRHR-Gs 复合物 Cryo-EM |
| `7V9M.pdb` | RCSB PDB | GHRH-GHRHR 变体复合物 Cryo-EM |

---

### `step0/` — 快速验证

| 文件 | 说明 |
|------|------|
| `YAD_tripeptide.pdb` | PyMOL fab 构建的 Tyr-Ala-Asp 三肽 |
| `YAD_aligned_to_DPP4.pdb` | pair_fit 对齐到 Diprotin A 的 YAD |
| `DPP4_with_YAD.pdb` | DPP-IV + YAD 复合物（用于几何分析） |
| `DPP4_with_YAD_only.pdb` | 去除 Diprotin A 后的版本 |
| `DPP4_potential.dx` | APBS 计算的静电势图 |
| `apbs.in` | APBS 输入文件 |
| `visualize_apbs.pml` | PyMOL 可视化脚本 |

---

### `step1/` — 结构准备

| 文件 | 说明 |
|------|------|
| `DPP4_clean.pdb` | DPP-IV 单体（chain B，去除溶剂/NAG） |
| `DPP4_with_diprotinA.pdb` | DPP-IV + Diprotin A（空间参考） |
| `GHRH_1-29_from_7CZ5.pdb` | 从 7CZ5 提取的 28 残基片段 |
| `GHRH_1-29.pdb` | 完整 29 残基（28 + 建模 Arg29） |
| `GHRH_1-29_DAla2.pdb` | D-Ala2 突变体（CB 镜像翻转） |

---

### `step2/` — Rosetta 对接

| 文件 | 说明 |
|------|------|
| `DPP4_GHRH_start.pdb` | 对接起始构象（GHRH N 端已预对齐） |
| `prepacked_DPP4_GHRH_start_0001.pdb` | Rosetta 预打包后的受体-配体复合物 |
| `GHRH_DPP4_dock_worker_*.silent` | 各 worker 的 silent 输出（共 16 个） |
| `GHRH_DPP4_dock_combined.silent` | 合并后的 silent 文件（对接完成后生成） |
| `parallel_dock.log` | 并行脚本主日志 |
| `worker_*.log` | 各 worker 的 Rosetta 日志 |
| `DPP4_prepack.silent` | 预打包阶段的 silent 输出 |

---

### `step3/` — MD 配置

| 文件 | 说明 |
|------|------|
| `ions.mdp` | 离子添加（genion）用空参数文件 |
| `em.mdp` | 能量最小化参数 |
| `nvt.mdp` | NVT 平衡（310 K，100 ps，位置限制） |
| `npt.mdp` | NPT 平衡（310 K，1 bar，100 ps，位置限制） |
| `md.mdp` | 生产模拟（200 ns，2 fs 步长，GPU 加速） |

---

### `results/` — 结果输出

| 文件 | 说明 |
|------|------|
| `step0_electrostatics.png` | Step 0 静电势可视化（2400×2400，300 dpi） |
| `step0_electrostatics.pse` | PyMOL 会话文件 |
| （待生成）| 对接最优构象聚类结果 |
| （待生成）| RMSD/RMSF/SASA 分析图 |
| （待生成）| 关键距离/角度时序图 |
| （待生成）| MM-PBSA 结合自由能结果 |

---

### `scripts/` — 脚本

| 文件 | 功能 |
|------|------|
| `prepare_structures_v2.py` | Step 1 结构预处理（PyMOL） |
| `build_GHRH_1-29.py` | 构建完整 GHRH(1-29) |
| `build_DAla_mutant.py` | 构建 D-Ala2 突变体 |
| `realign_YAD.py` | YAD 三肽对齐到 DPP-IV |
| `step0_analysis.py` | Step 0 几何参数分析 |
| `prepare_docking_start.py` | 准备 Rosetta 对接起始构象 |
| `run_flexpepdock_parallel.py` | 16 进程并行 FlexPepDock |
| `run_apbs.py` | APBS 静电势计算工作流 |
| `verify_openmm_dala.py` | OpenMM D-Ala 支持验证 |

---

## 文件命名约定

- `DPP4_*`：DPP-IV（二肽基肽酶-IV，即 DPP-IV/CD26）
- `GHRH_*`：生长激素释放激素（Growth Hormone-Releasing Hormone）
- `YAD`：Tyr-Ala-Asp 三肽
- `DAla2` / `D-Ala2`：第 2 位丙氨酸手性翻转为 D-型
- `*_start.pdb`：Rosetta 对接/模拟的起始构象
- `*_prepack.pdb` / `*_prepack.silent`：Rosetta 预打包输出

### `step3/` — MD 配置（更新 2026-05-15）

| 文件 | 说明 |
|------|------|
| `md.tpr` / `md*.trr` | **CJC-1295 (D-Ala2, 1-29)**, 130+ ns |
| `DAla2_production.tpr` / `DAla2_md.xtc` | **天然 GHRH (L-Ala, 1-29)**, dt=0.002 重启 |
| `short_peptide_md.*` | **D-Ala2 短肽 (1-10)**, 143+ ns |
| `short_peptide_LAla_md.*` | **L-Ala 短肽 (1-10)**, 新启动 |
| `lambda_00/` – `lambda_10/` | **FEP (D-Ala2 ↔ L-Ala2)**, 11 windows × 5 ns |
| `bar_output.txt` / `fep_summary.xvg` | FEP BAR 分析结果 |

重要：所有旧文件中 "wt"/"dala2" 标签与实际手性互换，详见 `CHIRALITY_LABELS.txt`。

---

*维护者：Claude Code*
*最后更新：2026-05-15 下午*
