# 第一阶段复现方案：DPP-IV 与 GHRH(1-29) "完美契合" 的计算机模拟

## 1. 目标与科学问题

**核心目标**：在硅基环境中复现 ConjuChem 团队早期的关键发现——天然 GHRH(1-29) 的第 2 位 L-丙氨酸（L-Ala2）与 DPP-IV 催化口袋在空间构象和静电势分布上形成"完美契合"，从而导致肽键在数分钟内被水解断裂。

**要回答的具体问题**：
1. L-Ala2 在空间上如何精确嵌入 DPP-IV 的 S1/S2 亚口袋？
2. GHRH 的 N 端（Tyr1-Ala2-Asp3）与 DPP-IV 的 Glu205/Glu206 双电 motif 形成怎样的静电互补？
3. 催化几何是否满足形成四面体过渡态的条件（Ser630 Oγ 到 Ala2 C=O 碳的距离 ~3.5 Å，角度接近 90°）？
4. 这种"契合"如何量化？（结合自由能、氢键网络、RMSD、静电势映射）

---

## 2. 生物化学背景

### 2.1 DPP-IV 的分子识别机制
- **DPP-IV**（二肽基肽酶-IV，CD26，Uniprot P27487）是 766 aa 的丝氨酸蛋白酶，以同源二聚体形式存在于血浆中。
- **催化三联体**：Ser630（亲核体）、Asp708（酸碱催化剂）、His740（广义碱）。
- **底物特异性**：专一性切割多肽 N 端倒数第二位为 **Pro** 或 **Ala** 的肽键（X-Pro↓ 或 X-Ala↓）。
- **关键功能区域**：
  - **S2 口袋（N 端识别区）**：Glu205、Glu206、Arg125。负责锚定底物 N 端的游离氨基（-NH₃⁺），形成盐桥。
  - **S1 口袋（疏水/催化核心区）**：Tyr631、Val656、Trp659、Tyr662、Tyr666、Val711。容纳 P1 位残基（此处为 Ala2）的侧链。
  - **氧负离子洞（Oxyanion hole）**：Tyr547、Ser631。稳定四面体过渡态的负电荷氧。
  - **S3 口袋**：Ser209、Phe357、Arg358。与 P2' 位残基（此处为 Asp3）相互作用。

### 2.2 GHRH(1-29) 的脆弱性
- **天然人源 GHRH(1-29) 序列**：
  ```
  Y-A-D-A-I-F-T-N-S-Y-R-K-V-L-G-Q-L-S-A-R-K-L-L-Q-D-I-M-S-R
  1 2 3 4 5 6 7 8 9 ...
  ```
- **致命弱点**：N 端第 2 位为 **L-Ala**，完全符合 DPP-IV 的底物偏好。DPP-IV 将切割 **Tyr1-Ala2** 与 **Asp3** 之间的肽键（即 Ala2↓Asp3）。

---

## 3. 计算方案总览

本方案采用 **"结构准备 → 柔性对接 → 分子动力学稳定 → 结合模式深度分析"** 的四步流水线。

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  Step 1: 结构准备 │ → │ Step 2: 分子对接 │ → │ Step 3: MD 模拟  │ → │ Step 4: 深度分析 │
└─────────────────┘    └─────────────────┘    └─────────────────┘    └─────────────────┘
       │                      │                      │                      │
       ▼                      ▼                      ▼                      ▼
   • DPP-IV 结构获取       • 定义结合盒             • 溶剂化/加离子          • 催化几何分析
   • GHRH(1-29)建模       • AutoDock/RosettaDock  • 能量最小化             • 静电势映射
   • 结构预处理           • 对接打分               • NVT/NPT 平衡           • 氢键/盐桥统计
   • 力场参数化           • 姿态聚类               • 生产模拟 (500ns+)      • 结合自由能计算
                                                                               • RMSF/SASA
```

---

## 4. 详细实施步骤

### Step 1: 大分子与多肽结构准备

#### 1.1 DPP-IV 受体结构获取与预处理
**推荐 PDB 结构**：
| PDB ID | 分辨率 | 特点 | 适用性 |
|--------|--------|------|--------|
| **1NU8** | 2.00 Å | 与底物类似物 Diprotin A (Ile-Pro-Ile) 共晶 | ⭐ 首选，含底物结合参考 |
| 1NU6 | 2.00 Å | apo 形式（无二聚体界面配体） | 可用 |
| 2ONC | 2.49 Å | 与抑制剂共晶 | 口袋构象可能关闭 |
| 4A5S | 2.20 Å | 与肽类配体共晶 | 可用 |

**预处理流程**（PyMOL / ChimeraX）：
```bash
# 下载 1NU8
wget https://files.rcsb.org/download/1NU8.pdb

# PyMOL 预处理命令
cmd.load("1NU8.pdb")
cmd.remove("solvent")          # 去除结晶水
cmd.remove("resn HOH")         # 去除水
cmd.remove("resn SO4")         # 去除硫酸根等结晶添加剂
# 保留一条链（DPP-IV 为二聚体，通常只取一个单体进行对接）
cmd.remove("chain B")          # 去除链 B（如果存在）
cmd.save("DPP4_clean.pdb")
```

**关键检查点**：
- 确认 **Ser630、Asp708、His740** 的侧链完整且几何合理。
- 检查 **Glu205、Glu206** 是否朝向活性口袋开口（它们是 N 端识别的关键）。
- 若使用 apo 结构（如 1NU6），建议以 1NU8 为参考进行结构叠加，确保口袋处于开放构象。

#### 1.2 GHRH(1-29) 多肽结构构建
由于 GHRH(1-29) 没有独立的晶体结构，需要从头建模：

**方案 A：基于 GHRHR 复合物提取（推荐）**
- 从 PDB **7CZ5** 或 **7V9M**（GHRH-GHRHR-Gs 复合物 Cryo-EM 结构）中提取 GHRH 的 N 端前 10-15 个残基。
- 在 PyMOL 中保留前 29 个残基，C 端加帽（NH₂ 或 NME）。

**方案 B：从头建模（备用）**
- 使用 **PEP-FOLD3/4** 在线服务器（https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/）提交 GHRH(1-29) 序列进行从头预测。
- 或使用 **AlphaFold3**（若有本地部署）。
- 或使用 **Rosetta** 的 `minirosetta` 进行肽段折叠。

**序列确认**（人源 GHRH(1-29)）：
```
>sp|P01286|GHRH_HUMAN (1-29)
YADAIFTNSYRKVLGQLSARKLLQDIMSR
```

**结构预处理**：
- N 端保持游离氨基（NH₃⁺，pH 7.4 下质子化）。
- C 端保持游离羧基（COO⁻）。
- 所有可电离侧链按生理 pH 设置：Asp/Glu 去质子化（-），Arg/Lys 质子化（+），Tyr 质子化，His 可视为默认质子化（HID/HIE）。

#### 1.3 力场参数化
- **DPP-IV**：使用标准 Amber ff14SB 或 CHARMM36m 力场。
- **GHRH(1-29)**：作为标准多肽，可直接被上述力场覆盖。
- **特殊处理**：无（本阶段不涉及非天然氨基酸或共价修饰）。

---

### Step 2: 分子对接（锁定"完美契合"的初始姿态）

目标：将 GHRH(1-29) 对接到 DPP-IV 活性位点，获得 L-Ala2 嵌入催化口袋的合理初始构象。

#### 2.1 对接策略选择

由于 GHRH(1-29) 是 29 个氨基酸的长肽，传统小分子对接工具（如 AutoDock Vina）的构象搜索空间不足。推荐以下策略：

**推荐方案：HADDOCK 肽-蛋白对接（首选）**
- **HADDOCK**（https://wenmr.science.uu.nl/haddock2.4/）专为肽-蛋白对接优化。
- 支持定义活性残基（active/passive residues）。
- 可将 **Ser630、Asp708、His740、Glu205、Glu206** 定义为活性位点残基。
- 支持柔性肽段采样（it1/itw 阶段）。

**备选方案：Rosetta FlexPepDock**
- 专为柔性肽段对接设计（详见文档第四阶段）。
- 适合后续做大规模构象采样。

**简化方案：AutoDock Vina（仅用于 N 端三肽片段验证）**
- 若计算资源有限，可先对接 **Tyr-Ala-Asp 三肽** 到 DPP-IV 活性位点，验证 Ala2 的结合模式。

#### 2.2 HADDOCK 对接参数设置
```
# 活性残基定义（DPP-IV 侧）
Active residues: 125, 205, 206, 547, 630, 631, 662, 708, 740
Passive residues: 209, 357, 358, 659, 666, 711

# GHRH(1-29) 侧（强制 N 端进入口袋）
Active residues: 1, 2, 3  (Tyr1, Ala2, Asp3)
Passive residues: 4-10

# 参数
Peptide sampling: flexible (C-terminally flexible)
Number of models: 1000 it0, 200 it1, 200 water
```

#### 2.3 对接结果筛选标准
从生成的 200 个水细化模型中，筛选出满足以下**全部**条件的姿态：

| 筛选指标 | 阈值 | 物理意义 |
|----------|------|----------|
| **HADDOCK score** | Top 10% | 整体结合亲和力 |
| **Ala2 Cα 到 S1 口袋中心距离** | < 5 Å | 嵌入催化核心区 |
| **Ser630 Oγ 到 Ala2 C=O 碳距离** | 2.5–4.5 Å | 亲核攻击距离 |
| **N-端 NH₃⁺ 到 Glu205/Glu206 距离** | < 4 Å | 盐桥形成 |
| **Ala2 C=O 氧到 Tyr547 OH 距离** | < 4 Å | 氧负离子洞稳定 |

**关键观察**：最优模型中，Ala2 的甲基侧链应紧密贴合由 Tyr631、Val656、Trp659、Tyr662、Tyr666 形成的疏水 S1 口袋，同时骨架羰基氧指向氧负离子洞。

---

### Step 3: 全原子分子动力学（MD）模拟

目标：验证对接获得的"契合"姿态在显式溶剂和生理温度下是否稳定，观察动态过程中的氢键网络和催化几何波动。

#### 3.1 模拟体系构建（GROMACS 流程）
```bash
# 1. 选择最佳对接模型，保存为 complex.pdb
# 2. 使用 Amber ff14SB 力场（通过 GROMACS 的 amber99sb-ildn.ff 或自建）
gmx pdb2gmx -f complex.pdb -o complex.gro -p topol.top \
    -i posre.itp -ff amber99sb-ildn -water tip3p -ignh

# 3. 定义模拟盒（十二面体，边缘距蛋白 1.2 nm）
gmx editconf -f complex.gro -o complex_box.gro -c -d 1.2 -bt dodecahedron

# 4. 溶剂化（TIP3P 水）
gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top

# 5. 加离子（0.15 M NaCl，中和体系净电荷）
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o complex_ions.gro -p topol.top -pname NA -nname CL \
    -neutral -conc 0.15
```

#### 3.2 能量最小化与平衡
```bash
# 能量最小化（最速下降法，Fmax < 1000 kJ/mol/nm）
gmx grompp -f em.mdp -c complex_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# NVT 平衡（310 K，100 ps，位置限制蛋白骨架）
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# NPT 平衡（310 K，1 bar，100 ps，位置限制）
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
```

#### 3.3 生产模拟
```bash
# 生产模拟（500 ns–1 μs，无位置限制，2 fs 步长）
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md -nb gpu -bonded gpu -pme gpu
```

**模拟参数建议**：
- **温度**：310 K（生理温度，使用 V-rescale 耦合器）
- **压力**：1 bar（Parrinello-Rahman 耦合器）
- **长程静电**：PME（ Particle Mesh Ewald）
- **截断**：范德华 1.0 nm，库仑 1.0 nm
- **约束**：LINCS（所有键）
- **输出频率**：每 10 ps 保存一帧坐标

---

### Step 4: 结合模式深度分析（揭示"完美契合"的证据）

这是复现的核心——用定量数据说话，证明 L-Ala2 的"脆弱性"源于几何与静电的双重完美匹配。

#### 4.1 结构稳定性与波动分析
```bash
# 4.1.1 体系整体 RMSD（反映复合物是否稳定）
gmx rms -s md.tpr -f md.xtc -o rmsd_complex.xvg
# 选择 "C-alpha" 或 "Backbone"

# 4.1.2 蛋白各残基 RMSF（反映哪些区域柔性大）
gmx rmsf -s md.tpr -f md.xtc -o rmsf_perres.xvg -res

# 4.1.3 多肽 RMSF（重点关注 Ala2 是否被"锁定"）
# 提取多肽链单独计算
gmx make_ndx -f md.tpr  # 创建多肽索引
gmx rmsf -s md.tpr -f md.xtc -n index.ndx -o rmsf_peptide.xvg
```
**预期结果**：GHRH 的 N 端（1-5 位）RMSF 值应显著低于 C 端（20-29 位），表明 N 端被 DPP-IV 牢牢锚定。

#### 4.2 催化几何动态监测（核心证据）
定义关键原子距离/角度的时序变化：

```bash
# 创建索引组（需手动选择原子）
gmx make_ndx -f md.tpr

# 1. Ser630 Oγ -- Ala2 C=O 碳距离（亲核攻击距离）
gmx distance -s md.tpr -f md.xtc -n index.ndx -o dist_Ser630_Og_Ala2_C.xvg \
    -select 'group "Ser630_OG" plus group "Ala2_C"'

# 2. Ser630 Oγ -- Ala2 C=O 氧距离（氧负离子洞距离）
gmx distance -s md.tpr -f md.xtc -n index.ndx -o dist_Ser630_Og_Ala2_O.xvg

# 3. N-端 NH₃⁺ -- Glu205 Oε 距离（盐桥）
gmx distance -s md.tpr -f md.xtc -n index.ndx -o dist_Nterm_Glu205.xvg

# 4. N-端 NH₃⁺ -- Glu206 Oε 距离（盐桥）
gmx distance -s md.tpr -f md.xtc -n index.ndx -o dist_Nterm_Glu206.xvg

# 5. 催化三联体角度：Ser630 Oγ - Ala2 C=O 碳 - Ala2 N（攻击角度，理想 ~90-110°）
gmx angle -s md.tpr -f md.xtc -n index.ndx -o angle_attack.xvg
```

**"完美契合"的量化标准**：
| 几何参数 | 理想范围 | 物理意义 |
|----------|----------|----------|
| d(Ser630 Oγ → Ala2 C=O C) | 2.5–4.0 Å | 亲核攻击准备态 |
| d(Ser630 Oγ → Ala2 C=O O) | 2.5–3.5 Å | 氧负离子洞预组织 |
| d(Tyr1 N → Glu205 Oε) | < 3.5 Å | N 端锚定盐桥 |
| d(Tyr1 N → Glu206 Oε) | < 3.5 Å | N 端锚定盐桥 |
| ∠(Ser630 Oγ–Ala2 C–Ala2 N) | 80–120° | 亲核攻击角度 |

#### 4.3 氢键与盐桥统计
```bash
# 统计 DPP-IV 与 GHRH(1-29) 之间的氢键
# 以及 GHRH 内部二级结构氢键
gmx hbond -s md.tpr -f md.xtc -num hbond_DPP4_GHRH.xvg -hbm hbond_matrix.xpm
```
**重点观察**：
- **Glu205/Glu206** 与 GHRH N 端之间的盐桥存续时间（occupancy）。
- **Ser630** 与 Ala2 羰基氧之间的氢键。
- **Tyr547** 与 Ala2 羰基氧之间的氢键（氧负离子洞）。

#### 4.4 结合自由能计算（MM-PBSA/GBSA）
```bash
# 使用 g_mmpbsa 或 MMPBSA.py 计算结合自由能
# 需要分离受体、配体轨迹
gmx trjconv -s md.tpr -f md.xtc -o receptor.xtc  # 提取 DPP-IV
gmx trjconv -s md.tpr -f md.xtc -o ligand.xtc     # 提取 GHRH(1-29)

# 运行 MM-PBSA（以 Amber 的 MMPBSA.py 为例）
MMPBSA.py -O -i mmpbsa.in -o RESULTS_MMPBSA.dat -sp complex.prmtop \
    -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
    -y md.xtc
```
**预期结果**：GHRH(1-29) 与 DPP-IV 的结合自由能应处于 **-8 至 -12 kcal/mol** 范围（与已知底物如 Diprotin A 相当），表明结合强而稳定。

#### 4.5 静电势表面映射（Electrostatic Potential Mapping）
这是文档强调的"静电势分布完美契合"的直接可视化证据。

**操作流程**：
1. 从 MD 轨迹中提取代表性的稳定构象（聚类中心）。
2. 使用 **APBS**（Adaptive Poisson-Boltzmann Solver）计算 DPP-IV 活性位点的静电势表面。
3. 使用 **PyMOL** 或 **VMD** 将 GHRH(1-29) 的 N 端叠加到静电势图上。

```python
# PyMOL 脚本示例
# 加载 DPP-IV 和 GHRH
cmd.load("DPP4_frame.pdb", "DPP4")
cmd.load("GHRH_frame.pdb", "GHRH")

# 计算 APBS 静电势（需安装 apbs-plugin）
cmd.apbs("DPP4")

# 显示表面并按静电势着色（红=负，蓝=正）
cmd.show("surface", "DPP4")
cmd.ramp_new("esp", "DPP4", [-5, 0, 5], "blue_white_red")
cmd.set("surface_color", "esp", "DPP4")

# 显示 GHRH N 端
cmd.show("sticks", "GHRH and resi 1-5")
cmd.color("yellow", "GHRH and resi 2")  # 高亮 Ala2
```

**预期可视化效果**：
- DPP-IV 的 S2 口袋（Glu205/Glu206 区域）呈现强烈的 **负电势（红色）**。
- GHRH 的 N 端游离氨基（-NH₃⁺）呈现 **正电势（蓝色）**。
- 两者形成鲜明的 **正负互补**，如同锁钥契合。
- S1 口袋内部为疏水区域（静电势接近中性，白色/绿色），与 Ala2 的疏水甲基侧链完美匹配。

#### 4.6 溶剂可及表面积（SASA）分析
```bash
# 计算 Ala2 的 SASA，判断其暴露/埋藏程度
gmx sasa -s md.tpr -f md.xtc -n index.ndx -o sasa_Ala2.xvg -surface "Protein" -output "Ala2"
```
**预期结果**：Ala2 的 SASA 应显著低于游离多肽中的 Ala2，表明其侧链和骨架被深埋于 DPP-IV 口袋内。

---

## 5. 对照实验设计（增强说服力）

为了证明"完美契合"是 L-Ala2 特有的，建议设置以下对照：

### 对照 A：D-Ala2 突变体（CJC-1295 的改造后形式）
- 将 GHRH(1-29) 的 Ala2 突变为 **D-Ala**。
- 重复对接 + 短 MD（50-100 ns）。
- **预期差异**：D-Ala2 的骨架手性翻转导致 Ramachandran 角度变化，羰基氧无法与氧负离子洞对齐，Ser630 攻击距离 > 5 Å 或角度严重偏离，整体结合自由能显著下降或复合物不稳定。

### 对照 B：Pro2 突变体
- 将 Ala2 突变为 **Pro**（另一种 DPP-IV 底物偏好残基）。
- 预期结合模式类似，但 Pro 的环状侧链可能引入不同构象限制。

### 对照 C：Gly2 突变体
- 将 Ala2 突变为 **Gly**（无侧链）。
- **预期差异**：S1 口袋内出现空腔，疏水互补丧失，RMSF 增大，结合能下降。

### 对照 D：无 GHRH 的 apo DPP-IV MD
- 模拟 apo DPP-IV 的口袋动态，证明口袋在没有底物时的预组织性。

---

## 6. 软件与计算资源需求

### 6.1 必需软件
| 软件 | 版本建议 | 用途 | 获取方式 |
|------|----------|------|----------|
| **GROMACS** | 2023.x / 2024.x | MD 模拟 | 官网 / conda |
| **PyMOL** | 2.5+ | 可视化、结构预处理 | 已安装 |
| **HADDOCK2.4/2.5** | 最新 | 肽-蛋白对接 | 在线/本地 |
| **AmberTools** | 23/24 | 力场参数、MM-PBSA | 官网（免费） |
| **APBS** | 3.4+ | 静电势计算 | conda/pip |
| **MDAnalysis** | 2.7+ | 轨迹分析（Python） | pip |
| **RDKit** | 2023+ | 小分子/肽处理 | pip/conda |

### 6.2 可选高级工具
| 软件 | 用途 |
|------|------|
| **Rosetta** | FlexPepDock 对接、序列设计 |
| **PLUMED** | 增强采样（元动力学、伞形采样） |
| **g_mmpbsa** | GROMACS 兼容的 MM-PBSA 计算 |
| **cpptraj** | Amber 轨迹高级分析 |

### 6.3 计算资源估算
| 任务 | 体系大小 | 建议时长 | GPU 需求 |
|------|----------|----------|----------|
| 能量最小化 + 平衡 | ~70,000 原子 | 1-2 小时 | 可选 |
| 生产 MD（500 ns） | ~70,000 原子 | 2-3 天 | **必需**（1x RTX 4090 / A100） |
| 对照模拟（4 x 100 ns） | ~70,000 原子 | 3-4 天 | 必需 |
| 对接（HADDOCK） | — | 数小时 | 无需 |
| 分析 | — | 半天 | 无需 |

---

## 7. 预期成果与交付物

### 7.1 定量数据表格
| 指标 | GHRH(1-29) WT | GHRH-D-Ala2 (对照A) | 差异说明 |
|------|---------------|---------------------|----------|
| HADDOCK Score | — | — | 结合亲和力打分 |
| 结合自由能 ΔG_bind | — | — | MM-PBSA 结果 |
| d(Ser630 Oγ → Ala2 C) 均值 | — | — | 亲核攻击距离 |
| d(N-term → Glu205/Glu206) 均值 | — | — | 盐桥距离 |
| Ala2 RMSF | — | — | 柔性/锁定程度 |
| Ala2 SASA | — | — | 埋藏程度 |
| H-bond occupancy (Ser630→Ala2) | — | — | 氢键存续率 |

### 7.2 关键可视化图形
1. **图 1**：DPP-IV 活性口袋表面静电势图，叠加 GHRH(1-29) N 端，展示正负互补。
2. **图 2**：MD 过程中关键距离（Ser630-Ala2, N-term-Glu205）的时间序列图。
3. **图 3**：催化三联体与 Ala2 的 2D/3D 相互作用示意图（LigPlot+/PyMOL）。
4. **图 4**：WT vs D-Ala2 的叠加结构对比，展示 D-Ala2 如何破坏几何对齐。
5. **图 5**：RMSF 和 RMSD 曲线，证明复合物稳定性。

### 7.3 结论性陈述（复现成功的标准）
> "在 500 ns 的全原子分子动力学模拟中，天然 GHRH(1-29) 的 N 端（Tyr1-Ala2-Asp3）稳定锚定于 DPP-IV 的 S2/S1 口袋。第 2 位 L-Ala 的甲基侧链深埋于由 Tyr631/Val656/Trp659/Tyr662 构成的疏水 S1 口袋（SASA 降低 XX%），其骨架羰基氧与氧负离子洞残基 Tyr547/Ser631 保持稳定的氢键作用（occupency > 80%）。催化残基 Ser630 的 Oγ 原子与 Ala2 的羰基碳平均距离为 X.XX Å，攻击角度为 XXX°，均处于丝氨酸蛋白酶形成四面体过渡态的理想几何范围内。静电势表面分析显示，带正电的 GHRH N 端与 Glu205/Glu206 的负电势口袋形成完美的静电互补。相比之下，D-Ala2 突变体中，由于手性翻转导致骨架偏移，Ser630 攻击距离增加至 X.XX Å，盐桥断裂，结合自由能下降 XX kcal/mol，无法形成有效的催化准备态。上述结果在原子层面定量复现了 L-Ala2 作为 DPP-IV 酶切最脆弱环节的'完美契合'机制。"

---

## 8. 风险评估与备选方案

| 风险 | 可能性 | 应对措施 |
|------|--------|----------|
| GHRH(1-29) 在 MD 中从 DPP-IV 口袋脱落 | 中 | 使用位置限制（restraint）在 N 端前 5 个残基，或采用伞形采样增强口袋内采样 |
| 长肽柔性过高，对接收敛性差 | 中 | 分段对接：先对接 N 端 5 肽，再逐步延伸 |
| 计算资源不足，无法运行 500 ns MD | 高 | 缩短至 100-200 ns，或使用 Martini 粗粒化力场加速采样 |
| 静电势计算复杂 | 低 | 使用 PyMOL 内置的 Vacuum Electrostatics 做定性展示 |
| D-Ala 力场参数缺失 | 低 | D-Ala 是标准残基，Amber/CHARMM 力场已覆盖 |

---

## 9. 下一步行动计划

1. **环境搭建**：安装 GROMACS、AmberTools、APBS、HADDOCK（或准备在线提交）。
2. **结构准备**：下载 1NU8，预处理 DPP-IV；构建 GHRH(1-29) 结构。
3. **对接测试**：先用 HADDOCK 在线版提交一个快速测试（it0=1000, it1=50）。
4. **小规模 MD 验证**：对最优对接模型运行 10 ns 测试 MD，检查体系稳定性。
5. **生产模拟**：正式运行 500 ns MD + 对照实验。
6. **数据分析与可视化**：提取数据，制作图表，撰写复现报告。

---

*方案版本：v1.0*
*撰写日期：2026-05-13*
*基于文档：CJC-1295 Computational Design Reproduction.md*
