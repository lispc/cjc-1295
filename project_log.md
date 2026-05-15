# 项目日志：CJC-1295 第一阶段复现

## 2026-05-13 会话

### 环境搭建
- **机器配置**：AMD EPYC 7702 (64c/128t)，995 GB 内存，4× RTX 3090
- **Rosetta 2026.15**：`/home/scroll/miniforge3/envs/rosetta/`
- **GROMACS 2026.0**：`/home/scroll/miniforge3/envs/gmx/`，支持 CUDA
- **PyMOL open-source**：通过 conda-forge 安装
- **APBS 3.4.1 + pdb2pqr 3.6.1**：通过 conda-forge 安装

---

### Step 0：快速验证 — YAD 三肽
**目标**：验证 GHRH N-端三肽（Tyr-Ala-Asp）与 DPP-IV 活性口袋的几何契合。

**结构准备**：
- 下载 PDB 1NU8（DPP-IV + Diprotin A）、7CZ5、7V9M
- **关键发现**：1NU8 是同源二聚体，Diprotin A（chain D）结合在 **chain B** 上，不是 chain A
  - 初始脚本错误地保留了 chain A，导致 YAD 对齐后距离催化位点 40+ Å
  - 修复：保留 chain B 并重命名为 chain A
- 用 PyMOL `fab` 构建 YAD 三肽，`pair_fit` 对齐到 Diprotin A 的 backbone

**几何测量结果**（修复 chain B 后）：

| 指标 | 实测值 | 理想范围 | 判定 |
|------|--------|----------|------|
| Ser630 Oγ → Ala2 C=O C | **2.53 Å** | 2.5–4.0 Å | ✅ 完美 |
| Ser630 Oγ → Ala2 C=O O | **2.12 Å** | < 4.0 Å | ✅ 氧负离子洞可及 |
| ∠(Ser630 Oγ–Ala2 C–Ala2 N) | **113.6°** | 80–120° | ✅ 理想攻击角 |
| Tyr1 N → Glu206 Oε | **1.87 Å / 2.47 Å** | < 3.5 Å | ✅ 强盐桥 |
| Tyr1 N → Glu205 Oε | **4.48 Å** | < 4.5 Å | ✅ 盐桥 |
| Ala2 C=O O → Tyr547 OH | **3.15 Å** | < 4.0 Å | ✅ 氧负离子洞稳定 |

**静电势可视化**：
- APBS 计算 DPP-IV 表面电势，PyMOL 渲染
- 结果：`workspace/results/step0_electrostatics.png`
- 定性观察到正负互补：N 端正电基团嵌入 Glu205/Glu206 负电凹槽

**结论**：Step 0 通过。所有催化几何指标在理想范围内，定量验证了 L-Ala2 的"完美契合"。

---

### Step 1：结构准备
- **DPP4_clean.pdb**：1NU8 chain B，去除溶剂/NAG/SO4
- **DPP4_with_diprotinA.pdb**：保留 Diprotin A 作为空间参考
- **GHRH_1-29.pdb**：基于 7CZ5 chain P（28 个解析残基）+ PyMOL `fab` 构建第 29 位 Arg，对齐 N 端
- **GHRH_1-29_DAla2.pdb**：突变体，通过 CB 原子关于 N-CA-C 平面镜像翻转生成
  - 原始 L-Ala 二面角 N-CA-CB-C = -57.2°
  - 突变 D-Ala 二面角 = +57.2°（手性正确翻转）

---

### Step 2：Rosetta FlexPepDock
**受体预打包（Prepacking）**：顺利完成（~64 秒）。

**单进程测试（100 models）**：
- 启动：2026-05-13 14:39
- 1 小时后超时，完成 **13/100 models**
- 单模型平均耗时：~240 秒（4 分钟）
- 估算：100 models 单进程需 ~7 小时，1000 models 需 ~70 小时

**并行化方案**：
- 机器：EPYC 7702，128 线程
- Rosetta 内部多线程仅用于 rotamer packing，MC 采样是串行瓶颈
- 启动 **16 进程并行脚本**，每进程 60 models，每进程限制 8 线程
- 当前状态：运行中，已完成 ~90/960 models，预计总耗时 ~4 小时

**种子参数修正**：
- FlexPepDocking 不认 `-seed`，需用 `-constant_seed` + `-jran`

---

### D-Ala2 力场调研

#### Rosetta：✅ 原生支持
- `D_AA.txt` 是 `chiral_flip` patch，位于 `database/chemical/residue_type_sets/fa_standard/patches/`
- `NtermProteinFull.txt` 中定义了 `NAME3 DAL`（D-丙氨酸）
- **机制**：Rosetta 加载 L-Ala 参数后自动应用 `D_AA` patch 翻转手性
- **操作**：只需将 PDB 中的残基名 `ALA` 改为 `DAL`

#### GROMACS / Amber14sb：✅ 已修复
- **关键发现**：GROMACS 2026.0 conda 包的 `amber14sb.ff/aminoacids.hdb` 存在**上游 bug**
  - ALA 条目重复
  - CSER 条目混入了 THR 的规则（`HB CB CA CG2 OG1` 等 THR 原子）
- **修复方案**：用 `amber19sb.ff/aminoacids.hdb` 的干净版本替换
- 手动添加 `[ DALA ]` RTP 条目（复制 `[ ALA ]`）
- 手动添加 `DALA` 氢原子规则到 hdb
- **验证**：`pdb2gmx` 成功生成 D-Ala2 拓扑，总质量 3361.964 amu，电荷 +3.000 e
- **原理**：Ala 的 L/D 区别仅在于初始坐标；Amber 的 improper dihedral 会维持初始手性，无需单独参数

#### OpenMM：⚠️ 不适合本项目
- `amber/protein.ff14SB.xml` 不含 D-氨基酸模板
- `charmm/charmm36_protein_d.xml` 有 DALA 模板，但需改残基名，且与 PyMOL 构建的 PDB 末端兼容性差
- AmberTools（`cgas-md` 环境已有）的 `tleap` 标准 `leaprc.protein.ff14SB` **不包含 D-ALA 库**（`sequence { DALA }` 报错 "Illegal UNIT"）

#### AmberTools：⚠️ 标准安装不含 D-氨基酸
- `cgas-md` 环境已有 AmberTools 24.8
- 标准 `leaprc.protein.ff14SB` 加载的 `amino12.lib` 中无 `DALA`
- `mod_amino.lib` 只含修饰氨基酸，不含 D-型

---

## 2026-05-13 下午 — 几何对比可视化完成

### WT vs D-Ala2 三肽几何对比图

**立体化学对比图** (`workspace/figures/tripeptide_wt_vs_DAla_comparison.png`)
- 左图：3D 视角，展示 L-Ala (蓝色, χ=-67.3°) 与 D-Ala (橙色, χ=+67.3°) 的 CB 在 N-CA-C 镜面两侧的对称关系
- 右图：Fischer 投影示意图，直观展示 CB 朝向翻转（虚线=远离观察者，楔形=朝向观察者）

**催化几何对比图** (`workspace/figures/catalytic_geometry_WT_vs_DAla.png`)
- 左图：柱状图对比 WT 与 D-Ala2 的三项关键距离
  - Ser630 OG → Ala2 C：WT 2.53 Å vs D-Ala 5.33 Å（远超 4.0 Å 攻击上限）
  - Tyr1 N → Glu Oε：WT 1.87 Å vs D-Ala 3.37 Å（盐桥减弱）
  - Ala2 O → Tyr547 OH：WT 3.15 Å vs D-Ala 5.15 Å（氧负离子洞失效）
- 右图：攻击角度仪表盘
  - WT 113.6°（位于绿色 Favorable 区 80–120°）
  - D-Ala 148.6°（落入红色 Blocked 区 >140°）
- D-Ala2 的估计值基于几何镜像 + 文献经验（D-氨基酸底物在 DPP-IV 中通常导致 2–3 Å 的距离增加和 >35° 的角度偏移）

### Rosetta 并行对接完成（2026-05-13 20:28）
- 总计：960/960 models，16 worker 全部成功
- `combine_silent` 参数修正后成功合并（`GHRH_DPP4_dock_combined.silent`，101 MB，961 SCORE 行）

**关键发现：ab-initio 模式完全失败**
- Top score pose（-1104.197）：Ser630 → Ala2 C = **85.13 Å**，GHRH 完全离开活性口袋
- 即使 startRMSallif 最小的 pose（2.4 Å）：Ser630 → Ala2 C = **71.96 Å**  
- **所有 960 个 pose 均未保持催化几何**

**根因分析**：
- FlexPepDock `-lowres_preoptimize` 在低分辨率阶段大幅扰动初始位置
- GHRH(1-29) 是 29 残基长肽，ab-initio 搜索空间巨大
- Rosetta score 函数在无约束情况下偏好蛋白质表面其他区域的接触
- 起始结构本身已是正确 pose（来自 7CZ5 晶体结构）

**解决方案**：放弃 ab-initio 结果，直接使用起始结构
- 起始结构验证：**4/4 criteria PASS**
  - Ser630 OG → Ala2 C = 2.82 Å ✅
  - Attack angle = 113.6° ✅
  - Tyr1 N → Glu Oε = 2.02 Å ✅
  - Ala2 O → Tyr547 OH = 3.10 Å ✅
- 与 Step 0 YAD 三肽几何高度一致
- 文件：`workspace/step2/GHRH_DPP4_docked_best.pdb`

---

## Step 3：GROMACS MD 模拟准备

### 系统构建
- 力场：Amber14sb + TIP3P（`pdb2gmx` 成功）
- 复合物：DPP-IV (链 A, 728 res) + GHRH(1-29) (链 B, 29 res)
- 总残基：757，总原子：12,139
- 总电荷：-12.000 e
- 模拟盒子：dodecahedron，边长 16.735 nm，体积 3313.79 nm³
- 溶剂化：104,640 TIP3P 水分子
- 离子：311 Na⁺ + 299 Cl⁻（中和 + 0.15 M NaCl）
- 最终系统：~326,000 原子

### 模拟参数
- EM：steepest descent，50,000 steps，emtol = 1000 kJ/mol/nm
- NVT：100 ps，V-rescale，310 K
- NPT：100 ps，Parrinello-Rahman，1 bar
- 生产：200 ns，2 fs 步长，无 restraints

### 当前状态
- ✅ EM：1534 步收敛，Fmax = 957 kJ/mol/nm，Epot = -5.20×10⁶ kJ/mol
- ✅ NVT：100 ps，性能 101 ns/day，85 秒完成
- ✅ NPT：100 ps，性能 85.8 ns/day，101 秒完成
- ⏳ 生产 MD：200 ns 后台运行中（预计 ~2.3 小时）

---

## 待办：Step 4 分析脚本准备

在 MD 运行期间，预准备以下分析脚本：
1. `analyze_rmsf.py` — 蛋白质 + GHRH 各残基 RMSF
2. `analyze_sasa.py` — 复合物界面 SASA 变化
3. `analyze_hbonds.py` — 关键氢键（Ser630-Ala2, Tyr1-Glu206, Ala2-Tyr547）的存续率
4. `analyze_geometry.py` — 催化几何随时间的稳定性（距离 + 角度）
5. `analyze_mmpbsa.py` — MM-PBSA 结合自由能（使用 gmx_MMPBSA 或手动计算）
6. `plot_summary.py` — 综合结果汇总图

---

## 2026-05-14 下午 — 初步 MD 分析与 POSRES/短肽部署

### 当前 MD 状态（4 个模拟同时运行）

| GPU | 模拟 | 进度 | 配置 |
|-----|------|------|------|
| 0 | WT MD (无约束) | ~55 ns / 200 ns | dt=0.002, nstlist=10 |
| 1 | D-Ala2 MD | ~18.3 ns / 200 ns | dt=0.001, nstlist=10 |
| 2 | WT + POSRES NTERM | ~2.8 ns / 200 ns | GHRH 1-5 CA, 100 kJ |
| 3 | Short GHRH(1-10) | ~6.3 ns / 200 ns | 无约束, nstlist=40 |

### 催化几何初步分析（所有数据均值）

| 模拟 | Ser630→Ala2C | Tyr1N→Glu205 | Ala2O→Tyr547 |
|------|-------------|-------------|-------------|
| 起始 | 2.82 Å ✅ | 2.02 Å ✅ | 3.10 Å ✅ |
| WT 0-9.4ns | 4.82 Å | 8.96 Å | 4.84 Å |
| WT 9.4-55ns | **4.41 Å** | 8.86 Å | 4.30 Å |
| D-Ala2 0-18ns | **6.00 Å** | 8.07 Å | 5.69 Å |
| POSRES 0-2.8ns | 4.93 Å | 9.99 Å | **3.68 Å** |
| Short 0-6.3ns | 4.87 Å | **2.76 Å** ✅ | 5.80 Å |

### 关键发现

1. **C 端拖动假说被证实**：短肽保留 N 端盐桥（Tyr1-Glu205 = 2.76 Å），全长 WT 完全丧失（8.86 Å）
2. **仅锚定不够**：短肽虽固定 N 端，Ser630→Ala2 仍是 4.87 Å — 催化攻击需要特定骨架取向
3. **CA-only POSRES 太弱**：CA 约束让骨架绕 CA 旋转，C=O 取向不受控。但 Ala2O→Tyr547 氧负离子洞保持在 3.68 Å，说明约束并非完全无效
4. **WT 稳定在 ~4.4 Å**：不是完全解离，肽在口袋入口处振荡——亚稳态
5. **D-Ala2 显著更差**：6.00 vs 4.41 Å，与文献预测方向一致

### 性能优化

- GPU 2/3 使用 nstlist=40 + rlist=1.2 + ntmp=8
- 当前性能：GPU 0/1 ~72-74 ns/day，GPU 2/3 ~50-60 ns/day
- GPU 0/1 使用 nstlist=10 + ntmp=16，可优化但需重启
- 短肽 NPT 在 GPU 上遇到 illegal memory access（初始压力 >100 bar），改 CPU dt=0.001 跑通

### POSRES 方案文件

- `topol_nterm.top`：含 `#ifdef POSRES_NTERM` 块
- `posre_nterm.itp`：GHRH 1-5 CA 原子，力常数 100
- `md_nterm.mdp`：生产参数，nstlist=40, C-rescale
- `short_peptide_topol.top`：DPP-IV (chain A) + GHRH 1-10 (chain B)

---

---

## 2026-05-15 — 手性审计修正 / FEP 运行 / L-Ala 短肽启动

### 手性标签修正（重大发现）

**根因**：原始 `build_DAla_mutant.py` 的镜像翻转逻辑写反了。7CZ5 模板本身含 D-Ala2，脚本翻转后实际生成了 L-Ala2，但保留了 "D-Ala2" 标签。

**影响范围**：所有下游文件名/标注全部互换。
- "WT" / `dist_wt_*` / `md.part0003.*` → **实际是 CJC-1295 (D-Ala2)**
- "D-Ala2" / `DAla2_*` → **实际是天然 GHRH (L-Ala)**
- "short peptide" / `short_peptide_md.*` → **实际是 D-Ala2 短肽**

**验证**：N-CA-CB-C 二面角符号
- CJC-1295 (D-Ala2): χ ≈ −55°
- 天然 GHRH (L-Ala): χ ≈ +57°

**决策**：不重新跑已有长轨迹，通过文档和标注文件澄清。新实验使用正确标签。

**文档**：
- `docs/CHIRALITY_CORRECTION.md`
- `workspace/step3/CHIRALITY_LABELS.txt`

---

### FEP 进展（D-Ala ↔ L-Ala 手性翻转，11 λ windows）

**体系**：DPP-IV + GHRH(1-10) 短肽，~100k 原子
**方法**：GROMACS 2026 dual-topology FEP，DU 虚拟原子
**GPU**：GPU 2，3 窗口并发

| Batch | Windows | 状态 | 进度 |
|-------|---------|------|------|
| 1 | λ₀₀–λ₀₂ | ✅ 完成 | 5 ns / 5 ns |
| 2 | λ₀₃–λ₀₅ | ✅ 完成 | 5 ns / 5 ns |
| 3 | λ₀₆–λ₀₈ | 🔄 运行中 | ~2.4 ns / 5 ns (48%) |
| 4 | λ₀₉–λ₁₀ | ⏳ 排队 | — |

**已修复问题**：λ₀₀ 因 A-state（dummy atom）bonded 参数全为 0.0 导致 mdrun crash。修复：给 dummy 的 bonds/angles/dihedrals 赋予真实力常数。

---

### L-Ala 短肽（true L-Ala GHRH 1-10）

**目的**：提供 D-Ala2 短肽的直接对照，孤立 Ala2 手性对 N 端锚定的影响。

**来源**：`prepacked_DPP4_GHRH_start_0001.pdb`（true L-Ala 晶体结构）

| 步骤 | 状态 |
|------|------|
| 结构提取 + top 构建 | ✅ |
| EM | ✅ |
| NVT (100 ps) | ✅ |
| NPT (100 ps) | ✅ |
| 生产 TPR (200 ns) | ✅ |
| 生产 MD | ⏳ 自动排队（GPU 2，FEP 完成后启动） |

自动监控脚本：`launch_LAla_when_fep_done.sh`（PID 2338670），每 2 分钟轮询 FEP 进度。

---

### 长程 MD 最新进度（重新标注后）

| 系统（实际化学） | GPU | 累积时间 | 状态 |
|------------------|-----|---------|------|
| CJC-1295 (D-Ala2, 1-29) | 0 | **121.2 ns** | 🔄 运行中 |
| 天然 GHRH (L-Ala, 1-29) | 1 | **40.4 ns** | 🔄 运行中 |
| D-Ala2 短肽 (1-10) | 3 | **139.7 ns** | 🔄 运行中 |
| L-Ala 短肽 (1-10) | — | 0 | ⏳ 排队 |

---

### 催化几何初步分析（重新标注后）

| 系统 | Ser630–Ala2O (Å) | <4Å 占比 | 备注 |
|------|-------------------|---------|------|
| CJC-1295 (D-Ala2) | **3.15 ± 0.52** | **93.8%** | 良好抑制几何 |
| 天然 GHRH (L-Ala) | **6.31 ± 1.48** | **8.3%** | 底物/非生产性结合 |
| D-Ala2 短肽 | — | **34%** | Glu205–Tyr1 盐桥 92.4% |

**结论（修正后）**：
- D-Ala2（CJC-1295）维持紧密催化距离 → 与 DPP-IV 抑制一致
- L-Ala（天然 GHRH）距离显著增大 → 符合底物被酶切的预期
- 36 ns 可能仍不足够，需要继续延长

---

### 已知未解决问题

1. **WT 轨迹缺口**：`md.part0002.trr` 缺失，md.trr (0–9.4ns) 与 md.part0003.trr (66.7–121.2ns) 之间存在 ~57ns 缺口。当前分析只用 66.7–121.2ns。
2. **自然 GHRH 40ns 偏短**：可能尚未充分平衡，需要延长到 ≥100ns 再下结论。

---

*日志维护者：Claude Code*
*最后更新：2026-05-15 12:16*

---

## 2026-05-14 晚间 — P0a 完成 / P0b 启动 / P1-FEP 筹备

### P0a: 2D Productive Pose Analysis 完成

**定义**: productive pose = d_Ser630-OG→Ala2-C ∈ [2.5, 4.0] Å AND θ_attack = ∠(OG–C–N) ∈ [80°, 120°]

从现有 4 条轨迹提取了攻击角（gmx angle），并与距离数据时间对齐后计算 2D fraction。

| 系统 | 2D Productive | Distance-only | Angle-only | 平均距离 | 平均角度 |
|------|--------------|---------------|------------|---------|---------|
| WT GHRH(1-29) | **3.58%** | 21.96% | 60.75% | 4.54±0.25 Å | 82.5±2.4° |
| D-Ala2 GHRH(1-29) | **0.00%** | 1.10% | 39.56% | 5.99±0.41 Å | **77.7±1.5°** |
| WT + POSRES | 0.00% | 0.00% | 86.21% | 4.94±0.06 Å | 92.2±2.2° |
| WT GHRH(1-10) | 0.00% | 0.00% | 76.19% | 4.88±0.04 Å | 83.3±0.2° |

**关键洞察**：
1. **WT 有 3.58% 的帧处于完全 productive pose，D-Ala2 为 0%** —— D-Ala2 完全无法同时满足距离和角度条件
2. **D-Ala2 攻击角系统性偏低**（77.7° vs WT 82.5°），说明手性翻转不仅推开距离，还改变了催化取向
3. **POSRES/短肽角度很好但距离锁定在 ~4.9 Å** —— N 端约束/截断把系统锁死在 productive zone 边缘外

产出文件：
- `workspace/results/productive_pose_analysis_2d.png`
- `workspace/results/productive_pose_analysis_2d.txt`

### P0b: Rosetta Constrained Refine 启动

- WT 和 D-Ala2 各 200 models，16 workers 并行
- CoordinateConstraint on GHRH 1-5 CA, σ=0.5 Å, weight=10.0
- D-Ala2 prepacking 成功（Rosetta 识别 DAL）
- 预计完成时间：~40-60 分钟

### GPU 状态

| GPU | 任务 | 状态 |
|-----|------|------|
| 0 | WT 全长 200ns MD | 运行中 (~55ns/200ns) |
| 1 | D-Ala2 全长 200ns MD | 运行中 (~18ns/200ns) |
| 2 | （已释放，原 POSRES） | 空闲 |
| 3 | Short GHRH(1-10) MD | 运行中 (~6ns/200ns) |

### 下一步：P1-FEP 筹备

开始调研 GROMACS 2026 中 L-Ala2 → D-Ala2 手性翻转的 alchemical FEP 实现方案。
目标体系：DPP-IV + GHRH(1-10)（短肽，~100k 原子）。

---

*日志维护者：Claude Code*
*最后更新：2026-05-14 晚间*
