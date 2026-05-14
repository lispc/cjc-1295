# 下一步策略建议：D-Ala2 几何契合度复现

**日期**：2026-05-14
**作者**：Claude (Opus 4.7, 1M context)
**触发问题**：当前无法稳定复现文献中"D-Ala2 破坏 DPP-IV/GHRH 几何契合度"的结论
**配套核对**：`active_task.md`、`project_log.md`、`workspace/DAla2_MD_status_report.md`、`workspace/step3/` 日志

---

## 一、为什么现在复现不出来

### 1.1 实验设计错位

文献的命题是 **"D-Ala2 能不能形成 productive catalytic pose"**（静态/动力学竞争）。
当前 200 ns 无约束 MD 在测的是 **"GHRH(1-29) 会不会从口袋脱落"**（解离）。这是两个不同的问题。

### 1.2 数据本身揭示的核心矛盾

来自 `project_log.md:198-206` 的几何均值表：

| 模拟 | Ser630→Ala2C | Tyr1N→Glu205 | Ala2O→Tyr547 |
|------|-------------|-------------|-------------|
| 起始 | 2.82 Å ✅ | 2.02 Å ✅ | 3.10 Å ✅ |
| WT 9.4-55ns | **4.41 Å** ❌ | 8.86 Å ❌ | 4.30 Å |
| D-Ala2 0-18ns | **6.00 Å** ❌ | 8.07 Å ❌ | 5.69 Å |
| Short GHRH(1-10) | 4.87 Å | **2.76 Å** ✅ | 5.80 Å |

**关键矛盾：WT 自己就已经不在催化态了**。没有 positive control 的"D-Ala2 比 WT 差"等于没结论——审稿人第一个问的就是"那 WT 怎么也飘了"。

### 1.3 三个证据指向同一个 root cause

1. **C 端拖动效应**（项目已验证）
   短肽 GHRH(1-10) 保留 N 端盐桥（2.76 Å），全长 WT 完全失去（8.86 Å）。29 残基里 24 个暴露在溶剂的 C 端尾巴在拉肽段。

2. **起始结构来自非催化复合物**
   7CZ5 是 GHRH/GHRHR 复合物，不是 DPP-IV/GHRH。GHRH 的 "productive bound state on DPP-IV" 在 PDB 里没有直接证据，只是基于 Diprotin A 类比 + 对齐拼出来的。当前的"完美契合"是几何构造的结果，不是观测。

3. **时间尺度错配**
   D-Ala2 vs L-Ala2 的差别可能就 1-2 kcal/mol（半衰期 ~7 天 vs ~7 分钟，因子 ~1000，对应 ΔΔG ~4 kcal/mol，含 enzyme turnover + 结合两部分）。无偏 MD 在 200 ns 内难以可靠区分这个量级，更何况 sampling window 远未覆盖 unbinding 事件。

### 1.4 当前数据的统计噪声水平

- WT vs D-Ala2 均值差 1.59 Å（4.41 vs 6.00）
- **没有报误差棒**。29 残基长肽的 Ser630→Ala2 距离自相关时间可能在 ns 量级，块平均后单 simulation 的标准误差可能 > 1 Å
- D-Ala2 使用 dt=0.001（WT 的一半），相同 wall-clock 下 sampling 时间也减半
- **当前的"D-Ala2 明显更差"声明在统计上站不住**

---

## 二、你真正应该问的问题（三选一）

| 命题 | 适合的方法 | 工作量 | 是否直接回答文献 |
|------|------------|--------|------------------|
| **(a)** "D-Ala2 的 productive pose 比 L-Ala2 不稳定" | PMF/umbrella sampling 沿 Ser630–Ala2C 距离，比较势阱深度 | 中（GPU 几天） | 间接 |
| **(b)** "D-Ala2 的结合自由能比 L-Ala2 差 ΔΔG" | Alchemical FEP（L-Ala ↔ D-Ala 在结合态 + 自由态） | 大（GPU 一周） | **直接** |
| **(c)** "D-Ala2 不形成 productive pose"（纯静态/几何） | Rosetta refine **加约束 + native**，比 I_sc + 催化几何 PASS 率 | 小（1 小时） | 部分 |

**说明**：文献的因果链是

```
催化几何丧失 → DPP-IV 不能切 → 半衰期变长
```

中间这一步是文献假设、不是直接观测。最容易反驳 / 支持的是 (c)（静态结构）和 (b)（自由能）。**(a) 和当前的 200 ns 不约束 MD 都是 indirect** —— 离 wet-lab 的可观测量（half-life、kcat/Km）有两步推理距离。

---

## 三、建议路线图（按性价比排序）

### P0a：把现有数据用对（~30 分钟）

**这是最便宜的步骤，可能直接就够用。**

1. **不要再用"均值距离"作主证据**。改成 **fraction of frames in catalytic pose**：
   ```
   pose_productive = (2.5 ≤ d_Ser630→Ala2C ≤ 4.0) ∧ (80° ≤ ∠ ≤ 120°)
   ```
   这个指标可能 WT 0.1% / D-Ala2 0.001% —— **三个数量级差，比 4.41 vs 6.00 Å 说服力强得多**。
   即使是 WT 0.05% / D-Ala2 0.005%，10× 比例差异在 paper 中就足够支撑命题。

2. **加 block averaging 算 mean ± SE**。1.6 Å 差异不带误差棒等于没差异。
   - 推荐用 `pyblock` 或手写：将轨迹切成 N 块，每块算均值，对块均值求 SE
   - 自相关时间用 `gmx analyze -ac` 估
   - 若块大小 > 5 × 自相关时间，SE 才可靠

3. **算 first passage time**：从起始 productive pose 第一次离开 (d > 4.0 OR θ ∉ [80,120]) 的时间。WT vs D-Ala2 这个 FPT 的对比比均值有意义得多。

4. **算 productive pose 内的驻留时间分布**：直方图 of dwell times。如果 WT 有长尾（偶尔会回到 productive），D-Ala2 没有，这才是 "D-Ala2 破坏几何契合度" 的动力学定义。

实现可以写在 `workspace/scripts/analyze_productive_pose.py`，复用现有的 `dist_*.xvg` 和 `ang_*.xvg`。

### P0b：补齐 Step 2 的 Rosetta 证据（~1 小时，包括运行）

上轮复核反复指出，目前仍未做。

```bash
FlexPepDocking \
    -database $DB \
    -s prepacked_DPP4_GHRH_start_0001.pdb \
    -native prepacked_DPP4_GHRH_start_0001.pdb \   # 关键：让 RMS 有意义
    -pep_refine \
    -constraints:cst_fa_file backbone_coord.cst \  # GHRH 1-5 CA 的 CoordinateConstraint, σ=0.5 Å
    -constraints:cst_fa_weight 10.0 \
    -ex1 -ex2aro -use_input_sc \
    -nstruct 200 \
    -constant_seed -jran $SEED \
    -out:file:silent GHRH_DPP4_refine_constrained.silent
```

WT 和 D-Ala2 各 200 个模型，对比：
- **ΔI_sc 中位数**（WT vs D-Ala2）
- **催化几何 PASS 率**（4/4 in [2.5, 4.0] Å 等）
- **接触模式分布**（CB → S1 口袋距离）

这是**论文里 Step 2 唯一能给出 "Rosetta 维度证据" 的方式**。

### P1：短肽 D-Ala2 对照 MD（~2 天 GPU）

当前缺 **D-Ala2 GHRH(1-10) 的对照**。短肽体系 N 端锚定保持得好（Tyr1-Glu205 = 2.76 Å），是目前唯一接近"真实催化前态"的体系。

```bash
gmx grompp -f md_short.mdp -c short_DAla2_npt.gro -p short_DAla2_topol.top -o short_DAla2_md
gmx mdrun -deffnm short_DAla2_md -nb gpu -pme gpu -bonded gpu -update gpu
```

预期信号比全长 MD 干净得多。如果短肽体系 WT vs D-Ala2 在 productive pose fraction 上出现明显差异，这就是论文的核心证据。

### P1：换问题做 alchemical FEP（~1 周 GPU，但定量答案）

- 体系：DPP-IV + GHRH(1-10) 复合物（C 端短，sampling 够）
- 转化：Ala2 的 χ1 从 L 到 D（手性 transform，dual topology 或 hybrid topology）
- 输出：ΔΔG_bind = ΔG_bind(D-Ala2) - ΔG_bind(L-Ala2)
- 工具：GROMACS BAR / `gmx_lambdadyn` / `alchemlyb` 后处理
- 注意：手性翻转不是常规的 vdW 关闭，需要 CB 的位置 (xyz) 转化路径，比常规 FEP 复杂；可考虑 NEQ-FEP (Crooks/Jarzynski)

**这是唯一一个直接对应"D-Ala2 改变结合"的定量量**，与文献的 kinetic 因子 ~1000（ΔΔG ~4 kcal/mol）可对照。

### P2：放弃全长 GHRH(1-29) 无约束 MD 作为主证据

- 让 GPU 0/1 跑完 200 ns，但**只作为 "我们也试了 unbiased MD，确实 C 端拖动" 的负结果展示**
- 主线证据走 P0a（重分析）+ P0b（Rosetta-with-native）+ P1（短肽 / FEP）
- POSRES (GPU 2) 路径**也建议改方向**：现在 CA-only 100 kJ 不够，但加到全骨架 500 kJ 又过强（成本 / 收益不好）；不如把 GPU 2 也用来跑短肽 D-Ala2

---

## 四、决策矩阵

| 你的需求 | 推荐路径 | 所需时间 |
|----------|----------|----------|
| **就想知道现有数据够不够** | P0a 重分析 | 30 分钟 |
| **写一篇可发表的论文，要 Rosetta 证据** | P0a + P0b | 半天 |
| **写一篇可发表的论文，要 MD 证据** | P0a + P1 (短肽 D-Ala2) | 3-4 天 |
| **要一个定量的 ΔΔG，能直接对照实验** | P0a + P0b + P1 (FEP) | 1-2 周 |
| **写 review / blog post，不要严格 paper-grade 证据** | P0a 重分析就够 | 30 分钟 |

---

## 五、关于"复现文献"这件事的元层面建议

文献中 D-Ala2 抗 DPP-IV 切割的证据**本身**来自：
1. **半衰期测量**（in vitro DPP-IV 降解 assay，最直接，但需 wet-lab）
2. **kcat/Km 测量**（酶动力学）
3. **静态结构论证**（D-amino acid 的 CB 朝 S1 口袋，几何不匹配——但这个论证强度低，只是"motivation"）

**没有任何文献用 unbiased MD "复现" 了这个观察**——因为 MD 的时间尺度根本不够看到酶催化（μs–ms）。所以"用 MD 复现文献"本身就是 ill-posed 的目标。

合理的替代目标（按强度递增）：
1. **复现 D-Ala2 的几何不匹配**（静态，已经做到，Step 0 + Step 2 起始结构）
2. **展示 D-Ala2 的 productive pose 在 MD 中比 WT 不稳定**（动力学比较，P0a + P1 可达）
3. **计算 D-Ala2 的结合自由能差**（热力学，P1 FEP 可达）
4. **计算 D-Ala2 的催化能垒变化**（QM/MM，远超当前 scope）

---

## 六、核心结论

> 当前复现不出来不是因为方法做错了具体步骤，而是因为**实验设计在回答错的问题**。无偏 200 ns MD 测的是"全长 GHRH 会不会脱落"，而文献命题是"D-Ala2 的催化几何是否不利"。两者交集很小。

**立即可做的最便宜动作**：P0a 重分析。30 分钟之内可以拿到 productive pose fraction 的对比。如果差异 ≥ 10×，文献结论被"动力学方式"重现了；如果差异 < 3×，需要补 P0b + P1。

**如果只能选一件事做**：选 P0b（Rosetta refine with `-native` + constraints）。它最便宜、最直接对应 Step 2 的预期产出、最容易写进 paper。

---

*文档维护者：Claude*
*与 `reviews/claude.md`（2026-05-13 / 2026-05-14 复核）配套阅读*
