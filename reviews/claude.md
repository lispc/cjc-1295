# 项目复核 — CJC-1295 第一阶段

**复核日期**：2026-05-13
**复核者**：Claude (Opus 4.7, 1M context)
**复核范围**：`project_log.md`、`AGENTS.md`、`active_task.md`、`workspace/scripts/*`、`workspace/step2/*` 日志与 silent 文件、`docs/Phase1_DPP4_GHRH_Perfect_Fit_Protocol.md`

---

## 一、Rosetta FlexPepDock "失败" 的真正原因（与补救方案）

### 1.1 现象复述

`project_log.md` 中将这次对接判定为 **"ab-initio 模式完全失败"**：960 个 pose 的 top 模型 Ser630 → Ala2 C 距离 = 85 Å，肽段完全离开口袋；最终方案是放弃所有对接结果，直接用起始结构（来自 7CZ5 对齐）作为 MD 起点。

这个**判定本身没错**，但**结论过早**——并不是 FlexPepDock "不行"，而是**协议选择错了**。本次对接的几个方法学错误，使它在原理上就不可能保持催化几何：

### 1.2 根因（协议层面）

| 问题 | 调用了什么 | 应该怎么做 |
|------|-----------|------------|
| **A. 误用 ab-initio 模式** | `-lowres_preoptimize` + `-pep_refine` | 起始结构来自 7CZ5 晶体对齐 + 手动放置，已经是接近-native 的 pose，应该用 **refine-only** 而不是 ab-initio |
| **B. 没有 native 参考** | 缺 `-native ...` | RMS* 列因此是相对于 *起始结构* 自身计算的（log 中 `[ WARNING ] No native supplied` 已点出），用最低 startRMSall 当 "近 native" 没有意义 |
| **C. 没有约束** | 无 `-constraints:cst_fa_file` | 在 lowres 阶段，质心模式下 GHRH 的 N 端会被 score 函数推向蛋白其他疏水接触面（你在日志里观察到的 85 Å 偏移正是这个机制） |
| **D. 没有 fragment 库** | 缺 `-frag3` / `-frag9` | 协议文档（`Phase1_..._Protocol.md` §2.3）里写了，但实际跑没用上。29 残基长肽在 lowres MC 中会失控，fragment 是抑制这种失控的关键 |
| **E. 锚定残基偏移** | 自动选 `Peptide anchor: 743`（即 GHRH 第 15 位） | 这意味着 fold tree 的固定点在肽中段；对于 N 端催化，应让锚点贴近 N 端（残基 1 或 2），否则 lowres MC 会把 N 端"甩"出口袋 |

可看到的证据：worker_00.log 中能看到 `protocols.TrialMover: Acceptance rate: 0.04`、`Line search failed even after resetting Hessian; aborting at iter#10` 反复出现——优化器在 lowres 阶段不收敛是常态，这本身就是上面 D + E 的副作用。

### 1.3 这次对接**不是没用**，至少有这些可保留的产物

- **prepacked 结构** `prepacked_DPP4_GHRH_start_0001.pdb`：受体侧链已 repacked，可继续作为 refine 起点（不要用 raw `DPP4_GHRH_start.pdb`）。
- **起始 pose 的 Rosetta 能量基线**：从 silent 文件可读出未扰动起始构象的 `score` / `I_sc` 等，作为 D-Ala2 比较的参照。

### 1.4 三条递进的补救方案（**强烈建议至少做方案 A**）

#### 方案 A：refine-only 重跑（最推荐，~30–60 分钟单进程，并行 ~5 分钟）

```bash
FlexPepDocking \
    -database $DB \
    -s prepacked_DPP4_GHRH_start_0001.pdb \
    -native prepacked_DPP4_GHRH_start_0001.pdb \   # 关键：让 RMS 有意义
    -pep_refine \                                   # 仅高分辨率精修，不再 lowres 全局搜索
    -ex1 -ex2aro -use_input_sc \
    -nstruct 200 \                                  # 200 个就够，refine 收敛快
    -constant_seed -jran $SEED \
    -out:file:silent GHRH_DPP4_refine.silent
```

去掉 `-lowres_preoptimize` 是核心改动。`-pep_refine` 只做小范围 backbone + 侧链精修，不会把肽段甩走。`-native` 让 RMS 列对得上 7CZ5 真实参考，可以用 RMS < 1 Å 的子集做能量分布。

#### 方案 B：约束式 ab-initio（如果想继续验证 sampling 鲁棒性）

写一个 `cat_geom.cst`，把 4 项关键催化几何作为软约束（HARMONIC 或 BOUNDED）：

```
AtomPair  OG  630A  C   2B   HARMONIC  2.8 0.5
AtomPair  OG  630A  O   2B   HARMONIC  2.5 0.5
AtomPair  N    1B  OE1 205A  HARMONIC  3.0 0.7
AtomPair  N    1B  OE1 206A  HARMONIC  3.0 0.7
AtomPair  O    2B  OH  547A  HARMONIC  3.2 0.5
Angle     OG  630A  C   2B  N  2B  HARMONIC  100  15  # 度，标准差 15°
```

附加：`-constraints:cst_fa_file cat_geom.cst -constraints:cst_fa_weight 5.0 -constraints:cst_file cat_geom.cst -constraints:cst_weight 5.0`（fa 和 cen 都加）。这样即便启用 `-lowres_preoptimize`，肽段也无法漂移。

#### 方案 C：score-only（最快，<1 分钟，专做能量打分）

```bash
FlexPepDocking -s start.pdb -flexpep_score_only -out:file:silent score_only.silent
```

只读起始结构，输出 `I_sc`、`I_bsa`、`I_hb`、`pep_sc` 等。这给你一个"起始 pose 在 Rosetta 评分下到底有多好"的数字，可作为 WT vs D-Ala2 的对比基线。

### 1.5 对当前 200 ns MD 的连带建议

既然 Rosetta 的 sampling 已被证明对长肽不友好，而**目前生产 MD 在无任何位置约束下跑 200 ns**，N 端从口袋滑出的风险是真的存在的（Step 0 的 YAD 三肽是 3 残基，整体被埋；GHRH 1-29 是 29 残基，C 端 24 个残基都暴露在溶剂里，会拖动 N 端）。

建议（按代价从低到高）：

1. **现在的 MD 跑完后先看 Ser630–Ala2 距离时序**。如果 100 ns 内距离稳定在 2.5–4 Å，那原方案没问题。
2. 如果出现脱落，加一个**前 50 ns 的弱 N 端 backbone 限制**重跑（mdp 里 `define = -DPOSRES_NTERM`，仅限制 GHRH 1–5 的 CA），50 ns 后释放，让构象自然平衡。
3. 进一步可考虑 **PLUMED 加 reaction coordinate restraint**（Ser630 OG–Ala2 C 距离），保持在 3 ± 0.5 Å。

`docs/Phase1_..._Protocol.md` §8 "风险评估" 第一项就是这个，但实际跑的时候没启用，建议立即评估。

---

## 二、其他做得不好的地方（按严重度排序）

### 🔴 严重

1. **`cc.sh` 在项目根目录，包含明文 OpenRouter API key**
   - 该文件未被 `.gitignore`（git status 显示 `??`，意味着只要 `git add .` 就会泄漏）
   - 立刻：把它移到家目录、加入 `.gitignore`、并**轮换该 key**（已经在磁盘上明文存在 ~7 小时，且我作为 review 工具读到了）
   - 永久：把环境变量移到 `~/.zshrc` / `~/.bashrc` 或 `direnv` 加载的 `.envrc`（`.envrc` 也要 gitignore）

2. **`active_task.md` 与现实严重不一致**
   - 文件里说 Task 1 "🔄 运行中"、Task 3 "等待中需 Task 1 完成"
   - 实际上 `project_log.md` 已记录 Task 1 完成（但失败）、Task 3 已经在跑、Task 4（D-Ala2）也在 EM 阶段
   - 这种文档漂移在多 agent / 多会话场景下很危险，下一个上来的 agent 会按错的状态行动。建议每次更新 log 时同步更新 task 状态，或干脆把两份文件合并

### 🟡 中等

3. **项目根目录散落 GROMACS/MC 中间产物**
   - `posre.itp`（GROMACS pdb2gmx 输出）、`io.mc`（APBS 的 MC shell 日志）、`#test_dala.top.1#`（emacs autosave，空文件）都不该在根目录，应在 `workspace/step3/` 或被清理。
   - `workspace/step3/` 里也有 `#topol.top.1#`、`#topol.top.2#`、`#DAla2_topol.top.1#` 等 emacs autosave —— 提示有人在运行 GROMACS 的同时手编辑了 topology，可能引入未追踪的差异，应核对 topol.top vs autosave 的内容。

4. **磁盘膨胀：step2/ 占 ~315 MB，多数已无用**
   - 16 个 `worker_*.silent`（每个 13 MB） + 16 个 worker log（每个 1.7 MB）+ `GHRH_DPP4_dock_combined.silent`（101 MB）
   - 既然结论是"放弃 ab-initio 结果"，建议保留 1 份 combined silent 用作历史记录，删除 16 个 worker silent；worker log 压缩到 `.tar.gz`。能省 ~250 MB。
   - 但**先做完方案 A 的 refine 重跑再清理**，refine 可能还会用到 prepacked 文件。

5. **D-Ala2 复合物用"CB 镜像"构建后未充分弛豫验证**
   - `build_DAla2_complex.py` 只镜像了 CB / HA / HB1-3 的坐标，没有在镜像后**重新对侧链做 minimization** —— 镜像把 CB 翻到 N-CA-C 平面另一侧后，与 S1 口袋 Tyr631/Trp659/Val656 之间的几何很可能产生轻微 clash（脚本最后做了一个静态 4 Å clash check，但没做 fix）
   - 建议：镜像后再做一次 `gmx mdrun -minimize` 或 PyMOL `cmd.minimize()`，或在 Rosetta 里 FastRelax + 邻居侧链 repack
   - DAla2_em 已经跑了 EM（看时间戳 21:19），可能已经隐式弛豫掉了，但**应该在 EM 前后输出 Ala2 的 Ramachandran (φ,ψ)** 验证 D-氨基酸的角度落在 D-鏡像区（φ≈+60, ψ≈+45），而不是 EM 把它"翻回" L 型——D/L 在共价拓扑里没区别，只在初始坐标 + improper dihedral 维持。这是 amber14sb 路线的隐性风险。

6. **`validate_docked_pose.py` 硬编码了 `("A", 630, "SER", "OG")` 等元组**
   - Rosetta 输出的 silent → PDB 中链 ID 经常变（特别是 `extract_pdbs` 后），脚本会 KeyError
   - 建议改成"先按残基号 + 残基名查找，链 ID 不限"

### 🟢 轻微

7. **`run_flexpepdock_parallel.py` 的 `combine_results()` 命令格式错** —— 已在最近会话被绕过修复（项目日志提到），但脚本本身仍是错的（`combine_silent` 需要 `-in:file:silent` 而非位置参数 + `-out:file:silent` 而非 `-out`）。如果不修，下次跑 D-Ala2 dock 还会再踩一次。

8. **协议文档 § 2.2 fragment 库** 未生成。如果以后要做正式的 ab-initio benchmark，仍需补这一步。

9. **MD 参数未在文档中列出 thermostat/barostat 详细配置**（只有"V-rescale 310K"），MD log 应在 `project_log.md` 里贴一个 `mdp` 摘要，方便复现。

10. **`workspace/scripts/`** 里有 `prepare_structures.py` 和 `prepare_structures_v2.py` 并存。"v2" 命名在生产代码里是技术债（哪个是当前正确的？日志里也没说），建议合并或在文件头注释中明确 deprecated 版本。

---

## 三、结论与建议优先级

| 优先级 | 行动 | 预期收益 | 工作量 |
|--------|------|----------|--------|
| **P0** | 轮换 OpenRouter key + gitignore `cc.sh` | 防泄漏 | 5 分钟 |
| **P0** | 跑方案 A（refine-only FlexPepDock）补 Rosetta 数据 | 让 Step 2 真正 PASS、给 paper 一个 Rosetta I_sc 数字、为 D-Ala2 提供对照 | 1 小时（含跑完） |
| **P1** | 同步更新 `active_task.md` 为当前真实状态 | 避免后续 agent 误操作 | 10 分钟 |
| **P1** | 给生产 MD 加上催化距离时序的实时监控（`gmx distance` 每 1ns 一次） | 提早发现脱落，避免白跑 200 ns | 20 分钟 |
| **P2** | 跑方案 C（score-only）对 WT 起始 pose 和 D-Ala2 起始 pose 都打分 | 给"完美契合"再加一个 Rosetta 维度的定量证据 | 10 分钟 |
| **P2** | 验证 D-Ala2 EM 前后 Ala2 的 (φ,ψ) | 保证 D 手性没在 EM 中被翻回 | 15 分钟 |
| **P2** | 清理 step2/ 中无用 silent + 整理根目录散文件 | 节省 ~250 MB，目录清爽 | 10 分钟 |
| **P3** | 修 `combine_silent` 命令、合并 `prepare_structures*.py`、改硬编码链 ID | 长期可维护性 | 30 分钟 |

**核心结论**：Rosetta 对接**没有从科学上失败**——只是用错了模式。当前"放弃 ab-initio 直接用起始结构"是务实选择，但还差一步：用 **refine-only + native reference** 把这个起始结构正式跑一遍，把 Rosetta 的能量学证据补齐。否则论文里 Step 2 就只有 "起始结构验证 4/4 PASS"，**没有任何 Rosetta sampling 的支撑**，这会让审稿人质疑"为什么用了 Rosetta 还不如不用"。

补救代价低（~1 小时），但能把 Step 2 从"妥协"升级回"完整证据链"。
