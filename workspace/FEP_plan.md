# P1-FEP 方案：L-Ala2 ↔ D-Ala2 手性翻转的结合自由能差

**目标**：计算 ΔΔG_bind = ΔG_bind(D-Ala2) - ΔG_bind(L-Ala2)
**体系**：DPP-IV + GHRH(1-10) 短肽（~100k 原子）
**预期信号**：文献报道 CJC-1295 半衰期延长 ~1000×，对应 ΔΔG ~4 kcal/mol
**工具**：GROMACS 2026.0 + Amber14sb

---

## 一、理论框架

### 热力学循环

```
                    DPP-IV + GHRH(L-Ala2)  ──ΔG_bind(L)──►  Complex(L-Ala2)
                          │                                    │
                     ΔG_mut(free)                        ΔG_mut(bound)
                          │                                    │
                          ▼                                    ▼
                    DPP-IV + GHRH(D-Ala2)  ──ΔG_bind(D)──►  Complex(D-Ala2)
```

$$
\Delta\Delta G_{bind} = \Delta G_{mut}(bound) - \Delta G_{mut}(free)
$$

只需要计算 **两个** FEP：
1. `ΔG_mut(bound)`：结合态中 Ala2 从 L→D 的自由能变化
2. `ΔG_mut(free)`：自由肽中 Ala2 从 L→D 的自由能变化

**优势**：受体（DPP-IV）不需要参与突变，简化计算。

---

## 二、三种实现路径

### 路径 A：Dual-Topology Alchemical FEP（标准方法）

**原理**：在 hybrid topology 中，Ala2 同时包含 L-构型和 D-构型的原子。
- 共享原子（N, CA, C, O, H, HA）：两个状态相同
- L-specific 原子（CB_L, HB1_L, HB2_L, HB3_L）：state A 完全相互作用，state B 为 dummy
- D-specific 原子（CB_D, HB1_D, HB2_D, HB3_D）：state A 为 dummy，state B 完全相互作用

**实现步骤**：
1. 用 `pdb2gmx` 分别生成 L-Ala2 和 D-Ala2 的拓扑
2. 手动合并拓扑，为 CB/HB 原子定义 A/B 状态参数
3. 创建 GROMACS FEP mdp（soft-core vdw, 10-15 个 lambda 窗口）
4. 每个 lambda 点平衡 5-10 ns，收集 dhdl
5. 用 `gmx bar` 做 BAR 分析

**难点**：
- 手动编辑拓扑极其繁琐（需要精确处理 bond/angle/dihedral 的 A/B 状态）
- dummy 原子与真实原子在空间上接近（同一个残基内），需要 careful soft-core 参数
- pmx 不支持 L↔D 手性翻转（pmx 只做常规氨基酸突变，如 Ala→Gly）

**工作量**：~1-2 周 GPU
**精度**：最高（如果收敛）

---

### 路径 B：Improper Dihedral λ-Driven 翻转（创新方法）

**原理**：利用 GROMACS 的 `bonded-lambdas` 和 `dihedral_restraints`。

Amber 力场通过 improper dihedral 维持手性。如果我们为 improper 定义 A/B 状态：
- State A (λ=0)：improper 相位 = L-手性值（如 -35° for N-CA-CB-C），k = 100 kJ/mol/rad²
- State B (λ=1)：improper 相位 = D-手性值（如 +35°），k = 100 kJ/mol/rad²

在 λ=0.5 时，两个 improper 同时存在但各半强度，系统自然寻找最低能路径（平面过渡态 → 翻转）。

**实现步骤**：
1. 用 `pdb2gmx` 生成 L-Ala2 拓扑
2. 在拓扑 `[ dihedrals ]` 部分，为 N-CA-CB-C improper 添加 B 状态参数
3. mdp：`bonded-lambdas` 控制 improper 切换，`sc-alpha` 启用 soft-core
4. 平衡 + 生产运行，收集 dhdl

**优势**：
- 不需要 dummy 原子，拓扑更简洁
- 物理上更接近真实的手性翻转过程

**难点**：
- 中间态（λ≈0.5）系统可能高度不稳定（平面构型能量很高）
- 需要大量 lambda 窗口（15-20）来确保收敛
- GROMACS 对 improper A/B 状态的支持需要验证

**工作量**：~1 周 GPU
**精度**：高（如果中间态稳定）

---

### 路径 C：Umbrella Sampling / PMF（替代方法）

**原理**：沿 χ1 二面角（N-CA-CB-C）做 umbrella sampling，计算 L-构型（χ1≈-60°）和 D-构型（χ1≈+60°）的相对势阱深度。

**实现步骤**：
1. 用 `gmx wham` 或 PLUMED 沿 χ1 做 US（10-15 个窗口，k=1000 kJ/mol/rad²）
2. 每个窗口平衡 5-10 ns
3. 计算 PMF，提取两个势阱的 ΔG

**优势**：
- 实现最简单，不需要修改拓扑
- 物理直观，中间态稳定可控
- 可以直接看到手性翻转的能垒高度

**难点**：
- 不是严格的 alchemical FEP（用户明确要求炼金术）
- 给出的是构象自由能差，不是结合自由能差
- 需要做结合态和自由态两套 US，工作量不小

**工作量**：~5-7 天 GPU
**精度**：中（取决于 US 窗口密度和采样）

---

## 三、推荐方案

**首选：路径 A（Dual-Topology FEP）**

理由：
1. 最标准的 alchemical 方法，审稿人认可
2. 与文献中的 amino acid mutation FEP 完全一致
3. `gmx bar` 后处理成熟

**但路径 A 的核心障碍是手动编辑拓扑。** 为解决此问题，推荐 **分步实现**：

### Phase 1（准备，~1 天）：
1. 用 Python 脚本自动生成 hybrid topology
2. 验证 topology（`gmx grompp` 通过，能量合理）
3. 单 lambda 点测试运行

### Phase 2（计算，~5-7 天 GPU）：
1. 结合态：15 个 lambda 窗口 × 10 ns = 150 ns 总采样
2. 自由肽态：15 个 lambda 窗口 × 10 ns = 150 ns 总采样
3. 用 `gmx bar` 分析

---

## 四、Hybrid Topology 自动生成思路

### 原子映射策略

Ala2 的原子组成（L-和 D-相同）：

| 原子 | 角色 | 处理方式 |
|------|------|----------|
| N, CA, C, O | 骨架 | 共享，A/B 相同 |
| H, HA | 骨架 H | 共享，A/B 相同 |
| CB | 手性中心 | **两套**：CB_L (state A), CB_D (state B) |
| HB1, HB2, HB3 | 甲基 H | **两套**：HB*_L, HB*_D |

拓扑中需要：
1. 添加 CB_D 和 HB*_D 为 dummy atoms（mass 保留，charge=0, type=DUM）
2. 为 CB_L 和 HB*_L 定义 B 状态为 DUM
3. 为 CB_D 和 HB*_D 定义 A 状态为 DUM
4. 重新分配 bond/angle/dihedral，确保 L-原子只在 L-原子间成键，D-原子只在 D-原子间成键

### 简化的拓扑生成脚本

写一个 Python 脚本：
1. 读取 L-Ala2 拓扑（`topol_L.top`）
2. 读取 D-Ala2 拓扑（`topol_D.top`）
3. 在 `[ atoms ]` 中为 CB/HB 添加 B 状态参数
4. 在 `[ bonds ]`、`[ angles ]`、`[ dihedrals ]` 中确保 A/B 状态正确
5. 输出 `topol_hybrid.top`

---

## 五、风险与对策

| 风险 | 对策 |
|------|------|
| Hybrid topology 生成错误 | 用 `gmx grompp` 验证；对比 L-和 D-单独拓扑的能量 |
| 中间 lambda 不稳定 | 启用 soft-core（sc-alpha=0.5）；减小 lambda 步长 |
| 采样不足 | 每个 lambda 10-20 ns；用 `alchemlyb` 检查自相关时间 |
| ΔG 不收敛 | 增加 lambda 窗口；检查 work 分布重叠 |
| 自由肽构象熵差异 | 自由肽用 restrain 维持 α-helix 或 random coil 采样 |

---

## 六、时间线估算

| 阶段 | 任务 | 时间 |
|------|------|------|
| Day 1 | Hybrid topology 脚本 + 验证 | 1 天 |
| Day 2-3 | 结合态 FEP（15 λ × 10 ns） | 2-3 天 GPU |
| Day 4-5 | 自由肽 FEP（15 λ × 10 ns） | 2-3 天 GPU |
| Day 6 | gmx bar 分析 + 误差估计 | 0.5 天 |
| Day 7 | 报告撰写 | 0.5 天 |

**总计：~1 周**

---

## 七、立即可以开始的工作

1. **写 hybrid topology 生成脚本**（Python，基于 `pdb2gmx` 输出）
2. **用短肽 GHRH(1-10) 体系测试**（已在 GPU 3 跑 MD，拓扑已准备好）
3. **单 lambda 点预测试**（验证 topology 不崩溃）

---

*文档维护者：Claude Code*
*日期：2026-05-14*
