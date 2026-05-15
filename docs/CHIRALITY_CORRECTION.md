# 手性标记修正说明 (Chirality Label Correction)

## 发现日期
2026-05-15

## 问题描述
项目中的所有 WT/D-Ala2 标记与实际的 Ala2 手性完全相反。

## 根因分析
1. 原始模板 `GHRH_1-29_from_7CZ5.pdb` 中的 Ala2 本身就是 **D-构型** (N-CA-CB-C χ₁ ≈ −55°)
2. `build_DAla_mutant.py` 从这个已经是 D-Ala 的结构出发，执行"镜像翻转"
3. 翻转后得到的文件被标记为 "D-Ala2"，但实际是 **L-Ala** (χ₁ ≈ +57°)
4. 未翻转的原始结构被标记为 "WT"，但实际是 **D-Ala2**

## 标记对照表

| 文件名/标签 | 原标记 | 实际手性 | 修正后标签 |
|------------|--------|----------|-----------|
| `GHRH_1-29.pdb` | WT / L-Ala | **D-Ala** | D-Ala2 模板 |
| `GHRH_1-29_DAla2.pdb` | D-Ala2 | **L-Ala** | L-Ala 模板 |
| `GHRH_DPP4_docked_best.pdb` | WT 对接 | **D-Ala** | D-Ala2 对接 |
| `GHRH_DPP4_docked_best_DAla2.pdb` | D-Ala2 对接 | **L-Ala** | L-Ala 对接 |
| `npt.gro` / `md.tpr` | WT 生产 | **D-Ala** | **CJC-1295 (D-Ala2)** |
| `DAla2_production.tpr` | D-Ala2 生产 | **L-Ala** | **天然 GHRH (L-Ala)** |
| `short_peptide_npt.gro` | WT 短肽 | **D-Ala** | **D-Ala2 短肽** |

## 模拟结论修正

### 全长肽 MD (~108 ns / ~36 ns)
- **CJC-1295 (D-Ala2)**：Ser630-Ala2O = 3.15 ± 0.52 Å，93.8% < 4 Å → 抑制剂维持良好结合姿态
- **天然 GHRH (L-Ala)**：Ser630-Ala2O = 6.31 ± 1.48 Å，8.3% < 4 Å → 底物在 MD 中偏离催化位点

### N-末端盐桥
- CJC-1295 (D-Ala2)：Glu205/Glu206-Tyr1 盐桥断裂（pivot motion）
- D-Ala2 短肽：Glu205-Tyr1 盐桥极好（92.4%）

### FEP (L-Ala ↔ D-Ala2)
- 实际计算方向：**D-Ala2 → L-Ala**（因为 λ=0 对应 D-Ala2，λ=1 对应 L-Ala）
- 报告值需取反：ΔG(L-Ala → D-Ala2) = −ΔG_FEP
- 热力学循环本身不受标记影响，只是方向相反

## 代码修正

### 分析脚本
所有 `dist_wt_*` / `dist_dala2_*` / `rmsd_wt_*` / `rmsd_dala2_*` 等输出文件：
- 原 `wt` → 实际为 `cjc1295` (D-Ala2)
- 原 `dala2` → 实际为 `native` (L-Ala)

### build_DAla_mutant.py
该脚本功能正确（镜像翻转），但输入应为 L-Ala、输出为 D-Ala2。
当前使用场景：输入已是 D-Ala，输出为 L-Ala。
**使用方式**：如需要生成真正的 D-Ala2，应从 L-Ala 结构出发。

## 后续工作建议

1. **不再重跑已完成模拟**（成本过高，108+36+112 ns）
2. **所有结论按修正后标签解释**
3. **FEP 结果取反即为正确 ΔG(L→D)**
4. **短肽分析**：现有短肽 = D-Ala2 短肽；如需 L-Ala 短肽对照，需重新生成

## 验证命令

```bash
# 检查任意 PDB/GRO 中 Ala2 的手性
python3 -c "
import math, numpy as np

def dihedral(p1, p2, p3, p4):
    b1 = np.array(p2) - np.array(p1)
    b2 = np.array(p3) - np.array(p2)
    b3 = np.array(p4) - np.array(p3)
    b2n = b2 / np.linalg.norm(b2)
    v = b1 - np.dot(b1, b2n) * b2n
    w = b3 - np.dot(b3, b2n) * b2n
    v = v / np.linalg.norm(v)
    w = w / np.linalg.norm(w)
    x = np.dot(v, w)
    y = np.dot(np.cross(b2n, v), w)
    return math.degrees(math.atan2(y, x))

# For PDB (chain B):
# N = [x,y,z] from 'ATOM ... N   ALA B   2 ...'
# CA = [x,y,z] from 'ATOM ... CA  ALA B   2 ...'
# CB = [x,y,z] from 'ATOM ... CB  ALA B   2 ...'
# C = [x,y,z] from 'ATOM ... C   ALA B   2 ...'
# chi1 = dihedral(N, CA, CB, C)
# chi1 > 0 → L-Ala, chi1 < 0 → D-Ala
"
```
