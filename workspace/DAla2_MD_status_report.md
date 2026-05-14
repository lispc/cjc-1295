# D-Ala2 MD Status Report

**Date**: 2026-05-13  
**WT MD Progress**: ~7.1 ns / 200 ns (step 3,543,200)

---

## Key Scientific Findings

### 1. D-Ala2 Handedness Does Not Destroy Catalytic Geometry (Static)

| Metric | Step 0 YAD | WT Start | D-Ala2 Start | Threshold |
|--------|-----------|----------|--------------|-----------|
| Ser630 OG → Ala2 C=O C | 2.53 Å | **2.82 Å** ✅ | **2.82 Å** ✅ | 2.0–4.0 Å |
| Attack angle ∠(OG–C–N) | 113.6° | **113.6°** ✅ | **113.6°** ✅ | 80–120° |
| Tyr1 N → Glu Oε | 1.87 Å | **2.02 Å** ✅ | **2.02 Å** ✅ | <3.5 Å |
| Ala2 O → Tyr547 OH | 3.15 Å | **3.10 Å** ✅ | **3.10 Å** ✅ | 2.5–4.0 Å |
| **Pass/Fail** | — | **4/4** | **4/4** | — |

**Conclusion**: The D-Ala2 mirror mutation (χ1: +126° → −116°) does **not** degrade the static catalytic geometry. Both WT and D-Ala2 starting poses are catalytically competent.

### 2. Rosetta Refine Destroys Catalytic Geometry

Rosetta `refine-only` (`-pep_refine` without `-lowres_preoptimize`) systematically displaces the peptide from the active site:

| Model | I_sc | Ser630→Ala2 C | Pass |
|-------|------|---------------|------|
| WT start | — | 2.82 Å | 4/4 |
| WT model 1 | 292 | 5.83 Å | 2/4 |
| WT model 4 (best) | 161 | 4.61 Å | 2/4 |
| D-Ala2 start | — | 2.82 Å | 4/4 |
| D-Ala2 model 3 (lowest I_sc) | **35** | **9.67 Å** | **0/4** |
| D-Ala2 model 4 | 280 | 3.81 Å | 2/4 |

**Critical insight**: Rosetta interface score (I_sc) is **anticorrelated** with catalytic competence. Lower I_sc = peptide pushed further from DPP-IV surface. The starting pose (unrefined) is the best catalytic configuration.

### 3. D-Ala2 EM Side Effects

- χ1 successfully flipped: +126° → −116°
- Backbone forced toward trans-linear: φ −146° → +174°, ψ +159° → +177°
- Catalytic geometry slightly degraded post-EM: Ser630→Ala2 C 3.47→3.62 Å, oxyanion hole 3.50→3.88 Å

---

## Technical Blocker: GROMACS 2026.0 Segfault

### Symptom
D-Ala2 MD crashes with **segmentation fault** at ~250–280 steps (0.5–0.6 ps) in both NVT and NPT ensembles, on both GPU and CPU.

### Diagnostic History
| Test | Steps | Result | Notes |
|------|-------|--------|-------|
| WT NPT | 200 | ✅ Success | Same mdp, same hardware |
| D-Ala2 NVT (GPU) | 50000 | ❌ Fail | Only step 0 |
| D-Ala2 NVT (CPU) | 100 | ✅ Success | With `-nb cpu` |
| D-Ala2 NVT (CPU) | 200 | ✅ Success | EM coords |
| D-Ala2 NVT (CPU) | 250 | ✅ Success | Threshold edge |
| D-Ala2 NVT (CPU) | 280 | ❌ Fail | Threshold crossed |
| D-Ala2 NPT (CPU) | 500 | ❌ Fail | All pressure coupling algorithms |
| D-Ala2 NPT (GPU) | 50000 | ❌ Fail | Segfault at step ~100 progress |
| D-Ala2 chunk from GRO | 200×3 | ✅ Success | Workaround found |

### Root Cause Analysis (Inconclusive)
- ❌ Not GPU-specific (CPU also fails)
- ❌ Not pressure-coupling-specific (NVT also fails >250 steps)
- ❌ Not coordinate-specific (EM coords fail; chunk restart from GRO succeeds)
- ❌ Not LINCS failure (no warning output)
- ❌ Not box/atom-count issue (gmx check validates TPR)
- ⚠️ Likely GROMACS 2026.0 bug triggered by specific atom configuration
- ⚠️ Suspect: `lincs-iter=1` + D-Ala2 improper dihedral strain → NaN propagation → segfault in force/PME kernel

### Workaround: Chunked Restart from GRO
MD can be sustained by restarting every 200 steps from a fresh `.gro` file:
```bash
gmx trjconv -f chunk.trr -s chunk.tpr -o next.gro -dump 0.4
gmx grompp -f md_nvt.mdp -c next.gro -p topol.top -o next.tpr
gmx mdrun -deffnm next -s next.tpr -nsteps 200
```
**Impracticality**: 200 ns requires ~500,000 chunks → ~60 days at CPU speed.

### OpenMM Alternative
- OpenMM 8.5.1 reads GROMACS top/gro successfully
- CUDA platform fails: `CUDA_ERROR_UNSUPPORTED_PTX_VERSION` (driver 580.95.05 vs CUDA 13.0)
- CPU platform works but EM of 324k atoms takes >3 minutes → production 200 ns would take months

---

## Recommendations

### Immediate (P0)
1. **Continue WT 200 ns MD** — on track, ~7.1 ns completed
2. **Monitor Ser630→Ala2 distance drift** in WT trajectory using `monitor_md_distances.py`

### Short-term (P1)
3. **For D-Ala2 MD**, options ranked by feasibility:
   - **Option A**: Rebuild D-Ala2 with proper D-氨基酸 force field (e.g., CHARMM with D-Ala parameters) and re-run GROMACS
   - **Option B**: Compile OpenMM from source with matching CUDA toolkit to enable GPU acceleration
   - **Option C**: Downgrade GROMACS to 2024.x or 2025.x (may avoid segfault)
   - **Option D**: Accept chunked GROMACS restart and run short D-Ala2 trajectory (e.g., 1–5 ns) for qualitative comparison

### Long-term (P2)
4. **Score-only Rosetta analysis** of WT vs D-Ala2 starting poses (WT completed, D-Ala2 blocked by file path issue)
5. **Clean up step2** bloated silent files (~315 MB)
6. **Publish methodology note** on Rosetta refine-only destroying near-native catalytic poses

---

## Files Updated
- `validate_docked_pose.py` — chain-agnostic, supports ALA/DALA lookup
- `monitor_md_distances.py` — MDAnalysis-based trajectory monitor
- `run_dala2_openmm.py` — OpenMM GROMACS-import script (CUDA incompatible)
