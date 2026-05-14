# D-Ala2 MD Status Report

**Date**: 2026-05-14
**WT MD Progress**: ~40 ns / 200 ns (GPU 0, dt=0.002)
**D-Ala2 MD Progress**: ~12.9 ns / 200 ns (GPU 1, dt=0.001)

---

## Key Scientific Findings

### 1. D-Ala2 Handedness: χ1 Flip Confirmed, Backbone Initially Distorted

> **Revision note (2026-05-14)**: Original Finding 1 was a tautology — the 4 catalytic metrics
> measured only backbone atoms (C, N, O, CA), none of which were moved by the CB mirror operation.
> The real finding is below.

**χ1 handeness flip**: +126° (L-Ala2) → −116° (D-Ala2) ✅ confirmed by `build_DAla2_complex.py`

**Backbone distortion (CRITICAL)**: The starting structure suffered from a `fab()` + `pair_fit()` bug
in `prepare_docking_start.py` that flattened the backbone to φ≈0°, ψ≈0° (Ramachandran forbidden).
This was **not** caused by the D-Ala2 mutation itself, but by the preparation pipeline.

| Stage | φ (Ala2) | ψ (Ala2) | χ1 | Status |
|-------|----------|----------|-----|--------|
| Native GHRH (7CZ5) | ~96° | ~41° | +126° | Reference |
| After `pair_fit()` (bug) | ~0° | ~0° | −116° | ❌ Forbidden |
| Post-EM | +174° | +177° | −115.68° | ❌ Trans-linear |
| After ~12.6 ns MD | −106° | −102° | — | ✅ Beta-sheet |

**Conclusion**: D-Ala2 mutation correctly flips χ1. The backbone flattening is a pipeline artifact,
not a physical consequence of D-alanination. The system self-corrected during MD to a valid
beta-sheet conformation.

### 2. Rosetta Refine Lacks Catalytic Constraints

Rosetta `refine-only` (`-pep_refine` without constraints) drifts the peptide away from
the catalytic pose. This is a **sampling protocol limitation**, not evidence that D-Ala2
destroys binding.

| Model | I_sc | rmsBB | Ser630→Ala2 C | Pass |
|-------|------|-------|---------------|------|
| WT start | — | 0.0 | 2.82 Å | 4/4 |
| WT model 1 | 292 | 1.6 Å | 5.83 Å | 2/4 |
| WT model 4 (best) | 161 | 3.8 Å | 4.61 Å | 2/4 |
| D-Ala2 start | — | 0.0 | 2.82 Å | 4/4 |
| D-Ala2 model 3 | 35 | 7.0 Å | 9.67 Å | 0/4 |
| D-Ala2 model 4 | 280 | 5.4 Å | 3.81 Å | 2/4 |

**Critical insight**: With only n=6–7 models, we can only state a **trend**: lower I_sc
correlates with greater backbone RMS drift (rmsBB 1.6–7 Å). This is expected for unconstrained
refinement of a 29-residue peptide. **The unrefined starting pose remains the best catalytic
configuration.** Future refine runs should use `-native` + CoordinateConstraint on GHRH 1–5 CA.

### 3. D-Ala2 EM Backbone Response

- χ1 successfully flipped: +126° → −116° ✅
- Backbone forced toward trans-linear by EM: φ −146° → +174°, ψ +159° → +177°
  - This is likely because the starting structure was already flattened by `pair_fit()`,
    and EM further pushed it toward the nearest local minimum (trans).
  - The improper dihedral in Amber14sb DALA (copied from ALA) may contribute, but
    the dominant driver is the unnatural starting geometry.
- Catalytic geometry slightly degraded post-EM: Ser630→Ala2 C 3.47→3.62 Å, oxyanion hole 3.50→3.88 Å

---

## Technical Blocker: GROMACS Segfault — RESOLVED

### Final Root Cause (2026-05-14)

**NOT** a GROMACS 2026.0 bug. **NOT** an improper dihedral issue.

**True root cause**: `prepare_docking_start.py` used `cmd.fab("YAD", ss=0)` to build a
linear peptide reference (φ=ψ=0), then `pair_fit()` forced GHRH backbone to match it.
Result: D-Ala2 starting structure had φ≈0°, ψ≈0° — Ramachandran forbidden.

- dt=0.002 → violent escape from forbidden region → NaN → segfault
- dt=0.001 → gentle relaxation → system escapes to beta-sheet (φ=-106°, ψ=-102°)
- Starting from 10+ ns checkpoint with dt=0.002 → stable (already in allowed region)

### Fix Applied
- `prepare_structures.py` / `prepare_structures_v2.py`: replaced `fab()` with native YAD extraction
- `prepare_docking_start.py`: added phi/psi validation with WARNING on forbidden region

### Current Status
- WT MD: ~40 ns / 200 ns, GPU 0, stable ✅
- D-Ala2 MD: ~12.9 ns / 200 ns, GPU 1, dt=0.001, stable ✅

---

## Recommendations

### Immediate (P0)
1. **Continue both MD runs to 200 ns** — both stable, no restart needed
2. **Monitor Ser630→Ala2 distance drift** in both trajectories

### Short-term (P1)
3. **Re-run Rosetta refine with constraints**:
   - `-native` parameter (currently missing, causing RMS self-reference)
   - CoordinateConstraint on GHRH 1–5 CA atoms
4. **Verify future structures** with `prepare_docking_start.py` phi/psi check

### Long-term (P2)
5. **Clean up step2** bloated silent files — ✅ done
6. **Publish methodology note** on pipeline backbone flattening bug

---

## Files Updated
- `prepare_structures.py` — replaced `fab()` with native YAD extraction
- `prepare_structures_v2.py` — same fix
- `prepare_docking_start.py` — added phi/psi/omega validation
- `validate_docked_pose.py` — chain-agnostic, supports ALA/DALA lookup
- `monitor_md_distances.py` — MDAnalysis-based trajectory monitor

*Report revised: 2026-05-14*
