#!/usr/bin/env python3
"""
Figure: D-Ala2 inhibits DPP-IV by locking the peptide in a non-productive pose.
The key finding: D-Ala2 binds tightly (close to Ser630) but at the wrong angle.
L-Ala can achieve productive geometry but rarely stays in the pocket.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch, Arc, Ellipse
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from scipy import stats

# ─── Load data ─────────────────────────────────────────────────
# NOTE: chirality corrected — "wt" files are actually D-Ala2 (CJC-1295)
# "dala2" files are actually L-Ala (natural GHRH)

def load_xvg(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

# D-Ala2 (CJC-1295) — from part 2 (9.4-53.5 ns of the original "WT" run)
dist_D = load_xvg("/home/scroll/personal/cjc-1295/workspace/step3/dist_wt_p2_d1.xvg")
ang_D = load_xvg("/home/scroll/personal/cjc-1295/workspace/step3/ang_wt_p2.xvg")

# L-Ala (natural GHRH)
dist_L = load_xvg("/home/scroll/personal/cjc-1295/workspace/step3/dist_dala2_d1.xvg")
ang_L = load_xvg("/home/scroll/personal/cjc-1295/workspace/step3/ang_dala2.xvg")

# Align time series — use the shorter one for both
n_D = min(len(dist_D), len(ang_D))
n_L = min(len(dist_L), len(ang_L))
dist_D, ang_D = dist_D[:n_D], ang_D[:n_D]
dist_L, ang_L = dist_L[:n_L], ang_L[:n_L]

# Convert nm to Å
dist_D_nm = dist_D[:, 1].copy()
dist_L_nm = dist_L[:, 1].copy()
dist_D_A = dist_D[:, 1] * 10
dist_L_A = dist_L[:, 1] * 10
ang_D_deg = ang_D[:, 1]
ang_L_deg = ang_L[:, 1]

# ─── Productive zone mask ──────────────────────────────────────
prod_mask_D = (dist_D_A >= 2.5) & (dist_D_A <= 4.0) & (ang_D_deg >= 80) & (ang_D_deg <= 120)
prod_mask_L = (dist_L_A >= 2.5) & (dist_L_A <= 4.0) & (ang_L_deg >= 80) & (ang_L_deg <= 120)

# ─── Setup figure ──────────────────────────────────────────────
fig = plt.figure(figsize=(16, 7))

# ── PANEL A: 2D scatter — Distance vs Angle ────────────────────
ax1 = fig.add_subplot(1, 2, 1)

# Productive zone
rect = Rectangle((2.5, 80), 1.5, 40, linewidth=0, facecolor='#4CAF50', alpha=0.10, zorder=1)
ax1.add_patch(rect)
ax1.text(3.25, 82, 'Productive\nZone', ha='center', va='bottom', fontsize=11,
         color='#2E7D32', fontweight='bold', alpha=0.7)

# D-Ala2 contour
ax1.scatter(dist_D_A[::5], ang_D_deg[::5], c='#E63946', alpha=0.18, s=4,
           rasterized=True, label='D-Ala2 (CJC-1295)', zorder=2)

# L-Ala contour
ax1.scatter(dist_L_A[::5], ang_L_deg[::5], c='#457B9D', alpha=0.18, s=4,
           rasterized=True, label='L-Ala2 (native GHRH)', zorder=2)

# Mean markers
ax1.scatter([np.mean(dist_D_A)], [np.mean(ang_D_deg)], c='#E63946', s=200, marker='X',
            edgecolors='white', linewidths=1.2, zorder=5)
ax1.scatter([np.mean(dist_L_A)], [np.mean(ang_L_deg)], c='#457B9D', s=200, marker='X',
            edgecolors='white', linewidths=1.2, zorder=5)

# Annotations
ax1.annotate('D-Ala²\nclose but\nwrong angle', xy=(np.mean(dist_D_A), np.mean(ang_D_deg)),
            xytext=(5.5, 65), fontsize=9, color='#C1121F', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#C1121F', lw=1.5),
            ha='center', va='top')

ax1.annotate('L-Ala²\nfar from\npocket', xy=(np.mean(dist_L_A), np.mean(ang_L_deg)),
            xytext=(8.0, 95), fontsize=9, color='#1D3557', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#1D3557', lw=1.5),
            ha='center', va='center')

ax1.set_xlabel('Ser630 Oγ → Ala2 C distance (Å)', fontsize=13)
ax1.set_ylabel('Attack angle ∠(Oγ–C–N) (°)', fontsize=13)
ax1.set_xlim(1.5, 9.5)
ax1.set_ylim(50, 130)
ax1.axhline(y=90, color='gray', linestyle=':', alpha=0.3)
ax1.legend(loc='upper left', fontsize=10, framealpha=0.8)
ax1.set_title('A. Catalytic Geometry Landscape', fontsize=14, fontweight='bold', loc='left')

# Stats text
stats_text = (
    f"D-Ala²: {np.mean(dist_D_A):.1f}±{np.std(dist_D_A):.1f} Å, "
    f"{np.mean(ang_D_deg):.1f}±{np.std(ang_D_deg):.1f}°\n"
    f"  Productive: {prod_mask_D.mean()*100:.1f}%\n"
    f"L-Ala²:  {np.mean(dist_L_A):.1f}±{np.std(dist_L_A):.1f} Å, "
    f"{np.mean(ang_L_deg):.1f}±{np.std(ang_L_deg):.1f}°\n"
    f"  Productive: {prod_mask_L.mean()*100:.1f}%"
)
ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=9,
         fontfamily='monospace', verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))

# ── PANEL B: Distance distribution + FEP inset ──────────────────
ax2 = fig.add_subplot(2, 2, 2)

bins = np.linspace(1.5, 10, 60)
ax2.hist(dist_D_A, bins=bins, color='#E63946', alpha=0.55, label='D-Ala² (CJC-1295)',
         edgecolor='#C1121F', linewidth=0.3)
ax2.hist(dist_L_A, bins=bins, color='#457B9D', alpha=0.55, label='L-Ala² (native GHRH)',
         edgecolor='#1D3557', linewidth=0.3)

# Productive zone shading
ax2.axvspan(2.5, 4.0, color='#4CAF50', alpha=0.12)

# Mean lines
ax2.axvline(np.mean(dist_D_A), color='#E63946', linestyle='--', linewidth=2)
ax2.axvline(np.mean(dist_L_A), color='#457B9D', linestyle='--', linewidth=2)

ax2.set_xlabel('Ser630 Oγ → Ala2 C distance (Å)', fontsize=12)
ax2.set_ylabel('Frame count', fontsize=12)
ax2.set_title('B. Distance Distribution', fontsize=13, fontweight='bold', loc='left')
ax2.legend(fontsize=9)

# ── PANEL C: Schematic of the mechanism ─────────────────────────
ax3 = fig.add_subplot(2, 2, 4)
ax3.set_xlim(0, 14)
ax3.set_ylim(0, 8)
ax3.axis('off')
ax3.set_title('C. Inhibition Mechanism', fontsize=13, fontweight='bold', loc='left')

# DPP-IV pocket (schematic)
pocket = FancyBboxPatch((0.8, 1.0), 5.5, 6.0, boxstyle="round,pad=0.3",
                         facecolor='#F0F0F0', edgecolor='#888888', linewidth=1.5)
ax3.add_patch(pocket)
ax3.text(3.55, 7.5, 'DPP-IV Catalytic Pocket', ha='center', fontsize=11,
         fontweight='bold', color='#333333')

# Ser630
ax3.plot(2.8, 3.5, 'o', color='#E76F51', markersize=18, markeredgecolor='white',
         markeredgewidth=1, zorder=5)
ax3.text(2.8, 2.4, 'Ser630\n(nucleophile)', ha='center', fontsize=9, color='#E76F51',
         fontweight='bold')

# Oxyanion hole
ax3.plot(4.5, 1.8, 's', color='#2A9D8F', markersize=14, markeredgecolor='white',
         markeredgewidth=1, zorder=4)
ax3.text(4.5, 1.2, 'Oxyanion\nhole', ha='center', fontsize=8, color='#2A9D8F')

# L-Ala substrate
l_ala_pos = (6.5, 5.0)
ax3.plot(l_ala_pos[0], l_ala_pos[1], 'o', color='#457B9D', markersize=22,
         markeredgecolor='white', markeredgewidth=1.5, zorder=5)
ax3.text(l_ala_pos[0], l_ala_pos[1], 'L', ha='center', va='center', fontsize=10,
         color='white', fontweight='bold')
ax3.text(l_ala_pos[0], l_ala_pos[1]-0.8, 'L-Ala²\n(can be cleaved)\nproductive angle\nbut drifts away',
         ha='center', fontsize=9, color='#457B9D', fontweight='bold')

# D-Ala inhibitor
d_ala_pos = (0.3, 5.5)  # outside pocket (left side misaligned)
# Actually, D-Ala is IN the pocket but at wrong angle. Let me put it inside.
d_ala_pos = (3.1, 5.8)
ax3.plot(d_ala_pos[0], d_ala_pos[1], 'o', color='#E63946', markersize=22,
         markeredgecolor='white', markeredgewidth=1.5, zorder=5)
ax3.text(d_ala_pos[0], d_ala_pos[1], 'D', ha='center', va='center', fontsize=10,
         color='white', fontweight='bold')
ax3.text(d_ala_pos[0], d_ala_pos[1]+0.8, 'D-Ala²\n(blocks cleavage)\ntightly bound\nwrong angle',
         ha='center', fontsize=9, color='#E63946', fontweight='bold')

# Attack arrows
# L-Ala: correct angle arrow
arrow_L = FancyArrowPatch((l_ala_pos[0]-0.5, l_ala_pos[1]-0.3),
                           (2.8+0.3, 3.5+0.3),
                           arrowstyle='->', color='#457B9D', lw=2, ls='--',
                           connectionstyle="arc3,rad=-0.2")
ax3.add_patch(arrow_L)
ax3.text(5.3, 3.9, 'correct\nangle', fontsize=8, color='#2E7D32', fontweight='bold')

# D-Ala: wrong angle arrow
arrow_D = FancyArrowPatch((d_ala_pos[0]-0.3, d_ala_pos[1]-0.2),
                           (2.8+0.1, 3.5-0.1),
                           arrowstyle='->', color='#E63946', lw=2, ls='--',
                           connectionstyle="arc3,rad=0.3")
ax3.add_patch(arrow_D)
ax3.text(2.0, 5.5, 'wrong\nangle', fontsize=8, color='#C1121F', fontweight='bold',
          ha='center')

# FEP result box
fep_box = FancyBboxPatch((9.0, 4.8), 4.5, 2.5, boxstyle="round,pad=0.2",
                          facecolor='white', edgecolor='#333333', linewidth=1.2, zorder=10)
ax3.add_patch(fep_box)
fep_text = (
    "FEP DG(D->L) = +0.83 kJ/mol\n"
    "<< kT (2.57 kJ/mol)\n\n"
    "Binding affinity ~ EQUAL\n"
    "Effect is purely geometric"
)
ax3.text(11.25, 6.05, fep_text, ha='center', va='center', fontsize=9,
         fontfamily='monospace', zorder=11)

# GHRH peptide cartoon
ax3.plot([7.2, 9.5], [5.0, 4.5], 'k-', lw=3, alpha=0.4)
ax3.text(9.8, 4.3, 'GHRH\nC-term', fontsize=7, color='gray')

ax3.plot([-0.5, 1.0], [5.8, 6.2], 'k-', lw=3, alpha=0.4)
ax3.text(-1.2, 6.3, 'GHRH\nN-term', fontsize=7, color='gray')

# FEP energy diagram inset
ax_energy = ax3.inset_axes([0.73, 0.05, 0.25, 0.35])
lambdas = np.arange(11)
dgs = np.array([0, -4.33, -4.15, -3.47, -2.17, -0.76, 1.00, 2.38, 3.27, 4.15, 4.93])
dg_cum = np.cumsum(dgs)
dg_cum -= dg_cum[0]  # start from 0

ax_energy.plot(lambdas, dg_cum, 'ko-', markersize=4, linewidth=1.5, color='#333333')
ax_energy.axhline(y=-0.83, color='#E63946', linestyle='--', linewidth=1, alpha=0.6)
ax_energy.fill_between(lambdas, dg_cum-0.3, dg_cum+0.3, alpha=0.1, color='#457B9D')
ax_energy.set_xlabel('λ (0=D-Ala, 1=L-Ala)', fontsize=7)
ax_energy.set_ylabel('ΔG (kJ/mol)', fontsize=7)
ax_energy.set_title('FEP Profile', fontsize=8, fontweight='bold')
ax_energy.tick_params(labelsize=6)
ax_energy.text(5, 1.5, 'ΔG = +0.83 ± 0.31', fontsize=7, color='#E63946',
               ha='center', fontweight='bold')

# ── Final layout ────────────────────────────────────────────────
plt.tight_layout(pad=2.5)
outpath = "/home/scroll/personal/cjc-1295/workspace/results/mechanism_figure.png"
plt.savefig(outpath, dpi=250, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
print(f"D-Ala2 productive: {prod_mask_D.mean()*100:.1f}%")
print(f"L-Ala2 productive:  {prod_mask_L.mean()*100:.1f}%")
print(f"D-Ala2 dist: {np.mean(dist_D_A):.2f} ± {np.std(dist_D_A):.2f} Å")
print(f"L-Ala2 dist:  {np.mean(dist_L_A):.2f} ± {np.std(dist_L_A):.2f} Å")
print(f"D-Ala2 angle: {np.mean(ang_D_deg):.1f} ± {np.std(ang_D_deg):.1f}°")
print(f"L-Ala2 angle:  {np.mean(ang_L_deg):.1f} ± {np.std(ang_L_deg):.1f}°")
