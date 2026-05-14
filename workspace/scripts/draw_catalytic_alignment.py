#!/usr/bin/env python3
"""
Draw side-by-side 3D comparison of WT vs D-Ala2 catalytic geometry.
Left: WT starting pose (perfect alignment, 2.82 Å)
Right: D-Ala2 after Rosetta refine drift (misaligned, 9.67 Å)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.patches as mpatches

def read_pdb_atoms(path):
    atoms = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21]
                resi = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atoms[(chain, resi, resn, name)] = np.array([x, y, z])
    return atoms

def get_atom(atoms, chain, resi, resn, name):
    key = (chain, resi, resn, name)
    if key in atoms:
        return atoms[key]
    # Fallback for D-Ala2
    alt = [(c, r, n, a) for (c, r, n, a) in atoms.keys() if r == resi and a == name]
    if alt:
        return atoms[alt[0]]
    return None

def draw_backbone(ax, atoms, chain, residues, color, linewidth=3, alpha=1.0):
    """Draw peptide backbone as a tube."""
    coords = []
    for resi, resn in residues:
        ca = get_atom(atoms, chain, resi, resn, "CA")
        if ca is not None:
            coords.append(ca)
    coords = np.array(coords)
    if len(coords) >= 2:
        for i in range(len(coords) - 1):
            ax.plot([coords[i,0], coords[i+1,0]],
                   [coords[i,1], coords[i+1,1]],
                   [coords[i,2], coords[i+1,2]],
                   color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='round')
    return coords

def draw_sphere(ax, center, radius, color, alpha=0.6, resolution=20):
    """Draw a sphere at center."""
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0)

def draw_distance_line(ax, p1, p2, color, label, offset=(0,0,0)):
    """Draw a dashed line between two points with a label."""
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
            color=color, linewidth=2, linestyle='--', alpha=0.8)
    mid = (p1 + p2) / 2 + np.array(offset)
    ax.text(mid[0], mid[1], mid[2], label, color=color, fontsize=11,
            fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=color, alpha=0.9))

# ---------------------------------------------------------------------------
# Load structures
# ---------------------------------------------------------------------------
WT = read_pdb_atoms("/home/scroll/personal/cjc-1295/workspace/step2/prepacked_DPP4_GHRH_start_0001.pdb")
DAL = read_pdb_atoms("/home/scroll/personal/cjc-1295/workspace/step2/GHRH_DPP4_docked_best_RosettaDAL_0003.pdb")

# WT key atoms
wt_ser630_og = get_atom(WT, 'A', 630, 'SER', 'OG')
wt_ser630_cb = get_atom(WT, 'A', 630, 'SER', 'CB')
wt_ala2_c = get_atom(WT, 'B', 2, 'ALA', 'C')
wt_ala2_ca = get_atom(WT, 'B', 2, 'ALA', 'CA')
wt_ala2_n = get_atom(WT, 'B', 2, 'ALA', 'N')
wt_ala2_o = get_atom(WT, 'B', 2, 'ALA', 'O')
wt_ala2_cb = get_atom(WT, 'B', 2, 'ALA', 'CB')
wt_tyr1_n = get_atom(WT, 'B', 1, 'TYR', 'N')
wt_tyr547_oh = get_atom(WT, 'A', 547, 'TYR', 'OH')
wt_glu205_oe1 = get_atom(WT, 'A', 205, 'GLU', 'OE1')
wt_glu206_oe1 = get_atom(WT, 'A', 206, 'GLU', 'OE1')

# D-Ala2 key atoms
dal_ser630_og = get_atom(DAL, 'A', 630, 'SER', 'OG')
dal_ser630_cb = get_atom(DAL, 'A', 630, 'SER', 'CB')
dal_ala2_c = get_atom(DAL, 'B', 2, 'DAL', 'C')
dal_ala2_ca = get_atom(DAL, 'B', 2, 'DAL', 'CA')
dal_ala2_n = get_atom(DAL, 'B', 2, 'DAL', 'N')
dal_ala2_o = get_atom(DAL, 'B', 2, 'DAL', 'O')
dal_ala2_cb = get_atom(DAL, 'B', 2, 'DAL', 'CB')
dal_tyr1_n = get_atom(DAL, 'B', 1, 'TYR', 'N')
dal_tyr547_oh = get_atom(DAL, 'A', 547, 'TYR', 'OH')
dal_glu205_oe1 = get_atom(DAL, 'A', 205, 'GLU', 'OE1')
dal_glu206_oe1 = get_atom(DAL, 'A', 206, 'GLU', 'OE1')

# ---------------------------------------------------------------------------
# Create figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(18, 9))
fig.patch.set_facecolor('white')

# Common view parameters
view_elev = 20
view_azim = -60
xlim = (62, 73)
ylim = (64, 72)
zlim = (63, 78)

# ========================================================================
# LEFT PANEL: WT — Perfect Alignment
# ========================================================================
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_facecolor('white')
ax1.set_title('WT GHRH(1-29) + DPP-IV\nPerfect Catalytic Alignment', 
              fontsize=14, fontweight='bold', pad=20, color='#2E7D32')

# Draw GHRH N-terminal backbone (residues 1-5)
wt_bb = draw_backbone(ax1, WT, 'B', [(1,'TYR'), (2,'ALA'), (3,'ASP'), (4,'ALA'), (5,'ILE')], 
                       color='#00BCD4', linewidth=5)

# Draw DPP-IV catalytic residues (simplified)
ax1.scatter([wt_ser630_cb[0]], [wt_ser630_cb[1]], [wt_ser630_cb[2]], 
           c='#FF5252', s=150, depthshade=False, edgecolors='black', linewidth=1.5)
ax1.text(wt_ser630_cb[0]+0.5, wt_ser630_cb[1]+0.5, wt_ser630_cb[2]+0.5,
        'Ser630', fontsize=9, color='#D32F2F', fontweight='bold')

# Ser630 OG — nucleophile (large red sphere)
draw_sphere(ax1, wt_ser630_og, 0.6, '#FF1744', alpha=0.5)
ax1.scatter([wt_ser630_og[0]], [wt_ser630_og[1]], [wt_ser630_og[2]],
           c='#FF1744', s=300, depthshade=False, edgecolors='darkred', linewidth=2, marker='o')
ax1.text(wt_ser630_og[0]-1.2, wt_ser630_og[1], wt_ser630_og[2]+0.8,
        'Ser630\nOG⁻', fontsize=10, color='#B71C1C', fontweight='bold',
        ha='right', va='bottom')

# Ala2 C=O (carbonyl carbon + oxygen)
ax1.scatter([wt_ala2_c[0]], [wt_ala2_c[1]], [wt_ala2_c[2]],
           c='#1976D2', s=250, depthshade=False, edgecolors='darkblue', linewidth=2, marker='o')
ax1.scatter([wt_ala2_o[0]], [wt_ala2_o[1]], [wt_ala2_o[2]],
           c='#FF5722', s=200, depthshade=False, edgecolors='darkred', linewidth=1.5, marker='o')
ax1.text(wt_ala2_c[0]+0.8, wt_ala2_c[1], wt_ala2_c[2]-0.5,
        'Ala2\nC=O', fontsize=10, color='#0D47A1', fontweight='bold',
        ha='left', va='top')

# Draw attack distance
draw_distance_line(ax1, wt_ser630_og, wt_ala2_c, '#2E7D32', '2.82 Å', offset=(0, 0.8, 0))

# Tyr1 N (salt bridge anchor)
ax1.scatter([wt_tyr1_n[0]], [wt_tyr1_n[1]], [wt_tyr1_n[2]],
           c='#4CAF50', s=200, depthshade=False, edgecolors='darkgreen', linewidth=1.5, marker='s')
ax1.text(wt_tyr1_n[0]-0.5, wt_tyr1_n[1]+0.5, wt_tyr1_n[2]+0.5,
        'Tyr1\nNH₃⁺', fontsize=9, color='#1B5E20', fontweight='bold', ha='right')

# Glu205/206 (salt bridge partners)
ax1.scatter([wt_glu205_oe1[0]], [wt_glu205_oe1[1]], [wt_glu205_oe1[2]],
           c='#8BC34A', s=120, depthshade=False, marker='^', alpha=0.7)
ax1.text(wt_glu205_oe1[0], wt_glu205_oe1[1]+0.5, wt_glu205_oe1[2],
        'Glu205', fontsize=8, color='#33691E', alpha=0.8)

# Tyr547 OH (oxyanion hole)
ax1.scatter([wt_tyr547_oh[0]], [wt_tyr547_oh[1]], [wt_tyr547_oh[2]],
           c='#FFC107', s=150, depthshade=False, edgecolors='#FF6F00', linewidth=1.5, marker='D')
ax1.text(wt_tyr547_oh[0], wt_tyr547_oh[1]+0.3, wt_tyr547_oh[2]+0.5,
        'Tyr547', fontsize=8, color='#E65100')

# Draw salt bridge distance
draw_distance_line(ax1, wt_tyr1_n, wt_glu205_oe1, '#4CAF50', '2.02 Å', offset=(0, -0.5, 0.5))

# Draw oxyanion hole distance
draw_distance_line(ax1, wt_ala2_o, wt_tyr547_oh, '#FFC107', '3.10 Å', offset=(0.5, 0, 0))

# Legend / annotation box
ax1.text2D(0.02, 0.98, 
          'OK Ser630 OG -> Ala2 C = 2.82 A\n'
          'OK Attack angle = 113.6 deg\n'
          'OK Tyr1 NH3+ -> Glu205 = 2.02 A\n'
          'OK Ala2 O -> Tyr547 = 3.10 A\n'
          'OK PASS: 4/4 criteria',
          transform=ax1.transAxes, fontsize=10, verticalalignment='top',
          bbox=dict(boxstyle='round,pad=0.5', facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2),
          family='monospace', color='#1B5E20')

ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_zlim(zlim)
ax1.view_init(elev=view_elev, azim=view_azim)
ax1.set_xlabel('X (Å)', fontsize=10)
ax1.set_ylabel('Y (Å)', fontsize=10)
ax1.set_zlabel('Z (Å)', fontsize=10)

# Hide grid for cleaner look
ax1.grid(False)
ax1.xaxis.pane.fill = False
ax1.yaxis.pane.fill = False
ax1.zaxis.pane.fill = False
ax1.xaxis.pane.set_edgecolor('gray')
ax1.yaxis.pane.set_edgecolor('gray')
ax1.zaxis.pane.set_edgecolor('gray')
ax1.xaxis.pane.set_alpha(0.1)
ax1.yaxis.pane.set_alpha(0.1)
ax1.zaxis.pane.set_alpha(0.1)

# ========================================================================
# RIGHT PANEL: D-Ala2 — Misaligned (Rosetta refine drift)
# ========================================================================
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_facecolor('white')
ax2.set_title('D-Ala2 GHRH(1-29) + DPP-IV\nCatalytic Geometry Destroyed', 
              fontsize=14, fontweight='bold', pad=20, color='#C62828')

# Draw GHRH N-terminal backbone
dal_bb = draw_backbone(ax2, DAL, 'B', [(1,'TYR'), (2,'DAL'), (3,'ASP'), (4,'ALA'), (5,'ILE')], 
                        color='#00BCD4', linewidth=5)

# Draw DPP-IV catalytic residues (same position as WT)
ax2.scatter([dal_ser630_cb[0]], [dal_ser630_cb[1]], [dal_ser630_cb[2]], 
           c='#FF5252', s=150, depthshade=False, edgecolors='black', linewidth=1.5)
ax2.text(dal_ser630_cb[0]+0.5, dal_ser630_cb[1]+0.5, dal_ser630_cb[2]+0.5,
        'Ser630', fontsize=9, color='#D32F2F', fontweight='bold')

# Ser630 OG — nucleophile
draw_sphere(ax2, dal_ser630_og, 0.6, '#FF1744', alpha=0.5)
ax2.scatter([dal_ser630_og[0]], [dal_ser630_og[1]], [dal_ser630_og[2]],
           c='#FF1744', s=300, depthshade=False, edgecolors='darkred', linewidth=2, marker='o')
ax2.text(dal_ser630_og[0]-1.2, dal_ser630_og[1], dal_ser630_og[2]+0.8,
        'Ser630\nOG⁻', fontsize=10, color='#B71C1C', fontweight='bold',
        ha='right', va='bottom')

# Ala2 C=O (now far away!)
ax2.scatter([dal_ala2_c[0]], [dal_ala2_c[1]], [dal_ala2_c[2]],
           c='#1976D2', s=250, depthshade=False, edgecolors='darkblue', linewidth=2, marker='o')
ax2.scatter([dal_ala2_o[0]], [dal_ala2_o[1]], [dal_ala2_o[2]],
           c='#FF5722', s=200, depthshade=False, edgecolors='darkred', linewidth=1.5, marker='o')
ax2.text(dal_ala2_c[0]+0.8, dal_ala2_c[1], dal_ala2_c[2]-0.5,
        'Ala2\nC=O', fontsize=10, color='#0D47A1', fontweight='bold',
        ha='left', va='top')

# Draw FAILED attack distance (long red dashed line)
draw_distance_line(ax2, dal_ser630_og, dal_ala2_c, '#C62828', '9.67 Å X', offset=(0, 0.8, 0))

# Tyr1 N (now far from Glu)
ax2.scatter([dal_tyr1_n[0]], [dal_tyr1_n[1]], [dal_tyr1_n[2]],
           c='#4CAF50', s=200, depthshade=False, edgecolors='darkgreen', linewidth=1.5, marker='s')
ax2.text(dal_tyr1_n[0]-0.5, dal_tyr1_n[1]+0.5, dal_tyr1_n[2]+0.5,
        'Tyr1\nNH₃⁺', fontsize=9, color='#1B5E20', fontweight='bold', ha='right')

# Glu205/206
ax2.scatter([dal_glu205_oe1[0]], [dal_glu205_oe1[1]], [dal_glu205_oe1[2]],
           c='#8BC34A', s=120, depthshade=False, marker='^', alpha=0.7)
ax2.text(dal_glu205_oe1[0], dal_glu205_oe1[1]+0.5, dal_glu205_oe1[2],
        'Glu205', fontsize=8, color='#33691E', alpha=0.8)

# Tyr547 OH
ax2.scatter([dal_tyr547_oh[0]], [dal_tyr547_oh[1]], [dal_tyr547_oh[2]],
           c='#FFC107', s=150, depthshade=False, edgecolors='#FF6F00', linewidth=1.5, marker='D')
ax2.text(dal_tyr547_oh[0], dal_tyr547_oh[1]+0.3, dal_tyr547_oh[2]+0.5,
        'Tyr547', fontsize=8, color='#E65100')

# Draw broken salt bridge
draw_distance_line(ax2, dal_tyr1_n, dal_glu205_oe1, '#C62828', '10.90 Å X', offset=(0, -0.5, 0.5))

# Draw broken oxyanion hole
draw_distance_line(ax2, dal_ala2_o, dal_tyr547_oh, '#C62828', '6.95 Å X', offset=(0.5, 0, 0))

# Add drift arrow showing peptide displacement
if wt_ala2_ca is not None and dal_ala2_ca is not None:
    ax2.quiver(wt_ala2_ca[0], wt_ala2_ca[1], wt_ala2_ca[2],
              dal_ala2_ca[0]-wt_ala2_ca[0], dal_ala2_ca[1]-wt_ala2_ca[1], dal_ala2_ca[2]-wt_ala2_ca[2],
              color='#C62828', arrow_length_ratio=0.15, linewidth=3, alpha=0.7)
    ax2.text((wt_ala2_ca[0]+dal_ala2_ca[0])/2, (wt_ala2_ca[1]+dal_ala2_ca[1])/2+1, 
            (wt_ala2_ca[2]+dal_ala2_ca[2])/2,
            'Peptide drift\n(Rosetta refine)', fontsize=9, color='#C62828', 
            fontweight='bold', ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#C62828', alpha=0.8))

# Legend / annotation box
ax2.text2D(0.02, 0.98, 
          'X  Ser630 OG -> Ala2 C = 9.67 A\n'
          'X  Attack angle = 126.1 deg (off)\n'
          'X  Tyr1 NH3+ -> Glu205 = 10.90 A\n'
          'X  Ala2 O -> Tyr547 = 6.95 A\n'
          'X  FAIL: 0/4 criteria',
          transform=ax2.transAxes, fontsize=10, verticalalignment='top',
          bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFEBEE', edgecolor='#C62828', linewidth=2),
          family='monospace', color='#B71C1C')

ax2.set_xlim(xlim)
ax2.set_ylim(ylim)
ax2.set_zlim(zlim)
ax2.view_init(elev=view_elev, azim=view_azim)
ax2.set_xlabel('X (Å)', fontsize=10)
ax2.set_ylabel('Y (Å)', fontsize=10)
ax2.set_zlabel('Z (Å)', fontsize=10)

ax2.grid(False)
ax2.xaxis.pane.fill = False
ax2.yaxis.pane.fill = False
ax2.zaxis.pane.fill = False
ax2.xaxis.pane.set_edgecolor('gray')
ax2.yaxis.pane.set_edgecolor('gray')
ax2.zaxis.pane.set_edgecolor('gray')
ax2.xaxis.pane.set_alpha(0.1)
ax2.yaxis.pane.set_alpha(0.1)
ax2.zaxis.pane.set_alpha(0.1)

# Overall title
fig.suptitle('Catalytic Geometry: WT (Perfect Alignment) vs D-Ala2 (Misaligned after Refine)',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.95])
outpath = "/home/scroll/personal/cjc-1295/workspace/figures/catalytic_alignment_WT_vs_DAla2.png"
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
