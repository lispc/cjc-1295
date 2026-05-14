#!/usr/bin/env python3
"""
Draw side-by-side 3D comparison of WT vs D-Ala2 from actual MD trajectories.
Uses the last frame of each production run (~40 ns WT, ~12.9 ns D-Ala2).
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_gro_atoms(path):
    atoms = {}
    with open(path) as f:
        lines = f.readlines()
    n_atoms = int(lines[1].strip())
    for line in lines[2:2+n_atoms]:
        if len(line) >= 44:
            resi_resn = line[:5].strip()
            resi = int(''.join(c for c in resi_resn if c.isdigit()))
            resn = ''.join(c for c in resi_resn if c.isalpha())
            atom_name = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            atoms[(resi, atom_name)] = np.array([x, y, z])
    return atoms

def get_atom(atoms, resi, name):
    return atoms.get((resi, name))

def draw_backbone(ax, atoms, residues, color, linewidth=4):
    coords = []
    for resi in residues:
        ca = get_atom(atoms, resi, 'CA')
        if ca is not None:
            coords.append(ca)
    coords = np.array(coords)
    if len(coords) >= 2:
        for i in range(len(coords)-1):
            ax.plot([coords[i,0], coords[i+1,0]],
                   [coords[i,1], coords[i+1,1]],
                   [coords[i,2], coords[i+1,2]],
                   color=color, linewidth=linewidth, solid_capstyle='round')
    return coords

def draw_sphere(ax, center, radius, color, alpha=0.5, res=15):
    u = np.linspace(0, 2*np.pi, res)
    v = np.linspace(0, np.pi, res)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0)

def draw_line(ax, p1, p2, color, label, offset=(0,0,0)):
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
            color=color, linewidth=2, linestyle='--', alpha=0.7)
    mid = (p1 + p2) / 2 + np.array(offset)
    ax.text(mid[0], mid[1], mid[2], label, color=color, fontsize=10,
            fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=color, alpha=0.9))

# Load MD frames
WT = read_gro_atoms("/tmp/wt_md_last.gro")
DAL = read_gro_atoms("/tmp/dala_md_last.gro")

# Extract key atoms
wt_ser_og = get_atom(WT, 630, 'OG')
wt_ala2_c = get_atom(WT, 2, 'C')
wt_ala2_ca = get_atom(WT, 2, 'CA')
wt_ala2_n = get_atom(WT, 2, 'N')
wt_ala2_o = get_atom(WT, 2, 'O')
wt_tyr1_n = get_atom(WT, 1, 'N')
wt_tyr547 = get_atom(WT, 547, 'OH')
wt_glu205 = get_atom(WT, 205, 'OE1')

dal_ser_og = get_atom(DAL, 630, 'OG')
dal_ala2_c = get_atom(DAL, 2, 'C')
dal_ala2_ca = get_atom(DAL, 2, 'CA')
dal_ala2_n = get_atom(DAL, 2, 'N')
dal_ala2_o = get_atom(DAL, 2, 'O')
dal_tyr1_n = get_atom(DAL, 1, 'N')
dal_tyr547 = get_atom(DAL, 547, 'OH')
dal_glu205 = get_atom(DAL, 205, 'OE1')

# Center both systems on Ser630 OG for comparison
center_wt = wt_ser_og.copy()
center_dal = dal_ser_og.copy()

for d in [WT, DAL]:
    for k in d:
        pass  # will center in plotting

fig = plt.figure(figsize=(18, 9))
fig.patch.set_facecolor('white')

def setup_axes(ax, title, title_color):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20, color=title_color)
    ax.set_facecolor('white')
    ax.grid(False)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('gray')
        pane.set_alpha(0.1)
    ax.set_xlabel('X (nm)', fontsize=10)
    ax.set_ylabel('Y (nm)', fontsize=10)
    ax.set_zlabel('Z (nm)', fontsize=10)

view_elev, view_azim = 20, -60

# ========================================================================
# LEFT: WT MD ~40 ns
# ========================================================================
ax1 = fig.add_subplot(121, projection='3d')
setup_axes(ax1, 'WT GHRH + DPP-IV\nAfter ~40 ns MD', '#2E7D32')

# Center on Ser630 OG
c = wt_ser_og
xlim = (c[0]-1.5, c[0]+1.5)
ylim = (c[1]-1.5, c[1]+1.5)
zlim = (c[2]-1.5, c[2]+1.5)

# Draw GHRH N-terminus (residues 1-8)
draw_backbone(ax1, WT, list(range(1, 9)), '#00BCD4', 5)

# Draw DPP-IV catalytic region (Ser630)
ax1.scatter([wt_ser_og[0]], [wt_ser_og[1]], [wt_ser_og[2]],
           c='#FF1744', s=400, depthshade=False, edgecolors='darkred', linewidth=2)
draw_sphere(ax1, wt_ser_og, 0.06, '#FF1744', 0.4)
ax1.text(wt_ser_og[0]-0.15, wt_ser_og[1], wt_ser_og[2]+0.1,
        'Ser630 OG\n(nucleophile)', fontsize=9, color='#B71C1C', fontweight='bold', ha='right')

# Ala2 C=O
ax1.scatter([wt_ala2_c[0]], [wt_ala2_c[1]], [wt_ala2_c[2]],
           c='#1976D2', s=300, depthshade=False, edgecolors='darkblue', linewidth=2)
ax1.scatter([wt_ala2_o[0]], [wt_ala2_o[1]], [wt_ala2_o[2]],
           c='#FF5722', s=200, depthshade=False, edgecolors='darkred', linewidth=1.5)
ax1.text(wt_ala2_c[0]+0.1, wt_ala2_c[1], wt_ala2_c[2]-0.05,
        'Ala2 C=O', fontsize=9, color='#0D47A1', fontweight='bold')

# Distances
draw_line(ax1, wt_ser_og, wt_ala2_c, '#C62828', '4.84 Å', offset=(0, 0, 0.1))

# Tyr1 N
draw_line(ax1, wt_tyr1_n, wt_glu205, '#C62828', '9.27 Å', offset=(0, 0, 0.15))

# Annotation
ax1.text2D(0.02, 0.98,
          'WT ~40 ns MD:\n'
          '  Ser630->Ala2 C = 4.84 Å\n'
          '  Attack angle = 49.5 deg\n'
          '  Tyr1->Glu205 = 9.27 Å\n'
          '  PASS: 0/4 (drifted)',
          transform=ax1.transAxes, fontsize=10, verticalalignment='top',
          bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF3E0', edgecolor='#E65100', linewidth=2),
          family='monospace', color='#BF360C')

ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_zlim(zlim)
ax1.view_init(elev=view_elev, azim=view_azim)

# ========================================================================
# RIGHT: D-Ala2 MD ~12.9 ns
# ========================================================================
ax2 = fig.add_subplot(122, projection='3d')
setup_axes(ax2, 'D-Ala2 GHRH + DPP-IV\nAfter ~12.9 ns MD', '#C62828')

c = dal_ser_og
xlim = (c[0]-1.5, c[0]+1.5)
ylim = (c[1]-1.5, c[1]+1.5)
zlim = (c[2]-1.5, c[2]+1.5)

# Draw GHRH N-terminus
draw_backbone(ax2, DAL, list(range(1, 9)), '#00BCD4', 5)

# Ser630 OG
ax2.scatter([dal_ser_og[0]], [dal_ser_og[1]], [dal_ser_og[2]],
           c='#FF1744', s=400, depthshade=False, edgecolors='darkred', linewidth=2)
draw_sphere(ax2, dal_ser_og, 0.06, '#FF1744', 0.4)
ax2.text(dal_ser_og[0]-0.15, dal_ser_og[1], dal_ser_og[2]+0.1,
        'Ser630 OG\n(nucleophile)', fontsize=9, color='#B71C1C', fontweight='bold', ha='right')

# Ala2 C=O
ax2.scatter([dal_ala2_c[0]], [dal_ala2_c[1]], [dal_ala2_c[2]],
           c='#1976D2', s=300, depthshade=False, edgecolors='darkblue', linewidth=2)
ax2.scatter([dal_ala2_o[0]], [dal_ala2_o[1]], [dal_ala2_o[2]],
           c='#FF5722', s=200, depthshade=False, edgecolors='darkred', linewidth=1.5)
ax2.text(dal_ala2_c[0]+0.1, dal_ala2_c[1], dal_ala2_c[2]-0.05,
        'Ala2 C=O', fontsize=9, color='#0D47A1', fontweight='bold')

# Distances
draw_line(ax2, dal_ser_og, dal_ala2_c, '#C62828', '6.40 Å', offset=(0, 0, 0.1))

# Tyr1 N
draw_line(ax2, dal_tyr1_n, dal_glu205, '#C62828', '6.08 Å', offset=(0, 0, 0.15))

# Drift arrow from WT position
if wt_ala2_ca is not None and dal_ala2_ca is not None:
    ax2.quiver(wt_ala2_ca[0], wt_ala2_ca[1], wt_ala2_ca[2],
              dal_ala2_ca[0]-wt_ala2_ca[0], dal_ala2_ca[1]-wt_ala2_ca[1], dal_ala2_ca[2]-wt_ala2_ca[2],
              color='#C62828', arrow_length_ratio=0.2, linewidth=2, alpha=0.6)
    mid = (wt_ala2_ca + dal_ala2_ca) / 2
    ax2.text(mid[0], mid[1]+0.1, mid[2],
            'Additional drift\nvs WT', fontsize=8, color='#C62828',
            ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#C62828', alpha=0.7))

# Annotation
ax2.text2D(0.02, 0.98,
          'D-Ala2 ~12.9 ns MD:\n'
          '  Ser630->Ala2 C = 6.40 Å\n'
          '  Attack angle = 86.6 deg\n'
          '  Tyr1->Glu205 = 6.08 Å\n'
          '  PASS: 1/4 (more drifted)',
          transform=ax2.transAxes, fontsize=10, verticalalignment='top',
          bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFEBEE', edgecolor='#C62828', linewidth=2),
          family='monospace', color='#B71C1C')

ax2.set_xlim(xlim)
ax2.set_ylim(ylim)
ax2.set_zlim(zlim)
ax2.view_init(elev=view_elev, azim=view_azim)

fig.suptitle('MD Trajectory Comparison: WT vs D-Ala2 Catalytic Geometry',
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.95])
outpath = "/home/scroll/personal/cjc-1295/workspace/figures/md_alignment_WT_vs_DAla2.png"
plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {outpath}")
plt.close()
