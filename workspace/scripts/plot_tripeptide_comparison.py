#!/usr/bin/env python3
"""
Generate WT vs D-Ala2 tripeptide geometric comparison figure.
Pure matplotlib — no PyMOL dependency.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ============================================================
# Read PDB atoms
# ============================================================
def read_pdb_atoms(path):
    atoms = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                name = line[12:16].strip()
                resn = line[17:20].strip()
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[(resi, resn, name)] = np.array([x, y, z])
    return atoms

pdb = read_pdb_atoms("/home/scroll/personal/cjc-1295/workspace/step0/DPP4_with_YAD.pdb")

# ============================================================
# Extract WT Ala2 backbone atoms
# ============================================================
ala2_atoms = {k[2]: v for k, v in pdb.items() if k[0] == 2 and k[1] == "ALA"}
N  = ala2_atoms["N"]
CA = ala2_atoms["CA"]
C  = ala2_atoms["C"]
O  = ala2_atoms["O"]
CB = ala2_atoms["CB"]
HA = ala2_atoms.get("HA", None)

# ============================================================
# Build D-Ala2 by mirroring CB/HA across N-CA-C plane
# ============================================================
def point_to_plane_mirror(p, a, b, c):
    """Mirror point p across plane defined by points a,b,c."""
    n = np.cross(b - a, c - a)
    n = n / np.linalg.norm(n)
    d = np.dot(p - a, n)
    return p - 2 * d * n

CB_D = point_to_plane_mirror(CB, N, CA, C)
HA_D = point_to_plane_mirror(HA, N, CA, C) if HA is not None else None

# Verify dihedral sign flip
n1 = np.cross(N - CA, CB - CA)
n2 = np.cross(CB - CA, C - CA)
m1 = n1 / np.linalg.norm(n1)
m2 = n2 / np.linalg.norm(n2)
wt_dihedral = np.degrees(np.arctan2(np.dot(np.cross(m1, m2), CB - CA), np.dot(m1, m2)))

n1d = np.cross(N - CA, CB_D - CA)
n2d = np.cross(CB_D - CA, C - CA)
m1d = n1d / np.linalg.norm(n1d)
m2d = n2d / np.linalg.norm(n2d)
da_dihedral = np.degrees(np.arctan2(np.dot(np.cross(m1d, m2d), CB_D - CA), np.dot(m1d, m2d)))

print(f"WT  N-CA-CB-C dihedral: {wt_dihedral:.1f}°")
print(f"D-Ala N-CA-CB-C dihedral: {da_dihedral:.1f}°")

# ============================================================
# Create figure: 3D view + 2D schematic
# ============================================================
fig = plt.figure(figsize=(14, 6))

# --- Left: 3D stereo view of Ala2 backbone ---
ax1 = fig.add_subplot(121, projection='3d')

# Draw N-CA-C backbone
for a, b, color, lw in [(N, CA, '#333333', 3), (CA, C, '#333333', 3)]:
    ax1.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color=color, lw=lw, zorder=5)

# Draw C=O
ax1.plot([C[0], O[0]], [C[1], O[1]], [C[2], O[2]], color='#e74c3c', lw=2.5, zorder=5)

# Draw WT CB (L-Ala) — pointing "up"
ax1.plot([CA[0], CB[0]], [CA[1], CB[1]], [CA[2], CB[2]], color='#3498db', lw=2.5, zorder=5)
if HA is not None:
    ax1.plot([CA[0], HA[0]], [CA[1], HA[1]], [CA[2], HA[2]], color='#3498db', lw=1.5, ls='--', zorder=4)
ax1.scatter([CB[0]], [CB[1]], [CB[2]], color='#3498db', s=80, zorder=6)

# Draw D-Ala CB — pointing "down"
ax1.plot([CA[0], CB_D[0]], [CA[1], CB_D[1]], [CA[2], CB_D[2]], color='#e67e22', lw=2.5, zorder=5)
if HA_D is not None:
    ax1.plot([CA[0], HA_D[0]], [CA[1], HA_D[1]], [CA[2], HA_D[2]], color='#e67e22', lw=1.5, ls='--', zorder=4)
ax1.scatter([CB_D[0]], [CB_D[1]], [CB_D[2]], color='#e67e22', s=80, zorder=6)

# Scatter backbone atoms
for pos, label, color in [(N, 'N', '#2ecc71'), (CA, 'CA', '#f1c40f'), (C, 'C', '#9b59b6'), (O, 'O', '#e74c3c')]:
    ax1.scatter([pos[0]], [pos[1]], [pos[2]], color=color, s=120, edgecolors='black', linewidths=1, zorder=7)
    ax1.text(pos[0]+0.15, pos[1]+0.15, pos[2]+0.15, label, fontsize=10, fontweight='bold')

# Mirror plane (N-CA-C triangle, semi-transparent)
verts = [N, CA, C]
plane = Poly3DCollection([verts], alpha=0.15, facecolor='gray', edgecolor='black', linewidths=1)
ax1.add_collection3d(plane)

ax1.set_title('Ala2 Stereochemistry: L-Ala (blue) vs D-Ala (orange)\nMirror plane = N–CA–C', fontsize=11, fontweight='bold')
ax1.set_xlabel('X (Å)')
ax1.set_ylabel('Y (Å)')
ax1.set_zlabel('Z (Å)')

# Equal aspect ratio
max_range = 2.5
center = CA
ax1.set_xlim(center[0]-max_range, center[0]+max_range)
ax1.set_ylim(center[1]-max_range, center[1]+max_range)
ax1.set_zlim(center[2]-max_range, center[2]+max_range)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='#3498db', lw=2.5, label=f'L-Ala (WT)  χ={wt_dihedral:.1f}°'),
    Line2D([0], [0], color='#e67e22', lw=2.5, label=f'D-Ala (mut) χ={da_dihedral:.1f}°'),
    Line2D([0], [0], color='gray', lw=1, ls='--', label='Mirror plane (N–CA–C)'),
]
ax1.legend(handles=legend_elements, loc='upper left', fontsize=9)

# --- Right: 2D schematic (looking down CA→N bond, standard stereochemistry view) ---
ax2 = fig.add_subplot(122)

# Project onto plane perpendicular to N→CA for schematic
# Simpler: use standard Fischer-like projection coordinates
# Place CA at origin, C to the right, N up, H down, CB left (L) or right (D)

# L-Ala standard view (Fischer-like, NH2 on top, COOH on right)
# In this projection: N up, C right, H down, CB left
ax2.set_xlim(-3, 3)
ax2.set_ylim(-3, 3)
ax2.set_aspect('equal')
ax2.axis('off')

# Central CA
circle_ca = plt.Circle((0, 0), 0.25, color='#f1c40f', ec='black', linewidth=2, zorder=5)
ax2.add_patch(circle_ca)
ax2.text(0, 0, 'CA', ha='center', va='center', fontsize=9, fontweight='bold', zorder=6)

# N (top)
ax2.plot([0, 0], [0, 1.8], color='#333333', lw=2.5, zorder=4)
circle_n = plt.Circle((0, 1.8), 0.25, color='#2ecc71', ec='black', linewidth=2, zorder=5)
ax2.add_patch(circle_n)
ax2.text(0, 1.8, 'N', ha='center', va='center', fontsize=9, fontweight='bold', zorder=6)

# C (right)  
ax2.plot([0, 1.8], [0, 0], color='#333333', lw=2.5, zorder=4)
circle_c = plt.Circle((1.8, 0), 0.25, color='#9b59b6', ec='black', linewidth=2, zorder=5)
ax2.add_patch(circle_c)
ax2.text(1.8, 0, 'C', ha='center', va='center', fontsize=9, fontweight='bold', zorder=6)

# C=O (up-right from C)
ax2.plot([1.8, 2.4], [0, 0.6], color='#e74c3c', lw=2, zorder=4)
circle_o = plt.Circle((2.4, 0.6), 0.22, color='#e74c3c', ec='black', linewidth=2, zorder=5)
ax2.add_patch(circle_o)
ax2.text(2.4, 0.6, 'O', ha='center', va='center', fontsize=9, fontweight='bold', color='white', zorder=6)

# H (bottom)
ax2.plot([0, 0], [0, -1.5], color='#333333', lw=2, zorder=4)
ax2.text(0, -1.7, 'H', ha='center', va='center', fontsize=9, fontweight='bold', zorder=6)

# L-Ala CB (left) — dashed to indicate "away from viewer" in standard Fischer
ax2.plot([0, -1.8], [0, 0], color='#3498db', lw=2.5, ls='--', zorder=4)
ax2.text(-1.8, 0, 'CH₃', ha='center', va='center', fontsize=10, color='#3498db', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#3498db', linewidth=1.5))

# D-Ala CB (right) — wedge to indicate "toward viewer"
# Use a filled triangle wedge pointing right
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches

# Draw D-Ala CB as solid wedge (toward viewer)
# Triangle wedge: base at CA, pointing right
wedge = mpatches.FancyBboxPatch((0.15, -0.15), 1.5, 0.3, boxstyle="square,pad=0",
                                 facecolor='#e67e22', edgecolor='#e67e22', linewidth=0, zorder=3)
# Better: draw a filled polygon
poly = plt.Polygon([[0, -0.12], [1.7, -0.25], [1.7, 0.25], [0, 0.12]], 
                    closed=True, facecolor='#e67e22', edgecolor='#d35400', linewidth=1.5, zorder=3)
ax2.add_patch(poly)
ax2.text(1.85, 0, 'CH₃', ha='center', va='center', fontsize=10, color='#e67e22', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#e67e22', linewidth=1.5))

# Titles and annotations
ax2.set_title('Fischer Projection: Ala2 Stereochemistry\n(Wedge = toward viewer, Dashed = away)', 
              fontsize=11, fontweight='bold')

# Annotation boxes
ax2.text(-2.7, 2.5, 'L-Alanine (WT)\n  – CB away from viewer\n  – N-CA-CB-C χ ≈ –60°\n  → Ser630 can attack C=O', 
         fontsize=9, color='#3498db', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#ebf5fb', edgecolor='#3498db', linewidth=1.5))

ax2.text(0.3, -2.5, 'D-Alanine (mutant)\n  – CB toward viewer\n  – N-CA-CB-C χ ≈ +60°\n  → C=O flipped away from Ser630', 
         fontsize=9, color='#e67e22', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef5e7', edgecolor='#e67e22', linewidth=1.5))

# Dashed separator line between the two projections
ax2.axvline(x=0, color='gray', linestyle=':', linewidth=1, alpha=0.5, zorder=1)

plt.tight_layout()
outpath = '/home/scroll/personal/cjc-1295/workspace/figures/tripeptide_wt_vs_DAla_comparison.png'
plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Figure saved to: {outpath}")
plt.close()

# ============================================================
# Also generate a second figure: catalytic geometry comparison
# ============================================================
fig2, axes = plt.subplots(1, 2, figsize=(13, 5))

# Read key distances from WT tripeptide
ser630_OG = pdb.get((2, "ALA", "O"), None)  # Actually we need DPP-IV Ser630
# Wait, the tripeptide PDB only has the tripeptide. We need DPP-IV atoms too.
# Actually the tripeptide was docked into DPP-IV, so the PDB should have both.
# Let me check which atoms are in the file.

res2_atoms = [(k, v) for k, v in pdb.items() if k[0] == 2]
print("Residue 2 atoms:", [k[2] for k, _ in res2_atoms[:10]])

# Check for DPP-IV atoms
dpp4_atoms = [k for k in pdb.keys() if k[1] not in ("TYR", "ALA", "ASP", "NME")]
print(f"Non-tripeptide residues: {len(dpp4_atoms)}")
if dpp4_atoms:
    print("Sample:", dpp4_atoms[:5])
