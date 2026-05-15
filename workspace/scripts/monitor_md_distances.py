
# ⚠️ CHIRALITY NOTE (2026-05-15):
# "WT" trajectory = CJC-1295 (D-Ala2), "D-Ala2" trajectory = native GHRH (L-Ala)
# See docs/CHIRALITY_CORRECTION.md

#!/usr/bin/env python3
"""
Monitor key catalytic distances during GROMACS MD trajectory.
Extracts Ser630 OG→Ala2 C (attack distance) and oxyanion hole (Tyr547 OH→Ala2 O)
from trajectory frames and plots time series.

Usage:
    python monitor_md_distances.py -x md.xtc -s npt.gro -o distance_plot.png

Requires: MDAnalysis, numpy, matplotlib
"""

import argparse
import sys

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
except ImportError:
    print("ERROR: MDAnalysis not installed. Install with: pip install MDAnalysis")
    sys.exit(1)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def compute_distances(u, selections):
    """
    selections: list of (name, sel1_str, sel2_str) tuples
    Returns: dict of name -> list of distances per frame
    """
    results = {name: [] for name, _, _ in selections}
    times = []

    for ts in u.trajectory:
        times.append(ts.time)
        for name, s1, s2 in selections:
            g1 = u.select_atoms(s1)
            g2 = u.select_atoms(s2)
            if len(g1) == 0 or len(g2) == 0:
                results[name].append(np.nan)
                continue
            d = distances.distance_array(g1.positions, g2.positions)[0, 0]
            results[name].append(d)

    return np.array(times), results


def plot_distances(times, results, output_path):
    fig, axes = plt.subplots(len(results), 1, figsize=(10, 3 * len(results)), sharex=True)
    if len(results) == 1:
        axes = [axes]

    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']
    target_lines = [3.5, 3.5, 3.0, 3.0]  # approximate catalytic thresholds

    for ax, (name, dists), color, target in zip(axes, results.items(), colors, target_lines):
        dists = np.array(dists)
        ax.plot(times / 1000, dists, color=color, linewidth=0.8, alpha=0.8)
        ax.axhline(y=target, color='gray', linestyle='--', linewidth=1, alpha=0.5,
                   label=f'target ~{target} Å')
        ax.fill_between(times / 1000, dists, alpha=0.15, color=color)

        # Stats
        mean_d = np.nanmean(dists)
        std_d = np.nanstd(dists)
        ax.axhline(y=mean_d, color=color, linestyle='-', linewidth=1.5, alpha=0.7,
                   label=f'mean = {mean_d:.2f} ± {std_d:.2f} Å')

        ax.set_ylabel(f'{name}\n(Å)', fontsize=10)
        ax.legend(loc='upper right', fontsize=8)
        ax.set_ylim(0, max(8, np.nanmax(dists) * 1.1))

    axes[-1].set_xlabel('Time (ns)', fontsize=11)
    fig.suptitle('Catalytic Geometry Monitoring During MD', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Plot saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Monitor catalytic distances in MD')
    parser.add_argument('-x', '--xtc', required=True, help='Trajectory file (.xtc)')
    parser.add_argument('-s', '--tpr', required=True, help='Structure file (.tpr or .gro)')
    parser.add_argument('-o', '--output', default='md_distance_monitor.png', help='Output plot')
    args = parser.parse_args()

    print(f"Loading trajectory: {args.xtc} + {args.tpr}")
    u = mda.Universe(args.tpr, args.xtc)

    # Selections for WT (or D-Ala2 with DALA residue name)
    # Try DALA first, fallback to ALA
    selections = [
        ("Ser630 OG → Ala2 C (attack)",
         "resname SER and resid 630 and name OG",
         "(resname ALA or resname DALA) and resid 2 and name C"),
        ("Tyr547 OH → Ala2 O (oxyanion)",
         "resname TYR and resid 547 and name OH",
         "(resname ALA or resname DALA) and resid 2 and name O"),
        ("Tyr1 N → Glu206 Oε (salt bridge)",
         "resname TYR and resid 1 and name N",
         "resname GLU and resid 206 and name OE*"),
    ]

    print("Computing distances...")
    times, results = compute_distances(u, selections)

    # Print summary
    print("\n=== Distance Summary ===")
    for name, dists in results.items():
        d = np.array(dists)
        print(f"{name}: mean={np.nanmean(d):.3f} Å, std={np.nanstd(d):.3f} Å, "
              f"min={np.nanmin(d):.3f} Å, max={np.nanmax(d):.3f} Å")

    # Fraction of frames within catalytic range
    attack = np.array(results["Ser630 OG → Ala2 C (attack)"])
    within_range = np.sum((attack >= 2.5) & (attack <= 4.0)) / len(attack) * 100
    print(f"\nFrames with attack distance 2.5–4.0 Å: {within_range:.1f}%")

    plot_distances(times, results, args.output)


if __name__ == '__main__':
    main()
