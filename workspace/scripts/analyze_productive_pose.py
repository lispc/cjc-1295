
# ⚠️ CHIRALITY NOTE (2026-05-15):
# "WT" = CJC-1295 (D-Ala2), "D-Ala2" = native GHRH (L-Ala)
# See docs/CHIRALITY_CORRECTION.md

#!/usr/bin/env python3
"""
P0a: Analyze productive pose fraction from existing MD distance + angle data.

Computes:
  1. Productive pose fraction (2D): d ∈ [2.5, 4.0] Å AND θ ∈ [80°, 120°]
  2. First passage time (FPT): time until first exit from productive zone
  3. Dwell time distribution: histogram of contiguous productive-state durations
  4. Block-averaged mean ± SE for all metrics

Uses existing .xvg files from gmx distance / gmx angle.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STEP3 = "/home/scroll/personal/cjc-1295/workspace/step3"
OUTDIR = "/home/scroll/personal/cjc-1295/workspace/results"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# XVG parser
# ---------------------------------------------------------------------------

def read_xvg(path):
    """Read a GROMACS xvg file, return (times, values) as numpy arrays."""
    t, v = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t.append(float(parts[0]))
                v.append(float(parts[1]))
    return np.array(t), np.array(v)


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

D_MIN, D_MAX = 0.25, 0.40       # nm = 2.5–4.0 Å
THETA_MIN, THETA_MAX = 80, 120  # degrees


def merge_xvgs(paths):
    """Read multiple xvg files and concatenate, filtering out non-existent."""
    all_t, all_v = [], []
    for p in paths:
        if os.path.exists(p):
            t, v = read_xvg(p)
            if len(t) > 0:
                all_t.append(t)
                all_v.append(v)
    if not all_t:
        return np.array([]), np.array([])
    t_cat = np.concatenate(all_t)
    v_cat = np.concatenate(all_v)
    idx = np.argsort(t_cat)
    return t_cat[idx], v_cat[idx]


def align_by_time(t1, v1, t2, v2):
    """Align two time series by nearest timestamp. Returns aligned arrays."""
    if len(t1) == 0 or len(t2) == 0:
        return np.array([]), np.array([]), np.array([])
    # Use t1 as reference, find nearest in t2
    v2_aligned = np.array([v2[np.argmin(np.abs(t2 - ti))] for ti in t1])
    return t1, v1, v2_aligned


def productive_mask(d, theta):
    """Boolean mask: True where both distance and angle are in productive zone."""
    d_ok = (d >= D_MIN) & (d <= D_MAX)
    theta_ok = (theta >= THETA_MIN) & (theta <= THETA_MAX)
    return d_ok & theta_ok


def compute_fpt(t, mask):
    """First passage time: time until first exit from productive zone."""
    if mask.all():
        return t[-1] - t[0]
    if mask[0]:
        first_exit = np.where(~mask)[0][0]
        return t[first_exit] - t[0]
    return 0.0


def dwell_times(t, mask):
    """Return list of contiguous dwell times inside productive zone (ps)."""
    if not mask.any():
        return []
    runs = []
    in_run = False
    start = 0
    for i, m in enumerate(mask):
        if m and not in_run:
            in_run = True
            start = i
        elif not m and in_run:
            in_run = False
            runs.append(t[i-1] - t[start])
    if in_run:
        runs.append(t[-1] - t[start])
    return runs


def block_average(data, n_blocks=5):
    """Split data into n_blocks, return mean of block means and SE."""
    n = len(data)
    if n < n_blocks:
        return np.mean(data), np.std(data, ddof=1) / np.sqrt(n) if n > 1 else 0.0
    block_size = n // n_blocks
    blocks = [data[i*block_size:(i+1)*block_size] for i in range(n_blocks)]
    block_means = [np.mean(b) for b in blocks]
    return np.mean(block_means), np.std(block_means, ddof=1) / np.sqrt(n_blocks)


# ---------------------------------------------------------------------------
# System definitions
# ---------------------------------------------------------------------------

SYSTEMS = {
    "WT_full": {
        "label": "WT GHRH(1-29)",
        "d_files": [
            f"{STEP3}/dist_wt_d1_part1.xvg",
            f"{STEP3}/dist_wt_p1_d1.xvg",
            f"{STEP3}/dist_wt_p2_d1.xvg",
        ],
        "a_files": [
            f"{STEP3}/ang_wt_p1.xvg",
            f"{STEP3}/ang_wt_p2.xvg",
        ],
        "color": "tab:blue",
    },
    "DAla2_full": {
        "label": "D-Ala2 GHRH(1-29)",
        "d_files": [f"{STEP3}/dist_dala2_d1.xvg"],
        "a_files": [f"{STEP3}/ang_dala2.xvg"],
        "color": "tab:orange",
    },
    "POSRES": {
        "label": "WT + POSRES (1-5 CA)",
        "d_files": [f"{STEP3}/dist_posres_d1.xvg"],
        "a_files": [f"{STEP3}/ang_posres.xvg"],
        "color": "tab:green",
    },
    "Short_WT": {
        "label": "WT GHRH(1-10)",
        "d_files": [f"{STEP3}/dist_short_d1.xvg"],
        "a_files": [f"{STEP3}/ang_short.xvg"],
        "color": "tab:purple",
    },
}


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main():
    results = {}

    print("=" * 70)
    print("P0a: Productive Pose Analysis (2D: distance + angle)")
    print("=" * 70)
    print(f"Productive zone: d ∈ [{D_MIN*10:.1f}, {D_MAX*10:.1f}] Å AND θ ∈ [{THETA_MIN}°, {THETA_MAX}°]")
    print()

    for sys_name, cfg in SYSTEMS.items():
        print(f">>> {cfg['label']}")

        t_d, d = merge_xvgs(cfg["d_files"])
        t_a, a = merge_xvgs(cfg["a_files"])

        if len(t_d) == 0:
            print(f"    WARNING: no distance data")
            continue
        if len(t_a) == 0:
            print(f"    WARNING: no angle data")
            continue

        # Align by time
        t, d_aligned, a_aligned = align_by_time(t_d, d, t_a, a)

        if len(t) == 0:
            print(f"    WARNING: could not align distance and angle data")
            continue

        # 2D productive pose fraction
        mask = productive_mask(d_aligned, a_aligned)
        frac = mask.mean() * 100

        # First passage time
        fpt = compute_fpt(t, mask)

        # Dwell times
        dwells = dwell_times(t, mask)
        n_dwells = len(dwells)
        max_dwell = max(dwells) / 1000.0 if dwells else 0.0
        total_dwell = sum(dwells) / 1000.0 if dwells else 0.0

        # Block averages
        d_ang = d_aligned * 10
        mean_d, se_d = block_average(d_ang, n_blocks=5)
        mean_a, se_a = block_average(a_aligned, n_blocks=5)

        # Distance-only and angle-only fractions for comparison
        d_only = ((d_aligned >= D_MIN) & (d_aligned <= D_MAX)).mean() * 100
        a_only = ((a_aligned >= THETA_MIN) & (a_aligned <= THETA_MAX)).mean() * 100

        results[sys_name] = {
            "label": cfg["label"],
            "color": cfg["color"],
            "t": t,
            "d": d_ang,
            "a": a_aligned,
            "frac_2d": frac,
            "frac_d": d_only,
            "frac_a": a_only,
            "fpt_ns": fpt / 1000.0,
            "n_dwells": n_dwells,
            "max_dwell_ns": max_dwell,
            "total_dwell_ns": total_dwell,
            "mean_d": mean_d,
            "se_d": se_d,
            "mean_a": mean_a,
            "se_a": se_a,
        }

        print(f"    Frames aligned: {len(t)}")
        print(f"    Time span: {t[0]/1000:.1f} – {t[-1]/1000:.1f} ns")
        print(f"    Distance-only productive: {d_only:.2f}%")
        print(f"    Angle-only productive:    {a_only:.2f}%")
        print(f"    2D productive fraction:   {frac:.2f}%")
        print(f"    First passage time:       {fpt/1000:.2f} ns")
        print(f"    Dwell events: {n_dwells}, max: {max_dwell:.2f} ns, total: {total_dwell:.2f} ns")
        print(f"    Ser630→Ala2C: {mean_d:.2f} ± {se_d:.2f} Å")
        print(f"    Attack angle:  {mean_a:.2f} ± {se_a:.2f}°")
        print()

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'System':<30} {'d (Å)':<14} {'θ (°)':<12} {'2D-Frac':<10} {'FPT(ns)':<8} {'Dwells':<8}")
    print("-" * 70)
    for k, r in results.items():
        print(f"{r['label']:<30} {r['mean_d']:.2f}±{r['se_d']:.2f}  {r['mean_a']:.1f}±{r['se_a']:.1f}   {r['frac_2d']:.2f}%     {r['fpt_ns']:.2f}     {r['n_dwells']}")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # (a) Distance time series
    ax = fig.add_subplot(gs[0, 0])
    for k, r in results.items():
        ax.plot(r["t"] / 1000.0, r["d"], label=r["label"], color=r["color"], alpha=0.8, lw=0.8)
    ax.axhspan(D_MIN*10, D_MAX*10, color="green", alpha=0.1)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ser630 OG → Ala2 C (Å)")
    ax.set_title("(a) Catalytic Distance")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(0, 12)

    # (b) Angle time series
    ax = fig.add_subplot(gs[0, 1])
    for k, r in results.items():
        ax.plot(r["t"] / 1000.0, r["a"], label=r["label"], color=r["color"], alpha=0.8, lw=0.8)
    ax.axhspan(THETA_MIN, THETA_MAX, color="green", alpha=0.1)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Attack angle (°)")
    ax.set_title("(b) Attack Angle ∠(OG–C–N)")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(0, 180)

    # (c) 2D productive fraction
    ax = fig.add_subplot(gs[0, 2])
    labels = [r["label"] for r in results.values()]
    fracs = [r["frac_2d"] for r in results.values()]
    colors = [r["color"] for r in results.values()]
    bars = ax.bar(range(len(labels)), fracs, color=colors, edgecolor="black")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("Productive Fraction (%)")
    ax.set_title("(c) 2D Productive Pose Fraction")
    ax.set_ylim(0, 105)
    for bar, frac in zip(bars, fracs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f"{frac:.1f}%", ha="center", va="bottom", fontsize=9)

    # (d) First passage time
    ax = fig.add_subplot(gs[1, 0])
    fpts = [r["fpt_ns"] for r in results.values()]
    bars = ax.bar(range(len(labels)), fpts, color=colors, edgecolor="black")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("First Passage Time (ns)")
    ax.set_title("(d) Time Until First Exit from 2D Zone")
    for bar, fpt in zip(bars, fpts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f"{fpt:.2f}", ha="center", va="bottom", fontsize=9)

    # (e) Mean distance with error bars
    ax = fig.add_subplot(gs[1, 1])
    means = [r["mean_d"] for r in results.values()]
    ses = [r["se_d"] for r in results.values()]
    x = np.arange(len(labels))
    ax.errorbar(x, means, yerr=ses, fmt="o", markersize=10, capsize=5, color="black", ecolor="black")
    for xi, mi, ci in zip(x, means, colors):
        ax.scatter(xi, mi, color=ci, s=150, zorder=3, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.axhspan(D_MIN*10, D_MAX*10, color="green", alpha=0.1)
    ax.set_ylabel("Mean Ser630→Ala2C (Å)")
    ax.set_title("(e) Block-Averaged Distance ± SE")
    ax.set_ylim(0, 12)

    # (f) Mean angle with error bars
    ax = fig.add_subplot(gs[1, 2])
    means_a = [r["mean_a"] for r in results.values()]
    ses_a = [r["se_a"] for r in results.values()]
    ax.errorbar(x, means_a, yerr=ses_a, fmt="o", markersize=10, capsize=5, color="black", ecolor="black")
    for xi, mi, ci in zip(x, means_a, colors):
        ax.scatter(xi, mi, color=ci, s=150, zorder=3, edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.axhspan(THETA_MIN, THETA_MAX, color="green", alpha=0.1)
    ax.set_ylabel("Mean Attack Angle (°)")
    ax.set_title("(f) Block-Averaged Angle ± SE")
    ax.set_ylim(0, 180)

    # (g) 2D scatter: distance vs angle (WT)
    ax = fig.add_subplot(gs[2, 0])
    if "WT_full" in results:
        r = results["WT_full"]
        mask = productive_mask(r["d"]/10, r["a"])
        ax.scatter(r["d"][~mask], r["a"][~mask], c="lightgray", s=5, alpha=0.5, label="Non-productive")
        ax.scatter(r["d"][mask], r["a"][mask], c="tab:blue", s=10, alpha=0.8, label="Productive")
        ax.axvspan(D_MIN*10, D_MAX*10, color="green", alpha=0.05)
        ax.axhspan(THETA_MIN, THETA_MAX, color="green", alpha=0.05)
        ax.set_xlabel("Ser630→Ala2C (Å)")
        ax.set_ylabel("Attack angle (°)")
        ax.set_title("(g) WT: Distance vs Angle")
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 180)
        ax.legend(fontsize=7)

    # (h) 2D scatter: distance vs angle (D-Ala2)
    ax = fig.add_subplot(gs[2, 1])
    if "DAla2_full" in results:
        r = results["DAla2_full"]
        mask = productive_mask(r["d"]/10, r["a"])
        ax.scatter(r["d"][~mask], r["a"][~mask], c="lightgray", s=5, alpha=0.5, label="Non-productive")
        ax.scatter(r["d"][mask], r["a"][mask], c="tab:orange", s=10, alpha=0.8, label="Productive")
        ax.axvspan(D_MIN*10, D_MAX*10, color="green", alpha=0.05)
        ax.axhspan(THETA_MIN, THETA_MAX, color="green", alpha=0.05)
        ax.set_xlabel("Ser630→Ala2C (Å)")
        ax.set_ylabel("Attack angle (°)")
        ax.set_title("(h) D-Ala2: Distance vs Angle")
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 180)
        ax.legend(fontsize=7)

    # (i) Fraction breakdown
    ax = fig.add_subplot(gs[2, 2])
    x = np.arange(len(labels))
    width = 0.25
    fracs_d = [r["frac_d"] for r in results.values()]
    fracs_a = [r["frac_a"] for r in results.values()]
    fracs_2d = [r["frac_2d"] for r in results.values()]
    ax.bar(x - width, fracs_d, width, label="Distance only", color="tab:cyan", edgecolor="black")
    ax.bar(x, fracs_a, width, label="Angle only", color="tab:olive", edgecolor="black")
    ax.bar(x + width, fracs_2d, width, label="2D (both)", color="tab:red", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right", fontsize=8)
    ax.set_ylabel("Fraction (%)")
    ax.set_title("(i) Productive Fraction Breakdown")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 105)

    plt.savefig(f"{OUTDIR}/productive_pose_analysis_2d.png", dpi=200, bbox_inches="tight")
    print(f"\nFigure saved: {OUTDIR}/productive_pose_analysis_2d.png")

    # -----------------------------------------------------------------------
    # Save numeric results to text
    # -----------------------------------------------------------------------
    txtpath = f"{OUTDIR}/productive_pose_analysis_2d.txt"
    with open(txtpath, "w") as f:
        f.write("P0a: Productive Pose Analysis Results (2D criteria)\n")
        f.write("=" * 70 + "\n\n")
        f.write("Productive pose defined as:\n")
        f.write(f"  d_Ser630-OG → Ala2-C ∈ [{D_MIN*10:.1f}, {D_MAX*10:.1f}] Å\n")
        f.write(f"  θ_attack = ∠(OG–C–N) ∈ [{THETA_MIN}°, {THETA_MAX}°]\n\n")
        for k, r in results.items():
            f.write(f"{r['label']}\n")
            f.write(f"  Time span: {r['t'][0]/1000:.1f} – {r['t'][-1]/1000:.1f} ns  ({len(r['t'])} frames)\n")
            f.write(f"  Ser630→Ala2C: {r['mean_d']:.2f} ± {r['se_d']:.2f} Å\n")
            f.write(f"  Attack angle:  {r['mean_a']:.1f} ± {r['se_a']:.1f}°\n")
            f.write(f"  Distance-only productive: {r['frac_d']:.2f}%\n")
            f.write(f"  Angle-only productive:    {r['frac_a']:.2f}%\n")
            f.write(f"  2D productive fraction:   {r['frac_2d']:.2f}%\n")
            f.write(f"  First passage time:       {r['fpt_ns']:.2f} ns\n")
            f.write(f"  Dwell events: {r['n_dwells']}, max: {r['max_dwell_ns']:.2f} ns, total: {r['total_dwell_ns']:.2f} ns\n\n")
    print(f"Text report saved: {txtpath}")


if __name__ == "__main__":
    main()
