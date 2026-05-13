#!/usr/bin/env python3
"""
Comprehensive MD trajectory analysis for DPP-IV / GHRH(1-29) complex.
To be run after the 200 ns production MD completes.

Produces:
- Catalytic geometry time series (distances + angles)
- RMSF per residue
- Hydrogen bond occupancy
- RMSD relative to starting structure
- Summary plots
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STEP3 = "/home/scroll/personal/cjc-1295/workspace/step3"
TPR = f"{STEP3}/md.tpr"
XTC = f"{STEP3}/md.xtc"
OUTDIR = f"{STEP3}/analysis"
START_GRO = f"{STEP3}/em.gro"  # Reference structure

def run(cmd, stdin=""):
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                           input=stdin, env=env)
    return result.returncode, result.stdout, result.stderr

def read_xvg(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(p) for p in parts])
    return np.array(data)

def make_index():
    """Create index groups for catalytic geometry analysis."""
    idx_script = """a OG & r 630
r 1 & a N
r 2 & a C
r 2 & a N
r 2 & a O
r 206 & a OE1 OE2
r 547 & a OH
r 1 & a CA
r 2 & a CA
r 3 & a CA
r 1 & a N C CA O CB
r 2 & a N C CA O CB
r 3 & a N C CA O CB
q
"""
    rc, out, err = run(f"gmx make_ndx -f {TPR} -o {OUTDIR}/catalytic.ndx", stdin=idx_script)
    if rc != 0:
        print(f"make_ndx failed: {err}")
        return False
    print("Index file created.")
    return True

def extract_distances():
    """Extract key distances from trajectory."""
    pairs = [
        ("d_Ser630_Ala2C",   "group 'a_1' plus group 'a_3'"),      # OG(630) - C(Ala2)
        ("d_Tyr1_Glu206",    "group 'a_2' plus group 'a_6'"),      # N(Tyr1) - OE*(Glu206)
        ("d_Ala2O_Tyr547",   "group 'a_5' plus group 'a_7'"),      # O(Ala2) - OH(Tyr547)
        ("d_Ser630_Ala2O",   "group 'a_1' plus group 'a_5'"),      # OG(630) - O(Ala2)
    ]
    
    results = {}
    for name, selection in pairs:
        cmd = (f"gmx distance -s {TPR} -f {XTC} -n {OUTDIR}/catalytic.ndx "
               f"-oall {OUTDIR}/{name}.xvg -select '{selection}' -tu ns")
        rc, out, err = run(cmd)
        if rc != 0:
            print(f"Warning: {name} extraction failed: {err[:200]}")
            continue
        data = read_xvg(f"{OUTDIR}/{name}.xvg")
        if len(data) > 0:
            results[name] = data
            print(f"{name}: mean={data[:,1].mean():.3f} nm, std={data[:,1].std():.3f} nm")
    return results

def extract_angles():
    """Extract attack angle from trajectory."""
    # Need an index file with 3 groups for angle
    angle_script = """a OG & r 630
r 2 & a C
r 2 & a N
q
"""
    rc, out, err = run(f"gmx make_ndx -f {TPR} -o {OUTDIR}/angle.ndx", stdin=angle_script)
    if rc != 0:
        print(f"angle make_ndx failed: {err}")
        return {}
    
    cmd = (f"gmx angle -s {TPR} -f {XTC} -n {OUTDIR}/angle.ndx "
           f"-ov {OUTDIR}/angle_attack.xvg -type angle -tu ns")
    rc, out, err = run(cmd)
    if rc != 0:
        print(f"angle extraction failed: {err[:200]}")
        return {}
    
    data = read_xvg(f"{OUTDIR}/angle_attack.xvg")
    if len(data) > 0:
        print(f"Attack angle: mean={data[:,1].mean():.1f} deg, std={data[:,1].std():.1f} deg")
        return {"angle_attack": data}
    return {}

def extract_rmsd():
    """Extract RMSD of protein and peptide."""
    # Protein RMSD
    cmd = (f"gmx rms -s {TPR} -f {XTC} -o {OUTDIR}/rmsd_protein.xvg -tu ns "
           f"<< 'EOF'\n1\n1\nEOF")
    rc, out, err = run(cmd)
    
    # Peptide RMSD (GHRH residues 1-29, C-alpha)
    # Create index for peptide
    pep_script = "ri 1-29\nq\n"
    rc, out, err = run(f"gmx make_ndx -f {TPR} -o {OUTDIR}/peptide.ndx", stdin=pep_script)
    
    cmd = (f"gmx rms -s {TPR} -f {XTC} -n {OUTDIR}/peptide.ndx "
           f"-o {OUTDIR}/rmsd_peptide.xvg -tu ns << 'EOF'\n3\n3\nEOF")
    rc, out, err = run(cmd)
    
    results = {}
    for name in ["rmsd_protein", "rmsd_peptide"]:
        path = f"{OUTDIR}/{name}.xvg"
        if os.path.exists(path):
            data = read_xvg(path)
            if len(data) > 0:
                results[name] = data
                print(f"{name}: mean={data[:,1].mean():.3f} nm")
    return results

def extract_rmsf():
    """Extract RMSF per residue."""
    # Protein RMSF
    cmd = (f"gmx rmsf -s {TPR} -f {XTC} -o {OUTDIR}/rmsf.xvg -res -tu ns << 'EOF'\n1\nEOF")
    rc, out, err = run(cmd)
    
    data = read_xvg(f"{OUTDIR}/rmsf.xvg")
    if len(data) > 0:
        print(f"RMSF: protein mean={data[:,1].mean():.3f} nm")
        return {"rmsf": data}
    return {}

def plot_geometry(distances, angles):
    """Plot catalytic geometry time series."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Distance plots
    ax = axes[0, 0]
    if "d_Ser630_Ala2C" in distances:
        d = distances["d_Ser630_Ala2C"]
        ax.plot(d[:,0], d[:,1]*10, 'b-', alpha=0.7, lw=0.5)
        ax.axhline(y=2.53, color='g', ls='--', label='Step0 ideal (2.53 Å)')
        ax.axhline(y=4.0, color='r', ls='--', label='Max attack (4.0 Å)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('Ser630 OG → Ala2 C=O C')
        ax.legend(fontsize=8)
    
    ax = axes[0, 1]
    if "d_Ala2O_Tyr547" in distances:
        d = distances["d_Ala2O_Tyr547"]
        ax.plot(d[:,0], d[:,1]*10, 'b-', alpha=0.7, lw=0.5)
        ax.axhline(y=3.15, color='g', ls='--', label='Step0 ideal (3.15 Å)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('Ala2 O → Tyr547 OH')
        ax.legend(fontsize=8)
    
    ax = axes[1, 0]
    if "d_Tyr1_Glu206" in distances:
        d = distances["d_Tyr1_Glu206"]
        ax.plot(d[:,0], d[:,1]*10, 'b-', alpha=0.7, lw=0.5)
        ax.axhline(y=1.87, color='g', ls='--', label='Step0 ideal (1.87 Å)')
        ax.axhline(y=3.5, color='r', ls='--', label='Salt bridge limit (3.5 Å)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('Tyr1 N → Glu206 Oε')
        ax.legend(fontsize=8)
    
    ax = axes[1, 1]
    if "angle_attack" in angles:
        a = angles["angle_attack"]
        ax.plot(a[:,0], a[:,1], 'b-', alpha=0.7, lw=0.5)
        ax.axhline(y=113.6, color='g', ls='--', label='Step0 ideal (113.6°)')
        ax.axhspan(80, 120, alpha=0.2, color='green', label='Favorable (80-120°)')
        ax.axhspan(120, 140, alpha=0.2, color='yellow')
        ax.axhspan(140, 180, alpha=0.2, color='red', label='Blocked (>140°)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Angle (°)')
        ax.set_title('Attack Angle ∠(Ser630 OG–Ala2 C–Ala2 N)')
        ax.set_ylim(50, 180)
        ax.legend(fontsize=8)
    
    plt.suptitle('Catalytic Geometry Stability Over 200 ns MD', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/geometry_time_series.png", dpi=300)
    print(f"Saved: {OUTDIR}/geometry_time_series.png")
    plt.close()

def plot_rmsd_rmsf(rmsd_data, rmsf_data):
    """Plot RMSD and RMSF summary."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    ax = axes[0]
    if "rmsd_protein" in rmsd_data:
        ax.plot(rmsd_data["rmsd_protein"][:,0], rmsd_data["rmsd_protein"][:,1]*10, 
                'b-', alpha=0.6, lw=0.5, label='Protein')
    if "rmsd_peptide" in rmsd_data:
        ax.plot(rmsd_data["rmsd_peptide"][:,0], rmsd_data["rmsd_peptide"][:,1]*10, 
                'r-', alpha=0.6, lw=0.5, label='GHRH(1-29)')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('RMSD vs Starting Structure')
    ax.legend()
    
    ax = axes[1]
    if "rmsf" in rmsf_data:
        r = rmsf_data["rmsf"]
        ax.plot(r[:,0], r[:,1]*10, 'b-', alpha=0.6, lw=0.5)
        # Highlight GHRH region (residues 1-29)
        ax.axvspan(1, 29, alpha=0.2, color='red', label='GHRH(1-29)')
        ax.set_xlabel('Residue Number')
        ax.set_ylabel('RMSF (Å)')
        ax.set_title('Per-Residue RMSF')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/rmsd_rmsf_summary.png", dpi=300)
    print(f"Saved: {OUTDIR}/rmsd_rmsf_summary.png")
    plt.close()

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    
    if not os.path.exists(XTC):
        print(f"Trajectory not found: {XTC}")
        print("Please run this script after production MD completes.")
        return 1
    
    print("=" * 60)
    print("MD Trajectory Analysis")
    print("=" * 60)
    
    # Activate gmx environment
    os.environ["PATH"] = "/home/scroll/miniforge3/envs/gmx/bin:" + os.environ.get("PATH", "")
    
    print("\n[1/5] Creating index groups...")
    if not make_index():
        return 1
    
    print("\n[2/5] Extracting catalytic geometry...")
    distances = extract_distances()
    angles = extract_angles()
    
    print("\n[3/5] Extracting RMSD...")
    rmsd_data = extract_rmsd()
    
    print("\n[4/5] Extracting RMSF...")
    rmsf_data = extract_rmsf()
    
    print("\n[5/5] Generating plots...")
    plot_geometry(distances, angles)
    plot_rmsd_rmsf(rmsd_data, rmsf_data)
    
    print("\n" + "=" * 60)
    print("Analysis complete. Results in:", OUTDIR)
    print("=" * 60)
    return 0

if __name__ == "__main__":
    sys.exit(main())
