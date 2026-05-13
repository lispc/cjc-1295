#!/usr/bin/env python3
"""
Analyze catalytic geometry stability over the MD trajectory.
Extract key distances and angles from the trajectory and plot their time evolution.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

STEP3 = "/home/scroll/personal/cjc-1295/workspace/step3"
TPR = f"{STEP3}/md.tpr"
XTC = f"{STEP3}/md.xtc"
OUTDIR = f"{STEP3}/analysis"

def run_gmx(command, stdin=""):
    """Run a gmx command and return stdout."""
    env = os.environ.copy()
    env["GMX_MAXBACKUP"] = "-1"
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True,
        input=stdin, env=env
    )
    return result.stdout, result.stderr, result.returncode

def extract_distance(traj, tpr, selection1, selection2, outfile):
    """Extract distance between two atom selections using gmx distance."""
    # Create index group
    idx_content = f"[sel1]\n{selection1}\n[sel2]\n{selection2}\n"
    with open(f"{outfile}.ndx", "w") as f:
        f.write(idx_content)
    
    stdout, stderr, rc = run_gmx(
        f"gmx distance -s {tpr} -f {traj} -n {outfile}.ndx -oall {outfile}.xvg -select 'group \"sel1\" plus group \"sel2\"' -tu ns"
    )
    if rc != 0:
        print(f"Warning: distance extraction failed: {stderr[:200]}")
        return None
    
    # Parse xvg
    data = []
    with open(f"{outfile}.xvg") as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

def extract_angle(traj, tpr, sel1, sel2, sel3, outfile):
    """Extract angle between three atom selections."""
    idx_content = f"[a1]\n{sel1}\n[a2]\n{sel2}\n[a3]\n{sel3}\n"
    with open(f"{outfile}.ndx", "w") as f:
        f.write(idx_content)
    
    stdout, stderr, rc = run_gmx(
        f"gmx angle -s {tpr} -f {traj} -n {outfile}.ndx -ov {outfile}.xvg -type angle -tu ns"
    )
    if rc != 0:
        print(f"Warning: angle extraction failed: {stderr[:200]}")
        return None
    
    data = []
    with open(f"{outfile}.xvg") as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

def read_xvg(path):
    """Read a simple xvg file."""
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    
    # Check if trajectory exists
    if not os.path.exists(XTC):
        print(f"Trajectory not found: {XTC}")
        print("This script should be run after the production MD completes.")
        return
    
    print("=" * 60)
    print("Catalytic Geometry Analysis")
    print("=" * 60)
    
    # Atom selections (using GROMACS index syntax)
    # Chain A = DPP-IV, Chain B = GHRH
    # Note: these are 1-based residue numbers as in the PDB
    selections = {
        "Ser630_OG": "2",   # Will need actual atom index — better to use make_ndx
        # Actually gmx distance/angle with make_ndx is more robust
    }
    
    # Use make_ndx to create proper index groups
    print("\nCreating index groups...")
    
    # For now, let's use residue-based selections with make_ndx
    # We need to know the actual atom numbers in the tpr
    # Let's use a simpler approach: gmx select or awk-based index
    
    # Create index file with the atoms we need
    ndx_script = """
a OG & r 630 & chain A
a C & r 2 & chain B
a N & r 2 & chain B
a O & r 2 & chain B
a N & r 1 & chain B
a OE* & r 206 & chain A
a OH & r 547 & chain A
q
"""
    with open(f"{OUTDIR}/make_ndx_input.txt", "w") as f:
        f.write(ndx_script)
    
    stdout, stderr, rc = run_gmx(
        f"gmx make_ndx -f {TPR} -o {OUTDIR}/catalytic.ndx",
        stdin=ndx_script
    )
    if rc != 0:
        print(f"make_ndx failed: {stderr}")
        return
    
    # The groups will be numbered. Let's check what numbers they got.
    print(stdout[-500:] if len(stdout) > 500 else stdout)
    
    print(f"\nIndex file created: {OUTDIR}/catalytic.ndx")
    print("\nTo extract distances/angles after MD completes, run:")
    print(f"  gmx distance -s {TPR} -f {XTC} -n {OUTDIR}/catalytic.ndx ...")

if __name__ == "__main__":
    main()
