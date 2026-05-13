#!/usr/bin/env python3
"""
Run APBS electrostatics calculation for DPP-IV with YAD overlay.
Workflow:
1. Use pdb2pqr to prepare PQR file for DPP-IV
2. Run APBS to calculate electrostatic potential
3. Use PyMOL to generate visualization
"""

import os, subprocess, sys

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP0 = f"{BASE}/step0"
STEP1 = f"{BASE}/step1"
RESULTS = f"{BASE}/results"

os.makedirs(RESULTS, exist_ok=True)

# Paths
PDB2PQR = "/home/scroll/miniforge3/bin/pdb2pqr"
APBS = "/home/scroll/miniforge3/bin/apbs"

# Step 1: Run pdb2pqr on DPP4_clean.pdb
print("Step 1: Running pdb2pqr on DPP4_clean.pdb...")
cmd_pdb2pqr = [
    PDB2PQR,
    "--ff=AMBER",
    "--keep-chain",
    f"{STEP1}/DPP4_clean.pdb",
    f"{STEP0}/DPP4_clean.pqr"
]
result = subprocess.run(cmd_pdb2pqr, capture_output=True, text=True)
if result.returncode != 0:
    print(f"pdb2pqr failed: {result.stderr}")
    sys.exit(1)
print("  pdb2pqr completed.")

# Step 2: Determine grid dimensions from PQR file
# Read PQR to get min/max coords
atoms = []
with open(f"{STEP0}/DPP4_clean.pqr") as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((x,y,z))

if not atoms:
    print("No atoms found in PQR file!")
    sys.exit(1)

xs = [a[0] for a in atoms]
ys = [a[1] for a in atoms]
zs = [a[2] for a in atoms]

minx, maxx = min(xs), max(xs)
miny, maxy = min(ys), max(ys)
minz, maxz = min(zs), max(zs)

# Center
cx = (minx + maxx) / 2
cy = (miny + maxy) / 2
cz = (minz + maxz) / 2

# Grid size (add padding)
padding = 20.0  # Angstroms
lenx = (maxx - minx) + 2*padding
leny = (maxy - miny) + 2*padding
lenz = (maxz - minz) + 2*padding

# Grid dimensions (must be odd for multigrid, but for mg-manual we can use any)
# Let's use reasonable grid spacing (0.5 A) and dimensions
grid_spacing = 0.5
npx = int(lenx / grid_spacing) + 1
npy = int(leny / grid_spacing) + 1
npz = int(lenz / grid_spacing) + 1

# Make sure dimensions are odd for mg-auto
if npx % 2 == 0: npx += 1
if npy % 2 == 0: npy += 1
if npz % 2 == 0: npz += 1

# For mg-auto, we need to specify cglen (coarse grid) and fglen (fine grid)
# Use a smaller fine grid around the protein
cglen = (lenx, leny, lenz)
fglen = (lenx*0.5, leny*0.5, lenz*0.5)

print(f"  Protein dimensions: {maxx-minx:.1f} x {maxy-miny:.1f} x {maxz-minz:.1f} Å")
print(f"  Grid center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")
print(f"  Grid lengths: cglen={cglen}, fglen={fglen}")

# Step 3: Write APBS input file
apbs_in = f"""read
    mol pqr {STEP0}/DPP4_clean.pqr
end

elec
    mg-auto
    dime {npx} {npy} {npz}
    cglen {cglen[0]:.2f} {cglen[1]:.2f} {cglen[2]:.2f}
    fglen {fglen[0]:.2f} {fglen[1]:.2f} {fglen[2]:.2f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    npbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    srfm mol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 310.0
    calcenergy total
    calcforce no
    write pot dx {STEP0}/DPP4_potential
end

quit
"""

with open(f"{STEP0}/apbs.in", "w") as f:
    f.write(apbs_in)

print("Step 2: Running APBS...")
result = subprocess.run([APBS, f"{STEP0}/apbs.in"], capture_output=True, text=True)
if result.returncode != 0:
    print(f"APBS failed: {result.stderr}")
    print(f"APBS stdout: {result.stdout[-500:]}")
    sys.exit(1)
print("  APBS completed.")

# Step 4: Generate PyMOL visualization script
pymol_script = f"""
load {STEP1}/DPP4_clean.pdb, DPP4
load {STEP0}/DPP4_potential.dx, pot
load {STEP0}/YAD_aligned_to_DPP4.pdb, YAD

# Color by chain for clarity
as cartoon, DPP4
show surface, DPP4
set transparency, 0.3, DPP4

# Set electrostatic potential colors
ramp_new esp, pot, [-5, 0, 5], [blue, white, red]
color esp, DPP4

# Show YAD as sticks
color yellow, YAD and resi 2
show sticks, YAD

# Save session
save {RESULTS}/step0_electrostatics.pse

# Save image
ray 2400, 2400
png {RESULTS}/step0_electrostatics.png, dpi=300

quit
"""

with open(f"{STEP0}/visualize_apbs.pml", "w") as f:
    f.write(pymol_script)

print("Step 3: Generating PyMOL visualization...")
result = subprocess.run(
    ["/home/scroll/miniforge3/bin/pymol", "-cq", f"{STEP0}/visualize_apbs.pml"],
    capture_output=True, text=True
)
if result.returncode != 0:
    print(f"PyMOL visualization failed: {result.stderr}")
else:
    print("  PyMOL visualization completed.")

print(f"\nResults saved to:")
print(f"  {RESULTS}/step0_electrostatics.pse")
print(f"  {RESULTS}/step0_electrostatics.png")
