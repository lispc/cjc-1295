#!/usr/bin/env python3
"""
Run D-Ala2 MD using OpenMM, reading GROMACS topology and coordinates directly.
This bypasses the GROMACS 2026.0 segfault issue.
"""

import sys
sys.path.insert(0, '/home/scroll/miniforge3/envs/base/lib/python3.13/site-packages')

from openmm import app, Platform, LangevinMiddleIntegrator, MonteCarloBarostat
import openmm.unit as unit

# Load GROMACS topology and coordinates
gmx_top_path = '/home/scroll/miniforge3/envs/gmx/share/gromacs/top'
top = app.GromacsTopFile(
    '/home/scroll/personal/cjc-1295/workspace/step3/DAla2_topol.top',
    periodicBoxVectors=app.GromacsGroFile(
        '/home/scroll/personal/cjc-1295/workspace/step3/DAla2_em.gro'
    ).getPeriodicBoxVectors(),
    includeDir=gmx_top_path,
)
gro = app.GromacsGroFile('/home/scroll/personal/cjc-1295/workspace/step3/DAla2_em.gro')

# Create system with PME, constraints, and dispersion correction
system = top.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)

# Add pressure coupling
system.addForce(MonteCarloBarostat(1.0*unit.bar, 310*unit.kelvin))

# Create integrator (Langevin for temperature control)
integrator = LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)

# Try CUDA first, fallback to CPU
for platform_name in ['CUDA', 'CPU']:
    try:
        platform = Platform.getPlatformByName(platform_name)
        if platform_name == 'CUDA':
            props = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': '1'}
        else:
            props = {}
        simulation = app.Simulation(top.topology, system, integrator, platform, props)
        print(f"Using platform: {platform_name}")
        break
    except Exception as e:
        print(f"{platform_name} failed: {e}")
        continue
else:
    raise RuntimeError("No platform available")

# Set positions and box vectors
simulation.context.setPositions(gro.positions)
simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())

# Energy minimization
print("Energy minimization...")
simulation.minimizeEnergy(maxIterations=1000)
state = simulation.context.getState(getEnergy=True)
print(f"Minimized potential energy: {state.getPotentialEnergy()}")

# Short equilibration (100 ps)
print("Equilibration (100 ps)...")
simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True, 
    density=True, speed=True, elapsedTime=True
))
simulation.reporters.append(app.DCDReporter(
    '/home/scroll/personal/cjc-1295/workspace/step3/DAla2_openmm_eq.dcd', 5000
))
simulation.step(50000)  # 100 ps

# Production (10 ns as a test)
print("\nProduction (10 ns)...")
simulation.reporters.clear()
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, 50000, step=True, potentialEnergy=True, temperature=True,
    density=True, speed=True, elapsedTime=True
))
simulation.reporters.append(app.DCDReporter(
    '/home/scroll/personal/cjc-1295/workspace/step3/DAla2_openmm_md.dcd', 50000
))
simulation.step(5000000)  # 10 ns

print("Done!")
