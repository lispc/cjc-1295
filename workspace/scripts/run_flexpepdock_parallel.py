#!/usr/bin/env python3
"""
Parallel Rosetta FlexPepDock runner.
Launches N independent processes, each with its own seed and thread limit.
"""

import os
import sys
import subprocess
import multiprocessing

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP2 = f"{BASE}/step2"
ROSETTA_BIN = "/home/scroll/miniforge3/envs/rosetta/bin"
DB = "/home/scroll/miniforge3/envs/rosetta/database"
INPUT_PDB = f"{STEP2}/prepacked_DPP4_GHRH_start_0001.pdb"

N_WORKERS = 16
MODELS_PER_WORKER = 60
THREADS_PER_WORKER = 8
BASE_SEED = 1542402926  # same base as original run

def run_worker(args):
    worker_id, nstruct, seed = args
    silent_out = f"{STEP2}/GHRH_DPP4_dock_worker_{worker_id:02d}.silent"
    log_out = f"{STEP2}/worker_{worker_id:02d}.log"

    cmd = [
        f"{ROSETTA_BIN}/FlexPepDocking",
        f"-database", DB,
        f"-s", INPUT_PDB,
        f"-lowres_preoptimize",
        f"-pep_refine",
        f"-ex1", f"-ex2aro",
        f"-use_input_sc",
        f"-nstruct", str(nstruct),
        f"-constant_seed",
        f"-jran", str(seed),
        f"-out:file:silent", silent_out,
    ]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(THREADS_PER_WORKER)

    print(f"[Worker {worker_id}] Starting {nstruct} models with seed {seed}...")
    with open(log_out, "w") as logfile:
        result = subprocess.run(cmd, env=env, stdout=logfile, stderr=subprocess.STDOUT)

    print(f"[Worker {worker_id}] Finished with exit code {result.returncode}")
    return worker_id, result.returncode

def combine_results():
    print("\n=== Combining silent files ===")
    silent_files = sorted([
        f for f in os.listdir(STEP2)
        if f.startswith("GHRH_DPP4_dock_worker_") and f.endswith(".silent")
    ])

    if not silent_files:
        print("No silent files found!")
        return

    cmd = [
        f"{ROSETTA_BIN}/combine_silent",
        "-out", f"{STEP2}/GHRH_DPP4_dock_combined.silent",
    ] + [f"{STEP2}/{f}" for f in silent_files]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Combined into {STEP2}/GHRH_DPP4_dock_combined.silent")
        # Count structures
        # Simple way: grep for a pattern
        total = sum(
            subprocess.run(
                ["grep", "-c", "reported success", f"{STEP2}/{f.replace('.silent', '.log')}"],
                capture_output=True, text=True
            ).stdout.strip() or "0"
            for f in silent_files
            if os.path.exists(f"{STEP2}/{f.replace('.silent', '.log')}")
        )
        print(f"Estimated models completed: {total}")
    else:
        print(f"combine_silent failed: {result.stderr}")

if __name__ == "__main__":
    tasks = [(i, MODELS_PER_WORKER, BASE_SEED + i * 1000) for i in range(N_WORKERS)]

    print(f"Launching {N_WORKERS} workers, {MODELS_PER_WORKER} models each, {THREADS_PER_WORKER} threads each")
    print(f"Total models: {N_WORKERS * MODELS_PER_WORKER}")
    print(f"Expected runtime: ~{(MODELS_PER_WORKER * 240) / 60:.0f} min per worker ≈ 4 hours total\n")

    # Use multiprocessing Pool
    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(run_worker, tasks)

    print("\n=== All workers finished ===")
    for wid, code in results:
        status = "OK" if code == 0 else f"FAILED ({code})"
        print(f"  Worker {wid}: {status}")

    combine_results()
