#!/usr/bin/env python3
"""
Parallel Rosetta FlexPepDock refine with CoordinateConstraints.
Runs WT and D-Ala2, 200 models each, 16 workers total.
"""
import os
import sys
import subprocess
import multiprocessing
import glob

ROSETTA_BIN = "/home/scroll/miniforge3/envs/rosetta/bin"
DB = "/home/scroll/miniforge3/envs/rosetta/database"
STEP2 = "/home/scroll/personal/cjc-1295/workspace/step2"
N_WORKERS = 16
N_STRUCT_PER_WORKER = 25  # 16 * 25 = 400 total; 8 workers per system

def run_worker(args):
    """Run a single FlexPepDocking worker."""
    system, worker_id, seed, nstruct = args

    if system == "WT":
        input_pdb = f"{STEP2}/prepacked_DPP4_GHRH_start_0001.pdb"
        native_pdb = f"{STEP2}/prepacked_DPP4_GHRH_start_0001.pdb"
        out_prefix = f"{STEP2}/refine_constrained_WT"
    else:
        input_pdb = f"{STEP2}/prepacked_DPP4_GHRH_DAla2_start_0001.pdb"
        native_pdb = f"{STEP2}/prepacked_DPP4_GHRH_DAla2_start_0001.pdb"
        out_prefix = f"{STEP2}/refine_constrained_DAla2"

    out_silent = f"{out_prefix}_w{worker_id:02d}.silent"

    cmd = [
        f"{ROSETTA_BIN}/FlexPepDocking",
        f"-database", DB,
        f"-s", input_pdb,
        f"-native", native_pdb,
        f"-pep_refine",
        f"-constraints:cst_fa_file", f"{STEP2}/backbone_coord.cst",
        f"-constraints:cst_fa_weight", "10.0",
        f"-ex1", "-ex2aro", "-use_input_sc",
        f"-nstruct", str(nstruct),
        f"-constant_seed", "-jran", str(seed),
        f"-out:file:silent", out_silent,
    ]

    log = f"{out_prefix}_w{worker_id:02d}.log"
    with open(log, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

    return (system, worker_id, result.returncode, out_silent)


def combine_silents(system):
    out_prefix = f"{STEP2}/refine_constrained_{system}"
    silent_files = sorted(glob.glob(f"{out_prefix}_w*.silent"))

    if not silent_files:
        print(f"No silent files found for {system}")
        return

    combined = f"{out_prefix}_combined.silent"
    cmd = [
        f"{ROSETTA_BIN}/combine_silent",
        "-in:file:silent", *silent_files,
        "-out:file:silent", combined,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Combined {len(silent_files)} files -> {combined}")
    else:
        print(f"combine_silent failed for {system}: {result.stderr}")


if __name__ == "__main__":
    tasks = []
    # WT: workers 0-7
    for i in range(8):
        tasks.append(("WT", i, 100000 + i, N_STRUCT_PER_WORKER))
    # D-Ala2: workers 8-15
    for i in range(8):
        tasks.append(("DAla2", i + 8, 200000 + i, N_STRUCT_PER_WORKER))

    print(f"Starting {len(tasks)} workers ({N_STRUCT_PER_WORKER} models each)")
    print(f"Total: {len(tasks) * N_STRUCT_PER_WORKER} models")

    with multiprocessing.Pool(N_WORKERS) as pool:
        results = pool.map(run_worker, tasks)

    print("\n--- Worker Results ---")
    for system, wid, rc, sf in results:
        status = "OK" if rc == 0 else f"ERR({rc})"
        print(f"  {system} worker {wid}: {status} -> {sf}")

    print("\n--- Combining silent files ---")
    combine_silents("WT")
    combine_silents("DAla2")
