#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch driver for cpmod_mp with automated multiple runs.

Workflow per case:
1. Generate ice_nucleating_particles.in and ice.in
2. Pre-run cleanup: remove old *.out files from pbe/
3. Run cpmod_mp
4. Copy pbe outputs to WORK_DIR/results/<case_label>_<timestamp>
5. Copy input files for reference
6. Post-run cleanup: remove everything from pbe/ except pbe.in
"""

import os
import subprocess
import shutil
import glob
from datetime import datetime
import INPS_WRITER as INPS

# ==================================================
# Paths and directories
# ==================================================
WORK_DIR = os.path.dirname(os.path.abspath(__file__))
CPMOD_DIR = os.path.abspath(os.path.join(WORK_DIR, ".."))
INPUT_DIR = os.path.join(CPMOD_DIR, "psr")
PBE_OUTPUT_DIR = os.path.join(CPMOD_DIR, "pbe")
RESULTS_ROOT = os.path.join(WORK_DIR, "results")
os.makedirs(RESULTS_ROOT, exist_ok=True)
FORTRAN_EXE = os.path.join(CPMOD_DIR, "cpmod_mp")

ICE_INP_FILE = os.path.join(INPUT_DIR, "ice_nucleating_particles.in")
ICE_FILE = os.path.join(INPUT_DIR, "ice.in")

# Safety checks
assert os.path.isfile(FORTRAN_EXE), f"Executable not found: {FORTRAN_EXE}"
assert os.path.isdir(INPUT_DIR), f"Input directory not found: {INPUT_DIR}"
assert os.path.isdir(PBE_OUTPUT_DIR), f"PBE directory not found: {PBE_OUTPUT_DIR}"

# ==================================================
# Batch run function
# ==================================================
def run_case(mode1, mode2, Nbins, threshold_diameter, ice_params, case_label=None):
    """
    Run cpmod_mp for given parameters and save outputs.
    """
    print(f"\n=== Running case: {case_label if case_label else 'unnamed'} ===")

    # --- Write input files ---
    print("Writing ice_nucleating_particles.in")
    INPS.write_ice_inp_file(
        output_filename=ICE_INP_FILE,
        Nbins=Nbins,
        mode1=mode1,
        mode2=mode2,
        threshold_diameter=threshold_diameter
    )

    print("Writing ice.in")
    INPS.write_ice_file(
        output_filename=ICE_FILE,
        **ice_params
    )

    # --- Pre-run cleanup: remove *.out files from pbe/ ---
    print("Pre-run cleanup: removing old .out files from pbe/")
    for out_file in glob.glob(os.path.join(PBE_OUTPUT_DIR, "*.out")):
        try:
            os.remove(out_file)
        except Exception as e:
            print(f"Could not remove {out_file}: {e}")

    # --- Run Fortran code ---
    print("Running cpmod_mp")
    subprocess.run([FORTRAN_EXE], cwd=CPMOD_DIR, check=True)
    print("cpmod_mp completed successfully")

    # --- Prepare case directory ---
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    case_name = f"{case_label}" if case_label else f"case_{timestamp}"
    CASE_DIR = os.path.join(RESULTS_ROOT, case_name)
    os.makedirs(CASE_DIR, exist_ok=False)

    # --- Copy outputs ---
    print(f"Copying outputs to {CASE_DIR}")
    shutil.copytree(PBE_OUTPUT_DIR, CASE_DIR, dirs_exist_ok=True)

    # --- Copy input files for reference ---
    print("Copying input files to case directory for reference")
    for input_file in [ICE_FILE, ICE_INP_FILE]:
        try:
            shutil.copy(input_file, CASE_DIR)
        except Exception as e:
            print(f"Could not copy {input_file}: {e}")

    # --- Post-run cleanup: keep only pbe.in ---
    KEEP_FILE = "pbe.in"
    for name in os.listdir(PBE_OUTPUT_DIR):
        path = os.path.join(PBE_OUTPUT_DIR, name)
        if name != KEEP_FILE:
            try:
                if os.path.isfile(path) or os.path.islink(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
            except Exception as e:
                print(f"Could not remove {path}: {e}")

    print(f" Case '{case_name}' completed successfully")
    return CASE_DIR

# ==================================================
# USER INPUTs for BATCH RUN MODE
# ==================================================
batch_params = [
    {
        "mode1": dict(GMD=2e-9, GSD=1.3, N0=2.1e12, rho=1800, kappa=0.5),
        "mode2": dict(GMD=35e-9, GSD=2.0, N0=1.05e11, rho=1550, kappa=0.005),
        "Nbins": 35,
        "threshold_diameter": 80e-9,
        "ice_params": dict(
            p_amb=18753.93,
            T_amb=225.15,
            mixing_slope=2.75,
            rho_left=1550.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=0.018,
            jet_velocity=96.4,
            jet_temperature=363.15,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=0.0,
        ),
        "label": "H-S_Nvpm2.1e12_Nnvpm1.05e11"
    },
    {
        "mode1": dict(GMD=2e-9, GSD=1.3, N0=2.1e14, rho=1800, kappa=0.5),
        "mode2": dict(GMD=35e-9, GSD=2.0, N0=1.05e11, rho=1550, kappa=0.005),
        "Nbins": 35,
        "threshold_diameter": 80e-9,
        "ice_params": dict(
            p_amb=18753.93,
            T_amb=225.15,
            mixing_slope=2.75,
            rho_left=1550.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=0.018,
            jet_velocity=96.4,
            jet_temperature=363.15,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=0.0,
        ),
        "label": "H-S_Nvpm2.1e14_Nnvpm1.05e11"
    },
    {
        "mode1": dict(GMD=2e-9, GSD=1.3, N0=2.1e14, rho=1800, kappa=0.5),
        "mode2": dict(GMD=35e-9, GSD=2.0, N0=0e11, rho=1550, kappa=0.005),
        "Nbins": 35,
        "threshold_diameter": 80e-9,
        "ice_params": dict(
            p_amb=18753.93,
            T_amb=225.15,
            mixing_slope=2.75,
            rho_left=1550.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=0.018,
            jet_velocity=96.4,
            jet_temperature=363.15,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=0.0,
        ),
        "label": "H-S_Nvpm2.1e14_Nnvpm0"
    },
    {
        "mode1": dict(GMD=2e-9, GSD=1.3, N0=0e14, rho=1800, kappa=0.5),
        "mode2": dict(GMD=35e-9, GSD=2.0, N0=1.05e11, rho=1550, kappa=0.005),
        "Nbins": 35,
        "threshold_diameter": 80e-9,
        "ice_params": dict(
            p_amb=18753.93,
            T_amb=225.15,
            mixing_slope=2.75,
            rho_left=1550.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=0.018,
            jet_velocity=96.4,
            jet_temperature=363.15,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=0.0,
        ),
        "label": "H-S_Nvpm0_Nnvpm1.05e11"
    }          

]

# Loop over all cases
for params in batch_params:
    run_case(
        mode1=params["mode1"],
        mode2=params["mode2"],
        Nbins=params["Nbins"],
        threshold_diameter=params["threshold_diameter"],
        ice_params=params["ice_params"],
        case_label=params.get("label")
    )

print("\n All batch runs completed successfully")
