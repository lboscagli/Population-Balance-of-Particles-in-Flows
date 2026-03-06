#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch driver for cpmod_mp with automated multiple runs.

Workflow per case:
1. Generate ice_nucleating_particles.in and ice.in
2. Pre-run cleanup: remove old *.out files from pbe/
3. Run cpmod_mp
4. Copy pbe outputs to WORK_DIR/results/<case_label>
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

assert os.path.isfile(FORTRAN_EXE)
assert os.path.isdir(INPUT_DIR)
assert os.path.isdir(PBE_OUTPUT_DIR)

# ==================================================
# Core runner (shared logic)
# ==================================================
def _run_cpmod(case_label):
    """Execute cpmod_mp and archive results."""
    # --- Pre-run cleanup ---
    for f in glob.glob(os.path.join(PBE_OUTPUT_DIR, "*.out")):
        os.remove(f)

    subprocess.run([FORTRAN_EXE], cwd=CPMOD_DIR, check=True)

    CASE_DIR = os.path.join(RESULTS_ROOT, case_label)
    os.makedirs(CASE_DIR, exist_ok=False)

    shutil.copytree(PBE_OUTPUT_DIR, CASE_DIR, dirs_exist_ok=True)

    for f in [ICE_FILE, ICE_INP_FILE]:
        shutil.copy(f, CASE_DIR)

    # --- Post-run cleanup ---
    for name in os.listdir(PBE_OUTPUT_DIR):
        if name != "pbe.in":
            path = os.path.join(PBE_OUTPUT_DIR, name)
            if os.path.isfile(path):
                os.remove(path)
            else:
                shutil.rmtree(path)

    return CASE_DIR

# ==================================================
# LEGACY: run case from number concentration
# ==================================================
def run_case(mode1, mode2, Nbins, threshold_diameter, ice_params, case_label):
    print(f"\n=== Running case (N-based): {case_label} ===")

    INPS.write_ice_inp_file(
        output_filename=ICE_INP_FILE,
        Nbins=Nbins,
        mode1=mode1,
        mode2=mode2,
        threshold_diameter=threshold_diameter,
    )

    INPS.write_ice_file(ICE_FILE, **ice_params)

    return _run_cpmod(case_label)

# ==================================================
# NEW: run case from emission index (EI)
# ==================================================
def run_case_N_from_EI(
    mode1,
    mode2,
    jet_params,
    Nbins,
    threshold_diameter,
    ice_params,
    case_label,
):
    print(f"\n=== Running case (EI-based): {case_label} ===")

    INPS.write_ice_inp_file_from_EI(
        output_filename=ICE_INP_FILE,
        Nbins=Nbins,
        mode1=mode1,
        mode2=mode2,
        jet_params=jet_params,
        ice_params=ice_params,
        threshold_diameter=threshold_diameter,
    )

    INPS.write_ice_file(ICE_FILE, **ice_params)

    return _run_cpmod(case_label)

# # ==================================================
# # USER INPUTS — N-based cases
# # ==================================================
# batch_params_N = [
#     dict(
#         mode1=dict(GMD=2e-9, GSD=1.3, N0=2.1e14, rho=1800, kappa=0.5),
#         mode2=dict(GMD=35e-9, GSD=2.0, N0=1.05e11, rho=1550, kappa=0.005),
#         Nbins=35,
#         threshold_diameter=80e-9,
#         ice_params=dict(
#             p_amb=18753.93,
#             T_amb=225.15,
#             mixing_slope=2.2,
#             rho_left=1550.0,
#             alpha=0.1,
#             jet_temp_model_flag=1,
#             jet_diameter=0.018,
#             jet_velocity=96.4,
#             jet_temperature=363.15,
#             kappa=0.005,
#             ss_consumption=True,
#             read_inp_distribution=True,
#             RHi_amb=0.0,
#         ),
#         label="JP25_Nbased",
#     )
# ]

# # ==================================================
# # Run all cases
# # ==================================================
# for p in batch_params_N:
#     run_case(**p)

# ==================================================
# USER INPUTS — EI-based cases
# ==================================================
batch_params_EI = [
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e16, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=23000,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=211,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=0.0,
        ),
        case_label="JP25_EIv1e17_EInv1e16_T211",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e15, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=23000,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=211,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e15_T211",
    ),    
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e14, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=211,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e14_T211",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e13, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=211,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e13_T211",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e12, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=211,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e12_T211",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e16, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=215,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e16_T215",
    ),    
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e15, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=215,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e15_T215",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e14, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=215,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e14_T215",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e13, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=215,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e13_T215",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e12, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=215,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e12_T215",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e16, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=220,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e16_T220",
    ),    
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e15, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=220,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e15_T220",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e14, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=220,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e14_T220",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e13, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=220,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e13_T220",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e12, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=220,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e12_T220",
    ), 
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e16, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=213,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e16_T213",
    ),     
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e15, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=213,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e15_T213",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e14, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=213,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e14_T213",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e13, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=213,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e13_T213",
    ),
    dict(
        mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
        mode2=dict(GMD=35e-9, GSD=2.0, EI=1e12, rho=1500, kappa=0.005),
        jet_params=dict(
            T_jet=600.0,
            pressure=22919.57,
            xh2o=0.0276,
            Phi=0.1613,
            f_st=0.08957,
        ),
        Nbins=35,
        threshold_diameter=50e-9,
        ice_params=dict(
            p_amb=22919.57,
            T_amb=213,
            mixing_slope=1.64,
            rho_left=1500.0,
            alpha=0.1,
            jet_temp_model_flag=1,
            jet_diameter=1.0,
            jet_velocity=418.85,
            jet_temperature=600.0,
            kappa=0.005,
            ss_consumption=True,
            read_inp_distribution=True,
            RHi_amb=1.0,
        ),
        case_label="JP25_EIv1e17_EInv1e12_T213",
    )
]

# ==================================================
# Run all cases
# ==================================================
for p in batch_params_EI:
    run_case_N_from_EI(**p)
    
    
    

print("\n All batch runs completed successfully")
