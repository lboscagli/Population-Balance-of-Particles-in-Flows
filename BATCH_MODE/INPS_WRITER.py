# -*- coding: utf-8 -*-
"""
INPS_WRITER.py

Utilities to generate ice_nucleating_particles.in for CPMOD

Supports:
- Number concentration input (N0)
- Emission index input (EI -> N0 conversion)

Author: lboscagl
"""

import numpy as np
import os
from scipy.optimize import root_scalar

# ============================================================
# -------------------- PHYSICS HELPERS -----------------------
# ============================================================

def es(T_c):
    """Saturation vapor pressure over liquid water [Pa]"""
    return 611.2 * np.exp(17.67 * T_c / (T_c + 243.5))


def rho_air(T, P, RH):
    """Moist air density [kg/m3]"""
    Rd = 8.314 / (28.96 / 1.0e3)
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / (Rd * Tv)


def d_to_vol(d):
    """Diameter [m] -> volume [m3]"""
    return (np.pi / 6.0) * d**3


def vol_to_d(v):
    """Volume [m3] -> diameter [m]"""
    return (6.0 * v / np.pi) ** (1 / 3)


def lognormal_n_v(v, N0, GMD, GSD):
    """Lognormal number distribution in volume space"""
    d = vol_to_d(v)
    n_d = (
        N0
        / (np.sqrt(2 * np.pi) * np.log(GSD))
        * (1 / d)
        * np.exp(-(np.log(d) - np.log(GMD)) ** 2 / (2 * np.log(GSD) ** 2))
    )
    jac = (1 / 3) * (6 / np.pi) ** (1 / 3) * v ** (-2 / 3)
    return n_d * jac


# ============================================================
# -------------------- CORE WRITER ----------------------------
# ============================================================

def _write_ice_file_from_modes(
    output_filename,
    Nbins,
    modes,
    threshold_diameter,
):
    """
    INTERNAL routine.

    modes: list of dicts, each containing:
        N0, GMD, GSD, rho, kappa
    """

    # ---------------- Bin grid ----------------
    GMDs = [m["GMD"] for m in modes]
    GSDs = [m["GSD"] for m in modes]

    vmin = d_to_vol(1.0e-9)#d_to_vol(min(GMDs) / max(GSDs))
    vmax = d_to_vol(max(GMDs) * max(GSDs) * 2e3)

    def objective(r):
        return (vmin / r) * r**Nbins - vmax

    sol = root_scalar(objective, bracket=[1.01, 100.0], method="bisect")
    if not sol.converged:
        raise RuntimeError("Could not determine bin grid.")

    r = sol.root
    edges = (vmin / r) * r ** np.arange(Nbins + 1)
    grid = 0.5 * (edges[:-1] + edges[1:])
    dv = edges[1:] - edges[:-1]

    # ---------------- Distributions ----------------
    n_modes = []
    for m in modes:
        n_arr = np.array([
            lognormal_n_v(grid[i], m["N0"], m["GMD"], m["GSD"]) * dv[i]
            for i in range(Nbins)
        ])
        n_modes.append(n_arr)

    n_total = np.sum(n_modes, axis=0)

    # ---------------- Threshold ----------------
    threshold_v = d_to_vol(threshold_diameter)
    threshold_v_min = d_to_vol(1.0e-9)
    for i in range(Nbins):
        if grid[i] > threshold_v:
            for nm in n_modes:
                nm[i] = 0.0
            n_total[i] = 0.0
        if grid[i] < threshold_v_min:
            for nm in n_modes:
                nm[i] = 0.0
            n_total[i] = 0.0            

    # ---------------- Renormalization ----------------
    for i, m in enumerate(modes):
        s = np.sum(n_modes[i])
        if s > 0:
            n_modes[i] *= m["N0"] / s

    n_total = np.sum(n_modes, axis=0)

    # ---------------- Bin properties ----------------
    rho = np.zeros(Nbins)
    kappa = np.zeros(Nbins)
    flag = ['.false.'] * Nbins

    for i in range(Nbins):
        dominant = np.argmax([nm[i] for nm in n_modes])
        rho[i] = modes[dominant]["rho"]
        kappa[i] = modes[dominant]["kappa"]
        if n_total[i] > 0 and grid[i] <= threshold_v:
            flag[i] = '.true.'

    # ---------------- Write file ----------------
    os.makedirs(os.path.dirname(output_filename) or ".", exist_ok=True)

    with open(output_filename, "w") as f:
        f.write(".true.\n")
        f.write(f"{Nbins}\n")
        for i in range(Nbins):
            f.write(f"{edges[i]:.16e}\n")
            f.write(
                f"{grid[i]:.16e} {dv[i]:.16e} "
                f"{rho[i]:.16e} {kappa[i]:.16e} "
                f"{n_total[i]:.16e} {flag[i]}\n"
            )
        f.write(f"{edges[-1]:.16e}\n")

    print(f"📄 ice_nucleating_particles.in written to {output_filename}")


# ============================================================
# -------------------- PUBLIC APIS ----------------------------
# ============================================================

def write_ice_file(
    output_filename,
    p_amb,
    T_amb,
    mixing_slope,
    rho_left,
    alpha,
    jet_temp_model_flag,
    jet_diameter,
    jet_velocity,
    jet_temperature,
    kappa,
    ss_consumption,
    read_inp_distribution,
    RHi_amb,
):
    """
    Writes ice.in input file for cpmod_mp
    """

    def fbool(val):
        return ".true." if val else ".false."

    with open(output_filename, "w") as f:
        f.write("ICE INPUT DATA\n\n")

        f.write(f"{p_amb:.6f}                             ! Ambient pressure [Pa]\n")
        f.write(f"{T_amb:.6f}                               ! Ambient temperature [K]\n")
        f.write(f"{mixing_slope:.6f}                       ! Mixing line slope [Pa/K]\n")
        f.write(f"{rho_left:.6f}                           ! Density of particles on left PSD [kg/m^3]\n")
        f.write(f"{alpha:.6f}                              ! alpha [-]\n")
        f.write(f"{jet_temp_model_flag:d}                                   ! Jet temperature model flag\n")
        f.write(f"{jet_diameter:.6f}                         ! Jet diameter [m]\n")
        f.write(f"{jet_velocity:.6f}                         ! Jet initial velocity [m/s]\n")
        f.write(f"{jet_temperature:.6f}                       ! Jet initial temperature [K]\n")
        f.write(f"{kappa:.6f}                               ! Hygroscopicity [-]\n")
        f.write(f"{fbool(ss_consumption):<36} ! Supersaturation consumption\n")
        f.write(f"{fbool(read_inp_distribution):<36} ! Read INPs from input file\n")
        f.write(f"{RHi_amb:.6f}                               ! Ambient ice relative humidity (RHi)\n")


    return None



def write_ice_inp_file(
    output_filename,
    Nbins,
    mode1,
    mode2,
    threshold_diameter,
):
    """
    Write ice input file using number concentrations (N0).

    modeX must contain:
        N0, GMD, GSD, rho, kappa
    """
    _write_ice_file_from_modes(
        output_filename,
        Nbins,
        [mode1, mode2],
        threshold_diameter,
    )


def write_ice_inp_file_from_EI(
    output_filename,
    Nbins,
    mode1,
    mode2,
    jet_params,
    ice_params,
    threshold_diameter,
):
    """
    Write ice input file using emission indices (EI).

    modeX must contain:
        EI, GMD, GSD, rho, kappa

    jet_params must contain:
        T_jet, pressure, xh2o, Phi, f_st
    """

    Tj = jet_params["T_jet"]
    P = jet_params["pressure"]
    xh2o = jet_params["xh2o"]
    Phi = jet_params["Phi"]
    f_st = jet_params["f_st"]

    pwater = xh2o * P
    psat = np.exp(
        54.842763 - 6763.22 / Tj - 4.21 * np.log(Tj)
        + 0.000367 * Tj
    )
    RH = np.clip(pwater / psat, 0.0, 1.0)
    rho_j = rho_air(Tj, P, RH)

    f_actual = Phi * f_st
    conversion = (rho_j * f_actual) / (1 + f_actual)
    
    #conversion =

    # Convert EI -> N0
    for m in (mode1, mode2):
        m["N0"] = m["EI"] * conversion

    _write_ice_file_from_modes(
        output_filename,
        Nbins,
        [mode1, mode2],
        threshold_diameter,
    )


# import INPS_WRITER as INPS

# INPS.write_ice_inp_file_from_EI(
#     output_filename="./psr/ice_nucleating_particles.in",
#     Nbins=55,
#     mode1=dict(GMD=2.5e-9, GSD=1.3, EI=1e17, rho=1800, kappa=0.5),
#     mode2=dict(GMD=35e-9, GSD=2.0, EI=1e15, rho=1500, kappa=0.005),
#     jet_params=dict(
#         T_jet=600.0,
#         pressure=22919.57,
#         xh2o=0.0276,
#         Phi=0.1613,
#         f_st=0.08957,
#     ),
#     threshold_diameter=50e-9,
# )
