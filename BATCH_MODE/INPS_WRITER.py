# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 12:45:25 2026

@author: lboscagl
"""

import numpy as np
import os
from scipy.optimize import root_scalar

def vol_to_d(vol):
    return ((6.0 * vol / np.pi) ** (1/3))

def d_to_vol(d):
    return (np.pi/6.0 * d**3)

def lognormal_n_v(v, N0, GMD, GSD):
    d = (6 * v / np.pi) ** (1/3)
    n_d = (N0 / (np.sqrt(2 * np.pi) * np.log(GSD))) * (1/d) * \
          np.exp(- (np.log(d) - np.log(GMD))**2 / (2 * np.log(GSD)**2))
    jac = (1/3) * (6/np.pi)**(1/3) * v**(-2/3)
    return n_d * jac


def write_ice_inp_file(
    output_filename,
    Nbins,
    mode1,
    mode2,
    threshold_diameter,
):
    """
    Write ice_nucleating_particles.in for the Fortran code.
    """

    # --- Unpack modes ---
    GMD1, GSD1, N01, rho1, kappa1 = (
        mode1["GMD"], mode1["GSD"], mode1["N0"], mode1["rho"], mode1["kappa"]
    )
    GMD2, GSD2, N02, rho2, kappa2 = (
        mode2["GMD"], mode2["GSD"], mode2["N0"], mode2["rho"], mode2["kappa"]
    )

    # --- Bin grid ---
    vmin = (np.pi/6) * (GMD1 / GSD1)**3
    vmax = (np.pi/6) * (GMD2 * GSD2 * 2e3)**3

    def objective(r):
        edges0 = vmin / r
        edges = edges0 * r**np.arange(Nbins + 1)
        return edges[-1] - vmax

    sol = root_scalar(objective, bracket=[1.01, 100.0], method='bisect')
    if not sol.converged:
        raise RuntimeError("Could not solve for geometric ratio r.")

    r = sol.root
    edges0 = vmin / r
    edges = edges0 * r**np.arange(Nbins + 1)
    grid = 0.5 * (edges[:-1] + edges[1:])
    dv = edges[1:] - edges[:-1]

    # --- Number distributions ---
    n1 = np.zeros(Nbins)
    n2 = np.zeros(Nbins)

    for i in range(Nbins):
        n1[i] = lognormal_n_v(grid[i], N01, GMD1, GSD1) * dv[i]
        n2[i] = lognormal_n_v(grid[i], N02, GMD2, GSD2) * dv[i]

    # --- Thresholding ---
    for i in range(Nbins):
        if vol_to_d(grid[i]) > threshold_diameter:
            n1[i] = 0.0
            n2[i] = 0.0

    # --- Renormalize ---
    if n1.sum() > 0:
        n1 *= N01 / n1.sum()
    if n2.sum() > 0:
        n2 *= N02 / n2.sum()

    n = n1 + n2

    # --- Assign properties ---
    rho = np.zeros(Nbins)
    kappa = np.zeros(Nbins)
    flag_nuclei = ['.false.'] * Nbins
    threshold_v = d_to_vol(threshold_diameter)

    for i in range(Nbins):
        if n1[i] > n2[i]:
            rho[i] = rho1
            kappa[i] = kappa1
        else:
            rho[i] = rho2
            kappa[i] = kappa2

        if n[i] > 0 and grid[i] <= threshold_v:
            flag_nuclei[i] = '.true.'

    # --- Write file ---
    os.makedirs(os.path.dirname(output_filename) or ".", exist_ok=True)

    with open(output_filename, "w") as f:
        f.write(".true.\n")
        f.write(f"{Nbins}\n")
        for i in range(Nbins):
            f.write(f"{edges[i]:.16e}\n")
            f.write(
                f"{grid[i]:.16e} {dv[i]:.16e} "
                f"{rho[i]:.16e} {kappa[i]:.16e} "
                f"{n[i]:.16e} {flag_nuclei[i]}\n"
            )
        f.write(f"{edges[-1]:.16e}\n")

    return {
        "edges": edges,
        "grid": grid,
        "n": n,
        "n1": n1,
        "n2": n2,
    }


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


# mode1 = dict(GMD=2e-9, GSD=1.3, N0=2.1e14, rho=1800, kappa=0.5)
# mode2 = dict(GMD=35e-9, GSD=2.0, N0=1.05e11, rho=1550, kappa=0.005)

# write_ice_inp_file(
#     output_filename="./ice_nucleating_particles.in",
#     Nbins=35,
#     mode1=mode1,
#     mode2=mode2,
#     threshold_diameter=80e-9
# )
