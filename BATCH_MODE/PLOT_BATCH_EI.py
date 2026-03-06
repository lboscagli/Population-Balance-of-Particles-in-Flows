# -*- coding: utf-8 -*-
"""
Summary post-processing across multiple result groups

x-axis: N_nvPM
y-axis: N_ice at final timestep
"""

import os
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.rcParams['text.usetex'] = True

# ---------------- Helper functions ----------------
def read_inps_file(filename):
    """
    Reads ice_nucleating_particles.in and returns a dictionary
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    try:
        Nbins = int(lines[0])
        idx = 1
    except ValueError:
        Nbins = int(lines[1])
        idx = 2

    kappa = np.zeros(Nbins)
    n = np.zeros(Nbins)

    for i in range(Nbins):
        idx += 1  # skip bin edge
        parts = lines[idx].split()
        idx += 1
        kappa[i] = float(parts[3])
        n[i] = float(parts[4])

    return {"kappa": kappa, "n": n}


def get_number_concentration_by_kappa(inps_dict, kappa_condition):
    """
    Sum number concentration for bins satisfying kappa_condition(kappa)
    """
    mask = kappa_condition(inps_dict["kappa"])
    return np.sum(inps_dict["n"][mask])


def concentration_to_emission_index(
    C_N_jet,        # number concentration at jet conditions [#/m^3]
    Q_fuel_mL_min,  # fuel volumetric flow [mL/min]
    rho_fuel,       # fuel density [kg/m^3]
    Q_exh_slpm,     # exhaust flow [slpm at T_std, p_std]
    T_jet,          # jet temperature [K]
    p_jet,          # jet pressure [Pa]
    T_std=273.15,   # standard temperature for slpm [K]
    p_std=101325.0  # standard pressure for slpm [Pa]
):
    """
    C_N at jet conditions -> EI.

    Returns:
        EI         : emission index [#/kg]
        N_dot      : particle emission rate [#/s]
        m_fuel_dot : fuel mass flow rate [kg/s]
        V_exh_jet  : exhaust volumetric flow at jet conditions [m^3/s]
    """

    # 1) Exhaust volumetric flow at jet conditions [m^3/s]
    V_exh_std = Q_exh_slpm * 1e-3 / 60.0       # [m^3/s at std]
    V_exh_jet = V_exh_std * (T_jet / T_std) * (p_std / p_jet)

    # 2) Particle emission rate [#/s] from C_N and V_exh_jet
    N_dot = C_N_jet * V_exh_jet

    # 3) Fuel mass flow [kg/s]
    V_fuel_dot = Q_fuel_mL_min * 1e-6 / 60.0   # [m^3/s]
    m_fuel_dot = V_fuel_dot * rho_fuel         # [kg/s]

    # 4) Emission index [#/kg]
    EI = N_dot / m_fuel_dot

    return EI, N_dot, m_fuel_dot, V_exh_jet


# ---------------- User settings ----------------
parent_dir = "./results"

legend_title ="$EI_{\mathrm{vPM}} = 1\\times10^{17}$"
groups = {
    "$T = 220\\,K$": [
        "JP25_EIv1e17_EInv1e12_T220",
        "JP25_EIv1e17_EInv1e13_T220",
        "JP25_EIv1e17_EInv1e14_T220",
        "JP25_EIv1e17_EInv1e15_T220",
        "JP25_EIv1e17_EInv1e16_T220",
    ],
    "$T = 215\\,K$": [
        "JP25_EIv1e17_EInv1e12_T215",
        "JP25_EIv1e17_EInv1e13_T215",
        "JP25_EIv1e17_EInv1e14_T215",
        "JP25_EIv1e17_EInv1e15_T215",
        "JP25_EIv1e17_EInv1e16_T215",
    ],
    "$T = 213\\,K$": [
        "JP25_EIv1e17_EInv1e12_T213",
        "JP25_EIv1e17_EInv1e13_T213",
        "JP25_EIv1e17_EInv1e14_T213",
        "JP25_EIv1e17_EInv1e15_T213",
        "JP25_EIv1e17_EInv1e16_T213",
    ],    
    "$T = 211\\,K$": [
        "JP25_EIv1e17_EInv1e12_T211",
        "JP25_EIv1e17_EInv1e13_T211",
        "JP25_EIv1e17_EInv1e14_T211",
        "JP25_EIv1e17_EInv1e15_T211",
        "JP25_EIv1e17_EInv1e16_T211",
    ],
}

# ---------------- Plot setup ----------------
colors = cm.cool(np.linspace(0, 1, len(groups)))
markers = ["o", "s", "D", "^", "v", "P", "X"]

fig, ax = plt.subplots()

# store all data for identity line range
N_nvpm_all = []
N_ice_all = []

# ---------------- Main loop ----------------
for i_group, (group_name, folder_list) in enumerate(groups.items()):
    N_nvpm_list = []
    N_ice_list = []

    print(f"\n--- Processing {group_name} ---")

    for folder in folder_list:
        post_dir = os.path.join(parent_dir, folder, "post-processing")
        stats_file = os.path.join(post_dir, "statistics.mat")
        inps_file = os.path.join(parent_dir, folder, "ice_nucleating_particles.in")

        if not os.path.exists(stats_file) or not os.path.exists(inps_file):
            print(f"WARNING: missing files in {folder}")
            continue

        data = loadmat(stats_file)
        moment_0_ice = data["moment_0_ice"].flatten()
        N_ice_final = moment_0_ice[-1]

        inps = read_inps_file(inps_file)
        N_nvpm = get_number_concentration_by_kappa(inps, lambda k: k < 0.1)

        N_nvpm_list.append(N_nvpm)
        N_ice_list.append(N_ice_final)

        N_nvpm_all.append(N_nvpm)
        N_ice_all.append(N_ice_final)

        print(f"{folder}: N_nvpm = {N_nvpm:.2e}, N_ice(final) = {N_ice_final:.2e}")

    if len(N_nvpm_list) == 0:
        continue

    N_nvpm_arr = np.array(N_nvpm_list)
    N_ice_arr = np.array(N_ice_list)

    sort_idx = np.argsort(N_nvpm_arr)

    ax.plot(
        N_nvpm_arr[sort_idx],
        N_ice_arr[sort_idx],
        marker=markers[i_group % len(markers)],
        linestyle="-",
        fillstyle="none",
        markersize=10,
        color=colors[i_group],
        label=group_name,
    )


# ---------------- Identity line ----------------
N_nvpm_all = np.array(N_nvpm_all)
N_ice_all = np.array(N_ice_all)

mask = (N_nvpm_all > 0) & (N_ice_all > 0)
xmin = min(N_nvpm_all[mask].min(), N_ice_all[mask].min())
xmax = max(N_nvpm_all[mask].max(), N_ice_all[mask].max())

x_ref = np.logspace(np.log10(xmin), np.log10(xmax), 300)

ax.plot(
    x_ref,
    x_ref,
    linestyle=":",
    color="black",
    linewidth=2.0,
    alpha=0.75,
    label=r"$N_{\mathrm{ice}} = N_{\mathrm{nvPM}}$",
)

# ---------------- Final formatting ----------------
ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$N_{\mathrm{nvPM}}\;[\#/m^3]$", fontsize=18)
ax.set_ylabel(r"$N_{\mathrm{ice}}\;[\#/m^3]$", fontsize=18)

ax.tick_params(labelsize=16)
ax.legend(fontsize=13,title=legend_title,title_fontsize=13)
ax.grid(True, which="both", ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig(parent_dir+"/NnvPM_vs_Nice_summary.png", dpi=600)
plt.show()

print("\nSummary post-processing completed.")
