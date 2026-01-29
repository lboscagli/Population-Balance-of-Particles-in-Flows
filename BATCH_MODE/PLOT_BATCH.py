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


# ---------------- User settings ----------------
parent_dir = "./results"

legend_title ="$N_{\mathrm{vPM}} = 2.1\\times10^{14}\\,\#/m^3$"
groups = {
    "$T = 225\\,K$": [
        "H-S_Nvpm2.1e14_Nnvpm1.05e8_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e9_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e10_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e11_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e12_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e13_T225",
        "H-S_Nvpm2.1e14_Nnvpm1.05e14_T225",
    ],
    "$T = 220\\,K$": [
        "H-S_Nvpm2.1e14_Nnvpm1.05e8_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e9_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e10_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e11_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e12_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e13_T220",
        "H-S_Nvpm2.1e14_Nnvpm1.05e14_T220",
    ],
    "$T = 215\\,K$": [
        "H-S_Nvpm2.1e14_Nnvpm1.05e8_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e9_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e10_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e11_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e12_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e13_T215",
        "H-S_Nvpm2.1e14_Nnvpm1.05e14_T215",
    ],
    "$T = 210\\,K$": [
        "H-S_Nvpm2.1e14_Nnvpm1.05e8_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e9_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e10_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e11_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e12_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e13_T210",
        "H-S_Nvpm2.1e14_Nnvpm1.05e14_T210",
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
plt.savefig("NnvPM_vs_Nice_summary.png", dpi=600)
plt.show()

print("\nSummary post-processing completed.")
