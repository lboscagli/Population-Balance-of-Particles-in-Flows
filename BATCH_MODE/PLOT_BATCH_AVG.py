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


# ---------------- User settings ----------------
parent_dir = "./results_LES_HS_streamlines_G04p75_new"

# analysis_name = 'Ta_225'
# legend_title ="$T=225\\,K$"
# groups = {
#     r"$\textnormal{H-S-1}$": [
#         "H-S-1_r0_Ta225",
#         "H-S-1_r4_Ta225",
#         "H-S-1_r5_Ta225"
#     ],
#     r"$\textnormal{H-S-2}$": [
#         "H-S-2_r0_Ta225",
#         "H-S-2_r4_Ta225",
#         "H-S-2_r5_Ta225"
#     ],
#     r"$\textnormal{H-S-3}$": [
#         "H-S-3_r0_Ta225",
#         "H-S-3_r4_Ta225",
#         "H-S-3_r5_Ta225"
#     ] 
# }

# analysis_name = 'Ta_225_CENTERLINE'
# legend_title ="$T=225\\,K$"
# groups = {
#     r"$\textnormal{H-S-1}$": [
#         "H-S-1_r0_Ta225"
#     ],
#     r"$\textnormal{H-S-2}$": [
#         "H-S-2_r0_Ta225"
#     ],
#     r"$\textnormal{H-S-3}$": [
#         "H-S-3_r0_Ta225"
#     ] 
# }

analysis_name = 'Ta_220'
legend_title ="$T=220\\,K$"
groups = {
    r"$\textnormal{H-S-1}$": [
        "H-S-1_r0_Ta220",
        "H-S-1_r4_Ta220",
        "H-S-1_r5_Ta220"
    ],
    r"$\textnormal{H-S-2}$": [
        "H-S-2_r0_Ta220",
        "H-S-2_r4_Ta220",
        "H-S-2_r5_Ta220"
    ],
    r"$\textnormal{H-S-3}$": [
        "H-S-3_r0_Ta220",
        "H-S-3_r4_Ta220",
        "H-S-3_r5_Ta220"
    ] 
}

# analysis_name = 'Ta_220_CENTERLINE'
# legend_title ="$T=220\\,K$"
# groups = {
#     r"$\textnormal{H-S-1}$": [
#         "H-S-1_r0_Ta220"
#     ],
#     r"$\textnormal{H-S-2}$": [
#         "H-S-2_r0_Ta220"
#     ],
#     r"$\textnormal{H-S-3}$": [
#         "H-S-3_r0_Ta220"
#     ] 
# }


# ---------------- Plot setup ----------------
fig, ax = plt.subplots()

# store all data for identity line range
N_nvpm_all = []
N_ice_all = []

# New: store fractions
vpm_fractions_all = []
nvpm_fractions_all = []

# ---------------- Main loop ----------------
group_names = []
vpm_avg_fractions = []
nvpm_avg_fractions = []

for i_group, (group_name, folder_list) in enumerate(groups.items()):
    
    markers = ['s','^','v'] #["o", "s", "D", "^", "v", "P", "X"]
    
    N_nvpm_list = []
    N_ice_list = []
    vpm_fraction_list = []
    nvpm_fraction_list = []

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

        # New: load fractions
        fraction_vpm = data.get('fraction_activated_vpm', np.array([0]))[0][0]
        fraction_nvpm = data.get('fraction_activated_nvpm', np.array([0]))[0][0]

        inps = read_inps_file(inps_file)
        N_nvpm = get_number_concentration_by_kappa(inps, lambda k: k < 0.1)

        N_nvpm_list.append(N_nvpm)
        N_ice_list.append(N_ice_final)
        vpm_fraction_list.append(fraction_vpm)
        nvpm_fraction_list.append(fraction_nvpm)

        N_nvpm_all.append(N_nvpm)
        N_ice_all.append(N_ice_final)
        vpm_fractions_all.append(fraction_vpm)
        nvpm_fractions_all.append(fraction_nvpm)

        print(f"{folder}: N_nvpm = {N_nvpm:.2e}, N_ice(final) = {N_ice_final:.2e}, frac_vpm = {fraction_vpm:.3f}, frac_nvpm = {fraction_nvpm:.3f}")

    if len(N_nvpm_list) == 0:
        continue

    N_nvpm_arr = np.array(N_nvpm_list)
    N_ice_arr = np.array(N_ice_list)
    vpm_fraction_arr = np.array(vpm_fraction_list)
    nvpm_fraction_arr = np.array(nvpm_fraction_list)
    
    #Average across the folder_list (this is equivalent to average along multiple streamlines for how the case has been designed)
    N_nvpm_avg = np.mean(N_nvpm_arr)
    N_ice_avg = np.mean(N_ice_arr)
    vpm_fraction_avg = np.mean(vpm_fraction_arr)
    nvpm_fraction_avg = np.mean(nvpm_fraction_arr)
    
    # Store for bar plots
    group_names.append(group_name)
    vpm_avg_fractions.append(vpm_fraction_avg * 100)  # as %
    nvpm_avg_fractions.append(nvpm_fraction_avg * 100)
    
    #sort_idx = np.argsort(N_nvpm_arr)

    ax.plot(
        N_nvpm_avg,
        N_ice_avg,
        marker=markers[i_group % len(markers)],
        linestyle="-",
        fillstyle="none",
        markersize=10,
        color='k',
        label=group_name,
    )


# ---------------- Identity line ----------------
N_nvpm_all = np.array(N_nvpm_all)
N_ice_all = np.array(N_ice_all)

mask = (N_nvpm_all > 0) & (N_ice_all > 0)
xmin = min(N_nvpm_all[mask].min(), N_ice_all[mask].min())
xmax = max(N_nvpm_all[mask].max(), N_ice_all[mask].max())

x_ref = np.logspace(np.log10(1e10), np.log10(1e13), 300)

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
ax.set_xlim(1e10,1e13)
ax.set_ylim(1e10,1e13)

ax.set_xlabel(r"$N_{\mathrm{nvPM}}\;[\#/m^3]$", fontsize=18)
ax.set_ylabel(r"$N_{\mathrm{ice}}\;[\#/m^3]$", fontsize=18)

ax.tick_params(labelsize=16)
ax.legend(fontsize=13,title=legend_title,title_fontsize=13)
ax.grid(True, which="both", ls="--", alpha=0.5)

plt.tight_layout()
plt.savefig(parent_dir+"/AVG_NnvPM_vs_Nice_summary_"+analysis_name+".png", dpi=600)
plt.show()

# ---------------- New: Bar plot for activated fractions ----------------
fig2, ax2 = plt.subplots()
x_pos = np.arange(len(group_names))
bar_width = 0.35

bars1 = ax2.bar(x_pos - bar_width/2, nvpm_avg_fractions, bar_width, color='red', alpha=0.7, hatch='//', label=r'$\textnormal{nvPM}$')
bars2 = ax2.bar(x_pos + bar_width/2, vpm_avg_fractions, bar_width, color='black', alpha=0.7, hatch='\\\\', label=r'$\textnormal{vPM}$')

ax2.set_xlabel(r'$\textnormal{Case}$', fontsize=18)
ax2.set_ylabel(r'$\textnormal{Activated Fraction [\%]}$', fontsize=18)
ax2.set_xticks(x_pos)
ax2.set_xticklabels(group_names, rotation=0,fontsize=16)
ax2.tick_params(labelsize=16)
ax2.legend(loc='best',fontsize=14)
plt.tight_layout()
plt.savefig(parent_dir+"/Activated_Fraction_"+analysis_name+".png", dpi=600)
plt.show()

print("\nSummary post-processing completed.")
