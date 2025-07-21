# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:37:28 2025

@author: lboscagl
"""

import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt

# === USER INPUT ===
Nbins = 55
bins_between_v01_v02 = 5  # exactly 5 bins between v01 and v02
r01 = 2e-9
r02 = 20e-9
n01 = 0.0
n02 = 4.5e11
rho1 = 1500
rho2 = 1550
kappa1 = 0.5
kappa2 = 0.005
output_filename = '../psr/ice_nucleating_particles.in'

# === Compute volumes from radii ===
v01 = (4.0 / 3.0) * np.pi * r01**3
v02 = (4.0 / 3.0) * np.pi * r02**3

# === Setup the right tail max volume ===
vmax = v02 * 2e3  # or some multiplier > v02 to cover the tail

# === Part 1: Create bins_between_v01_v02 bins between v01 and v02 (inclusive) ---
part1_bins = bins_between_v01_v02
part1_volumes = np.geomspace(v01, v02, part1_bins)

# === Part 2: Create remaining bins (Nbins - bins_between_v01_v02) from v02 to vmax using standard lognormal quantiles
part2_bins = Nbins - part1_bins
if part2_bins > 0:
    q = np.linspace(0, 1, part2_bins)
    # Avoid exactly 0 and 1 quantiles
    q = np.clip(q, 1e-6, 1 - 1e-6)

    log_mean = np.log(np.sqrt(v02 * vmax))
    log_sigma = np.log(vmax / v02) / 4
    dist_tail = lognorm(s=log_sigma, scale=np.exp(log_mean))
    part2_volumes = dist_tail.ppf(q)
else:
    part2_volumes = np.array([])

# === Combine two parts, remove duplicate v02 if present ===
if part2_volumes.size > 0 and np.isclose(part2_volumes[0], v02):
    part2_volumes = part2_volumes[1:]

grid = np.concatenate((part1_volumes, part2_volumes))

# Ensure grid length matches Nbins
if len(grid) > Nbins:
    grid = grid[:Nbins]
elif len(grid) < Nbins:
    last = grid[-1]
    extra_bins = Nbins - len(grid)
    extra_volumes = last * np.geomspace(1.1, 1.1**extra_bins, extra_bins)
    grid = np.concatenate((grid, extra_volumes))

# === Initialize arrays ===
rho = np.zeros(Nbins)
kappa = np.zeros(Nbins)
n = np.zeros(Nbins)

# === Assign values for v01 (first bin) ===
rho[0] = rho1
kappa[0] = kappa1
n[0] = n01

# === Assign values at bin closest to v02 ===
idx_v02 = (np.abs(grid - v02)).argmin()
rho[idx_v02] = rho2
kappa[idx_v02] = kappa2
n[idx_v02] = n02

# === Write output file ===
with open(output_filename, 'w') as f:
    f.write(f"{Nbins}\n")
    for i in range(Nbins):
        f.write(f"{grid[i]:.8e} {rho[i]:.8e} {kappa[i]:.8e} {n[i]:.8e}\n")

print(f"File '{output_filename}' written successfully.")


# === Plot for checks only ===
plt.semilogx(grid, np.arange(Nbins), marker='o')
plt.xlabel("Volume")
plt.ylabel("Bin index")
plt.title("Grid point distribution (lognormal, right-tail focused)")
plt.grid(True)
plt.show()

