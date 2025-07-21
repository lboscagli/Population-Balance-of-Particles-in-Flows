# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:37:28 2025

@author: lboscagl
"""

import numpy as np
from scipy.stats import lognorm

# === USER INPUT ===
Nbins = 50
r01 = 2E-9
r02 = 20E-9
rho1 = 1000
rho2 = 1550
kappa1 = 0.5
kappa2 = 0.005
output_filename = 'ice_nucleating_particles.in'

# === Convert radius to volume ===
v01 = 4/3*np.pi*r01**3
v02 = 4/3*np.pi*r02**3

# === Compute Lognormal Distribution Parameters ===
vmin = min(v01, v02) / 2.0
vmax = max(v01, v02) * 2.0

log_mean = np.log(np.sqrt(v01 * v02))
log_sigma = np.log(vmax / vmin) / 4

dist = lognorm(s=log_sigma, scale=np.exp(log_mean))

# === Generate lognormal grid excluding v01 and v02 ===
eps = 1e-3
quantiles = np.linspace(eps, 1 - eps, Nbins - 2)
grid_core = dist.ppf(quantiles)

# Add v01 and v02 explicitly to the grid
grid_full = np.concatenate((grid_core, [v01, v02]))
grid_sorted = np.sort(grid_full)

# === Initialize rho and kappa arrays ===
rho = np.zeros(Nbins)
kappa = np.zeros(Nbins)

# Set values exactly at v01 and v02
for i, v in enumerate(grid_sorted):
    if np.isclose(v, v01, rtol=1e-10):
        rho[i] = rho1
        kappa[i] = kappa1
    elif np.isclose(v, v02, rtol=1e-10):
        rho[i] = rho2
        kappa[i] = kappa2

# === Write to File ===
with open(output_filename, 'w') as f:
    f.write(f"{Nbins}\n")
    for i in range(Nbins):
        f.write(f"{grid_sorted[i]:.8e} {rho[i]:.8e} {kappa[i]:.8e}\n")

print(f"File '{output_filename}' with v01 and v02 included written successfully.")
