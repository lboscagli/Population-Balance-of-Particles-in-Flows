# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import root_scalar
# import os

# # === Utility Functions ===
# def es(T_c):
#     return 611.2 * np.exp(17.67 * T_c / (T_c + 243.5))

# def rho_air(T, P, RH):
#     Rd = 8.314 / (28.96 / 1.0e3)
#     qsat = RH * 0.622 * (es(T - 273.15) / P)
#     Tv = T * (1.0 + 0.61 * qsat)
#     return P / (Rd * Tv)

# # === USER INPUTS ===
# Nbins = 35
# r01 = 1e-9
# r02 = 20e-9
# EI01, EI02 = 1e17, 1e12
# rho1, rho2 = 1500, 1500
# kappa1, kappa2 = 0.5, 0.005
# output_filename = "../psr/ice_nucleating_particles.in"

# v01 = (4/3) * np.pi * r01**3
# v02 = (4/3) * np.pi * r02**3
# vmax = v02*1e10#(4/3) * np.pi * (1000e-9)**3  # User sets this!

# # === Jet parameters ===
# jet_V = 350.0
# diameter_jet = 1.0
# temp_jet = 500.0
# xh2o_jet = 0.059243
# f_st = 0.08957
# Phi = 0.38
# pressr = 20000.0
# temp_air = 220.0

# # === Environmental properties ===
# pwater = xh2o_jet * pressr
# p_water_sat_liq = np.exp(54.842763 - 6763.22 / temp_jet - 4.21 * np.log(temp_jet) +
#                           0.000367 * temp_jet + np.tanh(0.0415 * (temp_jet - 218.8)) *
#                           (53.878 - 1331.22 / temp_jet - 9.44523 * np.log(temp_jet) +
#                            0.014025 * temp_jet))
# RH = np.clip(pwater / p_water_sat_liq, 0.0, 1.0)
# rho_j = rho_air(temp_jet, pressr, RH)

# # === Find geometric ratio r so that grid[0]=v01 and edges[Nbins]=vmax ===
# def objective(r):
#     edges0 = 2*v01/(1 + r)
#     return edges0 * r**Nbins - vmax

# sol = root_scalar(objective, bracket=[1.01, 100.0], method='bisect')
# if not sol.converged:
#     raise RuntimeError("Could not solve for geometric ratio r.")

# r = sol.root
# edges0 = 2*v01/(1 + r)
# edges = edges0 * r**np.arange(Nbins+1)
# grid = 0.5 * (edges[:-1] + edges[1:])
# dv = edges[1:] - edges[:-1]

# # Check centering
# assert np.allclose(grid, 0.5 * (edges[:-1] + edges[1:]), rtol=1e-12), "Grid not centered between edges!"
# assert np.isclose(grid[0], v01, rtol=1e-12), "First grid point not equal to v01!"
# assert np.isclose(edges[-1], vmax, rtol=1e-12), "Last edge not equal to vmax!"

# # Find which grid point is closest to v02 in first 5 bins
# idx02 = np.argmin(np.abs(grid[:] - v02))
# idx01 = 0  # by construction

# # === Assign physical properties ===
# rho = np.zeros(Nbins)
# kappa = np.zeros(Nbins)
# n = np.zeros(Nbins)
# flag_nuclei = ['.false.'] * Nbins

# f_actual = Phi * f_st

# if EI01 > 0:
#     rho[idx01] = rho1
#     kappa[idx01] = kappa1
#     n[idx01] = EI01 * (rho_j * f_actual) / (1 + f_actual)
#     flag_nuclei[idx01] = '.true.'

# if EI02 > 0:
#     rho[idx02] = rho2
#     kappa[idx02] = kappa2
#     n[idx02] = EI02 * (rho_j * f_actual) / (1 + f_actual)
#     flag_nuclei[idx02] = '.true.'

# # === Write to output file ===
# os.makedirs(os.path.dirname(output_filename), exist_ok=True)
# with open(output_filename, 'w') as f:
#     f.write(f"{Nbins}\n")
#     for i in range(Nbins):
#         f.write(f"{edges[i]:.16e}\n")
#         f.write(f"{grid[i]:.16e} {dv[i]:.16e} {rho[i]:.16e} {kappa[i]:.16e} {n[i]:.16e} {flag_nuclei[i]}\n")
#     f.write(f"{edges[-1]:.16e}\n")

# print(f"📄 File '{output_filename}' written successfully.")
# print(f"v01 = {v01:.3e}, grid[0] = {grid[0]:.3e}")
# print(f"v02 = {v02:.3e}, closest grid = {grid[idx02]:.3e} (error {abs((grid[idx02]-v02)/v02)*100:.2f}%)")
# print(f"vmax = {vmax:.3e}, edges[-1] = {edges[-1]:.3e}")

# # === Plot for verification ===
# fig, ax = plt.subplots(figsize=(7,6))
# colors = ['red' if f == '.true.' else 'blue' for f in flag_nuclei]
# ax.scatter(grid, np.arange(Nbins), c=colors, marker='o', label='Grid midpoints')
# ax.scatter(edges, np.arange(Nbins+1), c='green', marker='^', label='Edges')
# ax.plot(grid, np.arange(Nbins), 'k--', lw=0.7)

# ax.set_xscale('log')
# ax.set_xlabel("Volume [m³]", fontsize=14)
# ax.set_ylabel("Bin index", fontsize=14)
# ax.grid(True, which='both')

# def vol_to_d(vol): return ((6.0 * vol / np.pi) ** (1/3)) * 1e9
# def d_to_vol(d): return (np.pi / 6.0) * (d * 1e-9)**3

# secax = ax.secondary_xaxis('top', functions=(vol_to_d, d_to_vol))
# secax.set_xlabel("Equivalent diameter [nm]", fontsize=14)
# ax.legend()
# plt.tight_layout()
# plt.savefig("pbe_grid_resolution.png", dpi=600)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import os

# === Utility Functions ===
def es(T_c):
    return 611.2 * np.exp(17.67 * T_c / (T_c + 243.5))

def rho_air(T, P, RH):
    Rd = 8.314 / (28.96 / 1.0e3)
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / (Rd * Tv)

# === USER INPUTS ===
Nbins = 35
r01 = 2e-9
r02 = 20e-9
EI01, EI02 = 1e17, 1e12
rho1, rho2 = 1500, 1500
kappa1, kappa2 = 0.5, 0.005
output_filename = "../psr/ice_nucleating_particles.in"

v01 = (4/3) * np.pi * r01**3
v02 = (4/3) * np.pi * r02**3
vmax = v02*1e10  # User sets this!

# === Jet parameters ===
jet_V = 350.0
diameter_jet = 1.0
temp_jet = 500.0
xh2o_jet = 0.059243
f_st = 0.08957
Phi = 0.38
pressr = 20000.0
temp_air = 220.0

# === Environmental properties ===
pwater = xh2o_jet * pressr
p_water_sat_liq = np.exp(54.842763 - 6763.22 / temp_jet - 4.21 * np.log(temp_jet) +
                          0.000367 * temp_jet + np.tanh(0.0415 * (temp_jet - 218.8)) *
                          (53.878 - 1331.22 / temp_jet - 9.44523 * np.log(temp_jet) +
                           0.014025 * temp_jet))
RH = np.clip(pwater / p_water_sat_liq, 0.0, 1.0)
rho_j = rho_air(temp_jet, pressr, RH)

# === Find geometric ratio r so that grid[1]=v01 and edges[Nbins]=vmax ===
def objective(r):
    edges0 = 2 * v01 / (r * (1 + r))
    edges = edges0 * r**np.arange(Nbins + 1)
    grid = 0.5 * (edges[:-1] + edges[1:])
    return edges[-1] - vmax

sol = root_scalar(objective, bracket=[1.01, 100.0], method='bisect')
if not sol.converged:
    raise RuntimeError("Could not solve for geometric ratio r.")

r = sol.root
edges0 = 2 * v01 / (r * (1 + r))
edges = edges0 * r**np.arange(Nbins + 1)
grid = 0.5 * (edges[:-1] + edges[1:])
dv = edges[1:] - edges[:-1]

# Check centering
assert np.allclose(grid, 0.5 * (edges[:-1] + edges[1:]), rtol=1e-12), "Grid not centered between edges!"
assert np.isclose(grid[1], v01, rtol=1e-12), "Second grid point not equal to v01!"
assert np.isclose(edges[-1], vmax, rtol=1e-12), "Last edge not equal to vmax!"

# Find which grid point is closest to v02 in first 5 bins
idx02 = np.argmin(np.abs(grid[:] - v02))
idx01 = 1  # by construction: v01 is now at grid[1]

# === Assign physical properties ===
rho = np.zeros(Nbins)
kappa = np.zeros(Nbins)
n = np.zeros(Nbins)
flag_nuclei = ['.false.'] * Nbins

f_actual = Phi * f_st

if EI01 > 0:
    rho[idx01] = rho1
    kappa[idx01] = kappa1
    n[idx01] = EI01 * (rho_j * f_actual) / (1 + f_actual)
    flag_nuclei[idx01] = '.true.'

if EI02 > 0:
    rho[idx02] = rho2
    kappa[idx02] = kappa2
    n[idx02] = EI02 * (rho_j * f_actual) / (1 + f_actual)
    flag_nuclei[idx02] = '.true.'

# === Write to output file ===
os.makedirs(os.path.dirname(output_filename), exist_ok=True)
with open(output_filename, 'w') as f:
    f.write(f"{Nbins}\n")
    for i in range(Nbins):
        f.write(f"{edges[i]:.16e}\n")
        f.write(f"{grid[i]:.16e} {dv[i]:.16e} {rho[i]:.16e} {kappa[i]:.16e} {n[i]:.16e} {flag_nuclei[i]}\n")
    f.write(f"{edges[-1]:.16e}\n")

print(f"📄 File '{output_filename}' written successfully.")
print(f"v01 = {v01:.3e}, grid[1] = {grid[1]:.3e}")
print(f"v02 = {v02:.3e}, closest grid = {grid[idx02]:.3e} (error {abs((grid[idx02]-v02)/v02)*100:.2f}%)")
print(f"vmax = {vmax:.3e}, edges[-1] = {edges[-1]:.3e}")

# === Plot for verification ===
fig, ax = plt.subplots(figsize=(7,6))
colors = ['red' if f == '.true.' else 'blue' for f in flag_nuclei]
ax.scatter(grid, np.arange(Nbins), c=colors, marker='o', label='Grid midpoints')
ax.scatter(edges, np.arange(Nbins+1), c='green', marker='^', label='Edges')
ax.plot(grid, np.arange(Nbins), 'k--', lw=0.7)

ax.set_xscale('log')
ax.set_xlabel("Volume [m³]", fontsize=14)
ax.set_ylabel("Bin index", fontsize=14)
ax.grid(True, which='both')

def vol_to_d(vol): return ((6.0 * vol / np.pi) ** (1/3)) * 1e9
def d_to_vol(d): return (np.pi / 6.0) * (d * 1e-9)**3

secax = ax.secondary_xaxis('top', functions=(vol_to_d, d_to_vol))
secax.set_xlabel("Equivalent diameter [nm]", fontsize=14)
ax.legend()
plt.tight_layout()
plt.savefig("pbe_grid_resolution.png", dpi=600)
plt.show()



