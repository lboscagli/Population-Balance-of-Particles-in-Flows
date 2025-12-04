import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import os

# === USER INPUTS ===
Nbins = 55

# Log-normal mode 1
GMD1 = 2e-9       # Geometric mean diameter [m]
GSD1 = 1.3
EI1  = 1e17        # Emission index [#/kg-fuel]
rho1 = 1500        # Density [kg/m³]
kappa1 = 0.5       # Hygroscopicity

# Log-normal mode 2
GMD2 = 35e-9      # Geometric mean diameter [m]
GSD2 = 2.0
EI2  = 1e12        # Emission index [#/kg-fuel]
rho2 = 1500        # Density [kg/m³]
kappa2 = 0.005     # Hygroscopicity

# Threshold for maximum diameter for nuclei (user can set)
threshold_diameter = 45e-9  # [m]

output_filename = "../psr/ice_nucleating_particles.in"

# === Jet & Fuel Parameters ===
jet_V = 350.0
diameter_jet = 1.0
temp_jet = 500.0
xh2o_jet = 0.059243
f_st = 0.08957
Phi = 0.38
pressr = 20000.0
temp_air = 220.0

def es(T_c):
    return 611.2 * np.exp(17.67 * T_c / (T_c + 243.5))

def rho_air(T, P, RH):
    Rd = 8.314 / (28.96 / 1.0e3)
    qsat = RH * 0.622 * (es(T - 273.15) / P)
    Tv = T * (1.0 + 0.61 * qsat)
    return P / (Rd * Tv)

pwater = xh2o_jet * pressr
p_water_sat_liq = np.exp(54.842763 - 6763.22 / temp_jet - 4.21 * np.log(temp_jet) +
                         0.000367 * temp_jet + np.tanh(0.0415 * (temp_jet - 218.8)) *
                         (53.878 - 1331.22 / temp_jet - 9.44523 * np.log(temp_jet) +
                          0.014025 * temp_jet))
RH = np.clip(pwater / p_water_sat_liq, 0.0, 1.0)
rho_j = rho_air(temp_jet, pressr, RH)
f_actual = Phi * f_st

# === Scale EI to Number Concentration ===
N01 = EI1 * (rho_j * f_actual) / (1 + f_actual)
N02 = EI2 * (rho_j * f_actual) / (1 + f_actual)

# === Bin grid setup ===
vmin = (np.pi/6) * (GMD1 / GSD1)**3
vmax = (np.pi/6) * (GMD2 * GSD2 * 2e3)**3  # Covers up to 5x GMD2

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

# === Calculate number concentration for each bin ===
n = np.zeros(Nbins)
n1_arr = np.zeros(Nbins)
n2_arr = np.zeros(Nbins)
for i in range(Nbins):
    v_mid = grid[i]
    dv_bin = dv[i]
    n1_arr[i] = lognormal_n_v(v_mid, N01, GMD1, GSD1) * dv_bin
    n2_arr[i] = lognormal_n_v(v_mid, N02, GMD2, GSD2) * dv_bin
    n[i] = n1_arr[i] + n2_arr[i]

# === Apply threshold to n1 and n2 based on diameter ===
threshold_v = d_to_vol(threshold_diameter)
for i in range(Nbins):
    d_eq = vol_to_d(grid[i])
    if d_eq > threshold_diameter:
        n1_arr[i] = 0.0
        n2_arr[i] = 0.0
        n[i] = 0.0

# === Renormalize n1 and n2 so that sum matches N01 and N02 ===
n1_sum = np.sum(n1_arr)
n2_sum = np.sum(n2_arr)
scale1 = N01 / n1_sum if n1_sum > 0 else 0.0
scale2 = N02 / n2_sum if n2_sum > 0 else 0.0
n1_arr *= scale1
n2_arr *= scale2
n = n1_arr + n2_arr

# === Assign physical properties per bin, re-flag nuclei ===
rho = np.zeros(Nbins)
kappa = np.zeros(Nbins)
flag_nuclei = ['.false.'] * Nbins
for i in range(Nbins):
    if n1_arr[i] > n2_arr[i]:
        rho[i] = rho1
        kappa[i] = kappa1
    else:
        rho[i] = rho2
        kappa[i] = kappa2
    if n[i] > 0 and grid[i] <= threshold_v:
        flag_nuclei[i] = '.true.'
    else:
        flag_nuclei[i] = '.false.'

# === Write to output file ===
os.makedirs(os.path.dirname(output_filename) or '.', exist_ok=True)
with open(output_filename, 'w') as f:
    f.write(".true.\n")
    f.write(f"{Nbins}\n")
    for i in range(Nbins):
        f.write(f"{edges[i]:.16e}\n")
        f.write(f"{grid[i]:.16e} {dv[i]:.16e} {rho[i]:.16e} {kappa[i]:.16e} {n[i]:.16e} {flag_nuclei[i]}\n")
    f.write(f"{edges[-1]:.16e}\n")

print(f"📄 File '{output_filename}' written successfully.")
print(f"GMD1 = {GMD1:.3e} m, GMD2 = {GMD2:.3e} m")
print(f"EI1 = {EI1:.2e} #/kg-fuel, EI2 = {EI2:.2e} #/kg-fuel")
print(f"Scaled N01 = {N01:.2e} #/m³, N02 = {N02:.2e} #/m³")
print(f"Sum after threshold and renormalization: n1 = {np.sum(n1_arr):.2e}, n2 = {np.sum(n2_arr):.2e}, total = {np.sum(n):.2e}")

# === Plot for verification ===
fig, ax = plt.subplots(figsize=(8,6))
colors = ['red' if f == '.true.' else 'blue' for f in flag_nuclei]
ax.scatter(grid, n, c=colors, marker='o', label='Total number conc. (bin mid)')
ax.plot(grid, n, 'k--', lw=0.7, label='Sum of lognormals')
ax.plot(grid, n1_arr, 'b:', lw=1.5, label='Mode 1')
ax.plot(grid, n2_arr, 'g:', lw=1.5, label='Mode 2')

# Add bin edges: vertical green lines
for edge in edges:
    ax.axvline(edge, color='green', linestyle=':', linewidth=1, alpha=0.7)

# Add edge markers at the bottom
ymin = n[n > 0].min() * 0.1 if np.any(n > 0) else 1e5
ax.scatter(edges, [ymin]*len(edges), c='green', marker='^', label='Bin edges')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Volume [$m^3$]", fontsize=14)
ax.set_ylabel("Number concentration [$\#/m^3$]", fontsize=14)
ax.set_ylim(1,max(n))
ax.grid(True, which='both')

def vol_to_d_nm(vol): return ((6.0 * vol / np.pi) ** (1/3)) * 1e9
def d_nm_to_vol(d): return (np.pi / 6.0) * (d * 1e-9)**3

secax = ax.secondary_xaxis('top', functions=(vol_to_d_nm, d_nm_to_vol))
secax.set_xlabel("Equivalent diameter [nm]", fontsize=14)
ax.legend()
plt.tight_layout()
plt.savefig("pbe_grid_lognormal.png", dpi=600)
plt.show()