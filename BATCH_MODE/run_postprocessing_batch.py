# -*- coding: utf-8 -*-
"""
Standalone post-processing for CPMOD batch results.

For each case in results/, generates plots and statistics in a
'post-processing' subfolder.
"""

import os
import numpy as np
from matplotlib import pyplot as plt, cm
from scipy.io import savemat

plt.rcParams['text.usetex'] = True

# ---------------- Constants ----------------
MW = {"H2O":18.01528, "CO2":44.0095, "O2":32.0, "N2":28.0134}

# ---------------- Helper functions ----------------

def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def water_flow_rate(Sv, pv_sat_liq, p_air, rho_air):
    X_H2O = Sv*pv_sat_liq/p_air
    X_O2 = (1.0 - X_H2O) * 0.233
    X_N2 = 1.0 - X_H2O - X_O2
    Mbar = X_H2O*MW["H2O"] + X_O2*MW["O2"] + X_N2*MW["N2"]
    Y_H2O = (X_H2O*MW["H2O"]) / Mbar
    C_H2O = rho_air * Y_H2O
    return C_H2O, Y_H2O

def Lw_to_dmdt(Loss_Sw, pv_sat_liq, p_air, Temperature):
    Rgas = 8314.0
    Mair_dry = 0.233*MW["O2"] + (1.0-0.233)*MW["N2"]
    Rv = Rgas/MW["H2O"]
    Rd = Rgas/Mair_dry
    rho_air_dry = p_air[-1]/Rd/Temperature[-1]
    dmh2odt_loss = Loss_Sw * rho_air_dry * Rd * pv_sat_liq / p_air / Rv
    return dmh2odt_loss

def parse_flag(val):
    val = val.strip().lower()
    return val in [".true.", "true", "1"]

def read_injection_file(filename):
    """Read ice_nucleating_particles.in file."""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    # Skip optional first line if not integer
    try:
        Nbins = int(lines[0])
        idx = 1
    except ValueError:
        Nbins = int(lines[1])
        idx = 2
    edges = np.zeros(Nbins+1)
    grid = np.zeros(Nbins)
    dv = np.zeros(Nbins)
    rho = np.zeros(Nbins)
    kappa = np.zeros(Nbins)
    n = np.zeros(Nbins)
    flag = np.zeros(Nbins, dtype=bool)
    for i in range(Nbins):
        edges[i] = float(lines[idx]); idx+=1
        parts = lines[idx].split(); idx+=1
        grid[i], dv[i], rho[i], kappa[i], n[i], f = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), parts[5]
        flag[i] = parse_flag(f)
    edges[-1] = float(lines[idx])
    return {"Nbins":Nbins, "edges":edges, "grid":grid, "dv":dv, "rho":rho, "kappa":kappa, "n":n, "flag":flag}

# ---------------- Main Post-processing ----------------

def postprocess_case(case_dir):
    """
    Process one case directory.
    All plots and .mat files are saved in 'post-processing' subfolder.
    """
    print(f"Processing case: {case_dir}")
    
    plot_dir = create_folder(os.path.join(case_dir, 'post-processing'))
    
    # Input files
    inp_file = os.path.join(case_dir, 'ice_nucleating_particles.in')
    ice_file = os.path.join(case_dir, 'ice.in')
    dic_INPs = read_injection_file(inp_file)
    
    # Find all timestep files
    timestep_files = sorted([f for f in os.listdir(case_dir) if f.startswith('psd_T')])
    
    # Temperature & main output
    temp_file = os.path.join(case_dir, 'ice_jet_temperature.out')
    if not os.path.exists(temp_file):
        print(f"Warning: {temp_file} not found. Skipping case.")
        return
    data_T = np.loadtxt(temp_file)
    time = data_T[:,0]; Temperature = data_T[:,1]; Density = data_T[:,2]
    Smw = data_T[:,7]; tau_g = data_T[:,4]; Smw_consumed = data_T[:,8]
    activation_binary = data_T[:,9]; moment_0 = data_T[:,10]; moment_1 = data_T[:,11]
    meansize = 2*((3/4/np.pi*data_T[:,12])**(1/3))
    Pressure = data_T[:,13]; particle_mass = data_T[:,14]; growth_rate_pbe = data_T[:,15]
    meansize_ice = 2*((3/4/np.pi*data_T[:,16])**(1/3))

    # Saturated water vapor
    p_water_sat_liq = np.exp(54.842763 - 6763.22 / Temperature - 4.21 * np.log(Temperature) + 0.000367 * Temperature +
                             np.tanh(0.0415*(Temperature-218.8))*(53.878 - 1331.22/Temperature - 9.44523*np.log(Temperature) + 0.014025*Temperature))
    p_water_sat_ice = np.exp(9.550426 - 5723.265/Temperature + 3.53068*np.log(Temperature) - 0.00728332*Temperature)
    P_v_consumed = Smw_consumed * p_water_sat_liq

    # Water consumption
    m_H2O_unconsumed, _ = water_flow_rate(Smw,p_water_sat_liq,Pressure,Density)
    m_H2O_consumed, _ = water_flow_rate(Smw_consumed,p_water_sat_liq,Pressure,Density)
    H2O_consumed = (m_H2O_unconsumed - m_H2O_consumed)
    H2O_consumed_perc = H2O_consumed/m_H2O_unconsumed * 100
    Sw_con_delta = Smw_consumed - Smw
    Loss_Sw = np.zeros(len(Sw_con_delta))
    Loss_Sw[1:] = (Sw_con_delta[1:] - Sw_con_delta[:-1])/(time[1:]-time[:-1])
    dmh2odt_loss = Lw_to_dmdt(Loss_Sw, p_water_sat_liq, Pressure, Temperature)

    # ---------------- Read PSD timestep files ----------------
    v_m, d_m, number_density, particle_type = [], [], [], []

    for t_filename in timestep_files:
        v_i, d_i, n_i, type_i, dv_i = [], [], [], [], []
        with open(os.path.join(case_dir, t_filename)) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 8: continue
                v_i.append(float(parts[0]))
                d_i.append(float(parts[1]))
                n_i.append(float(parts[2]))
                dv_i.append(float(parts[6]))
                type_i.append(float(parts[7]))
        v_m.append(np.array(v_i))
        d_m.append(np.array(d_i))
        number_density.append(np.array(n_i)*np.array(dv_i))
        particle_type.append(np.array(type_i))

    # ---------------- Save .mat ----------------
    dic = {'time':time,'Temperature':Temperature,'activation_binary':activation_binary,
           'P_v':P_v_consumed,'P_sat_liq':p_water_sat_liq,'P_sat_ice':p_water_sat_ice,
           'Saturation':Smw_consumed,'Mean_diameter':meansize,'Mean_diameter_ice':meansize_ice,
           'moment_0':moment_0,'moment_1':moment_1}
    savemat(os.path.join(plot_dir,'statistics.mat'), dic)

    # ---------------- Save growth rate .mat ----------------
    dic_h2o_growth = {'growth_rate_pbe': growth_rate_pbe, 'Temperature': Temperature}
    savemat(os.path.join(plot_dir,'h2o_consumption.mat'), dic_h2o_growth)



    print(f"Post-processing done for {case_dir}. Results saved in {plot_dir}.")


# ---------------- Run over all cases ----------------
if __name__ == "__main__":
    results_dir = 'results'  # <-- folder containing all case_* subdirectories
    case_dirs = [os.path.join(results_dir,f) for f in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir,f))]
    for case in case_dirs:
        postprocess_case(case)
