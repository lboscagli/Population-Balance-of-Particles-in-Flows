# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 09:11:23 2025

@author: lboscagl
"""

import os
import sys
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.io import loadmat,savemat

plt.rcParams['text.usetex'] = True

# Constants
MW = {
    "H2O": 18.01528,
    "CO2": 44.0095,
    "O2": 32.0,
    "N2": 28.0134,
}

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def create_folder(name):
    import os
    
    current_directory = os.getcwd()
    new_dir = os.path.join(current_directory, name) 
    
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    
    return new_dir

def water_flow_rate(Sv,pv_sat_liq,p_air,rho_air):
    
    #Determine mass flow rate of water vapor
    X_H2O = Sv*pv_sat_liq/p_air

    X_O2 = (1.0 - X_H2O) * 0.233
    X_N2 = (1.0 - X_H2O - X_O2) 
    Mbar = X_H2O*MW["H2O"] + X_O2*MW["O2"] + X_N2*MW["N2"]   # kg/kmol
    Y_H2O = (X_H2O * MW["H2O"]) / Mbar
    C_H2O = rho_air * Y_H2O #/ MW["H2O"]   # since V = 1 m3 
    #mdot_H2O = C_H2O      # because Q=1 -> mdot = rho*Y*Q = rho*Y
    
    return C_H2O, Y_H2O

def Lw_to_dmdt(Loss_Sw,pv_sat_liq,p_air,Temperature):
    Rgas = 8314.0 
    Mair_dry = 0.233*MW["O2"] + (1.0-0.233)*MW["N2"]
    Rv = Rgas/MW["H2O"]
    Rd = Rgas/Mair_dry
    rho_air_dry = p_air[-1]/Rd/Temperature[-1]
    dmh2odt_loss = Loss_Sw * rho_air_dry * Rd * pv_sat_liq / p_air / Rv
    
    return dmh2odt_loss


###############################################################################
""" USER SETTINGs """

path = '../pbe/'

case_name = r''#'$\rho_p [kg/m^3]$'

analysis_name = 'K15_EInv1e14_20nm_EIv1e17_2nm_Ta210'#'DEBUG_SESSION_EIvPM0e17_EInvPM1e14_K15'#'DEBUG_SESSION_EInvPM1e16_K15'#'DEBUG_SESSION_EInvPM1e17_K15'#'EXAMPLE_Soot_poor_v2'#'EXAMPLE_Soot_only'#'15_dd20nm_k0p005'#'N04p5e10_K96'#'_no_consumption'#'LES_XH2O_55bins_geometric_buthanol'#'TEST'#

NUMBER_DENSITY_BAR_PLOT = False
NUMBER_DENSITY_SCATTER_PLOT = False
SATURATION_ANIMATED_PLOTS = False

###############################################################################

plot_dir = create_folder('PLOTS_Time_history/'+analysis_name)

# In[0]: Read timestep directories
timestep_files=[]
Sampled_iterations=[]
for idx_t,tmp_file in enumerate(os.listdir(path)):
    if tmp_file.startswith('psd_T'):
        timestep_files.append(tmp_file)
        Sampled_iterations.append(float(tmp_file.replace('psd_T','').replace('.out','')))

#Reorder cronologically
timestep_files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
Sampled_iterations.sort()

colors = cm.cool(np.linspace(0, 1, len(Sampled_iterations)))

# In[0]: Read temperature profile
data_T = np.loadtxt(path+'ice_jet_temperature.out')
time = data_T[:,0]
t_step = np.linspace(0,len(time),num=len(time))
Temperature = data_T[:,1]
Density = data_T[:,2]
Smw = data_T[:,7] # saturation ratio of the mixing line assuming no particles gorwth
tau_g = data_T[:,4] # growth timescale
N_step_Save = len(time)//len(Sampled_iterations)
Smw_consumed = data_T[:,8] # saturation ratio with consumption
activation_binary = data_T[:,9]
moment_0 = data_T[:,10]
moment_1 = data_T[:,11]
meansize = 2*((3/4/np.pi*(data_T[:,12]))**(1/3))
Pressure = data_T[:,13]
particle_mass = data_T[:,14]
growth_rate_pbe = data_T[:,15]

## Compute saturated (relative to liquid) water vapour partial pressure 
p_water_sat_liq = np.exp(54.842763 - 6763.22 / Temperature - 4.21 * np.log(Temperature) + \
+ 0.000367 * Temperature + np.tanh(0.0415 * (Temperature - 218.8)) * (53.878 - 1331.22 / Temperature - 9.44523 * np.log(Temperature) + \
+ 0.014025 * Temperature))

## Compute saturated (relative to ice) water vapour partial pressure 
p_water_sat_ice = np.exp(9.550426 - 5723.265 / Temperature + 3.53068 * np.log(Temperature) + \
                  - 0.00728332 * Temperature) 

## Compute water vapor partial pressure from saturation ratio
P_v = Smw * p_water_sat_liq
P_v_consumed = Smw_consumed*p_water_sat_liq


#' Determine mass flow rate of water vapor - without consumption - this decreases simply due to dilution
m_H2O_unconsumed, Y_H2O_unconsumed = water_flow_rate(Smw,p_water_sat_liq,Pressure,Density)
mdot_H2O_unconsumed=np.zeros(len(m_H2O_unconsumed)) 
mdot_H2O_unconsumed[1:] = (m_H2O_unconsumed[1:]-m_H2O_unconsumed[:-1])/(time[1:]-time[:-1])

## Determine mass flow rate of water vapor - with consumption - this decreases due to dilution and H2O consumption due to particle growth
m_H2O_consumed, Y_H2O_consumed = water_flow_rate(Smw_consumed,p_water_sat_liq,Pressure,Density)
mdot_H2O_consumed=np.zeros(len(m_H2O_consumed)) 
mdot_H2O_consumed[1:] = (m_H2O_consumed[1:]-m_H2O_consumed[:-1])/(time[1:]-time[:-1])

#Percentage of consumed water vapor relative to the original one. This is a relative change, to eliminate the effect of dilution. So this is purely due to consumption 
H2O_consumed = (m_H2O_unconsumed - m_H2O_consumed) 
H2O_consumed_perc = H2O_consumed/m_H2O_unconsumed * 100

mdot_H2O_change = np.zeros(len(H2O_consumed))
mdot_H2O_change[1:] = (H2O_consumed[1:]-H2O_consumed[:-1])/(time[1:]-time[:-1])

##  particle mass increase due to water vapor uptake 
#M_p_change = np.zeros(len(particle_mass))
#M_p_change[1:] = (particle_mass[1:] - particle_mass[:-1])/(time[1]-time[0])*time[1:]
M_p_change = particle_mass - particle_mass[0]
mdot_p = np.zeros(len(particle_mass))
mdot_p[1:]=(particle_mass[1:]-particle_mass[:-1])/(time[1]-time[0])


#Rate of change of water vapor consumption
Sw_con_delta = Smw_consumed - Smw
Loss_Sw = np.zeros(len(Sw_con_delta))
Loss_Sw[1:] = (Sw_con_delta[1:]-Sw_con_delta[:-1])/(time[1:]-time[:-1])
dmh2odt_loss = Lw_to_dmdt(Loss_Sw, p_water_sat_liq, Pressure, Temperature)

#Plot water vapor mass flow balance - Ideally the c

# fig,ax=plt.subplots()

# plt.plot(Temperature,mdot_H2O_unconsumed,'k',label=r'w/o unconsumed')
# plt.plot(Temperature,mdot_H2O_consumed,'g-.',label=r'with consum.')
# #plt.plot(Temperature,mdot_p,'b-.',label=r'$\frac{d(\rho_p M_1)}{dt}$ - FD')
# plt.plot(Temperature,growth_rate_pbe,'orange',alpha=0.5,label=r'$\frac{d(\rho_p M_1)}{dt}$ - PBE')
# plt.plot(Temperature,mdot_H2O_consumed+growth_rate_pbe,'r:',label=r'with consum.+$\frac{d(\rho_p M_1)}{dt}$')

# ax.tick_params(labelsize=18)
# ax.set_xlabel(r'$T [K]$',fontsize=18)
# ax.set_ylabel(r'$\frac{\partial (\rho Y_{H_2O})}{\partial t} [kg/(s \, m^3)]$',fontsize=18)

# plt.xlim(min(Temperature),max(Temperature))

# plt.title('H2O mass balance', fontsize=14)

# plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

# plt.tight_layout()

# plt.savefig(plot_dir+'/H2O_balance.png',dpi=600)


fig,ax=plt.subplots()

#plt.plot(Temperature,-mdot_H2O_change,'k',label=r'consumed H2O flow rate')
plt.plot(Temperature,dmh2odt_loss,'k',alpha=0.5,label=r'consumed H2O flow rate')
#plt.plot(Temperature,mdot_p,'b:',label=r'$\frac{d(\rho_p M_1)}{dt}$ - FD')
plt.plot(Temperature,growth_rate_pbe,'r-.',label=r'$\frac{d(\rho_p M_1)}{dt}$ - PBE')
plt.plot(Temperature,growth_rate_pbe+dmh2odt_loss,'g:',label=r'sum')

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$T [K]$',fontsize=18)
ax.set_ylabel(r'$[kg/(s \, m^3)]$',fontsize=18)

plt.xlim(min(Temperature),270)#max(Temperature))

plt.title('H2O mass balance', fontsize=14)

plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/H2O_balance_v2.png',dpi=600)

fig,ax=plt.subplots()
ax1=ax.twinx()

ax.plot(Temperature,m_H2O_unconsumed,'k-.',alpha=0.2,label='w/o consumption')
ax.plot(Temperature,m_H2O_consumed,'k:',alpha=0.7,label='with consumption')

ax1.plot(Temperature,H2O_consumed_perc,'r')


ax.tick_params(labelsize=18)
ax1.tick_params(labelsize=18)
ax.set_xlabel(r'$T [K]$',fontsize=18)
ax.set_ylabel(r'$\rho Y_{H_2O} [kg/m^3]$',fontsize=18)
ax1.set_ylabel(r'$[\%]$',fontsize=18)

ax.set_xlim(min(Temperature),270)#max([max(d_m_i) for d_m_i in d_m])*1E9)
ax.set_ylim(0,m_H2O_unconsumed[Temperature<270].max())
ax1.set_ylim(-1,101)
plt.title('Consumed H2O', fontsize=14)

ax.legend(loc='center right',fontsize=12)

#Change color
ax1.tick_params(axis='y', colors='r')
ax1.yaxis.label.set_color('r')

plt.tight_layout()

plt.savefig(plot_dir+'/H2O_consumed_percentage.png',dpi=600)

fig,ax=plt.subplots()

plt.plot(Temperature,meansize*1e9,'k')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$T [K]$',fontsize=18)
ax.set_ylabel(r'$d_m^{mean} [nm]$',fontsize=18)


#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.ylim(0,max(meansize*1e9)+0.05*max(meansize*1E9))


plt.tight_layout()

plt.savefig(plot_dir+'/mean_diameter.png',dpi=600)


dic_h2o_growth = {'growth_rate_pbe':growth_rate_pbe,
                  'Temperature':Temperature}

savemat(plot_dir+'/h2o_consumption.mat',dic_h2o_growth)

# In[1]: Read PBE input file 
with open(path+'pbe.in', 'r') as file:
    #Read DNS bcset.f90 file
    lines = file.readlines()
    for line in lines:
        if 'Number of grid points' in line:#.split():
            Nbins=float(line.split()[0])
        if 'Minimum size in grid' in line:#.split(): 
            v_min=float(line.split()[0])
        if 'Maximum size in grid' in line:#.split(): 
            v_max=float(line.split()[0])

#volume of the bins assuming uniform spacing
dv_uniform = (v_max-v_min)/Nbins
    
# In[1]: Read timestep directories
v_m = []
d_m = []
mean_d_m = []
peak_d_m = []
number_density = []
Sum_nv = []
Mean_d_m = []
particle_type = []
saturation_loss = []
mass_loss = []
for t_filename in timestep_files:
    v_m_i = []
    d_m_i = []
    dv_i = []
    particle_type_i = []
    number_density_i = []
    saturation_loss_i = []
    mass_loss_i = []
    with open(path + t_filename) as f:
        for line in f:
            parts = line.split()
            # Safely parse the first 8 columns
            v_m_i.append(float(parts[0]) if len(parts) > 0 else 0.0)
            d_m_i.append(float(parts[1]) if len(parts) > 1 else 0.0)
            number_density_i.append(float(parts[2]) if len(parts) > 2 else 0.0)
            dv_i.append(float(parts[6]) if len(parts) > 6 else 0.0)
            particle_type_i.append(float(parts[7]) if len(parts) > 7 else 0.0)
            
            # Columns 8 and 9 might be missing, default to 0.0
            try:
                saturation_loss_i.append(float(parts[8]))
            except (IndexError, ValueError):
                saturation_loss_i.append(0.0)
                
            try:
                mass_loss_i.append(float(parts[9]))
            except (IndexError, ValueError):
                mass_loss_i.append(0.0)                            
               
    
    mask_inactive = np.array(particle_type_i)>0 
    v_m.append(np.array(v_m_i))
    d_m.append(np.array(d_m_i))
    number_density.append(np.array(number_density_i)*np.array(dv_i))
    idx_mean = find_nearest(number_density_i, np.mean(number_density_i))
    idx_peak = np.argmax(number_density_i)
    mean_d_m.append(d_m_i[idx_mean])
    peak_d_m.append(d_m_i[idx_peak])
    Sum_nv.append(np.sum(np.array(number_density_i)*np.array(dv_i)))
    Mean_d_m.append(((np.sum(np.array(v_m_i)[mask_inactive]*np.array(number_density_i)[mask_inactive]*np.array(dv_i)[mask_inactive])/np.sum(np.array(number_density_i)[mask_inactive]*np.array(dv_i)[mask_inactive]))*3/4/np.pi)**(1/3))
    particle_type.append(np.array(particle_type_i))
    saturation_loss.append(np.array(saturation_loss_i))
    mass_loss.append(np.array(mass_loss_i))


# dic = {'time':time[::N_step_Save],
#        'timesteps':timestep_files,
#        'Mean_diameter':mean_d_m,
#        'Peak_diameter':peak_d_m,
#        'Sum_nv':Sum_nv,
#        'Mean_d_m_nvw':Mean_d_m}

dic = {'time':time,
       'Temperature':Temperature,
       'activation_binary':activation_binary,
       'P_v':P_v_consumed,
       'P_sat_liq':p_water_sat_liq,
       'P_sat_ice':p_water_sat_ice,
       'Saturation':Smw_consumed,
       'Mean_diameter':meansize,
       'moment_0':moment_0,
       'moment_1':moment_1}

savemat(plot_dir+'/statistics.mat',dic)

        
    
fig,ax=plt.subplots()
for i, case in enumerate(Sampled_iterations):
    plt.semilogy(d_m[i]*1e9,number_density[i],marker='o',markersize=10,linestyle='solid',color=colors[i],fillstyle='none')#,label=case)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$d_m [nm]$',fontsize=18)
ax.set_ylabel(r'$n_v [\#/m^3]$',fontsize=18)

plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,10000)#max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.ylim(1,max([max(nv) for nv in number_density]))
plt.title('Particle number density\n per cell unit volume', fontsize=14)

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/number_density.png',dpi=600)


fig,ax=plt.subplots()

plt.plot(time,Temperature,'k',ls='solid') 
#plt.plot(time[::N_step_Save],Temperature[::N_step_Save],'ro',ls='none',fillstyle='none')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$t [s]$',fontsize=18)
ax.set_ylabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(time),max(time))

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Temperature_history.png',dpi=600)


fig,ax=plt.subplots()

plt.plot(Temperature[Temperature<270], p_water_sat_liq[Temperature<270], 'k',alpha=0.7,ls=':',label=r'$P_{v,sat}^{liq}$') 
plt.plot(Temperature[Temperature<270], p_water_sat_ice[Temperature<270], 'k',alpha=0.4,ls=':',label=r'$P_{v,sat}^{ice}$')   
plt.plot(Temperature[Temperature<270], P_v[Temperature<270], 'c',ls='-.',label=r'$\textnormal{mixing line}$')
plt.plot(Temperature[Temperature<270], P_v_consumed[Temperature<270], 'r',ls='-.',label=r'$\textnormal{mixing line}$')

#plt.plot(Temperature[1:], Smw_consumed[1:]*p_water_sat_liq[1:], 'c',ls='-.',label=r'$\textnormal{mixing line}$')
#plt.plot(Temperature[1:], P_v[1:], 'b--',ls='-.')#,label=r'$\textnormal{mixing line}$')  

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$P_{v} [Pa]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature),270)
plt.ylim(1,max([max(P_v[Temperature<270]),max(P_v[Temperature<270])]))

plt.legend(loc='best',fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Saturation_partial_pressure.png',dpi=600)


fig,ax=plt.subplots()

plt.plot(Temperature, Smw - 1, 'b',ls='solid',label=r'$\textnormal{no particles}$') 
plt.plot(Temperature, Smw_consumed -1, 'c',ls='dashed',label=r'$\textnormal{with consumption}$')   

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$s_{w} [-]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature),240)
plt.ylim(0,(max(Smw[Temperature<250])-1) + 0.1*(max(Smw[Temperature<250])-1))

plt.legend(loc='best',fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Supersaturation_consumption.png',dpi=600)

fig,ax=plt.subplots()

idx_growth = tau_g>0
plt.plot(time[idx_growth],tau_g[idx_growth],'k',ls='solid')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$t [s]$',fontsize=18)
ax.set_ylabel(r'$\tau_g [s]$',fontsize=18)

plt.xlim(min(time),max(time))

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Growth_timescale.png',dpi=600)


#Determine mixing timescale
tau_m = time[Temperature<Temperature[0]][0]

fig,ax=plt.subplots()

idx_growth = tau_g>0
plt.semilogy(time[idx_growth],tau_m/tau_g[idx_growth],'k',ls='solid')
plt.semilogy(time,np.ones(len(time)),'k',ls='dotted')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$t [s]$',fontsize=18)
ax.set_ylabel(r'$Da = \frac{\tau_m}{\tau_g} [-]$',fontsize=18)

plt.xlim(min(time),max(time))

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Damkohler_number.png',dpi=600)


if max(activation_binary)>1.0:
    Smw_activated = [val for i, val in enumerate(activation_binary) if val != 0]
else:    
    Smw_activated = 1.0#Smw_consumed[(activation_binary*Smw_consumed)>0][0]

fig,ax=plt.subplots()
 
ax1=ax.twinx()

ax.semilogy(time,moment_0,'k')#,label=r'$0^{th}-moment$')#meansize*1E9,'k')
ax1.plot(time,Smw_consumed,'r')#,label=r'$S_{mw}$')

if max(activation_binary)>1.0:
    ax1.plot(time,max(Smw_activated)*np.ones(len(time)),'r',ls='dashed',label=r'$S_{v,c}-vPM$')
    ax1.plot(time,min(Smw_activated)*np.ones(len(time)),'r',ls='dotted',label=r'$S_{v,c}-nvPM$')
else:
    ax1.plot(time,Smw_activated*np.ones(len(time)),'r',ls='dotted',label=r'$S_{v,c}$')

ax.set_xlabel(r'$t [s]$',fontsize=18)
ax.set_ylabel(r'$0^{th}-moment$',fontsize=18)
ax1.set_ylabel(r'$S_{mw}$',fontsize=18)
ax.set_xlim(0,max(time))
#ax.set_ylim(min(moment_0)-0.12*min(moment_0),max(moment_0)+0.12*max(moment_0))
ax1.plot(time,np.ones(len(time)),'r-.',label=r'$S_{v,liq}$')
ax1.set_ylim(0,2.5)
ax.tick_params(labelsize=18)
ax1.tick_params(labelsize=18)

ax1.legend(loc='lower right',fontsize=12)
#ax.legend(loc='center right',fontsize=12)

#Change color
ax1.tick_params(axis='y', colors='r')
ax1.yaxis.label.set_color('r')

plt.tight_layout()

plt.savefig(plot_dir+'/Supersaturation_and_number_density.png',dpi=600)


if NUMBER_DENSITY_SCATTER_PLOT:
    plot_hist_dir = create_folder(plot_dir+'/Scatter_plots')
 
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    
    # Define your custom labels (adjust as needed)
    category_labels = [r'$\textnormal{Inactive}$', r'$\textnormal{Water}$', r'$\textnormal{Ice}$']
    category_levels = [0, 1, 2]  # Must match values in particle_type[i]
    
    # Use a discrete colormap with 3 colors
    cmap = cm.get_cmap('Blues', 3)  # 3 discrete shades of blue
    norm = mcolors.BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5], ncolors=3)
    
    for i, case in enumerate(Sampled_iterations):
        # fig, ax = plt.subplots()
        
        # x_vals = d_m[i] * 1e9  # Convert to nm
        # y_vals = number_density[i]
        # z_vals = particle_type[i]  # Should contain only 0, 1, 2
        
        # # Scatter plot with discrete colors
        # sc = ax.scatter(x_vals, y_vals, c=z_vals, cmap=cmap, norm=norm,
        #                 s=100, edgecolors='k')
        
        # ax.set_yscale('log')  # Log scale for y-axis
        # ax.set_xscale('log')
        
        # ax.tick_params(labelsize=18)
        # ax.set_xlabel(r'$d_m\ [nm]$', fontsize=18)
        # ax.set_ylabel(r'$n_v\ [\#/m^3]$', fontsize=18)
        
        # ax.set_xlim(1, 1E4)#max([max(d_m_i) for d_m_i in d_m]) * 1e9)
        # ax.set_ylim(1, max([max(nv) for nv in number_density])+10*max([max(nv) for nv in number_density]))
        
        # plt.title(r'$T_{jet} \textnormal{[K]}=$'+str("${:.1f}$".format(Temperature[int(case)-1]))+'\n'+r'$t \textnormal{[s]}=$'+str("${:.4f}$".format(time[int(case)-1])), fontsize=14)
        
        # # Colorbar with custom category labels
        # cbar = plt.colorbar(sc, ax=ax, ticks=category_levels)
        # cbar.ax.set_yticklabels(category_labels)
        # cbar.set_label(r'$\textnormal{Particle type}$', fontsize=18)
        
        # plt.tight_layout()
        # plt.savefig(plot_hist_dir+'/number_density_it-'+str(int(case))+'.png', dpi=600)
        # plt.close()
        
        fig, ax = plt.subplots()
        
        x_vals = d_m[i] * 1e9  # Convert to nm
        y_vals = number_density[i]
        z_vals = particle_type[i]  # Should contain only 0, 1, 2
        
        # Scatter plot with discrete colors
        sc = ax.scatter(x_vals, y_vals, c=z_vals, cmap=cmap, norm=norm,
                        s=100, edgecolors='k')
        
        ax.set_yscale('log')  # Log scale for y-axis
        ax.set_xscale('log')
        
        ax.tick_params(labelsize=18)
        ax.set_xlabel(r'$d_m\ [nm]$', fontsize=18)
        ax.set_ylabel(r'$n_v\ [\#/m^3]$', fontsize=18)
        
        ax.set_xlim(1, 1E4)#max([max(d_m_i) for d_m_i in d_m]) * 1e9)
        ax.set_ylim(1, max([max(nv) for nv in number_density])+10*max([max(nv) for nv in number_density]))
        
        plt.title(r'$S_{mw} \textnormal{[-]}=$' + str("${:.2f}$".format(Smw_consumed[int(case) - 1])) + ',\t \t' +
                r'$T_{jet} \textnormal{[K]}=$'+str("${:.1f}$".format(Temperature[int(case)-1]))+'\n'+
                  r'$t \textnormal{[s]}=$'+str("${:.4f}$".format(time[int(case)-1])), fontsize=14)
        
        # Colorbar with custom category labels
        cbar = plt.colorbar(sc, ax=ax, ticks=category_levels)
        cbar.ax.set_yticklabels(category_labels)
        cbar.set_label(r'$\textnormal{Particle status}$', fontsize=18)
        
        plt.tight_layout()
        plt.savefig(plot_hist_dir+'/number_density_it-'+str(int(case))+'.png', dpi=600)
        plt.close()        

if NUMBER_DENSITY_BAR_PLOT:
    plot_hist_dir = create_folder(plot_dir+'/Bar_plots')

    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    # Define your custom labels (adjust as needed)
    category_labels = [r'$\textnormal{Inactive}$',
                       r'$\textnormal{Water}$',
                       r'$\textnormal{Ice}$']
    category_levels = [0, 1, 2]  # Must match values in particle_type[i]

    # Use discrete colormap with 3 colors
    category_colors = cm.get_cmap('Blues', 3)([0, 1, 2])  # Get 3 distinct shades

    for i, case in enumerate(Sampled_iterations):
        fig, ax = plt.subplots()

        x_vals = d_m[i] * 1e9  # Convert to nm
        y_vals = number_density[i]
        z_vals = particle_type[i]  # 0, 1, or 2

        # Sort by x for neater plotting (optional)
        sorted_indices = np.argsort(x_vals)
        x_vals = x_vals[sorted_indices]
        y_vals = y_vals[sorted_indices]
        z_vals = z_vals[sorted_indices]

        # Map each bar's color from particle type
        bar_colors = [category_colors[int(t)] for t in z_vals]

        # Plot bars instead of scatter
        ax.bar(x_vals, y_vals, width=10, color=bar_colors, edgecolor='k', align='center')

        ax.set_yscale('log')
        ax.set_xscale('log')

        ax.tick_params(labelsize=18)
        ax.set_xlabel(r'$d_m\ [nm]$', fontsize=18)
        ax.set_ylabel(r'$n_v\ [\#/m^3]$', fontsize=18)

        ax.set_xlim(1, 1E4)
        ax.set_ylim(1, max([max(nv) for nv in number_density]) + 10 * max([max(nv) for nv in number_density]))

        plt.title(
            r'$S_{mw} \textnormal{[-]}=$' + str("${:.2f}$".format(Smw_consumed[int(case) - 1])) + ',\t \t' +
            r'$T_{jet} \textnormal{[K]}=$' + str("${:.1f}$".format(Temperature[int(case) - 1])) + '\n' +
            r'$t \textnormal{[s]}=$' + str("${:.4f}$".format(time[int(case) - 1])),
            fontsize=14
        )

        # Custom legend instead of colorbar
        handles = [plt.Rectangle((0, 0), 1, 1, color=category_colors[j], label=category_labels[j])
                   for j in range(3)]
        ax.legend(handles=handles, title=r'$\textnormal{Particle status}$', fontsize=14, title_fontsize=14)

        plt.tight_layout()
        plt.savefig(plot_hist_dir + '/number_density_it-' + str(int(case)) + '.png', dpi=600)
        plt.close()


if SATURATION_ANIMATED_PLOTS:
    plot_hist_dir = create_folder(plot_dir+'/Saturation_plots')


    for i, case in enumerate(Sampled_iterations):
        fig,ax=plt.subplots()
         
        ax1=ax.twinx()
        
        ax.semilogy(time,moment_0,'k')#,label=r'$0^{th}-moment$')#meansize*1E9,'k')
        ax.semilogy(time[int(case) - 1],moment_0[int(case) - 1],'ko',fillstyle='none',markersize=10)#,label=r'$0^{th}-moment$')#meansize*1E9,'k')
        ax1.plot(time,Smw_consumed,'r')#,label=r'$S_{mw}$')
        ax1.plot(time[int(case) - 1],Smw_consumed[int(case) - 1],'rv',fillstyle='none',markersize=10)#,label=r'$S_{mw}$')
        
        if max(activation_binary)>1.0:
            if max(Smw_activated) != min(Smw_activated):
                ax1.plot(time,max(Smw_activated)*np.ones(len(time)),'r',ls='dashed',label=r'$S_{v,c}-vPM$')
            ax1.plot(time,min(Smw_activated)*np.ones(len(time)),'r',ls='dotted',label=r'$S_{v,c}-nvPM$')
        else:
            ax1.plot(time,Smw_activated*np.ones(len(time)),'r',ls='dotted',label=r'$S_{v,c}$')
        
        ax.set_xlabel(r'$t [s]$',fontsize=18)
        ax.set_ylabel(r'$0^{th}-moment$',fontsize=18)
        ax1.set_ylabel(r'$S_{mw}$',fontsize=18)
        ax.set_xlim(0,max(time))
        ax.set_ylim(min(moment_0)-0.12*min(moment_0),max(moment_0)+0.12*max(moment_0))
        ax1.plot(time,np.ones(len(time)),'r-.',label=r'$S_{v,liq}$')
        ax1.set_ylim(0,1.5)
        ax.tick_params(labelsize=18)
        ax1.tick_params(labelsize=18)
        
        ax1.legend(loc='lower right',fontsize=12)
        #ax.legend(loc='center right',fontsize=12)
        
        #Change color
        ax1.tick_params(axis='y', colors='r')
        ax1.yaxis.label.set_color('r')
        
        plt.tight_layout()
        
        plt.savefig(plot_hist_dir + '/Supersaturation_and_number_density_it-' + str(int(case)) + '.png',dpi=600)
        plt.close()
        


##Plot saturation loss as a function of temperature for various particle sizes
saturation_loss = np.array(saturation_loss)

#Filter out zero arrays (no growth)
saturation_loss_filtered = []
particle_diameter_filtered = []
for i in range(0,saturation_loss.shape[1]):
    if any(saturation_loss[:,i]>0):
        saturation_loss_filtered.append(saturation_loss[:,i])
        particle_diameter_filtered.append(d_m_i[i])

saturation_loss_filtered = np.array(saturation_loss_filtered).T
particle_diameter_filtered = np.array(particle_diameter_filtered)

#Reduce array size for plotting
bins_interval = 1
saturation_loss_selected = saturation_loss_filtered[:,::bins_interval]
particle_diameter_selected = particle_diameter_filtered[::bins_interval]

colors_bin = cm.seismic(np.linspace(0, 1, num=saturation_loss_selected.shape[1]))
Temperature_plot = [Temperature[int(i)-1] for i in Sampled_iterations]        


fig,ax=plt.subplots()
for i,dm in enumerate(particle_diameter_selected):
    if i==0 or i==len(particle_diameter_selected)-1:
        plt.semilogy(Temperature_plot,saturation_loss_selected[:,i],linestyle='solid',color=colors_bin[i],fillstyle='none',label="${:.1f}$".format(dm*1e9))    
    else:
        plt.semilogy(Temperature_plot,saturation_loss_selected[:,i],linestyle='solid',color=colors_bin[i],fillstyle='none')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$T \; [K]$',fontsize=18)
ax.set_ylabel(r'$\mathcal{L}_{w,i} \; [-]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,10000)#max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.xlim(220,270)
#plt.title('Particle number density\n per cell unit volume', fontsize=14)

plt.legend(loc='best',title=r'$d_m [nm]$',fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/saturation_loss.png',dpi=600)

###############################################################################
##Plot mass loss as a function of temperature for various particle sizes
mass_loss = np.array(mass_loss)

#Filter out zero arrays (no growth)
mass_loss_filtered = []
particle_diameter_filtered = []
for i in range(0,mass_loss.shape[1]):
    if any(mass_loss[:,i]>0):
        mass_loss_filtered.append(mass_loss[:,i])
        particle_diameter_filtered.append(d_m_i[i])

mass_loss_filtered = np.array(mass_loss_filtered).T
particle_diameter_filtered = np.array(particle_diameter_filtered)

#Reduce array size for plotting
bins_interval = 1
mass_loss_selected = mass_loss_filtered[:,::bins_interval]
particle_diameter_selected = particle_diameter_filtered[::bins_interval]

colors_bin = cm.seismic(np.linspace(0, 1, num=mass_loss_selected.shape[1]))
Temperature_plot = [Temperature[int(i)-1] for i in Sampled_iterations]        


fig,ax=plt.subplots()
for i,dm in enumerate(particle_diameter_selected):
    if i==0 or i==len(particle_diameter_selected)-1:
        plt.semilogy(Temperature_plot,mass_loss_selected[:,i],linestyle='solid',color=colors_bin[i],fillstyle='none',label="${:.1f}$".format(dm*1e9))    
    else:
        plt.semilogy(Temperature_plot,mass_loss_selected[:,i],linestyle='solid',color=colors_bin[i],fillstyle='none')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$T \; [K]$',fontsize=18)
ax.set_ylabel(r'$|\overline{\dot{\omega}_{V,i}}| \; [kg/(m^3s)]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,10000)#max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.xlim(220,270)
plt.ylim(1e-20,1)
#plt.title('Particle number density\n per cell unit volume', fontsize=14)

plt.legend(loc='best',title=r'$d_m [nm]$',fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/mass_loss.png',dpi=600)


fig,ax=plt.subplots()

plt.semilogy(Temperature_plot,np.sum(mass_loss_selected,axis=1)/MW["H2O"],ls='solid',fillstyle='none')  
plt.semilogy(Temperature,growth_rate_pbe/MW["H2O"],'r-.')  

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$T \; [K]$',fontsize=18)
ax.set_ylabel(r'$|\overline{\dot{\omega}_{V,i}}/MW_{H_2O}| \; [kmol/(m^3s)]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,10000)#max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.xlim(220,270)
plt.ylim(1e-20,1)
#plt.title('Particle number density\n per cell unit volume', fontsize=14)

#plt.legend(loc='best',title=r'$d_m [nm]$',fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/water_vapor_loss.png',dpi=600)

    
    
