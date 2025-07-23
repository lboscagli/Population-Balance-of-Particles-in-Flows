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

###############################################################################
""" USER SETTINGs """

path = '../pbe/'

case_name = r''#'$\rho_p [kg/m^3]$'

analysis_name = 'TEST'#'15_dd20nm_k0p005'#'N04p5e10_K96'#'_no_consumption'#'LES_XH2O_55bins_geometric_buthanol'#'TEST'#


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
for t_filename in timestep_files:
    v_m_i = []
    d_m_i = []
    dv_i = []
    particle_type_i = []
    number_density_i = []
    with open(path+t_filename) as f:
        for line in f.readlines():
            v_m_i.append(float(line.split()[0]))
            d_m_i.append(float(line.split()[1]))
            number_density_i.append(float(line.split()[2]))
            dv_i.append(float(line.split()[6]))
            particle_type_i.append(float(line.split()[7]))
        
    v_m.append(np.array(v_m_i))
    d_m.append(np.array(d_m_i))
    number_density.append(np.array(number_density_i)*np.array(dv_i))
    idx_mean = find_nearest(number_density_i, np.mean(number_density_i))
    idx_peak = np.argmax(number_density_i)
    mean_d_m.append(d_m_i[idx_mean])
    peak_d_m.append(d_m_i[idx_peak])
    Sum_nv.append(np.sum(np.array(number_density_i)*np.array(dv_i)))
    Mean_d_m.append(np.sum(np.array(d_m_i)*np.array(number_density_i)*np.array(dv_i))/np.sum(np.array(number_density_i)*np.array(dv_i)))
    particle_type.append(np.array(particle_type_i))


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

plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.ylim(1,max([max(nv) for nv in number_density]))
plt.title('Particle number density\n per cell unit volume', fontsize=14)

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/number_density.png',dpi=600)


fig,ax=plt.subplots()

plt.plot(time,Temperature,'k',ls='solid') 
plt.plot(time[::N_step_Save],Temperature[::N_step_Save],'ro',ls='none',fillstyle='none')    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$t [s]$',fontsize=18)
ax.set_ylabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(time),max(time))

#plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Temperature_history.png',dpi=600)


fig,ax=plt.subplots()

plt.plot(Temperature, p_water_sat_liq, 'k',ls='solid',label=r'$P_{v,sat}^{liq}$') 
plt.plot(Temperature, p_water_sat_ice, 'r',ls='dashed',label=r'$P_{v,sat}^{ice}$')   
plt.plot(Temperature[1:], P_v_consumed[1:], 'c',ls='-.',label=r'$\textnormal{mixing line}$')
#plt.plot(Temperature[1:], P_v[1:], 'b--',ls='-.')#,label=r'$\textnormal{mixing line}$')  

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$P_{v} [Pa]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature),260)
plt.ylim(1,max([max(p_water_sat_ice[Temperature<250]),max(p_water_sat_liq[Temperature<250])]))

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
    Smw_activated = Smw_consumed[(activation_binary*Smw_consumed)>0][0]

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

plt.savefig(plot_dir+'/Supersaturation_and_number_density.png',dpi=600)