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

path = './pbe/'

case_name = r''#'$\rho_p [kg/m^3]$'

analysis_name = 'RH1p6_55bins'


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
P_v = data_T[:,3] # water vapor partial pressure
tau_g = data_T[:,4] # growth timescale
N_step_Save = len(time)//len(Sampled_iterations)

## Compute saturated (relative to liquid) water vapour partial pressure 
p_water_sat_liq = np.exp(54.842763 - 6763.22 / Temperature - 4.21 * np.log(Temperature) + \
+ 0.000367 * Temperature + np.tanh(0.0415 * (Temperature - 218.8)) * (53.878 - 1331.22 / Temperature - 9.44523 * np.log(Temperature) + \
+ 0.014025 * Temperature))

## Compute saturated (relative to ice) water vapour partial pressure 
p_water_sat_ice = np.exp(9.550426 - 5723.265 / Temperature + 3.53068 * np.log(Temperature) + \
                  - 0.00728332 * Temperature) 

# In[1]: Read timestep directories
v_m = []
d_m = []
mean_d_m = []
peak_d_m = []
number_density = []
for t_filename in timestep_files:
    v_m_i = []
    d_m_i = []
    number_density_i = []
    with open(path+t_filename) as f:
        for line in f.readlines():
            v_m_i.append(float(line.split()[0]))
            d_m_i.append(float(line.split()[1]))
            number_density_i.append(float(line.split()[2]))
        
    v_m.append(np.array(v_m_i))
    d_m.append(np.array(d_m_i))
    number_density.append(np.array(number_density_i))
    idx_mean = find_nearest(number_density_i, np.mean(number_density_i))
    idx_peak = np.argmax(number_density_i)
    mean_d_m.append(d_m_i[idx_mean])
    peak_d_m.append(d_m_i[idx_peak])

dic = {'timesteps':timestep_files,
       'Mean_diameter':mean_d_m,
       'Peak_diameter':peak_d_m}

savemat(plot_dir+'/statistics.mat',dic)

        
    
fig,ax=plt.subplots()
for i, case in enumerate(Sampled_iterations):
    plt.semilogy(d_m[i]*1e9,number_density[i],marker='o',markersize=10,linestyle='solid',color=colors[i],fillstyle='none')#,label=case)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$d_m [nm]$',fontsize=18)
ax.set_ylabel(r'$n_v [\#/m^3]$',fontsize=18)

plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
#plt.xlim(0,2000)

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

plt.semilogy(Temperature, p_water_sat_liq, 'k',ls='solid',label=r'$P_{v,sat}^{liq}$') 
plt.semilogy(Temperature, p_water_sat_ice, 'r',ls='dashed',label=r'$P_{v,sat}^{ice}$')   
plt.semilogy(Temperature[1:], P_v[1:], 'b',ls='-.',label=r'$\textnormal{mixing line, } RH_{liq}=const$')  

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$P_{v} [Pa]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature),max(Temperature))
plt.ylim(1,max([max(p_water_sat_ice),max(p_water_sat_liq)]))

plt.legend(loc='best',fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/Saturation_partial_pressure.png',dpi=600)


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
