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

# paths = ['PLOTS_Time_history/RH1p6_25bins/',
#          'PLOTS_Time_history/RH1p6_35bins/',
#          'PLOTS_Time_history/RH1p6_45bins/',
#          'PLOTS_Time_history/RH1p6_55bins/',
#          'PLOTS_Time_history/RH1p6_100bins/',
#          'PLOTS_Time_history/RH1p6_200bins/']

# leg_names = ['$25$','$35$','$45$','$55$','$100$','$200$']
# leg_title = r'$\textnormal{\# bins}$'
# analysis_name = 'GRID_CONVERGENCE'

# paths = ['PLOTS_Time_history/LES_XH2O_35bins_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_45bins_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_55bins_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_85bins_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_100bins_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_200bins_buthanol/']

# leg_names = ['$35$','$45$','$55$','$85$','$100$','$200$']
# leg_title = r'$\textnormal{\# bins}$'
# analysis_name = 'GRID_CONVERGENCE_buthanol'

# paths = ['PLOTS_Time_history/LES_XH2O_35bins_uniform_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_55bins_uniform_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_100bins_uniform_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_200bins_uniform_buthanol/']

# leg_names = ['$35$','$55$','$100$','$200$']
# leg_title = r'$\textnormal{\# bins}$'
# analysis_name = 'GRID_CONVERGENCE_uniform_buthanol'

# paths = ['PLOTS_Time_history/LES_XH2O_35bins_geom_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_55bins_geom_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_100bins_geom_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_200bins_geom_buthanol/']

# leg_names = ['$35$','$55$','$100$','$200$']
# leg_title = r'$\textnormal{\# bins}$'
# analysis_name = 'GRID_CONVERGENCE_geometric_buthanol'

# paths = ['PLOTS_Time_history/LES_XH2O_35bins_composite_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_55bins_composite_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_100bins_composite_buthanol/',
#           'PLOTS_Time_history/LES_XH2O_200bins_composite_buthanol/']

# leg_names = ['$35$','$55$','$100$','$200$']
# leg_title = r'$\textnormal{\# bins}$'
# analysis_name = 'GRID_CONVERGENCE_composite_buthanol'

# paths = ['PLOTS_Time_history/N04p5e9_K96/',
#           'PLOTS_Time_history/N04p5e10_K96/',
#           'PLOTS_Time_history/N04p5e10_K96_no_consumption/',
#           'PLOTS_Time_history/N04p5e11_K96/']

# leg_names = ['$4.5e9$','$4.5e10$','$4.5e10$-no cons.','$4.5e11$']
# leg_title = r'$n_0 [\#/m^3]$'
# analysis_name = 'CONSUMPTION_EFFECT_K96'

# paths = ['PLOTS_Time_history/N04p5e9_K15_dd40nm_k0p005/',
#           'PLOTS_Time_history/N04p5e10_K15_dd40nm_k0p005/',
#           'PLOTS_Time_history/N04p5e11_K15_dd40nm_k0p005/',
#           'PLOTS_Time_history/N04p5e12_K15_dd40nm_k0p005/']

# leg_names = ['$4.5e9$','$4.5e10$','$4.5e11$','$4.5e12$']
# leg_title = r'$n_0 [\#/m^3]$'
# analysis_name = 'CONSUMPTION_EFFECT_K15_Dd40_K0p005'

# paths = ['PLOTS_Time_history/N04p5e12_K15_dd10nm_k0p005/',
#           'PLOTS_Time_history/N04p5e12_K15_dd20nm_k0p005/',
#           'PLOTS_Time_history/N04p5e12_K15_dd40nm_k0p005/']

# leg_names = ['$[10,0.005]$','$[20,0.005]$','$[40,0.005]$']
# leg_title = r'$[d_d,\kappa], [nm,-]$'
# analysis_name = 'CURVATURE_EFFECT_K15_n4P5e12_K0p005'

paths = ['PLOTS_Time_history/N04p5e12_K96/',
          'PLOTS_Time_history/N04p5e12_K15_dd40nm_k0p005/']

leg_names = ['$K96$','$K15$']
leg_title = r'$\textnormal{Ice kinetic model}$'
analysis_name = 'K96vsK15_n4P5e12_K0p005_dd40nm'

colors = cm.seismic(np.linspace(0.1, 1, len(paths))) 
markers = ['o','s','v','^','o','*'] 

###############################################################################

plot_dir = create_folder('PLOTS_'+analysis_name)


# In[1]: Read timestep directories
mean_d_m = []
sum_nv = []
time = []
Temperature = []
P_v = []
p_water_sat_liq = []
p_water_sat_ice = []
activation_binary = []
for path in paths:
    data = loadmat(path+'/statistics.mat')
    mean_d_m.append(data['Mean_diameter'][0])
    sum_nv.append(data['moment_0'][0])
    time.append(data['time'][0])
    Temperature.append(data['Temperature'][0])
    P_v.append(data['P_v'][0])
    p_water_sat_liq.append(data['P_sat_liq'][0])
    p_water_sat_ice.append(data['P_sat_ice'][0])
    activation_binary.append(data['activation_binary'][0])

#Compute saturation ratio and critical saturation ratio
S_v = [(P_v[i]/p_water_sat_liq[i]) for i in range(0,len(paths))]
S_vc = [S_v[i][(activation_binary[i]*S_v[i])>0][0] for i in range(0,len(paths))]

fig,ax=plt.subplots()
for i, ln in enumerate(leg_names):
    plt.plot(time[i][::1000],mean_d_m[i][::1000]*1e9,marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$\textnormal{time [s]}$',fontsize=18)
ax.set_ylabel(r'$d_m^{mean} [nm]$',fontsize=18)


#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.ylim(min([min(d_mi*1E9) for d_mi in mean_d_m]),max([max(d_mi*1E9) for d_mi in mean_d_m]))

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

plt.tight_layout()

plt.savefig(plot_dir+'/mean_diameter.png',dpi=600)

fig,ax=plt.subplots()
for i, ln in enumerate(leg_names):
    plt.semilogy(time[i][::1000],sum_nv[i][::1000],marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$\textnormal{time [s]}$',fontsize=18)
ax.set_ylabel(r'$n [\#/m^3]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
#plt.xlim(0,2000)

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

plt.tight_layout()

plt.savefig(plot_dir+'/sum_nv.png',dpi=600)


fig,ax=plt.subplots()

for i, ln in enumerate(leg_names):
    #plt.semilogy(time[i][::1000],sum_nv[i][::1000],marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    
    plt.plot(Temperature[i], P_v[i], color=colors[i],ls='solid',label=ln)

plt.plot(Temperature[-1], p_water_sat_liq[-1], 'k',ls='dotted',alpha=0.3)#,label=r'$P_{v,sat}^{liq}$') 
plt.plot(Temperature[-1], p_water_sat_ice[-1], 'k',ls='dotted',alpha=0.5)#,label=r'$P_{v,sat}^{ice}$')   

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$P_{v} [Pa]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature[-1]),260)
plt.ylim(1,max([max(p_water_sat_ice[-1][Temperature[-1]<250]),max(p_water_sat_liq[-1][Temperature[-1]<250])]))

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=1)

plt.tight_layout()

plt.savefig(plot_dir+'/Saturation_partial_pressure.png',dpi=600)


fig,ax=plt.subplots()


for i, ln in enumerate(leg_names):
    plt.plot(Temperature[i], S_v[i], color=colors[i],ls='solid',label=ln)
    plt.plot(Temperature[i], S_vc[i]*np.ones(len(Temperature[i])), color=colors[i],ls='dotted')


ax.tick_params(labelsize=18)
ax.set_ylabel(r'$S_{v} [-]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature[-1]),250)
plt.ylim(1,np.max(S_v)+0.05*np.max(S_v))

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=1)

plt.tight_layout()

plt.savefig(plot_dir+'/Supersaturation.png',dpi=600)


