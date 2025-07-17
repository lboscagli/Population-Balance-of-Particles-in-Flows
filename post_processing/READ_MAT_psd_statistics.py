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

paths = ['PLOTS_Time_history/N04p5e12_Depositional/',
          'PLOTS_Time_history/N04p5e13_Depositional/',
          'PLOTS_Time_history/N04p5e13_Depositional_no_consumption/',
          'PLOTS_Time_history/N04p5e14_Depositional/']

leg_names = ['$4.5e12$','$4.5e13$','$4.5e13$-no cons.','$4.5e14$']
leg_title = r'$n_0 [\#/m^3]$'
analysis_name = 'CONSUMPTION_EFFECT'

colors = cm.seismic(np.linspace(0.1, 1, len(paths))) 
markers = ['o','s','v','^','o','*'] 

###############################################################################

plot_dir = create_folder('PLOTS_'+analysis_name)


# In[1]: Read timestep directories
mean_d_m = []
peak_d_m = []
sum_nv = []
mean_d_m_nvw = []
time = []
for path in paths:
    data = loadmat(path+'/statistics.mat')
    mean_d_m.append(data['Mean_diameter'][0])
    peak_d_m.append(data['Peak_diameter'][0])
    sum_nv.append(data['Sum_nv'][0])
    mean_d_m_nvw.append(data['Mean_d_m_nvw'][0])
    tmp_time = data['timesteps']
    tmp_time = [float(tt.replace('psd_T','').replace('.out','')) for tt in tmp_time]
    time.append(tmp_time)

      
    
# fig,ax=plt.subplots()
# for i, ln in enumerate(leg_names):
#     plt.plot(time[i],mean_d_m[i]*1e9,marker='o',markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

# ax.tick_params(labelsize=18)
# ax.set_xlabel(r'$\textnormal{timestep [-]}$',fontsize=18)
# ax.set_ylabel(r'$d_m^{mean} [nm]$',fontsize=18)

# #plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
# #plt.xlim(0,2000)

# plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

# plt.tight_layout()

# plt.savefig(plot_dir+'/mean_diameter.png',dpi=600)

fig,ax=plt.subplots()
for i, ln in enumerate(leg_names):
    plt.plot(time[i],mean_d_m_nvw[i]*1e9,marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$\textnormal{timestep [-]}$',fontsize=18)
ax.set_ylabel(r'$d_m^{mean} [nm]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
#plt.ylim(0,1000)

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

plt.tight_layout()

plt.savefig(plot_dir+'/mean_diameter_nvw.png',dpi=600)

fig,ax=plt.subplots()
for i, ln in enumerate(leg_names):
    plt.semilogy(time[i],sum_nv[i],marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$\textnormal{timestep [-]}$',fontsize=18)
ax.set_ylabel(r'$n [\#/m^3]$',fontsize=18)

#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
#plt.xlim(0,2000)

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

plt.tight_layout()

plt.savefig(plot_dir+'/sum_nv.png',dpi=600)

# fig,ax=plt.subplots()
# for i, ln in enumerate(leg_names):
#     plt.plot(time[i],peak_d_m[i]*1e9,marker='o',markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

# ax.tick_params(labelsize=18)
# ax.set_xlabel(r'$\textnormal{timestep [-]}$',fontsize=18)
# ax.set_ylabel(r'$d_m^{peak} [nm]$',fontsize=18)

# #plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
# #plt.xlim(0,2000)

# plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=2)

# plt.tight_layout()

# plt.savefig(plot_dir+'/peak_diameter.png',dpi=600)
