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


# paths = ['PLOTS_Time_history/EXAMPLE_Soot_only/',
#          'PLOTS_Time_history/EXAMPLE_Soot_rich/',
#           'PLOTS_Time_history/EXAMPLE_Soot_poor/']

# leg_names = [r'$\textnormal{Soot-only}$',r'$\textnormal{Soot-rich}$',r'$\textnormal{Soot-poor}$']
# leg_title = r'$\textnormal{Regime}$'
# analysis_name = 'EXAMPLE_extended_K15'

paths = ['PLOTS_Time_history/EXAMPLE_Soot_poor/',
         'PLOTS_Time_history/EXAMPLE_Soot_poor_v2/']

leg_names = [r'$\textnormal{original}$',r'$\textnormal{Revised}$']
leg_title = r'$\textnormal{Implementation}$'
analysis_name = 'EXAMPLE_extended_K15_revision'

colors = ['k','r','b']#cm.brg(np.linspace(0.1, 1, len(paths))) 
markers = ['o','s','v','^','o','*'] 
ls = ['solid','dashed','dotted']
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

S_vc = []
for i in range(0,len(paths)):
    if max(activation_binary[i])>1.0:
        Smw_activated = [val for i, val in enumerate(activation_binary[i]) if val != 0]
    else:    
        Smw_activated = S_v[(activation_binary[i]*S_v[i])>0][0]
    S_vc.append(Smw_activated)

fig,ax=plt.subplots()
for i, ln in enumerate(leg_names):
    plt.plot(time[i][::1000],mean_d_m[i][::1000]*1e9,marker=markers[i],markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=ln)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$\textnormal{time [s]}$',fontsize=18)
ax.set_ylabel(r'$d_m^{mean} [nm]$',fontsize=18)


#plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)
plt.ylim(min([min(d_mi*1E9) for d_mi in mean_d_m]),max([max(d_mi*1E9) for d_mi in mean_d_m])+0.05*max([max(d_mi*1E9) for d_mi in mean_d_m]))

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
    plt.plot(Temperature[i], P_v[i], color=colors[i],ls=ls[i],label=ln)

plt.plot(Temperature[-1], p_water_sat_liq[-1], 'k',ls='solid',alpha=0.3)#,label=r'$P_{v,sat}^{liq}$') 
plt.plot(Temperature[-1], p_water_sat_ice[-1], 'k',ls='solid',alpha=0.5)#,label=r'$P_{v,sat}^{ice}$')   

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
    plt.plot(Temperature[i], S_v[i], color=colors[i],ls=ls[i],label=ln)
    
    plt.plot(Temperature[i], min(S_vc[i])*np.ones(len(Temperature[i])), color=colors[i],ls=ls[i])
    if max(S_vc[i]) != min(S_vc[i]):
        plt.plot(Temperature[i], max(S_vc[i])*np.ones(len(Temperature[i])), color=colors[i],ls=ls[i])

ax.tick_params(labelsize=18)
ax.set_ylabel(r'$S_{v} [-]$',fontsize=18)
ax.set_xlabel(r'$T_{jet} [K]$',fontsize=18)

plt.xlim(min(Temperature[-1]),250)
plt.ylim(1,np.max(S_v)+0.05*np.max(S_v))

plt.legend(loc='best',title=leg_title,fontsize=14,title_fontsize=14,ncols=1)

plt.tight_layout()

plt.savefig(plot_dir+'/Supersaturation.png',dpi=600)


