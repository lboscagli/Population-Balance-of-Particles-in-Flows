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

def create_folder(name):
    import os
    
    current_directory = os.getcwd()
    new_dir = os.path.join(current_directory, name) 
    
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    
    return new_dir

###############################################################################
""" USER SETTINGs """

paths = ['./pbe/psd_RH1p4.out','./pbe/psd_RH1p5.out','./pbe/psd_RH1p6.out']
cases = ['1.4','1.5','1.6']
colors = cm.cool(np.linspace(0, 1, len(cases)))
case_name = 'RH [-]'

analysis_name = 'RH_Effect_Tj234_const_Time0p0155'

# paths = ['./pbe/psd_25bins.out','./pbe/psd_35bins.out','./pbe/psd_45bins.out','./pbe/psd_55bins.out','./pbe/psd_65bins.out','./pbe/psd_100bins.out','./pbe/psd_150bins.out']
# cases = ['25','35','45','55','65','100','150']
# colors = cm.cool(np.linspace(0, 1, len(cases)))
# case_name = '# bins'

# analysis_name = 'discretization_RH1p0'


# paths = ['./pbe/psd_55bins_T1.out','./pbe/psd_55bins_T2.out','./pbe/psd_55bins_T3.out','./pbe/psd_55bins_T4.out']
# cases = ['1','2','3','4']
# colors = cm.cool(np.linspace(0, 1, len(cases)))
# case_name = 't [s]'

# analysis_name = 'RH1p0_55bins_time_evolution'


# paths = ['./pbe/psd_55bins_T4_RH1p0.out','./pbe/psd_55bins_T4_RH1p1.out','./pbe/psd_55bins_T4_RH1p2.out']
# cases = ['1.0','1.1','1.2']
# colors = cm.cool(np.linspace(0, 1, len(cases)))
# case_name = 'RH [-]'

# analysis_name = 'T4_55bins_RH_effect'


# paths = ['./pbe/psd_55bins_T4_RH1p0_alpha0p1.out','./pbe/psd_55bins_T4_RH1p0_alpha0p2.out','./pbe/psd_55bins_T4_RH1p0_alpha0p3.out']
# cases = ['0.1 (Karcher et al, 1996)','0.2','0.3']
# colors = cm.cool(np.linspace(0, 1, len(cases)))
# case_name = r'$\alpha$ [-]'

# analysis_name = 'T4_RH1p0_55bins_alpha_effect'

###############################################################################

plot_dir = create_folder('PLOTS/'+analysis_name)

v_m = []
d_m = []
number_density = []
for path in paths:
    v_m_i = []
    d_m_i = []
    number_density_i = []
    with open(path) as f:
        for line in f.readlines():
            v_m_i.append(float(line.split()[0]))
            d_m_i.append(float(line.split()[1]))
            number_density_i.append(float(line.split()[2]))
        
    v_m.append(np.array(v_m_i))
    d_m.append(np.array(d_m_i))
    number_density.append(np.array(number_density_i))
        
    
fig,ax=plt.subplots()
for i, case in enumerate(cases):
    plt.semilogy(d_m[i]*1e9,number_density[i],marker='o',markersize=10,linestyle='solid',color=colors[i],fillstyle='none',label=case)    

ax.tick_params(labelsize=18)
ax.set_xlabel(r'$d_m [nm]$',fontsize=18)
ax.set_ylabel(r'$n_v [\#/m^3]$',fontsize=18)

plt.xlim(min([min(d_m_i) for d_m_i in d_m])*1E9,max([max(d_m_i) for d_m_i in d_m])*1E9)

plt.legend(loc='best',title=case_name,fontsize=14,title_fontsize=14)

plt.tight_layout()

plt.savefig(plot_dir+'/number_density.png',dpi=600)