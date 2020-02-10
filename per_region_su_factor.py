# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 23:54:41 2020

@author: Libra
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

plot_figure=False
fstr=np.load('consec_tensor_comp_trial160_R11.npz')
[key for key in fstr.keys()]
su_factor=fstr['SU']
in_trial=fstr['in_trial']
reg_list=fstr['reg_list']


reg_set=np.array(list(set(reg_list.tolist())),dtype='|S10')
if plot_figure:
    fh=plt.figure(figsize=(11.34,40),dpi=300)
su_factors=[]
su_factors.append(list(set(reg_list.astype(np.str).tolist())))

for sub in range(11):
    
    per_region_su_factor=np.zeros_like(reg_set,dtype=np.float64)
    for i in range(reg_set.shape[0]):
        per_region_su_factor[i]=np.mean(su_factor[reg_list==reg_set[i],sub])
    su_factors.append(per_region_su_factor)    
    reg_idx=np.argsort(per_region_su_factor)
    if plot_figure:
        ax = plt.subplot(11, 1, sub+1)
        plt.bar(np.arange(reg_set.shape[0]),per_region_su_factor[reg_idx])
        ax.set_xticks(np.arange(reg_set.shape[0]))
        ax.set_xticklabels(reg_set[reg_idx].astype(np.str).tolist(),rotation=90,ha='center',va='top',fontsize=7.5)
        ax.set_xlim(-1,reg_set.shape[0])
        ax.set_ylabel('average neuron coefficient')

if plot_figure:
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fh.savefig("tca_coeff_su_factor.png", dpi=300, bbox_inches="tight")
    plt.show()
su_factors=[list(i) for i in zip(*su_factors)]
with open("tca_su_factors.csv", "w", newline="") as cf:
    cwriter = csv.writer(cf, dialect="excel")
    for row in su_factors:
        cwriter.writerow(row)