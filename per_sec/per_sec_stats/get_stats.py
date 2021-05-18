# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 08:45:24 2021

@author: Libra
"""
import numpy as np
from per_sec_stats import statsfun

def get_stats(trial_FR, trials, delay=6, debug=False, complete=False):

    '''
    Calculate FR statistics, e.g. selectivity
    '''
    if debug:
        trial_FR=trial_FR[:,:10,:]

    ### TODO: when variables are none
    if complete:
        s1sel=np.all((trials[:,7]==delay, trials[:,4]==4),axis=0)
        s2sel=np.all((trials[:,7]==delay, trials[:,4]==8),axis=0)

    else:
        s1sel=np.all((trials[:,9]>0, trials[:,8]>0, trials[:,7]==delay, trials[:,4]==4),axis=0)
        s2sel=np.all((trials[:,9]>0, trials[:,8]>0, trials[:,7]==delay, trials[:,4]==8),axis=0)


    per_sec_selectivity = np.zeros((trial_FR.shape[1],trial_FR.shape[0]))
    # per_sec_auc = np.zeros((trial_FR.shape[0], trial_FR.shape[1]))
    # per_sec_fr_s1 = np.zeros((trial_FR.shape[0], trial_FR.shape[1]))
    # per_sec_fr_s2 = np.zeros((trial_FR.shape[0], trial_FR.shape[1]))
    per_sec_wrs_p = np.zeros((trial_FR.shape[1],trial_FR.shape[0]))
    trial_counts=np.tile([np.sum(s1sel),np.sum(s2sel)],(trial_FR.shape[1],1))

    for su_idx in range(trial_FR.shape[1]):
        if su_idx % 50 == 0:
            print(su_idx)
        for bin_idx in range(0, trial_FR.shape[0]):
                mm=[np.mean(trial_FR[:,su_idx,x],axis=1) for x in [s1sel,s2sel]]

                # per_sec_fr_s1[:,su_idx]=mm[0]
                # per_sec_fr_s2[:,su_idx]=mm[1]
                per_sec_selectivity[su_idx,:]=[
                    (mm[0][x]-mm[1][x])/(mm[0][x]+mm[1][x])
                 if (mm[0][x]+mm[1][x])>0
                 else 0
                 for x in range(trial_FR.shape[0])]

                # AUC and WRS depends on external func.
                # TODO ufunc-> reduce?
                per_sec_wrs_p[su_idx,:] = [
                    statsfun.wrs(trial_FR[x,su_idx,s1sel],trial_FR[x,su_idx,s2sel])
                    if np.abs(per_sec_selectivity[su_idx,x])>0.05
                    else 1
                    for x in range(trial_FR.shape[0])]
    return {'per_sec_selectivity':per_sec_selectivity,
            # 'per_sec_fr_s1':per_sec_fr_s1,
            # 'per_sec_fr_s2':per_sec_fr_s2,
            'per_sec_wrs_p':per_sec_wrs_p,
            'trial_counts':trial_counts}
                # TODO AUC