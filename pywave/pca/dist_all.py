# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:45:14 2021

@author: Libra
"""

import numpy as np
from diststats import Dist_stats

def dist_all(repeats=2):
    fstr = np.load(f'PCA_comp_delay_diff.npz')
    comp_all = fstr['comp_all']
    comp_sust = fstr['comp_sust']
    comp_trans = fstr['comp_trans']

    dist_all = Dist_stats(comp_all)
    dist_sust = Dist_stats(comp_sust)
    dist_trans = Dist_stats(comp_trans)

    trans_fstr3 = np.load("sus_trans_pie_3.npz")
    sust3 = trans_fstr3["sust"]
    trans3 = trans_fstr3["transient"]

    trans_fstr6 = np.load("sus_trans_pie_6.npz")
    sust6 = trans_fstr6["sust"]
    trans6 = trans_fstr6["transient"]

    sust = np.logical_or(sust3, sust6)
    trans = np.logical_or(trans3, trans6)

    (features_per_su, reg_list, avails) = get_dataset(denovo=True, random_type=None)
    trans_avail = trans[avails]
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
    sust_avail = sust[avails]
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]

    coeff_all = dist_all.get_coeff(features_per_su)
    coeff_sust = dist_sust.get_coeff(fr_sust)
    coeff_trans = dist_trans.get_coeff(fr_trans)

    for r in range(repeats):
        (features_per_su, reg_list, _avails) = get_dataset(denovo=True, random_type='one_trial')
        fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
        fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]

        dist_all.append_data(features_per_su, coeff_all)
        dist_sust.append_data(fr_sust, coeff_sust)
        dist_trans.append_data(fr_trans, coeff_trans)

    dist_all_list = dist_all.get_data()
    dist_sust_list = dist_sust.get_data()
    dist_trans_list = dist_trans.get_data()

    np.savez(f'pca_dist_{repeats}.npz', dist_all=dist_all_list, dist_sust=dist_sust_list, dist_trans=dist_trans_list)

    return (dist_all_list, dist_sust_list, dist_trans_list)