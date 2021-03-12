# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:33:49 2021

@author: Libra
"""

def quickStats(delay=6):
    trans_fstr = np.load(f'sus_trans_pie_{delay}.npz')
    # list(trans6.keys())
    sust = trans_fstr['sust']
    trans = trans_fstr['transient']
    reg_arr = trans_fstr['reg_arr']

    reg_set = tuple(set(reg_arr.tolist()))

    count = []
    for one_reg in reg_set:
        sust_count = np.count_nonzero(np.logical_and(reg_arr == one_reg, sust))
        trans_count = np.count_nonzero(np.logical_and(reg_arr == one_reg, trans))
        count.append([one_reg, sust_count, trans_count])