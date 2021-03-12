# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:33:18 2021

@author: Libra
"""



def subgroup_equiv(delay, typeIdx):
    fstr6 = np.load('sus_trans_pie_6.npz')
    fstr3 = np.load('sus_trans_pie_3.npz')
    fstr_e3 = np.load('sus_trans_pie_early3in6.npz')
    fstr_l3 = np.load('sus_trans_pie_late3in6.npz')

    sample_only_3 = fstr3['sample_only']
    sample_only_6 = fstr6['sample_only']
    sample_only_e3 = fstr_e3['sample_only']
    sample_only_l3 = fstr_l3['sample_only']

    sus_3 = fstr3['sust']
    sus_6 = fstr6['sust']
    sus_e3 = fstr_e3['sust']
    sus_l3 = fstr_l3['sust']

    transient_3 = fstr3['transient']
    transient_6 = fstr6['transient']
    transient_e3 = fstr_e3['transient']
    transient_l3 = fstr_l3['transient']

    unclassified_3 = fstr3['unclassified']
    unclassified_6 = fstr6['unclassified']
    unclassified_e3 = fstr_e3['unclassified']
    unclassified_l3 = fstr_l3['unclassified']

    switched_3 = fstr3['switched']
    switched_6 = fstr6['switched']
    switched_e3 = fstr_e3['switched']
    switched_l3 = fstr_l3['switched']

    non_sel_mod_3 = fstr3['non_sel_mod']
    non_sel_mod_6 = fstr6['non_sel_mod']
    non_sel_mod_e3 = fstr_e3['non_sel_mod']
    non_sel_mod_l3 = fstr_l3['non_sel_mod']

    non_mod_3 = fstr3['non_mod']
    non_mod_6 = fstr6['non_mod']
    non_mod_e3 = fstr_e3['non_mod']
    non_mod_l3 = fstr_l3['non_mod']

    lists_3 = [sus_3, transient_3, switched_3, unclassified_3, sample_only_3, non_sel_mod_3, non_mod_3]
    lists_6 = [sus_6, transient_6, switched_6, unclassified_6, sample_only_6, non_sel_mod_6, non_mod_6]
    lists_e3 = [sus_e3, transient_e3, switched_e3, unclassified_e3, sample_only_e3, non_sel_mod_e3, non_mod_e3]
    lists_l3 = [sus_l3, transient_l3, switched_l3, unclassified_l3, sample_only_l3, non_sel_mod_l3, non_mod_l3]

    subgrp_dist = None
    if delay == 3:
        subgrp_dist = [np.count_nonzero(np.logical_and(lists_3[typeIdx], one_list)) for one_list in lists_6]
    elif delay == 6:
        subgrp_dist = [np.count_nonzero(np.logical_and(lists_6[typeIdx], one_list)) for one_list in lists_3]
    elif delay == 'early3in6':
        subgrp_dist = [np.count_nonzero(np.logical_and(lists_3[typeIdx], one_list)) for one_list in lists_e3]
    elif delay == 'late3in6':
        subgrp_dist = [np.count_nonzero(np.logical_and(lists_3[typeIdx], one_list)) for one_list in lists_e3]
    elif delay == 'early_late':
        subgrp_dist = [np.count_nonzero(np.logical_and(lists_l3[typeIdx], one_list)) for one_list in lists_e3]
    return (np.array(subgrp_dist), sum(subgrp_dist))