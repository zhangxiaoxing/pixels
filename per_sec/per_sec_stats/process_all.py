# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:28:16 2021

@author: Libra
"""

import h5py
import numpy as np
from per_sec_stats.plotpie import plotpie
from per_sec_stats.prepare_data import prepare_data

### entry point
def process_all(denovo=False, toPlot=False, toExport=False, delay=6, counterclock=False):


    if denovo:
        (per_sec_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list, bs_sel_list, raw_sel_list, cluster_id,
         path_list,auc_list,wrs_p,fr) = prepare_data(delay=delay)
        ### save raw data file
        per_sec_sel_arr = np.hstack(per_sec_list)
        non_sel_mod_arr = np.hstack(non_sel_mod_list)
        perfS1_arr = np.hstack(perfS1_list)
        perfS2_arr = np.hstack(perfS2_list)
        bs_sel = np.hstack(bs_sel_list)
        raw_sel_arr = np.hstack(raw_sel_list)
        auc_arr=np.hstack(auc_list)
        fr_arr=np.hstack(fr)
        wrs_p_arr = np.hstack(wrs_p)
        reg_arr = np.hstack(reg_list)
        clusterid_arr = np.hstack(cluster_id)
        path_arr = np.hstack(path_list)
        np.savez_compressed(f'per_sec_sel_{delay_num}.npz', per_sec_sel_arr=per_sec_sel_arr,
                            non_sel_mod_arr=non_sel_mod_arr,
                            # non_mod_arr=non_mod_arr,
                            perfS1_arr=perfS1_arr,
                            perfS2_arr=perfS2_arr,
                            reg_arr=reg_arr,
                            delay=delay,
                            bs_sel=bs_sel,
                            raw_sel_arr=raw_sel_arr,
                            clusterid_arr=clusterid_arr,
                            path_arr=path_arr,
                            auc_arr=auc_arr,
                            wrs_p_arr=wrs_p_arr,
                            fr_arr=fr_arr)
    else:
        ### load back saved raw data
        if not os.path.isfile(f"per_sec_sel_{delay_num}.npz"):
            print('missing data file!')
            return
        fstr = np.load(f'per_sec_sel_{delay_num}.npz', allow_pickle=True)
        per_sec_sel_arr = fstr['per_sec_sel_arr']
        non_sel_mod_arr = fstr['non_sel_mod_arr']
        if fstr['delay'] != delay_num:
            print('Delay duration mismatch')
        perfS1_arr = fstr['perfS1_arr']
        perfS2_arr = fstr['perfS2_arr']
        reg_arr = fstr['reg_arr']
        bs_sel = fstr['bs_sel']
        raw_sel_arr = fstr['raw_sel_arr']
        clusterid_arr = fstr['clusterid_arr']
        path_arr = fstr['path_arr']
        auc_arr=fstr['auc_arr']
        wrs_p_arr=fstr['wrs_p_arr']
        fr_arr=fstr['fr_arr']

    # rpt_workaround = ('M23_20191109_g0', '191018-DPA-Learning5_28_g1', '191226_64_learning6_g0_imec1_cleaned')
    # rpt = np.zeros_like(path_arr)

    if delay == 6:
        delay_bins = np.arange(1, 7)
        early_bins = np.arange(1, 4)
        late_bins = np.arange(4, 7)
    elif delay == 3:
        delay_bins = np.arange(1, 4)
        early_bins = np.arange(1, 3)
        late_bins = np.arange(2, 4)

    bs_count = np.count_nonzero(bs_sel[0,:])
    non_bs = np.logical_not(bs_sel[0,:])

    # any_sel = np.logical_and(non_bs, np.any(per_sec_sel_arr, axis=0))
    any_sel = np.any(per_sec_sel_arr, axis=0)
    any_sel_count = np.count_nonzero(any_sel)

    sample_only = np.logical_and(any_sel, np.logical_not(np.any(per_sec_sel_arr[delay_bins, :], axis=0)))

    sample_only_count = np.count_nonzero(sample_only)

    delay_sel = np.logical_and(any_sel, np.logical_not(sample_only))

    delay_sel_count = np.count_nonzero(delay_sel)

    # non_sel = np.logical_and(non_bs, np.logical_not(any_sel))
    non_sel = np.logical_not(any_sel)
    non_sel_count = np.count_nonzero(non_sel)

    non_sel_mod = np.logical_and(non_sel, np.any(non_sel_mod_arr[delay_bins, :], axis=0))
    non_sel_mod_count = np.count_nonzero(non_sel_mod)

    non_mod = np.logical_and(non_sel, np.logical_not(non_sel_mod))
    non_mod_count = np.count_nonzero(non_mod)

    sust = np.logical_and(delay_sel, np.logical_and(
        np.logical_xor(
            np.any(perfS1_arr[delay_bins, :], axis=0),
            np.any(perfS2_arr[delay_bins, :], axis=0)
        ), np.all(per_sec_sel_arr[delay_bins, :], axis=0)))
    sust_count = np.count_nonzero(sust)

    non_sust = np.logical_and(delay_sel, np.logical_not(sust))
    non_sust_count = np.count_nonzero(non_sust)

    cqtrans = np.logical_and(non_sust, CQ_transient.flatten())
    cqtrans_count = np.count_nonzero(cqtrans)

    switched = np.logical_and(cqtrans, np.logical_and(
        np.any(perfS1_arr[delay_bins, :], axis=0),
        np.any(perfS2_arr[delay_bins, :], axis=0)
    ))

    switched_count = np.count_nonzero(switched)

    transient = np.logical_and(cqtrans, np.logical_not(switched))
    early_in_6s = np.logical_and(transient, np.any(per_sec_sel_arr[early_bins, :], axis=0))
    late_in_6s = np.logical_and(transient, np.any(per_sec_sel_arr[late_bins, :], axis=0))

    transient_count = np.count_nonzero(transient)

    unclassified = np.logical_and(non_sust, np.logical_not(cqtrans))
    unclassified_count = np.count_nonzero(unclassified)

    if toPlot:
        plotpie()

    ### export list
    prefer_s = perfS1_arr + perfS2_arr * 2
    export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
    np.savez_compressed(f'sus_trans_pie_{delay}.npz', sust=sust, transient=transient,
                        switched=switched, unclassified=unclassified, sample_only=sample_only,
                        non_sel_mod=non_sel_mod, non_mod=non_mod, bs_sel=bs_sel, reg_arr=reg_arr,
                        raw_sel_arr=raw_sel_arr, auc_arr=auc_arr,fr_arr=fr_arr)

    if toExport:
        # np.savetxt(f'transient_{delay}.csv', export_arr, fmt='%d', delimiter=',',
        #            header='Sustained,transient,switched,unclassified,early_in_6s,late_in_6s,preferred)

        with h5py.File(f'transient_{delay}.hdf5', "w") as fw:
            fw.create_dataset('sus_trans', data=export_arr.astype('int8'))
            fw.create_dataset('raw_selectivity', data=raw_sel_arr.astype('float64'))
            fw.create_dataset('reg', data=reg_arr.astype('S10'))
            fw.create_dataset('path', data=path_arr.astype('S200'))
            fw.create_dataset('cluster_id', data=clusterid_arr.astype('uint16'))
            fw.create_dataset('auc', data=auc_arr.astype('float64'))
            fw.create_dataset('wrs_p',data=wrs_p_arr.astype('float64'))
            fw.create_dataset('fr',data=fr_arr.astype('float64'))
            fw.create_dataset('bs_sel',data=bs_sel.astype('float64'))


        # np.savetxt(f'transient_{delay}_reg.csv', reg_arr, fmt='%s', delimiter=',')

        all_sess_arr = np.vstack(
            (sust, transient, sample_only, non_sel_mod, non_mod, early_in_6s, late_in_6s))

        reg_set = list(set(reg_arr.tolist()))
        su_factors = []
        su_factors.append(reg_set)
        reg_count = []
        for reg in reg_set:
            reg_count.append(np.count_nonzero(reg_arr == reg))
        su_factors.append(reg_count)

        for feature in range(all_sess_arr.shape[0]):
            per_region_su_factor = np.zeros_like(reg_set, dtype=np.float64)
            for i in range(len(reg_set)):
                per_region_su_factor[i] = np.mean(all_sess_arr[feature, reg_arr == reg_set[i]])
            su_factors.append(per_region_su_factor.tolist())

        su_factors = [list(i) for i in zip(*su_factors)]

        ### export csv for matlab GLM
        with open("glm_coding_features_per_second.csv", "w", newline="") as cf:
            cwriter = csv.writer(cf, dialect="excel")
            for row in su_factors:
                cwriter.writerow(row)

    return export_arr