# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra


@author: Libra

Generate transient sustained percentage pie chart
Depends on intermediate data from CQTransient.py-> CQ_transient.m

For 1D and 2D decoding, use sus_transient_decoding.py

the output of this script is used for optogenetic mapping as of 4.19

"""

import os
import h5py
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import scipy.stats as stats
import su_region_align as align
from matplotlib import rcParams
from sklearn.metrics import roc_auc_score
import re


class per_sec_stats:
    def __init__(self):
        self.per_sec_sel = None  # 7 x SU , sample + 6s delay
        self.per_sec_sel_raw = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS1 = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS2 = None  # 7 x SU , sample + 6s delay
        self.per_sec_auc = None  # 7 x SU , sample + 6s delay
        self.per_sec_fr = None
        self.per_sec_wrs_p = None  # 7 x SU , sample + 6s delay
        self.non_sel_mod = None  # 7 x SU
        self.non_mod = None  # 7 x SU
        self.baseline_sel = None  # 1 x SU
        self.early_in_6s = None
        self.late_in_6s = None

        self.use_ranksum = True

    def bool_stats_test(self, A, B, bonf=1):
        ### TODO alternative use of perm_test
        if self.use_ranksum:
            try:
                # flatten is used instead of mean to control for variance
                # ranksum is supposedly insensitive to length of data
                (stat, p) = stats.mannwhitneyu(
                    A.flatten(),
                    B.flatten(),
                    alternative="two-sided",
                )
            except ValueError:
                p = 1

            try:
                auc = roc_auc_score(
                    np.concatenate((
                        np.zeros(A.size),
                        np.ones(B.size)
                    )),
                    np.concatenate((
                        A.flatten(),
                        B.flatten(),
                    ))
                )
            except ValueError:
                auc = 0.5

            if np.sum(A)+np.sum(B)==0:
                selectivity=0
            else:
                mma = np.mean(A)
                mmb = np.mean(B)
                selectivity=(mma-mmb)/(mma+mmb)

            return (p < (0.05 / bonf), p, selectivity,auc)




    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc

    def process_select_stats(self, trial_FR, trials, welltrain_window=None, correct_resp=None, delay=6):

        # trial_FR [bin, trial, SU]
        ### TODO: when variables are none

        firing_rate_6 = trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        firing_rate_3 = trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :]
        trial_perf_sel_6 = trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        trial_perf_sel_3 = trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :]

        trial_sel_left_6 = trial_perf_sel_6[:, 2] == 4
        trial_sel_left_3 = trial_perf_sel_3[:, 2] == 4

        self.per_sec_sel = np.zeros((delay+1, trial_FR.shape[2]))
        self.per_sec_sel_raw = np.zeros((delay+1, trial_FR.shape[2]))
        self.per_sec_auc = np.zeros((delay+1, trial_FR.shape[2]))
        self.per_sec_fr = np.zeros((delay+1, trial_FR.shape[2]))
        self.per_sec_wrs_p = np.zeros((delay + 1, trial_FR.shape[2]))

        self.per_sec_prefS1 = np.zeros_like(self.per_sec_sel)
        self.per_sec_prefS2 = np.zeros_like(self.per_sec_sel)
        self.non_sel_mod = np.zeros_like(self.per_sec_sel)
        self.baseline_sel = np.zeros((4, trial_FR.shape[2]))

        for su_idx in range(trial_FR.shape[2]):
            onesu_6 = np.squeeze(firing_rate_6[:, :, su_idx]).T
            onesu_3 = np.squeeze(firing_rate_3[:, :, su_idx]).T
            left_trials_6 = onesu_6[trial_sel_left_6, :]
            right_trials_6 = onesu_6[~trial_sel_left_6, :]

            left_trials_3 = onesu_3[trial_sel_left_3, :]
            right_trials_3 = onesu_3[~trial_sel_left_3, :]

            # SP = np.arange(12, 16) ED = np.arange(18, 28) LD = np.arange(28, 38)

            baseL = np.mean(trial_FR[:10, :, su_idx][:, (trials[:, 2] == 4) & welltrain_window & correct_resp],
                            axis=0)
            baseR = np.mean(trial_FR[:10, :, su_idx][:, (trials[:, 2] == 8) & welltrain_window & correct_resp],
                            axis=0)
            self.baseline_sel[:,su_idx] = self.bool_stats_test(baseL, baseR)

            if delay == 6:
                for bin_idx in range(0, 7):
                    bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)
                    spksum = (np.sum(left_trials_6[:, bins]) + np.sum(right_trials_6[:, bins]))

                    if spksum == 0:
                        self.per_sec_sel_raw[bin_idx, su_idx] = 0
                        self.per_sec_auc[bin_idx, su_idx] = 0.5
                        self.per_sec_fr[bin_idx, su_idx] = 0
                        self.per_sec_sel[bin_idx, su_idx] = False
                        self.per_sec_wrs_p[bin_idx, su_idx] = 1
                    else:
                        binstat=self.bool_stats_test(left_trials_6[:, bins],
                                             right_trials_6[:, bins], 7)
                        self.per_sec_sel_raw[bin_idx, su_idx] = binstat[2]
                        self.per_sec_sel[bin_idx, su_idx] = binstat[0]
                        self.per_sec_auc[bin_idx, su_idx] =binstat[3]
                        self.per_sec_fr[bin_idx, su_idx] =np.mean(np.concatenate((left_trials_6[:,bins],right_trials_6[:,bins])))
                        self.per_sec_wrs_p[bin_idx, su_idx] = binstat[1]

                    if self.per_sec_sel[bin_idx, su_idx]:
                        if np.mean(left_trials_6[:, bins]) > np.mean(right_trials_6[:, bins]):
                            self.per_sec_prefS1[bin_idx, su_idx] = 1
                        else:
                            self.per_sec_prefS2[bin_idx, su_idx] = 1

                    self.non_sel_mod[bin_idx, su_idx] = ((not self.per_sec_sel[bin_idx, su_idx]) and
                                                         self.bool_stats_test(onesu_6[:, bins],
                                                                              onesu_6[:, 6:10], 7)[0])
            elif delay == 3: # TODO update new stats
                for bin_idx in range(0, 4):
                    bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)

                    spksum = (np.sum(left_trials_6[:, bins]) - np.sum(right_trials_6[:, bins]))
                    if spksum == 0:
                        self.per_sec_sel_raw[bin_idx, su_idx] = 0
                        self.per_sec_sel[bin_idx, su_idx] = False
                    else:
                        self.per_sec_sel_raw[bin_idx, su_idx] = spksum / (
                                np.sum(left_trials_6[:, bins]) + np.sum(right_trials_6[:, bins]))
                        self.per_sec_sel[bin_idx, su_idx] = self.bool_stats_test(left_trials_3[:, bins],
                                                                                 right_trials_3[:, bins], 4)[0]

                    if self.per_sec_sel[bin_idx, su_idx]:
                        if np.mean(left_trials_3[:, bins]) > np.mean(right_trials_3[:, bins]):
                            self.per_sec_prefS1[bin_idx, su_idx] = 1
                        else:
                            self.per_sec_prefS2[bin_idx, su_idx] = 1

                    self.non_sel_mod[bin_idx, su_idx] = ((not self.per_sec_sel[bin_idx, su_idx]) and
                                                         self.bool_stats_test(onesu_3[:, bins],
                                                                              onesu_3[:, 6:10], 4))[0]
                # self.non_mod[bin_idx, su_idx] = not (
                #         self.per_sec_sel[bin_idx, su_idx] or self.non_sel_mod[bin_idx, su_idx])

    def getFeatures(self):
        return (self.per_sec_sel, self.non_sel_mod, self.per_sec_prefS1, self.per_sec_prefS2, self.baseline_sel,
                self.per_sec_sel_raw,self.per_sec_auc,self.per_sec_wrs_p,self.per_sec_fr)


### all brain region entry point

def prepare_data(delay=6):
    curr_stats = per_sec_stats()
    per_sec_sel_list = []
    non_sel_mod_list = []
    bs_sel_list = []
    perfS1_list = []
    perfS2_list = []
    raw_sel_list = []
    auc_list = []
    fr_list=[]
    wrs_p_list=[]
    # non_mod_list = []
    folder_list = []
    cluster_id_list = []
    reg_list = []
    dpath = align.get_root_path()
    for path in sorted(align.traverse(dpath)):
        print(path)
        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            continue
        done_read = False
        while not done_read:
            try:
                with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                    # print(list(ffr.keys()))
                    if not "SU_id" in ffr.keys():
                        done_read = True
                        print("missing su_id key in path ", path)
                        continue
                    dset = ffr["SU_id"]
                    SU_ids = np.array(dset, dtype="uint16")
                    dset = ffr["FR_All"]
                    trial_FR = np.array(dset, dtype="double")
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.process_select_stats(trial_FR, trials, welltrain_window, correct_resp, delay=delay)

        (per_sec_sel, non_sel_mod, perfS1, perfS2, bs_sel, raw_sel,auc,wrsp,fr) = curr_stats.getFeatures()
        per_sec_sel_list.append(per_sec_sel)
        non_sel_mod_list.append(non_sel_mod)
        bs_sel_list.append(bs_sel)
        # non_mod_list.append(non_sel_mod)
        perfS1_list.append(perfS1)
        perfS2_list.append(perfS2)
        raw_sel_list.append(raw_sel)
        auc_list.append(auc)
        fr_list.append(fr)
        wrs_p_list.append(wrsp)
        cluster_id_list.append(SU_ids)
        folder_list.append([re.search(r"(?<=DataSum\\)(.*)", path)[1]] * SU_ids.shape[1])

        for i in range(SU_ids.shape[1]):
            if suid_reg[0][i] != str(SU_ids[0, i]):
                print(f"unmatch cluster id and id_reg, {path}, {SU_ids[0, i]}")
                input("press Enter to continue")
        reg_list.append(suid_reg[1])
    return (
        per_sec_sel_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list, bs_sel_list, raw_sel_list,
        cluster_id_list, folder_list,auc_list,wrs_p_list,fr_list)


def prepare_data_sync():
    dpath = align.get_root_path()
    syncdata = []
    for path in align.traverse(dpath):
        print(path)
        SU_ids = []
        trials = None
        if not os.path.isfile(os.path.join(path, "su_id2reg.csv")):
            continue
        done_read = False
        while not done_read:
            try:
                with h5py.File(os.path.join(path, "FR_All.hdf5"), "r") as ffr:
                    # print(list(ffr.keys()))
                    if not "SU_id" in ffr.keys():
                        done_read = True
                        print("missing su_id key in path ", path)
                        continue
                    dset = ffr["SU_id"]
                    SU_ids = np.array(dset, dtype="uint16")
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue
        suid_reg = []
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue
        ### TODO: export per su per file alignment file

        for one_su in SU_ids.flatten():
            # breakpoint()
            # print(path+','+one_su.astype(np.str))
            syncdata.append([path, one_su])

    with open("su_list.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in syncdata:
            cwriter.writerow(row)

    return syncdata


def plot_features():
    # TODO: move figure plot here
    pass


### entry point
def process_all(denovo=False, toPlot=False, toExport=False, delay=6, counterclock=False):
    # output from CQ's algorithm of transient coding, both 3s and 6s obtained
    with h5py.File('CQ_transient.hdf5', 'r') as fr:
        if delay == 6:
            CQ_transient = np.array(fr["transient6"]).T
        elif delay == 3:
            CQ_transient = np.array(fr["transient3"]).T
        elif delay == 'early3in6':
            CQ_transient = np.array(fr['transient_early3in6']).T
        elif delay == 'late3in6':
            CQ_transient = np.array(fr['transient_late3in6']).T
        else:
            print("error delay time")
            sys.exit(-1)

    # per_sec_sel_arr = None
    # non_sel_mod_arr = None
    # perfS1_arr = None
    # perfS2_arr = None
    # reg_arr = None
    delay_num = delay
    if delay == 'early3in6' or delay == 'late3in6':
        delay_num = 6
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
    elif delay == 3 or delay == 'early3in6':
        delay_bins = np.arange(1, 4)
        early_bins = np.arange(1, 3)
        late_bins = np.arange(2, 4)
    elif delay == 'late3in6':
        delay_bins = np.arange(4, 7)
        early_bins = np.arange(4, 6)
        late_bins = np.arange(5, 7)

    bs_count = np.count_nonzero(bs_sel[0,:])
    non_bs = np.logical_not(bs_sel)

    any_sel = np.logical_and(non_bs, np.any(per_sec_sel_arr, axis=0))
    any_sel_count = np.count_nonzero(any_sel)

    sample_only = np.logical_and(any_sel, np.logical_not(np.any(per_sec_sel_arr[delay_bins, :], axis=0)))

    sample_only_count = np.count_nonzero(sample_only)

    delay_sel = np.logical_and(any_sel, np.logical_not(sample_only))

    delay_sel_count = np.count_nonzero(delay_sel)

    non_sel = np.logical_and(non_bs, np.logical_not(any_sel))
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
        rcParams['pdf.fonttype'] = 42
        rcParams['ps.fonttype'] = 42
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Arial']

        frac = [switched_count + unclassified_count + bs_count + sample_only_count, sust_count, transient_count,
                non_sel_mod_count + non_mod_count, ]
        print(np.sum(frac))
        explode = (0, 0.1, 0.1, 0)
        labels = ('unclassified', 'sustained', 'transient',
                  'non-selective')

        (fh, ax) = plt.subplots(1, 1, figsize=(12 / 2.54, 4 / 2.54), dpi=300)
        if counterclock:
            startangle = -60
        else:
            startangle = 240
        ax.pie(frac, explode=explode, labels=labels, autopct='%1.1f%%', shadow=False, startangle=startangle,
               counterclock=counterclock, colors=('grey', 'blue', 'red', 'black'))
        # ax.axis('equal')
        ax.set_xlim((-0.7, 0.7))
        fh.suptitle(f'{delay}s delay')
        plt.show()
        fh.savefig(f'sus_trans_pie_{delay}.pdf')

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


def bars():
    # load and align 6s 3s delay selectivity data

    # plot transient bars
    # 3s transient from
    norm = Normalize(vmin=0, vmax=0.6)
    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv(3, row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay trials')
    ax.set_ylabel('coding in 3s delay trials')
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_3vs6.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")

    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv('early3in6', row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay early half')
    ax.set_ylabel('coding in 3s delay trials')
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_3vse3.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")

    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv('early_late', row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay early half')
    ax.set_ylabel('coding in 6s delay late half')

    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_early_vs_late.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")

    # (fig, ax) = plt.subplots(7, 5, sharex=True, sharey=True, figsize=(12, 12), dpi=200)
    # # 3s equiv in 6s
    # for row in range(7):
    #     temp = subgroup_equiv(3, row)
    #     ax[row, 0].bar(range(len(temp)), temp)
    #
    # for row in range(7):
    #     temp = subgroup_equiv(6, row)
    #     ax[row, 1].bar(range(len(temp)), temp)
    #
    # for row in range(7):
    #     temp = subgroup_equiv('early3in6', row)
    #     ax[row, 2].bar(range(len(temp)), temp)
    #
    # for row in range(7):
    #     temp = subgroup_equiv('late3in6', row)
    #     ax[row, 3].bar(range(len(temp)), temp)
    #
    # for row in range(7):
    #     temp = subgroup_equiv('early_late', row)
    #     ax[row, 4].bar(range(len(temp)), temp)
    #
    # [ax[6, x].set_xticks(range(7)) for x in [0, 1, 2, 3, 4]]
    # [ax[6, x].set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
    #                            'nonselective_modulation', 'non-modulated'), rotation=45, va='top', ha='right') for x in
    #  [0, 1, 2, 3, 4]]
    # ax[6, 0].set_xlabel('3s sel. SU in 6s trials')
    # ax[6, 1].set_xlabel('6s sel. SU in 3s trials')
    # ax[6, 2].set_xlabel('6s early sel. in 3s trials')
    # ax[6, 3].set_xlabel('6s late sel. in 3s trials')
    # ax[6, 4].set_xlabel('6s late sel. in early sel.')
    #
    # ax[0, 0].set_ylabel('sustained')
    # ax[1, 0].set_ylabel('transient')
    # ax[2, 0].set_ylabel('switched')
    # ax[3, 0].set_ylabel('unclassified')
    # ax[4, 0].set_ylabel('only during sample')
    # ax[5, 0].set_ylabel('non-selective mod.')
    # ax[6, 0].set_ylabel('non-modulated')
    #
    # fig.savefig('redistribt_between_3s_and_6s.png', bbox_inches='tight')
    # plt.show()
    # 6s transient from
    # 6s transient to


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


if __name__ == "__main__":
    # prepare_data_sync()
    # delay can be 'early3in6','late3in6','3','6'
    process_all(denovo=True, toPlot=False, toExport=True, delay=6, counterclock=False)
    # process_all(denovo=False, toPlot=False, toExport=True, delay=3, counterclock=False)
    # process_all(denovo=False, toPlot=True, toExport=False, delay='early3in6')
    # process_all(denovo=False, toPlot=True, toExport=False, delay='late3in6')

    # bars()
