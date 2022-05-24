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
from pixelStats import baselineVector
import per_second_stats


class per_sec_stats:
    def __init__(self):
        self.per_sec_sel = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS1 = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS2 = None  # 7 x SU , sample + 6s delay
        self.non_sel_mod = None  # 7 x SU
        self.peak = None  # 7 x SU
        self.baseline_sel = None  # 1 x SU

    def bool_stats_test(self, A, B, bonf=1):
        try:
            # flatten is used instead of mean to control for variance
            # ranksum is supposedly insensitive to length of data
            (stat, p) = stats.mannwhitneyu(
                A.flatten(),
                B.flatten(),
                alternative="two-sided",
            )
            return p < (0.05 / bonf)

        except ValueError:
            return False

    def process_select_stats(self, trial_FR, trials, welltrain_window=None, correct_resp=None, delay=6):

        # trial_FR [bin, trial, SU]
        ### TODO: when variables are none

        firing_rate_6 = trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        firing_rate_3 = trial_FR[:, (trials[:, 5] == 3) & welltrain_window & correct_resp, :]
        trial_perf_sel_6 = trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        trial_perf_sel_3 = trials[(trials[:, 5] == 3) & welltrain_window & correct_resp, :]

        trial_sel_left_6 = trial_perf_sel_6[:, 2] == 4
        trial_sel_left_3 = trial_perf_sel_3[:, 2] == 4
        if delay == 6:
            self.per_sec_sel = np.zeros((7, trial_FR.shape[2]))
        else:
            self.per_sec_sel = np.zeros((4, trial_FR.shape[2]))

        self.per_sec_prefS1 = np.zeros_like(self.per_sec_sel)
        self.per_sec_prefS2 = np.zeros_like(self.per_sec_sel)
        self.non_sel_mod = np.zeros_like(self.per_sec_sel)
        self.baseline_sel = np.zeros(trial_FR.shape[2])
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
            self.baseline_sel[su_idx] = self.bool_stats_test(baseL, baseR)

            if delay == 6:
                for bin_idx in range(0, 7):
                    bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)
                    self.per_sec_sel[bin_idx, su_idx] = self.bool_stats_test(left_trials_6[:, bins],
                                                                             right_trials_6[:, bins], 7)
                    if self.per_sec_sel[bin_idx, su_idx]:
                        if np.mean(left_trials_6[:, bins]) > np.mean(right_trials_6[:, bins]):
                            self.per_sec_prefS1[bin_idx, su_idx] = 1
                        else:
                            self.per_sec_prefS2[bin_idx, su_idx] = 1

                    self.non_sel_mod[bin_idx, su_idx] = ((not self.per_sec_sel[bin_idx, su_idx]) and
                                                         self.bool_stats_test(onesu_6[:, bins],
                                                                              onesu_6[:, 6:10], 7))
            elif delay == 3:
                for bin_idx in range(0, 4):
                    bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)
                    self.per_sec_sel[bin_idx, su_idx] = self.bool_stats_test(left_trials_3[:, bins],
                                                                             right_trials_3[:, bins], 4)
                    if self.per_sec_sel[bin_idx, su_idx]:
                        if np.mean(left_trials_3[:, bins]) > np.mean(right_trials_3[:, bins]):
                            self.per_sec_prefS1[bin_idx, su_idx] = 1
                        else:
                            self.per_sec_prefS2[bin_idx, su_idx] = 1

                    self.non_sel_mod[bin_idx, su_idx] = ((not self.per_sec_sel[bin_idx, su_idx]) and
                                                         self.bool_stats_test(onesu_3[:, bins],
                                                                              onesu_3[:, 6:10], 4))
                # self.non_mod[bin_idx, su_idx] = not (
                #         self.per_sec_sel[bin_idx, su_idx] or self.non_sel_mod[bin_idx, su_idx])

    def getFeatures(self):
        return (
            self.per_sec_sel, self.non_sel_mod, self.per_sec_prefS1, self.per_sec_prefS2, self.baseline_sel, self.peak)


### all brain region entry point

def prepare_data(delay=6, reg_idx=1):
    curr_stats = per_sec_stats()
    per_sec_sel_list = []
    non_sel_mod_list = []
    bs_sel_list = []
    perfS1_list = []
    perfS2_list = []
    # non_mod_list = []
    reg_list = []
    dpath = align.get_root_path()
    for path in sorted(align.traverse(dpath)):
        print(path)
        # SU_ids = []
        trial_FR = None
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
                    # dset = ffr["SU_id"]
                    # SU_ids = np.array(dset, dtype="uint16")
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

        (per_sec_sel, non_sel_mod, perfS1, perfS2, bs_sel) = curr_stats.getFeatures()
        per_sec_sel_list.append(per_sec_sel)
        non_sel_mod_list.append(non_sel_mod)
        bs_sel_list.append(bs_sel)
        # non_mod_list.append(non_sel_mod)
        perfS1_list.append(perfS1)
        perfS2_list.append(perfS2)
        reg_list.extend(suid_reg[reg_idx])
    return (per_sec_sel_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list, bs_sel_list)


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
def process_all(denovo=False, toPlot=False, toExport=False, delay=6, reg_idx=1, counterclock=False):
    delay_num = delay
    if denovo:
        (per_sec_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list, bs_sel_list) = prepare_data(delay=delay,
                                                                                                         reg_idx=reg_idx)
        ### save raw data file
        reg_arr = np.array(reg_list)
        per_sec_sel_arr = np.hstack(per_sec_list)
        non_sel_mod_arr = np.hstack(non_sel_mod_list)
        perfS1_arr = np.hstack(perfS1_list)
        perfS2_arr = np.hstack(perfS2_list)
        bs_sel = np.hstack(bs_sel_list)
        np.savez_compressed(f'per_sec_sel_{delay_num}.npz', per_sec_sel_arr=per_sec_sel_arr,
                            non_sel_mod_arr=non_sel_mod_arr,
                            # non_mod_arr=non_mod_arr,
                            perfS1_arr=perfS1_arr,
                            perfS2_arr=perfS2_arr,
                            reg_arr=reg_arr,
                            delay=delay,
                            bs_sel=bs_sel)
    else:
        ### load back saved raw data
        if not os.path.isfile(f"per_sec_sel_{delay_num}.npz"):
            print('missing data file!')
            return
        fstr = np.load(f'per_sec_sel_{delay_num}.npz')
        per_sec_sel_arr = fstr['per_sec_sel_arr']
        non_sel_mod_arr = fstr['non_sel_mod_arr']
        if fstr['delay'] != delay_num:
            print('Delay duration mismatch')
        perfS1_arr = fstr['perfS1_arr']
        perfS2_arr = fstr['perfS2_arr']
        reg_arr = fstr['reg_arr']
        bs_sel = fstr['bs_sel']

    if delay == 6:
        delay_bins = np.arange(1, 7)
        early_bins = np.arange(1, 4)
        late_bins = np.arange(4, 7)
    elif delay == 3:
        delay_bins = np.arange(1, 4)
        early_bins = np.arange(1, 3)
        late_bins = np.arange(2, 4)

    bs_count = np.count_nonzero(bs_sel)
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

    # output from CQ's algorithm of transient coding, both 3s and 6s obtained
    with h5py.File(os.path.join('transient', 'CQ_transient.hdf5'), 'r') as fr:
        if delay == 6:
            CQ_transient = np.array(fr["transient6"]).T
        elif delay == 3:
            CQ_transient = np.array(fr["transient3"]).T
        else:
            print("error delay time")
            sys.exit(-1)
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
    ### select bins

    transient_count = np.count_nonzero(transient)

    unclassified = np.logical_and(non_sust, np.logical_not(cqtrans))
    unclassified_count = np.count_nonzero(unclassified)

    ### export list

    export_arr = np.vstack((sust, transient, switched, unclassified))
    np.savez_compressed(f'sus_trans_pie_{delay}.npz', sust=sust, transient=transient,
                        switched=switched, unclassified=unclassified, sample_only=sample_only,
                        non_sel_mod=non_sel_mod, non_mod=non_mod, bs_sel=bs_sel, reg_arr=reg_arr)

    if toExport:
        np.savetxt(f'transient_{delay}.csv', export_arr, fmt='%d', delimiter=',',
                   header='Sustained,transient,switched,unclassified')

        np.savetxt(f'transient_{delay}_reg.csv', reg_arr, fmt='%s', delimiter=',')

        all_sess_arr = np.vstack(
            (sust, transient, sample_only))

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


def gauss_average(x):
    return np.convolve(x, [0.06136, 0.24477, 0.38774, 0.24477, 0.06136], "same")


def find_peak(to_plot=False, delay=6, cpu=30, repeats=1000, n_trial=(20, 25), n_neuron=50):
    fstr = np.load("ctd.npz", allow_pickle=True)
    features_per_su = fstr["features_per_su"].tolist()
    # reg_list = fstr["reg_list"].tolist()
    sus_trans_flag = per_second_stats.process_all(denovo=False,
                                                  delay=delay)  # 33172 x 4, sust,trans,switch,unclassified
    # sus_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[0, :])[0]]
    trans_feat = [features_per_su[i] for i in np.nonzero(sus_trans_flag[1, :])[0]]

    keys = ["S1_3", "S2_3"] if delay == 3 else ["S1_6", "S2_6"]
    avail_sel = [(x[keys[0]].shape[1] >= n_trial[1] and x[keys[1]].shape[1] >= n_trial[1]) for x in trans_feat]
    su_index = np.nonzero(avail_sel)[0]
    su_selected_features = [trans_feat[i] for i in su_index]
    if delay == 3:
        bin_range = np.arange(16, 28)
    else:
        bin_range = np.arange(16, 40)
    all_prefered = []
    all_non_prefered = []

    for one_su in su_selected_features:
        (basemm, basestd) = baselineVector(np.hstack(tuple([x for x in one_su.values() if x.ndim > 1])))

        S1Norm = gauss_average(np.mean((one_su[keys[0]] - basemm) / basestd, axis=1))
        S2Norm = gauss_average(np.mean((one_su[keys[1]] - basemm) / basestd, axis=1))

        if np.max(S1Norm[bin_range]) > np.max(S2Norm[bin_range]):
            peak_idx = np.argmax(S1Norm[bin_range]) + 16

            prefered = S1Norm[peak_idx - 12:peak_idx + 13]
            nonpref = S2Norm[peak_idx - 12:peak_idx + 13]
            all_prefered.append(prefered)
            all_non_prefered.append(nonpref)
        else:
            peak_idx = np.argmax(S2Norm[bin_range]) + 16

            prefered = S2Norm[peak_idx - 12:peak_idx + 13]
            nonpref = S1Norm[peak_idx - 12:peak_idx + 13]
            all_prefered.append(prefered)
            all_non_prefered.append(nonpref)

    # np.savez_compressed('latent_inhibitory.npz', all_prefered=all_prefered, all_non_prefered=all_non_prefered)
    (fig, ax) = plt.subplots(1, 1)
    h0 = ax.plot(np.mean(np.vstack(all_prefered), axis=0), '-r')
    h1 = ax.plot(np.mean(np.vstack(all_non_prefered), axis=0), '-b')
    ax.legend([h0[0], h1[0]], ['prefered', 'non-prefered'])
    ax.set_xticks([0, 12, 24])
    ax.set_xticklabels([-3, 0, 3])
    ax.set_xlabel('time(s)')
    ax.set_ylabel('firing rate z-score')
    ax.set_title('peak firing rate')
    fig.savefig('latent_inhibation_avg_FR.png')
    plt.show()

    sidx = np.argsort(np.vstack(all_prefered)[:, 12])
    (fig, ax) = plt.subplots(1, 2, figsize=(3, 6), dpi=600)
    ax[0].imshow(np.vstack(all_prefered)[sidx, :], cmap="jet", aspect="auto", origin='lower', vmin=-2, vmax=2)
    ax[1].imshow(np.vstack(all_non_prefered)[sidx, :], cmap="jet", aspect="auto", origin='lower', vmin=-2, vmax=2)

    ax[0].set_xticks([0, 12, 24])
    ax[0].set_xticklabels([-3, 0, 3])
    ax[0].set_xlabel('time(s)')

    ax[1].set_xticks([0, 12, 24])
    ax[1].set_xticklabels([-3, 0, 3])
    ax[1].set_xlabel('time(s)')
    fig.suptitle('peak firing rate')
    fig.savefig('latent_inhibation_peak_FR.png',bbox_inches='tight')
    plt.show()


    all_prefered = []
    all_non_prefered = []
    for one_su in su_selected_features:
        (basemm, basestd) = baselineVector(np.hstack(tuple([x for x in one_su.values() if x.ndim > 1])))

        S1Norm = gauss_average(np.mean((one_su[keys[0]] - basemm) / basestd, axis=1))
        S2Norm = gauss_average(np.mean((one_su[keys[1]] - basemm) / basestd, axis=1))

        selectivity = S1Norm - S2Norm

        if np.max(selectivity[bin_range]) > (-np.min(selectivity[bin_range])):
            peak_idx = np.argmax(selectivity[bin_range]) + 16

            prefered = S1Norm[peak_idx - 12:peak_idx + 13]
            nonpref = S2Norm[peak_idx - 12:peak_idx + 13]
            all_prefered.append(prefered)
            all_non_prefered.append(nonpref)
        else:
            peak_idx = np.argmin(selectivity[bin_range]) + 16

            prefered = S2Norm[peak_idx - 12:peak_idx + 13]
            nonpref = S1Norm[peak_idx - 12:peak_idx + 13]
            all_prefered.append(prefered)
            all_non_prefered.append(nonpref)

        # selectivity = np.mean(one_su[keys[0]], axis=1) - np.mean(one_su[keys[1]], axis=1)
        #
        # selectivity = gauss_average(selectivity)[bin_range]

        # if np.max(selectivity) > (-np.min(selectivity)):
        #     peak_idx = np.argmin(selectivity) + 16

        #     prefered = (one_su[keys[1]][peak_idx - 12: peak_idx + 13, :] - basemm) / basestd
        #     nonpref = (one_su[keys[0]][peak_idx - 12: peak_idx + 13, :] - basemm) / basestd
        #     all_prefered.append(np.mean(prefered, axis=1))
        #     all_non_prefered.append(np.mean(nonpref, axis=1))
        # else:
        #     peak_idx = np.argmax(selectivity) + 16
        #     prefered = (one_su[keys[0]][peak_idx - 12: peak_idx + 13, :] - basemm) / basestd
        #     nonpref = (one_su[keys[1]][peak_idx - 12: peak_idx + 13, :] - basemm) / basestd
        #     all_prefered.append(np.mean(prefered, axis=1))
        #     all_non_prefered.append(np.mean(nonpref, axis=1))
    # np.savez_compressed('latent_inhibitory.npz', all_prefered=all_prefered, all_non_prefered=all_non_prefered)
    (fig, ax) = plt.subplots(1, 1)
    h0 = ax.plot(np.mean(np.vstack(all_prefered), axis=0), '-r')
    h1 = ax.plot(np.mean(np.vstack(all_non_prefered), axis=0), '-b')
    ax.legend([h0[0], h1[0]], ['prefered', 'non-prefered'])
    ax.set_xticks([0, 12, 24])
    ax.set_xticklabels([-3, 0, 3])
    ax.set_xlabel('time(s)')
    ax.set_ylabel('firing rate z-score')
    ax.set_title('peak selectivity')
    fig.savefig('latent_inhibation_avg_FR_sel.png')
    plt.show()

    sidx = np.argsort(np.vstack(all_prefered)[:, 12])
    (fig, ax) = plt.subplots(1, 2, figsize=(3, 6), dpi=600)
    ax[0].imshow(np.vstack(all_prefered)[sidx, :], cmap="jet", aspect="auto", origin='lower', vmin=-2, vmax=2)
    ax[1].imshow(np.vstack(all_non_prefered)[sidx, :], cmap="jet", aspect="auto", origin='lower', vmin=-2, vmax=2)

    ax[0].set_xticks([0, 12, 24])
    ax[0].set_xticklabels([-3, 0, 3])
    ax[0].set_xlabel('time(s)')

    ax[1].set_xticks([0, 12, 24])
    ax[1].set_xticklabels([-3, 0, 3])
    ax[1].set_xlabel('time(s)')
    fig.suptitle('peak selectivity')
    fig.savefig('latent_inhibation_peak_FR_sel.png',bbox_inches='tight')
    plt.show()

    return


if __name__ == "__main__":
    find_peak()
