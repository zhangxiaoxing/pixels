# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra
"""
import os
import h5py
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import su_region_align as align


class per_sec_stats:
    def __init__(self):
        self.per_sec_sel = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS1 = None  # 7 x SU , sample + 6s delay
        self.per_sec_prefS2 = None  # 7 x SU , sample + 6s delay
        self.non_sel_mod = None  # 7 x SU
        self.non_mod = None  # 7 x SU

        self.use_ranksum = True

    def bool_stats_test(self, A, B):
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
                return p < (0.05 / 7)

            except ValueError:
                return False

    def exact_mc_perm_test(self, xs, ys, nmc):
        n, k = len(xs), 0
        diff = np.abs(np.mean(xs) - np.mean(ys))
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
        return k / nmc

    def processGLMStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):

        # trial_FR [bin, trial, SU]
        ### TODO: when variables are none

        firing_rate_6 = trial_FR[:, (trials[:, 5] == 6) & welltrain_window & correct_resp, :]
        trial_perf_sel = trials[(trials[:, 5] == 6) & welltrain_window & correct_resp, :]

        trial_sel_left = trial_perf_sel[:, 2] == 4

        self.per_sec_sel = np.zeros((7, trial_FR.shape[2]))
        self.per_sec_prefS1 = np.zeros((7, trial_FR.shape[2]))
        self.per_sec_prefS2 = np.zeros((7, trial_FR.shape[2]))
        self.non_sel_mod = np.zeros_like(self.per_sec_sel)
        # self.non_mod = np.zeros_like(self.per_sec_sel)

        for su_idx in range(trial_FR.shape[2]):
            onesu = np.squeeze(firing_rate_6[:, :, su_idx]).T
            left_trials = onesu[trial_sel_left, :]
            right_trials = onesu[~trial_sel_left, :]

            # SP = np.arange(12, 16) ED = np.arange(18, 28) LD = np.arange(28, 38)
            for bin_idx in range(0, 7):
                bins = np.arange(bin_idx * 4 + 12, bin_idx * 4 + 16)
                self.per_sec_sel[bin_idx, su_idx] = self.bool_stats_test(left_trials[:, bins],
                                                                         right_trials[:, bins])
                if self.per_sec_sel[bin_idx, su_idx]:
                    if np.mean(left_trials[:, bins]) > np.mean(right_trials[:, bins]):
                        self.per_sec_prefS1[bin_idx, su_idx] = 1
                    else:
                        self.per_sec_prefS2[bin_idx, su_idx] = 1

                self.non_sel_mod[bin_idx, su_idx] = ((not self.per_sec_sel[bin_idx, su_idx]) and
                                                     self.bool_stats_test(onesu[:, bins],
                                                                          onesu[:, 6:10]))
                # self.non_mod[bin_idx, su_idx] = not (
                #         self.per_sec_sel[bin_idx, su_idx] or self.non_sel_mod[bin_idx, su_idx])

    def getFeatures(self):
        return (self.per_sec_sel, self.non_sel_mod, self.per_sec_prefS1, self.per_sec_prefS2)


### all brain region entry point

def prepare_data():
    curr_stats = per_sec_stats()
    per_sec_sel_list = []
    non_sel_mod_list = []
    perfS1_list = []
    perfS2_list = []
    # non_mod_list = []
    reg_list = []
    dpath = align.get_root_path()
    for path in align.traverse(dpath):
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
        suid_reg = []
        with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
            l = list(csv.reader(csvfile))[1:]
            suid_reg = [list(i) for i in zip(*l)]

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.processGLMStats(trial_FR, trials, welltrain_window, correct_resp)

        (per_sec_sel, non_sel_mod, perfS1, perfS2) = curr_stats.getFeatures()
        per_sec_sel_list.append(per_sec_sel)
        non_sel_mod_list.append(non_sel_mod)
        # non_mod_list.append(non_sel_mod)
        perfS1_list.append(perfS1)
        perfS2_list.append(perfS2)
        reg_list.extend(suid_reg[1])
    return (per_sec_sel_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list)


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


# %% main
def process_all(denovo=False, toPlot=False, toExport=False):
    per_sec_sel_arr = None
    non_sel_mod_arr = None
    perfS1_arr = None
    perfS2_arr = None
    reg_arr = None
    if denovo:
        ### save raw data file
        (per_sec_list, non_sel_mod_list, perfS1_list, perfS2_list, reg_list) = prepare_data()
        per_sec_sel_arr = np.hstack(per_sec_list)
        non_sel_mod_arr = np.hstack(non_sel_mod_list)
        perfS1_arr = np.hstack(perfS1_list)
        perfS2_arr = np.hstack(perfS2_list)

        np.savez_compressed('per_sec_sel.npz', per_sec_sel_arr=per_sec_sel_arr,
                            non_sel_mod_arr=non_sel_mod_arr,
                            # non_mod_arr=non_mod_arr,
                            perfS1_arr=perfS1_arr,
                            perfS2_arr=perfS2_arr,
                            reg_arr=reg_list)
    else:
        ### load back saved raw data
        if not os.path.isfile("per_sec_sel.npz"):
            print('missing data file!')
            return
        fstr = np.load('per_sec_sel.npz')
        per_sec_sel_arr = fstr['per_sec_sel_arr']
        non_sel_mod_arr = fstr['non_sel_mod_arr']
        perfS1_arr = fstr['perfS1_arr']
        perfS2_arr = fstr['perfS2_arr']
        reg_arr = fstr['reg_arr']

    sample_only = np.count_nonzero(per_sec_sel_arr[0, :].astype(np.bool) & ~np.any(per_sec_sel_arr[1:7, :], axis=0))
    delay_sel = np.count_nonzero(np.any(per_sec_sel_arr[1:7, :], axis=0))
    non_sel_mod = np.count_nonzero(np.any(non_sel_mod_arr[1:7, :], axis=0) & ~np.any(per_sec_sel_arr[1:7, :], axis=0))
    non_mod = per_sec_sel_arr.shape[1] - sample_only - delay_sel - non_sel_mod
    ### TODO: load transient permute result

    if not os.path.isfile("per_sec_sel.npz"):
        print('missing data file!')
        sys.exit(0)
    fstr = np.load('per_sec_sel.npz')
    per_sec_sel_arr = fstr['per_sec_sel_arr']
    non_sel_mod_arr = fstr['non_sel_mod_arr']
    perfS1_arr = fstr['perfS1_arr']
    perfS2_arr = fstr['perfS2_arr']
    # reg_arr = fstr['reg_arr']
    transient6 = None
    transient3 = None
    with h5py.File(os.path.join('transient', 'CQ_transient.hdf5'), 'r') as fr:
        transient6 = np.array(fr["transient6"]).T
        transient3 = np.array(fr["transient3"]).T

    frac = [delay_sel, sample_only, non_sel_mod, non_mod]
    explode = (0.1, 0.1, 0, 0)
    labels = ('selective during delay', 'selective only during sample', 'Non-selective modulation', 'Unmodulated')
    fh = None
    axes = None
    if toPlot:
        (fh, axes) = plt.subplots(1, 2, figsize=(10, 5), dpi=96)
        axes[0].pie(frac, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
        axes[0].axis('equal')

        explode = (0.1, 0, 0, 0)
        labels = ('sustained', 'transient', 'transient-switched', 'unclassified')

    switched = np.any(perfS1_arr[1:7, :], axis=0) & np.any(perfS2_arr[1:7, :], axis=0) & transient6.astype(
        np.bool).flatten()

    sust = np.count_nonzero(np.all(per_sec_sel_arr[1:7, :], axis=0) & ~switched)
    transient = np.count_nonzero(np.any(per_sec_sel_arr[1:7, :], axis=0) & (~switched)
                                 & (~np.all(per_sec_sel_arr[1:7, :], axis=0)) & transient6.astype(np.bool))
    switched_count = np.count_nonzero(switched)

    unclassified = np.count_nonzero(np.any(per_sec_sel_arr[1:7, :], axis=0) & (~switched)
                                    & (~np.all(per_sec_sel_arr[1:7, :], axis=0)) & ~transient6.astype(np.bool))
    if toPlot:
        axes[1].pie([sust, transient, switched_count, unclassified], explode=explode, labels=labels,
                    autopct=lambda p: '{:.1f}%'.format(p * delay_sel / switched.shape[0]),
                    radius=np.sqrt(delay_sel / per_sec_sel_arr.shape[1]), shadow=True)
        axes[1].axis('equal')
        axes[0].set_xlim((-1.25, 1.25))
        axes[1].set_xlim((-1.25, 1.25))
        fh.savefig('sus_trans_pie.png')
        plt.show()


    ### export list
    sust_list = np.all(per_sec_sel_arr[1:7, :], axis=0) & ~switched
    transient_list = (np.any(per_sec_sel_arr[1:7, :], axis=0) & (~switched)
                      & (~np.all(per_sec_sel_arr[1:7, :], axis=0)) & transient6.astype(np.bool).flatten())
    switched_list = switched
    unclassified_list = (np.any(per_sec_sel_arr[1:7, :], axis=0) & (~switched)
                         & (~np.all(per_sec_sel_arr[1:7, :], axis=0)) & ~transient6.astype(np.bool).flatten())

    export_arr = np.vstack((sust_list, transient_list, switched_list, unclassified_list)).T.astype(np.int8)
    if toExport:
        np.savetxt('transient.csv', export_arr, fmt='%d', delimiter=',',
                   header='Sustained,transient,switched,unclassified')
    return export_arr




if __name__ == "__main__":
    # prepare_data_sync()
    process_all(False,True,False)
