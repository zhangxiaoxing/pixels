# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 00:14:22 2020

@author: Libra

walk through all trials to identify potential unbalance
"""
import os
import h5py
import csv
import numpy as np
# from scipy.stats.stats import pearsonr
import scipy.stats as stats
import su_region_align as align


class PrevTrialStats:
    def __init__(self):
        # self.per_trial_base_avg = None
        self.stats = None
        # self.correct_base_avg = None
        # self.err_base_avg = None

    def process_baseline_stats(self, trials):
        self.stats = []
        for tidx in range(1, trials.shape[0]):
            self.stats.append(np.hstack((trials[tidx - 1, 2:], (trials[tidx, 2:]))))

    def getFeatures(self):
        return self.stats


### all brain region entry point
def prepare_baseline():
    curr_stats = PrevTrialStats()
    all_sess_list = []
    dpath = align.get_root_path()
    # counter = 0
    for path in align.traverse(dpath):
        # path = r'K:\neupix\DataSum\200104_75_learning4_g1\200104_75_learning4_g1_imec0_cleaned'
        print(path)
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
                    dset = ffr["Trials"]
                    trials = np.array(dset, dtype="double").T
                done_read = True
            except OSError:
                print("h5py read error handled")
        if trials is None:
            continue

        (perf_desc, perf_code, welltrain_window, correct_resp) = align.judgePerformance(trials)

        if perf_code != 3:
            continue

        curr_stats.process_baseline_stats(trials)

        onesession = curr_stats.getFeatures()
        all_sess_list.extend(onesession)
        # counter += 1
        # if counter > 2:
        #     break

    export = []
    export.extend(list(zip(*all_sess_list)))

    export = list(zip(*export))
    with open("baseline_stats.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in export:
            cwriter.writerow(row)

    return (all_sess_list)


def plot_features():
    # TODO: move figure plot here
    pass


# %% main
def process_all(denovo=False):
    ### save raw data file
    if denovo:
        (all_sess_list) = prepare_baseline()
        arr = np.array(all_sess_list)
        np.savez_compressed('prev_trial_stats.npz', all_sess_arr=arr)
    else:
        fstr = np.load('prev_trial_stats.npz')
        arr = fstr['all_sess_arr']

    t44 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 4] == 4))
    t48 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 4] == 8))
    t84 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 4] == 4))
    t88 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 4] == 8))

    s44 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 4] == 4))
    s48 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 4] == 8))
    s84 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 4] == 4))
    s88 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 4] == 8))

    correct = np.logical_xor(arr[:, 4] == arr[:, 5], arr[:, 6] > 0)

    samp_Cong = np.count_nonzero(np.logical_and(arr[:, 0] == arr[:, 4], correct)) / np.count_nonzero(
        arr[:, 0] == arr[:, 4])
    samp_non_Cong = np.count_nonzero(np.logical_and(arr[:, 0] != arr[:, 4], correct)) / np.count_nonzero(
        arr[:, 0] != arr[:, 4])

    test_Cong = np.count_nonzero(np.logical_and(arr[:, 1] == arr[:, 4], correct)) / np.count_nonzero(
        arr[:, 1] == arr[:, 4])
    test_non_Cong = np.count_nonzero(np.logical_and(arr[:, 1] != arr[:, 4], correct)) / np.count_nonzero(
        arr[:, 1] != arr[:, 4])

    print([s44, s48, s84, s88, t44, t48, t84, t88])
    print([samp_Cong, samp_non_Cong, test_Cong, test_non_Cong])


if __name__ == "__main__":
    process_all()
