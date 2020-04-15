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
import scipy.stats as stats
import su_region_align as align
import matplotlib.pyplot as plt
import statsmodels.stats.proportion as mdl


def exact_mc_perm_test(xs, ys, nmc):
    n, k = len(xs), 0
    diff = np.abs(np.mean(xs) - np.mean(ys))
    zs = np.concatenate([xs, ys])
    for j in range(nmc):
        np.random.shuffle(zs)
        k += diff < np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
    return k / nmc


def ratio(x):
    return x[0] / x[1]


class PrevTrialStats:
    def __init__(self):
        # self.per_trial_base_avg = None
        self.stats = None
        # self.correct_base_avg = None
        # self.err_base_avg = None

    def process_baseline_stats(self, trials, welltrain_window):
        self.stats = []
        for tidx in range(1, trials.shape[0]):
            if welltrain_window[tidx]:
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

        curr_stats.process_baseline_stats(trials, welltrain_window)

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

    s44 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 4] == 4))
    s48 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 4] == 8))
    s84 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 4] == 4))
    s88 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 4] == 8))

    st44 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 5] == 4))
    st48 = np.count_nonzero(np.logical_and(arr[:, 0] == 4, arr[:, 5] == 8))
    st84 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 5] == 4))
    st88 = np.count_nonzero(np.logical_and(arr[:, 0] == 8, arr[:, 5] == 8))

    ts44 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 4] == 4))
    ts48 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 4] == 8))
    ts84 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 4] == 4))
    ts88 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 4] == 8))

    t44 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 5] == 4))
    t48 = np.count_nonzero(np.logical_and(arr[:, 1] == 4, arr[:, 5] == 8))
    t84 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 5] == 4))
    t88 = np.count_nonzero(np.logical_and(arr[:, 1] == 8, arr[:, 5] == 8))

    d33 = np.count_nonzero(np.logical_and(arr[:, 3] < 5, arr[:, 7] < 5))
    d36 = np.count_nonzero(np.logical_and(arr[:, 3] < 5, arr[:, 7] > 5))
    d63 = np.count_nonzero(np.logical_and(arr[:, 3] > 5, arr[:, 7] < 5))
    d66 = np.count_nonzero(np.logical_and(arr[:, 3] > 5, arr[:, 7] > 5))

    correct = np.logical_xor(arr[:, 4] == arr[:, 5], arr[:, 6] > 0)

    samp_Cong = (np.count_nonzero(np.logical_and(arr[:, 0] == arr[:, 4], correct)), np.count_nonzero(
        arr[:, 0] == arr[:, 4]))
    samp_non_Cong = (np.count_nonzero(np.logical_and(arr[:, 0] != arr[:, 4], correct)), np.count_nonzero(
        arr[:, 0] != arr[:, 4]))

    samp_test_pair = (np.count_nonzero(np.logical_and(arr[:, 0] != arr[:, 5], correct)), np.count_nonzero(
        arr[:, 0] != arr[:, 5]))
    samp_test_non_pair = (np.count_nonzero(np.logical_and(arr[:, 0] == arr[:, 5], correct)), np.count_nonzero(
        arr[:, 0] == arr[:, 5]))

    test_sample_pair = (np.count_nonzero(np.logical_and(arr[:, 1] != arr[:, 4], correct)), np.count_nonzero(
        arr[:, 1] != arr[:, 4]))
    test_sample_non_pair = (np.count_nonzero(np.logical_and(arr[:, 1] == arr[:, 4], correct)), np.count_nonzero(
        arr[:, 1] == arr[:, 4]))

    test_Cong = (np.count_nonzero(np.logical_and(arr[:, 1] == arr[:, 5], correct)), np.count_nonzero(
        arr[:, 1] == arr[:, 5]))
    test_non_Cong = (np.count_nonzero(np.logical_and(arr[:, 1] != arr[:, 5], correct)), np.count_nonzero(
        arr[:, 1] != arr[:, 5]))

    trial_stats = [(s44 + s88) / (s44 + s84 + s48 + s88), (s48 + s84) / (s48 + s88 + s44 + s84),
                   (t44 + t88) / (t44 + t84 + t48 + t88), (t48 + t84) / (t48 + t88 + t44 + t84),
                   (st48 + st84) / (st48 + st84 + st44 + st88), (st44 + st88) / (st48 + st84 + st44 + st88),
                   (ts48 + ts84) / (ts48 + ts84 + ts44 + ts88), (ts44 + ts88) / (ts48 + ts84 + ts44 + ts88),
                   (d33 + d66) / (d33 + d36 + d63 + d66), (d36 + d63) / (d33 + d36 + d63 + d66),
                   ]

    perf_by_prev_trial = [ratio(samp_Cong), ratio(samp_non_Cong), ratio(test_Cong), ratio(test_non_Cong),
                          ratio(samp_test_pair), ratio(samp_test_non_pair), ratio(test_sample_pair),
                          ratio(test_sample_non_pair)]

    print(trial_stats)
    print(perf_by_prev_trial)

    # ps = stats.binom_test(np.count_nonzero(np.logical_and(arr[:, 0] != arr[:, 4], correct)), np.count_nonzero(
    #     arr[:, 0] != arr[:, 4]), samp_Cong, alternative='two-sided')
    # pst = stats.binom_test(np.count_nonzero(np.logical_and(arr[:, 0] == arr[:, 5], correct)), np.count_nonzero(
    #     arr[:, 0] == arr[:, 5]), samp_test_pair, alternative='two-sided')
    # pts = stats.binom_test(np.count_nonzero(np.logical_and(arr[:, 1] == arr[:, 4], correct)), np.count_nonzero(
    #     arr[:, 1] == arr[:, 4]), test_sample_pair, alternative='two-sided')
    # pt = stats.binom_test(np.count_nonzero(np.logical_and(arr[:, 1] != arr[:, 5], correct)), np.count_nonzero(
    #     arr[:, 1] != arr[:, 5]), test_Cong, alternative='two-sided')

    ps = stats.binom_test(*samp_non_Cong, ratio(samp_Cong), alternative='two-sided')
    pst = stats.binom_test(*samp_test_non_pair, ratio(samp_test_pair), alternative='two-sided')
    pts = stats.binom_test(*test_sample_non_pair, ratio(test_sample_pair), alternative='two-sided')
    pt = stats.binom_test(*test_non_Cong, ratio(test_Cong), alternative='two-sided')

    # np.count_nonzero(np.logical_and(arr[:, 0] == arr[:, 4], correct))
    # exact_mc_perm_test()

    print([ps, pt, pst, pts])

    ci = (np.vstack((mdl.proportion_confint(*samp_Cong),
                     mdl.proportion_confint(*samp_non_Cong),
                     mdl.proportion_confint(*test_Cong),
                     mdl.proportion_confint(*test_non_Cong),
                     mdl.proportion_confint(*samp_test_pair),
                     mdl.proportion_confint(*samp_test_non_pair),
                     mdl.proportion_confint(*test_sample_pair),
                     mdl.proportion_confint(*test_sample_non_pair),
                     )) - np.expand_dims(np.array(perf_by_prev_trial), axis=1)).T

    (fh, ax) = plt.subplots(1, 1, figsize=(6, 2), dpi=300)
    ax.bar([1, 2, 4, 5, 7, 8, 10, 11, 13, 14], trial_stats)
    ax.set_xticks([1, 2, 4, 5, 7, 8, 10, 11, 13, 14])
    ax.set_xticklabels(['same sample', 'oppo sample',
                           'same test', 'oppo test',
                           "S'-T pair", "S'-T non-pair",
                           "S-T' pair", "S-T' non-pair",
                           'same delay','diff delay'
                           ], rotation=30, va='top',ha='right')
    ax.set_ylabel('proportion in trials')
    fh.savefig('trial_design_stats1.png',bbox_inches='tight')
    plt.show()


    (fh, ax) = plt.subplots(1, 1, figsize=(6, 2), dpi=300)
    ax.bar([1, 2, 4, 5, 7, 8, 10, 11], perf_by_prev_trial)
    ax.errorbar([1, 2, 4, 5, 7, 8, 10, 11], perf_by_prev_trial, ci[0, :], fmt='k,', capsize=2, barsabove=True,lw=0.5)
    ax.set_xticks([1, 2, 4, 5, 7, 8, 10, 11])
    ax.set_xticklabels(['same sample', 'oppo. sample',
                           'same test', 'oppo. test',
                           "S'-T pair", "S'-T non-pair",
                           "S-T' pair", "S-T' non-pair"], rotation=30, va='top',ha='right')
    ax.set_ylabel('correct rate')
    ax.set_ylim([0.6, 1])
    fh.savefig('trial_design_stats2.png',bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    process_all(False)
