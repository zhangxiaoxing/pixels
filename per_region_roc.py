# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:27:55 2020

@author: Libra

unfinished as of 2/21/2020 due to change in data used in the study
"""

import os

# import sys
import csv
import h5py
import numpy as np
import scipy.stats as stats
# from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
# from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import selectivity as zpy
import sklearn.metrics as metrics


class auc_stats:
    def __init__(self):
        self.auc_SP = []
        self.auc_ED = []
        self.auc_LD = []

        self.selectivity_SP = []
        self.selectivity_ED = []
        self.selectivity_LD = []

    def selectivity(self, A, B):
        mA = np.mean(A)
        mB = np.mean(B)
        sums = np.mean(A) + np.mean(B)
        if sums > 0:
            return (mA - mB) / sums
        else:
            return 0

    def auroc(self, A, B):
        if np.mean(A) > np.mean(B):
            (A, B) = (B, A)
        try:
            return metrics.roc_auc_score(
                np.concatenate((
                    np.zeros(A.shape[1]),
                    np.ones(B.shape[1])
                )),
                np.concatenate((
                    np.mean(A, axis=0),
                    np.mean(B, axis=0),
                ))
            )
        except ValueError:
            return 0.5
        pass

    def one_class(self, trial_FR, trials, delay, sample, welltrain_window, correct_resp):
        rtn = trial_FR[:,
              np.all(np.vstack((trials[:, 5] == delay, trials[:, 2] == sample, welltrain_window, correct_resp,)),
                     axis=0, ), :, ]
        return rtn

    # unprocessed data entry point
    def process_auc_stats(self, trial_firing_rate, trials, welltrain_window=None, correct_resp=None):
        ### TODO: when variables are empty
        # [bin:trial:SU]
        # breakpoint()
        # auc contrast: sample vs pair

        # FR_D3_S1 = self.one_class(trial_firing_rate, trials, 3, 4, welltrain_window, correct_resp)
        # FR_D3_S2 = self.one_class(trial_firing_rate, trials, 3, 8, welltrain_window, correct_resp)
        FR_D6_S1 = self.one_class(trial_firing_rate, trials, 6, 4, welltrain_window, correct_resp)
        FR_D6_S2 = self.one_class(trial_firing_rate, trials, 6, 8, welltrain_window, correct_resp)
        # (bin, trial, su), e.g.(68, 39, 274)

        BS = np.arange(0, 10)
        SP = np.arange(12, 16)
        ED = np.arange(18, 28)
        LD = np.arange(28, 38)

        for su_idx in range(trial_firing_rate.shape[2]):
            self.selectivity_SP.append(self.selectivity(FR_D6_S1[SP, :, su_idx], FR_D6_S2[SP, :, su_idx]))
            self.selectivity_ED.append(self.selectivity(FR_D6_S1[ED, :, su_idx], FR_D6_S2[ED, :, su_idx]))
            self.selectivity_LD.append(self.selectivity(FR_D6_S1[LD, :, su_idx], FR_D6_S2[LD, :, su_idx]))
            self.auc_SP.append(self.auroc(FR_D6_S1[SP, :, su_idx], FR_D6_S2[SP, :, su_idx]))
            self.auc_ED.append(self.auroc(FR_D6_S1[ED, :, su_idx], FR_D6_S2[ED, :, su_idx]))
            self.auc_LD.append(self.auroc(FR_D6_S1[LD, :, su_idx], FR_D6_S2[LD, :, su_idx]))

    def get_features(self):
        return np.vstack((self.auc_SP,
                          self.auc_ED,
                          self.auc_LD,
                          self.selectivity_SP,
                          self.selectivity_ED,
                          self.selectivity_LD))


### TODO: Data logic entry point
def get_dataset(denovo):
    if denovo:
        features_per_su = []
        reg_list = []
        dpath = None
        if os.path.exists("/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"):
            dpath = "/gpfsdata/home/zhangxiaoxing/pixels/DataSum/"
        else:
            dpath = r"D:\neupix\DataSum"
        for path in zpy.traverse(dpath):
            print(path)
            # SU_ids = []
            trial_firing_rate = None
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
                        trial_firing_rate = np.array(dset, dtype="double")
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

            (perf_desc, perf_code, welltrain_window, correct_resp) = zpy.judgePerformance(trials)

            if perf_code != 3:
                continue

            curr_stats = auc_stats()
            curr_stats.process_auc_stats(trial_firing_rate, trials, welltrain_window, correct_resp)
            features_per_su.append(curr_stats.get_features())
            reg_list.extend(suid_reg[1])

            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

            ### save to npz file
        features_arr = np.hstack(tuple(features_per_su))
        reg_arr = np.array(reg_list)
        np.savez_compressed(
            "auc.npz", features_arr=features_arr, reg_arr=reg_arr
        )
    ### load from npz file
    else:
        if not os.path.isfile("auc.npz"):
            print('missing data file!')
            return
        fstr = np.load("auc.npz", allow_pickle=True)
        features_arr = fstr["features_arr"]
        reg_arr = fstr["reg_arr"]

    return (features_arr, reg_arr)


def process_all(denovo=False):
    (features_arr, reg_arr) = get_dataset(denovo)
    reg_set = list(set(reg_arr.tolist()))
    features_arr = np.abs(features_arr)

    # plot_figure = False
    # if plot_figure:
    #     reg_n = []
    #     for i in range(len(reg_set)):
    #         reg_n.append(f'{np.count_nonzero(reg_arr == reg_set[i])} SU {reg_set[i]}')
    #     feature_tag = ['sample selective during sample',
    #                    'sample selective during early delay',
    #                    'sample selective during late delay',
    #                    'sample selective ONLY during late delay',
    #                    'sample selective during decision making',
    #                    'test selective during decision making',
    #                    'pair selective during decision making',
    #                    'pair selective ONLY during reward',
    #                    'nonselective modulation during late delay',
    #                    'nonselective modulation during decision making',
    #                    'unmodulated']
    #     import matplotlib.pyplot as plt
    #     fh = plt.figure(figsize=(11.34, 40), dpi=300)
    su_factors = []
    su_factors.append(reg_set)

    reg_count = []
    for reg in reg_set:
        reg_count.append(np.count_nonzero(reg_arr == reg))
    su_factors.append(reg_count)

    for feature in range(features_arr.shape[0]):
        # 6 features, global mean, 5% 10% 15% top, expect 24 factors output

        # global mean
        per_region_su_factor_global = np.full(len(reg_set), np.nan, dtype=np.float64)
        per_region_su_factor_1pct = np.full(len(reg_set), np.nan, dtype=np.float64)
        per_region_su_factor_5pct = np.full(len(reg_set), np.nan, dtype=np.float64)
        per_region_su_factor_10pct = np.full(len(reg_set), np.nan, dtype=np.float64)
        for i in range(len(reg_set)):
            reg_sel = features_arr[feature, reg_arr == reg_set[i]]
            per_region_su_factor_global[i] = np.mean(reg_sel)
            if reg_count[i] >= 100:
                per_region_su_factor_1pct[i] = np.mean(
                    np.partition(reg_sel, np.ceil(reg_count[i] * 0.99).astype(np.int))[
                    np.ceil(reg_count[i] * 0.99).astype(np.int):])
            if reg_count[i] >= 20:
                per_region_su_factor_5pct[i] = np.mean(
                    np.partition(reg_sel, np.ceil(reg_count[i] * 0.95).astype(np.int))[
                    np.ceil(reg_count[i] * 0.95).astype(np.int):])
            if reg_count[i] >= 10:
                per_region_su_factor_10pct[i] = np.mean(
                    np.partition(reg_sel, np.ceil(reg_count[i] * 0.9).astype(np.int))[
                    np.ceil(reg_count[i] * 0.9).astype(np.int):])

        su_factors.append(per_region_su_factor_global.tolist())
        su_factors.append(per_region_su_factor_1pct.tolist())
        su_factors.append(per_region_su_factor_5pct.tolist())
        su_factors.append(per_region_su_factor_10pct.tolist())

        # if plot_figure:
        #     reg_idx = np.argsort(per_region_su_factor)
        #     ax = plt.subplot(11, 1, feature + 1)
        #     plt.bar(np.arange(len(reg_set)), per_region_su_factor[reg_idx])
        #     ax.set_xticks(np.arange(len(reg_set)))
        #     ax.set_xticklabels(np.array(reg_n)[reg_idx], rotation=90, ha='center', va='top', fontsize=7.5)
        #     ax.set_xlim(-1, len(reg_set))
        #     ax.set_ylabel(feature_tag[feature])

    # if plot_figure:
    #     plt.tight_layout(rect=[0, 0, 1, 0.95])
    #     fh.savefig("coding_feature.png", dpi=300, bbox_inches="tight")
    #     plt.show()

    su_factors = [list(i) for i in zip(*su_factors)]
    ### export csv for matlab GLM
    with open("auc_selectivity.csv", "w", newline="") as cf:
        cwriter = csv.writer(cf, dialect="excel")
        for row in su_factors:
            cwriter.writerow(row)


if __name__ == "__main__":
    process_all(True)



    # with open("reg.csv", "w", newline="") as cf:
    #     cwriter = csv.writer(cf, dialect="excel")
    #     for row in reg_list:
    #         cwriter.writerow(row)