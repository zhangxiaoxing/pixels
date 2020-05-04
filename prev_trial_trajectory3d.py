# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os

import sys
import csv
import h5py
import numpy as np

# import scipy.stats as stats
# from sklearn.svm import SVC
from sklearn.model_selection import KFold

# from mcc.mcc import MaximumCorrelationClassifier
# from sklearn.model_selection import cross_val_score
# from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import su_region_align as align

# import per_second_stats
# from multiprocessing import Pool
# from matplotlib import rcParams
from sklearn.decomposition import PCA
from mpl_toolkits import mplot3d
from pixelStats import baselineVector


class ctd_stats:
    def __init__(self):
        self.su_list = []

    def splitMean(self, FR, suidx, mm, std):
        sel = list(self.kf.split(np.arange(FR.shape[1])))[0]
        return ()
        pass

        ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None):
        shifted_trial = np.vstack((np.zeros(6), trials[:-1, :]))

        ### TODO: when variables are empty

        FR_D3_S1_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S2_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S1_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D3_S2_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 3, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        FR_D6_S1_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D6_S2_S1 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 4, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        FR_D6_S1_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 4, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]
        FR_D6_S2_S2 = trial_FR[:, np.all(np.vstack(
            (trials[:, 5] == 6, trials[:, 2] == 8, shifted_trial[:, 2] == 8, welltrain_window, correct_resp,)),
            axis=0, ), :, ]

        for su_idx in range(trial_FR.shape[2]):
            (mm, std) = baselineVector(
                np.squeeze(trial_FR[:, np.logical_and(welltrain_window, correct_resp), su_idx])
            )
            self.su_list.append(
                {
                    # "S1_S1_3": np.squeeze(FR_D3_S1_S1[:, :, su_idx]),
                    # "S2_S1_3": np.squeeze(FR_D3_S2_S1[:, :, su_idx]),
                    # "S1_S2_3": np.squeeze(FR_D3_S1_S2[:, :, su_idx]),
                    # "S2_S2_3": np.squeeze(FR_D3_S2_S2[:, :, su_idx]),
                    # "S1_S1_6": np.squeeze(FR_D6_S1_S1[:, :, su_idx]),
                    # "S2_S1_6": np.squeeze(FR_D6_S2_S1[:, :, su_idx]),
                    # "S1_S2_6": np.squeeze(FR_D6_S1_S2[:, :, su_idx]),
                    # "S2_S2_6": np.squeeze(FR_D6_S2_S2[:, :, su_idx]),
                    "S1_S1_3m": np.squeeze(np.mean((FR_D3_S1_S1[:, :, su_idx] - mm) / std, axis=1)),
                    "S1_S2_3m": np.squeeze(np.mean((FR_D3_S1_S2[:, :, su_idx] - mm) / std, axis=1)),
                    "S2_S1_3m": np.squeeze(np.mean((FR_D3_S2_S1[:, :, su_idx] - mm) / std, axis=1)),
                    "S2_S2_3m": np.squeeze(np.mean((FR_D3_S2_S2[:, :, su_idx] - mm) / std, axis=1)),
                    "S1_S1_6m": np.squeeze(np.mean((FR_D6_S1_S1[:, :, su_idx] - mm) / std, axis=1)),
                    "S1_S2_6m": np.squeeze(np.mean((FR_D6_S1_S2[:, :, su_idx] - mm) / std, axis=1)),
                    "S2_S1_6m": np.squeeze(np.mean((FR_D6_S2_S1[:, :, su_idx] - mm) / std, axis=1)),
                    "S2_S2_6m": np.squeeze(np.mean((FR_D6_S2_S2[:, :, su_idx] - mm) / std, axis=1)),
                }
            )

    def get_features(self):
        # breakpoint()
        return self.su_list


def get_dataset(denovo=False, delay=6):
    features_per_su = []
    reg_list = []
    if denovo:
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

            (_perf_desc, perf_code, welltrain_window, correct_resp,) = align.judgePerformance(trials)

            if perf_code != 3:
                continue

            suid_reg = []
            with open(os.path.join(path, "su_id2reg.csv")) as csvfile:
                l = list(csv.reader(csvfile))[1:]
                suid_reg = [list(i) for i in zip(*l)]

            currStats = ctd_stats()
            currStats.processCTDStats(trial_FR, trials, welltrain_window, correct_resp)
            features_per_su.extend(currStats.get_features())
            reg_list.extend(suid_reg[1])

            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

        ### save to npz file
        np.savez_compressed("ctd_prev.npz", features_per_su=features_per_su, reg_list=reg_list)

    ### load from npz file
    else:
        fstr = np.load("ctd_prev.npz", allow_pickle=True)
        features_per_su = fstr["features_per_su"].tolist()
        reg_list = fstr["reg_list"].tolist()

    return (features_per_su, reg_list)


def plotPCA3D_FR(fr,delay=6):
    pcamat = np.hstack(
        [
            np.squeeze(np.mean(np.vstack([x["S1_S1_3m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S2_S1_3m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S1_S2_3m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S2_S2_3m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S1_S1_6m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S2_S1_6m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S1_S2_6m"] for x in fr]),axis=1)),
            np.squeeze(np.mean(np.vstack([x["S2_S2_6m"] for x in fr]),axis=1)),
        ]
    )
    np.save(f"pcamat_{delay}_sust.npy", pcamat)

    pca = PCA(n_components=20)
    comp = pca.fit_transform(pcamat.T)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax = plt.axes(projection="3d")

    ax.plot3D(
        comp[(68 * 0 + 4): (68 * 0 + 40), 0], comp[(68 * 0 + 4): (68 * 0 + 40), 1],
        comp[(68 * 0 + 4): (68 * 0 + 40), 2], "r-", alpha=0.4
    )  # S1S1
    # S1S1
    ax.plot3D(
        comp[(68 * 1 + 4): (68 * 1 + 40), 0], comp[(68 * 1 + 4): (68 * 1 + 40), 1],
        comp[(68 * 1 + 4): (68 * 1 + 40), 2], "b-", alpha=0.4
    )
    ax.plot3D(
        comp[(68 * 4 + 4): (68 * 4 + 40), 0], comp[(68 * 4 + 4): (68 * 4 + 40), 1],
        comp[(68 * 4 + 4): (68 * 4 + 40), 2], "r-", alpha=0.4
    )  # S1S1
    # S1S1
    ax.plot3D(
        comp[(68 * 5 + 4): (68 * 5 + 40), 0], comp[(68 * 5 + 4): (68 * 5 + 40), 1],
        comp[(68 * 5 + 4): (68 * 5 + 40), 2], "b-", alpha=0.4
    )  # S2S1
    # ax.plot3D(
    #     comp[(68 * 6 + 4): (68 * 6 + 40), 0], comp[(68 * 6 + 4): (68 * 6 + 40), 1],
    #     comp[(68 * 6 + 4): (68 * 6 + 40), 2], "m-", alpha=0.4
    # )
    # ax.plot3D(
    #     comp[(68 * 7 + 4): (68 * 7 + 40), 0], comp[(68 * 7 + 4): (68 * 7 + 40), 1],
    #     comp[(68 * 7 + 4): (68 * 7 + 40), 2], "c-", alpha=0.4
    # )
    for bidx in [0, 1, 4, 5]:
        ax.plot3D(
            comp[(68 * bidx + 4): (68 * bidx + 12), 0], comp[(68 * bidx + 4): (68 * bidx + 12), 1],
            comp[(68 * bidx + 4): (68 * bidx + 12), 2], "k.", markersize=1, alpha=0.2
        )  # S1S1
    # S1S1
    # ax.plot3D(
    #     comp[(68 * 5 + 4): (68 * 5 + 12), 0], comp[(68 * 5 + 4): (68 * 5 + 12), 1],
    #     comp[(68 * 5 + 4): (68 * 5 + 12), 2], "k.", markersize=1, alpha=0.2
    # )  # S2S1
    # ax.plot3D(
    #     comp[(68 * 6 + 4): (68 * 6 + 12), 0], comp[(68 * 6 + 4): (68 * 6 + 12), 1],
    #     comp[(68 * 6 + 4): (68 * 6 + 12), 2], "k.", markersize=1, alpha=0.2
    # )
    # ax.plot3D(
    #     comp[(68 * 7 + 4): (68 * 7 + 12), 0], comp[(68 * 7 + 4): (68 * 7 + 12), 1],
    #     comp[(68 * 7 + 4): (68 * 7 + 12), 2], "k.", markersize=1, alpha=0.2
    # )

    plt.show()
    return pcamat


if __name__ == "__main__":
    delay = 6
    (features_per_su, reg_list) = get_dataset(False)
    trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # list(trans6.keys())
    sust = trans_fstr["sust"]
    trans = trans_fstr["transient"]
    # reg_arr = trans_fstr['reg_arr']

    fr_trans = [f for (f, t) in zip(features_per_su, trans) if t]
    fr_sust = [f for (f, t) in zip(features_per_su, sust) if t]
    plotPCA3D_FR(fr_sust)
    plotPCA3D_FR(fr_trans)
