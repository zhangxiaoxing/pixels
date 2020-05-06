# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:49:17 2020

@author: Libra
"""
import os
import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt
import su_region_align as align
from sklearn.decomposition import PCA
from pixelStats import baselineVector
from Arrow3D import Arrow3D
from matplotlib import rcParams


class ctd_stats:
    def __init__(self):
        self.su_list = []
        self.avail = []

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
            if std > 0 and np.all([x.shape[1] >= 2 for x in (
                    FR_D3_S1_S1,
                    FR_D3_S1_S2,
                    FR_D3_S2_S1,
                    FR_D3_S2_S2,
                    FR_D6_S1_S1,
                    FR_D6_S1_S2,
                    FR_D6_S2_S1,
                    FR_D6_S2_S2,
            )]):
                self.su_list.append(
                    {
                        "S1_S1_3m": np.squeeze(np.mean((FR_D3_S1_S1[:, :, su_idx] - mm) / std, axis=1)),
                        "S1_S2_3m": np.squeeze(np.mean((FR_D3_S1_S2[:, :, su_idx] - mm) / std, axis=1)),
                        "S2_S1_3m": np.squeeze(np.mean((FR_D3_S2_S1[:, :, su_idx] - mm) / std, axis=1)),
                        "S2_S2_3m": np.squeeze(np.mean((FR_D3_S2_S2[:, :, su_idx] - mm) / std, axis=1)),
                        "S1_S1_6m": np.squeeze(np.mean((FR_D6_S1_S1[:, :, su_idx] - mm) / std, axis=1)),
                        "S1_S2_6m": np.squeeze(np.mean((FR_D6_S1_S2[:, :, su_idx] - mm) / std, axis=1)),
                        "S2_S1_6m": np.squeeze(np.mean((FR_D6_S2_S1[:, :, su_idx] - mm) / std, axis=1)),
                        "S2_S2_6m": np.squeeze(np.mean((FR_D6_S2_S2[:, :, su_idx] - mm) / std, axis=1)),
                    })
                self.avail.append(True)
            else:
                self.avail.append(False)

    def get_features(self):
        # breakpoint()
        return (self.su_list, self.avail)


def get_dataset(denovo=False, delay=6):
    features_per_su = []
    avails = []
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
            (feat, avail) = currStats.get_features()
            features_per_su.extend(feat)
            avails.extend(avail)
            reg_list.extend(suid_reg[1])

            ### DEBUG in small subsets
            # if len(features_per_su)>50:
            #     break

        ### save to npz file
        np.savez_compressed("ctd_prev_stats.npz", features_per_su=features_per_su, reg_list=reg_list)

    ### load from npz file
    else:
        pass
        # fstr = np.load("ctd_prev_stats.npz", allow_pickle=True)
        # features_per_su = fstr["features_per_su"].tolist()
        # reg_list = fstr["reg_list"].tolist()

    return (features_per_su, reg_list,avails)


def plotPCA3D_FR(fr, delay=6):
    pcamat = np.hstack(
        [
            np.vstack([np.mean(np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
            np.vstack([np.mean(np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
            # np.vstack([x["S1_S2_3m"][8:28] for x in fr]),
            # np.vstack([x["S2_S2_3m"][8:28] for x in fr]),
            np.vstack([np.mean(np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
            np.vstack([np.mean(np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            # np.vstack([x["S1_S2_6m"][8:40] for x in fr]),
            # np.vstack([x["S2_S2_6m"][8:40] for x in fr]),
        ]
    )
    # np.save(f"pcamat_{delay}_sust.npy", pcamat)

    pca = PCA(n_components=20)
    comp = pca.fit_transform(pcamat.T)
    coeff = pca.components_
    ratio = pca.explained_variance_ratio_

    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax = plt.axes(projection="3d")

    for i in range(1, 20):
        arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                      [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                      lw=1, arrowstyle="->,head_length=0.1, head_width=0.05", color="r", alpha=0.5)
        ax.add_artist(arr)
        j = i + 20
        arr = Arrow3D([comp[j - 1, 0], comp[j, 0]], [comp[j - 1, 1], comp[j, 1]],
                      [comp[j - 1, 2], comp[j, 2]], mutation_scale=20,
                      lw=1, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", alpha=0.5)
        ax.add_artist(arr)

    for i in range(41, 72):
        # arr = Arrow3D([comp[i-1,0],comp[i,0]], [comp[i-1,1], comp[i,1]],
        #             [comp[i-1,2], comp[i,2]], mutation_scale=20,
        #             lw=1, arrowstyle="->,head_length=0.1, head_width=0.05", color="m",alpha=0.5)
        # ax.add_artist(arr)
        j = i + 32
        arr = Arrow3D([comp[j - 1, 0], comp[j, 0]], [comp[j - 1, 1], comp[j, 1]],
                      [comp[j - 1, 2], comp[j, 2]], mutation_scale=20,
                      lw=1, arrowstyle="->,head_length=0.1, head_width=0.05", color="c", alpha=0.5)
        ax.add_artist(arr)

    for bidx in [0, 20, 72]:
        ax.plot3D(
            comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
            comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
        )  # S1S1

        ax.plot3D(
            comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
            comp[bidx + 4:bidx + 8, 2], "k*", markersize=4, alpha=0.25
        )

    for bidx in [72]:
        ax.plot3D(
            comp[bidx + 20:bidx + 32, 0], comp[bidx + 20:bidx + 32, 1],
            comp[bidx + 20:bidx + 32, 2], "k+", markersize=3, alpha=0.25
        )

    plt.show()
    return (fig, ax)


def plot3d_all():
    delay = 3
    (features_per_su, reg_list, avails) = get_dataset(True)
    trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    sust = trans_fstr["sust"]
    trans = trans_fstr["transient"]

    trans_avail = trans[avails]
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
    sust_avail = sust[avails]
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
    # reg_arr = trans_fstr['reg_arr']

    fr_trans = [f for (f, t) in zip(features_per_su, trans) if t]
    fr_sust = [f for (f, t) in zip(features_per_su, sust) if t]

    (fig_all, ax_all) = plotPCA3D_FR(features_per_su)
    (fig_sust, ax_sust) = plotPCA3D_FR(fr_sust)
    (fig_trans, ax_trans) = plotPCA3D_FR(fr_trans)

    ax_all.set_ylim(-60, 65)
    ax_sust.set_zlim(-15, 20)
    ax_trans.set_ylim(-50, 50)
    ax_trans.set_zlim(-30, 30)

    for oneax in [ax_all, ax_sust, ax_trans]:
        oneax.set_xlabel('PC1')
        oneax.set_ylabel('PC2')
        oneax.set_zlabel('PC3')

    for onefig in [fig_all, fig_sust, fig_trans]:
        onefig.set_size_inches(65 / 25.4, 65 / 25.4)

    fig_all.savefig('traj_all.pdf', dpi=300, bbox_inches='tight')
    fig_sust.savefig('traj_sust.pdf', dpi=300, bbox_inches='tight')
    fig_trans.savefig('traj_trans.pdf', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']

    plot3d_all()
    # delay = 3
    # (features_per_su, reg_list) = get_dataset(True, rnd_half=True)
    # trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # # list(trans6.keys())
    # sust = trans_fstr["sust"]
    # trans = trans_fstr["transient"]
