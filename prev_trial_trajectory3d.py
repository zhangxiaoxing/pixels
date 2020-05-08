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
from scipy.spatial.distance import euclidean
from sklearn.model_selection import KFold
import scipy.stats as stats


class PCA_stats:
    def __init__(self):
        self.su_list = []
        self.avail = []
        self.kf = KFold(2, shuffle=True)

    def splitMean(self, FR, su_idx, mm, std, random_type=None):
        if random_type is None:
            out = np.squeeze(np.mean((FR[:, :, su_idx] - mm) / std, axis=1))
        elif random_type == 'half':
            sel = list(self.kf.split(np.arange(FR.shape[1])))[0]
            out = np.squeeze(np.mean((FR[:, sel[0], su_idx] - mm) / std, axis=1))
        else:
            sel = np.random.choice(FR.shape[1])
            out = np.squeeze((FR[:, sel, su_idx] - mm) / std)
        return out

        ### unprocessed data entry point

    def processCTDStats(self, trial_FR, trials, welltrain_window=None, correct_resp=None, random_type=None):
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
                        "S1_S1_3m": self.splitMean(FR_D3_S1_S1, su_idx, mm, std, random_type),
                        "S1_S2_3m": self.splitMean(FR_D3_S1_S2, su_idx, mm, std, random_type),
                        "S2_S1_3m": self.splitMean(FR_D3_S2_S1, su_idx, mm, std, random_type),
                        "S2_S2_3m": self.splitMean(FR_D3_S2_S2, su_idx, mm, std, random_type),
                        "S1_S1_6m": self.splitMean(FR_D6_S1_S1, su_idx, mm, std, random_type),
                        "S1_S2_6m": self.splitMean(FR_D6_S1_S2, su_idx, mm, std, random_type),
                        "S2_S1_6m": self.splitMean(FR_D6_S2_S1, su_idx, mm, std, random_type),
                        "S2_S2_6m": self.splitMean(FR_D6_S2_S2, su_idx, mm, std, random_type),
                    }
                )
                self.avail.append(True)
            else:
                self.avail.append(False)

    def get_features(self):
        # breakpoint()
        return (self.su_list, self.avail)


def get_dataset(denovo=False, random_type=None):
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

            currStats = PCA_stats()
            currStats.processCTDStats(trial_FR, trials, welltrain_window, correct_resp, random_type)
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

    return (features_per_su, reg_list, avails)


def plotPCA3D_FR(fr, plot_delay=6, to_plot=True, fig=None, ax=None, the_ref=True, coeff=None):
    if plot_delay == 'delay_diff':
        pcamat = np.hstack(
            [
                np.vstack([np.mean(np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
    else:
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
    if the_ref:
        pca = PCA(n_components=20)
        comp = pca.fit_transform(pcamat.T)
        coeff = pca.components_
        ratio = pca.explained_variance_ratio_
    else:
        pcamat_cent = pcamat - np.expand_dims(np.mean(pcamat, axis=1), axis=1)
        comp = np.matmul(coeff, pcamat_cent).T
    if the_ref:
        alpha = 1
        lw = 1.5
    else:
        alpha = 0.25
        lw = 0.5

    if to_plot:
        if plot_delay == 'delay_diff':

            for i in range(21, 40):
                arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                              [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                              lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", alpha=alpha)
                ax.add_artist(arr)
            for i in range(73, 104):
                arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                              [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                              lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="c", alpha=alpha)
                ax.add_artist(arr)
            if the_ref:
                for bidx in [0, 20]:
                    ax.plot3D(
                        comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
                        comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
                    )  # S1S1

                    ax.plot3D(
                        comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
                        comp[bidx + 4:bidx + 8, 2], "k+", markersize=4, alpha=0.25
                    )

        else:
            if plot_delay == 3:
                for i in range(1, 20):
                    arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                                  [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                                  lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="r", alpha=alpha)
                    ax.add_artist(arr)
                    j = i + 20
                    arr = Arrow3D([comp[j - 1, 0], comp[j, 0]], [comp[j - 1, 1], comp[j, 1]],
                                  [comp[j - 1, 2], comp[j, 2]], mutation_scale=20,
                                  lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", alpha=alpha)
                    ax.add_artist(arr)
                if the_ref:
                    for bidx in [0, 20]:
                        ax.plot3D(
                            comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
                            comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
                        )  # S1S1

                        ax.plot3D(
                            comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
                            comp[bidx + 4:bidx + 8, 2], "k+", markersize=4, alpha=0.25
                        )
            elif plot_delay == 6:
                for i in range(41, 72):
                    arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                                  [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                                  lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="m", alpha=alpha)
                    ax.add_artist(arr)
                    j = i + 32
                    arr = Arrow3D([comp[j - 1, 0], comp[j, 0]], [comp[j - 1, 1], comp[j, 1]],
                                  [comp[j - 1, 2], comp[j, 2]], mutation_scale=20,
                                  lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="c", alpha=alpha)
                    ax.add_artist(arr)
                if the_ref:
                    for bidx in [40, 72]:
                        ax.plot3D(
                            comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
                            comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
                        )  # S1S1

                        ax.plot3D(
                            comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
                            comp[bidx + 4:bidx + 8, 2], "k+", markersize=4, alpha=0.25
                        )
    return (comp, coeff)


def plot3d_all(plot_delay=3, repeats=5, to_plot=True):
    fig_all = plt.figure()
    ax_all = fig_all.add_subplot(111, projection='3d')

    fig_sust = plt.figure()
    ax_sust = fig_sust.add_subplot(111, projection='3d')

    fig_trans = plt.figure()
    ax_trans = fig_trans.add_subplot(111, projection='3d')

    (features_per_su, reg_list, avails) = get_dataset(denovo=True, random_type=None)
    if plot_delay == 'delay_diff':
        trans_fstr3 = np.load(f"sus_trans_pie_3.npz")
        sust3 = trans_fstr3["sust"]
        trans3 = trans_fstr3["transient"]

        trans_fstr6 = np.load(f"sus_trans_pie_6.npz")
        sust6 = trans_fstr6["sust"]
        trans6 = trans_fstr6["transient"]

        sust = np.logical_or(sust3, sust6)
        trans = np.logical_or(trans3, trans6)
    else:
        trans_fstr = np.load(f"sus_trans_pie_{plot_delay}.npz")
        sust = trans_fstr["sust"]
        trans = trans_fstr["transient"]

    trans_avail = trans[avails]
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
    sust_avail = sust[avails]
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
    # reg_arr = trans_fstr['reg_arr']

    (comp_all, coeff_all) = plotPCA3D_FR(features_per_su, plot_delay, fig=fig_all, ax=ax_all, the_ref=True)
    (comp_sust, coeff_sust) = plotPCA3D_FR(fr_sust, plot_delay, fig=fig_sust, ax=ax_sust, the_ref=True)
    (comp_trans, coeff_trans) = plotPCA3D_FR(fr_trans, plot_delay, fig=fig_trans, ax=ax_trans, the_ref=True)

    for r in range(repeats):
        (features_per_su, reg_list, _avails) = get_dataset(denovo=True, random_type='one_trial')
        if plot_delay == 'delay_diff':
            trans_fstr3 = np.load(f"sus_trans_pie_3.npz")
            sust3 = trans_fstr3["sust"]
            trans3 = trans_fstr3["transient"]

            trans_fstr6 = np.load(f"sus_trans_pie_6.npz")
            sust6 = trans_fstr6["sust"]
            trans6 = trans_fstr6["transient"]

            sust = np.logical_or(sust3, sust6)
            trans = np.logical_or(trans3, trans6)
        else:
            trans_fstr = np.load(f"sus_trans_pie_{plot_delay}.npz")
            sust = trans_fstr["sust"]
            trans = trans_fstr["transient"]

        trans_avail = trans[avails]
        fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
        sust_avail = sust[avails]
        fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]
        # reg_arr = trans_fstr['reg_arr']

        (_comp, _c) = plotPCA3D_FR(features_per_su, plot_delay, fig=fig_all, ax=ax_all, the_ref=False, coeff=coeff_all)
        (_comp, _c) = plotPCA3D_FR(fr_sust, plot_delay, fig=fig_sust, ax=ax_sust, the_ref=False, coeff=coeff_sust)
        (_comp, _c) = plotPCA3D_FR(fr_trans, plot_delay, fig=fig_trans, ax=ax_trans, the_ref=False, coeff=coeff_trans)

    np.savez(f'PCA_comp_{plot_delay}.npz', comp_all=comp_all, comp_sust=comp_sust, comp_trans=comp_trans)

    if plot_delay == 'delay_diff':

        ax_all.set_xlim(-10, 200)
        ax_all.set_ylim(-60, 60)
        ax_all.set_zlim(-40, 40)

        ax_sust.set_xlim(-15, 45)
        ax_sust.set_zlim(-12, 10)
        ax_sust.set_ylim(-15, 23)

        ax_trans.set_xlim(-35, 120)
        ax_trans.set_ylim(-45, 50)
        ax_trans.set_zlim(-30, 30)

    elif plot_delay == 3:
        ax_all.set_ylim(-60, 65)
        ax_sust.set_zlim(-12, 20)
        ax_sust.set_ylim(-20, 30)
        ax_trans.set_ylim(-35, 40)
        ax_trans.set_zlim(-24, 25)

    elif plot_delay == 6:
        ax_all.set_xlim(-70, 200)
        ax_all.set_ylim(-70, 80)
        ax_all.set_zlim(-50, 50)

        ax_sust.set_xlim(-15, 35)
        ax_sust.set_zlim(-12, 15)
        ax_sust.set_ylim(-10, 23)

        ax_trans.set_xlim(-35, 90)
        ax_trans.set_ylim(-45, 50)
        ax_trans.set_zlim(-30, 30)

    for oneax in [ax_all, ax_sust, ax_trans]:
        oneax.set_xlabel('PC1')
        oneax.set_ylabel('PC2')
        oneax.set_zlabel('PC3')

    for onefig in [fig_all, fig_sust, fig_trans]:
        onefig.set_size_inches(65 / 25.4, 65 / 25.4)
    if to_plot:
        fig_all.savefig(f'traj_all_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_sust.savefig(f'traj_sust_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        fig_trans.savefig(f'traj_trans_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        plt.show()
    return (fig_all, ax_all, fig_sust, ax_sust, fig_trans, ax_trans)


class Dist_stats:
    def __init__(self, ref_comp):
        self.ref_comp = ref_comp
        self.S1_S2_3_dist = []
        self.S1_S2_6_dist = []
        self.exit_dist_3_S1 = []
        self.exit_dist_3_S2 = []
        self.exit_dist_6_S1 = []
        self.exit_dist_6_S2 = []

        self.latedelay_dist_3_S1 = []
        self.latedelay_dist_3_S2 = []
        self.latedelay_dist_6_S1 = []
        self.latedelay_dist_6_S2 = []

        self.ref_dist_3_S1 = []
        self.ref_dist_3_S2 = []
        self.ref_dist_6_S1 = []
        self.ref_dist_6_S2 = []

        self.coeff = []

    def append_data(self, fr, coeff):
        pcamat = np.hstack(
            [
                np.vstack([np.mean(np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
        # pca = PCA(n_components=20)
        # comp = pca.fit_transform(pcamat.T)
        # coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        pcamat_cent = pcamat - np.expand_dims(np.mean(pcamat, axis=1), axis=1)
        comp = np.matmul(coeff, pcamat_cent).T

        S1_S2_3 = []
        deci_vec_S1_3 = []
        deci_vec_S2_3 = []
        late_vec_S1_3 = []
        late_vec_S2_3 = []
        ref_S1_3 = []
        ref_S2_3 = []
        for bin in range(20):  # time range
            S1_S2_3.append(euclidean(comp[bin, :3], comp[bin + 20, :3]))  # comp range
            deci_vec_S1_3.append(euclidean(comp[bin, :3], self.ref_comp[19, :3]))
            deci_vec_S2_3.append(euclidean(comp[bin + 20, :3], self.ref_comp[39, :3]))
            late_vec_S1_3.append(euclidean(comp[bin, :3], self.ref_comp[17, :3]))
            late_vec_S2_3.append(euclidean(comp[bin + 20, :3], self.ref_comp[37, :3]))
            ref_S1_3.append(euclidean(comp[bin, :3], self.ref_comp[bin, :3]))
            ref_S2_3.append(euclidean(comp[bin + 20, :3], self.ref_comp[bin + 20, :3]))
        self.S1_S2_3_dist.append(S1_S2_3)
        self.exit_dist_3_S1.append(deci_vec_S1_3)
        self.exit_dist_3_S2.append(deci_vec_S2_3)
        self.latedelay_dist_3_S1.append(late_vec_S1_3)
        self.latedelay_dist_3_S2.append(late_vec_S2_3)
        self.ref_dist_3_S1.append(ref_S1_3)
        self.ref_dist_3_S2.append(ref_S2_3)
        S1_S2_6 = []
        deci_vec_S1_6 = []
        deci_vec_S2_6 = []
        late_vec_S1_6 = []
        late_vec_S2_6 = []
        ref_S1_6 = []
        ref_S2_6 = []
        for bin in range(40, 72):  # time range
            S1_S2_6.append(euclidean(comp[bin, :3], comp[bin + 32, :3]))  # comp range
            deci_vec_S1_6.append(euclidean(comp[bin, :3], self.ref_comp[19, :3]))
            deci_vec_S2_6.append(euclidean(comp[bin + 32, :3], self.ref_comp[39, :3]))
            late_vec_S1_6.append(euclidean(comp[bin, :3], self.ref_comp[17, :3]))
            late_vec_S2_6.append(euclidean(comp[bin + 32, :3], self.ref_comp[37, :3]))
            ref_S1_6.append(euclidean(comp[bin, :3], self.ref_comp[bin, :3]))
            ref_S2_6.append(euclidean(comp[bin + 32, :3], self.ref_comp[bin + 32, :3]))

        self.S1_S2_6_dist.append(S1_S2_6)
        self.exit_dist_6_S1.append(deci_vec_S1_6)
        self.exit_dist_6_S2.append(deci_vec_S2_6)
        self.latedelay_dist_6_S1.append(late_vec_S1_6)
        self.latedelay_dist_6_S2.append(late_vec_S2_6)
        self.ref_dist_6_S1.append(ref_S1_6)
        self.ref_dist_6_S2.append(ref_S2_6)

    def get_coeff(self, fr):
        pcamat = np.hstack(
            [
                np.vstack([np.mean(np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
        pca = PCA(n_components=20)
        comp = pca.fit_transform(pcamat.T)
        coeff = pca.components_
        # ratio = pca.explained_variance_ratio_

        return coeff

    def get_data(self):
        return (self.S1_S2_3_dist, self.S1_S2_6_dist, self.exit_dist_3_S1,
                self.exit_dist_3_S2, self.exit_dist_6_S1, self.exit_dist_6_S2,
                self.ref_dist_3_S1, self.ref_dist_3_S2, self.ref_dist_6_S1, self.ref_dist_6_S2,
                self.latedelay_dist_3_S1, self.latedelay_dist_3_S2, self.latedelay_dist_6_S1, self.latedelay_dist_6_S2)


def dist_all(repeats=2):
    fstr = np.load(f'PCA_comp_delay_diff.npz')
    comp_all = fstr['comp_all']
    comp_sust = fstr['comp_sust']
    comp_trans = fstr['comp_trans']

    dist_all = Dist_stats(comp_all)
    dist_sust = Dist_stats(comp_sust)
    dist_trans = Dist_stats(comp_trans)

    trans_fstr3 = np.load("sus_trans_pie_3.npz")
    sust3 = trans_fstr3["sust"]
    trans3 = trans_fstr3["transient"]

    trans_fstr6 = np.load("sus_trans_pie_6.npz")
    sust6 = trans_fstr6["sust"]
    trans6 = trans_fstr6["transient"]

    sust = np.logical_or(sust3, sust6)
    trans = np.logical_or(trans3, trans6)

    (features_per_su, reg_list, avails) = get_dataset(denovo=True, random_type=None)
    trans_avail = trans[avails]
    fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
    sust_avail = sust[avails]
    fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]

    coeff_all = dist_all.get_coeff(features_per_su)
    coeff_sust = dist_sust.get_coeff(fr_sust)
    coeff_trans = dist_trans.get_coeff(fr_trans)

    for r in range(repeats):
        (features_per_su, reg_list, _avails) = get_dataset(denovo=True, random_type='one_trial')
        fr_trans = [f for (f, t) in zip(features_per_su, trans_avail) if t]
        fr_sust = [f for (f, t) in zip(features_per_su, sust_avail) if t]

        dist_all.append_data(features_per_su, coeff_all)
        dist_sust.append_data(fr_sust, coeff_sust)
        dist_trans.append_data(fr_trans, coeff_trans)

    dist_all_list = dist_all.get_data()
    dist_sust_list = dist_sust.get_data()
    dist_trans_list = dist_trans.get_data()

    np.savez(f'pca_dist_{repeats}.npz', dist_all=dist_all_list, dist_sust=dist_sust_list, dist_trans=dist_trans_list)

    return (dist_all_list, dist_sust_list, dist_trans_list)


def plot_dist_one(d_one, subgroup, repeats):
    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[0]), axis=0)
    sem = stats.sem(np.vstack(d_one[0]))
    plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="r", alpha=0.2, ls='None')
    ax.plot(mm, '-r', lw=1)

    mm = np.mean(np.vstack(d_one[1]), axis=0)
    sem = stats.sem(np.vstack(d_one[1]))
    plt.fill_between(np.arange(mm.shape[0]), mm - sem, mm + sem, color="b", alpha=0.2, ls='None')
    ax.plot(mm, '-b', lw=1)
    ax.set_ylabel('S1-S2 PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'Sample_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[2]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[2]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[4]), axis=0)
    ax.plot(np.arange(32), np.vstack(d_one[4]).T, '-b', alpha=0.1, lw=0.25)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.set_ylabel('S1-exit PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S1_exit_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[3]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[3]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[5]), axis=0)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.plot(np.arange(32), np.vstack(d_one[5]).T, '-b', alpha=0.1, lw=0.25)
    ax.set_ylabel('S2-exit PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S2_exit_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[6]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[6]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[8]), axis=0)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.plot(np.arange(32), np.vstack(d_one[8]).T, '-b', alpha=0.1, lw=0.25)
    ax.set_ylabel('S1-all trial mean distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S1_alltrial_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[7]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[7]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[9]), axis=0)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.plot(np.arange(32), np.vstack(d_one[9]).T, '-b', alpha=0.1, lw=0.25)
    ax.set_ylabel('S2-all trial mean PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S2_alltrial_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')
    plt.show()


    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[10]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[10]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[12]), axis=0)
    ax.plot(np.arange(32), np.vstack(d_one[12]).T, '-b', alpha=0.1, lw=0.25)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.set_ylabel('S1-late delay PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S1_latedelay_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[11]), axis=0)
    ax.plot(np.arange(20), np.vstack(d_one[11]).T, '-r', alpha=0.1, lw=0.25)
    ax.plot(mm, '-r', lw=1, alpha=0.5)

    mm = np.mean(np.vstack(d_one[13]), axis=0)
    ax.plot(mm, '-b', lw=1, alpha=0.5)
    ax.plot(np.arange(32), np.vstack(d_one[13]).T, '-b', alpha=0.1, lw=0.25)
    ax.set_ylabel('S2-late delay PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(f'S2_latedelay_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')




def plot_dist(repeats):
    fstr = np.load(f'pca_dist_{repeats}.npz', allow_pickle=True)
    dist_all = fstr['dist_all']
    dist_sust = fstr['dist_sust']
    dist_trans = fstr['dist_trans']

    plot_dist_one(dist_all, 'all', repeats)
    plot_dist_one(dist_sust, 'sust', repeats)
    plot_dist_one(dist_trans, 'trans', repeats)


if __name__ == "__main__":
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    # plot_delay = 'delay_diff'

    # # plot_dist()
    # (fig_all, ax_all, fig_sust, ax_sust, fig_trans, ax_trans) = plot3d_all(plot_delay=plot_delay, repeats=0)
    # (dist_all_list, dist_sust_list, dist_trans_list) = dist_all(repeats=100)
    plot_dist(repeats=100)
    # delay = 3
    # (features_per_su, reg_list) = get_dataset(True, rnd_half=True)
    # trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # # list(trans6.keys())
    # sust = trans_fstr["sust"]
    # trans = trans_fstr["transient"]
