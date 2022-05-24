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

from sklearn.model_selection import KFold
import scipy.stats as stats


def plotPCA3D_FR(fr, plot_delay=6, to_plot=True, fig=None, ax=None, the_ref=True, coeff=None):
    if plot_delay == 'delay_diff':
        pcamat = np.hstack(
            [
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
            ]
        )
    else:
        pcamat = np.hstack(
            [
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_3m"][8:28], x["S1_S2_3m"][8:28])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_3m"][8:28], x["S2_S2_3m"][8:28])), axis=0) for x in fr]),
                # np.vstack([x["S1_S2_3m"][8:28] for x in fr]),
                # np.vstack([x["S2_S2_3m"][8:28] for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S1_S1_6m"][8:40], x["S1_S2_6m"][8:40])), axis=0) for x in fr]),
                np.vstack([np.mean(
                    np.vstack((x["S2_S1_6m"][8:40], x["S2_S2_6m"][8:40])), axis=0) for x in fr]),
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

            # for i in range(36, 40):
            for i in range(16, 20):
                arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                              [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                              lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="b", alpha=alpha)
                ax.add_artist(arr)
            # for i in range(88, 104):
            for i in range(56, 72):
                arr = Arrow3D([comp[i - 1, 0], comp[i, 0]], [comp[i - 1, 1], comp[i, 1]],
                              [comp[i - 1, 2], comp[i, 2]], mutation_scale=20,
                              lw=lw, arrowstyle="->,head_length=0.1, head_width=0.05", color="c", alpha=alpha)
                ax.add_artist(arr)
            # if the_ref:
                # for bidx in [35, 40]:
                # ax.plot3D(
                #     comp[bidx:bidx + 4, 0], comp[bidx:bidx + 4, 1],
                #     comp[bidx:bidx + 4, 2], "k.", markersize=1, alpha=0.25
                # )  # S1S1

                # ax.plot3D(
                #     comp[bidx + 4:bidx + 8, 0], comp[bidx + 4:bidx + 8, 1],
                #     comp[bidx + 4:bidx + 8, 2], "k+", markersize=4, alpha=0.25
                # )

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
                            comp[bidx + 4:bidx + 8,
                                 0], comp[bidx + 4:bidx + 8, 1],
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
                            comp[bidx + 4:bidx + 8,
                                 0], comp[bidx + 4:bidx + 8, 1],
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
    (features_per_su, reg_list, avails) = get_dataset(
        denovo=True, random_type=None)
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

    (comp_all, coeff_all) = plotPCA3D_FR(features_per_su,
                                         plot_delay, fig=fig_all, ax=ax_all, the_ref=True)
    (comp_sust, coeff_sust) = plotPCA3D_FR(
        fr_sust, plot_delay, fig=fig_sust, ax=ax_sust, the_ref=True)
    (comp_trans, coeff_trans) = plotPCA3D_FR(fr_trans,
                                             plot_delay, fig=fig_trans, ax=ax_trans, the_ref=True)
    for r in range(repeats):
        (features_per_su, reg_list, _avails) = get_dataset(
            denovo=True, random_type='one_trial')
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
        (_comp, _c) = plotPCA3D_FR(features_per_su, plot_delay,
                                   fig=fig_all, ax=ax_all, the_ref=False, coeff=coeff_all)
        (_comp, _c) = plotPCA3D_FR(fr_sust, plot_delay,
                                   fig=fig_sust, ax=ax_sust, the_ref=False, coeff=coeff_sust)
        (_comp, _c) = plotPCA3D_FR(fr_trans, plot_delay,
                                   fig=fig_trans, ax=ax_trans, the_ref=False, coeff=coeff_trans)
    np.savez(f'PCA_comp_{plot_delay}.npz', comp_all=comp_all,
             comp_sust=comp_sust, comp_trans=comp_trans)
    if plot_delay == 'delay_diff':
        ax_all.set_xlim(-40, -20)
        ax_all.set_ylim(-50, 10)
        ax_all.set_zlim(5, 35)
        ax_sust.set_xlim(-12.5, 12.5)
        ax_sust.set_ylim(-20, 0)
        ax_sust.set_zlim(-10, 7)

        ax_trans.set_xlim(-35, -15)
        ax_trans.set_ylim(-45, 5)
        ax_trans.set_zlim(-30, 0)

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
        fig_all.savefig(f'traj_all_{plot_delay}.pdf',
                        dpi=300, bbox_inches='tight')
        fig_sust.savefig(f'traj_sust_{plot_delay}.pdf',
                         dpi=300, bbox_inches='tight')
        fig_trans.savefig(
            f'traj_trans_{plot_delay}.pdf', dpi=300, bbox_inches='tight')
        plt.show()
    return (fig_all, ax_all, fig_sust, ax_sust, fig_trans, ax_trans)



def plot_dist_one(d_one, subgroup, repeats):
    # plot time-distance 2d figure

    # S1-S2 distance

    (fig, ax) = plt.subplots(1, 1, figsize=(90 / 25.4, 90 / 25.4), dpi=300)
    mm = np.mean(np.vstack(d_one[0]), axis=0)
    sem = stats.sem(np.vstack(d_one[0]))
    plt.fill_between(
        np.arange(mm.shape[0]), mm - sem, mm + sem, color="r", alpha=0.2, ls='None')
    ax.plot(mm, '-r', lw=1)

    mm = np.mean(np.vstack(d_one[1]), axis=0)
    sem = stats.sem(np.vstack(d_one[1]))
    plt.fill_between(
        np.arange(mm.shape[0]), mm - sem, mm + sem, color="b", alpha=0.2, ls='None')
    ax.plot(mm, '-b', lw=1)
    ax.set_ylabel('S1-S2 PC distance')
    ax.set_xticks([3.5, 23.5])
    ax.set_xticklabels([0, 5])
    ax.set_xlabel('time (s)')
    [ax.axvline(x, color='k', ls=':') for x in [3.5, 7.5, 19.5]]
    fig.savefig(
        f'Sample_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    # S1-Exit
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
    fig.savefig(
        f'S1_exit_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    # S2-Exit
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
    fig.savefig(
        f'S2_exit_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    # S1-all trial mean distance
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
    fig.savefig(
        f'S1_alltrial_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    # S2-all trial mean distance
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
    fig.savefig(
        f'S2_alltrial_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')
    plt.show()

    # S1-late delay distance
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
    ax.set_yscale('log')
    fig.savefig(
        f'S1_latedelay_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')

    # S2-late delay distance
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
    ax.set_yscale('log')
    fig.savefig(
        f'S2_latedelay_PC_Dist_{subgroup}_{repeats}.pdf', bbox_inches='tight')


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
    # (fig_all, ax_all, fig_sust, ax_sust, fig_trans, ax_trans) = plot3d_all(plot_delay=plot_delay, repeats=100)

    # (dist_all_list, dist_sust_list, dist_trans_list) = dist_all(repeats=100)
    plot_dist(repeats=100)
    # delay = 3
    # (features_per_su, reg_list) = get_dataset(True, rnd_half=True)
    # trans_fstr = np.load(f"sus_trans_pie_{delay}.npz")
    # # list(trans6.keys())
    # sust = trans_fstr["sust"]
    # trans = trans_fstr["transient"]